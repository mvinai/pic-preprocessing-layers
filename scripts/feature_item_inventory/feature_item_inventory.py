#!/usr/bin/env python3
"""
Feature item inventory, Stage 1 of the feature sparsity and measurement-frequency audit (read-only, aggregate-only).
Implements `feature_sparsity_frequency_audit_spec.md` v0.2, section 5.1 (decisions of Nov 14, 2026).

Per raw table, per item (ITEMID; drug name pair; specimen; organism), TRAIN split only, cohort encounters only:
    rows in stay, rows with a numeric value, encounters with >= 1 row in stay / in the Layer B window /
    in the pre-ICU window [t0 - 6 h, t0), the unit strings seen, and the dictionary label.
Linkage (contract v12 H3): SUBJECT_ID and HADM_ID in the encounter's hadm_id_list, then the event time in [t0, OUT).
Event times: LABEVENTS release time = CHARTTIME + 1 h (decision 11 / A2.3); PRESCRIPTIONS STARTDATE; others CHARTTIME.
It also writes blank map templates (channel ITEMID map, drug-group map incl. NMB and AED, specimen map, organism map)
listing only items whose train encounter count is >= 10.

Reuse: ID parsing, the disclosure registry, complementary suppression and the recoverability solver come from the
duration audit v1.7 module (`audit_surgery_complexity_duration.py`), SHA-256-pinned and verified BEFORE import.

Usage (Windows, from the project root; this script and the v1.7 module in scripts/):
    python scripts/feature_item_inventory.py --raw-dir "C:/.../PICDataset" --split-file "C:/.../subject_splits.csv"

Outputs (data/interim/):
    feature_item_inventory.json / .md                    shared (suppressed); written ONLY if the release gate passes
    feature_channel_map_TEMPLATE.csv, drug_group_map_TEMPLATE.csv,
    specimen_map_TEMPLATE.csv, organism_map_TEMPLATE.csv shared templates (gated with the JSON)
    feature_item_inventory_raw_LOCAL.json                local only, never share (unsuppressed counts)
Exit codes: 0 ok; 2 input/hash/reconciliation failure; 3 release gate failure (no shared file left in place).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

# =============================================================================
# CONFIG BLOCK (spec v0.2, Nov 14, 2026)
# =============================================================================
SPEC_VERSION = "feature_item_inventory v1.0 (feature_sparsity_frequency_audit_spec v0.2, stage 1)"
V17_MODULE_FILE = "audit_surgery_complexity_duration.py"
V17_SPEC_VERSION = "surgery_complexity_duration_audit v1.7"
V17_MODULE_SHA256 = "1de2b17b2d894caff6a8362a9c72a458967691a8d2c7a533db863e87ad8dad94"  # public copy (cohort pins omitted); original private pin was ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0
V17_PINNED_SHA256 = {
    "canonical": "e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418",
    "harmonization": "ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d",
    "complexity_map": "5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5",
    "split": "60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31",
}
PINNED_SHA256 = {k: V17_PINNED_SHA256[k] for k in ("canonical", "harmonization", "split")}
MIN_CELL = 10
LAB_RELEASE_HOURS = 1.0         # decision 11 / A2.3: release = CHARTTIME + 1 h
LAYER_B_HOURS = 24.0
PRE_ICU_HOURS = 6.0             # pre-ICU window [t0 - 6 h, t0) (half-open; boundary at t0 to be confirmed)
LANDMARKS = (24, 48, 72, 120, 168, 240, 336)
EXPECTED_ELIGIBLE = {}  # public copy: landmark eligibility pins omitted
EXPECTED_RAW_ROWS = {}  # public copy: table-size pins omitted; see configs/cohort_pins.example.json
MAX_UNITS = 10
CHUNKSIZE = 1_000_000
SPLITS = ("train", "validation", "test")
DRUG_GROUPS = ("vasopressor", "inotrope", "antibiotic", "sedative_analgesic", "diuretic", "steroid", "nmb", "aed", "none")
SPECIMEN_GROUPS = ("blood", "urine", "respiratory", "csf", "other")
BLANK = "(blank)"

# Raw tables: time column, release offset, item key columns, unit and numeric candidates, dictionary.
TABLES = {
    "CHARTEVENTS": dict(time=("CHARTTIME",), offset=0.0, key="itemid", unit=("VALUEUOM",), numeric=("VALUENUM",), dictionary="D_ITEMS"),
    "LABEVENTS": dict(time=("CHARTTIME",), offset=LAB_RELEASE_HOURS, key="itemid", unit=("VALUEUOM",), numeric=("VALUENUM",), dictionary="D_LABITEMS"),
    "OUTPUTEVENTS": dict(time=("CHARTTIME",), offset=0.0, key="itemid", unit=("VALUEUOM",), numeric=("VALUE",), dictionary="D_ITEMS"),
    "INPUTEVENTS": dict(time=("CHARTTIME",), offset=0.0, key="itemid", unit=("AMOUNTUOM", "AMOUNT_UOM", "VALUEUOM"), numeric=("AMOUNT",), dictionary="D_ITEMS"),
    "PRESCRIPTIONS": dict(time=("STARTDATE", "STARTTIME"), offset=0.0, key="drug", unit=("DOSE_UNIT_RX",), numeric=(), dictionary=None),
    "MICROBIOLOGYEVENTS": dict(time=("CHARTTIME",), offset=0.0, key="micro", unit=(), numeric=(), dictionary=None),
}
# Inventories (an inventory = one keying of one raw table)
INVENTORIES = ("CHARTEVENTS", "LABEVENTS", "OUTPUTEVENTS", "INPUTEVENTS", "PRESCRIPTIONS",
               "MICROBIOLOGYEVENTS_specimen", "MICROBIOLOGYEVENTS_organism")
INV_TABLE = {i: i.split("_")[0] for i in INVENTORIES}
ROW_CLASSES = ("in_stay_train", "in_stay_validation", "in_stay_test", "cohort_admission_outside_stay",
               "no_cohort_admission", "cohort_admission_time_missing")
ITEM_COUNTS = ("rows", "rows_numeric", "encounters", "encounters_layer_b", "encounters_pre_icu_6h")
COUNT_KEYS = frozenset(set(ROW_CLASSES) | set(ITEM_COUNTS) | {"rows_total", "rows_not_in_dictionary", "train_encounters"})
NON_COUNT_TOP = ("disclosure", "suppression", "input_hashes", "v17_module", "raw_tables", "cohort_checks", "config", "search_hits")
TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}")

# Label / name search (keys for Vinai to verify; never an assignment). English and Chinese keywords.
SEARCH_TERMS = {
    "weight": r"weight|\bwt\b|体重",
    "height": r"height|身高|身长",
    "urine_output": r"urine|尿量|\bUO\b",
    "fluid_intake": r"intake|入量|输液量|fluid in",
    "gcs": r"\bgcs\b|glasgow|昏迷",
    "pupils": r"pupil|瞳孔",
    "fio2": r"fio2|fi02|inspired o|oxygen concentration|氧浓度",
    "ventilation": r"ventilat|\bpeep\b|tidal|呼吸机|潮气",
    "vasoactive": r"dopamine|dobutamine|epinephrine|adrenaline|norepinephrine|noradrenaline|milrinone|vasopressin|多巴|肾上腺素|米力农|血管加压素",
    "pain": r"pain|疼痛",
    "occult_blood": r"occult|潜血|隐血",
    "endocrine": r"cortisol|\btsh\b|thyroxine|free t4|\bft4\b|insulin|ketone|osmolal|皮质醇|促甲状腺|甲状腺素|胰岛素|酮体|渗透压",
    "sodium": r"sodium|\bna\b|钠",
    "nmb": r"vecuronium|rocuronium|atracurium|cisatracurium|pancuronium|succinylcholine|suxamethonium|维库溴铵|罗库溴铵|阿曲库铵|顺阿曲库铵|泮库溴铵|琥珀胆碱",
    "aed": r"levetiracetam|phenobarbital|phenytoin|valpro|carbamazepine|oxcarbazepine|topiramate|lamotrigine|lacosamide|左乙拉西坦|苯巴比妥|苯妥英|丙戊酸|卡马西平|奥卡西平|托吡酯|拉莫三嗪",
}
TEMPLATE_FILES = {
    "channel": "feature_channel_map_TEMPLATE.csv",
    "drug": "drug_group_map_TEMPLATE.csv",
    "specimen": "specimen_map_TEMPLATE.csv",
    "organism": "organism_map_TEMPLATE.csv",
}
TEMPLATE_COLUMNS = {
    "channel": ["table", "itemid", "label", "category", "units", "train_encounters", "train_rows",
                "channel", "table51_system", "triggered", "contract16", "unit_conversion", "plausibility_low",
                "plausibility_high", "notes"],
    "drug": ["drug_name", "drug_name_en", "units", "train_encounters", "group", "notes"],
    "specimen": ["spec_itemid", "spec_type_desc", "train_encounters", "specimen_group", "notes"],
    "organism": ["org_itemid", "org_name", "train_encounters", "canonical_name", "category", "contaminant_rule", "notes"],
}


class AuditError(Exception):
    """Input, hash or reconciliation failure (exit 2)."""


class HashMismatchError(AuditError):
    pass


class ReleaseGateError(Exception):
    """Release gate failure (exit 3)."""


# =============================================================================
# Pinned v1.7 module
# =============================================================================
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_v17(path: Path, expected_sha256: str = V17_MODULE_SHA256):
    path = Path(path)
    if not path.exists():
        raise AuditError(f"duration audit v1.7 module not found: {path} (copy {V17_MODULE_FILE} into the same folder)")
    if not re.fullmatch(r"[0-9a-f]{64}", str(expected_sha256)):
        raise HashMismatchError("pinned v1.7 module hash is not a full 64-character SHA-256")
    got = sha256_file(path)
    if got != expected_sha256:
        raise HashMismatchError(f"v1.7 module SHA-256 mismatch (expected {expected_sha256}, got {got}); "
                                "re-copy the reference file, do not edit it")
    name = f"_pinned_duration_audit_v17_{got[:12]}"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        sys.modules.pop(name, None)
        raise
    if getattr(mod, "SPEC_VERSION", None) != V17_SPEC_VERSION or getattr(mod, "PINNED_SHA256", None) != V17_PINNED_SHA256:
        raise AuditError("v1.7 module constants differ from this script's pinned values")
    return mod


def verify_hash(label: str, path: Path, expected: str) -> str:
    path = Path(path)
    if not path.exists():
        raise AuditError(f"{label} file not found: {path}")
    if not re.fullmatch(r"[0-9a-f]{64}", str(expected)):
        raise HashMismatchError(f"pinned {label} hash is not a full 64-character SHA-256")
    got = sha256_file(path)
    if got != expected:
        raise HashMismatchError(f"{label} SHA-256 mismatch (expected {expected}, got {got})")
    return got


# =============================================================================
# Config
# =============================================================================
@dataclass
class InvConfig:
    raw_dir: Path
    canonical_path: Path
    harmonization_path: Path
    split_path: Path
    out_dir: Path
    v17_module_path: Path
    pinned_sha256: dict = field(default_factory=lambda: dict(PINNED_SHA256))
    v17_module_sha256: str = V17_MODULE_SHA256
    harm_t0_col: str | None = None
    harm_out_col: str | None = None
    split_subject_col: str = "SUBJECT_ID"
    split_col: str = "split"
    expected_encounters: int | None = None  # public copy: pin omitted
    expected_subjects: int | None = None
    expected_split_encounters: dict = field(default_factory=dict)
    expected_eligible: dict = field(default_factory=lambda: dict(EXPECTED_ELIGIBLE))
    expected_raw_rows: dict = field(default_factory=lambda: dict(EXPECTED_RAW_ROWS))
    chunksize: int = CHUNKSIZE
    verbose: bool = True

    @property
    def shared_json(self) -> Path:
        return self.out_dir / "feature_item_inventory.json"

    @property
    def shared_md(self) -> Path:
        return self.out_dir / "feature_item_inventory.md"

    @property
    def local_json(self) -> Path:
        return self.out_dir / "feature_item_inventory_raw_LOCAL.json"

    def template(self, kind: str) -> Path:
        return self.out_dir / TEMPLATE_FILES[kind]


def log(cfg: InvConfig, msg: str) -> None:
    if cfg.verbose:
        print(msg, flush=True)


# =============================================================================
# Pure helpers (unit tested)
# =============================================================================
_NULLS = {"", "nan", "none", "null", "nat"}


def vnorm_id(s: pd.Series) -> pd.Series:
    """Vectorized equivalent of v1.7 norm_id for CSV text: trims, drops quotes, '123.0' -> '123'; other text unchanged."""
    s = s.astype(str).str.strip().str.strip("'\"").str.strip()
    s = s.where(~s.str.lower().isin(_NULLS), "")
    m = s.str.fullmatch(r"[+-]?\d+\.0*")
    if m.any():  # exactly as v1.7 norm_id: '0123.0' -> '123', but '0123' stays '0123'
        u = {v: str(int(v.split(".")[0])) for v in pd.unique(s[m])}
        s = s.where(~m, s.map(u))
    return s


def vclean_text(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
    return s.where(~s.str.lower().isin(_NULLS), "")


def fast_time(s: pd.Series) -> pd.Series:
    """Parse with the PIC format first; only the remainder goes through the slow mixed-format parser."""
    s = s.astype(str).str.strip()
    t = pd.to_datetime(s, format="%Y-%m-%d %H:%M:%S", errors="coerce")
    rest = t.isna() & ~s.str.lower().isin(_NULLS)
    if rest.any():
        try:
            t2 = pd.to_datetime(s[rest], errors="coerce", format="mixed")
        except (TypeError, ValueError):
            t2 = pd.to_datetime(s[rest], errors="coerce")
        t.loc[rest] = t2
    return t


def item_sort_key(k: str) -> tuple:
    return (0, int(k), "") if re.fullmatch(r"\d+", k) else (1, 0, k)


def resolve(header: list[str], candidates: Iterable[str], label: str, required: bool = True) -> str | None:
    up = {h.upper(): h for h in header}
    for c in candidates:
        if c.upper() in up:
            return up[c.upper()]
    if required:
        raise AuditError(f"{label}: none of the columns {list(candidates)} found; header has {header}")
    return None


def md_cell(v: Any) -> str:
    if isinstance(v, list):
        v = "; ".join(str(x) for x in v)
    return str(v).replace("|", "/").replace("\n", " ").strip()


# =============================================================================
# Inputs: cohort, split, landmark eligibility
# =============================================================================
def load_cohort(cfg: InvConfig, A) -> tuple[pd.DataFrame, dict]:
    hashes = {k: verify_hash(k, getattr(cfg, f"{k}_path"), cfg.pinned_sha256[k]) for k in ("canonical", "harmonization", "split")}
    log(cfg, "input hash checks: canonical OK, harmonization OK, split OK (exact 64-char)")
    checks: dict[str, Any] = {}
    ch = A.read_header(cfg.canonical_path)
    need = ["encounter_id", "SUBJECT_ID", "HADM_ID", "hadm_id_list"]
    miss = [c for c in need if c not in ch]
    if miss:
        raise AuditError(f"canonical CSV missing columns: {miss}")
    canon = A.read_csv_str(cfg.canonical_path, usecols=need)
    canon["encounter_id"] = canon["encounter_id"].map(lambda v: str(v).strip())
    if (canon["encounter_id"] == "").any() or canon["encounter_id"].duplicated().any():
        raise AuditError("canonical encounter_id empty or not unique")
    canon["SUBJECT_ID"] = canon["SUBJECT_ID"].map(A.norm_id)
    canon["HADM_ID"] = canon["HADM_ID"].map(A.norm_id)
    canon["hadm_list"] = canon["hadm_id_list"].map(A.parse_hadm_list)
    empty = canon["hadm_list"].map(len) == 0
    canon.loc[empty, "hadm_list"] = pd.Series([[h] if h else [] for h in canon.loc[empty, "HADM_ID"]],
                                              index=canon.index[empty], dtype=object)
    checks["hadm_list_empty_filled_from_HADM_ID"] = int(empty.sum())

    hh = A.read_header(cfg.harmonization_path)
    if "encounter_id" not in hh:
        raise AuditError("harmonization file has no encounter_id column")
    try:
        t0c = A.resolve_column(hh, cfg.harm_t0_col, A.HARM_T0_CANDIDATES, "T0", A.HARM_T0_FALLBACK)
        outc = A.resolve_column(hh, cfg.harm_out_col, A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK)
    except A.AuditError as e:
        raise AuditError(str(e)) from e
    harm = A.read_csv_str(cfg.harmonization_path, usecols=["encounter_id", t0c, outc])
    harm["encounter_id"] = harm["encounter_id"].map(lambda v: str(v).strip())
    if harm["encounter_id"].duplicated().any():
        raise AuditError("harmonization encounter_id not unique")
    if set(harm["encounter_id"]) != set(canon["encounter_id"]):
        raise AuditError("harmonization encounter_id set differs from canonical encounter_id set")
    harm = pd.DataFrame({"encounter_id": harm["encounter_id"], "t0": fast_time(harm[t0c]), "OUT": fast_time(harm[outc])})
    if harm["t0"].isna().any() or harm["OUT"].isna().any():
        raise AuditError("harmonized t0 or OUT missing/unparseable for some encounters")
    if (harm["OUT"] < harm["t0"]).any():
        raise AuditError("harmonized OUT earlier than t0 for some encounters")
    checks["harmonization_columns"] = {"t0": t0c, "OUT": outc}
    enc = canon.merge(harm, on="encounter_id", how="left")

    sh = A.read_header(cfg.split_path)
    for c in (cfg.split_subject_col, cfg.split_col):
        if c not in sh:
            raise AuditError(f"split file missing column {c!r}; header has {sh}")
    spl = A.read_csv_str(cfg.split_path, usecols=[cfg.split_subject_col, cfg.split_col])
    spl = pd.DataFrame({"SUBJECT_ID": spl[cfg.split_subject_col].map(A.norm_id),
                        "split": spl[cfg.split_col].map(lambda v: A.SPLIT_LABELS.get(str(v).strip().lower(), "?"))})
    if (spl["split"] == "?").any():
        raise AuditError("split file has unknown split labels")
    spl = spl.drop_duplicates()
    if spl["SUBJECT_ID"].duplicated().any():
        raise AuditError("a subject appears in more than one split")
    enc = enc.merge(spl, on="SUBJECT_ID", how="left")
    if enc["split"].isna().any():
        raise AuditError(f"{int(enc['split'].isna().sum())} cohort encounters have no split")
    n_enc, n_subj = len(enc), enc["SUBJECT_ID"].nunique()
    sc = {s: int((enc["split"] == s).sum()) for s in SPLITS}
    if n_enc != cfg.expected_encounters or n_subj != cfg.expected_subjects or sc != cfg.expected_split_encounters:
        raise AuditError(f"cohort/split totals differ from spec: encounters {n_enc}, subjects {n_subj}, splits {sc}")
    elig = {}
    for L in LANDMARKS:
        elig[L] = int((enc["OUT"] > enc["t0"] + pd.Timedelta(hours=L)).sum())
    exp = {int(k): int(v) for k, v in cfg.expected_eligible.items()}
    if elig != exp:
        raise AuditError(f"core-landmark eligible counts do not reproduce: got {elig}, expected {exp}")
    checks.update({"encounters": n_enc, "subjects": n_subj, "split_encounters": sc,
                   "eligible_at_landmark": {f"{L}h": elig[L] for L in LANDMARKS}})
    log(cfg, f"cohort ok ({n_enc} encounters, {n_subj} subjects, splits {sc}); landmark eligibility reproduces")
    enc = enc.reset_index(drop=True)
    enc["enc_i"] = np.arange(len(enc))
    return enc, {"hashes": hashes, "checks": checks}


def encounter_hadm_index(enc: pd.DataFrame) -> pd.DataFrame:
    rows = [(f"{s}|{h}", i) for s, lst, i in zip(enc["SUBJECT_ID"], enc["hadm_list"], enc["enc_i"]) for h in lst]
    eh = pd.DataFrame(rows, columns=["key", "enc_i"]).drop_duplicates()
    return eh.merge(enc[["enc_i", "t0", "OUT", "split"]], on="enc_i", how="left")


# =============================================================================
# Dictionaries
# =============================================================================
def load_dictionary(path: Path, A, table: str) -> tuple[dict, str]:
    if not Path(path).exists():
        raise AuditError(f"{table} not found: {path}")
    digest = sha256_file(path)
    hdr = A.read_header(path)
    idc = resolve(hdr, ("ITEMID",), table)
    labc = resolve(hdr, ("LABEL", "LABEL_EN"), table)
    catc = resolve(hdr, ("CATEGORY",), table, required=False)
    fluc = resolve(hdr, ("FLUID",), table, required=False)
    cols = [c for c in (idc, labc, catc, fluc) if c]
    d = A.read_csv_str(path, usecols=cols)
    out = {}
    for _, r in d.iterrows():
        k = A.norm_id(r[idc])
        if not k:
            continue
        out[k] = {"label": md_cell(r[labc]), "category": md_cell(r[catc]) if catc else "",
                  "fluid": md_cell(r[fluc]) if fluc else ""}
    return out, digest


# =============================================================================
# Inventory pass over one raw table (chunked)
# =============================================================================
class TableInventory:
    def __init__(self, name: str, keyings: list[str]):
        self.name = name
        self.keyings = keyings
        self.row_classes = {c: 0 for c in ROW_CLASSES}
        self.rows_total = 0
        self.table_enc: set = set()
        self.item_rows = {k: {} for k in keyings}
        self.item_numeric = {k: {} for k in keyings}
        self.pairs = {k: {w: [] for w in ("stay", "layer_b", "pre")} for k in keyings}
        self.units = {k: {} for k in keyings}
        self.fields = {k: {} for k in keyings}
        self.columns_used: list[str] = []

    @staticmethod
    def _add(d: dict, s: pd.Series) -> None:
        for k, v in s.items():
            d[k] = d.get(k, 0) + int(v)

    def process(self, ch: pd.DataFrame, eh: pd.DataFrame, offset_h: float) -> None:
        n = len(ch)
        self.rows_total += n
        ch = ch.reset_index(drop=True)
        ch["r"] = np.arange(n)
        ch["key"] = ch["subj"] + "|" + ch["hadm"]
        ok = (ch["subj"] != "") & (ch["hadm"] != "")
        m = ch.loc[ok, ["r", "key", "t"]].merge(eh, on="key", how="inner")
        cohort_r = pd.unique(m["r"])
        self.row_classes["no_cohort_admission"] += n - len(cohort_r)
        if offset_h:
            m["t"] = m["t"] + pd.Timedelta(hours=offset_h)
        tmiss_r = pd.unique(m.loc[m["t"].isna(), "r"])
        self.row_classes["cohort_admission_time_missing"] += len(tmiss_r)
        mv = m[m["t"].notna()]
        stay = mv[(mv["t"] >= mv["t0"]) & (mv["t"] < mv["OUT"])].drop_duplicates("r")
        for sp in SPLITS:
            self.row_classes[f"in_stay_{sp}"] += int((stay["split"] == sp).sum())
        self.row_classes["cohort_admission_outside_stay"] += len(cohort_r) - len(tmiss_r) - len(stay)
        hb = pd.Timedelta(hours=LAYER_B_HOURS)
        pre = mv[(mv["split"] == "train") & (mv["t"] >= mv["t0"] - pd.Timedelta(hours=PRE_ICU_HOURS)) & (mv["t"] < mv["t0"])]
        st = stay[stay["split"] == "train"]
        lb = st[st["t"] < np.minimum(st["OUT"], st["t0"] + hb)]
        self.table_enc.update(st["enc_i"].tolist())
        for kg in self.keyings:
            col = f"k_{kg}"
            s_items = ch.loc[st["r"].to_numpy(), [col, "unit", "num"]].assign(enc_i=st["enc_i"].to_numpy())
            self._add(self.item_rows[kg], s_items.groupby(col).size())
            self._add(self.item_numeric[kg], s_items[s_items["num"]].groupby(col).size())
            self.pairs[kg]["stay"].append(s_items[[col, "enc_i"]].drop_duplicates().to_numpy())
            self.pairs[kg]["layer_b"].append(pd.DataFrame({col: ch.loc[lb["r"].to_numpy(), col].to_numpy(),
                                                           "enc_i": lb["enc_i"].to_numpy()}).drop_duplicates().to_numpy())
            self.pairs[kg]["pre"].append(pd.DataFrame({col: ch.loc[pre["r"].to_numpy(), col].to_numpy(),
                                                       "enc_i": pre["enc_i"].to_numpy()}).drop_duplicates().to_numpy())
            uu = s_items[s_items["unit"] != ""][[col, "unit"]].drop_duplicates()
            for k, u in zip(uu[col], uu["unit"]):
                us = self.units[kg].setdefault(k, set())
                if len(us) < 50:
                    us.add(u)
            fc = [c for c in ch.columns if c.startswith(f"f_{kg}_")]
            if fc:
                ff = ch.loc[st["r"].to_numpy(), [col] + fc].drop_duplicates(col)
                for _, r in ff.iterrows():
                    self.fields[kg].setdefault(r[col], {c[len(f"f_{kg}_"):]: r[c] for c in fc})

    def finish(self) -> None:
        self.enc_counts = {}
        for kg in self.keyings:
            res = {}
            for w, parts in self.pairs[kg].items():
                parts = [p for p in parts if len(p)]
                if parts:
                    arr = pd.DataFrame(np.concatenate(parts), columns=["k", "e"]).drop_duplicates()
                    res[w] = arr.groupby("k").size().to_dict()
                else:
                    res[w] = {}
            self.enc_counts[kg] = res
        self.pairs = None


def read_table(cfg: InvConfig, A, name: str, inv: TableInventory, eh: pd.DataFrame) -> dict:
    spec = TABLES[name]
    path = Path(cfg.raw_dir) / f"{name}.csv"
    if not path.exists():
        raise AuditError(f"raw table not found: {path}")
    digest = sha256_file(path)
    hdr = A.read_header(path)
    sc = resolve(hdr, ("SUBJECT_ID",), name)
    hc = resolve(hdr, ("HADM_ID",), name)
    tc = resolve(hdr, spec["time"], name)
    uc = resolve(hdr, spec["unit"], name, required=False) if spec["unit"] else None
    nc = resolve(hdr, spec["numeric"], name, required=False) if spec["numeric"] else None
    keycols: dict[str, list[str]] = {}
    if spec["key"] == "itemid":
        ic = resolve(hdr, ("ITEMID",), name, required=False)
        keycols[name] = [ic] if ic else []
    elif spec["key"] == "drug":
        dn = resolve(hdr, ("DRUG_NAME",), name, required=False)
        de = resolve(hdr, ("DRUG_NAME_EN",), name, required=False)
        if not dn and not de:
            raise AuditError(f"{name}: neither DRUG_NAME nor DRUG_NAME_EN found")
        keycols[name] = [c for c in (dn, de)]
    else:
        keycols["MICROBIOLOGYEVENTS_specimen"] = [resolve(hdr, ("SPEC_ITEMID",), name, required=False),
                                                  resolve(hdr, ("SPEC_TYPE_DESC",), name)]
        keycols["MICROBIOLOGYEVENTS_organism"] = [resolve(hdr, ("ORG_ITEMID",), name, required=False),
                                                  resolve(hdr, ("ORG_NAME",), name)]
    used = sorted({c for c in [sc, hc, tc, uc, nc] + [c for v in keycols.values() for c in v] if c})
    inv.columns_used = used
    reader = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False, usecols=used,
                         encoding="utf-8-sig", chunksize=cfg.chunksize)
    for raw in reader:
        ch = pd.DataFrame({"subj": vnorm_id(raw[sc]), "hadm": vnorm_id(raw[hc]), "t": fast_time(raw[tc]),
                           "unit": vclean_text(raw[uc]) if uc else "",
                           "num": (pd.to_numeric(raw[nc].str.strip(), errors="coerce").notna() if nc
                                   else pd.Series(False, index=raw.index))})
        for kg, cols in keycols.items():
            if spec["key"] == "itemid":
                if cols:
                    k = vnorm_id(raw[cols[0]])
                else:  # no ITEMID column: key by unit
                    k = "unit:" + ch["unit"]
            elif spec["key"] == "drug":
                a = vclean_text(raw[cols[0]]) if cols[0] else pd.Series("", index=raw.index)
                b = vclean_text(raw[cols[1]]) if cols[1] else pd.Series("", index=raw.index)
                k = a + " | " + b
                ch[f"f_{kg}_drug_name"] = a.to_numpy()
                ch[f"f_{kg}_drug_name_en"] = b.to_numpy()
            else:
                iid = vnorm_id(raw[cols[0]]) if cols[0] else pd.Series("", index=raw.index)
                nm = vclean_text(raw[cols[1]])
                k = iid + " | " + nm
                pre = "spec" if kg.endswith("specimen") else "org"
                ch[f"f_{kg}_{pre}_itemid"] = iid.to_numpy()
                ch[f"f_{kg}_{pre}_name"] = nm.to_numpy()
            k = k.where(k.str.strip(" |") != "", BLANK)
            ch[f"k_{kg}"] = k.to_numpy()
        inv.process(ch, eh, spec["offset"])
    inv.finish()
    if sum(inv.row_classes.values()) != inv.rows_total:
        raise AuditError(f"{name}: row classes do not sum to the rows read")
    exp = cfg.expected_raw_rows.get(name)
    if exp is not None and inv.rows_total != exp:
        raise AuditError(f"{name}: {inv.rows_total} rows read, schema inventory says {exp}")
    log(cfg, f"{name}: {inv.rows_total} rows read, SHA-256 {digest[:12]}")
    return {"file": path.name, "sha256": digest, "rows_read": inv.rows_total, "columns_used": used}


# =============================================================================
# Registry construction
# =============================================================================
def build_registry(A, invs: dict, dicts: dict, n_train: int) -> tuple[Any, dict]:
    reg = A.Registry()
    P = A.P
    den = reg.count(P("denominators", "train_encounters"), n_train)
    meta: dict[str, dict] = {}
    for tname, inv in invs.items():
        tot = reg.count(P("tables", tname, "row_classes", "rows_total"), inv.rows_total)
        parts = [reg.count(P("tables", tname, "row_classes", c), inv.row_classes[c]) for c in ROW_CLASSES]
        reg.identity(tot, parts)
        te = reg.count(P("tables", tname, "train", "encounters"), len(inv.table_enc))
        reg.subset(te, den, "table_enc_in_cohort")
        rows_train = P("tables", tname, "row_classes", "in_stay_train")
        for kg in inv.keyings:
            items = sorted(set(inv.item_rows[kg]), key=item_sort_key)
            meta[kg] = {}
            dct = dicts.get(TABLES[tname]["dictionary"]) if TABLES[tname]["dictionary"] else None
            rcells, nd_cells = [], []
            for it in items:
                base = ("inventories", kg, "items", it)
                rc = reg.count(P(*base, "train", "rows"), inv.item_rows[kg][it])
                rcells.append(rc)
                rn = reg.count(P(*base, "train", "rows_numeric"), inv.item_numeric[kg].get(it, 0))
                reg.subset(rn, rc, "numeric")
                ec = reg.count(P(*base, "train", "encounters"), inv.enc_counts[kg]["stay"].get(it, 0))
                reg.subset(ec, te, "item_enc")
                lb = reg.count(P(*base, "train", "encounters_layer_b"), inv.enc_counts[kg]["layer_b"].get(it, 0))
                reg.subset(lb, ec, "layer_b")
                pr = reg.count(P(*base, "train", "encounters_pre_icu_6h"), inv.enc_counts[kg]["pre"].get(it, 0))
                reg.subset(pr, den, "pre_icu")
                info: dict[str, Any] = {"units": sorted(inv.units[kg].get(it, set()))[:MAX_UNITS]}
                if dct is not None:
                    d = dct.get(it)
                    info["in_dictionary"] = d is not None
                    info["label"] = d["label"] if d else ""
                    info["category"] = d["category"] if d else ""
                    if TABLES[tname]["dictionary"] == "D_LABITEMS":
                        info["fluid"] = d["fluid"] if d else ""
                    if d is None:
                        nd_cells.append(rc)
                else:
                    info.update({k: md_cell(v) for k, v in inv.fields[kg].get(it, {}).items()})
                meta[kg][it] = info
            reg.identity(rows_train, rcells)
            if dct is not None and nd_cells:  # only when some item is missing from the dictionary (a 0 = sum() identity would be recoverable)
                nd = reg.count(P("inventories", kg, "rows_not_in_dictionary"), sum(reg.raw[c] for c in nd_cells))
                reg.identity(nd, nd_cells)
    return reg, meta


# =============================================================================
# Document, rendering
# =============================================================================
def label_text(inv: str, info: dict) -> str:
    if inv == "PRESCRIPTIONS":
        return f"{info.get('drug_name', '')} {info.get('drug_name_en', '')}"
    if inv.endswith("specimen"):
        return info.get("spec_name", "")
    if inv.endswith("organism"):
        return info.get("org_name", "")
    return f"{info.get('label', '')} {info.get('category', '')} {info.get('fluid', '')}"


def search_hits(meta: dict) -> dict:
    out: dict[str, list] = {}
    for concept, pat in SEARCH_TERMS.items():
        rx = re.compile(pat, re.IGNORECASE)
        hits = []
        for inv in INVENTORIES:
            for it, info in meta.get(inv, {}).items():
                if rx.search(label_text(inv, info)):
                    hits.append(f"{inv}:{it}")
        out[concept] = hits
    return out


def build_document(A, reg, status: dict, supp: dict, meta: dict, hashes: dict, raw_tables: dict, checks: dict,
                   v17_sha: str) -> dict:
    doc: dict[str, Any] = {
        "spec": SPEC_VERSION, "release_gate": "pending",
        "input_hashes": hashes, "v17_module": {"file": V17_MODULE_FILE, "sha256": v17_sha},
        "raw_tables": raw_tables, "cohort_checks": checks,
        "config": {"lab_release_hours": LAB_RELEASE_HOURS, "layer_b_hours": LAYER_B_HOURS,
                   "pre_icu_window": f"[t0 - {PRE_ICU_HOURS:g} h, t0)", "split_published": "train",
                   "min_cell": MIN_CELL, "drug_groups": list(DRUG_GROUPS), "specimen_groups": list(SPECIMEN_GROUPS)},
        "search_hits": search_hits(meta),
        "denominators": {}, "tables": {}, "inventories": {},
    }
    for inv in INVENTORIES:
        if inv in meta:
            doc["inventories"].setdefault(inv, {})["items"] = {}
            for it, info in meta[inv].items():
                doc["inventories"][inv]["items"][it] = dict(info)
    for c in reg.published():
        st = status[c]
        A.set_ptr(doc, c, reg.raw[c] if st == "pub" else ("<10" if st == "P" else "S"))
    # units only where the item's rows are published as >= 10
    for inv, d in doc["inventories"].items():
        for it, o in d["items"].items():
            if not isinstance(o["train"]["rows"], int):
                o["units"] = []
    published = reg.published()
    doc["disclosure"] = {"count_fields": sorted(published),
                         "identities": [[t, list(p)] for t, p in reg.identities],
                         "hidden_cells": sorted(reg.hidden)}
    doc["suppression"] = dict(supp)
    return doc


def detect_count_fields(obj: Any, prefix: str = "") -> list[str]:
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if prefix == "" and k in NON_COUNT_TOP:
                continue
            ptr = prefix + "/" + str(k).replace("~", "~0").replace("/", "~1")
            if k in COUNT_KEYS and not isinstance(v, (dict, list)):
                out.append(ptr)
            else:
                out.extend(detect_count_fields(v, ptr))
    return out


def _f(v: Any) -> str:
    return "–" if v is None else str(v)


def render_markdown(doc: dict) -> str:
    L = ["# Feature item inventory (stage 1)", "",
         f"- spec: {doc['spec']}", f"- release gate: {doc['release_gate']}",
         "- aggregate only; train split; counts 0–9 shown as <10 (zero included, v1.7 rule), complementary suppression shown as S; items listed in key order",
         "- input hashes: " + ", ".join(f"{k} {v[:12]}" for k, v in doc["input_hashes"].items()),
         f"- v1.7 module: {doc['v17_module']['sha256'][:12]}", ""]
    c = doc["cohort_checks"]
    L += ["## Cohort checks", "", f"- encounters {c['encounters']}, subjects {c['subjects']}, splits " +
          ", ".join(f"{k} {v}" for k, v in c["split_encounters"].items()),
          "- eligible at landmark: " + ", ".join(f"{k} {v}" for k, v in c["eligible_at_landmark"].items()),
          f"- train encounters (denominator): {_f(doc['denominators'].get('train_encounters'))}", ""]
    L += ["## Raw tables read", "", "| table | file | sha256 | rows_read |", "|---|---|---|---|"]
    for t, o in doc["raw_tables"].items():
        L.append(f"| {t} | {o['file']} | {o['sha256']} | {o['rows_read']} |")
    L += ["", "## Row reconciliation (all splits)", "",
          "| table | rows_total | " + " | ".join(ROW_CLASSES) + " | train encounters |",
          "|---" * (len(ROW_CLASSES) + 3) + "|"]
    for t, o in doc["tables"].items():
        rc = o["row_classes"]
        L.append(f"| {t} | {_f(rc['rows_total'])} | " + " | ".join(_f(rc[k]) for k in ROW_CLASSES)
                 + f" | {_f(o['train']['encounters'])} |")
    L += ["", "## Search hits (labels and names; keys to verify)", ""]
    for k, v in doc["search_hits"].items():
        L.append(f"- {k}: " + (", ".join(v) if v else "none"))
    for inv, d in doc["inventories"].items():
        L += ["", f"## {inv}", ""]
        if "rows_not_in_dictionary" in d:
            L += [f"- rows of items not in the dictionary: {_f(d['rows_not_in_dictionary'])}", ""]
        L += ["| key | name | units | rows | rows_numeric | encounters | encounters_layer_b | encounters_pre_icu_6h |",
              "|---|---|---|---|---|---|---|---|"]
        for it, o in d["items"].items():
            tr = o["train"]
            L.append(f"| {md_cell(it)} | {md_cell(label_text(inv, o).strip())} | {md_cell(o['units'])} | "
                     + " | ".join(_f(tr[k]) for k in ITEM_COUNTS) + " |")
    s = doc["suppression"]
    L += ["", "## Suppression", "", f"- primary (<10): {s.get('primary_cells')}; complementary (S): {s.get('complementary_cells')}; "
          f"recoverable: {s.get('recoverable')}", ""]
    return "\n".join(L) + "\n"


def template_rows(doc: dict) -> dict[str, list[list]]:
    """Only items whose train encounter count is >= 10 (published value or complementary S)."""
    out: dict[str, list[list]] = {k: [] for k in TEMPLATE_FILES}
    for inv, d in doc["inventories"].items():
        for it, o in d["items"].items():
            e = o["train"]["encounters"]
            if e == "<10":
                continue
            u = "; ".join(o.get("units", []))
            if inv in ("CHARTEVENTS", "LABEVENTS", "OUTPUTEVENTS", "INPUTEVENTS"):
                out["channel"].append([inv, it, o.get("label", ""), o.get("category", ""), u, e, o["train"]["rows"]]
                                      + [""] * 8)
            elif inv == "PRESCRIPTIONS":
                out["drug"].append([o.get("drug_name", ""), o.get("drug_name_en", ""), u, e, "", ""])
            elif inv.endswith("specimen"):
                out["specimen"].append([o.get("spec_itemid", ""), o.get("spec_name", ""), e, "", ""])
            else:
                out["organism"].append([o.get("org_itemid", ""), o.get("org_name", ""), e, "", "", "", ""])
    return out


def render_templates(doc: dict) -> dict[str, str]:
    rows = template_rows(doc)
    res = {}
    for kind, cols in TEMPLATE_COLUMNS.items():
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows[kind])
        res[kind] = buf.getvalue()
    return res


# =============================================================================
# Release gate (disk reload) and publishing
# =============================================================================
def release_gate(A, json_path: Path, md_path: Path, tpl_paths: dict, expect_status: str,
                 sentinels: Iterable[str] = ()) -> list[str]:
    fails: list[str] = []
    try:
        doc = json.loads(Path(json_path).read_text(encoding="utf-8"))
        md = Path(md_path).read_text(encoding="utf-8")
        tpls = {k: Path(p).read_text(encoding="utf-8") for k, p in tpl_paths.items()}
    except Exception as e:  # noqa: BLE001
        return [f"cannot reload shared files: {e}"]
    if doc.get("spec") != SPEC_VERSION or doc.get("release_gate") != expect_status:
        fails.append("header fields (spec / release_gate) wrong")
    for p in [json_path, md_path, *tpl_paths.values()]:
        if "_raw_LOCAL" in Path(p).name:
            fails.append("a shared file name contains _raw_LOCAL")
    declared = set(doc.get("disclosure", {}).get("count_fields", []))
    detected = set(detect_count_fields(doc))
    if declared != detected:
        fails.append(f"count-field schema mismatch: {len(detected - declared)} undeclared, {len(declared - detected)} missing")
    known: dict[str, int] = {}
    unknown: dict[str, str] = {}
    n_p = n_s = 0
    for ptr in sorted(detected):
        try:
            v = A.get_ptr(doc, ptr)
        except KeyError:
            fails.append(f"declared count field missing at {ptr}")
            continue
        if isinstance(v, bool) or not (isinstance(v, int) or v in ("<10", "S")):
            fails.append(f"count field with an invalid value at {ptr}")
        elif isinstance(v, int):
            if v < MIN_CELL:
                fails.append(f"published count below 10 at {ptr}")
            known[ptr] = v
        elif v == "<10":
            unknown[ptr] = "P"; n_p += 1
        else:
            unknown[ptr] = "S"; n_s += 1
    for h in doc.get("disclosure", {}).get("hidden_cells", []):
        if not str(h).startswith("hidden:"):
            fails.append("hidden cell without the hidden: prefix")
        unknown[h] = "H"
    sup = doc.get("suppression", {})
    if sup.get("primary_cells") != n_p or sup.get("complementary_cells") != n_s or sup.get("recoverable") != 0:
        fails.append("suppression counters do not match the published cells")
    ids = [(t, tuple(p)) for t, p in doc.get("disclosure", {}).get("identities", [])]
    if not ids:
        fails.append("no disclosure identities published")
    if not fails:
        try:
            checked, unsafe, _ = A.recoverability_check(ids, known, unknown)
            if unsafe:
                fails.append(f"recoverability: {len(unsafe)} suppressed cell(s) recoverable or bounded inside [1, 9]")
            if checked != sup.get("recoverability_cells_checked"):
                fails.append("recoverability_cells_checked differs from the recomputed value")
        except A.ReleaseGateError as e:
            fails.append(f"recoverability solver: {e}")

    def walk(o: Any, path: str = "") -> None:
        if isinstance(o, dict):
            for k, v in o.items():
                if str(k).lower() in ("min", "max"):
                    fails.append(f"min/max published at {path}/{k}")
                walk(v, f"{path}/{k}")
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}/{i}")
        elif isinstance(o, str) and TIMESTAMP_RE.search(o):
            fails.append(f"timestamp-like string at {path}")
    walk({k: v for k, v in doc.items() if k != "disclosure"})
    for inv, d in doc.get("inventories", {}).items():
        for it, o in d.get("items", {}).items():
            if o.get("units") and not isinstance(o.get("train", {}).get("rows"), int):
                fails.append(f"units published for an item whose rows are suppressed ({inv})")
    try:
        if render_markdown(doc) != md:
            fails.append("markdown on disk is not the rendering of the shared JSON")
        rt = render_templates(doc)
        for k, text in tpls.items():
            if rt[k] != text:
                fails.append(f"template {k} on disk is not the rendering of the shared JSON")
    except Exception as e:  # noqa: BLE001
        fails.append(f"shared files cannot be re-rendered: {e}")
    for k, text in tpls.items():
        rdr = csv.DictReader(io.StringIO(text))
        for row in rdr:
            v = row.get("train_encounters", "")
            if v == "<10" or (v.isdigit() and int(v) < MIN_CELL):
                fails.append(f"template {k} lists an item with fewer than 10 train encounters")
                break
    text = json.dumps(doc, ensure_ascii=False) + md + "".join(tpls.values())
    for sv in sentinels:
        if sv and re.search(rf"(?<![0-9A-Za-z]){re.escape(str(sv))}(?![0-9A-Za-z])", text):
            fails.append("sentinel identifier found in shared output")
    headers = None
    for line in md.splitlines():
        if line.startswith("|") and not line.startswith("|---"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if headers is None:
                headers = cells
                continue
            for h, c in zip(headers, cells):
                if h in COUNT_KEYS and c not in ("<10", "S", "", "–") and not (c.isdigit() and int(c) >= MIN_CELL):
                    fails.append(f"markdown count column '{h}' has a disallowed value")
        elif not line.startswith("|---"):
            headers = None
    return fails


def _atomic_write(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8", newline="")
    os.replace(tmp, path)


def _remove(*paths: Path) -> None:
    for p in paths:
        try:
            Path(p).unlink()
        except FileNotFoundError:
            pass


def _pending(p: Path) -> Path:
    stem, ext = p.name.rsplit(".", 1)
    return p.with_name(f"{stem}.pending.{ext}")


def write_set(json_p: Path, md_p: Path, tpl: dict, doc: dict) -> None:
    _atomic_write(json_p, json.dumps(doc, indent=1, ensure_ascii=False))
    _atomic_write(md_p, render_markdown(doc))
    for k, text in render_templates(doc).items():
        _atomic_write(tpl[k], text)


def publish(A, cfg: InvConfig, doc: dict, sentinels: Iterable[str] = ()) -> None:
    ftpl = {k: cfg.template(k) for k in TEMPLATE_FILES}
    ptpl = {k: _pending(p) for k, p in ftpl.items()}
    pj, pm = _pending(cfg.shared_json), _pending(cfg.shared_md)
    allp = [cfg.shared_json, cfg.shared_md, *ftpl.values(), pj, pm, *ptpl.values()]
    _remove(*allp)
    write_set(pj, pm, ptpl, doc)
    fails = release_gate(A, pj, pm, ptpl, "pending", sentinels)
    if fails:
        _remove(*allp)
        raise ReleaseGateError("; ".join(fails[:20]))
    final = json.loads(pj.read_text(encoding="utf-8"))
    final["release_gate"] = "ok"
    write_set(cfg.shared_json, cfg.shared_md, ftpl, final)
    _remove(pj, pm, *ptpl.values())
    fails = release_gate(A, cfg.shared_json, cfg.shared_md, ftpl, "ok", sentinels)
    if fails:
        _remove(*allp)
        raise ReleaseGateError("; ".join(fails[:20]))


# =============================================================================
# Orchestration
# =============================================================================
def run_inventory(cfg: InvConfig, sentinels: Iterable[str] = ()) -> dict:
    A = load_v17(cfg.v17_module_path, cfg.v17_module_sha256)
    v17_sha = cfg.v17_module_sha256
    try:
        enc, info = load_cohort(cfg, A)
        eh = encounter_hadm_index(enc)
        n_train = int((enc["split"] == "train").sum())
        raw_tables: dict[str, dict] = {}
        dicts: dict[str, dict] = {}
        for dn in ("D_ITEMS", "D_LABITEMS"):
            dicts[dn], dg = load_dictionary(Path(cfg.raw_dir) / f"{dn}.csv", A, dn)
            raw_tables[dn] = {"file": f"{dn}.csv", "sha256": dg, "rows_read": len(dicts[dn]), "columns_used": []}
        invs: dict[str, TableInventory] = {}
        for t in TABLES:
            keyings = [k for k in INVENTORIES if INV_TABLE[k] == t]
            inv = TableInventory(t, keyings)
            raw_tables[t] = read_table(cfg, A, t, inv, eh)
            invs[t] = inv
        reg, meta = build_registry(A, invs, dicts, n_train)
        log(cfg, f"registry: {len(reg.raw)} cells, {len(reg.identities)} identities; suppressing ...")
        status, supp = A.apply_suppression(reg)
    except A.HashMismatchError as e:
        raise HashMismatchError(str(e)) from e
    except A.AuditError as e:
        raise AuditError(str(e)) from e
    except A.ReleaseGateError as e:
        raise ReleaseGateError(str(e)) from e
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    local = {"WARNING": "LOCAL ONLY - unsuppressed aggregates; never copy, upload or share",
             "raw_counts": reg.raw, "checks": info["checks"]}
    _atomic_write(cfg.local_json, json.dumps(local, indent=1, default=str, ensure_ascii=False))
    doc = build_document(A, reg, status, supp, meta, info["hashes"], raw_tables, info["checks"], v17_sha)
    publish(A, cfg, doc, sentinels)
    final = json.loads(cfg.shared_json.read_text(encoding="utf-8"))
    print_console_summary(cfg, final)
    return final


def print_console_summary(cfg: InvConfig, doc: dict) -> None:
    """Short; prints only published values."""
    if not cfg.verbose:
        return
    print(f"=== {doc['spec']} | release_gate={doc['release_gate']} ===")
    print("hashes OK: " + ", ".join(f"{k} {v[:8]}" for k, v in doc["input_hashes"].items()) + f", v1.7 module {doc['v17_module']['sha256'][:8]}")
    print("cohort and landmark eligibility reproduce: " + ", ".join(f"{k}={v}" for k, v in doc["cohort_checks"]["eligible_at_landmark"].items()))
    for inv, d in doc["inventories"].items():
        items = d["items"].values()
        n10 = sum(1 for o in items if o["train"]["encounters"] != "<10")
        print(f"{inv:28s} items={len(d['items'])} with>=10 train encounters={n10}")
    s = doc["suppression"]
    print(f"suppression: primary={s['primary_cells']} complementary={s['complementary_cells']} recoverable={s['recoverable']}")
    print("shared: " + ", ".join([cfg.shared_json.name, cfg.shared_md.name] + list(TEMPLATE_FILES.values())))
    print("release gate ok")


def build_config(argv: list[str] | None = None) -> InvConfig:
    here = Path(__file__).resolve()
    default_root = Path(os.environ.get("FOUNDATIONS_ROOT", str(here.parents[1])))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", default=str(default_root))
    ap.add_argument("--raw-dir", default=os.environ.get("PIC_RAW_DIR"))
    ap.add_argument("--split-file", default=os.environ.get("SPLIT_FILE_PATH"))
    ap.add_argument("--canonical", default=os.environ.get("CANONICAL_PATH"))
    ap.add_argument("--harmonization", default=os.environ.get("HARMONIZATION_PATH"))
    ap.add_argument("--out-dir", default=os.environ.get("AUDIT_OUT_DIR"))
    ap.add_argument("--harm-t0-col", default=os.environ.get("HARM_T0_COL"))
    ap.add_argument("--harm-out-col", default=os.environ.get("HARM_OUT_COL"))
    ap.add_argument("--split-col", default=os.environ.get("SPLIT_COL", "split"))
    ap.add_argument("--v17-module", default=os.environ.get("V17_MODULE_PATH", str(here.parent / V17_MODULE_FILE)))
    ap.add_argument("--chunksize", type=int, default=CHUNKSIZE)
    a = ap.parse_args(argv)
    if not a.raw_dir or not a.split_file:
        raise AuditError("--raw-dir and --split-file are required (or set PIC_RAW_DIR and SPLIT_FILE_PATH); no local default path is shipped")
    root = Path(a.project_root)
    interim = root / "data" / "interim"
    return InvConfig(
        raw_dir=Path(a.raw_dir),
        canonical_path=Path(a.canonical) if a.canonical else interim / "canonical_icu_encounters.csv",
        harmonization_path=Path(a.harmonization) if a.harmonization else interim / "encounter_harmonization_fields.csv",
        split_path=Path(a.split_file), out_dir=Path(a.out_dir) if a.out_dir else interim,
        v17_module_path=Path(a.v17_module), harm_t0_col=a.harm_t0_col, harm_out_col=a.harm_out_col,
        split_col=a.split_col, chunksize=a.chunksize)


def _force_utf8_console() -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass


def main(argv: list[str] | None = None) -> int:
    _force_utf8_console()
    cfg = build_config(argv)
    try:
        run_inventory(cfg)
    except ReleaseGateError as e:
        print(f"RELEASE GATE FAILED - no shared file left in place: {e}", file=sys.stderr)
        return 3
    except AuditError as e:
        print(f"AUDIT STOPPED: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
