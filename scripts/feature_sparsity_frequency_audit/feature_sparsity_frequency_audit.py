#!/usr/bin/env python3
"""
Feature sparsity and measurement-frequency audit, Stage 2 (read-only, aggregate-only).
Implements `feature_sparsity_frequency_audit_spec.md` v0.8 (sections 5.2, 5.3, 6, 7, 8, 9, 11, 12; decisions D1-D20).

What it does, in order (any failure stops the run; nothing shared is left behind):
  1. Hash gates (exact 64-character SHA-256, checked BEFORE a file is read): the duration audit v1.7 module, the
     canonical encounters file, the harmonization file, the split file, the surgery complexity map, the four v1 maps,
     the neuro admitting-diagnosis set v1, and the raw PIC tables whose hashes Stage 1 recorded (D_ITEMS, D_LABITEMS,
     CHARTEVENTS, LABEVENTS, OUTPUTEVENTS, INPUTEVENTS, PRESCRIPTIONS, MICROBIOLOGYEVENTS). PATIENTS, ADMISSIONS,
     ICUSTAYS and SURGERY_INFO were not hashed by Stage 1: their SHA-256 is recorded and their row counts are gated.
  2. Cohort, split and landmark eligibility must reproduce the frozen pins (values omitted from this public copy)
     (OUT > t0 + L, strict) reproduced exactly; in-hospital death counts omitted from this public copy.
  3. Map checks (spec 5.3): duplicate-item merge distributions (quantiles, KS, same-draw differences), HR merge (1003 kept
     over 1002 at the same CHARTTIME), MAP from same-CHARTTIME SBP/DBP pairs, P/F from the same draw (pSOFA cut points,
     no support requirement, D19), S/F only where no P/F exists (SpO2 <= 97 %, 5-minute pairing, D15), unit conversions
     (5257, 6525 per row; magnesium 'ug/dl' rows held out and compared), urine/intake charting intervals and fluid
     intake mL/kg/day, MICROBIOLOGYEVENTS column fill rates, the prescription in-stay / date-only start check, and the
     osmotic-agent resolution by the ADMISSIONS admitting diagnosis (neuro set v1, Q28 by 4 characters: Q28.2/Q28.3).
  4. Metrics (spec 6-8): Layer A missingness / prevalence / level support; Layers B and C metrics 1-12 per channel,
     drug group, specimen group and Layer C surgery; strata (age band, surgical status, first_careunit); outcome
     descriptors (train only) ONLY with --outcome-descriptors (the flag's state is written to the output).
  5. Disclosure: counts 0-9 -> '<10'; complementary 'S' until the v1.7 LP recoverability solver finds nothing
     recoverable; percentages, quantiles and KS only where every underlying count is >= 10 and published; never
     min/max, IDs or timestamps. Shared files are written only after the disk-reload release gate passes.

Usage (Windows, from the project root; this script and audit_surgery_complexity_duration.py both in scripts/):
    python scripts/feature_sparsity_frequency_audit.py --raw-dir "C:/.../PICDataset" --split-file "C:/.../subject_splits.csv" ^
        --harm-t0-col t0 --harm-out-col OUTTIME --maps-dir data/interim/maps --outcome-descriptors

Outputs (data/interim/):
    feature_sparsity_frequency_audit.json / .md         shared (suppressed); written ONLY if the release gate passes
    feature_sparsity_frequency_audit_raw_LOCAL.json     local only - NEVER share, copy or upload (unsuppressed)
Exit codes: 0 ok; 2 input/hash/reconciliation/acceptance failure; 3 release gate failure (no shared file left).
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import stats as sps

# =============================================================================
# CONFIG BLOCK (spec v0.7, Nov 14, 2026)
# =============================================================================
SPEC_VERSION = "feature_sparsity_frequency_audit v1.0.1 (feature_sparsity_frequency_audit_spec v0.8, stage 2)"
V17_MODULE_FILE = "audit_surgery_complexity_duration.py"
V17_SPEC_VERSION = "surgery_complexity_duration_audit v1.7"
V17_MODULE_SHA256 = "1de2b17b2d894caff6a8362a9c72a458967691a8d2c7a533db863e87ad8dad94"  # public copy (cohort pins omitted); original private pin was ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0
V17_PINNED_SHA256 = {
    "canonical": "e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418",
    "harmonization": "ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d",
    "complexity_map": "5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5",
    "split": "60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31",
}
PINNED_SHA256 = dict(V17_PINNED_SHA256)
# v1 maps and the neuro set (decided by Vinai Nov 14, 2026, D1-D19); exact SHA-256 and row counts
MAP_FILES = {
    "channel": "feature_channel_map_v1.csv",
    "drug": "drug_group_map_v1.csv",
    "specimen": "specimen_map_v1.csv",
    "organism": "organism_map_v1.csv",
    "neuro": "neuro_admitting_dx_icd10_v1.csv",
}
MAP_SHA256 = {
    "channel": "56e8bf72fcbcaac765a39c2e705aedef61323e94841dbdcc8b2657191af370b6",
    "drug": "c557c3a12a9e1bbaf4e95c54b75affd4b3310cfd498cb6a485a13ae4e8ff37d7",
    "specimen": "caed429de4f2d484b40785d31ecbe3375f308db58a332d108c247b5ce331e3a5",
    "organism": "087e5ba3f7263b55e78b1fc08e41e6165bfa2f8282e94f10a249eca1eb6817ce",
    "neuro": "3cc8f4acfe9df29037b60d903ad4204c40db5bc76923de8b64234ccb2553b9d9",
}
MAP_ROWS = {"channel": 685, "drug": 343, "specimen": 71, "organism": 101, "neuro": 67}
MAP_COLUMNS = {
    "channel": ["table", "itemid", "label", "category", "units", "train_encounters", "train_rows", "channel",
                "table51_system", "triggered", "contract16", "unit_conversion", "plausibility_low", "plausibility_high",
                "source_itemid", "merge_rule", "stage2_check", "status", "confidence", "decision", "notes"],
    "drug": ["drug_name", "drug_name_en", "units", "train_encounters", "group", "sensitivity_subgroup",
             "conditional_rule", "subclass", "confidence", "decision", "notes"],
    "specimen": ["spec_itemid", "spec_type_desc", "train_encounters", "specimen_group", "base_code", "test_type",
                 "susceptibility_addons", "result_leak_risk", "confidence", "decision", "notes"],
    "organism": ["org_itemid", "org_name", "train_encounters", "canonical_name", "category", "contaminant_rule",
                 "confidence", "decision", "notes"],
    "neuro": ["icd10_3char", "label", "domain", "status", "match_rule", "reason"],
}
# Raw tables: SHA-256 recorded by Stage 1 (feature_item_inventory.json raw_tables) -> exact gate
RAW_STAGE1_SHA256 = {
    "D_ITEMS": "a7b65cf5d61ca84f6ce5d07d3775b01116473176a980463d4fe9df1dd12e1591",
    "D_LABITEMS": "d4557986fe9c5475d0db6506fdf12801ecfb0e2e7a3c01caf6ac86bbab69ed60",
    "CHARTEVENTS": "fe209d47d15825f180624f2cefb9bc38f8513e1ac0b3709d57aef37ae88f328e",
    "LABEVENTS": "784ecb929a358ee1a18fc9beafbdf67f30a91646b780d9983714a5e9ce0018c5",
    "OUTPUTEVENTS": "4c49f3c6910b4e8aaa817a13acf16d7120e5ecdd53cb6675bfc7760fda3d4d0d",
    "INPUTEVENTS": "7d8aaebf51d82fdde32841a534f3d296911e92c3421885e40b379bc7c394f9c6",
    "PRESCRIPTIONS": "61fa43dcb1ed977f4219c902b01977da2d898a16a2c33640997ad8337d2eb2c6",
    "MICROBIOLOGYEVENTS": "e19cf9d41732f74911039acae6b2cd7bf5c5dd65179cbd31aa8634dd1e6bb45c",
}
EXPECTED_RAW_ROWS = {}  # public copy: table-size pins omitted; see configs/cohort_pins.example.json
# Not hashed by Stage 1: SHA-256 recorded (not gated); row count gated (open question for Vinai)
RAW_RECORD_ONLY = ("PATIENTS", "ADMISSIONS", "ICUSTAYS", "SURGERY_INFO")

MIN_CELL = 10
SPLITS = ("train", "validation", "test")
LANDMARKS = (24, 48, 72, 120, 168, 240, 336)
EXPECTED_ELIGIBLE = {}  # public copy: landmark eligibility pins omitted
EXPECTED_DIED = None  # public copy: pin omitted
MAX_H = 336.0
BIN_H = 6.0
N_BINS = 56
LAYER_B_H = 24.0
PRE_ICU_H = (6.0, 24.0)               # [t0 - 6 h, t0) primary, [t0 - 24 h, t0) sensitivity
LAB_RELEASE_H = 1.0                   # decision 11 / A2.3
LAB_OFFSETS = (0.0, 1.0, 2.0)         # availability sensitivities (labs only)
KEEP_LO_H = -26.0                     # earliest chart time kept (pre-ICU 24 h at the +0/+2 h offsets)
MED_KEEP_LO_H = -72.0                 # prescriptions kept from t0 - 72 h (courses that are active at t0)
DETAIL_LANDMARKS = (24, 72, 168)      # per-bin (d) availability and prefix-bin summaries (spec 9.1 landmarks)
STRATA_LANDMARK = 72
DETAIL_MIN_TRAIN = 100                # detailed outputs only for channels with >= 100 train encounters (9.1 drop rule)
SF_PAIR_MINUTES = 5.0
SF_SPO2_MAX = 97.0
VASO_CUTOFF_H = 12.0
DEFAULT_CUTOFF_H = 24.0
CUTOFF_FACTORS = (0.5, 1.0, 2.0)
ENDDATE_TOL_H = 6.0
MICRO_TA_H = {"gram_class": 24.0, "identification": 48.0, "susceptibility": 72.0, "negative_final": 120.0}
MICRO_GRID = (0.5, 1.0, 1.5)
ADM_SOURCE_X = (6, 12, 24, 48)
LEVEL_M = (20, 10)
AGE_BANDS = ("neonate", "infant", "toddler", "preschool", "school", "adolescent", "ge18", "age_invalid")
DOB_STATUS = ("missing", "unparseable", "negative", "valid")
DOB_NEG_BINS = ("lt1d", "1to30d", "gt30d")
MAX_INVALID_AGE_FRAC = 0.01  # D20 safety gate: 1% of the cohort
WEIGHT_SOURCES = ("charted", "imputed_age_sex", "imputed_age_band", "imputed_sex_only", "imputed_overall")
WEIGHT_IMPUTE_RULE = (
    "charted weight 1014 if present; else for a valid age band, the train median by age band x sex "
    "(n>=10), else by age band (n>=10); age_invalid is never imputed from an age band and falls back "
    "to the train median by sex (n>=10), else the overall train median (D20, A2.12)"
)
OUTCOME_GROUPS = ("0-2", "3-5", "6-7", "8")
OUTCOME_EDGES = (24, 48, 72, 120, 168, 240, 336)
MED_GROUPS = ("vasoactive", "antibiotic", "antifungal", "antiviral", "sedative_analgesic", "non_opioid_analgesic",
              "antiemetic", "diuretic", "osmotherapy", "electrolyte_replacement", "steroid", "nmb", "aed")
DRUG_MAP_GROUPS = set(MED_GROUPS) - {"osmotherapy"} | {"none", "osmotic_conditional"}
SPECIMEN_GROUPS_V1 = ("blood", "urine", "respiratory", "csf", "stool", "ureaplasma", "culture_site_unknown", "other")
# LOCF families (plan v12 Layer C table; contract 10.2-10.3); cap used for time since last (spec 9.1)
LOCF_4BIN = {"creatinine", "bun", "hemoglobin", "platelets", "wbc", "alt", "ast", "bilirubin_total", "bilirubin_direct",
             "albumin", "ggt", "inr", "ptt", "fibrinogen", "d_dimer", "anc", "lymphocytes_abs", "cystatin_c",
             "amylase", "lipase", "ldh"}
LOCF_1BIN = {"ph", "lactate", "bicarbonate", "sodium", "potassium", "glucose", "chloride", "calcium_total", "magnesium",
             "phosphate", "anion_gap", "pco2", "po2", "base_excess", "sao2_bloodgas"}
VITALS = ("heart_rate", "respiratory_rate", "map", "spo2", "temperature")
SCORE_INPUTS = {"pf_ratio", "sf_ratio", "platelets", "bilirubin_total", "map", "creatinine", "lactate",  # pSOFA
                "heart_rate", "sbp", "temperature", "ph", "bicarbonate", "pco2", "po2", "glucose", "potassium",
                "bun", "wbc", "inr", "ptt"}  # PRISM III
# Per-row and fixed unit conversions: itemid -> (exact map text, kind, args). Any other conversion text stops the run.
CONVERSIONS = {
    "1011": ("mmol/L -> mg/dL x18.016", "factor", 18.016, "mg/dL"),
    "5024": ("g/L -> g/dL /10", "factor", 0.1, "g/dL"),
    "5032": ("umol/L -> mg/dL /88.4", "factor", 1 / 88.4, "mg/dL"),
    "5041": ("umol/L -> mg/dL /88.4", "factor", 1 / 88.4, "mg/dL"),
    "5033": ("urea mmol/L -> BUN mg/dL x2.801", "factor", 2.801, "mg/dL"),
    "6477": ("urea mmol/L -> BUN mg/dL x2.801", "factor", 2.801, "mg/dL"),
    "5034": ("mmol/L -> mg/dL x4.008", "factor", 4.008, "mg/dL"),
    "5215": ("mmol/L -> mg/dL x4.008", "factor", 4.008, "mg/dL"),
    "5042": ("umol/L -> mg/dL /17.1", "factor", 1 / 17.1, "mg/dL"),
    "5075": ("umol/L -> mg/dL /17.1", "factor", 1 / 17.1, "mg/dL"),
    "5255": ("umol/L -> mg/dL /17.1", "factor", 1 / 17.1, "mg/dL"),
    "5047": ("mmol/L -> mg/dL x18.016", "factor", 18.016, "mg/dL"),
    "5223": ("mmol/L -> mg/dL x18.016", "factor", 18.016, "mg/dL"),
    "5072": ("mmol/L -> mg/dL x3.097", "factor", 3.097, "mg/dL"),
    "5099": ("g/L /10 -> g/dL", "factor", 0.1, "g/dL"),
    "5222": ("% -> fraction /100", "factor", 0.01, "fraction"),
    "5235": ("kPa -> mmHg x7.5 if any", "per_row_gate", {"kpa": 7.5, "mmhg": 1.0, "": 1.0}, "mmHg"),
    "5257": ("per-row VALUEUOM: g/L /10; g/dL x1 -> g/dL", "per_row_gate", {"g/l": 1.0, "g/dl": 10.0, "": 1.0}, "g/L"),
    "6525": ("per-row VALUEUOM: ng/mL x1; pg/mL /1000 -> ng/mL", "per_row_gate", {"ng/ml": 1.0, "pg/ml": 0.001, "": 1.0}, "ng/mL"),
    "6095": ("mIU/L = uU/mL", "factor", 1.0, "mIU/L"),
    "6156": ("ug/dL -> nmol/L x27.59", "factor", 27.59, "nmol/L"),
    "5059": ("mmol/L as is; 'ug/dl' rows: Stage 2 distribution check - convert if consistent after scaling, else drop those rows",
             "per_row_gate", {"mmol/l": 1.0, "": 1.0, "ug/dl": None}, "mmol/L"),
    "1008": ("aggregation chosen after the Stage 2 charting-interval check: daily total or mL/kg/h over the charting interval (imputed weight, A2.12)", "aggregation", None, "mL"),
    "1009": ("mL/kg/day: daily total of charted intake / weight (charted weight 1014, else imputed from age and gender, A2.12)", "aggregation", None, "mL"),
}
MG_SCALINGS = {"as_recorded": 1.0, "as_mg_per_dl": 0.4114, "as_ug_per_dl": 0.0004114}  # -> mmol/L
HR_KEEP, HR_DROP = "1003", "1002"
SBP_ITEM, DBP_ITEM, SPO2_ITEM, FIO2_ITEM, WEIGHT_ITEM = "1016", "1015", "1006", "5222", "1014"
URINE_ITEM, INTAKE_ITEM = "1008", "1009"
PO2_ITEMS, PF_LAB_ITEM = ("5239", "5244"), "5243"
CALCIUM_ITEMS, IONIZED_RANGE = ("5215", "5034"), (0.9, 1.6)
PSOFA_PF = ((400, 0), (300, 1), (200, 2), (100, 3), (-math.inf, 4))   # D19: >=400 0, 300-399 1, 200-299 2, 100-199 3, <100 4
PSOFA_SF = ((292, 0), (264, 1), (221, 2), (148, 3), (-math.inf, 4))   # D15: >=292 0, 264-291 1, 221-263 2, 148-220 3, <148 4
CHUNKSIZE = 1_000_000
BLANK = "(blank)"
NON_COUNT_TOP = ("disclosure", "suppression", "input_hashes", "v17_module", "raw_tables", "maps", "cohort_checks",
                 "config", "outcome_descriptors_flag", "spec", "release_gate", "provisional_classes", "notes")
COUNT_KEY_RE = re.compile(r"^n(_|$)")
TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}")
QKEYS = ("p10", "p25", "p50", "p75", "p90")


class AuditError(Exception):
    """Input, hash, reconciliation or acceptance failure (exit 2)."""


class HashMismatchError(AuditError):
    pass


class ReleaseGateError(Exception):
    """Release gate failure (exit 3)."""


# =============================================================================
# Hashes and the pinned v1.7 module
# =============================================================================
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


# =============================================================================
# Config
# =============================================================================
@dataclass
class S2Config:
    raw_dir: Path
    canonical_path: Path
    harmonization_path: Path
    split_path: Path
    complexity_map_path: Path
    maps_dir: Path
    out_dir: Path
    v17_module_path: Path
    outcome_descriptors: bool = False
    pinned_sha256: dict = field(default_factory=lambda: dict(PINNED_SHA256))
    map_sha256: dict = field(default_factory=lambda: dict(MAP_SHA256))
    map_rows: dict = field(default_factory=lambda: dict(MAP_ROWS))
    raw_sha256: dict = field(default_factory=lambda: dict(RAW_STAGE1_SHA256))
    v17_module_sha256: str = V17_MODULE_SHA256
    harm_t0_col: str | None = None
    harm_out_col: str | None = None
    split_subject_col: str = "SUBJECT_ID"
    split_col: str = "split"
    expected_encounters: int | None = None  # public copy: pin omitted
    expected_subjects: int | None = None
    expected_split_encounters: dict = field(default_factory=dict)
    expected_eligible: dict = field(default_factory=lambda: dict(EXPECTED_ELIGIBLE))
    expected_died: int | None = EXPECTED_DIED
    expected_raw_rows: dict = field(default_factory=lambda: dict(EXPECTED_RAW_ROWS))
    expected_surgery_rows: int | None = None  # public copy: pin omitted
    expected_complexity_map_rows: int = 1495
    surgery_expected_missing: dict | None = None
    chunksize: int = CHUNKSIZE
    verbose: bool = True

    @property
    def shared_json(self) -> Path:
        return self.out_dir / "feature_sparsity_frequency_audit.json"

    @property
    def shared_md(self) -> Path:
        return self.out_dir / "feature_sparsity_frequency_audit.md"

    @property
    def local_json(self) -> Path:
        return self.out_dir / "feature_sparsity_frequency_audit_raw_LOCAL.json"


def log(cfg: S2Config, msg: str) -> None:
    if cfg.verbose:
        print(msg, flush=True)


# =============================================================================
# Pure helpers (unit tested)
# =============================================================================
_NULLS = {"", "nan", "none", "null", "nat"}


def vnorm_id(s: pd.Series) -> pd.Series:
    """Vectorized v1.7 norm_id for CSV text: trims, drops quotes, '123.0' -> '123'; other text unchanged."""
    s = s.astype(str).str.strip().str.strip("'\"").str.strip()
    s = s.where(~s.str.lower().isin(_NULLS), "")
    m = s.str.fullmatch(r"[+-]?\d+\.0*")
    if m.any():
        u = {v: str(int(v.split(".")[0])) for v in pd.unique(s[m])}
        s = s.where(~m, s.map(u))
    return s


def vclean_text(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
    return s.where(~s.str.lower().isin(_NULLS), "")


def fast_time(s: pd.Series) -> pd.Series:
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


def has_clock(s: pd.Series) -> pd.Series:
    """True where the raw time string carries a clock time (HH:MM); False for date-only strings."""
    return s.astype(str).str.contains(r"\d{1,2}:\d{2}", regex=True)


def norm_unit(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.lower().str.replace("μ", "u").str.replace("µ", "u").str.replace(" ", "")
    return s.where(~s.isin(_NULLS), "")


def resolve(header: list[str], candidates: Iterable[str], label: str, required: bool = True) -> str | None:
    up = {h.upper(): h for h in header}
    for c in candidates:
        if c.upper() in up:
            return up[c.upper()]
    if required:
        raise AuditError(f"{label}: none of the columns {list(candidates)} found; header has {header}")
    return None


def hours(td: pd.Series | np.ndarray) -> np.ndarray:
    return (pd.to_timedelta(pd.Series(td)).dt.total_seconds() / 3600.0).to_numpy(dtype=float)


def r1(x: float) -> float:
    return float(round(float(x), 1))


def qdict(vals: np.ndarray) -> dict | None:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < MIN_CELL:
        return None
    p = np.percentile(vals, [10, 25, 50, 75, 90])
    return {k: r1(v) for k, v in zip(QKEYS, p)}


def psofa_score(v: np.ndarray, cuts: tuple) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    out = np.full(len(v), -1, dtype=int)
    for lo, sc in cuts:
        m = (out < 0) & (v >= lo)
        out[m] = sc
    return out


def norm_icd(code: Any) -> tuple[str, str]:
    """Return (normalized code without dots/spaces, status) where status is ok / missing / invalid / multiple.
    Several codes in one field: the first is used and the row is counted as 'multiple' (ambiguous)."""
    s = str(code if code is not None else "").strip().upper()
    if s.lower() in _NULLS:
        return "", "missing"
    parts = [p for p in re.split(r"[,;/|，；、\s]+", s) if p]
    first = re.sub(r"[^A-Z0-9]", "", parts[0]) if parts else ""
    if not re.match(r"^[A-Z][0-9]{2}", first):
        return first, "invalid"
    return first, ("multiple" if len(parts) > 1 else "ok")


def icd_domain(c3: str) -> str:
    """Draft Appendix A (plan v12) diagnosis-only domain with the draft defaults: chapter I, D00-D49, V-Y and Q00-Q07
    stay other_medical (open points for Vinai). Missing or invalid -> other_medical."""
    if not re.fullmatch(r"[A-Z][0-9]{2}", c3 or ""):
        return "other_medical"
    L, n = c3[0], int(c3[1:])
    if L == "Q" and 20 <= n <= 28:
        return "cardiac_congenital_heart"
    if L == "P":
        return "neonatal_perinatal"
    if L in ("S", "T"):
        return "surgical_trauma"
    if L in ("J", "A", "B"):
        return "respiratory_infection"
    if L == "G":
        return "neurological"
    if L == "C" or (L == "D" and 50 <= n <= 89):
        return "oncology_hematology"
    if L in ("K", "N"):
        return "gi_renal"
    return "other_medical"


def neuro_matcher(neuro: pd.DataFrame):
    """Neuro set v1: status include/partial; '3-character prefix' rows match the 3-character category; the Q28 row
    ('4-character prefix Q28.2 or Q28.3 only') matches only normalized codes starting Q282 or Q283."""
    three, four = set(), set()
    for _, r in neuro.iterrows():
        if r["status"] not in ("include", "partial"):
            continue
        rule = r["match_rule"].strip()
        c3 = r["icd10_3char"].strip().upper()
        if rule == "3-character prefix":
            three.add(c3)
        elif rule.startswith("4-character prefix"):
            codes = re.findall(r"([A-Z][0-9]{2})\.([0-9])", rule)
            if not codes:
                raise AuditError(f"neuro set: cannot parse the 4-character rule for {c3}: {rule!r}")
            four |= {a + b for a, b in codes}
        else:
            raise AuditError(f"neuro set: unknown match_rule {rule!r} for an included row ({c3})")

    def match(norm_code: str) -> bool:
        if not norm_code or not re.match(r"^[A-Z][0-9]{2}", norm_code):
            return False
        if norm_code[:3] in three:
            return True
        return any(norm_code.startswith(f) for f in four)
    match.three, match.four = three, four
    return match


def outcome_group(died: np.ndarray, out_h: np.ndarray, L: float) -> np.ndarray:
    """Outcome 1 group at landmark L (amendment A1.2): 8 = died in hospital; else the remaining ICU time
    R = OUT - (t0 + L) -> level 0 (<= 24 h) ... 7 (> 336 h); groups 0-2 / 3-5 / 6-7 / 8."""
    R = out_h - L
    lvl = np.searchsorted(np.array(OUTCOME_EDGES, dtype=float), R, side="left")  # R <= 24 -> 0; (24,48] -> 1 ...
    lvl = np.where(died, 8, lvl)
    grp = np.where(lvl <= 2, "0-2", np.where(lvl <= 5, "3-5", np.where(lvl <= 7, "6-7", "8")))
    return grp


def age_band(days: np.ndarray) -> np.ndarray:
    days = np.asarray(days, dtype=float)
    y = days / 365.25
    band = np.select([days < 28, y < 1, y < 3, y < 6, y < 12, y < 18],
                     ["neonate", "infant", "toddler", "preschool", "school", "adolescent"], "ge18")
    return np.where(np.isfinite(days) & (days >= 0), band, "age_invalid")


def classify_ages(t0: pd.Series, dob_parsed: pd.Series, dob_raw: pd.Series, has_row: np.ndarray
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Classify each encounter's DOB/age (D20). age_days is NaN unless the age is valid and >= 0.

    Status: missing (null/blank DOB, or no PATIENTS row); unparseable (DOB present but not a timestamp,
    or a parsed DOB whose age cannot be computed); negative (t0 < DOB); valid. Negative magnitude bins
    are lt1d (|age| < 1 day), 1to30d (1-30 days), gt30d (> 30 days); other rows have an empty bin.
    """
    raw = vclean_text(dob_raw.fillna("")).to_numpy()
    has = np.asarray(has_row, dtype=bool)
    missing = (~has) | (raw == "")
    parsed = dob_parsed
    unparseable = (~missing) & parsed.isna().to_numpy()
    age = (pd.to_datetime(t0) - parsed).dt.total_seconds().to_numpy() / 86400.0
    finite = np.isfinite(age)
    unparseable = unparseable | ((~missing) & (~unparseable) & (~finite))
    negative = (~missing) & (~unparseable) & finite & (age < 0)
    valid = (~missing) & (~unparseable) & (~negative)
    status = np.select([missing, unparseable, negative],
                       ["missing", "unparseable", "negative"], "valid")
    age_out = np.where(valid, age, np.nan)
    mag = np.where(negative, -age, np.nan)
    nbin = np.select([negative & (mag < 1.0), negative & (mag <= 30.0), negative],
                     ["lt1d", "1to30d", "gt30d"], "")
    return status, age_out, nbin


def invalid_age_limit(n_enc: int) -> int:
    """Largest allowed invalid-age count: 1% of the cohort."""
    return int(MAX_INVALID_AGE_FRAC * int(n_enc))


def check_invalid_age_rate(n_invalid: int, n_enc: int) -> None:
    limit = invalid_age_limit(n_enc)
    if n_invalid > limit:
        pct = 0.0 if n_enc == 0 else 100.0 * n_invalid / n_enc
        raise AuditError(
            f"invalid DOB/age for {int(n_invalid)} of {int(n_enc)} cohort encounters "
            f"({pct:.2f}% > {100 * MAX_INVALID_AGE_FRAC:.0f}% safety gate, {limit} encounters, D20)"
        )


def shift_of(hour: np.ndarray) -> np.ndarray:
    return np.where((hour >= 8) & (hour < 15), "day", np.where((hour >= 15) & (hour < 23), "evening", "night"))


def interval_band(h: np.ndarray) -> np.ndarray:
    return np.select([np.isnan(h), h < 24, h < 72, h < 168], ["missing", "lt24h", "24to72h", "72to168h"], "ge168h")


# =============================================================================
# Output builder on top of the v1.7 Registry (cells, identities, quantiles, ratios, KS)
# =============================================================================
class Out:
    """Every count is a registry cell whose JSON key starts with 'n_' (or is 'n'). Quantiles live at 'q_<x>' next to
    their own count 'n_<x>'; percentages at 'pct_<x>' from two cells; KS at 'ks_<x>' from the cells it needs.
    Non-count values (labels, booleans) never use those prefixes."""

    def __init__(self, A):
        self.A = A
        self.reg = A.Registry()
        self.values: dict[str, Any] = {}
        self.quants: list[tuple[str, str, dict | None]] = []
        self.ratios: list[tuple[str, str, str, float]] = []
        self.kss: list[tuple[str, list[str], float | None]] = []
        self.frozen: set[str] = set()  # already-public totals (cohort and split sizes): never complementary-suppressed

    def cell(self, path: tuple, value: int, parent: str | None = None, tag: str = "sub") -> str:
        if not COUNT_KEY_RE.match(str(path[-1])):
            raise AuditError(f"count key without the n_ prefix: {path}")
        cid = self.A.P(*path)
        self.reg.count(cid, int(value))
        if parent is not None:
            self.reg.subset(cid, parent, tag)
        return cid

    def ident(self, total: str, parts: Iterable[str]) -> None:
        self.reg.identity(total, list(parts))

    def partition(self, base: tuple, counts: dict[str, int], total: str | None = None, total_value: int | None = None) -> dict[str, str]:
        """Cells base/<n_level> for every level, plus an identity total = sum(levels)."""
        cids = {k: self.cell(base + (f"n_{k}",), v) for k, v in counts.items()}
        if total is None:
            total = self.cell(base + ("n_total",), sum(counts.values()) if total_value is None else total_value)
        self.ident(total, cids.values())
        return cids

    def q(self, base: tuple, name: str, vals: np.ndarray, parent: str | None = None, equal_to: str | None = None) -> str:
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        n = self.cell(base + (f"n_{name}",), len(vals), parent=parent)
        if equal_to is not None:
            self.reg.equal(n, equal_to)
        self.quants.append((self.A.P(*base, f"q_{name}"), n, qdict(vals)))
        return n

    def pct(self, base: tuple, name: str, num: str, den: str) -> None:
        self.ratios.append((self.A.P(*base, f"pct_{name}"), num, den, 100.0))

    def ks(self, base: tuple, name: str, needs: list[str], value: float | None) -> None:
        self.kss.append((self.A.P(*base, f"ks_{name}"), list(needs), value))

    def dep(self, base: tuple, name: str, needs: list[str], value: Any) -> None:
        """A derived value (share of hours, flag) published only if every count it rests on is published."""
        self.kss.append((self.A.P(*base, f"x_{name}"), list(needs), value))

    def val(self, path: tuple, v: Any) -> None:
        k = str(path[-1])
        if COUNT_KEY_RE.match(k) or k.startswith(("q_", "pct_", "ks_", "x_")):
            raise AuditError(f"non-count value under a reserved key: {path}")
        self.values[self.A.P(*path)] = v


def apply_suppression_frozen(A, reg, frozen: set[str], max_rounds: int = 1000) -> tuple[dict[str, str], dict]:
    """v1.7 apply_suppression (same primary rule, same greedy order, same final full-LP round) with one change:
    cells in `frozen` (cohort and split sizes that are already public) are never chosen for complementary
    suppression; they stay known to the recoverability solver, which is the conservative direction."""
    status = {c: ("P" if reg.raw[c] < A.MIN_CELL else "pub") for c in reg.published()}
    bad = [c for c in frozen if status.get(c) != "pub"]
    if bad:
        raise ReleaseGateError(f"a frozen (already public) total is under {A.MIN_CELL}: {bad[:3]}")
    by_cell: dict[str, list[int]] = {}
    for k, (t, parts) in enumerate(reg.identities):
        for c in (t,) + parts:
            by_cell.setdefault(c, []).append(k)
    cache: dict = {}
    proven: set[str] = set()
    pickable = lambda m: status.get(m) == "pub" and m not in frozen  # noqa: E731
    fast = True
    for _ in range(max_rounds):
        known = {c: reg.raw[c] for c, st in status.items() if st == "pub"}
        unknown = {c: st for c, st in status.items() if st != "pub"} | {h: "H" for h in reg.hidden}
        checked, unsafe, _ = A.recoverability_check(reg.identities, known, unknown, cache, propagation_only=fast,
                                                    proven_safe=proven)
        if not unsafe and fast:
            fast = False
            continue
        if not unsafe:
            return status, {"primary_cells": sum(1 for v in status.values() if v == "P"),
                            "complementary_cells": sum(1 for v in status.values() if v == "S"),
                            "recoverability_cells_checked": checked, "recoverable": 0}
        touched: set[int] = set()
        progressed = False
        for u in sorted(unsafe):
            if any(k in touched for k in by_cell.get(u, [])):
                continue
            best = None
            for k in by_cell.get(u, []):
                t, parts = reg.identities[k]
                members = (t,) + parts
                pubs = [m for m in members if pickable(m)]
                if not pubs:
                    continue
                n_unknown = sum(1 for m in members if status.get(m) != "pub")
                cand = min(pubs, key=lambda m: (len(by_cell.get(m, [])), reg.raw[m], m))
                score = (n_unknown, len(by_cell.get(cand, [])), reg.raw[cand], cand)
                if best is None or score < best[0]:
                    best = (score, k, cand)
            if best is None:  # breadth-first walk through identities linked by unknown cells
                seen_c, seen_k, frontier = {u}, set(), [u]
                while frontier and best is None:
                    cands, nxt = [], []
                    for c in frontier:
                        for k in by_cell.get(c, []):
                            if k in seen_k:
                                continue
                            seen_k.add(k)
                            t, parts = reg.identities[k]
                            for m in (t,) + parts:
                                if pickable(m):
                                    cands.append((reg.raw[m], m, k))
                                elif m not in seen_c and status.get(m) != "pub":
                                    seen_c.add(m)
                                    nxt.append(m)
                    if cands:
                        v, m, k = min(cands)
                        best = ((0, v, m), k, m)
                    frontier = nxt
            if best is None:
                continue
            status[best[2]] = "S"
            touched.add(best[1])
            progressed = True
        fast = True
        if not progressed:
            raise ReleaseGateError(f"cannot protect {len(unsafe)} cell(s): no publishable member left to suppress "
                                   f"(cell ids: {sorted(unsafe)[:5]})")
    raise ReleaseGateError("complementary suppression did not converge")


def fill_document(out: Out, status: dict, doc: dict) -> None:
    A = out.A
    for c in out.reg.published():
        st = status[c]
        A.set_ptr(doc, c, out.reg.raw[c] if st == "pub" else ("<10" if st == "P" else "S"))
    pub = {c for c, st in status.items() if st == "pub"}
    for ptr, n, qd in out.quants:
        A.set_ptr(doc, ptr, qd if (n in pub and out.reg.raw[n] >= MIN_CELL and qd is not None) else "suppressed")
    for ptr, num, den, scale in out.ratios:
        ok = num in pub and den in pub and out.reg.raw[den] > 0
        A.set_ptr(doc, ptr, r1(scale * out.reg.raw[num] / out.reg.raw[den]) if ok else "suppressed")
    for ptr, needs, v in out.kss:
        ok = v is not None and all(n in pub for n in needs)
        if ok and isinstance(v, float):
            v = float(round(v, 3))
        A.set_ptr(doc, ptr, v if ok else "suppressed")
    for ptr, v in out.values.items():
        A.set_ptr(doc, ptr, v)


# =============================================================================
# Inputs: cohort, split, landmark eligibility, harmonization extras
# =============================================================================
def load_cohort(cfg: S2Config, A) -> tuple[pd.DataFrame, dict]:
    hashes = {k: verify_hash(k, p, cfg.pinned_sha256[k]) for k, p in
              (("canonical", cfg.canonical_path), ("harmonization", cfg.harmonization_path),
               ("split", cfg.split_path), ("complexity_map", cfg.complexity_map_path))}
    log(cfg, "input hash checks: canonical OK, harmonization OK, split OK, complexity_map OK (exact 64-char)")
    checks: dict[str, Any] = {}
    ch = A.read_header(cfg.canonical_path)
    need = ["encounter_id", "SUBJECT_ID", "HADM_ID", "hadm_id_list"]
    miss = [c for c in need if c not in ch]
    if miss:
        raise AuditError(f"canonical CSV missing columns: {miss}")
    died_c = resolve(ch, ("Died",), "canonical Died", required=False)
    canon = A.read_csv_str(cfg.canonical_path, usecols=need + ([died_c] if died_c else []))
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
    if died_c:
        dv = canon[died_c].map(lambda v: str(v).strip().lower())
        bad = ~dv.isin({"0", "1", "0.0", "1.0", "true", "false"})
        if bad.any():
            raise AuditError(f"canonical Died has {int(bad.sum())} values that are not 0/1")
        canon["died"] = dv.isin({"1", "1.0", "true"})
    else:
        canon["died"] = False
    checks["died_column_present"] = bool(died_c)

    hh = A.read_header(cfg.harmonization_path)
    if "encounter_id" not in hh:
        raise AuditError("harmonization file has no encounter_id column")
    try:
        t0c = A.resolve_column(hh, cfg.harm_t0_col, A.HARM_T0_CANDIDATES, "T0", A.HARM_T0_FALLBACK)
        outc = A.resolve_column(hh, cfg.harm_out_col, A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK)
    except A.AuditError as e:
        raise AuditError(str(e)) from e
    extras = {"first_careunit": resolve(hh, ("first_careunit", "FIRST_CAREUNIT"), "h", required=False),
              "gap_count": resolve(hh, ("gap_count_between_stays", "n_gaps"), "h", required=False),
              "prior": resolve(hh, ("prior_encounters", "prior_encounters_n"), "h", required=False)}
    cols = ["encounter_id", t0c, outc] + [c for c in extras.values() if c and c not in (t0c, outc)]
    harm = A.read_csv_str(cfg.harmonization_path, usecols=list(dict.fromkeys(cols)))
    harm["encounter_id"] = harm["encounter_id"].map(lambda v: str(v).strip())
    if harm["encounter_id"].duplicated().any():
        raise AuditError("harmonization encounter_id not unique")
    if set(harm["encounter_id"]) != set(canon["encounter_id"]):
        raise AuditError("harmonization encounter_id set differs from canonical encounter_id set")
    hdf = pd.DataFrame({"encounter_id": harm["encounter_id"], "t0": fast_time(harm[t0c]), "OUT": fast_time(harm[outc])})
    for k, c in extras.items():
        hdf[f"h_{k}"] = harm[c].map(lambda v: str(v).strip()) if c else ""
    if hdf["t0"].isna().any() or hdf["OUT"].isna().any():
        raise AuditError("harmonized t0 or OUT missing/unparseable for some encounters")
    if (hdf["OUT"] < hdf["t0"]).any():
        raise AuditError("harmonized OUT earlier than t0 for some encounters")
    checks["harmonization_columns"] = {"t0": t0c, "OUT": outc, **{k: v for k, v in extras.items()}}
    enc = canon.merge(hdf, on="encounter_id", how="left")

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
    elig = {L: int((enc["OUT"] > enc["t0"] + pd.Timedelta(hours=L)).sum()) for L in LANDMARKS}
    exp = {int(k): int(v) for k, v in cfg.expected_eligible.items()}
    if elig != exp:
        raise AuditError(f"core-landmark eligible counts do not reproduce: got {elig}, expected {exp}")
    if cfg.expected_died is not None and died_c and int(enc["died"].sum()) != cfg.expected_died:
        raise AuditError(f"canonical Died total {int(enc['died'].sum())} != expected {cfg.expected_died}")
    if cfg.outcome_descriptors and not died_c:
        raise AuditError("--outcome-descriptors needs the canonical Died column, which is absent")
    checks.update({"encounters": n_enc, "subjects": n_subj, "split_encounters": sc,
                   "eligible_at_landmark": {f"{L}h": elig[L] for L in LANDMARKS},
                   "died_total_matches_expected": bool(died_c) and cfg.expected_died is not None})
    log(cfg, f"cohort ok ({n_enc} encounters, {n_subj} subjects, splits {sc}); landmark eligibility reproduces")
    enc = enc.sort_values("encounter_id", kind="mergesort").reset_index(drop=True)
    enc["enc_i"] = np.arange(len(enc))
    enc["out_h"] = hours(enc["OUT"] - enc["t0"])
    enc["W"] = np.minimum(enc["out_h"].to_numpy(), MAX_H)
    return enc, {"hashes": hashes, "checks": checks}


def encounter_hadm_index(enc: pd.DataFrame) -> pd.DataFrame:
    rows = [(f"{s}|{h}", i) for s, lst, i in zip(enc["SUBJECT_ID"], enc["hadm_list"], enc["enc_i"]) for h in lst]
    eh = pd.DataFrame(rows, columns=["key", "enc_i"]).drop_duplicates()
    return eh.merge(enc[["enc_i", "t0", "OUT", "split", "W"]], on="enc_i", how="left")


@dataclass
class Cohort:
    enc: pd.DataFrame
    eh: pd.DataFrame
    n: int
    split: np.ndarray
    out_h: np.ndarray
    W: np.ndarray
    nopen: np.ndarray
    gap_bins: np.ndarray
    gaps: dict
    elig: dict
    masks: dict

    @staticmethod
    def build(enc: pd.DataFrame, gaps: dict[int, list[tuple[float, float]]]) -> "Cohort":
        n = len(enc)
        out_h = enc["out_h"].to_numpy(dtype=float)
        W = enc["W"].to_numpy(dtype=float)
        nopen = np.minimum(N_BINS, np.ceil(out_h / BIN_H - 1e-12)).astype(int)
        nopen = np.maximum(nopen, 0)
        gb = np.zeros((n, N_BINS), dtype=bool)
        for i, gl in gaps.items():
            for a, b in gl:
                for k in range(N_BINS):
                    if a <= BIN_H * k and BIN_H * (k + 1) <= b:
                        gb[i, k] = True
        split = enc["split"].to_numpy()
        return Cohort(enc=enc, eh=encounter_hadm_index(enc), n=n, split=split, out_h=out_h, W=W, nopen=nopen,
                      gap_bins=gb, gaps=gaps, elig={L: out_h > L for L in LANDMARKS},
                      masks={s: split == s for s in SPLITS})


# =============================================================================
# Maps (v1) and the neuro set
# =============================================================================
def load_maps(cfg: S2Config, A) -> tuple[dict, dict]:
    maps, hashes = {}, {}
    for k, fn in MAP_FILES.items():
        p = Path(cfg.maps_dir) / fn
        hashes[k] = verify_hash(f"map {fn}", p, cfg.map_sha256[k])
        hdr = A.read_header(p)
        if hdr != MAP_COLUMNS[k]:
            raise AuditError(f"map {fn}: columns differ from the v1 schema: {hdr}")
        d = A.read_csv_str(p)
        d = d.apply(lambda s: s.map(lambda v: str(v).strip()))
        if len(d) != cfg.map_rows[k]:
            raise AuditError(f"map {fn}: {len(d)} rows, expected {cfg.map_rows[k]}")
        maps[k] = d
    log(cfg, "map hash checks: channel, drug, specimen, organism, neuro OK (exact 64-char)")
    return maps, hashes


def parse_channel_map(cm: pd.DataFrame) -> dict:
    """Channel definitions from the v1 channel map. Returns {'items': {(table, itemid): info}, 'channels': {...},
    'merge_main': {channel: main itemid}}. Unknown conversion text, gates that do not parse, or a merge group without
    a main item stop the run."""
    use = cm[cm["status"].isin(["assigned", "assigned_merge_pending_stage2", "derived"])].copy()
    items, chans, main, dgates = {}, {}, {}, {}
    for _, r in use.iterrows():
        ch = r["channel"]
        if not ch:
            raise AuditError(f"channel map row {r['table']}/{r['itemid']} is {r['status']} but has no channel")
        lo = float(r["plausibility_low"]) if r["plausibility_low"] else None
        hi = float(r["plausibility_high"]) if r["plausibility_high"] else None
        conv = r["unit_conversion"]
        if r["table"] != "DERIVED" and conv:
            if r["itemid"] not in CONVERSIONS or CONVERSIONS[r["itemid"]][0] != conv:
                raise AuditError(f"unrecognized unit_conversion for item {r['itemid']}: {conv!r}")
        info = {"table": r["table"], "itemid": r["itemid"], "channel": ch, "low": lo, "high": hi,
                "triggered": r["triggered"] == "triggered", "contract16": r["contract16"] == "likely",
                "status": r["status"], "system": r["table51_system"], "units": r["units"]}
        if r["table"] != "DERIVED":
            items[(r["table"], r["itemid"])] = info
        else:
            dgates[ch] = {"low": lo, "high": hi}
        c = chans.setdefault(ch, {"items": [], "tables": set(), "triggered": False, "contract16": False,
                                  "system": r["table51_system"], "derived": r["table"] == "DERIVED", "merge": False})
        c["items"].append(r["itemid"])
        c["tables"].add(r["table"])
        c["triggered"] |= info["triggered"]
        c["contract16"] |= info["contract16"]
        if r["status"] == "assigned_merge_pending_stage2":
            c["merge"] = True
            m = re.findall(r"main item (\d+)", r["merge_rule"]) or re.findall(r"same scale as (\d+)", r["merge_rule"])
            if not m:
                raise AuditError(f"merge group {ch}: no main item in merge_rule of {r['itemid']}")
            if ch in main and main[ch] != m[0]:
                raise AuditError(f"merge group {ch}: conflicting main items")
            main[ch] = m[0]
    for ch, c in chans.items():
        if c["merge"] and main.get(ch) not in c["items"]:
            raise AuditError(f"merge group {ch}: main item {main.get(ch)} is not one of its items")
        c["tables"] = sorted(c["tables"])
        c["source"] = ("derived" if c["derived"] else "chart" if c["tables"] == ["CHARTEVENTS"] else
                       "lab" if c["tables"] == ["LABEVENTS"] else "mixed")
        if c["source"] == "mixed":
            raise AuditError(f"channel {ch} mixes tables {c['tables']}")
    for need in ("heart_rate", "sbp", "dbp", "map", "spo2", "fio2", "po2", "pf_ratio", "sf_ratio", "weight",
                 "urine_output", "fluid_intake", "pf_ratio_lab_crosscheck"):
        if need not in chans:
            raise AuditError(f"channel map has no '{need}' channel")
    for need in ("map", "pf_ratio"):
        g = dgates.get(need, {})
        if g.get("low") is None or g.get("high") is None:
            raise AuditError(f"derived channel {need} has no plausibility gate in the channel map")
    return {"items": items, "channels": chans, "merge_main": main, "derived_gates": dgates}


def locf_family(ch: str, triggered: bool) -> tuple[str, float]:
    if triggered:
        return "triggered_no_carry", 24.0
    if ch in LOCF_4BIN:
        return "4_bins_24h", 24.0
    if ch in LOCF_1BIN or ch in ("pf_ratio", "sf_ratio", "fio2", "aa_gradient"):
        return "1_bin_6h", 6.0
    return "no_carry", 24.0


def parse_drug_map(dm: pd.DataFrame) -> dict:
    out = {}
    for _, r in dm.iterrows():
        g = r["group"]
        if g not in DRUG_MAP_GROUPS:
            raise AuditError(f"drug map: unknown group {g!r}")
        other = None
        if g == "osmotic_conditional":
            m = re.search(r"otherwise\s*->\s*([a-z_]+)", r["conditional_rule"])
            if not m or m.group(1) not in MED_GROUPS:
                raise AuditError(f"drug map: cannot parse the osmotic rule for {r['drug_name_en']!r}")
            if "neuro admitting dx -> osmotherapy" not in r["conditional_rule"] or "sensitivity: all -> osmotherapy" not in r["conditional_rule"]:
                raise AuditError(f"drug map: osmotic rule text differs from D3 for {r['drug_name_en']!r}")
            other = m.group(1)
        key = f"{r['drug_name']} | {r['drug_name_en']}"
        if key in out:
            raise AuditError(f"drug map: duplicate key {key!r}")
        out[key] = {"group": g, "otherwise": other, "subgroup": r["sensitivity_subgroup"]}
    return out


def parse_specimen_map(sm: pd.DataFrame) -> dict:
    out = {}
    for _, r in sm.iterrows():
        if r["specimen_group"] not in SPECIMEN_GROUPS_V1:
            raise AuditError(f"specimen map: unknown group {r['specimen_group']!r}")
        if r["base_code"] != r["spec_itemid"].split("+")[0]:
            raise AuditError(f"specimen map: base_code is not the first code of {r['spec_itemid']!r}")
        out[f"{r['spec_itemid']} | {r['spec_type_desc']}"] = {"group": r["specimen_group"], "base": r["base_code"],
                                                               "addons": r["result_leak_risk"] == "YES"}
    bases = {}
    for v in out.values():
        if bases.setdefault(v["base"], v["group"]) != v["group"]:
            raise AuditError(f"specimen map: base_code {v['base']} maps to two groups")
    return out


def organism_result_types(category: str) -> tuple[str, ...]:
    c = category.lower()
    if c.startswith(("gram-", "fungus", "atypical")):
        return ("gram_class", "identification")
    if c.startswith("smear finding"):
        return ("gram_class",)
    if c.startswith(("negative/no growth", "normal flora")):
        return ("negative_final",)
    if c.startswith(("blank", "specimen quality", "pending")):
        return ()
    raise AuditError(f"organism map: unknown category {category!r}")


def parse_organism_map(om: pd.DataFrame) -> dict:
    out = {}
    for _, r in om.iterrows():
        key = f"{r['org_itemid']} | {r['org_name']}"
        out[key] = {"category": r["category"], "types": organism_result_types(r["category"]),
                    "isolate": organism_result_types(r["category"]) == ("gram_class", "identification")}
    return out


# =============================================================================
# Raw tables: gate, read, link
# =============================================================================
def raw_path(cfg: S2Config, name: str) -> Path:
    p = Path(cfg.raw_dir) / f"{name}.csv"
    if not p.exists():
        raise AuditError(f"raw table not found: {p}")
    return p


def gate_raw(cfg: S2Config, name: str) -> dict:
    p = raw_path(cfg, name)
    if name in cfg.raw_sha256:
        return {"file": p.name, "sha256": verify_hash(f"raw {name}", p, cfg.raw_sha256[name]), "gate": "stage1_sha256"}
    return {"file": p.name, "sha256": sha256_file(p), "gate": "recorded_only_row_count_gated"}


ROW_CLASSES = ("in_stay_train", "in_stay_validation", "in_stay_test", "cohort_admission_outside_stay",
               "no_cohort_admission", "cohort_admission_time_missing")


def link_chunk(subj: pd.Series, hadm: pd.Series, t: pd.Series, eh: pd.DataFrame, offset_h: float,
               keep_lo: float, classes: dict, keep_all: bool = False) -> pd.DataFrame:
    """Contract H3 linkage: SUBJECT_ID and HADM_ID in hadm_id_list, then time. Row classes use the availability time
    (t + offset) against [t0, OUT). Returns rows kept for analysis: r (row position), enc_i, th = (t - t0) in hours
    (chart time, before any offset), for chart times in [keep_lo, W) where W = min(OUT, t0 + 336 h)."""
    n = len(subj)
    df = pd.DataFrame({"r": np.arange(n), "key": (subj + "|" + hadm).to_numpy(), "t": t.to_numpy()})
    ok = ((subj != "") & (hadm != "")).to_numpy()
    m = df[ok].merge(eh, on="key", how="inner")
    cohort_r = pd.unique(m["r"])
    classes["no_cohort_admission"] += n - len(cohort_r)
    tm = m["t"].isna()
    tmiss_r = pd.unique(m.loc[tm, "r"])
    classes["cohort_admission_time_missing"] += len(tmiss_r)
    mv = m[~tm].copy()
    mv["th"] = (mv["t"] - mv["t0"]).dt.total_seconds().to_numpy() / 3600.0
    ta = mv["th"] + offset_h
    mv["stay"] = (ta >= 0) & (ta < hours(mv["OUT"] - mv["t0"]))
    st = mv[mv["stay"]].drop_duplicates("r")
    for sp in SPLITS:
        classes[f"in_stay_{sp}"] += int((st["split"] == sp).sum())
    classes["cohort_admission_outside_stay"] += len(cohort_r) - len(tmiss_r) - len(st)
    keep = mv if keep_all else mv[(mv["th"] >= keep_lo) & (mv["th"] < mv["W"])]
    # one encounter per row: the stay that contains it, else the encounter it precedes most closely, else the latest
    keep = keep.assign(_pre=keep["th"] < 0, _abs=keep["th"].abs())
    keep = keep.sort_values(["r", "stay", "_pre", "_abs"], ascending=[True, False, False, True]).drop_duplicates("r")
    return keep[["r", "enc_i", "th", "stay"]].reset_index(drop=True)


def read_chunks(cfg: S2Config, A, name: str, cols: list[str]):
    return pd.read_csv(raw_path(cfg, name), dtype=str, keep_default_na=False, na_filter=False, usecols=cols,
                       encoding="utf-8-sig", chunksize=cfg.chunksize)


def finish_classes(cfg: S2Config, name: str, classes: dict, rows: int) -> None:
    if sum(classes.values()) != rows:
        raise AuditError(f"{name}: row classes do not sum to the rows read")
    exp = cfg.expected_raw_rows.get(name)
    if exp is not None and rows != exp:
        raise AuditError(f"{name}: {rows} rows read, schema inventory says {exp}")


def read_events(cfg: S2Config, A, co: Cohort, name: str, wanted_items: set[str], offset_h: float,
                stays: pd.DataFrame | None) -> tuple[pd.DataFrame, dict, dict]:
    """CHARTEVENTS / LABEVENTS: all rows are classified; mapped items in the analysis range are kept."""
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc = resolve(hdr, ("SUBJECT_ID",), name), resolve(hdr, ("HADM_ID",), name)
    ic, tc = resolve(hdr, ("ITEMID",), name), resolve(hdr, ("CHARTTIME",), name)
    vc = resolve(hdr, ("VALUE",), name, required=False)
    nc = resolve(hdr, ("VALUENUM",), name, required=False)
    uc = resolve(hdr, ("VALUEUOM",), name, required=False)
    icu = resolve(hdr, ("ICUSTAY_ID",), name, required=False)
    used = [c for c in (sc, hc, ic, tc, vc, nc, uc, icu) if c]
    classes = {c: 0 for c in ROW_CLASSES}
    pre_any = {w: np.zeros(co.n, dtype=bool) for w in PRE_ICU_H}
    rows, parts = 0, []
    icu_check = {"rows_in_stay_checked": 0, "rows_icustay_blank": 0, "rows_icustay_not_in_encounter": 0}
    for raw in read_chunks(cfg, A, name, used):
        rows += len(raw)
        subj, hadm, t = vnorm_id(raw[sc]), vnorm_id(raw[hc]), fast_time(raw[tc])
        k = link_chunk(subj, hadm, t, co.eh, offset_h, KEEP_LO_H, classes)
        if not len(k):
            continue
        ta = k["th"].to_numpy() + offset_h
        for w in PRE_ICU_H:
            sel = (ta >= -w) & (ta < 0)
            pre_any[w][k["enc_i"].to_numpy()[sel]] = True
        r = k["r"].to_numpy()
        item = vnorm_id(raw[ic].iloc[r]).to_numpy()
        if icu and stays is not None:
            ins = k["stay"].to_numpy()
            ids = vnorm_id(raw[icu].iloc[r[ins]]).to_numpy()
            icu_check["rows_in_stay_checked"] += int(ins.sum())
            icu_check["rows_icustay_blank"] += int((ids == "").sum())
            pairs = pd.DataFrame({"enc_i": k["enc_i"].to_numpy()[ins], "ICUSTAY_ID": ids})
            pairs = pairs[pairs["ICUSTAY_ID"] != ""]
            mm = pairs.merge(stays[["enc_i", "ICUSTAY_ID"]].drop_duplicates().assign(ok=True), how="left",
                             on=["enc_i", "ICUSTAY_ID"])
            icu_check["rows_icustay_not_in_encounter"] += int(mm["ok"].isna().sum())
        sel = np.isin(item, list(wanted_items))
        if not sel.any():
            continue
        rr = r[sel]
        vnum = pd.to_numeric(raw[nc].iloc[rr].str.strip(), errors="coerce").to_numpy() if nc else np.full(len(rr), np.nan)
        vtext = vclean_text(raw[vc].iloc[rr]).to_numpy() if vc else np.full(len(rr), "", dtype=object)
        if vc:
            alt = pd.to_numeric(pd.Series(vtext), errors="coerce").to_numpy()
            vnum = np.where(np.isnan(vnum), alt, vnum)
        parts.append(pd.DataFrame({"enc_i": k["enc_i"].to_numpy()[sel].astype(np.int32), "th": k["th"].to_numpy()[sel],
                                   "item": item[sel], "vnum": vnum, "vtext": vtext,
                                   "unit": norm_unit(raw[uc].iloc[rr]).to_numpy() if uc else ""}))
    finish_classes(cfg, name, classes, rows)
    ev = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int32), "th": pd.Series(dtype=float), "item": pd.Series(dtype=object),
         "vnum": pd.Series(dtype=float), "vtext": pd.Series(dtype=object), "unit": pd.Series(dtype=object)})
    log(cfg, f"{name}: {rows} rows read; {len(ev)} mapped-item rows kept in [t0 - 26 h, W)")
    return ev, {"row_classes": classes, "rows": rows, "pre_any": pre_any, "icustay_check": icu_check if icu else None}, {}


def assign_date_only_start(th: np.ndarray, date_only: np.ndarray) -> np.ndarray:
    """A date-only STARTDATE parses to its midnight (hours from t0). Assigned start (contract; plan v12 S3 (e)):
    the start of the first in-window bin on that date -> t0 on the t0 date, else the first bin start >= midnight.
    Dates before the t0 date keep their midnight (pre-ICU documentation)."""
    th = np.asarray(th, dtype=float).copy()
    do = np.asarray(date_only, dtype=bool)
    later = do & (th > 0)
    same = do & (th <= 0) & (th > -24)
    th[later] = np.ceil(th[later] / BIN_H - 1e-12) * BIN_H
    th[same] = 0.0
    return th


def read_prescriptions(cfg: S2Config, A, co: Cohort) -> tuple[pd.DataFrame, dict]:
    name = "PRESCRIPTIONS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc = resolve(hdr, ("SUBJECT_ID",), name), resolve(hdr, ("HADM_ID",), name)
    stc = resolve(hdr, ("STARTDATE", "STARTTIME"), name)
    enc_ = resolve(hdr, ("ENDDATE", "ENDTIME"), name, required=False)
    dn = resolve(hdr, ("DRUG_NAME",), name, required=False)
    de = resolve(hdr, ("DRUG_NAME_EN",), name, required=False)
    if not dn and not de:
        raise AuditError(f"{name}: neither DRUG_NAME nor DRUG_NAME_EN found")
    used = [c for c in (sc, hc, stc, enc_, dn, de) if c]
    classes = {c: 0 for c in ROW_CLASSES}
    pre_any = {w: np.zeros(co.n, dtype=bool) for w in PRE_ICU_H}
    rows, parts = 0, []
    chk = {k: 0 for k in ("cohort_rows", "in_stay", "before_t0_same_date", "before_t0_earlier_date", "at_or_after_out",
                          "time_missing", "clock_time", "date_only", "midnight_clock", "in_stay_date_only",
                          "in_stay_midnight_clock", "date_only_on_t0_date_before_t0")}
    t0_by = co.enc.set_index("enc_i")["t0"]
    for raw in read_chunks(cfg, A, name, used):
        rows += len(raw)
        subj, hadm = vnorm_id(raw[sc]), vnorm_id(raw[hc])
        s_raw = raw[stc]
        t = fast_time(s_raw)
        clock = has_clock(s_raw)
        k = link_chunk(subj, hadm, t, co.eh, 0.0, -1e9, classes, keep_all=True)  # every linked, timed row
        if not len(k):
            continue
        r = k["r"].to_numpy()
        e = k["enc_i"].to_numpy()
        th = k["th"].to_numpy()
        cl = clock.to_numpy()[r]
        tt = t.iloc[r]
        mid = cl & (tt.dt.hour.to_numpy() == 0) & (tt.dt.minute.to_numpy() == 0) & (tt.dt.second.to_numpy() == 0)
        t0d = t0_by.reindex(e).dt.normalize().to_numpy()
        same_date = tt.dt.normalize().to_numpy() == t0d
        out_h = co.out_h[e]
        chk["cohort_rows"] += len(k)
        chk["in_stay"] += int(((th >= 0) & (th < out_h)).sum())
        chk["before_t0_same_date"] += int(((th < 0) & same_date).sum())
        chk["before_t0_earlier_date"] += int(((th < 0) & ~same_date).sum())
        chk["at_or_after_out"] += int((th >= out_h).sum())
        chk["clock_time"] += int(cl.sum())
        chk["date_only"] += int((~cl).sum())
        chk["midnight_clock"] += int(mid.sum())
        ins = (th >= 0) & (th < out_h)
        chk["in_stay_date_only"] += int((ins & ~cl).sum())
        chk["in_stay_midnight_clock"] += int((ins & mid).sum())
        chk["date_only_on_t0_date_before_t0"] += int((~cl & same_date & (th < 0)).sum())
        for w in PRE_ICU_H:
            sel = (th >= -w) & (th < 0)
            pre_any[w][e[sel]] = True
        th = assign_date_only_start(th, ~cl)
        keep = (th >= MED_KEEP_LO_H) & (th < co.W[e])
        if not keep.any():
            continue
        rk = r[keep]
        a = vclean_text(raw[dn].iloc[rk]) if dn else pd.Series("", index=raw.index[rk])
        b = vclean_text(raw[de].iloc[rk]) if de else pd.Series("", index=raw.index[rk])
        key = (a + " | " + b).to_numpy()
        if enc_:
            e_raw = raw[enc_].iloc[rk]
            et = fast_time(e_raw)
            eclock = has_clock(e_raw).to_numpy()
            et = et + pd.to_timedelta(np.where(eclock | et.isna().to_numpy(), 0, 24), unit="h").to_numpy()  # date-only end -> next midnight
            eh_ = ((et - t0_by.reindex(e[keep]).to_numpy()).dt.total_seconds() / 3600.0).to_numpy()
        else:
            eh_ = np.full(len(rk), np.nan)
        parts.append(pd.DataFrame({"enc_i": e[keep].astype(np.int32), "S": th[keep], "E": eh_, "key": key,
                                   "date_only": ~cl[keep], "midnight": mid[keep]}))
    finish_classes(cfg, name, classes, rows)
    chk["time_missing"] = classes["cohort_admission_time_missing"]
    docs = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int32), "S": pd.Series(dtype=float), "E": pd.Series(dtype=float),
         "key": pd.Series(dtype=object), "date_only": pd.Series(dtype=bool), "midnight": pd.Series(dtype=bool)})
    log(cfg, f"{name}: {rows} rows read; {len(docs)} documentation rows kept in [t0 - 72 h, W)")
    return docs, {"row_classes": classes, "rows": rows, "pre_any": pre_any, "start_check": chk}


def read_micro(cfg: S2Config, A, co: Cohort) -> tuple[pd.DataFrame, dict]:
    name = "MICROBIOLOGYEVENTS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc = resolve(hdr, ("SUBJECT_ID",), name), resolve(hdr, ("HADM_ID",), name)
    tc = resolve(hdr, ("CHARTTIME",), name)
    si, sd = resolve(hdr, ("SPEC_ITEMID",), name, required=False), resolve(hdr, ("SPEC_TYPE_DESC",), name)
    oi, on = resolve(hdr, ("ORG_ITEMID",), name, required=False), resolve(hdr, ("ORG_NAME",), name)
    classes = {c: 0 for c in ROW_CLASSES}
    pre_any = {w: np.zeros(co.n, dtype=bool) for w in PRE_ICU_H}
    fill = {"rows_all": 0, "rows_in_stay_train": 0, "columns": {c: {"all": 0, "in_stay_train": 0} for c in hdr}}
    rows, parts = 0, []
    for raw in read_chunks(cfg, A, name, hdr):
        rows += len(raw)
        subj, hadm, t = vnorm_id(raw[sc]), vnorm_id(raw[hc]), fast_time(raw[tc])
        nonblank = {c: vclean_text(raw[c]) != "" for c in hdr}
        fill["rows_all"] += len(raw)
        for c in hdr:
            fill["columns"][c]["all"] += int(nonblank[c].sum())
        k = link_chunk(subj, hadm, t, co.eh, 0.0, KEEP_LO_H, classes)
        if not len(k):
            continue
        r, e, th = k["r"].to_numpy(), k["enc_i"].to_numpy(), k["th"].to_numpy()
        tr = k["stay"].to_numpy() & (co.split[e] == "train")
        # in-stay train rows (the window is W; rows beyond 336 h in stay are counted through the class totals only)
        fill["rows_in_stay_train"] += int(tr.sum())
        for c in hdr:
            fill["columns"][c]["in_stay_train"] += int(nonblank[c].to_numpy()[r[tr]].sum())
        for w in PRE_ICU_H:
            sel = (th >= -w) & (th < 0)
            pre_any[w][e[sel]] = True
        sid = vnorm_id(raw[si].iloc[r]) if si else pd.Series("", index=raw.index[r])
        sdesc = vclean_text(raw[sd].iloc[r])
        oid = vnorm_id(raw[oi].iloc[r]) if oi else pd.Series("", index=raw.index[r])
        oname = vclean_text(raw[on].iloc[r])
        parts.append(pd.DataFrame({"enc_i": e.astype(np.int32), "th": th, "spec_key": (sid + " | " + sdesc).to_numpy(),
                                   "spec_itemid": sid.to_numpy(), "org_key": (oid + " | " + oname).to_numpy()}))
    finish_classes(cfg, name, classes, rows)
    ev = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int32), "th": pd.Series(dtype=float), "spec_key": pd.Series(dtype=object),
         "spec_itemid": pd.Series(dtype=object), "org_key": pd.Series(dtype=object)})
    log(cfg, f"{name}: {rows} rows read; {len(ev)} rows kept in [t0 - 26 h, W)")
    return ev, {"row_classes": classes, "rows": rows, "pre_any": pre_any, "fill": fill}


def read_icustays(cfg: S2Config, A, enc: pd.DataFrame) -> tuple[pd.DataFrame, dict[int, list], dict]:
    """ICUSTAYS rebuild: stays of the encounter's subject and admissions within [t0, OUT] (1 min tolerance).
    Gaps between consecutive stays are the merge gaps (hours from t0); first_careunit = earliest stay's unit."""
    name = "ICUSTAYS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    cols = [resolve(hdr, (c,), name) for c in ("SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "INTIME", "OUTTIME")]
    fcu = resolve(hdr, ("FIRST_CAREUNIT",), name)
    d = A.read_csv_str(p, usecols=cols + [fcu])
    if len(d) != cfg.expected_raw_rows.get(name, len(d)):
        raise AuditError(f"{name}: {len(d)} rows, schema inventory says {cfg.expected_raw_rows.get(name)}")
    st = pd.DataFrame({"key": vnorm_id(d[cols[0]]) + "|" + vnorm_id(d[cols[1]]), "ICUSTAY_ID": vnorm_id(d[cols[2]]),
                       "IN": fast_time(d[cols[3]]), "OUTS": fast_time(d[cols[4]]), "unit": vclean_text(d[fcu])})
    eh = encounter_hadm_index(enc)
    m = st.merge(eh[["key", "enc_i", "t0", "OUT"]], on="key", how="inner")
    tol = pd.Timedelta(minutes=1)
    m = m[(m["IN"] >= m["t0"] - tol) & (m["OUTS"] <= m["OUT"] + tol)].sort_values(["enc_i", "IN"])
    m = m.drop_duplicates(["enc_i", "ICUSTAY_ID"])
    gaps: dict[int, list] = {}
    for i, g in m.groupby("enc_i", sort=False):
        ins, outs = g["IN"].to_numpy(), g["OUTS"].to_numpy()
        t0 = g["t0"].iloc[0]
        gl = []
        run_out = outs[0]
        for j in range(1, len(g)):
            if ins[j] > run_out:
                gl.append((float((run_out - t0) / np.timedelta64(1, "h")), float((ins[j] - t0) / np.timedelta64(1, "h"))))
            run_out = max(run_out, outs[j])
        if gl:
            gaps[int(i)] = gl
    first = m.drop_duplicates("enc_i", keep="first").set_index("enc_i")["unit"]
    info = {"encounters_without_a_matching_stay": int(len(enc) - m["enc_i"].nunique()),
            "encounters_with_gaps": len(gaps)}
    return m[["enc_i", "ICUSTAY_ID"]], gaps, {"first_unit": first, "info": info, "rows": len(d)}


# =============================================================================
# Channel values: gates, conversions, HR merge, MAP, P/F, S/F (spec 5.3; D6-D9, D13, D15, D19)
# =============================================================================
ROW_STATUS = ("measurement", "out_of_gate", "nonnumeric", "blank", "unit_pending", "unit_unrecognized")


def classify_item_rows(info: dict, df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Value in the gate unit (the map's recorded unit) and the row status. Gated channels need a numeric value inside
    [low, high]; ungated (qualitative) channels need a numeric value or a nonblank text."""
    it = info["itemid"]
    v = df["vnum"].to_numpy(dtype=float).copy()
    status = np.full(len(df), "measurement", dtype=object)
    conv = CONVERSIONS.get(it)
    if conv and conv[1] == "per_row_gate":
        u = df["unit"].astype(str).to_numpy()
        known = np.isin(u, list(conv[2]))
        f = np.array([conv[2].get(x, np.nan) if x in conv[2] and conv[2].get(x) is not None else np.nan for x in u])
        pending = known & np.isnan(f)
        status[~known] = "unit_unrecognized"
        status[pending] = "unit_pending"
        v = v * f
    has_num = ~np.isnan(v)
    free = status == "measurement"
    if info["low"] is not None or info["high"] is not None:
        lo = -np.inf if info["low"] is None else info["low"]
        hi = np.inf if info["high"] is None else info["high"]
        inside = has_num & (v >= lo) & (v <= hi)
        status[free & ~has_num] = "nonnumeric"
        status[free & has_num & ~inside] = "out_of_gate"
    else:
        txt = df["vtext"].astype(str).to_numpy() != ""
        status[free & ~has_num & ~txt] = "blank"
    return v, status


@dataclass
class ValueBuild:
    events: dict            # channel -> DataFrame(enc_i, t, tc, item, v)
    rowstat: pd.DataFrame   # enc_i, avail, channel, item, status, unit (all rows of mapped items)
    hr: dict
    map_pairs: dict
    pf: dict
    sf: dict


def _empty_events() -> pd.DataFrame:
    return pd.DataFrame({"enc_i": pd.Series(dtype=np.int32), "t": pd.Series(dtype=float), "tc": pd.Series(dtype=float),
                         "item": pd.Series(dtype=object), "v": pd.Series(dtype=float)})


def build_values(chart: pd.DataFrame, lab: pd.DataFrame, cm: dict) -> ValueBuild:
    rows = []
    for table, ev, off in (("CHARTEVENTS", chart, 0.0), ("LABEVENTS", lab, LAB_RELEASE_H)):
        for it, g in ev.groupby("item", sort=True):
            info = cm["items"].get((table, it))
            if info is None:
                continue
            v, st = classify_item_rows(info, g)
            rows.append(pd.DataFrame({"enc_i": g["enc_i"].to_numpy(), "tc": g["th"].to_numpy(), "t": g["th"].to_numpy() + off,
                                      "item": it, "channel": info["channel"], "v": v, "status": st,
                                      "vnum": g["vnum"].to_numpy(), "vtext": g["vtext"].to_numpy(), "unit": g["unit"].to_numpy(),
                                      "table": table}))
    allr = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["enc_i", "tc", "t", "item", "channel", "v", "status", "vnum", "vtext", "unit", "table"])
    meas = allr[allr["status"] == "measurement"]
    events: dict[str, pd.DataFrame] = {}
    # HR (D9): 1003 kept when both exist at the same CHARTTIME
    hrm = meas[meas["channel"] == "heart_rate"]
    k3 = hrm[hrm["item"] == HR_KEEP][["enc_i", "tc", "v"]].drop_duplicates(["enc_i", "tc"])
    k2 = hrm[hrm["item"] == HR_DROP]
    both = k2.merge(k3, on=["enc_i", "tc"], how="left", suffixes=("", "_1003"), indicator=True)
    drop2 = (both["_merge"] == "both").to_numpy()
    hr_keep = pd.concat([hrm[hrm["item"] != HR_DROP], k2[~drop2]], ignore_index=True)
    hr = {"pairs": both[drop2][["enc_i", "tc", "v", "v_1003"]].reset_index(drop=True), "dropped_1002": k2[drop2]}
    for ch, g in meas.groupby("channel", sort=True):
        if ch == "heart_rate":
            g = hr_keep
        events[ch] = g[["enc_i", "t", "tc", "item", "v"]].sort_values(["enc_i", "t"], kind="mergesort").reset_index(drop=True)
    # MAP (D8): same-CHARTTIME SBP/DBP pairs from numeric rows; drop DBP >= SBP or any value outside its gate
    gs = {k: cm["items"][("CHARTEVENTS", it)] for k, it in (("sbp", SBP_ITEM), ("dbp", DBP_ITEM))}
    num = allr[allr["table"] == "CHARTEVENTS"]
    sb = num[(num["item"] == SBP_ITEM) & num["vnum"].notna()].drop_duplicates(["enc_i", "tc"])
    db = num[(num["item"] == DBP_ITEM) & num["vnum"].notna()].drop_duplicates(["enc_i", "tc"])
    pr = sb[["enc_i", "tc", "t", "vnum"]].merge(db[["enc_i", "tc", "vnum"]], on=["enc_i", "tc"], suffixes=("_s", "_d"))
    s_, d_ = pr["vnum_s"].to_numpy(dtype=float), pr["vnum_d"].to_numpy(dtype=float)
    mapv = (s_ + 2 * d_) / 3.0
    dmap = cm["derived_gates"]["map"]
    ge = d_ >= s_
    in_s = (s_ >= gs["sbp"]["low"]) & (s_ <= gs["sbp"]["high"])
    in_d = (d_ >= gs["dbp"]["low"]) & (d_ <= gs["dbp"]["high"])
    in_m = (mapv >= dmap["low"]) & (mapv <= dmap["high"])
    reason = np.where(ge, "dbp_ge_sbp", np.where(~(in_s & in_d & in_m), "outside_gate", "valid"))
    pr = pr.assign(v=mapv, reason=reason)
    events["map"] = pr[pr["reason"] == "valid"].assign(item="derived:map")[["enc_i", "t", "tc", "item", "v"]] \
        .sort_values(["enc_i", "t"], kind="mergesort").reset_index(drop=True)
    map_pairs = {"sbp_rows": sb[["enc_i", "t"]], "dbp_rows": db[["enc_i", "t"]], "pairs": pr}
    # P/F (D7, D19): pO2 and FiO2 5222 from the same draw (encounter and CHARTTIME); 5239 preferred over 5244
    lm = meas[meas["table"] == "LABEVENTS"]
    po = lm[lm["item"].isin(PO2_ITEMS)].assign(_o=lambda d: (d["item"] != PO2_ITEMS[0]).astype(int))
    po = po.sort_values(["enc_i", "tc", "_o"]).drop_duplicates(["enc_i", "tc"])
    fi = lm[lm["item"] == FIO2_ITEM].drop_duplicates(["enc_i", "tc"])
    pf = po[["enc_i", "tc", "t", "v"]].merge(fi[["enc_i", "tc", "v"]], on=["enc_i", "tc"], suffixes=("_po2", "_fio2"))
    pf["v"] = pf["v_po2"] / (pf["v_fio2"] / 100.0)
    pfg = cm["derived_gates"]["pf_ratio"]
    pf["valid"] = (pf["v"] >= pfg["low"]) & (pf["v"] <= pfg["high"])
    lab_pf = lm[lm["item"] == PF_LAB_ITEM].drop_duplicates(["enc_i", "tc"])[["enc_i", "tc", "v"]]
    pf = pf.merge(lab_pf.rename(columns={"v": "v_lab"}), on=["enc_i", "tc"], how="left")
    events["pf_ratio"] = pf[pf["valid"]].assign(item="derived:pf_ratio")[["enc_i", "t", "tc", "item", "v"]] \
        .sort_values(["enc_i", "t"], kind="mergesort").reset_index(drop=True)
    # S/F (D15): FiO2 draws without a valid P/F; nearest SpO2 1006 within 5 min; only when that SpO2 <= 97 %
    has_pf = fi[["enc_i", "tc"]].merge(pf[pf["valid"]][["enc_i", "tc"]].assign(_pf=True), on=["enc_i", "tc"], how="left")
    nofi = fi[has_pf["_pf"].isna().to_numpy()][["enc_i", "tc", "t", "v"]].rename(columns={"v": "v_fio2"})
    sp = meas[(meas["table"] == "CHARTEVENTS") & (meas["item"] == SPO2_ITEM)][["enc_i", "tc", "v"]].rename(
        columns={"tc": "tc_spo2", "v": "v_spo2"})
    sf_rows = pd.DataFrame()
    if len(nofi):
        a = nofi.sort_values("tc").reset_index(drop=True)
        b = sp.sort_values("tc_spo2").reset_index(drop=True)
        a["enc_i"] = a["enc_i"].astype(np.int64)
        b["enc_i"] = b["enc_i"].astype(np.int64)
        sf_rows = pd.merge_asof(a, b, left_on="tc", right_on="tc_spo2", by="enc_i", direction="nearest",
                                tolerance=SF_PAIR_MINUTES / 60.0 + 1e-9)
        sf_rows["reason"] = np.where(sf_rows["v_spo2"].isna(), "no_spo2_within_5min",
                                     np.where(sf_rows["v_spo2"] > SF_SPO2_MAX, "spo2_gt_97", "sf"))
        sf_rows["v"] = sf_rows["v_spo2"] / (sf_rows["v_fio2"] / 100.0)
        sf_rows["t"] = np.maximum(sf_rows["t"].to_numpy(), sf_rows["tc_spo2"].fillna(-np.inf).to_numpy())
    else:
        sf_rows = pd.DataFrame(columns=["enc_i", "tc", "t", "v", "reason", "v_fio2", "v_spo2", "tc_spo2"])
    sfv = sf_rows[sf_rows["reason"] == "sf"] if len(sf_rows) else sf_rows
    events["sf_ratio"] = (sfv.assign(item="derived:sf_ratio")[["enc_i", "t", "tc", "item", "v"]]
                          .astype({"enc_i": np.int32}).sort_values(["enc_i", "t"], kind="mergesort").reset_index(drop=True)
                          if len(sfv) else _empty_events())
    for ch in cm["channels"]:
        events.setdefault(ch, _empty_events())
    return ValueBuild(events=events, rowstat=allr, hr=hr, map_pairs=map_pairs,
                      pf={"pairs": pf, "fio2_draws": fi[["enc_i", "tc", "t"]]}, sf={"rows": sf_rows})


# =============================================================================
# Generic per-encounter features of one event channel (metrics 1-8)
# =============================================================================
def _avail(co: Cohort, e: np.ndarray, t: np.ndarray) -> dict:
    out = {}
    for L in LANDMARKS:
        for kind, lo in (("full", 0.0), ("last24", L - 24.0), ("72h", max(0.0, L - 72.0))):
            m = (t >= max(lo, 0.0)) & (t < L)
            out[(kind, L)] = np.bincount(e[m], minlength=co.n) > 0
    return out


def event_features(co: Cohort, e: np.ndarray, t: np.ndarray, alt: dict | None = None) -> dict:
    """e: encounter index, t: availability time (hours from t0). Events outside [-24 h, W) are ignored except for the
    pre-ICU flags. alt: {offset: t_alt} availability times under the lab-release sensitivities."""
    n = co.n
    e = np.asarray(e, dtype=np.int64)
    t = np.asarray(t, dtype=float)
    o = np.lexsort((t, e))
    e, t = e[o], t[o]
    W = co.W[e]
    inw = (t >= 0) & (t < W)
    ei, ti = e[inw], t[inw]
    F: dict[str, Any] = {}
    F["ever_w"] = np.bincount(ei, minlength=n) > 0
    bend = np.minimum(co.out_h, LAYER_B_H)[e]
    F["ever_b"] = np.bincount(e[(t >= 0) & (t < bend)], minlength=n) > 0
    F["cnt"] = np.bincount(ei, minlength=n)
    first = np.full(n, np.nan)
    if len(ei):
        idx = np.r_[0, np.flatnonzero(np.diff(ei)) + 1]
        first[ei[idx]] = ti[idx]
    F["first"] = first
    for w in PRE_ICU_H:
        F[f"pre{int(w)}"] = np.bincount(e[(t >= -w) & (t < 0)], minlength=n) > 0
    F["avail"] = {1.0: _avail(co, e, t)}
    for off, ta in (alt or {}).items():
        ta = np.asarray(ta, dtype=float)[o]
        F["avail"][off] = _avail(co, e, ta)
    last = {}
    for L in LANDMARKS:
        m = (t >= 0) & (t < L)
        lt = np.full(n, -np.inf)
        np.maximum.at(lt, e[m], t[m])
        last[L] = np.where(np.isfinite(lt), L - lt, np.nan)
    F["last_age"] = last
    k = np.floor(ti / BIN_H).astype(np.int64)
    k = np.minimum(k, N_BINS - 1)
    F["occ"] = np.bincount(ei * N_BINS + k, minlength=n * N_BINS).reshape(n, N_BINS)
    same = ei[1:] == ei[:-1]
    d = (ti[1:] - ti[:-1])[same]
    de = ei[1:][same]
    tp, tn = ti[:-1][same], ti[1:][same]
    cross = np.zeros(len(d), dtype=bool)
    for i, gl in co.gaps.items():
        sel = de == i
        if sel.any():
            for a, b in gl:
                cross[sel] |= (tp[sel] <= a) & (tn[sel] >= b)
    F["int_e"], F["int_d"], F["int_cross"] = de, d, cross
    nz = d > 0
    med = np.full(n, np.nan)
    if nz.any():
        s = pd.Series(d[nz]).groupby(de[nz]).median()
        med[s.index.to_numpy()] = s.to_numpy()
    F["int_med"] = med
    return F


@dataclass
class Dens:
    enc: dict
    elig: dict
    w6: dict
    bin_full: dict
    strata: dict
    strata_levels: dict
    outcome: dict
    n_all: str = ""


def register_denominators(out: Out, co: Cohort, strata: dict[str, np.ndarray], outcome_on: bool, died: np.ndarray) -> Dens:
    enc, elig, w6, bf = {}, {}, {}, {}
    tot = out.cell(("denominators", "n_encounters"), co.n)
    out.frozen.add(tot)
    for s in SPLITS:
        M = co.masks[s]
        enc[s] = out.cell(("denominators", s, "n_encounters"), M.sum())
        out.frozen.add(enc[s])  # the split sizes are published in cohort_checks and in the earlier audits
        elig[s] = {L: out.cell(("denominators", s, "eligible", f"{L}h", "n_eligible"), (M & co.elig[L]).sum(), parent=enc[s])
                   for L in LANDMARKS}
        w6[s] = out.cell(("denominators", s, "n_open_window_ge_6h"), (M & (co.W >= BIN_H)).sum(), parent=enc[s])
        bf[s] = {k: out.cell(("denominators", s, "bins", str(k), "n_full_bin_open"), (M & (co.out_h >= BIN_H * (k + 1))).sum(),
                             parent=enc[s]) for k in range(N_BINS)}
    out.ident(tot, [enc[s] for s in SPLITS])
    st, lv = {}, {}
    for var, arr in strata.items():
        levels = sorted(pd.unique(arr).tolist(), key=str)
        lv[var] = levels
        for s in SPLITS:
            M = co.masks[s]
            cells = {}
            for v in levels:
                Mv = M & (arr == v)
                b = ("strata_denominators", var, s, str(v))
                c_enc = out.cell(b + ("n_encounters",), Mv.sum())
                cells[v] = {"enc": c_enc,
                            "elig": out.cell(b + ("n_eligible_72h",), (Mv & co.elig[STRATA_LANDMARK]).sum(), parent=c_enc),
                            "w6": out.cell(b + ("n_open_window_ge_6h",), (Mv & (co.W >= BIN_H)).sum(), parent=c_enc)}
            out.ident(enc[s], [c["enc"] for c in cells.values()])
            out.ident(elig[s][STRATA_LANDMARK], [c["elig"] for c in cells.values()])
            out.ident(w6[s], [c["w6"] for c in cells.values()])
            st[(var, s)] = cells
    oc: dict = {}
    if outcome_on:
        M = co.masks["train"]
        for L in LANDMARKS:
            g = outcome_group(died, co.out_h, L)
            cells = {}
            for grp in OUTCOME_GROUPS:
                cells[grp] = out.cell(("outcome_denominators", "train", "outcome1", f"{L}h", f"group_{grp}", "n_eligible"),
                                      (M & co.elig[L] & (g == grp)).sum())
            out.ident(elig["train"][L], cells.values())
            oc[("outcome1", L)] = (cells, g)
        mc = {}
        for lab, val in (("died", True), ("survived", False)):
            mc[lab] = {"enc": out.cell(("outcome_denominators", "train", "mortality", lab, "n_encounters"), (M & (died == val)).sum())}
            mc[lab]["w6"] = out.cell(("outcome_denominators", "train", "mortality", lab, "n_open_window_ge_6h"),
                                     (M & (died == val) & (co.W >= BIN_H)).sum(), parent=mc[lab]["enc"])
            for L in LANDMARKS:
                mc[lab][L] = out.cell(("outcome_denominators", "train", "mortality", lab, f"{L}h", "n_eligible"),
                                      (M & (died == val) & co.elig[L]).sum(), parent=mc[lab]["enc"])
        out.ident(enc["train"], [mc["died"]["enc"], mc["survived"]["enc"]])
        out.ident(w6["train"], [mc["died"]["w6"], mc["survived"]["w6"]])
        for L in LANDMARKS:
            out.ident(elig["train"][L], [mc["died"][L], mc["survived"][L]])
        oc["mortality"] = mc
    return Dens(enc=enc, elig=elig, w6=w6, bin_full=bf, strata=st, strata_levels=lv, outcome=oc, n_all=tot)


def write_event_metrics(out: Out, co: Cohort, dens: Dens, base: tuple, F: dict, *, detail: bool, cap_h: float,
                        lab_offsets: bool, strata: dict[str, np.ndarray], outcome_on: bool, died: np.ndarray) -> None:
    """Metrics 1-8 per split (detail outputs only for channels with >= 100 train encounters), strata (1, 2a at 72 h,
    3, 7) and the train-only outcome descriptors (1, 2a, 3)."""
    for s in SPLITS:
        M = co.masks[s]
        b = base + (s,)
        n_ew = out.cell(b + ("n_ever_336h",), (M & F["ever_w"]).sum(), parent=dens.enc[s])
        out.cell(b + ("n_ever_layer_b",), (M & F["ever_b"]).sum(), parent=dens.enc[s])
        out.pct(b, "ever_336h", n_ew, dens.enc[s])
        for w in PRE_ICU_H:
            out.cell(b + (f"n_pre_icu_{int(w)}h",), (M & F[f"pre{int(w)}"]).sum(), parent=dens.enc[s])
        offs = [1.0] + ([0.0, 2.0] if lab_offsets else [])
        for L in LANDMARKS:
            bl = b + ("landmarks", f"{L}h")
            E = M & co.elig[L]
            for off in offs:
                kinds = ("full", "last24", "72h") if (off == 1.0 or detail) else ("full",)
                for kind in kinds:
                    sfx = "" if off == 1.0 else f"_release{int(off)}h"
                    c = out.cell(bl + (f"n_avail_{kind}{sfx}",), (E & F["avail"][off][(kind, L)]).sum(), parent=dens.elig[s][L])
                    if off == 1.0 and kind == "full":
                        n_full = c
                        out.pct(bl, "avail_full", c, dens.elig[s][L])
            if not detail:
                continue
            A_ = E & F["avail"][1.0][("full", L)]
            ages = F["last_age"][L][A_]
            out.q(bl, "last_obs_age_hours", ages, equal_to=n_full)
            out.cell(bl + ("n_last_obs_within_cap",), int((ages <= cap_h).sum()), parent=n_full)
            if L in DETAIL_LANDMARKS:
                nb = int(L // BIN_H)
                occ = F["occ"][:, :nb]
                if s == "train":
                    for k in range(nb):
                        out.cell(bl + ("bins", str(k), "n_occupied"), (E & (occ[:, k] > 0)).sum(), parent=dens.elig[s][L])
                tot = out.cell(bl + ("n_prefix_bins",), int(E.sum()) * nb)
                out.cell(bl + ("n_prefix_bins_empty",), int((occ[E] == 0).sum()), parent=tot)
                totm = out.cell(bl + ("n_prefix_bins_measured",), int(A_.sum()) * nb)
                out.cell(bl + ("n_prefix_bins_measured_empty",), int((occ[A_] == 0).sum()), parent=totm)
                if L == STRATA_LANDMARK:
                    ng = ~co.gap_bins[:, :nb]
                    tg = out.cell(bl + ("n_prefix_bins_measured_excl_gap",), int(ng[A_].sum()))
                    out.cell(bl + ("n_prefix_bins_measured_empty_excl_gap",), int(((occ == 0) & ng)[A_].sum()), parent=tg)
        cw = F["cnt"] / np.maximum(co.W, 1e-9) * 24.0
        G6 = M & (co.W >= BIN_H)
        out.q(b, "rate_per_day_all", cw[G6], equal_to=dens.w6[s])
        out.q(b, "rate_per_day_measured", cw[G6 & F["ever_w"]], parent=n_ew)
        out.q(b, "first_hours", F["first"][M & F["ever_w"]], equal_to=n_ew)
        if not detail:
            continue
        sel = M[F["int_e"]]
        d = F["int_d"][sel]
        nzi = out.q(b, "interval_hours", d[d > 0])
        out.cell(b + ("n_zero_intervals",), int((d == 0).sum()))
        out.cell(b + ("n_intervals_crossing_gap",), int(F["int_cross"][sel].sum()), parent=nzi)
        out.q(b, "encounter_median_interval_hours", F["int_med"][M], parent=n_ew)
        opn = np.arange(N_BINS)[None, :] < co.nopen[:, None]
        occ = F["occ"]
        Mo = opn & M[:, None]
        n_open = out.cell(b + ("n_open_bins",), int(Mo.sum()))
        parts = [out.cell(b + ("n_open_bins_empty",), int((Mo & (occ == 0)).sum())),
                 out.cell(b + ("n_open_bins_1",), int((Mo & (occ == 1)).sum())),
                 out.cell(b + ("n_open_bins_2",), int((Mo & (occ == 2)).sum())),
                 out.cell(b + ("n_open_bins_3plus",), int((Mo & (occ >= 3)).sum()))]
        out.ident(n_open, parts)
        Mg = Mo & ~co.gap_bins
        n_og = out.cell(b + ("n_open_bins_excl_gap",), int(Mg.sum()))
        out.cell(b + ("n_open_bins_excl_gap_empty",), int((Mg & (occ == 0)).sum()), parent=n_og)
        if s == "train":
            for k in range(N_BINS):
                full = M & (co.out_h >= BIN_H * (k + 1))
                out.cell(b + ("occupancy_curve", str(k), "n_occupied"), int((full & (occ[:, k] > 0)).sum()), parent=dens.bin_full[s][k])
    if not detail:
        return
    opn = np.arange(N_BINS)[None, :] < co.nopen[:, None]
    for var, arr in strata.items():
        for s in SPLITS:
            for v, dc in dens.strata[(var, s)].items():
                Mv = co.masks[s] & (arr == v)
                b = ("strata", var, s, str(v)) + base[1:]
                out.cell(b + ("n_ever_layer_b",), (Mv & F["ever_b"]).sum(), parent=dc["enc"])
                out.cell(b + ("n_ever_336h",), (Mv & F["ever_w"]).sum(), parent=dc["enc"])
                out.cell(b + ("n_avail_full_72h",), (Mv & co.elig[STRATA_LANDMARK] & F["avail"][1.0][("full", STRATA_LANDMARK)]).sum(),
                         parent=dc["elig"])
                G6 = Mv & (co.W >= BIN_H)
                out.q(b, "rate_per_day_all", (F["cnt"] / np.maximum(co.W, 1e-9) * 24.0)[G6], equal_to=dc["w6"])
                Mo = opn & Mv[:, None]
                no = out.cell(b + ("n_open_bins",), int(Mo.sum()))
                out.cell(b + ("n_open_bins_empty",), int((Mo & (F["occ"] == 0)).sum()), parent=no)
    if outcome_on:
        M = co.masks["train"]
        rate = F["cnt"] / np.maximum(co.W, 1e-9) * 24.0
        for L in LANDMARKS:
            cells, g = dens.outcome[("outcome1", L)]
            for grp, den in cells.items():
                E = M & co.elig[L] & (g == grp)
                b = ("outcome_descriptors",) + base[1:] + ("outcome1", f"{L}h", f"group_{grp}")
                out.cell(b + ("n_ever_layer_b",), (E & F["ever_b"]).sum(), parent=den)
                out.cell(b + ("n_ever_336h",), (E & F["ever_w"]).sum(), parent=den)
                out.cell(b + ("n_avail_full",), (E & F["avail"][1.0][("full", L)]).sum(), parent=den)
                out.q(b, "rate_per_day_all", rate[E], equal_to=den)
        for lab, val in (("died", True), ("survived", False)):
            mc = dens.outcome["mortality"][lab]
            Md = M & (died == val)
            b = ("outcome_descriptors",) + base[1:] + ("mortality", lab)
            out.cell(b + ("n_ever_layer_b",), (Md & F["ever_b"]).sum(), parent=mc["enc"])
            out.cell(b + ("n_ever_336h",), (Md & F["ever_w"]).sum(), parent=mc["enc"])
            out.q(b, "rate_per_day_all", rate[Md & (co.W >= BIN_H)], equal_to=mc["w6"])
            for L in LANDMARKS:
                out.cell(b + (f"{L}h", "n_avail_full"), (Md & co.elig[L] & F["avail"][1.0][("full", L)]).sum(), parent=mc[L])


def provisional_class(co: Cohort, F: dict, triggered: bool, score_input: bool, is_vital: bool) -> dict:
    """Spec 9.1-9.3 approved defaults on TRAIN raw values (labels only; Vinai applies the rules)."""
    M = co.masks["train"]
    n = int(M.sum())
    ever = int((M & F["ever_w"]).sum())
    p = ever / n if n else 0.0
    E = M & co.elig[72]
    A_ = E & F["avail"][1.0][("full", 72)]
    occ = F["occ"][:, :12]
    nb = int(A_.sum()) * 12
    empty_share = float((occ[A_] == 0).sum() / nb) if nb else None
    ages = F["last_age"][72][A_]
    med_age = float(np.median(ages)) if len(ages) else None
    if p < 0.01 or ever < DETAIL_MIN_TRAIN:
        cls = "drop"
    elif p < 0.10 or (med_age is not None and med_age > 72):
        cls = "measured_indicator_only"
    elif empty_share is not None and (1 - empty_share) >= 0.20:
        cls = "value_mask_time_since_last"
    else:
        cls = "measured_indicator_only"
    if cls != "drop" and triggered:
        cls = "triggered_S9_form"
    elif cls == "measured_indicator_only" and score_input:
        cls = "value_mask_score_input_override"
    sparsity = None
    if empty_share is not None and cls != "drop":
        sparsity = "dense" if empty_share <= 0.20 else ("moderate" if empty_share <= 0.80 else "sparse")
    res = {"class": cls, "sparsity_class_72h": sparsity, "score_input": score_input, "triggered": triggered}
    if is_vital and cls != "drop":
        Mo = (np.arange(N_BINS)[None, :] < co.nopen[:, None]) & M[:, None] & (F["occ"] > 0)
        share3 = float((F["occ"][Mo] >= 3).mean()) if Mo.any() else 0.0
        d = F["int_d"][M[F["int_e"]]]
        d = d[d > 0]
        medi = float(np.median(d)) if len(d) else math.inf
        res["bin_width_6h_adequate"] = bool(share3 >= 0.5 and medi <= 2.0)
    return res


# =============================================================================
# Metric 9: co-measurement and panel structure (train only, aggregated)
# =============================================================================
PANEL_MIN = 100
PANEL_MAX = 30
VITAL_CO = ("heart_rate", "respiratory_rate", "sbp", "spo2", "temperature")
CO_BINS = 12


def _pair_block(out: Out, base: tuple, X: np.ndarray, names: list[str], unit: str) -> None:
    """X: bool matrix units x channels. Cells: total units, units with A, units with A and B (subset of A);
    P(B | A), P(A | B) as ratios and the Jaccard index as a dependent value."""
    tot = out.cell(base + (f"n_{unit}",), X.shape[0])
    cA = {}
    for j, a in enumerate(names):
        cA[a] = out.cell(base + ("channels", a, f"n_{unit}_with_channel"), int(X[:, j].sum()), parent=tot)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            nab = int((X[:, i] & X[:, j]).sum())
            pb = base + ("pairs", f"{a}__{b}")
            c = out.cell(pb + (f"n_{unit}_both",), nab, parent=cA[a])
            out.reg.subset(c, cA[b], "pair_b")
            out.pct(pb, f"b_given_a", c, cA[a])
            out.pct(pb, f"a_given_b", c, cA[b])
            na, nb = int(X[:, i].sum()), int(X[:, j].sum())
            den = na + nb - nab
            out.dep(pb, "jaccard", [c, cA[a], cA[b]], round(nab / den, 3) if den else None)


def write_comeasurement(out: Out, co: Cohort, vb: ValueBuild, cm: dict, occ12: dict[str, np.ndarray],
                        detail: set[str]) -> None:
    """Spec 7 metric 9 on train: lab draws (channels sharing encounter and CHARTTIME), the most frequent exact draw
    panels (support >= 100 train draws), vitals same-CHARTTIME co-recording, and bin-level co-occupancy pooled over
    bins 0-11. Only channels with >= 100 train encounters (detail) enter."""
    tr = co.masks["train"]

    def unit_matrix(chans: list[str]) -> tuple[np.ndarray, list[str], pd.DataFrame]:
        parts = []
        for ch in chans:
            ev = vb.events[ch]
            m = tr[ev["enc_i"].to_numpy()] & (ev["t"].to_numpy() >= 0) & (ev["t"].to_numpy() < co.W[ev["enc_i"].to_numpy()])
            parts.append(ev.loc[m, ["enc_i", "tc"]].drop_duplicates().assign(ch=ch))
        if not parts:
            return np.zeros((0, 0), dtype=bool), [], pd.DataFrame()
        d = pd.concat(parts, ignore_index=True)
        if not len(d):
            return np.zeros((0, len(chans)), dtype=bool), chans, d
        piv = pd.crosstab([d["enc_i"], d["tc"]], d["ch"]).reindex(columns=chans, fill_value=0) > 0
        return piv.to_numpy(), chans, piv

    labs = sorted(ch for ch, c in cm["channels"].items() if c["source"] == "lab" and ch in detail)
    X, names, piv = unit_matrix(labs)
    base = ("comeasurement", "train", "lab_draws")
    if len(names) >= 2:
        _pair_block(out, base, X, names, "draws")
        tot = out.A.P(*base, "n_draws")
        keys = pd.Series(["+".join(n for n, x in zip(names, row) if x) for row in X], dtype=object)
        vc = keys.value_counts()
        vc = vc[(vc >= PANEL_MIN) & vc.index.to_series().str.contains("+", regex=False).to_numpy()]
        vc = vc.sort_values(ascending=False, kind="mergesort")
        for k in sorted(vc.index[:PANEL_MAX]):
            out.cell(base + ("panels", k, "n_draws_exact_panel"), int(vc[k]), parent=tot)
        out.val(base + ("panel_rule",), f"exact channel set of a draw; sets of >= 2 channels with >= {PANEL_MIN} train draws; "
                                         f"at most {PANEL_MAX} (largest), listed by name")
    vit = [ch for ch in VITAL_CO if ch in detail]
    Xv, nv, _ = unit_matrix(vit)
    if len(nv) >= 2:
        _pair_block(out, ("comeasurement", "train", "vitals_same_charttime"), Xv, nv, "recordings")
    chans = sorted(occ12)
    if len(chans) >= 2:
        open_ = (np.arange(CO_BINS)[None, :] < co.nopen[:, None]) & tr[:, None]
        Xb = np.stack([occ12[c][open_] for c in chans], axis=1)
        _pair_block(out, ("comeasurement", "train", "bin_cooccupancy_bins_0_11"), Xb, chans, "open_bins")



# =============================================================================
# Medications: S3 state machine (plan v12 section 5; spec 7.10, 9.4; D1-D5, D14, D16, D18)
# =============================================================================
def group_cutoff(group: str) -> float:
    return VASO_CUTOFF_H if group == "vasoactive" else DEFAULT_CUTOFF_H


def resolve_drug_groups(docs: pd.DataFrame, dmap: dict, neuro_enc: np.ndarray) -> pd.DataFrame:
    """Primary group per row (osmotic agents resolved by the encounter's admitting diagnosis, D3) and the sensitivity
    group (all osmotic agents -> osmotherapy)."""
    info = docs["key"].map(dmap)
    base = info.map(lambda d: d["group"] if isinstance(d, dict) else "unmapped")
    other = info.map(lambda d: d["otherwise"] if isinstance(d, dict) else None)
    sub = info.map(lambda d: d["subgroup"] if isinstance(d, dict) else "")
    neuro = neuro_enc[docs["enc_i"].to_numpy()]
    osm = (base == "osmotic_conditional").to_numpy()
    prim = np.where(osm, np.where(neuro, "osmotherapy", other.fillna("").to_numpy()), base.to_numpy())
    sens = np.where(osm, "osmotherapy", base.to_numpy())
    return docs.assign(base_group=base.to_numpy(), group=prim, group_sens=sens, subgroup=sub.to_numpy(), osmotic=osm)


def build_courses(enc: np.ndarray, key: np.ndarray, S: np.ndarray, E: np.ndarray, cutoff: np.ndarray) -> pd.DataFrame:
    """S3 courses per (encounter, course key). Inputs sorted by (enc, key, S). A documentation joins the open course
    if S - last <= cutoff and no ENDDATE (after its own start) of the course is at or before S; otherwise a new course
    starts. Course stop = min(last + cutoff, earliest valid ENDDATE)."""
    n = len(S)
    rs, rl, re_, rc, rnd, re_i, rk = [], [], [], [], [], [], []
    i = 0
    while i < n:
        ge, gk, c = enc[i], key[i], cutoff[i]
        start = last = S[i]
        endE = E[i] if (E[i] == E[i] and E[i] > S[i]) else math.inf
        nd = 1
        j = i + 1
        while j < n and enc[j] == ge and key[j] == gk:
            sj = S[j]
            if sj - last <= c and sj < endE:
                last = sj
                ej = E[j]
                if ej == ej and ej > sj and ej < endE:
                    endE = ej
                nd += 1
            else:
                rs.append(start); rl.append(last); re_.append(endE); rc.append(c); rnd.append(nd); re_i.append(ge); rk.append(gk)
                start = last = sj
                ej = E[j]
                endE = ej if (ej == ej and ej > sj) else math.inf
                nd = 1
            j += 1
        rs.append(start); rl.append(last); re_.append(endE); rc.append(c); rnd.append(nd); re_i.append(ge); rk.append(gk)
        i = j
    cr = pd.DataFrame({"enc_i": np.array(re_i, dtype=np.int64), "key": np.array(rk, dtype=object),
                       "start": np.array(rs, dtype=float), "last": np.array(rl, dtype=float),
                       "endE": np.array(re_, dtype=float), "cutoff": np.array(rc, dtype=float), "n_docs": np.array(rnd)})
    cr["stop"] = np.minimum(cr["last"] + cr["cutoff"], cr["endE"])
    cr["end_type"] = np.where(cr["endE"] <= cr["last"] + cr["cutoff"], "enddate", "cutoff")
    cr["has_enddate"] = np.isfinite(cr["endE"])
    return cr


def courses_for(docs: pd.DataFrame, group_col: str, factor: float, as_of: float | None = None) -> pd.DataFrame:
    d = docs[docs[group_col].isin(MED_GROUPS)]
    if as_of is not None:  # as-of rebuild: only documentation at or before t; a future ENDDATE is unknown
        d = d[d["S"] <= as_of].copy()
        d.loc[d["E"] > as_of, "E"] = np.nan
    d = d.sort_values(["enc_i", "key", "S"], kind="mergesort")
    cut = d[group_col].map(group_cutoff).to_numpy(dtype=float) * factor
    cr = build_courses(d["enc_i"].to_numpy(), d["key"].to_numpy(), d["S"].to_numpy(dtype=float),
                       d["E"].to_numpy(dtype=float), cut)
    gk = d.drop_duplicates(["enc_i", "key"])[["enc_i", "key", group_col, "subgroup"]].rename(columns={group_col: "group"})
    gk = gk.astype({"enc_i": np.int64})
    cr = cr.merge(gk, on=["enc_i", "key"], how="left", validate="many_to_one")
    return cr


def active_at(cr: pd.DataFrame, L: float) -> np.ndarray:
    """active(t) = 1 if t - last doc <= cutoff and no ENDDATE at or before t (course-level evaluation)."""
    return ((cr["start"] <= L) & (L <= cr["last"] + cr["cutoff"]) & (L < cr["endE"])).to_numpy()


def bin_activity(cr: pd.DataFrame, n_enc: int, groups: Iterable[str], nbins: int) -> dict:
    """Layer C: group active in bin k if a course overlaps [6k, 6k + 6) (start < 6k + 6 and stop > 6k)."""
    out = {}
    for g in groups:
        c = cr[cr["group"] == g]
        m = np.zeros((n_enc, nbins), dtype=bool)
        for k in range(nbins):
            sel = (c["start"].to_numpy() < BIN_H * (k + 1)) & (c["stop"].to_numpy() > BIN_H * k)
            m[c["enc_i"].to_numpy()[sel], k] = True
        out[g] = m
    return out


def as_of_identity_check(docs: pd.DataFrame, n_enc: int, as_of: float = 48.0) -> int:
    """Acceptance 5: rebuilding with only documentation <= 48 h must give the same bins 0-7 as the full build."""
    nb = int(as_of // BIN_H)
    full = bin_activity(courses_for(docs, "group", 1.0), n_enc, MED_GROUPS, nb)
    trunc = bin_activity(courses_for(docs, "group", 1.0, as_of=as_of), n_enc, MED_GROUPS, nb)
    return int(sum((full[g] != trunc[g]).sum() for g in MED_GROUPS))


def max_concurrent(cr: pd.DataFrame, n_enc: int, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Per encounter, the maximum number of distinct course keys active at once inside [lo, hi)."""
    if not len(cr):
        return np.zeros(n_enc, dtype=int)
    e = cr["enc_i"].to_numpy()
    a = np.maximum(cr["start"].to_numpy(), lo[e])
    b = np.minimum(np.minimum(cr["stop"].to_numpy(), cr["last"].to_numpy() + cr["cutoff"].to_numpy()), hi[e])
    ok = b > a
    ev = pd.DataFrame({"e": np.r_[e[ok], e[ok]], "t": np.r_[a[ok], b[ok]],
                       "d": np.r_[np.ones(ok.sum(), dtype=int), -np.ones(ok.sum(), dtype=int)]})
    ev = ev.sort_values(["e", "t", "d"], kind="mergesort")
    ev["c"] = ev.groupby("e")["d"].cumsum()
    m = ev.groupby("e")["c"].max()
    res = np.zeros(n_enc, dtype=int)
    res[m.index.to_numpy()] = m.to_numpy()
    return res


def write_medications(out: Out, co: Cohort, dens: Dens, docs: pd.DataFrame, strata: dict, outcome_on: bool,
                      died: np.ndarray, local: dict) -> dict:
    res: dict[str, Any] = {"provisional": {}}
    insw = (docs["S"].to_numpy() >= 0) & (docs["S"].to_numpy() < co.W[docs["enc_i"].to_numpy()])
    crs = {f: courses_for(docs, "group", f) for f in CUTOFF_FACTORS}
    cr1 = crs[1.0]
    cr_sens = courses_for(docs, "group_sens", 1.0)
    for g in MED_GROUPS:
        dg = docs[(docs["group"] == g).to_numpy()]
        F = event_features(co, dg["enc_i"].to_numpy(), dg["S"].to_numpy())
        detail = int((co.masks["train"] & F["ever_w"]).sum()) >= DETAIL_MIN_TRAIN
        write_event_metrics(out, co, dens, ("drug_groups", g), F, detail=detail, cap_h=group_cutoff(g), lab_offsets=False,
                            strata=strata, outcome_on=outcome_on, died=died)
        res["provisional"][g] = {"class": "drop" if (int((co.masks["train"] & F["ever_w"]).sum()) < max(DETAIL_MIN_TRAIN, 0.01 * co.masks["train"].sum())) else "keep",
                                 "cutoff_hours": group_cutoff(g)}
        cg = {f: c[c["group"] == g] for f, c in crs.items()}
        flags = {}
        share_end = {}
        for s in SPLITS:
            M = co.masks[s]
            b = ("drug_groups", g, s, "medication")
            # documentation intervals within the same course key (spec 7.10)
            di = dg[insw[(docs["group"] == g).to_numpy()]]
            di = di[M[di["enc_i"].to_numpy()]].sort_values(["enc_i", "key", "S"], kind="mergesort")
            same = (di["enc_i"].to_numpy()[1:] == di["enc_i"].to_numpy()[:-1]) & (di["key"].to_numpy()[1:] == di["key"].to_numpy()[:-1])
            d = np.diff(di["S"].to_numpy())[same] if len(di) > 1 else np.array([])
            ni = out.q(b, "doc_interval_hours", d)
            for f in CUTOFF_FACTORS:
                out.cell(b + (f"n_doc_intervals_gt_cutoff_x{f:g}",), int((d > group_cutoff(g) * f).sum()), parent=ni)
            # ENDDATE per group (rows documented in [t0, W))
            rr = dg[insw[(docs["group"] == g).to_numpy()] & M[dg["enc_i"].to_numpy()]]
            S_, E_ = rr["S"].to_numpy(), rr["E"].to_numpy()
            tot = out.cell(b + ("enddate", "n_rows"), len(rr))
            parts = {"missing": int(np.isnan(E_).sum()), "after_start": int((E_ > S_).sum()),
                     "at_or_before_start": int((E_ <= S_).sum())}
            pc = out.partition(b + ("enddate",), parts, total=tot)
            share_end[s] = parts["after_start"] / len(rr) if len(rr) else None
            out.cell(b + ("enddate", "n_date_only_start"), int(rr["date_only"].sum()), parent=tot)
            out.cell(b + ("enddate", "n_midnight_clock_start"), int(rr["midnight"].sum()), parent=tot)
            # courses (start in [t0, W))
            for f, c in cg.items():
                cs = c[M[c["enc_i"].to_numpy()] & (c["start"].to_numpy() >= 0) & (c["start"].to_numpy() < co.W[c["enc_i"].to_numpy()])]
                bf = b + ("courses", f"cutoff_x{f:g}")
                nc = out.cell(bf + ("n_courses",), len(cs))
                dur = (cs["stop"] - cs["start"]).to_numpy()
                ce = out.q(bf, "duration_hours_cutoff_end", dur[cs["end_type"].to_numpy() == "cutoff"])
                ee = out.q(bf, "duration_hours_enddate_end", dur[cs["end_type"].to_numpy() == "enddate"])
                out.ident(nc, [ce, ee])
                if f != 1.0:
                    continue
                we = cs[cs["has_enddate"].to_numpy()]
                diff = (we["last"] + we["cutoff"] - we["endE"]).to_numpy()
                nw = out.q(bf, "cutoff_stop_minus_enddate_hours", diff, parent=nc)
                out.partition(bf + ("cutoff_vs_enddate",), {"cutoff_earlier": int((diff < -ENDDATE_TOL_H).sum()),
                                                            "cutoff_later": int((diff > ENDDATE_TOL_H).sum()),
                                                            "within_6h": int((np.abs(diff) <= ENDDATE_TOL_H).sum())}, total=nw)
                n24 = out.cell(bf + ("n_abs_difference_gt_24h",), int((np.abs(diff) > 24).sum()), parent=nw)
                if s == "train":
                    flags["enddate_differs_gt24h_in_gt25pct"] = (bool((np.abs(diff) > 24).mean() > 0.25) if len(diff) >= MIN_CELL else None, [nw, n24])
                for L in LANDMARKS:
                    E = M & co.elig[L]
                    ce_ = c[E[c["enc_i"].to_numpy()]]
                    ov = np.clip(np.minimum(ce_["stop"].to_numpy(), L) - np.maximum(ce_["start"].to_numpy(), 0.0), 0, None)
                    nov = out.cell(bf + ("landmarks", f"{L}h", "n_courses_overlapping_prefix"), int((ov > 0).sum()))
                    tot_h = float(ov.sum())
                    share = float(ov[~ce_["has_enddate"].to_numpy()].sum() / tot_h) if tot_h > 0 else None
                    out.dep(bf + ("landmarks", f"{L}h"), "share_active_hours_without_enddate", [nov], share)
            # active prevalence at landmarks (primary and the 0.5x / 2x cutoffs) and at t0 (pre-ICU)
            for f, c in cg.items():
                for L in LANDMARKS:
                    act = active_at(c, L)
                    a = np.zeros(co.n, dtype=bool)
                    a[c["enc_i"].to_numpy()[act]] = True
                    out.cell(b + ("active_prevalence", f"cutoff_x{f:g}", f"{L}h", "n_active"), (M & co.elig[L] & a).sum(),
                             parent=dens.elig[s][L])
            # metric 10 by stratum (decision 7): active prevalence at the strata landmark, cutoff x1
            c1 = cg[1.0]
            a72 = np.zeros(co.n, dtype=bool)
            a72[c1["enc_i"].to_numpy()[active_at(c1, STRATA_LANDMARK)]] = True
            tot72 = out.A.P(*b, "active_prevalence", "cutoff_x1", f"{STRATA_LANDMARK}h", "n_active")
            for var, arr in strata.items():
                cells = []
                for v in dens.strata_levels[var]:
                    cells.append(out.cell(("drug_groups", g, s, "strata", var, str(v), "n_active_72h"),
                                          (M & co.elig[STRATA_LANDMARK] & a72 & (arr == v)).sum(),
                                          parent=dens.strata[(var, s)][v]["elig"]))
                out.ident(tot72, cells)
            a0 = np.zeros(co.n, dtype=bool)
            a0[cg[1.0]["enc_i"].to_numpy()[active_at(cg[1.0], 0.0)]] = True
            out.cell(b + ("n_active_at_t0",), (M & a0).sum(), parent=dens.enc[s])
        vals = [v for v in share_end.values() if v is not None]
        flags["enddate_presence_differs_gt10pp_between_splits"] = (
            bool(max(vals) - min(vals) > 0.10) if len(vals) == len(SPLITS) else None,
            [out.A.P("drug_groups", g, s, "medication", "enddate", "n_rows") for s in SPLITS])
        for fk, (v, needs) in flags.items():
            out.dep(("drug_groups", g, "flags_9_4"), fk, needs, v)
    # vasoactive count (D1), Layer B max, subgroup sensitivity
    va = cr1[cr1["group"] == "vasoactive"]
    lo = np.zeros(co.n)
    hi = np.minimum(co.out_h, LAYER_B_H)
    mxb = max_concurrent(va, co.n, lo, hi)
    for s in SPLITS:
        M = co.masks[s]
        b = ("vasoactive_count", s)
        for L in LANDMARKS:
            act = active_at(va, L)
            cnt = np.bincount(va["enc_i"].to_numpy()[act], minlength=co.n)
            E = M & co.elig[L]
            out.partition(b + (f"{L}h",), {"0": int((E & (cnt == 0)).sum()), "1": int((E & (cnt == 1)).sum()),
                                           "2": int((E & (cnt == 2)).sum()), "3plus": int((E & (cnt >= 3)).sum())},
                          total=dens.elig[s][L])
            for sg in ("vasopressor", "inotrope"):
                a = np.zeros(co.n, dtype=bool)
                a[va["enc_i"].to_numpy()[act & (va["subgroup"].to_numpy() == sg)]] = True
                out.cell(b + (f"{L}h", "subgroup_sensitivity", f"n_{sg}_active"), (E & a).sum(), parent=dens.elig[s][L])
        out.partition(b + ("layer_b_max",), {"0": int((M & (mxb == 0)).sum()), "1": int((M & (mxb == 1)).sum()),
                                             "2": int((M & (mxb == 2)).sum()), "3plus": int((M & (mxb >= 3)).sum())},
                      total=dens.enc[s])
    # osmotic agents (D3): courses per resolved class, primary and sensitivity
    osm_keys = set(docs.loc[docs["osmotic"], "key"])
    for s in SPLITS:
        M = co.masks[s]
        for lab, c in (("primary", cr1), ("sensitivity_all_osmotherapy", cr_sens)):
            cs = c[c["key"].isin(osm_keys) & M[c["enc_i"].to_numpy()] & (c["start"].to_numpy() >= 0)
                   & (c["start"].to_numpy() < co.W[c["enc_i"].to_numpy()])]
            out.partition(("osmotic_resolution", s, lab), {g: int((cs["group"] == g).sum())
                                                           for g in ("osmotherapy", "diuretic", "electrolyte_replacement")})
    # unmapped / none documentation (counted, not analysed)
    for s in SPLITS:
        M = co.masks[s]
        dd = docs[insw & M[docs["enc_i"].to_numpy()]]
        bg = dd["base_group"].to_numpy()
        out.partition(("prescription_groups_in_window", s), {
            "analysed_groups": int(np.isin(dd["group"].to_numpy(), MED_GROUPS).sum()),
            "none": int((bg == "none").sum()), "unmapped": int((bg == "unmapped").sum())})
    # acceptance 5: as-of rebuild at 48 h equals the full build for bins 0-7
    mism = as_of_identity_check(docs, co.n, 48.0)
    if mism:
        raise AuditError(f"S3 as-of rebuild at 48 h differs from the full build in {mism} encounter-group-bins (bins 0-7)")
    res["as_of_48h_identity"] = "ok"
    local["medication_courses"] = {f"cutoff_x{f:g}": int(len(c)) for f, c in crs.items()}
    return res


# =============================================================================
# Microbiology: M1 collections (base_code only, D12) and M3 availability
# =============================================================================
def micro_collections(ev: pd.DataFrame, smap: dict, omap: dict) -> tuple[pd.DataFrame, dict]:
    si = ev["spec_key"].map(smap)
    mapped = si.notna().to_numpy()
    base = np.where(mapped, si.map(lambda d: d["base"] if isinstance(d, dict) else ""),
                    ev["spec_itemid"].astype(str).str.split("+").str[0].to_numpy())
    group = np.where(mapped, si.map(lambda d: d["group"] if isinstance(d, dict) else ""), "unmapped")
    addons = si.map(lambda d: bool(d["addons"]) if isinstance(d, dict) else False).to_numpy(dtype=bool)
    oi = ev["org_key"].map(omap)
    omapped = oi.notna().to_numpy()
    types = oi.map(lambda d: d["types"] if isinstance(d, dict) else ())
    rows = pd.DataFrame({"enc_i": ev["enc_i"].to_numpy(), "th": ev["th"].to_numpy(), "base": base, "group": group,
                         "addons": addons, "gram": types.map(lambda t: "gram_class" in t).to_numpy(dtype=bool),
                         "ident": types.map(lambda t: "identification" in t).to_numpy(dtype=bool),
                         "neg": types.map(lambda t: "negative_final" in t).to_numpy(dtype=bool)})
    col = rows.groupby(["enc_i", "th", "base"], sort=True).agg(group=("group", "first"), addons=("addons", "max"),
                                                                gram=("gram", "max"), ident=("ident", "max"),
                                                                neg=("neg", "max")).reset_index()
    col["susc"] = col["ident"] & col["addons"]
    return col, {"rows_spec_unmapped": int((~mapped).sum()), "rows_org_unmapped": int((~omapped).sum()),
                 "rows": len(ev), "spec_unmapped_mask": ~mapped, "org_unmapped_mask": ~omapped}


M3_TYPES = (("gram_class", "gram"), ("identification", "ident"), ("susceptibility", "susc"), ("negative_final", "neg"))


def write_micro(out: Out, co: Cohort, dens: Dens, ev: pd.DataFrame, col: pd.DataFrame, minfo: dict, strata: dict,
                outcome_on: bool, died: np.ndarray, fill: dict) -> dict:
    prov = {}
    for g in SPECIMEN_GROUPS_V1 + ("unmapped",):
        c = col[col["group"] == g]
        F = event_features(co, c["enc_i"].to_numpy(), c["th"].to_numpy())
        detail = int((co.masks["train"] & F["ever_w"]).sum()) >= DETAIL_MIN_TRAIN
        write_event_metrics(out, co, dens, ("specimen_groups", g), F, detail=detail, cap_h=24.0, lab_offsets=False,
                            strata=strata, outcome_on=outcome_on, died=died)
        prov[g] = {"detail": detail}
    e = col["enc_i"].to_numpy()
    th = col["th"].to_numpy()
    for s in SPLITS:
        M = co.masks[s][e]
        b = ("microbiology", s)
        out.cell(b + ("n_collections_in_window",), int((M & (th >= 0) & (th < co.W[e])).sum()))
        out.cell(b + ("n_collections_pre_icu_6h",), int((M & (th >= -6) & (th < 0)).sum()))
        out.cell(b + ("n_collections_pre_icu_24h",), int((M & (th >= -24) & (th < 0)).sum()))
        for L in LANDMARKS:
            E = co.masks[s] & co.elig[L]
            inL = E[e] & (th >= 0) & (th < L)
            for f in MICRO_GRID:
                bl = b + ("m3", f"{L}h", f"grid_x{f:g}")
                pend = np.zeros(len(col), dtype=bool)
                for tname, colname in M3_TYPES:
                    has = col[colname].to_numpy(dtype=bool)
                    av = th + MICRO_TA_H[tname] * f < L
                    nc = out.cell(bl + (tname, "n_collections"), int((inL & has).sum()))
                    out.cell(bl + (tname, "n_available"), int((inL & has & av).sum()), parent=nc)
                    pend |= inL & has & ~av
                differs = np.zeros(co.n, dtype=bool)
                differs[e[pend]] = True
                out.cell(bl + ("n_encounters_m3_differs_from_collection_time_results",), int((E & differs).sum()),
                         parent=dens.elig[s][L])
        # Layer B window (all encounters, window end min(OUT, 24 h))
        bend = np.minimum(co.out_h, LAYER_B_H)
        inB = co.masks[s][e] & (th >= 0) & (th < bend[e])
        for f in MICRO_GRID:
            pend = np.zeros(len(col), dtype=bool)
            for tname, colname in M3_TYPES:
                pend |= inB & col[colname].to_numpy(dtype=bool) & (th + MICRO_TA_H[tname] * f >= bend[e])
            differs = np.zeros(co.n, dtype=bool)
            differs[e[pend]] = True
            out.cell(b + ("m3_layer_b", f"grid_x{f:g}", "n_encounters_m3_differs_from_collection_time_results"),
                     int((co.masks[s] & differs).sum()), parent=dens.enc[s])
        tr = co.masks[s][ev["enc_i"].to_numpy()] & (ev["th"].to_numpy() >= 0) & (ev["th"].to_numpy() < co.W[ev["enc_i"].to_numpy()])
        nr = out.cell(b + ("n_rows_in_window",), int(tr.sum()))
        out.cell(b + ("n_rows_specimen_key_unmapped",), int((tr & minfo["spec_unmapped_mask"]).sum()), parent=nr)
        out.cell(b + ("n_rows_organism_key_unmapped",), int((tr & minfo["org_unmapped_mask"]).sum()), parent=nr)
    # D12: MICROBIOLOGYEVENTS column fill rates (all rows, and train rows in stay)
    tot_all = out.cell(("map_checks", "micro_fill_rates", "n_rows_all"), fill["rows_all"])
    tot_tr = out.cell(("map_checks", "micro_fill_rates", "n_rows_in_stay_train"), fill["rows_in_stay_train"])
    for cname, v in fill["columns"].items():
        b = ("map_checks", "micro_fill_rates", "columns", cname)
        a = out.cell(b + ("n_nonblank_all",), v["all"], parent=tot_all)
        t = out.cell(b + ("n_nonblank_in_stay_train",), v["in_stay_train"], parent=tot_tr)
        out.pct(b, "nonblank_all", a, tot_all)
        out.pct(b, "nonblank_in_stay_train", t, tot_tr)
    return prov


# =============================================================================
# Surgery (v1.7 linkage, A2.1): Layer C events, Layer A recent and remote
# =============================================================================
def surgery_frame(cfg: S2Config, A, enc: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    acfg = A.AuditConfig(canonical_path=cfg.canonical_path, harmonization_path=cfg.harmonization_path,
                         map_path=cfg.complexity_map_path, split_path=cfg.split_path,
                         surgery_info_path=raw_path(cfg, "SURGERY_INFO"), out_dir=cfg.out_dir,
                         admissions_path=raw_path(cfg, "ADMISSIONS"), harm_t0_col=cfg.harm_t0_col,
                         harm_out_col=cfg.harm_out_col, split_subject_col=cfg.split_subject_col, split_col=cfg.split_col,
                         expected_surgery_rows=cfg.expected_surgery_rows, expected_encounters=cfg.expected_encounters,
                         expected_subjects=cfg.expected_subjects, expected_split_encounters=dict(cfg.expected_split_encounters),
                         expected_map_rows=cfg.expected_complexity_map_rows, verbose=False)
    if cfg.surgery_expected_missing is not None:
        acfg.expected_missing = dict(cfg.surgery_expected_missing)
    acfg.pinned_sha256 = dict(cfg.pinned_sha256)
    inp = A.load_inputs(acfg)
    s, meta = A.build_frame(inp)
    idx = enc.set_index("encounter_id")
    s["enc_i"] = s["encounter_id"].map(idx["enc_i"])
    s["T_h"] = np.where(s["enc_i"].notna(), (s["T"] - s["t0"]).dt.total_seconds() / 3600.0, np.nan)
    s["B_h"] = np.where(s["enc_i"].notna(), (s["B"] - s["t0"]).dt.total_seconds() / 3600.0, np.nan)
    return s, {"systems": meta["systems"], "link_class_counts": s["link_class"].value_counts().to_dict(),
               "link_classes": list(A.LINK_CLASSES)}


def write_surgery(out: Out, co: Cohort, dens: Dens, s: pd.DataFrame, sinfo: dict, strata: dict, outcome_on: bool,
                  died: np.ndarray) -> dict:
    lc = {c: int((s["link_class"] == c).sum()) for c in sinfo["link_classes"]}
    if sum(lc.values()) != len(s):
        raise AuditError("surgery link classes do not sum to the SURGERY_INFO rows")
    out.partition(("surgery", "link_classes"), lc)
    lk = s[s["link_class"] == "linked"].copy()
    lk["enc_i"] = lk["enc_i"].astype(int)
    c = lk[lk["timing"] == "during_icu_layer_c"]
    e = c["enc_i"].to_numpy()
    T = np.minimum(c["T_h"].to_numpy(dtype=float), co.out_h[e] - 1e-9)
    F = event_features(co, e, T)
    write_event_metrics(out, co, dens, ("surgery", "layer_c"), F, detail=True, cap_h=24.0, lab_offsets=False,
                        strata=strata, outcome_on=outcome_on, died=died)
    Bc = np.clip(c["B_h"].to_numpy(dtype=float), 0.0, co.out_h[e] - 1e-9)
    Fb = event_features(co, e, Bc)
    write_event_metrics(out, co, dens, ("surgery", "layer_c_begin_time_sensitivity"), Fb, detail=False, cap_h=24.0,
                        lab_offsets=False, strata={}, outcome_on=False, died=died)
    in_gap = np.zeros(len(c), dtype=bool)
    for j, (i, t) in enumerate(zip(e, c["T_h"].to_numpy(dtype=float))):
        for a, b in co.gaps.get(int(i), []):
            if a < t < b:
                in_gap[j] = True
    # reuse (A2.1): Layer C for one encounter and remote Layer A (T < t0 - 24 h, same HADM) for a later encounter
    eh = co.eh[["key", "enc_i", "t0"]]
    cc = pd.DataFrame({"j": np.arange(len(c)), "key": (c["SUBJECT_ID"] + "|" + c["HADM_ID"]).to_numpy(), "T": c["T"].to_numpy(),
                       "own": e})
    mm = cc.merge(eh, on="key", how="inner")
    reuse_j = set(mm.loc[(mm["enc_i"] != mm["own"]) & (mm["T"] < mm["t0"] - pd.Timedelta(hours=24)), "j"])
    reuse = np.isin(np.arange(len(c)), list(reuse_j))
    systems = sorted(set(sinfo["systems"]) | {"Unclassified", "__unmapped__"})
    for sp in SPLITS:
        M = co.masks[sp][e]
        b = ("surgery", "layer_c_extra", sp)
        tot = out.cell(b + ("n_surgeries",), int(M.sum()))
        out.partition(b + ("by_system",), {sy: int((M & (c["body_system"].to_numpy() == sy)).sum()) for sy in systems}, total=tot)
        out.cell(b + ("n_in_gap",), int((M & in_gap).sum()), parent=tot)
        out.cell(b + ("n_reused_as_remote_layer_a_of_later_encounter",), int((M & reuse).sum()), parent=tot)
        per = np.bincount(e[M], minlength=co.n)
        out.q(b, "surgeries_per_encounter_with_any", per[co.masks[sp] & (per > 0)].astype(float))
        # Layer A recent / remote
        bA = ("surgery", "layer_a", sp)
        rec = lk[lk["timing"] == "recent_layer_a"]
        rem = lk[lk["timing"] == "remote_same_hadm"]
        for lab, d in (("recent", rec), ("remote", rem)):
            m = np.zeros(co.n, dtype=bool)
            m[d["enc_i"].to_numpy()] = True
            out.cell(bA + (f"n_encounters_{lab}",), int((co.masks[sp] & m).sum()), parent=dens.enc[sp])
            ds = d[co.masks[sp][d["enc_i"].to_numpy()]]
            t2 = out.cell(bA + (lab, "n_surgeries"), len(ds))
            out.partition(bA + (lab, "by_system"), {sy: int((ds["body_system"].to_numpy() == sy).sum()) for sy in systems}, total=t2)
            if lab == "remote":
                out.cell(bA + (lab, "n_ending_before_admittime"), int(ds["pre_admit"].sum()), parent=t2)
    out.val(("surgery", "notes", "phi1_band"), "not computed by this script: the frozen medians are not an input (open question)")
    return {"recent_enc": set(lk.loc[lk["timing"] == "recent_layer_a", "enc_i"])}


# =============================================================================
# Layer A inputs (PATIENTS, ADMISSIONS) and the Layer A report (spec 6)
# =============================================================================
def read_patients(cfg: S2Config, A, enc: pd.DataFrame) -> pd.DataFrame:
    name = "PATIENTS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, gc, dc = (resolve(hdr, (c,), name) for c in ("SUBJECT_ID", "GENDER", "DOB"))
    d = A.read_csv_str(p, usecols=[sc, gc, dc])
    if len(d) != cfg.expected_raw_rows.get(name, len(d)):
        raise AuditError(f"{name}: {len(d)} rows, schema inventory says {cfg.expected_raw_rows.get(name)}")
    pt = pd.DataFrame({"SUBJECT_ID": vnorm_id(d[sc]), "sex": vclean_text(d[gc]).str.upper(),
                       "DOB_raw": d[dc], "DOB": fast_time(d[dc])})
    if pt["SUBJECT_ID"].duplicated().any():
        raise AuditError("PATIENTS SUBJECT_ID not unique")
    m = enc[["enc_i", "SUBJECT_ID", "t0"]].merge(pt, on="SUBJECT_ID", how="left", indicator=True)
    if len(m) != len(enc):
        raise AuditError("PATIENTS merge changed the encounter count")
    m = m.sort_values("enc_i")
    status, age_days, nbin = classify_ages(m["t0"], m["DOB"], m["DOB_raw"], m["_merge"].eq("both").to_numpy())
    n_invalid = int((status != "valid").sum())
    check_invalid_age_rate(n_invalid, len(m))
    sex = np.where(np.isin(m["sex"].fillna("").to_numpy(), ["M", "F"]), m["sex"].fillna("").to_numpy(), "unrecognized")
    return pd.DataFrame({"enc_i": m["enc_i"].to_numpy(), "sex": sex, "age_days": age_days,
                         "dob_status": status, "dob_neg_bin": nbin})


def read_admissions(cfg: S2Config, A, enc: pd.DataFrame) -> pd.DataFrame:
    """Per encounter: the earliest linked admission (min ADMITTIME over hadm_id_list): ADMITTIME, ICD10_CODE_CN
    (preliminary admitting diagnosis) and EDREGTIME."""
    name = "ADMISSIONS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc, at = (resolve(hdr, (c,), name) for c in ("SUBJECT_ID", "HADM_ID", "ADMITTIME"))
    dxc = resolve(hdr, ("ICD10_CODE_CN", "ICD10_CODE"), name)
    edc = resolve(hdr, ("EDREGTIME",), name, required=False)
    d = A.read_csv_str(p, usecols=[c for c in (sc, hc, at, dxc, edc) if c])
    if len(d) != cfg.expected_raw_rows.get(name, len(d)):
        raise AuditError(f"{name}: {len(d)} rows, schema inventory says {cfg.expected_raw_rows.get(name)}")
    ad = pd.DataFrame({"key": vnorm_id(d[sc]) + "|" + vnorm_id(d[hc]), "ADMITTIME": fast_time(d[at]),
                       "dx": d[dxc].astype(str), "ED": fast_time(d[edc]) if edc else pd.NaT,
                       "ed_raw": vclean_text(d[edc]) if edc else ""})
    eh = encounter_hadm_index(enc)[["key", "enc_i"]]
    m = eh.merge(ad, on="key", how="inner").sort_values(["enc_i", "ADMITTIME"], na_position="last")
    first = m.drop_duplicates("enc_i", keep="first").set_index("enc_i")
    out = pd.DataFrame({"enc_i": enc["enc_i"].to_numpy()})
    for c in ("ADMITTIME", "dx", "ED", "ed_raw"):
        out[c] = out["enc_i"].map(first[c])
    out["edc_present"] = bool(edc)
    return out


def level_support(out: Out, base: tuple, co: Cohort, lvl: np.ndarray, subj: np.ndarray, dens: Dens, per_level: bool = True) -> None:
    """Categorical level support (spec 6): levels, levels with >= m train subjects (m = 20, 10), train rows pooled to
    other (levels under 20 train subjects), validation/test rows whose level is unseen in train, per-level counts."""
    tr = co.masks["train"]
    subj_per = pd.Series(subj[tr]).groupby(lvl[tr]).nunique()
    out.val(base + ("levels_total_train",), int(len(subj_per)))
    for m in LEVEL_M:
        out.val(base + (f"levels_ge{m}_train_subjects",), int((subj_per >= m).sum()))
    keep = set(subj_per[subj_per >= LEVEL_M[0]].index)
    out.cell(base + ("n_train_rows_pooled_to_other_m20",), int((tr & ~np.isin(lvl, list(keep))).sum()), parent=dens.enc["train"])
    seen = set(subj_per.index)
    for s in ("validation", "test"):
        out.cell(base + (f"n_{s}_rows_level_unseen_in_train",), int((co.masks[s] & ~np.isin(lvl, list(seen))).sum()),
                 parent=dens.enc[s])
    if per_level:
        levels = sorted(pd.unique(lvl).tolist(), key=str)
        for s in SPLITS:
            M = co.masks[s]
            out.partition(base + ("levels", s), {str(v): int((M & (lvl == v)).sum()) for v in levels}, total=dens.enc[s])


def write_dob_breakdown(out: Out, co: Cohort, dens: Dens, la: pd.DataFrame,
                        weight_source: np.ndarray | None) -> None:
    """Aggregate DOB/age report (D20). All counts go through the registry (primary <10, complementary S)."""
    status = la["dob_status"].to_numpy()
    nbin = la["dob_neg_bin"].to_numpy()
    fcu = co.enc["first_careunit"].astype(str).to_numpy()
    units = sorted(pd.unique(fcu).tolist(), key=str)
    b0 = ("layer_a", "dob_breakdown")
    out.val(b0 + ("safety_gate",),
            f"stop if invalid > {100 * MAX_INVALID_AGE_FRAC:.0f}% of the cohort "
            f"(1% of the cohort; the public copy omits the pinned cohort size)")
    out.val(b0 + ("weight_imputation_rule",), WEIGHT_IMPUTE_RULE)
    st_counts = {k: int((status == k).sum()) for k in DOB_STATUS}
    sc = out.partition(b0 + ("status",), st_counts, total=dens.n_all)
    n_inv = out.cell(b0 + ("n_invalid",), sum(st_counts[k] for k in ("missing", "unparseable", "negative")))
    out.ident(n_inv, [sc["missing"], sc["unparseable"], sc["negative"]])
    # age-continuous quantiles: valid ages only (NaN already dropped by Out.q)
    out.q(b0, "age_days", la["age_days"].to_numpy(), equal_to=sc["valid"])
    nbc = {k: int((nbin == k).sum()) for k in DOB_NEG_BINS}
    out.partition(b0 + ("negative_bins",), nbc, total=sc["negative"])
    inv_split = {}
    for s in SPLITS:
        M = co.masks[s]
        bs = b0 + ("by_split", s)
        stc = {k: int((M & (status == k)).sum()) for k in DOB_STATUS}
        scs = out.partition(bs + ("status",), stc, total=dens.enc[s])
        inv_split[s] = out.cell(bs + ("n_invalid",), sum(stc[k] for k in ("missing", "unparseable", "negative")))
        out.ident(inv_split[s], [scs["missing"], scs["unparseable"], scs["negative"]])
        out.q(bs, "age_days", la.loc[M, "age_days"].to_numpy(), equal_to=scs["valid"])
        out.partition(bs + ("negative_bins",), {k: int((M & (nbin == k)).sum()) for k in DOB_NEG_BINS},
                      total=scs["negative"])
    # invalid by first_careunit (overall) and by split x first_careunit
    unit_inv = {}
    for u in units:
        Mu = fcu == u
        bu = b0 + ("by_first_careunit", u)
        stc = {k: int((Mu & (status == k)).sum()) for k in ("missing", "unparseable", "negative")}
        scu = out.partition(bu + ("status",), stc)
        unit_inv[u] = out.cell(bu + ("n_invalid",), sum(stc.values()))
        out.ident(unit_inv[u], list(scu.values()))
    if unit_inv:
        out.ident(n_inv, list(unit_inv.values()))
    for s in SPLITS:
        M = co.masks[s]
        parts = []
        for u in units:
            bu = b0 + ("by_split_and_first_careunit", s, u)
            stc = {k: int((M & (fcu == u) & (status == k)).sum()) for k in ("missing", "unparseable", "negative")}
            scu = out.partition(bu + ("status",), stc)
            cid = out.cell(bu + ("n_invalid",), sum(stc.values()))
            out.ident(cid, list(scu.values()))
            parts.append(cid)
        if parts:
            out.ident(inv_split[s], parts)
    # weight fallback counts (cohort; plus age_invalid subset)
    if weight_source is not None:
        ws = np.asarray(weight_source, dtype=object)
        inv = status != "valid"
        bw = b0 + ("weight_imputation",)
        out.partition(bw, {k: int((ws == k).sum()) for k in WEIGHT_SOURCES})
        out.partition(bw + ("age_invalid",), {k: int((inv & (ws == k)).sum()) for k in WEIGHT_SOURCES},
                      total=n_inv)


def write_layer_a(out: Out, co: Cohort, dens: Dens, la: pd.DataFrame, pre_tables: dict, surg: pd.DataFrame,
                  neuro_match, weight_any: np.ndarray, gaps_harm_disagree: int | None, prior_disagree: int | None,
                  weight_source: np.ndarray | None = None) -> dict:
    enc = co.enc
    write_dob_breakdown(out, co, dens, la, weight_source)
    subj = enc["SUBJECT_ID"].to_numpy()
    t0 = enc["t0"]
    hrs = la["hours_admit_to_icu"].to_numpy(dtype=float)
    hour = t0.dt.hour.to_numpy()
    wkend = t0.dt.dayofweek.to_numpy() >= 5
    shift = shift_of(hour)
    day = np.where(shift == "day", "daytime", "evening_night")
    phi2 = np.where(wkend, "weekend_", "weekday_") + day
    band = interval_band(hrs)
    dxn = la["dx_norm"].to_numpy()
    dxs = la["dx_status"].to_numpy()
    c3 = np.array([x[:3] if st in ("ok", "multiple") else "" for x, st in zip(dxn, dxs)], dtype=object)
    chapter = np.array([x[0] if x else "missing" for x in c3], dtype=object)
    domain = np.array([icd_domain(x) for x in c3], dtype=object)
    cat = np.where(c3 == "", "missing_or_invalid", c3)
    neuro = la["neuro"].to_numpy(dtype=bool)
    # admission source (first match) at the X grid
    bmin = np.full(co.n, np.inf)
    if len(surg):
        lk = surg[(surg["link_class"] == "linked") & surg["enc_i"].notna()]
        np.minimum.at(bmin, lk["enc_i"].astype(int).to_numpy(), lk["B_h"].to_numpy(dtype=float))
    or_before = bmin < 0
    at_intime = np.zeros(co.n, dtype=bool)
    if len(surg):
        lk = surg[(surg["link_class"] == "linked") & surg["enc_i"].notna()]
        at_intime[lk.loc[lk["B_h"] == 0, "enc_i"].astype(int).to_numpy()] = True
    ed = la["ED"]
    ed_before = (ed < t0).fillna(False).to_numpy(dtype=bool)
    ed_class = np.select([la["ed_raw"].fillna("").to_numpy() == "", ed.isna().to_numpy(),
                          (ed == la["ADMITTIME"]).fillna(False).to_numpy(dtype=bool), ed_before],
                         ["missing", "unparseable", "equals_admittime", "before_intime_other"], "at_or_after_intime")
    prior = np.zeros(co.n, dtype=int)
    o = enc.sort_values(["SUBJECT_ID", "t0", "encounter_id"])
    prior[o["enc_i"].to_numpy()] = o.groupby("SUBJECT_ID").cumcount().to_numpy()
    for s in SPLITS:
        M = co.masks[s]
        b = ("layer_a", s)
        part = lambda name, arr, levels: out.partition(b + (name,), {str(v): int((M & (arr == v)).sum()) for v in levels},
                                                       total=dens.enc[s])
        part("sex", la["sex"].to_numpy(), ("M", "F", "unrecognized"))
        part("age_band", la["age_band"].to_numpy(), AGE_BANDS)
        part("hours_admit_to_icu_band", band, ("missing", "lt24h", "24to72h", "72to168h", "ge168h"))
        out.cell(b + ("n_hours_admit_to_icu_negative",), int((M & (hrs < 0)).sum()), parent=dens.enc[s])
        part("shift", shift, ("day", "evening", "night"))
        part("day_type", np.where(wkend, "weekend", "weekday"), ("weekday", "weekend"))
        part("timing_phi2", phi2, ("weekday_daytime", "weekday_evening_night", "weekend_daytime", "weekend_evening_night"))
        part("first_careunit", enc["first_careunit"].to_numpy(), sorted(pd.unique(enc["first_careunit"]).tolist()))
        part("diagnosis_code_status", dxs, ("ok", "multiple", "invalid", "missing"))
        part("diagnosis_domain_draft", domain, ("cardiac_congenital_heart", "neonatal_perinatal", "surgical_trauma",
                                                "respiratory_infection", "neurological", "oncology_hematology", "gi_renal",
                                                "other_medical"))
        out.cell(b + ("n_neuro_admitting_dx",), int((M & neuro).sum()), parent=dens.enc[s])
        for X in ADM_SOURCE_X:
            floor = hrs >= X
            src = np.select([floor, or_before, ed_before], ["floor", "or", "ed"], "other")
            part(f"admission_source_x{X}h", src, ("floor", "or", "ed", "other"))
            out.cell(b + (f"admission_source_x{X}h", "n_conflict_ed_registration_and_interval_ge_x"),
                     int((M & floor & ed_before).sum()), parent=dens.enc[s])
        out.cell(b + ("n_surgery_begin_exactly_at_intime",), int((M & at_intime).sum()), parent=dens.enc[s])
        part("edregtime_class", ed_class, ("missing", "unparseable", "equals_admittime", "before_intime_other", "at_or_after_intime"))
        part("prior_encounters", np.where(prior >= 3, "3plus", prior.astype(str)), ("0", "1", "2", "3plus"))
        for tname, pa in pre_tables.items():
            for w in PRE_ICU_H:
                out.cell(b + ("pre_icu_any_row", tname, f"n_{int(w)}h"), int((M & pa[w]).sum()), parent=dens.enc[s])
        out.cell(b + ("n_charted_weight",), int((M & weight_any).sum()), parent=dens.enc[s])
    level_support(out, ("layer_a", "level_support", "diagnosis_category"), co, cat, subj, dens)
    level_support(out, ("layer_a", "level_support", "diagnosis_chapter"), co, chapter, subj, dens)
    level_support(out, ("layer_a", "level_support", "diagnosis_domain_draft"), co, domain, subj, dens, per_level=False)
    level_support(out, ("layer_a", "level_support", "first_careunit"), co, enc["first_careunit"].to_numpy(), subj, dens, per_level=False)
    for X in ADM_SOURCE_X:
        src = np.select([hrs >= X, or_before, ed_before], ["floor", "or", "ed"], "other")
        level_support(out, ("layer_a", "level_support", f"admission_source_x{X}h"), co, src, subj, dens, per_level=False)
    if len(surg):
        rec = surg[(surg["link_class"] == "linked") & (surg["timing"] == "recent_layer_a")]
        sysarr = np.full(co.n, "none", dtype=object)
        rr = rec.sort_values(["enc_i", "T"]).drop_duplicates("enc_i")
        sysarr[rr["enc_i"].astype(int).to_numpy()] = rr["body_system"].to_numpy()
        level_support(out, ("layer_a", "level_support", "surgery_system_first_recent"), co, sysarr, subj, dens, per_level=False)
    if gaps_harm_disagree is not None:
        out.cell(("cohort_extra", "n_gap_presence_disagrees_with_harmonization"), gaps_harm_disagree)
    if prior_disagree is not None:
        out.cell(("cohort_extra", "n_prior_encounters_disagree_with_harmonization"), prior_disagree)
    return {"neuro": neuro}


# =============================================================================
# Map checks (spec 5.3)
# =============================================================================
def _train_window(co: Cohort, e: np.ndarray, t: np.ndarray, split: str = "train") -> np.ndarray:
    e = np.asarray(e, dtype=np.int64)
    return co.masks[split][e] & (np.asarray(t) >= 0) & (np.asarray(t) < co.W[e])


def _safe_label(v: str) -> str | None:
    v = str(v)
    if not v or len(v) > 20 or TIMESTAMP_RE.search(v) or re.search(r"\d{6,}", v):
        return None
    return v


def impute_weights(co: Cohort, la: pd.DataFrame, wt_first: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Charted weight 1014 (first value in [t0 - 24 h, W)), else TRAIN medians (A2.12, D20).

    Valid age bands: age band x sex (n>=10), else age band (n>=10), else sex (n>=10), else overall.
    age_invalid is never imputed from an age band; it uses sex-only, else overall.
    Medians are from train encounters with a charted weight; age-band medians exclude age_invalid.
    Returns (weight, source) with source in WEIGHT_SOURCES (empty if no value could be filled).
    """
    ab_all = la["age_band"].to_numpy()
    sx_all = la["sex"].to_numpy()
    tr = co.masks["train"] & ~np.isnan(wt_first)
    df = pd.DataFrame({"w": wt_first[tr], "ab": ab_all[tr], "sx": sx_all[tr]})
    df_age = df[df["ab"] != "age_invalid"]
    g2 = df_age.groupby(["ab", "sx"])["w"].agg(["median", "size"])
    g1 = df_age.groupby("ab")["w"].agg(["median", "size"])
    g_sex = df.groupby("sx")["w"].agg(["median", "size"])
    overall = float(df["w"].median()) if len(df) else np.nan
    w = wt_first.copy()
    source = np.where(np.isfinite(w), "charted", "").astype(object)
    for i in np.flatnonzero(~np.isfinite(w)):
        ab, sx = ab_all[i], sx_all[i]
        if ab != "age_invalid" and (ab, sx) in g2.index and g2.loc[(ab, sx), "size"] >= MIN_CELL:
            w[i] = g2.loc[(ab, sx), "median"]
            source[i] = "imputed_age_sex"
        elif ab != "age_invalid" and ab in g1.index and g1.loc[ab, "size"] >= MIN_CELL:
            w[i] = g1.loc[ab, "median"]
            source[i] = "imputed_age_band"
        elif sx in g_sex.index and g_sex.loc[sx, "size"] >= MIN_CELL:
            w[i] = g_sex.loc[sx, "median"]
            source[i] = "imputed_sex_only"
        else:
            w[i] = overall
            source[i] = "imputed_overall" if np.isfinite(overall) else ""
    return w, source


def write_map_checks(out: Out, co: Cohort, dens: Dens, vb: ValueBuild, cm: dict, la: pd.DataFrame) -> dict:
    A = out.A
    res: dict[str, Any] = {}
    rs = vb.rowstat
    rsel = _train_window(co, rs["enc_i"].to_numpy(), rs["t"].to_numpy()) if len(rs) else np.zeros(0, dtype=bool)
    r = rs[rsel]
    # row status per channel (gates: out-of-gate values are counted, never measurements)
    for ch in sorted(cm["channels"]):
        if cm["channels"][ch]["derived"]:
            continue
        g = r[r["channel"] == ch]
        out.partition(("map_checks", "row_status", ch), {st: int((g["status"] == st).sum()) for st in ROW_STATUS})
    # per-unit row counts (items with more than one unit string) - D13
    for (tb, it), info in sorted(cm["items"].items()):
        g = r[(r["table"] == tb) & (r["item"] == it)]
        units = sorted(pd.unique(g["unit"]).tolist())
        if len(units) > 1:
            out.partition(("map_checks", "units", f"{tb}_{it}"), {(u if u else BLANK): int((g["unit"] == u).sum()) for u in units})
    # duplicate merge groups (D6): quantiles, KS vs the main item, same-draw differences, text values
    meas = r[r["status"] == "measurement"]
    for ch, main in sorted(cm["merge_main"].items()):
        items = sorted(cm["channels"][ch]["items"])
        mv = meas[(meas["channel"] == ch)]
        vals = {it: mv.loc[mv["item"] == it, "v"].to_numpy(dtype=float) for it in items}
        ncell = {}
        for it in items:
            b = ("map_checks", "duplicates", ch, "items", it)
            ncell[it] = out.q(b, "values", vals[it])
            if ch == "calcium_total":
                vv = vals[it][np.isfinite(vals[it])]
                out.cell(b + ("n_in_ionized_range",), int(((vv >= IONIZED_RANGE[0]) & (vv <= IONIZED_RANGE[1])).sum()), parent=ncell[it])
            tx = mv.loc[mv["item"] == it, "vtext"].astype(str)
            if not np.isfinite(vals[it]).any() or ch in ("lipase", "stool_occult_blood"):
                vc = tx.value_counts()
                listed = {}
                for v, c in vc.items():
                    lab = _safe_label(v)
                    if lab is not None and c >= MIN_CELL and len(listed) < 8:
                        listed[lab] = int(c)
                listed["(other)"] = int(len(tx) - sum(listed.values()))
                out.partition(b + ("text_values",), listed)
        for it in items:
            if it == main:
                continue
            a, bb = vals[main][np.isfinite(vals[main])], vals[it][np.isfinite(vals[it])]
            ksv = float(sps.ks_2samp(a, bb).statistic) if len(a) >= MIN_CELL and len(bb) >= MIN_CELL else None
            out.ks(("map_checks", "duplicates", ch, "items", it), "vs_main_item", [ncell[main], ncell[it]], ksv)
            pm = mv[mv["item"] == main].drop_duplicates(["enc_i", "tc"])[["enc_i", "tc", "v"]]
            po = mv[mv["item"] == it].drop_duplicates(["enc_i", "tc"])[["enc_i", "tc", "v"]]
            pr = po.merge(pm, on=["enc_i", "tc"], suffixes=("_o", "_m"))
            out.q(("map_checks", "duplicates", ch, "items", it), "same_draw_difference_vs_main", (pr["v_o"] - pr["v_m"]).to_numpy())
        out.val(("map_checks", "duplicates", ch, "main_item"), main)
    # magnesium 'ug/dl' rows (D13): held out of the measurements; compared as recorded and after scaling
    mg = r[(r["table"] == "LABEVENTS") & (r["item"] == "5059")]
    mm = mg.loc[mg["status"] == "measurement", "v"].to_numpy(dtype=float)
    ug = mg.loc[(mg["unit"] == "ug/dl") & mg["vnum"].notna(), "vnum"].to_numpy(dtype=float)
    b = ("map_checks", "magnesium_ug_dl_rows")
    n_mm = out.q(b, "mmol_l_values", mm)
    for lab, f in MG_SCALINGS.items():
        n_u = out.q(b + (lab,), "values_mmol_l", ug * f)
        out.ks(b + (lab,), "vs_mmol_l_rows", [n_mm, n_u],
               float(sps.ks_2samp(mm, ug * f).statistic) if len(mm) >= MIN_CELL and len(ug) >= MIN_CELL else None)
    # HR (D9)
    hm = meas[meas["channel"] == "heart_rate"]
    b = ("map_checks", "heart_rate")
    n2 = out.q(b + ("1002",), "values", hm.loc[hm["item"] == HR_DROP, "v"].to_numpy(dtype=float))
    n3 = out.q(b + ("1003",), "values", hm.loc[hm["item"] == HR_KEEP, "v"].to_numpy(dtype=float))
    a2, a3 = hm.loc[hm["item"] == HR_DROP, "v"].to_numpy(dtype=float), hm.loc[hm["item"] == HR_KEEP, "v"].to_numpy(dtype=float)
    out.ks(b, "1002_vs_1003", [n2, n3], float(sps.ks_2samp(a2, a3).statistic) if len(a2) >= MIN_CELL and len(a3) >= MIN_CELL else None)
    hp = vb.hr["pairs"]
    hp = hp[_train_window(co, hp["enc_i"].to_numpy(), hp["tc"].to_numpy())]
    nd = out.q(b, "same_time_difference_1002_minus_1003", (hp["v"] - hp["v_1003"]).to_numpy(), parent=n2)
    out.val(b + ("rule",), "1003 kept when both exist at the same CHARTTIME; 1002 kept otherwise")
    # MAP (D8)
    b = ("map_checks", "map")
    mp = vb.map_pairs
    sbr = mp["sbp_rows"][_train_window(co, mp["sbp_rows"]["enc_i"].to_numpy(), mp["sbp_rows"]["t"].to_numpy())]
    dbr = mp["dbp_rows"][_train_window(co, mp["dbp_rows"]["enc_i"].to_numpy(), mp["dbp_rows"]["t"].to_numpy())]
    pr = mp["pairs"][_train_window(co, mp["pairs"]["enc_i"].to_numpy(), mp["pairs"]["t"].to_numpy())]
    ns = out.cell(b + ("n_sbp_rows_numeric",), len(sbr))
    out.cell(b + ("n_sbp_rows_with_same_time_dbp",), len(pr), parent=ns)
    ndb = out.cell(b + ("n_dbp_rows_numeric",), len(dbr))
    out.cell(b + ("n_dbp_rows_with_same_time_sbp",), len(pr), parent=ndb)
    out.partition(b + ("pairs",), {k: int((pr["reason"] == k).sum()) for k in ("valid", "dbp_ge_sbp", "outside_gate")})
    # P/F (D7, D19) and S/F (D15)
    pf = vb.pf["pairs"]
    pft = pf[_train_window(co, pf["enc_i"].to_numpy(), pf["t"].to_numpy())]
    b = ("map_checks", "pf_ratio")
    out.partition(b + ("same_draw_pairs",), {"valid": int(pft["valid"].sum()), "outside_gate": int((~pft["valid"]).sum())})
    wl = pft[pft["valid"] & pft["v_lab"].notna()]
    out.q(b, "difference_vs_analyser_5243", (wl["v"] - wl["v_lab"]).to_numpy())
    sc = psofa_score(pft.loc[pft["valid"], "v"].to_numpy(), PSOFA_PF)
    out.partition(b + ("psofa_points_d19",), {str(k): int((sc == k).sum()) for k in range(5)})
    sfr = vb.sf["rows"]
    b = ("map_checks", "sf_ratio")
    fd = vb.pf["fio2_draws"]
    fdt = fd[_train_window(co, fd["enc_i"].to_numpy(), fd["t"].to_numpy())]
    ntot = out.cell(b + ("n_fio2_draws",), len(fdt))
    sft = sfr[_train_window(co, sfr["enc_i"].to_numpy(), sfr["tc"].to_numpy() + LAB_RELEASE_H)] if len(sfr) else sfr
    n_with_pf = len(fdt) - len(sft)
    out.partition(b + ("fio2_draws",), {"with_valid_pf": n_with_pf,
                                        "no_spo2_within_5min": int((sft["reason"] == "no_spo2_within_5min").sum()) if len(sft) else 0,
                                        "spo2_gt_97": int((sft["reason"] == "spo2_gt_97").sum()) if len(sft) else 0,
                                        "sf_computed": int((sft["reason"] == "sf").sum()) if len(sft) else 0}, total=ntot)
    sfv = sft[sft["reason"] == "sf"]["v"].to_numpy(dtype=float) if len(sft) else np.array([])
    sc = psofa_score(sfv, PSOFA_SF)
    out.partition(b + ("psofa_points_d15",), {str(k): int((sc == k).sum()) for k in range(5)})
    out.val(("map_checks", "oxygenation_limitations"), "no ventilation item in PIC (support requirement not applied); "
            "blood gases not marked arterial (P/F may include venous or capillary samples)")
    # coverage per landmark: P/F, and S/F where no P/F exists (S/F coverage count)
    Fpf = event_features(co, vb.events["pf_ratio"]["enc_i"].to_numpy(), vb.events["pf_ratio"]["t"].to_numpy())
    Fsf = event_features(co, vb.events["sf_ratio"]["enc_i"].to_numpy(), vb.events["sf_ratio"]["t"].to_numpy())
    for s in SPLITS:
        for L in LANDMARKS:
            E = co.masks[s] & co.elig[L]
            ap, as_ = Fpf["avail"][1.0][("full", L)], Fsf["avail"][1.0][("full", L)]
            bl = ("map_checks", "oxygenation_coverage", s, f"{L}h")
            out.cell(bl + ("n_pf_available",), int((E & ap).sum()), parent=dens.elig[s][L])
            out.cell(bl + ("n_sf_available_without_pf",), int((E & as_ & ~ap).sum()), parent=dens.elig[s][L])
            out.cell(bl + ("n_pf_or_sf_available",), int((E & (ap | as_)).sum()), parent=dens.elig[s][L])
    # urine output 1008 and fluid intake 1009 (D10, D17)
    wev = vb.events["weight"]
    wsel = (wev["t"].to_numpy() >= -24) & (wev["t"].to_numpy() < co.W[wev["enc_i"].to_numpy()])
    wf = np.full(co.n, np.nan)
    ww = wev[wsel].sort_values(["enc_i", "t"]).drop_duplicates("enc_i")
    wf[ww["enc_i"].to_numpy()] = ww["v"].to_numpy(dtype=float)
    weight, wsrc = impute_weights(co, la, wf)
    res["weight_any"] = ~np.isnan(wf)
    res["weight_source"] = wsrc
    t0h = (co.enc["t0"].dt.hour + co.enc["t0"].dt.minute / 60.0).to_numpy()
    for ch in ("urine_output", "fluid_intake"):
        ev = vb.events[ch]
        sel = _train_window(co, ev["enc_i"].to_numpy(), ev["t"].to_numpy())
        ev = ev[sel].sort_values(["enc_i", "t"], kind="mergesort")
        e, t, v = ev["enc_i"].to_numpy(), ev["t"].to_numpy(), ev["v"].to_numpy(dtype=float)
        b = ("map_checks", "charting_interval", ch)
        same = e[1:] == e[:-1]
        gap = (t[1:] - t[:-1])[same]
        out.q(b, "gap_between_entries_hours", gap)
        hod = np.floor((t0h[e] + t) % 24).astype(int)
        out.partition(b + ("hour_of_day",), {f"{h:02d}": int((hod == h).sum()) for h in range(24)})
        out.partition(b + ("shift",), {k: int((shift_of(hod) == k).sum()) for k in ("day", "evening", "night")})
        day = np.floor(t / 24.0).astype(int)
        complete = (day + 1) * 24.0 <= co.W[e]
        dd = pd.DataFrame({"e": e[complete], "d": day[complete], "v": v[complete]}).groupby(["e", "d"])["v"].agg(["size", "sum"]).reset_index()
        out.partition(b + ("entries_per_complete_day",), {"1": int((dd["size"] == 1).sum()), "2to3": int(dd["size"].between(2, 3).sum()),
                                                          "4to11": int(dd["size"].between(4, 11).sum()), "12plus": int((dd["size"] >= 12).sum())})
        full_days = int(np.floor(co.W[co.masks["train"]] / 24.0).sum())
        tdays = out.cell(b + ("n_complete_days_train",), full_days)
        out.cell(b + ("n_complete_days_without_entry",), full_days - len(dd), parent=tdays)
        wk = weight[dd["e"].to_numpy()]
        out.q(b, "daily_total_ml_per_kg_per_day", dd["sum"].to_numpy() / wk)
        if ch == "urine_output":
            ge = e[1:][same]
            ok = gap > 0
            out.q(b, "ml_per_kg_per_hour_over_charting_interval", (v[1:][same][ok] / weight[ge[ok]] / gap[ok]))
        has = np.zeros(co.n, dtype=bool)
        has[e] = True
        out.partition(b + ("weight_source_encounters",),
                      {k: int((has & (wsrc == k)).sum()) for k in WEIGHT_SOURCES})
    out.val(("map_checks", "charting_interval", "rule"), "an empty bin or a day with no entry is never read as zero output or intake")
    return res


# =============================================================================
# Document, markdown, release gate, publish
# =============================================================================
def detect_count_fields(obj: Any, prefix: str = "") -> list[str]:
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if prefix == "" and k in NON_COUNT_TOP:
                continue
            ptr = prefix + "/" + str(k).replace("~", "~0").replace("/", "~1")
            if COUNT_KEY_RE.match(str(k)) and not isinstance(v, (dict, list)):
                out.append(ptr)
            else:
                out.extend(detect_count_fields(v, ptr))
    return out


def detect_prefixed(obj: Any, prefixes: tuple, prefix: str = "") -> list[str]:
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if prefix == "" and k in NON_COUNT_TOP:
                continue
            ptr = prefix + "/" + str(k).replace("~", "~0").replace("/", "~1")
            if str(k).startswith(prefixes):
                out.append(ptr)
            elif isinstance(v, dict):
                out.extend(detect_prefixed(v, prefixes, ptr))
    return out


def _f(v: Any) -> str:
    if isinstance(v, dict):
        return "/".join(_f(v.get(k)) for k in ("p25", "p50", "p75"))
    return "–" if v is None else str(v).replace("|", "/")


def _g(doc: dict, *path) -> Any:
    cur = doc
    for p in path:
        if not isinstance(cur, dict) or p not in cur:
            return None
        cur = cur[p]
    return cur


def render_markdown(doc: dict) -> str:
    L = ["# Feature sparsity and measurement-frequency audit (stage 2)", "",
         f"- spec: {doc['spec']}", f"- release gate: {doc['release_gate']}",
         f"- outcome descriptors (--outcome-descriptors): {doc['outcome_descriptors_flag']}",
         f"- suppression: primary {doc['suppression'].get('primary_cells')}, complementary {doc['suppression'].get('complementary_cells')}, "
         f"recoverable {doc['suppression'].get('recoverable')}", "",
         "Aggregate only. Counts under 10 are shown as <10; S marks complementary suppression. Train split shown; "
         "validation and test are in the JSON. Provisional classes apply the approved spec 9 defaults on train and "
         "are for Vinai's review, not decisions.", "", "## Input hashes", "", "| input | sha256 |", "|---|---|"]
    for k, v in sorted(doc["input_hashes"].items()):
        L.append(f"| {k} | {v} |")
    L += ["", "## Raw tables", "", "| table | gate | sha256 |", "|---|---|---|"]
    for k, v in sorted(doc["raw_tables"].items()):
        L.append(f"| {k} | {v['gate']} | {v['sha256']} |")
    for sec, title in (("channels", "Channels"), ("drug_groups", "Drug groups"), ("specimen_groups", "Specimen groups (M1)")):
        L += ["", f"## {title} (train)", "",
              "| name | n_ever_layer_b | n_ever_336h | pct_ever_336h | n_avail_full_24h | n_avail_full_72h | n_avail_full_168h | rate_per_day_all p25/p50/p75 | class |",
              "|---|---|---|---|---|---|---|---|---|"]
        for name in sorted(doc.get(sec, {})):
            t = doc[sec][name].get("train", {})
            cls = _g(doc, "provisional_classes", sec, name, "class")
            L.append("| " + " | ".join([name, _f(t.get("n_ever_layer_b")), _f(t.get("n_ever_336h")), _f(t.get("pct_ever_336h")),
                                         _f(_g(t, "landmarks", "24h", "n_avail_full")), _f(_g(t, "landmarks", "72h", "n_avail_full")),
                                         _f(_g(t, "landmarks", "168h", "n_avail_full")), _f(t.get("q_rate_per_day_all")),
                                         _f(cls)]) + " |")
    L += ["", "## Provisional class counts (channels)", ""]
    cc: dict[str, int] = {}
    for v in doc.get("provisional_classes", {}).get("channels", {}).values():
        cc[v["class"]] = cc.get(v["class"], 0) + 1
    for k in sorted(cc):
        L.append(f"- {k}: {cc[k]}")
    db = _g(doc, "layer_a", "dob_breakdown")
    if isinstance(db, dict):
        L += ["", "## Layer A DOB breakdown (D20)", "",
              "Invalid ages (missing, unparseable, or negative) are the age_invalid band. "
              "Counts under 10 are shown as <10; S marks complementary suppression. "
              "Age quantiles use valid ages only.", ""]
        st = db.get("status", {})
        L.append("- status: " + ", ".join(f"{k}={_f(st.get(f'n_{k}'))}" for k in DOB_STATUS) +
                 f"; invalid={_f(db.get('n_invalid'))}")
        nb = db.get("negative_bins", {})
        L.append("- negative bins: " + ", ".join(f"{k}={_f(nb.get(f'n_{k}'))}" for k in DOB_NEG_BINS))
        L.append(f"- age_days quantiles (valid only): {_f(db.get('q_age_days'))} (n={_f(db.get('n_age_days'))})")
        L.append(f"- safety gate: {_f(db.get('safety_gate'))}")
        L.append(f"- weight imputation: {_f(db.get('weight_imputation_rule'))}")
        wi = db.get("weight_imputation", {})
        if isinstance(wi, dict):
            L.append("- weight source (cohort): " +
                     ", ".join(f"{k}={_f(wi.get(f'n_{k}'))}" for k in WEIGHT_SOURCES))
            wia = wi.get("age_invalid", {})
            if isinstance(wia, dict):
                L.append("- weight source (age_invalid): " +
                         ", ".join(f"{k}={_f(wia.get(f'n_{k}'))}" for k in WEIGHT_SOURCES))
        L.append("- by split (invalid): " +
                 ", ".join(f"{s}={_f(_g(db, 'by_split', s, 'n_invalid'))}" for s in SPLITS))
        fcu = db.get("by_first_careunit") or {}
        if isinstance(fcu, dict) and fcu:
            units = sorted(k for k in fcu if isinstance(fcu[k], dict))
            L.append("- by first_careunit (invalid): " +
                     ", ".join(f"{u}={_f(_g(fcu, u, 'n_invalid'))}" for u in units))
    L += ["", "## Notes", ""]
    for k, v in sorted(doc.get("notes", {}).items()):
        L.append(f"- {k}: {v}")
    return "\n".join(L) + "\n"


def build_document(out: Out, status: dict, supp: dict, hashes: dict, raw_tables: dict, map_hashes: dict, checks: dict,
                   cfg: S2Config, provisional: dict, notes: dict) -> dict:
    doc: dict[str, Any] = {
        "spec": SPEC_VERSION, "release_gate": "pending", "outcome_descriptors_flag": bool(cfg.outcome_descriptors),
        "input_hashes": dict(hashes), "v17_module": {"file": V17_MODULE_FILE, "sha256": cfg.v17_module_sha256},
        "maps": {k: {"file": MAP_FILES[k], "sha256": v} for k, v in map_hashes.items()},
        "raw_tables": raw_tables, "cohort_checks": checks,
        "config": {"lab_release_hours": LAB_RELEASE_H, "lab_release_sensitivities_hours": [0.0, 2.0],
                   "layer_b_window": "[t0, min(OUT, t0 + 24 h))", "analysis_window": "[t0, min(OUT, t0 + 336 h))",
                   "pre_icu_windows": ["[t0 - 6 h, t0)", "[t0 - 24 h, t0)"], "bin_hours": BIN_H, "bins": N_BINS,
                   "landmarks_hours": list(LANDMARKS), "detail_min_train_encounters": DETAIL_MIN_TRAIN,
                   "strata_landmark_hours": STRATA_LANDMARK, "min_cell": MIN_CELL,
                   "cutoffs_hours": {g: group_cutoff(g) for g in MED_GROUPS}, "cutoff_factors": list(CUTOFF_FACTORS),
                   "micro_turnaround_hours": MICRO_TA_H, "micro_grid": list(MICRO_GRID),
                   "admission_source_x_grid_hours": list(ADM_SOURCE_X), "sf_pairing_minutes": SF_PAIR_MINUTES,
                   "sf_spo2_max": SF_SPO2_MAX, "psofa_pf_d19": ">=400 0; 300-399 1; 200-299 2; 100-199 3; <100 4",
                   "psofa_sf_d15": ">=292 0; 264-291 1; 221-263 2; 148-220 3; <148 4",
                   "dob_invalid_safety_gate_frac": MAX_INVALID_AGE_FRAC,
                   "dob_invalid_safety_gate_encounters_on_cohort": None,  # public copy: pin omitted
                   "weight_imputation_rule": WEIGHT_IMPUTE_RULE},
        "provisional_classes": provisional, "notes": notes,
    }
    fill_document(out, status, doc)
    published = out.reg.published()
    doc["disclosure"] = {"count_fields": sorted(published),
                         "identities": [[t, list(p)] for t, p in out.reg.identities],
                         "hidden_cells": sorted(out.reg.hidden),
                         "ratios": [[p, n, d] for p, n, d, _ in out.ratios],
                         "dependents": [[p, list(nd)] for p, nd, _ in out.kss],
                         "quantiles": [[p, n] for p, n, _ in out.quants]}
    doc["suppression"] = dict(supp)
    return doc


def release_gate(A, json_path: Path, md_path: Path, expect_status: str, sentinels: Iterable[str] = ()) -> list[str]:
    fails: list[str] = []
    try:
        doc = json.loads(Path(json_path).read_text(encoding="utf-8"))
        md = Path(md_path).read_text(encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        return [f"cannot reload shared files: {e}"]
    if doc.get("spec") != SPEC_VERSION or doc.get("release_gate") != expect_status:
        fails.append("header fields (spec / release_gate) wrong")
    if not isinstance(doc.get("outcome_descriptors_flag"), bool):
        fails.append("outcome_descriptors_flag not recorded")
    elif not doc["outcome_descriptors_flag"] and ("outcome_descriptors" in doc or "outcome_denominators" in doc):
        fails.append("outcome descriptors present although the flag is off")
    for p in (json_path, md_path):
        if "_raw_LOCAL" in Path(p).name:
            fails.append("a shared file name contains _raw_LOCAL")
    dis = doc.get("disclosure", {})
    declared = set(dis.get("count_fields", []))
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
    for h in dis.get("hidden_cells", []):
        if not str(h).startswith("hidden:"):
            fails.append("hidden cell without the hidden: prefix")
        unknown[h] = "H"
    sup = doc.get("suppression", {})
    if sup.get("primary_cells") != n_p or sup.get("complementary_cells") != n_s or sup.get("recoverable") != 0:
        fails.append("suppression counters do not match the published cells")
    ids = [(t, tuple(p)) for t, p in dis.get("identities", [])]
    if not ids:
        fails.append("no disclosure identities published")
    # quantiles: only next to their own published count >= 10
    qdecl = {p: n for p, n in dis.get("quantiles", [])}
    for ptr in detect_prefixed(doc, ("q_",)):
        v = A.get_ptr(doc, ptr)
        if ptr not in qdecl:
            fails.append(f"undeclared quantile at {ptr}")
            continue
        sib = ptr.rsplit("/", 1)[0] + "/n_" + ptr.rsplit("/", 1)[1][2:]
        if qdecl[ptr] != sib:
            fails.append(f"quantile not paired with its own count at {ptr}")
        if v == "suppressed":
            continue
        if not isinstance(v, dict) or set(v) != set(QKEYS):
            fails.append(f"quantile with disallowed keys at {ptr}")
        n = known.get(sib)
        if n is None or n < MIN_CELL:
            fails.append(f"quantile published without a published count >= 10 at {ptr}")
    rdecl = {p: (n, d) for p, n, d in dis.get("ratios", [])}
    for ptr in detect_prefixed(doc, ("pct_",)):
        v = A.get_ptr(doc, ptr)
        if ptr not in rdecl:
            fails.append(f"undeclared percentage at {ptr}")
            continue
        n, d = rdecl[ptr]
        if v == "suppressed":
            continue
        if n not in known or d not in known or known[d] == 0 or v != r1(100.0 * known[n] / known[d]):
            fails.append(f"percentage not reproducible from published counts at {ptr}")
    ddecl = {p: nd for p, nd in dis.get("dependents", [])}
    for ptr in detect_prefixed(doc, ("ks_", "x_")):
        v = A.get_ptr(doc, ptr)
        if ptr not in ddecl:
            fails.append(f"undeclared dependent value at {ptr}")
            continue
        if v != "suppressed" and not all(n in known for n in ddecl[ptr]):
            fails.append(f"dependent value published although a count it rests on is suppressed at {ptr}")
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
                if str(k).lower() in ("min", "max", "minimum", "maximum"):
                    fails.append(f"min/max published at {path}/{k}")
                if TIMESTAMP_RE.search(str(k)):
                    fails.append(f"timestamp-like key at {path}")
                walk(v, f"{path}/{k}")
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}/{i}")
        elif isinstance(o, str) and TIMESTAMP_RE.search(o):
            fails.append(f"timestamp-like string at {path}")
    walk({k: v for k, v in doc.items() if k != "disclosure"})
    try:
        if render_markdown(doc) != md:
            fails.append("markdown on disk is not the rendering of the shared JSON")
    except Exception as e:  # noqa: BLE001
        fails.append(f"shared markdown cannot be re-rendered: {e}")
    text = json.dumps(doc, ensure_ascii=False) + md
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
                if COUNT_KEY_RE.match(h) and c not in ("<10", "S", "–") and not (c.isdigit() and int(c) >= MIN_CELL):
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


def publish(A, cfg: S2Config, doc: dict, sentinels: Iterable[str] = ()) -> None:
    pj, pm = _pending(cfg.shared_json), _pending(cfg.shared_md)
    allp = [cfg.shared_json, cfg.shared_md, pj, pm]
    _remove(*allp)
    _atomic_write(pj, json.dumps(doc, indent=1, ensure_ascii=False))
    _atomic_write(pm, render_markdown(doc))
    fails = release_gate(A, pj, pm, "pending", sentinels)
    if fails:
        _remove(*allp)
        raise ReleaseGateError("; ".join(fails[:20]))
    final = json.loads(pj.read_text(encoding="utf-8"))
    final["release_gate"] = "ok"
    _atomic_write(cfg.shared_json, json.dumps(final, indent=1, ensure_ascii=False))
    _atomic_write(cfg.shared_md, render_markdown(final))
    _remove(pj, pm)
    fails = release_gate(A, cfg.shared_json, cfg.shared_md, "ok", sentinels)
    if fails:
        _remove(*allp)
        raise ReleaseGateError("; ".join(fails[:20]))


# =============================================================================
# Orchestration
# =============================================================================
def write_reconciliation(out: Out, name: str, info: dict) -> None:
    out.partition(("reconciliation", name), dict(info["row_classes"]))


def run_audit(cfg: S2Config, sentinels: Iterable[str] = ()) -> dict:
    A = load_v17(cfg.v17_module_path, cfg.v17_module_sha256)
    try:
        enc, info = load_cohort(cfg, A)
        maps, map_hashes = load_maps(cfg, A)
        cm = parse_channel_map(maps["channel"])
        dmap = parse_drug_map(maps["drug"])
        smap = parse_specimen_map(maps["specimen"])
        omap = parse_organism_map(maps["organism"])
        neuro_match = neuro_matcher(maps["neuro"])
        raw_tables = {}
        for name in list(RAW_STAGE1_SHA256) + list(RAW_RECORD_ONLY):
            raw_tables[name] = gate_raw(cfg, name)
        log(cfg, "raw table hash checks: " + ", ".join(f"{k} OK" for k in RAW_STAGE1_SHA256) +
            " (Stage 1 SHA-256); recorded only: " + ", ".join(RAW_RECORD_ONLY))
        checks = info["checks"]
        stays, gaps, ii = read_icustays(cfg, A, enc)
        hfc = enc["h_first_careunit"].astype(str)
        fu = enc["enc_i"].map(ii["first_unit"]).fillna("")
        enc["first_careunit"] = np.where(hfc != "", hfc, np.where(fu != "", fu, "unknown"))
        checks["icustays"] = dict(ii["info"])
        if (hfc != "").any():
            checks["first_careunit_harmonization_vs_icustays_disagreements"] = "see cohort_extra"
        co = Cohort.build(enc, gaps)
        gap_dis = None
        if (enc["h_gap_count"].astype(str) != "").any():
            hg = pd.to_numeric(enc["h_gap_count"], errors="coerce").fillna(0).to_numpy() > 0
            rg = np.zeros(co.n, dtype=bool)
            rg[list(gaps)] = True
            gap_dis = int((hg != rg).sum())
        prior_dis = None
        if (enc["h_prior"].astype(str) != "").any():
            o = enc.sort_values(["SUBJECT_ID", "t0", "encounter_id"])
            pr = np.zeros(co.n, dtype=int)
            pr[o["enc_i"].to_numpy()] = o.groupby("SUBJECT_ID").cumcount().to_numpy()
            prior_dis = int((pd.to_numeric(enc["h_prior"], errors="coerce").fillna(-1).to_numpy() != pr).sum())
        pt = read_patients(cfg, A, enc)
        ad = read_admissions(cfg, A, enc)
        la = pt.merge(ad, on="enc_i", how="left").sort_values("enc_i").reset_index(drop=True)
        la["age_band"] = age_band(la["age_days"].to_numpy())
        la["hours_admit_to_icu"] = hours(enc["t0"] - la["ADMITTIME"])
        dx = [norm_icd(x) for x in la["dx"].fillna("")]
        la["dx_norm"] = [a for a, _ in dx]
        la["dx_status"] = [b for _, b in dx]
        la["neuro"] = [neuro_match(a) if st in ("ok", "multiple") else False for a, st in dx]
        surg, sinfo = surgery_frame(cfg, A, enc)
        surgical = np.zeros(co.n, dtype=bool)
        rec = surg[(surg["link_class"] == "linked") & (surg["timing"] == "recent_layer_a")]
        surgical[rec["enc_i"].astype(int).to_numpy()] = True
        strata = {"age_band": la["age_band"].to_numpy(), "surgical_status": np.where(surgical, "recent_surgery", "no_recent_surgery"),
                  "first_careunit": enc["first_careunit"].to_numpy()}
        died = enc["died"].to_numpy(dtype=bool)
        out = Out(A)
        dens = register_denominators(out, co, strata, cfg.outcome_descriptors, died)
        chart_items = {it for (tb, it) in cm["items"] if tb == "CHARTEVENTS"}
        lab_items = {it for (tb, it) in cm["items"] if tb == "LABEVENTS"}
        chart, cinfo, _ = read_events(cfg, A, co, "CHARTEVENTS", chart_items, 0.0, stays)
        lab, linfo, _ = read_events(cfg, A, co, "LABEVENTS", lab_items, LAB_RELEASE_H, None)
        for name, inf in (("CHARTEVENTS", cinfo), ("LABEVENTS", linfo)):
            write_reconciliation(out, name, inf)
        if cinfo["icustay_check"] is not None:
            ic = cinfo["icustay_check"]
            t = out.cell(("reconciliation", "CHARTEVENTS_icustay_id", "n_rows_in_stay_kept_checked"), ic["rows_in_stay_checked"])
            out.cell(("reconciliation", "CHARTEVENTS_icustay_id", "n_rows_icustay_blank"), ic["rows_icustay_blank"], parent=t)
            out.cell(("reconciliation", "CHARTEVENTS_icustay_id", "n_rows_icustay_not_in_encounter"), ic["rows_icustay_not_in_encounter"], parent=t)
        vb = build_values(chart, lab, cm)
        del chart, lab
        mc = write_map_checks(out, co, dens, vb, cm, la)
        provisional: dict[str, dict] = {"channels": {}, "drug_groups": {}, "specimen_groups": {}}
        occ12: dict[str, np.ndarray] = {}
        detail_set: set[str] = set()
        for ch in sorted(cm["channels"]):
            c = cm["channels"][ch]
            ev = vb.events[ch]
            is_lab = c["source"] == "lab" or ch == "pf_ratio"
            alt = {off: ev["tc"].to_numpy() + off for off in (0.0, 2.0)} if is_lab else None
            F = event_features(co, ev["enc_i"].to_numpy(), ev["t"].to_numpy(), alt)
            detail = int((co.masks["train"] & F["ever_w"]).sum()) >= DETAIL_MIN_TRAIN
            fam, cap = locf_family(ch, c["triggered"])
            write_event_metrics(out, co, dens, ("channels", ch), F, detail=detail, cap_h=cap, lab_offsets=is_lab,
                                strata=strata, outcome_on=cfg.outcome_descriptors, died=died)
            if c["merge"]:
                mi = cm["merge_main"][ch]
                em = ev[ev["item"] == mi]
                Fm = event_features(co, em["enc_i"].to_numpy(), em["t"].to_numpy())
                for s in SPLITS:
                    b = ("channels", ch, s, "main_item_only_sensitivity")
                    out.cell(b + ("n_ever_336h",), (co.masks[s] & Fm["ever_w"]).sum(), parent=dens.enc[s])
                    for L in LANDMARKS:
                        out.cell(b + (f"{L}h", "n_avail_full"), (co.masks[s] & co.elig[L] & Fm["avail"][1.0][("full", L)]).sum(),
                                 parent=dens.elig[s][L])
            if detail:
                detail_set.add(ch)
                occ12[ch] = F["occ"][:, :CO_BINS] > 0
            pc = provisional_class(co, F, c["triggered"], ch in SCORE_INPUTS, ch in VITALS)
            pc.update({"source": c["source"], "locf_family": fam, "time_since_last_cap_hours": cap, "detail": detail,
                       "merge_pending": c["merge"], "system": c["system"]})
            if detail:
                over = {}
                for L in LANDMARKS:
                    E = co.masks["train"] & co.elig[L] & F["avail"][1.0][("full", L)]
                    ages = F["last_age"][L][E]
                    over[f"{L}h"] = bool((ages > cap).mean() > 0.5) if len(ages) >= MIN_CELL else None
                pc["locf_cap_exceeded_by_majority"] = over
            provisional["channels"][ch] = pc
        write_comeasurement(out, co, vb, cm, occ12, detail_set)
        del occ12
        labs_sparse = [ch for ch, v in provisional["channels"].items()
                       if v["source"] == "lab" and v.get("sparsity_class_72h") == "sparse"]
        docs, pinfo = read_prescriptions(cfg, A, co)
        write_reconciliation(out, "PRESCRIPTIONS", pinfo)
        sc = pinfo["start_check"]
        b = ("map_checks", "prescription_start_check")
        tot = out.cell(b + ("n_cohort_rows_timed",), sc["cohort_rows"])
        out.partition(b + ("position",), {k: sc[k] for k in ("in_stay", "before_t0_same_date", "before_t0_earlier_date",
                                                             "at_or_after_out")}, total=tot)
        out.partition(b + ("start_time_kind",), {k: sc[k] for k in ("clock_time", "date_only")}, total=tot)
        out.cell(b + ("n_midnight_clock",), sc["midnight_clock"], parent=tot)
        out.cell(b + ("n_in_stay_date_only",), sc["in_stay_date_only"], parent=tot)
        out.cell(b + ("n_in_stay_midnight_clock",), sc["in_stay_midnight_clock"], parent=tot)
        out.cell(b + ("n_date_only_on_t0_date_before_t0",), sc["date_only_on_t0_date_before_t0"], parent=tot)
        docs = resolve_drug_groups(docs, dmap, la["neuro"].to_numpy(dtype=bool))
        local: dict[str, Any] = {}
        mres = write_medications(out, co, dens, docs, strata, cfg.outcome_descriptors, died, local)
        provisional["drug_groups"] = mres["provisional"]
        mev, minfo_ = read_micro(cfg, A, co)
        write_reconciliation(out, "MICROBIOLOGYEVENTS", minfo_)
        col, minfo = micro_collections(mev, smap, omap)
        provisional["specimen_groups"] = write_micro(out, co, dens, mev, col, minfo, strata, cfg.outcome_descriptors, died,
                                                     minfo_["fill"])
        write_surgery(out, co, dens, surg, sinfo, strata, cfg.outcome_descriptors, died)
        pre_tables = {"CHARTEVENTS": cinfo["pre_any"], "LABEVENTS": linfo["pre_any"], "PRESCRIPTIONS": pinfo["pre_any"],
                      "MICROBIOLOGYEVENTS": minfo_["pre_any"]}
        write_layer_a(out, co, dens, la, pre_tables, surg, neuro_match, mc["weight_any"], gap_dis, prior_dis,
                      weight_source=mc.get("weight_source"))
        notes = {
            "as_of_48h_identity": mres["as_of_48h_identity"],
            "lab_channels_sparse_at_6h_train": len(labs_sparse),
            "lab_12h_24h_sensitivity_triggered": len(labs_sparse) >= 5,
            "merge_pending_channels": "primary metrics use all items of the merge group; main-item-only sensitivity reported",
            "magnesium_ug_dl_rows": "held out of the measurements pending Vinai's decision (map_checks/magnesium_ug_dl_rows)",
            "phi1_surgery_band": "not computed (the frozen medians are not an input of this script)",
            "outcome_descriptors": "train only, descriptive, never used for spec 9" if cfg.outcome_descriptors else "off (flag not given)",
        }
        thin = [k for k in range(N_BINS) if int((co.masks["train"] & (co.out_h >= BIN_H * (k + 1))).sum()) < 100]
        notes["thin_bins_train_lt100_from_bin"] = (min(thin) if thin else "none")
        log(cfg, f"registry: {len(out.reg.raw)} cells, {len(out.reg.identities)} identities; suppressing (LP) ...")
        status, supp = apply_suppression_frozen(A, out.reg, out.frozen)
        log(cfg, f"suppression: {supp['primary_cells']} primary, {supp['complementary_cells']} complementary, "
                 f"{supp['recoverable']} recoverable")
    except A.HashMismatchError as e:
        raise HashMismatchError(str(e)) from e
    except A.AuditError as e:
        raise AuditError(str(e)) from e
    except A.ReleaseGateError as e:
        raise ReleaseGateError(str(e)) from e
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    local.update({"WARNING": "LOCAL ONLY - unsuppressed aggregates; never copy, upload or share",
                  "raw_counts": out.reg.raw, "checks": checks})
    _atomic_write(cfg.local_json, json.dumps(local, indent=1, default=str, ensure_ascii=False))
    doc = build_document(out, status, supp, info["hashes"], raw_tables, map_hashes, checks, cfg, provisional, notes)
    publish(A, cfg, doc, sentinels)
    final = json.loads(cfg.shared_json.read_text(encoding="utf-8"))
    print_console_summary(cfg, final)
    return final


def print_console_summary(cfg: S2Config, doc: dict) -> None:
    """Short; prints only published values and labels."""
    if not cfg.verbose:
        return
    print(f"=== {doc['spec']} | release_gate={doc['release_gate']} | outcome_descriptors={doc['outcome_descriptors_flag']} ===")
    print("hashes OK: " + ", ".join(f"{k} {v[:8]}" for k, v in doc["input_hashes"].items()) +
          f", v1.7 module {doc['v17_module']['sha256'][:8]}")
    print("maps OK: " + ", ".join(f"{k} {v['sha256'][:8]}" for k, v in doc["maps"].items()))
    print("raw tables: " + ", ".join(f"{k} {v['sha256'][:8]} ({'gated' if v['gate'] == 'stage1_sha256' else 'recorded'})"
                                     for k, v in doc["raw_tables"].items()))
    print("cohort and landmark eligibility reproduce: " +
          ", ".join(f"{k}={v}" for k, v in doc["cohort_checks"]["eligible_at_landmark"].items()))
    cc: dict[str, int] = {}
    for v in doc["provisional_classes"]["channels"].values():
        cc[v["class"]] = cc.get(v["class"], 0) + 1
    print("provisional channel classes (train, spec 9 defaults): " + ", ".join(f"{k}={cc[k]}" for k in sorted(cc)))
    print(f"S3 as-of 48 h identity: {doc['notes']['as_of_48h_identity']}; sparse lab channels at 6 h: "
          f"{doc['notes']['lab_channels_sparse_at_6h_train']}")
    s = doc["suppression"]
    print(f"suppression: primary={s['primary_cells']} complementary={s['complementary_cells']} recoverable={s['recoverable']}")
    db = _g(doc, "layer_a", "dob_breakdown")
    if isinstance(db, dict):
        st = db.get("status") or {}
        print("DOB breakdown (D20): " +
              ", ".join(f"{k}={st.get(f'n_{k}')}" for k in DOB_STATUS) +
              f"; invalid={db.get('n_invalid')}; safety_gate={db.get('safety_gate')}")
    print(f"shared: {cfg.shared_json.name}, {cfg.shared_md.name}  (local only, never share: {cfg.local_json.name})")
    print("release gate ok")


def build_config(argv: list[str] | None = None) -> S2Config:
    here = Path(__file__).resolve()
    default_root = Path(os.environ.get("FOUNDATIONS_ROOT", str(here.parents[1])))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", default=str(default_root))
    ap.add_argument("--raw-dir", default=os.environ.get("PIC_RAW_DIR"))
    ap.add_argument("--split-file", default=os.environ.get("SPLIT_FILE_PATH"))
    ap.add_argument("--canonical", default=os.environ.get("CANONICAL_PATH"))
    ap.add_argument("--harmonization", default=os.environ.get("HARMONIZATION_PATH"))
    ap.add_argument("--complexity-map", default=os.environ.get("COMPLEXITY_MAP_PATH"))
    ap.add_argument("--maps-dir", default=os.environ.get("MAPS_DIR"))
    ap.add_argument("--out-dir", default=os.environ.get("AUDIT_OUT_DIR"))
    ap.add_argument("--harm-t0-col", default=os.environ.get("HARM_T0_COL"))
    ap.add_argument("--harm-out-col", default=os.environ.get("HARM_OUT_COL"))
    ap.add_argument("--split-col", default=os.environ.get("SPLIT_COL", "split"))
    ap.add_argument("--v17-module", default=os.environ.get("V17_MODULE_PATH", str(here.parent / V17_MODULE_FILE)))
    ap.add_argument("--outcome-descriptors", action="store_true",
                    help="compute the train-only outcome descriptors (spec 8); off by default; the state is recorded")
    ap.add_argument("--chunksize", type=int, default=CHUNKSIZE)
    a = ap.parse_args(argv)
    if not a.raw_dir or not a.split_file:
        raise AuditError("--raw-dir and --split-file are required (or set PIC_RAW_DIR and SPLIT_FILE_PATH); no local default path is shipped")
    root = Path(a.project_root)
    interim = root / "data" / "interim"
    return S2Config(
        raw_dir=Path(a.raw_dir),
        canonical_path=Path(a.canonical) if a.canonical else interim / "canonical_icu_encounters.csv",
        harmonization_path=Path(a.harmonization) if a.harmonization else interim / "encounter_harmonization_fields.csv",
        split_path=Path(a.split_file),
        complexity_map_path=Path(a.complexity_map) if a.complexity_map else root / "config" / "surgery_complexity_map_v1.csv",
        maps_dir=Path(a.maps_dir) if a.maps_dir else interim / "maps",
        out_dir=Path(a.out_dir) if a.out_dir else interim, v17_module_path=Path(a.v17_module),
        outcome_descriptors=bool(a.outcome_descriptors), harm_t0_col=a.harm_t0_col, harm_out_col=a.harm_out_col,
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
        run_audit(cfg)
    except ReleaseGateError as e:
        print(f"RELEASE GATE FAILED - no shared file left in place: {e}", file=sys.stderr)
        return 3
    except AuditError as e:
        print(f"AUDIT STOPPED: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
