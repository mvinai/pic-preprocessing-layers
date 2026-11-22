#!/usr/bin/env python3
"""
Surgery complexity vs duration audit -- clean reference implementation.

Implements `specs/surgery_complexity_duration_audit_spec.md` **v1.7** together with
`surgery_complexity_duration_audit_FIXES_v1.md` (Amendment A: signed D2; Amendment B: ROW_ID surgery key,
OPER_ID descriptor only, harmonization OUTTIME fallback for OUT; Amendment C: surgery timing by ANESTHESIA END,
Layer A recent window [t0 - 24 h, t0 + 6 h], Layer C after t0 + 6 h up to OUT incl. merge gaps).
Replaces Copilot's run-1 script. Read-only and aggregate-only.

Usage (Windows, from the project root):
    python scripts/audit_surgery_complexity_duration.py
    python scripts/audit_surgery_complexity_duration.py --raw-dir "C:/.../PICDataset" --split-file "C:/.../subject_splits.csv"
Every path can also be set by environment variable (see `build_config`).

Outputs (data/interim/):
    surgery_complexity_duration_audit_MAP_v1.json / .md   shared (suppressed); written ONLY if the release gate passes
    surgery_complexity_duration_audit_MAP_v1_raw_LOCAL.json   local only, never share
Exit codes: 0 ok; 2 input/hash/reconciliation failure; 3 release gate failure (no shared file left in place).
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import os
import re
import sys
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sps
from scipy.optimize import linprog

# =============================================================================
# CONFIG BLOCK (spec v1.6; criteria unchanged since v1.5). Thresholds are read from here, never hard-coded in logic.
# =============================================================================
SPEC_VERSION = "surgery_complexity_duration_audit v1.7"
MAP_LABEL = "MAP_v1"
MAP_STATUS = "reviewed_v1"
MAP_REVIEWED_DATE = "2026-11-10"
VERDICT_LABEL = "approved_criteria_2026-11-10_D2signed"

PINNED_SHA256 = {
    "canonical": "e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418",
    "harmonization": "ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d",
    "complexity_map": "5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5",
    "split": "60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31",
}

# Approved discrepancy criteria (Vinai, Nov 10, 2026; D2 signed per spec v1.5).
CRITERIA_SIGNED = True
CRITERIA_APPROVED_DATE = "2026-11-10"
D2_RULE = "signed"
TAU_B_SIGNED_MAX = 0.2          # D2 met if signed tau-b < TAU_B_SIGNED_MAX
D3_OVERLAP_MAX = 0.25           # D3 met if share > D3_OVERLAP_MAX
D4_MIRROR_OVERLAP_MAX = 0.25    # D4 (reported only)
D5_DISCORDANT_NAMES_MAX = 0.20  # D5 (reported only)

MIN_CELL = 10                   # primary suppression / eligibility / quantile threshold
SCREEN_LPS = 8                  # recoverability solver: screening LPs per component (speed only; exact LPs decide)
SCREEN_CAP = 1e7                # cap for open upper bounds in screening LPs only
CORR_MIN_TOTAL = 30             # correlation published only if n >= 30 and every group >= 10
D_GT_HOURS = 24.0
# v1.7 (Vinai, Nov 11, 2026): timing by anaesthesia end T = ANES_END_TIME (fallback SURGERY_END_TIME, counted).
RECENT_BEFORE_T0_HOURS = 24.0   # Layer A recent: t0 - 24 h <= T <= t0 + 6 h (Layer A ONLY, never also Layer C)
RECENT_AFTER_T0_HOURS = 6.0     # Layer C: t0 + 6 h < T <= OUT (merge gaps included); remote: T < t0 - 24 h

DEFAULT_RAW_DIR = os.environ.get("PIC_RAW_DIR")
DEFAULT_SPLIT_FILE = os.environ.get("SPLIT_FILE_PATH")

# Harmonization columns: the spec names them only as "t0" and "OUT" (last OUTTIME). The exact header is
# resolved from these candidates (exactly one must match) unless set explicitly (--harm-t0-col/--harm-out-col).
HARM_T0_CANDIDATES = ("t0", "T0", "t0_harmonized", "harmonized_t0", "first_intime", "FIRST_INTIME", "t0_intime")
HARM_OUT_CANDIDATES = ("OUT", "out", "last_outtime", "LAST_OUTTIME", "harmonized_out", "OUT_harmonized", "t_out", "last_OUTTIME")
# v1.6: the real encounter_harmonization_fields.csv carries t0 and OUTTIME (no OUT column). OUTTIME is a FALLBACK
# default: used only when no primary OUT candidate is present (an explicit --harm-out-col always wins).
HARM_OUT_FALLBACK = ("OUTTIME",)
HARM_T0_FALLBACK: tuple = ()
SURGERY_KEY = "ROW_ID"  # v1.6 (Vinai, Nov 11): OPER_ID is not unique within an admission; ROW_ID is the surgery key

SPLITS = ("train", "validation", "test")
ALL_SPLITS = SPLITS + ("pooled",)
SPLIT_LABELS = {"train": "train", "training": "train", "val": "validation", "valid": "validation",
                "validation": "validation", "dev": "validation", "test": "test", "testing": "test"}
BANDS = ("lt2h", "2to4h", "gt4h")
TIMING = ("recent_layer_a", "remote_same_hadm", "during_icu_layer_c")
LINK_CLASSES = ("linked", "post_icu_same_hadm", "other_hadm", "not_in_cohort", "invalid_time",
                "not_linkable_missing_encounter_time")
INVALID_REASONS = ("begin_after_end", "zero_duration", "missing_window", "fallback_nonpositive")

SURGERY_INFO_COLUMNS = ["ROW_ID", "SUBJECT_ID", "HADM_ID", "VISIT_ID", "OPER_ID", "ANES_START_TIME", "ANES_END_TIME",
                        "SURGERY_BEGIN_TIME", "SURGERY_END_TIME", "SURGERY_NAME", "ANES_METHOD", "SURGERY_POSITION"]
SURGERY_INFO_EXPECTED_MISSING = {c: 0 for c in SURGERY_INFO_COLUMNS} | {"ANES_START_TIME": 1, "ANES_END_TIME": 1}
CANONICAL_REQUIRED = ["encounter_id", "SUBJECT_ID", "hadm_id_list", "HADM_ID", "INTIME", "OUTTIME"]
MAP_COLUMNS = ["raw_name", "standardized_name", "body_system", "intensity", "I_group", "STAT", "STAT_group",
               "complexity_group", "review_flag", "n_components", "STAT_group_crosses"]

CARDIAC = "Cardiac"
UNCLASSIFIED = "Unclassified"
UNMAPPED = "__unmapped__"
G_ORDER = ["G1", "G2", "G3"]
I_ORDER = ["I1", "I2", "I3"]
STAT_ORDER = ["1", "2", "3", "4", "5"]
INTENSITY_ORDER = ["minor/diagnostic", "catheter/interventional", "minimally invasive", "major open"]
INTENSITY_TO_I = {"minor/diagnostic": "I1", "catheter/interventional": "I1", "minimally invasive": "I2", "major open": "I3"}
STAT_TO_G = {"1": "G1", "2": "G1", "3": "G2", "4": "G3", "5": "G3", "NC": "NC"}

# Shared-output count fields are recognised by key name (gate cross-checks against the declared list).
COUNT_KEYS = frozenset(
    {"n", "total", "count", "complement", "linked", "analysed", "unmapped", "excluded_unclassified",
     "valid_duration", "invalid_total", "duration_fallback_anes", "d_gt_24h", "straddling_t0", "recent_ending_after_t0",
     "timing_fallback_surgery_end", "pre_admit_time", "valid", "invalid", "surgery_info_rows", "valid_total", "n_multi_candidate",
     "resolved_rule1_during", "resolved_rule2_next_t0", "resolved_rule3_post_icu", "earlier_hadm", "later_hadm",
     "shared_oper_id_groups", "rows_in_shared_oper_id_groups", "possible_true_duplicate_groups",
     "rows_in_possible_true_duplicate_groups"}
    | set(BANDS) | set(TIMING) | set(LINK_CLASSES) | set(INVALID_REASONS))


class AuditError(Exception):
    """Input, hash or reconciliation failure (exit 2)."""


class HashMismatchError(AuditError):
    pass


class ReconciliationError(AuditError):
    pass


class ReleaseGateError(Exception):
    """Release gate failure (exit 3)."""


@dataclass
class AuditConfig:
    canonical_path: Path
    harmonization_path: Path
    map_path: Path
    split_path: Path
    surgery_info_path: Path
    out_dir: Path
    admissions_path: Path | None = None
    pinned_sha256: dict = field(default_factory=lambda: dict(PINNED_SHA256))
    harm_t0_col: str | None = None
    harm_out_col: str | None = None
    split_subject_col: str = "SUBJECT_ID"
    split_col: str = "split"
    expected_surgery_rows: int | None = None  # public copy: pin omitted
    expected_encounters: int | None = None  # public copy: pin omitted
    expected_subjects: int | None = None
    expected_split_encounters: dict = field(default_factory=dict)
    expected_map_rows: int = 1495
    expected_missing: dict = field(default_factory=lambda: dict(SURGERY_INFO_EXPECTED_MISSING))
    bootstrap_reps: int = 2000
    seed: int = 20261110
    verbose: bool = True

    @property
    def shared_json(self) -> Path:
        return self.out_dir / f"surgery_complexity_duration_audit_{MAP_LABEL}.json"

    @property
    def shared_md(self) -> Path:
        return self.out_dir / f"surgery_complexity_duration_audit_{MAP_LABEL}.md"

    @property
    def local_json(self) -> Path:
        return self.out_dir / f"surgery_complexity_duration_audit_{MAP_LABEL}_raw_LOCAL.json"


def criteria_config() -> dict:
    return {"signed": CRITERIA_SIGNED, "approved_date": CRITERIA_APPROVED_DATE, "decisive": ["D1", "D2", "D3"],
            "reported_only": ["D4", "D5"], "d2_rule_name": D2_RULE, "tau_b_signed_max": TAU_B_SIGNED_MAX,
            "d2_rule": "met if tau_b < tau_b_signed_max (signed)", "overlap_max": D3_OVERLAP_MAX,
            "mirror_overlap_max": D4_MIRROR_OVERLAP_MAX, "discordant_names_max": D5_DISCORDANT_NAMES_MAX}


def validate_criteria() -> None:
    """Acceptance check 10: refuse to run if the config differs from the approved values."""
    approved = {"CRITERIA_SIGNED": True, "D2_RULE": "signed", "TAU_B_SIGNED_MAX": 0.2, "D3_OVERLAP_MAX": 0.25,
                "D4_MIRROR_OVERLAP_MAX": 0.25, "D5_DISCORDANT_NAMES_MAX": 0.20, "CRITERIA_APPROVED_DATE": "2026-11-10"}
    current = {"CRITERIA_SIGNED": CRITERIA_SIGNED, "D2_RULE": D2_RULE, "TAU_B_SIGNED_MAX": TAU_B_SIGNED_MAX,
               "D3_OVERLAP_MAX": D3_OVERLAP_MAX, "D4_MIRROR_OVERLAP_MAX": D4_MIRROR_OVERLAP_MAX,
               "D5_DISCORDANT_NAMES_MAX": D5_DISCORDANT_NAMES_MAX, "CRITERIA_APPROVED_DATE": CRITERIA_APPROVED_DATE}
    bad = [k for k in approved if current[k] != approved[k]]
    if bad:
        raise AuditError(f"criteria config differs from approved spec v1.5 values: {bad}")


# =============================================================================
# Small utilities
# =============================================================================
def log(cfg: "AuditConfig | None", msg: str) -> None:
    if cfg is None or cfg.verbose:
        print(msg, flush=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_hash(label: str, path: Path, expected: str) -> str:
    """Hash is computed from the file bytes BEFORE the file is read."""
    path = Path(path)
    if not path.exists():
        raise AuditError(f"{label}: file not found: {path}")
    if not (isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected)):
        raise HashMismatchError(f"{label}: pinned hash is not a full 64-character SHA-256")
    got = sha256_file(path)
    if got != expected:
        raise HashMismatchError(f"{label}: SHA-256 mismatch (expected {expected}, got {got})")
    return got


def norm_id(x: Any) -> str:
    """Canonical string ID: '123', 123, 123.0, '123.0', ' 123 ' -> '123'. Missing -> ''."""
    if x is None:
        return ""
    if isinstance(x, (float, np.floating)):
        if math.isnan(x):
            return ""
        if float(x).is_integer():
            return str(int(x))
        return repr(float(x))
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    s = str(x).strip().strip("'\"").strip()
    if s.lower() in {"", "nan", "none", "null", "nat"}:
        return ""
    m = re.fullmatch(r"([+-]?\d+)\.0*", s)
    if m:
        return str(int(m.group(1)))
    return s


def parse_hadm_list(v: Any) -> list[str]:
    """Robust to Python lists, JSON '[1, 2]', Python-repr "['1', '2']" / '[1.0, 2.0]', and '1;2' / '1,2' / '1|2'."""
    if isinstance(v, (list, tuple, set, np.ndarray, pd.Series)):
        items = list(v)
    else:
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return []
        s = str(v).strip()
        if s.lower() in {"", "nan", "none", "null", "[]", "()"}:
            return []
        items = None
        if s[0] in "[(":
            for loader in (json.loads, ast.literal_eval):
                try:
                    obj = loader(s)
                except Exception:
                    continue
                items = list(obj) if isinstance(obj, (list, tuple, set)) else [obj]
                break
        if items is None:
            items = [t for t in re.split(r"[,;|\s]+", s.strip("[]()")) if t]
    out: list[str] = []
    for it in items:
        nid = norm_id(it)
        if nid and nid not in out:
            out.append(nid)
    return out


def norm_name(x: Any) -> str:
    """Same rule as the type map: trim, case-fold (lower), collapse whitespace."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return re.sub(r"\s+", " ", str(x).strip().lower())


def read_csv_str(path: Path, usecols: list[str] | None = None) -> pd.DataFrame:
    """Every CSV is read as text; parsing is explicit afterwards. Empty cells stay ''."""
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False, usecols=usecols,
                       encoding="utf-8-sig")


def read_header(path: Path) -> list[str]:
    return list(pd.read_csv(path, dtype=str, nrows=0, encoding="utf-8-sig").columns)


def parse_time(s: pd.Series) -> tuple[pd.Series, pd.Series, pd.Series]:
    """Returns (datetime series, empty mask, non-empty-but-unparseable mask)."""
    s = s.astype(object).map(lambda v: "" if v is None else str(v).strip())
    empty = s.isin(["", "nan", "NaN", "NaT", "None", "null", "NULL"])
    vals = s.where(~empty, None)
    try:
        t = pd.to_datetime(vals, errors="coerce", format="mixed")
    except (TypeError, ValueError):
        t = pd.to_datetime(vals, errors="coerce")
    unparseable = (~empty) & t.isna()
    return t, empty, unparseable


def band_of(d: float) -> str:
    if d < 2:
        return "lt2h"
    if d <= 4:
        return "2to4h"
    return "gt4h"


def P(*parts: Any) -> str:
    """JSON pointer (RFC 6901) for a cell id; body systems such as 'GI/abdominal' are escaped."""
    return "/" + "/".join(str(p).replace("~", "~0").replace("/", "~1") for p in parts)


def ptr_parts(ptr: str) -> list[str]:
    return [p.replace("~1", "/").replace("~0", "~") for p in ptr.split("/")[1:]]


def set_ptr(doc: dict, ptr: str, value: Any) -> None:
    parts = ptr_parts(ptr)
    cur = doc
    for p in parts[:-1]:
        cur = cur.setdefault(p, {})
    cur[parts[-1]] = value


def get_ptr(doc: dict, ptr: str) -> Any:
    cur: Any = doc
    for p in ptr_parts(ptr):
        if not isinstance(cur, dict) or p not in cur:
            raise KeyError(ptr)
        cur = cur[p]
    return cur


def r1(x: float) -> float:
    return float(round(float(x), 1))


def r3(x: float) -> float:
    return float(round(float(x), 3))


# =============================================================================
# Group mapping, durations, verdict (pure functions; unit tested)
# =============================================================================
def primary_group(body_system: str, stat: str, intensity: str, map_group: str) -> str:
    """Cardiac uses STAT (never intensity); others use the 3-group intensity collapse."""
    if body_system == CARDIAC:
        g = STAT_TO_G.get(stat)
        if g is None:
            raise AuditError("cardiac map row with STAT outside 1-5/NC")
    else:
        g = INTENSITY_TO_I.get(intensity)
        if g is None:
            raise AuditError("non-cardiac map row with unknown intensity")
    if map_group and map_group != g:
        raise AuditError("map complexity_group disagrees with STAT/intensity rule (map must not be edited)")
    return g


def compute_durations(B: pd.Series, E: pd.Series, AS: pd.Series, AE: pd.Series) -> pd.DataFrame:
    """D = E - B (hours). Fallback AE - AS only when B or E is missing/unparseable. Invalid reasons counted."""
    n = len(B)
    D = np.full(n, np.nan)
    reason = np.array([""] * n, dtype=object)
    fallback = np.zeros(n, dtype=bool)
    both = (B.notna() & E.notna()).to_numpy()
    d_be = ((E - B).dt.total_seconds() / 3600.0).to_numpy()
    d_an = ((AE - AS).dt.total_seconds() / 3600.0).to_numpy()
    an_ok = (AS.notna() & AE.notna()).to_numpy()
    for i in range(n):
        if both[i]:
            if d_be[i] > 0:
                D[i] = d_be[i]
            elif d_be[i] == 0:
                reason[i] = "zero_duration"
            else:
                reason[i] = "begin_after_end"
        elif an_ok[i]:
            if d_an[i] > 0:
                D[i] = d_an[i]
                fallback[i] = True
            else:
                reason[i] = "fallback_nonpositive"
        else:
            reason[i] = "missing_window"
    return pd.DataFrame({"D": D, "invalid_reason": reason, "fallback": fallback}, index=B.index)


def evaluate_criteria(ordered_medians: list[float] | None, tau_b: float | None, d3_share: float | None,
                      d4_share: float | None, d5_share: float | None, criteria: dict | None = None) -> dict:
    """Verdict semantics (spec v1.4/v1.5): 'met' = the discrepancy condition holds.
    D2 uses the SIGNED tau-b: met iff tau_b < tau_b_signed_max. NaN/None -> not_assessable (never not_met).
    big_discrepancy = Y iff any of D1-D3 is met; N if none met and at least one computable; else not_assessable."""
    c = criteria or criteria_config()

    def finite(v: Any) -> bool:
        return v is not None and not isinstance(v, bool) and isinstance(v, (int, float, np.floating, np.integer)) and bool(np.isfinite(v))

    if ordered_medians is None or len(ordered_medians) < 2 or not all(finite(m) for m in ordered_medians):
        d1 = "not_assessable"
    else:
        strictly_increasing = all(b > a for a, b in zip(ordered_medians, ordered_medians[1:]))
        d1 = "not_met" if strictly_increasing else "met"
    d2 = "not_assessable" if not finite(tau_b) else ("met" if tau_b < c["tau_b_signed_max"] else "not_met")
    d3 = "not_assessable" if not finite(d3_share) else ("met" if d3_share > c["overlap_max"] else "not_met")
    d4 = "not_assessable" if not finite(d4_share) else ("met" if d4_share > c["mirror_overlap_max"] else "not_met")
    d5 = "not_assessable" if not finite(d5_share) else ("met" if d5_share > c["discordant_names_max"] else "not_met")
    decisive = [d1, d2, d3]
    if "met" in decisive:
        big = "Y"
    elif "not_met" in decisive:
        big = "N"
    else:
        big = "not_assessable"
    return {"D1": d1, "D2": d2, "D3": d3, "D4": d4, "D5": d5, "big_discrepancy": big, "label": VERDICT_LABEL}


def kendall_tau_b_surgery_level(group_codes: np.ndarray, durations: np.ndarray) -> float:
    """Signed Kendall tau-b between ordinal group code and duration, one row per SURGERY (never group medians)."""
    x = np.asarray(group_codes, dtype=float)
    y = np.asarray(durations, dtype=float)
    if x.shape != y.shape or len(x) < 3 or np.isnan(x).any() or np.isnan(y).any():
        return float("nan")
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return float("nan")
    res = sps.kendalltau(x, y, variant="b")
    return float(res[0])


def spearman_rho(group_codes: np.ndarray, durations: np.ndarray) -> float:
    x = np.asarray(group_codes, dtype=float)
    y = np.asarray(durations, dtype=float)
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return float("nan")
    return float(sps.spearmanr(x, y)[0])


def prob_superiority(high: np.ndarray, low: np.ndarray) -> float:
    """P(D_high > D_low) + 0.5 P(tie)."""
    if len(high) == 0 or len(low) == 0:
        return float("nan")
    u = sps.mannwhitneyu(high, low, alternative="two-sided").statistic
    return float(u) / (len(high) * len(low))


class Bootstrapper:
    """Subject-cluster percentile bootstrap. Realised replicate counts are recorded (metadata honesty)."""

    def __init__(self, reps: int, seed: int):
        self.reps = int(reps)
        self.seed = int(seed)
        self.runs = 0
        self.realised: list[int] = []

    @staticmethod
    def cluster_weights(inv: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
        """Resample SUBJECTS with replacement; every surgery of a drawn subject comes along as a unit."""
        draw = rng.integers(0, k, size=k)
        return np.bincount(draw, minlength=k)[inv]

    def ci(self, subjects: np.ndarray, stat: Callable[[np.ndarray], float], key: str) -> dict:
        if self.reps <= 0:
            return {"ci": "insufficient", "reason": "bootstrap disabled", "realised_replicates": 0}
        uniq, inv = np.unique(np.asarray(subjects).astype(str), return_inverse=True)
        k = len(uniq)
        if k < 2:
            return {"ci": "insufficient", "reason": "fewer than 2 subjects", "realised_replicates": 0}
        rng = np.random.default_rng([self.seed, zlib.crc32(key.encode("utf-8"))])
        base = np.arange(len(inv))
        vals = []
        for _ in range(self.reps):
            w = self.cluster_weights(inv, k, rng)
            v = stat(np.repeat(base, w))
            if v is not None and np.isfinite(v):
                vals.append(v)
        self.runs += 1
        self.realised.append(len(vals))
        if len(vals) < 0.95 * self.reps:
            return {"ci": "insufficient", "reason": "too many degenerate replicates", "realised_replicates": len(vals)}
        lo, hi = np.percentile(vals, [2.5, 97.5])
        return {"ci": [r3(lo), r3(hi)], "realised_replicates": len(vals)}


# =============================================================================
# Inputs
# =============================================================================
def resolve_column(columns: list[str], explicit: str | None, candidates: Iterable[str], label: str,
                   fallback: Iterable[str] = ()) -> str:
    if explicit:
        if explicit not in columns:
            raise AuditError(f"harmonization column for {label} not found: {explicit!r}; header has {columns}")
        return explicit
    candidates = tuple(candidates)
    found = [c for c in candidates if c in columns]
    if not found:
        found = [c for c in fallback if c in columns]
        candidates = candidates + tuple(fallback)
    if len(found) != 1:
        raise AuditError(f"cannot resolve the harmonization {label} column (candidates {list(candidates)} matched "
                         f"{found}); header has {columns}. Set it with --harm-{label.lower()}-col.")
    return found[0]


def load_inputs(cfg: AuditConfig) -> dict:
    validate_criteria()
    hashes = {
        "canonical": verify_hash("canonical", cfg.canonical_path, cfg.pinned_sha256["canonical"]),
        "harmonization": verify_hash("harmonization", cfg.harmonization_path, cfg.pinned_sha256["harmonization"]),
        "complexity_map": verify_hash("complexity_map", cfg.map_path, cfg.pinned_sha256["complexity_map"]),
        "split": verify_hash("split", cfg.split_path, cfg.pinned_sha256["split"]),
    }
    log(cfg, "input hash checks: canonical OK, harmonization OK, complexity_map OK, split OK (all 4 exact 64-char)")
    info: dict[str, Any] = {"hashes": hashes, "checks": {}}

    # ---- SURGERY_INFO: exact headers, row count, missingness
    if not Path(cfg.surgery_info_path).exists():
        raise AuditError(f"SURGERY_INFO not found: {cfg.surgery_info_path}")
    header = read_header(cfg.surgery_info_path)
    missing_cols = [c for c in SURGERY_INFO_COLUMNS if c not in header]
    if missing_cols:
        raise AuditError(f"SURGERY_INFO missing columns: {missing_cols}")
    surg = read_csv_str(cfg.surgery_info_path, usecols=SURGERY_INFO_COLUMNS)
    if len(surg) != cfg.expected_surgery_rows:
        raise AuditError(f"SURGERY_INFO rows {len(surg)} != expected {cfg.expected_surgery_rows}")
    miss = {c: int(surg[c].map(lambda v: str(v).strip() == "").sum()) for c in SURGERY_INFO_COLUMNS}
    diffs = {c: (miss[c], cfg.expected_missing.get(c, 0)) for c in SURGERY_INFO_COLUMNS if miss[c] != cfg.expected_missing.get(c, 0)}
    if diffs:
        raise AuditError(f"SURGERY_INFO missingness differs from spec (observed, expected): {diffs}")
    rid = surg[SURGERY_KEY].map(norm_id)
    if (rid == "").any():
        raise AuditError(f"SURGERY_INFO {SURGERY_KEY} is blank on {int((rid == '').sum())} rows (surgery key must be nonblank)")
    if rid.duplicated().any():
        raise AuditError(f"SURGERY_INFO {SURGERY_KEY} is not unique ({int(rid.duplicated().sum())} repeated rows)")
    info["checks"]["surgery_key"] = SURGERY_KEY
    log(cfg, f"SURGERY_INFO headers ok (12/12), row count ok, missingness ok, {SURGERY_KEY} unique and nonblank")

    # ---- Canonical encounters
    cheader = read_header(cfg.canonical_path)
    miss_c = [c for c in CANONICAL_REQUIRED if c not in cheader]
    if miss_c:
        raise AuditError(f"canonical CSV missing columns: {miss_c}")
    canon = read_csv_str(cfg.canonical_path, usecols=CANONICAL_REQUIRED)
    canon["encounter_id"] = canon["encounter_id"].map(lambda v: str(v).strip())
    if (canon["encounter_id"] == "").any() or canon["encounter_id"].duplicated().any():
        raise AuditError("canonical encounter_id empty or not unique")
    canon["SUBJECT_ID"] = canon["SUBJECT_ID"].map(norm_id)
    canon["HADM_ID"] = canon["HADM_ID"].map(norm_id)
    canon["hadm_list"] = canon["hadm_id_list"].map(parse_hadm_list)  # parsed exactly once, from the raw text
    empty_list = canon["hadm_list"].map(len) == 0
    info["checks"]["hadm_list_empty_filled_from_HADM_ID"] = int(empty_list.sum())
    canon.loc[empty_list, "hadm_list"] = pd.Series([[h] if h else [] for h in canon.loc[empty_list, "HADM_ID"]],
                                                   index=canon.index[empty_list], dtype=object)
    info["checks"]["hadm_id_not_in_hadm_list"] = int(sum(h not in lst for h, lst in zip(canon["HADM_ID"], canon["hadm_list"])))
    info["checks"]["multi_hadm_encounters"] = int((canon["hadm_list"].map(len) > 1).sum())
    cin, _, _ = parse_time(canon["INTIME"])

    # ---- Harmonization (t0 / OUT)
    hheader = read_header(cfg.harmonization_path)
    if "encounter_id" not in hheader:
        raise AuditError("harmonization file has no encounter_id column")
    t0c = resolve_column(hheader, cfg.harm_t0_col, HARM_T0_CANDIDATES, "T0", HARM_T0_FALLBACK)
    outc = resolve_column(hheader, cfg.harm_out_col, HARM_OUT_CANDIDATES, "OUT", HARM_OUT_FALLBACK)
    span_col = "span_icu_hours" if "span_icu_hours" in hheader else None
    harm = read_csv_str(cfg.harmonization_path, usecols=["encounter_id", t0c, outc] + ([span_col] if span_col else []))
    harm["encounter_id"] = harm["encounter_id"].map(lambda v: str(v).strip())
    if harm["encounter_id"].duplicated().any():
        raise AuditError("harmonization encounter_id not unique")
    if set(harm["encounter_id"]) != set(canon["encounter_id"]):
        raise AuditError("harmonization encounter_id set differs from canonical encounter_id set")
    t0, _, _ = parse_time(harm[t0c])
    out, _, _ = parse_time(harm[outc])
    if span_col:  # OUT must be the last (max) OUTTIME = t0 + span_icu_hours (band-support spec; amendment decision 6)
        span = pd.to_numeric(harm[span_col].map(lambda v: str(v).strip()), errors="coerce")
        off = ((out - t0).dt.total_seconds() / 3600.0 - span).abs() > (1.0 / 60.0)
        info["checks"]["OUT_minus_t0_ne_span_icu_hours_gt_1min"] = int(off.sum())
    harm = pd.DataFrame({"encounter_id": harm["encounter_id"], "t0": t0, "OUT": out})
    info["checks"]["harmonization_columns"] = {"t0": t0c, "OUT": outc}
    info["checks"]["harmonized_time_missing_or_unparseable"] = int(harm["t0"].isna().sum() + harm["OUT"].isna().sum())
    if ((harm["OUT"] < harm["t0"]) & harm["OUT"].notna() & harm["t0"].notna()).any():
        raise AuditError("harmonized OUT earlier than t0 for some encounters")
    enc = canon[["encounter_id", "SUBJECT_ID", "HADM_ID", "hadm_list"]].assign(canon_intime=cin).merge(harm, on="encounter_id", how="left")
    info["checks"]["t0_ne_canonical_INTIME"] = int(((enc["t0"] != enc["canon_intime"]) & enc["t0"].notna()).sum())

    # ---- Split (subject-grouped)
    sheader = read_header(cfg.split_path)
    for c in (cfg.split_subject_col, cfg.split_col):
        if c not in sheader:
            raise AuditError(f"split file missing column {c!r}; header has {sheader}")
    spl = read_csv_str(cfg.split_path, usecols=[cfg.split_subject_col, cfg.split_col])
    spl = pd.DataFrame({"SUBJECT_ID": spl[cfg.split_subject_col].map(norm_id),
                        "split": spl[cfg.split_col].map(lambda v: SPLIT_LABELS.get(str(v).strip().lower(), "?"))})
    if (spl["split"] == "?").any():
        raise AuditError("split file has unknown split labels")
    spl = spl.drop_duplicates()
    if spl["SUBJECT_ID"].duplicated().any():
        raise AuditError("a subject appears in more than one split")
    enc = enc.merge(spl, on="SUBJECT_ID", how="left")
    if enc["split"].isna().any():
        raise AuditError(f"{int(enc['split'].isna().sum())} cohort encounters have no split")
    if "split" in cheader:  # the canonical file carries its own split column: it must agree with the split file
        cs = read_csv_str(cfg.canonical_path, usecols=["encounter_id", "split"])
        cs = pd.Series(cs["split"].map(lambda v: SPLIT_LABELS.get(str(v).strip().lower(), "?")).to_numpy(),
                       index=cs["encounter_id"].map(lambda v: str(v).strip()))
        disagree = int((enc["encounter_id"].map(cs) != enc["split"]).sum())
        info["checks"]["canonical_split_disagreements"] = disagree
        if disagree:
            raise AuditError(f"canonical split column disagrees with the split file for {disagree} encounters")
    n_enc, n_subj = len(enc), enc["SUBJECT_ID"].nunique()
    split_counts = {s: int((enc["split"] == s).sum()) for s in SPLITS}
    if n_enc != cfg.expected_encounters or n_subj != cfg.expected_subjects or split_counts != cfg.expected_split_encounters:
        raise AuditError(f"cohort/split totals differ from spec: encounters {n_enc}, subjects {n_subj}, splits {split_counts}")
    log(cfg, "cohort and split totals ok")

    # ---- Complexity map (name-level)
    mheader = read_header(cfg.map_path)
    if mheader != MAP_COLUMNS:
        raise AuditError(f"complexity map columns differ from spec: {mheader}")
    cmap = read_csv_str(cfg.map_path)
    if len(cmap) != cfg.expected_map_rows:
        raise AuditError(f"complexity map rows {len(cmap)} != {cfg.expected_map_rows}")
    cmap["key"] = cmap["raw_name"].map(norm_name)
    if cmap["key"].duplicated().any():
        raise AuditError("complexity map keys are not unique after normalization")

    # ---- Optional ADMISSIONS (only for the pre_admit_time data-quality count)
    adm = None
    if cfg.admissions_path is not None and Path(cfg.admissions_path).exists():
        a = read_csv_str(cfg.admissions_path, usecols=["SUBJECT_ID", "HADM_ID", "ADMITTIME"])
        at, _, _ = parse_time(a["ADMITTIME"])
        adm = pd.DataFrame({"SUBJECT_ID": a["SUBJECT_ID"].map(norm_id), "HADM_ID": a["HADM_ID"].map(norm_id), "ADMITTIME": at})
        adm = adm.drop_duplicates(["SUBJECT_ID", "HADM_ID"])
    info["checks"]["admissions_read"] = adm is not None
    return {"surg": surg, "enc": enc, "map": cmap, "adm": adm, "info": info}


# =============================================================================
# Linkage: one surgery, one encounter (spec v1.4; anaesthesia-end timing v1.7)
# =============================================================================
def link_surgeries(s: pd.DataFrame, enc: pd.DataFrame) -> pd.DataFrame:
    """s: sid, SUBJECT_ID, HADM_ID, T (anaesthesia end; fallback surgery end). enc: encounter_id, SUBJECT_ID,
    hadm_list, t0, OUT, split. Candidates = cohort encounters of the same subject whose hadm_id_list contains the
    surgery HADM_ID. Spec v1.7 (Vinai, Nov 11, 2026):
    Rule 1: t0 - 24 h <= T <= OUT for a candidate -> recent_layer_a if T <= t0 + 6 h (Layer A only), else
            during_icu_layer_c (merge gaps included). The windows of two encounters cannot overlap (>24 h apart).
    Rule 2: else the first candidate with t0 > T (the later encounter) -> remote_same_hadm (T < t0 - 24 h).
    Rule 3: else post_icu_same_hadm (counted, not analysed).
    Returns one row per sid."""
    sidx = s.set_index("sid")
    res = pd.DataFrame(index=sidx.index)
    res["link_class"] = ""
    res["timing"] = ""
    res["encounter_id"] = ""
    res["split"] = ""
    res["t0"] = pd.Series(pd.NaT, index=res.index, dtype="datetime64[ns]")
    res["OUT"] = pd.Series(pd.NaT, index=res.index, dtype="datetime64[ns]")
    res["rule"] = 0
    res["n_candidates"] = 0
    res["other_hadm_side"] = ""

    invalid = sidx["T"].isna()
    res.loc[invalid, "link_class"] = "invalid_time"
    notcoh = (~invalid) & (~sidx["SUBJECT_ID"].isin(set(enc["SUBJECT_ID"])))
    res.loc[notcoh, "link_class"] = "not_in_cohort"

    eh = enc[["encounter_id", "SUBJECT_ID", "hadm_list", "t0", "OUT", "split"]].explode("hadm_list")
    eh = eh.rename(columns={"hadm_list": "HADM_ID"}).dropna(subset=["HADM_ID"])
    eh = eh.drop_duplicates(["encounter_id", "HADM_ID"])
    pending = sidx.loc[(~invalid) & (~notcoh), ["SUBJECT_ID", "HADM_ID", "T"]].reset_index()
    cand = pending.merge(eh, on=["SUBJECT_ID", "HADM_ID"], how="inner")
    ncand = cand.groupby("sid").size()
    res.loc[ncand.index, "n_candidates"] = ncand.to_numpy()

    no_cand = pending.loc[~pending["sid"].isin(ncand.index)]
    if len(no_cand):
        first_t0 = enc.groupby("SUBJECT_ID")["t0"].min()
        ft = no_cand["SUBJECT_ID"].map(first_t0)
        side = np.where(no_cand["T"].to_numpy() < ft.to_numpy(), "earlier_hadm", "later_hadm")
        res.loc[no_cand["sid"].to_numpy(), "link_class"] = "other_hadm"
        res.loc[no_cand["sid"].to_numpy(), "other_hadm_side"] = side

    cand_ok = cand[cand["t0"].notna() & cand["OUT"].notna()]
    unlinkable = sorted(set(ncand.index) - set(cand_ok["sid"]))
    res.loc[unlinkable, "link_class"] = "not_linkable_missing_encounter_time"

    def assign(sub: pd.DataFrame, rule: int) -> list:
        sub = sub.sort_values(["sid", "t0", "encounter_id"]).drop_duplicates("sid", keep="first")
        idx = sub["sid"].to_numpy()
        res.loc[idx, "encounter_id"] = sub["encounter_id"].to_numpy()
        res.loc[idx, "split"] = sub["split"].to_numpy()
        res.loc[idx, "t0"] = sub["t0"].to_numpy()
        res.loc[idx, "OUT"] = sub["OUT"].to_numpy()
        res.loc[idx, "rule"] = rule
        res.loc[idx, "link_class"] = "linked"
        return list(idx)

    lo = cand_ok["t0"] - pd.Timedelta(hours=RECENT_BEFORE_T0_HOURS)
    done1 = assign(cand_ok[(lo <= cand_ok["T"]) & (cand_ok["T"] <= cand_ok["OUT"])], 1)
    if done1:
        after = (sidx.loc[done1, "T"] - res.loc[done1, "t0"]).dt.total_seconds() / 3600.0
        res.loc[done1, "timing"] = np.where(after.to_numpy() <= RECENT_AFTER_T0_HOURS, "recent_layer_a", "during_icu_layer_c")
    rest = cand_ok[~cand_ok["sid"].isin(done1)]
    done2 = assign(rest[rest["t0"] > rest["T"]], 2)
    res.loc[done2, "timing"] = "remote_same_hadm"  # T < t0 - 24 h necessarily (else rule 1 would have matched)
    rest3 = sorted(set(cand_ok["sid"]) - set(done1) - set(done2))
    res.loc[rest3, "link_class"] = "post_icu_same_hadm"
    res.loc[rest3, "rule"] = 3
    if (res["link_class"] == "").any():
        raise ReconciliationError("some surgeries were not assigned a link class")
    return res.reset_index()


# =============================================================================
# Disclosure registry, suppression, recoverability solver
# =============================================================================
class Registry:
    """Every published count is a cell (JSON pointer). Linear identities among cells are verified on RAW values
    (this is the reconciliation, done before suppression) and are later used by complementary suppression and by the
    recoverability solver. Hidden cells are true but unpublished quantities (e.g. 'primary n minus sensitivity n'),
    so differencing between published tables is also checked."""

    def __init__(self) -> None:
        self.raw: dict[str, int] = {}
        self.hidden: set[str] = set()
        self.identities: list[tuple[str, tuple[str, ...]]] = []

    def count(self, cid: str, value: int) -> str:
        if cid in self.raw:
            raise AuditError(f"duplicate cell id {cid}")
        v = int(value)
        if v < 0:
            raise ReconciliationError(f"negative count at {cid}")
        self.raw[cid] = v
        return cid

    def hide(self, cid: str, value: int) -> str:
        cid = "hidden:" + cid
        self.count(cid, value)
        self.hidden.add(cid)
        return cid

    def identity(self, total: str, parts: Iterable[str]) -> None:
        parts = tuple(parts)
        s = sum(self.raw[p] for p in parts)
        if self.raw[total] != s:
            raise ReconciliationError(f"reconciliation identity failed at {total}: {self.raw[total]} != {s} (sum of {len(parts)} cells)")
        self.identities.append((total, parts))

    def subset(self, child: str, parent: str, tag: str) -> None:
        diff = self.raw[parent] - self.raw[child]
        if diff < 0:
            raise ReconciliationError(f"subset larger than parent: {child} > {parent}")
        h = self.hide(f"{child}|complement|{tag}", diff)
        self.identity(parent, [child, h])

    def equal(self, a: str, b: str) -> None:
        self.identity(a, [b])

    def published(self) -> list[str]:
        return [c for c in self.raw if c not in self.hidden]


def _bounds(kind: str) -> tuple[float, float | None]:
    return {"P": (0.0, float(MIN_CELL - 1)), "S": (float(MIN_CELL), None), "H": (0.0, None)}[kind]


def recoverability_check(identities: list[tuple[str, tuple[str, ...]]], known: dict[str, int],
                         unknown: dict[str, str], cache: dict | None = None,
                         propagation_only: bool = False, proven_safe: set | None = None) -> tuple[int, list[str], dict]:
    """Bounds for every unknown cell given ALL published identities (exact propagation of single-unknown identities,
    then LP on each remaining connected component; continuous relaxation, rounded inward to integers).
    Unknown kinds: P = '<10' (0..9), S = complementary (>= 10), H = hidden/unpublished difference (>= 0).
    Unsafe: P exactly determined; S exactly determined; H bounded inside [1, 9] (a provably small hidden count)."""
    rels = []
    for total, parts in identities:
        unk: dict[str, float] = {}
        const = 0.0
        for i, c in [(total, 1.0)] + [(p, -1.0) for p in parts]:
            if i in unknown:
                unk[i] = unk.get(i, 0.0) + c
            else:
                if i not in known:
                    raise ReleaseGateError(f"identity references a cell with no published value: {i}")
                const += c * known[i]
        unk = {i: c for i, c in unk.items() if c != 0}
        if not unk:
            if abs(const) > 1e-9:
                raise ReleaseGateError(f"published cells violate an identity at {total}")
            continue
        rels.append([unk, -const])
    # 1) exact propagation: an identity with a single unknown determines it
    fixed: dict[str, float] = {}
    by_var: dict[str, list[int]] = {}
    for k, (unk, _) in enumerate(rels):
        for v in unk:
            by_var.setdefault(v, []).append(k)
    queue = [k for k, (unk, _) in enumerate(rels) if len(unk) == 1]
    while queue:
        k = queue.pop()
        unk, b = rels[k]
        if len(unk) != 1:
            continue
        (v, c), = unk.items()
        val = b / c
        if v in fixed:
            if abs(fixed[v] - val) > 1e-6:
                raise ReleaseGateError("published cells are inconsistent (propagation)")
        else:
            fixed[v] = val
        for k2 in by_var[v]:
            u2, b2 = rels[k2]
            if v in u2:
                rels[k2][1] = b2 - u2.pop(v) * val
                if len(u2) == 1:
                    queue.append(k2)
    bounds_out: dict[str, tuple[float, float]] = {v: (x, x) for v, x in fixed.items()}
    live = [(u, b) for u, b in rels if u]
    # 2) components of the remaining unknowns
    parent: dict[str, str] = {}

    def find(a: str) -> str:
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for u, _ in live:
        vs = list(u)
        for v in vs[1:]:
            ra, rb = find(vs[0]), find(v)
            if ra != rb:
                parent[rb] = ra
    comps: dict[str, list] = {}
    for u, b in live:
        comps.setdefault(find(next(iter(u))), []).append((u, b))
    # proven_safe (used only inside apply_suppression): suppression only ever turns published cells into S, which
    # relaxes the system, so a cell once proven safe stays safe. The release gate never passes it (full check).
    for crels in ([] if propagation_only else comps.values()):
        vars_ = sorted({v for u, _ in crels for v in u})
        if proven_safe is not None and all(v in proven_safe for v in vars_):
            continue
        key = (tuple((v, unknown[v]) for v in vars_), tuple(sorted((tuple(sorted(u.items())), round(b, 6)) for u, b in crels)))
        if cache is not None and key in cache:
            bounds_out.update(cache[key])
            continue
        col = {v: j for j, v in enumerate(vars_)}
        rows, cols, data, rhs = [], [], [], []
        for r, (u, b) in enumerate(crels):
            for v, c in u.items():
                rows.append(r); cols.append(col[v]); data.append(c)
            rhs.append(b)
        A = sparse.csr_matrix((data, (rows, cols)), shape=(len(crels), len(vars_)))
        bnds = [_bounds(unknown[v]) for v in vars_]
        res: dict[str, tuple[float, float]] = {}
        # Screening: feasible points from a few LPs (random objectives; open bounds capped so every LP is bounded).
        # A cell whose observed values already span two integers (or, for H, reach 0 or 10) is provably safe, so
        # only the remaining cells need exact min/max LPs. Every screening point is feasible for the uncapped system.
        nv = len(vars_)
        cap = [(lo_, hi_ if hi_ is not None else SCREEN_CAP) for lo_, hi_ in bnds]
        obs_lo, obs_hi = np.full(nv, np.inf), np.full(nv, -np.inf)
        rng = np.random.default_rng(nv)
        need = [j for j, v in enumerate(vars_) if proven_safe is None or v not in proven_safe]
        for t in range(SCREEN_LPS if len(need) > 2 * SCREEN_LPS else 0):
            cvec = rng.normal(size=nv) * (1 if t % 2 else -1)
            sol = linprog(cvec, A_eq=A, b_eq=rhs, bounds=cap, method="highs")
            if sol.status == 2:
                raise ReleaseGateError("recoverability LP infeasible; published cells inconsistent")
            if sol.status == 0:
                obs_lo, obs_hi = np.minimum(obs_lo, sol.x), np.maximum(obs_hi, sol.x)
        todo = []
        for j in need:
            v = vars_[j]
            a, b = obs_lo[j], obs_hi[j]
            if not np.isfinite(a):
                todo.append(j); continue
            if unknown[v] == "H":
                safe = b >= MIN_CELL - 1e-6 or a <= 1e-6
            else:
                safe = math.ceil(a - 1e-6) < math.floor(b + 1e-6)
            if safe:
                res[v] = (float(a), float(b))  # observed feasible range (sub-interval of the true bounds)
            else:
                todo.append(j)
        for j in todo:
            v = vars_[j]
            cvec = np.zeros(len(vars_)); cvec[j] = 1.0
            hi = linprog(-cvec, A_eq=A, b_eq=rhs, bounds=bnds, method="highs")
            if hi.status not in (0, 3):
                raise ReleaseGateError(f"recoverability LP infeasible (status {hi.status}); published cells inconsistent")
            hv = -hi.fun if hi.status == 0 else math.inf
            if unknown[v] == "H" and hv >= MIN_CELL - 1e-6:
                res[v] = (0.0, hv)  # cannot be bounded inside [1, 9]; min not needed
                continue
            lo = linprog(cvec, A_eq=A, b_eq=rhs, bounds=bnds, method="highs")
            if lo.status != 0:
                raise ReleaseGateError(f"recoverability LP infeasible (status {lo.status}); published cells inconsistent")
            res[v] = (lo.fun, hv)
        if cache is not None:
            cache[key] = res
        bounds_out.update(res)
        if proven_safe is not None:
            proven_safe.update(v for v, (lo_, hi_) in res.items() if not _unsafe_kind(unknown[v], lo_, hi_))
    unsafe = []
    for u, kind in unknown.items():
        dlo, dhi = _bounds(kind)
        lo, hi = bounds_out.get(u, (dlo, math.inf if dhi is None else dhi))
        if _unsafe_kind(kind, lo, hi):
            unsafe.append(u)
    return len(unknown), unsafe, bounds_out


def _unsafe_kind(kind: str, lo: float, hi: float) -> bool:
    """P / S: unsafe if the integer interval is a single value. H: unsafe if it lies inside [1, 9]."""
    li = math.ceil(lo - 1e-6)
    hi_i = math.floor(hi + 1e-6) if math.isfinite(hi) else math.inf
    if kind in ("P", "S"):
        return li == hi_i
    return li >= 1 and hi_i <= MIN_CELL - 1


def _neighbour_candidate(reg: "Registry", u: str, status: dict, by_cell: dict) -> tuple | None:
    """Greedy order elsewhere: fewest unknowns in the identity, then the published member that sits in the fewest
    identities (limits cascades), then the smallest value.
    When none of the identities of an unsafe cell has a published member, walk outwards through identities linked
    by unknown cells (breadth first) and return the smallest published cell at the nearest distance."""
    seen_cells, seen_ids = {u}, set()
    frontier = [u]
    while frontier:
        cands = []
        nxt = []
        for c in frontier:
            for k in by_cell.get(c, []):
                if k in seen_ids:
                    continue
                seen_ids.add(k)
                t, parts = reg.identities[k]
                for m in (t,) + parts:
                    if status.get(m) == "pub":
                        cands.append((reg.raw[m], m, k))
                    elif m not in seen_cells:
                        seen_cells.add(m)
                        nxt.append(m)
        if cands:
            v, m, k = min(cands)
            return ((0, v, m), k, m)
        frontier = nxt
    return None


def apply_suppression(reg: Registry, max_rounds: int = 1000) -> tuple[dict[str, str], dict]:
    """Primary: counts 0-9 -> '<10' (spec v1.5: 'counts of 1-9 and zero'). Complementary 'S' is added
    algorithmically, round by round, until the recoverability solver finds no recoverable cell."""
    status = {c: ("P" if reg.raw[c] < MIN_CELL else "pub") for c in reg.published()}
    by_cell: dict[str, list[int]] = {}
    for k, (t, parts) in enumerate(reg.identities):
        for c in (t,) + parts:
            by_cell.setdefault(c, []).append(k)
    cache: dict = {}
    proven: set[str] = set()
    fast = True  # cheap exact-propagation rounds first; the full LP round is the one that must come back clean
    for _ in range(max_rounds):
        known = {c: reg.raw[c] for c, st in status.items() if st == "pub"}
        unknown = {c: st for c, st in status.items() if st != "pub"} | {h: "H" for h in reg.hidden}
        checked, unsafe, _ = recoverability_check(reg.identities, known, unknown, cache, propagation_only=fast,
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
                continue  # an identity of this cell was already opened this round; re-check next round
            best = None
            for k in by_cell.get(u, []):
                t, parts = reg.identities[k]
                members = (t,) + parts
                pubs = [m for m in members if status.get(m) == "pub"]
                if not pubs:
                    continue
                n_unknown = sum(1 for m in members if m not in status or status[m] != "pub")
                cand = min(pubs, key=lambda m: (len(by_cell.get(m, [])), reg.raw[m], m))
                score = (n_unknown, len(by_cell.get(cand, [])), reg.raw[cand], cand)
                if best is None or score < best[0]:
                    best = (score, k, cand)
            if best is None:
                best = _neighbour_candidate(reg, u, status, by_cell)
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


# =============================================================================
# Analysis
# =============================================================================
def qstats(vals: np.ndarray) -> dict:
    """Quantiles only for n >= 10 (raw); never min/max."""
    if len(vals) < MIN_CELL:
        return {k: "<10" for k in ("median", "p25", "p75", "p10", "p90")}
    p10, p25, p50, p75, p90 = np.percentile(vals, [10, 25, 50, 75, 90])
    return {"median": r1(p50), "p25": r1(p25), "p75": r1(p75), "p10": r1(p10), "p90": r1(p90)}


def build_frame(inp: dict) -> tuple[pd.DataFrame, dict]:
    surg, enc, cmap, adm = inp["surg"], inp["enc"], inp["map"], inp["adm"]
    s = pd.DataFrame({"sid": np.arange(len(surg))})
    s["ROW_ID"] = surg[SURGERY_KEY].map(norm_id).to_numpy()
    if (s["ROW_ID"] == "").any() or s["ROW_ID"].duplicated().any():
        raise ReconciliationError(f"surgery key {SURGERY_KEY} blank or not unique in SURGERY_INFO")
    s["SUBJECT_ID"] = surg["SUBJECT_ID"].map(norm_id).to_numpy()
    s["HADM_ID"] = surg["HADM_ID"].map(norm_id).to_numpy()
    s["OPER_ID"] = surg["OPER_ID"].map(norm_id).to_numpy()  # descriptor only (v1.6): not unique within an admission
    B, _, _ = parse_time(surg["SURGERY_BEGIN_TIME"])
    E, _, _ = parse_time(surg["SURGERY_END_TIME"])
    AS, _, _ = parse_time(surg["ANES_START_TIME"])
    AE, _, _ = parse_time(surg["ANES_END_TIME"])
    s["B"], s["E"] = B.to_numpy(), E.to_numpy()
    dur = compute_durations(B, E, AS, AE)
    s["D"], s["invalid_reason"], s["fallback"] = dur["D"].to_numpy(), dur["invalid_reason"].to_numpy(), dur["fallback"].to_numpy()
    s["valid"] = s["D"].notna()
    s["gt24"] = s["valid"] & (s["D"] > D_GT_HOURS)
    s["name"] = surg["SURGERY_NAME"].map(norm_name).to_numpy()
    # v1.7 timing time: anaesthesia end, falling back to surgery end (counted)
    s["AE"] = AE.to_numpy()
    s["T"] = s["AE"].where(s["AE"].notna(), s["E"])
    s["timing_fallback"] = s["AE"].isna() & s["E"].notna()
    key_diag = surgery_key_diagnostics(s)

    link = link_surgeries(s[["sid", "SUBJECT_ID", "HADM_ID", "T"]], enc)
    s = s.merge(link, on="sid", how="left", validate="one_to_one")
    if len(s) != len(surg):
        raise ReconciliationError("linkage changed the number of surgery rows")

    m = cmap.set_index("key")
    s["mapped"] = s["name"].isin(m.index)

    def get(col: str) -> pd.Series:
        return s["name"].map(m[col]).fillna("").astype(str)
    s["body_system"] = np.where(s["mapped"], get("body_system"), UNMAPPED)
    s["intensity"], s["STAT"] = get("intensity"), get("STAT")
    s["review_flag"], s["crosses"] = get("review_flag"), get("STAT_group_crosses")
    s["n_components"] = pd.to_numeric(get("n_components"), errors="coerce")
    s["classified"] = s["mapped"] & (s["body_system"] != UNCLASSIFIED)
    cg, eg = [], []
    for bs, st, it, mg, cl in zip(s["body_system"], s["STAT"], s["intensity"], get("complexity_group"), s["classified"]):
        if not cl:
            cg.append(""); eg.append(""); continue
        cg.append(primary_group(bs, st, it, mg))
        eg.append(st if bs == CARDIAC else it)
    s["cgroup"], s["egroup"] = cg, eg
    s["straddle"] = ((s["link_class"] == "linked") & (s["B"] < s["t0"]) & s["E"].notna() & (s["E"] >= s["t0"])).fillna(False)
    s["recent_after_t0"] = ((s["timing"] == "recent_layer_a") & (s["T"] > s["t0"])).fillna(False)
    s["timing_fallback_linked"] = (s["link_class"] == "linked") & s["timing_fallback"]
    s["pre_admit"] = False
    if adm is not None:  # remote surgeries ending before ADMITTIME are kept and flagged
        at = s[["SUBJECT_ID", "HADM_ID"]].merge(adm, on=["SUBJECT_ID", "HADM_ID"], how="left")["ADMITTIME"]
        s["pre_admit"] = ((s["timing"] == "remote_same_hadm").to_numpy() & (s["T"].to_numpy() < at.to_numpy()))
    s["analysed"] = (s["link_class"] == "linked") & s["classified"]
    lk = s[s["link_class"] == "linked"]
    if lk.duplicated(["ROW_ID"]).any() or len(lk) > len(surg):
        raise ReconciliationError("a surgery key appears twice in the analysis set")
    systems = sorted(b for b in cmap["body_system"].unique() if b != UNCLASSIFIED)
    groups: dict[str, dict] = {}
    for bs in systems:
        sub = cmap[cmap["body_system"] == bs]
        if bs == CARDIAC:
            stats_present = set(sub["STAT"])
            present_p = {STAT_TO_G[x] for x in stats_present}
            groups[bs] = {"primary": [g for g in G_ORDER + ["NC"] if g in present_p],
                          "primary_ordinal": [g for g in G_ORDER if g in present_p],
                          "e": [x for x in STAT_ORDER + ["NC"] if x in stats_present],
                          "e_ordinal": [x for x in STAT_ORDER if x in stats_present]}
        else:
            present_i = set(sub["intensity"])
            prim = [g for g in I_ORDER if g in {INTENSITY_TO_I[i] for i in present_i}]
            e = [i for i in INTENSITY_ORDER if i in present_i]
            groups[bs] = {"primary": prim, "primary_ordinal": prim, "e": e, "e_ordinal": e}
    return s, {"systems": systems, "groups": groups, "key_diag": key_diag}


def surgery_key_diagnostics(s: pd.DataFrame) -> dict:
    """v1.6 key diagnostics (reported, nothing dropped): rows sharing (SUBJECT_ID, HADM_ID, OPER_ID), and possible
    true duplicates = rows sharing (SUBJECT_ID, HADM_ID, SURGERY_BEGIN_TIME, normalized SURGERY_NAME), begin time
    present. Raw sizes stay in the LOCAL file; the shared file gets suppressed counts only."""
    g = s.groupby(["SUBJECT_ID", "HADM_ID", "OPER_ID"], sort=False).size()
    shared = g[g >= 2]
    d = s[s["B"].notna()]
    g2 = d.groupby(["SUBJECT_ID", "HADM_ID", "B", "name"], sort=False).size()
    dup = g2[g2 >= 2]
    return {"shared_oper_id_groups": int(len(shared)), "rows_in_shared_oper_id_groups": int(shared.sum()),
            "shared_oper_id_group_sizes": {str(k): int(v) for k, v in shared.value_counts().sort_index().items()},
            "possible_true_duplicate_groups": int(len(dup)), "rows_in_possible_true_duplicate_groups": int(dup.sum())}


FAMILIES = {
    # key: (json section, row mask, splits, bands)
    "primary": (None, lambda d: np.ones(len(d), bool), ALL_SPLITS, True),
    "e_four_level_or_stat5": ("sensitivity", lambda d: np.ones(len(d), bool), ALL_SPLITS, True),
    "a_single_component": ("sensitivity", lambda d: (d["n_components"] == 1).to_numpy(), ("train", "pooled"), False),
    "b_no_cross_boundary": ("sensitivity", lambda d: (d["crosses"] != "Y").to_numpy(), ("train", "pooled"), False),
    "c_no_gt24_no_fallback": ("sensitivity", lambda d: (~d["gt24"] & ~d["fallback"]).to_numpy(), ("train", "pooled"), False),
    "d_during_icu": ("sensitivity", lambda d: (d["timing"] == "during_icu_layer_c").to_numpy(), ("train", "pooled"), False),
    "review_flag_Y": ("flagged_vs_unflagged", lambda d: (d["review_flag"] == "Y").to_numpy(), ("train", "pooled"), False),
    "review_flag_N": ("flagged_vs_unflagged", lambda d: (d["review_flag"] == "N").to_numpy(), ("train", "pooled"), False),
    "STAT_group_crosses_Y": ("flagged_vs_unflagged", lambda d: (d["crosses"] == "Y").to_numpy(), ("train", "pooled"), False),
    "STAT_group_crosses_N": ("flagged_vs_unflagged", lambda d: (d["crosses"] == "N").to_numpy(), ("train", "pooled"), False),
}
CARDIAC_ONLY_FAMILIES = ("b_no_cross_boundary", "STAT_group_crosses_Y", "STAT_group_crosses_N")


def fam_base(bs: str, fam: str) -> tuple:
    sec = FAMILIES[fam][0]
    return ("systems", bs) if sec is None else ("systems", bs, sec, fam)


def split_rows(df: pd.DataFrame, sp: str) -> pd.DataFrame:
    return df if sp == "pooled" else df[df["split"] == sp]


class Analysis:
    def __init__(self, cfg: AuditConfig, frame: pd.DataFrame, meta: dict, info: dict):
        self.cfg, self.f, self.meta, self.info = cfg, frame, meta, info
        self.reg = Registry()
        self.stats: dict[str, Any] = {}          # pointer -> non-count value (never suppressed as a count)
        self.local: dict[str, Any] = {}          # raw extras for the LOCAL file only
        self.boot = Bootstrapper(cfg.bootstrap_reps, cfg.seed)
        self.post: list[Callable[[dict, dict], None]] = []  # publishers that depend on the suppression status
        self.verdicts: dict[str, dict] = {}
        self.done: dict[str, set] = {}

    def st(self, ptr: str, value: Any) -> None:
        self.stats[ptr] = value

    def _xsplit(self, ids_for: Callable[[str], list[str]]) -> None:
        per = {sp: ids_for(sp) for sp in ALL_SPLITS}
        for j, pooled_id in enumerate(per["pooled"]):
            self.reg.identity(pooled_id, [per[sp][j] for sp in SPLITS])

    # ---------------------------------------------------------------- population & reconciliation
    def population(self) -> None:
        f, reg = self.f, self.reg
        systems = self.meta["systems"]
        total = reg.count(P("reconciliation", "surgery_info_rows"), len(f))
        reg.identity(total, [reg.count(P("reconciliation", c), int((f["link_class"] == c).sum())) for c in LINK_CLASSES])
        oh = f[f["link_class"] == "other_hadm"]
        reg.identity(P("reconciliation", "other_hadm"),
                     [reg.count(P("reconciliation", "other_hadm_by_side", k), int((oh["other_hadm_side"] == k).sum()))
                      for k in ("earlier_hadm", "later_hadm")])
        multi = f[f["n_candidates"] >= 2]
        nm = reg.count(P("reconciliation", "n_multi_candidate"), len(multi))
        rs = [reg.count(P("reconciliation", k), int((multi["rule"] == r).sum()))
              for k, r in (("resolved_rule1_during", 1), ("resolved_rule2_next_t0", 2), ("resolved_rule3_post_icu", 3))]
        reg.identity(nm, rs + [reg.hide("reconciliation/multi_candidate_unresolved", len(multi) - sum(reg.raw[x] for x in rs))])
        reg.subset(nm, total, "rows")
        kd = self.meta.get("key_diag", {})
        for k in ("shared_oper_id_groups", "rows_in_shared_oper_id_groups", "possible_true_duplicate_groups",
                  "rows_in_possible_true_duplicate_groups"):
            reg.count(P("reconciliation", "surgery_key_ROW_ID", k), int(kd.get(k, 0)))
        self.local["surgery_key_diagnostics"] = kd

        lk = f[f["link_class"] == "linked"]
        for sp in ALL_SPLITS:
            d = split_rows(lk, sp)
            base = ("population", sp)
            L = reg.count(P(*base, "linked"), len(d))
            tc = [reg.count(P(*base, "by_timing_class", c), int((d["timing"] == c).sum())) for c in TIMING]
            reg.identity(L, tc)
            an = d[d["analysed"]]
            A = reg.count(P(*base, "analysed"), len(an))
            U = reg.count(P(*base, "unmapped"), int((~d["mapped"]).sum()))
            X = reg.count(P(*base, "excluded_unclassified"), int((d["mapped"] & ~d["classified"]).sum()))
            reg.identity(L, [A, U, X])
            V = reg.count(P(*base, "valid_duration"), int(an["valid"].sum()))
            I = reg.count(P(*base, "invalid_duration", "invalid_total"), int((~an["valid"]).sum()))
            reg.identity(A, [V, I])
            reg.identity(I, [reg.count(P(*base, "invalid_duration", r), int((an["invalid_reason"] == r).sum())) for r in INVALID_REASONS])
            for key, mask, parent in (("duration_fallback_anes", an["fallback"], V), ("d_gt_24h", an["gt24"], V),
                                      ("straddling_t0", d["straddle"], L), ("recent_ending_after_t0", d["recent_after_t0"], tc[0]),
                                      ("timing_fallback_surgery_end", d["timing_fallback_linked"], L),
                                      ("pre_admit_time", d["pre_admit"], tc[1])):
                reg.subset(reg.count(P(*base, key), int(mask.sum())), parent, "pop")
            sysL, sysV, sysI = [], [], []
            sysT: dict[str, list] = {c: [] for c in TIMING}
            for bs in systems:
                ds = an[an["body_system"] == bs]
                l_ = reg.count(P(*base, "by_system", bs, "linked"), len(ds))
                v_ = reg.count(P(*base, "by_system", bs, "valid"), int(ds["valid"].sum()))
                i_ = reg.count(P(*base, "by_system", bs, "invalid"), int((~ds["valid"]).sum()))
                reg.identity(l_, [v_, i_])
                ts = [reg.count(P(*base, "by_system", bs, "by_timing_class", c), int((ds["timing"] == c).sum())) for c in TIMING]
                reg.identity(l_, ts)
                for c, t in zip(TIMING, ts):
                    sysT[c].append(t)
                sysL.append(l_); sysV.append(v_); sysI.append(i_)
            reg.identity(A, sysL); reg.identity(V, sysV); reg.identity(I, sysI)
            for c, t in zip(TIMING, tc):
                excl = int(((d["timing"] == c) & ~d["analysed"]).sum())
                reg.identity(t, sysT[c] + [reg.hide(f"population/{sp}/timing_not_analysed/{c}", excl)])
        self._xsplit(lambda sp: [P("population", sp, k) for k in ("linked", "analysed", "unmapped", "excluded_unclassified", "valid_duration",
                                                                       "duration_fallback_anes", "d_gt_24h", "straddling_t0", "recent_ending_after_t0",
                                                                       "timing_fallback_surgery_end", "pre_admit_time")]
                     + [P("population", sp, "by_timing_class", c) for c in TIMING]
                     + [P("population", sp, "invalid_duration", r) for r in INVALID_REASONS + ("invalid_total",)]
                     + [P("population", sp, "by_system", bs, k) for bs in systems for k in ("linked", "valid", "invalid")]
                     + [P("population", sp, "by_system", bs, "by_timing_class", c) for bs in systems for c in TIMING])
        reg.equal(P("reconciliation", "linked"), P("population", "pooled", "linked"))

    # ---------------------------------------------------------------- families (group stats, bands, margins)
    def family(self, bs: str, fam: str) -> None:
        _, maskf, splits, with_bands = FAMILIES[fam]
        reg, g = self.reg, self.meta["groups"][bs]
        is_e = fam == "e_four_level_or_stat5"
        col = "egroup" if is_e else "cgroup"
        groups = g["e"] if is_e else g["primary"]
        base = fam_base(bs, fam)
        d_all = self.f[self.f["analysed"] & (self.f["body_system"] == bs)]
        d = d_all[maskf(d_all)] if len(d_all) else d_all
        d = d[d["valid"]]
        for grp in groups:
            for sp in splits:
                rows = split_rows(d, sp)
                vals = rows.loc[rows[col] == grp, "D"].to_numpy(float)
                n = reg.count(P(*base, "groups", grp, sp, "n"), len(vals))
                if with_bands:
                    bl = [band_of(x) for x in vals]
                    reg.identity(n, [reg.count(P(*base, "groups", grp, sp, "bands", b), bl.count(b)) for b in BANDS])
                for k, v in qstats(vals).items():
                    self.st(P(*base, "groups", grp, sp, k), v)
            if splits == ALL_SPLITS:
                reg.identity(P(*base, "groups", grp, "pooled", "n"), [P(*base, "groups", grp, sp, "n") for sp in SPLITS])
                if with_bands:
                    for b in BANDS:
                        reg.identity(P(*base, "groups", grp, "pooled", "bands", b), [P(*base, "groups", grp, sp, "bands", b) for sp in SPLITS])
            else:
                reg.subset(P(*base, "groups", grp, "train", "n"), P(*base, "groups", grp, "pooled", "n"), "valtest")
        if fam == "primary":
            for sp in ALL_SPLITS:
                tot = reg.count(P("systems", bs, "crosstab_group_band", sp, "total"), int(len(split_rows(d, sp))))
                reg.identity(tot, [P(*base, "groups", grp, sp, "n") for grp in groups])
                reg.equal(tot, P("population", sp, "by_system", bs, "valid"))
                bt = []
                for b in BANDS:
                    ids = [P(*base, "groups", grp, sp, "bands", b) for grp in groups]
                    c = reg.count(P("systems", bs, "crosstab_group_band", sp, "band_totals", b), sum(reg.raw[i] for i in ids))
                    reg.identity(c, ids)
                    bt.append(c)
                reg.identity(tot, bt)
            for b in BANDS:
                reg.identity(P("systems", bs, "crosstab_group_band", "pooled", "band_totals", b),
                             [P("systems", bs, "crosstab_group_band", sp, "band_totals", b) for sp in SPLITS])
            reg.identity(P("systems", bs, "crosstab_group_band", "pooled", "total"),
                         [P("systems", bs, "crosstab_group_band", sp, "total") for sp in SPLITS])
        elif is_e:
            to_p = (lambda x: STAT_TO_G[x]) if bs == CARDIAC else (lambda x: INTENSITY_TO_I[x])
            for sp in ALL_SPLITS:
                for pg in g["primary"]:
                    members = [e for e in groups if to_p(e) == pg]
                    reg.identity(P("systems", bs, "groups", pg, sp, "n"), [P(*base, "groups", e, sp, "n") for e in members])
                    for b in BANDS:
                        reg.identity(P("systems", bs, "groups", pg, sp, "bands", b), [P(*base, "groups", e, sp, "bands", b) for e in members])
        else:
            for grp in groups:
                for sp in splits:
                    reg.subset(P(*base, "groups", grp, sp, "n"), P("systems", bs, "groups", grp, sp, "n"), fam)
            pair = {"review_flag_N": "review_flag_Y", "STAT_group_crosses_N": "STAT_group_crosses_Y"}.get(fam)
            if pair and pair in self.done.get(bs, set()):
                for grp in groups:
                    for sp in splits:
                        tot = P("systems", bs, "groups", grp, sp, "n")
                        a, b_ = P(*base, "groups", grp, sp, "n"), P(*fam_base(bs, pair), "groups", grp, sp, "n")
                        h = reg.hide(f"{bs}/{fam}/{grp}/{sp}/neither", reg.raw[tot] - reg.raw[a] - reg.raw[b_])
                        reg.identity(tot, [a, b_, h])
            if fam == "STAT_group_crosses_N" and "b_no_cross_boundary" in self.done.get(bs, set()):
                # b (cardiac names with crosses = Y excluded) contains crosses = N plus any blank flag.
                for grp in groups:
                    for sp in splits:
                        b_id = P(*fam_base(bs, "b_no_cross_boundary"), "groups", grp, sp, "n")
                        n_id = P(*base, "groups", grp, sp, "n")
                        reg.identity(b_id, [n_id, reg.hide(f"{bs}/b_minus_crossesN/{grp}/{sp}", reg.raw[b_id] - reg.raw[n_id])])
        self.done.setdefault(bs, set()).add(fam)
        self.comparisons(bs, fam, d, col, g["e_ordinal"] if is_e else g["primary_ordinal"], base)

    # ---------------------------------------------------------------- correlation, overlap, verdict inputs
    def comparisons(self, bs: str, fam: str, d: pd.DataFrame, col: str, ordinal: list, base: tuple) -> None:
        reg = self.reg
        tr = d[d["split"] == "train"]
        elig = [grp for grp in ordinal if int((tr[col] == grp).sum()) >= MIN_CELL]
        assessable = len(elig) >= 2
        reason = None if assessable else ("single_group" if len(ordinal) <= 1 else "groups_below_10")
        res: dict[str, Any] = {"eligible_groups_train": elig, "assessable": assessable, "reason": reason,
                               "tau_b": float("nan"), "spearman": float("nan"), "tau_ci": "insufficient"}
        code = {grp: i for i, grp in enumerate(ordinal)}
        corr_splits = ALL_SPLITS if fam == "primary" else ("train",)
        for sp in corr_splits:
            ds = split_rows(d, sp)
            inc = [grp for grp in ordinal if int((ds[col] == grp).sum()) >= MIN_CELL]
            sub = ds[ds[col].isin(inc)]
            ptr = P(*base, "correlation", sp)
            if len(inc) < 2 or len(sub) < CORR_MIN_TOTAL:
                self.st(ptr, {"kendall_tau_b": "insufficient", "spearman": "insufficient", "groups_included": inc,
                              "tau_b_from_crosstab": "insufficient"})
                continue
            x = sub[col].map(code).to_numpy(float)
            y = sub["D"].to_numpy(float)
            subj = sub["SUBJECT_ID"].to_numpy()
            tau = kendall_tau_b_surgery_level(x, y)
            rho = spearman_rho(x, y)
            key = f"{bs}|{fam}|{sp}"
            tci = self.boot.ci(subj, lambda idx, x=x, y=y: kendall_tau_b_surgery_level(x[idx], y[idx]), key + "|tau")
            rci = self.boot.ci(subj, lambda idx, x=x, y=y: spearman_rho(x[idx], y[idx]), key + "|rho")
            self.st(ptr, {"kendall_tau_b": {"est": r3(tau) if np.isfinite(tau) else "insufficient", "ci": tci["ci"],
                                            "ci_includes_0": bool(isinstance(tci["ci"], list) and tci["ci"][0] <= 0 <= tci["ci"][1])},
                          "spearman": {"est": r3(rho) if np.isfinite(rho) else "insufficient", "ci": rci["ci"]},
                          "abs_tau_b": r3(abs(tau)) if np.isfinite(tau) else "insufficient",
                          "groups_included": inc, "unit": "surgery-level; SUBJECT_ID-cluster bootstrap"})
            if FAMILIES[fam][3]:
                cell_ids = [P(*base, "groups", grp, sp, "bands", b) for grp in inc for b in BANDS] + [P(*base, "groups", grp, sp, "n") for grp in inc]
                bidx = np.array([BANDS.index(band_of(v)) for v in y], float)

                def pub_ct(status: dict, stats: dict, ptr=ptr, cell_ids=cell_ids, x=x, bidx=bidx) -> None:
                    if all(status.get(c) == "pub" for c in cell_ids):
                        t = kendall_tau_b_surgery_level(x, bidx)
                        stats[ptr] = dict(stats[ptr], tau_b_from_crosstab=r3(t) if np.isfinite(t) else "insufficient")
                    else:
                        stats[ptr] = dict(stats[ptr], tau_b_from_crosstab="not_published_suppressed_cells")
                self.post.append(pub_ct)
            else:
                self.stats[ptr]["tau_b_from_crosstab"] = "not_applicable_no_bands"
            if sp == "train":
                res["tau_b"], res["spearman"], res["tau_ci"] = tau, rho, tci["ci"]
        pairs = []
        if assessable:
            pairs.append((elig[0], elig[-1]))
            if bs == CARDIAC and fam == "primary":
                for a, b in (("G1", "G2"), ("G2", "G3")):
                    if a in elig and b in elig and (a, b) not in pairs:
                        pairs.append((a, b))
        med_tr = {grp: float(np.median(tr.loc[tr[col] == grp, "D"])) for grp in elig}
        for lo, hi in pairs:
            pk = f"{lo}_vs_{hi}"
            for sp in ("train", "pooled"):
                ds = split_rows(d, sp)
                lv = ds.loc[ds[col] == lo, "D"].to_numpy(float)
                hv = ds.loc[ds[col] == hi, "D"].to_numpy(float)
                obase = (*base, "overlap", pk, sp)
                for metric, arr, thr, cmp, grp in (("p_low_above_high_median", lv, med_tr[hi], np.greater, lo),
                                                   ("p_high_below_low_median", hv, med_tr[lo], np.less, hi)):
                    num = int(cmp(arr, thr).sum())
                    c1 = reg.count(P(*obase, metric, "count"), num)
                    c2 = reg.count(P(*obase, metric, "complement"), len(arr) - num)
                    reg.identity(P(*base, "groups", grp, sp, "n"), [c1, c2])
                    share = num / len(arr) if len(arr) else float("nan")
                    self.local[P(*obase, metric, "share_raw")] = share
                    if sp == "train" and (lo, hi) == pairs[0]:
                        res["d3_share" if metric.startswith("p_low") else "d4_share"] = share

                    def pub_share(status: dict, stats: dict, p=P(*obase, metric, "share"), c1=c1, c2=c2, share=share) -> None:
                        s1, s2 = status.get(c1), status.get(c2)
                        stats[p] = r3(share) if (s1 == "pub" and s2 == "pub") else "suppressed"
                    self.post.append(pub_share)
                if len(lv) >= MIN_CELL and len(hv) >= MIN_CELL:
                    allv = np.concatenate([lv, hv])
                    is_hi = np.r_[np.zeros(len(lv), bool), np.ones(len(hv), bool)]
                    subj = np.concatenate([ds.loc[ds[col] == lo, "SUBJECT_ID"].to_numpy(), ds.loc[ds[col] == hi, "SUBJECT_ID"].to_numpy()])
                    ci = self.boot.ci(subj, lambda idx, allv=allv, is_hi=is_hi: prob_superiority(allv[idx][is_hi[idx]], allv[idx][~is_hi[idx]]),
                                      f"{bs}|{fam}|{pk}|{sp}|ps")
                    self.st(P(*obase, "prob_superiority"), {"est": r3(prob_superiority(hv, lv)), "ci": ci["ci"]})
                else:
                    self.st(P(*obase, "prob_superiority"), "<10")
            self.st(P(*base, "overlap", pk, "thresholds_train_medians"), {lo: r1(med_tr[lo]), hi: r1(med_tr[hi])})
        res["ordered_medians_train"] = [med_tr[grp] for grp in elig]
        res["pair0"] = f"{pairs[0][0]}_vs_{pairs[0][1]}" if pairs else None
        res["base"] = base
        self.verdicts[f"{bs}|{fam}"] = res

    # ---------------------------------------------------------------- names (discordance), D5
    def names(self, bs: str) -> float | None:
        reg, g = self.reg, self.meta["groups"][bs]
        d = self.f[self.f["analysed"] & (self.f["body_system"] == bs) & self.f["valid"]]
        tr = d[d["split"] == "train"]
        gstat = {}
        for grp in g["primary"]:
            v = tr.loc[tr["cgroup"] == grp, "D"].to_numpy(float)
            gstat[grp] = np.percentile(v, [25, 50, 75]) if len(v) >= MIN_CELL else None
        counts = d.groupby("name").size()
        listed = sorted(counts[counts >= MIN_CELL].index)  # pooled analysis n >= 10 only
        n_det = n_disc = 0
        per_group_ids: dict[tuple, list] = {(grp, sp): [] for grp in g["primary"] for sp in ALL_SPLITS}
        ordn = g["primary_ordinal"]
        for nm in listed:
            dn = d[d["name"] == nm]
            grp = dn["cgroup"].iloc[0]
            base = ("systems", bs, "discordant_names", nm)
            for sp in ALL_SPLITS:
                per_group_ids[(grp, sp)].append(reg.count(P(*base, sp, "n"), len(split_rows(dn, sp))))
            reg.identity(P(*base, "pooled", "n"), [P(*base, sp, "n") for sp in SPLITS])
            med = float(np.median(dn["D"]))
            trn = dn.loc[dn["split"] == "train", "D"].to_numpy(float)
            entry: dict[str, Any] = {"group": grp, "stat_or_intensity": dn["egroup"].iloc[0], "median_pooled": r1(med),
                                     "median_train": r1(np.median(trn)) if len(trn) >= MIN_CELL else "<10"}
            q = gstat.get(grp)
            if q is None:
                entry.update({"discordant": "not_assessable", "direction": None, "crosses_adjacent_group_median": None})
            else:
                disc = bool(med < q[0] or med > q[2])
                direction = None if not disc else ("shorter" if med < q[0] else "longer")
                cross = None
                if disc and grp in ordn:
                    i = ordn.index(grp)
                    adj = None
                    if direction == "longer" and i + 1 < len(ordn):
                        adj = ordn[i + 1]
                    elif direction == "shorter" and i > 0:
                        adj = ordn[i - 1]
                    if adj is not None and gstat.get(adj) is not None:
                        cross = bool(med > gstat[adj][1]) if direction == "longer" else bool(med < gstat[adj][1])
                entry.update({"discordant": disc, "direction": direction, "crosses_adjacent_group_median": cross})
                if grp in ordn:
                    n_det += 1
                    n_disc += int(disc)
            for k, v in entry.items():
                self.st(P(*base, k), v)
        for (grp, sp), ids in per_group_ids.items():
            tot = P("systems", bs, "groups", grp, sp, "n")
            reg.identity(tot, ids + [reg.hide(f"{bs}/names_below_10/{grp}/{sp}", reg.raw[tot] - sum(reg.raw[i] for i in ids))])
        self.st(P("systems", bs, "discordance_summary"), {"names_listed_ordinal": n_det, "names_discordant": n_disc,
                                                          "d5_share": r3(n_disc / n_det) if n_det else "insufficient"})
        return (n_disc / n_det) if n_det else None

    # ---------------------------------------------------------------- across-system calibration (train)
    def across_system(self) -> None:
        reg, f = self.reg, self.f
        noncard = [bs for bs in self.meta["systems"] if bs != CARDIAC]
        d = f[f["analysed"] & f["valid"] & (f["split"] == "train") & (f["body_system"] != CARDIAC)]
        levels = [(lv, "intensity", "e_four_level_or_stat5", "e") for lv in INTENSITY_ORDER] + [(lv, "cgroup", "primary", "primary") for lv in I_ORDER]
        for lv, col, fam, gk in levels:
            sysl = [bs for bs in noncard if lv in self.meta["groups"][bs][gk]]
            if not sysl:
                continue
            dl = d[d[col] == lv]
            base = ("across_system", lv)
            pn = reg.count(P(*base, "pooled_train", "n"), len(dl))
            ids = []
            for bs in sysl:
                c = reg.count(P(*base, "by_system", bs, "n"), int((dl["body_system"] == bs).sum()))
                reg.equal(c, P(*fam_base(bs, fam), "groups", lv, "train", "n"))
                ids.append(c)
            reg.identity(pn, ids)
            pq = qstats(dl["D"].to_numpy(float))
            for k in ("median", "p25", "p75"):
                self.st(P(*base, "pooled_train", k), pq[k])
            inc = []
            for bs in sysl:
                v = dl.loc[dl["body_system"] == bs, "D"].to_numpy(float)
                q = qstats(v)
                for k in ("median", "p25", "p75"):
                    self.st(P(*base, "by_system", bs, k), q[k])
                if len(v) >= MIN_CELL and len(dl) >= MIN_CELL:
                    p25, p75 = np.percentile(dl["D"], [25, 75])
                    mv = float(np.median(v))
                    self.st(P(*base, "by_system", bs, "outside_pooled_iqr"), bool(mv < p25 or mv > p75))
                    inc.append(v)
                else:
                    self.st(P(*base, "by_system", bs, "outside_pooled_iqr"), "<10")
            if len(inc) >= 2:
                h = sps.kruskal(*inc).statistic
                self.st(P(*base, "kruskal_epsilon2"), r3(h / (sum(len(v) for v in inc) - 1)))
            else:
                self.st(P(*base, "kruskal_epsilon2"), "insufficient")
        if CARDIAC in self.meta["systems"]:
            tr = f[f["analysed"] & f["valid"] & (f["split"] == "train") & (f["body_system"] == CARDIAC)]
            self.st(P("across_system", "cardiac_reference_not_tested"),
                    {g_: qstats(tr.loc[tr["cgroup"] == g_, "D"].to_numpy(float))["median"]
                     for g_ in G_ORDER if g_ in self.meta["groups"][CARDIAC]["primary"]})

    # ---------------------------------------------------------------- descriptive secondary (all rows, no cohort)
    def descriptive(self) -> None:
        reg, f = self.reg, self.f
        v = f[f["valid"]]
        vt = reg.count(P("descriptive_all_rows", "valid_total"), len(v))
        reg.subset(vt, P("reconciliation", "surgery_info_rows"), "invalid_all_rows")
        sys_ids = []
        for bs in self.meta["systems"]:
            dv = v[v["classified"] & (v["body_system"] == bs)]
            st_ = reg.count(P("descriptive_all_rows", "systems", bs, "total"), len(dv))
            gids = []
            for grp in self.meta["groups"][bs]["primary"]:
                vals = dv.loc[dv["cgroup"] == grp, "D"].to_numpy(float)
                c = reg.count(P("descriptive_all_rows", "systems", bs, "groups", grp, "n"), len(vals))
                reg.subset(P("systems", bs, "groups", grp, "pooled", "n"), c, "not_analysed")
                for k, val in qstats(vals).items():
                    self.st(P("descriptive_all_rows", "systems", bs, "groups", grp, k), val)
                gids.append(c)
            reg.identity(st_, gids)
            sys_ids.append(st_)
        reg.identity(vt, sys_ids + [reg.hide("descriptive/unclassified_or_unmapped_valid", int((~v["classified"]).sum()))])

    # ---------------------------------------------------------------- run
    def run(self) -> None:
        self.population()
        d5 = {}
        for bs in self.meta["systems"]:
            for fam in FAMILIES:
                if fam in CARDIAC_ONLY_FAMILIES and bs != CARDIAC:
                    continue
                self.family(bs, fam)
            d5[bs] = self.names(bs)
        self.across_system()
        self.descriptive()
        self.final_verdicts(d5)

    def final_verdicts(self, d5: dict) -> None:
        per_sys = {}
        for bs in self.meta["systems"]:
            for fam in FAMILIES:
                key = f"{bs}|{fam}"
                if key not in self.verdicts:
                    continue
                r = self.verdicts[key]
                if r["assessable"]:
                    v = evaluate_criteria(r["ordered_medians_train"], r["tau_b"], r.get("d3_share"), r.get("d4_share"),
                                          d5.get(bs) if fam == "primary" else None)
                else:
                    v = {k: "not_assessable" for k in ("D1", "D2", "D3", "D4", "D5")} | {"big_discrepancy": "not_assessable", "label": VERDICT_LABEL}
                tau = r["tau_b"]
                vals = {"d1_medians_train": [r1(m) for m in r["ordered_medians_train"]], "groups": r["eligible_groups_train"],
                        "tau_b": r3(tau) if np.isfinite(tau) else "insufficient", "tau_b_ci": r["tau_ci"],
                        "spearman": r3(r["spearman"]) if np.isfinite(r["spearman"]) else "insufficient",
                        "abs_tau_b": r3(abs(tau)) if np.isfinite(tau) else "insufficient",
                        "d3_share": "not_assessable", "d4_share": "not_assessable",
                        "d5_share": r3(d5[bs]) if (fam == "primary" and d5.get(bs) is not None) else "insufficient"}
                self.local[f"/verdict_raw/{key}"] = {"d3_share": r.get("d3_share"), "d4_share": r.get("d4_share"),
                                                     "d5_share": d5.get(bs), "tau_b": tau}
                vptr = P("systems", bs, "verdict") if fam == "primary" else P(*fam_base(bs, fam), "criteria_descriptive")
                if r["assessable"] and r.get("pair0"):
                    def pub_d34(status: dict, stats: dict, vptr=vptr, ob=(*r["base"], "overlap", r["pair0"], "train")) -> None:
                        stats[vptr] = dict(stats[vptr], d3_share=stats[P(*ob, "p_low_above_high_median", "share")],
                                           d4_share=stats[P(*ob, "p_high_below_low_median", "share")])
                    self.post.append(pub_d34)
                if fam == "primary":
                    self.st(P("systems", bs, "grouping"), "STAT_group" if bs == CARDIAC else "intensity_3")
                    self.st(P("systems", bs, "assessable"), r["assessable"])
                    self.st(P("systems", bs, "not_assessable_reason"), r["reason"])
                    self.st(P("systems", bs, "verdict"), v | vals)
                    per_sys[bs] = (v, r)
                else:
                    self.st(P(*fam_base(bs, fam), "criteria_descriptive"), v | vals | {"note": "sensitivity; descriptive only"})
        assess = [bs for bs, (_, r) in per_sys.items() if r["assessable"]]
        f = self.f
        strat: dict[str, float] = {}
        for sp in ("train", "pooled"):
            parts = []
            for bs in assess:
                elig = per_sys[bs][1]["eligible_groups_train"]
                ordn = self.meta["groups"][bs]["primary_ordinal"]
                d = split_rows(f[f["analysed"] & f["valid"] & (f["body_system"] == bs)], sp)
                d = d[d["cgroup"].isin(elig)]
                parts.append(pd.DataFrame({"sys": bs, "x": d["cgroup"].map({g_: i for i, g_ in enumerate(ordn)}).astype(float).to_numpy(),
                                           "y": d["D"].astype(float).to_numpy(), "subj": d["SUBJECT_ID"].to_numpy()}))
            if not parts:
                self.st(P("pooled", "stratified_tau_b", sp), {"est": "insufficient", "ci": "insufficient", "systems": []})
                strat[sp] = float("nan")
                continue
            allp = pd.concat(parts, ignore_index=True)
            sysv, xv, yv = allp["sys"].to_numpy(), allp["x"].to_numpy(), allp["y"].to_numpy()

            def stratified(idx: np.ndarray, sysv=sysv, xv=xv, yv=yv) -> float:
                num = den = 0.0
                for bs in assess:
                    sel = idx[sysv[idx] == bs]
                    t = kendall_tau_b_surgery_level(xv[sel], yv[sel])
                    if not np.isfinite(t):
                        return float("nan")
                    num += len(sel) * t
                    den += len(sel)
                return num / den if den else float("nan")
            est = stratified(np.arange(len(allp)))
            ci = self.boot.ci(allp["subj"].to_numpy(), stratified, f"pooled|stratified|{sp}")
            strat[sp] = est
            self.st(P("pooled", "stratified_tau_b", sp), {"est": r3(est) if np.isfinite(est) else "insufficient", "ci": ci["ci"],
                                                          "systems": assess, "weighting": "n-weighted within-system tau-b (pairs only within a system)"})

        def anyrule(k: str) -> str:
            vals = [per_sys[bs][0][k] for bs in assess]
            return "met" if "met" in vals else ("not_met" if "not_met" in vals else "not_assessable")
        est = strat.get("train", float("nan"))
        pv = evaluate_criteria(None, est, None, None, None)
        pv.update({"D1": anyrule("D1"), "D3": anyrule("D3"), "D4": anyrule("D4"), "D5": anyrule("D5")})
        dec = [pv["D1"], pv["D2"], pv["D3"]]
        pv["big_discrepancy"] = "Y" if "met" in dec else ("N" if "not_met" in dec else "not_assessable")
        pv["tau_b"] = r3(est) if np.isfinite(est) else "insufficient"
        pv["systems_with_big_discrepancy"] = [bs for bs in assess if per_sys[bs][0]["big_discrepancy"] == "Y"]
        pv["rule"] = "D2 from the stratified train tau-b (signed); D1/D3/D4/D5 met if met in any assessable system"
        self.st(P("pooled", "verdict"), pv)


# =============================================================================
# Shared document, markdown, release gate
# =============================================================================
def build_document(an: Analysis, status: dict, supp: dict, hashes: dict) -> dict:
    doc: dict[str, Any] = {
        "spec": SPEC_VERSION,
        "input_hashes": dict(hashes),
        "map_status": MAP_STATUS, "map_label": MAP_LABEL, "map_reviewed_date": MAP_REVIEWED_DATE,
        "criteria": criteria_config(),
        "zero_policy": "counts 0-9 are published as '<10' (spec v1.5 Suppression: 'counts of 1-9 and zero')",
    }
    if an.boot.runs:  # metadata honesty: only if a bootstrap actually ran
        doc["bootstrap"] = {"replicates": an.boot.reps, "cluster": "SUBJECT_ID", "seed": an.boot.seed,
                            "ci": "percentile 95%", "bootstraps_run": an.boot.runs,
                            "realised_replicates_min": int(min(an.boot.realised)),
                            "realised_replicates_max": int(max(an.boot.realised))}
    stats = dict(an.stats)
    for fn in an.post:
        fn(status, stats)
    for ptr, v in stats.items():
        set_ptr(doc, ptr, v)
    for cid in an.reg.published():
        st = status[cid]
        set_ptr(doc, cid, an.reg.raw[cid] if st == "pub" else ("<10" if st == "P" else "S"))
    doc["suppression"] = dict(supp)
    doc["disclosure"] = {"count_fields": sorted(an.reg.published()),
                         "identities": [[t, list(p)] for t, p in an.reg.identities],
                         "hidden_cells": sorted(an.reg.hidden),
                         "markers": {"primary": "<10", "complementary": "S"}}
    doc["release_gate"] = "pending"
    return doc


def detect_count_fields(obj: Any, prefix: str = "") -> list[str]:
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if prefix == "" and k in ("disclosure", "suppression", "criteria", "bootstrap", "input_hashes"):
                continue
            p = prefix + "/" + str(k).replace("~", "~0").replace("/", "~1")
            if isinstance(v, (dict, list)):
                out += detect_count_fields(v, p)
            elif k in COUNT_KEYS:
                out.append(p)
    return out


def _fmt(v: Any) -> str:
    if v is None:
        return "–"
    if isinstance(v, bool):
        return "Y" if v else "N"
    if isinstance(v, float):
        return f"{v:.3f}".rstrip("0").rstrip(".") if abs(v) < 1 else f"{v:.1f}"
    if isinstance(v, list):
        return "[" + ", ".join(_fmt(x) for x in v) + "]"
    return str(v)


def _est_ci(o: Any) -> str:
    if isinstance(o, dict):
        return f"{_fmt(o.get('est'))} {_fmt(o.get('ci'))}"
    return _fmt(o)


def render_markdown(doc: dict) -> str:
    """Human-readable summary rendered ONLY from the shared JSON (the gate re-renders and compares)."""
    sup = doc["suppression"]
    L = [f"# Surgery complexity vs duration audit: {doc['map_label']} (reviewed {doc['map_reviewed_date']})", "",
         f"- Spec `{doc['spec']}`; verdict label `{VERDICT_LABEL}`; release gate: **{doc['release_gate']}**.",
         "- Input hashes (exact 64-character match): " + ", ".join(f"{k} `{v}`" for k, v in doc["input_hashes"].items()) + ".",
         "- Aggregate only. Counts of 0–9 are `<10`; `S` = complementary suppression. Quantiles only for groups with n ≥ 10; never min/max.",
         f"- Suppression: {sup['primary_cells']} primary cells, {sup['complementary_cells']} complementary cells; "
         f"{sup['recoverability_cells_checked']} cells checked by the recoverability solver; {sup['recoverable']} recoverable."]
    b = doc.get("bootstrap")
    if b:
        L.append(f"- Bootstrap: {b['replicates']} replicates, cluster {b['cluster']}, seed {b['seed']}, {b['ci']}; "
                 f"{b['bootstraps_run']} bootstraps run; realised replicates {b['realised_replicates_min']}–{b['realised_replicates_max']}.")
    rc = doc["reconciliation"]
    L += ["", "## 1. Reconciliation (all SURGERY_INFO rows)", "", "| class | count |", "|---|---|"]
    for k in ("surgery_info_rows",) + LINK_CLASSES:
        L.append(f"| {k} | {_fmt(rc[k])} |")
    L += ["", "| multi-candidate surgeries | count |", "|---|---|"]
    for k in ("n_multi_candidate", "resolved_rule1_during", "resolved_rule2_next_t0", "resolved_rule3_post_icu"):
        L.append(f"| {k} | {_fmt(rc[k])} |")
    L += ["", "## 1b. Population and validity by split", "",
          "| split | linked | recent_layer_a | remote_same_hadm | during_icu_layer_c | analysed | unmapped | excluded_unclassified | valid_duration | invalid_total | duration_fallback_anes | d_gt_24h |",
          "|---|" + "---|" * 11]
    for sp in ALL_SPLITS:
        p = doc["population"][sp]
        L.append(f"| {sp} | {_fmt(p['linked'])} | " + " | ".join(_fmt(p["by_timing_class"][c]) for c in TIMING)
                 + f" | {_fmt(p['analysed'])} | {_fmt(p['unmapped'])} | {_fmt(p['excluded_unclassified'])} | {_fmt(p['valid_duration'])} | "
                 f"{_fmt(p['invalid_duration']['invalid_total'])} | {_fmt(p['duration_fallback_anes'])} | {_fmt(p['d_gt_24h'])} |")
    L += ["", "## 2–5. Per system (primary grouping)", ""]
    for bs, s in doc["systems"].items():
        L += [f"### {bs} ({s['grouping']})", "",
              "Assessable: " + ("yes" if s["assessable"] else f"no ({s['not_assessable_reason']})"), "",
              "| group | split | n | median | p25 | p75 | p10 | p90 | lt2h | 2to4h | gt4h |", "|---|" + "---|" * 10]
        for grp, gs in s["groups"].items():
            for sp in ALL_SPLITS:
                r = gs[sp]
                L.append(f"| {grp} | {sp} | {_fmt(r['n'])} | " + " | ".join(_fmt(r[k]) for k in ("median", "p25", "p75", "p10", "p90"))
                         + " | " + " | ".join(_fmt(r["bands"][bb]) for bb in BANDS) + " |")
        for sp in ALL_SPLITS:
            ct = s["crosstab_group_band"][sp]
            L.append(f"| all groups | {sp} | {_fmt(ct['total'])} | | | | | | " + " | ".join(_fmt(ct["band_totals"][bb]) for bb in BANDS) + " |")
        L += ["", "| split | Kendall tau-b (signed) [95% CI] | Spearman rho [95% CI] | tau-b from crosstab |", "|---|---|---|---|"]
        for sp in ALL_SPLITS:
            c = s.get("correlation", {}).get(sp, {})
            L.append(f"| {sp} | {_est_ci(c.get('kendall_tau_b', 'insufficient'))} | {_est_ci(c.get('spearman', 'insufficient'))} | {_fmt(c.get('tau_b_from_crosstab', '–'))} |")
        if s.get("overlap"):
            L += ["", "| pair | split | share low above high train median | share high below low train median | prob. superiority [95% CI] |", "|---|---|---|---|---|"]
            for pk, ov in s["overlap"].items():
                for sp in ("train", "pooled"):
                    o = ov[sp]
                    L.append(f"| {pk} | {sp} | {_fmt(o['p_low_above_high_median']['share'])} | {_fmt(o['p_high_below_low_median']['share'])} | {_est_ci(o['prob_superiority'])} |")
        v = s["verdict"]
        L += ["", f"**Verdict** (`{v['label']}`): D1 {v['D1']}; D2 {v['D2']} (signed tau-b {_fmt(v['tau_b'])}, CI {_fmt(v['tau_b_ci'])}); "
              f"D3 {v['D3']} (share {_fmt(v['d3_share'])}); D4 {v['D4']} (reported only); D5 {v['D5']} (reported only) "
              f"→ **big_discrepancy = {v['big_discrepancy']}**. Train medians {_fmt(v['d1_medians_train'])} for {_fmt(v['groups'])}.", ""]
        dn = s.get("discordant_names", {})
        if dn:
            L += ["Names with at least 10 surgeries (pooled analysis population):", "",
                  "| name | group | n | median pooled | median train | discordant | direction | crosses adjacent median |", "|---|---|---|---|---|---|---|---|"]
            for nm, e in dn.items():
                L.append(f"| {nm} | {e['group']} | {_fmt(e['pooled']['n'])} | {_fmt(e['median_pooled'])} | {_fmt(e['median_train'])} | "
                         f"{_fmt(e['discordant'])} | {_fmt(e['direction'])} | {_fmt(e['crosses_adjacent_group_median'])} |")
            L.append("")
    L += ["## 9. Across-system check (train; non-cardiac systems)", "",
          "| level | pooled median (IQR) | system | median | outside pooled IQR | KW epsilon² |", "|---|---|---|---|---|---|"]
    for lv, a in doc.get("across_system", {}).items():
        if lv == "cardiac_reference_not_tested":
            continue
        pt = a["pooled_train"]
        for bs, e in a["by_system"].items():
            L.append(f"| {lv} | {_fmt(pt['median'])} ({_fmt(pt['p25'])}–{_fmt(pt['p75'])}) | {bs} | {_fmt(e['median'])} | {_fmt(e['outside_pooled_iqr'])} | {_fmt(a['kruskal_epsilon2'])} |")
    if "cardiac_reference_not_tested" in doc.get("across_system", {}):
        L += ["", f"Cardiac reference (not tested): train medians {_fmt(doc['across_system']['cardiac_reference_not_tested'])}."]
    pv = doc["pooled"]["verdict"]
    L += ["", "## 10. Verdicts", "", "| system | D1 | D2 | D3 | D4 (reported) | D5 (reported) | big_discrepancy |", "|---|---|---|---|---|---|---|"]
    for bs, s in doc["systems"].items():
        v = s["verdict"]
        L.append(f"| {bs} | {v['D1']} | {v['D2']} | {v['D3']} | {v['D4']} | {v['D5']} | {v['big_discrepancy']} |")
    L.append(f"| pooled (stratified tau-b {_fmt(pv['tau_b'])}) | {pv['D1']} | {pv['D2']} | {pv['D3']} | {pv['D4']} | {pv['D5']} | {pv['big_discrepancy']} |")
    na = [f"{bs} ({s['not_assessable_reason']})" for bs, s in doc["systems"].items() if not s["assessable"]]
    L += ["", "Not assessable: " + (", ".join(na) if na else "none") + ".", "",
          "Small cells: train sets every threshold and criterion; validation and test are reported for completeness, and pooled figures are "
          "descriptive only. Many split-level cells are `<10` or `S` by design. Sensitivity analyses (a–e), flagged-vs-unflagged blocks and the "
          "descriptive all-rows block are in the JSON.", ""]
    return "\n".join(L)


TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}")


def release_gate(json_path: Path, md_path: Path, expect_status: str, sentinels: Iterable[str] = ()) -> list[str]:
    """Disk-reload validator. Re-reads the shared JSON and markdown and independently checks: count-field schema,
    marker values (no published count 0-9), suppression counters, published identities, LP recoverability, quantiles
    of small groups, min/max, timestamps, sentinel IDs and the markdown. Returns failures (empty list = pass)."""
    fails: list[str] = []
    try:
        doc = json.loads(Path(json_path).read_text(encoding="utf-8"))
        md = Path(md_path).read_text(encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        return [f"cannot reload shared files: {e}"]
    if doc.get("spec") != SPEC_VERSION or doc.get("map_label") != MAP_LABEL or doc.get("release_gate") != expect_status:
        fails.append("header fields (spec / map_label / release_gate) wrong")
    if not Path(json_path).name.split(".json")[0].removesuffix(".pending").endswith(f"_{MAP_LABEL}"):
        fails.append("shared JSON name does not end in _MAP_v1")
    declared = set(doc.get("disclosure", {}).get("count_fields", []))
    detected = set(detect_count_fields(doc))
    if declared != detected:
        fails.append(f"count-field schema mismatch: {len(detected - declared)} undeclared, {len(declared - detected)} missing")
    known: dict[str, int] = {}
    unknown: dict[str, str] = {}
    n_p = n_s = 0
    for ptr in sorted(detected):
        v = get_ptr(doc, ptr)
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
            checked, unsafe, _ = recoverability_check(ids, known, unknown)
            if unsafe:
                fails.append(f"recoverability: {len(unsafe)} suppressed cell(s) recoverable or bounded inside [1, 9]")
            if checked != sup.get("recoverability_cells_checked"):
                fails.append("recoverability_cells_checked differs from the recomputed value")
        except ReleaseGateError as e:
            fails.append(f"recoverability solver: {e}")

    def walk(o: Any, path: str = "") -> None:
        if isinstance(o, dict):
            if o.get("n") == "<10":
                for q in ("median", "p25", "p75", "p10", "p90"):
                    if q in o and o[q] != "<10":
                        fails.append(f"quantile published for a group under 10 at {path}/{q}")
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
    try:
        if render_markdown(doc) != md:
            fails.append("markdown on disk is not the rendering of the shared JSON")
    except Exception as e:  # noqa: BLE001
        fails.append(f"markdown cannot be re-rendered: {e}")
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
                if h in COUNT_KEYS and c not in ("<10", "S", "", "–") and not (c.isdigit() and int(c) >= MIN_CELL):
                    fails.append(f"markdown count column '{h}' has a disallowed value")
        elif not line.startswith("|---"):
            headers = None
    return fails


def _atomic_write(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _remove(*paths: Path) -> None:
    for p in paths:
        try:
            Path(p).unlink()
        except FileNotFoundError:
            pass


def publish(cfg: AuditConfig, doc: dict, sentinels: Iterable[str] = ()) -> None:
    """Write pending files, run the disk-reload gate, then write the final files (release_gate = ok, set by the
    validator) and re-run the gate on them. On any failure every shared file is removed and ReleaseGateError raised."""
    pj = cfg.shared_json.with_name(cfg.shared_json.name.replace(".json", ".pending.json"))
    pm = cfg.shared_md.with_name(cfg.shared_md.name.replace(".md", ".pending.md"))
    _remove(cfg.shared_json, cfg.shared_md, pj, pm)
    _atomic_write(pj, json.dumps(doc, indent=1, ensure_ascii=False))
    _atomic_write(pm, render_markdown(doc))
    fails = release_gate(pj, pm, "pending", sentinels)
    if fails:
        _remove(pj, pm)
        raise ReleaseGateError("; ".join(fails[:20]))
    final = json.loads(pj.read_text(encoding="utf-8"))
    final["release_gate"] = "ok"
    _atomic_write(cfg.shared_json, json.dumps(final, indent=1, ensure_ascii=False))
    _atomic_write(cfg.shared_md, render_markdown(final))
    _remove(pj, pm)
    fails = release_gate(cfg.shared_json, cfg.shared_md, "ok", sentinels)
    if fails:
        _remove(cfg.shared_json, cfg.shared_md)
        raise ReleaseGateError("; ".join(fails[:20]))


# =============================================================================
# Orchestration
# =============================================================================
def run_audit(cfg: AuditConfig, sentinels: Iterable[str] = ()) -> dict:
    inp = load_inputs(cfg)
    frame, meta = build_frame(inp)
    rc = {c: int((frame["link_class"] == c).sum()) for c in LINK_CLASSES}
    if sum(rc.values()) != cfg.expected_surgery_rows or rc["linked"] > cfg.expected_surgery_rows:
        raise ReconciliationError("link classes do not sum exactly to the SURGERY_INFO row count")
    an = Analysis(cfg, frame, meta, inp["info"])
    an.run()  # every Registry.identity call is a raw-count reconciliation check and raises on failure
    log(cfg, f"reconciliation ok ({len(an.reg.identities)} identities verified on raw counts; link classes sum to SURGERY_INFO rows)")
    status, supp = apply_suppression(an.reg)
    log(cfg, f"suppression computed: {supp['primary_cells']} primary, {supp['complementary_cells']} complementary, "
             f"{supp['recoverability_cells_checked']} checked, {supp['recoverable']} recoverable")
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    local = {"WARNING": "LOCAL ONLY - unsuppressed aggregates; never copy, upload or share",
             "raw_counts": an.reg.raw, "extras": an.local, "checks": inp["info"]["checks"], "link_class_counts": rc}
    _atomic_write(cfg.local_json, json.dumps(local, indent=1, default=str))
    doc = build_document(an, status, supp, inp["info"]["hashes"])
    publish(cfg, doc, sentinels)
    final = json.loads(cfg.shared_json.read_text(encoding="utf-8"))
    print_console_summary(cfg, final, inp["info"]["checks"])
    return final


def print_console_summary(cfg: AuditConfig, doc: dict, checks: dict) -> None:
    """Console summary prints ONLY published (suppressed) values, so it can be pasted back."""
    if not cfg.verbose:
        return
    print(f"\n=== {doc['spec']} | {doc['map_label']} (reviewed {doc['map_reviewed_date']}) | release_gate={doc['release_gate']} ===")
    print("input hashes: " + ", ".join(f"{k} OK ({v[:12]}...)" for k, v in doc["input_hashes"].items()))
    print(f"harmonization columns used: t0={checks['harmonization_columns']['t0']!r}, OUT={checks['harmonization_columns']['OUT']!r}; "
          f"ADMISSIONS read: {checks['admissions_read']}; surgery key: {checks.get('surgery_key', SURGERY_KEY)}")
    if "OUT_minus_t0_ne_span_icu_hours_gt_1min" in checks:
        v = checks["OUT_minus_t0_ne_span_icu_hours_gt_1min"]
        print(f"OUT - t0 vs span_icu_hours (>1 min) disagreements: {'<10' if 0 < v < 10 else v}")
    kd = doc["reconciliation"].get("surgery_key_ROW_ID", {})
    if kd:
        print("surgery key ROW_ID (nothing dropped): " + ", ".join(f"{k}={v}" for k, v in kd.items()))
    rc = doc["reconciliation"]
    print("reconciliation: " + ", ".join(f"{k}={rc[k]}" for k in ("surgery_info_rows",) + LINK_CLASSES))
    print("multi-candidate: " + ", ".join(f"{k}={rc[k]}" for k in ("n_multi_candidate", "resolved_rule1_during", "resolved_rule2_next_t0", "resolved_rule3_post_icu")))
    p = doc["population"]["pooled"]
    print(f"pooled population: linked={p['linked']} analysed={p['analysed']} unmapped={p['unmapped']} excluded_unclassified={p['excluded_unclassified']} "
          f"invalid={p['invalid_duration']} fallback_anes={p['duration_fallback_anes']} d_gt_24h={p['d_gt_24h']}")
    for bs, s in doc["systems"].items():
        v = s["verdict"]
        ass = "assessable" if s["assessable"] else f"not_assessable ({s['not_assessable_reason']})"
        print(f"  {bs:34s} {ass:32s} D1={v['D1']:14s} D2={v['D2']:14s} D3={v['D3']:14s} D4={v['D4']:14s} D5={v['D5']:14s} "
              f"big={v['big_discrepancy']:14s} tau_b={v['tau_b']} ci={v['tau_b_ci']} d3_share={v['d3_share']} medians={v['d1_medians_train']}")
    pv = doc["pooled"]["verdict"]
    print(f"  {'POOLED':34s} {'':32s} D1={pv['D1']:14s} D2={pv['D2']:14s} D3={pv['D3']:14s} D4={pv['D4']:14s} D5={pv['D5']:14s} "
          f"big={pv['big_discrepancy']:14s} stratified_tau_b={pv['tau_b']}")
    s = doc["suppression"]
    print(f"suppression: primary={s['primary_cells']} complementary={s['complementary_cells']} "
          f"checked={s['recoverability_cells_checked']} recoverable={s['recoverable']}")
    b = doc.get("bootstrap")
    if b:
        print(f"bootstrap: replicates={b['replicates']} seed={b['seed']} runs={b['bootstraps_run']} realised={b['realised_replicates_min']}-{b['realised_replicates_max']}")
    print("release gate ok")


def build_config(argv: list[str] | None = None) -> AuditConfig:
    here = Path(__file__).resolve()
    default_root = Path(os.environ.get("FOUNDATIONS_ROOT", str(here.parents[1])))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", default=str(default_root))
    ap.add_argument("--raw-dir", default=os.environ.get("PIC_RAW_DIR", DEFAULT_RAW_DIR))
    ap.add_argument("--split-file", default=os.environ.get("SPLIT_FILE_PATH", DEFAULT_SPLIT_FILE))
    ap.add_argument("--canonical", default=os.environ.get("CANONICAL_PATH"))
    ap.add_argument("--harmonization", default=os.environ.get("HARMONIZATION_PATH"))
    ap.add_argument("--map", default=os.environ.get("COMPLEXITY_MAP_PATH"))
    ap.add_argument("--out-dir", default=os.environ.get("AUDIT_OUT_DIR"))
    ap.add_argument("--harm-t0-col", default=os.environ.get("HARM_T0_COL"))
    ap.add_argument("--harm-out-col", default=os.environ.get("HARM_OUT_COL"))
    ap.add_argument("--split-col", default=os.environ.get("SPLIT_COL", "split"))
    ap.add_argument("--bootstrap-reps", type=int, default=int(os.environ.get("BOOTSTRAP_REPS", "2000")))
    a = ap.parse_args(argv)
    if not a.raw_dir or not a.split_file:
        raise AuditError("--raw-dir and --split-file are required (or set PIC_RAW_DIR and SPLIT_FILE_PATH); no local default path is shipped")
    root = Path(a.project_root)
    interim = root / "data" / "interim"
    raw = Path(a.raw_dir)
    return AuditConfig(
        canonical_path=Path(a.canonical) if a.canonical else interim / "canonical_icu_encounters.csv",
        harmonization_path=Path(a.harmonization) if a.harmonization else interim / "encounter_harmonization_fields.csv",
        map_path=Path(a.map) if a.map else root / "config" / "surgery_complexity_map_v1.csv",
        split_path=Path(a.split_file), surgery_info_path=raw / "SURGERY_INFO.csv", admissions_path=raw / "ADMISSIONS.csv",
        out_dir=Path(a.out_dir) if a.out_dir else interim, harm_t0_col=a.harm_t0_col, harm_out_col=a.harm_out_col,
        split_col=a.split_col, bootstrap_reps=a.bootstrap_reps)


def _force_utf8_console() -> None:
    """Windows consoles default to cp1252; make stdout/stderr UTF-8 so symbols like ≥ never crash a run."""
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
