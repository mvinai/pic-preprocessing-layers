#!/usr/bin/env python3
"""
Stage 2 v1.0.2 release: differencing (complementary) suppression fix for the shared Stage 2 v1.0.1 output.
Implements `stage2_v1_0_2_suppression_spec.md` v1.0 (Nov 16, 2026).

Reads ONLY the released, already-suppressed Stage 2 v1.0.1 aggregate output
(feature_sparsity_frequency_audit.json and .md). It never reads patient data or raw PIC tables, never imports the
Stage 2 analysis script, never writes over its inputs, and never opens any *_raw_LOCAL* file.

What it does, in order (any failure stops the run; no release file is left behind):
  1. Input gates: file names (no _raw_LOCAL), optional exact SHA-256 pins, spec string and release_gate == ok, the
     declared count-field schema, every published count an integer >= 10 (else rule 1 suppresses it), every identity,
     link and frozen pointer resolvable.
  2. Input audit (what v1.0.1 exposes): interval / LP audit of every suppressed cell under all published identities
     plus the declared cross-table links, and the residual (differencing) check of every additive group.
  3. Rule 1 (primary): every count 0-9 stays '<10'.
     Rule 2 (complementary): while any cell fails the audit, suppress one more published cell per unsafe component,
     chosen by a short look-ahead (fewest extra suppressions; ties -> the group total, then the smallest value).
     Frozen (already public) cells are never chosen. Derived values (pct_, q_, ks_, x_) whose counts are no longer all
     published become 'suppressed'.
  4. Rule 3 (audit): every suppressed cell's feasible integer interval from all published equations and bounds must
     contain >= W integers (default W = 10) or contain values both < 10 and >= 10; '<10' cells must keep their full
     0-9 range; hidden (unpublished) differences pass when provably >= 10 or exactly 0. Residual check: whenever a
     group total and some of its parts are published, total minus published parts (the sum of the suppressed parts)
     must be >= 10.
  5. Rule 4 (release gate): pending files -> reload from disk -> all checks -> final names -> manifest (SHA-256 of
     inputs and outputs) and a summary report (counts of cells only, plus pointers of newly suppressed cells, never
     values). Any failure: every output written by this run is removed and the exit code is 3.

Usage (Windows, from the project root):
    python scripts\\stage2_v1_0_2_suppression\\stage2_v1_0_2_suppress.py ^
        --in-json data\\interim\\feature_sparsity_frequency_audit.json ^
        --in-md   data\\interim\\feature_sparsity_frequency_audit.md ^
        --links   scripts\\stage2_v1_0_2_suppression\\stage2_v1_0_2_links.csv ^
        --frozen  scripts\\stage2_v1_0_2_suppression\\stage2_v1_0_2_frozen.csv ^
        --out-dir data\\interim\\stage2_v1.0.2
    Add --check-only to audit an output without writing anything (exit 0 = passes the new gate, 3 = fails).

Outputs (in --out-dir, which must be new or empty and must not be an input folder):
    feature_sparsity_frequency_audit_v1.0.2.json / .md   shared (suppressed)
    stage2_v1.0.2_release_summary.md                     counts of cells by rule; pointers only, no values
    stage2_v1.0.2_manifest.json                          SHA-256 of inputs and outputs
Exit codes: 0 ok; 2 input/schema failure (nothing written); 3 audit/release gate failure (nothing left behind).
Dependencies: Python 3.10+, numpy, pandas, scipy (HiGHS). Deterministic: same inputs -> byte-identical outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog

SCRIPT_VERSION = "stage2_v1_0_2_suppress v1.3"
SPEC_DOC = "stage2_v1_0_2_suppression_spec v1.3"
EXPECT_INPUT_SPEC = "feature_sparsity_frequency_audit v1.0.1 (feature_sparsity_frequency_audit_spec v0.8, stage 2)"
OUTPUT_SPEC = ("feature_sparsity_frequency_audit v1.0.2 (v1.0.1 output with the differencing fix; "
               "stage2_v1_0_2_suppression_spec v1.0, stage 2)")
OUT_JSON = "feature_sparsity_frequency_audit_v1.0.2.json"
OUT_MD = "feature_sparsity_frequency_audit_v1.0.2.md"
OUT_SUMMARY = "stage2_v1.0.2_release_summary.md"
OUT_MANIFEST = "stage2_v1.0.2_manifest.json"
PROVENANCE_KEY = "release_v1_0_2"

MIN_CELL = 10
UNKNOWN = ("P", "S", "H", "D")    # statuses whose value a reader does not see (D: dropped by the --drops list)
SMALL_MAX = MIN_CELL - 1          # '<10' covers 0..9 (the Stage 2 primary rule suppresses zeros too)
INF = math.inf
BIG = 1.0e9                       # cap used only for screening LPs (points stay feasible for the uncapped system)
EPS = 1e-6
LOOKAHEAD_K = 4                   # candidates evaluated by look-ahead per unsafe component
LOOKAHEAD_CAP = 25                # max simulated suppressions per look-ahead
LOOKAHEAD_MAX_LINKS = 50          # a candidate in more groups than this is a last resort (never simulated)
MAX_ROUNDS = 2000
LARGE_COMPONENT = 300             # components with more unknown cells than this: per-group fixes, no look-ahead
LOOKAHEAD_MAX_RELS = 2000         # a look-ahead whose sub-system exceeds this many equations is abandoned
EXACT_SINGLES = 12                # large components: candidates tried one by one with an exact LP check (static order)
EXACT_TOP_VALUE = 5               # ... plus this many of the largest-valued candidates (most cover)
EXACT_RINGS = 3                   # candidate sets: cells in the failing cell's groups, then 1-2 groups further out
STALL_ROUNDS = 5                  # no-progress guard: rounds without a new low failing count before escalating
#                                   (the v1.0 look-ahead path legitimately goes 3 rounds without one on synthetic cases)

# Stage 2 v1.0.1 schema constants (copied from feature_sparsity_frequency_audit.py v1.0.1, sha256 3b5b36fb...)
SPLITS = ("train", "validation", "test")
DOB_STATUS = ("missing", "unparseable", "negative", "valid")
DOB_NEG_BINS = ("lt1d", "1to30d", "gt30d")
WEIGHT_SOURCES = ("charted", "imputed_age_sex", "imputed_age_band", "imputed_sex_only", "imputed_overall")
NON_COUNT_TOP = ("disclosure", "suppression", "input_hashes", "v17_module", "raw_tables", "maps", "cohort_checks",
                 "config", "outcome_descriptors_flag", "spec", "release_gate", "provisional_classes", "notes",
                 PROVENANCE_KEY)
COUNT_KEY_RE = re.compile(r"^n(_|$)")
TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}")
QKEYS = ("p10", "p25", "p50", "p75", "p90")


class InputError(Exception):
    """Exit 2: an input, schema or parameter problem. Nothing is written."""


class GateError(Exception):
    """Exit 3: the audit or the release gate failed. Every output of this run is removed."""


# =============================================================================
# Small helpers
# =============================================================================
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ptr_parts(ptr: str) -> list[str]:
    return [p.replace("~1", "/").replace("~0", "~") for p in ptr.split("/")[1:]]


def get_ptr(doc: Any, ptr: str) -> Any:
    cur = doc
    for p in ptr_parts(ptr):
        if not isinstance(cur, dict) or p not in cur:
            raise KeyError(ptr)
        cur = cur[p]
    return cur


def set_ptr(doc: dict, ptr: str, value: Any) -> None:
    parts = ptr_parts(ptr)
    cur = doc
    for p in parts[:-1]:
        cur = cur[p]
    if parts[-1] not in cur:
        raise KeyError(ptr)
    cur[parts[-1]] = value


def r1(x: float) -> float:
    return float(round(float(x), 1))


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


def dump_json(doc: dict) -> str:
    # same serialisation as the Stage 2 script's shared JSON
    return json.dumps(doc, indent=1, ensure_ascii=False)


def write_text(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8", newline="")
    os.replace(tmp, path)


# =============================================================================
# Markdown renderer: verbatim copy of render_markdown (and _f, _g) from feature_sparsity_frequency_audit.py v1.0.1
# (sha256 3b5b36fb...). The gate proves the copy is faithful: rendering the v1.0.1 JSON must reproduce the v1.0.1
# markdown byte for byte, otherwise the run stops (exit 3, "renderer drift").
# =============================================================================
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


# =============================================================================
# Inputs: links (cross-table equations) and frozen (already public) cells
# =============================================================================
def read_links(path: Path | None) -> list[tuple[str, tuple[str, ...], str]]:
    """CSV columns: total, parts, note. `parts` is ';'-separated. One part = an equality total == part."""
    if path is None:
        return []
    df = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    need = {"total", "parts", "note"}
    if not need.issubset(df.columns):
        raise InputError(f"links file {path.name} needs columns {sorted(need)}")
    out = []
    for i, row in df.iterrows():
        t = row["total"].strip()
        parts = tuple(p.strip() for p in row["parts"].split(";") if p.strip())
        if not t.startswith("/") or not parts or not all(p.startswith("/") for p in parts):
            raise InputError(f"links file {path.name}, row {i + 2}: pointers must start with '/'")
        out.append((t, parts, row["note"].strip()))
    return out


def read_frozen(path: Path | None) -> list[tuple[str, str]]:
    """CSV columns: pointer, reason. Cells whose values are already public elsewhere (never suppressed here)."""
    if path is None:
        return []
    df = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    if not {"pointer", "reason"}.issubset(df.columns):
        raise InputError(f"frozen file {path.name} needs columns ['pointer', 'reason']")
    return [(r["pointer"].strip(), r["reason"].strip()) for _, r in df.iterrows() if r["pointer"].strip()]


DROP_FORBIDDEN_TOP = ("disclosure", "suppression", "input_hashes", "raw_tables", "spec", "release_gate", "notes",
                      "outcome_descriptors_flag", PROVENANCE_KEY)


def read_drops(path: Path | None) -> list[tuple[str, str]]:
    """CSV columns: subtree, reason. Each row names one JSON subtree that is removed from the release as a whole
    (every field under it: counts, percentages, quantiles, everything). Totals above it are kept."""
    if path is None:
        return []
    df = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    if not {"subtree", "reason"}.issubset(df.columns):
        raise InputError(f"drops file {path.name} needs columns ['subtree', 'reason']")
    out = []
    for _, r in df.iterrows():
        t = r["subtree"].strip()
        if not t:
            continue
        t = t.rstrip("/")
        if not t.startswith("/") or not r["reason"].strip():
            raise InputError(f"drops file row needs a JSON pointer starting with '/' and a reason: {t!r}")
        out.append((t, r["reason"].strip()))
    if not out:
        raise InputError(f"drops file {path.name} lists no subtree")
    return out


def _under(ptr: str, sub: str) -> bool:
    return ptr == sub or ptr.startswith(sub + "/")


def apply_drops(L: Loaded, drops: list[tuple[str, str]], link_cells: set[str], frozen: set[str]) -> None:
    """Mark every count cell under a dropped subtree as D (unpublished; still an unknown in every equation, audited
    with the hidden-cell rule) and build the output registry without the dropped pointers. Fails on anything
    ambiguous: a missing or non-object subtree, overlapping rows, a subtree holding no count cell, a frozen cell or
    a links-file cell inside a dropped subtree."""
    if not drops:
        L.out_disclosure = L.doc.get("disclosure")
        return
    subs = [t for t, _ in drops]
    for i, t in enumerate(subs):
        top = ptr_parts(t)[0] if ptr_parts(t) else ""
        if top in DROP_FORBIDDEN_TOP or len(ptr_parts(t)) < 2:
            raise InputError(f"drops file subtree is not a droppable data subtree: {t}")
        try:
            v = get_ptr(L.doc, t)
        except KeyError:
            raise InputError(f"drops file subtree not found in the input JSON: {t}")
        if not isinstance(v, dict):
            raise InputError(f"drops file pointer is not a subtree (JSON object): {t}")
        for j, u in enumerate(subs):
            if i != j and _under(t, u):
                raise InputError(f"drops file rows overlap: {t} and {u}")

    def dropped_ptr(c: str) -> bool:
        return any(_under(c, t) for t in subs)
    dropped = [c for c in L.count_fields if dropped_ptr(c)]
    for t in subs:
        if not any(_under(c, t) for c in dropped):
            raise InputError(f"drops file subtree holds no count cell: {t}")
    for c in sorted(frozen):
        if dropped_ptr(c):
            raise InputError(f"a frozen (already public) cell lies inside a dropped subtree: {c}")
    for c in sorted(link_cells):
        if dropped_ptr(c):
            raise InputError(f"a links-file cell lies inside a dropped subtree: {c}")
    for c in dropped:
        L.status[c] = "D"
    L.dropped = dropped
    L.primary_new = [c for c in L.primary_new if c not in set(dropped)]
    L.drop_subtrees = list(drops)
    dset = set(dropped)

    def names_dropped(x: str) -> bool:   # a hidden cell named after a dropped pointer
        return any(t + "/" in x or x.endswith(t) for t in subs)
    dis = json.loads(json.dumps(L.doc["disclosure"]))
    dis["count_fields"] = [c for c in dis["count_fields"] if c not in dset]
    kept_ids = []
    for t, parts in dis["identities"]:
        cells = [t] + list(parts)
        if any(c in dset or (str(c).startswith("hidden:") and names_dropped(c)) for c in cells):
            continue
        kept_ids.append([t, list(parts)])
    dis["identities"] = kept_ids
    used = {c for t, parts in kept_ids for c in [t] + list(parts)}
    dis["hidden_cells"] = [h for h in dis["hidden_cells"] if h in used and not names_dropped(h)]
    for key in ("ratios", "quantiles", "dependents"):
        if key in dis:
            dis[key] = [e for e in dis[key] if not dropped_ptr(e[0])]
    L.out_disclosure = dis


# =============================================================================
# The equation system and the interval (LP) audit
# =============================================================================
@dataclass
class Params:
    width: int = 10              # protection width W (integers in the feasible interval)
    s_rule: str = "width"        # 'width' (default) or 'exact' (v1.7 rule for S cells: unsafe only if pinned)
    s_lower: int = MIN_CELL      # a priori lower bound of an 'S' cell (a reader knows S >= 10)


def apriori(kind: str, p: Params) -> tuple[float, float]:
    if kind == "P":
        return 0.0, float(SMALL_MAX)
    if kind == "S":
        return float(p.s_lower), INF
    return 0.0, INF  # H, and D (a dropped, unpublished cell: a reader knows only that it is >= 0)


def int_bounds(lo: float, hi: float) -> tuple[float, float]:
    lo_i = math.ceil(lo - EPS) if lo != -INF else -INF
    hi_i = math.floor(hi + EPS) if hi != INF else INF
    return lo_i, hi_i


def is_safe(kind: str, lo: float, hi: float, p: Params) -> bool:
    """Rule 3 on the integer interval [lo, hi] (hi may be inf)."""
    lo, hi = int_bounds(lo, hi)
    n_int = INF if hi == INF else hi - lo + 1
    straddle = lo <= SMALL_MAX and hi >= MIN_CELL
    if kind == "P":   # a priori 0..9: the published equations must not narrow it (>= min(W, 10) integers)
        return n_int >= min(p.width, SMALL_MAX + 1)
    if kind == "S":
        if p.s_rule == "exact":
            return n_int >= 2
        return straddle or n_int >= p.width
    # H: hidden, never-published difference; not small if provably >= 10 or exactly 0
    # D: a cell dropped from the release (--drops); unpublished, so it gets the same rule as H
    if lo >= MIN_CELL or hi <= 0:
        return True
    return straddle or n_int >= p.width


def single_equation_bounds(u: dict[str, float], b: float,
                           box: dict[str, tuple[float, float]]) -> dict[str, tuple[float, float]]:
    """Exact feasible interval of every variable of ONE equation sum(a_i x_i) = b over a box (the projection of a
    box cut by a hyperplane onto one axis is the interval below, intersected with that axis's box)."""
    terms = {}
    for v, a in u.items():
        lo, hi = box[v]
        x1, x2 = a * lo, (a * hi if hi != INF else (INF if a > 0 else -INF))
        terms[v] = (min(x1, x2), max(x1, x2))
    lo_sum = sum(t[0] for t in terms.values())
    hi_sum = sum(t[1] for t in terms.values())
    if lo_sum > b + 1e-6 or hi_sum < b - 1e-6:
        raise InputError("published cells are inconsistent (single-equation infeasible): check identities and links")
    out = {}
    for v, a in u.items():
        tlo, thi = terms[v]
        # range of the other terms
        olo = lo_sum - tlo if tlo != -INF else sum(t[0] for w, t in terms.items() if w != v)
        ohi = hi_sum - thi if thi != INF else sum(t[1] for w, t in terms.items() if w != v)
        r_lo, r_hi = b - ohi, b - olo          # range of a*x_v
        x_lo, x_hi = (r_lo / a, r_hi / a) if a > 0 else (r_hi / a, r_lo / a)
        blo, bhi = box[v]
        out[v] = (max(blo, x_lo), min(bhi, x_hi))
    return out


PROGRESS_EVERY = 250              # progress line every N multi-equation components solved
_T0 = [0.0]
_DEADLINE = [math.inf]


def log(msg: str) -> None:
    """Progress to stderr (counts and timings only, never values)."""
    import time
    print(f"[{time.monotonic() - _T0[0]:8.1f}s] {msg}", file=sys.stderr, flush=True)


def time_left() -> float:
    import time
    return _DEADLINE[0] - time.monotonic()


def check_time() -> None:
    if time_left() <= 0:
        raise TimeLimitError()


class TimeLimitError(Exception):
    """The --max-seconds limit was reached (fails closed: exit 3, no release file)."""


def solve_component(vars_: list[str], kinds: tuple, rows: tuple, p: "Params",
                    seconds_left: float = math.inf) -> tuple[dict[str, tuple[float, float]], int]:
    """Bounds of one multi-equation component (HiGHS). Screening LPs give feasible points; any cell whose observed
    range already passes rule 3 is proven safe (the rule is monotone in the interval). Then batched LPs push all
    still-pending cells down together and up together (one LP each), which proves most of the rest safe; only
    cells still pending after a batch round that proves nothing new get their own exact min / max LPs (and only
    the side not already at its a priori bound). Open upper bounds are capped at 1e9 (treated as unbounded).
    Returns (bounds, number of LPs solved)."""
    import time
    t_end = time.monotonic() + seconds_left
    n = len(vars_)
    r_idx, c_idx, data, rhs = [], [], [], []
    for r, (terms, b) in enumerate(rows):
        for j, a in terms:
            r_idx.append(r)
            c_idx.append(j)
            data.append(a)
        rhs.append(b)
    A = sparse.csr_matrix((data, (r_idx, c_idx)), shape=(len(rows), n))
    b_eq = np.array(rhs, dtype=float)
    bnd = [apriori(k, p) for k in kinds]
    bnd_cap = [(lo, BIG if hi == INF else hi) for lo, hi in bnd]
    n_lp = [0]

    def solve(c: np.ndarray) -> np.ndarray:
        left = t_end - time.monotonic()
        if left <= 0:
            raise TimeLimitError()
        n_lp[0] += 1
        res = linprog(c, A_eq=A, b_eq=b_eq, bounds=bnd_cap, method="highs",
                      options={"time_limit": float(min(left, 1e7))})
        if res.status == 2:
            raise InputError("published cells are inconsistent (LP infeasible): check identities and links")
        if res.status == 1 and (t_end - time.monotonic()) <= 1.0:
            raise TimeLimitError()
        if res.status != 0:
            raise InputError(f"LP solver failed (status {res.status}: {res.message})")
        return res.x

    isP = np.array([k == "P" for k in kinds], dtype=float)
    rng = np.random.default_rng(n + len(rows))
    pts = [solve(c) for c in (np.ones(n), -isP, -np.ones(n), rng.uniform(-1, 1, n), rng.uniform(-1, 1, n))]
    obs_lo = np.min(np.vstack(pts), axis=0)
    obs_hi = np.max(np.vstack(pts), axis=0)

    def pending_of() -> list[int]:
        return [j for j in range(n) if not is_safe(kinds[j], obs_lo[j], obs_hi[j], p)]
    pend = pending_of()
    while pend:                       # batched rounds
        e = np.zeros(n)
        e[pend] = 1.0
        for x in (solve(e), solve(-e)):
            np.minimum(obs_lo, x, out=obs_lo)
            np.maximum(obs_hi, x, out=obs_hi)
        new = pending_of()
        if len(new) == len(pend):
            pend = new
            break
        pend = new
    out: dict[str, tuple[float, float]] = {}
    pset = set(pend)
    for j in range(n):
        if j not in pset:
            out[vars_[j]] = (float(obs_lo[j]), float(obs_hi[j]))
    for j in pend:                    # exact bounds for the remaining cells
        lo0, hi0 = bnd[j]
        e = np.zeros(n)
        e[j] = 1.0
        lo = lo0 if obs_lo[j] <= lo0 + EPS else float(solve(e)[j])
        if hi0 != INF and obs_hi[j] >= hi0 - EPS:
            hi = hi0
        else:
            xm = solve(-e)[j]
            hi = INF if xm >= BIG / 2 else float(xm)
        out[vars_[j]] = (lo, hi)
    return out, n_lp[0]


@dataclass
class AuditResult:
    unsafe: dict[str, str] = field(default_factory=dict)       # cell -> reason ('interval' or 'residual')
    residual_groups: list[int] = field(default_factory=list)   # relation indices failing the residual check
    checked: int = 0                                            # unknown cells whose interval was audited
    exact_recoverable: int = 0                                  # v1.7-style: P or S cells pinned to one value
    kind_fail: dict[str, int] = field(default_factory=dict)    # failing cells by kind (P/S/H)
    comp_of: dict[str, int] = field(default_factory=dict)      # unknown cell -> component id
    comp_size: dict[int, int] = field(default_factory=dict)    # component id -> number of unknown cells


class System:
    """Cells: published (known value), constants (non-count JSON ints used by links), P ('<10'), S, H (hidden).
    Relations: total - sum(parts) = 0, from the published identities and the declared links."""

    def __init__(self, relations: list[tuple[str, tuple[str, ...], str]], values: dict[str, int],
                 kinds: dict[str, str], params: Params):
        self.params = params
        self.values = values                  # published / constant values (input values of published cells)
        self.base_kind = kinds                # cell -> 'pub' | 'const' | 'P' | 'S' | 'H' (input status)
        self.rels: list[list[tuple[str, int]]] = []
        self.rel_total: list[str] = []
        self.rel_source: list[str] = []
        for total, parts, src in relations:
            coef: dict[str, int] = defaultdict(int)
            coef[total] += 1
            for q in parts:
                coef[q] -= 1
            terms = sorted((c, k) for c, k in coef.items() if k != 0)
            if not terms:
                continue
            self.rels.append(terms)
            self.rel_total.append(total)
            self.rel_source.append(src)
        self.by_cell: dict[str, list[int]] = defaultdict(list)
        self.totals_of: dict[str, set[int]] = defaultdict(set)
        for k, terms in enumerate(self.rels):
            for c, _ in terms:
                self.by_cell[c].append(k)
            self.totals_of[self.rel_total[k]].add(k)
        self.cache: dict = {}
        self.lp_solves = 0
        self.workers = 1

    # ---- core audit ------------------------------------------------------------------------------------------
    def audit(self, status, rel_ids: Iterable[int] | None = None, label: str | None = None) -> AuditResult:
        p = self.params
        rel_ids = range(len(self.rels)) if rel_ids is None else sorted(set(rel_ids))
        res = AuditResult()
        live: dict[int, list] = {}
        unknown_cells: set[str] = set()
        for k in rel_ids:
            unk: dict[str, float] = {}
            const = 0.0
            for c, a in self.rels[k]:
                st = status[c]
                if st in UNKNOWN:
                    unk[c] = unk.get(c, 0.0) + a
                else:
                    const += a * self.values[c]
            if not unk:
                if abs(const) > EPS:
                    raise InputError(f"published cells violate an identity/link at {self.rel_total[k]} "
                                     f"({self.rel_source[k]})")
                continue
            unknown_cells.update(unk)
            live[k] = [unk, -const]
            # residual (differencing) check: every unknown member is a suppressed part (P/S), none hidden, and the
            # total is published -> total minus published parts = sum of the suppressed parts; must be >= 10
            if all(status[c] in ("P", "S") for c in unk) and all(a < 0 for a in unk.values()):
                if const < MIN_CELL - EPS:   # sum of the suppressed parts = const
                    res.residual_groups.append(k)
                    for c in unk:
                        res.unsafe.setdefault(c, "residual")
        # components over unknown cells (for 'one fix per component per round')
        parent: dict[str, str] = {}

        def find(a: str) -> str:
            parent.setdefault(a, a)
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        for unk, _ in live.values():
            vs = sorted(unk)
            for v in vs[1:]:
                ra, rb = find(vs[0]), find(v)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
        roots = sorted({find(c) for c in unknown_cells})
        rid = {r: i for i, r in enumerate(roots)}
        res.comp_of = {c: rid[find(c)] for c in unknown_cells}
        sizes: dict[int, int] = defaultdict(int)
        for c, i in res.comp_of.items():
            sizes[i] += 1
        res.comp_size = dict(sizes)
        # 1) exact propagation of single-unknown relations
        rel_list = {k: [dict(u), b] for k, (u, b) in live.items()}
        fixed: dict[str, float] = {}
        by_var: dict[str, list[int]] = defaultdict(list)
        for k, (u, _) in rel_list.items():
            for v in u:
                by_var[v].append(k)
        queue = sorted(k for k, (u, _) in rel_list.items() if len(u) == 1)
        while queue:
            k = queue.pop()
            u, b = rel_list[k]
            if len(u) != 1:
                continue
            (v, a), = u.items()
            val = b / a
            if v in fixed:
                if abs(fixed[v] - val) > 1e-6:
                    raise InputError("published cells are inconsistent (exact propagation)")
                continue
            lo0, hi0 = apriori(status[v], p)
            if val < lo0 - 1e-6 or val > hi0 + 1e-6:
                raise InputError("published cells are inconsistent with the suppression bounds (propagation)")
            fixed[v] = val
            for k2 in by_var[v]:
                u2, b2 = rel_list[k2]
                if v in u2:
                    rel_list[k2][1] = b2 - u2.pop(v) * val
                    if len(u2) == 1:
                        queue.append(k2)
                    elif not u2 and abs(rel_list[k2][1]) > 1e-6:
                        raise InputError("published cells are inconsistent (propagation residue)")
        bounds: dict[str, tuple[float, float]] = {v: (x, x) for v, x in fixed.items()}
        # 2) LP on each remaining component
        rem = [(u, b) for u, b in rel_list.values() if u]
        par2: dict[str, str] = {}

        def find2(a: str) -> str:
            par2.setdefault(a, a)
            while par2[a] != a:
                par2[a] = par2[par2[a]]
                a = par2[a]
            return a
        for u, _ in rem:
            vs = sorted(u)
            for v in vs[1:]:
                ra, rb = find2(vs[0]), find2(v)
                if ra != rb:
                    par2[max(ra, rb)] = min(ra, rb)
        comps: dict[str, list] = defaultdict(list)
        for u, b in rem:
            comps[find2(min(u))].append((u, b))
        bounds.update(self._solve_components([comps[r] for r in sorted(comps)], status, label))
        # 3) cells in no live relation keep their a priori bounds (safe by construction)
        for c in unknown_cells:
            if c not in bounds:
                bounds[c] = apriori(status[c], p)
        res.checked = len(bounds)
        for c in sorted(bounds):
            lo, hi = bounds[c]
            kind = status[c]
            li, hi_i = int_bounds(lo, hi)
            if kind in ("P", "S") and li == hi_i:
                res.exact_recoverable += 1
            if not is_safe(kind, lo, hi, p):
                res.unsafe[c] = "interval+residual" if res.unsafe.get(c) == "residual" else "interval"
        for c in res.unsafe:
            k = status[c]
            res.kind_fail[k] = res.kind_fail.get(k, 0) + 1
        return res

    def _solve_components(self, comps: list[list], status, label: str | None) -> dict[str, tuple[float, float]]:
        """Feasible [min, max] of every variable, component by component (components share no variable, so each
        is solved on its own). One-equation components: exact projection, no LP. Larger components: batched HiGHS
        LPs (solve_component). Results are cached by the component's exact signature (variables, kinds,
        equations, parameters), so a component that did not change between suppression rounds is not re-solved.
        With --workers > 1, uncached multi-equation components are solved in parallel processes; results are
        merged in a fixed order, so the output does not depend on the number of workers."""
        p = self.params
        out: dict[str, tuple[float, float]] = {}
        todo: list[tuple[tuple, dict]] = []
        n_single = n_multi = n_cached = 0
        largest = 0
        for crels in comps:
            check_time()
            if len(crels) == 1:
                u, b = crels[0]
                out.update(single_equation_bounds(u, b, {v: apriori(status[v], p) for v in u}))
                n_single += 1
                continue
            n_multi += 1
            vars_ = sorted({v for u, _ in crels for v in u})
            largest = max(largest, len(vars_))
            kinds = tuple(status[v] for v in vars_)
            col = {v: j for j, v in enumerate(vars_)}
            rows = tuple(sorted((tuple(sorted((col[v], a) for v, a in u.items())), round(b, 6)) for u, b in crels))
            key = (tuple(vars_), kinds, rows, p.width, p.s_rule, p.s_lower)
            hit = self.cache.get(key)
            if hit is not None:
                n_cached += 1
                out.update(hit)
                continue
            todo.append((key, {"vars": vars_, "kinds": kinds, "rows": rows}))
        if label:
            log(f"{label}: {len(comps)} components ({n_single} one-equation, {n_multi} multi-equation, "
                f"{n_cached} cached, {len(todo)} to solve); largest multi-equation component {largest} cells")
        results: list = []
        if self.workers > 1 and len(todo) > 1:
            from concurrent.futures import ProcessPoolExecutor
            order = sorted(range(len(todo)), key=lambda i: -len(todo[i][1]["vars"]))
            with ProcessPoolExecutor(max_workers=self.workers) as ex:
                futs = {i: ex.submit(solve_component, todo[i][1]["vars"], todo[i][1]["kinds"], todo[i][1]["rows"],
                                     p, time_left()) for i in order}
                done = 0
                res_by = {}
                for i in order:
                    res_by[i] = futs[i].result(timeout=max(1.0, time_left()))
                    done += 1
                    if label and (done % PROGRESS_EVERY == 0 or done == len(todo)):
                        log(f"{label}: solved {done}/{len(todo)} multi-equation components")
                results = [res_by[i] for i in range(len(todo))]
        else:
            biggest = 0
            for i, (_key, pl) in enumerate(todo):
                check_time()
                results.append(solve_component(pl["vars"], pl["kinds"], pl["rows"], p, time_left()))
                biggest = max(biggest, len(pl["vars"]))
                if label and ((i + 1) % PROGRESS_EVERY == 0 or i + 1 == len(todo)):
                    log(f"{label}: solved {i + 1}/{len(todo)} multi-equation components "
                        f"(largest so far {biggest} cells)")
        for (key, _pl), (r, n_lp) in zip(todo, results):
            self.lp_solves += n_lp
            self.cache[key] = r
            out.update(r)
        return out

    # ---- exact bounds of a few target cells (v1.2 candidate check) ---------------------------------------------
    def target_bounds(self, status, rel_ids: Iterable[int], targets: list[str]) -> dict[str, tuple[float, float]]:
        """Exact feasible [min, max] of each target cell over the given equations (exact propagation, then HiGHS
        min / max LPs on the target's remaining component). Used only to choose suppressions; every decision is
        re-checked by the full audit and, independently, by the release gate."""
        p = self.params
        rel_list: dict[int, list] = {}
        for k in sorted(set(rel_ids)):
            unk: dict[str, float] = {}
            const = 0.0
            for c, a in self.rels[k]:
                if status[c] in UNKNOWN:
                    unk[c] = unk.get(c, 0.0) + a
                else:
                    const += a * self.values[c]
            if unk:
                rel_list[k] = [unk, -const]
        by_var: dict[str, list[int]] = defaultdict(list)
        for k, (u, _) in rel_list.items():
            for v in u:
                by_var[v].append(k)
        fixed: dict[str, float] = {}
        queue = sorted(k for k, (u, _) in rel_list.items() if len(u) == 1)
        while queue:
            k = queue.pop()
            u, b = rel_list[k]
            if len(u) != 1:
                continue
            (v, a), = u.items()
            if v in fixed:
                continue
            fixed[v] = b / a
            for k2 in by_var[v]:
                u2, b2 = rel_list[k2]
                if v in u2:
                    rel_list[k2][1] = b2 - u2.pop(v) * fixed[v]
                    if len(u2) == 1:
                        queue.append(k2)
        out: dict[str, tuple[float, float]] = {}
        for t in targets:
            if t in fixed:
                out[t] = (fixed[t], fixed[t])
                continue
            # the target's remaining component
            seen_v, seen_k, stack = {t}, set(), [t]
            while stack:
                v = stack.pop()
                for k in by_var.get(v, []):
                    if k in seen_k or not rel_list[k][0]:
                        continue
                    seen_k.add(k)
                    for w in rel_list[k][0]:
                        if w not in seen_v:
                            seen_v.add(w)
                            stack.append(w)
            if not seen_k:
                out[t] = apriori(status[t], p)
                continue
            vars_ = sorted(seen_v)
            col = {v: j for j, v in enumerate(vars_)}
            r_idx, c_idx, data, rhs = [], [], [], []
            for r, k in enumerate(sorted(seen_k)):
                u, b = rel_list[k]
                for v, a in u.items():
                    r_idx.append(r)
                    c_idx.append(col[v])
                    data.append(a)
                rhs.append(b)
            A = sparse.csr_matrix((data, (r_idx, c_idx)), shape=(len(seen_k), len(vars_)))
            bnd = [apriori(status[v], p) for v in vars_]
            bnd_cap = [(lo, BIG if hi == INF else hi) for lo, hi in bnd]
            e = np.zeros(len(vars_))
            e[col[t]] = 1.0
            vals = []
            for sign in (1.0, -1.0):
                check_time()
                self.lp_solves += 1
                r = linprog(sign * e, A_eq=A, b_eq=np.array(rhs, dtype=float), bounds=bnd_cap, method="highs",
                            options={"time_limit": float(min(max(time_left(), 1.0), 1e7))})
                if r.status == 2:
                    raise InputError("published cells are inconsistent (LP infeasible): check identities and links")
                if r.status != 0:
                    if time_left() <= 1.0:
                        raise TimeLimitError()
                    raise InputError(f"LP solver failed (status {r.status}: {r.message})")
                vals.append(float(r.x[col[t]]))
            lo, hi = vals
            out[t] = (lo, INF if hi >= BIG / 2 else hi)
        return out

    def residual_fails(self, status, cell: str) -> bool:
        """The residual (differencing) check of rule 2 on the groups of one cell (same test as in audit)."""
        for k in self.by_cell.get(cell, []):
            unk: dict[str, float] = {}
            const = 0.0
            for c, a in self.rels[k]:
                if status[c] in UNKNOWN:
                    unk[c] = unk.get(c, 0.0) + a
                else:
                    const += a * self.values[c]
            if unk and all(status[c] in ("P", "S") for c in unk) and all(a < 0 for a in unk.values()):
                if const < MIN_CELL - EPS:
                    return True
        return False

    # ---- local sub-system (component reachable through unknown cells) -----------------------------------------
    def component_rels(self, seeds: Iterable[str], status) -> set[int]:
        seen_c, seen_k = set(), set()
        stack = [s for s in seeds]
        while stack:
            c = stack.pop()
            if c in seen_c:
                continue
            seen_c.add(c)
            for k in self.by_cell.get(c, []):
                if k in seen_k:
                    continue
                seen_k.add(k)
                for m, _ in self.rels[k]:
                    if m not in seen_c and status[m] in UNKNOWN:
                        stack.append(m)
        return seen_k


class Overlay(dict):
    """status view: overrides on top of a base mapping (cheap copies for look-ahead)."""

    def __init__(self, base, over=None):
        super().__init__(over or {})
        self.base = base

    def __missing__(self, key):
        return self.base[key]


# =============================================================================
# Rule 2: complementary suppression with look-ahead
# =============================================================================
class _TooBig(Exception):
    pass


class Suppressor:
    def __init__(self, system: System, status: dict, frozen: set[str]):
        self.sys = system
        self.status = status
        self.frozen = frozen
        self.rule_of: dict[str, str] = {}
        self.rounds = 0
        self.lookaheads = 0
        self.exact_checks = 0
        self.local_cache: dict = {}

    def pickable(self, m: str, st) -> bool:
        return st[m] == "pub" and m not in self.frozen

    def candidates(self, u: str, st) -> list[str]:
        S = self.sys
        cands = sorted({m for k in S.by_cell.get(u, []) for m, _ in S.rels[k] if self.pickable(m, st)})
        if cands:
            return cands
        # breadth-first through relations linked by unknown cells
        seen_c, seen_k, frontier = {u}, set(), [u]
        while frontier:
            found, nxt = set(), []
            for c in frontier:
                for k in S.by_cell.get(c, []):
                    if k in seen_k:
                        continue
                    seen_k.add(k)
                    for m, _ in S.rels[k]:
                        if self.pickable(m, st):
                            found.add(m)
                        elif m not in seen_c and st[m] in UNKNOWN:
                            seen_c.add(m)
                            nxt.append(m)
            if found:
                return sorted(found)
            frontier = sorted(nxt)
        return []

    def static_key(self, m: str, u: str) -> tuple:
        S = self.sys
        is_total = bool(S.totals_of.get(m, set()).intersection(S.by_cell.get(u, [])))
        return (len(S.by_cell.get(m, [])), 0 if is_total else 1, S.values[m], m)

    def role(self, m: str, u: str) -> str:
        S = self.sys
        ks = sorted(set(S.by_cell.get(u, [])).intersection(S.by_cell.get(m, [])))
        if not ks:
            return "complementary_linked"
        if any(S.rel_source[k] == "link" for k in ks) and not any(S.rel_source[k] == "identity" for k in ks):
            return "complementary_link"
        return "complementary_total" if any(S.rel_total[k] == m for k in ks) else "complementary_part"

    def local_unsafe(self, seeds: list[str], st) -> dict[str, str]:
        rels = self.sys.component_rels(seeds, st)
        if not rels:
            return {}
        if len(rels) > LOOKAHEAD_MAX_RELS:
            raise _TooBig()
        key = tuple(sorted(rels))
        cells = sorted({c for k in key for c, _ in self.sys.rels[k]})
        key = (key, tuple(st[c] for c in cells))
        hit = self.local_cache.get(key)
        if hit is None:
            hit = self.sys.audit(st, rels).unsafe
            self.local_cache[key] = hit
        return hit

    def simulate(self, u: str, m: str) -> float:
        """Number of suppressions needed (m included) until u's component passes, choosing statically after m."""
        self.lookaheads += 1
        try:
            return self._simulate(u, m)
        except _TooBig:
            return INF          # the candidate merges u into a very large component: never preferred by look-ahead

    def _simulate(self, u: str, m: str) -> float:
        st = Overlay(self.status, {m: "S"})
        steps = 1
        seeds = [u, m]
        for _ in range(LOOKAHEAD_CAP):
            uns = self.local_unsafe(seeds, st)
            if not uns:
                return steps
            w = sorted(uns)[0]
            cs = self.candidates(w, st)
            if not cs:
                return INF
            c = min(cs, key=lambda x: self.static_key(x, w))
            st[c] = "S"
            seeds.append(c)
            steps += 1
        return INF

    def choose(self, u: str, lookahead: bool = True) -> str | None:
        cs = self.candidates(u, self.status)
        if not cs:
            return None
        cs = sorted(cs, key=lambda x: self.static_key(x, u))[:LOOKAHEAD_K]
        if len(cs) == 1 or not lookahead:
            return cs[0]
        scored = sorted((self.simulate(u, m) if len(self.sys.by_cell.get(m, [])) <= LOOKAHEAD_MAX_LINKS else INF,
                         self.static_key(m, u), m) for m in cs)
        return scored[0][2]

    # ---- v1.2: exact candidate check for large components ---------------------------------------------------
    def ring(self, u: str, depth: int) -> list[tuple[int, str]]:
        """Publishable, non-frozen cells within `depth` groups of u, with their distance: depth 1 = cells in u's
        own groups; depth 2 = cells in the groups of the unknown cells that share a group with u; and so on."""
        S, st = self.sys, self.status
        out: dict[str, int] = {}
        seen_c, seen_k, frontier = {u}, set(), [u]
        for d in range(1, depth + 1):
            nxt = []
            for c in frontier:
                for k in S.by_cell.get(c, []):
                    if k in seen_k:
                        continue
                    seen_k.add(k)
                    for m, _ in S.rels[k]:
                        if self.pickable(m, st):
                            out.setdefault(m, d)
                        elif m not in seen_c and st[m] in UNKNOWN:
                            seen_c.add(m)
                            nxt.append(m)
            frontier = sorted(nxt)
        return sorted((d, m) for m, d in out.items())

    def makes_safe(self, u: str, adds: list[str]) -> bool:
        """Exact check: with `adds` suppressed, does u meet rule 3 and the residual check?"""
        self.exact_checks += 1
        st = Overlay(self.status, {m: "S" for m in adds})
        rels = self.sys.component_rels([u, *adds], st)
        lo, hi = self.sys.target_bounds(st, rels, [u])[u]
        return is_safe(st[u], lo, hi, self.sys.params) and not self.sys.residual_fails(st, u)

    def unprotectable(self, u: str) -> bool:
        """Proof that no suppression can protect u: free EVERY publishable, non-frozen cell in u's linked
        equations (dropping all other equations those cells appear in, which only widens intervals). If u still
        fails, any real suppression pattern (a subset of this, with more equations) fails too."""
        self.exact_checks += 1
        rels = self.sys.component_rels([u], self.status)
        free = {m for k in rels for m, _ in self.sys.rels[k] if self.pickable(m, self.status)}
        st = Overlay(self.status, {m: "S" for m in free})
        lo, hi = self.sys.target_bounds(st, rels, [u])[u]
        return not is_safe(st[u], lo, hi, self.sys.params) or self.sys.residual_fails(st, u)

    def choose_exact(self, u: str) -> list[str] | None:
        """Large components (v1.2). Every returned suppression set is checked exactly to make u pass:
        1) single candidates from u's own groups (the static order first, then the largest values);
        2) otherwise the smallest prefix (u's own groups first, larger values first) of all candidates within
           1, 2, ... EXACT_RINGS groups that makes u pass (adding suppressions only widens intervals, so the
           prefix search is exact). None if even all of them together do not."""
        r1 = [m for d, m in self.ring(u, 1)]
        if r1:
            order = sorted(r1, key=lambda x: self.static_key(x, u))
            singles = order[:EXACT_SINGLES]
            singles += [m for m in sorted(r1, key=lambda x: (-self.sys.values[x], x))[:EXACT_TOP_VALUE]
                        if m not in singles]
            for m in singles:
                if self.makes_safe(u, [m]):
                    return [m]
        for depth in range(1, EXACT_RINGS + 1):
            cand = self.ring(u, depth)
            if not cand:
                continue
            ordered = [m for d, v, m in sorted((d, -self.sys.values[m], m) for d, m in cand)]
            if not self.makes_safe(u, ordered):
                continue
            lo, hi = 1, len(ordered)
            while lo < hi:
                mid = (lo + hi) // 2
                if self.makes_safe(u, ordered[:mid]):
                    hi = mid
                else:
                    lo = mid + 1
            return ordered[:lo]
        return None

    def run(self) -> AuditResult:
        """Rounds of: audit everything (unchanged components come from the cache), then new suppressions.
        - Components up to LARGE_COMPONENT cells: one new suppression per failing component, chosen by
          look-ahead (unchanged from v1.0).
        - Larger components: one fix per failing GROUP, chosen by an exact LP check that it makes the failing cell
          pass (choose_exact); a cell that provably cannot be protected stops the run (fail closed).
        - No-progress guard: if the failing count has not gone below its best for STALL_ROUNDS rounds, every
          failing cell gets the stronger rule (all candidates within EXACT_RINGS groups suppressed together);
          if it still has not improved after another STALL_ROUNDS rounds, the run fails closed."""
        best = INF
        stall = 0
        for r in range(1, MAX_ROUNDS + 1):
            check_time()
            res = self.sys.audit(self.status, label=f"suppression round {r}: audit")
            if not res.unsafe:
                log(f"suppression round {r}: no failing cell; done")
                return res
            n_fail = len(res.unsafe)
            if n_fail < best:
                trend = "first round" if best == INF else f"improved from {int(best)}"
                best, stall = n_fail, 0
            else:
                stall += 1
                trend = f"NOT improved for {stall} round(s) (best {int(best)})"
            if stall >= 2 * STALL_ROUNDS:
                bad = sorted(res.unsafe)
                raise GateError(f"suppression not converging: {n_fail} cell(s) still failing and the failing count "
                                f"has not decreased in {stall} rounds, {STALL_ROUNDS} of them with the stronger "
                                f"rule (first pointers: {bad[:5]})")
            escalate = stall >= STALL_ROUNDS
            log(f"suppression round {r}: {n_fail} failing cells, {trend}"
                + ("; ESCALATING to the stronger rule" if escalate else ""))
            self.rounds += 1
            touched: set[int] = set()
            touched_rel: set[int] = set()
            progressed = 0
            unprot: list[str] = []
            for u in sorted(res.unsafe):
                check_time()
                if self.status[u] not in UNKNOWN:
                    continue
                comp = res.comp_of.get(u, -1)
                large = res.comp_size.get(comp, 0) > LARGE_COMPONENT
                if escalate:
                    picks = [m for d, m in self.ring(u, EXACT_RINGS)]
                    if not picks and self.unprotectable(u):
                        unprot.append(u)
                elif large:
                    if touched_rel.intersection(self.sys.by_cell.get(u, [])):
                        continue
                    picks = self.choose_exact(u)
                    if picks is None:
                        if self.unprotectable(u):
                            unprot.append(u)
                            continue
                        # protectable, but only from further out: the stronger rule now (all candidates within
                        # EXACT_RINGS groups), else the nearest candidate; the no-progress guard bounds this
                        picks = [m for d, m in self.ring(u, EXACT_RINGS)]
                        if not picks:
                            m = self.choose(u, lookahead=False)
                            picks = [] if m is None else [m]
                else:
                    if comp in touched:
                        continue
                    m = self.choose(u, lookahead=True)
                    picks = [] if m is None else [m]
                for m in picks:
                    if self.status[m] != "pub":
                        continue
                    self.rule_of[m] = self.role(m, u)
                    self.status[m] = "S"
                    progressed += 1
                    if large or escalate:
                        touched_rel.update(self.sys.by_cell.get(m, []))
                if picks:
                    if large or escalate:
                        touched_rel.update(self.sys.by_cell.get(u, []))
                    else:
                        touched.add(comp)
            log(f"suppression round {r}: {progressed} new suppressions (total new {len(self.rule_of)}; look-ahead "
                f"simulations {self.lookaheads}; exact candidate checks {self.exact_checks})")
            if unprot:
                raise GateError(f"cannot protect {len(unprot)} cell(s): even with every publishable, non-frozen cell "
                                f"in their linked equations suppressed, the protection rule cannot be met; they are "
                                f"fixed by frozen cells, links or constants (first pointers: {sorted(unprot)[:5]})")
            if not progressed and escalate:
                bad = sorted(res.unsafe)
                raise GateError(f"suppression not converging: {n_fail} cell(s) still failing after the stronger rule "
                                f"suppressed every candidate within {EXACT_RINGS} groups (first pointers: {bad[:5]})")
            if not progressed:
                bad = sorted(res.unsafe)
                raise GateError(f"cannot protect {len(bad)} cell(s): no publishable, non-frozen cell left to "
                                f"suppress in their groups (first pointers: {bad[:5]})")
        raise GateError("suppression not converging: round limit reached")


# =============================================================================
# Loading the released output
# =============================================================================
@dataclass
class Loaded:
    doc: dict
    md: str | None
    count_fields: list[str]
    hidden: list[str]
    identities: list[tuple[str, tuple[str, ...]]]
    ratios: list[tuple[str, str, str]]
    quantiles: list[tuple[str, str]]
    dependents: list[tuple[str, tuple[str, ...]]]
    status: dict[str, str]          # count cells and hidden cells
    values: dict[str, int]          # published count values (+ link constants)
    primary_new: list[str]
    dropped: list[str] = field(default_factory=list)              # dropped count-cell pointers (status D)
    drop_subtrees: list[tuple[str, str]] = field(default_factory=list)  # (subtree pointer, reason)
    out_disclosure: dict | None = None    # registry carried into the output (input registry minus the drops)


def check_name(p: Path, role: str) -> None:
    if "_raw_local" in p.name.lower():
        raise InputError(f"{role} '{p.name}' is a _raw_LOCAL file: never used as input or output")


def load_output(json_path: Path, md_path: Path | None, profile: str, expect_spec: str | None,
                allow_v102: bool = False) -> Loaded:
    check_name(json_path, "input JSON")
    if md_path is not None:
        check_name(md_path, "input markdown")
    try:
        doc = json.loads(json_path.read_text(encoding="utf-8"))
    except Exception as e:  # noqa: BLE001
        raise InputError(f"cannot read the input JSON: {e}")
    md = md_path.read_text(encoding="utf-8") if md_path is not None else None
    if not isinstance(doc, dict):
        raise InputError("input JSON is not an object")
    if doc.get("release_gate") != "ok":
        raise InputError("input JSON is not a released output (release_gate != 'ok')")
    if expect_spec is not None and doc.get("spec") != expect_spec:
        raise InputError(f"input spec string is not the expected one ({expect_spec!r})")
    if PROVENANCE_KEY in doc and not allow_v102:
        raise InputError("input already carries a v1.0.2 provenance block (it is not a v1.0.1 output)")
    dis = doc.get("disclosure")
    if not isinstance(dis, dict):
        raise InputError("input JSON has no disclosure block")
    for k in ("count_fields", "identities", "hidden_cells"):
        if k not in dis:
            raise InputError(f"disclosure block has no '{k}'")
    declared = list(dis["count_fields"])
    detected = detect_count_fields(doc)
    if set(declared) != set(detected) or len(declared) != len(set(declared)):
        raise InputError(f"count-field schema mismatch: {len(set(detected) - set(declared))} undeclared, "
                         f"{len(set(declared) - set(detected))} missing")
    status: dict[str, str] = {}
    values: dict[str, int] = {}
    primary_new: list[str] = []
    for ptr in declared:
        v = get_ptr(doc, ptr)
        if isinstance(v, bool):
            raise InputError(f"count field holds a boolean at {ptr}")
        if isinstance(v, int):
            if v < 0:
                raise InputError(f"negative count at {ptr}")
            if v < MIN_CELL:          # rule 1 safety net: should not occur in a released output
                status[ptr] = "P"
                primary_new.append(ptr)
            else:
                status[ptr] = "pub"
            values[ptr] = v
        elif v == "<10":
            status[ptr] = "P"
        elif v == "S":
            status[ptr] = "S"
        else:
            raise InputError(f"count field with an invalid value at {ptr}")
    hidden = list(dis["hidden_cells"])
    for h in hidden:
        if not str(h).startswith("hidden:"):
            raise InputError("hidden cell without the hidden: prefix")
        status[h] = "H"
    identities = [(t, tuple(p)) for t, p in dis["identities"]]
    if not identities:
        raise InputError("no disclosure identities published")
    for t, parts in identities:
        for c in (t,) + parts:
            if c not in status:
                raise InputError(f"identity references an undeclared cell: {c}")
    ratios = [(a, b, c) for a, b, c in dis.get("ratios", [])]
    quantiles = [(a, b) for a, b in dis.get("quantiles", [])]
    dependents = [(a, tuple(b)) for a, b in dis.get("dependents", [])]
    if profile == "stage2" and md is None:
        raise InputError("--in-md is required for the stage2 profile")
    return Loaded(doc, md, declared, hidden, identities, ratios, quantiles, dependents, status, values, primary_new)


def resolve_links(L: Loaded, links: list[tuple[str, tuple[str, ...], str]]) -> list[tuple[str, tuple[str, ...], str]]:
    rels = []
    for t, parts, _note in links:
        for c in (t,) + parts:
            if c in L.status:
                continue
            try:
                v = get_ptr(L.doc, c)
            except KeyError:
                raise InputError(f"links file pointer not found in the input JSON: {c}")
            if isinstance(v, bool) or not isinstance(v, int):
                raise InputError(f"links file pointer is neither a count field nor an integer: {c}")
            L.status[c] = "const"
            L.values[c] = v
        rels.append((t, parts, "link"))
    return rels


def build_system(L: Loaded, link_rels, params: Params) -> System:
    rels = [(t, p, "identity") for t, p in L.identities] + list(link_rels)
    return System(rels, L.values, dict(L.status), params)


# =============================================================================
# Applying the new statuses to the document (counts and derived values)
# =============================================================================
def apply_statuses(L: Loaded, status: dict) -> dict:
    doc = json.loads(json.dumps(L.doc))  # deep copy, key order kept
    pub = set()
    subs = [t for t, _ in L.drop_subtrees]

    def gone(x: str) -> bool:
        return any(_under(x, t) for t in subs)
    for ptr in L.count_fields:
        st = status[ptr]
        if st == "pub":
            pub.add(ptr)
            continue
        if st == "D":
            continue
        set_ptr(doc, ptr, "<10" if st == "P" else "S")
    derived = {"pct": 0, "q": 0, "dependent": 0}
    for ptr, num, den in L.ratios:
        if gone(ptr):
            continue
        cur = get_ptr(doc, ptr)
        if cur != "suppressed" and not (num in pub and den in pub and L.values[den] > 0):
            set_ptr(doc, ptr, "suppressed")
            derived["pct"] += 1
    for ptr, n in L.quantiles:
        if gone(ptr):
            continue
        cur = get_ptr(doc, ptr)
        if cur != "suppressed" and not (n in pub and L.values[n] >= MIN_CELL):
            set_ptr(doc, ptr, "suppressed")
            derived["q"] += 1
    for ptr, needs in L.dependents:
        if gone(ptr):
            continue
        cur = get_ptr(doc, ptr)
        if cur != "suppressed" and not all(x in pub for x in needs):
            set_ptr(doc, ptr, "suppressed")
            derived["dependent"] += 1
    for t in subs:                     # remove each dropped subtree as a whole; the registry loses its pointers
        parts = ptr_parts(t)
        cur = doc
        for q in parts[:-1]:
            cur = cur[q]
        del cur[parts[-1]]
    if subs:
        doc["disclosure"] = json.loads(json.dumps(L.out_disclosure))
    doc["_derived_counts"] = derived  # removed again in finish_document
    return doc


def finish_document(doc: dict, status: dict, L: Loaded, res: AuditResult, prov: dict,
                    profile: str) -> tuple[dict, dict]:
    derived = doc.pop("_derived_counts")
    n_p = sum(1 for c in L.count_fields if status[c] == "P")
    n_s = sum(1 for c in L.count_fields if status[c] == "S")
    sup = dict(doc.get("suppression", {}))
    sup.update({"primary_cells": n_p, "complementary_cells": n_s, "recoverable": 0,
                "recoverability_cells_checked": res.checked})
    doc["suppression"] = sup
    if profile == "stage2":
        doc["spec"] = OUTPUT_SPEC
        notes = doc.setdefault("notes", {})
        notes["differencing_fix_v1_0_2"] = (
            "v1.0.2 re-releases the v1.0.1 output with complementary suppression closed against differencing "
            "(totals minus published parts, within and across tables); every published value is unchanged from "
            "v1.0.1; cells newly marked S and derived values newly marked suppressed are listed by pointer in the "
            "release summary")
    if L.drop_subtrees:
        notes = doc.setdefault("notes", {})
        notes["dropped_subtrees_v1_0_2"] = (
            f"{len(L.drop_subtrees)} subtree(s) withheld from the v1.0.2 release as a whole (hash-pinned drop list; "
            f"the totals above them are kept): " + "; ".join(t for t, _ in L.drop_subtrees))
    prov = dict(prov)
    prov["derived_values_newly_suppressed"] = derived
    doc[PROVENANCE_KEY] = prov
    doc["release_gate"] = "ok"
    return doc, derived


# =============================================================================
# Release gate (run on files reloaded from disk)
# =============================================================================
def gate_document(doc: dict, md: str | None, L: Loaded, system: System, frozen: set[str], profile: str,
                  params: Params, render_ok_input: bool) -> tuple[list[str], AuditResult | None]:
    fails: list[str] = []
    dis = doc.get("disclosure", {})
    if dis != L.out_disclosure:
        fails.append("disclosure block differs from the input (it must be carried unchanged, minus dropped pointers)")
    for t, _ in L.drop_subtrees:
        try:
            get_ptr(doc, t)
            fails.append(f"dropped subtree still present in the output at {t}")
        except KeyError:
            pass
    declared = list(dis.get("count_fields", []))
    dset = set(L.dropped)
    if any(c in dset for c in declared):
        fails.append("a dropped cell is still declared in the output")
    if set(declared) != set(detect_count_fields(doc)):
        fails.append("count-field schema mismatch in the output")
    status: dict[str, str] = {}
    n_p = n_s = 0
    for ptr in declared:
        try:
            v = get_ptr(doc, ptr)
        except KeyError:
            fails.append(f"declared count field missing at {ptr}")
            continue
        old = L.status[ptr]
        if isinstance(v, bool) or not (isinstance(v, int) or v in ("<10", "S")):
            fails.append(f"invalid count value at {ptr}")
            continue
        if isinstance(v, int):
            if v < MIN_CELL:
                fails.append(f"published count below 10 at {ptr}")
            if old != "pub" or v != L.values.get(ptr):
                fails.append(f"published value differs from v1.0.1 or was suppressed there, at {ptr}")
            status[ptr] = "pub"
        elif v == "<10":
            if old != "P":
                fails.append(f"'<10' at a cell that was not '<10' in the input, at {ptr}")
            status[ptr] = "P"
            n_p += 1
        else:
            if old == "P":
                fails.append(f"'<10' cell relabelled S at {ptr}")
            status[ptr] = "S"
            n_s += 1
    for h in L.hidden:
        status[h] = "H"
    for c in L.dropped:
        status[c] = "D"
    for c, st in L.status.items():
        if st == "const":
            status[c] = "const"
            try:
                if get_ptr(doc, c) != L.values[c]:
                    fails.append(f"link constant changed at {c}")
            except KeyError:
                fails.append(f"link constant missing at {c}")
    for c in sorted(frozen):
        if status.get(c) != "pub":
            fails.append(f"frozen (already public) cell is not published at {c}")
    sup = doc.get("suppression", {})
    if sup.get("primary_cells") != n_p or sup.get("complementary_cells") != n_s or sup.get("recoverable") != 0:
        fails.append("suppression counters do not match the published cells")
    pubv = {c: L.values[c] for c, st in status.items() if st == "pub"}
    subs = [t for t, _ in L.drop_subtrees]

    def gone(x: str) -> bool:
        return any(_under(x, t) for t in subs)
    for ptr, num, den in L.ratios:
        if gone(ptr):
            continue
        v = get_ptr(doc, ptr)
        if v == "suppressed":
            continue
        if num not in pubv or den not in pubv or pubv[den] == 0 or v != r1(100.0 * pubv[num] / pubv[den]):
            fails.append(f"percentage not reproducible from published counts at {ptr}")
    for ptr in detect_prefixed(doc, ("pct_",)):
        if ptr not in {r[0] for r in L.ratios}:
            fails.append(f"undeclared percentage at {ptr}")
            break
    for ptr, n in L.quantiles:
        if gone(ptr):
            continue
        v = get_ptr(doc, ptr)
        if v == "suppressed":
            continue
        if not isinstance(v, dict) or set(v) != set(QKEYS) or n not in pubv or pubv[n] < MIN_CELL:
            fails.append(f"quantile published without its own published count >= 10 at {ptr}")
    for ptr, needs in L.dependents:
        if gone(ptr):
            continue
        v = get_ptr(doc, ptr)
        if v != "suppressed" and not all(x in pubv for x in needs):
            fails.append(f"dependent value published although a count it rests on is suppressed at {ptr}")

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
    if profile == "stage2":
        if not render_ok_input:
            fails.append("renderer drift: the bundled renderer does not reproduce the v1.0.1 markdown from the "
                         "v1.0.1 JSON")
        if md is None or render_markdown(doc) != md:
            fails.append("markdown on disk is not the rendering of the v1.0.2 JSON")
    res = None
    if not fails:
        system.cache.clear()
        res = system.audit(status, label="release gate audit (independent re-solve)")
        if res.unsafe:
            fails.append(f"interval audit: {len(res.unsafe)} suppressed cell(s) fail the protection rule "
                         f"(by kind: {dict(sorted(res.kind_fail.items()))}); residual check: "
                         f"{len(res.residual_groups)} group(s) whose total minus published parts is under 10")
    return fails, res


# =============================================================================
# Reports
# =============================================================================
def summary_markdown(info: dict) -> str:
    L = [f"# Stage 2 v1.0.2 release summary ({info['mode']})", "",
         "Counts of cells only. No count values appear in this file; newly suppressed cells are listed by JSON "
         "pointer.", "",
         f"- script: {SCRIPT_VERSION} (sha256 {info['script_sha256']})",
         f"- spec: {SPEC_DOC}",
         f"- input spec: {info['input_spec']}",
         f"- result: {info['result']}", "",
         "## Parameters", "",
         f"- protection width W: {info['params']['width']} (integers in a suppressed cell's feasible interval; "
         f"'<10' cells must keep 0-9)",
         f"- S-cell rule: {info['params']['s_rule']}; S a priori lower bound: {info['params']['s_lower']}",
         f"- links used: {info['n_links']} (file sha256 {info['links_sha256']})",
         f"- frozen cells: {info['n_frozen']} (file sha256 {info['frozen_sha256']})", ""]
    if info.get("drops"):
        d = info["drops"]
        L[-1:] = [f"- dropped subtrees: {len(d['subtrees'])} holding {d['n_count_cells']} count cells, withheld "
                  f"from the release as a whole (file sha256 {d['sha256']}; pointers below, no values)", ""]
        L += [f"  - {t}" for t in d["subtrees"]] + [""]
    L += [
         "## Input audit (the released v1.0.1 output as it stands)", "",
         f"- count cells: {info['n_count_cells']}; hidden differences: {info['n_hidden']}; identities: "
         f"{info['n_identities']}",
         f"- suppressed cells audited: {info['input_audit']['checked']}",
         f"- cells failing the protection rule: {info['input_audit']['fail_total']} "
         f"(P {info['input_audit']['fail'].get('P', 0)}, S {info['input_audit']['fail'].get('S', 0)}, "
         f"H {info['input_audit']['fail'].get('H', 0)}" +
         (f", D {info['input_audit']['fail'].get('D', 0)}" if info.get("drops") else "") + ")",
         f"- groups failing the residual (differencing) check: {info['input_audit']['residual_groups']}",
         f"- cells exactly recoverable (the v1.0.1 gate's own criterion): {info['input_audit']['exact']}", ""]
    if info["mode"] == "check-only":
        L += ["## Failing cells (pointers only)", ""]
        L += [f"- {p} ({r})" for p, r in info.get("failing", [])] or ["- none"]
        return "\n".join(L) + "\n"
    r = info["rules"]
    L += ["## Cells by rule (v1.0.2)", "",
          "| rule | cells |", "|---|---|",
          f"| 1 primary, carried from v1.0.1 ('<10') | {r['primary_carried']} |",
          f"| 1 primary, new (published count under 10 found in the input) | {r['primary_new']} |",
          f"| 2 complementary, carried from v1.0.1 (S) | {r['complementary_carried']} |",
          f"| 2 complementary, new: group total suppressed | {r['complementary_total']} |",
          f"| 2 complementary, new: part suppressed | {r['complementary_part']} |",
          f"| 2 complementary, new: cross-table link member suppressed | {r['complementary_link']} |",
          f"| 2 complementary, new: linked group (reached through other suppressed cells) | {r['complementary_linked']} |",
          f"| derived percentages newly suppressed | {r['derived_pct']} |",
          f"| derived quantile sets newly suppressed | {r['derived_q']} |",
          f"| derived KS / dependent values newly suppressed | {r['derived_dependent']} |", "",
          f"- suppression rounds: {info['rounds']}; look-ahead simulations: {info['lookaheads']}; LP solves: "
          f"{info['lp_solves']}",
          f"- output audit: {info['output_audit']['checked']} suppressed cells audited; failing "
          f"{info['output_audit']['fail_total']}; residual failures {info['output_audit']['residual_groups']}",
          f"- output counters: primary {info['out_primary']}, complementary {info['out_complementary']}", "",
          "## Newly suppressed cells (pointers only)", ""]
    L += [f"- {p} ({rule})" for p, rule in info["new_cells"]] or ["- none"]
    return "\n".join(L) + "\n"


# =============================================================================
# Orchestration
# =============================================================================
def run(args: argparse.Namespace) -> int:
    in_json = Path(args.in_json).resolve()
    in_md = Path(args.in_md).resolve() if args.in_md else None
    links_p = Path(args.links).resolve() if args.links else None
    frozen_p = Path(args.frozen).resolve() if args.frozen else None
    drops_p = Path(args.drops).resolve() if args.drops else None
    inputs = [("input_json", in_json)] + ([("input_md", in_md)] if in_md else []) + \
             ([("links", links_p)] if links_p else []) + ([("frozen", frozen_p)] if frozen_p else []) + \
             ([("drops", drops_p)] if drops_p else [])
    for role, pth in inputs:
        check_name(pth, role)
        if not pth.is_file():
            raise InputError(f"{role} not found: {pth}")
    in_hash = {role: sha256_file(pth) for role, pth in inputs}
    if args.expect_input_sha256 and in_hash["input_json"] != args.expect_input_sha256.lower():
        raise InputError("input JSON SHA-256 does not match --expect-input-sha256")
    if args.expect_md_sha256 and in_hash.get("input_md") != args.expect_md_sha256.lower():
        raise InputError("input markdown SHA-256 does not match --expect-md-sha256")
    params = Params(width=args.protection_width, s_rule=args.s_rule, s_lower=args.s_lower_bound)
    if params.width < 1:
        raise InputError("--protection-width must be >= 1")
    expect_spec = None if args.profile == "generic" else (args.expect_spec or EXPECT_INPUT_SPEC)
    if args.check_only:
        expect_spec = None  # check-only audits any registry-style output (v1.0.1, v1.0.2, follow-up)
    log(f"{SCRIPT_VERSION}: reading the input JSON ({in_json.stat().st_size:,} bytes) and checking its schema")
    L = load_output(in_json, in_md, args.profile, expect_spec, allow_v102=args.check_only)
    log(f"input loaded: {len(L.count_fields)} count cells, {len(L.hidden)} hidden differences, "
        f"{len(L.identities)} identities")
    links = read_links(links_p)
    frozen_list = read_frozen(frozen_p)
    link_rels = resolve_links(L, links)
    frozen = set()
    for ptr, _reason in frozen_list:
        if L.status.get(ptr) != "pub":
            raise InputError(f"frozen pointer is not a published count cell in the input: {ptr}")
        frozen.add(ptr)
    drops = read_drops(drops_p)
    apply_drops(L, drops, {c for t, parts, _ in link_rels for c in (t,) + parts}, frozen)
    if drops:
        log(f"drop list: {len(drops)} subtree(s) withheld from the release, {len(L.dropped)} count cells marked "
            f"unpublished (audited as hidden)")
    system = build_system(L, link_rels, params)
    system.workers = args.workers
    log(f"equation system built: {len(system.rels)} equations ({len(link_rels)} links), {len(frozen)} frozen cells")
    base_status = dict(L.status)
    in_res = system.audit(base_status, label="input audit")
    log(f"input audit done: {in_res.checked} suppressed cells audited, {len(in_res.unsafe)} failing, "
        f"{len(in_res.residual_groups)} groups failing the residual check")
    info: dict[str, Any] = {
        "mode": "check-only" if args.check_only else "release",
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "input_spec": L.doc.get("spec"),
        "params": {"width": params.width, "s_rule": params.s_rule, "s_lower": params.s_lower},
        "n_links": len(link_rels), "links_sha256": in_hash.get("links", "none"),
        "n_frozen": len(frozen), "frozen_sha256": in_hash.get("frozen", "none"),
        "n_count_cells": len(L.count_fields), "n_hidden": len(L.hidden), "n_identities": len(L.identities),
        **({"drops": {"sha256": in_hash["drops"], "subtrees": [t for t, _ in drops], "n_count_cells": len(L.dropped)}}
           if drops else {}),
        "input_audit": {"checked": in_res.checked, "fail_total": len(in_res.unsafe), "fail": in_res.kind_fail,
                        "residual_groups": len(in_res.residual_groups), "exact": in_res.exact_recoverable},
    }
    if args.check_only:
        failed = bool(in_res.unsafe) or bool(L.primary_new)
        info["result"] = "GATE FAIL (check-only)" if failed else "GATE OK (check-only)"
        info["failing"] = sorted(in_res.unsafe.items())
        print(summary_markdown(info), end="")
        if L.primary_new:
            print(f"published counts under 10: {len(L.primary_new)}")
        return 3 if failed else 0

    out_dir = Path(args.out_dir).resolve()
    created_dir = not out_dir.exists()
    in_dirs = {p.parent for _, p in inputs}
    if out_dir in in_dirs:
        raise InputError("--out-dir must be a new folder, not a folder holding an input")
    if out_dir.exists() and any(out_dir.iterdir()):
        raise InputError("--out-dir exists and is not empty")
    check_name(out_dir, "output folder")
    render_ok_input = True
    if args.profile == "stage2":
        render_ok_input = render_markdown(L.doc) == L.md
        if not render_ok_input:
            raise GateError("renderer drift: the bundled renderer does not reproduce the v1.0.1 markdown from the "
                            "v1.0.1 JSON; nothing written")

    status = dict(base_status)
    sup = Suppressor(system, status, frozen)
    log("complementary suppression: start")
    sup.run()
    log(f"complementary suppression: done ({len(sup.rule_of)} new cells, {sup.rounds} rounds)")
    doc = apply_statuses(L, status)
    new_cells = sorted((c, sup.rule_of[c]) for c in sup.rule_of)
    new_cells += sorted((c, "primary_new") for c in L.primary_new)
    prov = {"source_json": in_json.name, "source_json_sha256": in_hash["input_json"],
            "source_md_sha256": in_hash.get("input_md", "none"), "script": SCRIPT_VERSION,
            "script_sha256": info["script_sha256"], "spec": SPEC_DOC,
            "links_sha256": in_hash.get("links", "none"), "frozen_sha256": in_hash.get("frozen", "none"),
            "protection_width": params.width, "s_rule": params.s_rule, "s_lower_bound": params.s_lower,
            "newly_suppressed_cells": len(new_cells)}
    if drops:
        prov["drops_sha256"] = in_hash["drops"]
        prov["dropped_subtrees"] = [t for t, _ in drops]
        prov["dropped_count_cells"] = len(L.dropped)
    pre_res = system.audit(status)
    doc, derived = finish_document(doc, status, L, pre_res, prov, args.profile)

    written: list[Path] = []
    out_json, out_md = out_dir / OUT_JSON, out_dir / OUT_MD
    pend_json, pend_md = out_dir / (OUT_JSON + ".pending"), out_dir / (OUT_MD + ".pending")
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        written += [pend_json]
        write_text(pend_json, dump_json(doc))
        if args.profile == "stage2":
            written += [pend_md]
            write_text(pend_md, render_markdown(doc))
        # gate on the reloaded pending files
        doc_r = json.loads(pend_json.read_text(encoding="utf-8"))
        md_r = pend_md.read_text(encoding="utf-8") if args.profile == "stage2" else None
        log("release gate: reloaded pending files; checking")
        fails, out_res = gate_document(doc_r, md_r, L, system, frozen, args.profile, params, render_ok_input)
        log(f"release gate: {'ok' if not fails else 'FAILED'}")
        for role, pth in inputs:
            if sha256_file(pth) != in_hash[role]:
                fails.append(f"input {role} changed during the run")
        if fails:
            raise GateError("; ".join(fails[:20]))
        written += [out_json]
        os.replace(pend_json, out_json)
        if args.profile == "stage2":
            written += [out_md]
            os.replace(pend_md, out_md)
        rule_counts = {k: 0 for k in ("complementary_total", "complementary_part", "complementary_link",
                                      "complementary_linked")}
        for _, rr in sup.rule_of.items():
            rule_counts[rr] += 1
        info.update({
            "result": "RELEASE GATE OK",
            "rules": {"primary_carried": sum(1 for c in L.count_fields if L.status[c] == "P") - len(L.primary_new),
                      "primary_new": len(L.primary_new),
                      "complementary_carried": sum(1 for c in L.count_fields if L.status[c] == "S"),
                      **rule_counts, "derived_pct": derived["pct"], "derived_q": derived["q"],
                      "derived_dependent": derived["dependent"]},
            "rounds": sup.rounds, "lookaheads": sup.lookaheads, "lp_solves": system.lp_solves,
            "output_audit": {"checked": out_res.checked, "fail_total": len(out_res.unsafe),
                             "residual_groups": len(out_res.residual_groups)},
            "out_primary": doc_r["suppression"]["primary_cells"],
            "out_complementary": doc_r["suppression"]["complementary_cells"],
            "new_cells": new_cells,
        })
        summ = out_dir / OUT_SUMMARY
        written += [summ]
        write_text(summ, summary_markdown(info))
        outs = [out_json] + ([out_md] if args.profile == "stage2" else []) + [summ]
        manifest = {
            "script": SCRIPT_VERSION, "script_sha256": info["script_sha256"], "spec": SPEC_DOC,
            "parameters": {"protection_width": params.width, "s_rule": params.s_rule,
                           "s_lower_bound": params.s_lower, "profile": args.profile},
            "inputs": [{"role": role, "file": pth.name, "bytes": pth.stat().st_size, "sha256": in_hash[role]}
                       for role, pth in inputs],
            "outputs": [{"file": p.name, "bytes": p.stat().st_size, "sha256": sha256_file(p)} for p in outs],
            "release_gate": "ok",
        }
        if drops:
            manifest["dropped"] = [
                {"subtree": t, "reason": why, "count_cells": [c for c in L.dropped if _under(c, t)]}
                for t, why in drops]
        man = out_dir / OUT_MANIFEST
        written += [man]
        write_text(man, json.dumps(manifest, indent=1) + "\n")
        # final verification of the manifest against the files on disk
        m2 = json.loads(man.read_text(encoding="utf-8"))
        for o in m2["outputs"]:
            if sha256_file(out_dir / o["file"]) != o["sha256"]:
                raise GateError(f"manifest hash mismatch for {o['file']}")
    except BaseException:
        log("stopping: removing every file written by this run")
        for pth in written:
            try:
                pth.unlink()
            except FileNotFoundError:
                pass
        for pth in (pend_json, pend_md, out_dir / (OUT_JSON + ".pending.tmp")):
            try:
                pth.unlink()
            except FileNotFoundError:
                pass
        if created_dir:
            try:
                out_dir.rmdir()
            except OSError:
                pass
        raise
    log("release written")
    print(summary_markdown(info), end="")
    print()
    print(man.read_text(encoding="utf-8"), end="")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Stage 2 v1.0.2 differencing-fix release (aggregate output only).")
    ap.add_argument("--in-json", required=True, help="released Stage 2 v1.0.1 JSON (feature_sparsity_frequency_audit.json)")
    ap.add_argument("--in-md", help="released Stage 2 v1.0.1 markdown (required for --profile stage2)")
    ap.add_argument("--links", help="CSV of cross-table equations (total, parts, note)")
    ap.add_argument("--frozen", help="CSV of already-public cells that are never suppressed (pointer, reason)")
    ap.add_argument("--drops", help="CSV of JSON subtrees withheld from the release as a whole (subtree, reason)")
    ap.add_argument("--out-dir", help="new, empty output folder (not an input folder); not used with --check-only")
    ap.add_argument("--check-only", action="store_true", help="audit the input only; write nothing; exit 3 if it fails")
    ap.add_argument("--profile", choices=("stage2", "generic"), default="stage2",
                    help="stage2: Stage 2 JSON + markdown (default); generic: any registry-style JSON, no markdown")
    ap.add_argument("--protection-width", type=int, default=10,
                    help="W: a suppressed cell's feasible interval must hold >= W integers or straddle 9|10 (default 10)")
    ap.add_argument("--s-rule", choices=("width", "exact"), default="width",
                    help="S cells: 'width' applies W (default); 'exact' fails them only when pinned to one value")
    ap.add_argument("--s-lower-bound", type=int, default=MIN_CELL,
                    help="a priori lower bound of an S cell known to a reader (default 10)")
    ap.add_argument("--max-seconds", type=float, default=3600.0,
                    help="overall time limit in seconds (default 3600); on expiry the run fails closed (exit 3)")
    ap.add_argument("--workers", type=int, default=1,
                    help="processes for solving independent components in parallel (default 1; results identical)")
    ap.add_argument("--expect-spec", help="exact spec string required in the input (default: Stage 2 v1.0.1)")
    ap.add_argument("--expect-input-sha256", help="optional exact SHA-256 pin of the input JSON")
    ap.add_argument("--expect-md-sha256", help="optional exact SHA-256 pin of the input markdown")
    a = ap.parse_args(argv)
    if a.max_seconds <= 0 or a.workers < 1:
        ap.error("--max-seconds must be > 0 and --workers >= 1")
    if a.check_only and a.drops:
        ap.error("--drops is used only for a release run, not with --check-only")
    if not a.check_only and not a.out_dir:
        ap.error("--out-dir is required unless --check-only is given")
    return a


def main(argv: list[str] | None = None) -> int:
    try:
        args = parse_args(argv)
    except SystemExit as e:
        return int(e.code or 0) if isinstance(e.code, int) else 2
    import time
    _T0[0] = time.monotonic()
    _DEADLINE[0] = _T0[0] + args.max_seconds
    try:
        return run(args)
    except TimeLimitError:
        print(f"RELEASE GATE FAILED: time limit of {args.max_seconds:g} s (--max-seconds) reached before the run "
              f"finished", file=sys.stderr)
        print("No release file was left behind. Do not patch the script; send this message back.", file=sys.stderr)
        return 3
    except InputError as e:
        print(f"STOPPED (input): {e}", file=sys.stderr)
        print("Nothing was written. Do not patch the script; send this message back.", file=sys.stderr)
        return 2
    except GateError as e:
        print(f"RELEASE GATE FAILED: {e}", file=sys.stderr)
        print("No release file was left behind. Do not patch the script; send this message back.", file=sys.stderr)
        return 3


if __name__ == "__main__":
    sys.exit(main())
