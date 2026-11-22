#!/usr/bin/env python3
"""
Stage 2 v1.0.3 release: differencing (complementary) suppression of the Stage 2 v1.0.1 output, re-gated with the
cross-release findings of Stage 3 (L3 public table sizes, L5 age links, DOB links, the v1.0 stop-event count) and the
hidden-cell exact-bound fix. Module v1.5; implements `stage2_v1_0_3_suppression_spec.md` v1.5 (Nov 17, 2026).
v1.5 (after the real v1.4 run stalled in round 5's suppression selection, Nov 17): every v1.4 rule and guarantee is
kept; it changes only HOW complementary cells are chosen and adds bounds and a clean stop:
  * heartbeat: a background thread prints a line whenever 30 s pass without one (phase, position, LP solves);
  * bounded selection: candidate checks run on a local sub-system (<= 3 groups out, <= 1,500 equations) with a 20 s
    cap per LP, a 60 s budget per failing cell and a 900 s budget per round; nothing is ever accepted because a
    check timed out (a timed-out check counts as "not shown to fix the cell"); the full audit decides every round;
  * pin-ranked, merge-aware, self-safe choice: candidates in the equation that pins the failing cell (fewest
    unknowns) come first, then those that join no other suppressed component, then v1.4's static order; a set of
    picks (a chain of up to 4 cells) is accepted only if the cell it fixes AND every pick pass rule 3 (v1.4 picked
    cells that fixed one cell and were pinned themselves: 75 -> 373 -> 19 -> 256 failing);
  * no mass suppression: the fallback (10 cells) and the stronger rule (40 per failing cell) are capped, 400 per round;
  * convergence stop: more than --max-new-suppressions (2,500) new cells, a component above --max-component-cells
    (20,000), or no improvement for 10 rounds -> exit 3, nothing written, and a WITHHOLD PROPOSAL (pointers only,
    checked by a full audit) printed and, with --withhold-proposal FILE, written as a CSV that --withhold reads;
  * --withhold FILE (pointer, reason): single count cells withheld from the release (a cell-level alternative to
    the D4 drop of a whole subtree); same rules as --drops (never a frozen, public or links-file cell).
v1.4 keeps every v1.3 rule and adds:
  * exact bounds for hidden / dropped cells: screening points may prove such a cell safe only by width or by
    straddling 9|10 (or a screened minimum >= 10); "exactly 0" is accepted only from an exact LP maximum;
  * --public: cells whose value a reader already knows from outside the release (public table sizes, structural
    zeros, a count quoted elsewhere). They are known constants in every audit; their label in the release is not
    changed (team-internal input; pointers and values never printed);
  * --carry-json: the previous release (v1.0.2). Every cell it suppressed or withheld stays suppressed or withheld,
    so v1.0.2 and v1.0.3 side by side give a reader nothing that v1.0.3 alone does not.

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
    python scripts\\stage2_v1_0_3_suppression\\stage2_v1_0_3_suppress.py ^
        --in-json data\\interim\\feature_sparsity_frequency_audit.json ^
        --in-md   data\\interim\\feature_sparsity_frequency_audit.md ^
        --carry-json data\\interim\\stage2_v1.0.2\\feature_sparsity_frequency_audit_v1.0.2.json ^
        --links   scripts\\stage2_v1_0_3_suppression\\stage2_v1_0_3_links.csv ^
        --frozen  scripts\\stage2_v1_0_3_suppression\\stage2_v1_0_3_frozen.csv ^
        --public  scripts\\stage2_v1_0_3_suppression\\stage2_v1_0_3_public_TEAM_INTERNAL.csv ^
        --drops   scripts\\stage2_v1_0_3_suppression\\stage2_v1_0_3_drops.csv ^
        --out-dir data\\interim\\stage2_v1.0.3
    Add --check-only to audit an output without writing anything (exit 0 = passes the new gate, 3 = fails).

Outputs (in --out-dir, which must be new or empty and must not be an input folder):
    feature_sparsity_frequency_audit_v1.0.3.json / .md   shared (suppressed)
    stage2_v1.0.3_release_summary.md                     counts of cells by rule; pointers only, no values
    stage2_v1.0.3_manifest.json                          SHA-256 of inputs and outputs
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

SCRIPT_VERSION = "stage2_v1_0_3_suppress v1.5"
SPEC_DOC = "stage2_v1_0_3_suppression_spec v1.5"
EXPECT_INPUT_SPEC = "feature_sparsity_frequency_audit v1.0.1 (feature_sparsity_frequency_audit_spec v0.8, stage 2)"
OUTPUT_SPEC = ("feature_sparsity_frequency_audit v1.0.3 (v1.0.1 output with the differencing fix, re-gated with the "
               "cross-release links; stage2_v1_0_3_suppression_spec v1.5, stage 2)")
OUT_JSON = "feature_sparsity_frequency_audit_v1.0.3.json"
OUT_MD = "feature_sparsity_frequency_audit_v1.0.3.md"
OUT_SUMMARY = "stage2_v1.0.3_release_summary.md"
OUT_MANIFEST = "stage2_v1.0.3_manifest.json"
PROVENANCE_KEY = "release_v1_0_3"
CARRY_PROVENANCE_KEY = "release_v1_0_2"   # provenance block of the carried (previous) release
RELEASE_LABEL = "v1.0.3"

MIN_CELL = 10
UNKNOWN = ("P", "S", "H", "D")    # statuses whose value a reader does not see (D: dropped by the --drops list)
#                                   v1.4: status "known" = a --public cell (value known to a reader from outside the
#                                   release; a constant in every equation, whatever its label in the release)
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
# v1.5 selection bounds and convergence stop
HEARTBEAT_SECONDS = 30.0          # a heartbeat line whenever this long passes without any progress line
SELECT_LP_SECONDS = 20.0          # cap per candidate-check LP (a capped LP never counts as a fix)
SELECT_CELL_SECONDS = 60.0        # selection budget per failing cell and round (then: left for the next round)
SELECT_ROUND_SECONDS = 900.0      # selection budget per round (then: the round ends with the picks made so far)
LOCAL_RADIUS = 3                  # candidate checks use equations within this many groups of the failing cell
LOCAL_MAX_RELS = 1500             # ... and at most this many equations (fewer equations only widen intervals)
PREFIX_MAX = 40                   # a prefix (set) fix may hold at most this many cells
ESCALATE_MAX_CELL = 40            # capped "stronger rule": at most this many cells per failing cell
ESCALATE_MAX_ROUND = 400          # ... and per round
FALLBACK_MAX_CELL = 10            # large component, no checked fix found: the 10 nearest least-merging cells
CHAIN_DEPTH = 4                   # chain search: at most this many cells in one checked fix
CHAIN_WIDTH = 4                   # ... trying this many least-merging candidates per step
CHAIN_NODES = 60                  # ... and at most this many exact checks per failing cell
MAX_NEW_DEFAULT = 2500            # --max-new-suppressions
MAX_COMPONENT_DEFAULT = 20000     # --max-component-cells
STALL_ROUNDS = 5                  # no-progress guard: rounds without a new low failing count before escalating
#                                   (the v1.0 look-ahead path legitimately goes 3 rounds without one on synthetic cases)

# Stage 2 v1.0.1 schema constants (copied from feature_sparsity_frequency_audit.py v1.0.1, sha256 3b5b36fb...)
SPLITS = ("train", "validation", "test")
DOB_STATUS = ("missing", "unparseable", "negative", "valid")
DOB_NEG_BINS = ("lt1d", "1to30d", "gt30d")
WEIGHT_SOURCES = ("charted", "imputed_age_sex", "imputed_age_band", "imputed_sex_only", "imputed_overall")
NON_COUNT_TOP = ("disclosure", "suppression", "input_hashes", "v17_module", "raw_tables", "maps", "cohort_checks",
                 "config", "outcome_descriptors_flag", "spec", "release_gate", "provisional_classes", "notes",
                 PROVENANCE_KEY, CARRY_PROVENANCE_KEY)
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


def read_public(path: Path | None) -> list[tuple[str, int | None, str]]:
    """v1.4. CSV columns: pointer, value, reason. A cell whose value a reader already knows from outside this release
    (a public table size, a structural zero, a count quoted in a shared document). `value` empty: the input's
    published value is used (the cell must be published in the input); an integer: the asserted public value (must
    equal the input's value when the input publishes the cell); 'derived': a cell that the other public / published
    values fix exactly through one identity or link (e.g. valid DOBs = cohort - public invalid DOBs); the module
    derives the value itself, so the file never needs to hold it (a cell published in the input keeps its input
    value, as for an empty value); 'publish': the input's published value, and the cell is PUBLISHED in this release
    even where the carried release suppressed it (only for a value public by definition, e.g. a raw table size; an
    explicit decision, never a default). Team-internal: pointers and values are never printed
    or written to any output (the output records the file's SHA-256 and the number of rows only)."""
    if path is None:
        return []
    df = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    if not {"pointer", "value", "reason"}.issubset(df.columns):
        raise InputError(f"public file {path.name} needs columns ['pointer', 'value', 'reason']")
    out, seen = [], set()
    for i, r in df.iterrows():
        ptr, val, why = r["pointer"].strip(), r["value"].strip(), r["reason"].strip()
        if not ptr:
            continue
        if not ptr.startswith("/") or not why:
            raise InputError(f"public file {path.name}, row {i + 2}: needs a JSON pointer starting with '/' and a reason")
        if ptr in seen:
            raise InputError(f"public file {path.name}, row {i + 2}: duplicate pointer")
        seen.add(ptr)
        if val == "":
            v = None
        elif re.fullmatch(r"\d+", val):
            v = int(val)
        elif val in ("derived", "publish"):
            v = val
        else:
            raise InputError(f"public file {path.name}, row {i + 2}: value must be empty, a non-negative integer, "
                             f"'derived' or 'publish'")
        out.append((ptr, v, why))
    return out


def apply_public(L: "Loaded", public: list, frozen: set[str], link_rels=()) -> None:
    """v1.4: mark --public cells as known constants (status 'known'), keeping their release label. 'derived' rows are
    resolved last, each from one identity / link whose other cells are all published or public."""
    derived = [(p, v, w) for p, v, w in public if v == "derived"]
    for ptr, v, _why in public:
        if v == "derived":
            continue
        if v == "publish":
            row = public.index((ptr, v, _why)) + 2
            if ptr not in set(L.count_fields) or L.input_status.get(ptr) != "pub" or L.status.get(ptr) == "D":
                raise InputError(f"'publish' public cell is not a published count cell of the input (row {row})")
            if ptr in L.carried:
                L.carried.remove(ptr)
                L.republished.append(ptr)
            L.public_label[ptr] = "pub"
            L.status[ptr] = "known"
            continue
        st = L.status.get(ptr)
        if st is None or ptr not in set(L.count_fields):
            raise InputError(f"public pointer is not a count cell of the input (row {public.index((ptr, v, _why)) + 2})")
        if st == "D":
            raise InputError(f"public cell lies inside a dropped subtree (row {public.index((ptr, v, _why)) + 2})")
        if ptr in frozen:
            raise InputError(f"public cell is also in the frozen file (row {public.index((ptr, v, _why)) + 2})")
        if L.input_status.get(ptr, st) == "pub":      # published in the input (v1.0.1), even if carried as S
            if v is not None and v != L.values[ptr]:
                raise InputError(f"public value contradicts the published input value (row "
                                 f"{public.index((ptr, v, _why)) + 2}); nothing is printed")
            val = L.values[ptr]
        else:
            if v is None:
                raise InputError(f"public cell is suppressed in the input and the public file gives no value (row "
                                 f"{public.index((ptr, v, _why)) + 2})")
            if st == "P" and v > SMALL_MAX:
                raise InputError(f"public value contradicts the input's '<10' label (row {public.index((ptr, v, _why)) + 2})")
            if st == "S" and v < MIN_CELL:
                raise InputError(f"public value contradicts the input's 'S' label (row {public.index((ptr, v, _why)) + 2})")
            val = v
        L.public_label[ptr] = {"pub": "pub", "P": "P", "S": "S"}[st]
        L.status[ptr] = "known"
        L.values[ptr] = val
    if not derived:
        return
    eqs = [(t, tuple(ps)) for t, ps in L.identities] + [(t, tuple(ps)) for t, ps, _n in link_rels]
    known = lambda c: L.status.get(c) in ("pub", "known", "const") and c in L.values  # noqa: E731
    todo = {}
    for ptr, _v, _why in derived:
        row = public.index((ptr, _v, _why)) + 2
        st = L.status.get(ptr)
        if st is None or ptr not in set(L.count_fields):
            raise InputError(f"public pointer is not a count cell of the input (row {row})")
        if st == "D" or ptr in frozen:
            raise InputError(f"derived public cell is dropped or frozen (row {row})")
        if L.input_status.get(ptr, st) == "pub":    # published in the input: its value is known, as for an empty value
            L.public_label[ptr] = {"pub": "pub", "S": "S"}[st]
            L.status[ptr] = "known"
            continue
        todo[ptr] = (row, st)
    # a derived cell must be fixed by the equations ALONE (every unknown free, no bounds): exact LP min = max over
    # the cell's connected block of equations; the value never comes from the '<10' / 'S' bounds
    by_cell: dict[str, list[int]] = {}
    for i, (t, ps) in enumerate(eqs):
        for c in (t,) + ps:
            if not known(c):
                by_cell.setdefault(c, []).append(i)
    for c in list(todo):
        rows, cells, stack = set(), {c}, [c]
        while stack:
            for i in by_cell.get(stack.pop(), []):
                if i in rows:
                    continue
                rows.add(i)
                t, ps = eqs[i]
                for q in (t,) + ps:
                    if not known(q) and q not in cells:
                        cells.add(q)
                        stack.append(q)
        cl = sorted(cells)
        col = {q: j for j, q in enumerate(cl)}
        coef: dict = defaultdict(float)
        b = np.zeros(len(rows))
        for k, i in enumerate(sorted(rows)):
            t, ps = eqs[i]
            for q, a in [(t, 1.0)] + [(q, -1.0) for q in ps]:
                if known(q):
                    b[k] -= a * L.values[q]
                else:
                    coef[(k, col[q])] += a
        kk = list(coef)
        A = sparse.csr_matrix(([coef[x] for x in kk], ([x[0] for x in kk], [x[1] for x in kk])),
                              shape=(len(rows), len(cl)))
        e = np.zeros(len(cl))
        e[col[c]] = 1.0
        free = [(None, None)] * len(cl)
        lo = linprog(e, A_eq=A, b_eq=b, bounds=free, method="highs") if len(rows) else None
        hi = linprog(-e, A_eq=A, b_eq=b, bounds=free, method="highs") if len(rows) else None
        if lo is None or lo.status != 0 or hi.status != 0 or abs(lo.x[col[c]] - hi.x[col[c]]) > 1e-6:
            continue
        v = int(round(lo.x[col[c]]))
        row, st = todo.pop(c)
        if v < 0 or (st == "P" and v > SMALL_MAX) or (st == "S" and v < MIN_CELL):
            raise InputError(f"derived public value contradicts the input's label (row {row}); nothing is printed")
        L.public_label[c] = {"P": "P", "S": "S"}[st]
        L.status[c] = "known"
        L.values[c] = v
    if todo:
        raise InputError(f"{len(todo)} 'derived' public cell(s) are not fixed by the identities and links over published / "
                         f"public cells (first row {min(r for r, _ in todo.values())})")


def apply_carry(L: "Loaded", carry_doc: dict, in_sha: str) -> None:
    """v1.4: the previous release (v1.0.2) as released. Every count cell it suppressed ('<10' / 'S') or withheld
    (dropped subtree) stays suppressed / withheld; every value it published must equal the input's. Run after the
    drop list (a withheld cell must be dropped again) and before --public (a public cell keeps the carried label)."""
    if carry_doc.get("release_gate") != "ok":
        raise InputError("carry JSON is not a released output (release_gate != 'ok')")
    prov = carry_doc.get(CARRY_PROVENANCE_KEY)
    if not isinstance(prov, dict) or prov.get("source_json_sha256") != in_sha:
        raise InputError("carry JSON is not the v1.0.2 release built from this input (its provenance "
                         "source_json_sha256 differs from the input's SHA-256)")
    dis = carry_doc.get("disclosure")
    if not isinstance(dis, dict) or "count_fields" not in dis:
        raise InputError("carry JSON has no disclosure block")
    cfs = list(dis["count_fields"])
    in_set = set(L.count_fields)
    extra = [c for c in cfs if c not in in_set]
    if extra:
        raise InputError(f"carry JSON declares {len(extra)} count cell(s) that the input does not have")
    cset = set(cfs)
    n_s = 0
    for c in L.count_fields:
        st = L.status[c]
        if st == "D":
            continue
        if c not in cset:
            raise InputError(f"a cell withheld in the carried release is not withheld here (add its subtree to the "
                             f"drop list): {c}")
        try:
            v = get_ptr(carry_doc, c)
        except KeyError:
            raise InputError(f"carry JSON declares a count cell it does not hold: {c}")
        if isinstance(v, bool):
            raise InputError(f"carry JSON holds a boolean count at {c}")
        if isinstance(v, int):
            if st != "pub" or v != L.values[c]:
                raise InputError(f"carry JSON publishes a value that differs from the input at {c}")
        elif v == "<10":
            if st != "P":
                raise InputError(f"carry JSON labels '<10' a cell that is not '<10' in the input at {c}")
        elif v == "S":
            if st == "P":
                raise InputError(f"carry JSON labels S a '<10' cell at {c}")
            if st == "pub":
                L.status[c] = "S"
                L.carried.append(c)
                n_s += 1
        else:
            raise InputError(f"carry JSON holds an invalid count value at {c}")


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


def read_withhold(path: Path | None) -> list[tuple[str, str]]:
    """v1.5. CSV columns: pointer, reason. Each row names ONE count cell withheld from the release (removed from the
    JSON and the registry, audited as an unpublished cell, like a cell of a dropped subtree). The cell-level
    alternative to dropping a whole subtree (D4); usually the rows of a --withhold-proposal file, reviewed by Vinai."""
    if path is None:
        return []
    df = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    if not {"pointer", "reason"}.issubset(df.columns):
        raise InputError(f"withhold file {path.name} needs columns ['pointer', 'reason']")
    out = []
    for _, r in df.iterrows():
        t = r["pointer"].strip()
        if not t:
            continue
        if not t.startswith("/") or not r["reason"].strip():
            raise InputError(f"withhold file row needs a JSON pointer starting with '/' and a reason: {t!r}")
        out.append((t, r["reason"].strip()))
    if not out:
        raise InputError(f"withhold file {path.name} lists no cell")
    if len({t for t, _ in out}) != len(out):
        raise InputError(f"withhold file {path.name} lists a cell twice")
    return out


def _under(ptr: str, sub: str) -> bool:
    return ptr == sub or ptr.startswith(sub + "/")


def apply_drops(L: Loaded, drops: list[tuple[str, str]], link_cells: set[str], frozen: set[str],
                withhold: list[tuple[str, str]] = ()) -> None:
    """Mark every count cell under a dropped subtree as D (unpublished; still an unknown in every equation, audited
    with the hidden-cell rule) and build the output registry without the dropped pointers. Fails on anything
    ambiguous: a missing or non-object subtree, overlapping rows, a subtree holding no count cell, a frozen cell or
    a links-file cell inside a dropped subtree.
    v1.5: `withhold` rows name single count cells (leaves), handled exactly like one-cell subtrees."""
    if not drops and not withhold:
        L.out_disclosure = L.doc.get("disclosure")
        return
    cf = set(L.count_fields)
    for t, _ in withhold:
        if t not in cf:
            raise InputError(f"withhold file pointer is not a count cell of the input: {t}")
        if t in frozen:
            raise InputError(f"a frozen (already public) cell cannot be withheld: {t}")
        if t in link_cells:
            raise InputError(f"a links-file cell cannot be withheld (send it back to Vinai): {t}")
        for u, _ in drops:
            if _under(t, u):
                raise InputError(f"withhold file cell lies inside a dropped subtree: {t}")
    L.withheld = [t for t, _ in withhold]
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

    subs = subs + [t for t, _ in withhold]

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
    L.drop_subtrees = list(drops) + list(withhold)
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


def screen_safe(kind: str, lo: float, hi: float, p: Params) -> bool:
    """v1.4: what a set of FEASIBLE POINTS (an inner approximation [lo, hi] of the true interval) may prove. Rule 3 is
    monotone in the interval for P and S cells, and for H / D cells by width, by straddling 9|10 or by a minimum >= 10
    (a wider true interval then still has lo >= 10 or straddles). Only the H / D clause "exactly 0" (hi <= 0) is not
    monotone: points that all sit at 0 do not show that the cell is 0 (v1.3 accepted it; its true range may be 0-8).
    Such cells stay pending and get an exact LP maximum."""
    if kind in ("H", "D") and int_bounds(lo, hi)[1] <= 0:
        return False                         # needs the exact maximum (solve_component's per-cell LP)
    return is_safe(kind, lo, hi, p)


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
_LAST_LINE = [0.0]                # v1.5: time of the last line printed (progress or heartbeat)
_HB = {"phase": "starting", "detail": "", "lp": 0}
import threading as _threading   # noqa: E402
_PRINT_LOCK = _threading.Lock()


def log(msg: str) -> None:
    """Progress to stderr (counts and timings only, never values)."""
    import time
    with _PRINT_LOCK:
        now = time.monotonic()
        _LAST_LINE[0] = now
        print(f"[{now - _T0[0]:8.1f}s] {msg}", file=sys.stderr, flush=True)


def beat(phase: str | None = None, detail: str = "") -> None:
    """v1.5: record what the run is doing (read by the heartbeat thread). Counts and positions only."""
    if phase is not None:
        _HB["phase"] = phase
    _HB["detail"] = detail


def lp_count(n: int = 1) -> None:
    _HB["lp"] += n


class Heartbeat:
    """v1.5: a daemon thread that prints '... heartbeat: ...' whenever HEARTBEAT_SECONDS pass without a line. HiGHS
    releases the GIL while it solves, so the line appears even during one long LP. Never prints values."""

    def __init__(self, every: float = HEARTBEAT_SECONDS):
        self.every = every
        self.stop_ev = _threading.Event()
        self.t = _threading.Thread(target=self._loop, name="heartbeat", daemon=True)

    def _loop(self) -> None:
        import time
        while not self.stop_ev.wait(min(5.0, self.every / 3)):
            if time.monotonic() - _LAST_LINE[0] >= self.every:
                d = _HB["detail"]
                log(f"heartbeat: still working - {_HB['phase']}" + (f" ({d})" if d else "")
                    + f"; LP solves so far {_HB['lp']}; time left {max(0.0, time_left()):.0f} s")

    def __enter__(self):
        self.t.start()
        return self

    def __exit__(self, *exc):
        self.stop_ev.set()
        self.t.join(timeout=2)
        return False


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
        lp_count()
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

    def pending_of() -> list[int]:   # v1.4: screen_safe (points prove "exactly 0" for no hidden / dropped cell)
        return [j for j in range(n) if not screen_safe(kinds[j], obs_lo[j], obs_hi[j], p)]
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
                beat(label or "audit", f"component {i + 1}/{len(todo)}, {len(pl['vars'])} cells")
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
    def target_bounds(self, status, rel_ids: Iterable[int], targets: list[str],
                      lp_cap: float | None = None) -> dict[str, tuple[float, float] | None]:
        """Exact feasible [min, max] of each target cell over the given equations (exact propagation, then HiGHS
        min / max LPs on the target's remaining component). Used only to choose suppressions; every decision is
        re-checked by the full audit and, independently, by the release gate.
        v1.5: with lp_cap, each LP gets at most lp_cap seconds; a target whose LP hits the cap maps to None
        ("unknown"), which a caller never treats as safe."""
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
            capped = False
            for sign in (1.0, -1.0):
                check_time()
                self.lp_solves += 1
                lp_count()
                limit = max(time_left(), 1.0) if lp_cap is None else min(lp_cap, max(time_left(), 1.0))
                r = linprog(sign * e, A_eq=A, b_eq=np.array(rhs, dtype=float), bounds=bnd_cap, method="highs",
                            options={"time_limit": float(min(limit, 1e7))})
                if r.status == 2:
                    raise InputError("published cells are inconsistent (LP infeasible): check identities and links")
                if r.status != 0:
                    if time_left() <= 1.0:
                        raise TimeLimitError()
                    if r.status == 1 and lp_cap is not None:
                        capped = True          # v1.5: capped selection LP -> unknown (never "safe")
                        break
                    raise InputError(f"LP solver failed (status {r.status}: {r.message})")
                vals.append(float(r.x[col[t]]))
            if capped:
                out[t] = None
                continue
            lo, hi = vals
            out[t] = (lo, INF if hi >= BIG / 2 else hi)
        return out

    def local_rels(self, seeds: Iterable[str], status, radius: int = LOCAL_RADIUS,
                   max_rels: int = LOCAL_MAX_RELS) -> set[int]:
        """v1.5: the equations within `radius` groups of the seeds, walking through unknown cells only, at most
        max_rels of them (nearest first). A subset of the component's equations: intervals computed on it can only
        be WIDER than the true ones, so it is used to choose candidates (re-checked by the full audit) and to prove
        that a cell cannot be protected (a cell failing on fewer equations fails on all of them)."""
        seen_c = set(seeds)
        seen_k: set[int] = set()
        frontier = sorted(seen_c)
        for _ in range(radius):
            nxt = []
            for c in frontier:
                for k in self.by_cell.get(c, []):
                    if k in seen_k:
                        continue
                    if len(seen_k) >= max_rels:
                        return seen_k
                    seen_k.add(k)
                    for m, _ in self.rels[k]:
                        if m not in seen_c and status[m] in UNKNOWN:
                            seen_c.add(m)
                            nxt.append(m)
            frontier = sorted(nxt)
            if not frontier:
                break
        return seen_k

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
        self.res: AuditResult | None = None          # v1.5: the current round's audit (components for merge_key)
        self.max_new = MAX_NEW_DEFAULT
        self.max_component = MAX_COMPONENT_DEFAULT
        self.cell_deadline = INF
        self.fail_at_escalation = INF
        self.last_failing_adds: list[str] = []
        self.paths = {k: 0 for k in ("lookahead_small", "single", "chain", "prefix", "single_fixing_only_u",
                                     "fallback", "stronger_rule")}

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

    # ---- v1.5: merge-aware ordering and bounded, self-safe exact checks ------------------------------------------
    def merge_key(self, m: str, u: str, st=None) -> tuple:
        """v1.5 preference order of a candidate m for the failing cell u (smaller first):
        1) 'pin': the fewest unknown cells in an equation that m shares with u (1 = an equation in which u is the
           only unknown, i.e. one that fixes u; freeing it is what helps); 99 when m shares no equation with u;
        2) 'join': the cells of OTHER suppressed components that suppressing m would join to u's (suppressing a
           cell turns its groups into equations among unknowns; how the real run grew 7,738 -> 36,087 cells);
        3) v1.4's static key (fewest groups, a total first, smallest value, pointer)."""
        S, res = self.sys, self.res
        st = self.status if st is None else st
        cu = res.comp_of.get(u) if res is not None else None
        join: set[int] = set()
        pin = 99
        ku = set(S.by_cell.get(u, []))
        for k in S.by_cell.get(m, []):
            if k in ku:
                pin = min(pin, sum(1 for c, _ in S.rels[k] if st[c] in UNKNOWN))
            if res is None:
                continue
            for c, _ in S.rels[k]:
                if c != m and st[c] in UNKNOWN:
                    ci = res.comp_of.get(c)
                    if ci is not None and ci != cu:
                        join.add(ci)
        size = sum(res.comp_size.get(ci, 0) for ci in join) if res is not None else 0
        return (pin, size) + self.static_key(m, u)

    def check_fix(self, u: str, adds: list[str], check_adds: bool = True) -> tuple[bool, bool]:
        """v1.5 exact candidate check on the LOCAL sub-system (local_rels; capped LPs). Returns (u passes,
        u and every added cell pass). An LP that hits its cap gives 'not shown' (False), never 'safe'.
        The added cells that do not pass are left in self.last_failing_adds."""
        self.last_failing_adds = []
        self.exact_checks += 1
        p = self.sys.params
        st = Overlay(self.status, {m: "S" for m in adds})
        rels = self.sys.local_rels([u, *adds], st)
        targets = [u] + ([m for m in adds if m != u] if check_adds else [])
        b = self.sys.target_bounds(st, rels, targets, lp_cap=SELECT_LP_SECONDS)
        bu = b.get(u)
        ok_u = bu is not None and is_safe(st[u], bu[0], bu[1], p) and not self.sys.residual_fails(st, u)
        if not ok_u or not check_adds:
            return ok_u, False
        for m in targets[1:]:
            bm = b.get(m)
            if bm is None or not is_safe("S", bm[0], bm[1], p) or self.sys.residual_fails(st, m):
                self.last_failing_adds.append(m)
        return True, not self.last_failing_adds

    def makes_safe(self, u: str, adds: list[str]) -> bool:
        """With `adds` suppressed, does u meet rule 3 and the residual check (local sub-system, capped LPs)?"""
        return self.check_fix(u, adds, check_adds=False)[0]

    def unprotectable(self, u: str) -> bool:
        """Proof that no suppression can protect u: free EVERY publishable, non-frozen cell in the equations near u
        (v1.5: local_rels, i.e. a subset of u's component; dropping equations only widens intervals). If u still
        fails, any real suppression pattern (more equations, fewer free cells) fails too. A capped LP proves
        nothing (False)."""
        self.exact_checks += 1
        rels = self.sys.local_rels([u], self.status)
        free = {m for k in rels for m, _ in self.sys.rels[k] if self.pickable(m, self.status)}
        st = Overlay(self.status, {m: "S" for m in free})
        b = self.sys.target_bounds(st, rels, [u], lp_cap=SELECT_LP_SECONDS)[u]
        if b is None:
            return False
        return not is_safe(st[u], b[0], b[1], self.sys.params) or self.sys.residual_fails(st, u)

    def choose_exact(self, u: str, deadline: float | None = None):
        """Large components (v1.2, reworked in v1.5). Returns a list of cells to suppress, None (nothing found
        within the rings) or "timeout" (the per-cell budget ran out; u is left for the next round).
        1) single candidates from u's own groups in merge_key order (+ the largest values that join nothing):
           the first whose check leaves u AND itself passing wins; a single that fixes u but fails itself is
           paired with one of ITS own candidates (both checked); only then a single that fixes u alone;
        2) otherwise the smallest prefix (merge_key order, nearest groups first, at most PREFIX_MAX cells) of the
           candidates within 1, 2, ... EXACT_RINGS groups that makes u pass."""
        import time
        if deadline is None:
            deadline = self.cell_deadline

        def late() -> bool:
            return time.monotonic() > deadline
        r1 = [m for d, m in self.ring(u, 1)]
        alone: list[str] | None = None
        if r1:
            order = sorted(r1, key=lambda x: self.merge_key(x, u))
            singles = order[:EXACT_SINGLES]
            singles += [m for m in sorted(r1, key=lambda x: (-self.sys.values[x], x))
                        if m not in singles and self.merge_key(m, u)[1] == 0][:EXACT_TOP_VALUE]
            half: list[str] = []
            for m in singles:
                if late():
                    return "timeout"
                ok_u, ok_all = self.check_fix(u, [m])
                if ok_all:
                    self.paths["single"] += 1
                    return [m]
                if ok_u:
                    half.append(m)
            # chain search (v1.5): breadth-first over small sets, starting from the singles (those that fix u
            # first). A set that leaves u or one of its own cells failing is extended with the least-merging
            # candidates of the failing cell and of the set's last cell, up to CHAIN_DEPTH cells and CHAIN_NODES
            # checks per failing cell. A set is accepted only when u and every cell of the set pass.
            start = half + [m for m in singles if m not in half]
            frontier = [[m] for m in start[:CHAIN_WIDTH]]
            seen_sets = {tuple(x) for x in frontier}
            nodes = 0
            while frontier and nodes < CHAIN_NODES:
                adds = frontier.pop(0)
                if late():
                    return "timeout"
                if len(adds) > 1:
                    nodes += 1
                    ok_u, ok_all = self.check_fix(u, adds)
                    if ok_all:
                        self.paths["chain"] += 1
                        return adds
                else:
                    ok_u = adds[0] in half
                    if not ok_u:
                        self.last_failing_adds = []
                    else:
                        self.check_fix(u, adds)
                        nodes += 1
                if len(adds) >= CHAIN_DEPTH:
                    continue
                ws = ([u] if not ok_u else []) + list(self.last_failing_adds) + [adds[-1]]
                pool: set[str] = set()
                for w in ws[:2]:
                    pool.update(m2 for d, m2 in self.ring(w, 1))
                pool.discard(u)
                pool.difference_update(adds)
                st_a = Overlay(self.status, {m: "S" for m in adds})
                fix = ws[:2]
                for m2 in sorted(pool, key=lambda x: min(self.merge_key(x, w, st_a) for w in fix))[:CHAIN_WIDTH]:
                    cand = adds + [m2]
                    key = tuple(sorted(cand))
                    if key not in seen_sets:
                        seen_sets.add(key)
                        frontier.append(cand)
            if half:
                alone = [half[0]]
        for depth in range(1, EXACT_RINGS + 1):
            if late():
                return alone if alone is not None else "timeout"
            cand = self.ring(u, depth)
            if not cand:
                continue
            ordered = [m for d, _k, m in sorted((d, self.merge_key(m, u), m) for d, m in cand)][:PREFIX_MAX]
            if not self.makes_safe(u, ordered):
                continue
            lo, hi = 1, len(ordered)
            while lo < hi:
                if late():
                    return alone if alone is not None else "timeout"
                mid = (lo + hi) // 2
                if self.makes_safe(u, ordered[:mid]):
                    hi = mid
                else:
                    lo = mid + 1
            pref = ordered[:lo]
            if alone is not None and len(pref) > 1:
                self.paths["single_fixing_only_u"] += 1
                return alone
            self.paths["prefix"] += 1
            return pref
        if alone is not None:
            self.paths["single_fixing_only_u"] += 1
        return alone

    def capped_ring(self, u: str, budget: int) -> list[str]:
        """v1.5: the 'stronger rule' (candidates within EXACT_RINGS groups), nearest and least-merging first, at most
        ESCALATE_MAX_CELL cells and the round's remaining budget."""
        cand = self.ring(u, EXACT_RINGS)
        ordered = [m for d, _k, m in sorted((d, self.merge_key(m, u), m) for d, m in cand)]
        return ordered[:max(0, min(ESCALATE_MAX_CELL, budget))]

    def run(self) -> AuditResult:
        """Rounds of: audit everything (unchanged components come from the cache), then new suppressions.
        - Components up to LARGE_COMPONENT cells: one new suppression per failing component, chosen by
          look-ahead (unchanged from v1.0).
        - Larger components: one fix per failing GROUP, chosen by a bounded exact check (choose_exact, v1.5); a
          cell that provably cannot be protected stops the run (fail closed, D4).
        - No-progress guard: if the failing count has not gone below its best for STALL_ROUNDS rounds, failing
          cells get the capped stronger rule; after another STALL_ROUNDS rounds without improvement the run stops.
        - v1.5 convergence stop (NotConverging, exit 3, withhold proposal): more than max_new new suppressions, a
          component above max_component cells, or the guard above."""
        import time
        best = INF
        stall = 0
        for r in range(1, MAX_ROUNDS + 1):
            check_time()
            beat(f"suppression round {r}: audit", "")
            res = self.sys.audit(self.status, label=f"suppression round {r}: audit")
            self.res = res
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
            largest = max(res.comp_size.values(), default=0)
            if len(self.rule_of) > self.max_new:
                raise NotConverging(f"{len(self.rule_of)} new suppressions exceed the limit of {self.max_new} "
                                    f"(--max-new-suppressions); {n_fail} cell(s) still failing", res)
            if largest > self.max_component:
                raise NotConverging(f"a component of {largest} unknown cells exceeds the limit of "
                                    f"{self.max_component} (--max-component-cells); {n_fail} cell(s) still "
                                    f"failing", res)
            if stall >= 2 * STALL_ROUNDS:
                raise NotConverging(f"{n_fail} cell(s) still failing and the failing count has not decreased in "
                                    f"{stall} rounds, {STALL_ROUNDS} of them with the stronger rule", res)
            if stall > STALL_ROUNDS and n_fail > self.fail_at_escalation:
                raise NotConverging(f"the stronger rule made it worse: {n_fail} cell(s) failing now, "
                                    f"{self.fail_at_escalation} when it started", res)
            escalate = stall >= STALL_ROUNDS
            if stall == STALL_ROUNDS:
                self.fail_at_escalation = n_fail
            log(f"suppression round {r}: {n_fail} failing cells, {trend}"
                + ("; ESCALATING to the stronger rule" if escalate else ""))
            self.rounds += 1
            touched: set[int] = set()
            touched_rel: set[int] = set()
            progressed = 0
            skipped = 0
            esc_budget = ESCALATE_MAX_ROUND
            unprot: list[str] = []
            fail_list = sorted(res.unsafe)
            t_round = time.monotonic()
            last_line = t_round
            for i, u in enumerate(fail_list):
                check_time()
                now = time.monotonic()
                beat(f"suppression round {r}: selection",
                     f"failing cell {i + 1}/{len(fail_list)}, {progressed} picks, {self.exact_checks} exact checks")
                if now - last_line >= 60.0:
                    log(f"suppression round {r}: selection {i}/{len(fail_list)} failing cells handled, "
                        f"{progressed} new suppressions so far, {self.exact_checks} exact candidate checks")
                    last_line = now
                if now - t_round > SELECT_ROUND_SECONDS:
                    skipped += len(fail_list) - i
                    log(f"suppression round {r}: selection budget of {SELECT_ROUND_SECONDS:g} s used; "
                        f"{len(fail_list) - i} failing cell(s) left for the next round")
                    break
                if self.status[u] not in UNKNOWN:
                    continue
                comp = res.comp_of.get(u, -1)
                large = res.comp_size.get(comp, 0) > LARGE_COMPONENT
                if escalate:
                    picks = self.capped_ring(u, esc_budget)
                    esc_budget -= len(picks)
                    self.paths["stronger_rule"] += 1
                    if not picks and not self.ring(u, EXACT_RINGS) and self.unprotectable(u):
                        unprot.append(u)
                elif large:
                    if touched_rel.intersection(self.sys.by_cell.get(u, [])):
                        continue
                    self.cell_deadline = min(time.monotonic() + SELECT_CELL_SECONDS, t_round + SELECT_ROUND_SECONDS)
                    picks = self.choose_exact(u)
                    if picks == "timeout":
                        skipped += 1
                        continue
                    if picks is None:
                        if self.unprotectable(u):
                            unprot.append(u)
                            continue
                        # protectable, but only from further out: the capped stronger rule, else the nearest
                        # candidate; the no-progress guard and the v1.5 limits bound this
                        picks = self.capped_ring(u, min(esc_budget, FALLBACK_MAX_CELL))
                        esc_budget -= len(picks)
                        self.paths["fallback"] += 1
                        if not picks:
                            m = self.choose(u, lookahead=False)
                            picks = [] if m is None else [m]
                else:
                    if comp in touched:
                        continue
                    m = self.choose(u, lookahead=True)
                    picks = [] if m is None else [m]
                    self.paths["lookahead_small"] += 1
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
                f"simulations {self.lookaheads}; exact candidate checks {self.exact_checks}"
                + (f"; {skipped} failing cell(s) left for the next round (time budget)" if skipped else "")
                + "; choices " + ", ".join(f"{k} {v}" for k, v in self.paths.items() if v)
                + f"; selection {time.monotonic() - t_round:.0f} s)")
            if unprot:
                raise GateError(f"cannot protect {len(unprot)} cell(s): even with every publishable, non-frozen cell "
                                f"in their linked equations suppressed, the protection rule cannot be met; they are "
                                f"fixed by frozen cells, links or constants (first pointers: {sorted(unprot)[:5]})")
            if not progressed and escalate:
                raise NotConverging(f"{n_fail} cell(s) still failing after the stronger rule found no further "
                                    f"candidate within {EXACT_RINGS} groups", res)
            if not progressed and not skipped:
                bad = sorted(res.unsafe)
                raise GateError(f"cannot protect {len(bad)} cell(s): no publishable, non-frozen cell left to "
                                f"suppress in their groups (first pointers: {bad[:5]})")
        raise GateError("suppression not converging: round limit reached")


class NotConverging(GateError):
    """v1.5: the suppression rounds stopped without converging (exit 3; a withhold proposal is computed)."""

    def __init__(self, why: str, res: AuditResult):
        super().__init__(f"suppression not converging: {why}")
        self.res = res


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
    public_label: dict = field(default_factory=dict)  # v1.4: --public cell -> its release label ('pub' / 'P' / 'S')
    carried: list = field(default_factory=list)       # v1.4: cells published in the input, 'S' in the carried release
    input_status: dict = field(default_factory=dict)  # v1.4: statuses as read from the input (before carry / public)
    republished: list = field(default_factory=list)   # v1.4: 'publish' public cells that the carried release suppressed
    withheld: list = field(default_factory=list)      # v1.5: single cells withheld (--withhold)


def check_name(p: Path, role: str) -> None:
    if "_raw_local" in p.name.lower():
        raise InputError(f"{role} '{p.name}' is a _raw_LOCAL file: never used as input or output")


def load_output(json_path: Path, md_path: Path | None, profile: str, expect_spec: str | None,
                allow_released: bool = False) -> Loaded:
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
    if not allow_released and any(str(k).startswith("release_v1_0_") for k in doc):
        raise InputError("input already carries a release provenance block (it is not the v1.0.1 output; pass the "
                         "previous release with --carry-json)")
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
    return Loaded(doc, md, declared, hidden, identities, ratios, quantiles, dependents, status, values, primary_new,
                  input_status=dict(status))


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
        if st == "known":                  # v1.4: a --public cell keeps its label (value known from outside)
            st = L.public_label[ptr]
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
    lab = {c: (L.public_label[c] if status[c] == "known" else status[c]) for c in L.count_fields}
    n_p = sum(1 for c in L.count_fields if lab[c] == "P")
    n_s = sum(1 for c in L.count_fields if lab[c] == "S")
    sup = dict(doc.get("suppression", {}))
    sup.update({"primary_cells": n_p, "complementary_cells": n_s, "recoverable": 0,
                "recoverability_cells_checked": res.checked})
    doc["suppression"] = sup
    if profile == "stage2":
        doc["spec"] = OUTPUT_SPEC
        notes = doc.setdefault("notes", {})
        notes["differencing_fix_v1_0_3"] = (
            "v1.0.3 re-releases the v1.0.1 output with complementary suppression closed against differencing "
            "(totals minus published parts, within and across tables), now also under the declared age-band, DOB and "
            "care-unit links, public table sizes and other values already public; every cell suppressed or withheld "
            "in v1.0.2 stays suppressed or withheld; every published value is unchanged from v1.0.1; cells newly "
            "marked S and derived values newly marked suppressed are listed by pointer in the release summary")
    if L.drop_subtrees:
        notes = doc.setdefault("notes", {})
        notes["dropped_subtrees_v1_0_3"] = (
            f"{len(L.drop_subtrees)} subtree(s) withheld from the v1.0.3 release as a whole (hash-pinned drop list; "
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
        if old == "known":                 # v1.4: a --public cell must keep its label
            want = L.public_label[ptr]
            got = "pub" if isinstance(v, int) else ("P" if v == "<10" else "S")
            if got != want or (got == "pub" and v != L.values.get(ptr)):
                fails.append(f"public cell's label or value changed at {ptr}")
            old = want
        if isinstance(v, int):
            if v < MIN_CELL:
                fails.append(f"published count below 10 at {ptr}")
            if old != "pub" or v != L.values.get(ptr):
                fails.append(f"published value differs from v1.0.1 or was suppressed there (or in the carried "
                             f"release), at {ptr}")
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
    published_label = {c for c, st in status.items() if st == "pub"}
    for c in L.public_label:               # v1.4: known constants in the audit, whatever their label
        if c in status:
            status[c] = "known"
    for c in L.carried:                    # v1.4: explicit carry check (also implied by the input statuses)
        if status.get(c) == "pub":
            fails.append(f"a cell suppressed in the carried release is published at {c}")
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
    pubv = {c: L.values[c] for c in published_label}
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
            fails.append("markdown on disk is not the rendering of the v1.0.3 JSON")
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
    L = [f"# Stage 2 {RELEASE_LABEL} release summary ({info['mode']})", "",
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
         f"- frozen cells: {info['n_frozen']} (file sha256 {info['frozen_sha256']})",
         f"- public (already known) cells: {info.get('n_public', 0)} (file sha256 {info.get('public_sha256', 'none')}; "
         f"team-internal file, pointers and values not listed)",
         f"- carried release: {info.get('carry', 'none')}", ""]
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
    L += [f"## Cells by rule ({RELEASE_LABEL})", "",
          "| rule | cells |", "|---|---|",
          f"| 1 primary, carried from v1.0.1 ('<10') | {r['primary_carried']} |",
          f"| 1 primary, new (published count under 10 found in the input) | {r['primary_new']} |",
          f"| 2 complementary, carried from v1.0.1 (S) | {r['complementary_carried']} |",
          f"| 2 complementary, carried from the previous release (S there, published in v1.0.1) | "
          f"{r.get('complementary_carried_release', 0)} |",
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
def withhold_proposal(system: System, status: dict, res: AuditResult, L: Loaded, frozen: set[str],
                      link_cells: set[str]) -> tuple[dict[str, str], list[str], int]:
    """v1.5: when suppression does not converge, name single cells whose withholding (--withhold) removes the
    remaining failures, as the cell-level alternative to dropping subtrees (D4). Per failing cell u:
      * u labelled S (published in the input, suppressed): withhold u itself; withheld, it is audited as an
        unpublished cell (rule 3's hidden clause: passes when its equations prove it >= 10, or it is still wide);
      * u '<10' or hidden: withhold the published / S count cells of u's own groups (frees u's equations).
    Never a frozen, public or links-file cell (those failing cells are listed as 'blocked' for Vinai). The proposal
    is checked with a full audit (all equations) on top of the stopped state; returns (pointer -> reason, blocked
    failing cells, cells still failing with the proposal applied). Pointers only, never values."""
    cf = set(L.count_fields)

    def ok(c: str) -> bool:
        return c in cf and c not in frozen and c not in link_cells and c not in L.public_label \
            and status.get(c) in ("pub", "S")
    prop: dict[str, str] = {}
    blocked: list[str] = []
    for u in sorted(res.unsafe):
        if status.get(u) == "S" and ok(u):
            prop.setdefault(u, "v1.5 withhold proposal: suppressed cell fixed by its equations (audited as "
                               "unpublished once withheld)")
            continue
        near = sorted({c for k in system.by_cell.get(u, []) for c, _ in system.rels[k]
                       if c != u and ok(c)})
        if near:
            for c in near:
                prop.setdefault(c, "v1.5 withhold proposal: frees a '<10' / hidden cell in its group")
        else:
            blocked.append(u)
    st2 = dict(status)
    for c in prop:
        st2[c] = "D"
    beat("withhold proposal: full audit", f"{len(prop)} cells")
    res2 = system.audit(st2, label="withhold proposal check (full audit)")
    return prop, blocked, len(res2.unsafe)


def report_withhold_proposal(prop: dict[str, str], blocked: list[str], n_after: int, args) -> None:
    """Print the proposal (pointers only) and, with --withhold-proposal FILE, write it as a --withhold CSV."""
    print(f"WITHHOLD PROPOSAL: {len(prop)} cell(s); with them withheld, {n_after} cell(s) still fail the full "
          f"audit; {len(blocked)} failing cell(s) need a links-file / frozen / public cell (Vinai decides)",
          file=sys.stderr)
    for c in sorted(prop)[:20]:
        print(f"  withhold: {c}", file=sys.stderr)
    if len(prop) > 20:
        print(f"  ... and {len(prop) - 20} more" + ("" if args.withhold_proposal else
              " (use --withhold-proposal FILE for the full list)"), file=sys.stderr)
    for c in blocked[:10]:
        print(f"  blocked (needs a decision): {c}", file=sys.stderr)
    path = getattr(args, "withhold_proposal", None)
    if path:
        pth = Path(path).resolve()
        check_name(pth, "withhold proposal")
        if pth.exists():
            print(f"withhold proposal NOT written: {pth.name} already exists", file=sys.stderr)
            return
        if args.out_dir and Path(args.out_dir).resolve() in pth.parents:
            print("withhold proposal NOT written: it must not be inside --out-dir", file=sys.stderr)
            return
        lines = ["# v1.5 withhold proposal (TEAM-INTERNAL; pointers only). Review before use with --withhold.",
                 f"# with these cells withheld, {n_after} cell(s) still fail the full audit; blocked: {len(blocked)}",
                 "pointer,reason"] + [f"{c},{prop[c]}" for c in sorted(prop)]
        pth.parent.mkdir(parents=True, exist_ok=True)
        write_text(pth, "\n".join(lines) + "\n")
        print(f"withhold proposal written: {pth} ({len(prop)} row(s))", file=sys.stderr)


def run(args: argparse.Namespace) -> int:
    in_json = Path(args.in_json).resolve()
    in_md = Path(args.in_md).resolve() if args.in_md else None
    links_p = Path(args.links).resolve() if args.links else None
    frozen_p = Path(args.frozen).resolve() if args.frozen else None
    drops_p = Path(args.drops).resolve() if args.drops else None
    public_p = Path(args.public).resolve() if args.public else None
    withhold_p = Path(args.withhold).resolve() if getattr(args, "withhold", None) else None
    carry_p = Path(args.carry_json).resolve() if args.carry_json else None
    inputs = [("input_json", in_json)] + ([("input_md", in_md)] if in_md else []) + \
             ([("carry_json", carry_p)] if carry_p else []) + \
             ([("links", links_p)] if links_p else []) + ([("frozen", frozen_p)] if frozen_p else []) + \
             ([("public", public_p)] if public_p else []) + ([("drops", drops_p)] if drops_p else []) + \
             ([("withhold", withhold_p)] if withhold_p else [])
    for role, pth in inputs:
        check_name(pth, role)
        if not pth.is_file():
            raise InputError(f"{role} not found: {pth}")
    in_hash = {role: sha256_file(pth) for role, pth in inputs}
    if args.expect_input_sha256 and in_hash["input_json"] != args.expect_input_sha256.lower():
        raise InputError("input JSON SHA-256 does not match --expect-input-sha256")
    if args.expect_md_sha256 and in_hash.get("input_md") != args.expect_md_sha256.lower():
        raise InputError("input markdown SHA-256 does not match --expect-md-sha256")
    if args.expect_carry_sha256 and in_hash.get("carry_json") != args.expect_carry_sha256.lower():
        raise InputError("carry JSON SHA-256 does not match --expect-carry-sha256")
    params = Params(width=args.protection_width, s_rule=args.s_rule, s_lower=args.s_lower_bound)
    if params.width < 1:
        raise InputError("--protection-width must be >= 1")
    expect_spec = None if args.profile == "generic" else (args.expect_spec or EXPECT_INPUT_SPEC)
    if args.check_only:
        expect_spec = None  # check-only audits any registry-style output (v1.0.1, v1.0.2, follow-up)
    log(f"{SCRIPT_VERSION}: reading the input JSON ({in_json.stat().st_size:,} bytes) and checking its schema")
    L = load_output(in_json, in_md, args.profile, expect_spec, allow_released=args.check_only)
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
    withhold = read_withhold(withhold_p)
    link_cells = {c for t, parts, _ in link_rels for c in (t,) + parts}
    apply_drops(L, drops, link_cells, frozen, withhold)
    if drops or withhold:
        log(f"drop list: {len(drops)} subtree(s) and {len(withhold)} single cell(s) withheld from the release, "
            f"{len(L.dropped)} count cells marked unpublished (audited as hidden)")
    if carry_p is not None:                # v1.4: the previous release's suppressions and withheld subtrees stay
        try:
            carry_doc = json.loads(carry_p.read_text(encoding="utf-8"))
        except Exception as e:  # noqa: BLE001
            raise InputError(f"cannot read the carry JSON: {e}")
        if not isinstance(carry_doc, dict):
            raise InputError("carry JSON is not an object")
        apply_carry(L, carry_doc, in_hash["input_json"])
        for c in L.carried:
            if c in frozen:
                raise InputError(f"a frozen cell is suppressed in the carried release: {c}")
        log(f"carried release: {len(L.carried)} cell(s) published in the input stay S (suppressed there)")
    public = read_public(public_p)
    if withhold and {t for t, _ in withhold} & {ptr for ptr, _v, _r in public}:
        raise InputError("a --public cell cannot be withheld")
    apply_public(L, public, frozen, link_rels)
    if public:
        log(f"public cells: {len(public)} value(s) treated as known in every audit (labels unchanged"
            + (f"; {len(L.republished)} carried cell(s) published again: 'publish' rows)" if L.republished else ")"))
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
        "n_public": len(public), "public_sha256": in_hash.get("public", "none"),
        "carry": (f"{carry_p.name} (sha256 {in_hash['carry_json']}; {len(L.carried)} cell(s) carried as S"
                  + (f"; {len(L.republished)} carried cell(s) published again as public by definition ('publish' rows)"
                     if L.republished else "") + ")" if carry_p else "none"),
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
    sup.max_new = args.max_new_suppressions
    sup.max_component = args.max_component_cells
    log("complementary suppression: start")
    try:
        sup.run()
    except NotConverging as e:
        log(f"STOPPING: {e}")
        prop, blocked, n_after = withhold_proposal(system, status, e.res, L, frozen, link_cells)
        report_withhold_proposal(prop, blocked, n_after, args)
        raise GateError(f"{e} (new suppressions so far: {len(sup.rule_of)}); withhold proposal: {len(prop)} "
                        f"cell(s), {n_after} cell(s) still failing with them withheld, {len(blocked)} failing "
                        f"cell(s) that need a links-file / frozen / public cell (Vinai decides)")
    log(f"complementary suppression: done ({len(sup.rule_of)} new cells, {sup.rounds} rounds)")
    doc = apply_statuses(L, status)
    new_cells = sorted((c, sup.rule_of[c]) for c in sup.rule_of)
    new_cells += sorted((c, "primary_new") for c in L.primary_new)
    prov = {"source_json": in_json.name, "source_json_sha256": in_hash["input_json"],
            "source_md_sha256": in_hash.get("input_md", "none"), "script": SCRIPT_VERSION,
            "script_sha256": info["script_sha256"], "spec": SPEC_DOC,
            "links_sha256": in_hash.get("links", "none"), "frozen_sha256": in_hash.get("frozen", "none"),
            "public_sha256": in_hash.get("public", "none"), "n_public_cells": len(public),
            "carry_json": carry_p.name if carry_p else "none", "carry_json_sha256": in_hash.get("carry_json", "none"),
            "carried_suppressions": len(L.carried), "republished_public_cells": len(L.republished),
            "protection_width": params.width, "s_rule": params.s_rule, "s_lower_bound": params.s_lower,
            "newly_suppressed_cells": len(new_cells)}
    if withhold:
        prov["withhold_sha256"] = in_hash["withhold"]
        prov["withheld_cells"] = len(withhold)
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
            "rules": {"primary_carried": sum(1 for c in L.count_fields if L.input_status[c] == "P"
                                             and c not in set(L.dropped)) - len(L.primary_new),
                      "primary_new": len(L.primary_new),
                      "complementary_carried": sum(1 for c in L.count_fields if L.input_status[c] == "S"
                                                   and c not in set(L.dropped)),
                      "complementary_carried_release": len(L.carried),
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
        if withhold:
            manifest["withheld_cells"] = [{"pointer": t, "reason": why} for t, why in withhold]
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
    ap = argparse.ArgumentParser(description="Stage 2 v1.0.3 differencing-fix release (aggregate output only).")
    ap.add_argument("--in-json", required=True, help="released Stage 2 v1.0.1 JSON (feature_sparsity_frequency_audit.json)")
    ap.add_argument("--in-md", help="released Stage 2 v1.0.1 markdown (required for --profile stage2)")
    ap.add_argument("--links", help="CSV of cross-table equations (total, parts, note)")
    ap.add_argument("--frozen", help="CSV of already-public cells that are never suppressed (pointer, reason)")
    ap.add_argument("--drops", help="CSV of JSON subtrees withheld from the release as a whole (subtree, reason)")
    ap.add_argument("--public", help="v1.4: team-internal CSV of cells whose value a reader already knows from outside "
                                     "the release (pointer, value, reason); known constants in every audit")
    ap.add_argument("--carry-json", help="v1.4: the previous release (v1.0.2 JSON); its suppressions and withheld "
                                         "subtrees are kept")
    ap.add_argument("--expect-carry-sha256", help="optional exact SHA-256 pin of the carry JSON")
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
    ap.add_argument("--withhold", help="v1.5: CSV of single count cells withheld from the release (pointer, reason)")
    ap.add_argument("--withhold-proposal", help="v1.5: on a convergence stop, write the withhold proposal (pointers "
                                                "only) to this new CSV file (outside --out-dir)")
    ap.add_argument("--max-new-suppressions", type=int, default=MAX_NEW_DEFAULT,
                    help=f"v1.5: stop (exit 3, withhold proposal) above this many new suppressions "
                         f"(default {MAX_NEW_DEFAULT})")
    ap.add_argument("--max-component-cells", type=int, default=MAX_COMPONENT_DEFAULT,
                    help=f"v1.5: stop (exit 3, withhold proposal) when a component exceeds this many unknown "
                         f"cells (default {MAX_COMPONENT_DEFAULT})")
    ap.add_argument("--heartbeat-seconds", type=float, default=HEARTBEAT_SECONDS,
                    help=f"v1.5: heartbeat line after this many seconds without a line (default {HEARTBEAT_SECONDS:g})")
    ap.add_argument("--expect-spec", help="exact spec string required in the input (default: Stage 2 v1.0.1)")
    ap.add_argument("--expect-input-sha256", help="optional exact SHA-256 pin of the input JSON")
    ap.add_argument("--expect-md-sha256", help="optional exact SHA-256 pin of the input markdown")
    a = ap.parse_args(argv)
    if a.max_seconds <= 0 or a.workers < 1:
        ap.error("--max-seconds must be > 0 and --workers >= 1")
    if a.check_only and a.drops:
        ap.error("--drops is used only for a release run, not with --check-only")
    if a.check_only and (a.withhold or a.withhold_proposal):
        ap.error("--withhold / --withhold-proposal are used only for a release run, not with --check-only")
    if a.max_new_suppressions < 1 or a.max_component_cells < 1 or a.heartbeat_seconds <= 0:
        ap.error("--max-new-suppressions, --max-component-cells and --heartbeat-seconds must be > 0")
    if a.check_only and a.carry_json:
        ap.error("--carry-json is used only for a release run, not with --check-only")
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
    _LAST_LINE[0] = _T0[0]
    _HB.update(phase="starting", detail="", lp=0)
    try:
        with Heartbeat(args.heartbeat_seconds):
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
