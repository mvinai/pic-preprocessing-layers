#!/usr/bin/env python3
"""Stage 2 follow-up checks (feature_sparsity_frequency_audit_spec v0.9, section 5.4; decisions D24 and D29).

Read-only and aggregate-only. Two checks on the train split, inside the Stage 2 analysis window [t0, min(OUT, t0 + 336 h)),
with lab availability = chart time + 1 h (as Stage 2):

(i)  Bilirubin 5255 vs 5075 (D24): for all train encounters, NICU vs non-NICU (first care unit) and neonate vs
     not-neonate (age < 28 days at t0; age_invalid counts as not-neonate): per item n values, n encounters, value
     quantiles, per-encounter-median quantiles; per group KS 5255 vs 5075 (only when both n >= 10), encounters with
     both items, same-draw (same encounter and CHARTTIME) pairs and the same-draw difference 5255 - 5075.
(ii) FiO2 5222 format (D29): value bands (<= 1, (1, 21), [21, 100], > 100) of all train-window FiO2 rows;
     same-draw pO2 (5239 preferred over 5244, measurement rows) and numeric FiO2 pairs; P/F under the percent rule
     (P/F = pO2 / (FiO2 / 100)) split in gate / below / above the derived 10-700 gate (FiO2 <= 0 kept apart); how many
     out-of-gate pairs fall in the gate with FiO2 x100 or /100; P/F under the rule "FiO2 <= 1 is a proportion,
     otherwise a percentage"; pO2 and FiO2 quantiles for in-gate vs out-of-gate pairs; out-of-gate pairs by FiO2 band;
     and the Stage 2 v1.0.1 pair definition (FiO2 gated 21-100) as a reproduction check.

It reuses the pinned Stage 2 script (feature_sparsity_frequency_audit.py v1.0.1) and, through it, the duration audit
v1.7 module: cohort, split and harmonization gates, the event reader, the gates, the <10 primary plus complementary
suppression with the LP recoverability check, and the same pending -> gate -> final -> gate publication.
Every input is checked against an exact SHA-256 before use. Outputs: stage2_followup_checks.json and
stage2_followup_checks.md (both shareable, suppressed). No unsuppressed file is written.

Exit codes: 0 ok; 2 input / hash / cohort / reconciliation stop (AUDIT STOPPED); 3 release gate failed.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

SCRIPT_VERSION = "stage2_followup_checks v1.0 (feature_sparsity_frequency_audit_spec v0.9, section 5.4; D24, D29)"
STAGE2_FILE = "feature_sparsity_frequency_audit.py"
STAGE2_SHA256 = "ce03f97f22f93749c7d325fa8411326ca575945386134e7fcd30089e9180c902"  # public copy (cohort pins omitted); original private pin was 3b5b36fbd0f73420beaa74b18abf2afec4725cf86c5e8efe514be19120043213
STAGE2_SPEC_VERSION = "feature_sparsity_frequency_audit v1.0.1 (feature_sparsity_frequency_audit_spec v0.8, stage 2)"
CHANNEL_MAP_V2_FILE = "feature_channel_map_v2.csv"
CHANNEL_MAP_V2_SHA256 = "cbfdc68f11bc2f8d9f96150803a1339930da81eb091f481576f709debdd5ab5c"
CHANNEL_MAP_V2_ROWS = 686
# D36 (Nov 15, 2026): the four tables Stage 1 did not hash are now gated by these SHA-256 values (Stage 2 v1.0.1 record)
RAW_D36_SHA256 = {
    "PATIENTS": "85b016fcc3050e5b5b6200103ed74ad1798720c881e26c60f1769be9f336ae0d",
    "ADMISSIONS": "1451ae8ea0ff463269b348522047a4bc99efd2665c770e9c5cef0ffe2647e1ae",
    "ICUSTAYS": "449119e698acaaf948921c34375db535b9fd3f6eb6e996a86069da95e90eea1b",
    "SURGERY_INFO": "04450bff4592d676a7ebe2dad22d9e6c2a6f8edd07dff223b182eb0002a38d83",
}
LABEVENTS_STAGE1_SHA256 = "784ecb929a358ee1a18fc9beafbdf67f30a91646b780d9983714a5e9ce0018c5"
BILI_MAIN, BILI_ALT = "5075", "5255"
FIO2, PO2_MAIN, PO2_ALT = "5222", "5239", "5244"
WANTED = {BILI_MAIN, BILI_ALT, FIO2, PO2_MAIN, PO2_ALT}
EXPECTED_MAP_ROWS = {  # (table, itemid): (channel, status) in channel map v2
    ("LABEVENTS", BILI_MAIN): ("bilirubin_total", "assigned_merge_pending_stage2"),
    ("LABEVENTS", BILI_ALT): ("bilirubin_total", "assigned_merge_pending_stage2"),
    ("LABEVENTS", FIO2): ("fio2", "assigned"),
    ("LABEVENTS", PO2_MAIN): ("po2", "assigned_merged"),
    ("LABEVENTS", PO2_ALT): ("po2", "assigned_merged"),
    ("DERIVED", "derived:pf_ratio"): ("pf_ratio", "derived"),
}
EXPECTED_STAGE2_PAIRS = {"valid": 12149, "outside_gate": 4910}  # Stage 2 v1.0.1 map_checks/pf_ratio/same_draw_pairs (train)
FIO2_BANDS = ("le_1", "gt_1_lt_21", "21_to_100", "gt_100")
NICU = "NICU"
NON_COUNT_TOP = ("spec", "release_gate", "input_hashes", "stage2_module", "v17_module", "maps", "raw_tables",
                 "cohort_checks", "config", "notes", "disclosure", "suppression")
OUT_JSON = "stage2_followup_checks.json"
OUT_MD = "stage2_followup_checks.md"


class FollowupError(Exception):
    pass


# ============================================================================= pinned Stage 2 module
def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_stage2(path: Path, expected_sha256: str = STAGE2_SHA256):
    path = Path(path)
    if not path.exists():
        raise FollowupError(f"Stage 2 script not found: {path} (it must sit next to this script)")
    if not re.fullmatch(r"[0-9a-f]{64}", str(expected_sha256)):
        raise FollowupError("pinned Stage 2 script hash is not a full 64-character SHA-256")
    got = _sha256(path)
    if got != expected_sha256:
        raise FollowupError(f"Stage 2 script SHA-256 mismatch (expected {expected_sha256}, got {got}); re-copy it, do not edit it")
    name = f"_pinned_stage2_{got[:12]}"
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
    if getattr(mod, "SPEC_VERSION", None) != STAGE2_SPEC_VERSION:
        raise FollowupError("Stage 2 module SPEC_VERSION differs from the pinned value")
    return mod


# ============================================================================= config
@dataclass
class FUConfig:
    raw_dir: Path
    canonical_path: Path
    harmonization_path: Path
    split_path: Path
    complexity_map_path: Path
    maps_dir: Path
    out_dir: Path
    stage2_path: Path
    v17_module_path: Path
    harm_t0_col: str | None = None
    harm_out_col: str | None = None
    split_col: str = "split"
    chunksize: int = 1_000_000
    verbose: bool = True
    stage2_sha256: str = STAGE2_SHA256
    channel_map_v2_sha256: str = CHANNEL_MAP_V2_SHA256
    channel_map_v2_rows: int = CHANNEL_MAP_V2_ROWS
    raw_sha256: dict = field(default_factory=lambda: {"LABEVENTS": LABEVENTS_STAGE1_SHA256, **RAW_D36_SHA256})
    expected_stage2_pairs: dict = field(default_factory=lambda: dict(EXPECTED_STAGE2_PAIRS))
    s2_overrides: dict = field(default_factory=dict)  # tests only (synthetic pins and expected counts)

    @property
    def shared_json(self) -> Path:
        return self.out_dir / OUT_JSON

    @property
    def shared_md(self) -> Path:
        return self.out_dir / OUT_MD


def log(cfg: FUConfig, msg: str) -> None:
    if cfg.verbose:
        print(msg, flush=True)


def make_s2config(M, cfg: FUConfig):
    kw = dict(raw_dir=Path(cfg.raw_dir), canonical_path=Path(cfg.canonical_path),
              harmonization_path=Path(cfg.harmonization_path), split_path=Path(cfg.split_path),
              complexity_map_path=Path(cfg.complexity_map_path), maps_dir=Path(cfg.maps_dir), out_dir=Path(cfg.out_dir),
              v17_module_path=Path(cfg.v17_module_path), harm_t0_col=cfg.harm_t0_col, harm_out_col=cfg.harm_out_col,
              split_col=cfg.split_col, chunksize=cfg.chunksize, verbose=cfg.verbose)
    kw.update(cfg.s2_overrides)
    kw["raw_sha256"] = dict(cfg.raw_sha256)
    return M.S2Config(**kw)


# ============================================================================= channel map v2
def load_channel_map_v2(M, cfg: FUConfig) -> tuple[dict, str]:
    p = Path(cfg.maps_dir) / CHANNEL_MAP_V2_FILE
    h = M.verify_hash("feature_channel_map_v2", p, cfg.channel_map_v2_sha256)
    cm = pd.read_csv(p, dtype=str, keep_default_na=False, encoding="utf-8")
    if len(cm) != cfg.channel_map_v2_rows:
        raise FollowupError(f"channel map v2 has {len(cm)} rows, expected {cfg.channel_map_v2_rows}")
    need = {"table", "itemid", "channel", "plausibility_low", "plausibility_high", "status", "mask_group"}
    if need - set(cm.columns):
        raise FollowupError(f"channel map v2 missing columns: {sorted(need - set(cm.columns))}")
    gates = {}
    for key, (ch, st) in EXPECTED_MAP_ROWS.items():
        r = cm[(cm["table"].str.strip() == key[0]) & (cm["itemid"].str.strip() == key[1])]
        if len(r) != 1:
            raise FollowupError(f"channel map v2: expected one row for {key}, found {len(r)}")
        r = r.iloc[0]
        if r["channel"].strip() != ch or r["status"].strip() != st:
            raise FollowupError(f"channel map v2 row {key}: channel/status {r['channel']}/{r['status']}, expected {ch}/{st}")
        lo, hi = r["plausibility_low"].strip(), r["plausibility_high"].strip()
        gates[key[1]] = {"itemid": key[1], "low": float(lo) if lo else None, "high": float(hi) if hi else None}
    return gates, h


# ============================================================================= helpers
def _q_and_n(out, base: tuple, name: str, vals: np.ndarray, equal_to: str | None = None, parent: str | None = None) -> str:
    return out.q(base, name, np.asarray(vals, dtype=float), parent=parent, equal_to=equal_to)


def fio2_band(f: np.ndarray) -> np.ndarray:
    f = np.asarray(f, dtype=float)
    return np.select([f <= 1, f < 21, f <= 100], ["le_1", "gt_1_lt_21", "21_to_100"], "gt_100")


def pf_category(pf: np.ndarray, fio2: np.ndarray, lo: float, hi: float) -> np.ndarray:
    pf, fio2 = np.asarray(pf, dtype=float), np.asarray(fio2, dtype=float)
    return np.select([fio2 <= 0, pf < lo, pf > hi], ["fio2_le_0", "below_gate", "above_gate"], "in_gate")


def ks_stat(a: np.ndarray, b: np.ndarray, M) -> float | None:
    from scipy import stats as sps
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < M.MIN_CELL or len(b) < M.MIN_CELL:
        return None
    return float(sps.ks_2samp(a, b).statistic)


# ============================================================================= (i) bilirubin
def bilirubin_block(out, M, bili: pd.DataFrame, groupings: dict, train_n: dict) -> None:
    """bili: train-window measurement rows of 5075 / 5255 (enc_i, tc, item, v). train_n: frozen train group sizes."""
    def group_cells(base: tuple, emask: np.ndarray, n_enc: int) -> dict:
        c = {"n_encounters": out.cell(base + ("n_encounters",), n_enc)}
        sel = bili[emask[bili["enc_i"].to_numpy()]]
        vals = {}
        for it in (BILI_MAIN, BILI_ALT):
            g = sel[sel["item"] == it]
            b = base + ("items", it)
            c[f"v{it}"] = _q_and_n(out, b, "values", g["v"].to_numpy())
            med = g.groupby("enc_i")["v"].median()
            c[f"e{it}"] = _q_and_n(out, b, "encounter_medians", med.to_numpy(), parent=c["n_encounters"])
            vals[it] = (g["v"].to_numpy(), med.to_numpy(), set(med.index))
        both = vals[BILI_MAIN][2] & vals[BILI_ALT][2]
        either = vals[BILI_MAIN][2] | vals[BILI_ALT][2]
        c["both"] = out.cell(base + ("n_encounters_both_items",), len(both), parent=c[f"e{BILI_MAIN}"], tag="sub5075")
        out.reg.subset(c["both"], c[f"e{BILI_ALT}"], "sub5255")
        c["either"] = out.cell(base + ("n_encounters_either_item",), len(either), parent=c["n_encounters"])
        hsum = out.reg.hide(out.A.P(*base, "sum_item_encounters"), len(vals[BILI_MAIN][2]) + len(vals[BILI_ALT][2]))
        out.ident(hsum, [c[f"e{BILI_MAIN}"], c[f"e{BILI_ALT}"]])
        out.ident(hsum, [c["either"], c["both"]])
        a = sel[sel["item"] == BILI_MAIN].drop_duplicates(["enc_i", "tc"])
        b2 = sel[sel["item"] == BILI_ALT].drop_duplicates(["enc_i", "tc"])
        sd = a[["enc_i", "tc", "v"]].merge(b2[["enc_i", "tc", "v"]], on=["enc_i", "tc"], suffixes=("_5075", "_5255"))
        c["pairs"] = _q_and_n(out, base, "same_draw_pairs_difference_5255_minus_5075", (sd["v_5255"] - sd["v_5075"]).to_numpy())
        out.reg.subset(c["pairs"], c[f"v{BILI_MAIN}"], "pair5075")
        out.reg.subset(c["pairs"], c[f"v{BILI_ALT}"], "pair5255")
        out.ks(base, "values_5255_vs_5075", [c[f"v{BILI_MAIN}"], c[f"v{BILI_ALT}"]],
               ks_stat(vals[BILI_ALT][0], vals[BILI_MAIN][0], M))
        out.ks(base, "encounter_medians_5255_vs_5075", [c[f"e{BILI_MAIN}"], c[f"e{BILI_ALT}"]],
               ks_stat(vals[BILI_ALT][1], vals[BILI_MAIN][1], M))
        return c

    allm = groupings["all_train"]
    top = group_cells(("bilirubin", "all_train"), allm, int(allm.sum()))
    out.frozen.add(top["n_encounters"])
    for gname, parts in groupings["splits"].items():
        cs = {}
        for lvl, m in parts.items():
            cs[lvl] = group_cells(("bilirubin", gname, lvl), m, int(m.sum()))
            if lvl in train_n.get(gname, {}) and int(m.sum()) >= M.MIN_CELL:  # public in Stage 2 v1.0.1 (if >= 10)
                out.frozen.add(cs[lvl]["n_encounters"])
        for k in top:
            out.ident(top[k], [cs[l][k] for l in parts])


# ============================================================================= (ii) FiO2 format
def fio2_block(out, M, fi_rows: pd.DataFrame, pairs: pd.DataFrame, s2pairs: pd.DataFrame, pfg: dict,
               expected_pairs: dict) -> dict:
    lo, hi = pfg["low"], pfg["high"]
    b = ("fio2_format",)
    # rows
    num = fi_rows["vnum"].notna().to_numpy()
    rows = out.partition(b + ("rows_5222",), {"numeric": int(num.sum()), "nonnumeric": int((~num).sum())})
    bands = fio2_band(fi_rows.loc[num, "vnum"].to_numpy())
    out.partition(b + ("rows_5222", "numeric_by_band"), {k: int((bands == k).sum()) for k in FIO2_BANDS}, total=rows["numeric"])
    # pairs, percent rule (Stage 2 formula)
    f = pairs["fio2"].to_numpy(dtype=float)
    po2 = pairs["po2"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        pf_pct = np.where(f > 0, po2 / (f / 100.0), np.inf)
    cat = pf_category(pf_pct, f, lo, hi)
    bp = b + ("pairs_percent_rule",)
    ntot = out.cell(bp + ("n_pairs",), len(pairs))
    top = out.partition(bp, {"in_gate": int((cat == "in_gate").sum()), "out_of_gate": int((cat != "in_gate").sum())}, total=ntot)
    oog = cat != "in_gate"
    out.partition(bp + ("out_of_gate",), {k: int((cat == k).sum()) for k in ("below_gate", "above_gate", "fio2_le_0")},
                  total=top["out_of_gate"])
    pb = fio2_band(f)
    out.partition(bp + ("pairs_by_fio2_band",), {k: int((pb == k).sum()) for k in FIO2_BANDS}, total=ntot)
    out.partition(bp + ("out_of_gate", "by_fio2_band"), {k: int((pb[oog] == k).sum()) for k in FIO2_BANDS}, total=top["out_of_gate"])
    with np.errstate(divide="ignore", invalid="ignore"):
        x100 = np.where(f > 0, po2 / f, np.inf)              # FiO2 read as a proportion (x100 to percent)
        d100 = np.where(f > 0, po2 / (f / 10000.0), np.inf)  # FiO2 read as percent x100 (/100)
    ing = lambda v: (v >= lo) & (v <= hi)  # noqa: E731
    out.cell(bp + ("out_of_gate", "n_in_gate_if_fio2_times_100"), int((oog & ing(x100)).sum()), parent=top["out_of_gate"], tag="x100")
    out.cell(bp + ("out_of_gate", "n_in_gate_if_fio2_divided_by_100"), int((oog & ing(d100)).sum()), parent=top["out_of_gate"], tag="d100")
    for lvl, m in (("in_gate", ~oog), ("out_of_gate", oog)):
        _q_and_n(out, bp + (lvl,), "po2", po2[m], equal_to=top[lvl])
        _q_and_n(out, bp + (lvl,), "fio2", f[m], equal_to=top[lvl])
    _q_and_n(out, bp + ("in_gate",), "pf", pf_pct[~oog], equal_to=top["in_gate"])
    # proportion rule: FiO2 <= 1 is a proportion, otherwise a percentage
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(f <= 1, f, f / 100.0)
        pf_r = np.where(f > 0, po2 / frac, np.inf)
    cat_r = pf_category(pf_r, f, lo, hi)
    br = b + ("pairs_proportion_rule",)
    rr = out.partition(br, {"computable": int((f > 0).sum()), "fio2_le_0": int((f <= 0).sum())}, total=ntot)
    rr2 = out.partition(br + ("computable",), {k: int((cat_r == k).sum()) for k in ("in_gate", "below_gate", "above_gate")},
                        total=rr["computable"])
    _q_and_n(out, br + ("computable",), "pf", pf_r[f > 0], equal_to=rr["computable"])
    _q_and_n(out, br + ("computable", "in_gate"), "pf", pf_r[cat_r == "in_gate"], equal_to=rr2["in_gate"])
    # Stage 2 v1.0.1 definition (FiO2 measurement rows gated 21-100)
    bs = b + ("stage2_definition",)
    sv = s2pairs["valid"].to_numpy(dtype=bool)
    got = {"valid": int(sv.sum()), "outside_gate": int((~sv).sum())}
    out.partition(bs, got)
    out.val(bs + ("expected_from_stage2_v1_0_1",), dict(expected_pairs))
    out.val(bs + ("reproduces_stage2_v1_0_1",), bool(got == {k: int(v) for k, v in expected_pairs.items()}))
    return got


# ============================================================================= document, markdown, gate
def _g(doc: dict, *path) -> Any:
    cur = doc
    for p in path:
        if not isinstance(cur, dict) or p not in cur:
            return None
        cur = cur[p]
    return cur


def _qf(v: Any) -> str:
    if isinstance(v, dict):
        return " / ".join(str(v.get(k)) for k in ("p10", "p50", "p90"))
    return "–" if v is None else str(v)


def render_markdown(doc: dict) -> str:
    L = ["# Stage 2 follow-up checks (bilirubin 5255 vs 5075; FiO2 5222 format)", "",
         f"- spec: {doc['spec']}", f"- release gate: {doc['release_gate']}",
         "- train split, analysis window [t0, min(OUT, t0 + 336 h)), labs available at chart time + 1 h",
         "- counts 1-9 shown as <10; S = complementary suppression; quantiles are p10 / p50 / p90", "",
         "## Input hashes (all gated)", "", "| input | sha256 |", "|---|---|"]
    for k, v in doc["input_hashes"].items():
        L.append(f"| {k} | {v} |")
    L += ["", "## (i) Bilirubin 5255 vs 5075 (D24)", "",
          "| grouping | group | n_encounters | n_encounters_either_item | n_encounters_both_items | n_same_draw_pairs | "
          "same-draw diff 5255-5075 | KS values | KS encounter medians |", "|---|---|---|---|---|---|---|---|---|"]
    groups = [("all", "all_train", ("bilirubin", "all_train"))]
    for gname in ("by_first_careunit", "by_age"):
        for lvl in (_g(doc, "bilirubin", gname) or {}):
            groups.append((gname, lvl, ("bilirubin", gname, lvl)))
    for gname, lvl, path in groups:
        d = _g(doc, *path) or {}
        L.append(f"| {gname} | {lvl} | {d.get('n_encounters', '–')} | {d.get('n_encounters_either_item', '–')} | "
                 f"{d.get('n_encounters_both_items', '–')} | {d.get('n_same_draw_pairs_difference_5255_minus_5075', '–')} | "
                 f"{_qf(d.get('q_same_draw_pairs_difference_5255_minus_5075'))} | {d.get('ks_values_5255_vs_5075', '–')} | "
                 f"{d.get('ks_encounter_medians_5255_vs_5075', '–')} |")
    L += ["", "| grouping | group | item | n_values | values (umol/L) | n_encounter_medians | encounter medians (umol/L) |",
          "|---|---|---|---|---|---|---|"]
    for gname, lvl, path in groups:
        for it in (BILI_MAIN, BILI_ALT):
            d = _g(doc, *path, "items", it) or {}
            L.append(f"| {gname} | {lvl} | {it} | {d.get('n_values', '–')} | {_qf(d.get('q_values'))} | "
                     f"{d.get('n_encounter_medians', '–')} | {_qf(d.get('q_encounter_medians'))} |")
    fo = doc.get("fio2_format", {})
    L += ["", "## (ii) FiO2 5222 format (D29)", "", "### FiO2 rows (train window)", "",
          "| n_total | n_numeric | n_nonnumeric | n_le_1 | n_gt_1_lt_21 | n_21_to_100 | n_gt_100 |", "|---|---|---|---|---|---|---|"]
    r = fo.get("rows_5222", {})
    nb = r.get("numeric_by_band", {})
    L.append(f"| {r.get('n_total')} | {r.get('n_numeric')} | {r.get('n_nonnumeric')} | {nb.get('n_le_1')} | "
             f"{nb.get('n_gt_1_lt_21')} | {nb.get('n_21_to_100')} | {nb.get('n_gt_100')} |")
    p = fo.get("pairs_percent_rule", {})
    og = p.get("out_of_gate", {})
    L += ["", "### Same-draw pO2 / FiO2 pairs, percent rule P/F = pO2 / (FiO2 / 100), gate 10-700", "",
          "| n_pairs | n_in_gate | n_out_of_gate | n_below_gate | n_above_gate | n_fio2_le_0 | n_in_gate_if_fio2_times_100 | "
          "n_in_gate_if_fio2_divided_by_100 |", "|---|---|---|---|---|---|---|---|",
          f"| {p.get('n_pairs')} | {p.get('n_in_gate')} | {p.get('n_out_of_gate')} | {og.get('n_below_gate')} | "
          f"{og.get('n_above_gate')} | {og.get('n_fio2_le_0')} | {og.get('n_in_gate_if_fio2_times_100')} | "
          f"{og.get('n_in_gate_if_fio2_divided_by_100')} |", "",
          "| pairs | n_le_1 | n_gt_1_lt_21 | n_21_to_100 | n_gt_100 |", "|---|---|---|---|---|"]
    for lab, d in (("all pairs", p.get("pairs_by_fio2_band", {})), ("out of gate", og.get("by_fio2_band", {}))):
        L.append(f"| {lab} | {d.get('n_le_1')} | {d.get('n_gt_1_lt_21')} | {d.get('n_21_to_100')} | {d.get('n_gt_100')} |")
    L += ["", "| pairs | n_po2 | pO2 (mmHg) | n_fio2 | FiO2 (as recorded) |", "|---|---|---|---|---|"]
    for lvl in ("in_gate", "out_of_gate"):
        d = p.get(lvl, {})
        L.append(f"| {lvl} | {d.get('n_po2')} | {_qf(d.get('q_po2'))} | {d.get('n_fio2')} | {_qf(d.get('q_fio2'))} |")
    L.append(f"\n- P/F in gate (percent rule): n_pf {_g(p, 'in_gate', 'n_pf')}, quantiles {_qf(_g(p, 'in_gate', 'q_pf'))}")
    pr = fo.get("pairs_proportion_rule", {})
    c = pr.get("computable", {})
    L += ["", "### Proportion rule (FiO2 <= 1 is a proportion, otherwise a percentage)", "",
          "| n_computable | n_fio2_le_0 | n_in_gate | n_below_gate | n_above_gate | P/F all computable | P/F in gate |",
          "|---|---|---|---|---|---|---|",
          f"| {pr.get('n_computable')} | {pr.get('n_fio2_le_0')} | {c.get('n_in_gate')} | {c.get('n_below_gate')} | "
          f"{c.get('n_above_gate')} | {_qf(c.get('q_pf'))} | {_qf(_g(c, 'in_gate', 'q_pf'))} |"]
    s = fo.get("stage2_definition", {})
    L += ["", "### Stage 2 v1.0.1 pair definition (FiO2 gated 21-100)", "",
          "| n_total | n_valid | n_outside_gate | reproduces v1.0.1 |", "|---|---|---|---|",
          f"| {s.get('n_total')} | {s.get('n_valid')} | {s.get('n_outside_gate')} | {s.get('reproduces_stage2_v1_0_1')} |", "",
          "## Notes", ""]
    for k, v in doc.get("notes", {}).items():
        L.append(f"- {k}: {v}")
    L += ["", f"Suppression: {doc['suppression'].get('primary_cells')} primary, {doc['suppression'].get('complementary_cells')} "
              f"complementary, recoverable {doc['suppression'].get('recoverable')}.", ""]
    return "\n".join(L)


def detect_count_fields(M, obj: Any, prefix: str = "") -> list[str]:
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if prefix == "" and k in NON_COUNT_TOP:
                continue
            ptr = prefix + "/" + str(k).replace("~", "~0").replace("/", "~1")
            if M.COUNT_KEY_RE.match(str(k)) and not isinstance(v, (dict, list)):
                out.append(ptr)
            else:
                out.extend(detect_count_fields(M, v, ptr))
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


def release_gate(M, A, json_path: Path, md_path: Path, expect_status: str) -> list[str]:
    """Same checks as the Stage 2 release gate, on this script's document."""
    fails: list[str] = []
    try:
        doc = json.loads(Path(json_path).read_text(encoding="utf-8"))
        md = Path(md_path).read_text(encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        return [f"cannot reload shared files: {e}"]
    if doc.get("spec") != SCRIPT_VERSION or doc.get("release_gate") != expect_status:
        fails.append("header fields (spec / release_gate) wrong")
    for p in (json_path, md_path):
        if "_raw_LOCAL" in Path(p).name:
            fails.append("a shared file name contains _raw_LOCAL")
    dis = doc.get("disclosure", {})
    declared = set(dis.get("count_fields", []))
    detected = set(detect_count_fields(M, doc))
    if declared != detected:
        fails.append(f"count-field schema mismatch: {len(detected - declared)} undeclared, {len(declared - detected)} missing")
    known: dict[str, int] = {}
    unknown: dict[str, str] = {}
    n_p = n_s = 0
    for ptr in sorted(detected):
        v = A.get_ptr(doc, ptr)
        if isinstance(v, bool) or not (isinstance(v, int) or v in ("<10", "S")):
            fails.append(f"count field with an invalid value at {ptr}")
        elif isinstance(v, int):
            if v < M.MIN_CELL:
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
    qdecl = {p: n for p, n in dis.get("quantiles", [])}
    for ptr in detect_prefixed(doc, ("q_",)):
        v = A.get_ptr(doc, ptr)
        sib = ptr.rsplit("/", 1)[0] + "/n_" + ptr.rsplit("/", 1)[1][2:]
        if qdecl.get(ptr) != sib:
            fails.append(f"quantile not declared or not paired with its own count at {ptr}")
            continue
        if v == "suppressed":
            continue
        if not isinstance(v, dict) or set(v) != set(M.QKEYS):
            fails.append(f"quantile with disallowed keys at {ptr}")
        if known.get(sib, 0) < M.MIN_CELL:
            fails.append(f"quantile published without a published count >= 10 at {ptr}")
    ddecl = {p: nd for p, nd in dis.get("dependents", [])}
    for ptr in detect_prefixed(doc, ("ks_", "x_", "pct_")):
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
                if M.TIMESTAMP_RE.search(str(k)):
                    fails.append(f"timestamp-like key at {path}")
                walk(v, f"{path}/{k}")
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}/{i}")
        elif isinstance(o, str) and M.TIMESTAMP_RE.search(o):
            fails.append(f"timestamp-like string at {path}")
    walk({k: v for k, v in doc.items() if k != "disclosure"})
    try:
        if render_markdown(doc) != md:
            fails.append("markdown on disk is not the rendering of the shared JSON")
    except Exception as e:  # noqa: BLE001
        fails.append(f"shared markdown cannot be re-rendered: {e}")
    headers = None
    for line in md.splitlines():
        if line.startswith("|") and not line.startswith("|---"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if headers is None:
                headers = cells
                continue
            for h, c in zip(headers, cells):
                if M.COUNT_KEY_RE.match(h) and c not in ("<10", "S", "–", "None") and not (c.isdigit() and int(c) >= M.MIN_CELL):
                    fails.append(f"markdown count column '{h}' has a disallowed value")
        elif not line.startswith("|---"):
            headers = None
    return fails


def publish(M, A, cfg: FUConfig, doc: dict) -> None:
    pj, pm = M._pending(cfg.shared_json), M._pending(cfg.shared_md)
    allp = [cfg.shared_json, cfg.shared_md, pj, pm]
    M._remove(*allp)
    M._atomic_write(pj, json.dumps(doc, indent=1, ensure_ascii=False))
    M._atomic_write(pm, render_markdown(doc))
    fails = release_gate(M, A, pj, pm, "pending")
    if fails:
        M._remove(*allp)
        raise M.ReleaseGateError("; ".join(fails[:20]))
    final = json.loads(pj.read_text(encoding="utf-8"))
    final["release_gate"] = "ok"
    M._atomic_write(cfg.shared_json, json.dumps(final, indent=1, ensure_ascii=False))
    M._atomic_write(cfg.shared_md, render_markdown(final))
    M._remove(pj, pm)
    fails = release_gate(M, A, cfg.shared_json, cfg.shared_md, "ok")
    if fails:
        M._remove(*allp)
        raise M.ReleaseGateError("; ".join(fails[:20]))


# ============================================================================= orchestration
def run_followup(cfg: FUConfig) -> dict:
    M = load_stage2(cfg.stage2_path, cfg.stage2_sha256)
    s2 = make_s2config(M, cfg)
    A = M.load_v17(s2.v17_module_path, s2.v17_module_sha256)
    try:
        log(cfg, f"Stage 2 script OK ({cfg.stage2_sha256[:12]}), v1.7 module OK ({s2.v17_module_sha256[:12]})")
        enc, info = M.load_cohort(s2, A)
        gates, map_hash = load_channel_map_v2(M, cfg)
        log(cfg, f"channel map v2 OK ({map_hash[:12]}, {cfg.channel_map_v2_rows} rows)")
        raw_tables = {}
        for name in ("LABEVENTS", "PATIENTS", "ADMISSIONS", "ICUSTAYS", "SURGERY_INFO"):
            if name not in cfg.raw_sha256:
                raise FollowupError(f"no pinned SHA-256 for raw {name}")
            raw_tables[name] = M.gate_raw(s2, name)
            raw_tables[name]["gate"] = "stage1_sha256" if name == "LABEVENTS" else "d36_pinned_sha256"
        log(cfg, "raw table hash checks: " + ", ".join(f"{k} OK" for k in raw_tables) + " (exact SHA-256)")
        stays, gaps, ii = M.read_icustays(s2, A, enc)
        hfc = enc["h_first_careunit"].astype(str)
        fu = enc["enc_i"].map(ii["first_unit"]).fillna("")
        enc["first_careunit"] = np.where(hfc != "", hfc, np.where(fu != "", fu, "unknown"))
        co = M.Cohort.build(enc, gaps)
        pt = M.read_patients(s2, A, enc)
        band = M.age_band(pt.set_index("enc_i").loc[enc["enc_i"], "age_days"].to_numpy())
        lab, linfo, _ = M.read_events(s2, A, co, "LABEVENTS", set(WANTED), M.LAB_RELEASE_H, None)
        rows = []
        for it, g in lab.groupby("item", sort=True):
            v, st = M.classify_item_rows(gates[it], g)
            rows.append(pd.DataFrame({"enc_i": g["enc_i"].to_numpy(), "tc": g["th"].to_numpy(),
                                      "t": g["th"].to_numpy() + M.LAB_RELEASE_H, "item": it, "v": v, "status": st,
                                      "vnum": g["vnum"].to_numpy()}))
        allr = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
            columns=["enc_i", "tc", "t", "item", "v", "status", "vnum"])
        allr["enc_i"] = allr["enc_i"].astype(np.int64)
        tw = M._train_window(co, allr["enc_i"].to_numpy(), allr["t"].to_numpy(dtype=float))
        meas = allr["status"].to_numpy() == "measurement"
        out = M.Out(A)
        # (i) bilirubin
        train = co.masks["train"]
        fcu = enc["first_careunit"].astype(str).to_numpy()
        nicu = fcu == NICU
        neo = band == "neonate"
        groupings = {"all_train": train,
                     "splits": {"by_first_careunit": {"NICU": train & nicu, "non_NICU": train & ~nicu},
                                "by_age": {"neonate": train & neo, "not_neonate": train & ~neo}}}
        bili = allr[tw & meas & allr["item"].isin([BILI_MAIN, BILI_ALT]).to_numpy()][["enc_i", "tc", "item", "v"]]
        bilirubin_block(out, M, bili, groupings, {"by_first_careunit": {"NICU": 1}, "by_age": {"neonate": 1}})
        # (ii) FiO2: all train-window 5222 rows; same-draw pairs keyed on the pO2 row's availability time (as Stage 2)
        fi_all = allr[allr["item"] == FIO2]
        fi_rows = fi_all[M._train_window(co, fi_all["enc_i"].to_numpy(), fi_all["t"].to_numpy(dtype=float))]
        po = allr[meas & allr["item"].isin([PO2_MAIN, PO2_ALT]).to_numpy()].assign(
            _o=lambda d: (d["item"] != PO2_MAIN).astype(int)).sort_values(["enc_i", "tc", "_o"]).drop_duplicates(["enc_i", "tc"])
        fin = fi_all[fi_all["vnum"].notna()].drop_duplicates(["enc_i", "tc"])
        pairs = po[["enc_i", "tc", "t", "v"]].rename(columns={"v": "po2"}).merge(
            fin[["enc_i", "tc", "vnum"]].rename(columns={"vnum": "fio2"}), on=["enc_i", "tc"])
        pairs = pairs[M._train_window(co, pairs["enc_i"].to_numpy(), pairs["t"].to_numpy(dtype=float))]
        fim = fi_all[fi_all["status"] == "measurement"].drop_duplicates(["enc_i", "tc"])
        s2p = po[["enc_i", "tc", "t", "v"]].merge(fim[["enc_i", "tc", "v"]], on=["enc_i", "tc"], suffixes=("_po2", "_fio2"))
        s2p = s2p[M._train_window(co, s2p["enc_i"].to_numpy(), s2p["t"].to_numpy(dtype=float))]
        pfg = gates["derived:pf_ratio"]
        pfv = s2p["v_po2"].to_numpy(dtype=float) / (s2p["v_fio2"].to_numpy(dtype=float) / 100.0)
        s2p = s2p.assign(valid=(pfv >= pfg["low"]) & (pfv <= pfg["high"]))
        got = fio2_block(out, M, fi_rows, pairs, s2p, pfg, cfg.expected_stage2_pairs)
        log(cfg, f"registry: {len(out.reg.raw)} cells, {len(out.reg.identities)} identities; suppressing (LP) ...")
        status, supp = M.apply_suppression_frozen(A, out.reg, out.frozen)
        log(cfg, f"suppression: {supp['primary_cells']} primary, {supp['complementary_cells']} complementary, "
                 f"{supp['recoverable']} recoverable")
        hashes = dict(info["hashes"])
        hashes.update({"stage2_script": cfg.stage2_sha256, "v17_module": s2.v17_module_sha256, "feature_channel_map_v2": map_hash})
        hashes.update({f"raw_{k}": v["sha256"] for k, v in raw_tables.items()})
        notes = {
            "groups": "NICU = first care unit 'NICU' (harmonization first_careunit, else the earliest ICUSTAYS unit, as Stage 2); "
                      "neonate = age < 28 days at t0; age_invalid encounters are in not_neonate",
            "bilirubin_values": "measurement rows inside the channel-map v2 gate, in the recorded unit (umol/L), train window by availability time",
            "same_draw": "same encounter and CHARTTIME; one value per item per draw",
            "fio2_rows": "all train-window 5222 rows (any value); bands on the numeric value as recorded",
            "fio2_pairs": "pO2 measurement rows (5239 preferred over 5244) joined to the numeric FiO2 row of the same draw; train window by the pO2 availability time",
            "reproduction": "stage2_definition uses FiO2 measurement rows (gate 21-100) as Stage 2 v1.0.1; expected valid / outside from the v1.0.1 output",
            "lab_rows_read": int(linfo["rows"]),
        }
        doc: dict[str, Any] = {
            "spec": SCRIPT_VERSION, "release_gate": "pending", "input_hashes": hashes,
            "stage2_module": {"file": STAGE2_FILE, "sha256": cfg.stage2_sha256},
            "v17_module": {"file": M.V17_MODULE_FILE, "sha256": s2.v17_module_sha256},
            "maps": {"channel": {"file": CHANNEL_MAP_V2_FILE, "sha256": map_hash}},
            "raw_tables": raw_tables, "cohort_checks": info["checks"],
            "config": {"lab_release_hours": M.LAB_RELEASE_H, "analysis_window": "[t0, min(OUT, t0 + 336 h))",
                       "pf_gate": [pfg["low"], pfg["high"]], "min_cell": M.MIN_CELL,
                       "fio2_bands": ["<= 1", "(1, 21)", "[21, 100]", "> 100"]},
            "notes": notes,
        }
        M.fill_document(out, status, doc)
        doc["disclosure"] = {"count_fields": sorted(out.reg.published()),
                             "identities": [[t, list(p)] for t, p in out.reg.identities],
                             "hidden_cells": sorted(out.reg.hidden),
                             "dependents": [[p, list(nd)] for p, nd, _ in out.kss],
                             "quantiles": [[p, n] for p, n, _ in out.quants]}
        doc["suppression"] = dict(supp)
        log(cfg, f"Stage 2 pair definition reproduces v1.0.1: {got == {k: int(v) for k, v in cfg.expected_stage2_pairs.items()}}")
    except (A.HashMismatchError, M.HashMismatchError) as e:
        raise FollowupError(str(e)) from e
    except (A.ReleaseGateError, M.ReleaseGateError):
        raise
    except (A.AuditError, M.AuditError) as e:
        raise FollowupError(str(e)) from e
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    publish(M, A, cfg, doc)
    final = json.loads(cfg.shared_json.read_text(encoding="utf-8"))
    if cfg.verbose:
        print(cfg.shared_md.read_text(encoding="utf-8"))
        print(f"release gate ok: {cfg.shared_json.name}, {cfg.shared_md.name} written to {cfg.out_dir}")
    return final


def build_config(argv: list[str] | None = None) -> FUConfig:
    here = Path(__file__).resolve()
    default_root = Path(os.environ.get("FOUNDATIONS_ROOT", str(here.parents[1])))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", default=str(default_root))
    ap.add_argument("--raw-dir", default=os.environ.get("PIC_RAW_DIR"))
    ap.add_argument("--split-file", default=os.environ.get("SPLIT_FILE_PATH"))
    ap.add_argument("--canonical", default=None)
    ap.add_argument("--harmonization", default=None)
    ap.add_argument("--complexity-map", default=None)
    ap.add_argument("--maps-dir", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--harm-t0-col", default=None)
    ap.add_argument("--harm-out-col", default=None)
    ap.add_argument("--split-col", default="split")
    ap.add_argument("--stage2-script", default=str(here.parent / STAGE2_FILE))
    ap.add_argument("--v17-module", default=str(here.parent / "audit_surgery_complexity_duration.py"))
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    a = ap.parse_args(argv)
    if not a.raw_dir or not a.split_file:
        raise AuditError("--raw-dir and --split-file are required (or set PIC_RAW_DIR and SPLIT_FILE_PATH); no local default path is shipped")
    root = Path(a.project_root)
    interim = root / "data" / "interim"
    return FUConfig(
        raw_dir=Path(a.raw_dir),
        canonical_path=Path(a.canonical) if a.canonical else interim / "canonical_icu_encounters.csv",
        harmonization_path=Path(a.harmonization) if a.harmonization else interim / "encounter_harmonization_fields.csv",
        split_path=Path(a.split_file),
        complexity_map_path=Path(a.complexity_map) if a.complexity_map else root / "config" / "surgery_complexity_map_v1.csv",
        maps_dir=Path(a.maps_dir) if a.maps_dir else interim / "maps",
        out_dir=Path(a.out_dir) if a.out_dir else interim,
        stage2_path=Path(a.stage2_script), v17_module_path=Path(a.v17_module),
        harm_t0_col=a.harm_t0_col, harm_out_col=a.harm_out_col, split_col=a.split_col, chunksize=a.chunksize)


def main(argv: list[str] | None = None) -> int:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    cfg = build_config(argv)
    try:
        run_followup(cfg)
    except FollowupError as e:
        print(f"AUDIT STOPPED: {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001
        if type(e).__name__ == "ReleaseGateError":
            print(f"RELEASE GATE FAILED - no shared file left in place: {e}", file=sys.stderr)
            return 3
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
