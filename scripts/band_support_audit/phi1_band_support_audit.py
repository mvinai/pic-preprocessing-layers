#!/usr/bin/env python3
"""
phi1 duration-band support audit (read-only, aggregate-only). Implements `band_support_audit_spec.md` v1.0.

phi1 = surgery type x type-specific duration band (decision 16 and its follow-up, Nov 11-12, 2026):
    D = SURGERY_END_TIME - SURGERY_BEGIN_TIME (hours), falling back to ANES_END_TIME - ANES_START_TIME only when a
        surgery time is missing/unparseable (the duration audit v1.7 rule, reused unchanged). D missing or <= 0 is
        excluded and counted; D > 24 h is counted, not excluded.
    reference median (TRAIN split only, frozen): own type (standardized_name) if >= MIN_N train surgeries with a
        valid D, else system x complexity group, else system, else all surgeries (first level with >= MIN_N).
        MIN_N = 10 (primary), 20 (sensitivity).
    r = D / reference median; short r < 0.75, expected 0.75 <= r <= 1.33, long r > 1.33.
    Sensitivity: within-type train tertiles, for types with their own median only.

Reuse: the inputs, full-hash gates, encounter linkage, map handling, duration rule, disclosure registry,
complementary suppression and recoverability solver come from the duration audit v1.7 module
(`audit_surgery_complexity_duration.py`), which is SHA-256-pinned below and verified BEFORE it is imported.
Population = surgeries linked to a cohort encounter (v1.7 linkage; split from that encounter) with a mapped name.

Usage (Windows, from the project root; both scripts in scripts/):
    python scripts/phi1_band_support_audit.py
    python scripts/phi1_band_support_audit.py --raw-dir "C:/.../PICDataset" --split-file "C:/.../subject_splits.csv"

Outputs (data/interim/):
    phi1_band_support_audit_MAP_v1.json / .md                  shared (suppressed); written ONLY if the gate passes
    phi1_band_support_frozen_medians_MAP_v1.csv                shared frozen-medians table (gated with the JSON)
    phi1_band_support_audit_MAP_v1_raw_LOCAL.json              local only, never share
    phi1_band_support_frozen_medians_MAP_v1_raw_LOCAL.csv      local only, never share (full-precision medians)
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
# CONFIG BLOCK (approved Nov 11-12, 2026). Thresholds are read from here, never hard-coded in logic.
# =============================================================================
SPEC_VERSION = "phi1_band_support_audit v1.0"
MAP_LABEL = "MAP_v1"
V17_MODULE_FILE = "audit_surgery_complexity_duration.py"
V17_SPEC_VERSION = "surgery_complexity_duration_audit v1.7"
V17_MODULE_SHA256 = "1de2b17b2d894caff6a8362a9c72a458967691a8d2c7a533db863e87ad8dad94"  # public copy (cohort pins omitted); original private pin was ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0
PINNED_SHA256 = {
    "canonical": "e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418",
    "harmonization": "ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d",
    "complexity_map": "5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5",
    "split": "60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31",
}
PRIMARY_MIN_N = 10
SENSITIVITY_MIN_N = 20
THRESHOLDS = (PRIMARY_MIN_N, SENSITIVITY_MIN_N)
SHORT_BELOW = 0.75        # short: r < 0.75
LONG_ABOVE = 1.33         # long: r > 1.33; expected: 0.75 <= r <= 1.33 (inclusive)
RATIO_DECIMALS = 9        # r is rounded to 9 decimals before banding (removes float noise at the cut points)
MIN_CELL = 10
D_GT_HOURS = 24.0
SHARED_MEDIAN_DECIMALS = 2

LEVELS = ("own_type", "system_x_group", "system", "all")
BANDS = ("short", "expected", "long")
SPLITS = ("train", "validation", "test")
UNCLASSIFIED_GROUP = "unclassified"
INVALID_REASONS = ("begin_after_end", "zero_duration", "missing_window", "fallback_nonpositive")
LINK_CLASSES = ("linked", "post_icu_same_hadm", "other_hadm", "not_in_cohort", "invalid_time",
                "not_linkable_missing_encounter_time")
COUNT_KEYS = frozenset({"n", "short", "expected", "long", "types", "surgeries", "surgery_info_rows", "unmapped",
                        "mapped", "invalid_duration", "valid", "duration_fallback_anes", "d_gt_24h",
                        "types_in_more_than_one_level"} | set(LINK_CLASSES) | set(INVALID_REASONS))
NON_COUNT_TOP = ("disclosure", "suppression", "input_hashes", "config", "v17_module")
TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}")


def tkey(k: int) -> str:
    return f"min_n_{k}"


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
    """Verify the duration audit v1.7 file byte-for-byte, then import it under a private module name."""
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
    sys.modules[name] = mod  # dataclasses need the module registered before exec
    try:
        spec.loader.exec_module(mod)
    except Exception:
        sys.modules.pop(name, None)
        raise
    if getattr(mod, "SPEC_VERSION", None) != V17_SPEC_VERSION or getattr(mod, "PINNED_SHA256", None) != PINNED_SHA256:
        raise AuditError("v1.7 module constants differ from this script's pinned values")
    return mod


# =============================================================================
# Config
# =============================================================================
@dataclass
class BandConfig:
    canonical_path: Path
    harmonization_path: Path
    map_path: Path
    split_path: Path
    surgery_info_path: Path
    out_dir: Path
    v17_module_path: Path
    admissions_path: Path | None = None
    pinned_sha256: dict = field(default_factory=lambda: dict(PINNED_SHA256))
    v17_module_sha256: str = V17_MODULE_SHA256
    harm_t0_col: str | None = None
    harm_out_col: str | None = None
    split_subject_col: str = "SUBJECT_ID"
    split_col: str = "split"
    expected_surgery_rows: int | None = None  # public copy: pin omitted
    expected_encounters: int | None = None  # public copy: pin omitted
    expected_subjects: int | None = None
    expected_split_encounters: dict = field(default_factory=dict)
    expected_map_rows: int = 1495
    expected_missing: dict | None = None
    verbose: bool = True

    def _p(self, stem: str, ext: str) -> Path:
        return self.out_dir / f"{stem}_{MAP_LABEL}{ext}"

    @property
    def shared_json(self) -> Path:
        return self._p("phi1_band_support_audit", ".json")

    @property
    def shared_md(self) -> Path:
        return self._p("phi1_band_support_audit", ".md")

    @property
    def shared_csv(self) -> Path:
        return self._p("phi1_band_support_frozen_medians", ".csv")

    @property
    def local_json(self) -> Path:
        return self._p("phi1_band_support_audit", "_raw_LOCAL.json")

    @property
    def local_csv(self) -> Path:
        return self._p("phi1_band_support_frozen_medians", "_raw_LOCAL.csv")


def v17_config(cfg: BandConfig, A) -> Any:
    kw = {}
    if cfg.expected_missing is not None:
        kw["expected_missing"] = dict(cfg.expected_missing)
    return A.AuditConfig(
        canonical_path=Path(cfg.canonical_path), harmonization_path=Path(cfg.harmonization_path), map_path=Path(cfg.map_path),
        split_path=Path(cfg.split_path), surgery_info_path=Path(cfg.surgery_info_path), out_dir=Path(cfg.out_dir),
        admissions_path=cfg.admissions_path, pinned_sha256=dict(cfg.pinned_sha256), harm_t0_col=cfg.harm_t0_col,
        harm_out_col=cfg.harm_out_col, split_subject_col=cfg.split_subject_col, split_col=cfg.split_col,
        expected_surgery_rows=cfg.expected_surgery_rows, expected_encounters=cfg.expected_encounters,
        expected_subjects=cfg.expected_subjects, expected_split_encounters=dict(cfg.expected_split_encounters),
        expected_map_rows=cfg.expected_map_rows, bootstrap_reps=0, verbose=False, **kw)


# =============================================================================
# Pure functions (unit tested)
# =============================================================================
def band_of_ratio(r: float) -> str:
    if r < SHORT_BELOW:
        return "short"
    if r <= LONG_ABOVE:
        return "expected"
    return "long"


def ratio(d: float, m: float) -> float:
    return round(float(d) / float(m), RATIO_DECIMALS)


def tertile_band(d: float, q1: float, q2: float) -> str:
    if d < q1:
        return "short"
    if d <= q2:
        return "expected"
    return "long"


def key_label(key: tuple) -> str:
    lvl = key[0]
    if lvl == "own_type":
        return f"type:{key[1]}"
    if lvl == "system_x_group":
        return f"system_x_group:{key[1]} x {key[2]}"
    if lvl == "system":
        return f"system:{key[1]}"
    return "all"


@dataclass
class Frozen:
    min_n: int
    type_stats: dict      # type -> (n, median) for every type seen in train
    sg_stats: dict        # (system, group) -> (n, median)
    sys_stats: dict       # system -> (n, median)
    all_stats: tuple      # (n, median)

    @property
    def own_types(self) -> set:
        return {t for t, (n, _) in self.type_stats.items() if n >= self.min_n}

    def reference(self, t: str, bs: str, g: str) -> tuple[tuple, float]:
        """Fallback order: own type -> system x group -> system -> all (first level with n >= min_n)."""
        n, m = self.type_stats.get(t, (0, float("nan")))
        if n >= self.min_n:
            return ("own_type", t), m
        n, m = self.sg_stats.get((bs, g), (0, float("nan")))
        if n >= self.min_n:
            return ("system_x_group", bs, g), m
        n, m = self.sys_stats.get(bs, (0, float("nan")))
        if n >= self.min_n:
            return ("system", bs), m
        return ("all",), self.all_stats[1]

    def stats_of(self, key: tuple) -> tuple[int, float]:
        lvl = key[0]
        if lvl == "own_type":
            return self.type_stats[key[1]]
        if lvl == "system_x_group":
            return self.sg_stats[(key[1], key[2])]
        if lvl == "system":
            return self.sys_stats[key[1]]
        return self.all_stats


def freeze_medians(train: pd.DataFrame, min_n: int) -> Frozen:
    """train: banded TRAIN rows only (columns type, system, group, D). Medians never see validation/test rows."""
    if set(train.get("split", pd.Series(["train"] * len(train)))) - {"train"}:
        raise AuditError("freeze_medians received non-train rows")
    if len(train) < min_n:
        raise AuditError(f"fewer than {min_n} train surgeries with a valid duration; cannot freeze the 'all' median")

    def st(g):
        return {k: (int(len(v)), float(np.median(v.to_numpy(dtype=float)))) for k, v in g["D"]}
    return Frozen(min_n=min_n, type_stats=st(train.groupby("type")), sg_stats=st(train.groupby(["system", "group"])),
                  sys_stats=st(train.groupby("system")),
                  all_stats=(int(len(train)), float(np.median(train["D"].to_numpy(dtype=float)))))


def assign_bands(df: pd.DataFrame, fr: Frozen) -> pd.DataFrame:
    """df: banded rows of any split (type, system, group, D). Returns level, key, ref median, r, band."""
    keys, refs = [], []
    cache: dict = {}
    for t, bs, g in zip(df["type"], df["system"], df["group"]):
        k = (t, bs, g)
        if k not in cache:
            cache[k] = fr.reference(t, bs, g)
        keys.append(cache[k][0]); refs.append(cache[k][1])
    r = [ratio(d, m) for d, m in zip(df["D"], refs)]
    return pd.DataFrame({"level": [k[0] for k in keys], "key": keys, "ref": refs, "r": r,
                         "band": [band_of_ratio(x) for x in r]}, index=df.index)


def freeze_tertiles(train: pd.DataFrame, own_types: Iterable[str]) -> dict:
    own = set(own_types)
    out = {}
    for t, g in train[train["type"].isin(own)].groupby("type"):
        q1, q2 = np.quantile(g["D"].to_numpy(dtype=float), [1 / 3, 2 / 3])
        out[t] = (float(q1), float(q2))
    return out


def prepare(frame: pd.DataFrame, cmap: pd.DataFrame, unclassified: str) -> pd.DataFrame:
    """One row per SURGERY_INFO row (ROW_ID). Adds type (standardized_name), system, group, population flags."""
    std = cmap.set_index("key")["standardized_name"].map(lambda v: re.sub(r"\s+", " ", str(v).strip()))
    df = frame[["ROW_ID", "SUBJECT_ID", "link_class", "split", "mapped", "D", "valid", "fallback", "gt24",
                "invalid_reason", "body_system", "cgroup"]].copy()
    df["type"] = frame["name"].map(std).fillna("")
    if (df["mapped"] & (df["type"] == "")).any():
        raise AuditError("a mapped surgery has an empty standardized_name in the map")
    df["system"] = df["body_system"]
    df["group"] = np.where(df["body_system"] == unclassified, UNCLASSIFIED_GROUP, df["cgroup"])
    if (df["mapped"] & (df["group"] == "")).any():
        raise AuditError("a mapped surgery has no complexity group")
    df["in_pop"] = (df["link_class"] == "linked") & df["mapped"]
    df["banded"] = df["in_pop"] & df["valid"]
    if (df["banded"] & ~(df["D"] > 0)).any():
        raise AuditError("a banded surgery has a non-positive duration")
    names = set(df.loc[df["mapped"], "type"]) | set(df.loc[df["mapped"], "system"]) | set(df.loc[df["mapped"], "group"])
    bad = [x for x in names if x in COUNT_KEYS or x.startswith("pct_")]
    if bad:
        raise AuditError("a type/system/group name collides with a reserved output key")
    return df


# =============================================================================
# Analysis -> disclosure registry
# =============================================================================
class BandAnalysis:
    def __init__(self, A, df: pd.DataFrame, rc: dict) -> None:
        self.A, self.df, self.rc = A, df, rc
        self.reg = A.Registry()
        self.stats: dict[str, Any] = {}
        self.local: dict[str, Any] = {}
        self.frozen: dict[int, Frozen] = {}
        self.assigned: dict[int, pd.DataFrame] = {}
        self.tertiles: dict = {}
        self.own_cells: dict[tuple, dict] = {}   # (split, type) -> {field: cell id}

    def c(self, value: int, *parts: Any) -> str:
        return self.reg.count(self.A.P(*parts), int(value))

    def superset(self, total: str, members: list[str], tag: str) -> None:
        """total = sum(members) + hidden rest (protects differencing against the own-type table)."""
        rest = self.reg.raw[total] - sum(self.reg.raw[m] for m in members)
        if rest < 0:
            raise AuditError(f"reconciliation: members exceed total at {total}")
        h = self.reg.hide(f"{total}|rest|{tag}", rest)
        self.reg.identity(total, list(members) + [h])

    def run(self) -> None:
        df, reg = self.df, self.reg
        # ---- 1. reconciliation over all SURGERY_INFO rows (keyed by ROW_ID)
        tot = self.c(len(df), "reconciliation", "surgery_info_rows")
        parts = [self.c(self.rc[k], "reconciliation", k) for k in LINK_CLASSES]
        reg.identity(tot, parts)
        linked_split = []
        valid_cell = {}
        for sp in SPLITS:
            d = df[(df["link_class"] == "linked") & (df["split"] == sp)]
            b = ("population", sp)
            ln = self.c(len(d), *b, "linked")
            linked_split.append(ln)
            mp = self.c(int(d["mapped"].sum()), *b, "mapped")
            un = self.c(int((~d["mapped"]).sum()), *b, "unmapped")
            reg.identity(ln, [mp, un])
            m = d[d["mapped"]]
            va = self.c(int(m["valid"].sum()), *b, "valid")
            iv = self.c(int((~m["valid"]).sum()), *b, "invalid_duration")
            reg.identity(mp, [va, iv])
            reg.identity(iv, [self.c(int((m["invalid_reason"] == r).sum()), *b, "invalid_by_reason", r) for r in INVALID_REASONS])
            reg.subset(self.c(int((m["valid"] & m["fallback"]).sum()), *b, "duration_fallback_anes"), va, "fallback")
            reg.subset(self.c(int(m["gt24"].sum()), *b, "d_gt_24h"), va, "gt24")
            self.c(int(m.loc[m["valid"], "type"].nunique()), *b, "types")
            valid_cell[sp] = va
        reg.identity(parts[0], linked_split)
        self.valid_cell = valid_cell

        bd = df[df["banded"]]
        train = bd[bd["split"] == "train"][["type", "system", "group", "D", "split"]]
        for k in THRESHOLDS:
            self.frozen[k] = freeze_medians(train, k)
            self.assigned[k] = assign_bands(bd, self.frozen[k])
        own10 = sorted(self.frozen[PRIMARY_MIN_N].own_types)
        self.tertiles = freeze_tertiles(train, own10)
        self.local["checks_train_rows_used_for_medians"] = int(len(train))

        # ---- 2. own-type table (primary threshold): per split n and bands
        a10 = self.assigned[PRIMARY_MIN_N]
        for sp in SPLITS:
            sub = bd[bd["split"] == sp]
            asg = a10.loc[sub.index]
            for t in own10:
                rows = asg[sub["type"] == t]
                cells = {b: self.c(int((rows["band"] == b).sum()), "own_types", sp, t, b) for b in BANDS}
                n = self.c(len(rows), "own_types", sp, t, "n")
                reg.identity(n, [cells[b] for b in BANDS])
                cells["n"] = n
                self.own_cells[(sp, t)] = cells
                self.stats[self.A.P("own_types", sp, t, "any_band_below_10")] = "Y" if any(int((rows["band"] == b).sum()) < MIN_CELL for b in BANDS) else "N"

        # ---- 3. median source levels and bands, both thresholds
        for k in THRESHOLDS:
            asg_all = self.assigned[k]
            for sp in SPLITS:
                sub = bd[bd["split"] == sp]
                asg = asg_all.loc[sub.index]
                lv_cells = {}
                for lv in LEVELS:
                    msk = asg["level"] == lv
                    lv_cells[lv] = self.c(int(msk.sum()), "median_source", tkey(k), sp, lv, "surgeries")
                    self.c(int(sub.loc[msk, "type"].nunique()), "median_source", tkey(k), sp, lv, "types")
                reg.identity(valid_cell[sp], [lv_cells[lv] for lv in LEVELS])
                lv_per_type = pd.DataFrame({"t": sub["type"], "l": asg["level"]}).drop_duplicates().groupby("t").size()
                self.c(int((lv_per_type > 1).sum()), "median_source", tkey(k), sp, "types_in_more_than_one_level")
                if k == PRIMARY_MIN_N:
                    reg.identity(lv_cells["own_type"], [self.own_cells[(sp, t)]["n"] for t in own10])
                else:
                    reg.subset(lv_cells["own_type"], self.c_ptr("median_source", tkey(PRIMARY_MIN_N), sp, "own_type", "surgeries"), f"own{k}")
                # overall bands
                ob = {b: self.c(int((asg["band"] == b).sum()), "bands", tkey(k), sp, "overall", b) for b in BANDS}
                on = self.c(len(sub), "bands", tkey(k), sp, "overall", "n")
                reg.identity(on, [ob[b] for b in BANDS])
                reg.equal(on, valid_cell[sp])
                if k != PRIMARY_MIN_N:
                    continue
                for b in BANDS:
                    self.superset(ob[b], [self.own_cells[(sp, t)][b] for t in own10], f"overall_{b}")
                for dim, col in (("by_system", "system"), ("by_group", "group")):
                    tot_n, tot_b = [], {b: [] for b in BANDS}
                    for val in sorted(sub[col].unique()):
                        msk = sub[col] == val
                        bc = {b: self.c(int((asg.loc[msk, "band"] == b).sum()), "bands", tkey(k), sp, dim, val, b) for b in BANDS}
                        nc = self.c(int(msk.sum()), "bands", tkey(k), sp, dim, val, "n")
                        reg.identity(nc, [bc[b] for b in BANDS])
                        tot_n.append(nc)
                        for b in BANDS:
                            tot_b[b].append(bc[b])
                        inside = [t for t in own10 if self._inside(sp, t, col, val)]
                        for f in ("n",) + BANDS:
                            members = [self.own_cells[(sp, t)][f] for t in inside]
                            if members:
                                self.superset(bc[f] if f != "n" else nc, members, f"{dim}_{f}")
                    reg.identity(on, tot_n)
                    for b in BANDS:
                        reg.identity(ob[b], tot_b[b])

        # ---- 4. tertile sensitivity (own types only; train tertiles)
        for sp in SPLITS:
            sub = bd[(bd["split"] == sp) & bd["type"].isin(self.tertiles)]
            tb = [tertile_band(d, *self.tertiles[t]) for d, t in zip(sub["D"], sub["type"])]
            cells = [self.c(sum(1 for x in tb if x == b), "tertile_sensitivity", sp, b) for b in BANDS]
            n = self.c(len(sub), "tertile_sensitivity", sp, "n")
            reg.identity(n, cells)
            reg.equal(n, self.c_ptr("median_source", tkey(PRIMARY_MIN_N), sp, "own_type", "surgeries"))

        # ---- 5. types unseen in train (validation / test)
        seen = set(train["type"])
        for sp in SPLITS[1:]:
            sub = bd[bd["split"] == sp]
            uns = sub[~sub["type"].isin(seen)]
            sc = self.c(len(uns), "unseen_in_train", sp, "surgeries")
            reg.subset(sc, valid_cell[sp], "unseen")
            self.c(int(uns["type"].nunique()), "unseen_in_train", sp, "types")
            reg.equal(self.c(len(sub), "unseen_in_train", sp, "n"), valid_cell[sp])

        # ---- 6. frozen medians (keys used by at least one surgery, plus every own type)
        for k in THRESHOLDS:
            fr = self.frozen[k]
            used = {kk for kk in self.assigned[k]["key"]} | {("own_type", t) for t in fr.own_types}
            local_rows = []
            for key in sorted(used, key=lambda x: (LEVELS.index(x[0]), x[1:])):
                n, med = fr.stats_of(key)
                if n < k:
                    raise AuditError(f"frozen key below threshold: {key_label(key)}")
                lab = key_label(key)
                cid = self.c(n, "frozen_medians", tkey(k), lab, "n")
                self.stats[self.A.P("frozen_medians", tkey(k), lab, "level")] = key[0]
                self.stats[self.A.P("frozen_medians", tkey(k), lab, "median_hours")] = round(med, SHARED_MEDIAN_DECIMALS)
                local_rows.append({"threshold": k, "key": lab, "level": key[0], "train_n": n, "median_hours": med})
                if key[0] == "own_type":
                    reg.equal(cid, self.own_cells[("train", key[1])]["n"])
                elif key[0] == "all":
                    reg.equal(cid, valid_cell["train"])
                elif key[0] == "system":
                    reg.equal(cid, self.c_ptr("bands", tkey(PRIMARY_MIN_N), "train", "by_system", key[1], "n"))
                else:
                    inside = [t for t in own10 if self._inside("train", t, "system", key[1]) and self._inside("train", t, "group", key[2])]
                    members = [self.own_cells[("train", t)]["n"] for t in inside]
                    if members:
                        self.superset(cid, members, "sg_pool")
            self.local[f"frozen_medians_{tkey(k)}"] = local_rows
        self.local["tertile_cuts_hours"] = {t: list(v) for t, v in self.tertiles.items()}

    def c_ptr(self, *parts: Any) -> str:
        p = self.A.P(*parts)
        if p not in self.reg.raw:
            raise AuditError(f"internal: missing cell {p}")
        return p

    def _inside(self, sp: str, t: str, col: str, val: str) -> bool:
        """True if every banded surgery of type t in split sp has df[col] == val (and there is at least one)."""
        key = (sp, t, col)
        if not hasattr(self, "_inside_cache"):
            self._inside_cache = {}
        if key not in self._inside_cache:
            d = self.df[self.df["banded"] & (self.df["split"] == sp) & (self.df["type"] == t)]
            vals = set(d[col])
            self._inside_cache[key] = next(iter(vals)) if len(vals) == 1 else None
        return self._inside_cache[key] == val


# =============================================================================
# Document, percentages, rendering
# =============================================================================
def add_percentages(doc: dict) -> None:
    """Percentages are published only when the count AND its denominator are published integers."""
    def walk(o: Any) -> None:
        if isinstance(o, dict):
            n = o.get("n")
            for b in BANDS:
                if b in o:
                    v = o[b]
                    o[f"pct_{b}"] = round(100.0 * v / n, 1) if isinstance(v, int) and isinstance(n, int) and n > 0 else "suppressed"
            if "surgeries" in o and "n" in o and "types" in o:  # unseen_in_train block
                v = o["surgeries"]
                o["pct_surgeries"] = round(100.0 * v / n, 1) if isinstance(v, int) and isinstance(n, int) and n > 0 else "suppressed"
            for v in list(o.values()):
                walk(v)
    for key in ("bands", "own_types", "tertile_sensitivity", "unseen_in_train"):
        walk(doc.get(key, {}))


def build_document(A, an: BandAnalysis, status: dict, supp: dict, hashes: dict) -> dict:
    doc: dict[str, Any] = {
        "spec": SPEC_VERSION, "map_label": MAP_LABEL, "input_hashes": dict(hashes),
        "v17_module": {"file": V17_MODULE_FILE, "spec": V17_SPEC_VERSION, "sha256": V17_MODULE_SHA256},
        "config": {"primary_min_n": PRIMARY_MIN_N, "sensitivity_min_n": SENSITIVITY_MIN_N, "short_below": SHORT_BELOW,
                   "long_above": LONG_ABOVE, "expected_inclusive": True, "ratio_decimals": RATIO_DECIMALS,
                   "fallback_order": list(LEVELS), "duration": "SURGERY_END - SURGERY_BEGIN; fallback ANES_END - ANES_START",
                   "medians": "train split only, frozen", "d_gt_24h": "counted, not excluded",
                   "population": "linked to a cohort encounter (v1.7 linkage) and mapped; valid D > 0"},
        "zero_policy": "counts 0-9 are published as '<10'; 'S' = complementary suppression",
    }
    for ptr, v in an.stats.items():
        A.set_ptr(doc, ptr, v)
    for cid in an.reg.published():
        st = status[cid]
        A.set_ptr(doc, cid, an.reg.raw[cid] if st == "pub" else ("<10" if st == "P" else "S"))
    add_percentages(doc)
    doc["suppression"] = dict(supp)
    doc["disclosure"] = {"count_fields": sorted(an.reg.published()),
                         "identities": [[t, list(p)] for t, p in an.reg.identities],
                         "hidden_cells": sorted(an.reg.hidden), "markers": {"primary": "<10", "complementary": "S"}}
    doc["release_gate"] = "pending"
    return doc


def detect_count_fields(obj: Any, prefix: str = "") -> list[str]:
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if prefix == "" and k in NON_COUNT_TOP:
                continue
            p = prefix + "/" + str(k).replace("~", "~0").replace("/", "~1")
            if isinstance(v, (dict, list)):
                out += detect_count_fields(v, p)
            elif k in COUNT_KEYS:
                out.append(p)
    return out


def _f(v: Any) -> str:
    if v is None:
        return "–"
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".") if v != int(v) else f"{v:.1f}"
    return str(v).replace("|", "/")


def _band_row(label: str, o: dict) -> str:
    return "| " + " | ".join([_f(label), _f(o.get("n"))] + [_f(o.get(b)) for b in BANDS] + [_f(o.get(f"pct_{b}")) for b in BANDS]) + " |"


BAND_HDR = "| group | n | short | expected | long | pct_short | pct_expected | pct_long |\n|---|---|---|---|---|---|---|---|"


def render_markdown(doc: dict) -> str:
    L = [f"# phi1 duration-band support audit: {doc['map_label']}", "",
         f"- Spec `{doc['spec']}`; release gate: **{doc['release_gate']}**; duration audit module `{doc['v17_module']['spec']}` "
         f"(SHA-256 `{doc['v17_module']['sha256']}`).",
         "- Input hashes (exact 64-character match): " + ", ".join(f"{k} `{v}`" for k, v in doc["input_hashes"].items()) + ".",
         "- Aggregate only. Counts of 0–9 are `<10`; `S` = complementary suppression. Percentages only where the count and its denominator are published.",
         f"- Bands: r = D / frozen train median; short r < {SHORT_BELOW}, expected {SHORT_BELOW} ≤ r ≤ {LONG_ABOVE}, long r > {LONG_ABOVE}. "
         f"Fallback: own type → system × complexity group → system → all; min n {PRIMARY_MIN_N} (sensitivity {SENSITIVITY_MIN_N}).",
         "- Suppression: " + ", ".join(f"{k} {v}" for k, v in doc["suppression"].items()) + ".", "",
         "## 1. Reconciliation (all SURGERY_INFO rows, ROW_ID)", "", "| class | count |", "|---|---|"]
    for k, v in doc["reconciliation"].items():
        L.append(f"| {k} | {_f(v)} |")
    L += ["", "## 2. Population by split", "",
          "| split | linked | mapped | unmapped | valid | invalid_duration | duration_fallback_anes | d_gt_24h | types |",
          "|---|---|---|---|---|---|---|---|---|"]
    for sp, p in doc["population"].items():
        L.append("| " + " | ".join([sp] + [_f(p[k]) for k in ("linked", "mapped", "unmapped", "valid", "invalid_duration",
                                                                  "duration_fallback_anes", "d_gt_24h", "types")]) + " |")
    L += ["", "| split | " + " | ".join(INVALID_REASONS) + " |", "|---|" + "---|" * len(INVALID_REASONS)]
    for sp, p in doc["population"].items():
        L.append("| " + " | ".join([sp] + [_f(p["invalid_by_reason"][r]) for r in INVALID_REASONS]) + " |")
    L += ["", "## 3. Median source level", ""]
    for tk, by in doc["median_source"].items():
        L += [f"### {tk}", "", "| split | level | types | surgeries |", "|---|---|---|---|"]
        for sp, lv in by.items():
            for lvl in LEVELS:
                L.append(f"| {sp} | {lvl} | {_f(lv[lvl]['types'])} | {_f(lv[lvl]['surgeries'])} |")
            L.append(f"| {sp} | types_in_more_than_one_level | {_f(lv['types_in_more_than_one_level'])} | – |")
        L.append("")
    L += ["## 4. Band distribution", ""]
    for tk, by in doc["bands"].items():
        for sp, o in by.items():
            L += [f"### {tk}, {sp}", "", BAND_HDR, _band_row("overall", o["overall"])]
            for dim in ("by_system", "by_group"):
                for val, oo in o.get(dim, {}).items():
                    L.append(_band_row(f"{dim}: {val}", oo))
            L.append("")
    L += [f"## 5. Own-median types (min n {PRIMARY_MIN_N}): n per band by split", ""]
    for sp, types in doc["own_types"].items():
        L += [f"### {sp}", "", "| type | n | short | expected | long | pct_short | pct_expected | pct_long | any_band_below_10 |",
              "|---|---|---|---|---|---|---|---|---|"]
        for t, o in types.items():
            L.append(_band_row(t, o)[:-1] + f"| {_f(o['any_band_below_10'])} |")
        L.append("")
    L += ["## 6. Types unseen in train", "", "| split | types | surgeries | n | pct_surgeries |", "|---|---|---|---|---|"]
    for sp, o in doc["unseen_in_train"].items():
        L.append(f"| {sp} | {_f(o['types'])} | {_f(o['surgeries'])} | {_f(o['n'])} | {_f(o['pct_surgeries'])} |")
    L += ["", "## 7. Tertile sensitivity (own-median types; train tertiles)", "", BAND_HDR]
    for sp, o in doc["tertile_sensitivity"].items():
        L.append(_band_row(sp, o))
    L += ["", "## 8. Frozen medians", "", "See the shared CSV `phi1_band_support_frozen_medians_MAP_v1.csv` (same values).", ""]
    return "\n".join(L) + "\n"


CSV_COLUMNS = ["threshold", "key", "level", "train_n", "median_hours"]


def render_csv(doc: dict) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(CSV_COLUMNS)
    for tk, keys in doc["frozen_medians"].items():
        for lab, o in keys.items():
            w.writerow([tk.removeprefix("min_n_"), lab, o["level"], o["n"], o["median_hours"]])
    return buf.getvalue()


# =============================================================================
# Release gate (disk reload) and publishing
# =============================================================================
def release_gate(A, json_path: Path, md_path: Path, csv_path: Path, expect_status: str, sentinels: Iterable[str] = ()) -> list[str]:
    fails: list[str] = []
    try:
        doc = json.loads(Path(json_path).read_text(encoding="utf-8"))
        md = Path(md_path).read_text(encoding="utf-8")
        csv_text = Path(csv_path).read_text(encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        return [f"cannot reload shared files: {e}"]
    if doc.get("spec") != SPEC_VERSION or doc.get("map_label") != MAP_LABEL or doc.get("release_gate") != expect_status:
        fails.append("header fields (spec / map_label / release_gate) wrong")
    for p in (json_path, md_path, csv_path):
        stem = Path(p).name.rsplit(".", 1)[0].removesuffix(".pending")
        if not stem.endswith(f"_{MAP_LABEL}") or "_raw_LOCAL" in stem:
            fails.append("shared file name does not end in _MAP_v1")
    declared = set(doc.get("disclosure", {}).get("count_fields", []))
    detected = set(detect_count_fields(doc))
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
                if str(k).startswith("pct_") and v != "suppressed":
                    base = "surgeries" if k == "pct_surgeries" else k[4:]
                    num, den = o.get(base), o.get("n")
                    if not (isinstance(num, int) and isinstance(den, int) and den > 0 and v == round(100.0 * num / den, 1)):
                        fails.append(f"percentage published without its published count/denominator at {path}/{k}")
                walk(v, f"{path}/{k}")
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}/{i}")
        elif isinstance(o, str) and TIMESTAMP_RE.search(o):
            fails.append(f"timestamp-like string at {path}")
    walk({k: v for k, v in doc.items() if k != "disclosure"})
    for tk, keys in doc.get("frozen_medians", {}).items():
        for lab, o in keys.items():
            if o.get("n") == "<10":
                fails.append(f"frozen median published for a key under 10 ({tk})")
    try:
        if render_markdown(doc) != md:
            fails.append("markdown on disk is not the rendering of the shared JSON")
        if render_csv(doc) != csv_text:
            fails.append("frozen-medians CSV on disk is not the rendering of the shared JSON")
    except Exception as e:  # noqa: BLE001
        fails.append(f"shared files cannot be re-rendered: {e}")
    for row in csv.DictReader(io.StringIO(csv_text)):
        if row.get("train_n") == "<10" or (row.get("train_n", "").isdigit() and int(row["train_n"]) < MIN_CELL):
            fails.append("frozen-medians CSV has a key with n < 10")
    text = json.dumps(doc, ensure_ascii=False) + md + csv_text
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


def publish(A, cfg: BandConfig, doc: dict, sentinels: Iterable[str] = ()) -> None:
    finals = (cfg.shared_json, cfg.shared_md, cfg.shared_csv)
    pend = tuple(_pending(p) for p in finals)
    _remove(*finals, *pend)
    _atomic_write(pend[0], json.dumps(doc, indent=1, ensure_ascii=False))
    _atomic_write(pend[1], render_markdown(doc))
    _atomic_write(pend[2], render_csv(doc))
    fails = release_gate(A, *pend, "pending", sentinels)
    if fails:
        _remove(*pend)
        raise ReleaseGateError("; ".join(fails[:20]))
    final = json.loads(pend[0].read_text(encoding="utf-8"))
    final["release_gate"] = "ok"
    _atomic_write(cfg.shared_json, json.dumps(final, indent=1, ensure_ascii=False))
    _atomic_write(cfg.shared_md, render_markdown(final))
    _atomic_write(cfg.shared_csv, render_csv(final))
    _remove(*pend)
    fails = release_gate(A, *finals, "ok", sentinels)
    if fails:
        _remove(*finals)
        raise ReleaseGateError("; ".join(fails[:20]))


# =============================================================================
# Orchestration
# =============================================================================
def run_audit(cfg: BandConfig, sentinels: Iterable[str] = ()) -> dict:
    A = load_v17(cfg.v17_module_path, cfg.v17_module_sha256)
    try:
        acfg = v17_config(cfg, A)
        inp = A.load_inputs(acfg)                 # 4 full-hash gates (before reading), headers, counts, splits
        frame, _meta = A.build_frame(inp)         # ROW_ID key, v1.7 linkage, map, duration rule
        rc = {c: int((frame["link_class"] == c).sum()) for c in LINK_CLASSES}
        if sum(rc.values()) != cfg.expected_surgery_rows:
            raise AuditError("link classes do not sum to the SURGERY_INFO row count")
        df = prepare(frame, inp["map"], A.UNCLASSIFIED)
        an = BandAnalysis(A, df, rc)
        an.run()
        status, supp = A.apply_suppression(an.reg)
    except A.HashMismatchError as e:
        raise HashMismatchError(str(e)) from e
    except A.AuditError as e:
        raise AuditError(str(e)) from e
    except A.ReleaseGateError as e:
        raise ReleaseGateError(str(e)) from e
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    local = {"WARNING": "LOCAL ONLY - unsuppressed aggregates; never copy, upload or share",
             "raw_counts": an.reg.raw, "extras": an.local, "checks": inp["info"]["checks"], "link_class_counts": rc}
    _atomic_write(cfg.local_json, json.dumps(local, indent=1, default=str))
    rows = an.local[f"frozen_medians_{tkey(PRIMARY_MIN_N)}"] + an.local[f"frozen_medians_{tkey(SENSITIVITY_MIN_N)}"]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=CSV_COLUMNS, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    _atomic_write(cfg.local_csv, buf.getvalue())
    doc = build_document(A, an, status, supp, inp["info"]["hashes"])
    publish(A, cfg, doc, sentinels)
    final = json.loads(cfg.shared_json.read_text(encoding="utf-8"))
    print_console_summary(cfg, final, inp["info"]["checks"])
    return final


def _pct(v: Any) -> str:
    return f"{v}%" if isinstance(v, (int, float)) else str(v)


def print_console_summary(cfg: BandConfig, doc: dict, checks: dict) -> None:
    """Short; prints ONLY published (suppressed) values, so it can be pasted back."""
    if not cfg.verbose:
        return
    print(f"=== {doc['spec']} | {doc['map_label']} | release_gate={doc['release_gate']} ===")
    print("hashes OK: " + ", ".join(f"{k} {v[:8]}" for k, v in doc["input_hashes"].items()) + f", v1.7 module {doc['v17_module']['sha256'][:8]}")
    print(f"harmonization columns: t0={checks['harmonization_columns']['t0']!r}, OUT={checks['harmonization_columns']['OUT']!r}")
    for sp in SPLITS:
        p = doc["population"][sp]
        ms = doc["median_source"][tkey(PRIMARY_MIN_N)][sp]
        o = doc["bands"][tkey(PRIMARY_MIN_N)][sp]["overall"]
        print(f"{sp:10s} valid={p['valid']} invalid={p['invalid_duration']} gt24h={p['d_gt_24h']} | "
              + " ".join(f"{lv}={ms[lv]['surgeries']}" for lv in LEVELS)
              + " | " + " ".join(f"{b}={_pct(o.get(f'pct_{b}'))}" for b in BANDS))
    for sp in SPLITS[1:]:
        u = doc["unseen_in_train"][sp]
        print(f"unseen in train ({sp}): surgeries={u['surgeries']} ({_pct(u['pct_surgeries'])})")
    s = doc["suppression"]
    print(f"suppression: primary={s['primary_cells']} complementary={s['complementary_cells']} recoverable={s['recoverable']}")
    print(f"shared: {cfg.shared_json.name}, {cfg.shared_md.name}, {cfg.shared_csv.name}")
    print("release gate ok")


def build_config(argv: list[str] | None = None) -> BandConfig:
    here = Path(__file__).resolve()
    default_root = Path(os.environ.get("FOUNDATIONS_ROOT", str(here.parents[1])))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", default=str(default_root))
    ap.add_argument("--raw-dir", default=os.environ.get("PIC_RAW_DIR"))
    ap.add_argument("--split-file", default=os.environ.get("SPLIT_FILE_PATH"))
    ap.add_argument("--canonical", default=os.environ.get("CANONICAL_PATH"))
    ap.add_argument("--harmonization", default=os.environ.get("HARMONIZATION_PATH"))
    ap.add_argument("--map", default=os.environ.get("COMPLEXITY_MAP_PATH"))
    ap.add_argument("--out-dir", default=os.environ.get("AUDIT_OUT_DIR"))
    ap.add_argument("--harm-t0-col", default=os.environ.get("HARM_T0_COL"))
    ap.add_argument("--harm-out-col", default=os.environ.get("HARM_OUT_COL"))
    ap.add_argument("--split-col", default=os.environ.get("SPLIT_COL", "split"))
    ap.add_argument("--v17-module", default=os.environ.get("V17_MODULE_PATH", str(here.parent / V17_MODULE_FILE)))
    a = ap.parse_args(argv)
    if not a.raw_dir or not a.split_file:
        raise AuditError("--raw-dir and --split-file are required (or set PIC_RAW_DIR and SPLIT_FILE_PATH); no local default path is shipped")
    root = Path(a.project_root)
    interim = root / "data" / "interim"
    raw = Path(a.raw_dir)
    return BandConfig(
        canonical_path=Path(a.canonical) if a.canonical else interim / "canonical_icu_encounters.csv",
        harmonization_path=Path(a.harmonization) if a.harmonization else interim / "encounter_harmonization_fields.csv",
        map_path=Path(a.map) if a.map else root / "config" / "surgery_complexity_map_v1.csv",
        split_path=Path(a.split_file), surgery_info_path=raw / "SURGERY_INFO.csv", admissions_path=raw / "ADMISSIONS.csv",
        out_dir=Path(a.out_dir) if a.out_dir else interim, v17_module_path=Path(a.v17_module),
        harm_t0_col=a.harm_t0_col, harm_out_col=a.harm_out_col, split_col=a.split_col)


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
