#!/usr/bin/env python3
"""synthetic_merge_v105.py: SYNTHETIC reproducer of the v1.0.3 real-run stall (Nov 17, 2026). All numbers invented.

Base: synthetic_scale_v101 (the ~105k-cell v1.0.1-shaped output). Added on top, shaped like the v1.0.3 L5 / L6 links
(names and structure only, never real values):
  * strata denominators /strata_denominators/{var}/{split}/{L}h/{lvl}/n_encounters: each the parent (subset) of the
    ~180 channel counts of its stratum, so suppressing one joins ~180 hidden differences into one component;
  * a layer_a breakdown per (split, landmark, var): /layer_a/{split}/{L}h/{var}/n_{lvl} with n_total = sum of levels,
    protected v1.0.1-style (one S next to a '<10' level), and the three var totals published;
  * links (CSV, links-v3 format): layer_a level == strata denominator (L5-like), and layer_a totals of the three vars
    == the eligible denominator (L6-like). A layer_a S cell is then pinned by its published strata twin and the '<10'
    level next to it is recovered by subtraction: the input-audit failures that start the suppression rounds.
Usage: python synthetic_merge_v105.py OUT_DIR [--scale 1.0] [--p-small 0.12]
Writes feature_sparsity_frequency_audit.json / .md and links_merge_v105.csv."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import stage2_v1_0_3_suppress as T
import synthetic_scale_v101 as SC
import synthetic_v101_output as SY

SPLITS = ("train", "validation", "test")
LS = (24, 48, 72, 96, 120, 144, 168)
VARS = ("age_band", "surgical", "first_careunit")
LVLS = ("a", "b", "c")
_CFG = {"p_small": 0.12, "seed": 23, "links": [], "hub_frac": 0.0}


class MergeBuilder(SY.Builder):
    def finish(self) -> dict:
        if "/surgery/giant_0/n_total_courses" not in self.raw:
            return super().finish()          # the small embedded synthetic (synthetic_v101_output.build)
        rng = np.random.default_rng(_CFG["seed"])
        links = _CFG["links"]
        index: dict[tuple, list[str]] = {}
        for p in self.raw:
            if p.startswith("/channels/") and p.endswith("/n_avail_full_72h_release0h"):
                q = p.split("/")       # '', channels, chan, split, landmarks, Lh, strata, var, lvl, leaf
                index.setdefault((q[3], q[5], q[7], q[8]), []).append(p)
        hub = {(var, L): bool(rng.random() < _CFG["hub_frac"]) for var in VARS for L in LS}
        made: dict[tuple, str] = {}
        for s in SPLITS:
            for L in LS:
                elig = f"/denominators/{s}/eligible/{L}h/n_eligible"
                E = self.raw[elig]
                for var in VARS:
                    # children of each stratum (the channel counts already built by synthetic_scale_v101)
                    kids = {lvl: index.get((s, f"{L}h", var, lvl), []) for lvl in LVLS}
                    need = {lvl: max([self.raw[k] for k in kids[lvl]] or [0]) for lvl in LVLS}
                    small = rng.random() < _CFG["p_small"]
                    # level sizes: >= every child; one level may be small (< 10) only if it has no big child
                    vals = {}
                    for i, lvl in enumerate(LVLS):
                        vals[lvl] = need[lvl] + int(rng.integers(0, 40))
                    if small:
                        lv = LVLS[int(rng.integers(0, 3))]
                        if need[lv] < 10:
                            vals[lv] = int(rng.integers(max(need[lv], 1), 10))
                    tot = sum(vals.values())
                    # stratum denominators: parents of their channel children
                    strata = {}
                    for lvl in LVLS:
                        p = self.cell(f"/strata_denominators/{var}/{s}/{L}h/{lvl}/n_encounters", vals[lvl],
                                      force="S" if hub[(var, L)] and vals[lvl] >= 10 else None)
                        made[(var, L, lvl, s)] = p
                        strata[lvl] = p
                        for k in kids[lvl]:
                            h = f"hidden:{k}|complement|stratum"
                            self.hidden.append(h)
                            self.raw[h] = vals[lvl] - self.raw[k]
                            self.identities.append([p, [k, h]])
                    # layer_a breakdown, v1.0.1-protected on its own (one S next to a '<10' level)
                    smalls = [l for l in LVLS if vals[l] < 10]
                    bigs = sorted((vals[l], l) for l in LVLS if vals[l] >= 10)
                    force = {bigs[0][1]: "S"} if smalls and bigs else {}
                    ids = self.partition(f"/layer_a/{s}/{L}h/{var}", {lvl: vals[lvl] for lvl in LVLS}, force=force)
                    for lvl in LVLS:
                        links.append((ids[lvl], [strata[lvl]], "L5-like: layer_a level == strata denominator"))
                    # var totals (the population with a recorded value; up to the eligible count)
                    if tot <= E:
                        h = f"hidden:/layer_a/{s}/{L}h/{var}/n_total|complement|eligible"
                        self.hidden.append(h)
                        self.raw[h] = E - tot
                        self.identities.append([elig, [ids["total"], h]])
        for (var, L), on in sorted(hub.items()):
            if not on:
                continue
            for lvl in LVLS:          # L6-like: the overall stratum = sum of the three splits (published)
                parts = [made[(var, L, lvl, s)] for s in SPLITS]
                t = self.cell(f"/strata_denominators/{var}/overall/{L}h/{lvl}/n_encounters",
                              sum(self.raw[q] for q in parts))
                self.ident(t, parts)
        return super().finish()


def write(folder: Path, scale: float = 1.0, p_small: float = 0.12, seed: int = 23, hub_frac: float = 0.0) -> dict[str, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    _CFG.update(p_small=p_small, seed=seed, links=[], hub_frac=hub_frac)
    old = SC.SY.Builder
    SC.SY.Builder = MergeBuilder
    try:
        doc = SC.build(scale)
    finally:
        SC.SY.Builder = old
    pj = folder / "feature_sparsity_frequency_audit.json"
    pm = folder / "feature_sparsity_frequency_audit.md"
    pl = folder / "links_merge_v105.csv"
    pj.write_text(T.dump_json(doc), encoding="utf-8", newline="")
    pm.write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    with open(pl, "w", encoding="utf-8", newline="") as f:
        f.write("# SYNTHETIC links for synthetic_merge_v105 (L5 / L6-like)\ntotal,parts,note\n")
        for t, parts, note in _CFG["links"]:
            f.write(f"{t},{';'.join(parts)},{note}\n")
    return {"json": pj, "md": pm, "links": pl}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--p-small", type=float, default=0.12)
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--hub-frac", type=float, default=0.0)
    a = ap.parse_args()
    out = write(Path(a.out_dir), a.scale, a.p_small, a.seed, a.hub_frac)
    print({k: (str(v), v.stat().st_size) for k, v in out.items()})
