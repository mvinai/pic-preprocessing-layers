"""Large synthetic stand-in (about 57 MB) for the Stage 2 v1.0.1 shared JSON, for profiling only. All numbers invented.

The equation structure is shaped like the real output's (structure statistics, counts only, Nov 16):
~105k count fields, ~87k hidden differences, ~99k identities (mostly parent = child + hidden; ~11k one-part
equalities; a few hundred partitions of 4-14 parts; a few 500-part identities), ~16k '<10' and ~6k 'S' cells,
two giant linked components of ~7k unknown cells each (chains of S 'hub' cells, each the parent of ~200 counts),
~3k small multi-equation components and ~8k one-equation components. A few FiO2-style differencing exposures are
planted, some inside the giant components.

Usage: python synthetic_scale_v101.py OUT_DIR [--scale 1.0]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import stage2_v1_0_3_suppress as T
import synthetic_v101_output as SY

LONG = "landmarks/{L}h/strata/{var}/{lvl}"


def build(scale: float = 1.0, seed: int = 11) -> dict:
    rng = np.random.default_rng(seed)
    b = SY.Builder()
    splits = {"train": 800, "validation": 150, "test": 150}  # generic synthetic scale (not a real split size)
    tot = b.cell("/denominators/n_encounters", sum(splits.values()))
    enc = {s: b.cell(f"/denominators/{s}/n_encounters", n) for s, n in splits.items()}
    b.ident(tot, list(enc.values()))
    Ls = (24, 48, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288, 312, 336)
    elig = {}
    for s, n in splits.items():
        for L in Ls:
            e = b.cell(f"/denominators/{s}/eligible/{L}h/n_eligible", int(n * (0.98 - 0.05 * Ls.index(L))))
            b.subset(e, enc[s])
            elig[(s, L)] = e

    def val(par: int, p_small: float) -> int:
        if rng.random() < p_small:
            return int(rng.integers(0, 10))
        return int(rng.integers(10, max(11, par)))

    # 1) the bulk: channel x split x landmark x stratum counts, each a subset of its eligible denominator
    n_ch = int(180 * scale)
    vars_ = ("age_band", "surgical", "first_careunit")
    lvls = ("a", "b", "c")
    for ch in range(n_ch):
        for s in splits:
            for L in Ls[:7]:
                par = elig[(s, L)]
                for var in vars_:
                    for lvl in lvls:
                        base = f"/channels/chan_{ch:03d}/{s}/" + LONG.format(L=L, var=var, lvl=lvl)
                        v = val(b.raw[par], 0.15)
                        c = b.cell(f"{base}/n_avail_full_72h_release0h", min(v, b.raw[par]))
                        b.subset(c, par)
                        if rng.random() < 0.25:            # a one-part equality (same count published twice)
                            e2 = b.cell(f"{base}/n_avail_full_72h_release0h_copy", b.raw[c])
                            b.ident(e2, [c])
                        if rng.random() < 0.15:            # a nested subset under an unknown-able child
                            c2 = b.cell(f"{base}/n_ever_layer_b", int(b.raw[c] * rng.uniform(0.3, 1.0)))
                            b.subset(c2, c)
                pc = b.cell(f"/channels/chan_{ch:03d}/{s}/landmarks/{L}h/pct_avail_full", 0)
                del b.raw[pc]
                b.ratios.append((pc, f"/channels/chan_{ch:03d}/{s}/" + LONG.format(L=L, var=vars_[0], lvl="a")
                                 + "/n_avail_full_72h_release0h", par))
            q = f"/channels/chan_{ch:03d}/{s}/n_rate_per_day_all"
            b.cell(q, int(rng.integers(10, 900)))
            b.quantiles.append((f"/channels/chan_{ch:03d}/{s}/q_rate_per_day_all", q, SY.Q))
    # 2) partitions (row status, units, reasons): 4-14 parts, v1.0.1-style protection (one S next to small parts)
    n_part = int(800 * scale)
    for g in range(n_part):
        k = int(rng.choice([4, 5, 6, 7, 8, 14]))
        parts = {f"reason_{j}": (int(rng.integers(0, 10)) if rng.random() < 0.4 else int(rng.integers(10, 5000)))
                 for j in range(k - 1)}
        force = {}
        smalls = [x for x, v in parts.items() if v < 10]
        bigs = sorted((v, x) for x, v in parts.items() if v >= 10)
        if len(smalls) >= 1 and len(bigs) >= 2:
            force[bigs[0][1]] = "S"
        b.partition(f"/map_checks/row_status_like/group_{g:04d}", parts, force=force)
    # 3) two giant linked components: chains of S hubs, each the parent of ~200 counts
    for giant in range(2):
        top = b.cell(f"/surgery/giant_{giant}/n_total_courses", 2_000_000)
        parent = top
        for h in range(int(35 * scale) or 1):
            hub = b.cell(f"/surgery/giant_{giant}/hub_{h:02d}/n_courses", b.raw[parent] - int(rng.integers(1000, 5000)),
                         force="S")
            b.subset(hub, parent)
            for j in range(200):
                c = b.cell(f"/surgery/giant_{giant}/hub_{h:02d}/items/item_{j:03d}/n_courses_in_window",
                           val(min(b.raw[hub], 3000), 0.2))
                b.subset(c, hub)
                if j % 50 == 0:   # a planted FiO2-style exposure: published child = big part + two '<10'
                    big = b.raw[c] - 3 if b.raw[c] >= 13 else None
                    if big is not None and big >= 10:
                        b.partition(f"/surgery/giant_{giant}/hub_{h:02d}/items/item_{j:03d}/by_reason",
                                    {"main": big, "minor_a": 2, "minor_b": 1}, total_name="n_total")
            parent = hub
    # 3b) ~3,000 medium linked components: short chains of S hubs (heavy-tailed sizes, up to ~1,000 cells)
    for m in range(int(3000 * scale)):
        n_h = int(min(12, 1 + rng.pareto(1.5)))
        n_c = int(min(90, 2 + rng.pareto(1.2) * 3))
        parent = elig[("train", Ls[m % 14])]
        for h in range(n_h):
            hub = b.cell(f"/medications/group_{m:04d}/course_level_{h:02d}/n_encounters_with_course",
                         max(10, b.raw[parent] - int(rng.integers(1, 40))), force="S")
            b.subset(hub, parent)
            for j in range(n_c):
                c = b.cell(f"/medications/group_{m:04d}/course_level_{h:02d}/by_drug/drug_{j:02d}/n_encounters",
                           val(min(b.raw[hub], 400), 0.25))
                b.subset(c, hub)
            parent = hub
    # 4) three wide identities (500 parts)
    for w in range(3):
        parts = [b.cell(f"/microbiology/wide_{w}/organisms/org_{j:03d}/n_isolates",
                        int(rng.integers(0, 10)) if rng.random() < 0.3 else int(rng.integers(10, 200)))
                 for j in range(499)]
        t = b.cell(f"/microbiology/wide_{w}/n_isolates_total", sum(b.raw[p] for p in parts))
        b.ident(t, parts)
    # 5) planted FiO2 / units / S/F cases of the small synthetic (same pointers as the links file)
    small = SY.build(include_linked_table=True, include_second_cell=True)
    for ptr, v in small["raw"].items():
        if ptr.startswith("/map_checks/") and ptr not in b.raw:
            b.raw[ptr] = v
    for t, parts in small["doc"]["disclosure"]["identities"]:
        if t.startswith("/map_checks/"):
            b.identities.append([t, list(parts)])
    for c in ("/map_checks/row_status/aa_like/n_nonnumeric",):
        b.forced[c] = "S"
    b.doc["cohort_checks"] = {"encounters": b.raw[tot]}
    return b.finish()


def write(folder: Path, scale: float = 1.0) -> dict[str, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    doc = build(scale)
    pj = folder / "feature_sparsity_frequency_audit.json"
    pm = folder / "feature_sparsity_frequency_audit.md"
    pj.write_text(T.dump_json(doc), encoding="utf-8", newline="")
    pm.write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    return {"json": pj, "md": pm}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--scale", type=float, default=1.0)
    a = ap.parse_args()
    out = write(Path(a.out_dir), a.scale)
    print({k: (str(v), v.stat().st_size) for k, v in out.items()})
