"""Synthetic reproducer of the v1.1 suppression stall (all numbers invented; test data only).

One giant linked component (a chain of S 'hub' cells, each the parent of 120 counts, about half of them '<10'),
so the component is far above the 300-cell threshold where v1.1 switched to per-group fixes without look-ahead.
Two kinds of planted failing cells sit inside it:

- kind 'pinned' (the stall in the real run): a '<10' item that is also a part of a small identity whose total and
  other part are already public (frozen). No allowed suppression can open that identity. v1.1 found no candidate
  in the item's own groups, searched outward through the component and suppressed one unrelated published cell
  per group per round: the failing count stayed constant, every new cell was safe, and the giant component grew by
  one cell per round, without end (one-per-round walk through ~360 publishable cells).
- kind 'low': a '<10' item that is a part of an identity whose other parts are many counts of exactly 10 and one
  large count. v1.1's static choice (fewest groups, then smallest value) suppressed the 10s one per round; a
  suppressed 10 can only be 10 or more, so it never widens the item's interval, and the failing count grew.

Usage: python synthetic_stall_v101.py OUT_DIR [--kind pinned|low|both]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

import stage2_v1_0_3_suppress as T
import synthetic_v101_output as SY


def build(kind: str = "both", n_hubs: int = 3, n_items: int = 120, n_low: int = 40, seed: int = 5) -> tuple[dict, list]:
    rng = np.random.default_rng(seed)
    b = SY.Builder()
    frozen: list[tuple[str, str]] = []
    top = b.cell("/surgery/giant/n_total_courses", 1_000_000)
    parent = top
    hubs = []
    for h in range(n_hubs):
        hub = b.cell(f"/surgery/giant/hub_{h:02d}/n_courses", b.raw[parent] - int(rng.integers(1000, 5000)), force="S")
        b.subset(hub, parent)
        hubs.append(hub)
        for j in range(n_items):
            v = int(rng.integers(0, 10)) if rng.random() < 0.5 else int(rng.integers(10, 3000))
            c = b.cell(f"/surgery/giant/hub_{h:02d}/items/item_{j:03d}/n_courses_in_window", v)
            b.subset(c, hub)
        parent = hub
    planted = 0
    for g, hub in enumerate(hubs[:2]):
        if kind in ("pinned", "both"):
            u = b.cell(f"/surgery/giant/hub_{g:02d}/pinned_item/n_courses_in_window", 3)
            b.subset(u, hub)
            other = b.cell(f"/surgery/giant/hub_{g:02d}/pinned_item/n_other_public", 41)
            tot = b.cell(f"/surgery/giant/hub_{g:02d}/pinned_item/n_total_public", 44)
            b.ident(tot, [u, other])
            frozen += [(tot, "synthetic: already public elsewhere"), (other, "synthetic: already public elsewhere")]
            planted += 1
        if kind in ("low", "both"):
            u = b.cell(f"/surgery/giant/hub_{g:02d}/low_item/n_courses_in_window", 2)
            b.subset(u, hub)
            parts = [u] + [b.cell(f"/surgery/giant/hub_{g:02d}/low_item/by_site/site_{i:02d}/n", 10)
                           for i in range(n_low)]
            parts.append(b.cell(f"/surgery/giant/hub_{g:02d}/low_item/by_site/site_big/n", 900))
            tot = b.cell(f"/surgery/giant/hub_{g:02d}/low_item/n_total", sum(b.raw[p] for p in parts))
            b.ident(tot, parts)
            copy = b.cell(f"/surgery/giant/hub_{g:02d}/low_item/n_total_copy", b.raw[tot])
            b.ident(copy, [tot])   # the total is published twice (two groups): the static choice ranks it last
            planted += 1
    b.doc["cohort_checks"] = {"encounters": b.raw[top]}
    return b.finish(), frozen


def write(folder: Path, kind: str = "both", **kw) -> dict[str, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    doc, frozen = build(kind, **kw)
    pj = folder / "feature_sparsity_frequency_audit.json"
    pm = folder / "feature_sparsity_frequency_audit.md"
    pj.write_text(T.dump_json(doc), encoding="utf-8", newline="")
    pm.write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    pf = folder / "frozen.csv"
    with open(pf, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pointer", "reason"])
        w.writerows(frozen)
    return {"json": pj, "md": pm, "frozen": pf}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--kind", choices=("pinned", "low", "both"), default="both")
    a = ap.parse_args()
    print(write(Path(a.out_dir), a.kind))
