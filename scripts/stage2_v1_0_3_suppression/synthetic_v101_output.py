"""Synthetic stand-in for the released Stage 2 v1.0.1 aggregate output (JSON + markdown), for tests only.

Every number here is invented. The structure mimics feature_sparsity_frequency_audit.json v1.0.1: count fields
('n_' keys) holding an integer >= 10, '<10' or 'S'; a disclosure block with count_fields, identities (total = sum of
parts), hidden_cells (unpublished parent-minus-child differences), ratios, quantiles and dependents; and the
markdown rendered from the JSON. The planted cases reproduce the v1.0.1 differencing exposure:

  fio2-like row status   total and measurement rows published beside five '<10' exclusion cells (sum small)
  unit totals            two unit tables whose totals equal the fio2-like row-status total (declared links)
  S/F-like derived count a published draws total (frozen) = two published parts + two '<10' parts (sum small)
  linked table           2 x 2 table with row, column and grand margins; one '<10' cell (row AND column groups)
  second cell needed     a frozen total with parts 12, 13 and a '<10' part: one complementary cell is not enough
  safe group             an aa_gradient-like row with an S cell that already protects its '<10' cells
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import stage2_v1_0_3_suppress as T

SPEC_IN = T.EXPECT_INPUT_SPEC


class Builder:
    def __init__(self):
        self.doc: dict = {
            "spec": SPEC_IN, "release_gate": "ok", "outcome_descriptors_flag": True,
            "input_hashes": {"canonical": "0" * 64, "split": "1" * 64},
            "raw_tables": {"LABEVENTS": {"file": "LABEVENTS.csv", "sha256": "2" * 64, "gate": "stage1_sha256"}},
            "cohort_checks": {}, "provisional_classes": {"channels": {}}, "notes": {"synthetic": "test data only"},
        }
        self.raw: dict[str, int] = {}
        self.forced: dict[str, str] = {}
        self.identities: list = []
        self.hidden: list[str] = []
        self.ratios: list = []
        self.quantiles: list = []
        self.dependents: list = []
        self.pending_derived: list = []

    def cell(self, ptr: str, raw: int, force: str | None = None) -> str:
        self.raw[ptr] = raw
        if force:
            self.forced[ptr] = force
        return ptr

    def ident(self, total: str, parts: list[str]) -> None:
        assert self.raw[total] == sum(self.raw[p] for p in parts), (total, parts)
        self.identities.append([total, list(parts)])

    def subset(self, child: str, parent: str) -> None:
        h = f"hidden:{child}|complement|sub"
        self.hidden.append(h)
        self.raw[h] = self.raw[parent] - self.raw[child]
        assert self.raw[h] >= 0
        self.identities.append([parent, [child, h]])

    def partition(self, base: str, parts: dict[str, int], force: dict[str, str] | None = None,
                  total_name: str = "n_total") -> dict[str, str]:
        force = force or {}
        ids = {k: self.cell(f"{base}/n_{k}", v, force.get(k)) for k, v in parts.items()}
        t = self.cell(f"{base}/{total_name}", sum(parts.values()), force.get("total"))
        self.ident(t, list(ids.values()))
        ids["total"] = t
        return ids

    def status(self, ptr: str) -> str:
        if ptr in self.forced:
            return self.forced[ptr]
        return "P" if self.raw[ptr] < 10 else "pub"

    def finish(self) -> dict:
        doc = self.doc
        pub = set()
        n_p = n_s = 0
        for ptr, v in self.raw.items():
            if ptr.startswith("hidden:"):
                continue
            st = self.status(ptr)
            if st == "pub":
                pub.add(ptr)
                val: object = v
            elif st == "P":
                val, n_p = "<10", n_p + 1
            else:
                val, n_s = "S", n_s + 1
            parts = T.ptr_parts(ptr)
            cur = doc
            for p in parts[:-1]:
                cur = cur.setdefault(p, {})
            cur[parts[-1]] = val
        for ptr, num, den in self.ratios:
            ok = num in pub and den in pub
            T_set(doc, ptr, T.r1(100.0 * self.raw[num] / self.raw[den]) if ok else "suppressed")
        for ptr, n, q in self.quantiles:
            T_set(doc, ptr, q if (n in pub and self.raw[n] >= 10) else "suppressed")
        for ptr, needs, v in self.dependents:
            T_set(doc, ptr, v if all(x in pub for x in needs) else "suppressed")
        doc["disclosure"] = {
            "count_fields": T.detect_count_fields(doc),
            "identities": self.identities,
            "hidden_cells": sorted(self.hidden),
            "ratios": [[a, b, c] for a, b, c in self.ratios],
            "dependents": [[a, list(b)] for a, b, _ in self.dependents],
            "quantiles": [[a, b] for a, b, _ in self.quantiles],
        }
        doc["suppression"] = {"primary_cells": n_p, "complementary_cells": n_s,
                              "recoverability_cells_checked": 0, "recoverable": 0}
        return doc


def T_set(doc: dict, ptr: str, v) -> None:
    parts = T.ptr_parts(ptr)
    cur = doc
    for p in parts[:-1]:
        cur = cur.setdefault(p, {})
    cur[parts[-1]] = v


Q = {"p10": 0.0, "p25": 0.1, "p50": 0.4, "p75": 1.0, "p90": 2.0}


def build(include_linked_table: bool = True, include_second_cell: bool = True) -> dict:
    b = Builder()
    # denominators (frozen: cohort and split sizes are public)
    splits = {"train": 700, "validation": 150, "test": 150}
    tot = b.cell("/denominators/n_encounters", sum(splits.values()))
    enc = {s: b.cell(f"/denominators/{s}/n_encounters", n) for s, n in splits.items()}
    b.ident(tot, list(enc.values()))
    elig = {"train": 600, "validation": 120, "test": 125}
    for s in splits:
        e = b.cell(f"/denominators/{s}/eligible/24h/n_eligible", elig[s])
        b.subset(e, enc[s])
    b.doc["cohort_checks"] = {"encounters": 1000, "eligible_at_landmark": {"24h": sum(elig.values())}}
    # channels (rendered in the markdown)
    for name, ever, everb, av, cls in (("alpha", 640, 600, 500, "value_mask_time_since_last"),
                                       ("beta", 120, 90, 80, "measured_indicator_only")):
        base = f"/channels/{name}/train"
        ne = b.cell(f"{base}/n_ever_336h", ever)
        b.subset(ne, enc["train"])
        nb = b.cell(f"{base}/n_ever_layer_b", everb)
        b.subset(nb, ne)
        na = b.cell(f"{base}/landmarks/24h/n_avail_full", av)
        b.subset(na, f"/denominators/train/eligible/24h/n_eligible")
        b.ratios.append((f"{base}/pct_ever_336h", ne, enc["train"]))
        nr = b.cell(f"{base}/n_rate_per_day_all", ever)
        b.quantiles.append((f"{base}/q_rate_per_day_all", nr, Q))
        b.doc["provisional_classes"]["channels"][name] = {"class": cls}
    # fio2-like row status: every exclusion cell '<10', total minus measurement rows small
    fio2 = b.partition("/map_checks/row_status/fio2", {"measurement": 1707, "out_of_gate": 1, "nonnumeric": 0,
                                                       "blank": 0, "unit_pending": 0, "unit_unrecognized": 0})
    # aa_gradient-like row status: an S cell already protects the '<10' cells (must stay untouched)
    b.partition("/map_checks/row_status/aa_like", {"measurement": 777, "out_of_gate": 2, "nonnumeric": 22,
                                                   "blank": 0, "unit_pending": 1, "unit_unrecognized": 3},
                force={"nonnumeric": "S"})
    # unit tables whose totals equal the fio2-like total (links)
    b.partition("/map_checks/units/LABEVENTS_5222", {"(blank)": 14, "%": 1694})
    b.partition("/map_checks/units/LABEVENTS_5236", {"(blank)": 12, "mmhg": 1696})
    # S/F-like: draws total (frozen) = valid P/F (frozen) + no-SpO2 + two '<10'
    sf = b.cell("/map_checks/sf_ratio/n_fio2_draws", 1707)
    parts = {"with_valid_pf": 1214, "no_spo2_within_5min": 489, "spo2_gt_97": 3, "sf_computed": 1}
    ids = [b.cell(f"/map_checks/sf_ratio/fio2_draws/n_{k}", v) for k, v in parts.items()]
    b.ident(sf, ids)
    b.partition("/map_checks/pf_ratio/same_draw_pairs", {"valid": 1214, "outside_gate": 491})
    if include_linked_table:
        base = "/map_checks/linked_demo"
        c = {("r1", "c1"): 3, ("r1", "c2"): 50, ("r2", "c1"): 40, ("r2", "c2"): 60}
        cells = {k: b.cell(f"{base}/{k[0]}/n_{k[1]}", v) for k, v in c.items()}
        rows = {r: b.cell(f"{base}/{r}/n_total", sum(v for k, v in c.items() if k[0] == r)) for r in ("r1", "r2")}
        cols = {cc: b.cell(f"{base}/col_totals/n_{cc}", sum(v for k, v in c.items() if k[1] == cc)) for cc in ("c1", "c2")}
        grand = b.cell(f"{base}/n_grand_total", sum(c.values()))
        for r in rows:
            b.ident(rows[r], [cells[(r, "c1")], cells[(r, "c2")]])
        for cc in cols:
            b.ident(cols[cc], [cells[("r1", cc)], cells[("r2", cc)]])
        b.ident(grand, list(rows.values()))
        b.ident(grand, list(cols.values()))
        # derived values resting on table cells
        b.ratios.append((f"{base}/r1/pct_c2", cells[("r1", "c2")], rows["r1"]))
        nv = b.cell(f"{base}/r2/n_values", 60)
        b.ident(nv, [cells[("r2", "c2")]])
        b.quantiles.append((f"{base}/r2/q_values", nv, Q))
        b.dependents.append((f"{base}/x_share_c1", [cols["c1"], grand], 0.3))
    if include_second_cell:
        b.partition("/map_checks/map/pairs", {"valid": 12, "dbp_ge_sbp": 13, "outside_gate": 4})
    return {"doc": b.finish(), "raw": b.raw}


LINKS = [
    ("/map_checks/row_status/fio2/n_total", "/map_checks/units/LABEVENTS_5222/n_total", "fio2-like = item 5222 rows"),
    ("/map_checks/units/LABEVENTS_5236/n_total", "/map_checks/row_status/fio2/n_total", "same analyser draws"),
    ("/cohort_checks/eligible_at_landmark/24h", "/denominators/train/eligible/24h/n_eligible;"
     "/denominators/validation/eligible/24h/n_eligible;/denominators/test/eligible/24h/n_eligible",
     "cohort_checks total across splits"),
]
FROZEN = [
    ("/denominators/n_encounters", "cohort size public"),
    ("/denominators/train/n_encounters", "split size public"),
    ("/denominators/validation/n_encounters", "split size public"),
    ("/denominators/test/n_encounters", "split size public"),
    ("/map_checks/row_status/fio2/n_measurement", "public in the follow-up output"),
    ("/map_checks/sf_ratio/n_fio2_draws", "same value as the measurement rows"),
    ("/map_checks/sf_ratio/fio2_draws/n_with_valid_pf", "public in the follow-up output"),
    ("/map_checks/pf_ratio/same_draw_pairs/n_valid", "public in the follow-up output"),
    ("/map_checks/pf_ratio/same_draw_pairs/n_outside_gate", "public in the follow-up output"),
    ("/map_checks/pf_ratio/same_draw_pairs/n_total", "public in the follow-up output"),
    ("/map_checks/map/pairs/n_total", "synthetic: an already-public total (forces a second cell)"),
]


def write_case(folder: Path, **kw) -> dict[str, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    doc = build(**kw)["doc"]
    pj = folder / "feature_sparsity_frequency_audit.json"
    pm = folder / "feature_sparsity_frequency_audit.md"
    pj.write_text(T.dump_json(doc), encoding="utf-8", newline="")
    pm.write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    pl = folder / "links.csv"
    with open(pl, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["total", "parts", "note"])
        w.writerows(LINKS)
    pf = folder / "frozen.csv"
    with open(pf, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pointer", "reason"])
        frozen = FROZEN if kw.get("include_second_cell", True) else FROZEN[:-1]
        w.writerows(frozen)
    return {"json": pj, "md": pm, "links": pl, "frozen": pf}


if __name__ == "__main__":
    import sys
    out = write_case(Path(sys.argv[1] if len(sys.argv) > 1 else "synthetic_case"))
    print(json.dumps({k: str(v) for k, v in out.items()}, indent=1))
