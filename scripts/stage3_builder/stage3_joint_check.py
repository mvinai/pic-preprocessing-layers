#!/usr/bin/env python3
"""stage3_joint_check.py: independent joint (Stage 3 + Stage 2) disclosure check for a Stage 3 v1.0.2+
audit (v1.0.4: links file v3 and the optional constants file; test helper; SYNTHETIC FIXTURE USE ONLY).

Written separately from the builder and the pinned suppression module: it rebuilds the attacker's equation system
from the published files only (the Stage 3 audit JSON, the Stage 2 release JSON, the links file v2/v3) plus the Stage 3
audit layout a reader knows from the v1.0.1 release and the spec (withheld cells come back as unknowns with their
layout identities), then computes exact LP min / max (scipy HiGHS) of every unknown cell.
Reports every '<10' cell (Stage 3 or Stage 2) and every hidden / withheld cell that the joint system narrows beyond
a public-information baseline (Stage 2 + the public Stage 3 totals: cohort, landmarks, raw rows, Layer C surgeries).
Constants file (optional, v1.0.4): Stage 2 cells known from outside the release; integer rows are known values,
'derived' rows are cells the reader can compute exactly (not judged as findings).
Usage: stage3_joint_check.py <stage3_audit.json> <stage2_release.json> <stage3_links_v3.csv> [<constants.csv>]
Exit 0: nothing narrowed beyond the baseline; 1: findings."""
import csv
import json
import math
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import linprog
from scipy import sparse

S3_PUBLIC = ("/cohort/", "/landmarks/", "/raw_rows/", "/surgery/timing/n_layer_c")


def get(doc, ptr):
    d = doc
    for k in ptr.split("/")[1:]:
        k = k.replace("~1", "/").replace("~0", "~")
        if not isinstance(d, dict) or k not in d:
            raise KeyError(ptr)
        d = d[k]
    return d


def main(a3, a2, lk, consts=None):
    d3, d2 = json.load(open(a3, encoding="utf-8")), json.load(open(a2, encoding="utf-8"))
    kind, val, rels = {}, {}, []

    def cell(key, v):
        if isinstance(v, int) and not isinstance(v, bool):
            kind[key], val[key] = "K", v
        elif v == "<10":
            kind[key] = "P"
        elif v == "S":
            kind[key] = "S"
        else:
            kind[key] = "H"
    for p in d3["disclosure"]["count_fields"]:
        cell(p, get(d3, p))
    for h in d3["disclosure"]["hidden_cells"]:
        kind[h] = "H"
    for t, ps in d3["disclosure"]["identities"]:
        rels.append((t, list(ps)))
    # withheld cells and their layout identities (known from the v1.0.1 audit layout / spec), as unknowns
    def W(p):
        if p not in kind:
            kind[p] = "W"
        return p
    links = [r for r in csv.DictReader(l for l in open(lk, encoding="utf-8") if not l.startswith("#"))]
    ptrs = [c for r in links for c in [r["total"]] + [q for q in r["parts"].split(";") if q]]
    plaus = sorted({c.split("/")[2] for c in ptrs if c.startswith("/plausibility/")})
    for ch in plaus:
        b = f"/plausibility/{ch}"
        if f"{b}/train/n_rows_in_gate" in kind:
            continue
        # layout (spec §9, Nov 17, 2026): numeric = in gate + out of gate (+ zeros read as missing where the overall
        # block has a zero cell, i.e. respiratory rate); the train block is withheld as a whole
        zero = f"{b}/overall/n_rows_zero_as_missing" in kind
        names = ("in_gate", "out_of_gate", "numeric") + (("zero_as_missing",) if zero else ())
        for nm in names:
            hv = W(f"hidden:{ch}_{nm}_validation_plus_test")
            kind[hv] = "H"
            rels.append((f"{b}/overall/n_rows_{nm}", [W(f"{b}/train/n_rows_{nm}"), hv]))
        rels.append((f"{b}/train/n_rows_numeric", [f"{b}/train/n_rows_in_gate", f"{b}/train/n_rows_out_of_gate"]
                     + ([f"{b}/train/n_rows_zero_as_missing"] if zero else [])))
    if "/admission_interval/withheld/n_missing" not in kind:
        rels.append(("/cohort/n_encounters", ["/admission_interval/n_negative_clipped_to_0",
                                              W("/admission_interval/withheld/n_nonnegative"),
                                              W("/admission_interval/withheld/n_missing")]))
    for t in ("CHARTEVENTS", "LABEVENTS", "PRESCRIPTIONS", "MICROBIOLOGYEVENTS"):
        if f"/event_linkage/{t}/n_rows_in_stay" not in kind:
            rels.append((f"/raw_rows/{t}/n_rows", [W(f"/event_linkage/{t}/n_rows_in_stay"),
                                                   W(f"/event_linkage/{t}/n_rows_outside_every_window")]))
    if "/age/invalid_detail/n_missing" not in kind:
        det = [W(f"/age/invalid_detail/n_{c}") for c in ("missing", "unparseable", "negative")]
        rels.append(("/age/n_invalid", det))
        rels.append((det[2], [W(f"/age/invalid_detail/negative_magnitude/n_{b}") for b in ("lt1d", "1to30d", "gt30d")]))
    if "/weight_source/by_source/n_cn_ref" not in kind and "/weight_source/n_imputed" in kind:
        rels.append(("/weight_source/n_imputed", [W(f"/weight_source/by_source/n_{s}") for s in (
            "cn_ref", "who_ref", "train_tier_age_sex", "train_tier_age_band", "train_tier_sex_only",
            "train_tier_overall")]))
    # Stage 2 (whole released registry)
    dis2 = d2["disclosure"]
    for p in dis2["count_fields"]:
        cell("s2:" + p, get(d2, p))
    for h in dis2["hidden_cells"]:
        kind["s2:" + h] = "H"
    for t, ps in dis2["identities"]:
        for c in [t] + list(ps):
            kind.setdefault("s2:" + c, "H")
        rels.append(("s2:" + t, ["s2:" + c for c in ps]))
    # v1.0.4: constants known from outside the Stage 2 release
    skip = set()
    if consts:
        for r in csv.DictReader(l for l in open(consts, encoding="utf-8") if not l.startswith("#")):
            c, v = r["pointer"].strip(), r["value"].strip()
            if v == "derived":
                skip.add(c)
            else:
                kind[c], val[c] = "K", int(v)
    # declared links
    for j, r in enumerate(links):
        t = r["total"]
        ps = [q for q in r["parts"].split(";") if q]
        for c in [t] + ps:
            if c not in kind:
                W(c)
        if r["kind"] == "le":
            sl = f"slack:{j}"
            kind[sl] = "H"
            ps = ps + [sl]
        rels.append((t, ps))

    def solve(known_override):
        kk = dict(kind)
        kk.update(known_override)
        unk = sorted(c for c, k in kk.items() if k != "K")
        col = {c: i for i, c in enumerate(unk)}
        ri, ci, dv, rhs = [], [], [], []
        r = 0
        for t, ps in rels:
            coef = defaultdict(float)
            coef[t] += 1
            for q in ps:
                coef[q] -= 1
            b = 0.0
            row = {}
            for c, a in coef.items():
                if a == 0:
                    continue
                if kk[c] == "K":
                    b -= a * val[c]
                else:
                    row[c] = a
            if not row:
                if abs(b) > 1e-6:
                    raise SystemExit(f"inconsistent published values at {t}")
                continue
            for c, a in row.items():
                ri.append(r), ci.append(col[c]), dv.append(a)
            rhs.append(b)
            r += 1
        A = sparse.csr_matrix((dv, (ri, ci)), shape=(r, len(unk)))
        bnd = [(0, 9) if kk[c] == "P" else (10, None) if kk[c] == "S" else (0, None) for c in unk]
        used = set(ci)
        out = {}
        for c in unk:
            if col[c] not in used:
                out[c] = (bnd[col[c]][0], math.inf if bnd[col[c]][1] is None else bnd[col[c]][1])
                continue
            e = np.zeros(len(unk))
            e[col[c]] = 1
            lo = linprog(e, A_eq=A, b_eq=rhs, bounds=bnd, method="highs")
            hi = linprog(-e, A_eq=A, b_eq=rhs, bounds=bnd, method="highs")
            if lo.status == 2:
                raise SystemExit("infeasible: the published files contradict the declared links")
            l = math.ceil(lo.x[col[c]] - 1e-6)
            h = math.inf if hi.status == 3 or hi.x[col[c]] > 1e8 else math.floor(hi.x[col[c]] + 1e-6)
            out[c] = (l, h)
        return out
    joint = solve({})
    base_over = {c: "W" for c, k in kind.items()
                 if not c.startswith(("s2:", "slack:")) and k in ("K", "P", "S") and not c.startswith(S3_PUBLIC)}
    base = solve(base_over)

    def unsafe(k, lo, hi):
        n = math.inf if hi == math.inf else hi - lo + 1
        if k == "P":
            return n < 10
        if k == "S":
            return not (lo <= 9 and hi >= 10) and n < 10
        return not (lo >= 10 or hi <= 0 or (lo <= 9 and hi >= 10) or n >= 10)
    findings, compromised = [], 0
    for c, (lo, hi) in sorted(joint.items()):
        k = kind[c]
        if c in skip:
            continue
        if not unsafe(k, lo, hi):
            continue
        b = base.get(c, (0, math.inf))
        bl, bh = b
        if k == "P":
            bl, bh = max(bl, 0), min(bh, 9)
        if k == "S":
            bl = max(bl, 10)
        if lo > bl or hi < bh:
            findings.append((c, k, (lo, hi), b))
        else:
            compromised += 1
    n_p = sum(1 for k in kind.values() if k == "P")
    print(f"independent joint check: {len(kind)} cells, {len(rels)} relations, {n_p} '<10' cells; "
          f"{compromised} cell(s) narrow only as far as public information alone (Stage 2 + public totals) already "
          f"does; {len(findings)} cell(s) narrowed further by Stage 3")
    for c, k, j, b in findings:
        print(f"  NARROWED {k} {c}: joint {j} vs baseline {b}")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
