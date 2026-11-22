#!/usr/bin/env python3
"""stage2_independent_check.py: independent disclosure check of ONE Stage 2 release (test helper; SYNTHETIC USE).

Written separately from the suppression module (no import of it): it rebuilds a reader's equation system from the
published release JSON alone (its registry: count fields, identities, hidden cells), plus
  --links   declared equations (Stage 2 links CSV 'total,parts,note'; or a Stage 3 links file 'kind,total,parts,note',
            of which only the rows whose cells are all 's2:' Stage 2 cells are used: the Stage 2 internal equalities)
  --public  cells a reader knows from outside (pointer,value[,reason]; an empty value takes the value from --layout;
            'derived': fixed by the other public values, kept as an unknown and not rule-checked; 'publish': as empty)
  --layout  the input (v1.0.1-format) JSON: identities that mention withheld (dropped) cells come back with those
            cells as unknowns, as a reader who knows the layout would add them; also the source of empty public values
and computes the exact LP min / max (scipy HiGHS) of every unknown cell. A cell fails when
  '<10'   its interval is narrower than the full 0-9 range;    'S'  it is narrower than 10 integers and does not
  straddle 9|10;   hidden / withheld   it is not provably >= 10, provably 0, straddling 9|10 or >= 10 integers wide.
Prints counts by kind and the failing pointers (no values; --show-bounds adds the intervals: synthetic data only).
Exit 0: no cell fails; 1: at least one fails; 2: input error."""
import argparse
import csv
import json
import math
import sys
from collections import defaultdict

import numpy as np
from scipy import sparse
from scipy.optimize import linprog


def get(doc, ptr):
    d = doc
    for k in ptr.split("/")[1:]:
        k = k.replace("~1", "/").replace("~0", "~")
        if not isinstance(d, dict) or k not in d:
            raise KeyError(ptr)
        d = d[k]
    return d


def is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def read_links(path):
    rows = list(csv.DictReader(l for l in open(path, encoding="utf-8") if not l.startswith("#")))
    out = []
    for r in rows:
        t, ps = r["total"].strip(), [q.strip() for q in r["parts"].split(";") if q.strip()]
        kind = (r.get("kind") or "eq").strip()
        if "kind" in r:                       # Stage 3 links file: Stage 2 internal rows only
            if not all(c.startswith("s2:") for c in [t] + ps):
                continue
            t, ps = t[3:], [c[3:] for c in ps]
        out.append((kind, t, ps))
    return out


def unsafe(k, lo, hi):
    n = math.inf if hi == math.inf else hi - lo + 1
    if k == "P":
        return n < 10
    if k == "S":
        return not (lo <= 9 and hi >= 10) and n < 10
    return not (lo >= 10 or hi <= 0 or (lo <= 9 and hi >= 10) or n >= 10)


def check(release, links=None, public=None, layout=None, show=False, out=print):
    doc = json.load(open(release, encoding="utf-8"))
    lay = json.load(open(layout, encoding="utf-8")) if layout else None
    dis = doc["disclosure"]
    kind, val, rels = {}, {}, []
    for p in dis["count_fields"]:
        v = get(doc, p)
        if is_int(v):
            kind[p], val[p] = "K", v
        elif v == "<10":
            kind[p] = "P"
        elif v == "S":
            kind[p] = "S"
        else:
            raise SystemExit(f"input error: unexpected count label at a count field")
    for h in dis["hidden_cells"]:
        kind[h] = "H"
    for t, ps in dis["identities"]:
        rels.append((t, list(ps)))
    if lay is not None:                          # withheld cells (dropped subtrees) with their layout identities
        have = set(map(lambda r: (r[0], tuple(r[1])), dis["identities"]))
        for t, ps in lay["disclosure"]["identities"]:
            if (t, tuple(ps)) in have:
                continue
            cells = [t] + list(ps)
            if all(c in kind for c in cells):
                continue
            for c in cells:
                kind.setdefault(c, "W")
            rels.append((t, list(ps)))
        for h in lay["disclosure"]["hidden_cells"]:
            if h in kind and kind[h] == "W":
                kind[h] = "H"
    for j, (k, t, ps) in enumerate(read_links(links) if links else []):
        for c in [t] + ps:
            if c not in kind:
                raise SystemExit("input error: a links pointer is not a count field of the release")
        if k == "le":
            kind[f"slack:{j}"] = "H"
            ps = ps + [f"slack:{j}"]
        rels.append((t, ps))
    n_pub = 0
    if public:
        for r in csv.DictReader(l for l in open(public, encoding="utf-8") if not l.startswith("#")):
            p, v = r["pointer"].strip(), (r.get("value") or "").strip()
            if p not in kind:
                raise SystemExit("input error: a public pointer is not a count field of the release")
            if v == "derived":                   # fixed by the other public values: an unknown, not rule-checked
                kind[p] = "X"
                n_pub += 1
                continue
            if v in ("", "publish"):
                if lay is None:
                    raise SystemExit("input error: an empty public value needs --layout")
                v = get(lay, p)
                if not is_int(v):
                    raise SystemExit("input error: an empty public value is not published in the layout file")
            v = int(v)
            if kind[p] == "K" and val[p] != v:
                raise SystemExit("input error: a public value contradicts the release")
            kind[p], val[p] = "K", v
            n_pub += 1
    unk = sorted(c for c, k in kind.items() if k != "K")
    col = {c: i for i, c in enumerate(unk)}
    ri, ci, dv, rhs = [], [], [], []
    r = 0
    for t, ps in rels:
        coef = defaultdict(float)
        coef[t] += 1
        for q in ps:
            coef[q] -= 1
        b, row = 0.0, {}
        for c, a in coef.items():
            if a == 0:
                continue
            if kind[c] == "K":
                b -= a * val[c]
            else:
                row[c] = a
        if not row:
            if abs(b) > 1e-6:
                out("INCONSISTENT: published values contradict an identity or link")
                return 2
            continue
        for c, a in row.items():
            ri.append(r), ci.append(col[c]), dv.append(a)
        rhs.append(b)
        r += 1
    A = sparse.csr_matrix((dv, (ri, ci)), shape=(r, len(unk)))
    bnd = [(0, 9) if kind[c] == "P" else (10, None) if kind[c] == "S" else (0, None) for c in unk]
    # component decomposition (cells sharing no row are independent)
    parent = list(range(len(unk)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    A = A.tocsr()
    for i in range(r):
        cols = A.indices[A.indptr[i]:A.indptr[i + 1]]
        for c in cols[1:]:
            a, b2 = find(cols[0]), find(c)
            if a != b2:
                parent[a] = b2
    comp = defaultdict(list)
    for i in range(len(unk)):
        comp[find(i)].append(i)
    row_comp = defaultdict(list)
    for i in range(r):
        cols = A.indices[A.indptr[i]:A.indptr[i + 1]]
        if len(cols):
            row_comp[find(cols[0])].append(i)
    rhs = np.array(rhs)
    fails, by_kind = [], defaultdict(int)
    for root, cs in comp.items():
        rows = row_comp.get(root, [])
        for c in cs:
            by_kind[kind[unk[c]]] += 1
        if not rows:
            continue
        sub = A[rows][:, cs]
        sb = rhs[rows]
        bb = [bnd[c] for c in cs]
        for j, c in enumerate(cs):
            e = np.zeros(len(cs))
            e[j] = 1
            lo = linprog(e, A_eq=sub, b_eq=sb, bounds=bb, method="highs")
            if lo.status == 2:
                out("INFEASIBLE: the release contradicts the declared links / public values")
                return 2
            hi = linprog(-e, A_eq=sub, b_eq=sb, bounds=bb, method="highs")
            l_ = math.ceil(lo.x[j] - 1e-6)
            h_ = math.inf if hi.status == 3 or hi.x[j] > 1e8 else math.floor(hi.x[j] + 1e-6)
            if kind[unk[c]] != "X" and unsafe(kind[unk[c]], l_, h_):
                fails.append((unk[c], kind[unk[c]], (l_, h_)))
    out(f"independent Stage 2 check: {len(kind)} cells ({len(unk)} unknown: "
        + ", ".join(f"{k} {by_kind[k]}" for k in sorted(by_kind)) + f"), {r} relations, {n_pub} public cell(s); "
        f"{len(fails)} cell(s) fail")
    for c, k, iv in sorted(fails):
        out(f"  FAIL {k} {c}" + (f" {iv}" if show else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release")
    ap.add_argument("--links")
    ap.add_argument("--public")
    ap.add_argument("--layout")
    ap.add_argument("--show-bounds", action="store_true")
    a = ap.parse_args()
    sys.exit(check(a.release, a.links, a.public, a.layout, a.show_bounds))
