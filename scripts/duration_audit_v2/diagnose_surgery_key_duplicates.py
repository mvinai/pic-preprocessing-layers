"""Read-only, aggregate-only diagnostic for duplicate surgery keys in SURGERY_INFO.csv.

Prints counts only (no IDs, names, or times). Any count 1-9 prints as "<10".
Usage:  python diagnose_surgery_key_duplicates.py --surgery-info <path to SURGERY_INFO.csv>
"""
import argparse, re, sys
import pandas as pd

def norm_id(x):
    s = str(x).strip().strip("'\"").strip()
    if s.lower() in {"", "nan", "none", "null", "nat"}:
        return ""
    m = re.fullmatch(r"([+-]?\d+)\.0*", s)
    return m.group(1) if m else s

def c(n):
    n = int(n)
    return "<10" if 0 < n < 10 else str(n)

def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--surgery-info", required=True)
    a = ap.parse_args()
    df = pd.read_csv(a.surgery_info, dtype=str, keep_default_na=False, na_filter=False, encoding="utf-8-sig")
    key = ["SUBJECT_ID", "HADM_ID", "OPER_ID"]
    raw = df[key].copy()
    nrm = df[key].apply(lambda col: col.map(norm_id))
    print("rows:", c(len(df)))
    for k in key:
        print(f"blank {k} (after normalization):", c((nrm[k] == "").sum()))
        print(f"{k} values changed by normalization:", c((nrm[k] != raw[k].str.strip()).sum()))
    dup_raw = raw.duplicated(keep=False)
    dup_nrm = nrm.duplicated(keep=False)
    print("rows in a duplicated key group (raw strings):", c(dup_raw.sum()))
    print("rows in a duplicated key group (normalized):", c(dup_nrm.sum()))
    print("duplicated key groups (normalized):", c(nrm[dup_nrm].drop_duplicates().shape[0]))
    print("rows that would be dropped by keep-first:", c(nrm.duplicated(keep="first").sum()))
    blank_key = (nrm == "").any(axis=1)
    print("duplicated rows that have a blank key part:", c((dup_nrm & blank_key).sum()))
    sizes = nrm[dup_nrm].groupby(key).size()
    print("duplicate group sizes:", {f"size_{s}": c(n) for s, n in sizes.value_counts().sort_index().items()})
    other = [col for col in df.columns if col not in key and col != "ROW_ID"]
    d = df[dup_nrm].copy()
    d[key] = nrm[dup_nrm]
    g = d.groupby(key)
    exact = g[other].nunique().max(axis=1).eq(1)
    print("duplicate groups that are exact copies (all non-ROW_ID columns identical):", c(exact.sum()))
    print("duplicate groups that differ in at least one column:", c((~exact).sum()))
    nun = g[other].nunique()
    print("columns that differ within duplicate groups (number of groups):")
    for col in other:
        print(f"  {col}:", c((nun[col] > 1).sum()))
    if "VISIT_ID" in df.columns:
        k2 = nrm.assign(VISIT_ID=df["VISIT_ID"].map(norm_id))
        print("rows duplicated on (SUBJECT_ID, HADM_ID, VISIT_ID, OPER_ID):", c(k2.duplicated(keep=False).sum()))

if __name__ == "__main__":
    main()
