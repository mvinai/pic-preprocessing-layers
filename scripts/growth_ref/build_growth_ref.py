"""Build the frozen weight-for-age reference CSV for Stage 3 (growth_weight_for_age_ref_v1.csv).
Sources (downloaded 2026-11-17, files in raw/, SHA-256 in raw/SOURCE_SHA256SUMS.txt):
  cn_li2009: Zong XN, Li H. PLoS ONE 2013;8(3):e59569, Table 2 (JATS XML from journals.plos.org; CC BY 4.0) = the
             LMS parameters of the China 0-18 y growth standard (Li H, Ji CY, Zong XN, Zhang YQ. Zhonghua Er Ke Za Zhi
             2009;47(7):487-492), weight-for-age, exact ages 0-216 months.
  who:       WHO Child Growth Standards 2006, weight-for-age expanded tables by day 0-1856 (cdn.who.int), and the
             WHO Reference 2007 weight-for-age expanded tables by month 61-120 (cdn.who.int)."""
import csv, re, sys, xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
R = Path(__file__).resolve().parent / "raw"
DPM = 365.25 / 12.0     # 30.4375 days per month (WHO convention; the builder's age_months uses the same)


def num(s):
    s = s.strip().replace("\u2212", "-")
    s = re.sub(r"^-+", "-", s)          # source typo '−−0.21' (girls, 10 months; PDF shows the same double sign)
    return float(s)


rows = []
root = ET.parse(R / "zong_li_2013.xml").getroot()
tw = list(root.iter("table-wrap"))[1]
assert "Weight-for-age SD curves" in "".join(tw.find("caption").itertext())
for tr in list(tw.iter("tr"))[2:]:
    c = ["".join(x.itertext()).strip() for x in tr]
    mo = int(c[0])
    for sex, off in (("M", 1), ("F", 9)):
        L, M, S = num(c[off]), num(c[off + 1]), num(c[off + 2])
        med = num(c[off + 5])
        assert abs(med - M) < 1e-9, (mo, sex)
        rows.append(dict(reference="cn_li2009", component="China 2009 (Li 2009; Zong & Li 2013 Table 2)", sex=sex,
                         age_days=round(mo * DPM, 4), age_months=float(mo), L=L, M=M, S=S))
for sex, fn in (("M", "wfa-boys-zscore-expanded-tables.xlsx"), ("F", "wfa-girls-zscore-expanded-tables.xlsx")):
    d = pd.read_excel(R / fn)
    assert list(d.columns[:4]) == ["Day", "L", "M", "S"] and d["Day"].tolist() == list(range(1857))
    for _, r in d.iterrows():
        rows.append(dict(reference="who", component="WHO Child Growth Standards 2006 (expanded, by day)", sex=sex,
                         age_days=float(int(r["Day"])), age_months=round(int(r["Day"]) / DPM, 6),
                         L=float(r["L"]), M=float(r["M"]), S=float(r["S"])))
for sex, fn in (("M", "who2007_boys_z_exp.xlsx"), ("F", "who2007_girls_z_exp.xlsx")):
    d = pd.read_excel(R / fn)
    assert list(d.columns[:4]) == ["Month", "L", "M", "S"] and d["Month"].tolist() == list(range(61, 121))
    for _, r in d.iterrows():
        rows.append(dict(reference="who", component="WHO Reference 2007 (expanded, by month)", sex=sex,
                         age_days=round(int(r["Month"]) * DPM, 4), age_months=float(int(r["Month"])),
                         L=float(r["L"]), M=float(r["M"]), S=float(r["S"])))
out = Path(sys.argv[1])
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["reference", "component", "sex", "age_days", "age_months", "L", "M", "S"])
    for r in rows:
        w.writerow([r["reference"], r["component"], r["sex"], f"{r['age_days']:.4f}", f"{r['age_months']:.6f}",
                    f"{r['L']:.4f}", f"{r['M']:.4f}", f"{r['S']:.5f}"])
print(len(rows), "rows ->", out)
