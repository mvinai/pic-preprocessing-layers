"""Synthetic PIC-like inputs for testing feature_sparsity_frequency_audit.py (Stage 2). NO real data.

Builds a fake project: canonical encounters, harmonization fields, subject split, the surgery complexity map, small
v1-schema maps (channel, drug, specimen, organism, neuro set) and every raw table the script reads or hashes, with
planted boundary cases (release time at landmarks, bin edges, gaps, HR duplicates, SBP/DBP pairs, P/F and S/F draws,
per-row units, magnesium 'ug/dl' rows, date-only starts, ENDDATE variants, osmotic agents by diagnosis, susceptibility
add-on keys, Layer A/C surgeries). Returns the paths and the expected values the run must reproduce."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

FMT = "%Y-%m-%d %H:%M:%S"
SENTINEL_SUBJECT = "987654321"
LANDMARKS = (24, 48, 72, 120, 168, 240, 336)


def ts(x: pd.Timestamp) -> str:
    return x.strftime(FMT)


def sha(p: Path) -> str:
    h = hashlib.sha256()
    h.update(Path(p).read_bytes())
    return h.hexdigest()


# ---------------------------------------------------------------- synthetic v1-schema maps
CH_COLS = ["table", "itemid", "label", "category", "units", "train_encounters", "train_rows", "channel", "table51_system",
           "triggered", "contract16", "unit_conversion", "plausibility_low", "plausibility_high", "source_itemid",
           "merge_rule", "stage2_check", "status", "confidence", "decision", "notes"]
CONV = {
    "1011": "mmol/L -> mg/dL x18.016", "5041": "umol/L -> mg/dL /88.4", "5032": "umol/L -> mg/dL /88.4",
    "5099": "g/L /10 -> g/dL", "5257": "per-row VALUEUOM: g/L /10; g/dL x1 -> g/dL",
    "6525": "per-row VALUEUOM: ng/mL x1; pg/mL /1000 -> ng/mL", "5222": "% -> fraction /100",
    "5034": "mmol/L -> mg/dL x4.008", "5215": "mmol/L -> mg/dL x4.008",
    "5059": "mmol/L as is; 'ug/dl' rows: Stage 2 distribution check - convert if consistent after scaling, else drop those rows",
    "1008": "aggregation chosen after the Stage 2 charting-interval check: daily total or mL/kg/h over the charting interval (imputed weight, A2.12)",
    "1009": "mL/kg/day: daily total of charted intake / weight (charted weight 1014, else imputed from age and gender, A2.12)",
}
# (table, itemid, channel, units, low, high, triggered, status, merge_rule)
CHANNEL_ROWS = [
    ("CHARTEVENTS", "1001", "temperature", "°C", "30", "43", "", "assigned", ""),
    ("CHARTEVENTS", "1002", "heart_rate", "bpm", "30", "300", "", "assigned", "1002 Pulse and 1003 Heart Rate form one channel"),
    ("CHARTEVENTS", "1003", "heart_rate", "bpm", "30", "300", "", "assigned", "1002 Pulse and 1003 Heart Rate form one channel"),
    ("CHARTEVENTS", "1004", "respiratory_rate", "bpm", "5", "120", "", "assigned", ""),
    ("CHARTEVENTS", "1006", "spo2", "%", "30", "100", "", "assigned", ""),
    ("CHARTEVENTS", "1008", "urine_output", "ml", "0", "3000", "", "assigned", ""),
    ("CHARTEVENTS", "1009", "fluid_intake", "ml", "0", "5000", "", "assigned", ""),
    ("CHARTEVENTS", "1010", "", "ml", "", "", "", "dropped", ""),
    ("CHARTEVENTS", "1011", "glucose_bedside", "mmol/L", "0.5", "50", "", "assigned", ""),
    ("CHARTEVENTS", "1014", "weight", "kg", "0.4", "150", "", "assigned", ""),
    ("CHARTEVENTS", "1015", "dbp", "mmHg", "10", "200", "", "assigned", ""),
    ("CHARTEVENTS", "1016", "sbp", "mmHg", "30", "250", "", "assigned", ""),
    ("CHARTEVENTS", "1099", "", "", "", "", "", "not_used", ""),
    ("LABEVENTS", "5041", "creatinine", "μmol/L", "5", "2000", "", "assigned_merge_pending_stage2",
     "merge 5041+5032 into 'creatinine' (same unit (umol/L)); main item 5041; source_itemid kept"),
    ("LABEVENTS", "5032", "creatinine", "μmol/L", "5", "2000", "", "assigned_merge_pending_stage2",
     "merge 5041+5032 into 'creatinine' (same unit (umol/L)); main item 5041; source_itemid kept"),
    ("LABEVENTS", "5222", "fio2", "%", "21", "100", "", "assigned", ""),
    ("LABEVENTS", "5239", "po2", "mmHg", "10", "700", "", "assigned_merge_pending_stage2", "merge 5239+5244 into 'po2'; main item 5239"),
    ("LABEVENTS", "5244", "po2", "mmHg", "10", "700", "", "assigned_merge_pending_stage2", "merge 5239+5244 into 'po2'; main item 5239"),
    ("LABEVENTS", "5243", "pf_ratio_lab_crosscheck", "mmHg", "10", "700", "", "assigned", ""),
    ("LABEVENTS", "5099", "hemoglobin", "g/L", "20", "250", "", "assigned_merge_pending_stage2", "merge 5099+5257 into 'hemoglobin'; main item 5099"),
    ("LABEVENTS", "5257", "hemoglobin", "g/L; g/dL", "20", "250", "", "assigned_merge_pending_stage2", "merge 5099+5257 into 'hemoglobin'; main item 5099"),
    ("LABEVENTS", "6525", "troponin_t", "ng/mL; pg/ml", "", "", "triggered", "assigned", ""),
    ("LABEVENTS", "5059", "magnesium", "mmol/L; ug/dl", "0.1", "5", "", "assigned", ""),
    ("LABEVENTS", "5215", "calcium_total", "mmol/L", "0.5", "5", "", "assigned_merge_pending_stage2", "merge 5215+5034; main item 5215"),
    ("LABEVENTS", "5034", "calcium_total", "mmol/L", "0.5", "5", "", "assigned_merge_pending_stage2", "merge 5215+5034; main item 5215"),
    ("LABEVENTS", "6796", "lipase", "U/L", "", "", "", "assigned_merge_pending_stage2", "6473 merge into 'lipase' only if clearly on the same scale as 6796"),
    ("LABEVENTS", "6473", "lipase", "", "", "", "", "assigned_merge_pending_stage2", "6473 merge into 'lipase' only if clearly on the same scale as 6796"),
    ("LABEVENTS", "5230", "sodium", "mmol/L", "100", "190", "", "assigned", ""),
    ("LABEVENTS", "5227", "lactate", "mmol/L", "0.1", "30", "", "assigned", ""),
    ("LABEVENTS", "6080", "", "", "", "", "", "dropped", ""),
    ("DERIVED", "derived:map", "map", "mmHg", "15", "200", "", "derived", ""),
    ("DERIVED", "derived:pf_ratio", "pf_ratio", "mmHg", "10", "700", "", "derived", ""),
    ("DERIVED", "derived:sf_ratio", "sf_ratio", "ratio", "", "", "", "derived", ""),
]
DERIVED_CONV = {"derived:map": "MAP = (SBP + 2 x DBP) / 3", "derived:pf_ratio": "P/F = pO2 (mmHg) / (FiO2 % / 100)",
                "derived:sf_ratio": "S/F = SpO2 % / (FiO2 % / 100)"}
DRUGS = [  # drug_name, drug_name_en, group, sensitivity_subgroup, conditional_rule
    ("多巴胺", "Dopamine", "vasoactive", "inotrope", ""),
    ("肾上腺素", "Epinephrine", "vasoactive", "vasopressor", ""),
    ("米力农", "Milrinone", "vasoactive", "inotrope", ""),
    ("头孢曲松", "Ceftriaxone", "antibiotic", "", ""),
    ("咪达唑仑", "Midazolam", "sedative_analgesic", "", ""),
    ("20%甘露醇注射液", "20% Mannitol Injection", "osmotic_conditional", "",
     "neuro admitting dx -> osmotherapy; otherwise -> diuretic; sensitivity: all -> osmotherapy"),
    ("10%氯化钠注射液", "10% Sodium Chloride Injection", "osmotic_conditional", "",
     "neuro admitting dx -> osmotherapy; otherwise -> electrolyte_replacement (group D18); sensitivity: all -> osmotherapy"),
    ("葡萄糖注射液", "Glucose Injection", "none", "", ""),
    ("左乙拉西坦", "Levetiracetam", "aed", "", ""),
]
SPECIMENS = [  # spec_itemid, spec_type_desc, group, base, addons
    ("LIS0001", "血培养", "blood", "LIS0001", ""),
    ("LIS0001+LIS0204", "血培养,常规药敏定量试验(MIC)", "blood", "LIS0001", "YES"),
    ("LIS0002", "尿培养", "urine", "LIS0002", ""),
    ("LIS0003", "脑脊液培养", "csf", "LIS0003", ""),
    ("LIS0156", "一般细菌培养鉴定", "culture_site_unknown", "LIS0156", ""),
]
ORGS = [  # org_itemid, org_name, category
    ("", "", "blank (no organism recorded)"),
    ("MIC1001", "金黄色葡萄球菌", "Gram-positive coccus"),
    ("MIC1002", "无细菌生长", "negative/no growth"),
    ("MIC1003", "革兰阳性球菌涂片", "smear finding (not an isolate)"),
    ("MIC1008", "正常菌群生长", "normal flora / mixed"),
]
NEURO = [  # icd10_3char, status, match_rule
    ("S06", "include", "3-character prefix"), ("G00", "include", "3-character prefix"),
    ("A39", "partial", "3-character prefix"), ("Q28", "partial", "4-character prefix Q28.2 or Q28.3 only"),
    ("G40", "exclude", ""), ("G91", "exclude", ""),
]
CMAP_NAMES = [("vsd repair", "Cardiac", "major open", "2"), ("laparotomy", "GI/abdominal", "major open", ""),
              ("craniotomy", "Neurosurgical", "major open", ""), ("central line placement", "Vascular access/line placement",
                                                                   "minor/diagnostic", "")]


def extra_lab_rows(n: int) -> list:
    """Generic extra lab channels for the scale test (no conversion, simple gate)."""
    return [("LABEVENTS", str(8000 + i), f"lab_x{i:02d}", "u", "0", "1000", "", "assigned", "") for i in range(n)]


def write_maps(d: Path, extra_labs: int = 0) -> dict:
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for tb, it, ch, un, lo, hi, trig, st, mr in CHANNEL_ROWS + extra_lab_rows(extra_labs):
        rows.append({"table": tb, "itemid": it, "label": f"label {it}", "category": "", "units": un, "train_encounters": "100",
                     "train_rows": "1000", "channel": ch, "table51_system": "Synthetic", "triggered": trig,
                     "contract16": "likely" if it in ("1003", "5041") else "",
                     "unit_conversion": DERIVED_CONV.get(it, CONV.get(it, "")) if st != "not_used" else "",
                     "plausibility_low": lo, "plausibility_high": hi, "source_itemid": it, "merge_rule": mr, "stage2_check": "",
                     "status": st, "confidence": "high", "decision": "synthetic", "notes": ""})
    pd.DataFrame(rows, columns=CH_COLS).to_csv(d / "feature_channel_map_v1.csv", index=False)
    pd.DataFrame([{"drug_name": a, "drug_name_en": b, "units": "mg", "train_encounters": "50", "group": g,
                   "sensitivity_subgroup": sg, "conditional_rule": cr, "subclass": "", "confidence": "high",
                   "decision": "synthetic", "notes": ""} for a, b, g, sg, cr in DRUGS]).to_csv(d / "drug_group_map_v1.csv", index=False)
    pd.DataFrame([{"spec_itemid": a, "spec_type_desc": b, "train_encounters": "40", "specimen_group": g, "base_code": bc,
                   "test_type": "culture", "susceptibility_addons": "MIC" if ad else "", "result_leak_risk": ad,
                   "confidence": "high", "decision": "synthetic", "notes": ""} for a, b, g, bc, ad in SPECIMENS]
                 ).to_csv(d / "specimen_map_v1.csv", index=False)
    pd.DataFrame([{"org_itemid": a, "org_name": b, "train_encounters": "30", "canonical_name": b or "(blank organism)",
                   "category": c, "contaminant_rule": "", "confidence": "high", "decision": "", "notes": ""} for a, b, c in ORGS]
                 ).to_csv(d / "organism_map_v1.csv", index=False)
    pd.DataFrame([{"icd10_3char": c, "label": f"label {c}", "domain": "x", "status": s, "match_rule": m, "reason": "synthetic"}
                  for c, s, m in NEURO]).to_csv(d / "neuro_admitting_dx_icd10_v1.csv", index=False)
    names = {"channel": "feature_channel_map_v1.csv", "drug": "drug_group_map_v1.csv", "specimen": "specimen_map_v1.csv",
             "organism": "organism_map_v1.csv", "neuro": "neuro_admitting_dx_icd10_v1.csv"}
    return {k: sha(d / v) for k, v in names.items()}, {k: len(pd.read_csv(d / v, dtype=str)) for k, v in names.items()}


def build_cmap() -> pd.DataFrame:
    rows = []
    for nm, bs, it, st in CMAP_NAMES:
        ig = {"minor/diagnostic": "I1", "catheter/interventional": "I1", "minimally invasive": "I2", "major open": "I3"}[it]
        sg = {"1": "G1", "2": "G1", "3": "G2", "4": "G3", "5": "G3"}.get(st, "")
        rows.append({"raw_name": nm, "standardized_name": nm, "body_system": bs, "intensity": it, "I_group": ig, "STAT": st,
                     "STAT_group": sg, "complexity_group": sg if bs == "Cardiac" else ig, "review_flag": "N",
                     "n_components": "1", "STAT_group_crosses": "N" if bs == "Cardiac" else ""})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- project
def make_project(root: Path, n_subjects: int = 420, seed: int = 5, extra_labs: int = 0, vitals_step: float | None = None) -> dict:
    rng = np.random.default_rng(seed)
    root = Path(root)
    raw = root / "raw"
    interim = root / "data" / "interim"
    for p in (raw, interim, root / "config"):
        p.mkdir(parents=True, exist_ok=True)
    subjects = [str(300000 + i) for i in range(n_subjects - 1)] + [SENTINEL_SUBJECT]
    split_of = {s: ("train" if (u := rng.random()) < 0.6 else ("val" if u < 0.8 else "test")) for s in subjects}
    split_of[SENTINEL_SUBJECT] = "train"
    base = pd.Timestamp("2015-03-01 00:00:00")
    encs, stays, adms, pats = [], [], [], []
    hadm, icu = 700000, 900000
    for si, s in enumerate(subjects):
        dob = base - pd.Timedelta(days=int(rng.choice([5, 20, 200, 500, 1500, 3000, 5000, 6200])))
        pats.append({"ROW_ID": si, "SUBJECT_ID": s, "GENDER": ("" if si == 3 else ("M" if rng.random() < 0.55 else "F")),
                     "DOB": ts(dob), "DOD": ""})
        n_enc = 1 if rng.random() < 0.85 else 2
        t = base + pd.Timedelta(days=int(rng.integers(0, 900)), hours=int(rng.integers(0, 24)), minutes=int(rng.integers(0, 60)))
        for e in range(n_enc):
            hadm += 1
            hl = [str(hadm)]
            admit = t - pd.Timedelta(hours=float(rng.choice([2.0, 30.0, 100.0, 200.0])))
            u = rng.random()
            dx = rng.choice(["S06.001", "Q28.200", "Q28.000", "J18.900", "G40.100", "A39.000", "K35.800", "", "XYZ", "J18.9,G40.1"],
                            p=[0.1, 0.05, 0.05, 0.3, 0.1, 0.05, 0.15, 0.08, 0.04, 0.08])
            ed = "" if u < 0.3 else (ts(admit) if u < 0.5 else ts(admit - pd.Timedelta(hours=1)))
            adms.append({"ROW_ID": len(adms), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ADMITTIME": ts(admit),
                         "DISCHTIME": ts(t + pd.Timedelta(days=30)), "EDREGTIME": ed, "ICD10_CODE_CN": dx})
            t0 = t
            span = float(rng.choice([rng.uniform(4, 30), rng.uniform(30, 200), rng.uniform(200, 420)]))
            out = (t0 + pd.Timedelta(hours=span)).floor("s")
            unit = rng.choice(["PICU", "NICU", "SICU", "CICU"])
            gap = None
            if span > 60 and rng.random() < 0.12:  # two stays with a gap (merged encounter)
                a = (t0 + pd.Timedelta(hours=20)).floor("s")
                b = (a + pd.Timedelta(hours=float(rng.uniform(13, 23)))).floor("s")
                stays.append({"ROW_ID": len(stays), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": str(icu), "FIRST_CAREUNIT": unit,
                              "INTIME": ts(t0), "OUTTIME": ts(a)}); icu += 1
                stays.append({"ROW_ID": len(stays), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": str(icu), "FIRST_CAREUNIT": "PICU",
                              "INTIME": ts(b), "OUTTIME": ts(out)}); icu += 1
                gap = ((a - t0).total_seconds() / 3600, (b - t0).total_seconds() / 3600)
            else:
                stays.append({"ROW_ID": len(stays), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": str(icu), "FIRST_CAREUNIT": unit,
                              "INTIME": ts(t0), "OUTTIME": ts(out)}); icu += 1
            encs.append({"encounter_id": f"E{s}_{e + 1}", "SUBJECT_ID": s, "HADM_ID": str(hadm), "hadm_list": hl, "t0": t0,
                         "OUT": out, "split": {"val": "validation"}.get(split_of[s], split_of[s]), "gap": gap,
                         "died": int(rng.random() < 0.08), "dx": dx, "unit": unit})
            t = (out + pd.Timedelta(hours=float(rng.uniform(40, 400)))).floor("s")
    enc = pd.DataFrame(encs)
    canon = pd.DataFrame({"encounter_id": enc["encounter_id"], "SUBJECT_ID": enc["SUBJECT_ID"],
                          "HADM_ID": enc["HADM_ID"].map(lambda h: f"{h}.0"), "hadm_id_list": ["[" + ", ".join(h) + "]" for h in enc["hadm_list"]],
                          "INTIME": enc["t0"].map(ts), "OUTTIME": enc["OUT"].map(ts), "Died": enc["died"]})
    canon.to_csv(interim / "canonical_icu_encounters.csv", index=False)
    pd.DataFrame({"encounter_id": enc["encounter_id"], "SUBJECT_ID": enc["SUBJECT_ID"], "t0": enc["t0"].map(ts),
                  "INTIME": enc["t0"].map(ts), "OUTTIME": enc["OUT"].map(ts),
                  "span_icu_hours": ((enc["OUT"] - enc["t0"]).dt.total_seconds() / 3600).round(4),
                  "gap_count_between_stays": enc["gap"].map(lambda g: 1 if g else 0)}
                 ).to_csv(interim / "encounter_harmonization_fields.csv", index=False)
    pd.DataFrame({"SUBJECT_ID": subjects, "split": [split_of[s] for s in subjects]}).to_csv(root / "subject_splits.csv", index=False)
    build_cmap().to_csv(root / "config" / "surgery_complexity_map_v1.csv", index=False)
    map_sha, map_rows = write_maps(interim / "maps", extra_labs)
    extra_p = np.linspace(0.005, 0.9, extra_labs) if extra_labs else np.array([])
    pd.DataFrame(pats).to_csv(raw / "PATIENTS.csv", index=False)
    pd.DataFrame(adms).to_csv(raw / "ADMISSIONS.csv", index=False)
    pd.DataFrame(stays).to_csv(raw / "ICUSTAYS.csv", index=False)
    pd.DataFrame([{"ROW_ID": 1, "ITEMID": "1003", "LABEL": "Heart Rate"}]).to_csv(raw / "D_ITEMS.csv", index=False)
    pd.DataFrame([{"ROW_ID": 1, "ITEMID": "5041", "LABEL": "Creatinine"}]).to_csv(raw / "D_LABITEMS.csv", index=False)
    pd.DataFrame([{"SUBJECT_ID": subjects[0], "HADM_ID": "1", "ITEMID": "1", "CHARTTIME": "", "VALUE": "1"}]).to_csv(raw / "OUTPUTEVENTS.csv", index=False)
    pd.DataFrame([{"SUBJECT_ID": subjects[0], "HADM_ID": "1", "CHARTTIME": "", "AMOUNT": "1"}]).to_csv(raw / "INPUTEVENTS.csv", index=False)

    chart, lab, rx, micro, surg = [], [], [], [], []
    H = lambda t0, h: (t0 + pd.Timedelta(hours=float(h))).floor("s")
    for r in enc.itertuples():
        sp = (r.OUT - r.t0).total_seconds() / 3600
        h = r.HADM_ID
        # vitals: hourly-ish HR 1003, some 1002 at the same time and alone; RR; temp; SpO2; SBP/DBP
        for k in np.arange(-3, min(sp, 340), vitals_step or float(rng.choice([1.0, 2.0, 4.0]))):
            chart.append([r.SUBJECT_ID, h, "1003", ts(H(r.t0, k)), f"{rng.uniform(80, 160):.0f}", "bpm"])
            if rng.random() < 0.2:
                chart.append([r.SUBJECT_ID, h, "1002", ts(H(r.t0, k)), f"{rng.uniform(80, 160):.0f}", "bpm"])
            if rng.random() < 0.6:
                chart.append([r.SUBJECT_ID, h, "1004", ts(H(r.t0, k)), f"{rng.uniform(15, 50):.0f}", "bpm"])
            if rng.random() < 0.5:
                sv = rng.uniform(60, 120)
                dv = sv + 5 if rng.random() < 0.05 else rng.uniform(30, sv - 5)
                chart.append([r.SUBJECT_ID, h, "1016", ts(H(r.t0, k)), f"{sv:.0f}", "mmHg"])
                chart.append([r.SUBJECT_ID, h, "1015", ts(H(r.t0, k)), f"{dv:.0f}", "mmHg"])
            if rng.random() < 0.4:
                chart.append([r.SUBJECT_ID, h, "1006", ts(H(r.t0, k)), f"{rng.uniform(85, 100):.0f}", "%"])
        if rng.random() < 0.3:
            chart.append([r.SUBJECT_ID, h, "1002", ts(H(r.t0, 0.5)), "120", "bpm"])
        for k in np.arange(0, min(sp, 340), 4.0):
            chart.append([r.SUBJECT_ID, h, "1001", ts(H(r.t0, k)), f"{rng.uniform(35, 40):.1f}" if rng.random() > 0.02 else "55", "°C"])
        if rng.random() < 0.6:
            chart.append([r.SUBJECT_ID, h, "1014", ts(H(r.t0, 1)), f"{rng.uniform(3, 40):.1f}", "kg"])
        for k in np.arange(1, min(sp, 340), float(rng.choice([2.0, 6.0, 24.0]))):
            chart.append([r.SUBJECT_ID, h, "1008", ts(H(r.t0, k)), f"{rng.uniform(5, 200):.0f}", "ml"])
            chart.append([r.SUBJECT_ID, h, "1009", ts(H(r.t0, k + 0.5)), f"{rng.uniform(5, 300):.0f}", "ml"])
        # labs: creatinine (two items), sodium, lactate, hemoglobin (per-row units), troponin, magnesium, calcium, lipase
        for k in rng.uniform(-8, min(sp, 340), size=int(rng.integers(1, 8))):
            lab.append([r.SUBJECT_ID, h, rng.choice(["5041", "5041", "5032"]), ts(H(r.t0, k)), f"{rng.uniform(20, 120):.0f}", "μmol/L"])
            lab.append([r.SUBJECT_ID, h, "5230", ts(H(r.t0, k)), f"{rng.uniform(125, 150):.0f}", "mmol/L"])
        if rng.random() < 0.5:
            k = rng.uniform(0, min(sp, 100))
            lab.append([r.SUBJECT_ID, h, "5099", ts(H(r.t0, k)), f"{rng.uniform(80, 140):.0f}", "g/L"])
            lab.append([r.SUBJECT_ID, h, "5257", ts(H(r.t0, k)), f"{rng.uniform(8, 14):.1f}", "g/dL"])
            lab.append([r.SUBJECT_ID, h, "5257", ts(H(r.t0, k + 1)), f"{rng.uniform(80, 140):.0f}", "g/L"])
        draws = rng.uniform(-8, min(sp, 340), size=4)
        for i, p_ in enumerate(extra_p):
            if rng.random() < p_:
                for k in rng.choice(draws, size=int(rng.integers(1, 4))):
                    lab.append([r.SUBJECT_ID, h, str(8000 + i), ts(H(r.t0, k)), f"{rng.uniform(1, 100):.1f}", "u"])
        if rng.random() < 0.15:
            lab.append([r.SUBJECT_ID, h, "6525", ts(H(r.t0, rng.uniform(0, 30))), f"{rng.uniform(5, 90):.0f}", "pg/ml"])
            lab.append([r.SUBJECT_ID, h, "6525", ts(H(r.t0, rng.uniform(0, 30))), f"{rng.uniform(0.01, 0.2):.3f}", "ng/mL"])
        if rng.random() < 0.4:
            lab.append([r.SUBJECT_ID, h, "5059", ts(H(r.t0, rng.uniform(0, 48))), f"{rng.uniform(0.6, 1.2):.2f}", "mmol/L"])
            if rng.random() < 0.3:
                lab.append([r.SUBJECT_ID, h, "5059", ts(H(r.t0, rng.uniform(0, 48))), f"{rng.uniform(1.5, 2.8):.2f}", "ug/dl"])
        if rng.random() < 0.4:
            k = rng.uniform(0, 48)
            lab.append([r.SUBJECT_ID, h, "5215", ts(H(r.t0, k)), f"{rng.uniform(1.0, 1.4):.2f}", "mmol/L"])
            lab.append([r.SUBJECT_ID, h, "5034", ts(H(r.t0, k)), f"{rng.uniform(2.1, 2.6):.2f}", "mmol/L"])
        if rng.random() < 0.2:
            lab.append([r.SUBJECT_ID, h, "6796", ts(H(r.t0, rng.uniform(0, 48))), f"{rng.uniform(10, 60):.0f}", "U/L"])
            lab.append([r.SUBJECT_ID, h, "6473", ts(H(r.t0, rng.uniform(0, 48))), "阴性", ""])
        if rng.random() < 0.3:
            lab.append([r.SUBJECT_ID, h, "5227", ts(H(r.t0, rng.uniform(-5, 30))), f"{rng.uniform(0.5, 6):.1f}", "mmol/L"])
        # blood gases: FiO2 + pO2 (same draw), FiO2 alone with SpO2 nearby (S/F), analyser P/F
        for k in rng.uniform(0, min(sp, 120), size=int(rng.integers(0, 4))):
            tt = H(r.t0, k)
            fi = rng.choice([21, 30, 40, 60, 100])
            lab.append([r.SUBJECT_ID, h, "5222", ts(tt), str(fi), "%"])
            if rng.random() < 0.6:
                po = rng.uniform(40, 200)
                lab.append([r.SUBJECT_ID, h, rng.choice(["5239", "5244"]), ts(tt), f"{po:.0f}", "mmHg"])
                lab.append([r.SUBJECT_ID, h, "5243", ts(tt), f"{po / (fi / 100) + rng.normal(0, 5):.0f}", "mmHg"])
            else:
                chart.append([r.SUBJECT_ID, h, "1006", ts(tt + pd.Timedelta(minutes=int(rng.choice([2, 4, 9])))),
                              str(int(rng.choice([88, 92, 96, 99]))), "%"])
        # prescriptions: timed and date-only starts, ENDDATE variants
        for dn, de, g, sg, cr in DRUGS:
            if rng.random() < (0.35 if g != "osmotic_conditional" else 0.15):
                k = rng.uniform(-30, min(sp, 200))
                n_docs = int(rng.integers(1, 5))
                for j in range(n_docs):
                    tt = H(r.t0, k + j * float(rng.choice([6.0, 10.0, 30.0])))
                    u = rng.random()
                    start = tt.strftime("%Y-%m-%d") if u < 0.1 else ts(tt)
                    v = rng.random()
                    end = "" if v < 0.75 else (ts(tt + pd.Timedelta(hours=float(rng.uniform(2, 40)))) if v < 0.95 else ts(tt - pd.Timedelta(hours=1)))
                    rx.append([r.SUBJECT_ID, h, start, end, dn, de, "mg"])
        # microbiology: collections with several organism rows; add-on key at the same time as the base key
        for k in rng.uniform(-5, min(sp, 150), size=int(rng.integers(0, 3))):
            tt = ts(H(r.t0, k))
            sid, sdesc, _, _, _ = SPECIMENS[int(rng.integers(0, len(SPECIMENS)))]
            oid, oname, _ = ORGS[int(rng.integers(0, len(ORGS)))]
            micro.append([r.SUBJECT_ID, h, tt, sid, sdesc, oid, oname, "OXACILLIN", "S"])
            if sid == "LIS0001" and oid == "MIC1001":
                micro.append([r.SUBJECT_ID, h, tt, "LIS0001+LIS0204", "血培养,常规药敏定量试验(MIC)", oid, oname, "VANCOMYCIN", "S"])
        # surgery: recent (Layer A), during ICU (Layer C), remote
        u = rng.random()
        if u < 0.25:
            nm = CMAP_NAMES[int(rng.integers(0, len(CMAP_NAMES)))][0]
            if u < 0.1:
                b_, e_ = H(r.t0, -4), H(r.t0, -1)
            elif u < 0.2 and sp > 12:
                b_, e_ = H(r.t0, min(sp - 1, 10)), H(r.t0, min(sp - 0.5, 11))
            else:
                b_, e_ = H(r.t0, -60), H(r.t0, -58)
            surg.append([r.SUBJECT_ID, h, h, str(len(surg)), ts(b_ - pd.Timedelta(minutes=20)), ts(e_ + pd.Timedelta(minutes=10)),
                         ts(b_), ts(e_), nm, "general", "supine"])
    surg.append([enc["SUBJECT_ID"].iloc[0], enc["HADM_ID"].iloc[0], enc["HADM_ID"].iloc[0], "x1", "", "",
                 ts(enc["t0"].iloc[0] - pd.Timedelta(hours=3)), ts(enc["t0"].iloc[0] - pd.Timedelta(hours=2)), "laparotomy", "general", "supine"])
    # planted landmark / bin boundary labs on train long stays (sodium): released at exactly t0 + 24 h
    planted = enc[(enc["split"] == "train") & ((enc["OUT"] - enc["t0"]).dt.total_seconds() > 30 * 3600)].head(12)
    for r in planted.itertuples():
        lab.append([r.SUBJECT_ID, r.HADM_ID, "5227", ts(r.t0 + pd.Timedelta(hours=23)), "2.0", "mmol/L"])
    # noise rows
    e0 = enc.iloc[0]
    for _ in range(12):
        chart.append(["111", "999999", "1003", ts(e0.t0), "80", "bpm"])
        chart.append([e0.SUBJECT_ID, e0.HADM_ID, "1003", "", "80", "bpm"])
        lab.append([e0.SUBJECT_ID, "123", "5041", ts(e0.t0), "1", "μmol/L"])
        rx.append([e0.SUBJECT_ID, e0.HADM_ID, "", "", "头孢曲松", "Ceftriaxone", "g"])
        micro.append(["111", "1", ts(e0.t0), "LIS0001", "血培养", "", "", "", ""])
    ch = pd.DataFrame(chart, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"])
    ch["VALUENUM"] = pd.to_numeric(ch["VALUE"], errors="coerce").map(lambda v: "" if pd.isna(v) else f"{v:g}")
    ch["ICUSTAY_ID"] = ""
    ch.to_csv(raw / "CHARTEVENTS.csv", index=False)
    lb = pd.DataFrame(lab, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"])
    lb["VALUENUM"] = pd.to_numeric(lb["VALUE"], errors="coerce").map(lambda v: "" if pd.isna(v) else f"{v:g}")
    lb.to_csv(raw / "LABEVENTS.csv", index=False)
    pd.DataFrame(rx, columns=["SUBJECT_ID", "HADM_ID", "STARTDATE", "ENDDATE", "DRUG_NAME", "DRUG_NAME_EN", "DOSE_UNIT_RX"]).to_csv(
        raw / "PRESCRIPTIONS.csv", index=False)
    pd.DataFrame(micro, columns=["SUBJECT_ID", "HADM_ID", "CHARTTIME", "SPEC_ITEMID", "SPEC_TYPE_DESC", "ORG_ITEMID", "ORG_NAME",
                                 "AB_NAME", "INTERPRETATION"]).to_csv(raw / "MICROBIOLOGYEVENTS.csv", index=False)
    pd.DataFrame(surg, columns=["SUBJECT_ID", "HADM_ID", "VISIT_ID", "ROW_ID", "ANES_START_TIME", "ANES_END_TIME",
                                "SURGERY_BEGIN_TIME", "SURGERY_END_TIME", "SURGERY_NAME", "ANES_METHOD", "SURGERY_POSITION"]
                 ).assign(OPER_ID=lambda d: d["ROW_ID"])[
        ["ROW_ID", "SUBJECT_ID", "HADM_ID", "VISIT_ID", "OPER_ID", "ANES_START_TIME", "ANES_END_TIME", "SURGERY_BEGIN_TIME",
         "SURGERY_END_TIME", "SURGERY_NAME", "ANES_METHOD", "SURGERY_POSITION"]].to_csv(raw / "SURGERY_INFO.csv", index=False)
    split_counts = {s: int((enc["split"] == s).sum()) for s in ("train", "validation", "test")}
    elig = {L: int((enc["OUT"] > enc["t0"] + pd.Timedelta(hours=L)).sum()) for L in LANDMARKS}
    raw_rows = {"CHARTEVENTS": len(chart), "LABEVENTS": len(lab), "PRESCRIPTIONS": len(rx), "MICROBIOLOGYEVENTS": len(micro),
                "PATIENTS": len(pats), "ADMISSIONS": len(adms), "ICUSTAYS": len(stays), "SURGERY_INFO": len(surg)}
    raw_sha = {n: sha(raw / f"{n}.csv") for n in ("D_ITEMS", "D_LABITEMS", "CHARTEVENTS", "LABEVENTS", "OUTPUTEVENTS",
                                                 "INPUTEVENTS", "PRESCRIPTIONS", "MICROBIOLOGYEVENTS")}
    pinned = {"canonical": sha(interim / "canonical_icu_encounters.csv"),
              "harmonization": sha(interim / "encounter_harmonization_fields.csv"),
              "split": sha(root / "subject_splits.csv"), "complexity_map": sha(root / "config" / "surgery_complexity_map_v1.csv")}
    return {"root": root, "raw": raw, "interim": interim, "enc": enc, "split_counts": split_counts, "eligible": elig,
            "raw_rows": raw_rows, "raw_sha": raw_sha, "pinned": pinned, "map_sha": map_sha, "map_rows": map_rows,
            "expected_encounters": len(enc), "expected_subjects": enc["SUBJECT_ID"].nunique(),
            "died": int(enc["died"].sum()), "surgery_rows": len(surg), "planted": planted}


def make_config(M, proj: dict, out_dir: Path, v17_path: Path, **kw):
    root = proj["root"]
    return M.S2Config(
        raw_dir=proj["raw"], canonical_path=proj["interim"] / "canonical_icu_encounters.csv",
        harmonization_path=proj["interim"] / "encounter_harmonization_fields.csv", split_path=root / "subject_splits.csv",
        complexity_map_path=root / "config" / "surgery_complexity_map_v1.csv", maps_dir=proj["interim"] / "maps",
        out_dir=out_dir, v17_module_path=v17_path, pinned_sha256=dict(proj["pinned"]), map_sha256=dict(proj["map_sha"]),
        map_rows=dict(proj["map_rows"]), raw_sha256=dict(proj["raw_sha"]), harm_t0_col="t0", harm_out_col="OUTTIME",
        expected_encounters=proj["expected_encounters"], expected_subjects=proj["expected_subjects"],
        expected_split_encounters=dict(proj["split_counts"]), expected_eligible=dict(proj["eligible"]),
        expected_died=proj["died"], expected_raw_rows=dict(proj["raw_rows"]), expected_surgery_rows=proj["surgery_rows"],
        expected_complexity_map_rows=len(CMAP_NAMES), surgery_expected_missing={
            "ROW_ID": 0, "SUBJECT_ID": 0, "HADM_ID": 0, "VISIT_ID": 0, "OPER_ID": 0, "ANES_START_TIME": 1, "ANES_END_TIME": 1,
            "SURGERY_BEGIN_TIME": 0, "SURGERY_END_TIME": 0, "SURGERY_NAME": 0, "ANES_METHOD": 0, "SURGERY_POSITION": 0},
        verbose=False, **kw)
