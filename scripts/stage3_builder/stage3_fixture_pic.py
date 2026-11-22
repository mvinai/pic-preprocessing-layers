"""PIC-format fixture generator for the Stage 3 tests (SYNTHETIC DATA ONLY; spec §13, Stage 3 Q3).

Writes a small project (raw PIC-schema tables, canonical/harmonization/split files) with planted encounters whose
expected Layer A/B/C values are known exactly. The decided maps (channel map v2.3, drug, specimen, organism, neuro,
growth reference v1, surgery complexity map v1, frozen φ1 medians) are copied unchanged from maps/ (decision files, no patient
rows), so their real pins hold. Only the data-file pins and the reproduction counts are fixture values.
The team synthetic dataset (synthetic_svd_dataset) is not used or changed."""
from __future__ import annotations

import csv
import json
import hashlib
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MAPS = HERE / "maps"
MAP_NAMES = ["feature_channel_map_v2_3.csv", "drug_group_map_v1.csv", "specimen_map_v1.csv", "organism_map_v1.csv",
             "neuro_admitting_dx_icd10_v1.csv", "growth_weight_for_age_ref_v1.csv"]
BASE = pd.Timestamp("2016-01-01 08:00:00")
DX_CN = {"J18.900": "J18.9", "S06.001": "S06.0", "G40.100": "G40.1", "Q21.100": "Q21.1", "K35.800": "K35.8",
         "A39.000": "A39.0", "P07.300": "P07.3", "C91.000": "C91.0"}
SURG_NAMES = ["asd closure", "angiocardiography", "anoplasty", "acromphalus repair",
              "adenoidal hypertrophy melt+tonsillectomy", "abdominal drainge"]
DRUGS = {"vasoactive": ("盐酸多巴胺注射液", "Dopamine Hydrochloride Injection"),
         "antibiotic": ("利奈唑胺注射液", "Linezolid Injection"),
         "sedative_analgesic": ("丙泊酚注射液", "Propofol Injectable Emulsion"),
         "diuretic": ("呋塞米注射液", "Furosemide Injection"),
         "osmotic": ("20％甘露醇注射液", "20% Mannitol Injection"),
         "none": ("(4:1)葡萄糖氯化钠注射液", "(4:1)Glucose and Sodium Chloride Injection")}
QUAL_TOKENS = ["阴性", "阳性", "阴性(-)", "+", "-", "weak?"]
# Stage 3 review decisions (Nov 17, 2026): planted rows (deterministic; the random stream is unchanged)
UNPARSED_PLANTED = [("1:100", 12), ("床号12", 12), ("12345", 12), ("rare", 3),
                    ("seeattachedreportforfulldetailsplease", 11)]     # stool occult blood 5326, encounter P0
# v1.0.4 (Vinai, Nov 17, 2026): per-channel token additions, planted on encounter P2 (item, hours from t0, text)
QUAL_V104_PLANTED = [(5326, 2.0, "positive+"),      # stool_occult_blood -> positive (Layer B window)
                     (5368, 3.0, "Positive ++"),    # urine_nitrite -> positive (normal form: spaces removed, lower case)
                     (5359, 40.0, "positive++"),    # urine_leukocyte_esterase -> positive (Layer C only)
                     (5377, 4.0, "positive+"),      # csf_protein_pandy -> positive
                     (6809, 5.0, "<3.00"),          # resp_panel_mycoplasma -> negative, below cutoff (Layer B flag)
                     (6809, 30.0, "positive+")]     # not listed for mycoplasma -> stays unparsed
DX_PLANTED = {3: "M98000", 6: "M9800/3", 7: "M80000"}                 # ICD-O x2; dotless ICD-10 M80.0 in D_ICD
D_ICD_EXTRA = {"M80000": "M80.0"}
# v1.0.3: ICD-O morphology codes in the exact forms the real v1.0.2 audit listed as unmapped (bare 6-character codes),
# the clinical-edition form with its behaviour suffix (12 x M98000/3: under the v1.0.2 rule its first code M98000 would
# be listed in /diagnosis/unmapped_codes, n >= 10) and other written variants; first encounter of subjects 200-230.
DX_ICDO_REAL = ["M98000", "M95000", "M80001", "M90800", "M89700", "M95900", "M91200", "M89600", "M91700", "M89000",
                "M98010"]
DX_ICDO_SUFFIX = ["M98000/3"] * 12
DX_ICDO_VARIANTS = ["m9800/3", " M 98000 / 3 ", "M98000/3,C92.0", "\uff2d\uff19\uff15\uff10\uff10\uff10\uff0f\uff13",
                    "M9800.0", "M980003", "M91200/0", "M80001/1"]
DX_ICDO_PLANTED = {200 + k: c for k, c in enumerate(DX_ICDO_REAL + DX_ICDO_SUFFIX + DX_ICDO_VARIANTS)}


def ts(x: pd.Timestamp) -> str:
    return pd.Timestamp(x).strftime("%Y-%m-%d %H:%M:%S")


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _channel_items() -> pd.DataFrame:
    m = pd.read_csv(MAPS / "feature_channel_map_v2_3.csv", dtype=str, keep_default_na=False, encoding="utf-8-sig")
    use = m[m["status"].isin(["assigned", "assigned_merged", "assigned_merged_conditional",
                              "assigned_merged_indicator_only"]) & (m["channel"] != "") & (m["table51_system"] != "")]
    return use[use["table"].isin(["CHARTEVENTS", "LABEVENTS"])].reset_index(drop=True)


def _value(rng, r) -> tuple[str, str]:
    """A VALUE / VALUEUOM pair inside the item's gate (or a qualitative token)."""
    it, units = r["itemid"], r["units"]
    lo = float(r["plausibility_low"]) if r["plausibility_low"] else 0.5
    hi = float(r["plausibility_high"]) if r["plausibility_high"] else 50.0
    u = units.split(";")[0].strip()
    if r["channel"] in ("stool_occult_blood", "urine_nitrite", "urine_leukocyte_esterase", "csf_protein_pandy") \
            or r["channel"].startswith(("resp_panel_", "systemic_viral_")) or r["channel"] == "bacteremia_pcr":
        return str(rng.choice(QUAL_TOKENS)), u
    if it == "5257" and rng.random() < 0.4:
        return f"{rng.uniform(8, 16):.1f}", "g/dL"
    if it == "6525" and rng.random() < 0.4:
        return f"{rng.uniform(5, 200):.0f}", "pg/ml"
    if it == "5059" and rng.random() < 0.03:
        return f"{rng.uniform(1, 3):.2f}", "ug/dl"
    span = hi - lo
    v = rng.uniform(lo + 0.1 * span, lo + 0.6 * span)
    if rng.random() < 0.02 and r["plausibility_high"]:
        v = hi * 3 + 1          # out of gate
    return f"{v:.3g}", u


def make_project(root: Path, n_subjects: int = 300, seed: int = 20261116, padding_extremes: bool = False) -> dict:
    rng = np.random.default_rng(seed)
    root = Path(root)
    raw, interim = root / "raw", root / "data" / "interim"
    for p in (raw, interim / "maps", root / "config"):
        p.mkdir(parents=True, exist_ok=True)
    for fn in MAP_NAMES:
        shutil.copy(MAPS / fn, interim / "maps" / fn)
    shutil.copy(MAPS / "surgery_complexity_map_v1.csv", root / "config" / "surgery_complexity_map_v1.csv")
    shutil.copy(MAPS / "phi1_band_support_frozen_medians_MAP_v1.csv", interim / "phi1_band_support_frozen_medians_MAP_v1.csv")
    items = _channel_items()
    chart_items = items[items["table"] == "CHARTEVENTS"]
    lab_items = items[items["table"] == "LABEVENTS"]
    pats, adms, stays, encs, chart, lab, rx, micro, surg = [], [], [], [], [], [], [], [], []
    hadm, icu = 500000, 800000
    subjects = [str(100001 + i) for i in range(n_subjects)]
    split_of = {}
    planted = {}
    for si, s in enumerate(subjects):
        split_of[s] = "train" if si < 8 else ("train" if (u := rng.random()) < 0.65 else ("val" if u < 0.82 else "test"))
        age_days = int(rng.choice([10, 200, 700, 1500, 3000, 5000]))
        dob = BASE - pd.Timedelta(days=age_days)
        dob_s = ts(dob) if si != 5 else ""            # one missing DOB
        pats.append({"ROW_ID": si, "SUBJECT_ID": s, "GENDER": "M" if rng.random() < 0.55 else "F", "DOB": dob_s,
                     "DOD": "", "EXPIRE_FLAG": "0"})
        n_enc = 1 if (si < 8 or rng.random() < 0.85) else 2
        t = BASE + pd.Timedelta(days=int(rng.integers(0, 600)), hours=int(rng.integers(0, 24)))
        if si < 8:
            t = BASE + pd.Timedelta(days=10 * si)
        for e in range(n_enc):
            hadm += 1
            t0 = t.floor("s")
            planted_enc = si < 8 and e == 0
            if planted_enc:
                span = 400.0 if si in (0, 1, 2, 3) else 90.0
            else:
                span = float(rng.choice([rng.uniform(10, 60), rng.uniform(60, 200), rng.uniform(340, 420)],
                                        p=[0.35, 0.35, 0.3]))
            out = (t0 + pd.Timedelta(hours=span)).floor("s")
            unit = "NICU" if si == 4 else ("PICU" if si < 8 else str(rng.choice(["PICU", "NICU", "SICU", "CICU"])))
            admit = t0 - pd.Timedelta(hours=float(rng.choice([2.0, 30.0, 100.0])))
            if si == 0:
                admit = t0 + pd.Timedelta(minutes=30)       # negative interval (Q4: clipped at 0)
            dx = str(rng.choice(list(DX_CN) + ["", "XYZ", "Z99.999", "J18.900,G40.100"]))
            if si == 0:
                dx = "XYZ"
            if si == 1:
                dx = "Z99.999"
            if si == 2:
                dx = "G40.100"
            dx = DX_PLANTED.get(si, dx) if e == 0 else dx
            dx = DX_ICDO_PLANTED.get(si, dx) if e == 0 else dx            # v1.0.3 (random stream unchanged)
            ed = "" if rng.random() < 0.3 else (ts(admit) if rng.random() < 0.5 else ts(admit - pd.Timedelta(hours=1)))
            adms.append({"ROW_ID": len(adms), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ADMITTIME": ts(admit),
                         "DISCHTIME": ts(out + pd.Timedelta(days=3)), "DEATHTIME": "", "EDREGTIME": ed,
                         "ICD10_CODE_CN": dx, "DIAGNOSIS": "dx"})
            stays.append({"ROW_ID": len(stays), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": str(icu),
                          "FIRST_CAREUNIT": unit, "LAST_CAREUNIT": unit, "INTIME": ts(t0), "OUTTIME": ts(out),
                          "LOS": f"{span / 24:.2f}"})
            icu += 1
            eid = f"E{s}_{e + 1}"
            died = 1 if (not planted_enc and rng.random() < 0.1) else 0
            encs.append({"encounter_id": eid, "SUBJECT_ID": s, "HADM_ID": str(hadm), "hadm_id_list": f"[{hadm}]",
                         "INTIME": ts(t0), "OUTTIME": ts(out), "Died": str(died), "t0": t0, "OUT": out,
                         "split": {"val": "validation"}.get(split_of[s], split_of[s]), "unit": unit})
            icustay = str(icu - 1)

            def ch(itemid, hrs, val, uom="", icu_id=icustay):
                chart.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": icu_id, "ITEMID": str(itemid),
                              "CHARTTIME": ts(t0 + pd.Timedelta(hours=hrs)), "VALUE": str(val),
                              "VALUENUM": str(val) if _isnum(val) else "", "VALUEUOM": uom})

            def lb(itemid, hrs, val, uom=""):
                lab.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ITEMID": str(itemid),
                            "CHARTTIME": ts(t0 + pd.Timedelta(hours=hrs)), "VALUE": str(val),
                            "VALUENUM": str(val) if _isnum(val) else "", "VALUEUOM": uom})
            if planted_enc and si == 0:
                planted["P0"] = eid
                ch(1003, 1, 100, "bpm"); ch(1003, 2, 120, "bpm"); ch(1002, 2, 90, "bpm")   # bin 0 median 110
                ch(1003, 7, 500, "bpm")                                                  # out of gate -> missing
                lb(5041, 22.5, 88.4, "μmol/L"); lb(5041, 23.5, 176.8, "μmol/L")          # bins 3 and 4 (release +1 h)
                lb(5227, 3, 2.0, "mmol/L")                                               # lactate bin 0, cap 1 bin
                ch(1016, 5, 90, "mmHg"); ch(1015, 5, 60, "mmHg")                         # MAP 70 in bin 0
                lb(5230, -3, 140, "mmol/L"); lb(5230, -20, 130, "mmol/L")                # pre-ICU block only
                ch(1014, 1, 10, "kg")                                                    # weight 10 kg
                ch(1008, 6, 120, "ml"); ch(1008, 30, 24, "ml")                           # UO 2.0 and 0.1 mL/kg/h
                ch(1001, 400 + 5, 45, "°C")                                              # after OUT: never in C
                # decision 2: respiratory rate 0 = missing; 2 below and 150 above the 5-120 gate; 30 in gate
                ch(1004, 1, 0, "bpm"); ch(1004, 2, 2, "bpm"); ch(1004, 3, 150, "bpm"); ch(1004, 4, 30, "bpm")
                # WBC gate 0.1-200: 0.05 and 0 below (0 stays out of gate for WBC), 250 above, 12 in gate
                lb(5141, 2, 0.05); lb(5141, 3, 250); lb(5141, 4, 0); lb(5141, 4.5, 12)
                # decision 5: same-draw base excess (5211 main, 5249) and bicarbonate (5248 main, 5224) -> one value
                lb(5211, 5, -2.0, "mmol/L"); lb(5249, 5, -3.0, "mmol/L")
                lb(5248, 5, 20, "mmol/L"); lb(5224, 5, 22, "mmol/L")
                # decision 4: unparsed text strings (some listed, some not)
                for jx, (txt, times) in enumerate(UNPARSED_PLANTED):
                    for kx in range(times):
                        lb(5326, 30 + jx + 0.05 * kx, txt)
                rx.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": icustay,
                           "STARTDATE": ts(t0 + pd.Timedelta(hours=1)), "ENDDATE": "",
                           "DRUG_NAME": DRUGS["vasoactive"][0], "DRUG_NAME_EN": DRUGS["vasoactive"][1]})
                rx.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": icustay,
                           "STARTDATE": ts(t0 - pd.Timedelta(hours=3)), "ENDDATE": "",
                           "DRUG_NAME": DRUGS["antibiotic"][0], "DRUG_NAME_EN": DRUGS["antibiotic"][1]})
                surg.append(_surg(s, hadm, t0, -1.0, 2.5, -1.5, 3.0, "asd closure"))            # recent (T = t0 + 3 h)
                surg.append(_surg(s, hadm, t0, 26.0, 29.0, 25.5, 30.0, "angiocardiography"))   # Layer C, bin 5
            elif planted_enc and si == 1:
                planted["P1"] = eid
                surg.append(_surg(s, hadm, t0, -10.0, 22.0, -10.5, -8.5, "anoplasty"))         # 32 h -> anesthesia window
                # decision 7: drugs at t0 from [t0 - 6 h, t0); an ENDDATE after t0 is not used, one before t0 is
                for g_, st_, en_ in (("sedative_analgesic", -5.0, 30.0), ("diuretic", -4.0, -1.0), ("antibiotic", -10.0, None)):
                    rx.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": icustay,
                               "STARTDATE": ts(t0 + pd.Timedelta(hours=st_)),
                               "ENDDATE": ts(t0 + pd.Timedelta(hours=en_)) if en_ is not None else "",
                               "DRUG_NAME": DRUGS[g_][0], "DRUG_NAME_EN": DRUGS[g_][1]})
            elif planted_enc and si == 2:
                planted["P2"] = eid
                surg.append(_surg(s, hadm, t0, -5.0, -6.0, None, -4.0, "anoplasty"))           # order + no anesthesia start
                # v1.0.4 per-channel tokens: positive+ / positive++ (4 channels), '<3.00' on mycoplasma (below cutoff)
                for it_, hrs_, txt_ in QUAL_V104_PLANTED:
                    lb(it_, hrs_, txt_)
                rx.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": icustay,
                           "STARTDATE": ts(t0 + pd.Timedelta(hours=2)), "ENDDATE": "",
                           "DRUG_NAME": DRUGS["osmotic"][0], "DRUG_NAME_EN": DRUGS["osmotic"][1]})
            elif planted_enc and si == 3:
                planted["P3"] = eid
                surg.append(_surg(s, hadm, t0, 40.0, 42.0, 39.5, None, "anoplasty"))           # T falls back to end
                rx.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": icustay,
                           "STARTDATE": (t0 + pd.Timedelta(hours=30)).strftime("%Y-%m-%d"), "ENDDATE": "",
                           "DRUG_NAME": DRUGS["diuretic"][0], "DRUG_NAME_EN": DRUGS["diuretic"][1]})  # date-only
            elif planted_enc and si == 4:
                planted["P4_NICU"] = eid
                lb(5255, 2, 171, "μmol/L")                                                  # kept (NICU)
            elif planted_enc and si == 5:
                planted["P5"] = eid
                lb(5255, 2, 171, "μmol/L")                                                  # dropped (not NICU)
                # decision 9: linkage by HADM_ID only; VISIT_ID never links
                sv = _surg(s, hadm, t0, -50.0, -48.0, -50.5, -47.5, "anoplasty")         # remote, VISIT_ID differs
                sv["VISIT_ID"] = "999999"
                surg.append(sv)
                sh = _surg(s, "999999", t0, -3.0, -2.0, -3.5, -1.5, "anoplasty")          # other HADM, VISIT = HADM
                sh["VISIT_ID"] = str(hadm)
                surg.append(sh)
            else:
                # random bulk events
                nv = max(1, int(min(span, 336) // 6))
                for k in range(nv):
                    hrs = 6 * k + float(rng.uniform(0, 6))
                    if hrs >= span:
                        break
                    for _, r in chart_items.iterrows():
                        if r["itemid"] in ("1008", "1009", "1014") or rng.random() < 0.45:
                            continue
                        v, u = _value(rng, r)
                        ch(r["itemid"], hrs, v, u)
                    if k % 4 == 0:
                        ch(1008, hrs, f"{rng.uniform(50, 800):.0f}", "ml")
                        ch(1009, hrs, f"{rng.uniform(100, 1500):.0f}", "ml")
                if rng.random() < 0.6:
                    ch(1014, float(rng.uniform(-10, 10)), f"{rng.uniform(3, 60):.1f}", "kg")
                for _, r in lab_items.iterrows():
                    for _j in range(int(rng.integers(0, 3))):
                        hrs = float(rng.uniform(-8, min(span, 336)))
                        v, u = _value(rng, r)
                        lb(r["itemid"], hrs, v, u)
                for g in ("vasoactive", "antibiotic", "sedative_analgesic", "diuretic", "none"):
                    if rng.random() < 0.4:
                        st = float(rng.uniform(-5, min(span, 200)))
                        for d in range(int(rng.integers(1, 5))):
                            rx.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": icustay,
                                       "STARTDATE": ts(t0 + pd.Timedelta(hours=st + 10 * d)),
                                       "ENDDATE": ts(t0 + pd.Timedelta(hours=st + 10 * d + 48)) if rng.random() < 0.2 else "",
                                       "DRUG_NAME": DRUGS[g][0], "DRUG_NAME_EN": DRUGS[g][1]})
                if rng.random() < 0.5:
                    micro.append({"SUBJECT_ID": s, "HADM_ID": str(hadm), "CHARTTIME": ts(t0 + pd.Timedelta(hours=float(rng.uniform(0, min(span, 300))))),
                                  "SPEC_ITEMID": "LIS0144", "SPEC_TYPE_DESC": "一般细菌+淋球菌涂片检查", "ORG_ITEMID": "",
                                  "ORG_NAME": ""})
                if rng.random() < 0.35:
                    a = float(rng.choice([rng.uniform(-60, -30), rng.uniform(-20, 2), rng.uniform(10, min(span, 300))]))
                    dur = float(rng.uniform(0.5, 4))
                    surg.append(_surg(s, hadm, t0, a, a + dur, a - 0.3, a + dur + 0.4, str(rng.choice(SURG_NAMES))))
                if padding_extremes:
                    ch(1003, span + 3, 299, "bpm"); ch(1003, span + 50, 31, "bpm")
                    lb(5041, span + 2, 1999, "μmol/L")
            t = (out + pd.Timedelta(hours=float(rng.uniform(100, 500)))).floor("s")
    enc = pd.DataFrame(encs)
    # files
    _csv(raw / "PATIENTS.csv", pats)
    _csv(raw / "ADMISSIONS.csv", adms)
    _csv(raw / "ICUSTAYS.csv", stays)
    _csv(raw / "CHARTEVENTS.csv", [dict(ROW_ID=i, **r) for i, r in enumerate(chart)])
    _csv(raw / "LABEVENTS.csv", [dict(ROW_ID=i, **r) for i, r in enumerate(lab)])
    _csv(raw / "PRESCRIPTIONS.csv", [dict(ROW_ID=i, **r) for i, r in enumerate(rx)])
    _csv(raw / "MICROBIOLOGYEVENTS.csv", [dict(ROW_ID=i, **r) for i, r in enumerate(micro)])
    _csv(raw / "SURGERY_INFO.csv", [dict(ROW_ID=str(70000 + i), **r) for i, r in enumerate(surg)],
         cols=["ROW_ID", "SUBJECT_ID", "HADM_ID", "VISIT_ID", "OPER_ID", "ANES_START_TIME", "ANES_END_TIME",
               "SURGERY_BEGIN_TIME", "SURGERY_END_TIME", "SURGERY_NAME", "ANES_METHOD", "SURGERY_POSITION"])
    _csv(raw / "D_ITEMS.csv", [{"ROW_ID": i, "ITEMID": it, "LABEL": "x"} for i, it in enumerate(chart_items["itemid"])])
    _csv(raw / "D_LABITEMS.csv", [{"ROW_ID": i, "ITEMID": it, "LABEL": "x"} for i, it in enumerate(lab_items["itemid"])])
    _csv(raw / "D_ICD_DIAGNOSES.csv", [{"ROW_ID": i, "ICD10_CODE_CN": k, "ICD10_CODE": v, "TITLE_CN": "t", "TITLE": "t"}
                                       for i, (k, v) in enumerate({**DX_CN, **D_ICD_EXTRA}.items())])
    _csv(raw / "SURGERY_VITAL_SIGNS.csv", [{"ROW_ID": 1, "SUBJECT_ID": "100001", "VALUE": "1"}])
    _csv(interim / "canonical_icu_encounters.csv",
         enc[["encounter_id", "SUBJECT_ID", "HADM_ID", "hadm_id_list", "INTIME", "OUTTIME", "Died"]].to_dict("records"))
    _csv(interim / "encounter_harmonization_fields.csv",
         [{"encounter_id": r.encounter_id, "t0": ts(r.t0), "OUTTIME": ts(r.OUT), "first_careunit": r.unit}
          for r in enc.itertuples()])
    _csv(root / "subject_splits.csv", [{"SUBJECT_ID": s, "split": split_of[s]} for s in subjects])
    names = ("PATIENTS", "ADMISSIONS", "ICUSTAYS", "SURGERY_INFO", "CHARTEVENTS", "LABEVENTS", "PRESCRIPTIONS",
             "MICROBIOLOGYEVENTS", "D_ITEMS", "D_LABITEMS")
    rows = {"PATIENTS": len(pats), "ADMISSIONS": len(adms), "ICUSTAYS": len(stays), "SURGERY_INFO": len(surg),
            "CHARTEVENTS": len(chart), "LABEVENTS": len(lab), "PRESCRIPTIONS": len(rx), "MICROBIOLOGYEVENTS": len(micro),
            "D_ITEMS": None, "D_LABITEMS": None}
    raw_pins = {nm: (sha(raw / f"{nm}.csv"), rows[nm]) for nm in names}
    out_h = (enc["OUT"] - enc["t0"]).dt.total_seconds() / 3600
    landmarks = (24, 48, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288, 312, 336)
    return {"root": root, "raw": raw, "interim": interim, "enc": enc, "planted": planted, "raw_pins": raw_pins,
            "pins": {"canonical": sha(interim / "canonical_icu_encounters.csv"),
                     "harmonization": sha(interim / "encounter_harmonization_fields.csv"),
                     "split": sha(root / "subject_splits.csv")},
            "d_icd": sha(raw / "D_ICD_DIAGNOSES.csv"), "svs": sha(raw / "SURGERY_VITAL_SIGNS.csv"),
            "expect": {"encounters": len(enc), "subjects": enc["SUBJECT_ID"].nunique(), "died": int(enc["Died"].astype(int).sum()),
                       "splits": {sp: int((enc["split"] == sp).sum()) for sp in ("train", "validation", "test")},
                       "eligible": {L: int((out_h > L).sum()) for L in landmarks}}}


def _isnum(v) -> bool:
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def _surg(s, hadm, t0, b, e, a_s, a_e, name) -> dict:
    f = lambda h: "" if h is None else ts(t0 + pd.Timedelta(hours=h))  # noqa: E731
    return {"SUBJECT_ID": s, "HADM_ID": str(hadm), "VISIT_ID": str(hadm), "OPER_ID": "1", "ANES_START_TIME": f(a_s),
            "ANES_END_TIME": f(a_e), "SURGERY_BEGIN_TIME": f(b), "SURGERY_END_TIME": f(e), "SURGERY_NAME": name,
            "ANES_METHOD": "GA", "SURGERY_POSITION": "supine"}


def _csv(path: Path, rows: list, cols: list | None = None) -> None:
    df = pd.DataFrame(rows)
    if cols:
        df = df[cols]
    df.to_csv(path, index=False, encoding="utf-8", quoting=csv.QUOTE_MINIMAL)


# ----------------------------------------------------------------------------------------------------------------------
# Synthetic Stage 2 v1.0.3-format release (builder v1.0.4 cross-release tests; SYNTHETIC VALUES ONLY)
# ----------------------------------------------------------------------------------------------------------------------
S2_L1 = ("temperature heart_rate respiratory_rate spo2 urine_output fluid_intake glucose_bedside pain_score dbp sbp "
         "albumin alt ast bun bilirubin_direct ggt ldh magnesium phosphate bilirubin_total anc hemoglobin "
         "lymphocytes_abs platelets wbc ptt d_dimer fibrinogen inr base_excess chloride bicarbonate potassium lactate "
         "sodium pco2 ph sao2_bloodgas crp cystatin_c").split()
S2_SPLITS = ("train", "validation", "test")
S2_TABLES = ("CHARTEVENTS", "LABEVENTS", "PRESCRIPTIONS", "MICROBIOLOGYEVENTS")
S2_BANDS = ("neonate", "infant", "toddler", "preschool", "school", "adolescent", "ge18", "age_invalid")
S2_UNITS = ("CICU", "General ICU", "NICU", "PICU", "SICU")
S2_GE18 = tuple(f"/layer_a/{sp}/age_band/n_ge18" for sp in S2_SPLITS)
S2_STOP_EVENT = ("/layer_a/dob_breakdown/n_invalid", "/layer_a/dob_breakdown/status/n_valid",
                 "/layer_a/dob_breakdown/n_age_days")      # the stop-event count and the two cells it fixes


def _dist(v: int, n: int) -> list[int]:
    q = v // n
    return [v - q * (n - 1)] + [q] * (n - 1)


def stage2_internal_links(links_path: Path | None = None) -> list[tuple[str, list[str]]]:
    """The Stage 2 internal equalities of the links file (rows whose pointers are all 's2:'), without the prefix."""
    lp = Path(links_path) if links_path else HERE / "stage3_links_v3.csv"
    df = pd.read_csv(lp, dtype=str, keep_default_na=False, comment="#")
    out = []
    for _, r in df.iterrows():
        cells = [r["total"].strip()] + [q.strip() for q in r["parts"].split(";") if q.strip()]
        if r["kind"].strip() == "eq" and all(c.startswith("s2:/") for c in cells):
            out.append((cells[0][3:], [c[3:] for c in cells[1:]]))
    return out


def make_stage2_release(path: Path, s3: dict | None = None, splits: dict | None = None,
                        override: dict | None = None, statuses: dict | None = None, protect: bool = False,
                        known: tuple | list | None = None, links_path: Path | None = None) -> str:
    """Write a small Stage 2 v1.0.3-format release JSON consistent with the Stage 3 audit values `s3` (pointer ->
    true value, as built by the Stage 3 builder on the fixture) under stage3_links_v3.csv (every pointer of its Stage 2
    internal rows exists: age bands incl. ge18, strata, care units, DOB breakdown by split / care unit, recent
    surgery); s3=None writes the same structure with placeholder values. `override` sets values, `statuses` forces
    '<10' / 'S' / int on given pointers. Returns the SHA-256.
    protect=False: primary <10 and a simple residual complementary suppression (one 'S' per identity left with a
    single unknown), as in the v1.0.2 fixture (it can leave cells exposed: the compromised-cell tests).
    protect=True: complementary suppression by the pinned Stage 2 v1.0.3 module v1.5 (Suppressor) over the identities
    and the Stage 2 internal links, with the denominators and the four table sizes (D3) frozen and the `known` pointers (default: the three n_ge18
    structural zeros) as known values (labels kept), as the real Stage 2 v1.0.3 release was gated."""
    s3 = s3 or {}
    g = lambda p, d=100: int(s3.get(p, d))  # noqa: E731
    sp_n = splits or {"train": 60, "validation": 20, "test": 20}
    N = sum(sp_n.values())
    V: dict[str, int] = {}
    ids: list[tuple[str, list[str]]] = []
    hidden: dict[str, int] = {}

    def ident(t, parts):
        ids.append((t, list(parts)))
    for sp in S2_SPLITS:
        V[f"/denominators/{sp}/n_encounters"] = sp_n[sp]
    V["/denominators/n_encounters"] = N
    ident("/denominators/n_encounters", [f"/denominators/{sp}/n_encounters" for sp in S2_SPLITS])
    for ch in S2_L1:
        b = f"/map_checks/row_status/{ch}"
        V[f"{b}/n_measurement"] = g(f"/plausibility/{ch}/train/n_rows_in_gate")
        # Stage 2 counts a respiratory rate of 0 as out of gate; Stage 3 reads it as missing (decision 2)
        V[f"{b}/n_out_of_gate"] = (g(f"/plausibility/{ch}/train/n_rows_out_of_gate")
                                   + int(s3.get(f"/plausibility/{ch}/train/n_rows_zero_as_missing", 0)))
        for k, x in (("n_nonnumeric", 12), ("n_blank", 3), ("n_unit_pending", 0), ("n_unit_unrecognized", 0)):
            V[f"{b}/{k}"] = x
        parts = [f"{b}/{k}" for k in ("n_measurement", "n_out_of_gate", "n_nonnumeric", "n_blank", "n_unit_pending",
                                      "n_unit_unrecognized")]
        V[f"{b}/n_total"] = sum(V[p] for p in parts)
        ident(f"{b}/n_total", parts)
    for t in S2_TABLES:
        b = f"/reconciliation/{t}"
        ins = _dist(g(f"/event_linkage/{t}/n_rows_in_stay"), 3)
        outs = _dist(g(f"/event_linkage/{t}/n_rows_outside_every_window"), 3)
        parts = []
        for sp, x in zip(S2_SPLITS, ins):
            V[f"{b}/n_in_stay_{sp}"] = x
            parts.append(f"{b}/n_in_stay_{sp}")
        for k, x in zip(("n_cohort_admission_outside_stay", "n_no_cohort_admission", "n_cohort_admission_time_missing"), outs):
            V[f"{b}/{k}"] = x
            parts.append(f"{b}/{k}")
        V[f"{b}/n_total"] = sum(V[p] for p in parts)
        ident(f"{b}/n_total", parts)
    # ---- DOB breakdown: overall, by split, by care unit, by split x care unit (v1.0.3 structure)
    d = "/layer_a/dob_breakdown"
    m_t, u_t = g("/age/invalid_detail/n_missing", 4), g("/age/invalid_detail/n_unparseable", 4)
    nb_t = [g(f"/age/invalid_detail/negative_magnitude/{k}", 0) for k in ("n_lt1d", "n_1to30d", "n_gt30d")]
    neg_t = g("/age/invalid_detail/n_negative", sum(nb_t))
    if sum(nb_t) != neg_t:                                  # placeholder values: all negatives in the first bin
        nb_t = [neg_t, 0, 0]
    inv_t = m_t + u_t + neg_t
    inv = dict(zip(S2_SPLITS, _dist(inv_t, 3)))
    pool = ["missing"] * m_t + ["unparseable"] * u_t + ["lt1d"] * nb_t[0] + ["1to30d"] * nb_t[1] + ["gt30d"] * nb_t[2]
    dob, k0 = {}, 0
    for sp in S2_SPLITS:
        take = pool[k0:k0 + inv[sp]]
        k0 += inv[sp]
        dob[sp] = {"missing": take.count("missing"), "unparseable": take.count("unparseable"),
                   "bins": [take.count(x) for x in ("lt1d", "1to30d", "gt30d")]}

    def dob_block(base: str, den: str, n_all: int, c: dict) -> None:
        neg = sum(c["bins"])
        st = {"n_valid": n_all - c["missing"] - c["unparseable"] - neg, "n_missing": c["missing"],
              "n_unparseable": c["unparseable"], "n_negative": neg}
        for k, x in st.items():
            V[f"{base}/status/{k}"] = x
        ident(den, [f"{base}/status/{k}" for k in st])
        V[f"{base}/n_invalid"] = c["missing"] + c["unparseable"] + neg
        ident(f"{base}/n_invalid", [f"{base}/status/{k}" for k in ("n_missing", "n_unparseable", "n_negative")])
        V[f"{base}/n_age_days"] = st["n_valid"]
        ident(f"{base}/n_age_days", [f"{base}/status/n_valid"])
        for k, x in zip(("n_lt1d", "n_1to30d", "n_gt30d"), c["bins"]):
            V[f"{base}/negative_bins/{k}"] = x
        ident(f"{base}/status/n_negative", [f"{base}/negative_bins/{k}" for k in ("n_lt1d", "n_1to30d", "n_gt30d")])
    dob_block(d, "/denominators/n_encounters", N, {"missing": m_t, "unparseable": u_t, "bins": nb_t})
    for sp in S2_SPLITS:
        dob_block(f"{d}/by_split/{sp}", f"/denominators/{sp}/n_encounters", sp_n[sp], dob[sp])
    if s3 and "/age/n_valid" in s3 and int(s3["/age/n_valid"]) != V[f"{d}/status/n_valid"]:
        raise ValueError("fixture: Stage 3 /age/n_valid differs from cohort - invalid (synthetic values)")
    wk = ("n_charted", "n_imputed_age_sex", "n_imputed_age_band", "n_imputed_sex_only", "n_imputed_overall")
    for j, k in enumerate(wk):
        V[f"{d}/weight_imputation/age_invalid/{k}"] = V[f"{d}/n_invalid"] if j == 0 else 0
    ident(f"{d}/n_invalid", [f"{d}/weight_imputation/age_invalid/{k}" for k in wk])
    # care units per split (synthetic shares); invalid DOBs: negatives -> NICU, missing / unparseable -> PICU
    units = {}
    for sp in S2_SPLITS:
        n = sp_n[sp]
        sh = [int(n * x) for x in (0.1, 0.2, 0.3, 0.3)]
        units[sp] = dict(zip(S2_UNITS, sh + [n - sum(sh)]))
    ucnt = {u: {"missing": 0, "unparseable": 0, "negative": 0} for u in S2_UNITS}
    for sp in S2_SPLITS:
        for u in S2_UNITS:
            c = {"missing": dob[sp]["missing"] if u == "PICU" else 0,
                 "unparseable": dob[sp]["unparseable"] if u == "PICU" else 0,
                 "negative": sum(dob[sp]["bins"]) if u == "NICU" else 0}
            base = f"{d}/by_split_and_first_careunit/{sp}/{u}"
            for k, x in c.items():
                V[f"{base}/status/n_{k}"] = x
                ucnt[u][k] += x
            V[f"{base}/status/n_total"] = sum(c.values())
            ident(f"{base}/status/n_total", [f"{base}/status/n_{k}" for k in c])
            V[f"{base}/n_invalid"] = sum(c.values())
            ident(f"{base}/n_invalid", [f"{base}/status/n_{k}" for k in c])
    for u in S2_UNITS:
        base = f"{d}/by_first_careunit/{u}"
        for k, x in ucnt[u].items():
            V[f"{base}/status/n_{k}"] = x
        V[f"{base}/status/n_total"] = sum(ucnt[u].values())
        ident(f"{base}/status/n_total", [f"{base}/status/n_{k}" for k in ucnt[u]])
        V[f"{base}/n_invalid"] = sum(ucnt[u].values())
        ident(f"{base}/n_invalid", [f"{base}/status/n_{k}" for k in ucnt[u]])
    ident(f"{d}/n_invalid", [f"{d}/by_first_careunit/{u}/n_invalid" for u in S2_UNITS])
    miss = _dist(g("/admission_interval/withheld/n_missing", 0), 3)
    neg = _dist(g("/admission_interval/n_negative_clipped_to_0", 0), 3)
    dxm = _dist(g("/diagnosis/n_missing", 0), 3)
    lc = _dist(g("/surgery/timing/n_layer_c", 0), 3)
    for j, sp in enumerate(S2_SPLITS):
        den = f"/denominators/{sp}/n_encounters"
        # age bands (8 levels in layer_a; 7 strata levels: no encounter is >= 18 years, so ge18 is no strata level)
        rest = sp_n[sp] - inv[sp]
        bb = [int(rest * x) for x in (0.2, 0.3, 0.2, 0.12, 0.12)]
        bands = dict(zip(S2_BANDS, bb + [rest - sum(bb), 0, inv[sp]]))
        a = f"/layer_a/{sp}/age_band"
        for lv, x in bands.items():
            V[f"{a}/n_{lv}"] = x
        ident(den, [f"{a}/n_{lv}" for lv in S2_BANDS])
        for lv in S2_BANDS:
            if lv != "ge18":
                V[f"/strata_denominators/age_band/{sp}/{lv}/n_encounters"] = bands[lv]
        ident(den, [f"/strata_denominators/age_band/{sp}/{lv}/n_encounters" for lv in S2_BANDS if lv != "ge18"])
        fu = f"/layer_a/{sp}/first_careunit"
        for u, x in units[sp].items():
            V[f"{fu}/n_{u}"] = x
            V[f"/strata_denominators/first_careunit/{sp}/{u}/n_encounters"] = x
        ident(den, [f"{fu}/n_{u}" for u in S2_UNITS])
        ident(den, [f"/strata_denominators/first_careunit/{sp}/{u}/n_encounters" for u in S2_UNITS])
        rec = max(12, sp_n[sp] // 5)
        V[f"/surgery/layer_a/{sp}/n_encounters_recent"] = rec
        V[f"/strata_denominators/surgical_status/{sp}/recent_surgery/n_encounters"] = rec
        V[f"/strata_denominators/surgical_status/{sp}/no_recent_surgery/n_encounters"] = sp_n[sp] - rec
        ident(den, [f"/strata_denominators/surgical_status/{sp}/no_recent_surgery/n_encounters",
                    f"/strata_denominators/surgical_status/{sp}/recent_surgery/n_encounters"])
        h = f"/layer_a/{sp}/hours_admit_to_icu_band"
        V[f"{h}/n_missing"] = miss[j]
        V[f"{h}/n_lt24h"] = sp_n[sp] - miss[j]
        ident(den, [f"{h}/n_missing", f"{h}/n_lt24h"])
        V[f"/layer_a/{sp}/n_hours_admit_to_icu_negative"] = neg[j]
        hc = f"hidden:/layer_a/{sp}/n_hours_admit_to_icu_negative|complement|sub"
        hidden[hc] = sp_n[sp] - neg[j]
        ident(den, [f"/layer_a/{sp}/n_hours_admit_to_icu_negative", hc])
        c = f"/layer_a/{sp}/diagnosis_code_status"
        V[f"{c}/n_missing"] = dxm[j]
        V[f"{c}/n_ok"] = sp_n[sp] - dxm[j]
        ident(den, [f"{c}/n_ok", f"{c}/n_missing"])
        V[f"/surgery/layer_c_extra/{sp}/n_surgeries"] = lc[j]
    ident(f"{d}/n_invalid", [f"{d}/by_split/{sp}/n_invalid" for sp in S2_SPLITS])
    V["/surgery/link_classes/n_linked"] = g("/surgery/n_linked", 0)
    V["/surgery/link_classes/n_other_hadm"] = 25
    V["/surgery/link_classes/n_total"] = V["/surgery/link_classes/n_linked"] + 25
    ident("/surgery/link_classes/n_total", ["/surgery/link_classes/n_linked", "/surgery/link_classes/n_other_hadm"])
    V.update(override or {})
    st = {p: ("P" if v < 10 else "pub") for p, v in V.items()}
    st.update({p: {"<10": "P", "S": "S"}.get(x, "pub") for p, x in (statuses or {}).items()})
    frozen_tot = {p for p in V if p.startswith("/denominators/")}
    if protect:   # decided D3 (Stage 2 v1.0.3): the raw table sizes are public and published
        frozen_tot |= {f"/reconciliation/{t}/n_total" for t in S2_TABLES}
        st = _protect_stage2(V, st, ids, hidden, frozen_tot, S2_GE18 if known is None else tuple(known), links_path)
    else:
        # statuses: primary <10, then one complementary S per identity left with a single unknown (smallest published part)
        for _ in range(50):
            changed = False
            for t, parts in ids:
                cells = [t] + parts
                unk = [c for c in cells if c in hidden or st[c] != "pub"]
                if len(unk) == 1 and not unk[0] in hidden:
                    pick = sorted((V[c], c) for c in parts if st[c] == "pub" and c not in frozen_tot)
                    if pick:
                        st[pick[0][1]] = "S"
                        changed = True
            if not changed:
                break
    doc: dict = {"spec": "synthetic Stage 2 v1.0.3-format release (fixture)", "release_gate": "ok",
                 "release_v1_0_3": {"script": "synthetic fixture (not the module)", "protected": bool(protect)}}
    for p, v in V.items():
        cur = doc
        ks = [k.replace("~1", "/").replace("~0", "~") for k in p.split("/")[1:]]
        for k in ks[:-1]:
            cur = cur.setdefault(k, {})
        cur[ks[-1]] = v if st[p] == "pub" else ("<10" if st[p] == "P" else "S")
    doc["disclosure"] = {"count_fields": sorted(V), "identities": [[t, p] for t, p in ids],
                         "hidden_cells": sorted(hidden), "ratios": [], "quantiles": [], "dependents": []}
    Path(path).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    _LAST["V"] = dict(V)
    return sha(Path(path))


def _protect_stage2(V: dict, st: dict, ids: list, hidden: dict, frozen: set, known: tuple,
                    links_path: Path | None) -> dict:
    """Complementary suppression of the synthetic Stage 2 cells by the pinned module v1.5 (as Stage 2 v1.0.3 was
    gated): identities + the Stage 2 internal links, denominators frozen, `known` cells known (labels kept)."""
    import importlib.util
    import sys as _sys
    spec = importlib.util.spec_from_file_location("s2_fixture_module", HERE / "stage2_v1_0_3_suppress.py")
    SM = importlib.util.module_from_spec(spec)
    _sys.modules["s2_fixture_module"] = SM
    spec.loader.exec_module(SM)
    rels = [(t, tuple(p), "identity") for t, p in ids]
    rels += [(t, tuple(p), "link") for t, p in stage2_internal_links(links_path)]
    vals = dict(V)
    vals.update(hidden)
    status = dict(st)
    for h in hidden:
        status[h] = "H"
    label = {}
    for c in known:
        label[c] = status[c]
        status[c] = "const"
    params = SM.Params()
    kinds = dict(status)
    system = SM.System(rels, vals, kinds, params)
    _q = SM.log
    SM.log = lambda msg: None          # quiet in tests
    try:
        sup = SM.Suppressor(system, status, set(frozen))
        sup.run()
        res = system.audit(status, label="fixture check")
    finally:
        SM.log = _q
    if res.unsafe:
        raise RuntimeError(f"fixture: synthetic Stage 2 protection left {len(res.unsafe)} failing cell(s)")
    out = {p: status[p] for p in V}
    out.update(label)
    return out


_LAST: dict = {}


def make_stage2_constants(path: Path, stop_event: bool = True, values: dict | None = None,
                          extra_rows: list | None = None) -> str:
    """Synthetic constants file (the v1.0.4 --stage2-constants format: pointer,value,reason; s2: pointers) for the
    last release written by make_stage2_release: the three n_ge18 structural zeros and, with stop_event, the
    invalid-DOB total (its synthetic value) with the two cells it fixes as 'derived'. Returns the SHA-256."""
    V = values if values is not None else _LAST["V"]
    rows = [(f"s2:{p}", "0", "structural zero (synthetic)") for p in S2_GE18]
    if stop_event:
        rows += [("s2:" + S2_STOP_EVENT[0], str(V[S2_STOP_EVENT[0]]), "stop-event count (synthetic)"),
                 ("s2:" + S2_STOP_EVENT[1], "derived", "cohort minus the invalid count (synthetic)"),
                 ("s2:" + S2_STOP_EVENT[2], "derived", "equals the valid count (synthetic)")]
    rows += list(extra_rows or [])
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("# synthetic constants (fixture)\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["pointer", "value", "reason"])
        w.writerows(rows)
    return sha(Path(path))
