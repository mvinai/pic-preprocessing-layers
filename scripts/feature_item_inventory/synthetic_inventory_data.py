"""Synthetic PIC-like inputs for testing feature_item_inventory.py. NO real data.

Builds a fake project (canonical encounters, harmonization fields, subject split, D_ITEMS, D_LABITEMS and the six
raw event tables) with planted boundary cases, plus an independent slow reference of the expected counts."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

FMT = "%Y-%m-%d %H:%M:%S"
SENTINEL_SUBJECT = "987654321"
LANDMARKS = (24, 48, 72, 120, 168, 240, 336)

CHART_ITEMS = {"220045": ("Heart Rate", "Vitals", "bpm", 0.95), "224639": ("Weight", "General", "kg", 0.6),
               "226559": ("Urine Output", "Output", "ml", 0.5), "220739": ("GCS Eye", "Neuro", "", 0.05),
               "777001": ("Rare chart item", "Misc", "u", 0.01)}
CHART_NOT_IN_DICT = "888001"
LAB_ITEMS = {"50912": ("Creatinine", "Blood", "Chemistry", "mg/dL", 0.8), "50983": ("Sodium", "Blood", "Chemistry", "mmol/L", 0.8),
             "51001": ("Cortisol", "Blood", "Endocrine", "ug/dL", 0.03), "59999": ("Occult blood", "Stool", "Hematology", "", 0.02)}
PLANT_LAB = "55555"   # planted release-time boundary item
DRUGS = [("多巴胺", "Dopamine", "mg", 0.4), ("维库溴铵", "Vecuronium", "mg", 0.15), ("左乙拉西坦", "Levetiracetam", "mg", 0.12),
         ("头孢曲松", "Ceftriaxone", "g", 0.5), ("罕见药", "Rare drug/x", "mg", 0.01)]
SPECIMENS = [("70012", "BLOOD CULTURE", 0.4), ("70079", "URINE", 0.3), ("70081", "CSF", 0.02)]
ORGS = [("80023", "STAPHYLOCOCCUS AUREUS", 0.2), ("", "", 0.5), ("80155", "CANDIDA ALBICANS", 0.02)]


def ts(x: pd.Timestamp) -> str:
    return x.strftime(FMT)


def make_project(root: Path, n_subjects: int = 360, seed: int = 11) -> dict:
    rng = np.random.default_rng(seed)
    root = Path(root)
    raw = root / "raw"
    (root / "data" / "interim").mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)
    subjects = [str(200000 + i) for i in range(n_subjects - 1)] + [SENTINEL_SUBJECT]
    split_of = {}
    for i, s in enumerate(subjects):
        u = rng.random()
        split_of[s] = "train" if u < 0.6 else ("val" if u < 0.8 else "test")
    split_of[SENTINEL_SUBJECT] = "train"
    encs = []
    hadm = 600000
    base = pd.Timestamp("2016-01-01 00:00:00")
    for s in subjects:
        n_enc = 1 if rng.random() < 0.85 else 2
        t = base + pd.Timedelta(days=int(rng.integers(0, 1500)), hours=int(rng.integers(0, 24)), minutes=int(rng.integers(0, 60)))
        same = n_enc == 2 and rng.random() < 0.5
        first_h = None
        for e in range(n_enc):
            if same and e == 1:
                hl = [first_h]
            else:
                hadm += 1
                hl = [str(hadm)]
                if rng.random() < 0.08:
                    hadm += 1
                    hl.append(str(hadm))
            first_h = first_h or hl[0]
            t0 = t
            out = (t0 + pd.Timedelta(hours=float(rng.choice([rng.uniform(5, 30), rng.uniform(30, 500)])))).floor("s")
            encs.append({"encounter_id": f"{s}_{e + 1}", "SUBJECT_ID": s, "HADM_ID": hl[0], "hadm_list": hl,
                         "t0": t0, "OUT": out, "split": {"val": "validation"}.get(split_of[s], split_of[s])})
            t = (out + pd.Timedelta(hours=float(rng.uniform(30, 300)))).floor("s")
    enc = pd.DataFrame(encs)
    hl_fmt = ["[" + ", ".join(h) + "]" if i % 3 else str(h) for i, h in enumerate(enc["hadm_list"])]
    canon = pd.DataFrame({"encounter_id": enc["encounter_id"], "SUBJECT_ID": enc["SUBJECT_ID"],
                          "HADM_ID": enc["HADM_ID"].map(lambda h: f"{h}.0"), "hadm_id_list": hl_fmt,
                          "ICUSTAY_ID": range(1, len(enc) + 1), "INTIME": enc["t0"].map(ts), "OUTTIME": enc["OUT"].map(ts)})
    canon_p = root / "data" / "interim" / "canonical_icu_encounters.csv"
    canon.to_csv(canon_p, index=False)
    harm = pd.DataFrame({"encounter_id": enc["encounter_id"], "OUTTIME": enc["OUT"].map(ts), "t0": enc["t0"].map(ts),
                         "span_icu_hours": ((enc["OUT"] - enc["t0"]).dt.total_seconds() / 3600).round(4)})
    harm_p = root / "data" / "interim" / "encounter_harmonization_fields.csv"
    harm.to_csv(harm_p, index=False)
    split_p = root / "subject_splits.csv"
    pd.DataFrame({"SUBJECT_ID": subjects, "split": [split_of[s] for s in subjects]}).to_csv(split_p, index=False)

    pd.DataFrame([{"ROW_ID": i, "ITEMID": k, "LABEL": v[0], "CATEGORY": v[1], "UNITNAME": v[2]}
                  for i, (k, v) in enumerate(CHART_ITEMS.items())]
                 + [{"ROW_ID": 99, "ITEMID": "90001", "LABEL": "Infusion volume", "CATEGORY": "Input", "UNITNAME": "ml"}]
                 ).to_csv(raw / "D_ITEMS.csv", index=False)
    pd.DataFrame([{"ROW_ID": i, "ITEMID": k, "LABEL": v[0], "FLUID": v[1], "CATEGORY": v[2]} for i, (k, v) in enumerate(LAB_ITEMS.items())]
                 + [{"ROW_ID": 50, "ITEMID": PLANT_LAB, "LABEL": "Planted lab", "FLUID": "Blood", "CATEGORY": "Test"}]
                 ).to_csv(raw / "D_LABITEMS.csv", index=False)

    def times(t0, out, k):
        span = (out - t0).total_seconds() / 3600
        return [(t0 + pd.Timedelta(hours=float(h))).floor("s") for h in rng.uniform(-10, span + 10, size=k)]

    chart, lab, outp, inp, rx, micro = [], [], [], [], [], []
    for r in enc.itertuples():
        hadm_pick = r.hadm_list
        for it, (_, _, unit, p) in CHART_ITEMS.items():
            if rng.random() < p:
                for t in times(r.t0, r.OUT, int(rng.integers(1, 6))):
                    chart.append([r.SUBJECT_ID, rng.choice(hadm_pick), it, ts(t), "" if rng.random() < 0.1 else f"{rng.uniform(1, 200):.1f}", unit])
        if rng.random() < 0.05:
            chart.append([r.SUBJECT_ID, hadm_pick[0], CHART_NOT_IN_DICT, ts(r.t0 + pd.Timedelta(hours=1)), "1", ""])
        for it, (_, _, _, unit, p) in LAB_ITEMS.items():
            if rng.random() < p:
                for t in times(r.t0, r.OUT, int(rng.integers(1, 4))):
                    lab.append([r.SUBJECT_ID, rng.choice(hadm_pick), it, ts(t), f"{rng.uniform(0, 10):.2f}", unit])
        if rng.random() < 0.4:
            outp.append([r.SUBJECT_ID, hadm_pick[0], "226559", ts(times(r.t0, r.OUT, 1)[0]), f"{rng.uniform(1, 50):.0f}", "ml"])
        if rng.random() < 0.4:
            inp.append([r.SUBJECT_ID, hadm_pick[0], ts(times(r.t0, r.OUT, 1)[0]), f"{rng.uniform(1, 50):.0f}", "ml" if rng.random() < 0.7 else "mg"])
        for dn, de, unit, p in DRUGS:
            if rng.random() < p:
                for t in times(r.t0, r.OUT, int(rng.integers(1, 3))):
                    rx.append([r.SUBJECT_ID, hadm_pick[0], ts(t), "", dn, de, unit])
        for sid, sdesc, p in SPECIMENS:
            if rng.random() < p:
                oid, oname, _ = ORGS[int(rng.choice(len(ORGS), p=[0.3, 0.68, 0.02]))]
                micro.append([r.SUBJECT_ID, hadm_pick[0], ts(times(r.t0, r.OUT, 1)[0]), sid, sdesc, oid, oname, "OXACILLIN"])
    # planted release-time boundary labs (train, long stays): released exactly at t0, at t0 + 24 h, and after OUT
    planted = enc[(enc["split"] == "train") & ((enc["OUT"] - enc["t0"]).dt.total_seconds() > 48 * 3600)].head(12)
    for r in planted.itertuples():
        lab.append([r.SUBJECT_ID, r.hadm_list[0], PLANT_LAB, ts(r.t0 - pd.Timedelta(hours=1)), "1.0", "u"])          # release t0: stay + Layer B
        lab.append([r.SUBJECT_ID, r.hadm_list[0], PLANT_LAB, ts(r.t0 + pd.Timedelta(hours=23)), "1.0", "u"])         # release t0+24h: stay, not B
        lab.append([r.SUBJECT_ID, r.hadm_list[0], PLANT_LAB, ts(r.OUT - pd.Timedelta(minutes=30)), "1.0", "u"])      # release after OUT
        lab.append([r.SUBJECT_ID, r.hadm_list[0], PLANT_LAB, ts(r.t0 - pd.Timedelta(hours=7)), "1.0", "u"])          # release t0-6h: pre-ICU
        lab.append([r.SUBJECT_ID, r.hadm_list[0], PLANT_LAB, ts(r.t0 - pd.Timedelta(hours=7, seconds=1)), "1.0", "u"])  # just before pre window
    # noise: non-cohort subject, blank HADM, other HADM, missing time
    e0 = enc.iloc[0]
    for _ in range(15):
        chart.append(["111", "999999", "220045", ts(e0.t0), "80", "bpm"])
        chart.append([e0.SUBJECT_ID, "", "220045", ts(e0.t0 + pd.Timedelta(hours=1)), "80", "bpm"])
        chart.append([e0.SUBJECT_ID, "123", "220045", ts(e0.t0 + pd.Timedelta(hours=1)), "80", "bpm"])
        chart.append([e0.SUBJECT_ID, e0.hadm_list[0], "220045", "", "80", "bpm"])
        lab.append([e0.SUBJECT_ID, e0.hadm_list[0], "50912", "not a time", "1", "mg/dL"])
    pd.DataFrame(chart, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUENUM", "VALUEUOM"]).assign(
        ICUSTAY_ID="", VALUE="x", STORETIME="").to_csv(raw / "CHARTEVENTS.csv", index=False)
    pd.DataFrame(lab, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUENUM", "VALUEUOM"]).assign(
        VALUE="x", FLAG="").to_csv(raw / "LABEVENTS.csv", index=False)
    pd.DataFrame(outp, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"]).to_csv(raw / "OUTPUTEVENTS.csv", index=False)
    pd.DataFrame(inp, columns=["SUBJECT_ID", "HADM_ID", "CHARTTIME", "AMOUNT", "AMOUNTUOM"]).to_csv(raw / "INPUTEVENTS.csv", index=False)
    pd.DataFrame(rx, columns=["SUBJECT_ID", "HADM_ID", "STARTDATE", "ENDDATE", "DRUG_NAME", "DRUG_NAME_EN", "DOSE_UNIT_RX"]).to_csv(
        raw / "PRESCRIPTIONS.csv", index=False)
    pd.DataFrame(micro, columns=["SUBJECT_ID", "HADM_ID", "CHARTTIME", "SPEC_ITEMID", "SPEC_TYPE_DESC", "ORG_ITEMID", "ORG_NAME", "AB_NAME"]).to_csv(
        raw / "MICROBIOLOGYEVENTS.csv", index=False)
    split_counts = {s: int((enc["split"] == s).sum()) for s in ("train", "validation", "test")}
    elig = {L: int((enc["OUT"] > enc["t0"] + pd.Timedelta(hours=L)).sum()) for L in LANDMARKS}
    raw_rows = {"CHARTEVENTS": len(chart), "LABEVENTS": len(lab), "PRESCRIPTIONS": len(rx), "MICROBIOLOGYEVENTS": len(micro)}
    return {"root": root, "raw": raw, "canonical": canon_p, "harmonization": harm_p, "split": split_p, "enc": enc,
            "expected_encounters": len(enc), "expected_subjects": enc["SUBJECT_ID"].nunique(),
            "split_counts": split_counts, "eligible": elig, "raw_rows": raw_rows, "planted": planted}


def reference_counts(proj: dict, table: str, key_fn, time_col: str, offset_h: float = 0.0) -> dict:
    """Slow, independent reference: per item, train rows in stay and train encounters in stay / Layer B / pre-ICU."""
    enc = proj["enc"]
    df = pd.read_csv(proj["raw"] / f"{table}.csv", dtype=str, keep_default_na=False)
    by = {}
    for r in enc.itertuples():
        for h in r.hadm_list:
            by.setdefault((r.SUBJECT_ID, h), []).append(r)
    out: dict = {}
    classes = {"in_stay_train": 0, "in_stay_validation": 0, "in_stay_test": 0, "cohort_admission_outside_stay": 0,
               "no_cohort_admission": 0, "cohort_admission_time_missing": 0}
    for _, row in df.iterrows():
        h = row["HADM_ID"].strip()
        h = h[:-2] if h.endswith(".0") else h
        cands = by.get((row["SUBJECT_ID"].strip(), h), []) if h else []
        if not cands:
            classes["no_cohort_admission"] += 1
            continue
        t = pd.to_datetime(row[time_col], errors="coerce", format="%Y-%m-%d %H:%M:%S")
        if pd.isna(t):
            classes["cohort_admission_time_missing"] += 1
            continue
        t = t + pd.Timedelta(hours=offset_h)
        k = key_fn(row)
        hit = None
        for e in cands:
            if e.t0 <= t < e.OUT:
                hit = e
            if e.split == "train" and e.t0 - pd.Timedelta(hours=6) <= t < e.t0:
                out.setdefault(k, {"rows": 0, "enc": set(), "lb": set(), "pre": set()})["pre"].add(e.encounter_id)
        if hit is None:
            classes["cohort_admission_outside_stay"] += 1
            continue
        classes[f"in_stay_{hit.split}"] += 1
        if hit.split != "train":
            continue
        o = out.setdefault(k, {"rows": 0, "enc": set(), "lb": set(), "pre": set()})
        o["rows"] += 1
        o["enc"].add(hit.encounter_id)
        if t < min(hit.OUT, hit.t0 + pd.Timedelta(hours=24)):
            o["lb"].add(hit.encounter_id)
    res = {k: {"rows": v["rows"], "encounters": len(v["enc"]), "encounters_layer_b": len(v["lb"]),
               "encounters_pre_icu_6h": len(v["pre"])} for k, v in out.items() if v["rows"] > 0}
    return {"items": res, "classes": classes}
