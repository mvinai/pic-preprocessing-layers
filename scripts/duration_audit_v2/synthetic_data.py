"""Synthetic PIC-like inputs for testing audit_surgery_complexity_duration.py. NO real data.

Generates fake subjects/encounters/surgeries with the real column schemas, including multi-HADM encounters,
JSON / Python-repr / float-repr / delimited hadm_id_list strings, float-like IDs ('123.0'), missing times,
zero and negative durations, surgeries straddling INTIME, after OUTTIME, in other HADMs, for subjects not in the
cohort, small groups, and an inverse-order system (Hepatobiliary: I2 longer than I3)."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

SENTINEL_SUBJECT = "987654321"
FMT = "%Y-%m-%d %H:%M:%S"

# (raw_name, body_system, intensity, STAT, base_hours)
SYN_NAMES = [
    ("asd closure", "Cardiac", "major open", "1", 2.2), ("vsd repair", "Cardiac", "major open", "2", 2.8),
    ("pda ligation", "Cardiac", "major open", "1", 1.6), ("balloon valvotomy", "Cardiac", "catheter/interventional", "2", 1.9),
    ("tof repair", "Cardiac", "major open", "3", 4.0), ("avsd repair", "Cardiac", "major open", "3", 4.4),
    ("arterial switch", "Cardiac", "major open", "4", 5.6), ("norwood procedure", "Cardiac", "major open", "5", 6.5),
    ("cardiac exploration", "Cardiac", "major open", "NC", 3.0),
    ("gastroscopy", "GI/abdominal", "minor/diagnostic", "", 0.6), ("ng tube placement", "GI/abdominal", "catheter/interventional", "", 0.5),
    ("laparoscopic appendectomy", "GI/abdominal", "minimally invasive", "", 1.4), ("laparoscopic fundoplication", "GI/abdominal", "minimally invasive", "", 2.0),
    ("laparotomy", "GI/abdominal", "major open", "", 2.6), ("bowel resection", "GI/abdominal", "major open", "", 3.2),
    ("laparoscopic cholecystectomy", "Hepatobiliary", "minimally invasive", "", 4.2),
    ("kasai portoenterostomy", "Hepatobiliary", "major open", "", 2.1),
    ("craniotomy", "Neurosurgical", "major open", "", 3.8), ("vp shunt", "Neurosurgical", "major open", "", 1.5),
    ("endoscopic third ventriculostomy", "Neurosurgical", "minimally invasive", "", 1.2),
    ("central line placement", "Vascular access/line placement", "minor/diagnostic", "", 0.5),
    ("unknown procedure x", "Unclassified", "NA", "", 1.0),
]
# relative sampling weights per name (small groups: avsd/tof small, ETV tiny, cardiac exploration tiny)
WEIGHTS = {"tof repair": 0.6, "avsd repair": 0.5, "arterial switch": 0.7, "norwood procedure": 0.5, "cardiac exploration": 0.25,
           "endoscopic third ventriculostomy": 0.12, "unknown procedure x": 0.1, "ng tube placement": 0.5}


def build_map(names=SYN_NAMES, rng=None) -> pd.DataFrame:
    rng = rng or np.random.default_rng(1)
    rows = []
    for nm, bs, it, st, _ in names:
        if bs == "Unclassified":
            ig, sg, cg = "NA", "", "NA"
        elif bs == "Cardiac":
            sg = {"1": "G1", "2": "G1", "3": "G2", "4": "G3", "5": "G3", "NC": "NC"}[st]
            ig = {"minor/diagnostic": "I1", "catheter/interventional": "I1", "minimally invasive": "I2", "major open": "I3"}[it]
            cg = sg
        else:
            ig = {"minor/diagnostic": "I1", "catheter/interventional": "I1", "minimally invasive": "I2", "major open": "I3"}[it]
            sg, cg = "", ig
        rows.append({"raw_name": nm, "standardized_name": nm, "body_system": bs, "intensity": it, "I_group": ig,
                     "STAT": st, "STAT_group": sg, "complexity_group": cg, "review_flag": "Y" if rng.random() < 0.35 else "N",
                     "n_components": str(int(rng.integers(1, 4))), "STAT_group_crosses": ("Y" if rng.random() < 0.3 else "N") if bs == "Cardiac" else ""})
    return pd.DataFrame(rows)


def _ts(x: pd.Timestamp) -> str:
    return x.strftime(FMT)


def make_project(root: Path, n_subjects: int = 600, seed: int = 7, map_df: pd.DataFrame | None = None,
                 base_hours: dict | None = None, val_shift_hours: float = 0.0, shared_oper_every: int = 9,
                 true_dup_every: int = 41, harm_layout: str = "OUT") -> dict:
    """Writes the synthetic project under `root` and returns paths, hashes and expected counts for AuditConfig.

    shared_oper_every / true_dup_every: add a distinct surgery sharing an OPER_ID (or a possible true duplicate)
    every n-th surgery row (0 disables). harm_layout="OUT" writes encounter_id,t0,OUT,icu_time_hours (v1.5
    fixtures); harm_layout="real" writes the real header subset encounter_id,...,OUTTIME,...,t0,...,span_icu_hours
    with no OUT column."""
    rng = np.random.default_rng(seed)
    root = Path(root)
    (root / "data" / "interim").mkdir(parents=True, exist_ok=True)
    (root / "config").mkdir(parents=True, exist_ok=True)
    (root / "raw").mkdir(parents=True, exist_ok=True)
    cmap = map_df if map_df is not None else build_map(rng=np.random.default_rng(seed + 1))
    if base_hours is None:
        base_hours = {nm: h for nm, _, _, _, h in SYN_NAMES}
    names = list(cmap["raw_name"])
    w = np.array([WEIGHTS.get(n, 1.0) for n in names], float)
    w = w / w.sum()

    subjects = [str(100000 + i) for i in range(n_subjects - 1)] + [SENTINEL_SUBJECT]
    split_of = {}
    for s in subjects:
        u = rng.random()
        split_of[s] = "train" if u < 0.6 else ("val" if u < 0.8 else "test")
    enc_rows, harm_rows, adm_rows, surg_rows = [], [], [], []
    hadm_counter = 500000
    t_base = pd.Timestamp("2015-01-01 08:00:00")
    list_formats = ["json", "repr", "floatrepr", "delim", "plain"]
    for si, s in enumerate(subjects):
        n_enc = 1 if rng.random() < 0.8 else 2
        same_hadm = n_enc == 2 and rng.random() < 0.5
        t = t_base + pd.Timedelta(days=int(rng.integers(0, 2000)), hours=int(rng.integers(0, 24)))
        encs = []
        hadm_first = None
        for e in range(n_enc):
            if same_hadm and e == 1:
                hadms = [hadm_first]
            else:
                hadm_counter += 1
                hadms = [str(hadm_counter)]
                if rng.random() < 0.08:  # multi-HADM encounter
                    hadm_counter += 1
                    hadms.append(str(hadm_counter))
            if hadm_first is None:
                hadm_first = hadms[0]
            t0 = t
            out = t0 + pd.Timedelta(hours=float(rng.uniform(30, 400)))
            eid = f"{s}_{e + 1}"
            fmt = list_formats[(si + e) % len(list_formats)] if len(hadms) > 1 else ("plain" if (si % 3) else "json")
            if fmt == "json":
                hl = "[" + ", ".join(hadms) + "]"
            elif fmt == "repr":
                hl = str(hadms)
            elif fmt == "floatrepr":
                hl = "[" + ", ".join(f"{h}.0" for h in hadms) + "]"
            elif fmt == "delim":
                hl = ";".join(hadms)
            else:
                hl = hadms[0] if len(hadms) == 1 else ",".join(hadms)
            enc_rows.append({"encounter_id": eid, "SUBJECT_ID": s, "hadm_id_list": hl, "HADM_ID": hadms[0],
                             "INTIME": _ts(t0), "OUTTIME": _ts(out), "Died": "0"})
            span_h = f"{(out - t0).total_seconds() / 3600:.4f}"
            if harm_layout == "real":
                harm_rows.append({"encounter_id": eid, "SUBJECT_ID": s, "INTIME": _ts(t0), "OUTTIME": _ts(out), "Died": "0",
                                  "t0": _ts(t0), "icu_time_hours": span_h, "span_icu_hours": span_h})
            else:
                harm_rows.append({"encounter_id": eid, "t0": _ts(t0), "OUT": _ts(out), "icu_time_hours": f"{(out - t0).total_seconds() / 3600:.2f}"})
            for h in hadms:
                adm_rows.append({"SUBJECT_ID": s, "HADM_ID": h, "ADMITTIME": _ts(t0 - pd.Timedelta(days=6))})
            encs.append((eid, hadms, t0, out))
            t = out + pd.Timedelta(hours=float(rng.uniform(48, 500)))
        # surgeries
        oper = {}
        for _ in range(int(rng.integers(1, 4))):
            eid, hadms, t0, out = encs[int(rng.integers(0, len(encs)))]
            h = hadms[int(rng.integers(0, len(hadms)))]
            u = rng.random()
            if u < 0.62:
                B = t0 + (out - t0) * float(rng.uniform(0, 1))
            elif u < 0.74:
                B = t0 - pd.Timedelta(hours=float(rng.uniform(0.5, 23.5)))
            elif u < 0.84:
                B = t0 - pd.Timedelta(hours=float(rng.uniform(30, 120)))
            elif u < 0.88:
                B = encs[-1][3] + pd.Timedelta(hours=float(rng.uniform(2, 40)))  # after last OUT
                eid, hadms, t0, out = encs[-1]
                h = hadms[0]
            elif u < 0.91:
                hadm_counter += 1
                h = str(hadm_counter)  # other HADM, not in any encounter list
                adm_rows.append({"SUBJECT_ID": s, "HADM_ID": h, "ADMITTIME": _ts(B0 := t0 - pd.Timedelta(days=200))})
                B = B0 + pd.Timedelta(days=1)
            else:
                # v1.7: begins shortly after t0; short surgeries end by t0 + 6 h -> recent_layer_a (Layer A only)
                B = t0 + pd.Timedelta(hours=float(rng.uniform(0, 3)))
            nm = names[int(rng.choice(len(names), p=w))]
            hours = float(base_hours[nm]) * float(np.exp(rng.normal(0, 0.35)))
            if split_of[s] == "val":
                hours += val_shift_hours
            E = B + pd.Timedelta(hours=hours)
            AS = B - pd.Timedelta(minutes=float(rng.uniform(10, 50)))
            AE = E + pd.Timedelta(minutes=float(rng.uniform(5, 30)))
            oper[h] = oper.get(h, 0) + 1
            surg_rows.append({"SUBJECT_ID": s, "HADM_ID": h, "OPER_ID": str(oper[h]), "B": B, "E": E, "AS": AS, "AE": AE, "name": nm})
            # v1.6: as in real SURGERY_INFO, OPER_ID is not unique within an admission. Some admissions carry a second,
            # DISTINCT surgery (different times, name, anaesthesia window) under the same OPER_ID.
            if shared_oper_every and len(surg_rows) % shared_oper_every == 0:
                nm2 = names[int(rng.choice(len(names), p=w))]
                B2 = B + pd.Timedelta(hours=float(rng.uniform(3, 30)))
                E2 = B2 + pd.Timedelta(hours=float(base_hours[nm2]) * float(np.exp(rng.normal(0, 0.35))))
                surg_rows.append({"SUBJECT_ID": s, "HADM_ID": h, "OPER_ID": str(oper[h]), "B": B2, "E": E2,
                                  "AS": B2 - pd.Timedelta(minutes=25), "AE": E2 + pd.Timedelta(minutes=15), "name": nm2})
            # possible true duplicate: same SUBJECT/HADM/begin/name re-entered under the same OPER_ID (kept, counted)
            if true_dup_every and len(surg_rows) % true_dup_every == 0:
                surg_rows.append({"SUBJECT_ID": s, "HADM_ID": h, "OPER_ID": str(oper[h]), "B": B, "E": E + pd.Timedelta(minutes=5),
                                  "AS": AS, "AE": AE + pd.Timedelta(minutes=5), "name": nm})
    # not-in-cohort subjects
    for k in range(25):
        s = str(900000 + k)
        B = t_base + pd.Timedelta(days=int(rng.integers(0, 2000)))
        nm = names[int(rng.choice(len(names), p=w))]
        surg_rows.append({"SUBJECT_ID": s, "HADM_ID": str(800000 + k), "OPER_ID": "1", "B": B,
                          "E": B + pd.Timedelta(hours=base_hours[nm]), "AS": B - pd.Timedelta(minutes=20),
                          "AE": B + pd.Timedelta(hours=base_hours[nm], minutes=10), "name": nm})
    surg = pd.DataFrame(surg_rows)
    n = len(surg)
    out_rows = []
    for i, r in surg.iterrows():
        B, E, AS, AE = _ts(r["B"]), _ts(r["E"]), _ts(r["AS"]), _ts(r["AE"])
        if i % 97 == 3:
            B = ""                          # invalid_time
        elif i % 89 == 5:
            E = ""                          # fallback to anesthesia window
        elif i % 211 == 7:
            E, AS = "", ""                  # missing_window
        elif i % 83 == 9:
            E = B                           # zero duration
        elif i % 157 == 11:
            E = _ts(r["B"] - pd.Timedelta(minutes=30))  # B > E
        elif i % 53 == 13 and B:
            B = _ts(r["B"])                 # keep
        elif i % 101 == 17:
            AE = ""                         # v1.7: timing falls back to SURGERY_END_TIME (counted)
        elif i % 193 == 19:
            AE, E = "", ""                  # v1.7: no anaesthesia end and no surgery end -> invalid_time
        sid = r["SUBJECT_ID"] + (".0" if i % 7 == 0 else "")      # float-like IDs
        hid = r["HADM_ID"] + (".0" if i % 5 == 0 else "")
        nm = r["name"].upper() if i % 11 == 0 else ("  " + r["name"].replace(" ", "  ") + " ")  # case / whitespace noise
        out_rows.append({"ROW_ID": str(i + 1), "SUBJECT_ID": sid, "HADM_ID": hid, "VISIT_ID": "1", "OPER_ID": r["OPER_ID"],
                         "ANES_START_TIME": AS, "ANES_END_TIME": AE, "SURGERY_BEGIN_TIME": B, "SURGERY_END_TIME": E,
                         "SURGERY_NAME": nm, "ANES_METHOD": "general", "SURGERY_POSITION": "supine"})
    si_df = pd.DataFrame(out_rows)
    # one straddling surgery: begins before t0, ends after it
    paths = {
        "canonical": root / "data" / "interim" / "canonical_icu_encounters.csv",
        "harmonization": root / "data" / "interim" / "encounter_harmonization_fields.csv",
        "complexity_map": root / "config" / "surgery_complexity_map_v1.csv",
        "split": root / "raw" / "subject_splits.csv",
        "surgery_info": root / "raw" / "SURGERY_INFO.csv",
        "admissions": root / "raw" / "ADMISSIONS.csv",
    }
    pd.DataFrame(enc_rows).to_csv(paths["canonical"], index=False)
    pd.DataFrame(harm_rows).to_csv(paths["harmonization"], index=False)
    cmap.to_csv(paths["complexity_map"], index=False)
    pd.DataFrame({"SUBJECT_ID": subjects, "split": [split_of[s] for s in subjects]}).to_csv(paths["split"], index=False)
    si_df.to_csv(paths["surgery_info"], index=False)
    pd.DataFrame(adm_rows).drop_duplicates(["SUBJECT_ID", "HADM_ID"]).to_csv(paths["admissions"], index=False)
    hashes = {k: hashlib.sha256(paths[k].read_bytes()).hexdigest() for k in ("canonical", "harmonization", "complexity_map", "split")}
    enc_df = pd.DataFrame(enc_rows)
    sp = enc_df["SUBJECT_ID"].map(split_of).map({"train": "train", "val": "validation", "test": "test"})
    return {
        "paths": paths, "hashes": hashes, "n_surgeries": n, "n_encounters": len(enc_df), "n_subjects": enc_df["SUBJECT_ID"].nunique(),
        "split_encounters": {k: int((sp == k).sum()) for k in ("train", "validation", "test")}, "map_rows": len(cmap),
        "missing": {c: int((si_df[c] == "").sum()) for c in si_df.columns},
        "shared_oper_every": shared_oper_every, "true_dup_every": true_dup_every,
    }


def make_config(proj: dict, out_dir: Path, bootstrap_reps: int = 50, **kw):
    import audit_surgery_complexity_duration as A
    p = proj["paths"]
    return A.AuditConfig(
        canonical_path=p["canonical"], harmonization_path=p["harmonization"], map_path=p["complexity_map"], split_path=p["split"],
        surgery_info_path=p["surgery_info"], admissions_path=p["admissions"], out_dir=Path(out_dir),
        pinned_sha256=dict(proj["hashes"]), expected_surgery_rows=proj["n_surgeries"], expected_encounters=proj["n_encounters"],
        expected_subjects=proj["n_subjects"], expected_split_encounters=proj["split_encounters"], expected_map_rows=proj["map_rows"],
        expected_missing=proj["missing"], bootstrap_reps=bootstrap_reps, **kw)
