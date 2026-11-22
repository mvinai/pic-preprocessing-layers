#!/usr/bin/env python3
"""Generate one synthetic PIC-like dataset in three stages for SVD algorithm work.

Standalone generator. Does not run the Stage 2 audit and does not modify the
audit script, tests, or fixture. Imports v1 map writers / schema constants from
the fixture module and (read-only) conversion, merge, and age helpers from the
audit module.

Planted low-rank structure is generated INSIDE the raw writer: each encounter
has three independent N(0, 1) latent factors, and event values are drawn around
means shifted by those scores (clipped to clinical ranges, plus observation
noise). Stage 2 (normalized) is derived from raw; stage 3 (features) is derived
from normalized, so the structure appears in the encounter-by-feature matrix
naturally. Factor scores are saved only in 4_answer_key/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
AUDIT_DIR = Path(__file__).resolve().parents[2] / "scripts" / "feature_sparsity_frequency_audit"
sys.path.insert(0, str(AUDIT_DIR))

import feature_sparsity_frequency_audit as M  # noqa: E402
import synthetic_stage2_data as S  # noqa: E402

SEED = 20261115
N_SUBJECTS_DEFAULT = 15000
VITALS_STEP_DEFAULT = 4.0
EVENT_HORIZON_H = 36.0
INVALID_TARGET_FRAC = 0.006
WINDOW = (0.0, 24.0)

# Raw-scale generative loadings (mean shift per +1 SD of the named factor).
# Age also shifts HR / RR / SBP / DBP / weight through age-conditional norms.
BETA = {
    "lactate_sev": 1.05,          # mmol/L
    "creatinine_sev": 28.0,       # umol/L
    "hr_sev": 14.0,               # bpm
    "map_sev": -8.0,              # applied via SBP -8, DBP -5
    "sbp_sev": -8.0,
    "dbp_sev": -5.0,
    "vaso_logit_intercept": -1.15,
    "vaso_logit_sev": 1.55,
    "spo2_resp": -3.8,            # percentage points
    "po2_resp": -22.0,            # mmHg
    "fio2_resp": 16.0,            # %
    "rr_resp": 7.0,               # /min
    "weight_log_age": 0.10,       # residual log-weight on z_age, beyond the age curve
}
NOISE_SD = {
    "hr": 5.0, "rr": 2.5, "sbp": 6.0, "dbp": 5.0, "spo2": 1.4,
    "temp": 0.28, "lactate": 0.22, "creatinine": 8.0, "po2": 8.0, "fio2": 3.5,
    "weight_log": 0.06, "glucose": 0.8, "sodium": 3.0, "hemoglobin": 8.0,
    "magnesium": 0.12, "calcium": 0.15, "urine": 40.0, "intake": 50.0,
}

SVD_COLS = [
    "age_days", "weight_kg",
    "heart_rate_mean_0_24h", "respiratory_rate_mean_0_24h",
    "sbp_mean_0_24h", "dbp_mean_0_24h", "map_mean_0_24h",
    "spo2_mean_0_24h", "fio2_mean_0_24h", "po2_mean_0_24h",
    "pf_ratio_mean_0_24h", "sf_ratio_mean_0_24h",
    "lactate_mean_0_24h", "creatinine_mean_0_24h",
    "vasoactive_exposed_0_24h", "vasoactive_max_concurrent_0_24h",
    "psofa_resp_0_24h",
    "temperature_mean_0_24h", "glucose_bedside_mean_0_24h", "sodium_mean_0_24h",
    "hemoglobin_mean_0_24h", "magnesium_mean_0_24h", "calcium_total_mean_0_24h",
    "urine_output_mean_0_24h", "fluid_intake_mean_0_24h",
]
SIGNAL_COLS = set(SVD_COLS[:17])
NOISE_COLS = set(SVD_COLS[17:])
FACTOR_NAMES = ("z_severity", "z_respiratory", "z_age")

FMT = S.FMT


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ts(x: pd.Timestamp) -> str:
    return x.strftime(FMT)


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))


def age_years_from_z(z: np.ndarray) -> np.ndarray:
    return np.clip(9.0 + 6.0 * np.asarray(z, dtype=float), 0.05, 21.0)


def expected_weight_kg(age_y: np.ndarray) -> np.ndarray:
    a = np.asarray(age_y, dtype=float)
    return np.clip(3.3 + 2.4 * a + 0.12 * np.power(np.maximum(a, 0.0), 1.6), 2.0, 90.0)


def hr_norm(age_y):
    return np.clip(155.0 - 4.0 * age_y, 70.0, 180.0)


def rr_norm(age_y):
    return np.clip(42.0 - 1.35 * age_y, 14.0, 70.0)


def sbp_norm(age_y):
    return np.clip(75.0 + 2.3 * age_y, 50.0, 130.0)


def dbp_norm(age_y):
    return np.clip(45.0 + 1.35 * age_y, 25.0, 85.0)


def unit_from_age(age_y: float, rng: np.random.Generator) -> str:
    if age_y < 0.08:
        return "NICU" if rng.random() < 0.85 else "PICU"
    if age_y < 1.0:
        return rng.choice(["NICU", "PICU"], p=[0.35, 0.65])
    return rng.choice(["PICU", "SICU", "CICU"], p=[0.62, 0.20, 0.18])


# =============================================================================
# Stage 1: raw
# =============================================================================
def generate_raw(out_raw: Path, n_subjects: int, seed: int, vitals_step: float) -> dict:
    """Write PIC-schema raw tables with factor-shifted event values. Returns metadata + factor frame."""
    rng = np.random.default_rng(seed)
    out_raw.mkdir(parents=True, exist_ok=True)
    scratch = out_raw.parent / "_scratch_project"
    maps_dir = scratch / "maps"
    interim = scratch / "interim"
    cfg = scratch / "config"
    for p in (maps_dir, interim, cfg):
        p.mkdir(parents=True, exist_ok=True)

    map_sha, map_rows = S.write_maps(maps_dir, extra_labs=0)
    cmap = S.build_cmap()
    cmap.to_csv(out_raw / "surgery_complexity_map_v1.csv", index=False)

    subjects = [str(300000 + i) for i in range(n_subjects)]
    # z_age is patient-level (DOB is a patient attribute). Severity and respiratory are per encounter.
    z_age_of = {s: float(rng.normal()) for s in subjects}
    split_of = {}
    for s in subjects:
        u = rng.random()
        split_of[s] = "train" if u < 0.60 else ("validation" if u < 0.80 else "test")

    base = pd.Timestamp("2015-03-01 00:00:00")
    encs, stays, adms, pats = [], [], [], []
    hadm, icu = 700000, 900000
    H = lambda t0, h: (t0 + pd.Timedelta(hours=float(h))).floor("s")

    for si, s in enumerate(subjects):
        z_a = z_age_of[s]
        age_y = float(age_years_from_z(z_a))
        sex = "M" if rng.random() < 0.55 else "F"
        n_enc = 1 if rng.random() < 0.85 else 2
        t = base + pd.Timedelta(days=int(rng.integers(0, 900)),
                                hours=int(rng.integers(0, 24)),
                                minutes=int(rng.integers(0, 60)))
        dob = (t - pd.Timedelta(days=float(age_y * 365.25))).floor("s")
        pats.append({"ROW_ID": si, "SUBJECT_ID": s, "GENDER": sex, "DOB": ts(dob), "DOD": ""})
        for e in range(n_enc):
            hadm += 1
            hl = [str(hadm)]
            admit = t - pd.Timedelta(hours=float(rng.choice([2.0, 30.0, 100.0, 200.0])))
            dx = rng.choice(
                ["S06.001", "Q28.200", "Q28.000", "J18.900", "G40.100", "A39.000", "K35.800", "", "XYZ", "J18.9,G40.1"],
                p=[0.10, 0.05, 0.05, 0.30, 0.10, 0.05, 0.15, 0.08, 0.04, 0.08],
            )
            u = rng.random()
            ed = "" if u < 0.3 else (ts(admit) if u < 0.5 else ts(admit - pd.Timedelta(hours=1)))
            adms.append({"ROW_ID": len(adms), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ADMITTIME": ts(admit),
                         "DISCHTIME": ts(t + pd.Timedelta(days=30)), "EDREGTIME": ed, "ICD10_CODE_CN": dx})
            t0 = t
            bucket = rng.choice([0, 1, 2, 3], p=[0.25, 0.40, 0.25, 0.10])
            span = float({0: rng.uniform(20, 36), 1: rng.uniform(36, 72), 2: rng.uniform(72, 168), 3: rng.uniform(168, 280)}[bucket])
            out = (t0 + pd.Timedelta(hours=span)).floor("s")
            unit = unit_from_age(age_y + (t0 - dob).total_seconds() / (365.25 * 86400) - age_y, rng)
            # age at this encounter ≈ age_y plus time since first INTIME; for the first encounter it is age_y
            age_enc = age_y + (t0 - dob).total_seconds() / (365.25 * 86400.0)
            unit = unit_from_age(age_enc, rng)
            z_sev = float(rng.normal())
            z_resp = float(rng.normal())
            died = int(rng.random() < sigmoid(-2.4 + 0.55 * z_sev))  # weak, not used in the SVD matrix
            stays.append({"ROW_ID": len(stays), "SUBJECT_ID": s, "HADM_ID": str(hadm), "ICUSTAY_ID": str(icu),
                          "FIRST_CAREUNIT": unit, "INTIME": ts(t0), "OUTTIME": ts(out)})
            icu += 1
            encs.append({
                "encounter_id": f"E{s}_{e + 1}", "SUBJECT_ID": s, "HADM_ID": str(hadm), "hadm_list": hl,
                "t0": t0, "OUT": out, "split": split_of[s], "died": died, "dx": dx, "unit": unit,
                "z_severity": z_sev, "z_respiratory": z_resp, "z_age": z_a,
                "age_years_planted": age_enc, "sex": sex, "span_h": span,
            })
            t = (out + pd.Timedelta(hours=float(rng.uniform(40, 400)))).floor("s")

    enc = pd.DataFrame(encs)
    pats_df = pd.DataFrame(pats)

    # Invalid DOB < 1% of encounters, mix of missing / unparseable / negative bins.
    # Only singleton subjects so one DOB change = one encounter.
    counts = enc.groupby("SUBJECT_ID").size()
    singletons = [s for s in counts[counts == 1].index.tolist() if s in set(pats_df["SUBJECT_ID"])]
    n_enc = len(enc)
    cap = max(0, int(0.01 * n_enc) - 1)
    n_inv = int(round(INVALID_TARGET_FRAC * n_enc))
    n_inv = min(max(n_inv, min(10, cap)), cap)
    # 5 types, as even as possible
    types = (["missing", "unparseable", "neg_lt1d", "neg_1to30d", "neg_gt30d"] * ((n_inv // 5) + 1))[:n_inv]
    rng.shuffle(types)
    pick = rng.choice(np.array(singletons), size=n_inv, replace=False)
    pats_df = pats_df.set_index("SUBJECT_ID")
    enc_by_s = enc.set_index("SUBJECT_ID")
    unparseable_strings = ["not-a-date", "9999-99-99", "FOOBAR", "15/31/2010", "----"]
    for s, typ in zip(pick, types):
        t0 = enc_by_s.loc[s, "t0"]
        if typ == "missing":
            pats_df.at[s, "DOB"] = ""
        elif typ == "unparseable":
            pats_df.at[s, "DOB"] = unparseable_strings[int(rng.integers(0, len(unparseable_strings)))]
        elif typ == "neg_lt1d":
            pats_df.at[s, "DOB"] = ts(t0 + pd.Timedelta(hours=float(rng.uniform(2, 18))))
        elif typ == "neg_1to30d":
            pats_df.at[s, "DOB"] = ts(t0 + pd.Timedelta(days=float(rng.uniform(2, 25))))
        else:
            pats_df.at[s, "DOB"] = ts(t0 + pd.Timedelta(days=float(rng.uniform(40, 200))))
    pats_df = pats_df.reset_index()
    inv_subjects = set(pick)

    # Event generation around factor-shifted means.
    chart, lab, rx, micro, surg = [], [], [], [], []
    vaso_names = [d for d in S.DRUGS if d[2] == "vasoactive"]
    other_drugs = [d for d in S.DRUGS if d[2] != "vasoactive"]

    for r in enc.itertuples():
        sp = float(r.span_h)
        h = r.HADM_ID
        t0 = r.t0
        z_sev, z_resp, z_a = float(r.z_severity), float(r.z_respiratory), float(r.z_age)
        age_y = float(r.age_years_planted)
        horizon = min(sp, EVENT_HORIZON_H)

        hr_mu = float(hr_norm(age_y) + BETA["hr_sev"] * z_sev)
        rr_mu = float(rr_norm(age_y) + BETA["rr_resp"] * z_resp)
        sbp_mu = float(sbp_norm(age_y) + BETA["sbp_sev"] * z_sev)
        dbp_mu = float(dbp_norm(age_y) + BETA["dbp_sev"] * z_sev)
        spo2_mu = 97.2 + BETA["spo2_resp"] * z_resp
        fio2_mu = 30.0 + BETA["fio2_resp"] * z_resp
        po2_mu = 95.0 + BETA["po2_resp"] * z_resp
        lact_mu = 1.15 + BETA["lactate_sev"] * z_sev
        creat_mu = 42.0 + BETA["creatinine_sev"] * z_sev
        # encounter-level intercepts for noise channels (independent of z)
        temp_mu = 37.0 + rng.normal(0, 0.35)
        glu_mu = 5.4 + rng.normal(0, 0.6)
        na_mu = 138.0 + rng.normal(0, 2.0)
        hb_mu = 115.0 + rng.normal(0, 10.0)
        mg_mu = 0.85 + rng.normal(0, 0.08)
        ca_mu = 2.30 + rng.normal(0, 0.12)
        urine_mu = 40.0 + rng.normal(0, 15.0)
        intake_mu = 60.0 + rng.normal(0, 20.0)

        wt = float(expected_weight_kg(age_y) * math.exp(BETA["weight_log_age"] * z_a + rng.normal(0, NOISE_SD["weight_log"])))
        wt = float(np.clip(wt, 0.5, 140.0))
        chart.append([r.SUBJECT_ID, h, "1014", ts(H(t0, 1.0)), f"{wt:.2f}", "kg"])

        times = np.arange(-2.0, horizon, float(vitals_step))
        for k in times:
            hr = float(np.clip(hr_mu + rng.normal(0, NOISE_SD["hr"]), 30, 300))
            chart.append([r.SUBJECT_ID, h, "1003", ts(H(t0, k)), f"{hr:.0f}", "bpm"])
            if rng.random() < 0.18:
                hr2 = float(np.clip(hr_mu + rng.normal(0, NOISE_SD["hr"]), 30, 300))
                chart.append([r.SUBJECT_ID, h, "1002", ts(H(t0, k)), f"{hr2:.0f}", "bpm"])
            if rng.random() < 0.85:
                rr = float(np.clip(rr_mu + rng.normal(0, NOISE_SD["rr"]), 5, 120))
                chart.append([r.SUBJECT_ID, h, "1004", ts(H(t0, k)), f"{rr:.0f}", "bpm"])
            if rng.random() < 0.90:
                sv = float(np.clip(sbp_mu + rng.normal(0, NOISE_SD["sbp"]), 30, 250))
                dv = float(np.clip(dbp_mu + rng.normal(0, NOISE_SD["dbp"]), 10, 200))
                if dv >= sv - 2:
                    dv = max(10.0, sv - 8.0)
                chart.append([r.SUBJECT_ID, h, "1016", ts(H(t0, k)), f"{sv:.0f}", "mmHg"])
                chart.append([r.SUBJECT_ID, h, "1015", ts(H(t0, k)), f"{dv:.0f}", "mmHg"])
            if rng.random() < 0.90:
                spo = float(np.clip(spo2_mu + rng.normal(0, NOISE_SD["spo2"]), 30, 100))
                chart.append([r.SUBJECT_ID, h, "1006", ts(H(t0, k)), f"{spo:.0f}", "%"])
            if rng.random() < 0.55:
                tp = float(np.clip(temp_mu + rng.normal(0, NOISE_SD["temp"]), 30, 43))
                chart.append([r.SUBJECT_ID, h, "1001", ts(H(t0, k)), f"{tp:.1f}", "°C"])
            if rng.random() < 0.45:
                chart.append([r.SUBJECT_ID, h, "1008", ts(H(t0, k)), f"{float(np.clip(urine_mu + rng.normal(0, NOISE_SD['urine']), 0, 3000)):.0f}", "ml"])
            if rng.random() < 0.45:
                chart.append([r.SUBJECT_ID, h, "1009", ts(H(t0, k + 0.25)), f"{float(np.clip(intake_mu + rng.normal(0, NOISE_SD['intake']), 0, 5000)):.0f}", "ml"])
            if rng.random() < 0.25:
                glu = float(np.clip(glu_mu + rng.normal(0, NOISE_SD["glucose"]), 0.5, 50))
                chart.append([r.SUBJECT_ID, h, "1011", ts(H(t0, k)), f"{glu:.1f}", "mmol/L"])

        # Labs in a window that lands inside [0, 24) after the +1 h lab-release lag used in normalization.
        n_chem = int(rng.integers(2, 5))
        for k in rng.uniform(0.2, min(horizon, 22.5), size=n_chem):
            cr = float(np.clip(creat_mu + rng.normal(0, NOISE_SD["creatinine"]), 5, 2000))
            item = "5041" if rng.random() < 0.75 else "5032"
            lab.append([r.SUBJECT_ID, h, item, ts(H(t0, k)), f"{cr:.0f}", "μmol/L"])
            na = float(np.clip(na_mu + rng.normal(0, NOISE_SD["sodium"]), 100, 190))
            lab.append([r.SUBJECT_ID, h, "5230", ts(H(t0, k)), f"{na:.0f}", "mmol/L"])
        n_lac = int(rng.integers(2, 5))
        for k in rng.uniform(0.2, min(horizon, 22.5), size=n_lac):
            lac = float(np.clip(lact_mu + rng.normal(0, NOISE_SD["lactate"]), 0.1, 30))
            lab.append([r.SUBJECT_ID, h, "5227", ts(H(t0, k)), f"{lac:.2f}", "mmol/L"])
        if rng.random() < 0.7:
            k = float(rng.uniform(0.2, min(horizon, 22.5)))
            hb = float(np.clip(hb_mu + rng.normal(0, NOISE_SD["hemoglobin"]), 20, 250))
            lab.append([r.SUBJECT_ID, h, "5099", ts(H(t0, k)), f"{hb:.0f}", "g/L"])
            if rng.random() < 0.4:
                lab.append([r.SUBJECT_ID, h, "5257", ts(H(t0, k + 0.1)), f"{hb / 10.0:.1f}", "g/dL"])
        if rng.random() < 0.5:
            k = float(rng.uniform(0.2, min(horizon, 22.5)))
            lab.append([r.SUBJECT_ID, h, "5059", ts(H(t0, k)), f"{float(np.clip(mg_mu + rng.normal(0, NOISE_SD['magnesium']), 0.1, 5)):.2f}", "mmol/L"])
        if rng.random() < 0.5:
            k = float(rng.uniform(0.2, min(horizon, 22.5)))
            ca = float(np.clip(ca_mu + rng.normal(0, NOISE_SD["calcium"]), 0.5, 5))
            lab.append([r.SUBJECT_ID, h, "5215" if rng.random() < 0.7 else "5034", ts(H(t0, k)), f"{ca:.2f}", "mmol/L"])

        # Blood gas / FiO2 / pO2 / S/F partners. Mix of P/F draws and S/F (FiO2 + nearby SpO2) draws.
        n_gas = int(rng.integers(2, 5))
        for k in rng.uniform(0.2, min(horizon, 22.5), size=n_gas):
            fi = float(np.clip(fio2_mu + rng.normal(0, NOISE_SD["fio2"]), 21, 100))
            lab.append([r.SUBJECT_ID, h, "5222", ts(H(t0, k)), f"{fi:.0f}", "%"])
            if rng.random() < 0.62:
                po = float(np.clip(po2_mu + rng.normal(0, NOISE_SD["po2"]), 10, 700))
                po_item = "5239" if rng.random() < 0.7 else "5244"
                lab.append([r.SUBJECT_ID, h, po_item, ts(H(t0, k)), f"{po:.0f}", "mmHg"])
                # analyser P/F cross-check (noisy)
                lab.append([r.SUBJECT_ID, h, "5243", ts(H(t0, k)), f"{po / (fi / 100.0) + rng.normal(0, 8):.0f}", "mmHg"])
            else:
                spo = float(np.clip(spo2_mu + rng.normal(0, NOISE_SD["spo2"]), 30, 100))
                # keep many S/F-eligible (SpO2 <= 97) especially when z_resp is high
                if spo > 97 and z_resp > 0:
                    spo = float(rng.uniform(90, 97))
                delta_min = int(rng.choice([-4, -2, 0, 2, 3]))
                chart.append([r.SUBJECT_ID, h, "1006", ts(H(t0, k) + pd.Timedelta(minutes=delta_min)), f"{spo:.0f}", "%"])

        # Vasoactive exposure / concurrent count driven by severity; other drugs are noise.
        p_vaso = float(sigmoid(BETA["vaso_logit_intercept"] + BETA["vaso_logit_sev"] * z_sev))
        if rng.random() < p_vaso:
            # more agents and more overlap as severity rises
            n_agents = 1 + int(rng.random() < sigmoid(-0.6 + 1.1 * z_sev)) + int(rng.random() < sigmoid(-1.4 + 1.1 * z_sev))
            n_agents = int(np.clip(n_agents, 1, len(vaso_names)))
            chosen = [vaso_names[i] for i in rng.choice(len(vaso_names), size=n_agents, replace=False)]
            # shared window so concurrent count rises with n_agents
            center = float(rng.uniform(2.0, 14.0))
            for dn, de, g, sg, cr in chosen:
                start_h = center + float(rng.uniform(-3.0, 3.0))
                dur = float(rng.uniform(8.0, 20.0))
                stt = H(t0, start_h)
                ent = H(t0, start_h + dur)
                rx.append([r.SUBJECT_ID, h, ts(stt), ts(ent), dn, de, "mg"])
        for dn, de, g, sg, cr in other_drugs:
            if rng.random() < 0.28:
                k = float(rng.uniform(-6, min(horizon, 20)))
                stt = H(t0, k)
                if rng.random() < 0.12:
                    start = stt.strftime("%Y-%m-%d")
                else:
                    start = ts(stt)
                v = rng.random()
                if v < 0.70:
                    end = ""
                elif v < 0.95:
                    end = ts(H(t0, k + float(rng.uniform(4, 30))))
                else:
                    end = ts(H(t0, k - 1.0))
                rx.append([r.SUBJECT_ID, h, start, end, dn, de, "mg"])

        if rng.random() < 0.25:
            k = float(rng.uniform(-5, min(horizon, 20)))
            tt = ts(H(t0, k))
            sid, sdesc, _, _, _ = S.SPECIMENS[int(rng.integers(0, len(S.SPECIMENS)))]
            oid, oname, _ = S.ORGS[int(rng.integers(0, len(S.ORGS)))]
            micro.append([r.SUBJECT_ID, h, tt, sid, sdesc, oid, oname, "OXACILLIN", "S"])

        if rng.random() < 0.18:
            nm = S.CMAP_NAMES[int(rng.integers(0, len(S.CMAP_NAMES)))][0]
            b_ = H(t0, float(rng.choice([-8.0, -3.0, 6.0, -40.0])))
            e_ = b_ + pd.Timedelta(hours=float(rng.uniform(0.8, 3.0)))
            surg.append([r.SUBJECT_ID, h, h, str(len(surg)), ts(b_ - pd.Timedelta(minutes=20)), ts(e_ + pd.Timedelta(minutes=10)),
                         ts(b_), ts(e_), nm, "general", "supine"])

    # Write raw tables and cohort / split / map companions into 1_raw.
    canon = pd.DataFrame({
        "encounter_id": enc["encounter_id"], "SUBJECT_ID": enc["SUBJECT_ID"],
        "HADM_ID": enc["HADM_ID"].map(lambda x: f"{x}.0"),
        "hadm_id_list": ["[" + ", ".join(h) + "]" for h in enc["hadm_list"]],
        "INTIME": enc["t0"].map(ts), "OUTTIME": enc["OUT"].map(ts), "Died": enc["died"],
    })
    harm = pd.DataFrame({
        "encounter_id": enc["encounter_id"], "SUBJECT_ID": enc["SUBJECT_ID"],
        "t0": enc["t0"].map(ts), "INTIME": enc["t0"].map(ts), "OUTTIME": enc["OUT"].map(ts),
        "span_icu_hours": ((enc["OUT"] - enc["t0"]).dt.total_seconds() / 3600).round(4),
        "gap_count_between_stays": 0,
    })
    splits = pd.DataFrame({"SUBJECT_ID": subjects, "split": [split_of[s] for s in subjects]})

    pats_df.to_csv(out_raw / "PATIENTS.csv", index=False)
    pd.DataFrame(adms).to_csv(out_raw / "ADMISSIONS.csv", index=False)
    pd.DataFrame(stays).to_csv(out_raw / "ICUSTAYS.csv", index=False)
    canon.to_csv(out_raw / "canonical_icu_encounters.csv", index=False)
    harm.to_csv(out_raw / "encounter_harmonization_fields.csv", index=False)
    splits.to_csv(out_raw / "subject_splits.csv", index=False)

    ch = pd.DataFrame(chart, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"])
    ch["VALUENUM"] = pd.to_numeric(ch["VALUE"], errors="coerce").map(lambda v: "" if pd.isna(v) else f"{v:g}")
    ch["ICUSTAY_ID"] = ""
    ch.to_csv(out_raw / "CHARTEVENTS.csv", index=False)
    lb = pd.DataFrame(lab, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"])
    lb["VALUENUM"] = pd.to_numeric(lb["VALUE"], errors="coerce").map(lambda v: "" if pd.isna(v) else f"{v:g}")
    lb.to_csv(out_raw / "LABEVENTS.csv", index=False)
    pd.DataFrame(rx, columns=["SUBJECT_ID", "HADM_ID", "STARTDATE", "ENDDATE", "DRUG_NAME", "DRUG_NAME_EN", "DOSE_UNIT_RX"]).to_csv(
        out_raw / "PRESCRIPTIONS.csv", index=False)
    pd.DataFrame(micro, columns=["SUBJECT_ID", "HADM_ID", "CHARTTIME", "SPEC_ITEMID", "SPEC_TYPE_DESC",
                                 "ORG_ITEMID", "ORG_NAME", "AB_NAME", "INTERPRETATION"]).to_csv(
        out_raw / "MICROBIOLOGYEVENTS.csv", index=False)
    if surg:
        sd = pd.DataFrame(surg, columns=["SUBJECT_ID", "HADM_ID", "VISIT_ID", "ROW_ID", "ANES_START_TIME", "ANES_END_TIME",
                                         "SURGERY_BEGIN_TIME", "SURGERY_END_TIME", "SURGERY_NAME", "ANES_METHOD", "SURGERY_POSITION"])
        sd = sd.assign(OPER_ID=lambda d: d["ROW_ID"])[
            ["ROW_ID", "SUBJECT_ID", "HADM_ID", "VISIT_ID", "OPER_ID", "ANES_START_TIME", "ANES_END_TIME",
             "SURGERY_BEGIN_TIME", "SURGERY_END_TIME", "SURGERY_NAME", "ANES_METHOD", "SURGERY_POSITION"]]
    else:
        sd = pd.DataFrame(columns=["ROW_ID", "SUBJECT_ID", "HADM_ID", "VISIT_ID", "OPER_ID", "ANES_START_TIME", "ANES_END_TIME",
                                   "SURGERY_BEGIN_TIME", "SURGERY_END_TIME", "SURGERY_NAME", "ANES_METHOD", "SURGERY_POSITION"])
    sd.to_csv(out_raw / "SURGERY_INFO.csv", index=False)

    d_items = pd.DataFrame([{"ROW_ID": i, "ITEMID": r[1], "LABEL": f"label {r[1]}"} for i, r in enumerate(S.CHANNEL_ROWS) if r[0] == "CHARTEVENTS"])
    d_labs = pd.DataFrame([{"ROW_ID": i, "ITEMID": r[1], "LABEL": f"label {r[1]}"} for i, r in enumerate(S.CHANNEL_ROWS) if r[0] == "LABEVENTS"])
    d_items.to_csv(out_raw / "D_ITEMS.csv", index=False)
    d_labs.to_csv(out_raw / "D_LABITEMS.csv", index=False)
    # Minimal PIC-schema placeholders (same columns as the fixture).
    pd.DataFrame([{"SUBJECT_ID": subjects[0], "HADM_ID": "1", "ITEMID": "1", "CHARTTIME": "", "VALUE": "1"}]).to_csv(
        out_raw / "OUTPUTEVENTS.csv", index=False)
    pd.DataFrame([{"SUBJECT_ID": subjects[0], "HADM_ID": "1", "CHARTTIME": "", "AMOUNT": "1"}]).to_csv(
        out_raw / "INPUTEVENTS.csv", index=False)

    for name in ("feature_channel_map_v1.csv", "drug_group_map_v1.csv", "specimen_map_v1.csv",
                 "organism_map_v1.csv", "neuro_admitting_dx_icd10_v1.csv"):
        src = maps_dir / name
        (out_raw / name).write_bytes(src.read_bytes())

    factors = enc[["encounter_id", "SUBJECT_ID", "HADM_ID", "z_severity", "z_respiratory", "z_age"]].copy()
    raw_files = {p.name: {"rows": _count_rows(p), "bytes": p.stat().st_size, "sha256": sha256_file(p)}
                 for p in sorted(out_raw.glob("*.csv"))}
    raw_bytes = sum(v["bytes"] for v in raw_files.values())
    return {
        "enc": enc, "factors": factors, "n_subjects": n_subjects, "n_encounters": n_enc,
        "n_invalid_dob_subjects": int(len(inv_subjects)), "invalid_subjects": sorted(inv_subjects),
        "raw_files": raw_files, "raw_bytes": raw_bytes, "map_sha": map_sha, "map_rows": map_rows,
        "maps_dir": maps_dir,
    }


def _count_rows(p: Path) -> int:
    # header + data; return data rows
    n = 0
    with p.open("rb") as f:
        for i, _ in enumerate(f):
            n = i
    return n  # lines including header, minus 1


# =============================================================================
# Stage 2: normalized (from raw, v1 maps and D1-D20 rules)
# =============================================================================
def _prep_events(df: pd.DataFrame, time_col: str, lookup: dict) -> pd.DataFrame:
    out = df.copy()
    key = M.vnorm_id(out["SUBJECT_ID"]) + "|" + M.vnorm_id(out["HADM_ID"])
    out["enc_i"] = key.map(lookup)
    out = out[out["enc_i"].notna()].copy()
    out["enc_i"] = out["enc_i"].astype(np.int32)
    if "ITEMID" in out.columns:
        out["item"] = out["ITEMID"].astype(str)
    t = M.fast_time(out[time_col])
    return out, t


def normalize(out_raw: Path, out_norm: Path, enc_meta: pd.DataFrame) -> dict:
    out_norm.mkdir(parents=True, exist_ok=True)
    canon = pd.read_csv(out_raw / "canonical_icu_encounters.csv", dtype=str)
    pats = pd.read_csv(out_raw / "PATIENTS.csv", dtype=str)
    icu = pd.read_csv(out_raw / "ICUSTAYS.csv", dtype=str)
    adms = pd.read_csv(out_raw / "ADMISSIONS.csv", dtype=str)
    splits = pd.read_csv(out_raw / "subject_splits.csv", dtype=str)
    chart = pd.read_csv(out_raw / "CHARTEVENTS.csv", dtype=str)
    lab = pd.read_csv(out_raw / "LABEVENTS.csv", dtype=str)
    rx = pd.read_csv(out_raw / "PRESCRIPTIONS.csv", dtype=str)

    cm = M.parse_channel_map(pd.read_csv(out_raw / "feature_channel_map_v1.csv", dtype=str).fillna(""))
    dm = M.parse_drug_map(pd.read_csv(out_raw / "drug_group_map_v1.csv", dtype=str).fillna(""))
    neuro = M.neuro_matcher(pd.read_csv(out_raw / "neuro_admitting_dx_icd10_v1.csv", dtype=str).fillna(""))

    enc = canon.copy()
    enc["SUBJECT_ID"] = M.vnorm_id(enc["SUBJECT_ID"])
    enc["HADM_ID"] = M.vnorm_id(enc["HADM_ID"])
    enc["t0"] = M.fast_time(enc["INTIME"])
    enc["out"] = M.fast_time(enc["OUTTIME"])
    enc["enc_i"] = np.arange(len(enc), dtype=np.int32)
    lookup = dict(zip(enc["SUBJECT_ID"] + "|" + enc["HADM_ID"], enc["enc_i"]))

    pats = pats.copy()
    pats["SUBJECT_ID"] = M.vnorm_id(pats["SUBJECT_ID"])
    pats["GENDER"] = M.vclean_text(pats["GENDER"].fillna(""))
    # one PATIENTS row per subject
    pats_u = pats.drop_duplicates("SUBJECT_ID", keep="first")
    enc = enc.merge(pats_u[["SUBJECT_ID", "GENDER", "DOB"]], on="SUBJECT_ID", how="left")
    has_row = enc["DOB"].notna() | (enc["SUBJECT_ID"].isin(set(pats_u["SUBJECT_ID"])))
    # merge leaves DOB NaN if no patient row
    has_row = enc["SUBJECT_ID"].isin(set(pats_u["SUBJECT_ID"])).to_numpy()
    dob_raw = enc["DOB"].fillna("")
    dob_parsed = M.fast_time(dob_raw)
    status, age_days, nbin = M.classify_ages(enc["t0"], dob_parsed, dob_raw, has_row)
    enc["dob_status"] = status
    enc["age_days"] = age_days
    enc["neg_age_bin"] = nbin
    enc["age_band"] = M.age_band(age_days)
    enc["sex"] = enc["GENDER"].replace({"": "U"})

    icu = icu.copy()
    icu["SUBJECT_ID"] = M.vnorm_id(icu["SUBJECT_ID"])
    icu["HADM_ID"] = M.vnorm_id(icu["HADM_ID"])
    icu["INTIME"] = M.fast_time(icu["INTIME"])
    icu = icu.sort_values("INTIME").drop_duplicates(["SUBJECT_ID", "HADM_ID"], keep="first")
    enc = enc.merge(icu[["SUBJECT_ID", "HADM_ID", "FIRST_CAREUNIT", "ICUSTAY_ID"]], on=["SUBJECT_ID", "HADM_ID"], how="left")

    splits["SUBJECT_ID"] = M.vnorm_id(splits["SUBJECT_ID"])
    enc = enc.merge(splits, on="SUBJECT_ID", how="left")

    adms = adms.copy()
    adms["SUBJECT_ID"] = M.vnorm_id(adms["SUBJECT_ID"])
    adms["HADM_ID"] = M.vnorm_id(adms["HADM_ID"])
    adms["dx_norm"] = [M.norm_icd(x)[0] for x in adms.get("ICD10_CODE_CN", pd.Series([""] * len(adms)))]
    adms["neuro_dx"] = [bool(neuro(c)) for c in adms["dx_norm"]]
    enc = enc.merge(adms[["SUBJECT_ID", "HADM_ID", "dx_norm", "neuro_dx"]], on=["SUBJECT_ID", "HADM_ID"], how="left")
    enc["neuro_dx"] = enc["neuro_dx"].fillna(False).astype(bool)

    def to_value_frame(df, time_col):
        ev, t = _prep_events(df, time_col, lookup)
        ev["th"] = (t - ev["enc_i"].map(dict(zip(enc["enc_i"], enc["t0"])))).dt.total_seconds() / 3600.0
        ev = ev[t.notna() & ev["th"].notna()].copy()
        ev["vnum"] = pd.to_numeric(ev.get("VALUENUM", ev.get("VALUE")), errors="coerce")
        ev["vtext"] = M.vclean_text(ev["VALUE"].fillna("")) if "VALUE" in ev.columns else ""
        ev["unit"] = ev["VALUEUOM"].fillna("").astype(str).str.strip().str.lower() if "VALUEUOM" in ev.columns else ""
        return ev[["enc_i", "item", "th", "vnum", "vtext", "unit"]]

    chart_v = to_value_frame(chart, "CHARTTIME")
    lab_v = to_value_frame(lab, "CHARTTIME")
    vb = M.build_values(chart_v, lab_v, cm)

    # Long-format observations: mapped rows (all statuses) + derived valid MAP / P/F / S/F.
    parts = []
    if len(vb.rowstat):
        rs = vb.rowstat.copy()
        rs["source_itemid"] = rs["item"]
        rs["source_table"] = rs["table"]
        rs["value"] = rs["v"]
        rs["t_hours"] = rs["t"]
        rs["t_chart_hours"] = rs["tc"]
        parts.append(rs[["enc_i", "t_hours", "t_chart_hours", "channel", "value", "unit", "status",
                         "source_itemid", "source_table"]])
    # Overwrite / add derived channels already in events (valid only)
    for ch, item_label, unit in (("map", "derived:map", "mmHg"), ("pf_ratio", "derived:pf_ratio", "mmHg"),
                                 ("sf_ratio", "derived:sf_ratio", "ratio")):
        ev = vb.events.get(ch, M._empty_events())
        if ev is None or len(ev) == 0:
            continue
        tmp = ev.copy()
        tmp["channel"] = ch
        tmp["value"] = tmp["v"]
        tmp["t_hours"] = tmp["t"]
        tmp["t_chart_hours"] = tmp["tc"]
        tmp["unit"] = unit
        tmp["status"] = "measurement"
        tmp["source_itemid"] = item_label
        tmp["source_table"] = "DERIVED"
        parts.append(tmp[["enc_i", "t_hours", "t_chart_hours", "channel", "value", "unit", "status",
                          "source_itemid", "source_table"]])
    obs = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        columns=["enc_i", "t_hours", "t_chart_hours", "channel", "value", "unit", "status",
                 "source_itemid", "source_table"])
    obs = obs.merge(enc[["enc_i", "encounter_id"]], on="enc_i", how="left")
    # Channel units for clean measurement rows: recorded map unit (after per-row conversions only).
    unit_by_ch = {ch: info.get("units") or "" for ch, info in
                  ((cm["channels"][c]["items"][0] if False else c, {"units": ""}) for c in [])}
    # prefer first assigned item's units
    ch_units = {}
    for (table, itemid), info in cm["items"].items():
        ch_units.setdefault(info["channel"], info.get("units") or "")
    ch_units["map"] = "mmHg"
    ch_units["pf_ratio"] = "mmHg"
    ch_units["sf_ratio"] = "ratio"
    meas_mask = (obs["status"] == "measurement") & (obs["source_table"] != "DERIVED")
    obs.loc[meas_mask, "unit"] = obs.loc[meas_mask, "channel"].map(ch_units).fillna(obs.loc[meas_mask, "unit"])
    obs.loc[obs["source_table"] == "DERIVED", "unit"] = obs.loc[obs["source_table"] == "DERIVED", "channel"].map(ch_units)

    obs_out = obs.sort_values(["encounter_id", "channel", "t_hours"], kind="mergesort")[
        ["encounter_id", "t_hours", "t_chart_hours", "channel", "value", "unit", "status", "source_itemid", "source_table"]
    ]
    obs_out.to_csv(out_norm / "observations_long.csv", index=False)

    # Medications with resolved groups (osmotic_conditional -> osmotherapy or otherwise, D3).
    rx2, tstart = _prep_events(rx, "STARTDATE", lookup)
    tend = M.fast_time(rx2["ENDDATE"]) if "ENDDATE" in rx2.columns else pd.Series(pd.NaT, index=rx2.index)
    t0_map = dict(zip(enc["enc_i"], enc["t0"]))
    rx2["start_hours"] = (tstart - rx2["enc_i"].map(t0_map)).dt.total_seconds() / 3600.0
    rx2["end_hours"] = (tend - rx2["enc_i"].map(t0_map)).dt.total_seconds() / 3600.0
    keys = (rx2["DRUG_NAME"].fillna("").astype(str) + " | " + rx2["DRUG_NAME_EN"].fillna("").astype(str))
    neuro_of = dict(zip(enc["enc_i"], enc["neuro_dx"]))
    groups, reasons = [], []
    for k, ei, sh, eh in zip(keys, rx2["enc_i"], rx2["start_hours"], rx2["end_hours"]):
        info = dm.get(k)
        if info is None:
            groups.append("unmapped")
            reasons.append("unmapped")
            continue
        g = info["group"]
        if g == "osmotic_conditional":
            if bool(neuro_of.get(int(ei), False)):
                groups.append("osmotherapy")
                reasons.append("neuro_admitting_dx")
            else:
                groups.append(info["otherwise"] or "other")
                reasons.append("otherwise")
        else:
            groups.append(g)
            reasons.append("map")
    rx2["group"] = groups
    rx2["group_reason"] = reasons
    rx2["cutoff_hours"] = [M.group_cutoff(g) if g in M.MED_GROUPS else M.DEFAULT_CUTOFF_H for g in rx2["group"]]
    # missing or inverted end -> start + group cutoff
    end = rx2["end_hours"].to_numpy(dtype=float)
    start = rx2["start_hours"].to_numpy(dtype=float)
    cut = rx2["cutoff_hours"].to_numpy(dtype=float)
    bad_end = ~np.isfinite(end) | (end <= start)
    end = np.where(bad_end, start + cut, end)
    rx2["end_hours"] = end
    rx2["start_status"] = np.where(tstart.isna(), "missing_start", "ok")
    rx2 = rx2.merge(enc[["enc_i", "encounter_id"]], on="enc_i", how="left")
    meds = rx2[rx2["start_status"] == "ok"][
        ["encounter_id", "DRUG_NAME", "DRUG_NAME_EN", "group", "group_reason", "start_hours", "end_hours", "cutoff_hours"]
    ].rename(columns={"DRUG_NAME": "drug_name", "DRUG_NAME_EN": "drug_name_en"})
    meds.to_csv(out_norm / "medications_long.csv", index=False)

    enc_out = enc[["encounter_id", "SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "INTIME", "OUTTIME", "split",
                   "sex", "FIRST_CAREUNIT", "age_days", "age_band", "dob_status", "neg_age_bin",
                   "neuro_dx", "Died"]].rename(columns={"FIRST_CAREUNIT": "first_careunit", "Died": "died"})
    enc_out["los_hours"] = (enc["out"] - enc["t0"]).dt.total_seconds() / 3600.0
    enc_out.to_csv(out_norm / "encounters.csv", index=False)

    n_invalid = int((enc["dob_status"] != "valid").sum())
    return {
        "enc": enc, "obs": obs, "meds": meds, "cm": cm, "events": vb.events,
        "n_obs": int(len(obs_out)), "n_meds": int(len(meds)), "n_invalid_dob": n_invalid,
        "dob_status_counts": enc["dob_status"].value_counts().to_dict(),
        "neg_bin_counts": enc.loc[enc["dob_status"] == "negative", "neg_age_bin"].value_counts().to_dict(),
    }


# =============================================================================
# Stage 3: features from normalized (0-24 h clinical summaries)
# =============================================================================
def _window_meas(obs: pd.DataFrame, channel: str, lo=WINDOW[0], hi=WINDOW[1]) -> pd.DataFrame:
    d = obs[(obs["channel"] == channel) & (obs["status"] == "measurement") &
            obs["value"].notna() & (obs["t_hours"] >= lo) & (obs["t_hours"] < hi)]
    return d


def _summaries(obs: pd.DataFrame, n_enc: int, channels: list[str]) -> pd.DataFrame:
    rows = {"enc_i": np.arange(n_enc, dtype=np.int32)}
    meas = obs[(obs["status"] == "measurement") & obs["value"].notna() &
               (obs["t_hours"] >= WINDOW[0]) & (obs["t_hours"] < WINDOW[1])]
    for ch in channels:
        g = meas[meas["channel"] == ch]
        if len(g) == 0:
            for stat in ("last", "min", "max", "mean", "count"):
                rows[f"{ch}_{stat}_0_24h"] = np.zeros(n_enc) if stat == "count" else np.full(n_enc, np.nan)
            rows[f"{ch}_missing_0_24h"] = np.ones(n_enc, dtype=np.int8)
            continue
        # last = value at latest t_hours
        g = g.sort_values(["enc_i", "t_hours"], kind="mergesort")
        last = g.groupby("enc_i", sort=False)["value"].last()
        agg = g.groupby("enc_i")["value"].agg(["min", "max", "mean", "count"])
        idx = agg.index.to_numpy(dtype=int)
        for stat in ("min", "max", "mean"):
            arr = np.full(n_enc, np.nan)
            arr[idx] = agg[stat].to_numpy(dtype=float)
            rows[f"{ch}_{stat}_0_24h"] = arr
        arr = np.zeros(n_enc, dtype=float)
        arr[idx] = agg["count"].to_numpy(dtype=float)
        rows[f"{ch}_count_0_24h"] = arr
        arr = np.full(n_enc, np.nan)
        arr[last.index.to_numpy(dtype=int)] = last.to_numpy(dtype=float)
        rows[f"{ch}_last_0_24h"] = arr
        rows[f"{ch}_missing_0_24h"] = (arr != arr).astype(np.int8)  # noqa: will use count==0
        rows[f"{ch}_missing_0_24h"] = (rows[f"{ch}_count_0_24h"] == 0).astype(np.int8)
    return pd.DataFrame(rows)


def _psofa_resp(obs: pd.DataFrame, n_enc: int) -> tuple[np.ndarray, np.ndarray]:
    """Worst (max) pSOFA respiratory score in 0-24 h: P/F if any valid, else S/F (D15/D19 cut points)."""
    score = np.full(n_enc, np.nan)
    src = np.full(n_enc, "", dtype=object)
    pf = _window_meas(obs, "pf_ratio")
    if len(pf):
        sc = M.psofa_score(pf["value"].to_numpy(dtype=float), M.PSOFA_PF)
        tmp = pd.DataFrame({"enc_i": pf["enc_i"].to_numpy(dtype=int), "sc": sc})
        worst = tmp.groupby("enc_i")["sc"].max()
        score[worst.index.to_numpy(dtype=int)] = worst.to_numpy(dtype=float)
        src[worst.index.to_numpy(dtype=int)] = "pf_ratio"
    need = np.isnan(score)
    sf = _window_meas(obs, "sf_ratio")
    if len(sf) and need.any():
        sc = M.psofa_score(sf["value"].to_numpy(dtype=float), M.PSOFA_SF)
        tmp = pd.DataFrame({"enc_i": sf["enc_i"].to_numpy(dtype=int), "sc": sc})
        worst = tmp.groupby("enc_i")["sc"].max()
        ix = worst.index.to_numpy(dtype=int)
        use = need[ix]
        score[ix[use]] = worst.to_numpy(dtype=float)[use]
        src[ix[use]] = "sf_ratio"
    return score, src


def _concurrent_and_exposure(meds: pd.DataFrame, enc: pd.DataFrame) -> pd.DataFrame:
    n_enc = len(enc)
    ei_of = dict(zip(enc["encounter_id"], enc["enc_i"]))
    flags = {g: np.zeros(n_enc, dtype=np.int8) for g in M.MED_GROUPS}
    vaso_conc = np.zeros(n_enc, dtype=float)
    if len(meds) == 0:
        out = pd.DataFrame({"enc_i": np.arange(n_enc)})
        for g in M.MED_GROUPS:
            out[f"{g}_exposed_0_24h"] = flags[g]
        out["vasoactive_max_concurrent_0_24h"] = vaso_conc
        return out
    m = meds.copy()
    m["enc_i"] = m["encounter_id"].map(ei_of)
    m = m[m["enc_i"].notna()].copy()
    m["enc_i"] = m["enc_i"].astype(int)
    lo, hi = WINDOW
    # overlap with [0, 24)
    overlap = (m["start_hours"] < hi) & (m["end_hours"] > lo)
    m = m[overlap]
    for g, sub in m.groupby("group"):
        if g in flags:
            flags[g][sub["enc_i"].unique()] = 1
    # max concurrent vasoactives in the window: sweep unique start/end times
    vaso = m[m["group"] == "vasoactive"]
    for ei, sub in vaso.groupby("enc_i"):
        starts = np.clip(sub["start_hours"].to_numpy(dtype=float), lo, hi)
        ends = np.clip(sub["end_hours"].to_numpy(dtype=float), lo, hi)
        evs = [(t, 1) for t in starts] + [(t, -1) for t in ends]
        evs.sort(key=lambda x: (x[0], x[1]))  # ends (-1) before starts at the same instant? use ends first so zero-length don't count
        # actually at a shared boundary we want starts after ends if one stops when another starts... count at starts.
        evs.sort(key=lambda x: (x[0], x[1]))
        c = 0
        mx = 0
        for t, d in evs:
            c += d
            if c > mx:
                mx = c
        vaso_conc[int(ei)] = mx
    out = pd.DataFrame({"enc_i": np.arange(n_enc)})
    for g in M.MED_GROUPS:
        out[f"{g}_exposed_0_24h"] = flags[g]
    out["vasoactive_max_concurrent_0_24h"] = vaso_conc
    return out


def impute_weight(enc: pd.DataFrame, weight: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Charted weight if present; else train median by age-band x sex, with D20 fallback for age_invalid."""
    n = len(enc)
    src = np.where(np.isfinite(weight), "charted", "").astype(object)
    out = weight.copy()
    train = (enc["split"].to_numpy() == "train") & np.isfinite(weight)
    bands = enc["age_band"].to_numpy()
    sex = enc["sex"].to_numpy()
    # medians from train charted
    def med(mask):
        v = weight[mask]
        v = v[np.isfinite(v)]
        return float(np.median(v)) if len(v) else np.nan

    overall = med(train)
    sex_med = {s: med(train & (sex == s)) for s in pd.unique(sex)}
    band_sex = {}
    band_only = {}
    for b in pd.unique(bands):
        band_only[b] = med(train & (bands == b) & (b != "age_invalid"))
        for s in pd.unique(sex):
            band_sex[(b, s)] = med(train & (bands == b) & (sex == s) & (b != "age_invalid"))
    for i in range(n):
        if np.isfinite(out[i]):
            continue
        b, s = bands[i], sex[i]
        val, how = np.nan, "overall_median"
        if b != "age_invalid":
            val = band_sex.get((b, s), np.nan)
            how = "age_band_sex_median"
            if not np.isfinite(val):
                val = band_only.get(b, np.nan)
                how = "age_band_median"
        if not np.isfinite(val):
            val = sex_med.get(s, np.nan)
            how = "sex_median"
        if not np.isfinite(val):
            val = overall
            how = "overall_median"
        out[i] = val
        src[i] = how if np.isfinite(val) else "missing"
    return out, src


def engineer_features(out_norm: Path, out_feat: Path, norm: dict) -> pd.DataFrame:
    out_feat.mkdir(parents=True, exist_ok=True)
    enc = norm["enc"].copy()
    obs = norm["obs"]
    meds = norm["meds"]
    cm = norm["cm"]
    n_enc = len(enc)
    channels = [ch for ch, c in cm["channels"].items() if not c.get("derived") or ch in ("map", "pf_ratio", "sf_ratio")]
    # skip the lab P/F crosscheck from the main SVD list but still summarize it
    feat = _summaries(obs, n_enc, channels)
    psofa, psofa_src = _psofa_resp(obs, n_enc)
    feat["psofa_resp_0_24h"] = psofa
    feat["psofa_resp_source_0_24h"] = psofa_src
    feat["psofa_resp_missing_0_24h"] = np.isnan(psofa).astype(np.int8)
    drugs = _concurrent_and_exposure(meds, enc)
    feat = feat.merge(drugs, on="enc_i", how="left")

    # hypotension indicator (simple, not a full age-adjusted pSOFA CV score)
    map_w = _window_meas(obs, "map")
    map_lt70 = np.zeros(n_enc, dtype=np.int8)
    if len(map_w):
        hit = map_w.loc[map_w["value"] < 70, "enc_i"].astype(int).unique()
        map_lt70[hit] = 1
    feat["map_lt70_0_24h"] = map_lt70
    lac_w = _window_meas(obs, "lactate")
    for thr, name in ((2.0, "lactate_gt2_0_24h"), (4.0, "lactate_gt4_0_24h")):
        arr = np.zeros(n_enc, dtype=np.int8)
        if len(lac_w):
            arr[lac_w.loc[lac_w["value"] > thr, "enc_i"].astype(int).unique()] = 1
        feat[name] = arr

    wcol = feat["weight_last_0_24h"] if "weight_last_0_24h" in feat.columns else np.full(n_enc, np.nan)
    # prefer mean/last charted weight
    wt_chart = feat["weight_mean_0_24h"] if "weight_mean_0_24h" in feat.columns else wcol
    wt_imp, wt_src = impute_weight(enc, wt_chart.to_numpy(dtype=float))
    feat["weight_kg"] = wt_imp
    feat["weight_source"] = wt_src
    intake = feat["fluid_intake_mean_0_24h"] if "fluid_intake_mean_0_24h" in feat.columns else np.full(n_enc, np.nan)
    # daily total: sum of charted intake in the window (not the mean)
    intake_sum = np.zeros(n_enc, dtype=float)
    intake_has = np.zeros(n_enc, dtype=bool)
    inn = _window_meas(obs, "fluid_intake")
    if len(inn):
        sm = inn.groupby("enc_i")["value"].sum()
        ix = sm.index.to_numpy(dtype=int)
        intake_sum[ix] = sm.to_numpy(dtype=float)
        intake_has[ix] = True
    feat["fluid_intake_ml_0_24h"] = np.where(intake_has, intake_sum, np.nan)
    feat["fluid_intake_ml_per_kg_per_day"] = feat["fluid_intake_ml_0_24h"] / feat["weight_kg"]

    static = enc[["enc_i", "encounter_id", "SUBJECT_ID", "HADM_ID", "split", "sex", "FIRST_CAREUNIT",
                  "age_days", "age_band", "dob_status", "neg_age_bin", "neuro_dx", "Died"]].copy()
    static = static.rename(columns={"FIRST_CAREUNIT": "first_careunit", "Died": "died"})
    static["los_hours"] = (enc["out"] - enc["t0"]).dt.total_seconds() / 3600.0
    static["age_invalid"] = (enc["age_band"].to_numpy() == "age_invalid").astype(np.int8)
    out = static.merge(feat, on="enc_i", how="left").drop(columns=["enc_i"])
    out.to_csv(out_feat / "features_encounter_0_24h.csv", index=False)
    # Compact SVD analysis matrix (one summary per construct + noise means).
    svd_df = out[["encounter_id"] + [c for c in SVD_COLS if c in out.columns]].copy()
    svd_df.to_csv(out_feat / "features_svd_matrix.csv", index=False)
    return out


# =============================================================================
# SVD verification against planted factors
# =============================================================================
def standardize(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0, ddof=1)
    sd = np.where(sd < 1e-12, 1.0, sd)
    Xc = np.where(np.isfinite(X), X, mu)
    Z = (Xc - mu) / sd
    return Z, mu, sd


def induced_loadings(Xstd: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """OLS coefficients of standardized features on the three factor scores (p x 3)."""
    # Z is n x 3; add no intercept because X and Z are (approx) centered
    B, *_ = np.linalg.lstsq(Z, Xstd, rcond=None)
    return B.T


def orthonormalize(B: np.ndarray) -> np.ndarray:
    Q, _ = np.linalg.qr(B)
    # flip columns so they align in sign with original
    for j in range(Q.shape[1]):
        if np.dot(Q[:, j], B[:, j]) < 0:
            Q[:, j] *= -1
    return Q


def procrustes_align(V: np.ndarray, B: np.ndarray):
    """Orthogonal Procrustes: Q = argmin ||V Q - B||. Returns Q, aligned V, cosines, principal angles (deg)."""
    U, s, Vt = np.linalg.svd(V.T @ B, full_matrices=False)
    Q = U @ Vt
    aligned = V @ Q
    cos = []
    for j in range(V.shape[1]):
        a, b = aligned[:, j], B[:, j]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        cos.append(float(np.dot(a, b) / (na * nb)) if na * nb else 0.0)
    angles = np.degrees(np.arccos(np.clip(s, 0.0, 1.0)))
    return Q, aligned, cos, angles, s


def verify_svd(features: pd.DataFrame, factors: pd.DataFrame, out_key: Path) -> dict:
    out_key.mkdir(parents=True, exist_ok=True)
    df = features.merge(factors, on="encounter_id", how="inner", suffixes=("", "_y"))
    cols = [c for c in SVD_COLS if c in df.columns]
    X = df[cols].to_numpy(dtype=float)
    Z = df[list(FACTOR_NAMES)].to_numpy(dtype=float)
    Xstd, mu, sd = standardize(X)
    # center factors (should already be ~mean 0) and keep raw z in the answer key
    Zc = Z - Z.mean(axis=0, keepdims=True)
    U, S, Vh = np.linalg.svd(Xstd, full_matrices=False)
    ev = (S ** 2)
    ev_frac = ev / ev.sum()
    V3 = Vh[:3, :].T  # p x 3 right singular vectors
    B_ols = induced_loadings(Xstd, Zc)
    B_on = orthonormalize(B_ols)
    Q, aligned, cos, angles, s_proc = procrustes_align(V3, B_on)
    # also cosine of unaligned (best sign flip only) for each pair in order
    sign_cos = []
    for j in range(3):
        c = float(np.dot(V3[:, j], B_on[:, j]))
        sign_cos.append(c if abs(c) >= 0 else c)

    loadings = pd.DataFrame({
        "feature": cols,
        "role": ["signal" if c in SIGNAL_COLS else "noise" for c in cols],
        "mean_used_for_standardization": mu,
        "sd_used_for_standardization": sd,
        "loading_severity_ols": B_ols[:, 0],
        "loading_respiratory_ols": B_ols[:, 1],
        "loading_age_ols": B_ols[:, 2],
        "loading_severity_orthonormal": B_on[:, 0],
        "loading_respiratory_orthonormal": B_on[:, 1],
        "loading_age_orthonormal": B_on[:, 2],
        "right_sv1": V3[:, 0],
        "right_sv2": V3[:, 1],
        "right_sv3": V3[:, 2],
        "right_sv1_procrustes": aligned[:, 0],
        "right_sv2_procrustes": aligned[:, 1],
        "right_sv3_procrustes": aligned[:, 2],
    })
    loadings.to_csv(out_key / "planted_loadings.csv", index=False)
    fac_out = df[["encounter_id", "SUBJECT_ID", "HADM_ID"] + list(FACTOR_NAMES)].copy()
    fac_out.to_csv(out_key / "latent_factors.csv", index=False)

    top_k = 10
    report = {
        "n_encounters": int(len(df)),
        "n_features": int(len(cols)),
        "features": cols,
        "missing_imputation": "column mean of the standardized SVD matrix",
        "standardization": "column mean 0, sample SD 1 (ddof=1), after mean imputation",
        "singular_values_top10": [float(x) for x in S[:top_k]],
        "variance_explained_top10": [float(x) for x in ev_frac[:top_k]],
        "variance_explained_top3_sum": float(ev_frac[:3].sum()),
        "eigenvalues_of_gram_top10": [float(x) for x in ev[:top_k]],
        "gap_sv3_over_sv4": float(S[2] / S[3]) if len(S) > 3 and S[3] else None,
        "procrustes_cosine_similarity": {
            FACTOR_NAMES[j]: cos[j] for j in range(3)
        },
        "principal_angles_degrees": [float(x) for x in angles],
        "procrustes_singular_values": [float(x) for x in s_proc],
        "unaligned_cosine_same_order": {
            FACTOR_NAMES[j]: sign_cos[j] for j in range(3)
        },
        "factor_score_moments": {
            name: {"mean": float(Z[:, j].mean()), "sd": float(Z[:, j].std(ddof=1))}
            for j, name in enumerate(FACTOR_NAMES)
        },
        "factor_score_correlations": pd.DataFrame(Z, columns=list(FACTOR_NAMES)).corr().round(6).to_dict(),
    }
    (out_key / "svd_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


# =============================================================================
# Data dictionary, README, zips
# =============================================================================
FEATURE_DOCS = [
    ("encounter_id", "id", "Synthetic encounter identifier (E{SUBJECT_ID}_{k})."),
    ("SUBJECT_ID", "id", "Synthetic patient identifier."),
    ("HADM_ID", "id", "Synthetic hospital admission identifier (normalized)."),
    ("split", "static", "Subject-level split: train / validation / test (~60/20/20)."),
    ("sex", "static", "Sex from PATIENTS.GENDER (M/F; U if blank). Not planted; treated as noise."),
    ("first_careunit", "static", "FIRST_CAREUNIT of the first ICU stay. Mildly age-related (NICU more common at young ages) but not used in the SVD analysis matrix."),
    ("age_days", "signal-age", "INTIME - DOB in days. NaN and age_band=age_invalid when DOB is missing, unparseable, or yields a negative age (D20). Driven by z_age."),
    ("age_band", "signal-age", "neonate / infant / toddler / preschool / school / adolescent / ge18 / age_invalid (D20)."),
    ("dob_status", "static", "valid / missing / unparseable / negative."),
    ("neg_age_bin", "static", "For negative ages: lt1d, 1to30d, gt30d; otherwise blank."),
    ("age_invalid", "static", "1 if age_band is age_invalid."),
    ("neuro_dx", "static", "Admitting ICD-10 matches the v1 neuro set (D3), used only to resolve osmotic_conditional drugs."),
    ("died", "static", "In-hospital death flag (weakly related to severity in raw generation; not in the SVD analysis matrix)."),
    ("los_hours", "static", "OUTTIME - INTIME in hours."),
    ("weight_kg", "signal-age", "Charted weight (item 1014) in 0-24 h if present, else imputed (age-band x sex median on train; age_invalid falls back to sex or overall median, D20). Driven by z_age via an age-weight curve plus a residual log-weight loading."),
    ("weight_source", "static", "charted / age_band_sex_median / age_band_median / sex_median / overall_median / missing."),
    ("fluid_intake_ml_0_24h", "noise", "Sum of charted fluid_intake (item 1009, mL) in [0, 24) h. Not factor-driven."),
    ("fluid_intake_ml_per_kg_per_day", "noise", "fluid_intake_ml_0_24h / weight_kg."),
    ("psofa_resp_0_24h", "signal-respiratory", "Worst (max) pSOFA respiratory score in [0, 24) h. P/F preferred when present (cuts D19: >=400 0, 300-399 1, 200-299 2, 100-199 3, <100 4); otherwise S/F (cuts D15: >=292 0, 264-291 1, 221-263 2, 148-220 3, <148 4). No ventilation requirement (D19)."),
    ("psofa_resp_source_0_24h", "static", "pf_ratio or sf_ratio."),
    ("psofa_resp_missing_0_24h", "missingness", "1 if neither P/F nor S/F is available in the window."),
    ("map_lt70_0_24h", "signal-severity", "1 if any valid derived MAP in the window is < 70 mmHg. Simple hypotension flag, not a full age-adjusted pSOFA cardiovascular score."),
    ("lactate_gt2_0_24h", "signal-severity", "1 if any lactate in the window is > 2 mmol/L."),
    ("lactate_gt4_0_24h", "signal-severity", "1 if any lactate in the window is > 4 mmol/L."),
    ("vasoactive_exposed_0_24h", "signal-severity", "1 if any vasoactive course overlaps [0, 24) h (12 h default duration when ENDDATE is missing, D-group cutoff)."),
    ("vasoactive_max_concurrent_0_24h", "signal-severity", "Maximum number of vasoactive agents simultaneously active in [0, 24) h."),
]


def document_features(feat: pd.DataFrame, path: Path) -> pd.DataFrame:
    known = {k: (role, doc) for k, role, doc in FEATURE_DOCS}
    rows = []
    for col in feat.columns:
        if col in known:
            role, doc = known[col]
            unit = ""
        elif col.endswith("_missing_0_24h"):
            role, doc, unit = "missingness", f"1 if {col.replace('_missing_0_24h', '')} has no valid measurement in [0, 24) h.", ""
        elif col.endswith("_exposed_0_24h"):
            g = col.replace("_exposed_0_24h", "")
            role, doc, unit = "medication", f"1 if any {g} course overlaps [0, 24) h. Other than vasoactive, groups are not factor-driven (noise).", ""
        elif col.endswith("_last_0_24h"):
            ch = col[:-len("_last_0_24h")]
            role = "signal" if f"{ch}_mean_0_24h" in SIGNAL_COLS else "channel-summary"
            doc, unit = f"Last valid {ch} in [0, 24) h (by t_hours, labs lagged +1 h release).", ""
        elif col.endswith("_min_0_24h"):
            ch = col[:-len("_min_0_24h")]
            role = "signal" if f"{ch}_mean_0_24h" in SIGNAL_COLS else "channel-summary"
            doc, unit = f"Minimum valid {ch} in [0, 24) h.", ""
        elif col.endswith("_max_0_24h"):
            ch = col[:-len("_max_0_24h")]
            role = "signal" if f"{ch}_mean_0_24h" in SIGNAL_COLS else "channel-summary"
            doc, unit = f"Maximum valid {ch} in [0, 24) h.", ""
        elif col.endswith("_mean_0_24h"):
            ch = col[:-len("_mean_0_24h")]
            role = "signal" if col in SIGNAL_COLS else ("noise" if col in NOISE_COLS else "channel-summary")
            doc, unit = f"Mean valid {ch} in [0, 24) h.", ""
        elif col.endswith("_count_0_24h"):
            ch = col[:-len("_count_0_24h")]
            role, doc, unit = "channel-summary", f"Count of valid {ch} measurements in [0, 24) h.", "count"
        else:
            role, doc, unit = "other", "See README.", ""
        rows.append({"feature": col, "role": role, "unit_or_coding": unit, "description": doc,
                     "in_svd_matrix": "yes" if col in SVD_COLS else "no"})
    dd = pd.DataFrame(rows)
    dd.to_csv(path, index=False)
    return dd


def zip_dir(src: Path, dest_zip: Path) -> None:
    dest_zip.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(dest_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(src.rglob("*")):
            if p.is_file() and p.suffix != ".zip":
                z.write(p, arcname=p.relative_to(src))


def write_readme(root: Path, meta: dict, report: dict, gen_sha: str) -> None:
    inv = meta["dob_status_counts"]
    n_enc = meta["n_encounters"]
    n_sub = meta["n_subjects"]
    n_inv = meta["n_invalid_dob"]
    raw_mb = meta["raw_bytes"] / (1024 * 1024)
    sv = report["singular_values_top10"]
    ve = report["variance_explained_top10"]
    cos = report["procrustes_cosine_similarity"]
    ang = report["principal_angles_degrees"]
    gap = report["gap_sv3_over_sv4"]

    def fmt_sv(xs):
        return ", ".join(f"{x:.3f}" for x in xs)

    def fmt_pct(xs):
        return ", ".join(f"{100 * x:.2f}%" for x in xs)

    text = f"""# Synthetic dataset for optimizing SVD algorithms

One synthetic PIC-like dataset in three related stages, plus a held-out answer key of planted latent factors. Built for algorithm work on feature matrices (SVD and related decompositions). **Not real patient data.**

- Generator: `make_synthetic_svd_dataset.py` (SHA-256 `{gen_sha}`)
- Seed: `{meta['seed']}`
- Subjects: **{n_sub:,}**
- Cohort encounters: **{n_enc:,}**
- Invalid DOB / age-invalid encounters: **{n_inv:,}** ({100 * n_inv / n_enc:.2f}% of encounters; safety cap is 1%). Breakdown: missing={inv.get('missing', 0)}, unparseable={inv.get('unparseable', 0)}, negative={inv.get('negative', 0)}, valid={inv.get('valid', 0)}. Negative bins: {meta.get('neg_bin_counts', {})}.
- Raw CSV total: **{raw_mb:.1f} MB** (event density thinned to a {meta['vitals_step']} h vital step and a {EVENT_HORIZON_H:.0f} h event horizon so the raw stage stays under ~300 MB).

Stages are stacked: normalized is derived from raw; features are derived from normalized. Do not treat this dataset as an input to the Stage 2 audit pipeline (the audit script was not run on it and was not modified).

## Folder layout

- `1_raw/` — PIC-schema tables, v1 maps, canonical encounters, harmonization fields, subject split, surgery complexity map. Also `1_raw.zip`.
- `2_normalized/` — `encounters.csv`, `observations_long.csv` (channel, time, clean value/unit, status), `medications_long.csv` (resolved drug groups). Also `2_normalized.zip`.
- `3_features/` — `features_encounter_0_24h.csv` (full encounter × feature matrix), `features_svd_matrix.csv` (the recommended standardized-SVD analysis columns, unstandardized values), `data_dictionary.csv`. Also `3_features.zip`.
- `4_answer_key/` — **held-out** `latent_factors.csv` (true per-encounter factor scores), `planted_loadings.csv` (induced loadings of the standardized SVD matrix on those scores, plus recovered right singular vectors), `svd_verification.json`.

## Planted latent factors

Three independent standard normal scores are drawn **inside the raw generator**. `z_age` is patient-level (DOB is a patient attribute) and copied to every encounter of that patient. `z_severity` and `z_respiratory` are per-encounter.

Event values are drawn around means shifted by the factor scores, then clipped to the channel plausibility gates, plus observation noise. Other channels (temperature, bedside glucose, sodium, hemoglobin, magnesium, calcium, urine output, fluid intake, non-vasoactive drugs) get an encounter-level intercept independent of the three factors, so they are noise for SVD.

| Factor | What it does at the raw level (approximate loadings) |
|---|---|
| **z_severity** (illness) | Raises lactate ({BETA['lactate_sev']:+.2f} mmol/L per SD), creatinine ({BETA['creatinine_sev']:+.1f} μmol/L), heart rate ({BETA['hr_sev']:+.1f} bpm); lowers SBP ({BETA['sbp_sev']:+.1f} mmHg) and DBP ({BETA['dbp_sev']:+.1f} mmHg), hence MAP; raises the logit of vasoactive exposure (intercept {BETA['vaso_logit_intercept']:+.2f}, slope {BETA['vaso_logit_sev']:+.2f}) and the chance of multiple overlapping agents. Platelets and bilirubin are not in the v1 synthetic channel map, so those pSOFA components are absent. |
| **z_respiratory** | Lowers SpO2 ({BETA['spo2_resp']:+.1f} pp) and PaO2 ({BETA['po2_resp']:+.1f} mmHg); raises FiO2 ({BETA['fio2_resp']:+.1f} pp), which lowers derived P/F and S/F; raises respiratory rate ({BETA['rr_resp']:+.1f} /min). P/F uses same-draw pO2 + FiO2; S/F uses FiO2 without a valid P/F, nearest SpO2 within 5 min, only if SpO2 ≤ 97% (D15). |
| **z_age** | Sets DOB so age_years ≈ clip(9 + 6 z_age, 0.05, 21). Weight follows a pediatric age-weight curve times exp({BETA['weight_log_age']:.2f} z_age + N(0, {NOISE_SD['weight_log']:.2f})). Age-conditional norms: higher HR and RR, lower SBP/DBP/weight at younger ages (e.g. HR ≈ clip(155 − 4 age_years, 70, 180)). |

Observation noise (independent, per event): HR SD {NOISE_SD['hr']}, RR {NOISE_SD['rr']}, SBP {NOISE_SD['sbp']}, DBP {NOISE_SD['dbp']}, SpO2 {NOISE_SD['spo2']}, lactate {NOISE_SD['lactate']} mmol/L, creatinine {NOISE_SD['creatinine']} μmol/L, PaO2 {NOISE_SD['po2']} mmHg, FiO2 {NOISE_SD['fio2']} pp.

Invalid-DOB encounters keep their hidden `z_age` in the answer key, but `age_days` is missing and `age_band` is `age_invalid` (D20). They are < 1% of encounters (a mix of missing, unparseable, and negative ages in bins <1 d, 1-30 d, and >30 d).

## Normalization (stage 2)

Applied from the v1 maps and D1-D20 rules, using the audit module helpers (read-only):

- Item → channel mapping, plausibility gates, per-row unit conversions where the audit stream applies them.
- HR merge (D9): item 1003 kept when both 1002 and 1003 exist at the same CHARTTIME.
- MAP (D8): (SBP + 2×DBP)/3 at the same CHARTTIME; drop DBP ≥ SBP or pairs outside the SBP/DBP/MAP gates.
- P/F (D7, D19): pO2 and FiO2 5222 from the same draw; item 5239 preferred over 5244.
- S/F (D15): FiO2 draws without a valid P/F; nearest SpO2 1006 within 5 minutes; only if that SpO2 ≤ 97%.
- Labs use a +1 h release lag (`t_hours` is the analysis time).
- Drug groups from the v1 drug map; `osmotic_conditional` → osmotherapy if the admitting ICD matches the v1 neuro set, else the mapped otherwise-group (D3). Missing/inverted ENDDATE → start + group cutoff (12 h vasoactive, 24 h otherwise).
- Age = INTIME − DOB with the `age_invalid` band (D20).

`observations_long.csv` is long-format (one row per encounter / channel / time / source item) including out-of-gate and unit-status rows (`status`). Use `status == measurement` for clean values. Units on measurement rows are the channel-map recorded units (creatinine stays μmol/L in this stream; FiO2 stays %).

## Feature matrix (stage 3)

`features_encounter_0_24h.csv` has one row per encounter and, for each mapped channel (plus derived MAP, P/F, S/F): last, min, max, mean, count, and a missingness indicator over [0, 24) h, plus the clinical features listed in the data dictionary (pSOFA respiratory component, vasoactive flag and concurrent count, drug-group exposure flags, fluid intake mL/kg/day, age band, sex, care unit).

`features_svd_matrix.csv` is the recommended analysis matrix for SVD: the 0-24 h **means** of the planted and noise channels, plus age, weight, vasoactive exposure and concurrent count, and pSOFA respiratory. Last/min/max/count of the same channel are near-duplicates of the mean and would inflate within-channel singular values if they were left in the decomposition. All columns of this matrix are documented in `data_dictionary.csv` (`in_svd_matrix = yes`).

Standardization used for the verification SVD: column-mean imputation, then mean 0 and sample SD 1.

## SVD verification (numpy)

Run on the standardized stage-3 analysis matrix (`features_svd_matrix.csv`) against the held-out scores in `4_answer_key/latent_factors.csv`. Planted loading vectors are the **induced** loadings: OLS of the standardized features on the (centered) factor scores, then column-orthonormalized. Recovered right singular vectors are aligned to that basis by orthogonal Procrustes. Principal angles are from the singular values of V₃ᵀ B.

- n = {report['n_encounters']:,} encounters, p = {report['n_features']} features
- Top 10 singular values: {fmt_sv(sv)}
- Variance explained (top 10): {fmt_pct(ve)}
- Top-3 variance explained (sum): {100 * report['variance_explained_top3_sum']:.2f}%
- Gap σ₃/σ₄: {gap:.3f}
- Procrustes cosine similarity (severity, respiratory, age): {cos['z_severity']:.4f}, {cos['z_respiratory']:.4f}, {cos['z_age']:.4f}
- Principal angles (degrees): {ang[0]:.3f}, {ang[1]:.3f}, {ang[2]:.3f}

The answer-key files are the ground truth for recovery experiments. Do not use them as features.

## How the stages were generated

```
python3 make_synthetic_svd_dataset.py --n-subjects {n_sub} --seed {meta['seed']} --vitals-step {meta['vitals_step']}
```

Fixed seed, no edits to the Stage 2 audit script, tests, or fixture. The fixture module is imported only for v1 map writers and schema constants.

## Row counts

{meta['rowcount_md']}
"""
    (root / "README_synthetic_svd_levels.md").write_text(text)


def rowcount_markdown(meta: dict) -> str:
    lines = ["### 1_raw", "", "| file | rows | bytes |", "|---|---:|---:|"]
    for name, info in meta["raw_files"].items():
        lines.append(f"| `{name}` | {info['rows']:,} | {info['bytes']:,} |")
    lines += ["", "### 2_normalized", "",
              f"- `encounters.csv`: {meta['n_encounters']:,} rows (1 per encounter)",
              f"- `observations_long.csv`: {meta['n_obs']:,} rows",
              f"- `medications_long.csv`: {meta['n_meds']:,} rows",
              "", "### 3_features", "",
              f"- `features_encounter_0_24h.csv`: {meta['n_encounters']:,} rows (1 per encounter)",
              f"- `features_svd_matrix.csv`: {meta['n_encounters']:,} rows, {len(SVD_COLS)} analysis columns + encounter_id"]
    return "\n".join(lines)


def write_manifest(root: Path, meta: dict, report: dict, gen_sha: str) -> None:
    doc = {
        "generator": "make_synthetic_svd_dataset.py",
        "generator_sha256": gen_sha,
        "seed": meta["seed"],
        "n_subjects": meta["n_subjects"],
        "n_encounters": meta["n_encounters"],
        "n_invalid_dob": meta["n_invalid_dob"],
        "dob_status_counts": meta["dob_status_counts"],
        "neg_bin_counts": meta["neg_bin_counts"],
        "vitals_step_hours": meta["vitals_step"],
        "event_horizon_hours": EVENT_HORIZON_H,
        "raw_bytes": meta["raw_bytes"],
        "raw_files": meta["raw_files"],
        "n_observations_long": meta["n_obs"],
        "n_medications_long": meta["n_meds"],
        "svd_verification": {k: report[k] for k in (
            "n_encounters", "n_features", "singular_values_top10", "variance_explained_top10",
            "variance_explained_top3_sum", "gap_sv3_over_sv4", "procrustes_cosine_similarity",
            "principal_angles_degrees")},
        "beta_raw_scale": BETA,
        "noise_sd": NOISE_SD,
        "svd_columns": SVD_COLS,
    }
    (root / "manifest.json").write_text(json.dumps(doc, indent=2) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Build the three-stage synthetic SVD dataset.")
    ap.add_argument("--out-dir", type=Path, default=HERE)
    ap.add_argument("--n-subjects", type=int, default=N_SUBJECTS_DEFAULT)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--vitals-step", type=float, default=VITALS_STEP_DEFAULT)
    ap.add_argument("--skip-zip", action="store_true")
    args = ap.parse_args()
    root = args.out_dir.resolve()
    raw_dir = root / "1_raw"
    norm_dir = root / "2_normalized"
    feat_dir = root / "3_features"
    key_dir = root / "4_answer_key"
    for p in (raw_dir, norm_dir, feat_dir, key_dir):
        p.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] raw: n_subjects={args.n_subjects} seed={args.seed} vitals_step={args.vitals_step}", flush=True)
    raw_meta = generate_raw(raw_dir, args.n_subjects, args.seed, args.vitals_step)
    print(f"      encounters={raw_meta['n_encounters']} raw_MB={raw_meta['raw_bytes'] / (1024 * 1024):.1f}", flush=True)
    if raw_meta["raw_bytes"] > 300 * 1024 * 1024:
        print("WARNING: raw CSVs exceed ~300 MB; re-run with a larger --vitals-step to thin events.", flush=True)

    print("[2/4] normalized", flush=True)
    norm = normalize(raw_dir, norm_dir, raw_meta["enc"])
    print(f"      obs={norm['n_obs']} meds={norm['n_meds']} invalid_dob={norm['n_invalid_dob']} {norm['dob_status_counts']}", flush=True)
    if norm["n_invalid_dob"] > int(0.01 * raw_meta["n_encounters"]):
        raise SystemExit(f"invalid DOB fraction exceeds 1% ({norm['n_invalid_dob']} / {raw_meta['n_encounters']})")

    print("[3/4] features", flush=True)
    feat = engineer_features(norm_dir, feat_dir, norm)
    document_features(feat, feat_dir / "data_dictionary.csv")
    document_features(feat, root / "data_dictionary.csv")

    print("[4/4] SVD verification + answer key", flush=True)
    report = verify_svd(feat, raw_meta["factors"], key_dir)
    print("      singular values top 10:", [round(x, 3) for x in report["singular_values_top10"]], flush=True)
    print("      var explained top 10 %:", [round(100 * x, 2) for x in report["variance_explained_top10"]], flush=True)
    print("      Procrustes cosines:", report["procrustes_cosine_similarity"], flush=True)
    print("      principal angles (deg):", [round(x, 3) for x in report["principal_angles_degrees"]], flush=True)
    print("      gap sv3/sv4:", report["gap_sv3_over_sv4"], flush=True)

    gen_path = Path(__file__).resolve()
    gen_sha = sha256_file(gen_path)
    meta = {
        "seed": args.seed,
        "vitals_step": args.vitals_step,
        "n_subjects": raw_meta["n_subjects"],
        "n_encounters": raw_meta["n_encounters"],
        "n_invalid_dob": norm["n_invalid_dob"],
        "dob_status_counts": norm["dob_status_counts"],
        "neg_bin_counts": norm["neg_bin_counts"],
        "raw_files": raw_meta["raw_files"],
        "raw_bytes": raw_meta["raw_bytes"],
        "n_obs": norm["n_obs"],
        "n_meds": norm["n_meds"],
    }
    meta["rowcount_md"] = rowcount_markdown(meta)
    write_readme(root, meta, report, gen_sha)
    write_manifest(root, meta, report, gen_sha)
    # also copy the generator next to the outputs
    if gen_path != root / "make_synthetic_svd_dataset.py":
        (root / "make_synthetic_svd_dataset.py").write_bytes(gen_path.read_bytes())

    scratch = raw_dir.parent / "_scratch_project"
    if scratch.exists():
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)

    if not args.skip_zip:
        print("zipping stages", flush=True)
        zip_dir(raw_dir, raw_dir / "1_raw.zip")
        zip_dir(norm_dir, norm_dir / "2_normalized.zip")
        zip_dir(feat_dir, feat_dir / "3_features.zip")
        zip_dir(key_dir, key_dir / "4_answer_key.zip")

    print("done", flush=True)
    print(json.dumps({
        "n_subjects": meta["n_subjects"], "n_encounters": meta["n_encounters"],
        "n_invalid_dob": meta["n_invalid_dob"], "raw_mb": round(meta["raw_bytes"] / (1024 * 1024), 2),
        "singular_values_top10": report["singular_values_top10"],
        "variance_explained_top10": report["variance_explained_top10"],
        "gap_sv3_over_sv4": report["gap_sv3_over_sv4"],
        "procrustes_cosine_similarity": report["procrustes_cosine_similarity"],
        "principal_angles_degrees": report["principal_angles_degrees"],
        "generator_sha256": gen_sha,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
