#!/usr/bin/env python3
"""Build a three-layer (A/B/C) landmark-forward tensor dataset from a companion
extended raw stay, following contract v12 + v12.1 amendment A1.2–A1.7.

Does not modify stages 1–3, the three planted static factor scores, or existing
answer-key files. Static scores are read from 4_answer_key/latent_factors.csv.
Extended events are generated under 5_layered_tensor/ and driven by
severity_t(e,t) = s_e + d_e(t) (and a milder respiratory trajectory).
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
AUDIT_DIR = Path(__file__).resolve().parents[2] / "scripts" / "feature_sparsity_frequency_audit"
sys.path.insert(0, str(AUDIT_DIR))
sys.path.insert(0, str(HERE))

import feature_sparsity_frequency_audit as M  # noqa: E402
import synthetic_stage2_data as S  # noqa: E402
from make_synthetic_svd_dataset import (  # noqa: E402
    BETA, NOISE_SD, SVD_COLS, expected_weight_kg, hr_norm, rr_norm, sbp_norm, dbp_norm,
    normalize, engineer_features, document_features,
)

LAYERED_SEED = 20261115
LANDMARKS = tuple(range(24, 337, 24))  # 14 landmarks: 24..336
CORE_LANDMARKS = (24, 48, 72, 120, 168, 240, 336)
BIN_H = 6.0
HORIZON_H = 336.0
N_BINS = int(HORIZON_H / BIN_H)  # 56
BIN_EDGES = np.arange(N_BINS + 1, dtype=np.float64) * BIN_H
BIN_MIDS = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])
PRE_ICU_H = 6.0
SHAPE_CLASSES = ("stable", "improving", "worsening", "spike_recovery", "late_deterioration")
AMP_LOW, AMP_HIGH = 1.10, 2.40
RESP_AMP_SCALE = 0.40  # milder respiratory trajectory, same class
NOISE_STD_LAYER = 0.05  # not used in raw; raw uses NOISE_SD from the static generator
DEATH_TARGET = 0.075
FMT = S.FMT

MEAN_CHANNELS = [
    "heart_rate", "respiratory_rate", "sbp", "dbp", "map", "spo2", "fio2", "po2",
    "pf_ratio", "sf_ratio", "lactate", "creatinine", "temperature", "glucose_bedside",
    "sodium", "hemoglobin", "magnesium", "calcium_total", "urine_output", "fluid_intake",
]
# time-since last cap: 24 h default (A2.13) when the channel has no LOCF cap;
# triggered / locf families from the audit module if available.
TSL_CAP_H = 24.0


def ts(x: pd.Timestamp) -> str:
    return x.strftime(FMT)


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30.0, 30.0)))


def shape_raw(kind: str, t_mid: np.ndarray, onset: float | None = None) -> np.ndarray:
    """Unnormalized shape on a time grid in hours, u = t/336.

    stable:               0
    improving:            exp(-2 u)
    worsening:            logistic centered at 168 h
    spike_recovery:       Gaussian bump, peak at `onset` (hours), width 18 h
    late_deterioration:   0 until `onset`, then a logistic rise
    """
    t_mid = np.asarray(t_mid, dtype=float)
    u = t_mid / HORIZON_H
    if kind == "stable":
        return np.zeros_like(u)
    if kind == "improving":
        return np.exp(-2.0 * u)
    if kind == "worsening":
        return 1.0 / (1.0 + np.exp(-8.0 * (u - 0.5)))
    if kind == "spike_recovery":
        peak = 48.0 if onset is None else float(onset)
        return np.exp(-0.5 * ((t_mid - peak) / 18.0) ** 2)
    if kind == "late_deterioration":
        on = 120.0 if onset is None else float(onset)
        return 1.0 / (1.0 + np.exp(-6.0 * ((t_mid - on) / 48.0)))
    raise ValueError(kind)


def unit_norm_center(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    if not np.any(np.abs(v) > 0):
        return v
    v = v - v.mean()
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def plant_trajectories(n: int, rng: np.random.Generator) -> dict:
    cls = rng.choice(np.array(SHAPE_CLASSES), size=n, replace=True)
    amp = np.zeros(n, dtype=float)
    onset = np.full(n, np.nan)
    D = np.zeros((n, N_BINS), dtype=float)
    non = cls != "stable"
    amp[non] = rng.uniform(AMP_LOW, AMP_HIGH, size=int(non.sum()))
    # onsets
    sp = cls == "spike_recovery"
    onset[sp] = rng.uniform(12.0, 96.0, size=int(sp.sum()))
    lt = cls == "late_deterioration"
    onset[lt] = rng.uniform(48.0, 200.0, size=int(lt.sum()))
    for i in range(n):
        k = cls[i]
        if k == "stable":
            continue
        on = None if not np.isfinite(onset[i]) else float(onset[i])
        D[i] = amp[i] * unit_norm_center(shape_raw(k, BIN_MIDS, onset=on))
    return {"class": cls, "amplitude": amp, "onset_h": onset, "D": D}


def d_at_time(D: np.ndarray, t_h: np.ndarray) -> np.ndarray:
    """Piecewise-constant d_e(t) on 6 h bins (value of the containing bin). t outside [0,336) uses the edge bin."""
    idx = np.floor(np.clip(t_h, 0, HORIZON_H - 1e-9) / BIN_H).astype(int)
    # D is (n, 56); t_h can be scalar or (n,)
    if np.ndim(t_h) == 0:
        return D[:, idx] if D.ndim == 2 and np.ndim(idx) == 0 else D[idx]
    return D[np.arange(len(t_h)), idx] if D.ndim == 2 and len(t_h) == len(D) else D[:, idx]


def admission_source(edreg, admit) -> str:
    if edreg is pd.NaT or (isinstance(edreg, float) and not np.isfinite(edreg)):
        return "unknown"
    if pd.isna(edreg):
        return "unknown"
    return "ed" if edreg <= admit else "unknown"


# =============================================================================
# Extended raw event generation, driven by severity_t / respiratory_t
# =============================================================================
def generate_extended_raw(root: Path, ext_raw: Path, plant: dict, rng: np.random.Generator) -> dict:
    ext_raw.mkdir(parents=True, exist_ok=True)
    # Copy identity / map / split tables (new copies; existing 1_raw is not touched).
    src = root / "1_raw"
    for name in (
        "PATIENTS.csv", "ADMISSIONS.csv", "ICUSTAYS.csv", "SURGERY_INFO.csv",
        "canonical_icu_encounters.csv", "encounter_harmonization_fields.csv",
        "subject_splits.csv", "D_ITEMS.csv", "D_LABITEMS.csv",
        "OUTPUTEVENTS.csv", "INPUTEVENTS.csv", "MICROBIOLOGYEVENTS.csv",
        "feature_channel_map_v1.csv", "drug_group_map_v1.csv", "specimen_map_v1.csv",
        "organism_map_v1.csv", "neuro_admitting_dx_icd10_v1.csv",
        "surgery_complexity_map_v1.csv",
    ):
        shutil.copy2(src / name, ext_raw / name)

    factors = pd.read_csv(root / "4_answer_key" / "latent_factors.csv", dtype={"encounter_id": str})
    enc = pd.read_csv(root / "2_normalized" / "encounters.csv", dtype={"encounter_id": str})
    enc = enc.merge(factors, on="encounter_id", how="left", suffixes=("", "_fac"))
    # keep factor file order for D alignment: reindex enc to factors order
    enc = enc.set_index("encounter_id").loc[factors["encounter_id"]].reset_index()
    assert len(enc) == len(factors)
    n = len(enc)
    D = plant["D"]
    # age in years for norms
    age_y = enc["age_days"].to_numpy(dtype=float) / 365.25
    age_y = np.where(np.isfinite(age_y) & (age_y > 0), age_y, 9.0)  # age_invalid: center of the z_age map
    z_sev = enc["z_severity"].to_numpy(dtype=float)
    z_resp = enc["z_respiratory"].to_numpy(dtype=float)
    z_age = enc["z_age"].to_numpy(dtype=float)
    los = enc["los_hours"].to_numpy(dtype=float)
    ids = enc["encounter_id"].to_numpy()
    sids = enc["SUBJECT_ID"].astype(str).to_numpy()
    hadms = enc["HADM_ID"].astype(str).to_numpy()
    t0 = pd.to_datetime(enc["INTIME"])
    sex = enc["sex"].to_numpy()

    chart, lab, rx = [], [], []
    H = lambda t0i, h: (t0i + pd.Timedelta(hours=float(h))).floor("s")
    vaso_names = [d for d in S.DRUGS if d[2] == "vasoactive"]
    other_drugs = [d for d in S.DRUGS if d[2] != "vasoactive"]

    # Precompute weights (charted near t0)
    wt = expected_weight_kg(age_y) * np.exp(BETA["weight_log_age"] * z_age + rng.normal(0, NOISE_SD["weight_log"], size=n))
    wt = np.clip(wt, 0.5, 140.0)

    print("  generating extended events ...", flush=True)
    for i in range(n):
        if i % 4000 == 0:
            print(f"    encounter {i}/{n}", flush=True)
        t0i = t0.iloc[i]
        sp = float(los[i])
        sid, h, ay = sids[i], hadms[i], float(age_y[i])
        # Encounter-level noise intercepts (independent of factors)
        temp_mu = 37.0 + rng.normal(0, 0.35)
        glu_mu = 5.4 + rng.normal(0, 0.6)
        na_mu = 138.0 + rng.normal(0, 2.0)
        hb_mu = 115.0 + rng.normal(0, 10.0)
        mg_mu = 0.85 + rng.normal(0, 0.08)
        ca_mu = 2.30 + rng.normal(0, 0.12)
        urine_mu = 40.0 + rng.normal(0, 15.0)
        intake_mu = 60.0 + rng.normal(0, 20.0)

        chart.append([sid, h, "1014", ts(H(t0i, 1.0)), f"{float(wt[i]):.2f}", "kg"])

        # Time grid: every 2 h on [-6, 24), every 4 h on [24, min(LOS, 336))
        early = np.arange(-PRE_ICU_H, min(sp, 24.0), 2.0)
        late = np.arange(24.0, min(sp, HORIZON_H), 4.0) if sp > 24 else np.array([])
        times = np.unique(np.concatenate([early, late])) if len(late) else early
        # d and respiratory deviation at each time (bin-constant)
        t_clip = np.clip(times, 0.0, HORIZON_H - 1e-9)
        t_idx = np.floor(t_clip / BIN_H).astype(int)
        t_idx = np.clip(t_idx, 0, N_BINS - 1)
        d_t = np.where(times >= 0, D[i, t_idx], 0.0)  # pre-t0 uses static s_e only
        sev_t = z_sev[i] + d_t
        resp_t = z_resp[i] + RESP_AMP_SCALE * d_t  # milder, same shape

        hr_mu = hr_norm(ay) + BETA["hr_sev"] * sev_t
        rr_mu = rr_norm(ay) + BETA["rr_resp"] * resp_t
        sbp_mu = sbp_norm(ay) + BETA["sbp_sev"] * sev_t
        dbp_mu = dbp_norm(ay) + BETA["dbp_sev"] * sev_t
        spo2_mu = 97.2 + BETA["spo2_resp"] * resp_t
        fio2_mu = 30.0 + BETA["fio2_resp"] * resp_t
        po2_mu = 95.0 + BETA["po2_resp"] * resp_t
        lact_mu = 1.15 + BETA["lactate_sev"] * sev_t
        creat_mu = 42.0 + BETA["creatinine_sev"] * sev_t

        for k, sev, resp, hrm, rrm, sbm, dbm, spm, fim, pom, lam, crm in zip(
            times, sev_t, resp_t, hr_mu, rr_mu, sbp_mu, dbp_mu, spo2_mu, fio2_mu, po2_mu, lact_mu, creat_mu
        ):
            hr = float(np.clip(hrm + rng.normal(0, NOISE_SD["hr"]), 30, 300))
            chart.append([sid, h, "1003", ts(H(t0i, k)), f"{hr:.0f}", "bpm"])
            if rng.random() < 0.15:
                chart.append([sid, h, "1002", ts(H(t0i, k)), f"{float(np.clip(hrm + rng.normal(0, NOISE_SD['hr']), 30, 300)):.0f}", "bpm"])
            if rng.random() < 0.85:
                chart.append([sid, h, "1004", ts(H(t0i, k)), f"{float(np.clip(rrm + rng.normal(0, NOISE_SD['rr']), 5, 120)):.0f}", "bpm"])
            if rng.random() < 0.88:
                sv = float(np.clip(sbm + rng.normal(0, NOISE_SD["sbp"]), 30, 250))
                dv = float(np.clip(dbm + rng.normal(0, NOISE_SD["dbp"]), 10, 200))
                if dv >= sv - 2:
                    dv = max(10.0, sv - 8.0)
                chart.append([sid, h, "1016", ts(H(t0i, k)), f"{sv:.0f}", "mmHg"])
                chart.append([sid, h, "1015", ts(H(t0i, k)), f"{dv:.0f}", "mmHg"])
            if rng.random() < 0.88:
                chart.append([sid, h, "1006", ts(H(t0i, k)), f"{float(np.clip(spm + rng.normal(0, NOISE_SD['spo2']), 30, 100)):.0f}", "%"])
            if rng.random() < 0.50:
                chart.append([sid, h, "1001", ts(H(t0i, k)), f"{float(np.clip(temp_mu + rng.normal(0, NOISE_SD['temp']), 30, 43)):.1f}", "°C"])
            if rng.random() < 0.40:
                chart.append([sid, h, "1008", ts(H(t0i, k)), f"{float(np.clip(urine_mu + rng.normal(0, NOISE_SD['urine']), 0, 3000)):.0f}", "ml"])
            if rng.random() < 0.40:
                chart.append([sid, h, "1009", ts(H(t0i, k + 0.2)), f"{float(np.clip(intake_mu + rng.normal(0, NOISE_SD['intake']), 0, 5000)):.0f}", "ml"])
            if rng.random() < 0.20:
                chart.append([sid, h, "1011", ts(H(t0i, k)), f"{float(np.clip(glu_mu + rng.normal(0, NOISE_SD['glucose']), 0.5, 50)):.1f}", "mmol/L"])

        # Labs: a couple in the pre-ICU 6 h window (for Layer A), several in 0-24 h, then ~daily to OUT
        lab_times = list(rng.uniform(-PRE_ICU_H + 0.2, -0.2, size=int(rng.integers(0, 2))))
        lab_times += list(rng.uniform(0.2, min(sp, 22.5), size=int(rng.integers(2, 5))))
        day = 24.0
        while day < min(sp, HORIZON_H) - 1:
            lab_times.append(float(day + rng.uniform(0.0, 6.0)))
            day += float(rng.uniform(18.0, 30.0))
        for k in lab_times:
            if k >= sp:
                continue
            t_idx = int(np.clip(np.floor(np.clip(k, 0, HORIZON_H - 1e-9) / BIN_H), 0, N_BINS - 1))
            dtk = D[i, t_idx] if k >= 0 else 0.0
            sev = z_sev[i] + dtk
            resp = z_resp[i] + RESP_AMP_SCALE * dtk
            cr = float(np.clip(42.0 + BETA["creatinine_sev"] * sev + rng.normal(0, NOISE_SD["creatinine"]), 5, 2000))
            lab.append([sid, h, "5041" if rng.random() < 0.75 else "5032", ts(H(t0i, k)), f"{cr:.0f}", "μmol/L"])
            lab.append([sid, h, "5230", ts(H(t0i, k)), f"{float(np.clip(na_mu + rng.normal(0, NOISE_SD['sodium']), 100, 190)):.0f}", "mmol/L"])
            lac = float(np.clip(1.15 + BETA["lactate_sev"] * sev + rng.normal(0, NOISE_SD["lactate"]), 0.1, 30))
            lab.append([sid, h, "5227", ts(H(t0i, k)), f"{lac:.2f}", "mmol/L"])
            if rng.random() < 0.55:
                lab.append([sid, h, "5099", ts(H(t0i, k)), f"{float(np.clip(hb_mu + rng.normal(0, NOISE_SD['hemoglobin']), 20, 250)):.0f}", "g/L"])
            if rng.random() < 0.40:
                lab.append([sid, h, "5059", ts(H(t0i, k)), f"{float(np.clip(mg_mu + rng.normal(0, NOISE_SD['magnesium']), 0.1, 5)):.2f}", "mmol/L"])
            if rng.random() < 0.40:
                lab.append([sid, h, "5215", ts(H(t0i, k)), f"{float(np.clip(ca_mu + rng.normal(0, NOISE_SD['calcium']), 0.5, 5)):.2f}", "mmol/L"])
            # blood gas / S/F partners
            if rng.random() < (0.70 if k < 24 else 0.45):
                fi = float(np.clip(30.0 + BETA["fio2_resp"] * resp + rng.normal(0, NOISE_SD["fio2"]), 21, 100))
                lab.append([sid, h, "5222", ts(H(t0i, k)), f"{fi:.0f}", "%"])
                if rng.random() < 0.62:
                    po = float(np.clip(95.0 + BETA["po2_resp"] * resp + rng.normal(0, NOISE_SD["po2"]), 10, 700))
                    lab.append([sid, h, "5239" if rng.random() < 0.7 else "5244", ts(H(t0i, k)), f"{po:.0f}", "mmHg"])
                    lab.append([sid, h, "5243", ts(H(t0i, k)), f"{po / (fi / 100.0) + rng.normal(0, 8):.0f}", "mmHg"])
                else:
                    spo = float(np.clip(97.2 + BETA["spo2_resp"] * resp + rng.normal(0, NOISE_SD["spo2"]), 30, 100))
                    if spo > 97 and resp > 0:
                        spo = float(rng.uniform(90, 97))
                    chart.append([sid, h, "1006", ts(H(t0i, k) + pd.Timedelta(minutes=int(rng.choice([-3, 0, 2])))), f"{spo:.0f}", "%"])

        # Vasoactive exposure tracks severity_t over the stay: start a course when a mid-stay
        # draw exceeds a logit threshold; number of agents rises with severity.
        # Evaluate at a few checkpoints.
        checks = np.arange(1.0, min(sp, HORIZON_H), 8.0)
        if len(checks) == 0:
            checks = np.array([min(6.0, max(1.0, sp / 2))])
        cidx = np.clip(np.floor(np.clip(checks, 0, HORIZON_H - 1e-9) / BIN_H).astype(int), 0, N_BINS - 1)
        sev_c = z_sev[i] + D[i, cidx]
        p = sigmoid(BETA["vaso_logit_intercept"] + BETA["vaso_logit_sev"] * sev_c)
        started = False
        for ck, pk, svk in zip(checks, p, sev_c):
            if rng.random() < pk:
                n_agents = 1 + int(rng.random() < sigmoid(-0.6 + 1.1 * svk)) + int(rng.random() < sigmoid(-1.4 + 1.1 * svk))
                n_agents = int(np.clip(n_agents, 1, len(vaso_names)))
                chosen = [vaso_names[j] for j in rng.choice(len(vaso_names), size=n_agents, replace=False)]
                dur = float(rng.uniform(8.0, 24.0))
                for dn, de, g, sg, cr in chosen:
                    stt = H(t0i, float(ck))
                    ent = H(t0i, float(min(sp, ck + dur)))
                    rx.append([sid, h, ts(stt), ts(ent), dn, de, "mg"])
                started = True
                break
        for dn, de, g, sg, cr in other_drugs:
            if rng.random() < 0.22:
                k = float(rng.uniform(-2.0, max(1.0, min(sp, 80.0))))
                stt = H(t0i, k)
                start = ts(stt) if rng.random() > 0.1 else stt.strftime("%Y-%m-%d")
                end = ts(H(t0i, k + float(rng.uniform(6, 36)))) if rng.random() < 0.35 else ""
                rx.append([sid, h, start, end, dn, de, "mg"])

    ch = pd.DataFrame(chart, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"])
    ch["VALUENUM"] = pd.to_numeric(ch["VALUE"], errors="coerce").map(lambda v: "" if pd.isna(v) else f"{v:g}")
    ch["ICUSTAY_ID"] = ""
    ch.to_csv(ext_raw / "CHARTEVENTS.csv", index=False)
    lb = pd.DataFrame(lab, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"])
    lb["VALUENUM"] = pd.to_numeric(lb["VALUE"], errors="coerce").map(lambda v: "" if pd.isna(v) else f"{v:g}")
    lb.to_csv(ext_raw / "LABEVENTS.csv", index=False)
    pd.DataFrame(rx, columns=["SUBJECT_ID", "HADM_ID", "STARTDATE", "ENDDATE", "DRUG_NAME", "DRUG_NAME_EN", "DOSE_UNIT_RX"]).to_csv(
        ext_raw / "PRESCRIPTIONS.csv", index=False)
    print(f"    chart={len(ch)} lab={len(lb)} rx={len(rx)}", flush=True)
    return {"enc": enc, "n": n, "ids": ids, "wt": wt, "chart_n": len(ch), "lab_n": len(lb), "rx_n": len(rx)}


# =============================================================================
# Outcome model (outcomes_source; does not modify ADMISSIONS / ICUSTAYS)
# =============================================================================
def generate_outcomes(enc: pd.DataFrame, plant: dict, rng: np.random.Generator) -> pd.DataFrame:
    """Documented generative model. ICU INTIME/OUTTIME stay the existing ones.

    Death, post-ICU discharge, and same-HADM readmission are written only here.
    Outcome 1 levels 0–7 for survivors are therefore fixed by existing LOS;
    death (level 8 / Outcome 2 level 4) and Outcome 2 levels 0–3 depend on the
    latent factors and trajectories via the logits below.

    logit P(died in hospital) = a0 + 0.85 s_e + 0.55 end_severity + 0.25 z_resp
                                + 0.45 1{late_deterioration} + 0.15 mean(d_e)
    a0 is calibrated so the cohort mortality is about 7.5% (PIC-like 7–8%).
    Among deaths, P(death in ICU, at last OUTTIME) = sigmoid(-0.2 + 0.4 end_severity),
    otherwise death is after ICU discharge.

    logit P(same-HADM ICU readmission >24 h after OUT) =
        b0 + 0.55 s_e + 0.40 mean(d_e) + 0.20 z_resp,  b0 chosen for ~8% of survivors.
    post_icu_days = clip(exp(1.15 + 0.28 s_e + 0.20 end_severity + 0.15 z_resp + N(0,0.55)), 0.4, 45)
    """
    n = len(enc)
    z_sev = enc["z_severity"].to_numpy(dtype=float)
    z_resp = enc["z_respiratory"].to_numpy(dtype=float)
    los = enc["los_hours"].to_numpy(dtype=float)
    D = plant["D"]
    # end-of-stay severity and mean d over bins that start before OUT
    end_sev = np.empty(n)
    mean_d = np.empty(n)
    for i in range(n):
        last = int(np.clip(np.floor((min(los[i], HORIZON_H) - 1e-9) / BIN_H), 0, N_BINS - 1))
        end_sev[i] = z_sev[i] + D[i, last]
        valid = max(1, last + 1)
        mean_d[i] = D[i, :valid].mean()
    late = (plant["class"] == "late_deterioration").astype(float)
    lin = 0.85 * z_sev + 0.55 * end_sev + 0.25 * z_resp + 0.45 * late + 0.15 * mean_d
    # calibrate intercept
    def mort(a0):
        return float(sigmoid(a0 + lin).mean())
    lo, hi = -6.0, 2.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if mort(mid) > DEATH_TARGET:
            hi = mid
        else:
            lo = mid
    a0 = 0.5 * (lo + hi)
    p_death = sigmoid(a0 + lin)
    died = rng.random(n) < p_death
    # ICU vs post-ICU death
    p_icu_death = sigmoid(-0.2 + 0.4 * end_sev)
    death_in_icu = died & (rng.random(n) < p_icu_death)
    t0 = pd.to_datetime(enc["INTIME"])
    out = pd.to_datetime(enc["OUTTIME"])
    death_time = np.array([pd.NaT] * n, dtype=object)
    death_loc = np.full(n, "", dtype=object)
    for i in range(n):
        if not died[i]:
            continue
        if death_in_icu[i]:
            # death at last OUTTIME (stay bounds unchanged; eligible iff existing OUT > t0+L)
            death_time[i] = out.iloc[i]
            death_loc[i] = "icu_at_outtime"
        else:
            extra = float(np.clip(rng.lognormal(mean=2.2, sigma=0.8), 2.0, 400.0))  # hours after ICU out
            death_time[i] = out.iloc[i] + pd.Timedelta(hours=extra)
            death_loc[i] = "post_icu"

    # Readmission among non-ICU-deaths
    lin_r = 0.55 * z_sev + 0.40 * mean_d + 0.20 * z_resp
    def rate(b0):
        return float(sigmoid(b0 + lin_r)[~death_in_icu].mean()) if (~death_in_icu).any() else 0.0
    lo, hi = -6.0, 2.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if rate(mid) > 0.08:
            hi = mid
        else:
            lo = mid
    b0 = 0.5 * (lo + hi)
    p_re = sigmoid(b0 + lin_r)
    readmit = (~death_in_icu) & (rng.random(n) < p_re)
    # if they die post-ICU, keep readmission only if it occurs before death
    readmit_in = np.array([pd.NaT] * n, dtype=object)
    readmit_out = np.array([pd.NaT] * n, dtype=object)
    readmit_gap_h = np.full(n, np.nan)
    for i in range(n):
        if not readmit[i]:
            continue
        gap = 24.0 + float(rng.lognormal(mean=3.0, sigma=0.7))  # hours after OUT, always >24
        gap = float(np.clip(gap, 24.01, 400.0))
        rin = out.iloc[i] + pd.Timedelta(hours=gap)
        if died[i] and death_loc[i] == "post_icu" and pd.notna(death_time[i]) and rin >= death_time[i]:
            readmit[i] = False
            continue
        rout = rin + pd.Timedelta(hours=float(np.clip(rng.lognormal(3.5, 0.6), 6.0, 200.0)))
        if died[i] and death_loc[i] == "post_icu" and pd.notna(death_time[i]) and rout > death_time[i]:
            rout = death_time[i]
        readmit_in[i] = rin
        readmit_out[i] = rout
        readmit_gap_h[i] = gap

    # Post-ICU hospital stay
    post_days = np.clip(
        np.exp(1.15 + 0.28 * z_sev + 0.20 * end_sev + 0.15 * z_resp + rng.normal(0, 0.55, size=n)),
        0.4, 45.0,
    )
    disch = np.array([pd.NaT] * n, dtype=object)
    post_h = np.full(n, np.nan)
    for i in range(n):
        if death_in_icu[i]:
            disch[i] = out.iloc[i]  # hospital death at ICU out
            post_h[i] = 0.0
            continue
        hours = float(post_days[i] * 24.0)
        if readmit[i] and pd.notna(readmit_out[i]):
            cand = readmit_out[i] + pd.Timedelta(hours=float(rng.uniform(12, 120)))
            base = out.iloc[i] + pd.Timedelta(hours=hours)
            dt = max(cand, base)
        else:
            dt = out.iloc[i] + pd.Timedelta(hours=hours)
        if died[i] and death_loc[i] == "post_icu" and pd.notna(death_time[i]):
            if dt > death_time[i]:
                dt = death_time[i]
        disch[i] = dt
        post_h[i] = (dt - out.iloc[i]).total_seconds() / 3600.0

    out_df = pd.DataFrame({
        "encounter_id": enc["encounter_id"],
        "SUBJECT_ID": enc["SUBJECT_ID"],
        "HADM_ID": enc["HADM_ID"],
        "split": enc["split"],
        "INTIME": enc["INTIME"],
        "OUTTIME": enc["OUTTIME"],
        "los_hours": los,
        "z_severity": z_sev,
        "z_respiratory": z_resp,
        "end_severity": end_sev,
        "mean_d": mean_d,
        "shape_class": plant["class"],
        "died_in_hospital": died.astype(int),
        "death_location": death_loc,
        "death_time": ["" if pd.isna(x) else ts(pd.Timestamp(x)) for x in death_time],
        "readmission": readmit.astype(int),
        "readmit_intime": ["" if pd.isna(x) else ts(pd.Timestamp(x)) for x in readmit_in],
        "readmit_outtime": ["" if pd.isna(x) else ts(pd.Timestamp(x)) for x in readmit_out],
        "readmit_gap_hours": readmit_gap_h,
        "hospital_dischtime": ["" if pd.isna(x) else ts(pd.Timestamp(x)) for x in disch],
        "post_icu_hours": post_h,
        "logit_death_a0": a0,
        "logit_readmit_b0": b0,
    })
    # stash intercepts on the frame for the answer key
    out_df.attrs["a0"] = a0
    out_df.attrs["b0"] = b0
    out_df.attrs["death_rate"] = float(died.mean())
    out_df.attrs["readmit_rate"] = float(readmit.mean())
    out_df.attrs["lin_death"] = lin
    return out_df


def outcome1_level(remaining_h: float, died_after_L: bool) -> int:
    if died_after_L:
        return 8
    if not np.isfinite(remaining_h):
        return -1
    # searchsorted on (24,48,72,120,168,240,336), side=left: R<=24 -> 0, ..., R>336 -> 7
    return int(np.searchsorted(np.array(M.OUTCOME_EDGES, dtype=float), remaining_h, side="left"))


def outcome2_level(died_after_L: bool, readmit: bool, post_icu_hours: float) -> int:
    if died_after_L:
        return 4
    if readmit:
        return 3
    if not np.isfinite(post_icu_hours):
        return -1  # unassignable
    days = post_icu_hours / 24.0
    if days >= 14:
        return 2
    if days >= 7:
        return 1
    return 0


def inherited_3_0_level(total_icu_h: float, died_after_L: bool) -> int:
    """Total ICU stay from t0, same cut points, death = 8. Collapses at late landmarks (A1.1)."""
    if died_after_L:
        return 8
    return int(np.searchsorted(np.array(M.OUTCOME_EDGES, dtype=float), total_icu_h, side="left"))


def build_landmark_outcomes(enc: pd.DataFrame, oc: pd.DataFrame) -> pd.DataFrame:
    t0 = pd.to_datetime(enc["INTIME"])
    outt = pd.to_datetime(enc["OUTTIME"])
    los = enc["los_hours"].to_numpy(dtype=float)
    death = pd.to_datetime(oc["death_time"], errors="coerce")
    died = oc["died_in_hospital"].to_numpy(dtype=int)
    readmit = oc["readmission"].to_numpy(dtype=int)
    post_h = oc["post_icu_hours"].to_numpy(dtype=float)
    rows = []
    for i in range(len(enc)):
        dt = death.iloc[i]
        for L in LANDMARKS:
            # Eligibility: still in ICU at L (existing OUTTIME > t0+L). If an ICU death
            # is recorded at OUTTIME, this is the same rule. A death at or before t0+L
            # would also exclude, but ICU deaths are at OUTTIME so it coincides.
            if not (los[i] > L):
                continue
            if pd.notna(dt) and dt <= t0.iloc[i] + pd.Timedelta(hours=float(L)):
                continue
            died_after = bool(died[i]) and (pd.isna(dt) or dt > t0.iloc[i] + pd.Timedelta(hours=float(L)))
            rem = los[i] - L
            o1 = outcome1_level(rem, died_after)
            o2 = outcome2_level(died_after, bool(readmit[i]), post_h[i])
            inh = inherited_3_0_level(los[i], died_after)
            rows.append({
                "encounter_id": enc["encounter_id"].iloc[i],
                "SUBJECT_ID": enc["SUBJECT_ID"].iloc[i],
                "split": enc["split"].iloc[i],
                "landmark_h": L,
                "eligible": 1,
                "n_bins": int(L / BIN_H),
                "remaining_icu_hours": rem,
                "outcome1_remaining_icu": o1,
                "outcome2_post_icu": o2,
                "in_hospital_mortality": int(died_after),
                "outcome_inherited_3_0_ordinal": inh if L in (24, 48, 72) else np.nan,
                "outcome2_unassignable": int(o2 < 0),
            })
    return pd.DataFrame(rows)


# =============================================================================
# Layers
# =============================================================================
def build_layer_a(root: Path, ext_norm: Path, enc: pd.DataFrame, feat_b: pd.DataFrame) -> pd.DataFrame:
    adms = pd.read_csv(root / "1_raw" / "ADMISSIONS.csv", dtype=str)
    surg = pd.read_csv(root / "1_raw" / "SURGERY_INFO.csv", dtype=str)
    adms["SUBJECT_ID"] = M.vnorm_id(adms["SUBJECT_ID"])
    adms["HADM_ID"] = M.vnorm_id(adms["HADM_ID"])
    adms["ADMITTIME"] = M.fast_time(adms["ADMITTIME"])
    adms["EDREGTIME"] = M.fast_time(adms["EDREGTIME"])
    adms["dx_norm"] = [M.norm_icd(x)[0] for x in adms.get("ICD10_CODE_CN", pd.Series([""] * len(adms)))]
    adms["dx_3char"] = adms["dx_norm"].str[:3]
    adms["dx_domain"] = [M.icd_domain(c) for c in adms["dx_3char"]]
    adms["admission_source"] = [
        admission_source(ed, ad) for ed, ad in zip(adms["EDREGTIME"], adms["ADMITTIME"])
    ]

    e = enc.copy()
    e["SUBJECT_ID"] = M.vnorm_id(e["SUBJECT_ID"].astype(str))
    e["HADM_ID"] = M.vnorm_id(e["HADM_ID"].astype(str))
    e["t0"] = M.fast_time(e["INTIME"])
    e = e.merge(adms[["SUBJECT_ID", "HADM_ID", "dx_norm", "dx_3char", "dx_domain", "admission_source", "ADMITTIME"]],
                on=["SUBJECT_ID", "HADM_ID"], how="left")

    # Surgery timing uses anesthesia end (A2.1). Layer A recent: T in [t0-24h, t0+6h].
    surg["SUBJECT_ID"] = M.vnorm_id(surg["SUBJECT_ID"])
    surg["HADM_ID"] = M.vnorm_id(surg["HADM_ID"])
    # A2.1: T = ANES_END_TIME; fall back to SURGERY_END_TIME then BEGIN
    T = M.fast_time(surg["ANES_END_TIME"])
    T = T.fillna(M.fast_time(surg["SURGERY_END_TIME"])).fillna(M.fast_time(surg["SURGERY_BEGIN_TIME"]))
    surg = surg.assign(T=T).dropna(subset=["T"])
    e = e.merge(surg[["SUBJECT_ID", "HADM_ID", "T", "SURGERY_NAME"]], on=["SUBJECT_ID", "HADM_ID"], how="left")
    # multiple surgeries: flags
    def flags(g):
        t0 = g["t0"].iloc[0]
        Ts = g["T"]
        recent = ((Ts >= t0 - pd.Timedelta(hours=24)) & (Ts <= t0 + pd.Timedelta(hours=6))).any()
        remote = (Ts < t0 - pd.Timedelta(hours=24)).any()
        return pd.Series({"recent_surgery": int(bool(recent)), "remote_surgery_same_hadm": int(bool(remote))})
    # if merge duplicated rows, aggregate
    key = ["encounter_id"]
    extra = e.groupby(key, sort=False).apply(flags, include_groups=False).reset_index()
    e0 = e.drop_duplicates("encounter_id")
    a = e0.merge(extra, on="encounter_id", how="left")
    a["recent_surgery"] = a["recent_surgery"].fillna(0).astype(int)
    a["remote_surgery_same_hadm"] = a["remote_surgery_same_hadm"].fillna(0).astype(int)

    # Pre-ICU lab block: last released value in [t0-6h, t0) from extended observations
    obs = pd.read_csv(ext_norm / "observations_long.csv",
                      usecols=["encounter_id", "t_hours", "channel", "value", "status"],
                      dtype={"encounter_id": str, "channel": str, "status": str})
    obs = obs[(obs["status"] == "measurement") & obs["value"].notna()].copy()
    obs["t_hours"] = pd.to_numeric(obs["t_hours"], errors="coerce")
    pre = obs[(obs["t_hours"] >= -PRE_ICU_H) & (obs["t_hours"] < 0)]
    lab_ch = ["lactate", "creatinine", "sodium", "hemoglobin", "magnesium", "calcium_total", "po2", "fio2", "pf_ratio"]
    last = (pre[pre["channel"].isin(lab_ch)]
            .sort_values(["encounter_id", "channel", "t_hours"])
            .groupby(["encounter_id", "channel"], sort=False)["value"].last()
            .unstack("channel"))
    last = last.reindex(columns=lab_ch)
    last.columns = [f"pre_icu_{c}_last" for c in last.columns]
    mask = last.notna().astype(np.int8)
    mask.columns = [c + "_observed" for c in last.columns]
    pre_block = last.join(mask)
    pre_block.index.name = "encounter_id"
    pre_block = pre_block.reset_index()

    wcol = feat_b[["encounter_id", "weight_kg", "weight_source"]] if "weight_source" in feat_b.columns else feat_b[["encounter_id", "weight_kg"]]
    if "age_band" in a.columns:
        # name may be duplicated after the surgery merge; take the first
        band = a["age_band"]
        if isinstance(band, pd.DataFrame):
            band = band.iloc[:, 0]
        a = a.loc[:, ~a.columns.duplicated()].copy()
        a["age_band"] = band.to_numpy()
        a["age_invalid"] = (a["age_band"].astype(str) == "age_invalid").astype(np.int8)
    cols = ["encounter_id", "SUBJECT_ID", "HADM_ID", "split", "sex", "first_careunit",
            "age_days", "age_band", "dob_status", "age_invalid",
            "dx_3char", "dx_domain", "admission_source", "recent_surgery", "remote_surgery_same_hadm",
            "neuro_dx", "los_hours"]
    keep = [c for c in cols if c in a.columns]
    out = a[keep].drop_duplicates("encounter_id")
    if "age_invalid" not in out.columns:
        out["age_invalid"] = (out["age_band"].astype(str) == "age_invalid").astype(np.int8)
    out = out.merge(wcol, on="encounter_id", how="left")
    out = out.merge(pre_block, on="encounter_id", how="left")
    for c in [f"pre_icu_{ch}_last_observed" for ch in lab_ch]:
        if c in out.columns:
            out[c] = out[c].fillna(0).astype(np.int8)
    return out


def build_layer_c(ext_norm: Path, enc: pd.DataFrame) -> dict:
    print("  Layer C binning ...", flush=True)
    obs = pd.read_csv(ext_norm / "observations_long.csv",
                      usecols=["encounter_id", "t_hours", "channel", "value", "status"],
                      dtype={"encounter_id": str, "channel": str, "status": str})
    meds = pd.read_csv(ext_norm / "medications_long.csv", dtype={"encounter_id": str})
    ids = enc["encounter_id"].astype(str).to_numpy()
    id_to_i = {e: i for i, e in enumerate(ids)}
    n = len(ids)
    los = enc["los_hours"].to_numpy(dtype=float)

    # feature layout
    feat_names = [f"{ch}_mean" for ch in MEAN_CHANNELS]
    feat_names += [f"{ch}_time_since_last_h" for ch in MEAN_CHANNELS]
    feat_names += [f"{ch}_last_carried" for ch in MEAN_CHANNELS]
    feat_names += [f"{g}_active" for g in M.MED_GROUPS]
    feat_names += ["vasoactive_max_concurrent", "psofa_resp"]
    F = len(feat_names)
    name_i = {k: j for j, k in enumerate(feat_names)}

    X = np.full((n, N_BINS, F), np.nan, dtype=np.float32)
    mask_obs = np.zeros((n, N_BINS, F), dtype=np.uint8)
    # valid bin: bin start < LOS (still in ICU when the bin opens)
    mask_valid = (BIN_EDGES[:-1][None, :] < los[:, None]).astype(np.uint8)

    obs = obs[(obs["status"] == "measurement") & obs["value"].notna()].copy()
    obs["t_hours"] = pd.to_numeric(obs["t_hours"], errors="coerce")
    obs = obs[obs["t_hours"].notna() & (obs["t_hours"] >= 0) & (obs["t_hours"] < HORIZON_H)]
    obs["t_step"] = np.floor(obs["t_hours"].to_numpy() / BIN_H).astype(np.int16)
    obs = obs[(obs["t_step"] >= 0) & (obs["t_step"] < N_BINS)]
    obs["enc_i"] = obs["encounter_id"].map(id_to_i)
    obs = obs[obs["enc_i"].notna()].copy()
    obs["enc_i"] = obs["enc_i"].astype(np.int32)

    # means
    need = set(MEAN_CHANNELS)
    meas = obs[obs["channel"].isin(need)]
    g = meas.groupby(["enc_i", "t_step", "channel"], sort=False)["value"].mean()
    for (ei, ts, ch), val in zip(g.index, g.to_numpy(dtype=float)):
        j = name_i[f"{ch}_mean"]
        ei, ts = int(ei), int(ts)
        X[ei, ts, j] = val
        mask_obs[ei, ts, j] = 1

    # last observation time per channel for TSL / last-carried
    last_t = meas.groupby(["enc_i", "t_step", "channel"], sort=False)["t_hours"].max()
    # walk bins in time for each (encounter, channel)
    # build dict (ei, ch) -> list of (step, t_hours)
    from collections import defaultdict
    ev = defaultdict(list)
    for (ei, ts, ch), tlast in zip(last_t.index, last_t.to_numpy(dtype=float)):
        ev[(int(ei), ch)].append((int(ts), float(tlast)))
    for (ei, ch), items in ev.items():
        items.sort()
        prev_t = None
        it = 0
        for b in range(N_BINS):
            observed = False
            t_obs = None
            while it < len(items) and items[it][0] == b:
                observed = True
                t_obs = items[it][1]
                it += 1
            if observed:
                prev_t = t_obs
                X[ei, b, name_i[f"{ch}_last_carried"]] = 0.0
                mask_obs[ei, b, name_i[f"{ch}_last_carried"]] = 1
                X[ei, b, name_i[f"{ch}_time_since_last_h"]] = 0.0
                mask_obs[ei, b, name_i[f"{ch}_time_since_last_h"]] = 1
            elif prev_t is not None:
                # no in-bin value; time since last, last-carried = 1
                tsl = float(BIN_MIDS[b] - prev_t)
                tsl = float(min(max(tsl, 0.0), TSL_CAP_H))
                X[ei, b, name_i[f"{ch}_last_carried"]] = 1.0
                mask_obs[ei, b, name_i[f"{ch}_last_carried"]] = 1
                X[ei, b, name_i[f"{ch}_time_since_last_h"]] = tsl
                mask_obs[ei, b, name_i[f"{ch}_time_since_last_h"]] = 1

    # pSOFA resp per bin
    psofa = np.full((n, N_BINS), np.nan, dtype=np.float32)
    pf = obs[obs["channel"] == "pf_ratio"]
    if len(pf):
        sc = M.psofa_score(pf["value"].to_numpy(dtype=float), M.PSOFA_PF)
        tmp = pd.DataFrame({"enc_i": pf["enc_i"].to_numpy(dtype=int), "t_step": pf["t_step"].to_numpy(dtype=int), "sc": sc})
        for (ei, ts), val in tmp.groupby(["enc_i", "t_step"])["sc"].max().items():
            psofa[int(ei), int(ts)] = val
    sf = obs[obs["channel"] == "sf_ratio"]
    if len(sf):
        sc = M.psofa_score(sf["value"].to_numpy(dtype=float), M.PSOFA_SF)
        tmp = pd.DataFrame({"enc_i": sf["enc_i"].to_numpy(dtype=int), "t_step": sf["t_step"].to_numpy(dtype=int), "sc": sc})
        for (ei, ts), val in tmp.groupby(["enc_i", "t_step"])["sc"].max().items():
            if not np.isfinite(psofa[int(ei), int(ts)]):
                psofa[int(ei), int(ts)] = val
    jps = name_i["psofa_resp"]
    X[:, :, jps] = psofa
    mask_obs[:, :, jps] = np.isfinite(psofa).astype(np.uint8)

    # drug-group active flags + vaso concurrent, always defined on valid bins (0 if none)
    meds = meds.copy()
    meds["enc_i"] = meds["encounter_id"].map(id_to_i)
    meds = meds[meds["enc_i"].notna()].copy()
    meds["enc_i"] = meds["enc_i"].astype(int)
    starts = meds["start_hours"].to_numpy(dtype=float)
    ends = meds["end_hours"].to_numpy(dtype=float)
    groups = meds["group"].astype(str).to_numpy()
    eis = meds["enc_i"].to_numpy(dtype=int)
    for b in range(N_BINS):
        lo, hi = BIN_EDGES[b], BIN_EDGES[b + 1]
        hit = (starts < hi) & (ends > lo)
        if not hit.any():
            # still mark group flags as observed (0) on valid bins
            continue
        for gname in M.MED_GROUPS:
            j = name_i[f"{gname}_active"]
            sel = hit & (groups == gname)
            if sel.any():
                X[np.unique(eis[sel]), b, j] = 1.0
        # vaso concurrent
        vsel = hit & (groups == "vasoactive")
        if vsel.any():
            sub_e, sub_s, sub_n = eis[vsel], np.clip(starts[vsel], lo, hi), np.clip(ends[vsel], lo, hi)
            order = np.argsort(sub_e)
            sub_e, sub_s, sub_n = sub_e[order], sub_s[order], sub_n[order]
            i0 = 0
            jv = name_i["vasoactive_max_concurrent"]
            while i0 < len(sub_e):
                i1 = i0 + 1
                while i1 < len(sub_e) and sub_e[i1] == sub_e[i0]:
                    i1 += 1
                evs = [(sub_s[k], 1) for k in range(i0, i1)] + [(sub_n[k], -1) for k in range(i0, i1)]
                evs.sort(key=lambda x: (x[0], x[1]))
                c = mx = 0
                for _, dlt in evs:
                    c += dlt
                    if c > mx:
                        mx = c
                X[int(sub_e[i0]), b, jv] = float(mx)
                i0 = i1
    # fill 0 for drug flags / vaso concurrent on valid bins that are still NaN
    for gname in M.MED_GROUPS:
        j = name_i[f"{gname}_active"]
        miss = ~np.isfinite(X[:, :, j])
        X[:, :, j] = np.where(miss, 0.0, X[:, :, j])
        mask_obs[:, :, j] = mask_valid  # defined whenever the stay covers the bin
    jv = name_i["vasoactive_max_concurrent"]
    miss = ~np.isfinite(X[:, :, jv])
    X[:, :, jv] = np.where(miss, 0.0, X[:, :, jv])
    mask_obs[:, :, jv] = mask_valid

    # features on invalid bins stay NaN (beyond OUTTIME)
    invalid = mask_valid == 0
    X[invalid] = np.nan
    mask_obs[invalid] = 0

    return {
        "X": X, "mask_obs": mask_obs, "mask_valid": mask_valid,
        "feature_names": np.array(feat_names, dtype="U"),
        "encounter_ids": ids.astype("U"),
        "bin_edges": BIN_EDGES, "bin_midpoints": BIN_MIDS,
    }


# =============================================================================
# Validation
# =============================================================================
def principal_angles_deg(A, B):
    Qa, _ = np.linalg.qr(A, mode="reduced")
    Qb, _ = np.linalg.qr(B, mode="reduced")
    s = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
    return np.degrees(np.arccos(np.clip(s, 0.0, 1.0)))


def unfold(X, mode):
    return np.moveaxis(X, mode, 0).reshape(X.shape[mode], -1)


def r2_score(y, yhat):
    y = np.asarray(y, dtype=float).ravel()
    yhat = np.asarray(yhat, dtype=float).ravel()
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float("nan") if ss_tot <= 0 else float(1.0 - ss_res / ss_tot)


def adjusted_rand_index(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    n = len(y_true)
    t_codes, t_inv = np.unique(y_true, return_inverse=True)
    p_codes, p_inv = np.unique(y_pred, return_inverse=True)
    ct = np.zeros((len(t_codes), len(p_codes)), dtype=np.int64)
    for i in range(n):
        ct[t_inv[i], p_inv[i]] += 1
    def comb2(x):
        return x * (x - 1) / 2.0
    sum_ct, sum_r, sum_c = comb2(ct).sum(), comb2(ct.sum(1)).sum(), comb2(ct.sum(0)).sum()
    comb_n = comb2(n)
    expected = (sum_r * sum_c) / comb_n if comb_n else 0.0
    max_ = 0.5 * (sum_r + sum_c)
    if max_ == expected:
        return 1.0 if sum_ct == expected else 0.0
    return float((sum_ct - expected) / (max_ - expected))


def kmeans(X, k, rng, n_init=20, n_iter=80):
    n, p = X.shape
    best = None
    for _ in range(n_init):
        centers = np.empty((k, p), dtype=float)
        centers[0] = X[rng.integers(0, n)]
        d2 = np.full(n, np.inf)
        for j in range(1, k):
            diff = X - centers[j - 1]
            d2 = np.minimum(d2, np.einsum("ij,ij->i", diff, diff))
            pvec = d2 / d2.sum() if d2.sum() > 0 else np.full(n, 1.0 / n)
            centers[j] = X[rng.choice(n, p=pvec)]
        labels = np.zeros(n, dtype=int)
        for it in range(n_iter):
            dist = ((X[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
            new = dist.argmin(axis=1)
            if it and np.array_equal(new, labels):
                break
            labels = new
            for j in range(k):
                m = labels == j
                centers[j] = X[m].mean(axis=0) if m.any() else X[rng.integers(0, n)]
        inertia = float(((X - centers[labels]) ** 2).sum())
        if best is None or inertia < best[2]:
            best = (labels.copy(), centers.copy(), inertia)
    return best


def validate_layer_c(lc: dict, plant: dict, enc: pd.DataFrame, rng: np.random.Generator,
                    loadings_path: Path | None = None) -> dict:
    """HOSVD on complete windows (LOS > W) so padding cannot leak duration.

    Validation uses the severity-related *mean* channels (plus vasoactive / pSOFA),
    not the full 75-feature tensor, so time-since-last / missingness patterns do not
    dominate mode-2. Trajectories are recovered by projecting onto the static
    severity OLS loading (mapped to Layer C names), then an intercept+scale OLS
    of demeaned d_e on that projection.
    """
    names = list(lc["feature_names"])
    name_i = {k: j for j, k in enumerate(names)}
    signal = [
        "heart_rate_mean", "sbp_mean", "dbp_mean", "map_mean",
        "lactate_mean", "creatinine_mean",
        "respiratory_rate_mean", "spo2_mean", "fio2_mean", "po2_mean",
        "pf_ratio_mean", "sf_ratio_mean",
        "vasoactive_active", "vasoactive_max_concurrent", "psofa_resp",
    ]
    cols = [c for c in signal if c in name_i]
    idx = [name_i[c] for c in cols]
    # map static OLS severity loadings onto these names
    beta = np.zeros(len(cols), dtype=float)
    default_sign = {
        "heart_rate_mean": 0.44, "sbp_mean": -0.40, "dbp_mean": -0.36, "map_mean": -0.38,
        "lactate_mean": 0.96, "creatinine_mean": 0.96,
        "respiratory_rate_mean": 0.0, "spo2_mean": 0.0, "fio2_mean": 0.0, "po2_mean": 0.0,
        "pf_ratio_mean": 0.0, "sf_ratio_mean": 0.0,
        "vasoactive_active": 0.52, "vasoactive_max_concurrent": 0.57, "psofa_resp": 0.01,
    }
    if loadings_path and Path(loadings_path).exists():
        ld = pd.read_csv(loadings_path).set_index("feature")
        rename = {c: c.replace("_mean", "_mean_0_24h") if c.endswith("_mean")
                  else ("vasoactive_exposed_0_24h" if c == "vasoactive_active"
                        else ("vasoactive_max_concurrent_0_24h" if c == "vasoactive_max_concurrent"
                              else ("psofa_resp_0_24h" if c == "psofa_resp" else c)))
                  for c in cols}
        for j, c in enumerate(cols):
            src = rename[c]
            if src in ld.index:
                beta[j] = float(ld.loc[src, "loading_severity_ols"])
            else:
                beta[j] = default_sign.get(c, 0.0)
    else:
        for j, c in enumerate(cols):
            beta[j] = default_sign.get(c, 0.0)
    if np.linalg.norm(beta) == 0:
        beta[0] = 1.0

    los = enc["los_hours"].to_numpy(dtype=float)
    cls_all = plant["class"]
    D_all = plant["D"]

    def run_window(hours: int, label: str) -> dict:
        Tsub = int(hours / BIN_H)
        keep = los > hours
        n_k = int(keep.sum())
        if n_k < 50:
            return {"n_encounters_in_window": n_k, "note": "too few encounters"}
        X = lc["X"][keep][:, :Tsub, :][:, :, idx].astype(np.float64)
        m = lc["mask_obs"][keep][:, :Tsub, :][:, :, idx].astype(bool)
        # standardize observed cells, impute 0
        for j in range(X.shape[-1]):
            v = X[:, :, j][m[:, :, j]]
            if len(v) >= 5:
                mu, s = v.mean(), v.std(ddof=1)
                s = s if s > 1e-8 else 1.0
                X[:, :, j] = np.where(m[:, :, j], (X[:, :, j] - mu) / s, 0.0)
            else:
                X[:, :, j] = 0.0
        Us, Ss = [], []
        for mode in range(3):
            U, S, _ = np.linalg.svd(unfold(X, mode), full_matrices=False)
            Us.append(U); Ss.append(S)
        U_time, sv_time = Us[1], Ss[1]
        mids = BIN_MIDS[:Tsub]
        cols_p, kinds_used = [], []
        default_on = {"spike_recovery": 48.0, "late_deterioration": 60.0}
        for kind in SHAPE_CLASSES:
            if kind == "stable":
                continue
            v = unit_norm_center(shape_raw(kind, mids, onset=default_on.get(kind)))
            if np.linalg.norm(v) > 0:
                cols_p.append(v); kinds_used.append(kind)
        planted = np.column_stack(cols_p)
        angles = {}
        for name, rec in {"mode2_0_4": U_time[:, :min(4, U_time.shape[1])],
                          "mode2_1_5": U_time[:, 1:min(5, U_time.shape[1])]}.items():
            k = min(rec.shape[1], planted.shape[1])
            ang = principal_angles_deg(planted, rec[:, :max(k, 1)])
            angles[name] = [float(a) for a in ang]
            angles[name + "_max"] = float(np.max(ang)) if len(ang) else None
        one = np.ones((Tsub, 1))
        def orth(M):
            return M - one @ (np.linalg.lstsq(one, M, rcond=None)[0])
        rec = U_time[:, 1:min(5, U_time.shape[1])]
        if rec.shape[1] >= 1:
            ang_c = principal_angles_deg(orth(planted), orth(rec))
            angles["mode2_1_5_orth_to_constant"] = [float(a) for a in ang_c]
            angles["mode2_1_5_orth_to_constant_max"] = float(np.max(ang_c))

        # project onto severity loading, recover d_e with intercept+scale OLS
        b = beta / np.linalg.norm(beta)
        proj = X @ b
        proj_c = proj - proj.mean(axis=1, keepdims=True)
        D = D_all[keep][:, :Tsub]
        D_c = D - D.mean(axis=1, keepdims=True)
        # OLS: vec(D_c) ~ a + b * vec(proj_c)
        y = D_c.ravel(); x = proj_c.ravel()
        A = np.column_stack([np.ones(len(x)), x])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        yhat = A @ coef
        r2 = r2_score(y, yhat)
        corr = float(np.corrcoef(y, x)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else float("nan")
        cls = cls_all[keep]
        r2_e = [r2_score(D_c[i], coef[0] + coef[1] * proj_c[i]) for i in range(len(D)) if cls[i] != "stable"]
        labels, _, _ = kmeans(proj_c, k=len(SHAPE_CLASSES), rng=rng, n_init=25)
        ari = adjusted_rand_index(cls, labels)
        # exclude late_deterioration at short windows (onset often after the window)
        mask_nl = cls != "late_deterioration"
        ari_nl = adjusted_rand_index(cls[mask_nl], labels[mask_nl]) if mask_nl.sum() > 50 else None
        return {
            "label": label,
            "n_encounters_in_window": n_k,
            "T_window": Tsub,
            "features_used": cols,
            "mode2_time_singular_values": [float(v) for v in sv_time],
            "mode2_variance_explained": [float(v) for v in (sv_time ** 2) / (sv_time ** 2).sum()],
            "shape_kinds_in_basis": kinds_used,
            "principal_angles_degrees": angles,
            "trajectory_R2_ols_pooled": r2,
            "trajectory_correlation": corr,
            "trajectory_R2_per_encounter_mean_nonstable": float(np.nanmean(r2_e)) if r2_e else None,
            "kmeans_ARI": ari,
            "kmeans_ARI_excluding_late_deterioration": ari_nl,
        }

    w72 = run_window(72, "72h complete window (LOS>72, first 12 bins)")
    w120 = run_window(120, "120h complete window (LOS>120, first 20 bins)")
    w168 = run_window(168, "168h complete window (LOS>168, first 28 bins)")
    return {
        "note": "HOSVD and trajectory recovery use severity-related Layer C mean channels on complete stay windows so mask_valid padding cannot leak LOS.",
        "window_72h": w72, "window_120h": w120, "window_168h": w168,
        # convenience aliases from the 72 h window for older report fields
        "n_encounters_in_window": w72.get("n_encounters_in_window"),
        "T_window": w72.get("T_window"),
        "mode2_time_singular_values": w72.get("mode2_time_singular_values", []),
        "principal_angles_degrees": w72.get("principal_angles_degrees", {}),
        "trajectory_R2_pooled": w72.get("trajectory_R2_ols_pooled"),
        "trajectory_R2_per_encounter_mean_nonstable": w72.get("trajectory_R2_per_encounter_mean_nonstable"),
        "kmeans_ARI": w72.get("kmeans_ARI"),
    }


def los_and_eligibility_report(enc: pd.DataFrame, lm: pd.DataFrame) -> dict:
    los = enc["los_hours"]
    q = {k: float(los.quantile(p)) for k, p in [("min", 0), ("p25", .25), ("p50", .5), ("p75", .75),
                                               ("p90", .9), ("p99", .99), ("max", 1)]}
    q["mean"] = float(los.mean())
    elig = {}
    for L in LANDMARKS:
        n = int((los > L).sum())
        elig[str(L)] = {"n": n, "pct": float(100.0 * n / len(enc))}
    # also from landmark file (death filter, should match)
    elig_file = lm.groupby("landmark_h").size().to_dict() if len(lm) else {}
    return {"n_encounters": int(len(enc)), "los_hours": q,
            "eligibility_outtime_gt_t0_plus_L": elig,
            "eligibility_landmark_index_rows": {str(k): int(v) for k, v in elig_file.items()},
            "note": "Existing ICUSTAYS/canonical OUTTIME were not modified. Generated stays have max LOS ≈ 280 h, so eligibility at 288, 312 and 336 h is zero, and Outcome 1 level 7 (remaining ICU time > 336 h) is empty at every landmark."}


def outcome_counts(lm: pd.DataFrame) -> dict:
    out = {}
    for L in LANDMARKS:
        sub = lm[lm["landmark_h"] == L]
        rec = {"n_eligible": int(len(sub)), "mortality": float(sub["in_hospital_mortality"].mean()) if len(sub) else None}
        for split, g in (("all", sub),) + tuple((s, sub[sub["split"] == s]) for s in ("train", "validation", "test")):
            if split != "all":
                rec.setdefault("by_split", {})[split] = {
                    "n": int(len(g)),
                    "outcome1": {str(k): int(v) for k, v in g["outcome1_remaining_icu"].value_counts().sort_index().items()} if len(g) else {},
                    "outcome2": {str(k): int(v) for k, v in g["outcome2_post_icu"].value_counts().sort_index().items()} if len(g) else {},
                    "mortality": float(g["in_hospital_mortality"].mean()) if len(g) else None,
                }
        rec["outcome1"] = {str(k): int(v) for k, v in sub["outcome1_remaining_icu"].value_counts().sort_index().items()} if len(sub) else {}
        rec["outcome2"] = {str(k): int(v) for k, v in sub["outcome2_post_icu"].value_counts().sort_index().items()} if len(sub) else {}
        if L in (24, 48, 72) and len(sub):
            rec["inherited_3_0"] = {str(k): int(v) for k, v in sub["outcome_inherited_3_0_ordinal"].value_counts().sort_index().items()}
        out[str(L)] = rec
    return out


def write_layer_c_long_gz(path: Path, lc: dict, max_rows_hint: int = 20_000_000) -> bool:
    """Observed cells only (mask_obs=1). Skip if it would be huge; npz is the primary format."""
    import gzip
    X, mo, ids, names = lc["X"], lc["mask_obs"], lc["encounter_ids"], lc["feature_names"]
    n_obs = int(mo.sum())
    if n_obs > max_rows_hint:
        return False
    with gzip.open(path, "wt", compresslevel=6) as g:
        g.write("encounter_id,t_step,t_start_h,t_end_h,t_mid_h,feature,value\n")
        n, t, f = X.shape
        for i in range(n):
            eid = ids[i]
            for b in range(t):
                if not mo[i, b].any():
                    continue
                lo, hi, mid = BIN_EDGES[b], BIN_EDGES[b + 1], BIN_MIDS[b]
                prefix = f"{eid},{b},{lo:g},{hi:g},{mid:g},"
                for j in range(f):
                    if mo[i, b, j]:
                        g.write(f"{prefix}{names[j]},{float(X[i, b, j]):.8g}\n")
    return True


def write_readme(root: Path, report: dict) -> None:
    text = f"""# Layered tensor dataset (A / B / C) for landmark-forward models

Companion to the three-stage synthetic SVD dataset. This folder follows the **three-layer process of the PIC Foundations pipeline** (contract v12 + v12.1 amendment A1.2–A1.7). It does **not** replace stages 1–3 or the three planted static factors; those files were not modified.

- Generator: `make_synthetic_svd_layered_tensor.py`
- Seed: {report['seed']}
- Subjects / encounters: {report['n_subjects']:,} / {report['n_encounters']:,} (same IDs as the existing dataset)
- Static factor scores: reused from `4_answer_key/latent_factors.csv` (unchanged)

## Why late landmarks are empty

Existing `ICUSTAYS` / canonical `OUTTIME` were **not** changed (same INTIME/OUTTIME/DOB as `1_raw`). The already-generated stays have

- LOS hours: min {report['los']['los_hours']['min']:.1f}, p25 {report['los']['los_hours']['p25']:.1f}, median {report['los']['los_hours']['p50']:.1f}, p75 {report['los']['los_hours']['p75']:.1f}, p90 {report['los']['los_hours']['p90']:.1f}, p99 {report['los']['los_hours']['p99']:.1f}, **max {report['los']['los_hours']['max']:.1f}**

so **eligibility at 288, 312 and 336 h is zero**, and Outcome 1 level 7 (remaining ICU time > 336 h) cannot occur. Core landmark 336 h therefore has no stacked rows. This is documented, not repaired by editing existing stays.

Eligibility (`last OUTTIME > t0 + L`):

{report['elig_md']}

## Folder layout

- `5_layered_tensor/1_extended_raw/` — companion raw events over each encounter's **whole existing stay** (plus a 6 h pre-ICU lookback). PATIENTS, ADMISSIONS, ICUSTAYS, canonical, maps and splits are **copies** of the originals (IDs and timestamps unchanged). New CHARTEVENTS / LABEVENTS / PRESCRIPTIONS are driven by `severity_t`. `outcomes_source.csv` holds death, same-HADM readmission and hospital discharge times used for A1.2 outcomes; existing ADMISSIONS.DISCHTIME and canonical `Died` are left as they were.
- `5_layered_tensor/2_normalized_extended/` — same v1 maps and D1–D20 rules as `2_normalized`.
- `5_layered_tensor/3_layers/` — Layer A CSV, Layer B CSV, Layer C npz, landmark index, outcomes-by-landmark, optional long gz.
- `4_answer_key/layered_tensor/` — **new** ground-truth files only (trajectories, shape basis, outcome coefficients, validation). Existing answer-key CSVs and zip were not modified.

## Time-varying ground truth (drives extended raw)

`severity_t(e, t) = s_e + d_e(t)`, with `s_e` the existing static `z_severity`.

`d_e(t) = amplitude_e * shape_{{k(e)}}(t)` on 6 h bin midpoints, k ∈ {{stable, improving, worsening, spike_recovery, late_deterioration}} (i.i.d. uniform), amplitude ~ Uniform({AMP_LOW}, {AMP_HIGH}) for non-stable classes (0 if stable).

| Class | Raw shape (then mean-centered and L2-normalized on the 56 midpoints) |
|---|---|
| stable | 0 |
| improving | exp(-2 u), u = t/336 |
| worsening | logistic 1/(1+exp(-8(u-0.5))) |
| spike_recovery | Gaussian bump, peak ~ Uniform(12, 96) h, width 18 h |
| late_deterioration | 0 then logistic rise after onset ~ Uniform(48, 200) h |

Respiratory: `resp_t(e,t) = z_respiratory + {RESP_AMP_SCALE} * d_e(t)` (same shape, milder). Age remains a static patient-level factor (DOB unchanged). Event values use the same raw-scale loadings as the static generator (`BETA` in `make_synthetic_svd_dataset.py`), with means shifted by `severity_t` / `resp_t` at the event time (pre-t0 uses the static scores only). Vasoactive starts are also tied to `severity_t`.

## Layers (A1.6)

**Layer A** — one row per encounter, fixed at t0. Age = INTIME − DOB with `age_invalid` (D20), sex, first care unit, admitting diagnosis domain from `ADMISSIONS.ICD10_CODE_CN` (3-character → `icd_domain`), admission source derived from `EDREGTIME` (ed / unknown), recent surgery (anesthesia-end in [t0 − 24 h, t0 + 6 h], A2.1), remote same-HADM surgery, weight (charted or imputed, A2.12/D20), and the Layer A pre-ICU laboratory block: last released value in **[t0 − 6 h, t0)** plus a mask (A2.11). Events before t0 never enter Layer B or C.

**Layer B** — one row per encounter, computed once on **[t0, t0 + 24 h)** from the extended normalized stay, held fixed at every landmark. Summaries match stage 3: last/min/max/mean/count, pSOFA respiratory (P/F then S/F, D19/D15), MAP, lactate, creatinine, vasoactive flag and concurrent count, drug-group flags, fluid intake mL/kg/day, missingness. Because the extended raw is a new companion stay (full horizon, trajectory-driven), Layer B is **not** a copy of `3_features/`; that file is unchanged.

**Layer C** — tensor `X_C` with shape **(N, 56, F_C)** = ({report['n_encounters']}, {N_BINS}, {report['F_C']}). 6 h bins on [0, 336) h; bin *b* is `[6b, 6b+6)`. Observation time is `t_hours` (labs already include the +1 h release lag). Per-bin features: channel means, time-since-last (capped at 24 h, A2.13) and last-carried indicators, drug-group active flags, vasoactive concurrent count, pSOFA respiratory. **Variable-length handling:** `mask_valid[e, b] = 1` iff the bin starts before that encounter's last OUTTIME (`6b < LOS`); `mask_obs[e, b, f] = 1` iff a value was observed (or a flag is defined) in a valid bin. Invalid bins are NaN in `X_C` (no duration leak through padding).

**Landmark index** (`landmark_index.csv`): one row per (encounter, L) with `OUTTIME > t0 + L`, L ∈ {{24, 48, …, 336}}, and `n_bins = L/6`. Slice the sequence as `X_C[:, :L//6, :]`.

**72 h fixed-window sensitivity (A1.6):** at landmark L use the last 12 bins before L, i.e. `X_C[:, max(0, L//6 - 12):L//6, :]`. Encounters with L ≤ 72 h contribute their full prefix.

## Outcomes (A1.2–A1.4)

Defined in `outcomes_by_landmark.csv` from `outcomes_source.csv`, **not** from existing ADMISSIONS.DISCHTIME or canonical `Died` (those columns were left unchanged; agreement is reported in the validation JSON).

- **Eligibility:** last OUTTIME > t0 + L (and death, when present, after t0 + L).
- **Outcome 1** (primary, 9 levels): 8 = died in hospital any time after L; 7 remaining ICU > 336 h; 6 (240, 336]; 5 (168, 240]; 4 (120, 168]; 3 (72, 120]; 2 (48, 72]; 1 (24, 48]; 0 ≤ 24 h. Remaining time is existing OUTTIME − (t0 + L).
- **Outcome 2** (key secondary, 5 levels): 4 death after L; 3 same-hospitalization ICU readmission > 24 h after OUTTIME; 2 no readmission, post-ICU stay ≥ 14 d; 1 [7, 14) d; 0 < 7 d. Post-ICU stay = hospital_dischtime − OUTTIME. Unassignable (no discharge time and no readmission/death) is flagged; it should be zero.
- **In-hospital mortality** (binary) at every landmark.
- **Inherited 3.0 ordinal-8** at 24 / 48 / 72 h only: **total** ICU stay from t0, same cut points, death = 8. This is the definition that collapses at late landmarks (A1.1). Documented as codes 0–7 (stay bands) + 8 (death).

Death takes precedence. Outcome 1 for readmitted patients is still measured to the **index** encounter's last OUTTIME.

**Generative model** (coefficients in the answer key): see `outcomes_source` docstring in the generator. Mortality is calibrated to about 7.5%. Outcome 1 levels 0–7 for survivors are fixed by the already-generated LOS (independent of this outcome model, because OUTTIME was not changed). Death and Outcome 2 are factor- and trajectory-dependent.

Cohort mortality in `outcomes_source`: {report.get('death_rate', float('nan')):.2%}. Readmission rate: {report.get('readmit_rate', float('nan')):.2%}.

{report['outcome_md']}

## Validation (HOSVD / class / outcomes)

HOSVD was run on the **complete 72 h window** (encounters with LOS > 72 h, first 12 bins), standardized on observed cells, missing features zero-filled *inside that complete window only*, so stay-length padding cannot leak.

{report['val_md']}

## Formats

- `layer_c.npz`: `X_C` float32 (N, 56, F_C), `mask_obs` uint8, `mask_valid` uint8 (N, 56), `encounter_ids`, `bin_edges`, `bin_midpoints`, `feature_names`.
- Layer A / B CSVs; `landmark_index.csv`; `outcomes_by_landmark.csv`; `outcomes_source.csv` (in 1_extended_raw).
- Long Layer C CSV.gz if the observed-cell count permits.
- Ground truth: `4_answer_key/layered_tensor/`.

Splits are the existing `subject_splits.csv`.
"""
    (root / "README_layered_tensor.md").write_text(text)
    (root / "5_layered_tensor" / "README_layered_tensor.md").write_text(text)


def append_main_readme(root: Path, report: dict) -> None:
    path = root / "README_synthetic_svd_levels.md"
    old = path.read_text()
    marker = "\n## Layered tensor extension (A / B / C, amendment A1.2–A1.7)\n"
    if marker.strip() in old:
        # replace previous append if re-run
        old = old.split(marker)[0].rstrip() + "\n"
    add = marker + f"""
A **new** companion dataset, `5_layered_tensor/`, follows the three-layer PIC pipeline (contract v12 + v12.1 A1.2–A1.7). Stages 1–3, the three static factor scores, and the original answer-key files in `4_answer_key/` were **not** modified.

See **[README_layered_tensor.md](README_layered_tensor.md)** for layers, binning, masking, outcome definitions, the LOS/eligibility limitation (existing stays end by ~{report['los']['los_hours']['max']:.0f} h, so landmarks 288–336 h have zero eligibility), ground truth, and validation numbers.

- Layer C tensor: N={report['n_encounters']}, T=56 (6 h bins, 0–336 h), F={report['F_C']}; `mask_valid` is 0 after each encounter's existing OUTTIME; `mask_obs` is per bin × feature.
- Landmark index slices `X_C[:, :L//6, :]`. 72 h sensitivity: `X_C[:, max(0, L//6-12):L//6, :]`.
- Outcomes are in `5_layered_tensor/3_layers/outcomes_by_landmark.csv`, generated from `outcomes_source.csv` (existing ADMISSIONS/ICUSTAYS timestamps were not overwritten).
- Trajectory and outcome-model ground truth: `4_answer_key/layered_tensor/` (new files only).
"""
    path.write_text(old.rstrip() + "\n" + add)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=HERE)
    ap.add_argument("--seed", type=int, default=LAYERED_SEED)
    ap.add_argument("--start-from", choices=("raw", "normalize", "layers", "validate"), default="raw",
                    help="Resume after a previous run: skip earlier stages if their outputs already exist.")
    args = ap.parse_args()
    root = args.root.resolve()
    ext_root = root / "5_layered_tensor"
    ext_raw = ext_root / "1_extended_raw"
    ext_norm = ext_root / "2_normalized_extended"
    ext_lay = ext_root / "3_layers"
    key = root / "4_answer_key" / "layered_tensor"
    for p in (ext_raw, ext_norm, ext_lay, key):
        p.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    factors = pd.read_csv(root / "4_answer_key" / "latent_factors.csv", dtype={"encounter_id": str})
    n = len(factors)
    print(f"[1/6] plant trajectories n={n} seed={args.seed}", flush=True)
    plant = plant_trajectories(n, rng)
    if args.start_from == "validate":
        print("[validate] reloading layers / outcomes and refreshing validation + README", flush=True)
        enc = pd.read_csv(ext_norm / "encounters.csv", dtype={"encounter_id": str})
        enc = enc.set_index("encounter_id").loc[factors["encounter_id"]].reset_index()
        if "los_hours" not in enc.columns:
            enc["los_hours"] = (pd.to_datetime(enc["OUTTIME"]) - pd.to_datetime(enc["INTIME"])).dt.total_seconds() / 3600.0
        for col in ("z_severity", "z_respiratory", "z_age"):
            enc[col] = factors[col].to_numpy()
        bundled = np.load(ext_lay / "layer_c.npz", allow_pickle=False)
        lc = {k: bundled[k] for k in bundled.files}
        # numpy string arrays
        if "X_C" in lc:
            lc["X"] = lc["X_C"]
        oc = pd.read_csv(ext_raw / "outcomes_source.csv", dtype={"encounter_id": str})
        lm = pd.read_csv(ext_lay / "outcomes_by_landmark.csv", dtype={"encounter_id": str})
        val = validate_layer_c(lc, plant, enc, np.random.default_rng(args.seed + 4),
                               loadings_path=root / "4_answer_key" / "planted_loadings.csv")
        los_rep = los_and_eligibility_report(enc, lm)
        counts = outcome_counts(lm)
        print("      72h R2", val["window_72h"].get("trajectory_R2_ols_pooled"),
              "corr", val["window_72h"].get("trajectory_correlation"),
              "ARI", val["window_72h"].get("kmeans_ARI"), flush=True)
        print("      120h R2", val["window_120h"].get("trajectory_R2_ols_pooled"),
              "corr", val["window_120h"].get("trajectory_correlation"),
              "ARI", val["window_120h"].get("kmeans_ARI"), flush=True)
        print("      168h R2", val["window_168h"].get("trajectory_R2_ols_pooled"),
              "corr", val["window_168h"].get("trajectory_correlation"),
              "ARI", val["window_168h"].get("kmeans_ARI"), flush=True)
        print("      angles 72h orth", val["window_72h"].get("principal_angles_degrees", {}).get("mode2_1_5_orth_to_constant"), flush=True)
        report_json = json.loads((key / "layered_tensor_validation.json").read_text()) if (key / "layered_tensor_validation.json").exists() else {}
        report_json.update({"validation": val, "los": los_rep, "outcomes": counts,
                            "death_rate": float(oc["died_in_hospital"].mean()),
                            "readmit_rate": float(oc["readmission"].mean()),
                            "n_encounters": n, "F_C": int(lc["X"].shape[-1])})
        (key / "layered_tensor_validation.json").write_text(json.dumps(report_json, indent=2) + "\n")
        (ext_lay / "layered_tensor_validation.json").write_text(json.dumps(report_json, indent=2) + "\n")
        elig_lines = ["| L (h) | n eligible | % of 17,239 |", "|---:|---:|---:|"]
        for L in LANDMARKS:
            e = los_rep["eligibility_outtime_gt_t0_plus_L"][str(L)]
            elig_lines.append(f"| {L} | {e['n']:,} | {e['pct']:.2f}% |")
        out_lines = ["### Outcome counts at core landmarks", "",
                     "| L | n | mortality | Outcome 1 (8/7/6/5/4/3/2/1/0) | Outcome 2 (4/3/2/1/0) |",
                     "|---:|---:|---:|---|---|"]
        for L in CORE_LANDMARKS:
            c = counts[str(L)]
            def fmt(d, keys):
                return "/".join(str(d.get(str(k), 0)) for k in keys)
            mort = "n/a" if c["mortality"] is None or c["n_eligible"] == 0 else f"{100*c['mortality']:.1f}%"
            out_lines.append(
                f"| {L} | {c['n_eligible']:,} | {mort} | {fmt(c['outcome1'], range(8,-1,-1))} | {fmt(c['outcome2'], range(4,-1,-1))} |"
            )
        def fmt_win(w):
            if not w or w.get("n_encounters_in_window", 0) < 50:
                return f"- {w.get('label','window')}: too few encounters\n"
            ang = w.get("principal_angles_degrees", {})
            return (
                f"- {w['label']}: n={w['n_encounters_in_window']:,}. "
                f"Mode-2 SVs: {', '.join(f'{x:.2f}' for x in w['mode2_time_singular_values'][:8])}. "
                f"Principal angles (deg) vs planted shapes, U2[:,1:5] ⊥ constant: "
                f"{', '.join(f'{a:.2f}' for a in ang.get('mode2_1_5_orth_to_constant', []))} "
                f"(max {ang.get('mode2_1_5_orth_to_constant_max', float('nan')):.2f}). "
                f"Trajectory OLS R²={w.get('trajectory_R2_ols_pooled'):.3f}, corr={w.get('trajectory_correlation'):.3f}. "
                f"k-means ARI={w.get('kmeans_ARI'):.3f} (excl. late_deterioration: {w.get('kmeans_ARI_excluding_late_deterioration')}).\n"
            )
        val_md = fmt_win(val["window_72h"]) + fmt_win(val["window_120h"]) + fmt_win(val["window_168h"])
        report = {
            "seed": args.seed, "n_subjects": int(enc["SUBJECT_ID"].nunique()), "n_encounters": n,
            "F_C": int(lc["X"].shape[-1]), "los": los_rep,
            "death_rate": float(oc["died_in_hospital"].mean()),
            "readmit_rate": float(oc["readmission"].mean()),
            "elig_md": "\n".join(elig_lines), "outcome_md": "\n".join(out_lines), "val_md": val_md,
        }
        write_readme(root, report)
        append_main_readme(root, report)
        print("validate refresh done", flush=True)
        return

    if args.start_from == "raw":
        print("[2/6] extended raw (whole existing stay, severity_t-driven)", flush=True)
        raw_meta = generate_extended_raw(root, ext_raw, plant, rng)
    else:
        print("[2/6] skipping extended raw (--start-from %s)" % args.start_from, flush=True)
        enc0 = pd.read_csv(root / "2_normalized" / "encounters.csv", dtype={"encounter_id": str})
        enc0 = enc0.merge(factors, on="encounter_id", how="left", suffixes=("", "_fac"))
        enc0 = enc0.set_index("encounter_id").loc[factors["encounter_id"]].reset_index()
        raw_meta = {"enc": enc0, "n": n, "ids": factors["encounter_id"].to_numpy(),
                    "chart_n": None, "lab_n": None, "rx_n": None}

    if args.start_from in ("raw", "normalize"):
        if args.start_from == "normalize":
            print("[3/6] normalize extended raw (v1 maps, D1–D20)", flush=True)
        else:
            print("[3/6] normalize extended raw (v1 maps, D1–D20)", flush=True)
        norm = normalize(ext_raw, ext_norm, raw_meta["enc"])
        print(f"      obs={norm['n_obs']} meds={norm['n_meds']} invalid_dob={norm['n_invalid_dob']}", flush=True)
    else:
        print("[3/6] skipping normalize; loading 2_normalized_extended", flush=True)
        # Reconstruct the norm dict needed by engineer_features / later steps
        obs = pd.read_csv(ext_norm / "observations_long.csv", dtype={"encounter_id": str})
        meds = pd.read_csv(ext_norm / "medications_long.csv", dtype={"encounter_id": str})
        enc_n = pd.read_csv(ext_norm / "encounters.csv", dtype={"encounter_id": str})
        # intern form expected by engineer_features: enc_i, t0, out, split, ...
        enc_n = enc_n.set_index("encounter_id").loc[factors["encounter_id"]].reset_index()
        enc_n["t0"] = M.fast_time(enc_n["INTIME"])
        enc_n["out"] = M.fast_time(enc_n["OUTTIME"])
        enc_n["enc_i"] = np.arange(len(enc_n), dtype=np.int32)
        if "FIRST_CAREUNIT" not in enc_n.columns and "first_careunit" in enc_n.columns:
            enc_n["FIRST_CAREUNIT"] = enc_n["first_careunit"]
        if "Died" not in enc_n.columns and "died" in enc_n.columns:
            enc_n["Died"] = enc_n["died"]
        # observations need enc_i
        ei_map = dict(zip(enc_n["encounter_id"], enc_n["enc_i"]))
        if "enc_i" not in obs.columns:
            obs["enc_i"] = obs["encounter_id"].map(ei_map)
        cm = M.parse_channel_map(pd.read_csv(ext_raw / "feature_channel_map_v1.csv", dtype=str).fillna(""))
        norm = {"enc": enc_n, "obs": obs, "meds": meds, "cm": cm,
                "n_obs": len(obs), "n_meds": len(meds),
                "n_invalid_dob": int((enc_n["dob_status"] != "valid").sum()) if "dob_status" in enc_n.columns else 0,
                "dob_status_counts": enc_n["dob_status"].value_counts().to_dict() if "dob_status" in enc_n.columns else {}}
        raw_meta["chart_n"] = raw_meta["lab_n"] = raw_meta["rx_n"] = "see existing CSVs"

    print("[4/6] Layer B (first 24 h) + Layer A + Layer C + outcomes", flush=True)
    feat_b = engineer_features(ext_norm, ext_lay / "_layer_b_scratch", norm)
    # keep Layer B as a clean CSV at 3_layers
    feat_b.to_csv(ext_lay / "layer_b.csv", index=False)
    shutil.rmtree(ext_lay / "_layer_b_scratch", ignore_errors=True)
    document_features(feat_b, ext_lay / "layer_b_data_dictionary.csv")

    enc = norm["enc"].copy()
    if "FIRST_CAREUNIT" in enc.columns and "first_careunit" not in enc.columns:
        enc["first_careunit"] = enc["FIRST_CAREUNIT"]
    if "los_hours" not in enc.columns:
        enc["los_hours"] = (pd.to_datetime(enc["OUTTIME"] if "OUTTIME" in enc.columns else enc["out"])
                            - pd.to_datetime(enc["INTIME"] if "INTIME" in enc.columns else enc["t0"])).dt.total_seconds() / 3600.0
    # align enc to factor order and attach static factor scores (existing, unchanged)
    enc = enc.set_index("encounter_id").loc[factors["encounter_id"]].reset_index()
    for col in ("z_severity", "z_respiratory", "z_age", "SUBJECT_ID", "HADM_ID"):
        if col not in enc.columns and col in factors.columns:
            enc[col] = factors[col].to_numpy()
        elif col in ("z_severity", "z_respiratory", "z_age"):
            enc[col] = factors[col].to_numpy()
    layer_a = build_layer_a(root, ext_norm, enc, feat_b)
    # preserve factor order
    layer_a = layer_a.set_index("encounter_id").loc[factors["encounter_id"]].reset_index()
    layer_a.to_csv(ext_lay / "layer_a.csv", index=False)

    lc = build_layer_c(ext_norm, enc)
    np.savez_compressed(
        ext_lay / "layer_c.npz",
        X_C=lc["X"], mask_obs=lc["mask_obs"], mask_valid=lc["mask_valid"],
        encounter_ids=lc["encounter_ids"], bin_edges=lc["bin_edges"],
        bin_midpoints=lc["bin_midpoints"], feature_names=lc["feature_names"],
    )
    fill = float((lc["mask_obs"].sum() / max(1, lc["mask_valid"].sum() * lc["X"].shape[-1])))
    fill_valid_bins = float(lc["mask_valid"].mean())
    print(f"      Layer C shape={lc['X'].shape} valid-bin fill={fill_valid_bins:.3f} obs-fill-in-valid≈{fill:.3f}", flush=True)

    oc = generate_outcomes(enc, plant, np.random.default_rng(args.seed + 3))
    oc.to_csv(ext_raw / "outcomes_source.csv", index=False)
    lm = build_landmark_outcomes(enc, oc)
    lm.to_csv(ext_lay / "outcomes_by_landmark.csv", index=False)
    idx = lm[["encounter_id", "SUBJECT_ID", "split", "landmark_h", "eligible", "n_bins"]].copy()
    idx.to_csv(ext_lay / "landmark_index.csv", index=False)
    print(f"      mortality={oc.attrs['death_rate']:.3f} readmit={oc.attrs['readmit_rate']:.3f} landmark_rows={len(lm)}", flush=True)

    print("[5/6] validate + answer key", flush=True)
    val = validate_layer_c(lc, plant, enc, np.random.default_rng(args.seed + 4),
                           loadings_path=root / "4_answer_key" / "planted_loadings.csv")
    los_rep = los_and_eligibility_report(enc, lm)
    counts = outcome_counts(lm)
    print("      mode-2 SVs:", [round(x, 2) for x in val["mode2_time_singular_values"][:8]], flush=True)
    print("      angles 1-5 orth:", val["principal_angles_degrees"].get("mode2_1_5_orth_to_constant"), flush=True)
    print("      R2:", val["trajectory_R2_pooled"], "ARI:", val["kmeans_ARI"], flush=True)

    # answer key (NEW files only)
    traj = pd.DataFrame({
        "encounter_id": factors["encounter_id"],
        "SUBJECT_ID": factors["SUBJECT_ID"],
        "HADM_ID": factors["HADM_ID"],
        "z_severity": factors["z_severity"],
        "z_respiratory": factors["z_respiratory"],
        "z_age": factors["z_age"],
        "shape_class": plant["class"],
        "amplitude": plant["amplitude"],
        "onset_h": plant["onset_h"],
    })
    for b in range(N_BINS):
        traj[f"d_t{b}"] = plant["D"][:, b]
    traj.to_csv(key / "trajectory_scores.csv", index=False)
    sh = pd.DataFrame({"t_step": np.arange(N_BINS), "t_start_h": BIN_EDGES[:-1],
                       "t_end_h": BIN_EDGES[1:], "t_mid_h": BIN_MIDS})
    for kind in SHAPE_CLASSES:
        # canonical shape at default onset (per-encounter onsets live in trajectory_scores)
        on = {"spike_recovery": 48.0, "late_deterioration": 120.0}.get(kind)
        sh[f"shape_{kind}"] = unit_norm_center(shape_raw(kind, BIN_MIDS, onset=on)) if kind != "stable" else 0.0
    sh.to_csv(key / "shape_basis.csv", index=False)
    coefs = {
        "death_logit": {
            "a0_calibrated": oc.attrs["a0"],
            "z_severity": 0.85, "end_severity": 0.55, "z_respiratory": 0.25,
            "late_deterioration": 0.45, "mean_d": 0.15,
            "target_mortality": DEATH_TARGET, "realized_mortality": oc.attrs["death_rate"],
        },
        "readmit_logit": {
            "b0_calibrated": oc.attrs["b0"],
            "z_severity": 0.55, "mean_d": 0.40, "z_respiratory": 0.20,
            "target_rate_among_non_icu_deaths": 0.08, "realized_rate": oc.attrs["readmit_rate"],
        },
        "post_icu_days": {
            "formula": "clip(exp(1.15 + 0.28 s_e + 0.20 end_severity + 0.15 z_resp + N(0,0.55)), 0.4, 45)",
        },
        "respiratory_trajectory": f"resp_t = z_respiratory + {RESP_AMP_SCALE} * d_e(t) (same class/shape, milder)",
        "static_factors_source": "4_answer_key/latent_factors.csv (unchanged)",
        "note_outcome1_survivor_levels": "Levels 0-7 are determined by existing OUTTIME (not redrawn). Death and Outcome 2 are generated here.",
    }
    (key / "outcome_model_coefficients.json").write_text(json.dumps(coefs, indent=2) + "\n")

    report_json = {
        "seed": args.seed, "n_subjects": int(enc["SUBJECT_ID"].nunique()),
        "n_encounters": n, "F_C": int(lc["X"].shape[-1]), "T": N_BINS,
        "layer_c_shape": list(lc["X"].shape),
        "layer_c_dtype": "float32",
        "valid_bin_rate": fill_valid_bins,
        "observed_fill_within_valid": fill,
        "extended_raw_rows": {"CHARTEVENTS": raw_meta["chart_n"], "LABEVENTS": raw_meta["lab_n"],
                              "PRESCRIPTIONS": raw_meta["rx_n"]},
        "normalized": {"n_obs": norm["n_obs"], "n_meds": norm["n_meds"], "n_invalid_dob": norm["n_invalid_dob"]},
        "los": los_rep, "outcomes": counts, "validation": val,
        "death_rate": oc.attrs["death_rate"], "readmit_rate": oc.attrs["readmit_rate"],
        "agreement_canonical_died": float((oc["died_in_hospital"].to_numpy() == enc["Died"].astype(int).to_numpy()).mean()) if "Died" in enc.columns else None,
        "layer_a_rows": int(len(layer_a)), "layer_b_rows": int(len(feat_b)),
        "landmark_index_rows": int(len(lm)),
        "canonical_died_note": "canonical/encounters Died is from the original generator and was not changed; layered outcomes use outcomes_source only.",
    }
    (key / "layered_tensor_validation.json").write_text(json.dumps(report_json, indent=2) + "\n")
    (ext_lay / "layered_tensor_validation.json").write_text(json.dumps(report_json, indent=2) + "\n")

    print("[6/6] long CSV + README", flush=True)
    wrote_long = write_layer_c_long_gz(ext_lay / "layer_c_long.csv.gz", lc)
    print("      long csv written:" , wrote_long, flush=True)

    # markdown snippets for the README
    elig_lines = ["| L (h) | n eligible | % of 17,239 |", "|---:|---:|---:|"]
    for L in LANDMARKS:
        e = los_rep["eligibility_outtime_gt_t0_plus_L"][str(L)]
        elig_lines.append(f"| {L} | {e['n']:,} | {e['pct']:.2f}% |")
    out_lines = ["### Outcome counts at core landmarks", "",
                 "| L | n | mortality | Outcome 1 (8/7/6/5/4/3/2/1/0) | Outcome 2 (4/3/2/1/0) |",
                 "|---:|---:|---:|---|---|"]
    for L in CORE_LANDMARKS:
        c = counts[str(L)]
        def fmt(d, keys):
            return "/".join(str(d.get(str(k), 0)) for k in keys)
        mort = "n/a" if c["mortality"] is None or c["n_eligible"] == 0 else f"{100*c['mortality']:.1f}%"
        out_lines.append(
            f"| {L} | {c['n_eligible']:,} | {mort} | {fmt(c['outcome1'], range(8,-1,-1))} | {fmt(c['outcome2'], range(4,-1,-1))} |"
        )
    ang = val["principal_angles_degrees"]
    val_md = (
        f"- n in 72 h window: {val['n_encounters_in_window']:,}\n"
        f"- Mode-2 singular values: {', '.join(f'{x:.2f}' for x in val['mode2_time_singular_values'])}\n"
        f"- Principal angles (deg) vs planted shape span, U2[:,1:5] after removing a constant: "
        f"{', '.join(f'{a:.2f}' for a in ang.get('mode2_1_5_orth_to_constant', []))} "
        f"(max {ang.get('mode2_1_5_orth_to_constant_max', float('nan')):.2f})\n"
        f"- Trajectory R² (pooled, demeaned d_e vs projected Layer C): {val['trajectory_R2_pooled']:.3f}; "
        f"mean per-encounter R² (non-stable): {val['trajectory_R2_per_encounter_mean_nonstable']}\n"
        f"- k-means ARI vs true shape class: {val['kmeans_ARI']:.3f}\n"
    )
    report = {
        "seed": args.seed, "n_subjects": report_json["n_subjects"], "n_encounters": n,
        "F_C": report_json["F_C"], "los": los_rep, "death_rate": oc.attrs["death_rate"],
        "readmit_rate": oc.attrs["readmit_rate"],
        "elig_md": "\n".join(elig_lines), "outcome_md": "\n".join(out_lines), "val_md": val_md,
    }
    write_readme(root, report)
    append_main_readme(root, report)
    print("done", flush=True)
    print(json.dumps({
        "N": n, "T": N_BINS, "F_C": report_json["F_C"],
        "valid_bin_rate": fill_valid_bins, "obs_fill": fill,
        "death_rate": oc.attrs["death_rate"], "readmit_rate": oc.attrs["readmit_rate"],
        "landmark_rows": int(len(lm)),
        "R2": val["trajectory_R2_pooled"], "ARI": val["kmeans_ARI"],
        "angles_orth": ang.get("mode2_1_5_orth_to_constant"),
        "elig_336": los_rep["eligibility_outtime_gt_t0_plus_L"]["336"],
        "elig_240": los_rep["eligibility_outtime_gt_t0_plus_L"]["240"],
        "elig_24": los_rep["eligibility_outtime_gt_t0_plus_L"]["24"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
