"""make_synthetic.py -- SYNTHETIC stand-ins for the PIC Foundations X_eng / X_raw matrices.

NOT PIC DATA. Every value is drawn from hand-specified distributions written below.
No real patient-level data or fitted real statistics were used. Only the column
schema (names/meanings of the 59 frozen columns) was taken from the methods draft.

Usage:  python make_synthetic.py [outdir]
"""
import sys, os
import numpy as np
import pandas as pd

SEED = 6310
N_ENC = 2000  # default size for a public stand-in; not a PIC cohort size
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(SEED)
n = N_ENC
expit = lambda x: 1 / (1 + np.exp(-x))

# ---------------- subjects, encounters, split ----------------
# ~85% of subjects have 1 encounter; the rest have 2-4 (repeat ICU stays)
subj = []
sid = 0
while len(subj) < n:
    sid += 1
    k = rng.choice([1, 2, 3, 4], p=[0.85, 0.11, 0.03, 0.01])
    subj += [sid] * k
subj = np.array(subj[:n])
subject_ids = np.unique(subj)
u = rng.random(len(subject_ids))
split_of_subject = np.where(u < 0.70, "train", np.where(u < 0.85, "val", "test"))
split = pd.Series(split_of_subject, index=subject_ids)[subj].to_numpy()

# ---------------- demographics / timing ----------------
sex = rng.choice(["M", "F"], size=n, p=[0.56, 0.44]).astype(object)
sex_missing = (rng.random(n) < 0.01).astype(int)
sex[sex_missing == 1] = np.nan
# age: mixture heavy in infants, clipped to [0, 17.99]
comp = rng.choice(3, size=n, p=[0.45, 0.30, 0.25])
age = np.where(comp == 0, rng.exponential(0.25, n),
      np.where(comp == 1, rng.uniform(1, 5, n), rng.uniform(5, 18, n)))
age = np.clip(age, 0, 17.99).round(2)
age_band = pd.cut(age, [-0.01, 0.0833, 1, 5, 12, 18],
                  labels=["neonate_lt1m", "infant_1m_1y", "child_1_4", "child_5_11", "adolescent_12_17"]).astype(str)
hours_admit = np.round(rng.lognormal(0.5, 1.2, n), 2)
hours_admit_missing = (rng.random(n) < 0.04).astype(int)
hours_admit[hours_admit_missing == 1] = np.nan
hour = rng.choice(24, size=n, p=(lambda w: w / w.sum())(1 + 0.8 * np.sin((np.arange(24) - 8) / 24 * 2 * np.pi) + 0.3))
is_night = ((hour >= 19) | (hour < 7)).astype(int)
dow = rng.choice(7, size=n, p=[0.16, 0.16, 0.16, 0.16, 0.16, 0.10, 0.10])  # 0=Mon

# clinical domain drives unit, chapter, surgery
domains = ["cardiac", "respiratory", "neurologic", "infectious", "surgical_other", "neonatal", "other"]
dom_p = np.array([0.22, 0.20, 0.10, 0.12, 0.14, 0.10, 0.12])
domain = rng.choice(domains, size=n, p=dom_p)
domain = np.where((age < 0.0833) & (rng.random(n) < 0.6), "neonatal", domain)
chapter_map = {"cardiac": "IX_circulatory", "respiratory": "X_respiratory", "neurologic": "VI_nervous",
               "infectious": "I_infectious", "surgical_other": "XIX_injury", "neonatal": "XVI_perinatal",
               "other": "XVIII_symptoms"}
chapter = np.array([chapter_map[d] for d in domain], dtype=object)
swap = rng.random(n) < 0.15  # chapter not a pure function of domain
chapter[swap] = rng.choice(list(chapter_map.values()) + ["XVII_congenital", "IV_endocrine"], size=swap.sum())
chapter_missing = (rng.random(n) < 0.05).astype(int)
chapter[chapter_missing == 1] = np.nan
unit = np.where(domain == "cardiac", rng.choice(["CICU", "PICU"], n, p=[0.8, 0.2]),
       np.where(domain == "neonatal", rng.choice(["NICU", "PICU"], n, p=[0.85, 0.15]),
       np.where(domain == "surgical_other", rng.choice(["SICU", "PICU"], n, p=[0.7, 0.3]),
                rng.choice(["PICU", "General_ICU"], n, p=[0.75, 0.25]))))
p_surg = np.select([domain == "cardiac", domain == "surgical_other"], [0.55, 0.65], 0.08)
surgery = (rng.random(n) < p_surg).astype(int)
postop = np.where(rng.random(n) < 0.92, surgery, 1 - surgery)  # ~92% agreement

# ---------------- latent factors ----------------
sev = rng.normal(0, 1, n) + 0.3 * (domain == "infectious") + 0.2 * surgery          # global acuity
resp = rng.normal(0, 1, n) + 0.6 * (domain == "respiratory")                         # respiratory/acid-base
renal = rng.normal(0, 1, n) + 0.3 * (domain == "cardiac")                            # renal/metabolic
heme = rng.normal(0, 1, n)                                                           # hematologic

# age-normative vitals (hand-specified smooth curves)
hr_mu = 150 - 45 * np.log1p(age) / np.log(19) - 20 * age / 18; hr_sd = 15 - 4 * age / 18
map_mu = 45 + 30 * np.log1p(age) / np.log(19); map_sd = 7 + 2 * age / 18
rr_mu = 45 - 22 * np.log1p(age) / np.log(19) - 6 * age / 18; rr_sd = 8 - 4 * age / 18

def obs(mu, sd, load):
    return mu + sd * (load + rng.normal(0, 1, n) * 0.8)

hr = obs(hr_mu, hr_sd, 0.7 * sev + 0.2 * resp)
mapv = obs(map_mu, map_sd, -0.6 * sev)
rr = obs(rr_mu, rr_sd, 0.4 * sev + 0.6 * resp)
temp = 37.1 + 0.6 * (0.3 * sev + 0.5 * (domain == "infectious") + rng.normal(0, 1, n))
spo2 = np.clip(97 - 3 * np.clip(0.5 * resp + 0.3 * sev + 0.6 * rng.normal(0, 1, n), 0, None) - 8 * (domain == "cardiac") * (rng.random(n) < 0.3), 60, 100)
ph = 7.39 - 0.05 * (0.5 * sev + 0.5 * resp + 0.7 * rng.normal(0, 1, n))
lactate = np.exp(0.3 + 0.35 * (0.8 * sev + 0.6 * rng.normal(0, 1, n)))
bicarb = 23 - 2.5 * (0.5 * sev + 0.4 * renal + 0.7 * rng.normal(0, 1, n))
cr_base = 0.25 + 0.6 * age / 18          # mg/dL, age-dependent baseline
creat = cr_base * np.exp(0.35 * (0.8 * renal + 0.3 * sev + 0.6 * rng.normal(0, 1, n)))
bun = np.exp(2.4 + 0.35 * (0.7 * renal + 0.3 * sev + 0.6 * rng.normal(0, 1, n)))
sodium = 139 + 3 * (0.3 * renal + rng.normal(0, 1, n))
potassium = 4.2 + 0.45 * (0.4 * renal + rng.normal(0, 1, n))
glucose = np.exp(4.65 + 0.25 * (0.4 * sev + rng.normal(0, 1, n)))
platelets = np.exp(5.5 + 0.4 * (-0.5 * sev + 0.6 * heme + 0.6 * rng.normal(0, 1, n)))
wbc = np.exp(2.3 + 0.35 * (0.4 * sev + 0.3 * heme + 0.7 * rng.normal(0, 1, n)))
hgb = 12 + 1.5 * (0.6 * heme - 0.3 * sev + 0.7 * rng.normal(0, 1, n))
hgb = np.where(age < 0.1, hgb + 3, hgb)

# first-24h median (raw) and "worst" value (eng); worst = median shifted by |noise|
def worst(med, scale, direction):
    return med + direction * np.abs(rng.normal(0, scale, n)) * (1 + 0.3 * np.clip(sev, 0, None))

# ---------------- missingness (MNAR: depends on acuity + panels co-missing) ----------------
def measured(base, sev_coef):
    return rng.random(n) < expit(base + sev_coef * sev)
vit_meas = {v: measured(b, 0.3) for v, b in [("hr", 5.0), ("map", 2.5), ("rr", 4.0), ("temp", 3.5), ("spo2", 4.0)]}
gas = measured(0.0, 1.0)      # blood gas panel: pH, lactate, bicarbonate
chem = measured(1.2, 0.8)     # chemistry: creatinine, BUN, Na, K, glucose
cbc = measured(1.5, 0.7)      # CBC: platelets, WBC, hemoglobin
def within(panel, p_keep):   # some component-level dropout inside a measured panel
    return panel & (rng.random(n) < p_keep)
meas = dict(vit_meas)
meas.update(ph=within(gas, 0.97), lactate=within(gas, 0.80), bicarbonate=within(gas, 0.95),
            creatinine=within(chem, 0.97), bun=within(chem, 0.95), sodium=within(chem, 0.98),
            potassium=within(chem, 0.97), glucose=within(chem, 0.90),
            platelets=within(cbc, 0.98), wbc=within(cbc, 0.98), hemoglobin=within(cbc, 0.98))
miss = {k: (~v).astype(int) for k, v in meas.items()}
# missing fraction over a broader 20-item physiology inventory (16 here + 4 unlisted items),
# so it is near-, but not exactly, a linear combination of the 16 masks
extra_missing = (rng.random((n, 4)) < expit(0.0 - 0.8 * sev)[:, None]).sum(1)
miss_frac = (sum(miss.values()) + extra_missing) / 20.0

def m(x, key):
    x = np.asarray(x, dtype=float).copy(); x[miss[key] == 1] = np.nan; return x

# ---------------- engineered features ----------------
hr_z = (hr - hr_mu) / hr_sd   # median HR z-score (age-normed)
map_z = (mapv - map_mu) / map_sd
rr_z = (rr - rr_mu) / rr_sd
temp_min = worst(temp, 0.4, -1); temp_max = worst(temp, 0.4, +1)
spo2_med = spo2
spo2_min = np.clip(worst(spo2, 2.0, -1), 50, 100)
ph_min = worst(ph, 0.03, -1); lact_max = worst(lactate, 0.5, +1); bic_min = worst(bicarb, 1.5, -1)
cr_max = worst(creat, 0.05, +1); bun_max = worst(bun, 2.0, +1)
na_min = worst(sodium, 1.5, -1); k_min = worst(potassium, 0.3, -1); glu_min = worst(glucose, 10, -1)
plt_min = worst(platelets, 25, -1); wbc_min = worst(wbc, 1.0, -1); hgb_min = worst(hgb, 0.6, -1)
ratio = cr_max / cr_base
aki = np.select([ratio >= 3, ratio >= 2, ratio >= 1.5], [3, 2, 1], 0)
vaso = (rng.random(n) < expit(-3.2 + 1.2 * sev + 1.0 * (domain == "cardiac") + 0.8 * surgery)).astype(int)

eng = pd.DataFrame({
    "sex": sex, "sex_missing": sex_missing, "age_years_at_intime": age, "age_band": age_band,
    "hours_admit_to_icu_intime": hours_admit, "hours_admit_to_icu_intime_missing": hours_admit_missing,
    "icu_intime_hour_of_day": hour, "icu_intime_is_night": is_night, "icu_intime_day_of_week": dow,
    "first_careunit": unit, "preliminary_diagnosis_chapter": chapter,
    "preliminary_diagnosis_chapter_missing": chapter_missing, "admission_clinical_domain": domain,
    "early_surgery_flag": surgery,
    "hr_zscore_median": m(hr_z, "hr"), "hr_zscore_missing": miss["hr"],
    "hr_tachycardia_flag": ((hr_z > 2) & (miss["hr"] == 0)).astype(int),
    "map_zscore_median": m(map_z, "map"), "map_zscore_missing": miss["map"],
    "map_hypotension_flag": ((map_z < -2) & (miss["map"] == 0)).astype(int),
    "rr_zscore_median": m(rr_z, "rr"), "rr_zscore_missing": miss["rr"],
    "rr_tachypnea_flag": ((rr_z > 2) & (miss["rr"] == 0)).astype(int),
    "temp_c_min": m(temp_min, "temp"), "temp_missing": miss["temp"],
    "temp_hypothermia_flag": ((temp_min < 36.0) & (miss["temp"] == 0)).astype(int),
    "temp_hyperthermia_flag": ((temp_max > 38.5) & (miss["temp"] == 0)).astype(int),
    "spo2_median": m(spo2_med, "spo2"), "spo2_missing": miss["spo2"],
    "spo2_hypoxia_flag": ((spo2_min < 90) & (miss["spo2"] == 0)).astype(int),
    "ph_min": m(ph_min, "ph"), "ph_missing": miss["ph"],
    "ph_acidosis_flag": ((ph_min < 7.30) & (miss["ph"] == 0)).astype(int),
    "lactate_max": m(lact_max, "lactate"), "lactate_missing": miss["lactate"],
    "lactate_elevated_flag": ((lact_max > 2.0) & (miss["lactate"] == 0)).astype(int),
    "bicarbonate_min": m(bic_min, "bicarbonate"), "bicarbonate_missing": miss["bicarbonate"],
    "creatinine_mgdl_max": m(cr_max, "creatinine"), "creatinine_missing": miss["creatinine"],
    "creatinine_aki_stage": np.where(miss["creatinine"] == 1, 0, aki),
    "bun_mgdl_max": m(bun_max, "bun"), "bun_missing": miss["bun"],
    "sodium_min": m(na_min, "sodium"), "sodium_missing": miss["sodium"],
    "potassium_min": m(k_min, "potassium"), "potassium_missing": miss["potassium"],
    "glucose_min": m(glu_min, "glucose"), "glucose_missing": miss["glucose"],
    "platelets_min": m(plt_min, "platelets"), "platelets_missing": miss["platelets"],
    "platelets_thrombocytopenia_flag": ((plt_min < 100) & (miss["platelets"] == 0)).astype(int),
    "wbc_min": m(wbc_min, "wbc"), "wbc_missing": miss["wbc"],
    "hemoglobin_min": m(hgb_min, "hemoglobin"), "hemoglobin_missing": miss["hemoglobin"],
    "early_vasopressor_inotrope_flag": vaso, "early_postop_surgery_flag": postop,
    "early_24h_physio_missing_fraction": miss_frac,
})
assert eng.shape[1] == 59

# ---------------- synthetic outcome (in-hospital death stand-in, ~4%) ----------------
lin = -3.9 + 0.9 * sev + 0.3 * resp + 0.25 * renal + 0.6 * vaso + 0.4 * (age < 0.0833) \
      + 0.3 * (miss["ph"] == 0) + 0.5 * (aki >= 2)
outcome = (rng.random(n) < expit(lin)).astype(int)

ids = pd.DataFrame({"subject_id": subj, "encounter_id": np.arange(1, n + 1), "split": split,
                    "outcome_synthetic": outcome})

# ---------------- minimally processed matrix ----------------
raw = pd.DataFrame({
    "sex": sex, "sex_missing": sex_missing, "age_years_at_intime": age,
    "hours_admit_to_icu_intime": hours_admit, "hours_admit_to_icu_intime_missing": hours_admit_missing,
    "icu_intime_hour_of_day": hour, "icu_intime_day_of_week": dow, "first_careunit": unit,
    "preliminary_diagnosis_chapter": chapter, "preliminary_diagnosis_chapter_missing": chapter_missing,
    "hr_median_bpm": m(hr, "hr"), "map_median_mmhg": m(mapv, "map"), "rr_median_bpm": m(rr, "rr"),
    "temp_median_c": m(temp, "temp"), "spo2_median_pct": m(spo2, "spo2"), "ph_median": m(ph, "ph"),
    "lactate_median_mmoll": m(lactate, "lactate"), "bicarbonate_median_mmoll": m(bicarb, "bicarbonate"),
    "creatinine_median_mgdl": m(creat, "creatinine"), "bun_median_mgdl": m(bun, "bun"),
    "sodium_median_mmoll": m(sodium, "sodium"), "potassium_median_mmoll": m(potassium, "potassium"),
    "glucose_median_mgdl": m(glucose, "glucose"), "platelets_median_k_ul": m(platelets, "platelets"),
    "wbc_median_k_ul": m(wbc, "wbc"), "hemoglobin_median_gdl": m(hgb, "hemoglobin"),
    **{f"{k}_missing": v for k, v in miss.items()},
    "vasopressor_inotrope_any_24h": vaso, "surgery_any_24h": surgery,
})

def rnd(df):
    return df.round({c: 4 for c in df.columns if df[c].dtype == float})
pd.concat([ids, rnd(eng)], axis=1).to_csv(os.path.join(OUT, "synthetic_X_eng.csv"), index=False)
pd.concat([ids, rnd(raw)], axis=1).to_csv(os.path.join(OUT, "synthetic_X_raw.csv"), index=False)
print("eng", eng.shape, "raw", raw.shape, "subjects", len(subject_ids),
      "outcome rate %.3f" % outcome.mean(), "split", pd.Series(split).value_counts().to_dict())
