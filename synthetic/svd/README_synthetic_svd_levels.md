# Synthetic dataset for optimizing SVD algorithms

One synthetic PIC-like dataset in three related stages, plus a held-out answer key of planted latent factors. Built for algorithm work on feature matrices (SVD and related decompositions). **Not real patient data.**

- Generator: `make_synthetic_svd_dataset.py` (SHA-256 `5240655d149a4afba0a1803a6e798d56eb65f0ae5c6ffb5a9f8a2ef1106864b8`)
- Seed: `20261115`
- Subjects: **15,000**
- Cohort encounters: **17,239**
- Invalid DOB / age-invalid encounters: **103** (0.60% of encounters; safety cap is 1%). Breakdown: missing=21, unparseable=21, negative=61, valid=17136. Negative bins: {'lt1d': 21, 'gt30d': 20, '1to30d': 20}.
- Raw CSV total: **77.7 MB** (event density thinned to a 4.0 h vital step and a 36 h event horizon so the raw stage stays under ~300 MB).

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
| **z_severity** (illness) | Raises lactate (+1.05 mmol/L per SD), creatinine (+28.0 μmol/L), heart rate (+14.0 bpm); lowers SBP (-8.0 mmHg) and DBP (-5.0 mmHg), hence MAP; raises the logit of vasoactive exposure (intercept -1.15, slope +1.55) and the chance of multiple overlapping agents. Platelets and bilirubin are not in the v1 synthetic channel map, so those pSOFA components are absent. |
| **z_respiratory** | Lowers SpO2 (-3.8 pp) and PaO2 (-22.0 mmHg); raises FiO2 (+16.0 pp), which lowers derived P/F and S/F; raises respiratory rate (+7.0 /min). P/F uses same-draw pO2 + FiO2; S/F uses FiO2 without a valid P/F, nearest SpO2 within 5 min, only if SpO2 ≤ 97% (D15). |
| **z_age** | Sets DOB so age_years ≈ clip(9 + 6 z_age, 0.05, 21). Weight follows a pediatric age-weight curve times exp(0.10 z_age + N(0, 0.06)). Age-conditional norms: higher HR and RR, lower SBP/DBP/weight at younger ages (e.g. HR ≈ clip(155 − 4 age_years, 70, 180)). |

Observation noise (independent, per event): HR SD 5.0, RR 2.5, SBP 6.0, DBP 5.0, SpO2 1.4, lactate 0.22 mmol/L, creatinine 8.0 μmol/L, PaO2 8.0 mmHg, FiO2 3.5 pp.

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

- n = 17,239 encounters, p = 25 features
- Top 10 singular values: 335.355, 286.330, 227.657, 131.439, 130.894, 129.927, 128.755, 121.666, 118.752, 109.744
- Variance explained (top 10): 28.73%, 20.94%, 13.24%, 4.41%, 4.38%, 4.31%, 4.23%, 3.78%, 3.60%, 3.08%
- Top-3 variance explained (sum): 62.90%
- Gap σ₃/σ₄: 1.732
- Procrustes cosine similarity (severity, respiratory, age): 0.9768, 0.9993, 0.9908
- Principal angles (degrees): 1.739, 2.244, 14.492

The answer-key files are the ground truth for recovery experiments. Do not use them as features.

## How the stages were generated

```
python3 make_synthetic_svd_dataset.py --n-subjects 15000 --seed 20261115 --vitals-step 4.0
```

Fixed seed, no edits to the Stage 2 audit script, tests, or fixture. The fixture module is imported only for v1 map writers and schema constants.

## Row counts

### 1_raw

| file | rows | bytes |
|---|---:|---:|
| `ADMISSIONS.csv` | 17,239 | 1,404,291 |
| `CHARTEVENTS.csv` | 1,090,738 | 55,373,029 |
| `D_ITEMS.csv` | 13 | 257 |
| `D_LABITEMS.csv` | 17 | 343 |
| `ICUSTAYS.csv` | 17,239 | 1,230,165 |
| `INPUTEVENTS.csv` | 1 | 48 |
| `LABEVENTS.csv` | 305,414 | 15,978,772 |
| `MICROBIOLOGYEVENTS.csv` | 4,278 | 424,492 |
| `OUTPUTEVENTS.csv` | 1 | 56 |
| `PATIENTS.csv` | 15,000 | 528,309 |
| `PRESCRIPTIONS.csv` | 39,487 | 3,081,227 |
| `SURGERY_INFO.csv` | 3,072 | 428,419 |
| `canonical_icu_encounters.csv` | 17,239 | 1,327,468 |
| `drug_group_map_v1.csv` | 9 | 970 |
| `encounter_harmonization_fields.csv` | 17,239 | 1,503,273 |
| `feature_channel_map_v1.csv` | 33 | 5,238 |
| `neuro_admitting_dx_icd10_v1.csv` | 6 | 352 |
| `organism_map_v1.csv` | 5 | 495 |
| `specimen_map_v1.csv` | 5 | 529 |
| `subject_splits.csv` | 15,000 | 206,701 |
| `surgery_complexity_map_v1.csv` | 4 | 418 |

### 2_normalized

- `encounters.csv`: 17,239 rows (1 per encounter)
- `observations_long.csv`: 1,586,573 rows
- `medications_long.csv`: 39,487 rows

### 3_features

- `features_encounter_0_24h.csv`: 17,239 rows (1 per encounter)
- `features_svd_matrix.csv`: 17,239 rows, 25 analysis columns + encounter_id

## Layered tensor extension (A / B / C, amendment A1.2–A1.7)

A **new** companion dataset, `5_layered_tensor/`, follows the three-layer PIC pipeline (contract v12 + v12.1 A1.2–A1.7). Stages 1–3, the three static factor scores, and the original answer-key files in `4_answer_key/` were **not** modified.

See **[README_layered_tensor.md](README_layered_tensor.md)** for layers, binning, masking, outcome definitions, the LOS/eligibility limitation (existing stays end by ~280 h, so landmarks 288–336 h have zero eligibility), ground truth, and validation numbers.

- Layer C tensor: N=17239, T=56 (6 h bins, 0–336 h), F=75; `mask_valid` is 0 after each encounter's existing OUTTIME; `mask_obs` is per bin × feature.
- Landmark index slices `X_C[:, :L//6, :]`. 72 h sensitivity: `X_C[:, max(0, L//6-12):L//6, :]`.
- Outcomes are in `5_layered_tensor/3_layers/outcomes_by_landmark.csv`, generated from `outcomes_source.csv` (existing ADMISSIONS/ICUSTAYS timestamps were not overwritten).
- Trajectory and outcome-model ground truth: `4_answer_key/layered_tensor/` (new files only).
