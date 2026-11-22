# Layered tensor dataset (A / B / C) for landmark-forward models

Companion to the three-stage synthetic SVD dataset. This folder follows the **three-layer process of the PIC Foundations pipeline** (contract v12 + v12.1 amendment A1.2–A1.7). It does **not** replace stages 1–3 or the three planted static factors; those files were not modified.

- Generator: `make_synthetic_svd_layered_tensor.py`
- Seed: 20261115
- Subjects / encounters: 15,000 / 17,239 (same IDs as the existing dataset)
- Static factor scores: reused from `4_answer_key/latent_factors.csv` (unchanged)

## Why late landmarks are empty

Existing `ICUSTAYS` / canonical `OUTTIME` were **not** changed (same INTIME/OUTTIME/DOB as `1_raw`). The already-generated stays have

- LOS hours: min 20.0, p25 36.2, median 58.6, p75 110.2, p90 167.6, p99 267.8, **max 280.0**

so **eligibility at 288, 312 and 336 h is zero**, and Outcome 1 level 7 (remaining ICU time > 336 h) cannot occur. Core landmark 336 h therefore has no stacked rows. This is documented, not repaired by editing existing stays.

Eligibility (`last OUTTIME > t0 + L`):

| L (h) | n eligible | % of 17,239 |
|---:|---:|---:|
| 24 | 16,192 | 93.93% |
| 48 | 10,688 | 62.00% |
| 72 | 6,038 | 35.03% |
| 96 | 4,986 | 28.92% |
| 120 | 3,857 | 22.37% |
| 144 | 2,838 | 16.46% |
| 168 | 1,706 | 9.90% |
| 192 | 1,340 | 7.77% |
| 216 | 969 | 5.62% |
| 240 | 587 | 3.41% |
| 264 | 226 | 1.31% |
| 288 | 0 | 0.00% |
| 312 | 0 | 0.00% |
| 336 | 0 | 0.00% |

## Folder layout

- `5_layered_tensor/1_extended_raw/` — companion raw events over each encounter's **whole existing stay** (plus a 6 h pre-ICU lookback). PATIENTS, ADMISSIONS, ICUSTAYS, canonical, maps and splits are **copies** of the originals (IDs and timestamps unchanged). New CHARTEVENTS / LABEVENTS / PRESCRIPTIONS are driven by `severity_t`. `outcomes_source.csv` holds death, same-HADM readmission and hospital discharge times used for A1.2 outcomes; existing ADMISSIONS.DISCHTIME and canonical `Died` are left as they were.
- `5_layered_tensor/2_normalized_extended/` — same v1 maps and D1–D20 rules as `2_normalized`.
- `5_layered_tensor/3_layers/` — Layer A CSV, Layer B CSV, Layer C npz, landmark index, outcomes-by-landmark, optional long gz.
- `4_answer_key/layered_tensor/` — **new** ground-truth files only (trajectories, shape basis, outcome coefficients, validation). Existing answer-key CSVs and zip were not modified.

## Time-varying ground truth (drives extended raw)

`severity_t(e, t) = s_e + d_e(t)`, with `s_e` the existing static `z_severity`.

`d_e(t) = amplitude_e * shape_{k(e)}(t)` on 6 h bin midpoints, k ∈ {stable, improving, worsening, spike_recovery, late_deterioration} (i.i.d. uniform), amplitude ~ Uniform(1.1, 2.4) for non-stable classes (0 if stable).

| Class | Raw shape (then mean-centered and L2-normalized on the 56 midpoints) |
|---|---|
| stable | 0 |
| improving | exp(-2 u), u = t/336 |
| worsening | logistic 1/(1+exp(-8(u-0.5))) |
| spike_recovery | Gaussian bump, peak ~ Uniform(12, 96) h, width 18 h |
| late_deterioration | 0 then logistic rise after onset ~ Uniform(48, 200) h |

Respiratory: `resp_t(e,t) = z_respiratory + 0.4 * d_e(t)` (same shape, milder). Age remains a static patient-level factor (DOB unchanged). Event values use the same raw-scale loadings as the static generator (`BETA` in `make_synthetic_svd_dataset.py`), with means shifted by `severity_t` / `resp_t` at the event time (pre-t0 uses the static scores only). Vasoactive starts are also tied to `severity_t`.

## Layers (A1.6)

**Layer A** — one row per encounter, fixed at t0. Age = INTIME − DOB with `age_invalid` (D20), sex, first care unit, admitting diagnosis domain from `ADMISSIONS.ICD10_CODE_CN` (3-character → `icd_domain`), admission source derived from `EDREGTIME` (ed / unknown), recent surgery (anesthesia-end in [t0 − 24 h, t0 + 6 h], A2.1), remote same-HADM surgery, weight (charted or imputed, A2.12/D20), and the Layer A pre-ICU laboratory block: last released value in **[t0 − 6 h, t0)** plus a mask (A2.11). Events before t0 never enter Layer B or C.

**Layer B** — one row per encounter, computed once on **[t0, t0 + 24 h)** from the extended normalized stay, held fixed at every landmark. Summaries match stage 3: last/min/max/mean/count, pSOFA respiratory (P/F then S/F, D19/D15), MAP, lactate, creatinine, vasoactive flag and concurrent count, drug-group flags, fluid intake mL/kg/day, missingness. Because the extended raw is a new companion stay (full horizon, trajectory-driven), Layer B is **not** a copy of `3_features/`; that file is unchanged.

**Layer C** — tensor `X_C` with shape **(N, 56, F_C)** = (17239, 56, 75). 6 h bins on [0, 336) h; bin *b* is `[6b, 6b+6)`. Observation time is `t_hours` (labs already include the +1 h release lag). Per-bin features: channel means, time-since-last (capped at 24 h, A2.13) and last-carried indicators, drug-group active flags, vasoactive concurrent count, pSOFA respiratory. **Variable-length handling:** `mask_valid[e, b] = 1` iff the bin starts before that encounter's last OUTTIME (`6b < LOS`); `mask_obs[e, b, f] = 1` iff a value was observed (or a flag is defined) in a valid bin. Invalid bins are NaN in `X_C` (no duration leak through padding).

**Landmark index** (`landmark_index.csv`): one row per (encounter, L) with `OUTTIME > t0 + L`, L ∈ {24, 48, …, 336}, and `n_bins = L/6`. Slice the sequence as `X_C[:, :L//6, :]`.

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

Cohort mortality in `outcomes_source`: 7.66%. Readmission rate: 6.97%.

### Outcome counts at core landmarks

| L | n | mortality | Outcome 1 (8/7/6/5/4/3/2/1/0) | Outcome 2 (4/3/2/1/0) |
|---:|---:|---:|---|---|
| 24 | 16,192 | 7.6% | 1233/0/210/1022/1396/1986/954/4296/5095 | 1233/1125/227/1427/12180 |
| 48 | 10,688 | 7.7% | 824/0/0/889/682/2001/1042/954/4296 | 824/749/150/942/8023 |
| 72 | 6,038 | 7.8% | 470/0/0/544/688/1396/944/1042/954 | 470/425/85/512/4546 |
| 120 | 3,857 | 7.4% | 285/0/0/0/544/688/339/1057/944 | 285/282/56/333/2901 |
| 168 | 1,706 | 7.9% | 135/0/0/0/0/544/345/343/339 | 135/125/26/148/1272 |
| 240 | 587 | 7.3% | 43/0/0/0/0/0/0/210/334 | 43/39/9/47/449 |
| 336 | 0 | n/a | 0/0/0/0/0/0/0/0/0 | 0/0/0/0/0 |

## Validation (HOSVD / class / outcomes)

HOSVD was run on the **complete 72 h window** (encounters with LOS > 72 h, first 12 bins), standardized on observed cells, missing features zero-filled *inside that complete window only*, so stay-length padding cannot leak.

- 72h complete window (LOS>72, first 12 bins): n=6,038. Mode-2 SVs: 653.67, 282.81, 193.78, 149.74, 136.53, 124.80, 124.14, 121.21. Principal angles (deg) vs planted shapes, U2[:,1:5] ⊥ constant: 6.35, 13.85, 39.35, 89.17 (max 89.17). Trajectory OLS R²=0.032, corr=0.179. k-means ARI=0.016 (excl. late_deterioration: 0.016856578123998604).
- 120h complete window (LOS>120, first 20 bins): n=3,857. Mode-2 SVs: 647.40, 285.40, 181.11, 143.08, 120.17, 115.16, 106.94, 103.50. Principal angles (deg) vs planted shapes, U2[:,1:5] ⊥ constant: 7.98, 9.62, 31.45, 69.29 (max 69.29). Trajectory OLS R²=0.050, corr=0.224. k-means ARI=0.018 (excl. late_deterioration: 0.02035729741648729).
- 168h complete window (LOS>168, first 28 bins): n=1,706. Mode-2 SVs: 501.11, 220.01, 137.24, 107.00, 89.41, 82.13, 77.07, 72.15. Principal angles (deg) vs planted shapes, U2[:,1:5] ⊥ constant: 7.60, 10.55, 36.62, 79.18 (max 79.18). Trajectory OLS R²=0.057, corr=0.239. k-means ARI=0.016 (excl. late_deterioration: 0.015517637204164665).


## Formats

- `layer_c.npz`: `X_C` float32 (N, 56, F_C), `mask_obs` uint8, `mask_valid` uint8 (N, 56), `encounter_ids`, `bin_edges`, `bin_midpoints`, `feature_names`.
- Layer A / B CSVs; `landmark_index.csv`; `outcomes_by_landmark.csv`; `outcomes_source.csv` (in 1_extended_raw).
- Long Layer C CSV.gz if the observed-cell count permits.
- Ground truth: `4_answer_key/layered_tensor/`.

Splits are the existing `subject_splits.csv`.

**Empirical planted time-subspace (additional).** Because spike and late-deterioration onsets vary by encounter, a fairer check is the principal angles between mode-2 of Layer C and the SVD of the actual planted `d_e(t)` matrix on the same window (not the four canonical curves at a single default onset). In the 72 h window the planted trajectories are essentially rank-2 (96.9% of `d_e` variance in two components). Angles vs mode-2 after dropping a constant: **19.5° and 45.4°**. At 120 h (rank 3): **11.4°, 31.2°, 80.8°**. At 168 h: **10.8°, 46.7°, 84.9°**. So the leading time-varying direction is partially recovered; later shape directions are not cleanly separated from other temporal structure in Layer C (observation noise, missing labs after 24 h, and variable onsets). Class ARI stays near chance. These numbers are reported honestly; the ground-truth files still allow algorithm work with a known `d_e(t)`.
