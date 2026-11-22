# Synthetic stand-in data (v1)

> **⚠️ WARNING: THIS IS SYNTHETIC DATA. IT IS NOT PIC DATA AND CONTAINS NO REAL PATIENT INFORMATION. Every value was drawn from hand-specified distributions in `make_synthetic.py`. No real records or statistics fitted to real data were used. ANY RESULT COMPUTED ON THESE FILES IS MEANINGLESS. Use them only to write and test code.**

## Files
| File | Contents |
|---|---|
| `synthetic_X_eng.csv` | 4 ID/meta columns + the 59 frozen engineered columns ([N_encounters omitted] rows) |
| `synthetic_X_raw.csv` | the same 4 ID/meta columns + 44 minimally processed columns (same rows, same order, same split) |
| `make_synthetic.py` | generator (`python make_synthetic.py [outdir]`), fixed seed 6310, fully reproducible |

ID/meta columns, which must **never** be used as features: `subject_id`, `encounter_id`, `split` (train/val/test, assigned by subject, about 70/15/15), and `outcome_synthetic` (binary, about 5%, a stand-in for the outcome).

## Rules for team code
1. **Read every path from one config variable** at the top of your script, for example:
   ```python
   DATA = {"eng": "synthetic_X_eng.csv", "raw": "synthetic_X_raw.csv"}   # Vinai swaps in the real paths
   ```
   Don't hard-code paths anywhere else. Vinai will run the same code on the real files.
2. **Fit preprocessing on the training rows only** (`split == "train"`), then apply it unchanged to val/test and to any bootstrap or CV resample:
   - Imputation: training-set median for continuous columns. Missing values are left as NaN in both files, and the `*_missing` masks are already columns.
   - One-hot encoding of categoricals **with one reference level dropped**. Fix the level set from train. Treat a NaN category as its own level, or rely on its mask.
     - X_eng categoricals: `sex, age_band, icu_intime_day_of_week, first_careunit, preliminary_diagnosis_chapter, admission_clinical_domain`
     - X_raw categoricals: `sex, icu_intime_day_of_week, first_careunit, preliminary_diagnosis_chapter`
   - Standardization (center and scale to unit variance) using train means and SDs. Drop zero-variance columns.
3. Split by `subject_id`, never by row. Repeat encounters of one subject are always in the same split.
4. Don't tune anything on the test split.

## Data dictionary: X_eng (59)
Type key: C = continuous, B = binary flag, M = missingness mask (1 = missing), K = categorical, O = ordinal.

**Layer A (14)**

| Column | Type | Synthetic definition |
|---|---|---|
| sex | K | M/F (NaN when masked) |
| sex_missing | M | about 1% |
| age_years_at_intime | C | years, 0–18, infant-heavy mixture |
| age_band | K | neonate_lt1m, infant_1m_1y, child_1_4, child_5_11, adolescent_12_17 (deterministic from age) |
| hours_admit_to_icu_intime | C | lognormal hours |
| hours_admit_to_icu_intime_missing | M | about 4% |
| icu_intime_hour_of_day | C/K | integer 0–23 |
| icu_intime_is_night | B | 1 if hour ≥ 19 or < 7 (deterministic) |
| icu_intime_day_of_week | K | integer 0–6 (0 = Monday) |
| first_careunit | K | PICU, NICU, CICU, SICU, General_ICU |
| preliminary_diagnosis_chapter | K | ICD-10-style chapter label, mostly determined by domain (15% swapped) |
| preliminary_diagnosis_chapter_missing | M | about 5% |
| admission_clinical_domain | K | cardiac, respiratory, neurologic, infectious, surgical_other, neonatal, other |
| early_surgery_flag | B | probability depends on domain |

**Layer B (45).** Each physiologic variable has its value, a mask (`*_missing`), and, where named, a threshold flag. Flags are 0 when the value is missing.

| Variable | Value column (C) | Flag (B / O) and rule |
|---|---|---|
| HR | hr_zscore_median (age-normed z of the median) | hr_tachycardia_flag: z > 2 |
| MAP | map_zscore_median | map_hypotension_flag: z < −2 |
| RR | rr_zscore_median | rr_tachypnea_flag: z > 2 |
| Temp | temp_c_min (°C) | temp_hypothermia_flag: min < 36; temp_hyperthermia_flag: a hidden 24-h max > 38.5 |
| SpO2 | spo2_median (%) | spo2_hypoxia_flag: a hidden 24-h min < 90 |
| pH | ph_min | ph_acidosis_flag: < 7.30 |
| Lactate | lactate_max (mmol/L) | lactate_elevated_flag: > 2 |
| Bicarbonate | bicarbonate_min (mmol/L) | (none) |
| Creatinine | creatinine_mgdl_max | creatinine_aki_stage (O, 0–3): max / age-baseline ≥ 1.5 / 2 / 3 |
| BUN | bun_mgdl_max | (none) |
| Na, K, glucose | sodium_min, potassium_min, glucose_min | (none) |
| Platelets | platelets_min (10³/µL) | platelets_thrombocytopenia_flag: < 100 |
| WBC, Hgb | wbc_min (10³/µL), hemoglobin_min (g/dL) | (none) |

Masks: `hr_zscore_missing, map_zscore_missing, rr_zscore_missing, temp_missing, spo2_missing, ph_missing, lactate_missing, bicarbonate_missing, creatinine_missing, bun_missing, sodium_missing, potassium_missing, glucose_missing, platelets_missing, wbc_missing, hemoglobin_missing` (M).
Other columns: `early_vasopressor_inotrope_flag` (B), `early_postop_surgery_flag` (B), `early_24h_physio_missing_fraction` (C, 0–1).

## Data dictionary: X_raw (44, minimally processed)
- Demographics, timing, unit and diagnosis: `sex, sex_missing, age_years_at_intime, hours_admit_to_icu_intime, hours_admit_to_icu_intime_missing, icu_intime_hour_of_day, icu_intime_day_of_week, first_careunit, preliminary_diagnosis_chapter, preliminary_diagnosis_chapter_missing`. There is no age_band, is_night or clinical domain, because those are derived.
- 16 first-24-hour **medians in native units**: `hr_median_bpm, map_median_mmhg, rr_median_bpm, temp_median_c, spo2_median_pct, ph_median, lactate_median_mmoll, bicarbonate_median_mmoll, creatinine_median_mgdl, bun_median_mgdl, sodium_median_mmoll, potassium_median_mmoll, glucose_median_mgdl, platelets_median_k_ul, wbc_median_k_ul, hemoglobin_median_gdl`.
- 16 masks, `hr_missing … hemoglobin_missing`, identical in value to the X_eng masks.
- Raw exposure indicators: `vasopressor_inotrope_any_24h` and `surgery_any_24h`.

## How the structure was built (all hand-specified)
- **Latent factors.** Four independent N(0,1) factors drive the labs and vitals with hand-chosen loadings plus noise:
  - acuity `sev`: HR, MAP, RR, lactate, pH, bicarbonate, platelets, glucose, WBC
  - respiratory `resp`: RR, SpO2, pH
  - renal/metabolic `renal`: creatinine, BUN, bicarbonate, Na, K
  - hematologic `heme`: platelets, WBC, Hgb
- **Medians vs. worst values.** X_raw stores the latent "median". X_eng stores a worst value (min or max) equal to the median shifted by |noise|, with larger shifts for sicker patients. Vitals in X_eng are age-normed z-scores computed with hand-written normative curves. This gives the age × value interaction that a linear factorization of X_raw cannot reproduce.
- **Built-in redundancies, for the diagnostics to find:**
  - each value vs. its threshold flag (deterministic threshold rules)
  - `age_years_at_intime` vs. `age_band` (deterministic)
  - hour vs. `icu_intime_is_night` (deterministic)
  - `early_surgery_flag` vs. `early_postop_surgery_flag` (92% agreement)
  - `early_24h_physio_missing_fraction` = (16 masks + 4 unlisted hidden items) / 20, so it is a near-exact but not exact linear combination of the masks
  - co-missing lab panels: a blood-gas panel (pH, lactate, bicarbonate), a chemistry panel (creatinine, BUN, Na, K, glucose) and a CBC (platelets, WBC, Hgb). Each panel is measured or not as a block, with small within-panel dropout.
- **Missing not at random.** The probability that a panel or vital is measured rises with acuity, and the outcome depends on acuity and on whether pH was measured.
- **Outcome.** A logistic function of the latent factors, vasopressor use, neonatal age, AKI stage ≥ 2 and blood-gas measurement. About 5% positive.

## Sanity check (preprocessing only, run by the generator's author)
Training rows only, after one-hot encoding with a dropped reference level (NaN category as its own level), train-median imputation and train standardization:
- X_eng: 9,505 × 82, κ₂ ≈ 27.5
- X_raw: 9,505 × 59, κ₂ ≈ 17.3

Both are finite and moderate. The real data will differ.
