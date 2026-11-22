# PIC dataset audit for Foundations Plan v12: `pic_dataset_audit_v12.py`

Aggregate-only audit that settles the bracketed v12 audit items (spec: `revised_plan_v12.md`, `v12_change_list.md`). It writes counts, percentages and quantiles only. It never writes SUBJECT_ID, HADM_ID, ICUSTAY_ID, encounter_id, or any row-level value.

**Status:** written and checked for syntax (`py_compile` and pyflakes are clean). It has **not** been run against PIC, and the self-tests have not been run yet. Run the self-tests first.

## How to run (locally, credentialed environment, VS Code terminal)

1. Put the script in the project, for example `<project>/audit/pic_dataset_audit_v12.py`. The project root is found by walking up to a folder that has `configs/` or `data/interim/`, or you can pass `--project-root`.
2. Check the `CONFIG` block at the top:
   - **Paths:** `COHORT_FILE` (default `data/interim/canonical_icu_encounters.csv`, per contract §1.1), `COHORT_COLUMNS`, `SPLIT_FILE` (`data/processed/subject_splits.csv`), `DATA_SOURCE_KEYS`.
   - **Channels:** ITEMID lists in `CHANNELS`, which come from contract §11 where known; the other Table 5.1 items are filled from the dictionary search.
   - **Keywords:** keyword lists in `CHANNELS`, `PATHOGEN_FAMILIES`, `DRUG_GROUPS`.
   - **Other choices:** plausibility ranges (in native units; most are DRAFT), `X_GRID_HOURS`, and the Appendix A domain choices.
3. Run the commands:
   ```
   python audit/pic_dataset_audit_v12.py --self-test     # synthetic tests only; reads no data
   pytest audit/pic_dataset_audit_v12.py                 # same tests via pytest
   python audit/pic_dataset_audit_v12.py                 # full audit (A,B,C,D,E)
   python audit/pic_dataset_audit_v12.py --sections A,B  # a subset
   python audit/pic_dataset_audit_v12.py --raw-root "D:/PICDataset" --chunksize 250000
   ```
4. Outputs: `data/interim/pic_dataset_audit_v12.json` (full) and `data/interim/pic_dataset_audit_v12.md` (readable summary).

**Raw root discovery:** first `PIC_ROOT` env / `--raw-root`; then `configs/data_source.json` (any key in `DATA_SOURCE_KEYS`, then any string value that is a folder); then `data/raw`, `../PICDataset`, and similar. The script searches up to 2 levels down for `ICUSTAYS.csv` or `ICUSTAYS.csv.gz`.

**Tables:** `.csv` and `.csv.gz` are both read. File and column names are matched regardless of case. CHARTEVENTS, LABEVENTS, OUTPUTEVENTS, INPUTEVENTS, PRESCRIPTIONS and MICROBIOLOGYEVENTS are read in chunks, loading only the columns needed (`usecols`).

**Column-role matching:** roles are matched on whole words, so "age" no longer matches `STAGE` or `DOSAGE`, "dt" no longer matches `WIDTH` or `ADT_FLAG`, and "los" no longer matches `CLOSE`. Self-test `test_role_match_word_boundaries` covers this.

**Privacy:**
- Counts from 1 to 9 are shown as "<10", and a percentage built on such a count is shown as "suppressed".
- Quantiles need n ≥ 10 and never include min or max.
- Complementary suppression (hiding a cell that could be worked out from a total) is **not** applied, so review the `.md` file before sharing.

**Linking events to encounters:** an event belongs to an encounter by SUBJECT_ID and time, if INTIME − 24 h ≤ t < min(OUTTIME, INTIME + 168 h). The windows don't overlap because encounters are more than 24 h apart (contract §4.2). Hours before INTIME are used only for weight and are counted separately.

## Sections and what each one resolves

| Section | What it reports | v12 bracket / open item it resolves |
|---|---|---|
| Cohort | N encounters and subjects, split counts vs [N_encounters omitted] and the frozen 70/15/15 split (sizes omitted from this public copy); source of the cohort and split | plan §2 cohort |
| A1 interval | INTIME − min(linked ADMITTIME): missing, negatives by size, quantiles, v12 bands, share ≥ X for X ∈ {6, 12, 24, 48, 72, 168} h, exact boundary hits | [handling of negative values: audit], [X set at audit], check A(5) |
| A2 ED | Which tables have an ED registration column; EDREGTIME present, before / at / after INTIME; lag | [EDREGTIME source table], ED rule of the admission source |
| A3 surgery | Begin time before / exactly at / after INTIME (<24 h, ≥24 h); counts within ±1/15/60 min (timestamp rounding); name or code columns | surgery at exactly INTIME → Layer A only (A(13), B(13)); [does surgery at INTIME count as OR: decide]; [surgery name / code field] |
| A4 diagnosis missing | ICD10_CODE_CN missing overall and by split; free-text DIAGNOSIS present when the code is missing | [diagnosis missingness indicator] |
| A5 sex | GENDER values; missing count; values other than M/F | [count of empty sex] (no sex_missing in v12) |
| A6 ICD mapping | Codes not in the dictionary; mapped codes that aren't valid ICD-10; CN 3-char vs ICD-10 3-char agreement; categories with fewer than 20 train patients (with sensitivity at 10/30/50); val/test encounters with unseen categories; first letter vs WHO chapter cross-tab (D00–D48 vs D50–D89, H00–H59 vs H60–H95); letters that span domains (Q, D, I) and rare categories per letter×domain combination | Chinese-ICD mapping audit (A(14)); [m = 20]; [chapter = letter or WHO]; pooling within letter×domain [confirm]; Appendix A I and D00–D49 |
| A7 care unit | first_careunit counts; contract domain × diagnosis-only domain cross-tab; chapter I and Q20–Q28 by care unit | first_careunit; check A(8) |
| A8 prior ICU | Number of earlier cohort encounters and earlier raw ICU stays | [prior ICU admissions: include or omit] |
| A9 admission source | floor / OR / ED / other at each X; a variant where surgery at INTIME counts as OR; conflicts (ED registration and interval ≥ X; surgery and interval ≥ X) | [X set at audit], check A(12) |
| A10 timing | v12 L1 shifts, weekday/weekend, L2 categories, exact 08:00 / 15:00 / 23:00 hits | check A(5) |
| B dictionary | Keyword search of D_LABITEMS and D_ITEMS (all label-like fields, English and Chinese) for every Table 5.1 channel, plus GCS, pupils, FiO2, ventilation, weight, height, cortisol, TSH, free T4, insulin, ketones, osmolality, stool occult blood, urine sodium, and the pathogen families. Each item is flagged configured / single candidate / ambiguous / not found, with a fluid check. Serum and urine sodium are listed side by side | [ITEMID list per channel]; GCS / pupils / FiO2 / ventilation / weight search (B(25), stage 07a); endocrine labs; stool occult blood; mis-mapped urine sodium; troponin I |
| C coverage | Per channel and per candidate ITEMID: share of encounters covered 0–24 h (all, train, and window ≥24 h); share with ≥3 values in 24 h; coverage in each 6 h bin to 168 h (denominator = open bins); units; numeric vs text values; records below / above the range; value quantiles; any-channel coverage per system | 0–24 h coverage on this cohort (the plan's coverage figures came from Study C); [GI standalone]; variability ≥3 values (B(21)); unit table (B(14)) |
| C panels | Labs or vitals recorded at the same (subject, CHARTTIME): share recorded alongside ≥2 / ≥5 other channels, top co-recorded partners | data-quality flag "panel-driven coverage"; deciding which tests count as triggered |
| C HR gaps | Gaps between HR records (0–168 h and 0–24 h), share >1/2/4/6/12/24 h, hours to first HR | data-quality flag "HR 76% suggests CHARTEVENTS gaps" |
| C weight / fluids | Weight in −24 to 24 h and at any time; weight + ≥2 urine records (mL/kg/h feasible); weight + intake + output (% fluid overload feasible); urine output and intake records per encounter, gaps, 24 h totals, units | [weight source and missingness] (B(26)); [urine output rate definition] (S2, conflict #5) |
| D release timing | Time columns beyond CHARTTIME in LABEVENTS, CHARTEVENTS (STORETIME), OUTPUTEVENTS, INPUTEVENTS and MICROBIOLOGYEVENTS, with lag from CHARTTIME (missing / negative / zero / quantiles); pathogen-family rows and encounters in MICROBIOLOGYEVENTS | S9 (conflict #8); pathogen positive-flag timing (conflict #26) |
| E drugs | PRESCRIPTIONS name and dose columns; per group (vasoactive, NMB, antiepileptic, dual-use benzodiazepines): agents matched, top names, encounters 0–24 h and 0–168 h, share of dose fields filled, dose units, ENDDATE missing; overlaps between groups; dose or rate columns in other tables; share of STARTDATE that is date-only | vasoactive dose for pSOFA (conflict #25); [NMB and antiepileptic drug groups: decide] (conflict #7); ENDDATE missingness |

## Bracketed in the script (to confirm on first run)

- **Cohort file columns:** `encounter_id`, `SUBJECT_ID`, `INTIME`, `OUTTIME`, `HADM_ID`, `hadm_id_list` (JSON array). The split column is assumed to be named `split` in `subject_splits.csv`. If the canonical file is missing, the cohort is rebuilt in memory with the contract §4.2 24 h collapse; this is flagged in the notes and nothing is written to disk.
- **`configs/data_source.json`:** the key that holds the raw root is unknown.
- **EDREGTIME:** assumed to be in ADMISSIONS; taken from the earliest linked admission (`ED_SOURCE_ROW`).
- **SURGERY_INFO:** name and procedure-code columns are found by role, so check what was picked. If HADM_ID is missing, surgeries are linked by time.
- **Chapter definition:** the default is the first letter; the WHO-chapter alternative is reported side by side.
- **Plausibility ranges:** contract gates were converted to native units (creatinine µmol/L, BUN mmol/L, Hb g/L). All other ranges are DRAFT.
- **Auto-filled channels:** channels without contract ITEMIDs use up to 8 dictionary candidates (`origin: auto`). The channel-level union can mix serum and urine or CSF items, so choose from the per-ITEMID rows and record the choice in `CONFIG`.
- **Pathogen family names in PIC:** especially "bacteremia PCR".
- **INPUTEVENTS:** may have no ITEMID or drug; vasoactive rates may not exist anywhere.
- **Panel groups:** a group split across two read chunks is counted twice (minor).
