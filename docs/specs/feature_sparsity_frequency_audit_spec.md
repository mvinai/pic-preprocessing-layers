# Spec: Feature sparsity and measurement-frequency audit, Layers A, B and C (read-only, aggregate-only)
**Version 0.9.2 DRAFT, Nov 15, 2026.** v0.1 (Nov 12, 2026) listed 12 open questions; **Vinai answered all of them on Nov 14, 2026** (§15, "Decisions"), and confirmed the five follow-up points the same day (§15, "Follow-up confirmations"). The Stage 1 item inventory has been **run by Vinai (Nov 14, 2026; release gate ok)**. Its four map templates were filled as proposals, reviewed, and **decided by Vinai (decisions D1–D13, §15, "Map decisions")**. The resulting **v1 maps are hash-pinned in §4.2**, and the Stage 2 checks they require are listed in §5.3. Vinai confirmed the neuro admitting-diagnosis set the same day (v1, pinned in §4.2; D3), and then settled the remaining loose ends (D14–D18, §15): the new-group cutoffs, the S/F cut points, promethazine/chlorpromazine, fluid status and a new electrolyte_replacement group. The channel and drug-group maps were rebuilt and re-pinned. Vinai then decided that P/F is scored like S/F, on the pSOFA cut points without the respiratory-support requirement (D19). The Stage 2 script (`feature_sparsity_frequency_audit.py` v1.0) was written and hash-gates the pinned inputs. After the v1.0 run stopped on 42 encounters with missing or invalid DOB, Vinai decided (D20, Nov 15, 2026) to continue with an `age_invalid` band and an aggregate DOB breakdown; the script is now v1.0.1. Vinai ran v1.0.1 (release gate ok), reviewed its results and decided **D21–D37** the same day (§15, "Stage 2 review decisions"): the channel map is now **v2** (§4.2), the four raw tables that Stage 1 did not hash are now **gated** (§4.3, D36), and two questions wait on an aggregate follow-up script (§5.4): whether bilirubin 5255 is a NICU-specific total bilirubin (D24) and how FiO2 5222 is recorded (D29). Vinai then confirmed the four remaining map details the same day (v0.9.1: hemoglobin unit placement, troponin T and ionized-calcium gates, blood-gas mask membership); channel map v2 was revised in place before any run used it and re-pinned. The follow-up checks then ran (release gate ok), and Vinai took the final decisions the same day (v0.9.2): **P/F and the pSOFA respiratory component are dropped** (D29; pSOFA is a **5-component score without the respiratory component** (D29): cardiovascular, coagulation, hepatic, renal and neurological), and bilirubin 5255 is merged into `bilirubin_total` **only in NICU** (D24). They are recorded in **channel map v2.1** (§4.2).

All §9 thresholds are **approved as starting defaults** (decision 10). Vinai reviews the per-channel Stage 2 results before any disposition is locked. The script applies no decision: it reports, and Vinai decides on the train split.

## 1. Purpose
Before Layer A/B/C features are built under contract v12.1, measure for every candidate feature how often it is measured, when, and how regularly. Three questions are answered per feature:
1. **Disposition:** drop it, keep it as a measured indicator only, or keep its value plus mask plus time since last measurement (§9.1).
2. **Layer C bin width:** whether 6 h bins and the contract LOCF caps fit how often the channel is measured (§9.2).
3. **Imputation:** which sparsity class the channel falls into, and therefore which imputation arms may touch it (§9.3).

It also quantifies how the 78% missing medication `ENDDATE` affects exposure duration (§9.4). It also completes the contract checklist §9 [AUDIT] item "coverage recomputed per encounter for 0 to 24 h on the [N_encounters omitted] cohort".

This is a design audit only. It builds no features, fits no models, changes no frozen file, and writes no patient-level output that is shared.

## 2. Sources this spec is grounded in
- Contract v12 (`contract/foundations_contract_v12.md`): H2 windows, H3 event linking and event times, H7 fields, H8 table notes, the Layer A/B/C level tables and the banned list (I7).
- Contract v12.1 amendment (`contract/v12_1_amendment_landmark_forward_outcome_DRAFT.md`):
  - A1.4 landmark grid and core landmarks;
  - A1.6 feature layers at landmark L;
  - decision 11 and A2.3 (labs at release time);
  - A2.1 surgery timing;
  - A2.4 Layer C bins (4 to 56 bins);
  - A2.7 surgery φ1 band.
- Revised plan v12:
  - §4–§6 for the feature lists;
  - Table 5.1 (ten organ systems);
  - the medication state machine (S3), the microbiology timing hierarchy (S8, M1–M3) and the lab release rule (S9);
  - the LOCF table;
  - the v11 imputation comparison arms (r), (a), (b), (c);
  - S2 (measurement counts are audit descriptors only).
- Preprocessing log: Stage 2 (pre-ICU lookback, 6 h primary / 24 h sensitivity), Stage 3 (gap handling), Stage 4 (landmark eligibility counts), Step 5.1 (lab release delay), Step 5.2 (surgery).
- The PIC_Analysis_3_0 Layer A feature spec and modeling contract §3 (16 Layer C channels, 6 event columns, missingness flags kept by default).
- The Stage 1 item-inventory outputs (Nov 14, 2026; aggregate, release gate ok) and the v1 maps decided from them (§4.2).
- Earlier audit specs for conventions: the duration audit v1.7, the band-support audit v1.0 and the landmark-forward band-support audit.

No feature is invented here. Channels named in the plan but not yet mapped to ITEMIDs are audited only through the pinned channel map (§4.2). Anything not in that map is reported in the item inventory, not as a feature.

## 3. Design rules (inherited)
- **Cohort and split:** [N_encounters omitted] canonical encounters, subject-grouped split (the frozen 70/15/15 split (sizes omitted from this public copy)). Never rebuilt.
- **Everything is stratified by split.** Decisions are made on **train only**; validation and test are descriptive. Every train-derived threshold is reported with the val/test values beside it, never refit on them.
- **Clock:** t0 = harmonized first INTIME; `OUT` = last (maximum) OUTTIME. One continuous encounter clock, including merge gaps (Stage 3).
- **Event linking (H3):** `SUBJECT_ID` and `HADM_ID ∈ hadm_id_list`, then time. ICUSTAY_ID disagreements are counted; the time-based assignment wins. Events outside every window are counted per table.
- **Event times (H3.4, A2.1, A2.3):**
  - CHARTEVENTS, INPUTEVENTS, OUTPUTEVENTS: `CHARTTIME`;
  - LABEVENTS: **release time = CHARTTIME + 1 h** (primary) everywhere: Layer B window, Layer C bins, landmark availability and the pre-ICU lookback (including the approved Layer A pre-ICU lab block, decision 11). The + 0 h and + 2 h sensitivities are reported only for the availability metrics (§7.2, §7.6);
  - MICROBIOLOGYEVENTS: `CHARTTIME` = collection (M1); M3 availability = CHARTTIME + turnaround;
  - PRESCRIPTIONS: `STARTDATE` (and `ENDDATE`) through the S3 state machine;
  - SURGERY_INFO: timing time T = `ANES_END_TIME`, falling back to `SURGERY_END_TIME` (A2.1). **The Layer C surgery event is placed in the bin of T** (decision 8, Nov 14, 2026; the begin-time bin is a sensitivity analysis).
  - `STORETIME` is never used.
- **Windows:**
  - Layer A: at t0, with the approved pre-ICU lookback [t0 − 6 h, t0) primary and [t0 − 24 h, t0) as a sensitivity analysis, capped by admission ID;
  - Layer B: [t0, min(OUT, t0 + 24 h));
  - Layer C: 6 h bins k = [t0 + 6k, t0 + 6k + 6), k = 0..55, truncated at OUT. An observation at exactly t0 + 24 h is in bin 4.
- **Landmarks (A1.4):** L ∈ {24, 48, 72, 120, 168, 240, 336} h are the core landmarks reported here, at prefix bins 4, 8, 12, 20, 28, 40 and 56.
  - Eligible at L means `OUT` > t0 + L (strict).
  - The script must reproduce the published eligible counts exactly or stop: overall **[landmark eligible counts omitted]**, and the published split counts.
  - The other 7 grid landmarks (every 24 h) are reported only in the per-bin occupancy curve (§7.7).
- **Intraoperative data:** `SURGERY_VITAL_SIGNS` is not a Layer B/C vital channel (H8). It is out of scope here except as a surgery descriptor, which the surgery audits already cover.
- **Excluded tables:** `OR_EXAM_REPORTS` is excluded from all matrices. `EMR_SYMPTOMS` is not a planned feature source; it gets an inventory line only.
- **Banned:** no LOS-derived quantity, OUTTIME, DISCHTIME or outcome field is a feature (I7). This audit uses `OUT` only to define windows and eligibility, as the landmark audit did.

## 4. Inputs (read-only)
### 4.1 Pinned derived files (exact full SHA-256, checked before the file is read; any mismatch stops the run)
- canonical encounters `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418` (gated in Stage 1 too, because `hadm_id_list` comes from it; confirmed Nov 14, 2026);
- harmonization fields `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`;
- subject split `60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31`;
- surgery complexity map v1 `5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5` (surgery system descriptors);
- optional: the frozen φ1 medians `phi1_band_support_frozen_medians_MAP_v1.csv` `e73c8b5a8d780b51edfd0935abdb05335f57e3a0cfc0c4c52eaaa85888db0047`, used only for the surgery band descriptor.

### 4.2 Pinned maps (v1, decided by Vinai Nov 14, 2026, D1–D19; channel map v2, Nov 15, 2026, D21–D37; full SHA-256 checked before use)
The four Stage 1 templates were filled as proposals from the shared, aggregate Stage 1 outputs only, reviewed, and every group-level call was **decided by Vinai** (D1–D19, §15). Copilot never assigns items, drugs, specimens or organisms to groups: Stage 2 only reads these pinned maps.

| Map | File | SHA-256 |
|---|---|---|
| Channel map | `feature_channel_map_v1.csv` | `56e8bf72fcbcaac765a39c2e705aedef61323e94841dbdcc8b2657191af370b6` |
| Drug-group map | `drug_group_map_v1.csv` | `6fa7d8a935f03ad71b4a7ecd8276f9a54f3a679f4cfc08c8b4504eda83455ead` |
| Specimen map | `specimen_map_v1.csv` | `177b165ad8e5888a3999583a4949e93908be45af34cfa9c6e4937c40c08666f1` |
| Organism map | `organism_map_v1.csv` | `21f28448cb31c654c5312ca97883a34fe5cb6b269e8323973eb8b899bb3583d2` |
| Neuro admitting-diagnosis set (confirmed by Vinai Nov 14, 2026) | `neuro_admitting_dx_icd10_v1.csv` | `3cc8f4acfe9df29037b60d903ad4204c40db5bc76923de8b64234ccb2553b9d9` |
| **Channel map v2 (D21–D37, Nov 15, 2026)** | `feature_channel_map_v2.csv` | `cbfdc68f11bc2f8d9f96150803a1339930da81eb091f481576f709debdd5ab5c` |
| **Channel map v2.1 (D24, D29 final; current)** | `feature_channel_map_v2_1.csv` | `7d1c060a228ab0772f35ae27e43a30971bcc13c36c98fcce2e369640b91bc899` |

- **Which pin applies.** Stage 2 v1.0.1 ran on channel map v1 and keeps it as its pin (its outputs are released). Channel map v2 is the pin of the Stage 2 follow-up checks (§5.4). **Channel map v2.1 (`feature_channel_map_v2_1.csv`, SHA-256 `7d1c060a228ab0772f35ae27e43a30971bcc13c36c98fcce2e369640b91bc899`) is the current pin for every later step (Stage 3 feature building).** The drug-group, specimen, organism and neuro files are unchanged (v1 pins stand). `maps_v2_sha256.json` records the v2.1, v2 and other pins.
- **Channel map versions.** The follow-up run consumed the first v2 revision (`9872e688…`); v2 was then revised in place (`cbfdc68f…`, the four confirmed details), which the follow-up script now pins. The final D24 / D29 decisions are a new file, v2.1, so no pinned file is overwritten again and the script does not need re-pinning. v2.1 changes against v2: P/F (`derived:pf_ratio`), FiO2 5222 and pO2 5244 set to `not_used` (5244 also leaves the `blood_gas` mask group, now 15 items); pO2 = 5239 only; SaO2 5252 stays its own channel; bilirubin 5075 main, 5255 `assigned_merged_conditional` (NICU only, source flag).
- **Channel map v2** (686 rows; built from v1 by a script that checks the v1 hash first): every D21–D37 change is written in the row's `decision` and `notes`.
  - New column `mask_group`: one shared mask per panel (D33): `cytokines` (4 items), `teg` (9), `urine_dipstick` (4), `blood_gas` (16 blood-gas analyser items; ionized calcium 5215 and blood-gas hemoglobin 5257 keep their own masks).
  - New statuses: `assigned_merged` (merge decided, D22) and `assigned_merged_indicator_only` (merged for the measured indicator; values not used, D26). `assigned_merge_pending_stage2` now remains only for bilirubin 5075/5255 (D24).
  - New channel `calcium_ionized` (5215, mmol/L, D21) and new derived row `kdigo_uo_stage` (D30). Dropped to `not_used`: creatinine 5032, glucose 5047, ferritin 6161 (D23), anion gap 5212/5213 (D25), lipase 6473 (D26), S/F (D28), analyser P/F 5243 (D33).

- **Channel map** (685 rows): table, ITEMID → channel, Table 5.1 system, unit conversion, plausibility gate, flags `triggered` and `contract16`, plus `source_itemid` (provenance, kept per value), `merge_rule`, `stage2_check`, `status` and `decision`.
  - Status values: `assigned`; `assigned_merge_pending_stage2` (merged only if the Stage 2 distribution check agrees, D6); `derived` (MAP, P/F, S/F; D7, D8); `dropped`; `not_used` (with a reason).
  - The 16 contract channels must reproduce the PIC_Analysis_3_0 contract §11 ITEMIDs. That list is still not among the Foundations documents, so `contract16` is marked "likely" for 18 items until it is checked.
- **Drug-group map** (343 names): `DRUG_NAME` / `DRUG_NAME_EN` → **vasoactive** (one group, D1; `sensitivity_subgroup` vasopressor / inotrope for the sensitivity split only), **antibiotic**, **antifungal**, **antiviral** (D2), **sedative_analgesic**, **non_opioid_analgesic**, **antiemetic** (D4), **diuretic**, **osmotic_conditional** (D3: resolved per encounter to osmotherapy, diuretic or electrolyte_replacement by the admitting diagnosis), **electrolyte_replacement** (D18: potassium chloride, calcium and magnesium preparations, plus concentrated / 10% saline when the admitting diagnosis is not neuro), **steroid**, **NMB**, **AED** (D5), or `none`.
  - Vasoactive exposure is a flag plus the **count of concurrently active vasoactive agents**: PRESCRIPTIONS has no infusion rates (D1).
  - Cutoffs: 12 h vasoactive; 24 h antibiotics; 24 h sedatives/analgesics, diuretics and steroids (decision 3); NMB and AED 24 h (confirmed Nov 14, 2026); each with a 0.5× / 2× sensitivity.
  - **24 h, with the same 0.5× / 2× sensitivity (confirmed by Vinai Nov 14, 2026, D14)** for the new groups: antifungal, antiviral, non-opioid analgesic, antiemetic, the osmotic outcomes (osmotherapy, and the diuretic / electrolyte-replacement resolutions) and electrolyte_replacement (D18).
  - Promethazine and chlorpromazine are antiemetics only; there is no sedation rule for pethidine co-administration (D16).
- **Specimen map** (71 keys): the key is a test-panel code. **M1 uses `base_code` only**: 46 keys carry susceptibility add-ons (MIC, ESBL, MRSA and qualitative susceptibility tests) that are ordered only after a positive culture, so they would leak the result (D12). Groups: blood, urine, respiratory, CSF, stool, **ureaplasma**, **culture_site_unknown** (D12) and other.
- **Organism map and M3 turnarounds** (101 keys, keyed by name because ORG_ITEMID is reused across names): canonical name, Gram/category and contaminant rule. The blank-key row means "no organism recorded", not a confirmed negative. Turnarounds: gram class 24 h, identification/positivity 48 h, susceptibility 72 h, negative final 120 h; grid 0.5× / 1× / 1.5×; still bracketed in the plan.

### 4.3 Raw PIC tables (resolve_column with fallbacks, UTF-8 handling, as in v1.7)
PATIENTS, ADMISSIONS, ICUSTAYS, CHARTEVENTS, LABEVENTS, D_ITEMS, D_LABITEMS, PRESCRIPTIONS, INPUTEVENTS, OUTPUTEVENTS, MICROBIOLOGYEVENTS, SURGERY_INFO.
- Stage 1 **records** each raw file's SHA-256 and row count (decision 9).
- Stage 2 **pins and gates** them: the Stage 1 hashes are copied into the Stage 2 script, and any mismatch stops the run.
- The row counts must equal the schema-audit inventory where one exists: CHARTEVENTS [N_chartevents omitted]; LABEVENTS [N_labevents omitted]; PRESCRIPTIONS [N_prescriptions omitted]; MICROBIOLOGYEVENTS [N_microbiologyevents omitted]. Otherwise the run stops.
- **D36 (Nov 15, 2026).** Stage 1 did not hash PATIENTS, ADMISSIONS, ICUSTAYS or SURGERY_INFO; Stage 2 v1.0.1 recorded their SHA-256 and gated only their row counts ([N_patients omitted] / [N_admissions omitted] / [N_icustays omitted] / 7,488). From the follow-up checks (§5.4) on, these hashes are **gated inputs** (any mismatch stops the run), and the row-count gates stay:
  - PATIENTS `85b016fcc3050e5b5b6200103ed74ad1798720c881e26c60f1769be9f336ae0d`;
  - ADMISSIONS `1451ae8ea0ff463269b348522047a4bc99efd2665c770e9c5cef0ffe2647e1ae`;
  - ICUSTAYS `449119e698acaaf948921c34375db535b9fd3f6eb6e996a86069da95e90eea1b`;
  - SURGERY_INFO `04450bff4592d676a7ebe2dad22d9e6c2a6f8edd07dff223b182eb0002a38d83`.

## 5. Two-stage run (mirrors the surgery item-level audit)
### 5.1 Stage 1: item inventory (approved to run now, decision 2; script written Nov 14, 2026)
Script `scripts/feature_item_inventory/feature_item_inventory.py` v1.0 (SHA-256 `86388f7b628a4ae64d3cf59467213c591962ff94c7867e0133a37657f1d9aa41`; 61 synthetic-data tests pass), run by Copilot per its README.
- **Gates:** canonical, harmonization and split full SHA-256 (§4.1) and the duration audit v1.7 module, which is verified and then imported for the shared ID parsing, disclosure registry, suppression and recoverability solver (as in the φ1 band-support audit).
  - The cohort ([N_encounters omitted] encounters, [N_subjects omitted] subjects, the frozen 70/15/15 split (sizes omitted from this public copy)) and the core-landmark eligible counts ([landmark eligible counts omitted]) must reproduce, or the run stops (exit 2).
  - The raw row counts of CHARTEVENTS, LABEVENTS, PRESCRIPTIONS and MICROBIOLOGYEVENTS must equal the schema inventory (§4.3).
- **Tables:** CHARTEVENTS, LABEVENTS, OUTPUTEVENTS, INPUTEVENTS (keyed by ITEMID, joined to D_ITEMS or D_LABITEMS for label, category and fluid), PRESCRIPTIONS (keyed by `DRUG_NAME` + `DRUG_NAME_EN`), and MICROBIOLOGYEVENTS twice (specimen: `SPEC_ITEMID` + `SPEC_TYPE_DESC`; organism: `ORG_ITEMID` + `ORG_NAME`).
  - OUTPUTEVENTS and INPUTEVENTS are inventoried **for documentation only** (decision 4).
  - If a table has no ITEMID column, it is keyed by its unit column.
- **Linkage (H3):** `SUBJECT_ID` and `HADM_ID ∈ hadm_id_list`, then the event time inside [t0, OUT). Labs use the release time CHARTTIME + 1 h; prescriptions `STARTDATE`; the others `CHARTTIME`.
  - Each table reconciles: all rows = in stay (train / validation / test) + cohort admission outside the stay + no cohort admission + cohort admission with a missing time.
- **Per item, train split only:**
  - rows in stay;
  - rows with a numeric value (CHARTEVENTS and LABEVENTS `VALUENUM`; OUTPUTEVENTS `VALUE`; INPUTEVENTS `AMOUNT`; values themselves are never read out);
  - encounters with ≥ 1 row in stay, in the Layer B window, and in the pre-ICU window [t0 − 6 h, t0);
  - the unit strings seen (published only when the item's row count is itself published);
  - whether the item is in the dictionary.
- **Search hits:** dictionary labels and drug names are matched (English and Chinese keywords) for weight, height, urine output, fluid intake, GCS, pupils, FiO2, ventilation, vasoactive drugs, pain, occult blood, endocrine labs, sodium, NMB and AED. These are keys for Vinai to verify, not assignments.
- **Map templates** (blank decision columns, shared CSVs): channel ITEMID map, drug-group map (template groups vasopressor, inotrope, antibiotic, sedative_analgesic, diuretic, steroid, nmb, aed, none; superseded by the v1 groups, §4.2), specimen map, organism map.
  - **Only items whose train encounter count is ≥ 10** (published or complementary-suppressed) are listed. Items under 10 stay in the JSON as `<10`.
- **Disclosure:** §11 rules; items are listed in key order, never by count. Outputs `data/interim/feature_item_inventory.json` / `.md` plus the four templates; the unsuppressed `feature_item_inventory_raw_LOCAL.json` is never shared.
- Exit codes: 0 ok; 2 input, hash or reconciliation failure; 3 release gate.
- It also answers the outstanding item-audit searches: GCS, pupils, FiO2, ventilation, vasoactive dose, stool occult blood, endocrine labs, urine-sodium mis-mapping. Weight, urine-output and fluid-intake items in CHARTEVENTS are inventoried as well (decision 4).

### 5.2 Stage 2: sparsity and frequency
With the four v1 maps and the neuro set v1 hash-pinned (§4.2), and the raw-table hashes from Stage 1 pinned (decision 9), compute the §5.3 map checks first, then §6–§8 for every mapped channel and drug group.

### 5.3 Stage 2 map checks (decisions D1–D13, Nov 14, 2026)
All of these are train-split aggregates with the §11 suppression. Validation and test are shown beside train, never used to decide. No value-level or patient-level listing is shared.

**Outcome (Nov 15, 2026).** Stage 2 v1.0.1 ran these checks (release gate ok). Vinai's decisions on them are D21–D37 (§15), written into channel map v2 (§4.2). Two items stay open until the follow-up checks (§5.4): the bilirubin merge (D24) and P/F (D29).
1. **Duplicate lab items (D6).** Each merge group in the channel map (creatinine, urea, ALT, total bilirubin, CRP, total calcium, anion gap, base excess, pH, pCO2, pO2, glucose, bicarbonate, hemoglobin, ferritin, TSH):
   - per item, after unit conversion: p10/p25/p50/p75/p90 and the KS distance between the items;
   - merge only if they agree; otherwise keep the main item and drop the other;
   - `source_itemid` stays as a provenance column.
   - The calcium check also tests whether 5215 is ionized rather than total calcium.
   - **Unitless duplicates** (lipase 6473; stool occult blood 5341 and 6010): merged only if clearly on the same scale as the main item (value sets and quantiles), else dropped.
2. **Heart rate (D9).** 1002 Pulse and 1003 Heart Rate: per-item quantiles, plus agreement at the same CHARTTIME; 1003 is kept when both exist at the same time.
3. **MAP (D8).** MAP = (SBP + 2 × DBP) / 3 from **same-CHARTTIME** SBP/DBP pairs only. Report:
   - the share of SBP and DBP rows that have a pair;
   - the share of pairs dropped as implausible (DBP ≥ SBP; either value outside its gate).
4. **Oxygenation (D7). Superseded by D29 (Nov 15, 2026): P/F and the pSOFA respiratory component are dropped; the text below is kept for the v1.0.1 record.**
   - **P/F** = pO2 / (FiO2 5222 / 100), from the same draw (same encounter and CHARTTIME), with an oxygenation mask. It is scored on the **pSOFA P/F cut points** (≥ 400 = 0, 300–399 = 1, 200–299 = 2, 100–199 = 3, < 100 = 4), **without the respiratory-support requirement**, exactly like S/F (D19, Vinai, Nov 14, 2026).
   - Report the agreement with the analyser P/F 5243 at the same draw (quantiles of the difference).
   - Report P/F coverage per landmark.
   - **S/F (dropped, D28, Nov 15, 2026; text kept for the v1.0.1 record)** only where no P/F exists: SpO2 1006 within 5 min of FiO2 5222, and only when SpO2 ≤ 97%. It is scored on the **pSOFA S/F cut points** (Matics and Sanchez-Pinto 2017; confirmed by Vinai Nov 14, 2026, D15): ≥ 292 = 0, 264–291 = 1, 221–263 = 2, 148–220 = 3, < 148 = 4, **without the respiratory-support requirement**, and never converted to P/F.
   - Report how often S/F adds coverage where P/F is absent, per landmark and split.
   - Limitations (D15, D19), reported with both scores: PIC has **no ventilation item**, so the respiratory-support requirement of the top pSOFA points cannot be applied; and blood gases are **not marked arterial**, so P/F may include venous or capillary samples.
5. **Urine output and fluid input (D10, D17).** For CHARTEVENTS 1008 and 1009 (fluid output 1010 is dropped), report:
   - the distribution of gaps between entries;
   - the hour-of-day pattern of entries (per shift vs once a day).
   - Fluid intake 1009 is fixed as **mL/kg/day** (daily total / weight; charted weight, else weight imputed from age and gender, A2.12). It is the **surrogate for fluid status**; fluid balance is dropped because PIC does not hold enough information for it (D17, an acknowledged limitation).
   - Vinai then picks the urine-output aggregation: daily totals, or mL/kg/h over the charting interval (imputed weight, A2.12).
   - An empty bin or a day with no entry is never read as zero output or intake.
6. **Units (D13).**
   - Troponin T 6525 is converted per row to ng/mL (pg/mL ÷ 1000).
   - Blood-gas hemoglobin 5257 is converted per row to g/L for the 20–250 g/L gate (g/dL rows × 10); the **final feature unit is g/dL** (g/L ÷ 10) for both 5099 and 5257 (D27). The v1.0.1 script converts correctly to g/L but labels the step "→ g/dL"; it stays frozen, and the label is corrected in channel map v2 and in the Stage 3 builder.
   - Per-unit row counts are published only where ≥ 10.
   - Magnesium 5059 "ug/dl" rows: compare their distribution with the mmol/L rows, as recorded and after scaling. Convert them if one scaling agrees; otherwise drop them.
7. **Microbiology fields (D12).** List the MICROBIOLOGYEVENTS columns with their fill rates (aggregate), to look for a site field for culture_site_unknown keys. The `base_code` rule against susceptibility leakage stays in force.
8. **Osmotic agents (D3).** Each encounter's osmotic-agent courses are resolved by the Layer A **primary admitting diagnosis** (ADMISSIONS `ICD10_CODE_CN`, 3-character; an admission diagnosis, never a discharge code, so there is no leakage):
   - a neuro diagnosis makes mannitol, glycerol-fructose and concentrated / 10% sodium chloride **osmotherapy**;
   - otherwise mannitol and glycerol-fructose are **diuretic**, and concentrated saline is **electrolyte_replacement** (the D18 group);
   - sensitivity analysis: all of them are osmotherapy.
   - Report courses per resolved class.
   - **Neuro set:** `neuro_admitting_dx_icd10_v1.csv` (§4.2), confirmed by Vinai Nov 14, 2026. A diagnosis is neuro if its code matches a row with status `include` or `partial` under that row's `match_rule`.
     - Include (3-character prefix): A17, A80, A81, A83–A89, G00–G09 (CNS infection); S06, S07, T90 (head injury); I60–I64, I67, I69, G45, G46, P52 (stroke/bleed); G92, G93, P91, G36, G37 (encephalopathy).
     - Partial (kept): A39, S02, S09, R40 (whole 3-character category) and Q28 (Q28.2 / Q28.3 only, brain AVM).
     - Excluded: seizures (G40, G41, R56, P90), CNS tumours (C70–C72, D32, D33, D43; C75, D35, D44) and hydrocephalus / malformations (G91, Q01, Q03, Q04, Q05, Q07), because they "don't typically require osmotherapy"; the neuromuscular / chronic codes G61, G70, G12, G71 and G80 stay excluded.
9. **Vasoactive count (D1).** The distribution of the number of concurrently active vasoactive agents, at each core landmark and in the Layer B window. The vasopressor / inotrope subgroup split is reported as a sensitivity analysis only.
10. **Serology (D11).** Serology items (IgM, CMV / EBV / HSV antibodies, Mycoplasma antibodies) are not time-varying channels and are not reported in §7.

### 5.4 Stage 2 follow-up checks (D24, D29; Nov 15, 2026)
Script `stage2_followup_checks.py` v1.0 (SHA-256 `e955974d8e9fd85414482549cc14e019d922215505088220ae9cfba990992f62`; 9 synthetic-data tests pass). It is read-only and aggregate-only, and it is run by Copilot unedited.
- **Gates:** the Stage 2 script v1.0.1 (imported, not run) and through it the v1.7 module; the §4.1 pins with the cohort, split and eligibility totals; channel map v2 (§4.2); LABEVENTS (Stage 1 hash) and the four D36 hashes (§4.3).
- **Scope:** train split; analysis window [t0, min(OUT, t0 + 336 h)); labs at release time (CHARTTIME + 1 h).
- **(i) Bilirubin 5255 vs 5075 (D24).** Groups: all train; NICU vs non-NICU (first care unit, as in Stage 2); neonate vs not-neonate (age < 28 days at t0; `age_invalid` in not-neonate).
  - Per item per group: n values, n encounters, value quantiles and per-encounter-median quantiles.
  - Per group: the KS distance 5255 vs 5075 (values and encounter medians; only when both n ≥ 10), encounters with either / both items, same-draw pairs and the same-draw difference.
  - Decision rule for Vinai: merge if the items agree within NICU (and within non-NICU); otherwise keep 5075 and treat 5255 as a separate or dropped item.
- **(ii) FiO2 5222 format (D29).**
  - FiO2 rows by band: ≤ 1, (1, 21), [21, 100], > 100, and nonnumeric.
  - Same-draw pO2 / FiO2 pairs: P/F = pO2 / (FiO2 / 100) in gate, below 10, above 700 (FiO2 ≤ 0 apart); out-of-gate pairs that fall in the gate with FiO2 × 100 or ÷ 100; P/F under the rule "FiO2 ≤ 1 is a proportion, otherwise a percentage".
  - pO2 and FiO2 quantiles for in-gate vs out-of-gate pairs; pairs and out-of-gate pairs by FiO2 band.
  - The Stage 2 pair definition (FiO2 gated 21–100) is recomputed and compared with the v1.0.1 valid and outside-gate pair counts (counts not quoted here; V29).
- **Disclosure and outputs:** §11 rules through the Stage 2 machinery; `data/interim/stage2_followup_checks.json` / `.md` only. No unsuppressed file is written.
- **Result (Nov 15, 2026; run with script `a5619ebb…` on map v2 `9872e688…`; release gate ok; output JSON SHA-256 `9280ef6f…`; aggregate, train).**
  - **Bilirubin:** NICU KS 0.08 on values / 0.083 on encounter medians (5075 p10 / p50 / p90 22.1 / 83.7 / 193.2 µmol/L; 5255 18.0 / 98.0 / 204.0); non-NICU KS 0.405 / 0.298 (median 10.0 vs 40.0); neonates KS 0.095. Same-draw pairs <10 in every group. In NICU 5255 has 12,229 values against 3,034 for 5075, and 473 NICU encounters have 5255 only; outside NICU 98 encounters have 5255 only.
  - **FiO2:** the numeric train rows lie in 21–100 (percent) apart from bands that are each <10. The Stage 2 pair definition reproduces exactly. With FiO2 21–100 and pO2 ≥ 10, P/F cannot fall below 10, so every out-of-gate pair is above 700 (reasoning); about half of the in-gate pairs have P/F ≥ 400 (Stage 2). Counts are not quoted (V29: they would expose suppressed cells by subtraction). The analyser P/F 5243 agrees with the derived P/F (p10–p90 −1.0 to 1.0), so the cause is the FiO2 entered at the analyser (likely left at 21%), not the arithmetic or the unit.
  - The pair tables of part (ii) are complementary-suppressed (they share identities with a few non-21–100 pairs); the Stage 2 definition carries the result.

## 6. Layer A (static): missingness, prevalence and level support
Per feature (plan v12 §4, contract v12 Layer A table, amendment A2.1), per split:

| Feature | Report |
|---|---|
| sex | share empty or unrecognized (no `sex_missing` column exists in v12; a nonzero count goes to Vinai) |
| age at INTIME, age_band | invalid DOB or age (missing, including no PATIENTS row; unparseable; or negative, t0 < DOB) → the `age_invalid` band, with `age_days` unset; aggregate breakdown in `layer_a.dob_breakdown` (status, negative magnitude bins <1 day / 1–30 days / >30 days, by split and first_careunit); age-continuous quantiles use valid ages only; a safety gate stops the run if invalid ages exceed 1% of the cohort (D20). Band counts include `age_invalid`. |
| hours_admit_to_icu_intime | missing; negative values (count); the four interval bands |
| admission timing | hour, weekday, the three shifts, weekday/weekend and the four φ2 categories (never missing; partition check) |
| first_careunit | level counts (never missing) |
| diagnosis (ICD10_CODE_CN → three-character ICD-10 → chapter → domain) | share missing; share unmapped or ambiguous; the domain map's coverage |
| surgery, Layer A recent (T in [t0 − 24 h, t0 + 6 h]) | prevalence; type and system support; φ1 band (3-level ordinal) counts from the frozen medians |
| surgery, remote (same HADM, T < t0 − 24 h) | prevalence; before-ADMITTIME flag count |
| admission source (inferred) | share per first-match rule; conflicting-evidence count; **EDREGTIME filler share** (equal to ADMITTIME or another placeholder; H8 [AUDIT]) |
| prior_encounters (candidate) | distribution (0, 1, 2, ≥ 3) |
| pre-ICU lookback block | encounters with ≥ 1 event in [t0 − 6 h, t0) and [t0 − 24 h, t0) per table and per mapped channel (labs at release time); medication groups active at t0 |
| **pre-ICU lab block (approved, decision 11)** | per lab channel: the share of encounters with a last released value (CHARTTIME + 1 h) in [t0 − 6 h, t0) (primary) and in [t0 − 24 h, t0) (sensitivity; half-open at t0, confirmed Nov 14, 2026), i.e. the value-present share for the block's value + mask columns |
| weight | share of encounters with a charted weight item (Stage 1 inventory). **Weight is treated as missing and is imputed for now from age and gender**; the future approach is a regression on age, gender, unit and primary diagnosis (decision 4) |

**Stage 2 review decisions for Layer A (Nov 15, 2026):**
- **Age (D34, D37).** `age_invalid` stays in the primary analysis. Two sensitivity analyses: excluding `age_invalid` encounters, and treating negative ages under 1 day as age 0. Descriptive age percentiles are reported from valid-age encounters only, with n rounded to the nearest 10.
- **first_careunit (D35).** Kept in Layer A. Chart-channel masks largely reflect the care unit (the general ICU has almost no charted vitals), which is stated as a limitation; a unit-stratified sensitivity analysis is added.

**Categorical level support** (diagnosis category, chapter, domain, surgery type, first_careunit, admission source):
- the number of levels;
- levels with ≥ m train **subjects**, at m = 20 (the plan's pooling threshold) and at m = 10;
- the share of train rows pooled to `other` (within chapter × domain for diagnosis);
- the share of validation/test rows whose level is unseen in train;
- per-level counts suppressed under 10.

## 7. Layers B and C (time-varying): metrics
Sources and units of measurement:

| Source | Unit | Notes |
|---|---|---|
| Vitals (CHARTEVENTS; monitor and charted) | one charted value | HR (1003 and 1002 merged, D9), RR, SBP, DBP, MAP (derived from same-time SBP/DBP pairs, D8), temperature, SpO2, bedside glucose (1011); v12 additions pain score, urine output, fluid intake (mL/kg/day, the fluid-status surrogate; no fluid balance, D17). **Fluid intake (and urine output) come from CHARTEVENTS** (decision 4: CHARTEVENTS is the best source; OUTPUTEVENTS is not well documented and INPUTEVENTS is not complete). Values outside the plausibility gate are counted separately and are not measurements. |
| Labs (LABEVENTS) | one released result (release = CHARTTIME + 1 h) | the 16-channel contract labs plus the Table 5.1 extension; triggered tests flagged |
| Medications (PRESCRIPTIONS) | one documentation event, then the S3 course | per group of the v1 drug map (§4.2; D1–D5): vasoactive (exposure flag + count of concurrent agents), antibiotic, antifungal, antiviral, sedative_analgesic, non_opioid_analgesic, antiemetic, diuretic, osmotic (resolved by admitting diagnosis), steroid, NMB, AED |
| Microbiology (MICROBIOLOGYEVENTS) | one collection (M1) | per specimen group; M3 result availability per result type and turnaround grid |
| Surgery (SURGERY_INFO, Layer C event) | one surgery with t0 + 6 h < T ≤ OUT | by system (map v1); **placed in the bin of T** (anesthesia end, A2.1; decision 8, Nov 14, 2026), with the begin-time bin as a sensitivity analysis. This supersedes the "start bin" of amendment decision 12 |

For each channel (and each drug group, specimen group, and surgery), per split:

1. **Ever measured:** the share of encounters with ≥ 1 measurement in the Layer B window, and in [t0, min(OUT, t0 + 336 h)).
2. **Landmark availability:** at each core landmark L, among the encounters eligible at L, the share with ≥ 1 value available before t0 + L:
   - (a) in the full prefix [t0, t0 + L) (**primary for §9.1**, decision 5);
   - (b) in the last 24 h, [t0 + L − 24, t0 + L) (**sensitivity**, decision 5);
   - (c) in the fixed 72 h window before L (the A1.6 sensitivity);
   - (d) in each prefix bin k < L/6.
   - Labs are also reported at + 0 h and + 2 h release.
3. **Measurements per patient-day:** the count in [t0, min(OUT, t0 + 336 h)) divided by the open window in days, as median and IQR, among measured encounters and among all encounters. Encounters with an open window under 6 h are excluded from this metric and counted. This is a descriptor only; per S2, counts and rates are never features.
4. **Inter-measurement interval (h):** between consecutive measurements of the same channel in the same encounter.
   - Reported as the pooled p10/p25/p50/p75/p90, and as the distribution of the per-encounter median interval.
   - Same-timestamp duplicates are counted separately (zero intervals).
   - Intervals that cross a merge gap are flagged.
5. **Time from t0 to first measurement (h):** quantiles among the encounters measured within [t0, OUT). The share whose first measurement falls in the pre-ICU lookback is reported separately.
6. **Last-observation age at L:** among the encounters eligible at L with ≥ 1 value before t0 + L, the quantiles of t0 + L − (time of the last available value).
   - Also reported: the share whose age is within the contract LOCF cap (1 bin / 6 h, or 4 bins / 24 h, by family; triggered tests and vitals have no carry).
7. **Empty bins:** the share of open bins with no measurement.
   - Reported over all open bins in [t0, OUT) (to bin 55), and over each core landmark's prefix (bins whose end is ≤ L).
   - Plus the per-bin occupancy curve: for k = 0..55, among encounters with OUT ≥ t0 + 6(k + 1), the share with ≥ 1 value in bin k.
   - Bins inside merge gaps (101 encounters) are **included as primary and excluded as a sensitivity analysis** (decision 12).
8. **Within-bin count:** among occupied open bins, the share with exactly 1, 2, or ≥ 3 values. This tells whether a bin median, min and max are informative.
9. **Co-measurement and panel structure (aggregated):**
   - Labs: a draw is the set of channels sharing (encounter, CHARTTIME). For each channel pair, report P(B drawn | A drawn) and the Jaccard index over train draws, with any cell < 10 suppressed.
   - The top panels (frequent channel sets drawn together, support ≥ 100 train draws), as channel lists with counts.
   - Vitals: same-CHARTTIME co-recording of HR / RR / BP / SpO2 / temperature.
   - Bin-level co-occupancy of channel pairs, at k = 0..11.
   - All of these are aggregates. No draw-level listing is shared.
10. **Medications (the v1 drug groups, §4.2):**
    - documentation events per patient-day and inter-documentation intervals per group (v1 groups), including the share of intervals longer than the group cutoff (12 h vasoactives; 24 h antibiotics, sedatives/analgesics, diuretics and steroids; NMB and AED 24 h, confirmed Nov 14, 2026; the new v1 groups and electrolyte_replacement 24 h, confirmed Nov 14, 2026, D14/D18, §4.2), and at 0.5× and 2× the cutoff. This share is the rate at which the state machine splits courses.
    - `ENDDATE` per group: missing; present and after start; present and ≤ start; date-only `STARTDATE` share.
    - Exposure duration per course under the S3 rule, split into cutoff-implied ends and ENDDATE ends.
    - Among courses with a valid ENDDATE, the difference between the cutoff-implied stop and ENDDATE (h): quantiles, plus the shares where the cutoff ends the course earlier, later, or within ± 6 h.
    - The share of active course-hours at each core landmark that are cutoff-implied.
    - Active prevalence at each core landmark per group.
11. **Microbiology:** M1 collection metrics 1–7 per specimen group, plus the share of pre-ICU collections.
    - M3: the share of collections whose result is available before t0 + L at each core landmark, per result type × turnaround grid.
    - The share of encounters whose M3 features differ from M1 at 24 h (Layer B) and at each landmark.
    - M2 is look-ahead and is not reported.
12. **Surgery (Layer C event):**
    - the share of encounters with ≥ 1 Layer C surgery, and per landmark prefix;
    - surgeries per encounter;
    - time from t0 to the first Layer C surgery;
    - the bins with a surgery;
    - by system;
    - in_gap descriptor count;
    - the reuse count (Layer C for one encounter and remote Layer A for a later one, A2.1).

**Strata** (decision 7; train, validation and test each):
- age_band (6 levels);
- surgical status = Layer A recent surgery (T in [t0 − 24 h, t0 + 6 h]) yes/no;
- `first_careunit`.

Strata are reported only for metrics 1, 2, 3, 7 and 10. A stratum's cell is suppressed under 10 like any other.

## 8. Informative-missingness descriptors (**approved, decision 1, Nov 14, 2026**)
Computed on **train only**:
- metrics 1, 2(a) and 3 by **Outcome 1 group** at each core landmark (levels 0–2 / 3–5 / 6–7 / 8);
- and by **in-hospital mortality**.
- They are **descriptive only** and are **never used** for the §9 disposition rules (using them would select features by outcome).

The Stage 2 script keeps this section behind the explicit `--outcome-descriptors` flag; the flag's state is recorded in the output. The Copilot command for the approved run includes the flag.

## 9. Decision criteria (approved Nov 14, 2026 as starting defaults, decision 10; applied by Vinai on train only after he reviews the per-channel results)
### 9.1 Feature disposition (per channel; Layer B at the 24 h window; Layer C at the core landmarks 24, 72 and 168 h; full prefix primary, last 24 h sensitivity)
| Class | Rule (train; approved default) | Representation |
|---|---|---|
| **Drop** | ever measured in [t0, min(OUT, t0 + 336 h)) in < 1% of train encounters, **or** < 100 train encounters | none (logged) |
| **Measured indicator only** | ever measured 1–10%, **or** (Layer C) median last-observation age at L above 72 h among measured encounters | the ordered/measured flag plus mask (Layer B: measured in window; Layer C: measured in bin, plus ever measured before L); no value column |
| **Value + mask + time since last** | ever measured ≥ 10% **and** (Layer C) ≥ 20% of open prefix bins occupied among measured encounters at L = 72 h | the value per the layer's level rules, with `*_missing`, `*_carried_forward` (C), and time since last measurement (C) |

Overrides and constraints:
- **Triggered tests** (Table 5.1 triggered column) keep their S9 form (ordered flag, worst released value, mask) whenever they are not dropped. The indicator-only class does not apply to them.
- **Score inputs:** components of modified pSOFA / PRISM III are kept as value + mask if not dropped, even in the indicator class, and flagged. pSOFA is a **5-component score without the respiratory component** (D29): cardiovascular, coagulation, hepatic, renal and neurological; the hepatic component carries the neonatal note of D24.
- **Time since last measurement (decision 6):** **included** for value-class Layer C channels, **capped at the channel's LOCF cap** (24 h for the 4-bin family, 6 h for the 1-bin family; **24 h for channels without an LOCF cap**, confirmed Nov 14, 2026). A sensitivity analysis drops it. It is known at t, so it is not leakage.
- **Counts:** measurement counts and rates stay audit-only (S2). The Layer C count column remains its own open decision; this audit informs it.
- **Medication groups and NMB/AED:** a group with train prevalence below 1% of encounters is dropped above φ1 (the plan's "[sparsity cutoff]", set at 1% as an approved starting default).

### 9.2 Layer C bin-width adequacy (6 h bins are approved; the audit tests fit, not the geometry)
- **Vitals:** 6 h is adequate if, for HR, RR, MAP, SpO2 and temperature, ≥ 50% of occupied open train bins have ≥ 3 values (bin median, min and max are meaningful) and the median inter-measurement interval is ≤ 2 h.
  - If a vital fails, its per-bin summary is reduced to median plus mask (no min/max/range).
- **Labs:** a channel is "sparse at 6 h" if, among measured encounters, ≥ 80% of open prefix bins are empty at L = 72 h.
  - Sparse channels keep 6 h bins with their LOCF cap plus time since last.
  - A 12 h or 24 h aggregation is run only as a sensitivity analysis, and only if ≥ 5 Table 5.1 channels are sparse.
- **LOCF caps:** report the share of last-observation ages beyond the cap at each landmark. If > 50% of measured encounters exceed the cap, flag the channel for Vinai. The contract caps stay unchanged unless he decides otherwise.
- **Long horizon:** if fewer than 100 train encounters remain open in a bin (≈ late bins at 336 h), report that bin range as thin.

### 9.3 Imputation per sparsity class (v11 arms; masks are kept in every arm)
Classes are defined at L = 72 h on train by the share of empty open prefix bins among measured encounters:

| Class | Rule (approved default) | Imputation |
|---|---|---|
| Dense | empty ≤ 20% | contract LOCF; then every arm, (r), (a), (b) and (c), is compared as planned |
| Moderate | 20–80% | contract LOCF; arm (a) (age-specific normal fill) and arm (r) are primary candidates; arms (b)/(c) are compared only if the channel's observed fraction in the train prefix matrix is ≥ 20% |
| Sparse | > 80% | no value imputation beyond LOCF; the cell gets the arm (r) constant only for matrix algebra; mask plus time since last carry the information; excluded from low-rank completion, like triggered tests |
| Indicator only / triggered | §9.1 | never imputed as if measured |

Everything is fit on train, time-respecting at each landmark (prefix-only for (b)/(c)). Masked-known-value RMSE stays secondary, because informative missingness breaks its assumption. This follows the Ali, Jiang et al. comparison (AMIA 2026, session S36) as adopted in plan v11/v12.

### 9.4 Medication ENDDATE missingness
- **Primary stays the S3 state machine:** cutoff-implied stops, with an ENDDATE at or before t ending a course early (the assumption confirmed Oct 29, 2026).
- **ENDDATE presence is never a feature** (a future ENDDATE is not known at t).
- **Flag a group** (approved default) if, among its courses with a valid ENDDATE, the cutoff-implied stop and ENDDATE differ by > 24 h in > 25% of courses, **or** if ENDDATE presence differs by > 10 percentage points between splits.
  - For flagged groups, `hours_on_course` and `hrs_active` are treated as low-reliability: kept as φ0/φ1 columns, but primary interpretation rests on `active`, `started` and `ever_received`.
  - A sensitivity analysis rescales the group cutoff (0.5× and 2×).
- The methods text must state that most active flags are cutoff-implied (78.07% missing ENDDATE).

## 10. Outputs
- Shared: `data/interim/feature_sparsity_frequency_audit.json` and `.md` (aggregate only). Stage 1: `feature_item_inventory.json` / `.md`.
- Local only: `*_raw_LOCAL.*`, unsuppressed, never shared.
- Console: a short summary. It lists hashes ok, eligible counts reproduced, the number of channels per class, suppression counts and the release-gate status.

## 11. Disclosure (same as the earlier audits)
- Counts 1–9 are published as `<10`. Complementary `S` suppression is added until the LP recoverability check finds nothing recoverable, across splits, strata and declared identities (e.g. stratum sums = split totals; per-bin occupancy ≤ eligible).
- Percentages and quantiles are published only where the underlying count is ≥ 10 and published.
- Min/max, IDs, timestamps and draw-level or patient-level rows are never published. Quantiles are limited to p10–p90.
- Shared files are written only after the disk-reload release gate passes.

## 12. Acceptance checks (the script stops on failure)
1. Any hash mismatch: pinned inputs, the v1 maps (§4.2) or, from the follow-up checks on, channel map v2, the neuro admitting-diagnosis set v1 (exact SHA-256 `3cc8f4acfe9df29037b60d903ad4204c40db5bc76923de8b64234ccb2553b9d9`), or a gated raw table (including the four D36 tables from the follow-up checks on).
2. The cohort is not [N_encounters omitted], or the split totals or the core-landmark eligible counts do not reproduce exactly.
3. Event-linking reconciliation per table: assigned + out of window + other admission + unlinked = rows.
4. The lab release offset is not applied: a unit check on synthetic rows at the bin and landmark boundaries (CHARTTIME + 1 h = t0 + L is not available at L; an observation at exactly t0 + 24 h is in bin 4).
5. The S3 state machine fails its unit tests, or the as-of rebuild at 48 h differs from the full build for bins 0–7.
6. Any metric computed from validation or test feeds a train threshold.
7. A release-gate failure.
8. Invalid ages exceed 1% of the cohort (137 encounters on [N_encounters omitted]) (D20 safety gate). Missing or invalid DOB below that threshold is not a stop.

## 13. Tests (synthetic data only; written with the script)
- Window and boundary inclusivity (bin 4 at t0 + 24 h; strict landmark eligibility; half-open bins).
- Lab release at + 0 / + 1 / + 2 h.
- The pre-ICU lookback at 6 h and 24 h, capped by admission ID.
- Gap bins.
- Inter-measurement intervals with duplicates.
- Last-observation age vs LOCF caps.
- Co-measurement on planted panels.
- The state machine with and without ENDDATE (including ENDDATE ≤ start and date-only starts).
- M3 availability over the grid.
- Surgery classes per A2.1.
- Suppression, the recoverability check and release-gate tamper cases.
- Hash-gate failures.
- `--outcome-descriptors` off by default.
- Map rules (D1–D13), on planted synthetic rows:
  - HR merge: 1003 is kept over 1002 at the same CHARTTIME;
  - MAP: same-CHARTTIME pairs only; DBP ≥ SBP and out-of-range pairs dropped;
  - P/F: same draw only, scored on the pSOFA P/F cut points without the respiratory-support requirement (D19);
  - S/F: within 5 min only, SpO2 ≤ 97% only, and only where no P/F exists, scored on the pSOFA S/F cut points without the respiratory-support requirement (D15);
  - per-row unit conversion (troponin T, blood-gas hemoglobin);
  - merge groups keep `source_itemid`;
  - osmotic resolution by neuro vs non-neuro admitting diagnosis, plus the all-osmotherapy sensitivity;
  - the concurrent vasoactive-agent count;
  - susceptibility add-ons never reach M1 (`base_code` only);
  - an empty urine-output bin is never zero;
  - D20: missing / unparseable / negative DOB, the `age_invalid` band, weight sex-only / overall fallback, the 1% safety gate, and small-cell suppression of `layer_a.dob_breakdown`.

## 14. Copilot run rules (same as before)
Copy the script into the project `scripts/` folder, verify its SHA-256, run it unedited, and stop on any failure without patching. Report only the shared JSON/MD and the console summary.

## 15. Decisions (Vinai, Nov 14, 2026; answers to the v0.1 open questions)
1. **Outcome-stratified descriptors:** yes, on train only, by Outcome 1 group (0–2 / 3–5 / 6–7 / 8) and in-hospital mortality; descriptive only and never used for §9; still behind `--outcome-descriptors`, whose state is recorded in the output (§8).
2. **Stage 1 item inventory:** run now (§5.1).
3. **NMB and AED** are separate drug groups. Sedatives/analgesics, diuretics and steroids use a 24 h cutoff, with 0.5× / 2× sensitivity.
4. **Fluid intake:** the best source is CHARTEVENTS; OUTPUTEVENTS is not well documented; INPUTEVENTS is not complete. **Weight** is missing and will be imputed for now based on age and gender; the best approach (future) is a regression using age, gender, unit and primary diagnosis. Stage 1 still inventories urine-output, intake and weight items in CHARTEVENTS, and reports OUTPUTEVENTS / INPUTEVENTS row and item counts, aggregate-only, for documentation.
5. **Availability window for §9.1:** full prefix primary; last 24 h before the landmark as a sensitivity analysis.
6. **Time since last measurement:** included for value-class Layer C channels, capped at the channel's LOCF cap; sensitivity analysis without it.
7. **Strata:** age_band, surgical status from Layer A recent surgery, and `first_careunit`.
8. **Layer C surgery event bin** = the bin of T (anesthesia end, A2.1); the begin-time bin is a sensitivity analysis. This supersedes the "start bin" of amendment decision 12.
9. **Raw-table hashes:** recorded in Stage 1; pinned and gated in Stage 2.
10. **§9 thresholds:** all approved as starting defaults; Vinai reviews the per-channel results before locking.
11. **Layer A pre-ICU lab block: approved** (H3.5): the last released value (CHARTTIME + 1 h) in [t0 − 6 h, t0) plus a mask (half-open at t0, confirmed Nov 14, 2026); 24 h lookback as a sensitivity analysis.
12. **Gap bins:** included in the empty-bin share as primary; excluded as a sensitivity analysis.

### Follow-up confirmations (Vinai, Nov 14, 2026)
The points left open by the decisions above are now settled:
1. **NMB and AED cutoffs:** 24 h, with the same 0.5× / 2× sensitivity (§4.2).
2. **Pre-ICU lab window:** [t0 − 6 h, t0), half-open at t0 because Layer B starts at t0 (invariant I6); the 24 h sensitivity is [t0 − 24 h, t0). The Stage 1 inventory already uses [t0 − 6 h, t0).
3. **Time-since-last cap for channels without an LOCF cap** (vitals such as HR, MAP, RR, temperature, SpO2; SBP, DBP, pain, urine output, fluid intake; fluid balance was later dropped, D17): **24 h** (not the 6 h proposed in v0.2). Channels with an LOCF cap keep their own cap.
4. **Urine output source:** CHARTEVENTS, as for fluid intake; OUTPUTEVENTS is documented only.
5. **Canonical-encounters hash gate** (`e5feeb22…`) in Stage 1: approved.

The amendment's A2.7 item "whether the expected duration appears in the Layer C start bin" is **closed as superseded** by decision 8 (surgery event in the bin of T; amendment A2.10).

### Map decisions (Vinai, Nov 14, 2026; D1–D19, recorded in the v1 maps and amendment A2.14)
1. **D1 Vasoactives:** PRESCRIPTIONS has no infusion rates, so dopamine, adrenaline, norepinephrine, phenylephrine, dobutamine, milrinone and isoprenaline form **one vasoactive group**: an exposure flag plus the count of concurrently active agents. The vasopressor / inotrope split is a sensitivity analysis only.
2. **D2:** antifungals and antivirals are **two separate groups**. They are not antibiotics and not "other".
3. **D3 Osmotic agents:** classified by the Layer A **primary admitting diagnosis** (ADMISSIONS, a true admission diagnosis, not a discharge code; no leakage).
   - With a neuro diagnosis: mannitol, glycerol-fructose and concentrated / 10% sodium chloride count as osmotherapy.
   - Otherwise: mannitol and glycerol-fructose count as diuretics, and concentrated saline as electrolyte replacement.
   - Sensitivity analysis: all of them are osmotherapy.
   - The neuro ICD-10 set v1 (3-character categories) was **confirmed by Vinai on Nov 14, 2026**: CNS infection, head injury, stroke/bleed and encephalopathy, with their partial codes, are neuro. Seizures, CNS tumours and hydrocephalus / malformations are excluded because they "don't typically require osmotherapy". Pinned in §4.2.
4. **D4:** promethazine and chlorpromazine form an **antiemetic** group, together with ondansetron, metoclopramide and the other antiemetics present. Non-opioid analgesics (paracetamol, ibuprofen, metamizole, etc.) form their own group, outside the sedatives.
5. **D5:** clonazepam and nitrazepam → AED; diazepam and midazolam stay sedatives.
6. **D6 Duplicate lab items:**
   - merge into one channel when the units match or convert **and** the distributions agree (Stage 2 check, §5.3);
   - unitless duplicates merge only if clearly on the same scale, else they are dropped;
   - the source ITEMID is kept as a provenance column.
7. **D7 Respiratory score:**
   - P/F from PaO2 plus blood-gas FiO2 5222 at the same draw, with lab P/F 5243 as a cross-check and an oxygenation mask;
   - without a P/F: S/F when SpO2 and FiO2 are within 5 min, only if SpO2 ≤ 97%, mapped to its own cut points, never converted to P/F;
   - Stage 2 counts how often S/F adds coverage.
8. **D8:** MAP = (SBP + 2 × DBP) / 3, from same-charting-time pairs only, dropping implausible pairs.
9. **D9:** Pulse 1002 and HR 1003 merge into one heart-rate channel, preferring 1003 at the same time; the source is kept as provenance, with a Stage 2 distribution check.
10. **D10:**
    - fluid output 1010 is **dropped** (documentation inconsistent);
    - only fluid input 1009 and urine output 1008 (CHARTEVENTS) are used;
    - the aggregation is chosen after the Stage 2 charting-interval check, and empty bins are never read as zero.
11. **D11:** serology (IgM; CMV / EBV / HSV antibodies, etc.) is dropped from the time-varying channels.
12. **D12:**
    - site-less cultures and smears (LIS0156 variants, LIS05088, LIS0146, LIS0144, LIS0258) → `culture_site_unknown`;
    - Ureaplasma LIS031292 → its own category;
    - Stage 2 lists the MICROBIOLOGYEVENTS columns and their fill rates;
    - the `base_code` rule against susceptibility leakage stays.
13. **D13 Units:**
    - troponin T 6525 is converted per row to ng/mL;
    - blood-gas Hb 5257 is converted per row to g/dL;
    - magnesium "ug/dl" rows get a Stage 2 distribution check: converted if consistent after scaling, else dropped.
14. **D14 Cutoffs for the new groups (confirmed):** antifungal, antiviral, antiemetic, non-opioid analgesic and the osmotic outcomes use **24 h**, with the 0.5× / 2× sensitivity.
15. **D15 S/F cut points:** the pSOFA S/F cut points (≥ 292 = 0, 264–291 = 1, 221–263 = 2, 148–220 = 3, < 148 = 4), scored **without the respiratory-support requirement**. Documented limitations: PIC has no ventilation item, and blood gases are not marked arterial.
16. **D16:** chlorpromazine and promethazine are **antiemetics only**; there is no sedation rule for pethidine co-administration.
17. **D17 Fluid status:** fluid balance is **dropped** (not enough information in PIC). Fluid intake in **mL/kg/day** (imputed weight where weight is missing, A2.12) is the surrogate for fluid status, an acknowledged limitation.
18. **D18 Electrolyte replacement:** a new **electrolyte_replacement** drug group: potassium chloride, calcium and magnesium preparations (moved out of `none`), plus concentrated / 10% saline when the admitting diagnosis is not neuro (D3). Cutoff 24 h, with the 0.5× / 2× sensitivity. Calcium acetate (a phosphate binder) and oral magnesium sulfate powder (a laxative / topical compress) are not replacement preparations and stay in `none` (accepted by Vinai, Nov 14, 2026).
19. **D19 P/F scoring (superseded by D29, Nov 15, 2026: P/F dropped):** P/F is scored on the pSOFA P/F cut points (≥ 400 = 0, 300–399 = 1, 200–299 = 2, 100–199 = 3, < 100 = 4) **without the respiratory-support requirement**, the same as S/F (D15). The same two limitations apply: no ventilation item, and blood gases not marked arterial.
20. **D20 Age-invalid band and DOB breakdown (Nov 15, 2026):** after the v1.0 Stage 2 run stopped on 42 cohort encounters with missing or invalid DOB, continue with an `age_invalid` band and report the aggregate DOB breakdown in v1.0.1. Status is missing (null/blank DOB, or no PATIENTS row), unparseable, negative (t0 < DOB), or valid. Negative ages are binned by magnitude (<1 day, 1–30 days, >30 days). Age-continuous quantiles use valid ages only. Weight is not imputed from the `age_invalid` band; those encounters use the train sex-only median (n ≥ 10) or, if that is unavailable, the overall train median. Shared outputs are aggregate only, with the usual <10 / complementary-S suppression. A safety gate still stops the run if invalid ages exceed 1% of the cohort.

### Stage 2 review decisions (Vinai, Nov 15, 2026; D21–D37, recorded in channel map v2 and amendment A2.16)
Taken after reviewing the Stage 2 v1.0.1 outputs (aggregate, release gate ok). D24 and D29 wait on the follow-up checks (§5.4).
21. **D21 Ionized calcium:** 5215 is ionized calcium (87,025 of 90,852 train values in the ionized range; KS 0.968 against 5034). It becomes a new channel `calcium_ionized` in mmol/L, with no conversion, plausibility gate 0.3–3.0 mmol/L (confirmed Nov 15, 2026) and its own mask. `calcium_total` is 5034 only (mmol/L → mg/dL × 4.008).
22. **D22 Merges:** base excess 5211 + 5249 (KS 0.011), bicarbonate 5248 + 5224 (KS 0.061), BUN 5033 + 6477 (KS 0.101), CRP 5626 + 5821 (KS 0.034), hemoglobin 5099 + 5257 (KS 0.114; a per-value source flag, CBC or blood gas, is kept as a feature) and stool occult blood 5326 + 5341 + 6010 (qualitative) are merged. pCO2 5235/5236, pH 5237/5238 and pO2 5239/5244 are deduplicated by draw (one value per encounter and CHARTTIME, main item kept; their same-draw difference is 0 at p10–p90).
23. **D23 Main item only:** creatinine 5041, glucose 5223 and ferritin 6459 are kept; 5032 (KS 0.411), 5047 (KS 0.604) and 6161 (KS 0.561) are dropped.
24. **D24 Bilirubin (final, Nov 15, 2026):** one `bilirubin_total` channel with 5075 as the main item. 5255 is merged, with a per-value source flag, **only when first_careunit = NICU**, and dropped elsewhere (98 non-NICU encounters lose bilirubin). Evidence (§5.4 result): the items agree in NICU (KS 0.08 / 0.083) and in neonates (0.095), disagree outside NICU (0.405 / 0.298), are never on the same draw, and 5255 is the only bilirubin for 473 NICU encounters. Sensitivity analyses: 5075 only; merged everywhere with the flag; the neonate (< 28 days) rule instead of NICU. **This departs from the §5.4 rule as written** (agreement required in both groups); approved by the PI. **Neonatal note for pSOFA hepatic scoring:** the neonatal median total bilirubin is about 90–100 µmol/L (≈ 5–6 mg/dL, about 2 hepatic points) and largely reflects physiological jaundice; a neonatal sensitivity analysis (neonates excluded from, or scored separately in, the hepatic component) is reported.
25. **D25 Anion gap dropped:** it is a linear combination of sodium, chloride and bicarbonate, so it adds no information to the matrix (5212/5213 also disagree, KS 0.472).
26. **D26 ALT and TSH; lipase 6473:** ALT 5026 + 5195 (KS 0.427) and TSH 6193 + 6095 (KS 0.168) are merged **for the measured indicator only**; their values are not used as features. Lipase 6473 is dropped as a likely mislabelled serum lipemic index (no unit; values p10–p90 0–1 against 7–40 U/L for 6796; counts that match the hemolytic and icteric index items 6471/6472; KS 0.929). Lipase is 6796 only.
27. **D27 Hemoglobin unit:** the final unit is g/dL. Blood-gas 5257 is gated in g/L (g/dL rows × 10) and then reported as g/dL (÷ 10), like 5099. The Stage 2 v1.0.1 label ("→ g/dL" on a g/L conversion) is corrected in channel map v2 and the Stage 3 builder; v1.0.1 stays frozen (confirmed by Vinai, Nov 15, 2026). 5257 keeps its own mask.
28. **D28 S/F dropped.** Limitation recorded: PIC has no charted FiO2, oxygen-device or ventilation item; the only FiO2 is blood-gas 5222, which is recorded only with a pO2 on the same draw, so S/F would add fewer than 10 train draws.
29. **D29 P/F (final, Nov 15, 2026):** **P/F and the pSOFA respiratory component are dropped**, a documented limitation. The analyser FiO2 5222 is likely under-recorded (left at 21%): it is in percent (its train rows lie in 21–100 apart from a suppressed remainder under 10), yet every out-of-gate pair is above 700, about half of the in-gate pairs have P/F ≥ 400, and the analyser's own P/F 5243 agrees with the derived one. pO2 5239 and SaO2 5252 stay as separate channels; pO2 5244 is dropped as an exact same-draw duplicate of 5239; the analyser P/F 5243 stays dropped (D33); FiO2 5222 is dropped as a model channel (not kept indicator-only, because it is recorded only with a same-draw pO2, so its presence duplicates the blood-gas mask). pSOFA is a **5-component score without the respiratory component** (D29): cardiovascular, coagulation, hepatic, renal and neurological. D19 (P/F cut points) is superseded.
30. **D30 Urine output in mL/kg/h:** volume / weight / hours since the previous entry, capped at 24 h; the first entry uses the hours since ICU admission (t0). The value goes in the bin that contains its charting time, never back-filled; it is carried forward at most 24 h; days without an entry stay missing. A derived **KDIGO urine-output stage** uses the 24 h criteria: stage 3 if < 0.3 mL/kg/h for 24 h or anuria ≥ 12 h; a daily average < 0.5 mL/kg/h is a looser stand-in for stages 1–2. The 6 h and 12 h criteria cannot be assessed from once-daily charting (about 98% of entries between 05:00 and 07:59). Fluid intake stays mL/kg/day.
31. **D31 Magnesium:** the 'ug/dl' rows (fewer than 10 in train) are dropped.
32. **D32 Troponin T:** the per-row pg/mL ÷ 1000 conversion stays, and a plausibility gate is added (0–50 ng/mL, confirmed by Vinai Nov 15, 2026).
33. **D33 Panel masks:** one shared mask per panel (cytokines, TEG, urine dipstick, blood gas; column `mask_group` in channel map v2), because the items of a panel are measured together. The blood-gas mask covers 16 analyser items; ionized calcium 5215 and hemoglobin 5257 keep their own masks (Vinai, Nov 15, 2026). The analyser P/F 5243 (`pf_ratio_lab_crosscheck`) is dropped; it duplicates the derived P/F (difference p10–p90 −1.0 to 1.0).
34. **D34 DOB:** `age_invalid` stays in the primary analysis. Sensitivity analyses: exclude `age_invalid` encounters; and treat negative ages under 1 day as age 0.
35. **D35 first_careunit:** kept in Layer A. The write-up states that chart-channel masks largely reflect the care unit. A unit-stratified sensitivity analysis is added.
36. **D36 Raw-table hashes:** PATIENTS, ADMISSIONS, ICUSTAYS and SURGERY_INFO are pinned as gated inputs with the SHA-256 values in §4.3.
37. **D37 Age percentiles:** descriptive age percentiles come from valid-age encounters only, with n rounded to the nearest 10.

## 16. Contract and inventory gaps found while drafting
- **No channel-level ITEMID map exists yet.** The stage 07a lab-dictionary audit has not run; contract checklist §9 [AUDIT] and the Table 5.1 [DECISION] are open. The PIC_Analysis_3_0 contract §11 ITEMID list for the 16 channels is not among the Foundations docs.
- **Maps:** resolved Nov 14, 2026. The four v1 maps are decided and pinned (§4.2). Still open:
  - the M3 turnarounds (bracketed);
  - the remaining low-confidence rows. Channels: osmolality 5229, urine ketones 5356, urine leukocyte esterase 5359, CSF Pandy 5377, CSF RBC 5382, insulin 6172, troponin I 6192, Mycoplasma RNA 6806 (odd unit), influenza A/B RNA 6966 / 6967. Drug: posterior pituitary injection (left out of the vasoactive group). Specimen (medium): the general smears LIS035635 / LIS0536 and anaerobic culture LIS0567 (placed in `culture_site_unknown` by the D12 principle), and the Cryptococcus smear LIS0400 (CSF).
- **No ventilation item and no blood-gas sample type** in the inventory. The pSOFA respiratory points that need respiratory support, and arterial-only P/F, cannot be applied as published; the score stays "modified". Since Nov 15, 2026 pSOFA has no respiratory component at all (D28, D29).
- **Fluid balance:** dropped (D17); fluid intake in mL/kg/day is the fluid-status surrogate, an acknowledged limitation.
- **Oxygenation (Nov 15, 2026):** no charted FiO2, oxygen-device or ventilation item exists; the only FiO2 is blood-gas 5222. S/F is dropped (D28). P/F and the pSOFA respiratory component are dropped (D29): the analyser FiO2 is likely under-recorded.
- **Urine output (Nov 15, 2026):** charting is once daily, so only the 24 h KDIGO urine-output criteria can be applied (D30).
- **The 3.0 Layer B and Layer C feature specs are not among the Foundations docs.** Only the 3.0 Layer A spec and modeling contract were found. The 3.0 Layer A spec is partly superseded (`sex_missing`, the night flag, SURGERY_BEGIN-based early surgery).
- **Layer C horizon:** the frozen 3.0 Layer C file covers 28 bins (168 h). A2.4 extends Layer C to 56 bins (336 h), but no extension spec covers bins 28–55 yet. Contract v12 H2 and the checklist landmark line (24/48/72 h) predate A1.4/A2.4.
- **Surgery event bin:** resolved Nov 14, 2026 (decision 8; amendment A2.10): the bin of T, with the begin-time bin as a sensitivity analysis.
- **Lab timing text:** plan v12 S9 still says routine labs keep the CHARTTIME rule. Decision 11 / A2.3 supersede it (release + 1 h everywhere); the plan text is stale.
- **Urine output / fluid intake source and weight:** resolved Nov 14, 2026 (decision 4; amendment A2.12): CHARTEVENTS; weight imputed from age and gender for now.
- **Informative missingness:** resolved Nov 14, 2026 (decision 1; amendment A2.13): train-only, descriptive.
- **Pre-ICU lab block:** approved Nov 14, 2026 (decision 11; amendment A2.11).
- **Layer C count column:** still open (carried over from v12). The plan's Layer C φ0 lists "count", which conflicts with S2 unless it is dropped.

## Changelog
- **v0.9.2 (Nov 15, 2026):**
  - Recorded the final D24 and D29 (§15) with the follow-up results as evidence (§5.4 "Result"). P/F and the pSOFA respiratory component are dropped; pSOFA has 5 components (§5.3 item 4, §9 score inputs, §15, §16). Bilirubin 5255 is merged only in NICU, with sensitivities and a neonatal hepatic-scoring note.
  - New **channel map v2.1** `feature_channel_map_v2_1.csv` (686 rows) SHA-256 `7d1c060a228ab0772f35ae27e43a30971bcc13c36c98fcce2e369640b91bc899`, the current pin; v2 (`cbfdc68f…`) is unchanged and remains the follow-up script's pin. The follow-up script is not re-pinned and not rerun.
- **v0.9.1 (Nov 15, 2026):**
  - Vinai confirmed the four open map details: (1) the final hemoglobin g/dL is set in channel map v2 and Stage 3, and Stage 2 v1.0.1 stays unchanged (D27, confirmed); (2) troponin T gate 0–50 ng/mL (D32; replaces the 0–100 placeholder); (3) ionized calcium 5215 gate 0.3–3.0 mmol/L (D21); (4) the blood-gas shared mask has 16 items, with 5215 and 5257 keeping their own masks (D33).
  - Channel map v2 was revised **in place** (same file name, no run had used the first v2 hash) and re-pinned: SHA-256 `cbfdc68f11bc2f8d9f96150803a1339930da81eb091f481576f709debdd5ab5c` (686 rows). The superseded first v2 hash `9872e688d4d6d06792eda934c909f149f00b69f9bcb3955297c049c04c3b9eaa` must not be used.
  - `stage2_followup_checks.py` re-pinned to the revised map: SHA-256 `e955974d8e9fd85414482549cc14e019d922215505088220ae9cfba990992f62` (logic unchanged; 9 synthetic tests pass).
- **v0.9 (Nov 15, 2026):**
  - Recorded the Stage 2 review decisions D21–D37 (§15). Changed: §4.2 (channel map v2 pin, which pin applies, v2 changes), §4.3 (D36 raw-table hash gates), §5.3 (outcome note, S/F dropped, hemoglobin unit), the new §5.4 (follow-up checks), §6 (D34, D35, D37), §12 and §16.
  - Channel map v2 `feature_channel_map_v2.csv` (686 rows) SHA-256 `cbfdc68f11bc2f8d9f96150803a1339930da81eb091f481576f709debdd5ab5c`; the other four map pins are unchanged.
  - New script `stage2_followup_checks.py` v1.0, SHA-256 `e955974d8e9fd85414482549cc14e019d922215505088220ae9cfba990992f62`. Stage 2 v1.0.1 is unchanged and stays pinned to the v1 maps.
- **v0.8 (Nov 15, 2026):**
  - Recorded D20 (§6, §15): invalid DOB or age is no longer fatal at any count. Those encounters go to an `age_invalid` band; the audit reports `layer_a.dob_breakdown` (status, negative bins, by split and first_careunit). Age-continuous quantiles exclude invalid ages. Weight imputation for `age_invalid` uses sex-only, else overall. A safety gate stops the run if invalid ages exceed 1% of the cohort. The Stage 2 script is v1.0.1. No map or hash gate changed.
- **v0.7 (Nov 14, 2026):**
  - Recorded D19 (§5.3 item 4, §15): P/F is scored on the pSOFA P/F cut points without the respiratory-support requirement, like S/F.
  - Recorded that Vinai accepted keeping calcium acetate granules and oral magnesium sulfate powder in `none` (D18).
  - No map changed; the map hashes in §4.2 are unchanged.
  - Later the same day: the Stage 2 script (`feature_sparsity_frequency_audit.py` v1.0) was written. It implements this v0.7 spec (D1–D19). No map or decision changed.
- **v0.6 (Nov 14, 2026):**
  - Recorded Vinai's decisions D14–D18 (§15): 24 h cutoffs (0.5× / 2×) for the new groups are confirmed; S/F uses the pSOFA cut points without the respiratory-support requirement, with both limitations documented; promethazine and chlorpromazine are antiemetics only; fluid balance is dropped and fluid intake in mL/kg/day is the fluid-status surrogate; a new electrolyte_replacement drug group.
  - Rebuilt and re-pinned the channel map (`56e8bf72fcbc…`) and the drug-group map (`6fa7d8a935f0…`) (§4.2). The specimen map, organism map and neuro set are unchanged.
  - Updated §4.2, §5.3 items 4, 5 and 8, §7, §15 and §16.
- **v0.5 (Nov 14, 2026):**
  - Pinned the neuro admitting-diagnosis set v1 (`neuro_admitting_dx_icd10_v1.csv`, `3cc8f4acfe9d…`), as confirmed by Vinai:
    - CNS infection, head injury, stroke/bleed and encephalopathy are included, with their partial codes (Q28 brain AVM only);
    - seizures, CNS tumours and hydrocephalus / malformations are excluded.
  - Replaced the pending-draft gate (§12) with an exact hash gate.
  - Updated §4.2, §5.3 item 8, §15 D3 and §16. No other change.
- **v0.4 (Nov 14, 2026):**
  - Recorded Vinai's map decisions D1–D13 (§15, "Map decisions").
  - Pinned the v1 channel, drug-group, specimen and organism maps by full SHA-256 (§4.2); the neuro admitting-diagnosis set is drafted and pending.
  - Drug groups updated: one vasoactive group; antifungal, antiviral, antiemetic and non-opioid analgesic groups; an osmotic group resolved by admitting diagnosis; clonazepam and nitrazepam moved to AED.
  - Added §5.3 (Stage 2 map checks): duplicate-item distributions, HR merge, MAP pairs, P/F and S/F coverage, the urine-output and intake charting interval, units, microbiology fill rates, osmotic resolution and the vasoactive count.
  - Updated §7 (sources), §12 (hash gate incl. the neuro set), §13 (tests) and §16.
  - No change to the Stage 1 script; the Stage 2 script is not written yet.
- **v0.3 (Nov 14, 2026):** recorded Vinai's follow-up confirmations (§15): NMB and AED cutoffs 24 h with 0.5× / 2× sensitivity; pre-ICU lab window [t0 − 6 h, t0) (sensitivity [t0 − 24 h, t0)); time-since-last cap 24 h for channels without an LOCF cap (replacing the 6 h proposal); urine output from CHARTEVENTS; the Stage 1 canonical-encounters hash gate approved. The A2.7 "start bin" item is closed as superseded. "Not resolved by the decisions" is replaced by "Follow-up confirmations". No change to the Stage 1 script.
- **v0.2 (Nov 14, 2026):** recorded Vinai's answers to all 12 open questions (§15). The [OPEN] and [PROPOSAL] markers are resolved: §9 thresholds approved as starting defaults; §8 outcome descriptors approved (train only, descriptive); NMB and AED groups and the 24 h cutoffs; CHARTEVENTS as the fluid-intake source and the weight rule; full prefix primary for §9.1; time since last measurement included (capped); strata fixed; surgery event in the bin of T; raw hashes recorded in Stage 1 and gated in Stage 2; Layer A pre-ICU lab block approved; gap bins included as primary. §5.1 now describes the written Stage 1 script. Added "Not resolved by the decisions".
- **v0.1 DRAFT (Nov 12, 2026):** first draft with 12 open questions.
