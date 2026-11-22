# Feature sparsity and measurement-frequency audit, Stage 2 (Copilot run instructions)

This is Stage 2 of the feature sparsity and measurement-frequency audit. It is read-only and aggregate-only.
- Spec: `feature_sparsity_frequency_audit_spec.md` v0.8 (§5.2, §5.3, §6–§9, §11–§13; decisions D1–D20).
- Script: `feature_sparsity_frequency_audit.py` (v1.0.1), SHA-256 `3b5b36fbd0f73420beaa74b18abf2afec4725cf86c5e8efe514be19120043213`.
- Tests: `test_feature_sparsity_frequency_audit.py` with `synthetic_stage2_data.py` (synthetic data only).

## Change in v1.0.1 (D20, Nov 15, 2026)
v1.0 stopped when 42 cohort encounters had a missing or invalid DOB. Vinai decided to continue with an `age_invalid` band and an aggregate DOB breakdown.
- Invalid ages (missing DOB, including no PATIENTS row; unparseable DOB; or t0 < DOB) are the `age_invalid` band; `age_days` is unset. Negative ages are binned <1 day / 1–30 days / >30 days.
- `layer_a.dob_breakdown` reports status, negative bins, and invalid counts by split and first care unit, with the usual <10 / complementary-`S` suppression. Age quantiles use valid ages only.
- Weight for `age_invalid` is not imputed from an age band: sex-only train median (n ≥ 10), else the overall train median.
- A safety gate still stops the run if invalid ages exceed 1% of the cohort.
- Every other hash gate is unchanged. Run the script unedited, as below.

## What it does
- **Hash gates.** Every gate is an exact 64-character SHA-256 match, checked before the file is used. Any mismatch stops the run with exit 2.
  - The duration audit v1.7 module `audit_surgery_complexity_duration.py`: `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`.
  - Canonical encounters `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418`.
  - Harmonization `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`.
  - Split `60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31`.
  - Surgery complexity map v1 `5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5`.
  - The v1 maps and the neuro set (in `data/interim/maps/`):

    | file | SHA-256 | rows |
    |---|---|---|
    | `feature_channel_map_v1.csv` | `56e8bf72fcbcaac765a39c2e705aedef61323e94841dbdcc8b2657191af370b6` | 685 |
    | `drug_group_map_v1.csv` | `6fa7d8a935f03ad71b4a7ecd8276f9a54f3a679f4cfc08c8b4504eda83455ead` | 343 |
    | `specimen_map_v1.csv` | `177b165ad8e5888a3999583a4949e93908be45af34cfa9c6e4937c40c08666f1` | 71 |
    | `organism_map_v1.csv` | `21f28448cb31c654c5312ca97883a34fe5cb6b269e8323973eb8b899bb3583d2` | 101 |
    | `neuro_admitting_dx_icd10_v1.csv` | `3cc8f4acfe9df29037b60d903ad4204c40db5bc76923de8b64234ccb2553b9d9` | 67 |

  - The raw PIC tables, using the hashes that Stage 1 recorded:

    | table | SHA-256 |
    |---|---|
    | D_ITEMS | `a7b65cf5d61ca84f6ce5d07d3775b01116473176a980463d4fe9df1dd12e1591` |
    | D_LABITEMS | `d4557986fe9c5475d0db6506fdf12801ecfb0e2e7a3c01caf6ac86bbab69ed60` |
    | CHARTEVENTS | `fe209d47d15825f180624f2cefb9bc38f8513e1ac0b3709d57aef37ae88f328e` |
    | LABEVENTS | `784ecb929a358ee1a18fc9beafbdf67f30a91646b780d9983714a5e9ce0018c5` |
    | OUTPUTEVENTS | `4c49f3c6910b4e8aaa817a13acf16d7120e5ecdd53cb6675bfc7760fda3d4d0d` |
    | INPUTEVENTS | `7d8aaebf51d82fdde32841a534f3d296911e92c3421885e40b379bc7c394f9c6` |
    | PRESCRIPTIONS | `61fa43dcb1ed977f4219c902b01977da2d898a16a2c33640997ad8337d2eb2c6` |
    | MICROBIOLOGYEVENTS | `e19cf9d41732f74911039acae6b2cd7bf5c5dd65179cbd31aa8634dd1e6bb45c` |

  - Stage 1 did not hash PATIENTS, ADMISSIONS, ICUSTAYS or SURGERY_INFO. Their SHA-256 is recorded in the output, and their row counts are gated: [N_patients omitted] / [N_admissions omitted] / [N_icustays omitted] / 7,488.
- **DOB / age (D20).** Missing, unparseable or negative ages become `age_invalid` (not a stop unless they exceed 1% of the cohort). See `layer_a.dob_breakdown` in the shared output.
- **Cohort.** The run must reproduce:
  - the frozen canonical cohort (counts omitted from this public copy);
  - splits of the frozen 70/15/15 split (sizes omitted from this public copy);
  - eligibility (OUT > t0 + L, strict) of [landmark eligible counts omitted] at the core landmarks;
  - in-hospital deaths (count omitted from this public copy).
- **Map rules (D1–D19) and the §5.3 checks.** All outputs are aggregate.
  - **Channels:**
    - channels come from the channel map only;
    - plausibility gates are applied, and out-of-gate rows are counted separately;
    - the HR merge keeps 1003 over 1002 at the same CHARTTIME;
    - MAP comes from same-CHARTTIME SBP/DBP pairs; pairs with DBP ≥ SBP or out of range are dropped.
  - **Oxygenation:**
    - P/F uses the same draw (5239 preferred) and is scored on the pSOFA cut points without the support requirement (D19);
    - S/F is computed only where there is no P/F, with SpO2 ≤ 97% and the nearest SpO2 within 5 minutes (D15);
    - the P/F-or-S/F coverage count is reported.
  - **Units:**
    - per-row conversions for 5257 and 6525;
    - magnesium `ug/dl` rows are held out and compared under three scalings;
    - any unrecognized conversion text in the map stops the run.
  - **Duplicate-item merge checks:** quantiles, KS and same-draw differences, plus the main-item-only sensitivity.
  - **Urine and intake:** the charting-interval check; intake is reported as mL/kg/day.
  - **MICROBIOLOGYEVENTS:** column fill rates.
  - **Prescriptions:** the in-stay / date-only start check.
  - **Osmotic agents:** classified by the ADMISSIONS admitting diagnosis using neuro set v1. Q28 is matched on 4 characters (Q28.2 / Q28.3 only). The all-osmotherapy sensitivity is also reported.
  - **Microbiology collections:** built on `base_code` only, so susceptibility add-ons never reach M1.
- **Metrics:** Layer A (§6); metrics 1–12 for Layers B and C per channel, drug group, specimen group and surgery (§7); the strata.
  - Outcome descriptors (§8, train only) are computed **only with `--outcome-descriptors`**. The flag's state is written to the output.
  - The run also reports the provisional §9 classes as labels for Vinai's review. They are not decisions.
- **Disclosure:**
  - counts 0–9 are published as `<10`;
  - complementary suppression `S` is applied until the v1.7 LP recoverability solver finds nothing recoverable;
  - quantiles, percentages and KS are published only where every count they rest on is published and ≥ 10;
  - no min/max, IDs or timestamps are published.

## Run order for Copilot
1. Put the files in place. Do **not** edit any file.
   - Copy `feature_sparsity_frequency_audit.py` into the project `scripts/` folder, next to `audit_surgery_complexity_duration.py` (v1.7) and `feature_item_inventory.py`.
   - Copy the five v1 CSVs from the shared `item_inventory_outputs/` folder into `data/interim/maps/`: `feature_channel_map_v1.csv`, `drug_group_map_v1.csv`, `specimen_map_v1.csv`, `organism_map_v1.csv` and `neuro_admitting_dx_icd10_v1.csv`.
2. Verify the hashes before running. In PowerShell, from the project root:
   ```
   Get-FileHash scripts\feature_sparsity_frequency_audit.py -Algorithm SHA256
   Get-FileHash scripts\audit_surgery_complexity_duration.py -Algorithm SHA256
   Get-FileHash data\interim\maps\*.csv -Algorithm SHA256
   ```
   - The script must equal `3b5b36fbd0f73420beaa74b18abf2afec4725cf86c5e8efe514be19120043213`.
   - The module must equal `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`.
   - The maps must equal the table above.
   - If any of them differs, stop and re-copy. Do not patch.
3. Run it unedited, from the project root, as one line:
   ```
   python scripts\feature_sparsity_frequency_audit.py --raw-dir "<PIC_RAW_DIR>" --split-file "<SPLIT_FILE>" --harm-t0-col t0 --harm-out-col OUTTIME --maps-dir data\interim\maps --outcome-descriptors
   ```
   - It reads about 14 million raw rows in 1-million-row chunks.
   - Suppression then runs an LP recoverability check over the whole output, so expect a long run (tens of minutes). Let it finish.
4. **Stop on any failure, without patching.**
   - Do not edit the script, the maps or any input, and do not rerun with changed options.
   - Exit 2 (`AUDIT STOPPED`) means an input, hash, cohort, eligibility, row-count, reconciliation or acceptance failure.
   - Exit 3 (`RELEASE GATE FAILED`) means the release gate failed, and no shared file is left in place.
   - In both cases, send back the full console text and the `Get-FileHash` output only.
5. On success (the console ends with `release gate ok`), send back:
   - the `Get-FileHash` output from step 2;
   - the full console text, which prints only published (suppressed) values and the input hashes;
   - the two shared files listed below.

## Outputs (`data/interim/`)
| file | share? |
|---|---|
| `feature_sparsity_frequency_audit.json` | yes (aggregate, suppressed) |
| `feature_sparsity_frequency_audit.md` | yes (rendered from the JSON) |
| `feature_sparsity_frequency_audit_raw_LOCAL.json` | **never** (unsuppressed counts). Do not open, copy, attach, upload or paste it anywhere. |

The shared files are written only after the disk-reload release gate passes. They are first written as pending files, then as final files, and the gate runs on both. The gate checks:
- the header, the flag state and the count schema;
- no published count from 0 to 9;
- the suppression counters and the identities;
- a re-run of LP recoverability;
- quantile, percentage and dependent rules;
- no min/max, no timestamps and no sentinel identifiers;
- that the md equals the rendering of the JSON;
- that there are no outcome keys when the flag is off.

## Notes
- Only aggregate counts are produced. Metric 9 (co-measurement) and the per-bin outputs are train only. Validation and test appear in the per-split metrics, the strata and the reconciliation.
- Tests use synthetic data only and never open real patient data:
  ```
  python -m pytest -q test_feature_sparsity_frequency_audit.py
  ```
  - They need `synthetic_stage2_data.py` next to the test file, and the v1.7 module next to the script or in `../duration_audit_v2/`.
  - They take about 6 minutes (81 passed, 1 skipped in the v1.0.1 suite).
  - One test parses the real v1 maps (decision files, no patient rows) only if `data/interim/maps/` exists; otherwise it is skipped.
  - The tests do not need to be run against PIC.
