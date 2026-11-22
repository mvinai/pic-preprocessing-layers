# Feature item inventory, Stage 1 (Copilot run instructions)

This is Stage 1 of the feature sparsity and measurement-frequency audit. It is read-only and aggregate-only.
- Spec: `feature_sparsity_frequency_audit_spec.md` v0.2, §5.1.
- Decisions: Vinai, Nov 14, 2026.
- Script: `feature_item_inventory.py` (v1.0).

## What it does
- **Raw tables read:** `CHARTEVENTS`, `LABEVENTS`, `PRESCRIPTIONS`, `MICROBIOLOGYEVENTS`, `OUTPUTEVENTS`, `INPUTEVENTS`, plus the `D_ITEMS` and `D_LABITEMS` dictionaries.
  - Every file read has its SHA-256 recorded in the output. The hashes are recorded here and gated in Stage 2.
  - It stops if a table is missing.
  - It also stops if the rows read differ from the schema inventory: CHARTEVENTS [N_chartevents omitted]; LABEVENTS [N_labevents omitted]; PRESCRIPTIONS [N_prescriptions omitted]; MICROBIOLOGYEVENTS [N_microbiologyevents omitted].
- **Item keys:**
  - `ITEMID`, with the label, category and fluid joined from the dictionaries;
  - drug name pair, `DRUG_NAME | DRUG_NAME_EN`;
  - specimen, `SPEC_ITEMID | SPEC_TYPE_DESC`;
  - organism, `ORG_ITEMID | ORG_NAME`;
  - if `INPUTEVENTS` has no `ITEMID` column, it is keyed by unit (`unit:<AMOUNTUOM>`).
- **Linkage (contract v12 H3):** `SUBJECT_ID` and `HADM_ID` must be in the encounter's `hadm_id_list`, and the event time must be in [t0, OUT).
- **Event times:**
  - LABEVENTS uses release time = CHARTTIME + 1 h;
  - PRESCRIPTIONS uses STARTDATE;
  - all other tables use CHARTTIME.
- **Per item, train split only, cohort encounters only:**
  - rows in stay;
  - rows with a numeric value;
  - encounters with ≥ 1 row in the stay;
  - encounters with ≥ 1 row in the Layer B window [t0, min(OUT, t0+24 h));
  - encounters with ≥ 1 row in the pre-ICU window [t0−6 h, t0);
  - the unit strings seen.
- **Row reconciliation per table, all splits:** each row read is classified as in stay (train, validation or test), cohort admission but outside the stay, no cohort admission, or time missing. The classes must sum to the rows read.
- **Search hits:** label and name keyword hits, in English and Chinese, for weight, height, urine output, fluid intake, GCS, pupils, FiO2, ventilation, vasoactives, pain, occult blood, endocrine, sodium, NMB and AED. These are keys for Vinai to verify, never assignments.
- **Blank map templates:** channel ITEMID map, drug-group map (groups include `nmb` and `aed`), specimen map and organism map.
  - They list only items with ≥ 10 train encounters.
  - The mapping columns are left blank for Vinai to fill.
- **Reused unchanged from the duration audit v1.7** (`audit_surgery_complexity_duration.py`, SHA-256 `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`): ID parsing, the disclosure registry, primary and complementary suppression, and the recoverability solver.
  - The script verifies that file's SHA-256 **before** importing it.
- **Input hash gates (exact 64-character match):**
  - harmonization `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`;
  - split `60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31`;
  - canonical encounters `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418` (needed for `hadm_id_list`).
- **Cohort check:** the cohort and core-landmark eligibility (OUT > t0 + L at 24 / 48 / 72 / 120 / 168 / 240 / 336 h) must reproduce the frozen pins (values omitted from this public copy). Otherwise the script exits 2.

## Run order for Copilot
1. Copy `feature_item_inventory.py` into the project `scripts/` folder, next to the existing `audit_surgery_complexity_duration.py` (v1.7). Do **not** edit either file.
2. Verify both SHA-256 values before running (PowerShell, from the project root):
   ```
   Get-FileHash scripts\feature_item_inventory.py -Algorithm SHA256
   Get-FileHash scripts\audit_surgery_complexity_duration.py -Algorithm SHA256
   ```
   - The first must equal `86388f7b628a4ae64d3cf59467213c591962ff94c7867e0133a37657f1d9aa41`.
   - The second must equal `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`.
   - If either differs, stop and re-copy. Do not patch.
3. Run it unedited, from the project root:
   ```
   python scripts\feature_item_inventory.py --raw-dir "<PIC_RAW_DIR>" --split-file "<SPLIT_FILE>" --harm-t0-col t0 --harm-out-col OUTTIME
   ```
   It reads about 14 million raw rows in 1-million-row chunks; expect several minutes.
4. **Stop on any failure, without patching.**
   - Exit 2 (`AUDIT STOPPED`) means an input, hash, cohort, eligibility, row-count or reconciliation failure.
   - Exit 3 (`RELEASE GATE FAILED`) means the release gate failed, and no shared file is left in place.
   - In both cases, send back the console text only.
5. On success (the last console line is `release gate ok`), send back:
   - the short console summary;
   - the shared files listed below.

## Outputs (`data/interim/`)
| file | share? |
|---|---|
| `feature_item_inventory.json` | yes (aggregate, suppressed) |
| `feature_item_inventory.md` | yes (rendered from the JSON) |
| `feature_channel_map_TEMPLATE.csv` | yes (items with ≥ 10 train encounters; blank map columns) |
| `drug_group_map_TEMPLATE.csv` | yes (same rule; `group` blank; allowed groups listed in the JSON `config.drug_groups`) |
| `specimen_map_TEMPLATE.csv` | yes (same rule) |
| `organism_map_TEMPLATE.csv` | yes (same rule) |
| `feature_item_inventory_raw_LOCAL.json` | **never** (unsuppressed counts) |

The shared files are written only after the disk-reload release gate passes, first as pending files, then as final files, with the gate run on both. The gate checks:
- header and schema: every count field is declared;
- no published count from 0 to 9; counts 0–9 are shown as `<10`, and complementary suppression is shown as `S`;
- suppression counters and published identities;
- that LP recoverability finds nothing recoverable;
- that units are published only when the item's rows are published;
- no min/max and no timestamps;
- that the md and templates equal the rendering of the JSON;
- no template row with fewer than 10 train encounters;
- no `_raw_LOCAL` file among the shared files.

Items are listed in key order, never by count.

## Notes
- Only train-split counts are published per item. Validation and test enter only the row reconciliation.
- Tests: `test_feature_item_inventory.py` (61 tests) uses synthetic data only, from `synthetic_inventory_data.py`.
  - It runs where `../duration_audit_v2/` holds the v1.7 module:
    ```
    python -m pytest -q test_feature_item_inventory.py
    ```
  - It never opens real patient data and does not need to be run against PIC.
