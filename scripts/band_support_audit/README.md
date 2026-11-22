# phi1 duration-band support audit (Copilot run instructions)

Read-only, aggregate-only audit of the approved phi1 duration bands (decision 16 and its follow-up, Nov 11–12, 2026).
Spec: `band_support_audit_spec.md` (v1.0).

## What it does
- **Duration:** D = SURGERY_END − SURGERY_BEGIN (hours). If a surgery time is missing, D falls back to ANES_END − ANES_START (counted).
  - D missing or ≤ 0 is excluded and counted by reason.
  - D > 24 h is counted, **not** excluded.
- **Reference median:** fit on the **train split only** and frozen.
  - Fallback order: own type (`standardized_name`) with ≥ 10 train surgeries → system × complexity group → system → all surgeries. The first level with ≥ 10 is used.
  - The threshold is 10 (primary) and 20 (sensitivity).
  - Types unseen in train go straight to the fallback.
- **Bands:** r = D / reference median.
  - short: r < 0.75
  - expected: 0.75 ≤ r ≤ 1.33
  - long: r > 1.33
  - Sensitivity: within-type train tertiles, for own-median types only.
- **Reused unchanged from the duration audit v1.7** (`audit_surgery_complexity_duration.py`, SHA-256 `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`):
  - the inputs and the four full-hash gates;
  - the surgery key, `ROW_ID` (all 7,488 rows reconciled);
  - encounter linkage and split;
  - map handling and the duration rule;
  - primary and complementary suppression, and the recoverability solver.
- **Pinned import:** the script verifies that file's SHA-256 **before** importing it.
- **Population:** surgeries linked to a cohort encounter (v1.7 linkage) with a mapped name.

## Run order for Copilot
1. Copy `phi1_band_support_audit.py` into the project `scripts/` folder, next to the existing `audit_surgery_complexity_duration.py` (v1.7). Do **not** edit either file.
2. Verify both SHA-256 values before running (PowerShell, from the project root):
   ```
   Get-FileHash scripts\phi1_band_support_audit.py -Algorithm SHA256
   Get-FileHash scripts\audit_surgery_complexity_duration.py -Algorithm SHA256
   ```
   - The first must equal `2d02154ff1275d53a6ccb4ce947fc8b78cf8e6dcd8b381875febee5bc99821f6`.
   - The second must equal `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`.
   - If either differs, stop and re-copy. Do not patch.
3. Run it unedited, from the project root:
   ```
   python scripts\phi1_band_support_audit.py --raw-dir "<PIC_RAW_DIR>" --split-file "<SPLIT_FILE>" --harm-t0-col t0 --harm-out-col OUTTIME
   ```
4. **Stop on any failure, without patching.**
   - Exit 2 (`AUDIT STOPPED`) means an input, hash or reconciliation failure.
   - Exit 3 (`RELEASE GATE FAILED`) means the release gate failed, and no shared file is left in place.
   - In both cases, send back the console text only.
5. On success (the last console line is `release gate ok`), send back:
   - the short console summary;
   - the three shared files.

## Outputs (`data/interim/`)
| file | share? |
|---|---|
| `phi1_band_support_audit_MAP_v1.json` | yes (aggregate, suppressed) |
| `phi1_band_support_audit_MAP_v1.md` | yes (rendered from the JSON) |
| `phi1_band_support_frozen_medians_MAP_v1.csv` | yes (`threshold, key, level, train_n, median_hours`; train_n suppressed where needed; no key with n < 10) |
| `phi1_band_support_audit_MAP_v1_raw_LOCAL.json` | **never** (unsuppressed counts) |
| `phi1_band_support_frozen_medians_MAP_v1_raw_LOCAL.csv` | **never** (full-precision medians and raw train n) |

The shared files are written only after the disk-reload release gate passes. The gate checks:
- no published count from 0 to 9;
- suppression counters and published identities;
- that LP recoverability finds nothing recoverable;
- that percentages appear only with their published count and denominator;
- no min/max, no timestamps;
- that the md and CSV equal the rendering of the JSON.

## Reported per split (train, validation, test)
- types and surgeries at each median-source level (min n 10 and 20);
- band counts and percentages: overall, by system and by complexity group (min n 10), plus overall at min n 20;
- per own-median type: n per band, suppressed under 10, and an `any_band_below_10` flag;
- validation and test: surgeries and types whose type is unseen in train (share in %);
- the tertile sensitivity band distribution.

## Notes
- Complementary suppression can hide a headline total. This happens when the structure makes a small cell exactly derivable (e.g., every fallback level is 0). The protection is deliberate.
- Tests: `test_phi1_band_support_audit.py` (48 tests, synthetic data only) runs where `../duration_audit_v2/` holds the v1.7 module and `synthetic_data.py`:
  ```
  python -m pytest -q test_phi1_band_support_audit.py
  ```
  It does not need to be run against PIC.
