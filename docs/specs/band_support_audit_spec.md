# Spec: phi1 duration-band support audit (read-only, aggregate-only)
Version 1.0, Nov 12, 2026. A results section was added Nov 12, 2026 after Vinai's MAP_v1 run; the script is unchanged. It implements the approved phi1 band definition: decision 16 and its follow-up (Nov 11, 2026), and the duration window (Nov 12, 2026). See the duration-audit spec v1.7, post-audit decision 4, and the contract v12.1 amendment A2.7.

## Purpose
Show how well the approved type-specific duration bands are supported in each split before they are frozen for modeling:
- how many surgeries get their own type median versus a fallback;
- the band distribution;
- which own-median types have thin bands;
- how often validation and test types are unseen in train.

## Inputs and conventions (unchanged from the duration audit v1.7, MAP_v1)
- **Files:** canonical encounters, harmonization fields, complexity map v1 and subject splits, gated by exact full SHA-256:
  - map `5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5`;
  - canonical `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418`;
  - harmonization `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`;
  - split `60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31`.
- **Shared code:** the duration audit v1.7 module (SHA-256 `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`) is verified, then imported. It supplies:
  - input checks, UTF-8 reading and column resolution with fallbacks;
  - ROW_ID keying (all SURGERY_INFO rows reconciled) and encounter linkage (split from the linked encounter);
  - map handling and the duration rule;
  - the disclosure registry, primary and complementary suppression, and the recoverability solver.
- **Population:** surgeries linked to a cohort encounter, with a mapped name.
  - Type = `standardized_name`.
  - System = `body_system`.
  - Complexity group: cardiac STAT 1–2 / 3 / 4–5 / not classified (G1 / G2 / G3 / NC); non-cardiac I1 / I2 / I3; Unclassified system → `unclassified`.
- **Duration** D = SURGERY_END − SURGERY_BEGIN (h). If a surgery time is missing, D falls back to ANES_END − ANES_START (counted).
  - Excluded and counted by reason: D missing or ≤ 0.
  - D > 24 h is counted, not excluded (an open item).

## Computation (train split only for medians; frozen)
1. **Reference median per surgery.** Use the own type median if the type has ≥ 10 train surgeries with a valid D. Otherwise fall back, in order, to system × complexity group, then system, then all surgeries, taking the first level with ≥ 10.
   - Sensitivity: the same with a threshold of 20.
   - Validation and test use the frozen train medians.
   - Types unseen in train go straight to the fallback.
2. **Bands.** r = D / reference median, rounded to 9 decimals.
   - short: r < 0.75
   - expected: 0.75 ≤ r ≤ 1.33 (inclusive)
   - long: r > 1.33
3. **Tertile sensitivity.** For own-median types only, use within-type train tertile cuts q1 and q2:
   - short: D < q1
   - expected: q1 ≤ D ≤ q2
   - long: D > q2

## Outputs (per split: train, validation, test)
- **Reconciliation:** link classes sum to all SURGERY_INFO rows. Per split, the population and invalid-duration reasons, the anesthesia-fallback count and the > 24 h count.
- **Median-source level:** types and surgeries at each level, at thresholds 10 and 20, plus the types split across levels.
- **Bands:** counts and percentages overall, by system and by complexity group (threshold 10), and overall at threshold 20.
- **Own-median types:** n per band, with `any_band_below_10` = Y/N.
- **Unseen in train:** validation and test surgeries and types whose type is unseen in train, with the share.
- **Tertile sensitivity:** the band distribution.
- **Frozen-medians table:** a separate shareable CSV with `threshold, key, level, train_n, median_hours`. Every key has train n ≥ 10. train_n is suppressed (`S`) where differencing would reveal a small count. Full-precision medians stay in the `_raw_LOCAL` CSV.

## Disclosure
- **Aggregate only.** Counts from 0 to 9 are published as `<10`. Complementary `S` is added until the LP recoverability check finds nothing recoverable.
- **Identities:** all published identities are declared, including:
  - pool and container totals versus the own-type rows, with hidden remainders;
  - threshold-20 own-level counts as a subset of threshold-10 counts.
- **Percentages:** published only where the count and its denominator are published.
- **Never published:** min/max, IDs and timestamps.
- **Raw outputs** are written as `_raw_LOCAL` and are never shared.
- **Shared JSON/MD/CSV** are written only after the disk-reload release gate passes.

## Acceptance checks (the script stops on failure)
1. Any input-hash or v1.7 module-hash mismatch. The hash is checked before the file is read or imported.
2. v1.7 input checks: headers, row counts, missingness, cohort and split totals, ROW_ID unique.
3. Link classes not summing to the SURGERY_INFO rows, or any declared identity failing on raw counts.
4. A frozen key with train n below the threshold, or a median computed from non-train rows.
5. A release-gate failure.

## Tests (synthetic data only)
48 pytest tests cover:
- the fallback order;
- thresholds 10 and 20, including exactly at n;
- band boundary inclusivity at 0.75 and 1.33, including after the ratio is computed;
- the surgery-to-anesthesia fallback and invalid durations;
- > 24 h kept;
- train-only medians, with no validation/test leakage;
- unseen types;
- tertiles;
- suppression and the release gate, with eight tamper cases;
- hash-gate failures stopping the run;
- a short console summary.

## Results (MAP_v1 run by Vinai, Nov 12, 2026; aggregate only)
- **Run:** script `2d02154f…` with v1.7 module `ce894b20…`. Exit 0, release gate ok. Suppression: 755 primary and 293 complementary cells, 0 recoverable.
- **Train surgeries with a valid duration:** 5,141. The all-surgeries train median is 1.53 h.
- **Bands (short / expected / long):**

  | split | short | expected | long |
  |---|---|---|---|
  | train | 22.3% | 50.1% | 27.5% |
  | validation | 21.8% | 49.5% | 28.7% |
  | test | 23.8% | 49.7% | 26.6% |

- **min n = 10:** 73 of 1,083 train types have their own median, covering 3,250 surgeries (63%). 990 types (1,859 surgeries) fall back to system × complexity group. The rest (system / all levels) are suppressed. Frozen keys: 73 own types, 18 system × group pools, 7 systems, all.
- **min n = 20 (sensitivity):** 41 types with their own median, covering 2,803 surgeries. Frozen keys: 41 / 17 / 6 / 1.
- **Types unseen in train:** 12.9% of validation and 12.5% of test surgeries.
- **Tertile sensitivity (own-median types):** about 33% / 35% / 32%.
- **Sparsity:** many own-median types have at least one band cell under 10.
- **Data quality:** durations > 24 h and ordering violations are each < 10 in every split.
- **Frozen medians file:** `phi1_band_support_frozen_medians_MAP_v1.csv`, SHA-256 `e73c8b5a8d780b51edfd0935abdb05335f57e3a0cfc0c4c52eaaa85888db0047`.

**Decisions (Vinai, Nov 12, 2026):**
- **Freeze** the train medians at **min n = 10**. The frozen medians file is `phi1_band_support_frozen_medians_MAP_v1.csv`, SHA-256 `e73c8b5a8d780b51edfd0935abdb05335f57e3a0cfc0c4c52eaaa85888db0047`.
- **Coding:** the duration band is coded as a **3-level ordinal (short / expected / long)** alongside the type and system descriptors. It is **not** a full type × band cross-product, because within-type band cells are sparse.
- **Type naming (resolved):** the type median uses `standardized_name`, with the fallback hierarchy for rare types. Collapsing combination names to their primary component is a **sensitivity analysis only**.

## Open (not part of this audit)
- Handling of durations > 24 h and ordering violations. Each is < 10 in every split.
- Whether the expected duration appears in the Layer C start bin.
- *Resolved Nov 12, 2026:* the type median uses `standardized_name`, with combination-name collapsing as a sensitivity analysis only.
