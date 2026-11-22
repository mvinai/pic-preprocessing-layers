# Copilot spec: Landmark-forward band-support audit (read-only, aggregate-only)
Version 2.3, Nov 10, 2026 (5-level Outcome 2; disclosure-control acceptance items; see Changelog). Implements contract v12.1 draft amendment A1.8 (`contract/v12_1_amendment_landmark_forward_outcome_DRAFT.md`). Style and guardrails follow `specs/landmark_eligibility_audit_spec.md` (v2).

## Purpose
Before amendment A1 is frozen, measure how well each level of the two landmark-forward ordinal outcomes (A1.2) is supported at each of the 14 landmarks (A1.4) and each split. Vinai uses the result to decide level merges and the policy for readmission stays that are also index encounters. Both are post-audit decisions, and the script applies neither. This is a design audit only: it builds no features, fits no models, and changes no frozen files.

## Primary decision target
- **Primary criterion: support on both sides of every cutpoint of Outcome 1 (remaining ICU time, 9 levels), in the test split, at the core landmarks.**
- Outcome 2 (post-ICU course, 5 levels as of v2.3) is the key secondary outcome. It is audited the same way.
- Mortality is context only: Outcome 1 level 8 is the same as Outcome 2 level 4.

## Inputs (read only; do not modify)
1. **Frozen canonical encounter CSV**, `data/interim/canonical_icu_encounters.csv` (frozen 3.0 cohort).
   - Compute its SHA-256 **before reading** and compare it with `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418`. On mismatch, raise an error and stop.
   - Required columns: `encounter_id`, `SUBJECT_ID`, `hadm_id_list`, `HADM_ID`, `INTIME`, `OUTTIME`.
   - If any of these is missing, raise an error naming the column; never print values.
2. **Harmonization fields**, `data/interim/encounter_harmonization_fields.csv` (reviewed version accepted by Vinai on Nov 10, 2026; SHA-256 `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`). The full SHA-256 must match this value exactly. A prefix match is not enough.
   - Use the harmonized **t0** (first INTIME) and **last OUTTIME** for each `encounter_id`.
   - Use these as given. **Do not rebuild the 24 h stay merge.**
3. **Death (inherited, primary):** `Died` from the existing pipeline label file (`data/processed/encounter_labels.csv`), the recorded ICU deaths (count omitted from this public copy). Do not redefine it.
   - `death_attribution_flag` (from the harmonization fields) feeds one sensitivity count only (see Outputs §7).
4. **Split:** the inherited split (`data/processed/subject_splits.csv`, or the split column already joined in the label file).
   - Expected encounter counts: the frozen subject-grouped 70/15/15 split (sizes omitted from this public copy); the frozen canonical cohort (counts omitted from this public copy) overall.
   - Raise an error on any mismatch, or if any subject appears in more than one split.
5. **Raw PIC, two tables only, read with `usecols`.** No other raw table is read.
   - `ICUSTAYS`: `SUBJECT_ID`, `HADM_ID`, `ICUSTAY_ID`, `INTIME`, `OUTTIME`. Used to find the stay that provides the last OUTTIME and to find ICU readmissions.
   - `ADMISSIONS`: `SUBJECT_ID`, `HADM_ID`, `DISCHTIME`. Used for hospital discharge time.
   - `DEATHTIME` is **not** read (see Limitations).
   - Apply the Stage 1 valid-stay rule to raw stays (drop `OUTTIME <= INTIME`; expect [N_valid_stays omitted] valid stays) and count anything it drops.

## Definitions
All times are in hours. Let L be the landmark, t0 the harmonized first INTIME, and `OUT` the harmonized last OUTTIME of the index encounter.

- **Encounter in-time and out-time (intentional):** the encounter in-time t0 is the earliest ICU stay INTIME in the encounter, and `OUT` is the last OUTTIME across all ICU stays in the encounter. Both come from the harmonization fields and are never changed by this audit.
- **Landmarks:** L ∈ {24, 48, 72, …, 336}, every 24 h (14 landmarks).
- **Core landmarks** (used for reporting and the merge flag): 24, 48, 72, 120, 168, 240, 336.
- **Eligible at L:** `OUT > t0 + L` (strictly).
- **Stacked rows:** one per (encounter, eligible landmark).

### Outcome 1: remaining ICU time (primary), levels 0–8
- Death takes precedence: `Died = 1` → **level 8**, whatever the remaining time.
- For survivors, `r = OUT − (t0 + L)`. Because the patient is eligible, r > 0.
- Levels and their exact boundaries:

| Level | Rule (survivors) |
|---|---|
| 0 | 0 < r ≤ 24 |
| 1 | 24 < r ≤ 48 |
| 2 | 48 < r ≤ 72 |
| 3 | 72 < r ≤ 120 |
| 4 | 120 < r ≤ 168 |
| 5 | 168 < r ≤ 240 |
| 6 | 240 < r ≤ 336 |
| 7 | r > 336 |

- **Clock:** r runs on the continuous encounter clock and includes within-encounter gaps of 24 h or less (Stage 3).
- **Readmissions:** for readmitted patients, r ends at the index encounter's `OUT`.

### Outcome 2: post-ICU course (key secondary), levels 0–4
Five levels, decided by Vinai on Nov 10, 2026 after the full band-support audit. The readmission timing sub-levels and the 14–30 d / > 30 d stay levels were merged; see the Changelog. Assign the level by the first rule that applies, in this order:

1. **Death:** `Died = 1` → **level 4**.
2. **Index HADM:** the `HADM_ID` of the raw valid stay whose OUTTIME equals `OUT`, matched within the encounter's `hadm_id_list`.
   - **Tie-break:** if the last stays share the same OUTTIME (`OUT`) but carry different HADM IDs, use the HADM with the latest DISCHTIME, and count those encounters (`index_hadm_ambiguous`).
   - The tie-break selects only two things: the DISCHTIME used for the post-ICU stay (Outcome 2), and the hospitalization searched for readmissions. **It does not change the encounter in-time (t0) or out-time (`OUT`).**
   - If no stay matches, count it (`index_hadm_unmatched`); the encounter's Outcome 2 is then unassignable.
3. **Readmission:** the **earliest** raw valid ICU stay of the same `SUBJECT_ID` with the same index `HADM_ID` and `INTIME > OUT`. Let `g` be that stay's INTIME minus `OUT`.
   - If `g ≤ 24`, count it in `return_within_24h`. A non-zero count is an error: report the count only, then stop. It should not exist after harmonization.
   - `g > 24` → **level 3** (ICU readmission more than 24 h after OUTTIME). This takes precedence over the stay levels 0–2.
   - The timing sub-bands (24 < g ≤ 48, 48 < g ≤ 72, g > 72) are **descriptive only**. They are reported in `readmission_check` and are not outcome levels.
4. **No readmission:** `p = DISCHTIME(index HADM) − OUT`, in days (hours / 24).
   - `p < 7` → **level 0**; `7 ≤ p < 14` → **level 1**; `p ≥ 14` → **level 2** (post-ICU stay over 14 d).
   - Boundaries: "7–14 d" is 7 to under 14 d. "Over 14 d" includes exactly 14 d, which is the union of the former 14–30 d and > 30 d levels under the v2 boundary convention.
5. **Unassignable (no readmission and no usable discharge time):** these encounters are reported on a **separate line**, `outcome2_unassignable`, broken down by reason (`dischtime_missing`, `dischtime_before_out`, `index_hadm_unmatched`). They are **never** silently dropped or forced into a level.
   - Identity (an acceptance check): assignable levels 0–4 + unassignable = total eligible, at every landmark × split.
   - Unassignables are **excluded** from the Outcome 2 level-support denominators (cutpoint support, weakest cutpoint) and from the merge-review flags.
   - The expected count is zero; the full audit found zero. If it is not zero, the count is reported, not raised as an error.
   - Sensitivity: report separately what the level counts would be if `dischtime_before_out` cases within 24 h were set to p = 0 (level 0). Earlier audits documented timestamp-rounding exceptions in which DISCHTIME falls before OUTTIME.

Outcome 2 does not depend on L except through eligibility. It is computed once per encounter and tabulated per landmark.

### Cutpoints
- Outcome 1: k = 0..7. Outcome 2: k = 0..3.
- **Support at k** is the pair `n(y ≤ k)` and `n(y > k)`. Every k is included, **k = 0 too**. There is no trivial cutpoint at the top level.
- Unassignable Outcome 2 rows are excluded from its support, its weakest cutpoint and its merge flags, and are reported on their own line next to it.
- **Weakest cutpoint** (per landmark × split × outcome): the k that minimizes `min(n(y ≤ k), n(y > k))`. Ties go to the lowest k. Report k, the smaller side, and which side it is.
- **No structural levels:** under these definitions every level is possible at every landmark. Nothing is hardcoded as structural or forced to 0. A level whose raw count is 0 is flagged (`zero_levels`) and suppressed like any other small cell.

## Run order: pilot first, then the full matrix
1. **Pilot run (required before the full run):**
   - Scope: core landmarks only (24, 48, 72, 120, 168, 240, 336 h), both outcomes, **train split only**.
   - Validate:
     - the boundary unit tests (acceptance check 12);
     - the readmission counts: `return_within_24h = 0`, the counts by readmission sub-level, and the `readmit_index_overlap` counts;
     - the raw reconciliation identities, including assignable + unassignable = eligible;
     - the train eligible count at each core landmark, which must be ≤ the overall counts in acceptance check 4 and non-increasing.
   - Output: the same aggregate-only JSON/markdown, with the suffix `_pilot`. Suppression applies as usual.
   - **If any check fails, stop and report the failing check (counts only). Do not start the full run.**
2. **Full run:** only after the pilot passes. All 14 landmarks × train / validation / test / overall, with every acceptance check.

## Outputs
Files:
- JSON: `data/interim/landmark_forward_band_support_audit.json`
- Markdown summary: `data/interim/landmark_forward_band_support_audit.md`

Suggested script: `scripts/audit_landmark_forward_band_support.py`. The top-level JSON carries the survivor-bias note (see that section) and the input hashes.

For each of the 14 landmarks × {train, validation, test, overall}:
1. **Eligible encounters.**
2. **Outcome 1 level counts** (0–8).
3. **Outcome 1 support** for k = 0..7 (both sides).
4. **Outcome 2 level counts** (0–4), with `outcome2_unassignable` (by reason) as a separate line.
5. **Outcome 2 support** for k = 0..3 (both sides). Internal only unless jointly suppressed; see acceptance check 21.
6. **Weakest cutpoint** for each outcome (k, smaller side).
7. **Deaths** (Outcome 1 level 8 / Outcome 2 level 4) as context. Also, as a one-line sensitivity count, deaths with `death_attribution_flag = 1` among the eligible.
8. **Merge-review flag**, at core landmarks only, test split:
   - `merge_review_needed = true` if any cutpoint has fewer than 30 on either side.
   - List the failing k and the side. The threshold is configurable; the amendment says "about 30".
   - This is a **flag for Vinai**. The script never merges levels.
   - The flag is computed from raw counts. Only the boolean and the k list are written.
9. **Per-row `reconciled`** boolean (see Suppression).

Cohort-level items (by split and overall):

10. **Readmission check:**
    - `return_within_24h` (must be 0).
    - Readmissions in the same HADM (= Outcome 2 level 3), with the descriptive timing sub-bands (24–48 / 48–72 / > 72 h).
    - Readmissions in a different HADM. These are informational and not used in Outcome 2.
11. **Readmission overlap:** the counts defined under "Named decision: `readmit_index_overlap`" below.
12. **Data-quality counts**, reported as counts and never silently dropped:
    - Index HADMs with missing DISCHTIME.
    - Index HADMs with DISCHTIME before `OUT`, split into ≤ 1 h, 1–24 h, and > 24 h before.
    - Raw ICU stays with a blank `HADM_ID`: overall, among cohort subjects, and among stays after an index `OUT`. These stays cannot be judged as same-HADM readmissions.
    - `index_hadm_ambiguous` and `index_hadm_unmatched`.
    - Raw stays dropped by the valid-stay rule.
    - Encounters with more than one HADM in `hadm_id_list` (expected 52).
    - Encounters with an ICU span (`OUT` − t0) over 60 days (> 1,440 h), by split and overall, suppressed if < 10. The reviewed file has a maximum span of 28,003.70 h, which needs checking; the count is required in the full run.
13. **Stacked rows:**
    - Total stacked rows by split and overall (sum of eligible over the 14 landmarks).
    - Rows per patient (sum over the patient's encounters), as median and IQR by split. Quantiles only; never min or max.
14. **Markdown summary:**
    - One table per outcome: rows = core landmarks; columns = eligible (test), weakest test cutpoint (k, smaller side), `merge_review_needed`.
    - The readmission and data-quality counts.
    - The stacked-row table.
    - The survivor-bias note.

## Named decision: `readmit_index_overlap`
- **Decision variable:** `readmit_index_overlap = keep_both` (default; Vinai, Nov 10, 2026). An encounter whose readmission stay is itself an index encounter in the frozen cohort keeps its Outcome 2 readmission level, and the readmission encounter also keeps its own rows in the stack.
- **Rationale:**
  - The split is by patient, so both encounters land in the same split. There is no leakage across splits.
  - The cluster bootstrap resamples by patient, so the dependence between the two encounters is carried into the confidence intervals.
  - Layer A prior-admission and recency features use only information from before each encounter's own INTIME, so the readmission encounter's features never see the index encounter's outcome.
- **What the audit reports** (by split and overall; counts suppressed when < 10, like every other cell):
  - The number of encounters whose readmission stay is itself an index encounter. Match on the same `SUBJECT_ID`, with the readmission INTIME falling in [t0, `OUT`] of another cohort encounter.
  - Readmission stays not found in the cohort.
  - The number of overlapping pairs in which **both** encounters are eligible at each core landmark (both contribute rows to the stack).
  - A check that both encounters of every pair are in the same split.
- **Leakage note:** keeping both is safe only if Layer A recency features are computed strictly from data before each encounter's INTIME. This is recorded as acceptance check 17 and must also be enforced when Layer A is built.

## Late-censoring and survivor bias
Eligibility at L means still being in the ICU at L, so later landmarks describe patients who have survived with an ongoing stay. They are not the full cohort. PIC stays are complete, so there is no censoring. The JSON carries this note at the top level.

## Suppression and disclosure (PhysioNet DUA)
- **Aggregate counts only.** No row-level output, no IDs, no example rows.
- **Raw-count reconciliation, done before suppression and not written out.** For every landmark × split × outcome, raise an error unless all of the following hold:
  - level counts (Outcome 2: plus unassignable) sum to the eligible count;
  - `n(y ≤ k) + n(y > k)` equals the eligible count (Outcome 2: the assignable count) for every k;
  - `n(y ≤ k)` equals the cumulative sum of level counts;
  - train + validation + test equals overall for every cell;
  - Outcome 1 level 8 equals Outcome 2 level 4.
- **In plain language:**
  - *Primary suppression* hides any small count, from 1 to 9, so that no small group of patients can be singled out.
  - *Secondary (complementary) suppression* hides some additional, larger counts. Without it, a hidden small count could be worked out by subtracting the published cells from a published total. For example, if the levels of a row sum to the published eligible count and only one level is hidden, that level equals the total minus the others.
- **Mandatory: raw-count reconciliation** (above), **primary suppression**, and **secondary suppression**. The audit is not valid without all three.
- **Primary suppression (mandatory):** any count from 1 to 9 becomes `"<10"`. Zero is also shown as `"<10"`: nothing is structural, so no zero is published as an exact 0. The design-relevant zeros are listed in `zero_levels`, which records only that the raw count was 0.
- **Secondary / complementary suppression (mandatory):** repeat until nothing changes:
  1. *Margins.* Within each published set that sums to a published total, if exactly one cell is suppressed, also suppress the next smallest unsuppressed cell in that set. The sets are:
     - levels within a split, which sum to eligible;
     - splits within a level, which sum to overall;
     - the pair `n(y ≤ k)` and `n(y > k)`, which sums to eligible (Outcome 2: assignable).
  2. *Cumulative differencing.* `level k = n(y ≤ k) − n(y ≤ k−1)`, so if level k is suppressed, the support pairs at k−1 and k must not both be published. Suppress the pair whose smaller side is smaller.
  3. *Pairs.* If either side of a support pair is suppressed, suppress both.
- **Optional defensive final verification (not the main logic):** after steps 1–3, the script may also check every suppressed cell with a linear program (min and max of the cell subject to all published cells, the constraints above and non-negativity, e.g. `scipy.optimize.linprog`). If a cell is exactly recoverable (min = max) or its feasible interval lies inside [1, 9], suppress more and repeat; if it still fails, raise an error. Steps 1–3 must be correct on their own; this check is only a safety net.
- Suppressed rows will not sum, by design. Each row carries `reconciled: true`, confirming that the raw counts reconciled.
- **Logs:** error messages and logs must not print subject, encounter or admission IDs, or timestamps of individual rows. They print counts only.

## Acceptance checks (the script must raise on failure)
1. Canonical SHA-256 mismatch, or a harmonization-fields SHA-256 that is not exactly `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`. Both checks compare the full 64-character hash and run before anything else.
2. Missing required columns in any input, named by column without values.
3. Cohort totals other than [N_encounters omitted] encounters / [N_subjects omitted] subjects / the frozen 70/15/15 split (sizes omitted from this public copy); or any subject in two splits.
4. Eligible counts that don't reproduce the earlier eligibility audit. The new code must reproduce [landmark eligible counts omitted] overall, and 1,406 / 1,097 / 947 / 723 / 601 / 473 / 355 in test, at 24 / 48 / 72 / 120 / 168 / 240 / 336 h. Eligibility must also be monotonically non-increasing across all 14 landmarks.
5. Deaths among the eligible that don't reproduce 767 / 644 / 564 / 450 / 379 / 296 / 235 at the core landmarks.
6. `return_within_24h > 0`.
7. Any raw reconciliation failure (see Suppression).
8. Any eligible survivor with r ≤ 0.
9. Any cutpoint list that doesn't start at k = 0 and end at the second-highest level. A unit test asserts k = 0..7 and 0..3.
10. Any level forced to 0, or any hardcoded structural-level list. A unit test on synthetic data with a zero-count level must show it in `zero_levels` and suppressed as `"<10"`, not dropped.
11. Secondary suppression not applied, or (if the optional LP verification is enabled) a verification failure.
12. Boundary unit tests on synthetic rows:
    - r = 24 → level 0; 24.0001 → level 1; 336 → level 6; 336.0001 → level 7.
    - g = 24 → `return_within_24h` (stop); 24.0001 → level 3; 48, 72 and 72.0001 → level 3 (the sub-band is recorded descriptively).
    - p = 6.999 d → level 0; 7 d → level 1; 13.999 d → level 1; 14 d → level 2; 30.0001 d → level 2.
    - Death overrides readmission and discharge time.
13. Readmission overlap pairs found in different splits.
14. Any raw table read beyond the ICUSTAYS and ADMISSIONS columns listed above.
15. Outcome 2 identity: assignable levels 0–4 + `outcome2_unassignable` ≠ total eligible at any landmark × split. A non-zero unassignable count is reported, not raised. Unassignables must not appear in any support denominator or merge flag (unit test).
16. Full run started without a passing pilot run (the full run checks for the pilot's pass record).
17. Leakage guard for `readmit_index_overlap = keep_both`: Layer A prior-admission and recency features use only information before each encounter's own INTIME. The audit records this as a required check for the Layer A build, and verifies here that every overlap pair is in one split (check 13).
18. Readmission-level identity, at every landmark × split:
    - Outcome 2 level 3 equals the same-HADM readmission count in `readmission_check`, and the descriptive timing sub-bands sum to it.
    - Readmission (level 3) takes precedence over post-ICU stay levels 0–2. A readmitted survivor is never dropped or given a stay level.
    - With check 15, levels 0–4 + unassignable = eligible.
    - A unit test on synthetic data includes readmitted survivors and asserts they land in level 3.
19. Per-landmark counter reset:
    - Every `readmission_check` and `data_quality` counter is computed fresh for each landmark, from that landmark's eligible set only, and never accumulates across landmarks.
    - Counts restricted to the eligible set must not increase from one landmark to the next (eligible sets are nested).
    - Encounter-level counts that do not depend on L (e.g., multiple HADMs among the eligible) must not exceed the eligible count.
    - A unit test runs two landmarks on synthetic data and asserts that the second landmark's counters are not the sum of both.
20. **Cross-split complementary suppression.** Overall is published next to train / validation / test, so overall − (the other splits) must not recover a suppressed split cell. In every set tied by a published total (levels within a split; train + validation + test = overall for each level, each eligible count and each count item), the number of suppressed cells is 0 or at least 2. The full audit found about 145 cells recoverable this way or by check 21's differencing.
21. **Cumulative support counts.** Cutpoint support pairs `n(y ≤ k)`, `n(y > k)` can recover suppressed level cells by differencing. Shared files must do one of the following:
    - (default) omit the support counts, keeping only the weakest cutpoint k, the `merge_review_needed` flag and the failing-k list, with the counts kept in the local raw file; or
    - suppress them jointly with the level cells, so that no `n(y ≤ k) − n(y ≤ k−1)` or `eligible − n(y ≤ k)` reveals a suppressed level.
22. **Recoverability self-check (mandatory before any shared file is written).** For every suppressed cell, apply every published linear identity: level sums, split sums, support pairs and cumulative differences, Outcome 1 level 8 = Outcome 2 level 4 = deaths, and readmission level 3 = the readmission count. Raise an error if any suppressed cell is exactly recoverable or bounded inside [1, 9]. The JSON records the number of suppressed cells checked and the number recoverable, which must be 0. The LP check remains an optional extra.
23. **ICU span count present.** The JSON must contain the count of encounters with an ICU span over 1,440 h (60 d), by split and overall, suppressed if < 10.
24. **Overlap count present.** The JSON must contain the `readmit_index_overlap` count: encounters whose readmission stay is itself an index encounter, by split and overall, suppressed if < 10.
- **Release gate:** no shared file is copied, uploaded or shared until checks 20–24 pass. Unsuppressed files stay local.

## Post-audit decisions (Vinai; not automated)
- Level merges per landmark. The merge rule: merge a level with its neighbor if any test-split cutpoint has fewer than about 30 on either side at a core landmark. Merges are prespecified and applied identically across models.
  - Outcome 1: no merges. In the full audit, all 9 levels have at least 12 in every split and landmark.
  - Outcome 2: decided as 5 levels (v2.3). The merge flag stays informational.
- Readmission stays that are also index encounters: default `readmit_index_overlap = keep_both` (see the named decision). Vinai may change it after seeing the overlap counts.
- Whether `dischtime_before_out` cases within 24 h are set to p = 0.

## Limitations
- `DEATHTIME` is not read. "Death any time after L" relies on the inherited encounter-level `Died`. The audit cannot confirm that no eligible death occurred before t0 + L, though an eligible patient was still in the ICU at t0 + L. The JSON states this.
- "Same hospitalization" means the HADM of the stay that ends at `OUT`. For the 52 multi-HADM encounters this is a definition, not a data fact.

## Do not
- Do not modify the frozen CSV, the harmonization fields file, the labels, or the split file.
- Do not rebuild the stay merge or the outcome labels of 3.0.
- Do not build features, train models, or use the test split for anything beyond these counts.
- Do not print or save any row-level data.

## Changelog
- **v2.3 (Nov 10, 2026):** after the full audit, Vinai set Outcome 2 to 5 levels (4 death; 3 ICU readmission > 24 h; 2 post-ICU stay over 14 d; 1 7–14 d; 0 under 7 d). Readmission timing is kept as descriptive sub-bands. Cutpoints are now k = 0..3, and the boundary tests and checks 9, 12, 15 and 18 were updated. Added acceptance items 20–24 and a release gate. The full audit's shared outputs had about 145 suppressed cells recoverable by cumulative differencing or overall − splits, and lacked the > 1,440 h span count and the overlap count.
- **v2.2 (Nov 10, 2026):** after the pilot review, added acceptance checks 18 (readmission-level identity and precedence) and 19 (per-landmark counter reset, and subset counters that don't increase across landmarks). The pilot found readmitted survivors dropped from Outcome 2 and counters accumulating across landmarks.
- **v2.1 (Nov 10, 2026):** the harmonization-fields reference changed from the `4b083d18` prefix to the full reviewed SHA-256 `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`, which Vinai accepted after a reproduction check. The gate is now an exact full-hash match. Added the count of encounters with an ICU span over 60 days to the data-quality counts.
- **v2 (Nov 10, 2026):** five clarifying edits from Copilot review and Vinai's answers.
  1. **Index HADM tie-break:** encounter in-time = earliest stay INTIME, out-time = last OUTTIME across all stays. When the last stays share the OUTTIME but have different HADMs, the HADM with the latest DISCHTIME is used, for Outcome 2's DISCHTIME and the readmission search only.
  2. **Unassignable bookkeeping:** Outcome 2 unassignables are a separate line; assignable + unassignable = eligible is an acceptance check; unassignables are excluded from support denominators and merge flags; zero is expected.
  3. **Suppression:** primary and secondary (complementary) suppression and raw-count reconciliation are mandatory; the LP recoverability check is now optional defensive verification; a plain-language explanation was added.
  4. **Pilot run:** required before the full matrix (core landmarks, both outcomes, train only); stop and report if a check fails.
  5. **Named decision `readmit_index_overlap = keep_both`:** rationale and leakage check added; the readmission-overlap counts moved under it.
  - Acceptance checks 15–17 added.
- **v1 (Nov 10, 2026):** initial spec.
