# Copilot spec: Landmark eligibility audit (read-only, aggregate-only)
Version 2, Nov 10, 2026 (incorporates Vinai's review comments)

## Purpose
Decide which landmarks are defensible for the composite ordinal-8 endpoint, while preserving the frozen cohort and the read-only audit guardrails. Candidates: 24, 48, 72, 120, 168, 240 and 336 h. Also decide whether Layer C should extend to 336 h (6 h bins, 56 bins). This is a design audit only. It creates no features and changes no frozen files.

## Primary decision target
- **The primary landmark viability criterion is support for the ordinal-8 outcome, not mortality alone.**
- The ordinal-8 outcome captures both death and the burden of ICU and hospital stay, so a landmark can be clinically useful even when its mortality rate is low.
- Mortality counts are a secondary sanity check, not the only gate.

## Inputs (read only; do not modify)
- Frozen canonical encounter CSV (SHA-256 `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418`). Check the hash before reading and raise an error on mismatch.
- `data/interim/encounter_harmonization_fields.csv` (t0, ICU span, `icu_time_hours` / `sum_icu_stay_hours`, `death_attribution_flag`).
- The existing ordinal-8 label and split columns as they are already defined in the pipeline. Do not redefine the outcome. If the ordinal-8 column can't be found, raise an error; do not rebuild it.

## Definitions
- Landmarks L ∈ {24, 48, 72, 120, 168, 240, 336} hours after t0.
- **Eligible at L:** the patient is still in the ICU at t0 + L, meaning the encounter's last OUTTIME is strictly after t0 + L, and the patient is alive at t0 + L. If time of death is available, exclude deaths at or before t0 + L. If it isn't, state that in the output.
- Bins: 6 h from t0. Landmark L corresponds to bin L/6: 4, 8, 12, 20, 28, 40, 56.
- Keep the existing 24/48/72 h eligibility logic. The new code must reproduce [N_eligible omitted], [N_eligible omitted] and [N_eligible omitted] exactly, or raise an error.

## Outputs: `data/interim/landmark_eligibility_audit.json`
For each landmark L:
1. **Cohort size:** eligible encounters, overall and by split (train, validation, test).
2. **Ordinal-8 distribution (primary):** counts of each level 1–8, overall and by split.
3. **Support by threshold (primary):** for each cutpoint k = 1..7, report count(y ≤ k) and count(y > k), overall and by split. A landmark can have an acceptable total but poor support at one end of the scale, and this table shows it.
4. **Structural impossibility:** list every ordinal level that cannot occur at landmark L by construction, with a one-sentence reason each. Label these "structural," not as a data issue or a small cell.
   - Example: the ordinal-8 outcome is partly built from ICU and hospital length of stay. A patient still in the ICU at L has an ICU stay longer than L. Any level that requires an ICU stay of L or less (or a hospital stay shorter than that) is therefore impossible at that landmark. At 24 h this may rule out the shortest-stay level; at 336 h it may rule out every short- and medium-stay level.
   - Derive impossibility from the outcome's definition (thresholds and rules), not from observed zeros. A level that is possible but has zero observed cases is a small cell, not structural.
5. **Mortality (secondary sanity check):** in-hospital deaths among the eligible, overall and by split.
6. **Composite summary:** a one-line summary per landmark giving eligible count, test-split minimum threshold support (the smaller side of the weakest non-structural cutpoint), and death count, in that order. This keeps the ordinal-8 support as the headline and mortality as context.
7. **Layer C bin occupancy:** for bins 1–56, the number of eligible-at-24 h encounters still in the ICU at the end of each bin. One row per bin, overall only.
8. **Remaining ICU time:** median and IQR of remaining ICU time after L among the eligible, in hours.

## Late-censoring and survivor bias
- The eligible cohort at L is defined by still being in the ICU at L. Later landmarks are therefore **conditional on survival and an ongoing ICU stay.**
- Results at 120, 168, 240 and 336 h inform the design for long-stay patients. They are not population-wide estimates.
- The JSON should carry this note at the top level.

## Candidate landmark selection rule (applied by Vinai after the audit; not automated)
- **Default core set: 24, 48 and 72 h.** These are the primary design landmarks.
- **120, 168, 240 and 336 h are sensitivity or late-trajectory candidates**, not primary design landmarks.
- Keep a landmark based on:
  1. stable cohort size,
  2. stable ordinal-8 support across splits,
  3. adequate coverage of non-structural thresholds, and
  4. minimal survivor-bias distortion.
- **A landmark is weak if any non-structural threshold in the test split has fewer than about 30 cases on either side.**
- If a level is sparse at a given landmark, consider merging it with its neighbor for that landmark before dropping the landmark.
- Layer C decomposition is fit per landmark, on that landmark's eligible encounters.

## Suppression and disclosure (PhysioNet DUA)
- Aggregate counts only. No row-level output, no IDs, no example rows.
- **Internal consistency check (before suppression, not written to output):** for every landmark and split, the ordinal level counts must sum to the eligible count, and count(y ≤ k) + count(y > k) must equal the eligible count for every k. Raise an error if not.
- **Then suppress every row of the ordinal distribution and threshold report before any comparison or output:**
  - Report any count from 1 to 9 as `"<10"`.
  - Show 0 only for structurally impossible levels; any other zero is reported as `"<10"`.
  - Apply secondary suppression: if a suppressed cell could be recovered by subtraction from a published total (for example train + validation + test = overall, or level counts summing to the eligible count, or y ≤ k vs y > k), also suppress the next smallest cell in that row or column.
  - Suppressed rows will not sum to the eligible count, by design. Instead, output a boolean per row confirming that the unsuppressed counts reconciled.
- Error messages must not print subject, encounter or admission IDs.

## Checks the script must raise on
- Canonical hash mismatch.
- Missing ordinal-8, split, t0 or OUTTIME columns.
- Any landmark count that is not monotonically non-increasing from 24 to 336 h.
- 24/48/72 h counts that differ from [N_eligible omitted] / [N_eligible omitted] / [N_eligible omitted].
- Any internal consistency failure described above.

## Do not
- Do not modify the frozen CSV or the harmonization fields file.
- Do not build features, train models or touch the test-set modeling outcome beyond these counts.
- Do not print or save any row-level data.
