# Copilot spec: Surgery complexity vs duration audit (read-only, aggregate-only)
Version 1.7, Nov 11, 2026. **Surgery timing is now by anesthesia end** (Vinai's decision, Nov 11, 2026; fix list Amendment C): Layer A recent = anesthesia end in [t0 − 24 h, t0 + 6 h]; Layer C = after t0 + 6 h up to `OUT`, merge gaps included; remote = before t0 − 24 h. v1.6: **the surgery key is `ROW_ID`** (Vinai's decision, Nov 11, 2026: `OPER_ID` is not unique within an admission; all 7,488 rows are kept; see `surgery_complexity_duration_audit_FIXES_v1.md`, Amendment B). v1.5: **D2 uses the signed reading** (Vinai's decision: D2 is met if τ-b < 0.2, signed). v1.4 clarified the spec after the first run's code review (see `surgery_complexity_duration_audit_FIXES_v1.md`, Amendment A). The map and hashes are unchanged. **The final reviewed map's hash is pinned** (map v1, reviewed by Vinai Nov 10, 2026; see Changelog). The discrepancy criteria D1–D3 were approved by Vinai on Nov 10, 2026.

**Post-run documentation update (Nov 11, 2026). The spec version stays at v1.7, and the script and its hash are unchanged.**
- The MAP_v1 run on the v1.7 script passed its release gate. Results are in "Results: MAP_v1 run (Nov 11, 2026)".
- Post-audit decisions 3 and 4 are **decided**: φ1 = surgery type × type-specific duration band.
- **φ1 band defaults are approved (Vinai, Nov 11, 2026; decision 16 follow-up):**
  - ratio to the type's train median, with cuts < 0.75 / 0.75–1.33 / > 1.33;
  - minimum n = 10 train surgeries per type;
  - fallback type → system × complexity group → system → all.
- **Band-support audit done and medians frozen (Nov 12, 2026):**
  - freeze at min n = 10, with frozen-medians CSV SHA-256 `e73c8b5a…0047`;
  - the band is a 3-level ordinal alongside the type and system descriptors;
  - combination-name collapsing is a sensitivity analysis only.
- Two minor items stay open; see "Open items for φ1 bands".

**Status: run complete (Nov 11, 2026).** The final run was authorized by Vinai on Nov 10, 2026.
- The complexity map is the **final reviewed map v1**, `surgery_complexity_map_v1.csv`, with SHA-256 `5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5`.
- Vinai reviewed the type map by system and the Cardiac_STAT sheet and accepted them as-is, with no edits. Flagged items are accepted at their default assignment.
- Every output is labelled `MAP_v1` (reviewed 2026-11-10).

Addendum to `specs/surgery_item_level_audit_spec.md` (v1.2). It uses that spec's inputs, encounter linkage, timing classes, windows, split, suppression and release gate unchanged, and adds one question: **do the planned complexity groups agree with actual surgery duration?** For cardiac surgery the groups are STAT groups. For every other body system they are the type map's `intensity` levels. Style follows v1.2.

## Purpose
Vinai plans to use both a complexity grouping and the duration of surgery (φ1 = type × duration band). This audit measures, per body system and pooled, how far actual duration agrees with the complexity grouping:

- if they agree closely, duration adds little beyond the group;
- if they disagree, duration carries separate information, or the grouping needs review.

It is descriptive only. **The audit sets no bands, merges no groups and changes no map.** Disagreement is judged against pre-specified criteria that Vinai **approved on Nov 10, 2026** (see "Discrepancy criteria"). STAT is a mortality-risk score, not a duration score, so a discrepancy is a finding to interpret, not an error.

## Inputs (read only; do not modify)
1. **Everything in v1.2 Inputs §1–4**, with the same exact-hash checks: the canonical CSV `e5feeb22…6418` and harmonization `ffd1b8ae…209d6d` (full 64-character match), the split, and the `SURGERY_INFO.csv` headers and row count (7,488). `SURGERY_VITAL_SIGNS` is **not** needed.
   - **(v1.6) Surgery key:** `SURGERY_INFO.ROW_ID` (normalized like the other IDs: trimmed, `123.0` → `123`). The script asserts it is nonblank and unique on all 7,488 rows and stops otherwise. `OPER_ID` is kept as a descriptor only. Real-data diagnostic (Copilot, Nov 11, aggregate): 889 rows fall in 389 (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) groups of size ≥ 2, with 0 exact copies. Every group differs in anesthesia and surgery times, so these are distinct surgeries.
   - **(v1.6) Harmonization columns:** t0 is read from `t0`. `OUT` is read from a column named `OUT` (or one of the listed aliases) when present; otherwise from `OUTTIME`, which is how the real `encounter_harmonization_fields.csv` names it. `--harm-t0-col` / `--harm-out-col` override both. When `span_icu_hours` is present, the script counts encounters with |`OUT` − t0 − `span_icu_hours`| > 1 min. This confirms that the column holds the **last** (max) OUTTIME. The count is reported and does not stop the run; any nonzero count goes to Vinai.
2. **Complexity map** (name-level only; no counts; not data). Copy it to `config/surgery_complexity_map_v1.csv`.
   - **Pinned (final, reviewed):** map v1, with `COMPLEXITY_MAP_SHA256 = "5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5"`.
   - It was built from `surgery_type_map_DRAFT.csv` and `cardiac_stat_names.csv` as Vinai reviewed and approved them on Nov 10, 2026, with no reviewer edits. The `reviewer_*` columns are empty, so every name keeps its default assignment, including the review-flagged names.
   - It has 1,495 rows, one per normalized raw name, covering every body system. Columns:
     - `raw_name`: the join key, already normalized (trim, case-fold, collapse whitespace);
     - `standardized_name`, `body_system`;
     - `intensity` (the 4 labels, or `NA` for Unclassified);
     - `I_group` (I1/I2/I3; `NA` for Unclassified);
     - `STAT` (1–5 or `NC`) and `STAT_group` (G1/G2/G3/NC), for cardiac rows only (blank otherwise);
     - `complexity_group`: G1/G2/G3/NC for cardiac and I1/I2/I3 otherwise. It is `NA` for the 2 Unclassified names, which are excluded and counted;
     - `review_flag` (Y/N). For cardiac rows it is Y if either the type-map flag or the STAT flag is set;
     - `n_components` (for sensitivity analysis a);
     - `STAT_group_crosses` (Y/N; cardiac only; for sensitivity analysis b).
   - **On a hash mismatch, the script raises an error and stops.** The JSON records `map_status = "reviewed_v1"`, `map_label = "MAP_v1"`, `map_reviewed_date = "2026-11-10"` and the full hash.
   - Any later change to the map needs a new map version, a new hash and a new spec version.
   - **Copilot does not assign or change any group, STAT category or intensity.**
3. **Join:** `SURGERY_INFO.SURGERY_NAME` is normalized (trim, case-fold, collapse whitespace; the same rule as the type map) and joined to `raw_name`. Map keys must be unique (asserted). Expected coverage is all 7,488 rows, since the map covers all 1,495 normalized names. The 1,496 raw names in `SURGERY_INFO` collapse to 1,495 after normalization: one pair differs only by case. Unmatched rows are counted (expected 0) and reported as `unmapped`.

### Groups actually present (context; name-list counts, not cohort output)
These counts come from Vinai's name list (no identifiers), before cohort linkage. Cells of 1–9 are shown as `<10`, and "–" means no names; `S` marks complementary suppression (GI/abdominal minor/diagnostic, which would otherwise reveal the catheter cell from the I1 total below).

**Intensity labels in the type map: four, not three:** `minor/diagnostic`, `catheter/interventional`, `minimally invasive` (laparoscopic, thoracoscopic or endoscopic) and `major open`. Two names (2 surgeries) are `Unclassified` with a blank intensity and are excluded. Surgeries per system × intensity:

| body_system | minor/diagnostic | catheter/interventional | minimally invasive | major open |
|---|---|---|---|---|
| Cardiac | – | 38 | – | 2,289 |
| Craniofacial/ENT/oral | – | – | <10 | 223 |
| Diagnostic/imaging-interventional | 283 | 155 | – | – |
| GI/abdominal | S | <10 | 491 | 1,293 |
| Hepatobiliary | – | <10 | 171 | 107 |
| Neurosurgical | – | – | 12 | 1,284 |
| Orthopedic/spine | <10 | – | – | 26 |
| Other | 154 | – | <10 | 10 |
| Thoracic/airway | – | <10 | 37 | 491 |
| Urologic/renal | <10 | <10 | 77 | 166 |
| Vascular access/line placement | 153 | – | – | – |

**Decided three-group collapse for non-cardiac systems** (Vinai, Nov 10, 2026): I1 = minor/diagnostic + catheter/interventional; I2 = minimally invasive; I3 = major open. Surgeries per system × I-group, and whether the system is **comparable** (at least 2 groups with ≥ 10 surgeries on the name list):

| body_system | I1 | I2 | I3 | comparable (eligible groups) |
|---|---|---|---|---|
| Craniofacial/ENT/oral | – | <10 | 223 | N |
| Diagnostic/imaging-interventional | 438 | – | – | N (single group) |
| GI/abdominal | 13 | 491 | 1,293 | Y (I1/I2/I3; I1 borderline) |
| Hepatobiliary | <10 | 171 | 107 | Y (I2/I3) |
| Neurosurgical | – | 12 | 1,284 | Y (I2/I3; I2 borderline) |
| Orthopedic/spine | <10 | – | 26 | N |
| Other | 154 | <10 | 10 | Y (I1/I3; I3 borderline) |
| Thoracic/airway | <10 | 37 | 491 | Y (I2/I3) |
| Urologic/renal | <10 | 77 | 166 | Y (I2/I3) |
| Vascular access/line placement | 153 | – | – | N (single group) |

Cardiac is grouped by STAT (below) and is comparable on G1/G2/G3. Cohort linkage and splitting reduce every count, so the borderline groups (13, 12, 10) may become ineligible at run time.

**Cardiac STAT groups** (map v1, reviewed; unchanged from the draft): G1 (STAT 1–2) 2,184; G2 (STAT 3) 36; G3 (STAT 4–5) 58; NC (not STAT-classified) 49. Of the 179 review-flagged cardiac names, 48 (207 surgeries) have a plausible STAT range that crosses a group boundary.

## Definitions
- **Analysis population (primary):** `SURGERY_INFO` rows linked to a cohort encounter by the v1.2 encounter link, in timing classes `recent_layer_a`, `remote_same_hadm` or `during_icu_layer_c` (v1.7 names). A surgery's split is its encounter's split.
  - **Timing time T (v1.7):** `ANES_END_TIME`. If it is missing, `SURGERY_END_TIME` is used and counted (`timing_fallback_surgery_end`); if both are missing the row is `invalid_time`. `SURGERY_BEGIN_TIME` no longer decides timing. A missing B makes only the duration invalid.
  - **One surgery, one encounter (v1.4).** Each surgery (one `SURGERY_INFO` row, key `ROW_ID`; v1.6) enters the analysis **at most once**. Rows that share (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) are separate surgeries and each one is linked on its own. Candidate encounters are the cohort encounters of the same subject whose `hadm_id_list` contains the surgery's `HADM_ID`. Times are the harmonized t0 and `OUT`, not raw `INTIME`/`OUTTIME`. Assignment rule, applied in order:
    1. If t0 − 24 h ≤ T ≤ `OUT` for a candidate, the surgery belongs to that encounter:
       - `recent_layer_a` if T ≤ t0 + 6 h. These surgeries are **Layer A only** and are never also Layer C events.
       - `during_icu_layer_c` if t0 + 6 h < T ≤ `OUT`. This **includes surgeries in merge gaps** between merged stays; there is no separate gap class.
       - Two encounters' windows cannot overlap, because encounters are more than 24 h apart.
    2. Otherwise it is assigned to the **later encounter**, the first candidate with t0 > T. It is then `remote_same_hadm` (T < t0 − 24 h). This covers surgery between encounters.
    3. Otherwise (T after the last candidate's `OUT`) it is `post_icu_same_hadm`: counted, and not analysed.
    - **Also reported:**
      - `recent_ending_after_t0`: recent surgeries whose T is in (t0, t0 + 6 h];
      - `straddling_t0`: B < t0 ≤ E;
      - `pre_admit_time`: remote surgeries with T before `ADMITTIME`. These are kept and flagged.
    - The B = t0 edge question from v1.2 is moot.
    - **Report:** `n_multi_candidate`, the number of surgeries with 2 or more candidate encounters, and how many of them were resolved by each rule.
    - **Note:** this rule applies only to this audit's unit of analysis (one surgery = one duration). **For features (approved Nov 11, 2026), one surgery may be a Layer C event for one encounter and Layer A history for a later encounter in the same admission.** That reuse is out of scope here.
  - **Reconciliation (v1.4; v1.6 wording):** linked (unique `ROW_ID`s) + `post_icu_same_hadm` + `other_hadm` + `not_in_cohort` + `invalid_time` + `not_linkable_missing_encounter_time` = 7,488. The last class covers rows whose candidate encounter lacks a usable t0/`OUT`; it is expected to be 0 and was already in the reference script. The script raises an error if linked exceeds 7,488 or if any `ROW_ID` appears twice in the analysis set.
  - **Key diagnostics (v1.6; reported, never dropped):** `reconciliation.surgery_key_ROW_ID` publishes four counts, suppressed like every other count:
    - `shared_oper_id_groups` and `rows_in_shared_oper_id_groups`: (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) groups of size ≥ 2;
    - `possible_true_duplicate_groups` and `rows_in_possible_true_duplicate_groups`: rows sharing `SUBJECT_ID`, `HADM_ID`, `SURGERY_BEGIN_TIME` and the normalized `SURGERY_NAME`.
    Group sizes go to the LOCAL file only. Possible true duplicates stay in the analysis. Whether to collapse them is a decision for Vinai after the run.
  - **Descriptive secondary:** all 7,488 `SURGERY_INFO` rows with a valid duration, pooled (no split), with no cohort restriction.
- **Duration** D = `SURGERY_END_TIME` − `SURGERY_BEGIN_TIME` (the v1.2 surgery window [B, E]), in hours.
  - **Fallback:** if B or E is missing or unparseable (expected 0), use the anesthesia window AE − AS, and count it (`duration_fallback_anes`).
  - **Invalid:** B > E or B = E (ordering violations, v1.2), or no usable window. These are counted and excluded, never silently dropped.
  - D > 24 h is counted and kept in the primary analysis; it is excluded in a sensitivity analysis.
- **Bands (v1.2, unchanged):** under 2 h (D < 2), 2–4 h (2 ≤ D ≤ 4), over 4 h (D > 4).
- **`complexity_group`** (ordinal low → high, from the map):
  - **Cardiac:** `STAT_group` G1 = STAT 1–2, G2 = STAT 3, G3 = STAT 4–5. `NC` is reported descriptively but is **not ordinal**, so it is excluded from correlations, monotonicity and overlap. Cardiac rows use STAT, not intensity. Cardiac intensity (catheter vs major open) is reported as context only.
  - **Every other system:** intensity collapsed to three groups (**decided by Vinai, Nov 10, 2026**): I1 = `minor/diagnostic` + `catheter/interventional`; I2 = `minimally invasive`; I3 = `major open`.
    - **Sensitivity:** the four intensity levels as an ordinal scale, minor/diagnostic < catheter/interventional < minimally invasive < major open (kept by Vinai as the sensitivity version).
- **Groups below 10 are reported, not dropped (v1.4).** Every group that exists in the map for a system appears in every split block (train, validation, test, pooled). When n < 10, its n and statistics are published as `"<10"`, including n = 0, as in the earlier specs. Cardiac `NC` is reported descriptively, with the same rule. The Unclassified surgeries are reported only as the `excluded_unclassified` count, not as a system.
- **Group eligibility:** a group enters a within-system comparison only if it has **at least 10 surgeries with a valid D** in the population being analysed.
  - A system is **assessable** if at least 2 groups are eligible. Otherwise it is reported as `not_assessable`, with the reason (`single_group` or `groups_below_10`).
  - Expected not assessable from the name list (I-group table above): Diagnostic/imaging-interventional and Vascular access (single group, I1), Orthopedic/spine (27 surgeries; one eligible group), Craniofacial/ENT/oral (I2 < 10). Borderline: GI/abdominal I1 (13), Neurosurgical I2 (12), Other I3 (10). Cohort linkage and splitting will shrink every count further, so eligibility is decided at run time.
  - **Diagnostic/imaging-interventional** is comparable only in the four-level sensitivity analysis (minor/diagnostic 283 vs catheter/interventional 155 on the name list).
- **Reference statistics (thresholds)** are computed on **train only**: group medians and IQRs used by the overlap metrics and by the discordance flag. Validation and test never set anything.

## Small cells and reporting level
- Cardiac G2 (~36) and G3 (~58) are small before linkage and splitting. The same holds for several intensity groups (e.g., GI/abdominal I1 = 13, Neurosurgical I2 = 12, Other I3 = 10, Orthopedic/spine 27 in total). Split-level cells will often be `<10` or `S`.
- **Recommendation:** read the audit at two levels:
  - **train**, which is used for every reference statistic and for the criteria;
  - **pooled** across splits, for descriptive reading only.
- Validation and test are reported separately for completeness. **No modeling decision is fit on validation or test**, and the pooled figures are descriptive, not a basis for thresholds.

## Outputs
Files:
- Shared: `data/interim/surgery_complexity_duration_audit_MAP_v1.json` / `.md` (aggregate-only, suppressed). The suffix follows the map label.
- Local only: `surgery_complexity_duration_audit_MAP_v1_raw_LOCAL.json`. It is never copied, uploaded or shared.
- Suggested script: `scripts/audit_surgery_complexity_duration.py`.

Unless stated otherwise, every item is reported **per body system** and **pooled over systems**, for **train / validation / test / pooled**.

1. **Population and validity:**
   - linked surgeries by system and timing class;
   - `unmapped`, `duration_fallback_anes`, invalid durations by reason, D > 24 h;
   - the excluded `Unclassified` surgeries.
2. **Per-group duration**, per group (cardiac: G1/G2/G3/NC and STAT 1–5; other systems: I1/I2/I3 and the four raw intensity levels):
   - n, median, IQR (p25, p75), p10 and p90, rounded to 0.1 h. Never min or max. Only for groups with n ≥ 10.
   - Counts per band (under 2 h, 2–4 h, over 4 h).
3. **Cross-tab group × band** per system and split, with row and column margins (suppressed), plus the row percentages computed from the published counts.
4. **Rank correlation** between the ordinal group and D, for assessable systems:
   - **Spearman ρ** and **Kendall τ-b** (tie-corrected). Cardiac uses G1/G2/G3; the STAT 1–5 version is a sensitivity analysis. Other systems use I1/I2/I3; the four-level intensity version is a sensitivity analysis.
   - **95% CI by cluster bootstrap:** resample `SUBJECT_ID` with replacement; 2,000 replicates; percentile interval; fixed seed recorded in the JSON. Computed locally on the raw durations. Only the point estimate and CI are published, and only when n ≥ 30 with every included group ≥ 10. Otherwise the result is `"insufficient"`.
   - Also τ-b on the group × band cross-tab. It is recomputable from the published table, and is published only if the table has no suppressed cell.
   - **Pooled over systems:** a **stratified** τ-b, in which pairs are compared only within a system (the n-weighted average of within-system τ-b, with a cluster bootstrap CI). Groups are never compared across systems here.
5. **Overlap metrics**, for each assessable system, between its lowest and highest eligible groups (cardiac: G1 vs G3, also G1 vs G2 and G2 vs G3):
   - **P_low_above_high_median:** the share of low-group surgeries with D above the high group's **train** median. For cardiac, this is the share of G1 surgeries longer than the G3 median.
   - **P_high_below_low_median:** the share of high-group surgeries with D below the low group's train median. For cardiac, the share of G3 shorter than the G1 median.
   - The **probability of superiority** P(D_high > D_low) + ½ P(tie), with a cluster bootstrap CI.
   - A share counts as a count for suppression: it is published only when its numerator and complement are each ≥ 10 (raw). Otherwise it is `"<10"`, or `"S"` when complementary.
6. **Per-name discordance list:** names with **at least 10 surgeries in the pooled analysis population** (names below 10 are never printed). Per name:
   - `name`, `body_system`, `STAT` (cardiac) or `intensity`, `complexity_group`, n (pooled), median D (pooled; also the train median if train n ≥ 10);
   - **`discordant` flag:** the name's median falls outside its group's **train** IQR within its system;
   - `direction` (`shorter` / `longer`), and whether the median also crosses an adjacent group's train median.
   - Split counts per name are suppressed under the usual rules.
7. **Flagged vs unflagged:** items 2–5 repeated separately for names with `review_flag = Y` and `N`. For cardiac, also for `STAT_group_crosses = Y` vs `N`.
8. **Sensitivity analyses** (items 2–5 repeated):
   - (a) multi-component names excluded (`n_components > 1`);
   - (b) cardiac names with `STAT_group_crosses = Y` excluded;
   - (c) D > 24 h and fallback durations excluded;
   - (d) `during_icu_layer_c` surgeries only (v1.7 class);
   - (e) the four-level intensity scale and STAT 1–5 in place of the 3-group versions.
9. **Across-system check (intensity calibration):** does one intensity level mean similar durations in different systems?
   - Per raw intensity level (and per I1/I2/I3): the median and IQR of D per system (systems with n ≥ 10), and the **pooled** level median and IQR (train).
   - **Flag** a system whose level median lies outside the pooled level's train IQR. Example: the `major open` median in Neurosurgical vs GI/abdominal vs Thoracic.
   - Also the Kruskal–Wallis ε² across systems within the level (published only when every included system has n ≥ 10).
   - For reference, cardiac G1/G2/G3 medians are listed next to the `major open` medians. They are not tested against them.
10. **Discrepancy verdict** per assessable system and pooled: each criterion's value, met or not, and `big_discrepancy` (Y/N), labelled `"approved_criteria_2026-11-10_D2signed"`.
11. **Markdown summary:** the tables for items 1–5 and 9, the discordance list, the verdicts, a not-assessable list with reasons, and a short note on small cells.

## Discrepancy criteria (pre-specified; APPROVED Nov 10, 2026)
**Status (v1.5): APPROVED by Vinai (Nov 10, 2026); D2 signed reading decided Nov 10, 2026.** D1–D3 are decisive with the cutoffs below (signed τ-b < 0.2; 25%). D4 and D5 are **reported only and non-decisive**.

Evaluated on **train** for each assessable system (cardiac on G1/G2/G3). The pooled-split values are shown alongside for description only. A system has a **big discrepancy** if **any** of D1–D3 holds.

**Verdict semantics (v1.4 clarification).**
- Each D is a *discrepancy condition*. `met` means the discrepancy condition **holds**, and `not_met` means it does not.
- `big_discrepancy = "Y"` **if and only if** at least one of D1–D3 is `met`. It is `"N"` if all three that can be computed are `not_met`, and `not_assessable` if none can be computed.
- Unit test (truth table): a synthetic system with strictly increasing medians, τ-b = 0.5 and a 10% overlap gives D1–D3 all `not_met` and `big_discrepancy = "N"`.

The conditions:

- **D1, non-monotone medians:** the group medians are not strictly increasing across eligible groups (cardiac: not G1 < G2 < G3; others: not I1 < I2 < I3 over the groups present). Ties count as non-monotone, so D1 is `met` when two medians are equal.
- **D2, weak or inverse rank association (signed; Vinai, Nov 10, 2026):** the signed Kendall τ-b is **< 0.2** (point estimate), with τ-b computed explicitly (`scipy.stats.kendalltau(..., variant="b")`).
  - D2 is `met` when the association is weak (−0.2 < τ-b < 0.2) **or inverse** (τ-b ≤ −0.2, meaning higher complexity groups have shorter durations). It is `not_met` only when τ-b ≥ 0.2.
  - This replaces the v1.2 wording |τ-b| < 0.2, which did not flag a strong inverse association.
  - Config: `TAU_B_SIGNED_MAX = 0.2`.
  - The signed τ-b and its subject-cluster bootstrap CI are always reported. If the CI includes 0, this is noted. |τ-b| may be published as a **descriptive field only** (`abs_tau_b`). It never enters a verdict.
- **D3, overlap:** more than 25% of the lowest group's surgeries are longer than the highest group's train median (cardiac: over 25% of G1 above the G3 median).
- **Supplementary, reported only and non-decisive (approved as such):** D4, the mirror of D3: more than 25% of the highest group shorter than the lowest group's train median. D5: more than 20% of listed names (n ≥ 10) in a system are `discordant`.

A criterion that cannot be computed (insufficient n) is `not_assessable` and does not count as met. The approved thresholds (D2 signed τ-b < 0.2, `TAU_B_SIGNED_MAX = 0.2`; D3 25%; and, for the reported-only items, D4 25% and D5 20%) are fixed in the script config with `CRITERIA_SIGNED = true` and `CRITERIA_APPROVED_DATE = "2026-11-10"`. They are not tuned after the audit is run; any change needs a new spec version.

## Suppression and disclosure (PhysioNet DUA)
The same as v1.2 "Suppression and disclosure". In brief:
- **Aggregate only:** no IDs, timestamps or individual durations.
- **Reconciliation:** raw counts reconcile before suppression and are not written out.
- **Primary suppression:** counts of 1–9 and zero are published as `"<10"`.
- **Complementary suppression:** published as `"S"`. It covers margins (group × band, system totals, systems → pooled), **cross-split** (train + validation + test = pooled), and shares treated as counts.
- **Quantiles:** only for groups with n ≥ 10. Never min or max.
- **Names:** only names with at least 10 surgeries are printed.
- **Recoverability self-check:** mandatory.
- **Release gate:** the disk-reload release gate prints `release gate ok`. Nothing is shared before it.

**Reconciliation identities** (raise an error on failure):
- groups + `NC` + excluded = system total;
- systems + `Unclassified` + `unmapped` = linked surgeries;
- valid + invalid = population;
- band counts = valid durations per group;
- train + validation + test = pooled for every cell.

## Output JSON schema (shared file)
```json
{
  "spec": "surgery_complexity_duration_audit v1.7",
  "input_hashes": {"canonical": "<64 hex>", "harmonization": "<64 hex>", "complexity_map": "<64 hex>", "split": "<64 hex>"},
  "reconciliation": {"surgery_info_rows": "", "linked": "", "post_icu_same_hadm": "", "other_hadm": "", "not_in_cohort": "", "invalid_time": "", "not_linkable_missing_encounter_time": "", "n_multi_candidate": "", "surgery_key_ROW_ID": {"shared_oper_id_groups": "", "rows_in_shared_oper_id_groups": "", "possible_true_duplicate_groups": "", "rows_in_possible_true_duplicate_groups": ""}},
  "map_status": "reviewed_v1",
  "map_label": "MAP_v1",
  "map_reviewed_date": "2026-11-10",
  "criteria": {"signed": true, "approved_date": "2026-11-10", "decisive": ["D1", "D2", "D3"], "reported_only": ["D4", "D5"], "tau_b_signed_max": 0.2, "d2_rule": "met if tau_b < tau_b_signed_max (signed)", "overlap_max": 0.25, "mirror_overlap_max": 0.25, "discordant_names_max": 0.20},
  "bootstrap": {"replicates": 2000, "cluster": "SUBJECT_ID", "seed": 0, "ci": "percentile 95%"},
  "population": {"<split>": {"linked": "n|<10|S", "by_timing_class": {"recent_layer_a": "", "remote_same_hadm": "", "during_icu_layer_c": ""}, "unmapped": "", "invalid_duration": {"reason": ""}, "duration_fallback_anes": "", "d_gt_24h": "", "straddling_t0": "", "recent_ending_after_t0": "", "timing_fallback_surgery_end": "", "pre_admit_time": ""}},
  "systems": {
    "<body_system>": {
      "grouping": "STAT_group | intensity_3",
      "assessable": true, "not_assessable_reason": null,
      "groups": {"<group>": {"<split>": {"n": "", "median": 0.0, "p25": 0.0, "p75": 0.0, "p10": 0.0, "p90": 0.0,
                  "bands": {"lt2h": "", "2to4h": "", "gt4h": ""}}}},
      "crosstab_group_band": {"<split>": {"<group>": {"lt2h": "", "2to4h": "", "gt4h": "", "total": ""}}},
      "correlation": {"<split>": {"spearman": {"est": 0.0, "ci": [0.0, 0.0]}, "kendall_tau_b": {"est": 0.0, "ci": [0.0, 0.0]}, "tau_b_from_crosstab": 0.0}},
      "overlap": {"<low>_vs_<high>": {"p_low_above_high_median": "", "p_high_below_low_median": "", "prob_superiority": {"est": 0.0, "ci": [0.0, 0.0]}}},
      "discordant_names": [{"name": "", "group": "", "stat_or_intensity": "", "n": 0, "median": 0.0, "train_median": 0.0, "discordant": true, "direction": "longer"}],
      "flagged_vs_unflagged": {"Y": {}, "N": {}},
      "sensitivity": {"a_single_component": {}, "b_no_cross_boundary": {}, "c_no_gt24_no_fallback": {}, "d_during_icu": {}, "e_four_level_or_stat5": {}},
      "verdict": {"D1": "met|not_met|not_assessable", "D2": "", "D3": "", "D4": "", "D5": "", "big_discrepancy": "Y|N|not_assessable", "label": "approved_criteria_2026-11-10_D2signed"}
    }
  },
  "pooled": {"stratified_tau_b": {"est": 0.0, "ci": [0.0, 0.0]}, "verdict": {}},
  "across_system": {"<intensity_level>": {"pooled_train": {"median": 0.0, "p25": 0.0, "p75": 0.0}, "by_system": {"<body_system>": {"n": "", "median": 0.0, "outside_pooled_iqr": false}}, "kruskal_epsilon2": 0.0}},
  "suppression": {"primary_cells": 0, "complementary_cells": 0, "recoverability_cells_checked": 0, "recoverable": 0},
  "release_gate": "ok"
}
```
Rules for the schema:
- Every count field is an integer ≥ 10, `"<10"` or `"S"`. This includes band counts, n per group, timing-class and data-quality counters, and the overlap numerators.
- Each system block carries `split_stats` for `train`, `validation`, `test` **and `pooled`**, with every map group present (see Definitions).
- `verdict` also carries the numeric values behind each criterion: `d1_medians_train`, `tau_b` (signed), `tau_b_ci`, `spearman`, `d3_share`, `d4_share` and `d5_share` (suppressed as shares-as-counts where needed), and optionally `abs_tau_b` (descriptive only). The v1.4 dual-reading fields `D2_signed`, `big_discrepancy_signed` and `d2_sign_disagreement` are **removed**.
- `suppression` counters are **computed, never hard-coded**. `release_gate` is written only by the disk-reload validator.
- Statistics that can't be computed are `"<10"` (group under 10) or `"insufficient"`.

## Acceptance checks (the script must raise on failure)
1. v1.2 checks 1–3 and 5 (hashes, headers, row counts, cohort and split totals, reconciliation) pass. `COMPLEXITY_MAP_SHA256` equals `5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5` exactly (all 64 characters). Map keys are unique. Every shared output (JSON, markdown, console summary) carries the label `MAP_v1` (reviewed 2026-11-10), and the output file names end in `_MAP_v1`.
2. `unmapped` = 0, or it is reported. The join uses the normalized name only, never an ID.
3. **Group mapping unit tests:**
   - STAT 1 and 2 → G1; 3 → G2; 4 and 5 → G3; `NC` → NC. NC is excluded from ρ, τ-b, D1 and overlap.
   - minor/diagnostic and catheter/interventional → I1; minimally invasive → I2; major open → I3.
   - Cardiac rows never use intensity for the group.
4. **Duration unit tests:**
   - D = E − B;
   - a missing B falls back to AE − AS and increments `duration_fallback_anes`;
   - B > E and B = E are invalid and counted;
   - band boundaries: 1.999 → under 2 h; 2 → 2–4 h; 4 → 2–4 h; 4.0001 → over 4 h.
5. **Train-only references.** Unit test: changing validation or test durations leaves every reference median, every IQR, every overlap threshold and every D1–D5 verdict on train unchanged.
6. **τ-b and ρ correctness:**
   - unit tests against a reference implementation on synthetic data with ties;
   - the bootstrap resamples subjects, not surgeries (unit test: a subject with 2 surgeries is drawn as a unit);
   - the seed is recorded and the result reproducible.
7. **Eligibility.** No statistic is published for a group under 10. Systems with fewer than 2 eligible groups are marked `not_assessable` with a reason, and no criterion is evaluated for them.
8. No quantile is computed on n < 10, and no min or max is published. No individual duration, ID or timestamp appears anywhere (sentinel-ID unit test).
9. The discordance list contains only names with at least 10 surgeries in the pooled population.
10. Every assessable system and the pooled level carry all of D1–D5 and a `big_discrepancy` verdict labelled `approved_criteria_2026-11-10_D2signed`. `big_discrepancy` depends on D1–D3 only.
    - A unit test confirms that changing D4 or D5 alone never changes the verdict.
    - The script raises an error if `CRITERIA_SIGNED` is not true, or if the approved values differ: `TAU_B_SIGNED_MAX = 0.2`, `D3 = 0.25`, `D4 = 0.25`, `D5 = 0.20`, and `D2_RULE = "signed"`.
    - **D2 unit tests (v1.5):** τ-b = 0.5 → `not_met`; τ-b = 0.2 → `not_met` (the boundary is not a discrepancy); τ-b = 0.19 → `met`; τ-b = 0.0 → `met`; τ-b = −0.5 → `met`. The test also asserts that no code path computes the D2 verdict from `abs(tau_b)`. The thresholds are read from config, not hard-coded in the logic.
11. Items 7 (flagged vs unflagged), 8 (sensitivities a–e) and 9 (across-system check) are present.
12. Primary and complementary suppression (margins, cross-split, shares), the recoverability self-check (0 recoverable) and the disk-reload release gate (`release gate ok`) pass.
13. No feature, band, merge, cluster or model output is written. No map file is modified.
14. **(v1.4; key = `ROW_ID` in v1.6) Linkage.** Each surgery (`ROW_ID`) appears at most once in the analysis set; linked ≤ 7,488; the reconciliation identity in Definitions holds; `n_multi_candidate` is reported. Unit test: two encounters in one HADM with one surgery during the first encounter gives exactly 1 linked row (`during_icu_layer_c` in v1.7, first encounter).
15. **(v1.4) Verdict semantics.** The truth-table unit test (Discrepancy criteria) passes, and `big_discrepancy` depends only on D1–D3 being `met`.
16. **(v1.4) Release gate is real.** The validator reloads the shared JSON and markdown from disk and fails if any published integer from 1 to 9 (or 0) appears in a count field, if any suppression counter was not computed, or if the recoverability check finds a recoverable cell. Unit test: injecting a band count of 5 makes the gate fail. The console prints `reconciliation ok` and `release gate ok` only after the corresponding checks have actually run and passed.
17. **(v1.4) No silent omission.** Every map group appears in every split block, as `"<10"` when n < 10.
18. **(v1.4) Metadata honesty.** The `bootstrap` block appears only if a bootstrap was actually run, and it records the realised replicate count.
19. **(v1.6) Surgery key.**
    - `ROW_ID` is nonblank and unique on all 7,488 rows; otherwise the run stops (exit 2).
    - Rows sharing an `OPER_ID` are all kept. Link classes sum to 7,488.
    - The four key-diagnostic counts are published (suppressed).
    - Unit tests on synthetic data cover:
      - a duplicate or blank `ROW_ID` (including a float-like `123.0` repeat) stops the run;
      - distinct surgeries sharing an `OPER_ID` are each linked;
      - the diagnostics equal an independent count;
      - a fixture with no shared `OPER_ID`s gives 0.
20. **(v1.6) Harmonization defaults.**
    - With the real header (`t0`, `OUTTIME`, no `OUT`), the script resolves t0 = `t0` and `OUT` = `OUTTIME` without flags. A primary `OUT` column takes precedence over `OUTTIME`, and an explicit flag beats both.
    - The `OUT` − t0 vs `span_icu_hours` disagreement count is reported. A unit test with a planted mismatch checks it.
21. **(v1.7) Anesthesia-end timing.** Boundary unit tests (anesthesia end relative to t0; e1 = t0..t0 + 48 h, e2 = t0 + 100..t0 + 150 h, same HADM):
    - t0 − 24 → recent; t0 − 24.01 → remote; t0 → recent; t0 + 6 → recent; t0 + 6.01 → Layer C; `OUT` → Layer C;
    - between encounters → the later encounter: 20 h before e2 → recent; 40 h and 51 h before e2 → remote;
    - after the last `OUT` → `post_icu_same_hadm`.

    Further tests:
    - no linked recent surgery is also Layer C, and the Layer C sensitivity (d) excludes it;
    - a missing `ANES_END_TIME` falls back to `SURGERY_END_TIME` (counted); both missing → `invalid_time`;
    - a missing B does not block linkage;
    - link classes still sum to 7,488.

## Results: MAP_v1 run (Nov 11, 2026; aggregate only, as reported to Vinai)
- **Run:**
  - Copilot ran the v1.7 script, SHA-256 `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`. This is identical to the reference v1.7 file.
  - Release gate: **ok**. Reconciliation: **ok**.
  - Suppression: 1,137 primary and 880 complementary cells; 2,703 cells checked by the solver; **0 recoverable**.
  - Bootstrap: 2,000 replicates, clustered on SUBJECT_ID.
  - Input hashes match the pinned values: canonical `e5feeb22…`, harmonization `ffd1b8ae…`, map `5c993b88…`, split `60196792…`.
- **Earlier rerun stopped:** an earlier rerun failed at the `OUT`-column gate because a stale pre-v1.7 copy of the script was still in Copilot's `scripts/`. Re-copying the v1.7 file fixed it, with no code change.
- **Verdicts** (`approved_criteria_2026-11-10_D2signed`; train sets every criterion):

  | system | grouping | train medians (h) | signed τ-b, train [95% CI] | big_discrepancy |
  |---|---|---|---|---|
  | Cardiac | STAT G1 / G2 / G3 | 1.9 / 3.9 / 3.1 (non-monotone) | 0.177 [0.142, 0.213] | **Y** |
  | GI/abdominal | I2 / I3 | 1.5 / 1.4 | 0.025 [−0.025, 0.079] | **Y** |
  | Hepatobiliary | I2 / I3 | 2.9 / 2.0 (reversed) | −0.262 [−0.389, −0.133] | **Y** |
  | Thoracic/airway | I2 / I3 | 1.3 / 1.0 | −0.066 [−0.15, 0.025] | **Y** |
  | Urologic/renal | I2 / I3 | 2.2 / 2.6 (monotone but weak) | 0.173 [0.068, 0.277] | **Y** |
  | Pooled (stratified, n-weighted within-system τ-b) | — | — | 0.076 [0.048, 0.103] | **Y** |

  - **Not assessable:**
    - Craniofacial/ENT/oral, Neurosurgical, Orthopedic/spine and Other: groups below 10;
    - Diagnostic/imaging-interventional and Vascular access/line placement: single group.
  - **Sensitivity analyses:** the verdicts are robust under (a) single-component names, (b) no cross-boundary cardiac items, (c) no > 24 h durations and no duration fallback, and (e) the four-level intensity scale or STAT 5. Sensitivity (d), during-ICU / Layer C surgeries only, is **not assessable**.
    - **Source:** the per-system sensitivity verdicts in the audit JSON (`systems.<system>.sensitivity.<a–e>.criteria_descriptive`), verified in review.
    - In `criteria_descriptive`, Cardiac, GI/abdominal, Hepatobiliary, Thoracic/airway and Urologic/renal keep `big_discrepancy = Y` under a, c and e. Sensitivity b applies to Cardiac only, and is also Y. Sensitivity d is `not_assessable` in every system.
    - Descriptive only: Diagnostic/imaging-interventional, not assessable in the primary analysis, becomes assessable under e (four-level scale) and gives N. This does not change any primary verdict.
- **Interpretation (Vinai):** STAT is a **mortality-based** risk score, and the intensity scale reflects the **surgical approach** (minimally invasive vs open) rather than complexity. Neither works as a proxy for operative duration.
  - Supporting detail from the name table (not Vinai's wording): Cardiac G1 spans short device closures (e.g. ASD closure, train median 0.6 h) and longer open repairs (e.g. tetralogy of Fallot repair, 3.0 h).
- **Harmonization consistency check (same day; real data, aggregate only):** no conflicts. Surgery timing classes (v1.7):

  | class | surgeries |
  |---|---|
  | recent (Layer A) | 6,008 |
  | remote | 850 |
  | Layer C | 486 |
  | post-ICU same HADM | 144 |
  | **total** | **7,488** |

  - 3,099 recent surgeries end in (t0, t0 + 6 h].
  - `in_gap` is < 10.
  - Details are in `docs/encounter_harmonization_consistency_review.md` and the preprocessing log, Step 5.2.

## Post-audit decisions (Vinai; not automated)
1. ~~Sign off the criteria~~ **Approved (Nov 10, 2026):** D1–D3 decisive (signed τ-b < 0.2, decided Nov 10, 2026; 25%); D4–D5 reported only.
2. ~~Confirm the intensity collapse~~ **Decided (Nov 10, 2026):** I1 = minor/diagnostic + catheter/interventional, I2 = minimally invasive, I3 = major open; the four-level scale is kept as a sensitivity analysis.
3. ~~Interpret each system's verdict~~ **Decided (Vinai, Nov 11, 2026): type-specific bands.** Every assessable system has big_discrepancy = Y, and so does the pooled analysis. Complexity groups are therefore **not** used as duration proxies.
   - STAT and intensity are kept as **separate descriptors**.
   - STAT is the cardiac **risk marker** (it is mortality-based).
   - Intensity describes the surgical approach.
   - Neither is used to stand in for duration.
4. ~~Decide whether any finding changes the φ1 definition~~ **Decided (Vinai, Nov 11, 2026): φ1 = surgery type × type-specific duration band.**
   - Each surgery's duration is banded relative to **its own type's train-split median**: shorter or longer than expected for that operation.
   - Rare types use a **pooled fallback by system × complexity group** (STAT group for cardiac, intensity group otherwise).
   - Medians and cuts are **fit on the train split only and frozen before modeling**.
   - **Approved defaults (Vinai, Nov 11, 2026; decision 16 follow-up):**
     - **Bands:** r = duration D / the train-split median of the surgery's own type.
       - **short:** r < 0.75;
       - **expected:** 0.75 ≤ r ≤ 1.33 (inclusive at both ends);
       - **long:** r > 1.33.
       - **Sensitivity:** within-type train tertiles.
     - **Own-type median:** only if the type has **≥ 10 train surgeries**. Sensitivity: n ≥ 20.
     - **Fallback order:** own type → system × complexity group (STAT group for cardiac, intensity I1/I2/I3 otherwise) → system → all surgeries.
       - Types not seen in train go straight to the fallback.
       - Medians come from the train split only and are frozen before modeling; the frozen level of each type is fixed.
       - Report counts at each fallback level by split. Every published cell must be ≥ 10 or suppressed.
     - **STAT and intensity:** separate descriptors. STAT is the cardiac risk marker. Neither is a duration proxy.
   - **Duration window (decided, Vinai, Nov 12, 2026):** the φ1 bands use **D = `SURGERY_END_TIME` − `SURGERY_BEGIN_TIME`**, falling back to the anesthesia window AE − AS (counted, `duration_fallback_anes`).
     - This is the same definition as in Definitions and the item-level spec, and the same one the audit's train medians use (e.g. Cardiac G1 1.9 h).
     - Both the type medians and each surgery's ratio are computed on this window.
   - This replaces the universal under 2 h / 2–4 h / over 4 h bands as the φ1 definition (item-level spec post-audit decision 4). The universal bands stay a descriptive and sensitivity level.
   - **Frozen and coded (Vinai, Nov 12, 2026, after the band-support audit):**
     - **Freeze** the train medians at **min n = 10**. The frozen medians file is `phi1_band_support_frozen_medians_MAP_v1.csv`, SHA-256 `e73c8b5a8d780b51edfd0935abdb05335f57e3a0cfc0c4c52eaaa85888db0047`.
     - **Coding:** the duration band is coded as a **3-level ordinal (short / expected / long)** alongside the type and system descriptors. It is **not** a full type × band cross-product, because within-type band cells are sparse.
     - **Type naming (resolved):** the type median uses `standardized_name`, with the fallback hierarchy for rare types. Collapsing combination names to their primary component is a **sensitivity analysis only**.

### Band-support audit results (MAP_v1, Nov 12, 2026; aggregate only)
Full results are in `band_support_audit_spec.md`, "Results".
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

### Open items for φ1 bands (remaining; minor)
The cut points, the rare-type threshold and the fallback hierarchy were **decided Nov 11, 2026** (decision 4 above):
- ratio cuts < 0.75 / 0.75–1.33 / > 1.33, with within-type tertiles as a sensitivity;
- min n = 10, with n ≥ 20 as a sensitivity;
- fallback type → system × complexity group → system → all.

These remain open:
1. Handling of durations > 24 h and of ordering violations. Each is < 10 in every split (band-support audit). They keep the v1.2 item-level handling until decided.
2. Whether the expected duration (the type median) is shown in the Layer C start bin, as in decision 12 of the contract amendment.

The type-median name was **resolved Nov 12, 2026**: `standardized_name`, with combination-name collapsing as a sensitivity analysis only.

The duration window was **decided Nov 12, 2026**: the surgery window, falling back to the anesthesia window (decision 4).

## Do not
- Do not assign or edit STAT categories, STAT groups, intensity or body systems.
- Do not fit bands, thresholds or groupings on validation or test data.
- Do not print names with fewer than 10 surgeries, row-level data, IDs or timestamps.
- Do not build features or train models.

## Changelog
- **v1.7, band-support results and frozen medians (Nov 12, 2026). Documentation only: the version string, script (`ce894b20…`) and hashes are unchanged.**
  - Added a summary of the band-support audit results.
  - Recorded the decisions: freeze at min n = 10 (frozen CSV SHA-256 `e73c8b5a…0047`); the band is a 3-level ordinal, not a type × band cross-product; combination-name collapsing is a sensitivity analysis only.
  - Removed the type-name open item.
- **v1.7, duration window decided (Nov 12, 2026). Documentation only: the version string, script (`ce894b20…`) and hashes are unchanged.**
  - Vinai confirmed that the φ1 bands use D = `SURGERY_END_TIME` − `SURGERY_BEGIN_TIME`, falling back to the anesthesia window. This is the same window as Definitions and the audit medians.
  - The Nov 11 "to confirm" flag is removed; it arose from a reversed window in the brief.
  - Open items are now only the three minor ones.
- **v1.7, φ1 band defaults (Nov 11, 2026). Documentation only: the version string, script (`ce894b20…`) and hashes are unchanged.**
  - Recorded Vinai's approved defaults: ratio-to-type-median cuts < 0.75 / 0.75–1.33 (inclusive) / > 1.33, with tertiles as a sensitivity; min n = 10 (sensitivity 20); fallback type → system × complexity group → system → all; types unseen in train go straight to the fallback; train-only, frozen.
  - The open items are reduced to the minor ones, plus a flagged duration-window question.
  - Corrected the robustness attribution: it now comes from the JSON `criteria_descriptive` sensitivity verdicts, verified in review, not from "Vinai's reading".
- **v1.7, post-run documentation update (Nov 11, 2026). Documentation only: the version, script (`ce894b20…`) and hashes are unchanged.**
  - Added the MAP_v1 run results: release gate ok; 1,137 / 880 suppressed cells, 0 recoverable; big_discrepancy = Y for Cardiac, GI, Hepatobiliary, Thoracic, Urologic and pooled; 6 systems not assessable; sensitivity (d) not assessable.
  - Recorded the harmonization-check timing classes.
  - Post-audit decisions 3 and 4 decided: φ1 = type × type-specific duration band, with a pooled system × complexity-group fallback. STAT and intensity are descriptors only, and STAT is the cardiac risk marker.
  - Added the open items: band cuts (ratio < 0.75 / 0.75–1.33 / > 1.33 vs within-type tertiles), min n = 10, the fallback hierarchy.
  - Recorded the stale-copy rerun stop.
- **v1.7 (Nov 11, 2026):** Vinai approved the harmonization-review decisions. **Surgery timing is by anesthesia end T**: `ANES_END_TIME`, falling back to `SURGERY_END_TIME`, counted.
  - Classes:
    - recent (Layer A only): T in [t0 − 24 h, t0 + 6 h];
    - Layer C: t0 + 6 h < T ≤ `OUT`, merge gaps included;
    - remote: T < t0 − 24 h, same HADM;
    - between encounters: the later encounter.
  - Class names are now `recent_layer_a`, `remote_same_hadm` and `during_icu_layer_c`.
  - `begin_equals_t0` is removed (moot). `recent_ending_after_t0` and `timing_fallback_surgery_end` are added. `pre_admit_time` now uses T.
  - Reuse of a surgery as later-encounter Layer A history is stated explicitly (features only).
  - Acceptance check 21 added.
  - Unchanged: `ROW_ID` key, reconciliation to 7,488, map, hashes, criteria, duration definition (still B→E, with the anesthesia fallback) and suppression.
- **v1.6 (Nov 11, 2026):** Copilot's real run stopped on duplicate (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) keys. In the aggregate diagnostic, 889 rows fall in 389 groups, with 0 exact copies; all groups differ in anesthesia and surgery times, and most differ in name. **Vinai decided (Nov 11): use `ROW_ID` as the surgery key and keep all 7,488 rows.**
  - `ROW_ID` is asserted unique and nonblank. `OPER_ID` becomes a descriptor.
  - New key diagnostics (shared-`OPER_ID` groups; possible true duplicates by subject/HADM/begin/name) are reported, not dropped.
  - `OUTTIME` is a recognized fallback name for the harmonized `OUT`. The real harmonization file has `t0` and `OUTTIME` and no `OUT`, and Copilot's run needed `--harm-out-col OUTTIME`.
  - Reconciliation wording now lists `not_linkable_missing_encounter_time`. The JSON example lists 4 input hashes (split included), matching the script and handoff.
  - Acceptance checks 19–20 added.
  - Map, hashes, criteria, thresholds and suppression are unchanged.
  - The item-level spec v1.2's (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) key and vitals join are **not** changed here. They are flagged for decision in `specs/encounter_harmonization_consistency_review.md`.
- **v1.5 (Nov 10, 2026):** Vinai decided **D2 uses the signed reading**: met if τ-b < 0.2 (weak or inverse association), replacing |τ-b| < 0.2.
  - Config is now `TAU_B_SIGNED_MAX = 0.2` and `D2_RULE = "signed"` (replacing `tau_b_min`). The verdict label is `approved_criteria_2026-11-10_D2signed`.
  - The v1.4 dual reporting (`D2_signed`, `big_discrepancy_signed`, `d2_sign_disagreement`) is removed. |τ-b| is kept only as the descriptive `abs_tau_b`.
  - Acceptance check 10 is updated (approved values and D2 unit tests).
  - The map, hashes and the other thresholds are unchanged.
- **v1.4 (Nov 10, 2026):** clarifications after the code review of the first run (`surgery_complexity_duration_audit_FIXES_v1.md`). The map, hashes and approved thresholds are unchanged.
  - Explicit verdict semantics: `met` means the discrepancy holds; `big_discrepancy = Y` if and only if any of D1–D3 is met. Ties count as non-monotone for D1.
  - D2 uses explicit τ-b, reports the signed value, and runs a signed reading (τ-b < 0.2) in parallel as a **proposed clarification pending Vinai**. The verdict keeps the approved |τ-b| wording until Vinai confirms.
  - A one-surgery-one-encounter assignment rule, using harmonized t0/`OUT` and `hadm_id_list`, and a reconciliation to 7,488.
  - Groups below 10 are published as `"<10"`, not dropped. The pooled split block is required.
  - Schema: numeric criterion values; suppression counters computed, not hard-coded.
  - Acceptance checks 14–18 added.
- **v1.3 (Nov 10, 2026):** the final reviewed map is pinned.
  - Vinai reviewed the type map by system and the Cardiac_STAT sheet and approved them as-is, with no edits; flagged items are accepted at their default assignment.
  - Map v1 is `surgery_complexity_map_v1.csv`, with SHA-256 `5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5`.
  - Inputs §2 is updated to the map's actual columns (`raw_name` join key, `I_group`, `NA` for Unclassified).
  - Outputs are labelled `MAP_v1` (reviewed 2026-11-10).
  - Final run authorized.
- **v1.2 (Nov 10, 2026):** Vinai approved the discrepancy criteria. D1 (non-monotone medians), D2 (|Kendall τ-b| < 0.2) and D3 (> 25% of the lowest group above the highest group's train median) are decisive. D4 and D5 are reported only and non-decisive. The verdict label is now `approved_criteria_2026-11-10`, the JSON `criteria` block is updated, acceptance check 10 is tightened (D1–D3 only, and approved values enforced), and post-audit decision 1 is marked done. Status: ready for Copilot.
- **v1.1 (Nov 10, 2026):** non-cardiac intensity collapse decided by Vinai: I1 = minor/diagnostic + catheter/interventional (was I1 = minor/diagnostic only), I2 = minimally invasive (was catheter + minimally invasive), I3 = major open; four levels kept as sensitivity. Added the system × I-group name-list table with the comparable-systems column (and complementary suppression in the four-level table); updated the expected not-assessable list (Diagnostic/imaging-interventional becomes single-group under the collapse), the mapping unit test (acceptance check 3) and post-audit decision 2. Discrepancy criteria D1–D3 unchanged and still proposed pending sign-off.
- **v1 (Nov 10, 2026):** initial spec, from Vinai's request (t445u) to audit STAT-based cardiac grouping (G1 = STAT 1–2, G2 = 3, G3 = 4–5) against actual surgery duration. Extended (t446u) to every body system using the type map's intensity levels, with per-system and pooled reporting and an across-system intensity calibration check. Supersedes the planned `cardiac_stat_duration_audit_spec` name.
