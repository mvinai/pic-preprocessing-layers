# Encounter harmonization: consistency review (documents and code)

> **Resolved Nov 11, 2026 (approved decisions, Vinai):** the open decisions in this review are settled. Surgery timing is by anesthesia end. Layer A recent is [t0 − 24 h, t0 + 6 h]; Layer C is (t0 + 6 h, OUTTIME] including gaps. Remote surgeries before ADMITTIME are kept and flagged, and a surgery between encounters goes to the later encounter. Vitals are assigned to a `ROW_ID` by anesthesia window (±30 min sensitivity). Labs enter at release time everywhere. `icu_time_hours` = span − gaps. See contract v12.1 amendment A2, duration spec v1.7 and item-level spec v1.3. The review text below is kept as written.

> **Real-data check results (Nov 11, 2026; Copilot ran `check_encounter_harmonization_consistency.py`, aggregate only): no conflicts. Both [VERIFY] items are resolved.**
> - **Rebuild [§11]:** all canonical encounters rebuild exactly from ICUSTAYS.
> - **OUTTIME [§4, §11] (finding 4, VERIFY resolved):** harmonized `OUTTIME` = the **maximum** stay OUTTIME in every encounter, not the last stay's.
> - **`icu_los_hours` [§4] (finding 5, VERIFY resolved):** `icu_los_hours` = `span_icu_hours` (t0 to OUTTIME, gaps included) in all canonical encounters.
> - **`icu_time_hours`:** = span minus gaps. It differs from the span in 100 encounters. The 40 nested stays explain why `icu_time_hours` ≠ the sum of stay durations, which confirms the log Step 1.5 definition (contract amendment A2.6).
> - **Surgery timing [§12]** (anesthesia end; one surgery, one encounter):
>
>   | class | surgeries |
>   |---|---|
>   | recent Layer A | 6,008 |
>   | remote | 850 |
>   | Layer C | 486 |
>   | post-ICU same HADM | 144 |
>   | **total** | **7,488** |
>
>   - 3,099 recent surgeries end in (t0, t0 + 6 h].
>   - `in_gap` < 10 (finding 8: negligible).
> - **Vitals [§13]:**
>
>   | rows | count |
>   |---|---|
>   | total `SURGERY_VITAL_SIGNS` rows | 1,944,187 |
>   | rows with no `SURGERY_INFO` key | 17,640 |
>   | keyed rows | 1,926,547 |
>   | of which assigned (strict anesthesia window) | 1,884,816 (97.8%) |
>   | of which unassigned (strict window) | 41,731 |
>   | unassigned at ±30 min | 10,196 |
>   | ambiguous | 0 |
>
>   **Primary rule (Vinai, Nov 11, 2026): the strict window**, with ±30 min as the sensitivity analysis.
> - **Admit time:** 26 encounters have a hospital admit time more than 1 h after t0. Read this alongside the remote-surgery "before ADMITTIME" flag, since both reflect admit-timestamp artifacts.
Nov 11, 2026. Requested by Vinai ("make sure all the harmonization related to encounters we did up until doesn't cause any conflict").

**Scope.** This is a document and code review only. No real data was read; the companion check script was tested on synthetic data. Approved specs and the log were **not** edited for this review. The duration-audit spec, fix list, handoff and log were updated separately under Vinai's ROW_ID decision (Nov 11; spec v1.6 / Amendment B).

**Sources reviewed:**
- the contract and its drafts:
  - `contract/foundations_contract_v12.md` (H1–H10) and its review checklist;
  - `contract/encounter_harmonization_section.md` (earlier draft, superseded);
  - `contract/v12_1_amendment_landmark_forward_outcome_DRAFT.md`;
- `preprocessing_log/foundations_preprocessing_log.md`, Steps 1.3–1.6, 2.x, 3.x, 4.6–4.10, 5.1 and 5.2;
- the specs:
  - `specs/landmark_eligibility_audit_spec.md`;
  - `specs/landmark_forward_band_support_audit_spec.md` (v2.3);
  - `specs/surgery_item_level_audit_spec.md` (v1.2);
  - `specs/surgery_complexity_duration_audit_spec.md` (v1.5 → v1.6), with its fix list and handoff;
- `scripts/duration_audit_v2/audit_surgery_complexity_duration.py`;
- the real headers Vinai supplied for the harmonization file (`ffd1b8ae…`) and the canonical file (`e5feeb22…`).

**Classes:**
- **CONFLICT:** two sources disagree, and following both is impossible.
- **NAMING:** same meaning, different names.
- **UNVERIFIABLE:** documents agree or are silent, and only the data can confirm. Each has a proposed aggregate-only check (counts 1–9 → `<10`).

Companion check: `scripts/check_encounter_harmonization_consistency.py`. It is read-only and aggregate-only. It verifies both hashes, prints counts only (no IDs or timestamps), and uses `<10` plus complementary `S` on partitions. 9 synthetic tests pass under the default locale and under `LC_ALL=C`. Section numbers [§n] below refer to its output.

---

## 1. Ranked findings

| # | finding | class | impact | where |
|---|---|---|---|---|
| 1 | **Surgery key and vitals join on (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) are invalid on real data.** 889 `SURGERY_INFO` rows fall in 389 groups of distinct surgeries (Copilot diagnostic, Nov 11). The item-level spec v1.2 keys every surgery and the `SURGERY_VITAL_SIGNS` join on this triple, expects 0 duplicates, and its monitor-time dedup groups by it. | CONFLICT (spec vs data) | **High**: the item-level run will fail acceptance checks 17/24, and vitals from different surgeries in a group would be pooled and cross-deduplicated. | item-level spec lines 111, 133, 180, 185, 306, 313. Duration audit already fixed (v1.6, `ROW_ID`). |
| 2 | **Where during-ICU surgery goes.** Contract v12 (H3, line 140): B ≤ INTIME → Layer A; B > INTIME → **Layer B**. Step 5.2 and amendment decision 12: during ICU → **Layer C**; Layer B is physiology only. | CONFLICT (contract vs log/amendment) | High for feature building. Resolved in substance by Step 5.2 ("replaces the v12 H3.4 wording"), but contract v12 still says Layer B until v12.1 is merged. | contract line 140; log Step 5.2; item-level spec line 121 |
| 3 | **B = t0 edge.** The contract puts B ≤ INTIME in Layer A. The item-level spec, duration spec and script classify B = t0 as `during_icu` (Layer C) and count it as `begin_equals_t0`. | CONFLICT (open edge, already flagged in item-level spec §decisions) | Low–medium: few rows, but it decides which layer the feature lands in. [§12] counts it. | contract line 140; item-level line 121 |
| 4 | **`OUTTIME` semantics (max vs last stay) are not proven in the data.** Every window uses the encounter end: Layer B/C end, landmark eligibility, `during_icu`, Outcome 2 `r` and readmission gap g. The log says the running end is the max OUTTIME and that the "latest-discharge" and "last-stay" rules gave identical *encounters* (membership). For a nested last stay, though, the *value* differs: last-stay OUTTIME < max OUTTIME. Contract H1 line 104 is still marked [VERIFY]. The real harmonization file has no `OUT` column, so the duration script reads `OUTTIME` as `OUT`. | UNVERIFIABLE | **High** if wrong: 40 nested stays could shorten encounters, changing eligibility and outcomes. | [§4] `abs(OUTTIME − t0 − span_icu_hours)`; [§11] rebuild from ICUSTAYS: `harm_OUTTIME != rebuilt max OUTTIME` (must be 0) and `… != last stay` (shows how many encounters are exposed). The duration script now prints the same span check. |
| 5 | **`icu_time_hours` definition.** Contract H7: sum of raw stay durations. Log Step 1.5: span minus within-encounter gaps ("so overlapping stays aren't double-counted"). These differ for the 88 overlapping stays. Contract H7 also leaves `icu_los_hours` as [VERIFY] (span or ICU time). | CONFLICT (docs) + UNVERIFIABLE (which one the file holds) | Low for modelling (banned as features), medium for descriptive tables and the > 1,440 h span check. | [§4] `abs(span − icu_time − gap_hours_total)` and `icu_los_hours_c_matches`; [§11] `abs(icu_time_hours − sum of stay durations)` |
| 6 | **One surgery feeding two encounters.** In a same-HADM ICU readmission, a surgery during encounter 1 is a Layer C event for encounter 1 and also a remote same-HADM Layer A feature for encounter 2 (Step 5.2). This is legitimate as history, but it is stated nowhere. The duration audit's one-surgery-one-encounter rule is audit-only (spec line 84); the item-level spec v1.2 has no such rule, so its timing-class counts can exceed 7,488. | CONFLICT (silent; needs wording) | Medium: double counting in item-level tables; ambiguity for feature builders. | [§12] `surgeries_linked_to_2+_encounters_of_the_same_HADM`; [§8] `next_encounter_in_same_HADM` |
| 7 | **Remote-surgery rule differs between the specs.** Item-level spec: `remote_same_hadm` requires B ≥ ADMITTIME, else `pre_admit_time`; the wording does not say whether those rows are excluded. Duration script: `remote_same_hadm` regardless, with `pre_admit_time` as a flag *within* it. ADMITTIME is after t0 in 2,582 encounters, by a median of about 1 s, a rounding artifact (log Step 1.5). | CONFLICT (minor) | Low: a handful of rows. | item-level line 116–117; duration script `pre_admit` |
| 8 | **Surgery inside a merge gap.** Encounters absorb gaps of 24 h or less: 112 gaps totalling 709.97 h. A surgery with t0 ≤ B ≤ OUTTIME that falls in a gap is `during_icu` / Layer C by the encounter clock, although the patient was out of the ICU. That conflicts with Step 5.2's premise that Layer C is ICU-documented time, with OR time as a gap. | UNVERIFIABLE | Low–medium: possibly the typical "out to theatre, back within 24 h" case. | [§12] `of_which_B_in_a_merge_gap` (needs `--icustays`) |
| 9 | **Lab release delay vs windows.** Amendment decision 11 says labs count from CHARTTIME + 1 h, and states it for Layer C ("lab results count only once released"). Two gaps follow: (a) whether Layer B (`[t0, t0+24 h)`) windows by CHARTTIME or by release time is not stated; (b) a lab collected in [t0 − 1 h, t0) is released after t0. It is outside the pre-ICU lookback by availability and outside Layer B by CHARTTIME, so it lands in no layer. | CONFLICT (gap in the rules) | Medium for Layer B lab features. | amendment lines 67, 118–120; contract H2 and H3 item 5 (pre-ICU lookback, line 141) |
| 10 | **Layer C length.** Contract H2: 28 bins / 168 h. Amendment A1.4: 4–56 bins (landmarks to 336 h). | CONFLICT, resolved by the amendment | None once v12.1 is merged. | contract line 121; amendment line 67 |
| 11 | **Death attribution wording.** Log: deaths "not on the patient's last encounter" (counts omitted from this public copy). Contract H4.3: "last encounter of the admission". A death followed by a later encounter in a *different* HADM is impossible, so the two coincide unless the data are wrong. | NAMING (verify) | Low. | [§7] both rules vs the flag (both must be 0) |
| 12 | **Split held in three places:** the split file, the canonical `split` and the harmonization `split`. The duration script checks canonical vs the split file (and stops on disagreement); the harmonization `split` is unchecked. The split-file hash `60196792…` comes from the run-1 script and is not independently verified. | UNVERIFIABLE | Low (expected 0). | [§9] |
| 13 | **Landmark eligibility.** Contract, amendment, eligibility spec and band-support spec all use OUTTIME > t0 + L **strictly**. The eligibility spec also says "alive", which is implied because death ends the stay. Counts are [N_eligible omitted] / [N_eligible omitted] / [N_eligible omitted] everywhere. The `in_icu_at_Lh` flags should equal the strict rule computed from OUTTIME. | consistent; UNVERIFIABLE for the flags | Low. | [§6] flag vs strict and vs ≥ |
| 14 | **Readmission rules vs the encounter merge.** A same-HADM ICU return more than 24 h after OUT is a readmission (Outcome 2 level 3; `readmit_index_overlap = keep_both`). Returns of 24 h or less are merged by construction, so `return_within_24h` must be 0. This is consistent with the surgery rules: outcomes look forward and surgery features look back. | consistent | None. | [§8] `next_encounter_starts_within_24h_of_OUTTIME` = 0 |
| 15 | **Duration spec documentation gaps:** the reconciliation omitted `not_linkable_missing_encounter_time`; the JSON example listed 3 hashes, not 4; the `OUT` column name was absent from the real file. | CONFLICT (spec vs script) | **Fixed** in spec v1.6 / Amendment B; the script falls back to `OUTTIME`. | — |
| 16 | **Spans over 1,440 h** (354; max 28,003.70 h) were flagged for checking "span against actual ICU time". | UNVERIFIABLE | Low (descriptor), but it may reveal bad OUTTIMEs. | [§4] `span_gt_1440h`, `icu_time < span` |

### Naming-only mismatches (docs vs real headers)
| documented (contract H7 / drafts / specs) | real header |
|---|---|
| `icu_span_hours`; "ICU length of stay" | `span_icu_hours` (also `icu_los_hours`, inherited; **resolved Nov 11, 2026: `icu_los_hours` = `span_icu_hours` in all [N_encounters omitted]**) |
| `sum_icu_stay_hours` (log), `icu_time_hours` (contract) | `icu_time_hours` |
| `n_gaps`, `total_gap_hours` | `gap_count_between_stays`, `gap_hours_total` |
| `age_at_t0` (contract), `age_at_t0_days` (superseded draft) | `age_years_at_intime` (units are years) |
| `prior_encounters_n` | `prior_encounters` |
| `OUT` (band-support, item-level and duration specs) | `OUTTIME` (in both files; there is no `OUT` column) |
| harmonized t0 | `t0` (also equal to `INTIME`, per log Step 1.5) |
| `SURGERY_END_TIME`, `age_at_t0_days`, last-encounter death as primary (`encounter_harmonization_section.md`) | superseded by contract v12: E is used only for duration; age is in years; inherited death is primary |

Recommendation: when v12.1 is merged, add a name-mapping table to contract H7. Do not rename the frozen or accepted files, because their hashes are pinned.

---

## 2. Item-level spec v1.2: surgery key and vitals join (flag; spec not edited)
Under the ROW_ID decision, (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) no longer identifies a surgery for 389 groups (889 rows, about 11.9%). `SURGERY_VITAL_SIGNS` has no `ROW_ID` link to `SURGERY_INFO`; it carries only `SUBJECT_ID`, `HADM_ID`, `VISIT_ID`, `OPER_ID`, `ITEM_NO` and `MONITOR_TIME`. So for those groups, monitoring rows cannot be attached to a surgery by key.

**Options (Vinai to choose; none applied):**
- **A. Window assignment (recommended primary).**
  - Within each shared group, assign every vitals row to the `ROW_ID` whose anesthesia window [`ANES_START_TIME`, `ANES_END_TIME`] contains `MONITOR_TIME`. Unique-key groups keep the direct join.
  - Count rows assigned to exactly one surgery, rows that are ambiguous (two or more windows contain the time, i.e. overlapping surgeries) and rows that are unassigned (no window).
  - Ambiguous and unassigned rows get no intraoperative vitals features in the primary analysis.
- **B. Tolerance sensitivity:** the same as A with ±30 min around the window, or with the surgery window [B, E]. Report how many rows move.
- **C. Restrict:** intraoperative-vitals features (φ2 cluster inputs) only for unique-key surgeries. This is simple, but it drops about 12% of surgeries, probably the multi-procedure ones, which risks bias.
- **D. Session level:** treat a shared group as one operative session. Pool its vitals and treat the surgery type as multi-label. This contradicts the ROW_ID decision (distinct surgeries with distinct times), so it is not recommended.

**Also needed in v1.3 of the item-level spec, if A is chosen:**
- key = `ROW_ID`;
- acceptance 17/24 reworded: the `OPER_ID` uniqueness check becomes a reported count, not a failure;
- the monitor-time duplicate rule (same key + ITEMID + time) becomes per assigned `ROW_ID`;
- the timing-class tables either adopt the one-surgery-one-encounter rule or state that rows can count twice (finding 6).

**Proposed aggregate check before deciding** (Copilot; the companion script [§13], with `--surgery-vitals`):
- vitals rows in shared groups by: assigned to exactly one ROW_ID, ambiguous, unassigned. Shown for the anesthesia window and for ±30 min, with complementary `S`;
- as a baseline, unique-key vitals rows that fall outside their surgery's anesthesia window.

---

## 3. Recommended fixes needing Vinai's decision
1. **Item-level key and vitals join:** choose A / B / C (section 2). This is required before the item-level run.
2. **B = t0:** Layer A (contract) or `during_icu` / Layer C (specs)? Recommendation: during ICU (Layer C), since the patient is in the ICU at t0. Record it in v12.1, and keep `begin_equals_t0` as a count.
3. **Merge contract v12.1** so that it carries:
   - surgery timing per Step 5.2 (H3.4 replaced);
   - Layer C as 4–56 bins;
   - the name-mapping table;
   - `icu_time_hours` = span − gaps, if [§4]/[§11] confirm it; otherwise correct the log.
4. **Surgery reuse across encounters:** state explicitly that one surgery may be a Layer C event for one encounter and remote-history Layer A for a later encounter in the same HADM. Give the item-level counts both per-encounter and per-ROW_ID.
5. **Gap surgeries:** if [§12] shows ≥ 10, decide whether a surgery in a merge gap is Layer C (encounter clock) or a separate "between-stays" event. Recommendation: Layer C with a `in_gap` flag.
6. **Lab release:**
   - state whether Layer B windows use CHARTTIME or release time;
   - state where a lab collected in [t0 − 1 h, t0) belongs.
   - Recommendation: windows by release time everywhere. A pre-t0 collection released after t0 enters Layer B, which matches "available at" semantics.
7. **Remote surgery before ADMITTIME:** keep it in `remote_same_hadm` with a flag (duration script) or exclude it (item-level). Recommendation: keep it and flag it (a rounding artifact).
8. **Run the companion check once (Copilot) and attach its aggregate output to the log.** It confirms findings 4, 5, 8, 11, 12, 13 and 16 in a single pass.

## 4. Proposed Copilot run (aggregate only)
```
python scripts/check_encounter_harmonization_consistency.py ^
  --icustays "<PIC_RAW_DIR>\ICUSTAYS.csv" ^
  --surgery-info "<PIC_RAW_DIR>\SURGERY_INFO.csv" ^
  --surgery-vitals "<PIC_RAW_DIR>\SURGERY_VITAL_SIGNS.csv" ^
  --json-out data\interim\harmonization_consistency_check.json
```
- Both hashes are checked before anything is read; on a mismatch it exits with code 2.
- Output is counts only, with 1–9 printed as `<10`; documented totals are shown as value / documented / matches.

**Expected if everything is consistent:**
- all `disagreeing_rows_per_shared_column` = 0;
- `harm_OUTTIME != rebuilt max OUTTIME` = 0;
- `abs(OUTTIME − t0 − span_icu_hours) > 1 min` = 0;
- the strict landmark-flag mismatch = 0;
- both death-rule mismatches = 0;
- split mismatches = 0;
- `next_encounter_starts_within_24h` = 0.

Send back the console output and the JSON. Neither contains IDs or timestamps.
