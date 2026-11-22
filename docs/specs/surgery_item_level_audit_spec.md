# Copilot spec: Surgery item-level audit (read-only, aggregate-only)
Version 1.3, Nov 11, 2026. **Approved decisions, Nov 11, 2026 (Vinai)** change three things:
- the surgery key becomes `ROW_ID`;
- timing is by **anesthesia end**, with one surgery per encounter;
- `SURGERY_VITAL_SIGNS` rows are assigned to a `ROW_ID` by anesthesia window.

See the Changelog, the duration-audit spec v1.7 and `docs/encounter_harmonization_consistency_review.md`. v1.2 (Nov 10, 2026) aligned the spec with the public PIC documentation: D_ITEMS-based SV map, ID semantics, monitor-time duplicates. Implements the surgery decision recorded as v12.1 draft amendment decision 12 (`contract/v12_1_amendment_landmark_forward_outcome_DRAFT.md`) and preprocessing log Step 5.2. Style, suppression and release gate follow `specs/landmark_forward_band_support_audit_spec.md` (v2.3).

**Public documentation used** (v1.2): the PIC documentation pages `surgery_vital_sign`, `d_items` and `or_exam_report` (pic-doc.nbscn.org; source in github.com/Healthink/PIC, `pic-doc/src/components/Pages/`), and the official demo `demo/bloodpressure_mortality_demo.ipynb` in the same repository. These are public documentation, not patient data.

## Purpose
Before the surgery features are built, measure what PIC's surgery data can support:

- the surgery names, and (after Vinai's manual map) their types;
- timing relative to the encounter;
- the anesthesia and surgery windows and their data quality;
- the record links;
- duration-band support;
- the training-split median duration by type;
- coverage of the intraoperative monitor table (`SURGERY_VITAL_SIGNS`);
- availability of pre-incision physiology;
- per-type support for the φ2 clusters.

Vinai uses the result to finalize the type map, the SV ITEMID map, the pre-incision source, the duration bands and the cluster inputs. These are post-audit decisions, and the script applies none of them. **Copilot never assigns surgery types.** This is a design audit only. It builds no features, fits no clusters or models, and changes no frozen files.

## Decisions being audited (for context; do not change)
- **Timing (approved Nov 11, 2026; replaces the B-based rule)**, by **anesthesia end** T = `ANES_END_TIME` relative to the encounter's harmonized t0 (`INTIME`) and `OUT` (`OUTTIME`). If T is missing, `SURGERY_END_TIME` is used, and this is counted.
  - T in [t0 − 24 h, t0 + 6 h] → Layer A recent-surgery flag. **Layer A only: never also a Layer C event.**
  - t0 + 6 h < T ≤ `OUT` → Layer C event, **including surgeries in merge gaps** between merged stays. There is no separate gap class; `in_gap` is a descriptor count only.
  - Same HADM and T < t0 − 24 h → Layer A comorbidity-type feature (remote surgery). Remote surgeries before `ADMITTIME` are kept and flagged.
  - Between two encounters → assigned to the **later** encounter and classed by the rules above, relative to its t0.
  - **One surgery may be a Layer C event for one encounter and Layer A history for a later encounter in the same admission.** In this audit, each surgery is still counted once (one surgery, one encounter); the reuse is reported as a count.
  - Surgery in earlier hospitalizations is **not** used; the Layer A prior-admission features cover it.
- **Levels, fine to coarse:**
  - φ0: the actual surgery, with rare names pooled.
  - φ1: type × duration band. Proposed bands are under 2 h, 2–4 h and over 4 h.
  - φ2: clusters, fit on train only from type plus pre-incision and intraoperative summaries, and frozen before modeling.
  - Type alone is a sensitivity level.
- **Layer C timing:**
  - Start bin: type, expected duration (train median for the type) and pre-incision physiology.
  - End bin: actual duration band and intraoperative summaries.
- **Layer C uses ICU-documented measurements only.** Operating-room (OR) time is a gap for ICU vitals in Layer C. **Labs (approved Nov 11, 2026)** enter every window at their release time, CHARTTIME + 1 h; Layer B windows use the same rule. Labs during surgery are **not** excluded from Layer C. The OR-window lab count below is descriptive only.

## Inputs (read only; do not modify)
1. **Frozen canonical encounter CSV**, `data/interim/canonical_icu_encounters.csv`.
   - Compute its SHA-256 **before reading** and compare it exactly (all 64 characters) with `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418`. On mismatch, raise an error and stop.
   - Required columns: `encounter_id`, `SUBJECT_ID`, `hadm_id_list`, `HADM_ID`, `INTIME`, `OUTTIME`.
2. **Harmonization fields**, `data/interim/encounter_harmonization_fields.csv`. The SHA-256 must be exactly `ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d`.
   - Use t0 and `OUT` (last OUTTIME) as given. Do not rebuild the 24 h stay merge.
3. **Split:** the inherited subject-grouped split. Expected the frozen 70/15/15 split (sizes omitted from this public copy) encounters; the frozen canonical cohort (counts omitted from this public copy) overall.
4. **Raw PIC tables, read with `usecols`** (chunked for the large tables). These are the actual headers (raw metadata, Nov 10, 2026). The script asserts that they are present by exact name.
   - **`SURGERY_INFO.csv`** (expected 7,488 rows): `ROW_ID`, `SUBJECT_ID`, `HADM_ID`, `VISIT_ID`, `OPER_ID`, `ANES_START_TIME`, `ANES_END_TIME`, `SURGERY_BEGIN_TIME`, `SURGERY_END_TIME`, `SURGERY_NAME`, `ANES_METHOD`, `SURGERY_POSITION`.
     - Expected missingness: `ANES_START_TIME` and `ANES_END_TIME` are each blank in 1 row; every other column is complete.
     - `SURGERY_NAME` has 1,496 distinct values. **There is no native type column.**
     - `ROW_ID` is not used.
   - **`SURGERY_VITAL_SIGNS.csv`** (expected 1,944,187 rows): `ROW_ID`, `SUBJECT_ID`, `HADM_ID`, `VISIT_ID`, `OPER_ID`, `ITEM_NO`, `MONITOR_TIME`, `ITEMID`, `VALUE`.
     - There are no missing values. `ITEMID` has 20 distinct codes (`SV1`, `SV2`, …).
     - **Column names:** the docs name the time column `MONITORTIME`, but the local header is `MONITOR_TIME`. Resolve it from the header by exact match against either name, and record which one was found. The docs type `VALUE` as INT; the audit checks whether it is integer-valued (Outputs §8).
     - **Row count:** the docs give 1,216,011 rows, which is the v1.0 count. The local v1.1.0 file has 1,944,187 rows, because v1.1.0 added SpO2, EtCO2 and ECG ST change. The audit expects 1,944,187 and also reports rows by ITEMID, so the v1.1.0 additions are visible.
     - **Documented semantics:**
       - readings are taken every 5 minutes during surgery;
       - `ITEM_NO` is the order of the reading within the surgery;
       - `OPER_ID` is the order of the surgery within the hospitalization;
       - `VISIT_ID` is how many times the patient has been admitted.
     - **Links (docs):** PATIENTS on `SUBJECT_ID`, ADMISSIONS on `HADM_ID`, and **D_ITEMS on `ITEMID`**.
   - **`ADMISSIONS`:** `SUBJECT_ID`, `HADM_ID`, `ADMITTIME`, `DISCHTIME`.
   - **`ICUSTAYS`:** `SUBJECT_ID`, `HADM_ID`, `ICUSTAY_ID`, `INTIME`, `OUTTIME`. Apply the Stage 1 valid-stay rule.
   - **`CHARTEVENTS`:** `HADM_ID`, `ICUSTAY_ID`, `ITEMID`, `CHARTTIME`, `VALUE`. Use the v12 contract ITEMIDs for HR, RR, SBP/DBP/MAP, temperature and SpO2. Values are read only to test numeric presence, and no value is written.
   - **`LABEVENTS`:** `HADM_ID`, `CHARTTIME` only (for the OR-window lab count).
   - **`PRESCRIPTIONS`:** `HADM_ID`, `STARTDATE`, `ENDDATE` and the drug-name columns, restricted to the v12 vasopressor group (for intraoperative vasopressor availability).
   - **`D_ITEMS`** (dictionary table, not patient data): `ROW_ID`, `ITEMID`, `LABEL_CN`, `LABEL`, `LINKSTO`, `CATEGORY`, `UNITNAME`.
     - The docs say it links to SURGERY_VITAL_SIGNS on ITEMID, and that heart rate is both ITEMID 1003 (CHARTEVENTS) and SV1 (SURGERY_VITAL_SIGNS).
     - The official demo treats **SV3 as systolic BP**.
     - D_ITEMS is the **primary source** of the SV map (see Inputs §5).
5. **Config maps (written by Vinai, not data):**
   - `config/surgery_type_map_DRAFT.csv`: see "Two-stage gate".
   - `config/surgery_vital_itemid_map_DRAFT.csv`: columns `ITEMID`, `label_cn`, `label`, `linksto`, `category`, `unit`, `plausible_low`, `plausible_high`, `source` (D_ITEMS / PIC docs / inferred), `confirmed` (Y/N), `notes`.
     - **Primary source: D_ITEMS.** Copilot pre-fills the map by joining the 20 observed SV ITEMIDs to D_ITEMS on `ITEMID` and copying the dictionary fields (`LABEL_CN`, `LABEL`, `LINKSTO`, `CATEGORY`, `UNITNAME`), with `source = D_ITEMS` and `confirmed = N`.
       - D_ITEMS is a dictionary, so its labels may be printed in full.
       - SV codes with no D_ITEMS row are listed with empty labels and `source = none`.
     - Documented anchors are noted in `notes`: SV1 = heart rate (same concept as CHARTEVENTS 1003; PIC d_items docs), and SV3 = systolic BP (official demo).
     - `plausible_low` and `plausible_high` stay empty for Vinai.
     - **Fallback and consistency check:** the per-ITEMID value distributions (Outputs §8) are used to infer the meaning of any unlabeled code, and to flag labels that don't fit their value ranges. Copilot never changes a label; it only flags.
     - **Vinai confirms** each row (`confirmed = Y`).
   - The SHA-256 of each map is pinned once Vinai approves it (`TYPE_MAP_SHA256`, `SV_MAP_SHA256` in the script's config block).

## Two-stage gate
**Stage 1, name inventory** (can run while the maps are stubs):
- Copilot may create `config/surgery_type_map_DRAFT.csv` as a **stub** with the columns `SURGERY_NAME`, `n_train`, `n_overall`, `standardized_name`, `surgery_type`, `notes`.
  - One row per raw name with at least 10 surgeries overall, sorted by `n_overall`.
  - Names with fewer than 10 overall are **not listed**. They form one pooled row, `RARE`, which gives only the number of rare names and their pooled counts.
  - `standardized_name` (for merging spelling and wording variants, e.g., VSD repair vs VSD closure, and for translating if needed), `surgery_type` and `notes` are left **empty**.
  - `n_train` follows primary and complementary suppression (`<10`, `S`) like any other cell.
- Stage 1 writes only the stub and the name-inventory section (Outputs §2) to `data/interim/surgery_item_level_audit_stage1.json` / `.md`.
- The stub and the Stage 1 files must pass the same recoverability self-check and the disk-reload release gate.
- **Copilot must not fill `standardized_name` or `surgery_type`, propose types, or use keyword rules to assign them.**

**Stage 2, full audit:**
- Runs only after Vinai has filled the type map manually and its SHA-256 is pinned. The script raises an error if the hash is unpinned or doesn't match.
- Every listed name must have a non-empty `surgery_type`. The `RARE` row's surgeries are reported as type `unmapped_rare`, and their share is reported.
- The SV map may be pre-filled from D_ITEMS in Stage 1, with `confirmed = N`.
- **Sections that depend on the SV ITEMID map** (intraoperative summaries such as lowest MAP, and the intraoperative part of cluster support) run only after the SV map is approved and its hash pinned. Until then, they write `"pending_sv_map"`, and only the per-ITEMID distributions (Outputs §8) are reported.

## Definitions
All times are in hours unless stated. Let t0 = harmonized first INTIME and `OUT` = harmonized last OUTTIME of the encounter.

- **Times:**
  - AS = `ANES_START_TIME`, AE = `ANES_END_TIME`;
  - B = `SURGERY_BEGIN_TIME`, E = `SURGERY_END_TIME`.
- **Windows:**
  - **Anesthesia window** [AS, AE] defines the **OR gap in Layer C** and the **OR-window lab check** in LABEVENTS. If AS or AE is missing, fall back to B or E, and count every fallback (`anes_fallback_start`, `anes_fallback_end`).
  - **Surgery window** [B, E] defines the **duration** D = E − B and the duration bands.
  - **Induction-to-incision** = B − AS. **Emergence** = AE − E.
  - **Ordering check:** AS ≤ B ≤ E ≤ AE. Count each violation separately: AS > B, B > E (including E = B), E > AE. Count any surgery with at least one violation.
- **Encounter link:**
  - A surgery belongs to an encounter if its `SUBJECT_ID` matches and its `HADM_ID` is in the encounter's `hadm_id_list`. Then assign the timing class.
  - `VISIT_ID` is **not** an HADM equivalent. The docs define it as how many times the patient has been admitted (an admission count), so it is never used to link. It is audited only as an admission-order consistency check (Outputs §5).
  - **v1.3: a surgery's key is `SURGERY_INFO.ROW_ID`** (Vinai, Nov 11, 2026). It is asserted nonblank and unique on all 7,488 rows, and all rows are kept.
    - `OPER_ID` is documented as the order of the surgery within the hospitalization, but it is **not unique** there. Copilot's diagnostic found 389 (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) groups covering 889 rows; these are distinct surgeries with different times. `OPER_ID` is a descriptor and the `SURGERY_VITAL_SIGNS` group key only.
    - Report shared-`OPER_ID` groups and possible true duplicates (same subject, HADM, begin time and name). Nothing is dropped.
  - **One surgery, one encounter (v1.3):** each `ROW_ID` is assigned to at most one encounter by the timing rule below, so the timing classes reconcile to 7,488.
  - Never link by time alone.
- **Timing class (v1.3: by anesthesia end T** = AE, falling back to E, counted as `timing_fallback_surgery_end`). Classes are mutually exclusive, with one encounter per surgery, assigned in this order:
  1. `invalid_time`: T missing (AE and E both missing or unparseable). Counted and excluded.
  2. **Within an encounter window:** t0 − 24 ≤ T ≤ `OUT` for a candidate encounter.
     - `recent_layer_a`: T ≤ t0 + 6 (Layer A only). Sub-count `recent_ending_after_t0`: T in (t0, t0 + 6].
     - `during_icu_layer_c`: t0 + 6 < T ≤ `OUT`. Sub-counts:
       - `in_stay`: T inside a raw valid ICU stay [INTIME, OUTTIME];
       - `in_gap`: T inside a within-encounter gap of 24 h or less. This is a descriptor only and stays in Layer C.
  3. `remote_same_hadm`: otherwise, assigned to the **later** encounter (the first with t0 > T). T < t0 − 24, same HADM.
     - T before ADMITTIME is kept and flagged as `pre_admit_time`.
  4. Also report `layer_c_or_recent_reused_as_later_layer_a_history`: assigned surgeries whose HADM also belongs to a later encounter of the subject.
  5. `post_icu_same_hadm`: T > the last candidate's `OUT`, same HADM. Not used as a feature, since it is outcome-period. Counted only.
  6. `other_hadm`: the surgery links to a cohort subject but to a HADM outside the encounter's `hadm_id_list`. Excluded by decision. Report it split into `earlier_hadm` and `later_hadm`.
- **Straddling surgeries:** B < t0 ≤ E, or AS < t0 ≤ AE. Counted. Their class follows T (v1.3).
- **Bands:** under 2 h (D < 2), 2–4 h (2 ≤ D ≤ 4), over 4 h (D > 4).
  - **Candidate cut grid:** 1, 1.5, 2, 3, 4, 6 h. These are counts only. Grid counts are cumulative; see the suppression rules for cumulative differencing.
- **Rare-name pooling:** a raw name is rare if it has fewer than `RARE_MIN` surgeries. Report `RARE_MIN` = 10, 20 and 30, both in train and overall.
  - After the type map is filled, also report the pooling on `standardized_name`.
- **Pre-incision physiology**, two candidate sources:
  - **(a) OR monitor:** `SURGERY_VITAL_SIGNS` readings in [AS, B), any ITEMID. After the SV map is approved, also report by mapped variable.
  - **(b) Last ICU vitals before AS:** the last ICU CHARTEVENTS value of HR, SBP/MAP, SpO2, RR or temperature in [AS − w, AS), for lookback w = **2 h** and **6 h**. Report per item and for "all of HR + BP + SpO2".
    - Source (b) is meaningful only for `during_icu_layer_c` surgeries, since the patient must be in the ICU before going to the OR. For the other classes, it is reported as not applicable.
- **Vitals → surgery assignment (v1.3, approved Nov 11, 2026).** `SURGERY_VITAL_SIGNS` carries (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) but no `ROW_ID`. Each vitals row is assigned to the `SURGERY_INFO` `ROW_ID` **within its key group** whose anesthesia window [AS, AE] contains `MONITOR_TIME`. A missing AS or AE falls back to B or E, as for the windows above, and is counted.
  - Count, by timing class and split, separately for unique-key groups and shared-`OPER_ID` groups:
    - `assigned` (exactly one window);
    - `ambiguous` (two or more windows contain the time; overlapping surgeries);
    - `unassigned` (no window).
  - Ambiguous and unassigned rows are **excluded** from every per-surgery intraoperative statistic and counted.
  - **Sensitivity:** the same assignment with the window widened by ±30 min. Report how many rows change status.
  - Rows whose key matches no `SURGERY_INFO` row are counted as before.
- **Monitor-time duplicates:** two or more `SURGERY_VITAL_SIGNS` rows **assigned to the same `ROW_ID`** (v1.3; was the `OPER_ID` key), with the same `ITEMID` and the same monitor time.
  - *Identical duplicates:* all `VALUE`s equal.
  - *Conflicting duplicates:* `VALUE`s differ.
  - **Dedup rule** (applied before every per-ITEMID statistic, coverage and pre-incision count): keep the row with the **highest `ITEM_NO`** (the latest-recorded reading). Conflicts are counted and reported, never silently resolved.
  - The official demo drops duplicate monitor times within a surgery and ITEMID, keeping the first row. The spec instead keeps the highest `ITEM_NO`, which is Vinai's rule.
  - The pre-dedup counts are reported alongside.
- **Intraoperative vasopressor availability:** a vasopressor-group PRESCRIPTIONS record with STARTDATE in [AS, AE]. Availability only.

## Outputs
Files:
- Stage 1: `data/interim/surgery_item_level_audit_stage1.json` / `.md`, plus the type-map stub.
- Stage 2, shared: `data/interim/surgery_item_level_audit.json` / `.md` (aggregate-only, suppressed).
- Local only: `data/interim/surgery_item_level_audit_raw_LOCAL.json` (unsuppressed aggregates). It is never copied, uploaded or shared.

Suggested script: `scripts/audit_surgery_item_level.py`, with `--stage 1|2`. The JSON top level carries the input and map hashes, the asserted headers, the stage and the suppression summary.

All count items are reported by **train / validation / test / overall** unless stated otherwise. A surgery's split is its encounter's split. Unlinked surgeries are reported overall only.

1. **Table presence and size:**
   - the row counts, which must equal 7,488 and 1,944,187;
   - per-column missing counts, which must match Inputs §4;
   - surgeries linked to cohort encounters, and surgeries not linked (with reasons: HADM not in cohort, subject not in cohort, HADM not in the encounter's list);
   - surgeries per encounter (quantiles: median, IQR, p90; never min or max), and the share of encounters with any surgery in each timing class.
2. **Name inventory (Stage 1 and Stage 2):**
   - The number of distinct names (raw, and after trim, case-fold and whitespace normalization). Expected raw: 1,496.
   - **Long-tail profile:** the number of names covering 50%, 80% and 90% of surgeries (overall and train). Also the number of names, and the share of surgeries, at or above `RARE_MIN` = 10 / 20 / 30.
   - Names with at least 10 overall, with counts by split (the stub carries train and overall).
   - Stage 2 adds: counts by `surgery_type` and by `standardized_name`, the `unmapped_rare` share, and the number of raw names merged per standardized name (quantiles).
3. **ANES_METHOD and SURGERY_POSITION:** each one's distinct values and counts, suppressed.
   - Values with fewer than 10 overall are pooled (`other`) and not listed.
   - Also their cross-tab with timing class.
   - Both are known at the start of surgery and are **candidate inputs** for the type map and the clusters (Vinai decides).
4. **Timing, windows and validity:**
   - Timing-class counts and sub-counts from Definitions, including:
     - `recent_ending_after_t0`, `in_stay`/`in_gap`, `timing_fallback_surgery_end`;
     - straddling, `pre_admit_time`, `earlier_hadm` and `later_hadm`;
     - the later-encounter reuse count.
     - For `recent_layer_a` and `remote_same_hadm`, give quantiles of t0 − T.
   - Missing or unparseable AS, AE, B and E; `anes_fallback_start`, `anes_fallback_end`.
   - Ordering violations: AS > B, B > E, E > AE, any violation.
   - D > 24 h (implausible; counted, not trimmed).
   - B or E outside [ADMITTIME, DISCHTIME].
   - Quantiles (median, IQR, p10, p90) of induction-to-incision (B − AS) and emergence (AE − E), by timing class, for valid ordered records.
5. **Linkage:**
   - Surgeries within cohort encounters by `HADM_ID` + timing: counts by timing class and split (this is §4's table, reconciled to linked surgeries).
   - **VISIT_ID as an admission count** (consistency check only; not a link):
     - how often `VISIT_ID` equals the rank of the HADM by ADMITTIME within the subject (in the full ADMISSIONS table);
     - how often it is constant within (`SUBJECT_ID`, `HADM_ID`);
     - how often it agrees between `SURGERY_INFO` and `SURGERY_VITAL_SIGNS` for the same surgery key.
     - Equality with `HADM_ID` is **not** tested.
   - **`OPER_ID` semantics:**
     - shared (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) groups and their rows in `SURGERY_INFO`. This is reported, not expected to be 0: the Nov 11 diagnostic found 389 groups and 889 rows. Also report possible true duplicates;
     - how often `OPER_ID` runs 1..n consecutively within each hospitalization;
     - how often the `OPER_ID` order agrees with the order of B within the hospitalization;
     - the number of `OPER_ID` values reused across hospitalizations, which is expected and informational only.
     - Also identical name + B + E within an encounter.
   - **Monitor join (v1.3):** surgeries (`ROW_ID`) with at least 1 assigned `SURGERY_VITAL_SIGNS` row, using the key-group-plus-anesthesia-window assignment in Definitions. Report the assigned / ambiguous / unassigned counts (primary and ±30 min) and:
     - the number of monitor rows matching no surgery;
     - the number of monitor surgery keys absent from `SURGERY_INFO`;
     - how often `ITEM_NO` increases with monitor time within a surgery;
     - all by timing class and split.
     - `OPER_ID` alone, or `OPER_ID` + `SUBJECT_ID`, is **not** a valid key.
   - IDs are never printed; only counts and rates.
6. **Duration**, for valid ordered B and E:
   - quantiles (median, IQR, p10, p90; never min or max), by timing class, and in Stage 2 by type (types with at least 10 in the group);
   - band counts for the proposed bands, by type × split (Stage 2), and counts at each cut of the candidate grid;
   - **Band support flags:** type × band cells with train < 10, and type × band cells with validation or test < 10. These are flags for Vinai, not merges.
7. **Training-split median duration by type** (Stage 2; the "expected duration" for the start bin):
   - Computed from **train only**, and only for types with at least 10 valid train durations. Types below that use a pooled fallback (the train median of all surgeries); report the number of such types.
   - For context, report the median absolute difference between actual duration and the type median in validation. Validation is not used to set anything.
8. **SV monitor table:**
   - **Per ITEMID** (all 20, reported whether or not the SV map exists):
     - row count and the number of surgeries with at least 1 reading;
     - the share of numeric `VALUE`s;
     - percentiles p1, p5, p25, p50, p75, p95 and p99 (computed only when n ≥ 10; never min or max);
     - shares ≤ 0, and the shares above 50, 100, 200 and 300. These are generic plausible-range flags; Vinai sets the real ranges in the SV map.
     - **D_ITEMS labels** (`LABEL`, `LABEL_CN`, `CATEGORY`, `UNITNAME`, `LINKSTO`), printed next to each code. These are dictionary fields.
     - **Integer check:** the share of integer-valued readings. The docs type `VALUE` as INT, so non-integer values are flagged.
     - **Label–range consistency flag:** for codes with a D_ITEMS label, flag the label if its p5–p95 falls outside a broad reference range for that concept. Broad ranges are given in the script config, e.g., heart rate 40–250 /min, SBP 30–200 mmHg, SpO2 50–100%, EtCO2 10–80 mmHg; other concepts are unchecked.
       - SV1 is also compared with the CHARTEVENTS 1003 percentiles among the same subjects (aggregate percentiles only).
       - Flags are for Vinai; nothing is relabeled.
     - **Duplicate monitor times:** identical and conflicting counts per ITEMID (before dedup), and the share of surgeries with any conflict.
     - the sampling interval after dedup (median and IQR of within-surgery gaps), compared with the documented 5 minutes: the share of gaps equal to 5 min, under 5 and over 5;
     - readings per surgery (quantiles);
     - `ITEM_NO` value counts per ITEMID (pooled below 10).
     - The number of SV codes found in D_ITEMS, and the number missing from it.
   - **Coverage per surgery:**
     - the fraction of the surgery window [B, E] with readings, defined as 5-minute slots containing at least 1 reading divided by total slots (quantiles);
     - the same fraction for the anesthesia window.
   - **Monitor readings outside the anesthesia window:** counts and shares before AS and after AE, and the share of surgeries with any such reading.
   - **After the SV map is approved:** per mapped variable (e.g., MAP, HR, SpO2, EtCO2), coverage, and the availability of the candidate intraoperative summaries (lowest MAP, desaturation, vasopressor availability, duration). Before approval, write `"pending_sv_map"`.
9. **Pre-incision physiology availability:**
   - Source (a): the share of surgeries with at least 1 reading in [AS, B), and readings in that window (quantiles). By mapped variable after SV-map approval.
   - Source (b): the share with each item, and with HR + BP + SpO2, at w = 2 h and 6 h. Also the median time from the last ICU value to AS (quantiles). For `during_icu_layer_c` surgeries only.
   - Report both by timing class and split, plus the overlap: the share with (a) only, (b) only, both and neither.
10. **OR-gap check**, using the anesthesia window [AS, AE] (with fallback), for `during_icu_layer_c` surgeries:
    - The share of surgeries with **zero** ICU CHARTEVENTS rows in (AS, AE). Expected to be high. Also counts of rows.
    - **LABEVENTS rows in (AS, AE).** LABEVENTS covers the whole hospital stay, so labs drawn in the OR could enter Layer C by time. This is a **flag for Vinai**, who decides whether OR-window labs are excluded from Layer C.
11. **Cluster input support** (Stage 2, for φ2, in train only):
    - Per type, the number of surgeries with complete inputs: type, a pre-incision source ((a) or (b); report each), and, after SV-map approval, the intraoperative BP + HR + SpO2 variables.
    - Flag types with fewer than 30 complete surgeries; these would be pooled in the cluster fit.
    - ANES_METHOD and SURGERY_POSITION completeness, as candidate inputs.
    - Validation and test complete counts are reported for context only.
    - **No clustering is run.**
12. **Landmark context**, at core landmarks 24, 48, 72, 120, 168, 240 and 336 h, among encounters eligible at L (`OUT` > t0 + L):
    - encounters with at least 1 `during_icu_layer_c` surgery that started (B) before t0 + L;
    - surgeries **in progress** at t0 + L (B < t0 + L ≤ E). For these, the end-bin information is not yet revealed;
    - anesthesia windows spanning t0 + L;
    - surgeries that ended before t0 + L;
    - `during_icu_layer_c` surgeries whose anesthesia window spans more than one 6 h bin.
    - Eligible counts must reproduce the band-support audit: [landmark eligible counts omitted] overall at the core landmarks.
13. **Markdown summary:**
    - timing classes;
    - the long-tail profile;
    - ANES_METHOD and SURGERY_POSITION;
    - window validity and ordering;
    - linkage;
    - duration and bands (Stage 2: by type), and the train median by type;
    - the SV per-ITEMID table;
    - coverage, pre-incision availability, the OR-gap check, cluster support and landmark context.

## Suppression and disclosure (PhysioNet DUA)
- **Aggregate counts only.** No row-level output, no IDs, no timestamps of individual rows, no individual values, no example rows.
- **Raw-count reconciliation, done before suppression and not written out.** Raise an error unless all of the following hold:
  - timing classes plus `invalid_time` sum to linked surgeries;
  - linked plus unlinked surgeries equal the `SURGERY_INFO` rows (7,488);
  - band counts plus invalid durations equal the valid-B surgeries in each group;
  - type counts (Stage 2) sum to linked surgeries, with `unmapped_rare` included;
  - ANES_METHOD and SURGERY_POSITION value counts sum to their totals;
  - per-ITEMID rows sum to the `SURGERY_VITAL_SIGNS` rows (1,944,187);
  - train + validation + test equals overall for every cell.
- **In plain language:**
  - *Primary suppression* hides any small count, from 1 to 9.
  - *Complementary (secondary) suppression* hides some larger counts as well, so a hidden count can't be worked out by subtraction from published totals, from other splits, or from cumulative counts.
- **Primary suppression (mandatory):** any count from 1 to 9 is published as `"<10"`. Zero is also published as `"<10"`; nothing is structural.
- **Complementary suppression (mandatory):** published as `"S"`. It must cover:
  - **margins**, including the type-map stub (listed names + `RARE` = total, for train and overall);
  - **cross-split** (train + validation + test = overall);
  - **cumulative differencing** (the cut-grid counts, the long-tail profile and any `n(≤ k)` style count);
  - **shares:** a published share with a published denominator is treated as its count.
- **Quantiles and medians:**
  - Computed only on groups of at least 10 (raw). Smaller groups publish `"<10"` in place of the statistic.
  - Never min or max.
  - Durations and gaps are rounded to 0.1 h. SV value percentiles are rounded to the item's native precision.
- **Names and categories:** only names, ANES_METHOD values and SURGERY_POSITION values with at least 10 surgeries overall are printed. Pooled values are never listed.
- **Recoverability self-check (mandatory):** for every suppressed cell, apply every published linear identity listed above. Raise an error if any suppressed cell is exactly recoverable or bounded inside [1, 9]. The JSON records the number of cells checked and the number recoverable (must be 0). An LP check is an optional extra.
- **Release gate (mandatory, both stages):**
  - After writing the shared JSON and, in Stage 1, the type-map stub, the validator **reloads the published files from disk** and re-runs reconciliation of the published cells, the marker checks and the recoverability self-check.
  - The markers must be exactly `"<10"` (primary) and `"S"` (complementary). A complementary cell must not be labeled `"<10"`.
  - It prints `release gate ok` only if every check passes.
  - No shared file is copied, uploaded or shared before this.
- **Logs:** never print subject, encounter, admission, visit or operation IDs, or individual timestamps or values. Counts only.

## Acceptance checks (the script must raise on failure)
1. Canonical or harmonization SHA-256 that is not an exact full-hash match.
2. Any expected header missing from `SURGERY_INFO` or `SURGERY_VITAL_SIGNS` (Inputs §4), named by column only. The monitor-time column must resolve to exactly one of `MONITOR_TIME` or `MONITORTIME`. Row counts other than 7,488 and 1,944,187, or missingness that differs from Inputs §4. Differences are reported as counts before raising.
3. Cohort totals or split counts that differ from Inputs §3, or any subject in two splits.
4. Core-landmark eligible counts that don't reproduce the band-support audit (Outputs §12).
5. Any raw reconciliation failure.
6. Timing classes that are not mutually exclusive and exhaustive over valid-T linked surgeries, or any `ROW_ID` assigned to two encounters, or classes that do not reconcile to 7,488. Unit test.
7. Boundary unit tests on synthetic rows (v1.3; T = anesthesia end):
   - T = t0 − 24 → `recent_layer_a`; t0 − 24.0001 → `remote_same_hadm`; t0 + 6 → `recent_layer_a` (and never Layer C); t0 + 6.0001 → `during_icu_layer_c`; `OUT` → `during_icu_layer_c`; `OUT` + 0.0001 → `post_icu_same_hadm`.
   - T in a merge gap → `during_icu_layer_c`, with `in_gap` incremented.
   - T between two encounters of one HADM → the later encounter.
   - AE missing → E is used, with the fallback counted.
   - Other HADM → `other_hadm`.
   - D = 1.999 → under 2 h; 2 → 2–4 h; 4 → 2–4 h; 4.0001 → over 4 h.
8. Any surgery from another hospitalization counted as remote. Unit test: a surgery in an earlier HADM goes to `earlier_hadm`, never to `remote_same_hadm`.
9. **Window roles.** The OR gap and the OR-window lab check must use [AS, AE], and duration and bands must use [B, E]. Unit tests:
   - a synthetic surgery with AS < B and E < AE puts a lab at AS + 0.01 inside the OR window, but outside the surgery window;
   - a synthetic surgery with a missing AS falls back to B, and `anes_fallback_start` increments.
10. **Ordering.** Each violation type is counted, and a unit test with one record per violation type gives counts of exactly 1 each. Violating records are counted, not silently dropped.
11. Training medians computed with validation or test rows. Unit test: changing validation or test durations leaves every published median unchanged.
12. Any quantile computed on a group under 10, or any min or max published.
13. Any printed name, ANES_METHOD value or SURGERY_POSITION value with fewer than 10 overall; or any ID, timestamp or individual value in the JSON, markdown, stub or logs. Unit test with a sentinel ID.
14. **Type-assignment guard:**
    - Stage 1 writes the stub with `standardized_name` and `surgery_type` empty, and contains no type-assignment or keyword-typing code.
    - Stage 2 raises an error if `TYPE_MAP_SHA256` is unpinned or doesn't match, or if any listed name has an empty `surgery_type`.
15. **SV-map guard.** No intraoperative summary (such as lowest MAP) and no mapped-variable coverage is computed unless `SV_MAP_SHA256` is pinned and matches. Otherwise those fields are `"pending_sv_map"`. The per-ITEMID distributions (Outputs §8) are always present for all 20 ITEMIDs.
16. **Pre-incision sources.** Both (a) and (b) are reported, with (b) at 2 h and 6 h, and with the (a)/(b) overlap table. (b) is reported as not applicable outside `during_icu_layer_c`.
17. **Linkage items present:** the VISIT_ID admission-count check, the `ROW_ID` uniqueness assertion, the shared-`OPER_ID` group counts, the vitals assignment counts (assigned, ambiguous, unassigned; primary and ±30 min), and coverage fractions.
18. **Long-tail profile present:** the names covering 50%, 80% and 90%, and pooling at 10, 20 and 30.
19. Primary suppression, complementary suppression (margins including the stub, cross-split, cumulative) or the recoverability self-check missing or failing.
20. Release gate not passed in either stage. The disk-reload validator must report `release gate ok`.
21. A raw `_LOCAL` file written outside `data/interim/`, or referenced in the shared JSON by anything other than its name.
22. Any feature, cluster or model output written.
23. **SV map from D_ITEMS.** The map stub contains all 20 observed SV codes, with D_ITEMS fields copied verbatim where they exist, `source` set, and `confirmed = N`. No label is edited by code. Label–range flags and the SV1 vs CHARTEVENTS 1003 comparison are present.
24. **Keys (v1.3).** The script raises in any of these cases:
    - `ROW_ID` is not unique or is blank;
    - any join uses `OPER_ID` without `SUBJECT_ID` + `HADM_ID`;
    - `VISIT_ID` is used as a link key;
    - a vitals row is attached to a surgery without the anesthesia-window assignment.

    Unit tests:
    - two hospitalizations that each have `OPER_ID = 1` stay separate;
    - two distinct surgeries sharing (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) stay separate;
    - their vitals rows go to the right `ROW_ID`;
    - an overlap gives `ambiguous`;
    - a reading outside both windows gives `unassigned`.
25. **Monitor-time duplicates.** Identical and conflicting duplicate counts are present. The dedup keeps the highest `ITEM_NO`, and a unit test with a conflicting pair keeps the higher-`ITEM_NO` value and increments the conflict count. Every downstream statistic uses the deduplicated rows.
26. **VALUE integer check** is present per ITEMID.

## Post-audit decisions (Vinai; not automated)
1. **Surgery type map:** `standardized_name` merges and `surgery_type` for every listed name; how `unmapped_rare` surgeries are handled; whether ANES_METHOD or SURGERY_POSITION informs the type. Then pin the hash.
2. **SV ITEMID map:** confirm the D_ITEMS labels (anchors: SV1 heart rate, SV3 systolic BP); resolve any label–range flags; label any codes missing from D_ITEMS using the distributions; set plausible ranges. Then pin the hash. Intraoperative summaries are built only after this.
3. **Pre-incision source:** (a) OR monitor, (b) last ICU vitals (2 h or 6 h), or both with a priority rule.
4. **Final duration bands:** keep under 2 h / 2–4 h / over 4 h, or choose from the cut grid, possibly type-specific if support differs. Also the handling of implausible durations (> 24 h) and ordering violations.
5. **Cluster inputs (φ2):** type, pre-incision summaries, intraoperative summaries, and whether to include ANES_METHOD or SURGERY_POSITION; how to pool thin types.
- Also: the rare-name threshold for φ0 (10, 20 or 30), and how to handle conflicting monitor-time duplicates if the highest-`ITEM_NO` rule looks wrong.
  - Two items are **decided Nov 11, 2026**: the B = t0 edge is moot (anesthesia-end timing), and OR-window labs are **not** excluded from Layer C.

## Limitations
- Surgery type depends entirely on Vinai's manual map. Rare names (fewer than 10 overall) are not listed, so they cannot be typed individually and stay `unmapped_rare` unless Vinai adds a rule.
- The SV labels come from D_ITEMS and are unconfirmed until Vinai approves the map. Distribution-based inference for unlabeled codes is a judgment, not a data fact.
- Intraoperative vasopressor need may not be recorded at the dose level. The audit reports availability only.
- "Same hospitalization" uses `HADM_ID`. For the 52 multi-HADM encounters, any HADM in `hadm_id_list` counts as the same.

## Out of scope: future candidate sources
- **`OR_EXAM_REPORTS`** (docs: 183,809 rows; `SUBJECT_ID`, `HADM_ID`, `EXAMTIME`, `REPORTTIME`, `EXAM_ITEM_TYPE_NAME`, `EXAM_ITEM_NAME`, `EXAM_PART_NAME`): hospital-wide exams and imaging performed during the stay. The exam notes are not included. Despite the "OR" prefix, it is not a surgery table.
  - It is **out of scope** for this audit.
  - It is a future candidate **Layer C event source**, released at `REPORTTIME` (not `EXAMTIME`), to avoid look-ahead. It would need its own item-level audit spec.

## Do not
- Do not modify the frozen CSV, the harmonization fields, the labels or the split file.
- Do not assign surgery types or SV labels, build features, fit clusters or train models.
- Do not use validation or test data to set any threshold, median or band.
- Do not print or save any row-level data, outside the local-only raw file of aggregate counts.

## Changelog
- **v1.3 (Nov 11, 2026): approved decisions (Vinai).**
  - **Key:** `ROW_ID`. `OPER_ID` is not unique within an admission (389 groups, 889 rows). All rows are kept.
  - **Timing by anesthesia end** (fallback: surgery end, counted):
    - recent = [t0 − 24 h, t0 + 6 h], Layer A only;
    - Layer C = (t0 + 6 h, `OUT`], merge gaps included, with `in_gap` as a descriptor;
    - remote < t0 − 24 h, same HADM, kept and flagged before `ADMITTIME`;
    - between encounters → the later encounter;
    - one surgery per encounter in the audit;
    - reuse as Layer A history of a later same-HADM encounter is stated and counted.
  - `begin_equals_t0` is removed (moot).
  - **Vitals join:** each row goes to the `ROW_ID` in its (`SUBJECT_ID`, `HADM_ID`, `OPER_ID`) group whose [AS, AE] contains `MONITOR_TIME`. Ambiguous and unassigned rows are counted and excluded; ±30 min is the sensitivity. Monitor-time dedup is per assigned `ROW_ID`.
  - **Labs:** release time (CHARTTIME + 1 h) everywhere, including Layer B. OR-window labs are not excluded from Layer C.
  - Acceptance checks 6, 7, 17 and 24 are updated.
- **v1.2 (Nov 10, 2026):** aligned with the public PIC documentation and the official demo.
  - **SV map:** built primarily by joining D_ITEMS on ITEMID (dictionary labels printed; anchors SV1 = heart rate, SV3 = systolic BP). Value distributions become the fallback and a consistency check. This replaces v1.1's check of only whether SV codes appear.
  - **IDs:** `VISIT_ID` is an admission count, not an HADM equivalent, and is never a link key. `OPER_ID` is unique only within `SUBJECT_ID` + `HADM_ID`. `ITEM_NO` is the reading order.
  - **Monitor-time duplicates:** identical and conflicting counts, and a dedup rule that keeps the highest `ITEM_NO`.
  - **Headers:** `MONITOR_TIME` vs the documented `MONITORTIME` is resolved from the header, and `VALUE` is checked for integer values.
  - **Row count:** noted (docs 1,216,011 for v1.0; local v1.1.0 has 1,944,187 after the SpO2, EtCO2 and ECG ST additions).
  - **OR_EXAM_REPORTS:** listed as an out-of-scope future Layer C source, released at REPORTTIME.
  - Acceptance checks 23–26 added.
- **v1.1 (Nov 10, 2026):** updated to the actual `SURGERY_INFO` and `SURGERY_VITAL_SIGNS` headers, row counts and missingness, and removed assumed columns and the D_ITEMS labeling assumption. Changes:
  - **Two-stage gate:** a Stage 1 name inventory with a type-map stub; Stage 2 runs after Vinai fills the map and its hash is pinned. Copilot never assigns types.
  - **Windows:** anesthesia times define the OR gap and the OR-window lab check; surgery times define duration. Fallback counts, induction-to-incision and emergence distributions, and ordering checks added.
  - **Pre-incision physiology:** two sources, (a) OR monitor in [AS, B) and (b) last ICU vitals at 2 h and 6 h.
  - **ANES_METHOD and SURGERY_POSITION:** counts added.
  - **SV ITEMID map:** stub, per-ITEMID distributions, and a guard on intraoperative summaries.
  - **Linkage:** expanded.
  - **Name cardinality:** long-tail profile added.
  - Acceptance checks 9–10 and 14–18 added; the post-audit decision list updated.
- **v1 (Nov 10, 2026):** initial spec, following Vinai's final surgery decisions (log Step 5.2; amendment decision 12).
