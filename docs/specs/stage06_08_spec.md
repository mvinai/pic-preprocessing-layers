---
title: "Stages 06–08 builder spec (stage06_08_spec v1.2.1)"
subtitle: "PIC Foundations — plan v13, contract v13 — Layer A levels, Layer B levels, train-only preprocessing"
date: "Nov 18, 2026"
---

# 1. Summary

One script, `06_08_build_levels_preprocess.py` (v1.2.1), runs stages 06, 07a/07 and 08 of the plan v13 stage map (§13.2) in one pass on the Stage 3 v1.0.4 release, plus one test file and four pinned published reference tables (decision D2) that sit next to the script. **Copilot runs: 3** (tests, precheck, build). Expected real build time 3–6 min on the laptop (allow 10; hard limit `--max-seconds 3600`).

```
06_08_build_levels_preprocess.py       SHA-256
f55cdaf5d7dd51470a9bd9494b1ed37730fcc29c7948de29eb8461d33daec8e7
test_06_08_build_levels_preprocess.py  SHA-256
617f3ca258d70801fb2d248f7b5f62d8efa9b5c738c4fcda90f2d6dc549fbdc9
fleming2011_hr_centiles.csv            SHA-256
4ff3c424dd398ded5b7e9a7a70de45f3b44fff8f9fe456a2a3be0e61d331f1ad
fleming2011_rr_centiles.csv            SHA-256
ec54810ee50120ad3ea5c69c90a5f7cbd6c3ef09de560d102786f8de7c94e552
roberts2020_map_percentiles.csv        SHA-256
421c965f262335e2445674a4fe78c9fa756b10d28b64a206eb3c38b57f1ef2f2
podium2022_baseline_scr.csv            SHA-256
e56de3e8f9340ada0d163e520c58dbf657a21fddfc179ca0213f184ca79190e7
```

Box testing: 43 tests pass on the synthetic PIC-format Stage 3 fixture (the Stage 3 v1.0.4 test recipe; 55 s on the box). Scale run: synthetic 4,646-encounter Stage 3 fixture (one third of the real cohort): build 27 s on the box, peak memory 1.2 GB, 130 of 130 acceptance checks pass (all three admission-source levels occur; A12 and A13 pass; reference tables verified; D2 columns present; a variant with injected blank categories and drug names containing "tab_"/"mab_" also passes 130 of 130), release gate ok, 43 files (96 MB). No real data was opened or run on the box.

**What each stage produces** (folder `data\processed\ladder\stage06_08\`, 43 files, published atomically from a `stage06_08.pending-<pid>` folder):

* **Stage 06 (Layer A levels).** `layer_a_a0.csv`, `layer_a_a1.csv`, `layer_a_a2.csv` (+ `.meta.json` each): one row per canonical encounter (count omitted from this public copy), `encounter_id`, `split`, then the level's features before preprocessing (empty cells where a value is not defined, each covered by a mask/flag column).
* **Stage 07 (Layer B levels; 07a absorbed).** `layer_b_b0.csv`, `layer_b_b1.csv`, `layer_b_b2s.csv`, `layer_b_b2d.csv`, `layer_b_b2p.csv` (+ `.meta.json`).
* **Stage 08 (train-only preprocessing).** `stage08_<level>.npz` for the 8 tabular levels (`X_imputed`, `X_z` primary set with one reference category dropped per block, `X_full_z` full one-hot set, `columns`, `kinds`, `split`, `encounter_id`, `eligible_24h`); `stage08_layer_c_L<L>.npz` for the 14 landmarks (`X_raw_z` with NaN where missing, `X_filled_z` from the Stage 3 landmark-specific filled copy, `missing`, `carried_forward`, `imputed`, `channels`); `stage08_fitted_transforms_INTERNAL.json` (all train statistics; team-internal).
* **Audit.** `stage06_08_audit.json` / `.md` (shareable: widths, column and category counts, check results, hashes; **no patient counts at all**, so no cell needs suppression — the release gate refuses any number outside a structural whitelist and re-checks the file after reloading it from disk); `stage06_08_category_report_INTERNAL.csv` (train-subject counts per category; team-internal, never shared); `stage06_08_manifest.json` (SHA-256 of every file).

# 2. Inputs and hash gates (precheck; nothing parsed before these pass)

**v1.2.1 (Nov 18, 2026; fixes after the first real run, which stopped safely with nothing written on two checks):** (1) **B8** tested the tokens `org_`, `ab_`, `dilution`, `interpret` as raw substrings of every column name, so real drug names in the data-driven `b0_drug__<name>_<hash>` columns (e.g. "…_tab_…", "…mab_…") tripped it; the tokens are now tested only on field stems — at the start of the prefix-stripped stem or after an underscore — and for one-hot / drug category columns only on the block part before `__`, never on category text; the audit lists, by block name only, the blocks whose category text contains a token. (2) **A6** failed because real Layer A has blank category values (a blank row is all-zero in its block); following D11 (missing as its own indicator), every one-hot block whose source has blanks (for `recent_`/`remote_` blocks, blanks among active rows) gets a `<block>_missing` mask in every level holding the block (`admit_interval_band` reuses the existing `admit_interval_missing`); A6 now checks block sum + mask = 1 (or the active indicator) and fails only on an uncovered row. **Sex is exempt** (D35 removed `sex_missing`, a banned column): a blank sex stays an all-zero row, accepted only on blank-source rows. Masks constant in train are dropped by the stage 08 zero-variance rule; A4 compares the sex and first-care-unit blocks plus their masks across a0–a2. Layer B has no one-hot block (drug columns are multi-hot flags; no B check sums them). (3) **Diagnostics:** a failed acceptance check now prints, per check, the level and column or block names involved (never row or patient counts). On blank-free synthetic data every level file and tensor is byte-identical to v1.2.0.

**v1.2.0 (Nov 18, 2026):** decision D2 — published age norms replace the empirical train age-band z for heart rate and respiratory rate (Fleming 2011), mean arterial pressure (Roberts 2020) and creatinine (ratio to the PODIUM 2022 baseline); the four reference tables are pinned files checked first in the precheck (any difference → `EXIT:2`, nothing read or written); Haque–Zaritsky 2007 is a team-internal MAP sensitivity only. Stage 3 (including its KDIGO staging) is unchanged.

**v1.1.0 (Nov 18, 2026):** decision D1 — the single-X admission-source rule is replaced by Vinai's split rule (surgical / floor / transfer_ed, fixed settings 48 h and 24 h recorded in the script, audit, metadata and manifest); `--adm-src-x-hours` is removed; check A12 now verifies the first-match rule and A13 verifies that the admission-source block is not a copy of the interval band.

* Reference tables (decision D2, §9) — `scripts\stage06_08\fleming2011_hr_centiles.csv`, `fleming2011_rr_centiles.csv`, `roberts2020_map_percentiles.csv`, `podium2022_baseline_scr.csv`: full SHA-256 pinned in the script (values in §1 and Step 0), checked before anything else, plus a shape check (14 / 13 / 21 / 10 rows; 13 centile columns for Fleming and Roberts).
* `data\processed\ladder\stage3\stage3_manifest.json` — pinned on the command line by its full SHA-256 (`--expect-stage3-manifest-sha256`, value printed by the precheck). The manifest must name builder `03_build_real_layers v1.0.4` with SHA-256 `d60fbdc0bda5a01b1237c17e7c55fbc646dcec327132de0d4ff5f151ee80ffd3`, release gate ok, every acceptance check pass (only `train_median_available_all_channels` may be `noted`), channel map v2.3 `05c1f917…97b1fe` and φ1 medians `e73c8b5a…db0047` (full hashes in the script).
* Every one of the 46 Stage 3 files is re-hashed and must equal its manifest entry; no extra or missing file is allowed in the folder.
* Canonical encounter file (`data\interim\canonical_icu_encounters.csv`; columns `encounter_id`, `SUBJECT_ID` only — `Died` is never read), harmonization file (`data\interim\encounter_harmonization_fields.csv`; `t0` and `OUT`/`OUTTIME` only, for the stay length computed exactly as Stage 3 does) and split file (`..\PIC_Analysis_3_0\data\processed\subject_splits.csv`) must equal the hashes recorded in the Stage 3 manifest; canonical and split must also equal the plan pins `e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418` and `60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31`.
* After loading: one row per canonical encounter; Layer A base 297 and Layer B base 667 columns; split per subject equals the Layer A base split; the frozen 70/15/15 split (sizes omitted from this public copy); the stay lengths reproduce the 14 landmark eligible counts (values omitted from this public copy); Layer C raw rows in Layer A order.
* The output folder must be empty or absent and no `stage06_08.pending-*` folder may exist.

The precheck prints the SHA-256 of each reference table, of the manifest, of each Stage 3 file used, and of canonical, harmonization and split.

# 3. Stage 06 — Layer A levels (plan §4)

Common to all levels: sex and first care unit one-hot (identical in every level, check A4). Rare categories (fewer than **m = 20 distinct training subjects**, counted on train only) are pooled; diagnosis codes pool within their (chapter, domain) cell, then to a global `other`.

* **a0 (φ0):** `age_months` + `age_invalid`; admission-to-ICU hours + missing flag; ICU admission hour as sin/cos; day-of-week one-hot; primary diagnosis 3-character code (with `missing`); recent and remote surgery standardized names; recent surgery T→t0 hours; pre-ICU 6 h labs (value + mask); medication groups active at t0.
* **a1 (φ1):** age bands; interval bands; shift + weekend; diagnosis chapter; early-surgery flag, recent body system, complexity group, φ1 duration band (ordinal + mask), remote body system; **admission source** (decision D1, Nov 18, 2026), one-hot with levels surgical / floor / transfer_ed, **first match wins**: (1) *surgical* = a recent surgery ended ≤ t0 (Stage 3 `adm_src_or_component`) **and** hospital-to-ICU interval < 48 h; (2) *floor* = not surgical **and** interval ≥ 24 h; (3) *transfer_ed* = everyone else (ED and outside transfer cannot be separated in PIC: EDREGTIME equals ADMITTIME almost always, Stage 3 decision 9). The interval is Stage 3's `adm_src_interval_h` (hospital-to-ICU hours, clipped at 0). **Missing interval** follows the existing handling: a missing interval never counts as a long interval, so it satisfies "< 48 h" and fails "≥ 24 h" — surgical when a qualifying surgery exists, otherwise transfer_ed (the interval missing flag stays in a0/a1). The two settings are fixed, recorded constants (surgical window 48 h, floor cut-off 24 h), not command-line options; pre-6 h labs as train age-band z-scores, except **creatinine**, which becomes the ratio to the PODIUM 2022 age- and sex-specific baseline (`a1_pre6_creatinine_ratio_podium` + `_missing` mask; decision D2).
* **a2 (φ2):** 4-level timing; diagnosis domain; surgery body system only (clusters deferred, §8); pre-6 h labs with organ-block names (creatinine again as the PODIUM ratio, `a2_pre6_renal__creatinine_ratio_podium` + mask).

Checks: A1 rows, A3 no banned column, A4 fixed fields identical, A5 bands partition, A6 one-hot blocks sum to 1 with the block's `_missing` mask (v1.2.1), A7 vocab train-only (tested), A9 early-surgery flag, A10 pooling (group sums), A12 admission source (three levels partition the cohort and every row follows the first-match rule), A13 duplicate-column check (no admission-source indicator equals an interval-band indicator, the interval ≥ 24 h indicator, or the complement of either), A16 weight not a feature.

# 4. Stage 07 — Layer B levels (plan §5.5–5.8; window [0, min(24 h, stay)))

* **b0 (φ0):** vitals median / min / max / range; labs first and worst (worst = direction table D7; two-sided labs take the value farther from the train age-band median); triggered tests worst; qualitative results; panel masks; indicator-only flags (ALT, TSH); hemoglobin and bilirubin source flags; drug-name multi-hot (pooled at m = 20).
* **b1 (φ1):** vital z (median / min / max) — HR and RR against Fleming 2011, MAP against Roberts 2020, each with a `<vital>_agez_missing` mask (value missing or age outside the table's cover); other vitals train age-band z; HR / RR / MAP high/low flags on the published z; variability = SD of published-z 6 h bin medians over bins 0–3 (≥ 3 observed bins, else empty + mask); creatinine first and worst as the PODIUM ratio (`creatinine_first_ratio_podium`, `creatinine_worst_ratio_podium` + `creatinine_ratio_podium_missing`); other labs z (first, worst; train age band); triggered worst z; urine output and fluid intake (raw); KDIGO UO stage (provisional); BUN/Cr ratio; AKI stage; 13 medication groups × 7 features + max concurrent vasoactives; microbiology M1 (any culture, culture count, hours to first culture, specimen families).
* **b2s / b2d / b2p (φ2):** b1 regrouped into organ blocks: single home (b2s), duplicated with `__dup__` copies (b2d), and the scored variant (b2p: four pSOFA components + input masks + unscored systems + drug block). Modified PRISM III in all three. Drug groups under 1 % train prevalence dropped above φ1.

Checks: B1, B3 (+ no raw vitals outside b0), B5, B6, B8 (M1 only; field stems, v1.2.1), B11, B13, B20, B21, B23, B24 (PRISM identical; four pSOFA components), B25 (no P/F, S/F, FiO2, A-a gradient, anion gap), B26 (urine output never back-filled), every empty cell covered by a mask/flag.

# 5. Stage 08 — train-only preprocessing (plan §7)

Per tabular level, fitted on train rows only and applied unchanged to validation and test: winsorize continuous columns at train p1 / p99 → log1p if non-negative and train skewness > 2 → train-median fill (the mask/flag already carries missingness) → drop zero-variance columns → standardize with train mean / SD. `X_z` drops one reference category per one-hot block (the most frequent in train); `X_full_z` keeps all. B10: the train rank deficiency of b2d minus that of b2s must equal the number of retained duplicates. Layer C: per landmark, per channel, winsorize / log / mean / SD from the observed train cells of that landmark's prefix; `X_raw_z` keeps NaN, `X_filled_z` standardizes the Stage 3 filled copy (one filled copy per landmark with landmark-specific train medians; never a whole-grid fill); a channel with no train variation at a landmark is dropped and listed (e.g. `kdigo_uo_stage` at 24 h). Layer C is never flattened.

# 6. Decisions applied (each listed in the audit)

| # | decision | applied |
|---|---|---|
| D1 | admission source (decided by Vinai, Nov 18) | split rule, first match: surgical (surgery ended ≤ t0 and interval < 48 h) / floor (interval ≥ 24 h) / transfer_ed; fixed 48 h and 24 h; missing interval never counts as long |
| D2 | age-norming (decided by Vinai, Nov 18) | HR, RR: Fleming 2011 normal z; MAP: Roberts 2020 centile-interpolated z (capped 1st/99th), interpolated by exact age; creatinine: ratio to the PODIUM 2022 baseline (age, sex); other labs and vitals: train age-band z; outside cover → empty + mask; details §9 |
| D2u | creatinine units | Stage 3 already delivers mg/dL (PIC µmol/L ÷ 88.4, as in Stage 3); not converted again; run stops if the train median is ≥ 10 mg/dL |
| D2s | Haque–Zaritsky 2007 | MAP sensitivity only, ages 1–17 y, team-internal file; not a feature |
| D3 | vital variability | SD of z-scored bin medians, bins 0–3, ≥ 3 bins, HR / RR / MAP |
| D4 | worst direction | table D7 in the script; two-sided rule vs the train age-band median |
| D5 | pSOFA | 4 components; cardiovascular 2 if any vasoactive; neurologic not scored (no GCS); no total; missing → 0 |
| D6 | PRISM III | non-neurologic points; PT not scored; pO2 may be venous; missing → 0 |
| D7 | b2p composition | pSOFA components + input masks + unscored systems + drug block (+ PRISM) |
| D8 | φ2 surgery | body system (clusters deferred) |
| D9 | hour encoding | sin / cos |
| D10 | prior encounters | omitted |
| D11 | diagnosis missingness | the `missing` category |
| D12 | log rule / empty-by-design fill | log1p if ≥ 0 and skew > 2; train median |
| D13 | Layer C standardization | pooled per channel over observed train cells |
| D14 | triggered-test "ordered" flag | ≡ 1 − mask, so only the mask is kept |

# 7. Open items carried from plan v13 (§18) — not built

* φ2 surgery clusters (inputs and algorithm undefined in the plan) → body system used.
* Microbiology M2 / M3 (Stage 3 has no organism fields) → M1 only.
* φ1-stop variant and the neonatal hepatic sensitivity analysis.
* a2o.
* Published age norms for labs other than creatinine (train age-band z kept; see D2).

**Limitation / future work (creatinine baselines).** Stages 06–08 use the PODIUM 2022 baselines for the creatinine ratio, but Stage 3's KDIGO stage (Layer B `aki_stage`, Layer C) is unchanged and still uses its own baselines. A later Stage 3 rebuild should swap the PODIUM baselines into Stage 3 KDIGO staging so both agree.

**Limitation (age cover).** PIC has no gestational age, so infants under 37 weeks post-menstrual age cannot be identified; the tables are applied from birth by postnatal age. Age missing or invalid, age ≥ 18 y (Fleming, Roberts after its last curve point), or sex not M/F (PODIUM) → empty + mask.

**Limitation / future work (admission source).** `ADMISSION_DEPARTMENT` in ADMISSIONS was never audited. A future Stage 3 audit of that field could separate emergency-department from ward admissions (and possibly from outside transfers), which the current transfer_ed level cannot do.

# 8. Copilot steps (Windows, project root `<PROJECT_ROOT>`)

**Rules.** Run from the project root. Do not edit any file. Run the commands exactly as written (the only allowed change is pasting the manifest hash into Step 3). Stop at the first failure (any `EXIT:` other than 0, or a hash that differs) and send back the full console text only; do not patch, retry with other options, or delete anything. Do not open any input or output data file; the only output file you may open is `data\processed\ladder\stage06_08\stage06_08_audit.md`. Never open or paste a file whose name contains `_INTERNAL` or `_raw_LOCAL`. **Hang rule:** the script prints a progress line or a heartbeat at least every 30 s; if nothing new appears for 5 minutes, press Ctrl+C once and send the console text.

**Step 0 — pause cloud file syncing and copy the package.** *Pause cloud file syncing first* (taskbar cloud-sync icon → Settings/gear → Pause syncing → 2 hours). A real run failed last night with WinError 5 because syncing was not paused. Keep it paused until Step 3 has finished. List leftovers (list only, do not open):

```
Get-ChildItem data\processed\ladder -Force -Filter 'stage06_08*' | Select-Object Name,Mode
```

If anything is listed, stop and report it. Then copy the folder `stage06_08\` unchanged (binary copy) to `scripts\stage06_08\` — it now holds six files: the script, the tests and the four reference tables named below (`.csv`; do not open, edit or re-save them in Excel) — and check the hashes (each must equal the value shown; the Stage 3 files are the unchanged v1.0.4 package the tests reuse):

```
Get-FileHash scripts\stage06_08\06_08_build_levels_preprocess.py -Algorithm SHA256
# f55cdaf5d7dd51470a9bd9494b1ed37730fcc29c7948de29eb8461d33daec8e7
Get-FileHash scripts\stage06_08\test_06_08_build_levels_preprocess.py -Algorithm SHA256
# 617f3ca258d70801fb2d248f7b5f62d8efa9b5c738c4fcda90f2d6dc549fbdc9
Get-FileHash scripts\stage06_08\fleming2011_hr_centiles.csv -Algorithm SHA256
# 4ff3c424dd398ded5b7e9a7a70de45f3b44fff8f9fe456a2a3be0e61d331f1ad
Get-FileHash scripts\stage06_08\fleming2011_rr_centiles.csv -Algorithm SHA256
# ec54810ee50120ad3ea5c69c90a5f7cbd6c3ef09de560d102786f8de7c94e552
Get-FileHash scripts\stage06_08\roberts2020_map_percentiles.csv -Algorithm SHA256
# 421c965f262335e2445674a4fe78c9fa756b10d28b64a206eb3c38b57f1ef2f2
Get-FileHash scripts\stage06_08\podium2022_baseline_scr.csv -Algorithm SHA256
# e56de3e8f9340ada0d163e520c58dbf657a21fddfc179ca0213f184ca79190e7
Get-FileHash scripts\stage3_builder\03_build_real_layers.py -Algorithm SHA256
# d60fbdc0bda5a01b1237c17e7c55fbc646dcec327132de0d4ff5f151ee80ffd3
Get-FileHash scripts\stage3_builder\test_03_build_real_layers.py -Algorithm SHA256
# 84529d8e10bf82f957ab5d3c69e170318ee7593cefb99f6efb26a64d53a4831e
Get-FileHash scripts\stage3_builder\stage3_fixture_pic.py -Algorithm SHA256
# dab2442235d72cbaf44a2f422e847940822dbc7b5d15af41af6aef70b96e5c11
```

**Step 1 — tests (synthetic data only, about 1–3 min).** Expect `43 passed` and `EXIT:0`.

```
$ErrorActionPreference='Continue'; cmd /c "python -m pytest -q scripts\stage06_08\test_06_08_build_levels_preprocess.py 2>&1"; Write-Host "EXIT:$LASTEXITCODE"
```

**Step 2 — precheck (hashes only; writes nothing; under 1 min).** Expect `EXIT:0` and `result: PRECHECK OK (reference tables, Stage 3 manifest, every Stage 3 file, canonical, harmonization, split); nothing written`, with four `reference <name> sha256 …` lines (equal to Step 0), `stage3_manifest.json sha256 <64 hex characters>` and the other hashes.

```
$ErrorActionPreference='Continue'; cmd /c "python scripts\stage06_08\06_08_build_levels_preprocess.py --precheck-only 2>&1"; Write-Host "EXIT:$LASTEXITCODE"
```

**Step 3 — build (about 3–6 min on the laptop (allow 10; hard limit `--max-seconds 3600`)).** Replace `<MANIFEST_SHA>` with the 64-character value after `stage3_manifest.json sha256` in Step 2 (no other change). There is no admission-source option any more (decision D1 is fixed in the script). Expect `EXIT:0` and a last line `result: STAGES 06-08 OK; RELEASE GATE OK`.

```
$ErrorActionPreference='Continue'; cmd /c "python scripts\stage06_08\06_08_build_levels_preprocess.py --expect-stage3-manifest-sha256 <MANIFEST_SHA> 2>&1"; Write-Host "EXIT:$LASTEXITCODE"
```

On `EXIT:2` nothing was written. On `EXIT:3` the console names any pending folder left behind; do not open or delete it — report it. Then resume cloud file syncing.

**Step 4 — report.** Send back the full console text of Steps 0–3 and the contents of `data\processed\ladder\stage06_08\stage06_08_audit.md` only.

# 9. Decision D2 — published age norms (details) and references

* **Heart rate, respiratory rate — Fleming 2011 [1].** Normal z with mean = the published median and SD = (p90 − p10) / (2 × 1.28155), by the published age row containing the exact age. HR: the "Birth" row is used only for age < 1 day; from 1 day the 0–3 month row. RR: no separate birth row is used; birth uses 0–3 months. Cover 0 ≤ age < 18 y.
* **Mean arterial pressure — Roberts 2020 [2] (main table).** At the exact age (days) each of the 13 published centile values is linearly interpolated between the table's ages (flat before the first and after the last point up to 18 y); z is linearly interpolated between the centiles' normal quantiles and capped at the 1st / 99th (|z| ≤ 2.3263).
* **Creatinine — PODIUM 2022 [4], Table 2.** Ratio = serum creatinine (mg/dL) ÷ the baseline for the age band (lower bound inclusive) and sex (M → boys, F → girls). Used in a1/a2 pre-ICU labs and in b1 and above (first, worst). Units: Stage 3 already converts PIC µmol/L ÷ 88.4 to mg/dL; Stages 06–08 do not convert again and stop if the train median is ≥ 10 mg/dL.
* **Haque & Zaritsky 2007 [3].** MAP sensitivity only, ages 1–17 y: z = (MAP − (1.5 × age in years + 55)) / (15 / 1.645), compared with the Roberts z on train rows in the team-internal fitted-transforms file (correlation and mean / SD of the difference; no counts). Not a feature.
* **Other labs and other vitals** keep the empirical train age-band mean / SD.
* Provenance of every table (source URL, download date, SHA-256 of the raw download): `reference_norms\SOURCES.md` in the delivery folder.

**References.**

1. Fleming S, Thompson M, Stevens R, et al. Normal ranges of heart rate and respiratory rate in children from birth to 18 years of age: a systematic review of observational studies. *Lancet.* 2011;377(9770):1011-1018.
2. Roberts JS, Yanay O, Barry D. Age-based percentiles of measured mean arterial pressure in pediatric patients in a hospital setting. *Pediatr Crit Care Med.* 2020;21(9):e759-e768.
3. Haque IU, Zaritsky AL. Analysis of the evidence for the lower limit of systolic and mean arterial pressure in children. *Pediatr Crit Care Med.* 2007;8(2):138-144.
4. Fitzgerald JC, Basu RK, Fuhrman DY, et al. Renal dysfunction criteria in critically ill children: the PODIUM consensus conference. *Pediatrics.* 2022;149(1 Suppl 1):S66-S73.
