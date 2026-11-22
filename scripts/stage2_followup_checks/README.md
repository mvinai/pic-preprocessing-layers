# Stage 2 follow-up checks (Copilot run instructions)

Read-only, aggregate-only follow-up to the Stage 2 audit. Spec: `feature_sparsity_frequency_audit_spec.md` v0.9, §5.4 (decisions D24 and D29).
- Script: `stage2_followup_checks.py` v1.0, SHA-256 `e955974d8e9fd85414482549cc14e019d922215505088220ae9cfba990992f62` (41,003 bytes).
- Tests: `test_stage2_followup_checks.py` (synthetic data only; needs `synthetic_stage2_data.py` from the Stage 2 folder).

## What it does
- **(i) Bilirubin 5255 vs 5075 (D24).** Train split, analysis window, lab availability = chart time + 1 h.
  - Groups: all train; NICU vs non-NICU (first care unit, as Stage 2); neonate vs not-neonate (age < 28 days; `age_invalid` counts as not-neonate).
  - Per item per group: n values, n encounters, value p10/p25/p50/p75/p90, per-encounter-median quantiles.
  - Per group: KS 5255 vs 5075 on values and on encounter medians (only when both n ≥ 10); encounters with either / both items; same-draw pairs and the same-draw difference 5255 − 5075.
- **(ii) FiO2 5222 format (D29).**
  - Train-window FiO2 rows by band: ≤ 1, (1, 21), [21, 100], > 100, plus nonnumeric.
  - Same-draw pO2 (5239 preferred) / numeric FiO2 pairs. Percent rule P/F = pO2 / (FiO2/100): in gate, below 10, above 700, FiO2 ≤ 0.
  - How many out-of-gate pairs fall in the gate with FiO2 ×100 or ÷100.
  - P/F under "FiO2 ≤ 1 is a proportion, otherwise a percentage".
  - pO2 and FiO2 quantiles for in-gate vs out-of-gate pairs; pairs and out-of-gate pairs by FiO2 band.
  - The Stage 2 v1.0.1 pair definition (FiO2 gated 21–100) as a reproduction check against the v1.0.1 valid and outside-gate pair counts (pinned in the script; not quoted here, V29).
- **Gates (exact SHA-256, checked before use; any mismatch stops with exit 2):**
  - `feature_sparsity_frequency_audit.py` v1.0.1 `3b5b36fbd0f73420beaa74b18abf2afec4725cf86c5e8efe514be19120043213` (imported, not run);
  - `audit_surgery_complexity_duration.py` v1.7 `ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0`;
  - canonical, harmonization, split and complexity map (the Stage 2 pins); cohort, split and eligibility totals must reproduce;
  - `data\interim\maps\feature_channel_map_v2.csv` `cbfdc68f11bc2f8d9f96150803a1339930da81eb091f481576f709debdd5ab5c` (686 rows);
  - LABEVENTS (Stage 1 hash) and, per D36, PATIENTS `85b016fc…9f336ae0d`, ADMISSIONS `1451ae8e…2647e1ae`, ICUSTAYS `449119e6…95e90eea1b`, SURGERY_INFO `04450bff…0002a38d83` (full values in the script).
- **Disclosure:** the Stage 2 machinery. Counts 1–9 → `<10`; complementary `S` until the LP recoverability check finds nothing recoverable; quantiles and KS only beside published counts ≥ 10; no min/max, IDs or timestamps. Pending → gate → final → gate.

## Outputs (`data\interim\`)
| file | share? |
|---|---|
| `stage2_followup_checks.json` | **TEAM-INTERNAL — never share** |
| `stage2_followup_checks.md` | rendered from the JSON; treat as team-internal until separately reviewed |

**`stage2_followup_checks.json` is TEAM-INTERNAL and is never shared.** The v1.3 suppression check-only (Step 5, November 16, 2026) failed the release gate on 3 FiO2 `rows_5222` numeric-band cells (pointers only, no values): `/fio2_format/rows_5222/numeric_by_band/n_gt_100`, `/fio2_format/rows_5222/numeric_by_band/n_gt_1_lt_21`, `/fio2_format/rows_5222/numeric_by_band/n_le_1`. Decision (Vinai, November 16, 2026). This does not affect the Stage 2 v1.0.2 release, which passed. Do not modify the pinned script `stage2_followup_checks.py`.

No unsuppressed (`_raw_LOCAL`) file is written. Do not open, copy or attach any existing `*_raw_LOCAL*` file.

## Run order for Copilot
See the precheck and the run command in the hand-off message (they are the same as below). Run from the project root
`<PROJECT_ROOT>`. Do not edit any file.
Stop on any failure (exit 2 `AUDIT STOPPED`, exit 3 `RELEASE GATE FAILED`) and send back the full console text only.
