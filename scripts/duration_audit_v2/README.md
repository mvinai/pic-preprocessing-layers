# duration_audit_v2: surgery complexity vs duration audit (clean reference implementation)

This folder implements `specs/surgery_complexity_duration_audit_spec.md` **v1.7** and `surgery_complexity_duration_audit_FIXES_v1.md`, including:
- **Amendment A** (signed D2);
- **Amendment B** (Nov 11, 2026: `ROW_ID` surgery key; `OUTTIME` accepted as the harmonized OUT);
- **Amendment C** (Nov 11, 2026: surgery timing by anesthesia end). It replaces Copilot's run-1 script. Do not reuse anything from that script.

| file | purpose |
|---|---|
| `audit_surgery_complexity_duration.py` | the audit (a single file; the config block is at the top) |
| `test_audit_surgery_complexity_duration.py` | pytest suite (synthetic data only) |
| `synthetic_data.py` | synthetic PIC-like generator used by the tests (no real data). Since v1.6 it plants distinct surgeries that share an `OPER_ID`, plus possible true duplicates. |
| `diagnose_surgery_key_duplicates.py` | aggregate-only duplicate-key diagnostic, kept for the record; not needed for the run |

## 0. Run order for Copilot (Nov 11, 2026)
1. **pytest** (synthetic only):
   - in this folder: `python -m pytest -q test_audit_surgery_complexity_duration.py` → `128 passed`;
   - in `scripts/`: `python -m pytest -q test_check_encounter_harmonization_consistency.py` → `9 passed`.
2. **Harmonization check** (read-only, aggregate-only): `python scripts/check_encounter_harmonization_consistency.py --icustays ... --surgery-info ... --surgery-vitals ... --json-out data/interim/harmonization_consistency_check.json`. For the full command, see `docs/encounter_harmonization_consistency_review.md` §4. Send back its console output.
3. **Duration audit:** §1 below. Send back the console output and the shared JSON/MD.

## 1. How Copilot runs it (Windows)

1. Copy `audit_surgery_complexity_duration.py` into the project's `scripts/` folder. The project root is taken as the parent of `scripts/`, the same layout as run 1. Inputs are then found at:
   - `data/interim/canonical_icu_encounters.csv`
   - `data/interim/encounter_harmonization_fields.csv`
   - `config/surgery_complexity_map_v1.csv`
   - `<PIC_RAW_DIR>\SURGERY_INFO.csv` (and `ADMISSIONS.csv`, optional)
   - `<SPLIT_FILE>`
2. Install the dependencies (Python ≥ 3.9): `pip install pandas numpy scipy`. To run the tests as well: `pip install pytest`.
3. Run the tests first. They use synthetic data only:
   `python -m pytest -q test_audit_surgery_complexity_duration.py` (copy `synthetic_data.py` next to it). Expected result: `128 passed`. The suite also passes under `LC_ALL=C`, and console and file I/O are UTF-8 on Windows.
4. Run the audit from the project root:
   ```
   python scripts/audit_surgery_complexity_duration.py
   ```
   Any path can be overridden with a CLI flag or an environment variable:

   | flag | env var |
   |---|---|
   | `--project-root` | `FOUNDATIONS_ROOT` |
   | `--raw-dir` | `PIC_RAW_DIR` |
   | `--split-file` | `SPLIT_FILE_PATH` |
   | `--canonical` | `CANONICAL_PATH` |
   | `--harmonization` | `HARMONIZATION_PATH` |
   | `--map` | `COMPLEXITY_MAP_PATH` |
   | `--out-dir` | `AUDIT_OUT_DIR` |
   | `--split-col` | `SPLIT_COL` (default `split`) |
   | `--bootstrap-reps` | `BOOTSTRAP_REPS` (default and spec value 2000; do not lower it for the real run) |

5. **Harmonized ICU times.** The script uses t0 (harmonized first INTIME) and OUT (harmonized last OUTTIME) from `encounter_harmonization_fields.csv`. Each column is resolved from these candidates, and exactly one must match:
   - t0: `t0, T0, t0_harmonized, harmonized_t0, first_intime, FIRST_INTIME, t0_intime`
   - OUT: `OUT, out, last_outtime, LAST_OUTTIME, harmonized_out, OUT_harmonized, t_out, last_OUTTIME`. **If none of these is present, it falls back to `OUTTIME` (v1.6).**

   The real file has `t0` and `OUTTIME` and no `OUT`, so **no flags are needed**: the console shows `t0='t0', OUT='OUTTIME'`. `--harm-t0-col t0 --harm-out-col OUTTIME` is equivalent and still accepted, and an explicit flag always wins. If `span_icu_hours` is present, the console prints `OUT - t0 vs span_icu_hours (>1 min) disagreements: N`. Expect 0; a nonzero value means the column is not the last (max) OUTTIME, so report it. If the names ever differ, the run stops with `cannot resolve the harmonization T0/OUT column ... header has [...]`. Rerun with the flags and **report the names used**.
6. **Surgery key (v1.6).** Each `SURGERY_INFO` row is one surgery, keyed by `ROW_ID`. The run stops (exit 2) if `ROW_ID` is blank or repeated, including `123` vs `123.0`.
   - `OPER_ID` is only a descriptor. In PIC it is not unique within an admission: 389 (SUBJECT_ID, HADM_ID, OPER_ID) groups cover 889 rows, and these are distinct surgeries.
   - All 7,488 rows are kept.
   - The console line `surgery key ROW_ID (nothing dropped): ...` reports the shared-OPER_ID groups and rows, and the possible true duplicates (same SUBJECT_ID, HADM_ID, SURGERY_BEGIN_TIME and normalized name). Nothing is removed.
7. Do not edit the config block: the pinned hashes, the criteria (`TAU_B_SIGNED_MAX = 0.2`, signed) or the thresholds. The script refuses to run if the criteria differ from the approved values. A hash mismatch stops the run before any file is read.

Runtime on the real cohort is roughly 3–6 minutes. Most of it is the 2,000-replicate SUBJECT_ID-cluster bootstrap and the recoverability solver.

## 2. Exit codes and outputs

| exit | meaning |
|---|---|
| 0 | ok; shared JSON + MD written after passing the release gate |
| 2 | input / hash / header / cohort / reconciliation failure (`AUDIT STOPPED: ...`); nothing shared is written |
| 3 | release gate failed (`RELEASE GATE FAILED - no shared file left in place: ...`); all shared files are deleted |

The outputs go to `data/interim/`:
- `surgery_complexity_duration_audit_MAP_v1.json`: **shared**, suppressed, with `release_gate: "ok"`.
- `surgery_complexity_duration_audit_MAP_v1.md`: **shared**, a human-readable rendering of the JSON. The gate checks that it is exactly that rendering.
- `surgery_complexity_duration_audit_MAP_v1_raw_LOCAL.json`: **LOCAL ONLY**. It holds unsuppressed counts. Never upload, paste or share it.

## 3. Expected console output (shape; numbers are from a synthetic run)

```
input hash checks: canonical OK, harmonization OK, complexity_map OK, split OK (all 4 exact 64-char)
SURGERY_INFO headers ok (12/12), row count ok, missingness ok, ROW_ID unique and nonblank
cohort and split totals ok
reconciliation ok (NNNN identities verified on raw counts; link classes sum to SURGERY_INFO rows)
suppression computed: NNN primary, NNN complementary, NNNN checked, 0 recoverable

=== surgery_complexity_duration_audit v1.7 | MAP_v1 (reviewed 2026-11-10) | release_gate=ok ===
input hashes: canonical OK (e5feeb22b545...), harmonization OK (ffd1b8ae6a13...), complexity_map OK (5c993b88adfc...), split OK (60196792457a...)
harmonization columns used: t0='t0', OUT='OUTTIME'; ADMISSIONS read: True; surgery key: ROW_ID
OUT - t0 vs span_icu_hours (>1 min) disagreements: 0
surgery key ROW_ID (nothing dropped): shared_oper_id_groups=..., rows_in_shared_oper_id_groups=..., possible_true_duplicate_groups=..., rows_in_possible_true_duplicate_groups=...
reconciliation: surgery_info_rows=[omitted], linked=..., post_icu_same_hadm=..., other_hadm=..., not_in_cohort=..., invalid_time=..., not_linkable_missing_encounter_time=...
multi-candidate: n_multi_candidate=..., resolved_rule1_during=..., resolved_rule2_next_t0=..., resolved_rule3_post_icu=...
pooled population: linked=... analysed=... unmapped=... excluded_unclassified=... invalid={...} fallback_anes=... d_gt_24h=...
  Cardiac        assessable   D1=not_met D2=not_met D3=not_met D4=not_met D5=not_met big=N  tau_b=0.617 ci=[0.59, 0.64] d3_share=... medians=[...]
  <one line per body system> ...
  POOLED                      D1=...     D2=...     D3=...     D4=...     D5=...     big=...  stratified_tau_b=...
suppression: primary=... complementary=... checked=... recoverable=0
bootstrap: replicates=2000 seed=20261110 runs=... realised=...-...
release gate ok
```

Count cells print as an integer ≥ 10, as `<10` (primary suppression, including zero) or as `S` (complementary suppression). The console prints **only** published (suppressed) values, so it is safe to paste back.

## 4. What to send back

1. The full console output, copied verbatim, including any error and the exit code.
2. `surgery_complexity_duration_audit_MAP_v1.json` and `surgery_complexity_duration_audit_MAP_v1.md` (shared files only).
3. The pytest summary line.
4. If you used `--harm-t0-col/--harm-out-col` or any path override, the exact command line.

**Never send** the `_raw_LOCAL.json`, any row-level data, IDs or timestamps.

## 5. Design notes (summary)

- **Hashes.** All four pinned SHA-256 values (canonical, harmonization, complexity map, split) are checked on the file bytes before parsing.
- **Parsing.** Every CSV is read as `dtype=str` with NA detection off, so the map's `NA` stays a literal. IDs are normalized (`"123.0"` → `"123"`). `hadm_id_list` is parsed once and handles JSON, Python repr, float repr and delimited strings.
- **Surgery key.** `ROW_ID` (unique, nonblank); `OPER_ID` is a descriptor. Rows sharing an OPER_ID are linked independently.
- **Linkage (v1.7, anesthesia-end timing).** T = `ANES_END_TIME`, falling back to `SURGERY_END_TIME` (counted as `timing_fallback_surgery_end`). One surgery is linked to one encounter, using candidate encounters for the same subject whose `hadm_id_list` contains the surgery's HADM_ID:
  1. If t0 − 24 h ≤ T ≤ OUT, the surgery belongs to that encounter:
     - `recent_layer_a` if T ≤ t0 + 6 h. It is **Layer A only, never also Layer C**.
     - otherwise `during_icu_layer_c`, including merge gaps.
  2. Otherwise it goes to the later encounter (first with t0 > T) as `remote_same_hadm`.
  3. Otherwise it is `post_icu_same_hadm`.

  Also counted: `recent_ending_after_t0`, `straddling_t0`, and `pre_admit_time` (remote surgeries that end before ADMITTIME; kept and flagged).

  The other classes are `other_hadm`, `not_in_cohort`, `invalid_time` and `not_linkable_missing_encounter_time`. All classes must sum exactly to the SURGERY_INFO rows (7,488), or the run stops.
- **Statistics.** All are computed at the surgery level:
  - The **signed** Kendall tau-b (`scipy.stats.kendalltau(..., variant="b")`) on train surgeries, with a 2,000-replicate SUBJECT_ID-cluster bootstrap (fixed seed).
  - NaN tau → `not_assessable`, never `not_met`. Groups < 10 and NC are reported, never dropped.
  - D1–D5 as in the spec; `big_discrepancy` = Y iff D1, D2 or D3 is met. Label: `approved_criteria_2026-11-10_D2signed`.
  - Sensitivities a–e, flagged vs unflagged, the across-system check, the pooled stratified tau-b and the descriptive all-rows block.
- **Disclosure.**
  - Suppression applies only to declared count fields. Counts 0–9 → `<10`.
  - Complementary `S` cells are chosen by an algorithm. A linear-programming solver checks every published identity (row, column, split, pooled and table-to-table differences). No suppressed cell may be exactly recoverable, and no hidden difference may be bounded inside [1, 9].
  - The release gate re-reads the JSON/MD from disk and repeats every check independently.
