# Spec: Stage 2 v1.0.3 release, cross-release leak fix (module v1.5; read-only, aggregate-only)

Nov 17, 2026. PIC Foundations feature matrices. Built and tested on the box with synthetic data only; nothing here was run on PIC data, the real v1.0.1 / v1.0.2 values were never opened (pointer names and label classes only), and no `_raw_LOCAL` file was used.

## 0. Change notes

### v1.5 (Nov 17, 2026): the stalled release run — convergence first, then a heartbeat

**What happened (real run, module v1.4, aggregate console lines only).** The input audit (largest component 7,738 cells, 75 failing cells) was fast. The suppression rounds did not converge: the failing count went 75 → 373 → 19 → 256 while the new suppressions grew to 8,849 and the largest component grew 7,738 → 16,143 → 24,758 → 32,955 → 36,087 cells. Round 5's audit took about 8 s; then the run printed nothing for more than 20 minutes and was stopped (no exit code, nothing written).

**Cause (reproduced on the box with a synthetic mirror; §5).**

1. *The silence.* For a failing cell in a component above 300 cells, v1.4 checks candidates with exact LPs over the **whole** component (uncapped, two LPs per check, several checks per failing cell), and the selection loop printed no line until the round was over. With a 36k-cell component and 256 failing cells this is hours, silently. (The audit itself solves each component once and is fast.)
2. *The non-convergence (the real problem).* v1.4's large-component choice (a) tried the largest values first (strata denominators, totals; they sit in many groups and in the links, so each pick joins components), (b) on no single fix suppressed the whole 3-group ring, (c) after 5 rounds without progress ("ESCALATING") suppressed whole rings for every failing cell, and (d) never checked that a pick could itself be protected: a new `S` cell pinned by a link to a published twin (`layer_a/…/n_<band>` = `strata_denominators/…/n_encounters`) fails in the next round. That is the oscillation and the balloon: each round fixes the old failures and creates new ones in a bigger component.

**What v1.5 changes (selection only; rules, audit and gate unchanged).**

- *Self-safe picks.* A candidate set is accepted only if the failing cell **and every newly suppressed cell** pass (checked exactly). A failing pick is paired with one of its own candidates (a breadth-first chain of at most 4 cells, 60 checks per failing cell). The synthetic stuck case — an `S` stratum pinned by its linked twin — is fixed by the twin plus the twin's partition total.
- *Ordering that fixes and does not merge.* Candidates are ranked (1) by the equation they share with the failing cell that has the fewest unknowns (the equation that pins it), (2) by the size of the other suppressed components the pick would join (0 first), (3) by v1.4's static key (fewest groups, a total first, smallest value). The largest-value shortcut only takes values that join nothing.
- *Bounded runtime.* Candidate checks use a local sub-system (equations within 3 groups of the cells, at most 1,500; fewer equations only widen intervals, so a pick is re-checked by the next full audit, and a proof that a cell cannot be protected stays valid), each LP capped at 20 s (a capped LP counts as "not shown", never as safe), 60 s per failing cell, 900 s per round (cells left over go to the next round, logged). The no-progress fallback and the stronger rule are capped (10 / 40 cells per failing cell, 400 per round).
- *Convergence stop (new exit-3 case).* The run stops cleanly, writing nothing in the output folder, when the new suppressions exceed `--max-new-suppressions` (default 2,500), when a component exceeds `--max-component-cells` (default 20,000), when the failing count has not decreased for 10 rounds, or when the stronger rule makes it worse. It then prints `STOPPING: suppression not converging: …` and a **withhold proposal**: single cells whose withholding removes the remaining failures (a failing `S` cell itself; for a failing `<10` / hidden cell, the published or `S` cells of its own groups), checked with a full audit, pointers only; with `--withhold-proposal FILE` the full list is written as a CSV (outside the output folder; refused if the file exists).
- *`--withhold FILE` (pointer, reason).* Single count cells withheld from the release — removed from the JSON, the registry and the markdown and audited as unpublished cells, exactly like a one-cell dropped subtree. The cell-level alternative to D4's subtree drop (decision D8). Refused for a frozen, `--public` or links-file cell, a pointer that is not a count cell, or a cell inside a dropped subtree. The manifest and provenance record the file's SHA-256 and the number of cells.
- *Heartbeat.* A line `heartbeat: still working - <phase> (<position>); LP solves so far N; time left T s` whenever 30 s pass without a line (`--heartbeat-seconds`), also during a single long LP (HiGHS releases the interpreter lock); a progress line every 60 s inside selection; each round's line ends with `choices …; selection N s`. Counts and positions only, never a pointer or value.

**What does not change.** Rule 3 and the residual check; the exact-bound fix (an exactly-0 hidden cell only from an exact LP); carried suppressions, public cells, frozen cells, links and drops (D1–D7); the independent gate re-solve; "cannot protect" (provably unprotectable cells, exit 3, D4); outputs and their names. When v1.4 converges on small components (≤ 300 cells, look-ahead path) v1.5 makes the same choices.

### v1.4 (Stage 2 v1.0.3, first draft of this spec): Stage 2 internal equalities, public values, the carried release, and the hidden-cell bound fix

**Why.** The Stage 3 builder v1.0.3 audits its release jointly with the published Stage 2 v1.0.2 release (Stage 3 spec §16). Its baseline (Stage 2 v1.0.2 plus public information alone, before any Stage 3 cell) already fails 15 Stage 2 cells and 25 hidden Stage 2 differences. Stage 3 cannot repair those; they are flagged for Stage 2 v1.0.3:

- `/layer_a/dob_breakdown/by_split/train/{n_age_days, n_invalid, status/n_valid}`
- `/layer_a/train/age_band/{n_age_invalid, n_ge18}`
- `/reconciliation/{CHARTEVENTS, MICROBIOLOGYEVENTS}/{n_cohort_admission_time_missing, n_no_cohort_admission, n_total}`
- `/strata_denominators/age_band/train/{adolescent, age_invalid}/{n_eligible_72h, n_encounters}`

The causes, and what v1.4 does about each:

1. **L5, age links (train-age leak).** The v1.0.1 registry does not declare that `layer_a/<split>/age_band/n_<band>` and `strata_denominators/age_band/<split>/<band>/n_encounters` count the same encounters (same `age_band()` array), nor that DOB status *not valid* is the `age_invalid` band. A reader can. With the published bands and the strata children (for example `heart_rate/n_ever_layer_b`, seen in every encounter, whose hidden complement is 0), the train invalid-age, adolescent and ≥ 18 cells follow. **Fix:** the links file v3 declares these equalities (§2), so the audit protects against them.
2. **≥ 18 years is a structural zero.** The strata age-band levels are the bands present in the whole cohort; `ge18` is not one of them, so no encounter in any split is ≥ 18 years. The published level list shows this. **Fix:** `n_ge18` = 0 in all three splits is a public value (decision D1).
3. **L3, public table sizes (event-linkage leak).** `reconciliation/<table>/n_total` is the raw table size, which is public (PIC documentation; Stage 3 `/raw_rows`). v1.0.2 suppressed it for CHARTEVENTS and MICROBIOLOGYEVENTS as a complementary cell, which protects nothing when the value is known. **Fix:** the four totals are public values in the audit, so the `<10` not-in-stay classes are protected with the total known.
4. **Invalid-DOB count.** The shared plan v13 quotes the v1.0 run-stop count of invalid-DOB cohort encounters. A reader may equate it with v1.0.1 `layer_a/dob_breakdown/n_invalid` (decision D2). **Fix:** that total is a public value in the audit (team-internal file; the number is never written into any output, log or document), and so are the two cells it fixes exactly with the frozen cohort size (overall valid DOBs and `n_age_days`, rows marked `derived`; the module derives them itself).
5. **Hidden-cell exact-bound bug (module v1.3).** The multi-equation solver screens each component with a few LPs whose solutions are feasible points: an *inner* approximation of each cell's interval. Rule 3 is monotone in the interval (a wider true interval still passes) except one clause: a hidden or dropped cell is safe when it is *exactly 0*. v1.3 accepted "exactly 0" from the screening points alone; when all points happen to sit at 0, a hidden difference whose true range is, say, 0–3 was reported safe. **Fix:** `screen_safe()` never accepts the "exactly 0" clause from points; such cells get the exact per-cell LP maximum. All other clauses still use the points (they are sound). Counterexample (in the tests): H, H, H, S, S cells with v0 + v1 − v3 + v4 = 14 and v0 + v1 + v2 + v3 = 13; v1.3 gives v2 = (0, 0) "safe", v1.4 the exact (0, 3), which fails. The scaled synthetic gives the same outputs with v1.3 and v1.4 (§5), so the fix costs only a few extra LPs.

**What does not change.** Rules 1–3, the residual check, the release gate and every v1.0.2 decision (the NMB drop, the frozen cells) are the same. No published value changes. v1.0.3 suppresses everything v1.0.2 suppressed (`--carry-json`, §3) and withholds the same subtree, except the CHARTEVENTS and MICROBIOLOGYEVENTS table sizes, which v1.0.2 suppressed and v1.0.3 publishes again with their v1.0.1 values (public by definition; decision D3).

**New options in v1.4** (all optional; the v1.0.2 command still gives the v1.0.2 behaviour, apart from the bug fix):

- `--public FILE` (columns `pointer, value, reason`): cells whose value a reader already knows. `value` empty: the input's published value; an integer: the asserted value; `derived`: the module derives it from the other known values through the identities and links (exact LP with all unknowns free; never from the `<10` / `S` bounds); `publish`: the input's value, and the cell is published even where the carried release suppressed it (used for the two table sizes, decision D3). Public cells keep their release label; in the audit they are known constants. The file is team-internal: outputs record only its SHA-256 and row count, never a pointer or value.
- `--carry-json FILE` (+ `--expect-carry-sha256`): the previous release. It must carry `release_gate: ok` and a `release_v1_0_2.source_json_sha256` equal to the input's SHA-256; every value it publishes must equal the input's, every `<10` must be `<10` in the input, and every cell it suppressed or withheld stays so (the gate checks it again). Not allowed with `--check-only`.
- Output names, provenance key and notes are v1.0.3: `release_v1_0_3`, `notes.differencing_fix_v1_0_3`, `notes.dropped_subtrees_v1_0_3`.

### Decisions applied (Vinai, Nov 17, 2026) and records

All seven open items of the v1.4 draft are decided (§7). What changed in the package: the public file publishes the CHARTEVENTS and MICROBIOLOGYEVENTS table sizes again (`publish` rows, D3), so the release command needs no extra option; the "cannot protect" rule (D4), the handling of the public file (D6) and the records note with the v1.0.2 deletion request (D7, §8) are written down. The module, the links, the frozen cells and the drop list are unchanged.

### Earlier versions

v1.3 (hash-pinned drop list), v1.2 (convergence guard), v1.1 (performance) and v1.0: see the v1.0.2 spec (`stage2_v1_0_2_suppression_spec.md`, v1.3). Their behaviour is kept; the 54 v1.3 tests are part of this suite, and the regression tests compare v1.4 with the delivered v1.3 (`reference/stage2_v1_0_2_suppress_v1_3.py`, `df77234b…`).

## 1. Purpose

Re-release Stage 2 so that no suppressed Stage 2 cell is exposed by Stage 2 alone plus public information (the Stage 3 baseline), and so that Stage 3 v1.0.4, re-pinned to this release, finds no compromised Stage 2 cell. The input is the shared v1.0.1 output (values) with the v1.0.2 release carried (its suppressions stay).

## 2. Inputs

| input | file | SHA-256 |
|---|---|---|
| Stage 2 v1.0.1 JSON (values) | `data\interim\feature_sparsity_frequency_audit.json` | `0a2123f48cb6bfbdcf2c0c98d8a886ea23cc9e4b05cf4bd605c31a7e7e2556d0` |
| Stage 2 v1.0.1 markdown | `data\interim\feature_sparsity_frequency_audit.md` | `9fe84830ece57ac4c03151db40907eb866d096f59f91e1ddd2b26aaedcc912ad` |
| Stage 2 v1.0.2 release (carried) | `data\interim\stage2_v1.0.2\feature_sparsity_frequency_audit_v1.0.2.json` | `a801e41ceb70a6bc2e56a5214a785bce31e35911536da3aadffb2281f162ac18` |
| module v1.5 | `stage2_v1_0_3_suppress.py` (137,025 bytes) | `63d903d1ceef80a61f257a08aff4acd4a8c7fd5c2361d1e36de1c0cf603dc751` |
| links v3 (89 rows) | `stage2_v1_0_3_links.csv` | `a9c563fecd0891689f96460fe0bc5604e47acb9fd514f22776f700cd85850ed9` |
| frozen (11 cells, = v1.0.2) | `stage2_v1_0_3_frozen.csv` | `81b38a2e94a784471aa11fa0169567ec20a237489f6df8cc1260c33af97a75d2` |
| drop list (1 subtree, = v1.0.2) | `stage2_v1_0_3_drops.csv` | `960e2c6560d4902b6ffdb980fcbb82cacf9589002c50cbc294b3880ca932923a` |
| public values (10 rows, TEAM-INTERNAL) | `stage2_v1_0_3_public_TEAM_INTERNAL.csv` | `a004a3331321fab55dc179ccc7d09965c04fa4e33dbd38f85bccc9e552b0b3e1` |
| tests | `test_stage2_v1_0_3_suppress.py` | `4ea6059c435313532f8d9b87959c05680068b87ede7a1814423b3ee42e5fae5e` |
| module v1.4 (regression tests only) | `reference\stage2_v1_0_3_suppress_v1_4.py` | `1acb7eb325e8f0a950f52e371c9a3679e1a1eaf7a9ce0ec2b88f33a4b7dd4670` |
| merge reproducer (tests) | `synthetic_merge_v105.py` | `f3ca3aa59b6461272af62e263edda5b3b80e202dcd06099a79b465f9f62c3a0f` |

**Links v3** = the 13 v1.0.2 rows unchanged, plus 76 Stage 2 internal equalities generated by `make_stage2_links_v1_0_3.py` from the v1.0.1 **count-field names** only (all pointers resolve in the v1.0.1 registry; checked on the box, names only):

| family | rows | equality |
|---|---|---|
| L5a | 21 | `layer_a/<split>/age_band/n_<level>` = `strata_denominators/age_band/<split>/<level>/n_encounters` (7 levels × 3 splits) |
| L5b | 3 | `layer_a/dob_breakdown/by_split/<split>/n_invalid` = `layer_a/<split>/age_band/n_age_invalid` |
| L5c | 15 | `layer_a/<split>/first_careunit/n_<unit>` = `strata_denominators/first_careunit/<split>/<unit>/n_encounters` |
| L6 | 34 | DOB breakdown overall = sum of the three splits (9 cells), by care unit = sum over splits (5 units × 5 cells) |
| L7 | 3 | `surgery/layer_a/<split>/n_encounters_recent` = `strata_denominators/surgical_status/<split>/recent_surgery/n_encounters` (both: linked surgery timed `recent_layer_a`) |

Each is read from the Stage 2 v1.0.1 script (inspected, not run). A link is a modelling claim: if the released values contradict one, the run stops (exit 2, LP infeasible) and nothing is written.

**Public file (team-internal; decisions D1, D2, D3, D6).** 3 × `n_ge18` = 0 (D1); `reconciliation/<table>/n_total`: CHARTEVENTS and MICROBIOLOGYEVENTS `publish` (published again, D3), LABEVENTS and PRESCRIPTIONS empty (value from the input; already published); `layer_a/dob_breakdown/n_invalid`, the plan-v13 stop-event value (D2: same definition as v1.0.1); `status/n_valid` and `n_age_days` overall, `derived` (follow from D2). **Handling (D6):** the file stays in the project folder only (confirmed not shared); it is never copied to a shared or synced team folder, attached to an email or chat, pasted into a console report, or quoted; its value must not be combined with v1.0.1 cells in any document. Outputs record only its SHA-256 and row count.

## 3. Rules (unchanged, plus v1.4; v1.5 changes selection only)

Rule 3 as in v1.0.2: a `<10` cell must keep the whole 0–9 range; an `S` cell at least 10 integers or a range straddling 9|10; a hidden or withheld difference must be provably ≥ 10, provably exactly 0, straddle 9|10, or span ≥ 10 integers. The residual check and the gate's independent re-solve are unchanged. v1.4 adds:

- public cells are known constants in every audit (input audit, suppression, gate), with their label unchanged in the release;
- carried cells (published in v1.0.1, suppressed in v1.0.2) start as `S`; the gate fails if one is published, except the two `publish` rows (D3), which the gate checks are published with their v1.0.1 values;
- "exactly 0" for a hidden or withheld cell is accepted only from an exact LP maximum (`screen_safe`).

## 4. Outputs (`data\interim\stage2_v1.0.3\`, a new folder)

| file | share? |
|---|---|
| `feature_sparsity_frequency_audit_v1.0.3.json` | yes, after Vinai's review; replaces v1.0.2 |
| `feature_sparsity_frequency_audit_v1.0.3.md` | yes, rendered from the JSON |
| `stage2_v1.0.3_release_summary.md` | yes; counts of cells by rule and the pointers of newly suppressed cells; no count values, no public-file content |
| `stage2_v1.0.3_manifest.json` | yes; SHA-256 and sizes of inputs (incl. the carried release and the public file) and outputs |

v1.0.1 and v1.0.2 files are never modified; `data\interim\stage2_v1.0.2\` is only read. The script refuses a non-empty or input folder as output, and refuses a released (v1.0.2 / v1.0.3) JSON as the `--in-json` input outside `--check-only`.

## 5. Tests and synthetic verification (box)

`test_stage2_v1_0_3_suppress.py`: **91 passed (about 6 min on the box)** (Python 3.13.5, numpy 2.5.5, pandas 3.0.6, scipy 1.18.1, pytest 9.2.0). The 54 v1.3 tests (regression against the delivered v1.3), 26 v1.4 tests and 11 v1.5 tests. v1.5:

- heartbeat lines during a silent phase (counts only, no pointer); option checks;
- a capped selection LP is never "safe" and never proves "cannot protect";
- every exact pick passes the full audit; pin-ranked candidates;
- the merge reproducer (`synthetic_merge_v105.py`, small, large-component path forced): v1.5 converges with a strictly decreasing failing count and the largest component at most twice its input size, while v1.4 merges it into a component more than twice v1.5's largest;
- convergence stops (`--max-new-suppressions`, `--max-component-cells`): exit 3, `STOPPING`, `WITHHOLD PROPOSAL`, the proposal file written outside the output folder (never over an existing file), no frozen cell in it, nothing in the output folder;
- `--withhold`: leaf cells leave the JSON and registry, the gate and `--check-only` pass; frozen, links-file and unknown pointers refused (exit 2); the proposal used as `--withhold` gives a release that passes.

v1.4 tests:

- the hidden-cell counterexample (v1.3 vs v1.4) and the `screen_safe` clauses;
- the links generator (21/3/15/34/3 rows on the mirror; shipped links v3 = 13 v1.0.2 rows + 76) and the shape of the real public file (row kinds only; no value printed);
- on the synthetic mirror (`synthetic_s2_v103_mirror.py`, which reproduces the real v1.0.1 label pattern of the flagged cells): the v1.0.2 mirror release fails the independent Stage 2 check with links v3 + public values (all six train families among the failures); the v1.0.3 mirror release passes it, and passes the Stage 3 view (links v2 internal rows + public table sizes); suppressions ⊇ v1.0.2; public labels unchanged; no public pointer or value in the summary, manifest or log; a tampered release (train invalid-age count published again) is caught;
- the decided release publishes the two table sizes (`publish` rows, D3) with their v1.0.1 values; a control with them kept `S` also passes alone; nine bad public files stop with exit 2 and write nothing; a `derived` row without the invalid total stops; the invalid total without the `derived` rows stops with `cannot protect` (the two cells are fixed exactly by public information);
- the carried release must be the v1.0.2 release of this input (hash, values, gate flag; not with `--check-only`; SHA pin);
- end to end with the Stage 3 builder (box only; skipped on the Copilot machine): `joint_check_synthetic.py`.

`stage2_independent_check.py` is a separate checker (no import of the module): it rebuilds a reader's system from the release JSON, the links, the public values and the input layout (withheld cells come back as unknowns) and computes the exact LP interval of every unknown cell.

**Synthetic results** (`demo/mirror_check_console.txt`, `demo/joint_check_console.txt`, `demo/merge_stress.txt`):

*Merge stress mirror (`synthetic_merge_v105.py --hub-frac 0.35`, all numbers invented): strata denominators linked to Layer A age / care-unit / surgical levels, forced-S hub strata with published overall rows — the structure that merged the real run's components. Same machine, `--max-seconds` set, `demo/merge_stress.txt`.*

| synthetic input | module | failing cells by round | new suppressions | largest component (input → end) | longest silence | wall clock | result |
|---|---|---|---|---|---|---|---|
| scale 0.5 (68.6k count cells, 189 links) | v1.4 | 535 → 816 → 21 → 2 → 0 | 4,120 | 4,124 → 48,954 | 271 s | 9 min 26 s | gate ok |
| scale 0.5 | **v1.5** | 535 → 133 → 50 → 23 → 0 | **667** | 4,124 → 4,124 | < 35 s | 2 min 10 s | gate ok |
| scale 1.0 (133.7k count cells, 189 links) | v1.4 | 1,003 → 1,695 → 22 → 3 → 0 | 8,856 | 8,461 → 78,528 | **1,004 s** (round 2 selection) | > 23 min | time limit (1,400 s) in the gate |
| scale 1.0 | **v1.5** | 1,003 → 139 → 68 → 39 → 15 → 5 → 0 | **960** | 8,461 → 15,796 | < 35 s | 3 min 33 s | gate ok |

v1.4 reproduces the real pattern (failing count up after round 1, a mass of new cells, a component growing several-fold, a selection phase silent for many minutes). v1.5's output passes `--check-only` with both v1.5 and v1.4 (0 failing). Small components (≤ 300 cells) take the unchanged look-ahead path: on the small merge mirrors (scale 0.05 / 0.1) v1.4 and v1.5 give the same round lines.

*Mirror (Stage 2 alone; `mirror_check_synthetic.py`, about 50 s):*

| release | independent check, links v3 + public | Stage 3 view (links v2 internal rows + table sizes) |
|---|---|---|
| v1.0.2 mirror (module v1.3) | exit 1: 163 cells fail, incl. 10 of the 15 flagged pointers | exit 1: 35 fail, incl. 11 of the 15 |
| **v1.0.3 mirror, decided D3: table sizes published** (module v1.4) | exit 0: 0 fail | exit 0: 0 fail |
| v1.0.3 mirror, control: table sizes kept `S` | exit 0: 0 fail | exit 0: 0 fail |

The remaining flagged pointers are not failures in the mirror by design: the table sizes are public values there (known, not audited), and the mirror's train `n_eligible_72h` values leave room. The decided mirror v1.0.3 run: gate OK in 16 s; 62 cells carried from the v1.0.2 mirror, of which 2 (the table sizes) are published again and 60 stay `S`; 10 public cells.

*Joint with the Stage 3 builder v1.0.3 (`joint_check_synthetic.py`, about 3 min):* Stage 3 fixture, 464 Stage 3 audit values captured; a Stage 2 mirror consistent with them under links v2; v1.0.2 by module v1.3, v1.0.3 by module v1.4 (the fixture's invalid-DOB total is `<10`, so it is not treated as public here).

| Stage 2 release | alone (both views) | Stage 3 gate | compromised Stage 2 cells | independent joint check |
|---|---|---|---|---|
| v1.0.2 mirror | fails (36 cells with links v3) | OK | 4 (MICROBIOLOGYEVENTS reconciliation, incl. `n_total`) | 0 narrowed |
| **v1.0.3 mirror, decided D3: table sizes published** | 0 fail | **OK** | **0** (**0 hidden**) | **0 / 0** |
| v1.0.3 mirror, control: table sizes kept `S` | 0 fail | OK | 1 (`/reconciliation/MICROBIOLOGYEVENTS/n_total`) | 0 narrowed |

*Scale (65 MB synthetic, `synthetic_scale_v101.py --scale 1.0`, no public / carry):* v1.3 105 s, v1.4 132 s (60,142 LP solves vs 52,278: the extra exact maxima of zero-screened hidden cells); byte-for-byte the same count cells (0 differing; 26,832 primary and 7,569 complementary in both). On a laptop expect roughly 5–15 min for the real run; the 2 h limit is a margin for the extra links and public values.

## 6. Copilot run instructions

Copy `stage2_v1_0_3_suppress.py`, `stage2_v1_0_3_links.csv`, `stage2_v1_0_3_frozen.csv`, `stage2_v1_0_3_drops.csv`, `stage2_v1_0_3_public_TEAM_INTERNAL.csv`, `test_stage2_v1_0_3_suppress.py`, `synthetic_v101_output.py`, `synthetic_scale_v101.py`, `synthetic_stall_v101.py`, `synthetic_s2_v103_mirror.py`, `synthetic_merge_v105.py`, `make_stage2_links_v1_0_3.py`, `stage2_independent_check.py`, `stage2_v1_0_2_links.csv` (read by one test) and the subfolder `reference\` (`stage2_v1_0_2_suppress_v1_0.py`, `stage2_v1_0_2_suppress_v1_3.py`, `stage2_v1_0_3_suppress_v1_4.py`; regression tests only) to the folder `scripts\stage2_v1_0_3_suppression\`. Keep the public file out of any shared or synced folder other than this project folder (D6). Run from the project root `<PROJECT_ROOT>` in PowerShell.

**Do not edit any file. Run each command exactly as written, one at a time. Stop at the first unexpected result and send back the full console text only. Do not patch, retry with other options, add a wrapper or logging, or open any `_raw_LOCAL` file.** Every Python command uses the same form: `$ErrorActionPreference='Continue'; cmd /c "python … 2>&1"; Write-Host "EXIT:$LASTEXITCODE"`.

**Step 1 – precheck.** All 11 hashes must match (case-insensitive), and neither the output folder nor the proposal file may exist (`False`, `False`). If anything differs, stop. (The stopped v1.4 run wrote nothing, so `data\interim\stage2_v1.0.3` should still be absent.)

```powershell
Get-FileHash scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_suppress.py -Algorithm SHA256
#   63d903d1ceef80a61f257a08aff4acd4a8c7fd5c2361d1e36de1c0cf603dc751
Get-FileHash scripts\stage2_v1_0_3_suppression\test_stage2_v1_0_3_suppress.py -Algorithm SHA256
#   4ea6059c435313532f8d9b87959c05680068b87ede7a1814423b3ee42e5fae5e
Get-FileHash scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_links.csv -Algorithm SHA256
#   a9c563fecd0891689f96460fe0bc5604e47acb9fd514f22776f700cd85850ed9
Get-FileHash scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_frozen.csv -Algorithm SHA256
#   81b38a2e94a784471aa11fa0169567ec20a237489f6df8cc1260c33af97a75d2
Get-FileHash scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_drops.csv -Algorithm SHA256
#   960e2c6560d4902b6ffdb980fcbb82cacf9589002c50cbc294b3880ca932923a
Get-FileHash scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_public_TEAM_INTERNAL.csv -Algorithm SHA256
#   a004a3331321fab55dc179ccc7d09965c04fa4e33dbd38f85bccc9e552b0b3e1
Get-FileHash scripts\stage2_v1_0_3_suppression\reference\stage2_v1_0_2_suppress_v1_3.py -Algorithm SHA256
#   df77234b6c68e1ca9b542f6f387f861141649fdd32e88f5f23ea009acac46d4f
Get-FileHash scripts\stage2_v1_0_3_suppression\reference\stage2_v1_0_3_suppress_v1_4.py -Algorithm SHA256
#   1acb7eb325e8f0a950f52e371c9a3679e1a1eaf7a9ce0ec2b88f33a4b7dd4670
Get-FileHash data\interim\feature_sparsity_frequency_audit.json -Algorithm SHA256
#   0a2123f48cb6bfbdcf2c0c98d8a886ea23cc9e4b05cf4bd605c31a7e7e2556d0
Get-FileHash data\interim\feature_sparsity_frequency_audit.md -Algorithm SHA256
#   9fe84830ece57ac4c03151db40907eb866d096f59f91e1ddd2b26aaedcc912ad
Get-FileHash data\interim\stage2_v1.0.2\feature_sparsity_frequency_audit_v1.0.2.json -Algorithm SHA256
#   a801e41ceb70a6bc2e56a5214a785bce31e35911536da3aadffb2281f162ac18
Test-Path data\interim\stage2_v1.0.3
#   False
Test-Path data\interim\stage2_v1.0.3_withhold_proposal_TEAM_INTERNAL.csv
#   False
```

**Step 2 – tests (synthetic only).** Expect `90 passed, 1 skipped` and `EXIT:0`, in about 8–15 min (about 7 min on the box). The suite prints only dots while it runs; stop it only if it passes 30 min.

```powershell
$ErrorActionPreference='Continue'; cmd /c "python -m pytest -q scripts\stage2_v1_0_3_suppression\test_stage2_v1_0_3_suppress.py 2>&1"; Write-Host "EXIT:$LASTEXITCODE"
```

**Step 3 – release run.** One line; copy it from the .md (in the PDF the ↪ marks only show line wrapping). Expect `result: RELEASE GATE OK` and `EXIT:0`.

```powershell
$ErrorActionPreference='Continue'; cmd /c "python scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_suppress.py --in-json data\interim\feature_sparsity_frequency_audit.json --in-md data\interim\feature_sparsity_frequency_audit.md --links scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_links.csv --frozen scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_frozen.csv --drops scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_drops.csv --public scripts\stage2_v1_0_3_suppression\stage2_v1_0_3_public_TEAM_INTERNAL.csv --carry-json data\interim\stage2_v1.0.2\feature_sparsity_frequency_audit_v1.0.2.json --expect-input-sha256 0a2123f48cb6bfbdcf2c0c98d8a886ea23cc9e4b05cf4bd605c31a7e7e2556d0 --expect-md-sha256 9fe84830ece57ac4c03151db40907eb866d096f59f91e1ddd2b26aaedcc912ad --expect-carry-sha256 a801e41ceb70a6bc2e56a5214a785bce31e35911536da3aadffb2281f162ac18 --max-seconds 5400 --withhold-proposal data\interim\stage2_v1.0.3_withhold_proposal_TEAM_INTERNAL.csv --out-dir data\interim\stage2_v1.0.3 2>&1"; Write-Host "EXIT:$LASTEXITCODE"
```

**Expected output lines** (counts only; `…` = a number not known before the run):

```text
[     0.1s] stage2_v1_0_3_suppress v1.5: reading the input JSON (56,866,492 bytes) and checking its schema
[     1.0s] input loaded: … count cells, … hidden differences, … identities
[     1.5s] drop list: 1 subtree(s) and 0 single cell(s) withheld from the release, 2 count cells marked unpublished (audited as hidden)
[     2.5s] carried release: 148 cell(s) published in the input stay S (suppressed there)
[     2.6s] public cells: 10 value(s) treated as known in every audit (labels unchanged; 2 carried cell(s) published again: 'publish' rows)
[     2.7s] equation system built: … equations (89 links), 11 frozen cells
[     …   ] input audit: … components (…); largest multi-equation component … cells
[     …   ] input audit done: … suppressed cells audited, … failing, … groups failing the residual check
[     …   ] complementary suppression: start
[     …   ] suppression round 1: audit: … components (…); largest multi-equation component … cells
[     …   ] suppression round 1: … failing cells, first round
[     …   ] heartbeat: still working - suppression round 1: selection (failing cell …/…, … picks, … exact checks); LP solves so far …; time left … s
[     …   ] suppression round 1: … new suppressions (total new …; look-ahead simulations …; exact candidate checks …; choices …; selection … s)
   … the same for each further round; the failing count should go down round by round …
[     …   ] complementary suppression: done (… new cells, … rounds)
[     …   ] release gate: reloaded pending files; checking
[     …   ] release gate audit (independent re-solve): …
[     …   ] release gate: ok
[     …   ] release written
# Stage 2 v1.0.3 release summary (release)
- script: stage2_v1_0_3_suppress v1.5 (sha256 63d903d1ceef80a61f257a08aff4acd4a8c7fd5c2361d1e36de1c0cf603dc751)
- spec: stage2_v1_0_3_suppression_spec v1.5
- result: RELEASE GATE OK
- links used: 89 (file sha256 a9c563fecd0891689f96460fe0bc5604e47acb9fd514f22776f700cd85850ed9)
- frozen cells: 11 (file sha256 81b38a2e94a784471aa11fa0169567ec20a237489f6df8cc1260c33af97a75d2)
- public (already known) cells: 10 (file sha256 a004a3331321fab55dc179ccc7d09965c04fa4e33dbd38f85bccc9e552b0b3e1; team-internal file, pointers and values not listed)
- carried release: feature_sparsity_frequency_audit_v1.0.2.json (sha256 a801e41ceb70a6bc2e56a5214a785bce31e35911536da3aadffb2281f162ac18; 146 cell(s) carried as S; 2 carried cell(s) published again as public by definition ('publish' rows))
- dropped subtrees: 1 holding 2 count cells, withheld from the release as a whole (file sha256 960e2c6560d4902b6ffdb980fcbb82cacf9589002c50cbc294b3880ca932923a; pointers below, no values)
| 2 complementary, carried from the previous release (S there, published in v1.0.1) | 146 |
EXIT:0
```

The input audit is expected to fail some cells (the new links and public values expose cells that v1.0.2 left published); the suppression rounds then protect them. **Runtime.** The stopped v1.4 run reached round 5 in about 6.6 min; on the box the larger synthetic merge mirror (134k count cells, §5) takes 3.5 min with v1.5. Expect **about 10–30 min** on the laptop, rarely up to an hour; the module stops itself at `--max-seconds 5400` (90 min, exit 3, nothing written).

**Timeout guidance for Copilot.** v1.5 prints a line at least every 30 s (a progress line or a `heartbeat:` line). Watch only for silence: if **no new line appears for 3 minutes**, the process is hung — stop it (Ctrl+C), and send back the full console text. Do not stop a run that is still printing heartbeats, and do not wait longer than **100 minutes** wall clock in total (the 90 min module limit plus the gate). Nothing is written to `data\interim\stage2_v1.0.3\` unless the gate passes, so a stopped run leaves nothing behind; do not delete or edit anything, and do not re-run with other options.

Possible stops (send the text back, change nothing): exit 2 `STOPPED (input): …` (a hash, a public row, the carried release, or a link contradicted by the released values — LP infeasible); exit 3 `RELEASE GATE FAILED: cannot protect N cell(s) … (first pointers: …)` (decision D4: send the pointers back; by default the subtrees holding them are withheld in a re-run with new drop-list rows made on the box, never an edit on the Copilot side); exit 3 `time limit of 5400 s`; exit 3 `STOPPING: suppression not converging: …` followed by `WITHHOLD PROPOSAL: N cell(s) …` and `RELEASE GATE FAILED: suppression not converging …` (v1.5 convergence stop): send back the console text and the proposal file `data\interim\stage2_v1.0.3_withhold_proposal_TEAM_INTERNAL.csv` (pointers only, no values; team-internal). Vinai reviews it under D8; the box then issues either a `--withhold` file or new drop-list rows with new hashes, and a fresh run into a new, empty folder.

**Step 4 – paste back** the console text of Step 3 (counts and pointers only) and the contents of `data\interim\stage2_v1.0.3\stage2_v1.0.3_manifest.json`. Do not paste the v1.0.3 JSON or markdown, any count value, or the public file. Do not share the v1.0.3 files until Vinai has reviewed the summary and Stage 3 v1.0.4 has passed jointly with it.

## 7. Decisions (Vinai, Nov 17, 2026)

- **D1 ≥ 18 structural zero — public.** `n_ge18` = 0 in all three splits is a known value (derivable from the release; PIC is pediatric).
- **D2 invalid-DOB count — public.** The v1.0 stop event used the same invalid-DOB definition as v1.0.1 (missing, unparseable or negative DOB), so the stop-event row stays, with its two `derived` cells (overall valid DOBs, `n_age_days`).
- **D3 table sizes — published.** CHARTEVENTS and MICROBIOLOGYEVENTS `reconciliation/<table>/n_total` are published again in v1.0.3 with their v1.0.1 values (`publish` rows in the public file; the default release files and command). Their values are public (raw table sizes); keeping them `S` protected nothing and left them on the Stage 3 compromised list.
- **D4 cannot protect — withhold by default.** If the real run stops with `RELEASE GATE FAILED: cannot protect N cell(s) … (first pointers: …)`, Copilot sends the console text back unchanged. Vinai decides from the pointers; **the default is to withhold the smallest JSON subtree holding each named cell** (whole subtree, totals above it kept), as for the v1.0.2 NMB drop: new rows in `stage2_v1_0_3_drops.csv` with the reason and date, new hashes from the box, the spec and Copilot steps re-issued, and a fresh run into a new, empty output folder. A subtree that holds a frozen or links-file cell cannot be dropped (the run stops with exit 2); that case goes back to Vinai. Nothing is ever edited on the Copilot side.
- **D5 uncertain links — left out.** Weight `n_charted` (Layer A) vs the DOB-breakdown weight imputation, and strata eligibility vs landmark denominators, are not links (different variables / definitions unconfirmed).
- **D6 public file — project folder only.** The project folder is confirmed not shared. The file is never placed in any shared or synced team folder, email, chat or report (see §2).
- **D7 v1.0.2 holders — asked to delete.** Once v1.0.3 is released, everyone who holds v1.0.2 is asked to delete it; the records note is §8.
- **D8 convergence stop — caps kept, withhold fallback (decided by Vinai, Nov 17, 2026).** (a) Keep the caps `--max-new-suppressions 2500` and `--max-component-cells 20000` (the real input's largest component is about 7.7k cells; on the synthetic mirror v1.5 at most doubled it, v1.4 grew it 5.8×); (b) on a convergence stop, accept the cell-level withhold (`--withhold`) for the rows marked "suppressed cell fixed by its equations" (an already-suppressed cell; withholding it removes only its label), and review the rows "frees a '<10' / hidden cell" against D4's subtree drop, taking whichever withholds fewer published cells; (c) a failing cell listed as `blocked` (it needs a frozen, public or links-file cell) goes to D4 (subtree) — a link or a public value is never removed to make room.

## 8. Records note (for the log)

- **What was exposed.** The Stage 2 v1.0.2 release (Nov 16–17, 2026) protected each table on its own and across the 13 declared links. A reader who also makes the undeclared Stage 2 internal equalities (age bands and DOB breakdown = strata denominators, L5) and uses public values (the raw table sizes; the ≥ 18 structural zero; the invalid-DOB count quoted in the plan) can narrow 15 suppressed Stage 2 cells and 25 hidden differences below the protection rule: train invalid-age, adolescent and ≥ 18 cells, and the small CHARTEVENTS / MICROBIOLOGYEVENTS reconciliation classes. Found by the Stage 3 builder v1.0.3 joint audit (its public-information baseline), Nov 17, 2026. No single cell is quoted here, and no count is.
- **Fix.** Stage 2 v1.0.3 (module v1.4): the equalities are declared links, the public values are known in the audit, every v1.0.2 suppression is carried, the two table sizes are published again (D3), and the hidden-cell exact-bound bug of module v1.3 is fixed. The gate re-solves everything; on synthetic data the v1.0.3 release passes alone and jointly with Stage 3.
- **Action (D7).** When v1.0.3 is released (gate OK, reviewed by Vinai, Stage 3 v1.0.4 passing jointly), the shared copy is replaced by v1.0.3 and everyone who received v1.0.2 (teammates, shared folders, Copilot's workspace, assistant copies) is asked to delete it. The log records the date v1.0.2 was pulled and replaced, who had it, and each confirmation of deletion. v1.0.1 was already withdrawn (v1.0.2 records note). Whether PhysioNet needs to be told is Vinai's decision under the PIC DUA text.
- **Rule going forward.** A release gate audits the release with every equality a reader can make from the variable definitions (not only the registry), with public values as known, and jointly with every other released stage.

## 9. Stage 3 v1.0.4 (plan; not built)

- **Re-pin Stage 2:** `STAGE2_RELEASE_SHA256` = the real v1.0.3 JSON hash from the Step 4 manifest; label "Stage 2 v1.0.3"; the release path `data\interim\stage2_v1.0.3\feature_sparsity_frequency_audit_v1.0.3.json`.
- **Re-pin the suppression module** to v1.5 (`stage2_v1_0_3_suppress.py`, `63d903d1…`): the hidden-cell bound fix applies to Stage 3's own audit too. Adapt the provenance / note renaming (module keys are now `release_v1_0_3`, `dropped_subtrees_v1_0_3`) and the output names to v1.0.4.
- **Links v3 (`stage3_links_v3.csv`):** links v2 plus the 76 Stage 2 internal rows (as `s2:` rows), so the joint system holds the same Stage 2 equations as Stage 2's own audit; keep the six L5 rows.
- **Public constants in the baseline:** a team-internal constants input (same format as the Stage 2 public file, `s2:` pointers): `n_ge18` = 0, the D2 value and its derived cells; table sizes are already public through `/raw_rows`, and Stage 2 v1.0.3 publishes them (D3), so no exemption is needed.
- **Drops:** keep `stage3_drops_v1.csv` (20 subtrees, incl. the withheld event-linkage cells) unchanged in v1.0.4; review each after the joint run shows what is publishable.
- **Token mappings (decided):** `positive+` and `positive++` → positive for `stool_occult_blood`, `urine_leukocyte_esterase`, `urine_nitrite`, `csf_protein_pandy` (channel-specific, not global); `resp_panel_mycoplasma` `<3.00` → negative with a new below-cutoff flag column (Layer B schema and audit cells change). Not trivial: needs a per-channel token table, the flag column, tests, and the unparsed-string list re-run.
- **Expected:** release gate OK, `flagged_for_stage2_v1_0_4` empty, `0 hidden Stage 2 difference(s)`, independent joint check 0 / 0.
