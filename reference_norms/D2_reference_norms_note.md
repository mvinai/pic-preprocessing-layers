# D2 prep: published age norms for HR / RR / MAP (no script change yet)

Nov 18, 2026. Files originally under `stage06_08/reference_norms/` (provenance and hashes: `SOURCES.md`). Every number below comes from a downloaded file or is computed from one (`derived_spread_sd_DERIVED.csv`); none was typed from memory.

## 1. What was obtained

| source | signal | age bands | statistics | status |
|---|---|---|---|---|
| Roberts, Yanay, Barry, PCCM 2020 (OSF upqtv, `Roberts_et_al_2020_MAP_percentiles.csv` = the paper's Table 2 lookup) | MAP (cuff) | 21 bands: 37 wk PMA–30 d, 1–3 mo, 3–6 mo, 6–12 mo, then yearly 1–2 y … 17–18 y (upper bounds in days from the R code: 30, 90, 180, 365.25, then multiples of 365.25 to 6575) | 13 percentiles: 1, 5, 10, 20, 30, …, 90, 95, 99 (integers) | **obtained** |
| same OSF project | MAP by birth sex | same bands × female / male | 9 percentiles (1, 5, 10, 25, 50, 75, 90, 95, 99) | obtained, not recommended (below) |
| same OSF project | SBP / MAP / DBP | same bands | measured 5th percentile only ("hypotension thresholds") | obtained |
| Fleming et al., Lancet 2011, Web Appendix Web Table 4 | RR | 13 bands: 0–3 m, 3–6 m, 6–9 m, 9–12 m, 12–18 m, 18–24 m, 2–3 y, 3–4 y, 4–6 y, 6–8 y, 8–12 y, 12–15 y, 15–18 y | centiles 1, 10, 25, 50, 75, 90, 99 (integers) | **obtained** |
| Fleming et al., Web Table 5 | HR | 14 bands: "Birth" (immediate neonatal period) + the same 13 | same 7 centiles | **obtained** |
| Haque & Zaritsky, PCCM 2007 | MAP | ages 1–17 y | 50th = 1.5 × age(y) + 55; 5th = 1.5 × age(y) + 40 | formula (from the brief); also coded identically in the downloaded Roberts R code (`HZ_50`, `HZ_5`; its `HZ_95` = +70 is Roberts' own extrapolation, not Haque's) |

How the tables were made (from the downloaded documents): Roberts — quantile regression (natural spline, 5 df) on single cuff MAP readings sampled from about 85,000 hospitalized children (Cerner EMR, Feb 2012–Jun 2016), each row being the fitted curve evaluated at one age (15, 60, 135, 275, 548 … 6392 days), rounded to integers. Fleming — kernel regression centile curves from 69 studies of healthy children; each band value is the mean of the centile curve over the band, rounded to an integer. The Fleming paper states it found no skewness and **assumed a normal distribution at each age**.

Discrepancy to note: the Roberts sex-specific file has higher medians than the main table for the same bands (e.g. neonatal p50 62 / 61 vs 56). The R code reads it from a different input (`MAP_birth_sex_bootstrapped.csv`); the reason is not documented in what was downloaded. Use the main table.

## 2. Turning each table into a fixed z-score table

**Option A — normal (mean, SD) per band.** mean = median (p50); SD from a centile spread, e.g. SD = (p90 − p10) / 2.5631 (or (p75 − p25) / 1.3490, (p99 − p1) / 4.6527). For Fleming this is consistent with the authors' normality assumption, and the three spreads agree to within about 0.5 units in most bands (`derived_spread_sd_DERIVED.csv`). Caveat: integer rounding is coarse where the spread is small (older-child RR: implied SD about 2.2–2.4 breaths/min, so ±0.5 rounding is about 0.2 SD) — averaging the three spread-based SDs reduces this.

**Option B — centile-based z (no normality).** Place the observed value on the band's centile function by piecewise-linear interpolation between the published centiles, convert the percentile p to z = Φ⁻¹(p), and cap at the outermost centiles (|z| ≤ 2.33 for p1/p99) or extrapolate the outer segment linearly. Better for MAP: the Roberts percentiles are right-skewed in infancy (p90 − p50 exceeds p50 − p10 by up to 4 mmHg), though nearly symmetric after age 3.

**Option C — smooth over age.** Both tables are step functions by band. To avoid jumps at band edges, interpolate each centile linearly in age between band reference ages (Roberts gives the exact evaluation age of each row; for Fleming use band midpoints). Lookup by each encounter's age in days/months (Stage 3 has `age_months`), not by the coarser PIC age bands.

**Suggested combination (for Vinai's D2 decision):** HR and RR — Fleming, Option A with SD from (p90 − p10) (or the mean of the three spreads), optionally with C; MAP — Roberts main table, Option B with C; Haque–Zaritsky only as a sensitivity check for 1–17 y (it implies a constant SD of (55 − 40) / 1.645 ≈ 9.1 mmHg and has no infant range). Labs stay on the empirical train age-band z (no published lab table was obtained), so D2 would become "published norms for HR/RR/MAP, empirical for labs".

## 3. Gaps and caveats

* **MAP under 1 year:** not covered by Haque–Zaritsky (1–17 y); **covered by Roberts** (37 wk PMA–30 d, 1–3, 3–6, 6–12 mo). Preterm infants below 37 wk PMA are not covered, and PIC has no gestational age, so every neonate ≤ 30 days would use the first Roberts row.
* **RR at birth:** Fleming RR starts at 0–3 m (no "Birth" row); HR has a "Birth" row for the immediate neonatal period that overlaps 0–3 m — a rule is needed (suggest: use 0–3 m for all ages 0–3 months and ignore "Birth", or use "Birth" only for the first days of life).
* **Age ≥ 18 y (PIC `ge18`) and `age_invalid`:** no table covers ≥ 18 y (Fleming ends at 18 y, Roberts at 6575 days). Options: carry the last band forward, or fall back to the empirical z; `age_invalid` stays empty with its mask.
* **Reference populations differ from a PICU cohort:** Fleming = healthy children (mostly developed countries); Roberts = hospitalized children (cuff MAP, one US EMR). PIC children are Chinese PICU patients, and PIC MAP may be invasive or cuff. The z-scores therefore measure deviation from a normal / ward reference, not from the PICU train distribution; stage 08 re-standardizes them on train anyway, so they act as an age adjustment.
* **Rounding:** all published values are integers.
* **Licensing / citation:** Roberts data under ODbL 1.0 (attribution and share-alike for the database); cite Roberts et al. 2020 and the OSF project, Fleming et al. 2011 and Haque & Zaritsky 2007.
* **Not obtained:** published age-specific lab reference cut-offs; the Roberts raw data zip (not needed).
