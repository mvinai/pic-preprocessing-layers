# growth_weight_for_age_ref_v1.csv: provenance (Stage 3 weight normalizer, spec §10)

Frozen November 17, 2026. SHA-256 of the CSV:
`8e53ca2d4fd9f053c2fa56b153e90b4ea63d090804add22414678fd80b88b54e` (3,926 rows + header; pinned in
`03_build_real_layers.py` as `PINS["growth_ref"]`). Columns: `reference, component, sex, age_days, age_months, L, M, S`.
The builder uses `M` (the LMS median = P50) only; L and S are kept for later z-scores. Ages: `age_days` = months × 30.4375
(WHO convention, same as the builder's `age_months`); WHO 2006 rows are by day as published.

## reference = cn_li2009 (primary): China 0–18 y growth standard, weight-for-age

- Standard: Li H, Ji CY, Zong XN, Zhang YQ. Height and weight standardized growth charts for Chinese children and
  adolescents aged 0 to 18 years. *Zhonghua Er Ke Za Zhi* 2009;47(7):487–492 (PMID 19951507). Recommended as the China
  national growth standard; data: 2005 National Growth Survey of Children under 7 in nine cities and the 2005 Chinese
  National Survey on Students' Constitution and Health (94,302 healthy urban children; Beijing, Shanghai, Harbin, Xi'an,
  Nanjing, Wuhan, Guangzhou, Fuzhou, Kunming); LMS smoothing; exclusions included preterm birth and birth weight < 2,500 g.
- Numbers: the same authors' open-access English publication of the LMS parameters, Zong XN, Li H. Construction of a
  new growth references for China based on urban Chinese children: comparison with the WHO growth standards. *PLoS ONE*
  2013;8(3):e59569, doi:10.1371/journal.pone.0059569, **Table 2** "Weight-for-age SD curves (in kg) for Chinese boys and
  girls, 0–18 years" (L, M, S, −2 SD … +2 SD; exact ages 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 21, 24, 30 … 216 months;
  46 ages × 2 sexes). Licence CC BY 4.0.
- Extracted from the article's JATS XML
  (`https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0059569&type=manuscript`, SHA-256
  `2a6cff2ea7df34652e695e72e9873752ed958ebf55b87a15f110ab27c920e5c0`); every M equals the table's Median column; all
  92 M/Median cells re-checked against the printed PDF
  (`…&type=printable`, SHA-256 `e20dbf8dda1edfefa92b092b2670d19e5ab3b5ba936291faf0e2277be29da6f2`): 0 mismatches.
  Source typo: girls 10 months L printed "−−0.21" (XML and PDF); stored as −0.21 (does not affect M).
- Cross-check (not used by the builder): WS/T 423-2022 *Growth standard for children under 7 years of age* (National
  Health Commission, 2022-09-19; `https://www.nhc.gov.cn/cms-search/downFiles/e38068f0a62d4a1eb1bd451414444ec1.pdf`,
  SHA-256 `1fa2f1f0b299028063fc82b58920a2ea058a04aa2bd7e5bf3e17e64e77c067e5`), Tables A.1/A.2 P50 (one decimal; 2015
  survey). Over 40 parsed ages 0–81 months, China 2009 M − WS/T 423-2022 P50 = mean −0.08 kg (range −0.35 to +0.05) in
  boys and −0.08 kg (−0.39 to +0.05) in girls; the largest gap is at birth (boys 3.32 vs 3.5 kg; girls 3.21 vs 3.3 kg).

## reference = who (fallback where China does not cover the age; also the WHO sensitivity arm)

- WHO Child Growth Standards (2006), weight-for-age, expanded tables by day 0–1856 (z-scores):
  `https://cdn.who.int/media/docs/default-source/child-growth/child-growth-standards/indicators/weight-for-age/expanded-tables/wfa-boys-zscore-expanded-tables.xlsx`
  (SHA-256 `b5b4748c6bfa5230e2eddafa1767629c349178b08d457f400b59422b8bfef86c`) and `…/wfa-girls-zscore-expanded-tables.xlsx`
  (`ee3ae12cb96c6c5541cdf43665c03ce6c984f877859a183a5f6104eb06a49a6e`).
- WHO Reference 2007 (5–19 y), weight-for-age 5–10 y, expanded tables by month 61–120 (z-scores), linked from
  `https://www.who.int/tools/growth-reference-data-for-5to19-years/indicators/weight-for-age-5to10-years`:
  `…/weight-for-age-(5-10-years)/hfa-boys-z-who-2007-exp_0ff9c43c-8cc0-4c23-9fc6-81290675e08b.xlsx`
  (`a6ed0d9f3cfa209747afe49f7503b7cc24a1be6b6a8e7222e85a36b7850cb99f`) and
  `…/hfa-girls-z-who-2007-exp_7ea58763-36a2-436d-bef0-7fcfbadd2820.xlsx`
  (`d747e068fcef4238cbaf1c201dd7788f266d9903685f4495360614eed0153562`) (WHO's file names say "hfa"; the sheets are
  weight-for-age, M 18.5057 kg at 61 months for boys as in the WHO simplified field table). WHO publishes no
  weight-for-age reference beyond 10 years.
- CDC 2000 (`wtage.csv`, for > 10 y) was **not** used: cdc.gov returned "Access Denied" (HTTP 403) from the box and the
  file is not on ftp.cdc.gov; it is also a US reference. With China covering 0–216 months, the primary chain never
  needed it.

## Build

`build_growth_ref.py` (see `scripts/growth_ref/build_growth_ref.py` and the provenance note; source hashes in
`raw/SOURCE_SHA256SUMS.txt`) writes the CSV deterministically from the downloaded files (no hand-typed numbers).
