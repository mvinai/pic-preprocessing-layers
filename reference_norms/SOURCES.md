# Reference-norm downloads (D2 prep) — provenance

Retrieved Wed Nov 18, 2026, 06:29–06:31 CST (America/Chicago) on the box with `curl`. Raw files are unmodified in `raw/`.

| raw file | source URL | SHA-256 |
|---|---|---|
| `osf_upqtv_providers.json` | https://api.osf.io/v2/nodes/upqtv/files/ (OSF API listing; the folder tree was walked via `.../files/osfstorage/`) | `f55c0bce9d08f8077e1f8c15d2f59ba8d65a971c0434a07b5ca777eab76aeba6` |
| `roberts_Roberts_et_al_2020_MAP_percentiles.csv` | https://osf.io/download/cgm6q/ | `38219710fc64000f61d8c5c88236d406fa6f95a1ced16a34349a9a6a924d2316` |
| `roberts_Roberts_et_al_2020_MAP_thresholds_birth_sex.csv` | https://osf.io/download/8ms96/ | `5b27f21494b9ad16f294af883ba805786270bfe057456f7adab799a6d9fda63a` |
| `roberts_Roberts_et_al_2020_hypotension_thresholds.csv` | https://osf.io/download/r62gw/ | `9ef195256dbe448dbbf690a65c7be2ab66bb44578890650d3fbec8555d2024a1` |
| `roberts_Roberts_et_al_2020_MAP_thresholds_code.R` | https://osf.io/download/3btvf/ | `113abee04fc5f8c0cbc0ecefb63d0e237940ca93d517d345f18960893c64939d` |
| `roberts_Roberts_Data_Dictionary_Inpatient_MAP.pdf` | https://osf.io/download/mb62v/ | `24b9a53088f364b2a6204d5979ffc53b6ff592d2cc836ac9c9309af6c16c854d` |
| `roberts_MAP_valid_data_thresholds_graphs.pdf` | https://osf.io/download/yk53b/ | `db9e25995bff5e34273c297c07aeb59267efcf276c6295bb912c9a4bbce5ee8a` |
| `roberts_license.txt` | https://osf.io/download/8ydxq/ (data: ODbL 1.0; contents: DbCL 1.0) | `b0086e1efbe6fc223e3c1bafd82ee4e5daf2970b06e54c2982bf305962e15d32` |
| `fleming_susz_Fleming2011_Lancet_Webappendix.pdf` | https://www.susz.me.uk/pubs/Fleming2011_Lancet_Webappendix.pdf (first author's site) | `01c896d1faf1936f9d313de5813c5542ee7dc729fccc28463a43cc10fd2e9cc0` |
| `fleming_ora_webappendix.pdf` | https://ora.ox.ac.uk/objects/uuid:d7afbf70-40f7-4027-bcbb-b729f59298ac/files/m96a2b837ffdff6168f77d3a7e67caf4d (Oxford Research Archive) — byte-identical to the file above | `01c896d1faf1936f9d313de5813c5542ee7dc729fccc28463a43cc10fd2e9cc0` |
| `fleming_PMC3789232_efetch.xml` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=3789232&rettype=xml (PMC author manuscript text; contains no tables) | `7c092cfaffecae145b3ac4e305e2e718a1dbd1b71b2d902093c2f6db72c09775` |

Not downloaded: `Roberts_et_al_2020_MAP_thresholds_data.zip` (7.2 MB, https://osf.io/download/hyvpb/; the bootstrapped single-reading data — not needed for the lookup table). The PMC copy of the Fleming Web Appendix (`NIHMS54714-supplement-Web_Appendix.pdf`) returned 403 / a browser-check page, so the identical author and ORA copies were used. Europe PMC full-text XML returned HTTP 500.

Derived files (produced by `extract_tables.py` from the raw files above; no value typed by hand):

| file | SHA-256 |
|---|---|
| `fleming2011_rr_centiles.csv` (Web Table 4) | `ec54810ee50120ad3ea5c69c90a5f7cbd6c3ef09de560d102786f8de7c94e552` |
| `fleming2011_hr_centiles.csv` (Web Table 5) | `4ff3c424dd398ded5b7e9a7a70de45f3b44fff8f9fe456a2a3be0e61d331f1ad` |
| `roberts2020_map_percentiles.csv` (OSF Table 2 lookup + age bounds parsed from the R code) | `421c965f262335e2445674a4fe78c9fa756b10d28b64a206eb3c38b57f1ef2f2` |
| `derived_spread_sd_DERIVED.csv` (normal-theory SDs implied by centile spreads; illustrative) | `1f74a877801e4ecde337bb72a178b9c04d3ec52d51696c05f39edf1013dae4e8` |

## PODIUM 2022 (added 2026-11-18 06:42 CST)
Fitzgerald JC, Basu RK, Fuhrman DY, et al. Renal Dysfunction Criteria in Critically Ill Children: The PODIUM Consensus Conference. Pediatrics. 2022;149(1 Suppl 1):S66-S73. doi:10.1542/peds.2021-052888J. Table 2 (baseline SCr by age and sex), transcribed from raw/podium_renal_2022.pdf (SHA-256 aa9fdb127f2dc6ad3729d1f45902bbe175b6b04a7c1573e76604d5ab1bf26753, link supplied by Vinai) into podium2022_baseline_scr.csv (SHA-256 e56de3e8f9340ada0d163e520c58dbf657a21fddfc179ca0213f184ca79190e7).
