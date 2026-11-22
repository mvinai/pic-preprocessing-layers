# PIC preprocessing and three-layer feature build

This repository holds the **preprocessing pipeline** that turns the credentialed [Paediatric Intensive Care (PIC) database](https://physionet.org/content/picdb/) into the three-layer A / B / C feature representation used by the Foundations analysis.

A companion repository, [pic-foundations-analysis](https://github.com/mvinai/pic-foundations-analysis), holds stages 09–18 (spectral diagnostics, models, attribution, and the desktop sensitivity runs).

**Authors.** Vinai (principal investigator), with contributions from Habib, Vinay and Yashwant. Scripts were written with AI assistance and validated on synthetic data before any real PIC run.

This repository is released under the MIT License. Copyright (c) 2026 Vinai Modem. See `LICENSE`.

## What this is not

This is not a copy of PIC. The database is credentialed health data under a PhysioNet data-use agreement (DUA). Nothing in this repository is PIC row-level data, a real-run output, or a cohort or split count. To work with the real database you must apply for access yourself and keep the data off this repository.

## Accessing PIC

1. Create a [PhysioNet](https://physionet.org/) account and complete credentialing.
2. Request access to the [Paediatric Intensive Care database](https://physionet.org/content/picdb/) and accept its DUA.
3. Download the tables locally (for example `PATIENTS`, `ADMISSIONS`, `ICUSTAYS`, `CHARTEVENTS`, `LABEVENTS`, `PRESCRIPTIONS`, `MICROBIOLOGYEVENTS`, `SURGERY_INFO`, and the dictionaries). Do not commit them.
4. Restore the omitted reproduction pins from your local project into `configs/cohort_pins.json` (see `configs/cohort_pins.example.json` and [`docs/PUBLIC_REDACTIONS.md`](docs/PUBLIC_REDACTIONS.md)). SHA-256 hashes of the pinned inputs are already in the scripts.

## How the pipeline fits together

| Stage | Script | Role |
|---|---|---|
| Stage 1 audit | `scripts/audit/pic_dataset_audit_v12.py` and the audits under `scripts/band_support_audit`, `scripts/duration_audit_v2`, `scripts/feature_item_inventory`, `scripts/feature_sparsity_frequency_audit`, `scripts/stage2_followup_checks`, plus `scripts/check_encounter_harmonization_consistency.py` | Aggregate-only audits of coverage, maps, surgery timing, and encounter harmonization. |
| Stage 2 | `scripts/stage2_v1_0_3_suppression/stage2_v1_0_3_suppress.py` (v1.0.3, module v1.5) | Cell suppression and disclosure control for the Stage 2 frequency audit. |
| Stage 3 | `scripts/stage3_builder/03_build_real_layers.py` (v1.0.4) | Builds the Layer A and B bases and the Layer C raw grid from PIC and the frozen maps. |
| Stages 06–08 | `scripts/stage06_08/06_08_build_levels_preprocess.py` (v1.2.1) | Layer A / B levels and train-only preprocessing of Layer C, including reference vital-sign and growth norms. |

Maps and reference tables used by Stage 3 and Stages 06–08 also live under `scripts/stage3_builder/maps/` and `scripts/stage06_08/`, with copies of the released complexity and growth tables under `maps/` and `reference_norms/`.

## Environment

The released runs used:

- Python 3.14
- numpy 2.4.2
- pandas 3.0.1
- scipy 1.17.x
- scikit-learn 1.8.0
- matplotlib 3.10.8
- torch 2.10 (not required for the preprocessing tests)

See `requirements.txt`. A nearby CPU-only stack is enough to run the synthetic tests.

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Running the tests (synthetic data only)

The tests build their own temporary synthetic data. They do not read PIC.

From the repository root, with `PYTHONNOUSERSITE=1` and `python -s` if you want the same isolated invocation used for the released runs:

```bash
export PYTHONNOUSERSITE=1

# Stage 1 audits
python -s -m pytest -q -p no:cacheprovider scripts/audit/pic_dataset_audit_v12.py
python -s -m pytest -q -p no:cacheprovider scripts/band_support_audit
python -s -m pytest -q -p no:cacheprovider scripts/duration_audit_v2
python -s -m pytest -q -p no:cacheprovider scripts/feature_item_inventory
python -s -m pytest -q -p no:cacheprovider scripts/feature_sparsity_frequency_audit
python -s -m pytest -q -p no:cacheprovider scripts/stage2_followup_checks
python -s -m pytest -q -p no:cacheprovider scripts/test_check_encounter_harmonization_consistency.py

# Stage 2 suppression (v1.0.3 / module v1.5)
python -s -m pytest -q -p no:cacheprovider scripts/stage2_v1_0_3_suppression/test_stage2_v1_0_3_suppress.py

# Stage 3 builder (v1.0.4)
python -s -m pytest -q -p no:cacheprovider scripts/stage3_builder/test_03_build_real_layers.py

# Stages 06–08 (v1.2.1)
python -s -m pytest -q -p no:cacheprovider scripts/stage06_08/test_06_08_build_levels_preprocess.py
```

Or, in one pass: `python -s -m pytest -q -p no:cacheprovider scripts`.

## Synthetic datasets

`synthetic/` contains generators and, where they are small enough, the generated tables:

- `synthetic/v1/` — a small tabular stand-in generator (`make_synthetic.py`, default 2,000 encounters). The previously generated full-size tables are not shipped.
- `synthetic/planted_v1/` — planted-spectrum checks.
- `synthetic/svd/` — generators and READMEs for a larger layered synthetic tensor. The generated arrays are not shipped (they are well over the size limit). Recreate them by running the two `make_synthetic_svd_*.py` scripts as described in those READMEs.
- `synthetic/SYNTHETIC_DATASET_GUIDE.md` — how the stand-ins relate to Layers A / B / C.

## Documentation

Final specifications (markdown only) are under `docs/`. `docs/PUBLIC_REDACTIONS.md` lists what was removed for this public copy and why.

## Acknowledgements

The synthetic data generators in `synthetic/` were written with Grok Bot (Research Associate), an AI research assistant, under Vinai's direction and review.

## Citation

See `CITATION.cff`.
