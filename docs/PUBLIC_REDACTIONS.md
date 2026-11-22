# What was omitted from this public copy

This repository is a public, DUA-safe copy of the analysis code. It is not a byte-for-byte dump of the private project folder.

## Not included

- Any PIC (Paediatric Intensive Care) row-level data, real-run outputs, `*_INTERNAL` files, `*_raw_LOCAL` files, or `*_TEAM` result tables from real runs.
- Cohort, split, landmark-eligibility, death, and raw-table row-count pins. In the original project these are hard-coded reproduction gates. Here the corresponding defaults are empty or `None`. Restore them locally from `configs/cohort_pins.example.json` (copy to `configs/cohort_pins.json`) before a real PIC run. SHA-256 hashes of pinned inputs are kept: they do not reveal counts.
- Per-item `train_encounters` / `train_rows` values in the channel, drug, organism and specimen maps. The mapping columns themselves are kept. The frozen φ1-medians map is shipped unchanged (including `train_n` and `median_hours`) so its SHA-256 still matches the pin the builder checks; some `train_n` cells were already suppressed as `S`. This is noted as an open question.
- The Stage 2 team-internal public-values file and the Stage 3 Stage-2-constants file. Tests that read the former are skipped when it is absent; Stage 3 tests build a synthetic constants file themselves.
- Internal planning notes and report drafts. Specs and the analysis plan/contract are included as cleaned markdown.

## How this differs from the laptop-run packages

The script SHA-256 values recorded in the Copilot instruction notes are the hashes of the **private** packages, which still contain the pins. The published `.py` files differ by this redaction (and, for one Stage 2 test, by a skip marker and a broader search path for the Stage 3 builder). Use the Copilot notes as a record of how the real runs were launched, not as a hash check against this public tree.

Source PDFs and spreadsheets of WHO / published vital-sign standards are not included (third-party copyright). The derived CSVs and their provenance notes are included.

## Personal path redaction

Local Windows paths (the project folder, the raw PIC directory, the split-file path, the desktop venv interpreter, and the username they contained) were replaced with placeholders: `<PROJECT_ROOT>`, `<PIC_RAW_DIR>`, `<SPLIT_FILE>`, `<VENV_PYTHON>`, and `<VENV_DIR>`. CLI defaults that used to hard-code those paths now read `PIC_RAW_DIR` and `SPLIT_FILE_PATH` from the environment, with no local fallback. Downstream scripts in this tree that hash-pin a sibling public script were retargeted to the published hash of that script (the hashes quoted in Copilot notes and some specifications remain the recorded hashes of the private packages). Script hashes below are of the files in this public tree after that redaction.

Changed scripts and their SHA-256 in this tree:

- `scripts/band_support_audit/phi1_band_support_audit.py` — `b127bbb37f791fb01e99cecaa161838e05f8fcc36865ca7b57c7956f85ba9804`
- `scripts/check_encounter_harmonization_consistency.py` — `60d55b888d92cc945495f094f61836d2f12811856d66e5d683125e2997f19e59`
- `scripts/duration_audit_v2/audit_surgery_complexity_duration.py` — `1de2b17b2d894caff6a8362a9c72a458967691a8d2c7a533db863e87ad8dad94`
- `scripts/duration_audit_v2/test_audit_surgery_complexity_duration.py` — `669fa3d54fc68958c97126b4064791047422ef47ae35ccf49d3badca5907272f`
- `scripts/feature_item_inventory/feature_item_inventory.py` — `6a11c2bdd935faca7d8f136718405a9d7f81339ff94c2598e01d09983309caea`
- `scripts/feature_item_inventory/test_feature_item_inventory.py` — `e1f6498e72040eecbba50e8cf0e3be13dd2de3088aed58a84ffa72597613e5e4`
- `scripts/feature_sparsity_frequency_audit/feature_sparsity_frequency_audit.py` — `ce03f97f22f93749c7d325fa8411326ca575945386134e7fcd30089e9180c902`
- `scripts/feature_sparsity_frequency_audit/test_feature_sparsity_frequency_audit.py` — `40498d72a28d86fed370b2d3dde4e82237054e8d0248fdd1d9f9a9055c710b5a`
- `scripts/stage06_08/06_08_build_levels_preprocess.py` — `cbadeb163bf875c8d40b12077879e638963e3793dd1931de6941835d6dc7c278`
- `scripts/stage2_followup_checks/audit_surgery_complexity_duration.py` — `1de2b17b2d894caff6a8362a9c72a458967691a8d2c7a533db863e87ad8dad94`
- `scripts/stage2_followup_checks/feature_sparsity_frequency_audit.py` — `ce03f97f22f93749c7d325fa8411326ca575945386134e7fcd30089e9180c902`
- `scripts/stage2_followup_checks/stage2_followup_checks.py` — `4675faf9b329c71e28a01ddce1bd63e36f4b541087cec664930925e465623e0a`
- `scripts/stage3_builder/03_build_real_layers.py` — `5fb7629678d4f4370ebd5066fb0762fb8ce3e38e4daef3d46392face842149bb`
- `scripts/stage3_builder/audit_surgery_complexity_duration.py` — `1de2b17b2d894caff6a8362a9c72a458967691a8d2c7a533db863e87ad8dad94`
- `scripts/stage3_builder/test_03_build_real_layers.py` — `027925abcac4323e6b7b2ffde48017009c6a0e6929aced1a4fa9f29ece6cf71e`
