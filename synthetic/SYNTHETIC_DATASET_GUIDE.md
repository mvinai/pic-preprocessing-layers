---
title: "Synthetic Dataset Resource Guide"
subtitle: "PIC Foundations feature matrices (Layers A / B / C)"
author: "Prepared for teammates working on SVD and landmark models"
date: "November 16, 2026"
---

# What this folder is

This is a **synthetic**, PIC-like dataset built on **15 November 2026** (generator seed `20261115`) for algorithm work: **15,000** patients and **17,239** ICU encounters (train 10,328 / validation 3,405 / test 3,506, split by `SUBJECT_ID`). It is **not real PIC data** and must not be treated as such.


# Start here — is `6_tensor_ready_layers` enough?

**Almost, but not by itself if you are fitting or scoring models.**

`6_tensor_ready_layers/` is the analysis-ready package: Layers A and B as matrices, Layer C as a **3-way tensor that is never flattened**, plus masks, a filled/imputed copy (levels 2 and 3), feature dictionaries, train-split statistics, and `encounter_index.csv` (IDs, subject-level `split`, `los_hours`). You do **not** need the raw PIC-schema tables, the v1 channel map, the root `data_dictionary.csv`, or `subject_splits.csv` — those last two are already represented inside folder 6.

You **do** still need one file from outside folder 6 for Objectives 1 and 2:

| Still needed | Where | Why |
|---|---|---|
| **`outcomes_by_landmark.csv`** | `5_layered_tensor/3_layers/` | Outcome 1 (0–8), Outcome 2 (0–4), in-hospital mortality, and eligibility at each landmark. Not stored in folder 6. |
| `landmark_index.csv` (optional) | same folder | Same eligible `(encounter, landmark)` rows without the outcome columns. Redundant if you have the outcomes file. |
| `4_answer_key/` (optional) | dataset root | **Held-out** planted factor scores and Layer C trajectories. Use only to check whether an SVD recovered the planted structure. **Never use as features.** |

**Recommended starting paths (inside folder 6):** `3_features/layer_{a,b,c}.npz` for the class analyses (clinically engineered, already imputed and standardized). Use `2_normalized/` if you want unit-harmonized z-scored channels without the extra clinical features. Use `1_raw/` only if you are studying preprocessing itself — values are native units with NaNs and **no filled copy**.

# Folder tree

Sizes below are from a complete local copy (~1.1 GB). A smaller packaged copy (~452 MB) leaves out some uncompressed long CSVs in favour of `.zip` / `.csv.gz`; the analysis files you need are all there.

```
synthetic data for optimizing SVD algorithms/
├── README_synthetic_svd_levels.md      # stages 1–4, planted static factors
├── README_layered_tensor.md            # Layers A/B/C, outcomes, LOS limit
├── data_dictionary.csv                 # for the older 0–24 h SVD matrix only
├── manifest.json
├── 1_raw/          ~91 MB   PIC-schema tables, v1 maps, splits  (event data)
├── 2_normalized/   ~142 MB  encounters, observations_long, medications
├── 3_features/     ~24 MB   0–24 h feature tables (older static SVD matrix)
├── 4_answer_key/   ~20 MB   held-out planted factors + trajectory scores
├── 5_layered_tensor/  ~665 MB (packaged copy ~214 MB)
│   ├── 1_extended_raw/         full-stay companion events + outcomes_source
│   ├── 2_normalized_extended/
│   └── 3_layers/               layer_a/b.csv, source layer_c.npz,
│                               outcomes_by_landmark.csv, landmark_index.csv
└── 6_tensor_ready_layers/      ~85 MB   ← START HERE
    ├── README_tensor_ready_layers.md
    ├── encounter_index.csv     # 17,239 rows; row order = every npz
    ├── manifest.csv            # paths, shapes, dtypes, sizes
    ├── validation.json
    ├── 1_raw/{layer_a,b,c}.npz + feature_dictionary.csv
    ├── 2_normalized/{layer_a,b,c}.npz + feature_dictionary_and_train_stats.csv
    └── 3_features/{layer_a,b,c}.npz + feature_dictionary_and_train_stats.csv
```

# What's inside folder 6

All 17,239 encounters are in the **same order** in `encounter_index.csv` and in every npz (`encounter_ids`, `splits`). Join on `encounter_id`. Layer C bins are 6 h wide from t0 to 336 h (`bin_edges`, `bin_midpoints` in the C files). `mask_valid[e, b] = 1` iff `6*b < los_hours[e]`. `mask_obs` marks charted / defined values; `mask_imputed` marks valid cells that were filled.

| Level | Layer A (static) | Layer B (first 24 h) | Layer C (tensor) |
|---|---|---|---|
| `1_raw` | (17239, 38) | (17239, 31) | (17239, 56, 31) |
| `2_normalized` | (17239, 49) | (17239, 35) | (17239, 56, 34) |
| `3_features` | (17239, 46) | (17239, 73) | (17239, 56, 76) |

Which array to load: **`1_raw`** — `X` only (NaNs kept; no filled copy; feature names are item IDs / raw drug names). **`2_normalized`** — `X` (with NaN) and `X_filled` (LOCF on C, then 0 = train mean on the z-scale). **`3_features`** — `X` (filled, original units), `X_standardized`, and `X_with_nan`; indicators stay 0/1.

Layer C at level 3 is ordered as 35 value channels (20 means + 13 drug-group flags + `vasoactive_max_concurrent` + `psofa_resp`), then 21 observed indicators, then 20 time-since-last columns. Full names, roles, and train mean / SD / median are in each level's feature dictionary CSV (and as `feature_names` / `train_*` inside the npz). Numerics are `float32`; masks are `uint8`.

**Landmark slice (do not unfold C for the intended tensor SVD).** An encounter is eligible at landmark *L* in {24, 48, …, 336} hours if `los_hours > L`. Take the first `n_bins = L // 6` time slices of the tensor (for *L* = 72 that is 12 bins). The 72 h sensitivity window uses only the last 12 bins before *L* (the full prefix when *L* ≤ 72). Existing stays end by ~280 h, so eligibility at 288 / 312 / 336 h is **zero** (this was not patched).

# Which file for which task

| Task | Use these files |
|---|---|
| **Objective 1** — primary A/B analysis at the 24 h landmark, Outcome 1 | `6_.../3_features/layer_a.npz` + `layer_b.npz` (`X_standardized`); `outcomes_by_landmark.csv` filtered to `landmark_h == 24` and `eligible == 1` (16,192 rows). |
| **Objective 2** — A+B vs A+B+C, refit at core landmarks 24/48/72/120/168/240/336 h (336 h will be empty); report all 14 | Same A/B files; `3_features/layer_c.npz` sliced as a tensor at each $L$; `outcomes_by_landmark.csv` for that $L$. Concatenate A and B (and a tensor reduction of C), do not unfold C into the primary contrast. |
| **Objective 3** — tensor SVD (Tucker / HOSVD) of Layer C | `3_features/layer_c.npz` → `X_standardized` (or `2_normalized/layer_c.npz` → `X_filled`). Restrict to a complete-stay window (e.g. `los_hours > 72`, first 12 bins) or otherwise honour `mask_valid`. Padded zeros leak length of stay if you treat them as data. |
| **Flattened Layer C baseline** (labeled baseline only) | There is **no pre-built flattened-C file**. After slicing, `X_L.reshape(n, -1)`. Do not confuse this with `3_features/features_svd_matrix.csv`, which is an older 0–24 h **static** 25-column matrix, not flattened C. |
| Recovery check (optional) | `4_answer_key/latent_factors.csv` (`z_severity`, `z_respiratory`, `z_age`); `4_answer_key/layered_tensor/trajectory_scores.csv` and `shape_basis.csv`. |

# Loading snippets (Python)

Tested against the 15 November 2026 build. Paths assume the dataset root is the current working directory.

```python
import numpy as np
import pandas as pd

# IDs, subject-level splits, and LOS — already inside folder 6
idx = pd.read_csv("6_tensor_ready_layers/encounter_index.csv")  # 17239 x 11

# Layer B, matrix SVD (Objective 1 features)
b = np.load("6_tensor_ready_layers/3_features/layer_b.npz", allow_pickle=True)
X_B = b["X_standardized"]          # (17239, 73), no NaN
U, S, Vt = np.linalg.svd(X_B, full_matrices=False)
# leading singular values: 606.07, 350.25, 293.34, ...

# Layer C, keep it a tensor. Complete 72 h window (LOS > 72, first 12 bins)
c = np.load("6_tensor_ready_layers/3_features/layer_c.npz", allow_pickle=True)
X = c["X_standardized"]            # (17239, 56, 76)
assert np.array_equal(idx["encounter_id"].to_numpy(), c["encounter_ids"])
W = X[idx["los_hours"].to_numpy() > 72][:, :12, :]   # (6038, 12, 76)

def unfold(X, mode):
    return np.moveaxis(X, mode, 0).reshape(X.shape[mode], -1)

_, s_time, _ = np.linalg.svd(unfold(W, 1), full_matrices=False)
# mode-2 (time) top 5: 1387.09, 778.46, 569.84, 450.68, 426.67

# Landmark L = 120 h, still a tensor; join outcomes by encounter_id
L = 120
elig = idx["los_hours"].to_numpy() > L
X_L = X[elig][:, : L // 6, :]      # (3857, 20, 76)
X_flat = X_L.reshape(X_L.shape[0], -1)   # labeled baseline only: (3857, 1520)

y = pd.read_csv("5_layered_tensor/3_layers/outcomes_by_landmark.csv")
y_L = (y[(y["landmark_h"] == L) & (y["eligible"] == 1)]
         .set_index("encounter_id")
         .loc[idx.loc[elig, "encounter_id"]])
# y_L["outcome1_remaining_icu"] is the primary ordinal label (0–8)
```

# Caveats

- **Synthetic only.** Event values, timestamps, and outcomes were generated; this is not PIC and is not a substitute for a later hospital dataset.
- **Predates decision D29 (and map v2.2).** The tensors still include P/F, S/F, FiO2, and the pSOFA respiratory component, which the analysis plan later dropped (analyser FiO2 likely under-recorded). Item 5210 (Aa gradient) is not in these tensors. **Keep the files as they are**; drop those oxygenation / pSOFA-respiratory columns in analysis code if you are following the locked plan.
- **Channel map is v1.** Folder 6's dictionaries are enough to interpret columns. Do not rebuild features from `1_raw/` or `5_layered_tensor/` unless you are studying the generator itself.
- **1,047 encounters never appear in the outcomes file** (LOS ≤ 24 h, ineligible at every landmark). Eligibility is `los_hours > L` in `encounter_index.csv`.
- **Do not unfold Layer C** for the intended Objective 3 analysis, and do not treat post-discharge padded zeros as observations.
- **Do not use `4_answer_key/` as model features.** Those files are the answer key.
- Outcome definitions and the planted trajectory classes are documented in `README_layered_tensor.md`. Inherited ordinal-8 is present only at 24 / 48 / 72 h and is not a primary endpoint.

Questions about column roles: start with `6_tensor_ready_layers/3_features/feature_dictionary_and_train_stats.csv` and `README_tensor_ready_layers.md`.
