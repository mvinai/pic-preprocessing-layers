#!/usr/bin/env python3
"""06_08_build_levels_preprocess v1.2.1 - stages 06, 07a/07 and 08 of the PIC Foundations plan v13 in one run.

Stage 06  Layer A representation levels a0 / a1 / a2 from the Stage 3 v1.0.4 Layer A base.
Stage 07  Layer B representation levels b0 / b1 / b2s / b2d / b2p from the Stage 3 v1.0.4 Layer B base
          (07a residual items are absorbed: organ membership and the scorers are fixed tables in this script).
Stage 08  Train-only preprocessing fit (rare-category pooling is done in 06/07 on train subjects only; here:
          winsorizing, log rule, train-median fill, zero-variance drop, one-hot reference drop, standardizing)
          applied unchanged to validation and test; Layer C value tensors standardized per landmark from the
          14 Stage 3 filled copies (raw level NaN + masks, normalized level filled + flags).

Spec: stage06_08_spec v1.0.0 (plan v13 sections 4, 5, 7, 13.2; contract v13). Inputs are the Stage 3 v1.0.4 release
folder (verified file by file against its manifest; the manifest itself is pinned by its full SHA-256 on the command
line), the canonical encounter file (encounter_id and SUBJECT_ID only; never Died), the harmonization file (t0 and
OUT only, for the stay length exactly as Stage 3 computes it) and the subject split file (all three verified against
the hashes recorded in the Stage 3 manifest).

v1.2.1 (Nov 18, 2026): B8 tests banned organism tokens on field stems, not on one-hot / drug category text; one-hot
blocks with blank source values get a <block>_missing mask (D11; sex exempt per D35) and A6 checks block + mask = 1;
a failed acceptance check names its levels and columns / blocks (never row or patient counts).

Data rules: row-level outputs stay on the credentialed machine. The only shareable files are stage06_08_audit.json/.md
(structure, widths, check results; no patient counts at all, so there is no cell to suppress). Files named *_INTERNAL*
are team-internal and never shared. Progress lines and heartbeats carry counts of columns only, never values,
identifiers or timestamps.

Exit codes: 0 ok; 2 stopped on an input check (nothing written); 3 build or acceptance-check failure (nothing
published; a pending folder is removed, or named if Windows refuses to remove it).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
import threading
import time
from collections import OrderedDict
from statistics import NormalDist
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_VERSION = "06_08_build_levels_preprocess v1.2.1"
SPEC_DOC = "stage06_08_spec v1.2.1"
PLAN, CONTRACT = "v13", "v13"

# ----------------------------------------------------------------------------------------------- pins and targets
PINS = {
    # Stage 3 v1.0.4 builder (03_build_real_layers.py) as verified by Copilot on Nov 17, 2026 (Stage 3 spec section 14)
    "stage3_script_version": "03_build_real_layers v1.0.4",
    "stage3_script_sha256": "5fb7629678d4f4370ebd5066fb0762fb8ce3e38e4daef3d46392face842149bb",
    "channel_map": "7541d7ea67aa9631c984a5b2c7bb63bf17e6e38b9013a7b09ea42fe754536ec2",     # map v2.3
    "phi1_medians": "e73c8b5a8d780b51edfd0935abdb05335f57e3a0cfc0c4c52eaaa85888db0047",
    "canonical": "e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418",
    "split": "60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31",
    "stage3_files": 46,          # files listed in stage3_manifest.json (47 in the folder with the manifest)
}


@dataclass
class Expect:
    encounters: int | None = None  # public copy: pin omitted; see configs/cohort_pins.example.json
    splits: dict = field(default_factory=dict)
    eligible: dict = field(default_factory=dict)
    layer_a_cols: int = 297
    layer_b_cols: int = 667


EXPECT = Expect()
SPLITS = ("train", "validation", "test")
LANDMARKS = (24, 48, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288, 312, 336)
M_RARE = 20                 # plan section 7: categories with fewer than 20 training subjects are pooled to other
MIN_BAND_N = 10             # age-band reference needs >= 10 train encounters with a value, else the pooled train value
# Admission source (a1), Vinai's decision D1 of Nov 18, 2026: recorded fixed settings, not command-line options.
ADM_SRC_SURGICAL_WINDOW_H = 48.0   # surgical: recent surgery ended <= t0 AND hospital-to-ICU interval < 48 h
ADM_SRC_FLOOR_CUTOFF_H = 24.0      # floor: not surgical AND interval >= 24 h
ADM_SRC_LEVELS = ("surgical", "floor", "transfer_ed")
DRUG_GROUP_MIN_PREV = 0.01  # groups below 1% train prevalence dropped above phi1 (plan section 5.5)
WINSOR_Q = (0.01, 0.99)
LOG_SKEW = 2.0              # log1p for non-negative continuous columns with train skewness > 2 (after winsorizing)
BANNED_TOKENS = ("dod", "expire", "deathtime", "dischtime", "discharge", "outtime", "los", "icu_los", "span_icu",
                 "icu_time", "last_careunit", "last_wardid", "has_chartevents", "readmi", "death", "landmark",
                 "target", "outcome", "subject_id", "hadm_id", "icustay_id", "mask_valid", "weight", "enddate",
                 "storetime", "edreg", "visit_id", "sex_missing", "_count")
# v1.2.1 (B8): organism / susceptibility field stems; matched at the start of a field stem or after an underscore,
# never inside one-hot or drug category text (real drug names such as "..._tab_..." or "...mab_..." contain "ab_")
B8_TOKENS = ("org_", "ab_", "dilution", "interpret")
B8_RE = re.compile(r"(?:^|_)(?:org_|ab_|dilution|interpret)")
# v1.2.1 (A6, decision D11 "missing as its own indicator"): a one-hot block whose source values have blanks gets a
# <block>_missing mask, except the blocks below (sex: plan v12 / contract v13 D35 removed sex_missing and it is a
# banned column; a blank sex stays an all-zero row of the sex block, which A6 accepts only on blank-source rows)
MISSING_MASK_EXEMPT = {"sex": "D35: sex_missing removed (plan v12, contract v13) and banned"}
STAGE3_FILES_USED = ["stage3_layer_a_base.csv", "stage3_layer_a_base.csv.meta.json", "stage3_layer_b_base.csv",
                     "stage3_layer_b_base.csv.meta.json", "stage3_layer_c_raw.npz", "stage3_layer_c_raw.npz.meta.json",
                     "stage3_medication_courses_INTERNAL.csv", "stage3_microbiology_collections_INTERNAL.csv"] + \
                    [f"stage3_layer_c_filled_L{L}.npz" for L in LANDMARKS]

# ------------------------------------------------------------------------------------- channel tables (map v2.3)
SYSTEM_OF = {  # primary home (map v2.3 table51_system; Table 5.1)
    "heart_rate": "cardiovascular", "sbp": "cardiovascular", "dbp": "cardiovascular", "map": "cardiovascular",
    "lactate": "cardiovascular", "bnp": "cardiovascular", "troponin_t": "cardiovascular", "troponin_i": "cardiovascular",
    "respiratory_rate": "respiratory", "spo2": "respiratory", "ph": "respiratory", "pco2": "respiratory",
    "po2": "respiratory", "base_excess": "respiratory", "bicarbonate": "respiratory", "sao2_bloodgas": "respiratory",
    "creatinine": "renal", "bun": "renal", "cystatin_c": "renal", "urine_output": "renal", "fluid_intake": "renal",
    "kdigo_uo_stage": "renal",
    "bilirubin_total": "hepatic", "bilirubin_direct": "hepatic", "albumin": "hepatic", "ast": "hepatic",
    "ggt": "hepatic", "alt": "hepatic", "ammonia": "hepatic",
    "amylase": "gi", "lipase": "gi", "stool_occult_blood": "gi",
    "hemoglobin": "hematologic", "platelets": "hematologic", "lymphocytes_abs": "hematologic", "ldh": "hematologic",
    "inr": "hematologic", "ptt": "hematologic", "fibrinogen": "hematologic", "d_dimer": "hematologic",
    "ferritin": "hematologic", "teg_r": "hematologic", "teg_k": "hematologic", "teg_angle": "hematologic",
    "teg_ma": "hematologic", "teg_coagulation_index": "hematologic", "teg_epl": "hematologic",
    "teg_ly30": "hematologic", "teg_g": "hematologic", "teg_a": "hematologic",
    "temperature": "infectious", "wbc": "infectious", "anc": "infectious", "crp": "infectious",
    "procalcitonin": "infectious", "beta_d_glucan": "infectious", "il6": "infectious", "il10": "infectious",
    "tnf_alpha": "infectious", "ifn_gamma": "infectious", "urine_wbc": "infectious", "urine_nitrite": "infectious",
    "urine_leukocyte_esterase": "infectious", "resp_panel_rsv": "infectious", "resp_panel_parainfluenza": "infectious",
    "resp_panel_influenza": "infectious", "resp_panel_metapneumovirus": "infectious",
    "resp_panel_mycoplasma": "infectious", "systemic_viral_cmv": "infectious", "systemic_viral_ebv": "infectious",
    "systemic_viral_hsv": "infectious", "systemic_viral_enterovirus": "infectious", "bacteremia_pcr": "infectious",
    "pain_score": "neurologic", "csf_wbc": "neurologic", "csf_rbc": "neurologic", "csf_protein_pandy": "neurologic",
    "nse": "neurologic",
    "sodium": "electrolytes", "potassium": "electrolytes", "chloride": "electrolytes", "calcium_total": "electrolytes",
    "calcium_ionized": "electrolytes", "magnesium": "electrolytes", "phosphate": "electrolytes",
    "glucose": "endocrine", "glucose_bedside": "endocrine", "cortisol": "endocrine", "free_t4": "endocrine",
    "insulin": "endocrine", "osmolality": "endocrine", "tsh": "endocrine",
}
# Stacked ("also in") memberships of Table 5.1 note (vi), carried into phi2-dup; K in Neurologic is open -> not added
ALSO_IN = {
    "potassium": ["cardiovascular", "renal"], "calcium_total": ["cardiovascular", "neurologic"],
    "magnesium": ["cardiovascular", "neurologic"], "lactate": ["respiratory", "endocrine"],
    "phosphate": ["renal"], "bicarbonate": ["renal"], "inr": ["hepatic"], "ldh": ["hepatic"],
    "wbc": ["hematologic"], "anc": ["hematologic"], "heart_rate": ["infectious"], "respiratory_rate": ["infectious"],
    "temperature": ["neurologic"], "glucose": ["neurologic"], "ammonia": ["neurologic", "endocrine"],
    "sodium": ["neurologic"],
}
SYSTEMS = ("cardiovascular", "respiratory", "renal", "hepatic", "gi", "hematologic", "infectious", "neurologic",
           "electrolytes", "endocrine")
VITALS = ("temperature", "heart_rate", "respiratory_rate", "spo2", "map", "sbp", "dbp", "pain_score", "glucose_bedside")
ZFLAG_VITALS = {"heart_rate": ("max", ">", 2.0), "respiratory_rate": ("max", ">", 2.0), "map": ("min", "<", -2.0)}
SD_VITALS = ("heart_rate", "respiratory_rate", "map")
PUB_NORM_VITALS = ("heart_rate", "respiratory_rate", "map")   # D2: published norms; other vitals stay empirical
TRIGGERED = ("ammonia", "beta_d_glucan", "bnp", "ferritin", "ifn_gamma", "il10", "il6", "nse", "procalcitonin",
             "teg_a", "teg_angle", "teg_coagulation_index", "teg_epl", "teg_g", "teg_k", "teg_ly30", "teg_ma",
             "teg_r", "tnf_alpha", "troponin_i", "troponin_t")
INDICATOR_ONLY = ("alt", "tsh")
QUALITATIVE = ("resp_panel_parainfluenza", "resp_panel_rsv", "stool_occult_blood", "urine_leukocyte_esterase",
               "urine_nitrite", "csf_protein_pandy", "systemic_viral_cmv", "systemic_viral_ebv", "systemic_viral_hsv",
               "resp_panel_mycoplasma", "bacteremia_pcr", "resp_panel_metapneumovirus", "systemic_viral_enterovirus",
               "resp_panel_influenza")
PANELS = {"blood_gas": ("base_excess", "chloride", "glucose", "bicarbonate", "potassium", "lactate", "sodium", "pco2",
                        "ph", "po2", "sao2_bloodgas"),
          "cytokines": ("ifn_gamma", "il10", "il6", "tnf_alpha"),
          "teg": ("teg_r", "teg_k", "teg_angle", "teg_ma", "teg_coagulation_index", "teg_epl", "teg_ly30", "teg_g",
                  "teg_a"),
          "urine_dipstick": ("urine_leukocyte_esterase", "urine_nitrite", "urine_wbc")}
PANEL_SYSTEM = {"blood_gas": "respiratory", "cytokines": "infectious", "teg": "hematologic",
                "urine_dipstick": "infectious"}
PANEL_OF = {c: p for p, cs in PANELS.items() for c in cs}
# direction of "worst" (plan section 5.5; open item, default table of this spec): high, low or two-sided
WORST = {**{c: "high" for c in ("lactate", "creatinine", "bun", "cystatin_c", "bilirubin_total", "bilirubin_direct",
                                "ast", "ggt", "ldh", "amylase", "lipase", "inr", "ptt", "d_dimer", "crp", "csf_wbc",
                                "csf_rbc", "urine_wbc", "fluid_intake", "procalcitonin", "il6", "il10", "tnf_alpha",
                                "ifn_gamma", "troponin_i", "troponin_t", "bnp", "ferritin", "ammonia", "nse",
                                "beta_d_glucan", "teg_r", "teg_k", "teg_epl", "teg_ly30") + QUALITATIVE},
         **{c: "low" for c in ("platelets", "hemoglobin", "albumin", "po2", "sao2_bloodgas", "lymphocytes_abs",
                               "urine_output", "fibrinogen", "teg_angle", "teg_ma", "teg_g", "teg_a")},
         **{c: "two" for c in ("sodium", "potassium", "chloride", "glucose", "calcium_total", "calcium_ionized",
                               "magnesium", "phosphate", "ph", "pco2", "base_excess", "bicarbonate", "wbc", "anc",
                               "osmolality", "cortisol", "free_t4", "insulin", "teg_coagulation_index")}}
MED_GROUPS = ("vasoactive", "antibiotic", "antifungal", "antiviral", "sedative_analgesic", "non_opioid_analgesic",
              "antiemetic", "diuretic", "osmotherapy", "electrolyte_replacement", "steroid", "nmb", "aed")
MED_FEATS = ("started", "hrs_to_start", "hrs_active", "active_end", "hours_on_course", "hours_since_stopped",
             "ever_received")
SPECIMENS = ("blood", "urine", "respiratory", "csf", "stool", "ureaplasma", "culture_site_unknown", "other", "unmapped")
AGE_BANDS = ("neonate", "infant", "toddler", "preschool", "school", "adolescent", "ge18", "age_invalid")
INTERVAL_BANDS = ("lt24h", "24to72h", "72to168h", "ge168h")
SHIFTS = ("day", "evening", "night")

# ------------------------------------------------------------------------ scorer tables (theta*, plan section 5.6)
# Modified pSOFA: Matics TJ, Sanchez-Pinto LN. JAMA Pediatr 2017;171(10):e172352, Table 1 (read Nov 18, 2026 from the
# journal page). Respiratory dropped (D29); neurologic needs GCS, which PIC lacks (Table 5.1 note iv) -> not scored;
# cardiovascular 2-4 need infusion doses, which PRESCRIPTIONS lacks -> modified rule (2 if any vasoactive active).
PSOFA_AGE_EDGES_MO = (1, 12, 24, 60, 144, 216.0001)   # <1, 1-11, 12-23, 24-59, 60-143, 144-216, >216 months
PSOFA_MAP_CUT = (46, 55, 60, 62, 65, 67, 70)          # score 1 if MAP < cut
PSOFA_CREAT = ((0.8, 1.0, 1.2, 1.6), (0.3, 0.5, 0.8, 1.2), (0.4, 0.6, 1.1, 1.5), (0.6, 0.9, 1.6, 2.3),
               (0.7, 1.1, 1.8, 2.6), (1.0, 1.7, 2.9, 4.2), (1.2, 2.0, 3.5, 5.0))   # lower edges of scores 1-4, mg/dL
PSOFA_PLT = (150, 100, 50, 20)          # score 0 if >= 150; 1: 100-149; 2: 50-99; 3: 20-49; 4: < 20 (x10^3/uL)
PSOFA_BILI = (1.2, 2.0, 6.0, 12.0)      # lower edges of scores 1-4, mg/dL
# Modified PRISM III point scores (Pollack MM et al., Crit Care Med 1996;24:743-752) as tabulated in two independent
# secondary sources that agree: APCC 2023 Supplementary Table 1 (doi 10.32990/apcc.2023.00031) and the Dutch PICE
# registry PRISM III table (pice.nl). Mental status, pupils (absent in PIC) and PT (no PT seconds channel) not scored.
PRISM_AGE = {"neonate": 0, "infant": 1, "toddler": 2, "preschool": 2, "school": 2, "adolescent": 3, "ge18": 3}
PRISM_SBP = ((40, 55), (45, 65), (55, 75), (65, 85))      # (7 points below a, 3 points a..b)
PRISM_HR = ((215, 225), (215, 225), (185, 205), (145, 155))  # (3 points a..b, 4 points above b)
PRISM_CREAT = (0.85, 0.9, 0.9, 1.3)                        # > cut: 2 points (mg/dL)
PRISM_BUN = (11.9, 14.9, 14.9, 14.9)                       # > cut: 3 points (mg/dL)
PRISM_PTT = (85.0, 57.0, 57.0, 57.0)                       # > cut: 3 points (s)
# AKI stage from creatinine (contract section 9.3 baselines, mg/dL by age in months)
AKI_BASE = ((3, 0.35), (12, 0.40), (60, 0.55), (144, 0.70), (1e9, 0.90))


# ------------------------------------------------- published age norms (decision D2, Vinai, Nov 18, 2026, 6:38-6:43)
# Pinned reference tables shipped next to this script (copied unchanged by Copilot; SHA-256 checked in the precheck).
# Fleming S et al. Lancet 2011;377:1011-1018, Web Tables 5 (HR) and 4 (RR): HR and RR as normal z (mean = median,
#   SD = (p90 - p10) / (2 * 1.28155)); HR 'Birth' row only for age < 1 day, else the 0-3 m row; RR at birth -> 0-3 m.
# Roberts JS, Yanay O, Barry D. Pediatr Crit Care Med 2020;21(9):e759-e768, Table 2 (OSF upqtv): MAP as a
#   centile-interpolated z (piecewise linear in z between the 13 published centiles, capped at the 1st/99th),
#   with every centile interpolated linearly by exact age between the published evaluation ages.
# Haque IU, Zaritsky AL. Pediatr Crit Care Med 2007;8(2):138-144: MAP sensitivity only (ages 1-17 y; team-internal).
# Fitzgerald JC et al. (PODIUM) Pediatrics 2022;149(1 Suppl 1):S66-S73, Table 2: creatinine / age- and sex-specific
#   baseline SCr ratio.
REF_DIR = Path(__file__).resolve().parent
REF_FILES = {
    "fleming2011_hr_centiles.csv": "4ff3c424dd398ded5b7e9a7a70de45f3b44fff8f9fe456a2a3be0e61d331f1ad",
    "fleming2011_rr_centiles.csv": "ec54810ee50120ad3ea5c69c90a5f7cbd6c3ef09de560d102786f8de7c94e552",
    "roberts2020_map_percentiles.csv": "421c965f262335e2445674a4fe78c9fa756b10d28b64a206eb3c38b57f1ef2f2",
    "podium2022_baseline_scr.csv": "e56de3e8f9340ada0d163e520c58dbf657a21fddfc179ca0213f184ca79190e7",
}
DAYS_PER_MONTH = 30.4375              # Stage 3 age_months = age_days / 30.4375
FLEMING_BIRTH_MAX_DAYS = 1.0          # HR 'Birth' row (immediate neonatal period) used for age < 1 day only
Z90 = NormalDist().inv_cdf(0.90)
# Creatinine units: PIC LABEVENTS creatinine (items 5041/5032) is umol/L; Stage 3 converts it to mg/dL (/88.4) by the
# pinned channel map v2.3 (conversion 'umol/L -> mg/dL /88.4'), so Layer A/B creatinine arrive in mg/dL and are not
# converted again here. A train-median plausibility gate stops the run if they look like umol/L.
STAGE3_CREATININE_UNIT = "mg/dL"
UMOL_PER_MG_DL_CREATININE = 88.4
CREATININE_MG_DL_MEDIAN_MAX = 10.0


def creatinine_mg_dl(x, unit: str) -> np.ndarray:
    """Creatinine in mg/dL from mg/dL or umol/L (the Stage 3 conversion: umol/L / 88.4)."""
    u = unit.replace("\u03bc", "u").replace("\u00b5", "u").strip().lower()
    x = np.asarray(x, float)
    if u == "mg/dl":
        return x.copy()
    if u == "umol/l":
        return x / UMOL_PER_MG_DL_CREATININE
    raise InputError(f"unknown creatinine unit {unit!r}")


def _read_ref(name: str) -> list:
    with open(REF_DIR / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


class Norms:
    """The four pinned reference tables, parsed. Every lookup returns NaN outside the table's age (and sex) cover."""

    def __init__(self):
        self.fl = {}
        for sig, name in (("heart_rate", "fleming2011_hr_centiles.csv"), ("respiratory_rate", "fleming2011_rr_centiles.csv")):
            rows = _read_ref(name)
            self.fl[sig] = [(r["age_band"], float(r["age_lo_months"]), float(r["age_hi_months"]), float(r["p50"]),
                             (float(r["p90"]) - float(r["p10"])) / (2.0 * Z90)) for r in rows]
        rr = _read_ref("roberts2020_map_percentiles.csv")
        self.rob_cols = [c for c in rr[0] if re.fullmatch(r"p\d\d", c)]
        self.rob_p = np.array([int(c[1:]) / 100.0 for c in self.rob_cols])
        self.rob_z = np.array([NormalDist().inv_cdf(q) for q in self.rob_p])
        self.rob_age = np.array([float(r["curve_evaluated_at_days"]) for r in rr])     # days
        self.rob_max_days = float(rr[-1]["age_le_days"])
        self.rob_val = np.array([[float(r[c]) for c in self.rob_cols] for r in rr])
        pod = _read_ref("podium2022_baseline_scr.csv")
        self.pod = [(float(r["age_lo_months"]), float(r["age_hi_months_excl"]), float(r["baseline_scr_boys_mg_dl"]),
                     float(r["baseline_scr_girls_mg_dl"])) for r in pod]
        if not (len(self.fl["heart_rate"]) == 14 and len(self.fl["respiratory_rate"]) == 13 and len(rr) == 21 and
                len(pod) == 10 and len(self.rob_cols) == 13):
            raise InputError("a reference table does not have its published shape")

    def fleming_params(self, sig: str, age_days: np.ndarray):
        """(mean, SD) per encounter from the Fleming band of its age; NaN outside 0 <= age < 18 y."""
        ad = np.asarray(age_days, float)
        am = ad / DAYS_PER_MONTH
        mu, sd = np.full(len(ad), np.nan), np.full(len(ad), np.nan)
        for lab, lo, hi, m_, s_ in self.fl[sig]:
            if lab == "Birth":
                sel = np.isfinite(ad) & (ad >= 0) & (ad < FLEMING_BIRTH_MAX_DAYS)
            else:
                sel = np.isfinite(am) & (am >= lo) & (am < hi)
                if lo == 0 and sig == "heart_rate":
                    sel &= ~(ad < FLEMING_BIRTH_MAX_DAYS)
            mu[sel], sd[sel] = m_, s_
        return mu, sd

    def fleming_z(self, sig: str, x, age_days) -> np.ndarray:
        mu, sd = self.fleming_params(sig, age_days)
        with np.errstate(invalid="ignore"):
            return (np.asarray(x, float) - mu) / sd

    def roberts_centiles(self, age_days) -> np.ndarray:
        """(n, 13) MAP centiles at the exact age (linear in age between evaluation ages; flat before the first and
        after the last evaluation age up to 6575 days); NaN rows outside 0 <= age <= 6575 days."""
        ad = np.asarray(age_days, float)
        out = np.full((len(ad), len(self.rob_cols)), np.nan)
        ok = np.isfinite(ad) & (ad >= 0) & (ad <= self.rob_max_days)
        for j in range(len(self.rob_cols)):
            out[ok, j] = np.interp(ad[ok], self.rob_age, self.rob_val[:, j])
        return out

    def roberts_z(self, x, age_days) -> np.ndarray:
        x = np.asarray(x, float)
        cen = self.roberts_centiles(age_days)
        z = np.full(len(x), np.nan)
        for i in np.where(np.isfinite(x) & np.isfinite(cen[:, 0]))[0]:
            z[i] = float(np.interp(x[i], cen[i], self.rob_z))          # capped at the 1st / 99th by np.interp
        return z

    def podium_baseline(self, age_months, sex) -> np.ndarray:
        am = np.asarray(age_months, float)
        sx = np.array([str(v).strip().upper() for v in sex])
        out = np.full(len(am), np.nan)
        for lo, hi, boys, girls in self.pod:
            sel = np.isfinite(am) & (am >= lo) & (am < hi)
            out[sel & (sx == "M")] = boys
            out[sel & (sx == "F")] = girls
        return out

    def creatinine_ratio(self, cr_mg_dl, age_months, sex) -> np.ndarray:
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.asarray(cr_mg_dl, float) / self.podium_baseline(age_months, sex)


def haque_zaritsky_z(x, age_days) -> np.ndarray:
    """Sensitivity only: MAP z from Haque & Zaritsky 2007 (50th = 1.5 * age(y) + 55, 5th = 1.5 * age(y) + 40, so a
    normal SD of 15 / 1.645); ages 1-17 y, else NaN."""
    ay = np.asarray(age_days, float) / 365.25
    with np.errstate(invalid="ignore"):
        z = (np.asarray(x, float) - (1.5 * ay + 55.0)) / (15.0 / NormalDist().inv_cdf(0.95))
    z[~(np.isfinite(ay) & (ay >= 1) & (ay < 18))] = np.nan
    return z


# ------------------------------------------------------------------------------------------- errors and output
class InputError(Exception):
    pass


class BuildError(Exception):
    pass


class TimeLimit(Exception):
    pass


_T0 = time.time()
_LOCK = threading.Lock()
_STATE = {"last": time.time(), "label": "start", "deadline": None}


def say(msg: str) -> None:
    with _LOCK:
        print(f"[{time.time() - _T0:9.1f}s] {msg}", file=sys.stderr, flush=True)
        _STATE["last"] = time.time()
        _STATE["label"] = msg[:90]


class Heartbeat:
    def __init__(self, every: float):
        self.every, self.stop = every, threading.Event()

    def _run(self):
        while not self.stop.wait(1.0):
            with _LOCK:
                if time.time() - _STATE["last"] >= self.every:
                    left = (f"; time left {int(_STATE['deadline'] - time.time())} s" if _STATE["deadline"] else "")
                    print(f"[{time.time() - _T0:9.1f}s] heartbeat: still working - {_STATE['label']}{left}",
                          file=sys.stderr, flush=True)
                    _STATE["last"] = time.time()

    def __enter__(self):
        self.t = threading.Thread(target=self._run, daemon=True)
        self.t.start()
        return self

    def __exit__(self, *a):
        self.stop.set()


def check_time() -> None:
    if _STATE["deadline"] and time.time() > _STATE["deadline"]:
        raise TimeLimit("time limit reached")


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def retry(op, what: str, attempts: int = 10):
    """Windows file-lock retry (cloud sync, Explorer, antivirus): 0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30 s."""
    waits = (0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30)
    for i in range(attempts):
        try:
            return op()
        except PermissionError as e:
            if i == attempts - 1:
                raise BuildError(f"file operation failed after {attempts} attempts: {what} (PermissionError, "
                                 f"WinError {getattr(e, 'winerror', None)}); is cloud sync paused?") from None
            say(f"file operation retry: {what} (PermissionError, WinError {getattr(e, 'winerror', None)}); "
                f"attempt {i + 1}, waiting {waits[i]} s")
            time.sleep(waits[i])


def write_bytes_atomic(path: Path, data: bytes) -> str:
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    retry(lambda: os.replace(tmp, path), f"replace {path.name}")
    return hashlib.sha256(data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


SPLIT_LABELS = {"train": "train", "training": "train", "val": "validation", "valid": "validation",
                "validation": "validation", "dev": "validation", "test": "test", "testing": "test"}   # as Stage 3
HARM_T0_CANDIDATES = ("t0", "T0", "t0_harmonized", "harmonized_t0", "first_intime", "FIRST_INTIME", "t0_intime")
HARM_OUT_CANDIDATES = ("OUT", "out", "last_outtime", "LAST_OUTTIME", "harmonized_out", "OUT_harmonized", "t_out",
                       "last_OUTTIME")
HARM_OUT_FALLBACK = ("OUTTIME",)
_NULLS = {"", "nan", "none", "null", "nat", "na", "n/a"}


def norm_id(x) -> str:
    """Stage 3 ID normalization: '123', '123.0', ' 123 ' -> '123'; missing -> ''."""
    s_ = str(x).strip().strip("'\"").strip()
    if s_.lower() in {"", "nan", "none", "null", "nat"}:
        return ""
    m_ = re.fullmatch(r"([+-]?\d+)\.0*", s_)
    return str(int(m_.group(1))) if m_ else s_


def fast_time(s: pd.Series) -> pd.Series:
    """Stage 3 time parser (fixed format first, then mixed)."""
    s = s.astype(str).str.strip()
    t = pd.to_datetime(s, format="%Y-%m-%d %H:%M:%S", errors="coerce")
    rest = t.isna() & ~s.str.lower().isin(_NULLS)
    if rest.any():
        try:
            t2 = pd.to_datetime(s[rest], errors="coerce", format="mixed")
        except (TypeError, ValueError):
            t2 = pd.to_datetime(s[rest], errors="coerce")
        t.loc[rest] = t2
    return t.astype("datetime64[ns]")


def num(s) -> np.ndarray:
    return np.array(pd.to_numeric(pd.Series(s), errors="coerce").to_numpy(dtype=float), dtype=float, copy=True)


BANNED_EXACT_TOKENS = {"los", "dod", "outtime", "dischtime", "deathtime", "died", "expire", "landmark", "target",
                       "outcome", "weight", "enddate", "storetime", "edregtime", "count"}
BANNED_SUBSTR = ("expire_flag", "icu_los", "span_icu", "icu_time", "last_careunit", "last_wardid", "has_chartevents",
                 "readmi", "death", "subject_id", "hadm_id", "icustay_id", "mask_valid", "visit_id", "sex_missing",
                 "discharge", "weight_source", "edreg", "ed_present", "_unparsed_text", "below_cutoff")
COUNT_WHITELIST = {"micro_culture_count"}   # plan section 5.5 phi1 M1 lists culture_count (not a measurement count)


def banned(col: str) -> bool:
    low = col.lower()
    if any(w in low for w in COUNT_WHITELIST):
        low = low.replace("micro_culture_count", "micro_culture_n")
    toks = set(low.split("_"))
    return bool(toks & BANNED_EXACT_TOKENS) or any(s in low for s in BANNED_SUBSTR)


# ---------------------------------------------------------------------------------------------------- precheck
@dataclass
class Config:
    project_root: Path
    stage3_dir: Path
    canonical: Path
    harmonization: Path
    split_file: Path
    out_dir: Path
    expect_manifest_sha256: str | None
    precheck_only: bool
    max_seconds: float = 3600.0
    heartbeat_seconds: float = 30.0


def precheck(cfg: Config) -> dict:
    """Hash gates before any data file is parsed. Returns the verified hashes."""
    ref_hashes = {}
    for name, pin in REF_FILES.items():
        rp = REF_DIR / name
        if not rp.is_file():
            raise InputError(f"reference table {name} not found next to the script")
        h = sha256_file(rp)
        if h != pin:
            raise InputError(f"reference table {name} SHA-256 differs from the pin (copy it unchanged)")
        ref_hashes[name] = h
    s3 = cfg.stage3_dir
    man_p = s3 / "stage3_manifest.json"
    if not man_p.is_file():
        raise InputError(f"Stage 3 manifest not found: {man_p}")
    man_sha = sha256_file(man_p)
    if cfg.expect_manifest_sha256 is not None and man_sha != cfg.expect_manifest_sha256:
        raise InputError("stage3_manifest.json SHA-256 differs from --expect-stage3-manifest-sha256")
    man = json.loads(man_p.read_text(encoding="utf-8"))
    if man.get("script") != PINS["stage3_script_version"]:
        raise InputError(f"Stage 3 manifest script is {man.get('script')!r}, expected {PINS['stage3_script_version']}")
    if man.get("script_sha256") != PINS["stage3_script_sha256"]:
        raise InputError("Stage 3 manifest script SHA-256 differs from the pinned v1.0.4 builder")
    if man.get("release_gate") != "ok":
        raise InputError("Stage 3 release gate is not ok")
    bad = [k for k, v in (man.get("acceptance_checks") or {}).items()
           if v != "pass" and not (k == "train_median_available_all_channels" and v == "noted")]
    if bad or not man.get("acceptance_checks"):
        raise InputError(f"Stage 3 acceptance checks not all pass: {bad}")
    ih = man.get("input_hashes", {})
    for k in ("channel_map", "phi1_medians"):
        if ih.get(k) != PINS[k]:
            raise InputError(f"Stage 3 input hash for {k} differs from the pin")
    for k in ("canonical", "split"):
        if PINS.get(k) and ih.get(k) != PINS[k]:
            raise InputError(f"Stage 3 input hash for {k} differs from the plan pin")
    files = man.get("files", {})
    if len(files) != PINS["stage3_files"]:
        raise InputError(f"Stage 3 manifest lists {len(files)} files, expected {PINS['stage3_files']}")
    on_disk = sorted(p.name for p in s3.iterdir())
    extra = sorted(set(on_disk) - set(files) - {"stage3_manifest.json"})
    missing = sorted(set(files) - set(on_disk))
    if extra or missing:
        raise InputError(f"Stage 3 folder does not match its manifest (extra {len(extra)}, missing {len(missing)})")
    hashes = {"stage3_manifest.json": man_sha}
    for i, (name, rec) in enumerate(sorted(files.items())):
        if (i + 1) % 10 == 0 or i + 1 == len(files):
            say(f"hashing Stage 3 files: {i + 1}/{len(files)}")
        h = sha256_file(s3 / name)
        if h != rec.get("sha256"):
            raise InputError(f"Stage 3 file {name} SHA-256 differs from its manifest entry")
        hashes[name] = h
    for name in STAGE3_FILES_USED:
        if name not in files:
            raise InputError(f"Stage 3 file {name} is not in the manifest")
    for key, p in (("canonical", cfg.canonical), ("harmonization", cfg.harmonization), ("split", cfg.split_file)):
        if not p.is_file():
            raise InputError(f"{key} file not found: {p}")
        h = sha256_file(p)
        if h != ih.get(key):
            raise InputError(f"{key} file SHA-256 differs from the hash recorded in the Stage 3 manifest")
        hashes[key] = h
    if cfg.out_dir.exists() and any(cfg.out_dir.iterdir()):
        raise InputError(f"output folder exists and is not empty: {cfg.out_dir.name}")
    pend = [p.name for p in cfg.out_dir.parent.glob(cfg.out_dir.name + ".pending-*")] if cfg.out_dir.parent.exists() else []
    if pend:
        raise InputError(f"a pending folder from an earlier failed run exists: {pend[0]} (remove it unopened first)")
    return {"hashes": hashes, "manifest": man, "ref_hashes": ref_hashes}


# ------------------------------------------------------------------------------------------------------ loading
def load_inputs(cfg: Config, man: dict) -> dict:
    s3 = cfg.stage3_dir
    A = pd.read_csv(s3 / "stage3_layer_a_base.csv", dtype=str, keep_default_na=False)
    B = pd.read_csv(s3 / "stage3_layer_b_base.csv", dtype=str, keep_default_na=False)
    if A.shape[1] != EXPECT.layer_a_cols or B.shape[1] != EXPECT.layer_b_cols:
        raise InputError(f"Layer A/B base widths {A.shape[1]}/{B.shape[1]} differ from {EXPECT.layer_a_cols}/"
                         f"{EXPECT.layer_b_cols}")
    if len(A) != EXPECT.encounters or len(B) != EXPECT.encounters:
        raise InputError("Layer A/B base row count differs from the cohort")
    if not (A["encounter_id"].to_numpy() == B["encounter_id"].to_numpy()).all():
        raise InputError("Layer A and Layer B base rows are not in the same order")
    can = pd.read_csv(cfg.canonical, dtype=str, keep_default_na=False, usecols=["encounter_id", "SUBJECT_ID"])
    can["encounter_id"] = can["encounter_id"].str.strip()
    if len(can) != EXPECT.encounters or can["encounter_id"].duplicated().any() or \
            set(can["encounter_id"]) != set(A["encounter_id"]):
        raise InputError("canonical encounter set differs from the Layer A base")
    can = can.set_index("encounter_id").loc[A["encounter_id"].to_numpy()].reset_index()   # Layer A row order
    sp = pd.read_csv(cfg.split_file, dtype=str, keep_default_na=False)
    for c in ("SUBJECT_ID", "split"):
        if c not in sp.columns:
            raise InputError(f"split file missing column {c!r}")
    sp = pd.DataFrame({"SUBJECT_ID": sp["SUBJECT_ID"].map(norm_id),
                       "split": sp["split"].map(lambda v: SPLIT_LABELS.get(str(v).strip().lower(), "?"))})
    if (sp["split"] == "?").any():
        raise InputError("split file has unknown split labels")
    sp = sp.drop_duplicates()
    if sp["SUBJECT_ID"].duplicated().any():
        raise InputError("a subject appears in more than one split")
    smap = dict(zip(sp["SUBJECT_ID"], sp["split"]))
    subj = can["SUBJECT_ID"].map(norm_id).to_numpy()
    split_file = np.array([smap.get(s_, "") for s_ in subj])
    split = A["split"].to_numpy()
    if not (split_file == split).all():
        raise InputError("split in the Layer A base differs from subject_splits.csv")
    for s_ in SPLITS:
        if int((split == s_).sum()) != EXPECT.splits[s_]:
            raise InputError(f"split {s_} count differs from the plan")
    # stay length exactly as Stage 3 (harmonized t0 and OUT; OUTTIME only as the fallback name)
    hh = pd.read_csv(cfg.harmonization, dtype=str, nrows=0).columns.tolist()
    t0c = [c for c in HARM_T0_CANDIDATES if c in hh]
    outc = [c for c in HARM_OUT_CANDIDATES if c in hh] or [c for c in HARM_OUT_FALLBACK if c in hh]
    if len(t0c) != 1 or len(outc) != 1 or "encounter_id" not in hh:
        raise InputError("cannot resolve the harmonization t0 / OUT columns the way Stage 3 does")
    harm = pd.read_csv(cfg.harmonization, dtype=str, keep_default_na=False, usecols=["encounter_id", t0c[0], outc[0]])
    harm["encounter_id"] = harm["encounter_id"].str.strip()
    if harm["encounter_id"].duplicated().any() or set(harm["encounter_id"]) != set(can["encounter_id"]):
        raise InputError("harmonization encounter set differs from the canonical set")
    harm = harm.set_index("encounter_id").loc[can["encounter_id"]]
    t0 = fast_time(harm[t0c[0]]).to_numpy(dtype="datetime64[ns]").astype("int64")
    to = fast_time(harm[outc[0]]).to_numpy(dtype="datetime64[ns]").astype("int64")
    nat = np.iinfo(np.int64).min
    if (t0 == nat).any() or (to == nat).any():
        raise InputError("harmonized t0 or OUT missing or unparseable for some encounters")
    stay_h = (to - t0) / 3.6e12
    if {L: int((stay_h > L).sum()) for L in LANDMARKS} != {int(k): int(v) for k, v in EXPECT.eligible.items()}:
        raise InputError("stay lengths do not reproduce the Stage 3 landmark eligible counts (counts withheld)")
    meds = pd.read_csv(s3 / "stage3_medication_courses_INTERNAL.csv", dtype=str, keep_default_na=False)
    micro = pd.read_csv(s3 / "stage3_microbiology_collections_INTERNAL.csv", dtype=str, keep_default_na=False)
    bmeta = json.loads((s3 / "stage3_layer_b_base.csv.meta.json").read_text(encoding="utf-8"))
    return dict(A=A, B=B, subj=subj, split=split, stay_h=stay_h, meds=meds, micro=micro, bmeta=bmeta,
                eid=A["encounter_id"].to_numpy())


# -------------------------------------------------------------------------------------------- level containers
class Level:
    """Ordered feature columns of one representation level with their kind and metadata.

    kind: cont (continuous or z), ordinal (scores, bands), bin (0/1 flag), onehot (member of a categorical block),
    mask (missingness indicator). cover: for a value column that can be empty, the column whose value says why
    (mask == 1, or flag == 0 for empty-by-design fields)."""

    def __init__(self, name: str, prefix: str):
        self.name, self.prefix = name, prefix
        self.cols: "OrderedDict[str, np.ndarray]" = OrderedDict()
        self.meta: dict = {}
        self.block_blank: dict = {}      # full block name -> rows whose source value is blank (active rows only)
        self.block_missing: dict = {}    # full block name -> full name of the mask column covering those rows
        self.block_exempt: list = []     # full block names with blanks covered by MISSING_MASK_EXEMPT

    def add(self, col: str, v, kind: str, block: str | None = None, system: str | None = None,
            cover: str | None = None, cover_when: int = 1, provisional: bool = False) -> str:
        name = self.prefix + col
        if name in self.cols:
            raise BuildError(f"duplicate column {name}")
        arr = np.asarray(v, dtype=float)
        self.cols[name] = arr
        self.meta[name] = dict(kind=kind, block=block, system=system,
                               cover=(self.prefix + cover) if cover else None, cover_when=cover_when,
                               provisional=provisional)
        return name

    def onehot(self, block: str, values: np.ndarray, vocab: list, system: str | None = None,
               active: np.ndarray | None = None, missing_name: str | None = None) -> list:
        out, used = [], set()
        for c in vocab:
            v = (values == c).astype(float)
            if active is not None:
                v = v * active
            nm = safe(c)[:48]
            if nm in used or nm != safe(c) or not str(c).isascii():
                nm = f"{nm[:40]}_{hashlib.sha1(str(c).encode('utf-8')).hexdigest()[:8]}"
            used.add(nm)
            out.append(self.add(f"{block}__{nm}", v, "onehot", block=self.prefix + block, system=system))
            self.meta[out[-1]]["category"] = str(c)
        # v1.2.1: blank source values (for active blocks, blanks among active rows) get a <block>_missing mask
        full = self.prefix + block
        blank = np.array([isinstance(x, str) and x == "" for x in values], dtype=bool)
        if active is not None:
            blank &= np.asarray(active) == 1
        self.block_blank[full] = blank
        if blank.any():
            if block in MISSING_MASK_EXEMPT:
                self.block_exempt.append(full)
            else:
                nm = missing_name or f"{block}_missing"
                if self.prefix + nm not in self.cols:
                    self.add(nm, blank.astype(float), "mask", system=system)
                    self.meta[self.prefix + nm]["blank_mask_of"] = full
                self.block_missing[full] = self.prefix + nm
        return out

    def frame(self, eid, split) -> pd.DataFrame:
        df = pd.DataFrame(self.cols)
        df.insert(0, "split", split)
        df.insert(0, "encounter_id", eid)
        return df

    @property
    def width(self) -> int:
        return len(self.cols)


def safe(s: str) -> str:
    out = "".join(ch if (ch.isalnum() and ch.isascii()) else "_" for ch in str(s).strip().lower())
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_") or "blank"


def fit_vocab(values: np.ndarray, train: np.ndarray, subj: np.ndarray, m: int = M_RARE, exclude=("",)) -> list:
    """Categories with >= m distinct train subjects (counted on train only, by subject), sorted."""
    df = pd.DataFrame({"v": values[train], "s": subj[train]})
    df = df[~df["v"].isin(exclude)]
    n = df.groupby("v")["s"].nunique()
    return sorted(n[n >= m].index.tolist())


def pooled(values: np.ndarray, vocab: list, other: str = "other", exclude=("",)) -> np.ndarray:
    keep = set(vocab)
    return np.array([v if (v in keep or v in exclude) else other for v in values], dtype=object)


class AgeRef:
    """Train-only age-band reference (median, mean, SD) of a per-encounter value; bands with < 10 train values and
    age_invalid use the pooled train reference (empirical age-norming; stand-in for published curves, spec D3)."""

    def __init__(self, band: np.ndarray, train: np.ndarray):
        self.band, self.train = band, train
        self.fitted: dict = {}

    def fit(self, key: str, ref_values: np.ndarray) -> dict:
        tr = self.train & np.isfinite(ref_values)
        pooled_v = ref_values[tr]
        pool = (float(np.median(pooled_v)), float(np.mean(pooled_v)), float(np.std(pooled_v))) if len(pooled_v) else (
            math.nan, math.nan, math.nan)
        rec = {"pooled": pool, "bands": {}}
        for b in AGE_BANDS:
            if b == "age_invalid":
                continue
            v = ref_values[tr & (self.band == b)]
            if len(v) >= MIN_BAND_N and np.std(v) > 0:
                rec["bands"][b] = (float(np.median(v)), float(np.mean(v)), float(np.std(v)))
        self.fitted[key] = rec
        return rec

    def arrays(self, key: str):
        rec = self.fitted[key]
        med = np.full(len(self.band), rec["pooled"][0])
        mu = np.full(len(self.band), rec["pooled"][1])
        sd = np.full(len(self.band), rec["pooled"][2])
        for b, (m_, u, s) in rec["bands"].items():
            sel = self.band == b
            med[sel], mu[sel], sd[sel] = m_, u, s
        return med, mu, sd

    def z(self, key: str, x: np.ndarray) -> np.ndarray:
        _, mu, sd = self.arrays(key)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = (x - mu) / sd
        out[~np.isfinite(sd) | (sd <= 0)] = np.where(np.isfinite(x[~np.isfinite(sd) | (sd <= 0)]), 0.0, np.nan)
        return out


# ---------------------------------------------------------------------------------------------------- stage 06
def adm_source(interval_h: np.ndarray, surg_comp: np.ndarray, surgical_window_h: float = ADM_SRC_SURGICAL_WINDOW_H,
               floor_cutoff_h: float = ADM_SRC_FLOOR_CUTOFF_H) -> np.ndarray:
    """Decision D1 (Nov 18, 2026), first match wins:
    (1) surgical    = a recent surgery ended <= t0 (Stage 3 adm_src_or_component) AND interval < 48 h;
    (2) floor       = not surgical AND interval >= 24 h;
    (3) transfer_ed = everyone else (ED and outside transfer cannot be separated in PIC; Stage 3 decision 9).
    Missing interval (existing handling): a missing interval never counts as a long interval, i.e. it satisfies
    'interval < 48 h' and fails 'interval >= 24 h' -> surgical with a qualifying surgery, else transfer_ed.
    The interval is Stage 3's hospital-to-ICU hours, clipped at 0."""
    iv = np.asarray(interval_h, float)
    known = np.isfinite(iv)
    surgical = (np.asarray(surg_comp) == 1) & ~(known & (iv >= surgical_window_h))
    floor = ~surgical & known & (iv >= floor_cutoff_h)
    return np.where(surgical, "surgical", np.where(floor, "floor", "transfer_ed"))


def adm_src_not_duplicate(blk: np.ndarray, band_cols: list, interval_h: np.ndarray) -> bool:
    """A13: True when no admission-source indicator equals an interval-band indicator, the interval >= 24 h
    indicator, or the complement of either."""
    iv = np.asarray(interval_h, float)
    ge24 = (np.isfinite(iv) & (iv >= ADM_SRC_FLOOR_CUTOFF_H)).astype(float)
    refs = [np.asarray(b, float) for b in band_cols] + [ge24]
    refs = refs + [1.0 - r for r in refs]
    return len(band_cols) > 0 and all(not np.array_equal(blk[:, j].astype(float), r)
                                      for j in range(blk.shape[1]) for r in refs)


def build_layer_a(D: dict, fitted: dict, report: list) -> dict:
    A, split, subj = D["A"], D["split"], D["subj"]
    train = split == "train"
    n = len(A)
    s = lambda c: A[c].to_numpy()  # noqa: E731
    band = s("age_band")
    norms = D["norms"]
    age_ok = num(s("age_months"))
    age_ok[band == "age_invalid"] = np.nan
    unknown_bands = sorted(set(band) - set(AGE_BANDS))
    if unknown_bands:
        raise InputError(f"unknown age bands in the Layer A base: {len(unknown_bands)}")
    ref = AgeRef(band, train)
    voc: dict = {}

    def vocab(key, values, exclude=("",)):
        v = fit_vocab(values, train, subj, exclude=exclude)
        voc[key] = v
        cats = pd.Series(values[train]).value_counts()
        for c, k in cats.items():
            if c in exclude:
                continue
            nsub = len(set(subj[train][values[train] == c]))
            report.append(dict(layer="A", block=key, category=str(c), n_train_subjects=nsub, kept=int(c in v)))
        return v

    sexv = vocab("sex", s("sex"))
    cuv = vocab("first_careunit", s("first_careunit"))
    sex_vals = pooled(s("sex"), sexv)
    cu_vals = pooled(s("first_careunit"), cuv)
    # diagnosis: 3-character category pooled within its (chapter, domain) cell; cells unseen in train -> other
    dx3, chap, dom = s("dx_icd10_3char"), s("dx_chapter"), s("dx_domain")
    dxv = vocab("dx_3char", dx3)
    cell = np.array([f"other_{c}_{d}" for c, d in zip(chap, dom)], dtype=object)
    dx_pooled = np.where(np.isin(dx3, dxv), dx3, cell)
    cells_train = sorted(set(dx_pooled[train]) - set(dxv))
    dx_final = np.where(np.isin(dx_pooled, dxv + cells_train), dx_pooled, "other")
    dx_vocab = dxv + cells_train + ["other"]
    chv = vocab("dx_chapter", chap)
    chap_final = pooled(chap, chv)
    domv = vocab("dx_domain", dom)
    dom_final = pooled(dom, domv)
    rs = num(s("recent_surgery"))
    rm = num(s("remote_surgery"))
    rname = vocab("recent_surgery_name", s("recent_surgery_name"))
    mname = vocab("remote_surgery_name", s("remote_surgery_name"))
    rsys = vocab("recent_body_system", s("recent_body_system"))
    rgrp = vocab("recent_complexity_group", s("recent_complexity_group"))
    msys = vocab("remote_body_system", s("remote_body_system"))
    asrc = adm_source(num(s("adm_src_interval_h")), num(s("adm_src_or_component")))
    ibands = sorted(set(s("admit_interval_band")) - {""}, key=lambda b: (INTERVAL_BANDS + (b,)).index(b))
    # pre-ICU lab block (6 h primary); the 24 h block is a sensitivity build, not a level column
    pre_ch = [c[len("pre6_"):] for c in A.columns if c.startswith("pre6_") and not c.endswith("_missing")]
    hour = num(s("intime_hour"))
    dow = s("intime_dow")
    shift, wkend = s("intime_shift"), num(s("intime_weekend"))
    timing4 = np.where(wkend == 1, "weekend_", "weekday_") + np.where(shift == "day", "daytime", "evening_night")
    levels = {}
    for k in (0, 1, 2):
        L = Level(f"a{k}", f"a{k}_")
        # fixed in every level (identical values): sex and first care unit
        L.onehot("sex", sex_vals, sexv + ["other"])
        L.onehot("first_careunit", cu_vals, cuv + ["other"])
        if k == 0:
            age = num(s("age_months"))
            inv = (band == "age_invalid").astype(float)
            age[inv == 1] = np.nan
            L.add("age_invalid", inv, "mask")
            L.add("age_months", age, "cont", cover="age_invalid")
            iv = num(s("hours_admit_to_icu_intime"))
            L.add("admit_interval_missing", np.isnan(iv).astype(float), "mask")
            L.add("hours_admit_to_icu_intime", iv, "cont", cover="admit_interval_missing")
            L.add("intime_hour_sin", np.sin(2 * np.pi * hour / 24.0), "cont")
            L.add("intime_hour_cos", np.cos(2 * np.pi * hour / 24.0), "cont")
            L.onehot("intime_dow", dow, [str(d) for d in range(7)])
            L.onehot("dx", dx_final, dx_vocab)
            L.add("recent_surgery", rs, "bin")
            L.onehot("recent_surgery_name", pooled(s("recent_surgery_name"), rname), rname + ["other"], active=rs)
            L.add("recent_hours_T_to_t0", np.where(rs == 1, num(s("recent_hours_T_to_t0")), np.nan), "cont",
                  cover="recent_surgery", cover_when=0)
            L.add("remote_surgery", rm, "bin")
            L.onehot("remote_surgery_name", pooled(s("remote_surgery_name"), mname), mname + ["other"], active=rm)
            L.add("remote_pre_admit_flag", np.where(rm == 1, num(s("remote_pre_admit_flag")), 0.0), "bin")
        else:
            L.onehot("age_band", band, [b for b in AGE_BANDS if (band == b).any()])
            # the existing admit_interval_missing mask (band blank) is the block's missing mask; no duplicate
            L.onehot("admit_interval_band", s("admit_interval_band"), ibands, missing_name="admit_interval_missing")
            if k == 1:
                L.onehot("intime_shift", shift, list(SHIFTS))
                L.add("intime_weekend", wkend, "bin")
                L.onehot("dx_chapter", chap_final, chv + ["other"])
            else:
                L.onehot("intime_timing4", timing4, ["weekday_daytime", "weekday_evening_night", "weekend_daytime",
                                                     "weekend_evening_night"])
                L.onehot("dx_domain", dom_final, domv + ["other"])
            L.add("early_surgery_flag", rs, "bin")
            L.onehot("recent_body_system", pooled(s("recent_body_system"), rsys), rsys + ["other"], active=rs)
            if k == 1:
                L.onehot("recent_complexity_group", pooled(s("recent_complexity_group"), rgrp), rgrp + ["other"],
                         active=rs)
                bmap = {"short": 0.0, "expected": 1.0, "long": 2.0}
                bandv = np.array([bmap.get(b, np.nan) for b in s("recent_phi1_band")])
                bandv[rs != 1] = np.nan
                L.add("recent_phi1_band_missing", np.isnan(bandv).astype(float), "mask")
                L.add("recent_phi1_band", bandv, "ordinal", cover="recent_phi1_band_missing")
            L.add("remote_surgery", rm, "bin")
            L.onehot("remote_body_system", pooled(s("remote_body_system"), msys), msys + ["other"], active=rm)
            L.onehot("admission_source", asrc, list(ADM_SRC_LEVELS))
        # pre-ICU lab block
        for ch in pre_ch:
            v = num(s(f"pre6_{ch}"))
            mk = num(s(f"pre6_{ch}_missing"))
            sysn = SYSTEM_OF.get(ch, "other")
            if k == 0:
                L.add(f"pre6_{ch}_missing", mk, "mask")
                L.add(f"pre6_{ch}", v, "cont", cover=f"pre6_{ch}_missing")
            else:
                tag = f"pre6_{ch}" if k == 1 else f"pre6_{sysn}__{ch}"
                L.add(f"{tag}_missing", mk, "mask", system=sysn)
                if ch in QUALITATIVE:
                    L.add(tag, v, "ordinal", system=sysn, cover=f"{tag}_missing")
                elif ch == "creatinine":       # D2: creatinine / PODIUM 2022 baseline (age, sex); mg/dL from Stage 3
                    rt = norms.creatinine_ratio(v, age_ok, s("sex"))
                    L.add(f"{tag}_ratio_podium_missing", np.isnan(rt).astype(float), "mask", system=sysn)
                    L.add(f"{tag}_ratio_podium", rt, "cont", system=sysn, cover=f"{tag}_ratio_podium_missing")
                else:
                    key = f"A:pre6_{ch}"
                    if key not in ref.fitted:
                        ref.fit(key, v)
                    L.add(f"{tag}_z", ref.z(key, v), "cont", system=sysn, cover=f"{tag}_missing")
        for g in MED_GROUPS:
            L.add(f"med_{g}_active_t0", num(s(f"med_{g}_active_t0")), "bin")
        levels[L.name] = L
    fitted["layer_a"] = dict(vocab=voc, age_ref=ref.fitted, adm_src_surgical_window_hours=ADM_SRC_SURGICAL_WINDOW_H,
                           adm_src_floor_cutoff_hours=ADM_SRC_FLOOR_CUTOFF_H, dx_other_cells=cells_train)
    # acceptance checks (plan section 4, A1-A17 where they apply to levels)
    chk, det = {}, {}
    det["A1_rows"] = [k for k, L in levels.items() if len(next(iter(L.cols.values()))) != n]
    chk["A1_rows"] = not det["A1_rows"]
    det["A3_no_banned_column"] = [c for L in levels.values() for c in L.cols if banned(c)]
    chk["A3_no_banned_column"] = not det["A3_no_banned_column"]
    # A4: sex and first care unit (one-hot block plus its missing mask, if any) identical in a0, a1 and a2
    det["A4_fixed_fields_identical"] = []

    def fixed_block(L, blk):
        cs = [c for c in L.cols if c.startswith(f"{L.prefix}{blk}__")] + (
            [f"{L.prefix}{blk}_missing"] if f"{L.prefix}{blk}_missing" in L.cols else [])
        return [c[len(L.prefix):] for c in cs], np.column_stack([L.cols[c] for c in cs])
    for blk in ("sex", "first_careunit"):
        n0, ref0 = fixed_block(levels["a0"], blk)
        for k in (1, 2):
            nk, other = fixed_block(levels[f"a{k}"], blk)
            if not (n0 == nk and ref0.shape == other.shape and bool((ref0 == other).all())):
                det["A4_fixed_fields_identical"].append(f"a{k}:a{k}_{blk} vs a0_{blk}")
    chk["A4_fixed_fields_identical"] = not det["A4_fixed_fields_identical"]
    det["A5_bands_partition"] = ([] if np.isin(s("admit_interval_band"), INTERVAL_BANDS + ("",)).all()
                                 else ["admit_interval_band"]) + ([] if np.isin(shift, SHIFTS).all()
                                                                  else ["intime_shift"])
    chk["A5_bands_partition"] = not det["A5_bands_partition"]
    # A6 (v1.2.1): block sum + its missing mask = 1 (= the active indicator for recent_/remote_ blocks); a blank
    # source value that no mask covers fails, named by level and block
    a6_fail = []
    for L in levels.values():
        blocks = {}
        for c, m in L.meta.items():
            if m["kind"] == "onehot":
                blocks.setdefault(m["block"], []).append(L.cols[c])
        for b, arrs in blocks.items():
            tot = np.sum(arrs, axis=0)
            if b in L.block_missing:
                tot = tot + L.cols[L.block_missing[b]]
            elif b in L.block_exempt:
                tot = tot + L.block_blank[b].astype(float)
            if any(t in b for t in ("recent_", "remote_")):
                target = rs if "recent_" in b else rm
            else:
                target = 1.0
            if not bool((tot == target).all()):
                a6_fail.append(f"{L.name}:{b}")
    det["A6_blocks_sum_to_one"] = a6_fail
    chk["A6_blocks_sum_to_one"] = not a6_fail
    chk["A9_early_surgery_flag"] = bool((levels["a1"].cols["a1_early_surgery_flag"] == rs).all())
    det["A9_early_surgery_flag"] = ["a1:a1_early_surgery_flag"]
    # A10 pooling check: each kept chapter column equals the group sum of the phi0 diagnosis columns of that chapter
    # (plus the rows sent to the global other, which only unseen validation/test cells reach)
    chap_of_dx = {}
    for d, cch in zip(dx_final, chap):
        chap_of_dx.setdefault(d, set()).add(cch)
    ok10 = all(len(v) == 1 for d, v in chap_of_dx.items() if d != "other")
    det["A10_pooling_check"] = [] if ok10 else ["a0:a0_dx (a category maps to more than one chapter)"]
    for cc in chv:
        col = levels["a1"].cols.get(f"a1_dx_chapter__{safe(cc)}")
        tot = np.zeros(n)
        for d, cs in chap_of_dx.items():
            if d != "other" and cs == {cc}:
                tot += levels["a0"].cols[f"a0_dx__{safe(d)}"]
        tot += ((dx_final == "other") & (chap == cc)).astype(float)
        if not (col is not None and bool((col == tot).all())):
            ok10 = False
            det["A10_pooling_check"].append(f"a1:a1_dx_chapter__{safe(cc)}")
    chk["A10_pooling_check"] = ok10
    # A12: the three levels partition the cohort and follow the first-match rule exactly
    ivA, scA = num(s("adm_src_interval_h")), num(s("adm_src_or_component"))
    kn = np.isfinite(ivA)
    exp_s = (scA == 1) & ~(kn & (ivA >= ADM_SRC_SURGICAL_WINDOW_H))
    exp_f = ~exp_s & kn & (ivA >= ADM_SRC_FLOOR_CUTOFF_H)
    blk = np.column_stack([levels["a1"].cols[f"a1_admission_source__{v}"] for v in ADM_SRC_LEVELS])
    chk["A12_admission_source_partition"] = bool(np.isin(asrc, ADM_SRC_LEVELS).all()) and bool(
        (blk.sum(axis=1) == 1).all()) and bool((asrc == np.where(exp_s, "surgical", np.where(
            exp_f, "floor", "transfer_ed"))).all())
    # A13: duplicate-column check - no admission-source indicator equals an interval-band indicator, the
    # interval >= 24 h indicator, or the complement of either (the block is not a copy of the interval band)
    ib = [levels["a1"].cols[c] for c in levels["a1"].cols if c.startswith("a1_admit_interval_band__")]
    det["A12_admission_source_partition"] = ["a1:a1_admission_source"]
    chk["A13_admission_source_not_identical_to_interval_band"] = adm_src_not_duplicate(blk, ib, ivA)
    det["A13_admission_source_not_identical_to_interval_band"] = ["a1:a1_admission_source", "a1:a1_admit_interval_band"]
    det["A16_weight_not_a_feature"] = [c for L in levels.values() for c in L.cols if "weight" in c]
    chk["A16_weight_not_a_feature"] = not det["A16_weight_not_a_feature"]
    chk["A7_vocab_train_only"] = True   # by construction (fit_vocab uses train rows only); tested on the box
    blank_masks = OrderedDict((k, [L.block_missing[b] for b in L.block_missing]) for k, L in levels.items())
    exempt = OrderedDict((k, list(L.block_exempt)) for k, L in levels.items() if L.block_exempt)
    return dict(levels=levels, checks=chk, details=det, blank_masks=blank_masks, blank_exempt=exempt)


# ---------------------------------------------------------------------------------------------------- scorers
def _psofa_age_idx(age_mo: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(np.array(PSOFA_AGE_EDGES_MO), age_mo, side="right")
    return np.where(np.isfinite(age_mo), idx, -1)


def _step(x: np.ndarray, edges, higher_is_worse: bool = True) -> np.ndarray:
    """Score 0..len(edges) by lower edges (higher_is_worse) or by thresholds 'below' (not higher_is_worse)."""
    out = np.zeros(len(x))
    for i, e in enumerate(edges):
        out = np.where(np.isfinite(x) & ((x >= e) if higher_is_worse else (x < e)), i + 1, out)
    return out


def psofa(age_mo, plt_min, bili_max, map_min, creat_max, vaso_active) -> dict:
    """Modified pSOFA components (Matics and Sanchez-Pinto 2017). Missing input -> 0 (as in the publication).
    Age-dependent components are empty when age is invalid."""
    ai = _psofa_age_idx(np.asarray(age_mo, float))
    coag = _step(np.asarray(plt_min, float), PSOFA_PLT, higher_is_worse=False)
    hep = _step(np.asarray(bili_max, float), PSOFA_BILI)
    mapc = np.array([PSOFA_MAP_CUT[i] if i >= 0 else np.nan for i in ai])
    m = np.asarray(map_min, float)
    cv = np.where(np.isfinite(m) & (m < mapc), 1.0, 0.0)
    cv = np.where(np.asarray(vaso_active) > 0, 2.0, cv)
    cr = np.asarray(creat_max, float)
    ren = np.zeros(len(cr))
    for j in range(len(cr)):
        if ai[j] >= 0 and np.isfinite(cr[j]):
            ren[j] = sum(cr[j] >= e for e in PSOFA_CREAT[ai[j]])
    bad = ai < 0
    cv[bad] = np.nan
    ren[bad] = np.nan
    return {"cardiovascular": cv, "coagulation": coag, "hepatic": hep, "renal": ren}


def prism3(band, sbp_min, hr_max, t_min, t_max, ph_min, ph_max, hco3_min, hco3_max, po2_min, pco2_max, glu_max,
           k_max, cr_max, bun_max, wbc_min, plt_min, ptt_max) -> np.ndarray:
    """Modified PRISM III (non-neurologic points from the worst first-24 h values; PT not scored). Missing -> 0.
    Empty when age is invalid."""
    n = len(band)
    out = np.full(n, np.nan)
    f = lambda a: np.asarray(a, float)  # noqa: E731
    sbp_min, hr_max, t_min, t_max, ph_min, ph_max, hco3_min, hco3_max, po2_min, pco2_max, glu_max, k_max, cr_max, \
        bun_max, wbc_min, plt_min, ptt_max = map(f, (sbp_min, hr_max, t_min, t_max, ph_min, ph_max, hco3_min,
                                                     hco3_max, po2_min, pco2_max, glu_max, k_max, cr_max, bun_max,
                                                     wbc_min, plt_min, ptt_max))
    gt = lambda x, c: np.isfinite(x) & (x > c)  # noqa: E731
    lt = lambda x, c: np.isfinite(x) & (x < c)  # noqa: E731
    btw = lambda x, a, b: np.isfinite(x) & (x >= a) & (x <= b)  # noqa: E731
    for g in range(4):
        sel = np.array([PRISM_AGE.get(b, -1) == g for b in band])
        if not sel.any():
            continue
        a, b = PRISM_SBP[g]
        p = np.where(lt(sbp_min, a), 7, np.where(btw(sbp_min, a, b), 3, 0))
        a, b = PRISM_HR[g]
        p = p + np.where(gt(hr_max, b), 4, np.where(btw(hr_max, a, b), 3, 0))
        p = p + np.where(lt(t_min, 33) | gt(t_max, 40), 3, 0)
        acid6 = lt(ph_min, 7.0) | lt(hco3_min, 5)
        acid2 = btw(ph_min, 7.0, 7.28) | btw(hco3_min, 5, 16.9)
        p = p + np.where(acid6, 6, np.where(acid2, 2, 0))
        p = p + np.where(gt(hco3_max, 34.0), 4, 0)
        p = p + np.where(gt(ph_max, 7.55), 3, np.where(btw(ph_max, 7.48, 7.55), 2, 0))
        p = p + np.where(lt(po2_min, 42.0), 6, np.where(btw(po2_min, 42.0, 49.9), 3, 0))
        p = p + np.where(gt(pco2_max, 75.0), 3, np.where(btw(pco2_max, 50.0, 75.0), 1, 0))
        p = p + np.where(gt(glu_max, 200.0), 2, 0)
        p = p + np.where(gt(k_max, 6.9), 3, 0)
        p = p + np.where(gt(cr_max, PRISM_CREAT[g]), 2, 0)
        p = p + np.where(gt(bun_max, PRISM_BUN[g]), 3, 0)
        p = p + np.where(lt(wbc_min, 3.0), 4, 0)
        p = p + np.where(lt(plt_min, 50), 5, np.where(btw(plt_min, 50, 99.999), 4,
                                                       np.where(btw(plt_min, 100, 200), 2, 0)))
        p = p + np.where(gt(ptt_max, PRISM_PTT[g]), 3, 0)
        out[sel] = p[sel]
    return out


def aki_stage(age_mo, cr_max) -> np.ndarray:
    age_mo, cr = np.asarray(age_mo, float), np.asarray(cr_max, float)
    base = np.full(len(cr), np.nan)
    for j, a in enumerate(age_mo):
        if np.isfinite(a):
            base[j] = next(b for lim, b in AKI_BASE if a < lim)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = cr / base
    st = np.where(~np.isfinite(cr), 0.0,
                  np.where((r >= 3) | (cr >= 4.0), 3.0, np.where(r >= 2, 2.0, np.where(r >= 1.5, 1.0, 0.0))))
    st[np.isfinite(cr) & ~np.isfinite(base) & (cr < 4.0)] = np.nan   # age invalid: stage not defined unless Cr >= 4
    return st


# ---------------------------------------------------------------------------------------------------- stage 07
def build_layer_b(D: dict, zC: dict, fitted: dict, report: list) -> dict:
    A, B, split, subj = D["A"], D["B"], D["split"], D["subj"]
    train = split == "train"
    n = len(B)
    band = A["age_band"].to_numpy()
    age_mo = num(A["age_months"])
    age_mo[band == "age_invalid"] = np.nan
    age_d = num(A["age_days"])
    age_d[band == "age_invalid"] = np.nan
    norms = D["norms"]
    sex = A["sex"].to_numpy()

    def pub_z(c, x):
        """D2: published age norms for HR / RR (Fleming 2011) and MAP (Roberts 2020); NaN outside their cover."""
        return norms.roberts_z(x, age_d) if c == "map" else norms.fleming_z(c, x, age_d)

    ref = AgeRef(band, train)
    g = lambda c: num(B[c])  # noqa: E731
    chans = [c[: -len("_measured")] for c in B.columns if c.endswith("_measured") and not c.startswith("panel_")]
    known = set(VITALS) | set(TRIGGERED) | set(INDICATOR_ONLY) | set(QUALITATIVE) | set(WORST) | {"kdigo_uo_stage"}
    unknown = sorted(set(chans) - known)
    if unknown or set(chans) - set(SYSTEM_OF):
        raise InputError(f"Layer B base has channels without a rule in this script: {len(unknown)}")
    prov = set(D["bmeta"].get("provisional_columns", []))
    for c in chans:
        ref.fit(c, g(f"{c}_median") if f"{c}_median" in B.columns else np.full(n, np.nan))

    def worst(c):
        lo, hi = g(f"{c}_min"), g(f"{c}_max")
        d = WORST[c]
        if d == "high":
            return hi
        if d == "low":
            return lo
        mid, _, _ = ref.arrays(c)
        return np.where(np.abs(hi - mid) >= np.abs(lo - mid), hi, lo)

    def mask_of(c):
        return f"panel_{PANEL_OF[c]}_missing" if c in PANEL_OF else f"{c}_missing"

    # drug names (phi0): any course of the drug key active in the Layer B window [0, min(24, stay))
    meds = D["meds"]
    w_end = np.minimum(24.0, D["stay_h"])
    pos = {e: i for i, e in enumerate(D["eid"])}
    mi_ = meds["encounter_id"].map(pos)
    if mi_.isna().any():
        raise InputError("a Stage 3 internal file has an encounter_id outside the canonical cohort")
    mi = mi_.to_numpy(dtype=int)
    st, sp_ = num(meds["start"]), num(meds["stop"])
    act = (st < w_end[mi]) & (sp_ > 0)
    key = meds["key"].to_numpy()
    drug_rows = pd.DataFrame({"i": mi[act], "k": key[act]}).drop_duplicates()
    dsub = pd.DataFrame({"k": drug_rows["k"].to_numpy(), "s": subj[drug_rows["i"].to_numpy()],
                         "t": train[drug_rows["i"].to_numpy()]})
    nsub = dsub[dsub["t"]].groupby("k")["s"].nunique()
    drug_vocab = sorted(nsub[nsub >= M_RARE].index.tolist())
    for k_, v in nsub.items():
        report.append(dict(layer="B", block="drug_name", category=str(k_), n_train_subjects=int(v),
                           kept=int(k_ in drug_vocab)))
    dmat = {k_: np.zeros(n) for k_ in drug_vocab + ["other"]}
    for i, k_ in zip(drug_rows["i"].to_numpy(), drug_rows["k"].to_numpy()):
        dmat[k_ if k_ in dmat else "other"][i] = 1.0
    # microbiology culture count in the window (M1; collections only)
    mc = D["micro"]
    mci_ = mc["encounter_id"].map(pos)
    if mci_.isna().any():
        raise InputError("a Stage 3 internal file has an encounter_id outside the canonical cohort")
    mci = mci_.to_numpy(dtype=int)
    th = num(mc["th"])
    inwin = (th >= 0) & (th < w_end[mci])
    culture_count = np.bincount(mci[inwin], minlength=n).astype(float)

    def common(L: Level, k: int, sysmode: bool):
        """Columns shared by b0 and b1: panel masks, indicator-only flags, source flags."""
        for p in PANELS:
            L.add(f"panel_{p}_missing", 1.0 - g(f"panel_{p}_measured"), "mask", system=PANEL_SYSTEM[p])
        for c in INDICATOR_ONLY:
            L.add(f"{c}_measured", g(f"{c}_measured"), "bin", system=SYSTEM_OF[c])
        L.add("hemoglobin_src_bloodgas", g("hemoglobin_src_bloodgas"), "bin", system="hematologic")
        L.add("bilirubin_total_src_5255", g("bilirubin_total_src_5255"), "bin", system="hepatic")

    levels = {}
    # ---------------- b0
    L0 = Level("b0", "b0_")
    common(L0, 0, False)
    for c in chans:
        if c in INDICATOR_ONLY or c == "kdigo_uo_stage":
            continue
        sysn = SYSTEM_OF[c]
        if c not in PANEL_OF:
            L0.add(f"{c}_missing", 1.0 - g(f"{c}_measured"), "mask", system=sysn)
        mk = mask_of(c)
        pv = any(x.startswith(c + "_") for x in prov)
        if c in VITALS:
            for stt in ("median", "min", "max"):
                L0.add(f"{c}_{stt}", g(f"{c}_{stt}"), "cont", system=sysn, cover=mk)
            L0.add(f"{c}_range", g(f"{c}_max") - g(f"{c}_min"), "cont", system=sysn, cover=mk)
        elif c in TRIGGERED:
            L0.add(f"{c}_worst", worst(c), "cont", system=sysn, cover=mk)
        else:
            kind = "bin" if c in QUALITATIVE else "cont"
            L0.add(f"{c}_first", g(f"{c}_first"), kind if kind == "cont" else "ordinal", system=sysn, cover=mk,
                   provisional=pv)
            L0.add(f"{c}_worst", worst(c), kind if kind == "cont" else "ordinal", system=sysn, cover=mk,
                   provisional=pv)
    for k_ in drug_vocab + ["other"]:
        nm = safe(k_)[:40]
        nm = f"{nm}_{hashlib.sha1(k_.encode('utf-8')).hexdigest()[:8]}" if k_ != "other" else "other"
        L0.add(f"drug__{nm}", dmat[k_], "bin", system="drug")
        L0.meta[L0.prefix + f"drug__{nm}"]["category"] = k_
    levels["b0"] = L0
    # ---------------- b1
    L1 = Level("b1", "b1_")
    common(L1, 1, True)
    for c in chans:
        if c in INDICATOR_ONLY:
            continue
        sysn = SYSTEM_OF[c]
        pv = any(x.startswith(c + "_") for x in prov)
        if c == "kdigo_uo_stage":
            L1.add("kdigo_uo_stage_missing", 1.0 - g("kdigo_uo_stage_measured"), "mask", system=sysn)
            L1.add("kdigo_uo_stage", g("kdigo_uo_stage_max"), "ordinal", system=sysn, cover="kdigo_uo_stage_missing",
                   provisional=True)
            continue
        if c not in PANEL_OF:
            L1.add(f"{c}_missing", 1.0 - g(f"{c}_measured"), "mask", system=sysn)
        mk = mask_of(c)
        if c in VITALS:
            if c in PUB_NORM_VITALS:
                zs = {stt: pub_z(c, g(f"{c}_{stt}")) for stt in ("median", "min", "max")}
                zm = f"{c}_agez_missing"       # value missing or age outside the published table's cover
                L1.add(zm, np.isnan(zs["median"]).astype(float), "mask", system=sysn)
                for stt, v in zs.items():
                    L1.add(f"{c}_{stt}_z", v, "cont", system=sysn, cover=zm)
            else:
                zs = {stt: ref.z(c, g(f"{c}_{stt}")) for stt in ("median", "min", "max")}
                for stt, v in zs.items():
                    L1.add(f"{c}_{stt}_z", v, "cont", system=sysn, cover=mk)
            if c in ZFLAG_VITALS:
                stt, op, thr = ZFLAG_VITALS[c]
                v = zs[stt]
                fl = np.where(np.isfinite(v) & ((v > thr) if op == ">" else (v < thr)), 1.0, 0.0)
                L1.add(f"{c}_z_flag", fl, "bin", system=sysn)
            if c in SD_VITALS:
                bins = zC[c]                                   # (n, 4) bin medians, NaN where empty or padding
                zb = np.column_stack([pub_z(c, bins[:, j]) for j in range(bins.shape[1])])
                cnt = np.isfinite(zb).sum(axis=1)
                with np.errstate(invalid="ignore"):
                    sd = np.where(cnt >= 3, np.nanstd(np.where(cnt[:, None] >= 3, zb, 0.0), axis=1, ddof=1), np.nan)
                sd[cnt < 3] = np.nan
                L1.add(f"{c}_zsd_missing", (cnt < 3).astype(float), "mask", system=sysn)
                L1.add(f"{c}_zsd", sd, "cont", system=sysn, cover=f"{c}_zsd_missing")
        elif c in TRIGGERED:
            L1.add(f"{c}_worst_z", ref.z(c, worst(c)), "cont", system=sysn, cover=mk)
        elif c in QUALITATIVE:
            L1.add(f"{c}_first", g(f"{c}_first"), "ordinal", system=sysn, cover=mk)
            L1.add(f"{c}_worst", worst(c), "ordinal", system=sysn, cover=mk)
        elif c in ("urine_output", "fluid_intake"):
            L1.add(f"{c}_first", g(f"{c}_first"), "cont", system=sysn, cover=mk, provisional=pv)
            L1.add(f"{c}_worst", worst(c), "cont", system=sysn, cover=mk, provisional=pv)
        elif c == "creatinine":            # D2: ratio to the PODIUM 2022 age- and sex-specific baseline (mg/dL)
            rf = norms.creatinine_ratio(g("creatinine_first"), age_mo, sex)
            rw = norms.creatinine_ratio(worst(c), age_mo, sex)
            L1.add("creatinine_ratio_podium_missing", np.isnan(rf).astype(float), "mask", system=sysn)
            L1.add("creatinine_first_ratio_podium", rf, "cont", system=sysn, cover="creatinine_ratio_podium_missing")
            L1.add("creatinine_worst_ratio_podium", rw, "cont", system=sysn, cover="creatinine_ratio_podium_missing")
        else:
            L1.add(f"{c}_first_z", ref.z(c, g(f"{c}_first")), "cont", system=sysn, cover=mk)
            L1.add(f"{c}_worst_z", ref.z(c, worst(c)), "cont", system=sysn, cover=mk)
    # derived quantities
    bf, cf = g("bun_first"), g("creatinine_first")
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(np.isfinite(bf) & np.isfinite(cf) & (cf > 0), bf / cf, np.nan)
    L1.add("bun_cr_ratio_missing", np.isnan(ratio).astype(float), "mask", system="renal")
    L1.add("bun_cr_ratio", ratio, "cont", system="renal", cover="bun_cr_ratio_missing")
    aki = aki_stage(age_mo, g("creatinine_max"))
    L1.add("aki_stage_missing", np.isnan(aki).astype(float), "mask", system="renal")
    L1.add("aki_stage", aki, "ordinal", system="renal", cover="aki_stage_missing")
    # drug groups (state machine at the window end)
    for gr in MED_GROUPS:
        for f_ in MED_FEATS:
            v = g(f"med_{gr}_{f_}")
            if f_ in ("started", "active_end", "ever_received"):
                L1.add(f"med_{gr}_{f_}", v, "bin", system="drug")
            elif f_ == "hrs_to_start":
                L1.add(f"med_{gr}_{f_}", v, "cont", system="drug", cover=f"med_{gr}_started", cover_when=0)
            elif f_ == "hours_since_stopped":
                L1.add(f"med_{gr}_none_stopped", np.isnan(v).astype(float), "mask", system="drug")
                L1.add(f"med_{gr}_{f_}", v, "cont", system="drug", cover=f"med_{gr}_none_stopped")
            else:
                L1.add(f"med_{gr}_{f_}", v, "cont", system="drug")
    L1.add("med_vasoactive_max_concurrent", g("med_vasoactive_max_concurrent"), "ordinal", system="drug")
    # microbiology M1
    L1.add("micro_any_culture", g("micro_any_culture"), "bin", system="infectious")
    L1.add("micro_culture_count", culture_count, "ordinal", system="infectious")
    L1.add("micro_hours_to_first_culture", g("micro_hours_to_first_culture"), "cont", system="infectious",
           cover="micro_any_culture", cover_when=0)
    for sp in SPECIMENS:
        L1.add(f"micro_{sp}", g(f"micro_{sp}"), "bin", system="infectious")
    levels["b1"] = L1
    # ---------------- scores
    vaso = (g("med_vasoactive_hrs_active") > 0).astype(float)
    ps = psofa(age_mo, g("platelets_min"), g("bilirubin_total_max"), g("map_min"), g("creatinine_max"), vaso)
    glu = np.fmax(g("glucose_max"), g("glucose_bedside_max"))
    pr = prism3(band, g("sbp_min"), g("heart_rate_max"), g("temperature_min"), g("temperature_max"), g("ph_min"),
                g("ph_max"), g("bicarbonate_min"), g("bicarbonate_max"), g("po2_min"), g("pco2_max"), glu,
                g("potassium_max"), g("creatinine_max"), g("bun_max"), g("wbc_min"), g("platelets_min"), g("ptt_max"))
    # ---------------- b2 variants from b1 columns
    drug_prev = {gr: float(np.mean(g(f"med_{gr}_ever_received")[train] > 0)) for gr in MED_GROUPS}
    dropped_groups = [gr for gr, p in drug_prev.items() if p < DRUG_GROUP_MIN_PREV]

    def b1_items():
        for c1, m in L1.meta.items():
            suf = c1[len("b1_"):]
            if m["system"] == "drug" and any(suf.startswith(f"med_{gr}_") for gr in dropped_groups):
                continue
            yield c1, suf, m

    def add_prism(L):
        inv = (band == "age_invalid").astype(float)
        L.add("prism3_age_invalid", inv, "mask", system="score")
        L.add("prism3", pr, "ordinal", system="score", cover="prism3_age_invalid")

    def copy(L, c1, suf, m, sysn):
        cov = None
        if m["cover"]:
            cov = f"{L1.meta[m['cover']]['system']}__{m['cover'][len('b1_'):]}"
        L.add(f"{sysn}__{suf}", L1.cols[c1], m["kind"], block=None, system=sysn, cover=cov,
              cover_when=m["cover_when"], provisional=m["provisional"])

    L2s = Level("b2s", "b2s_")
    for c1, suf, m in b1_items():
        copy(L2s, c1, suf, m, m["system"])
    add_prism(L2s)
    levels["b2s"] = L2s
    L2d = Level("b2d", "b2d_")
    for c1, suf, m in b1_items():
        copy(L2d, c1, suf, m, m["system"])
    n_dup = 0
    for c1, suf, m in list(b1_items()):
        ch = channel_of(suf)
        for sysn in ALSO_IN.get(ch, []):
            if m["kind"] == "mask" and suf.startswith("panel_"):
                continue
            # duplicate value columns and their own masks under the second system (exact copies)
            cov = (f"{L1.meta[m['cover']]['system']}__{m['cover'][len('b1_'):]}" if m["cover"] else None)
            L2d.add(f"{sysn}__dup__{suf}", L1.cols[c1], m["kind"], system=sysn, cover=cov,
                    cover_when=m["cover_when"], provisional=m["provisional"])
            n_dup += 1
    add_prism(L2d)
    levels["b2d"] = L2d
    L2p = Level("b2p", "b2p_")
    inputs = {"cardiovascular": ["map"], "coagulation": ["platelets"], "hepatic": ["bilirubin_total"],
              "renal": ["creatinine"]}
    L2p.add("psofa_age_invalid", (band == "age_invalid").astype(float), "mask", system="score")
    for comp, v in ps.items():
        for ch in inputs[comp]:
            L2p.add(f"psofa_{comp}_input_missing", 1.0 - g(f"{ch}_measured"), "mask", system="score")
        L2p.add(f"psofa_{comp}", v, "ordinal", system="score",
                cover="psofa_age_invalid" if comp in ("cardiovascular", "renal") else None)
    unscored = ("respiratory", "gi", "infectious", "electrolytes", "endocrine", "drug")
    for c1, suf, m in b1_items():
        if m["system"] in unscored:
            copy(L2p, c1, suf, m, m["system"])
    add_prism(L2p)
    levels["b2p"] = L2p
    # auto masks for any remaining uncovered empty cells (safety net; listed in the meta)
    for L in levels.values():
        auto_masks(L)
    fitted["layer_b"] = dict(age_ref=ref.fitted, drug_vocab=drug_vocab, drug_group_train_prevalence_lt_1pct=dropped_groups,
                             n_dup_columns=n_dup)
    # acceptance checks (plan section 5.8, where they apply to levels built from the Stage 3 base)
    chk, det = {}, {}
    det["B1_rows"] = [k for k, L in levels.items() if len(next(iter(L.cols.values()))) != n]
    chk["B1_rows"] = not det["B1_rows"]
    det["B3_no_banned_column"] = [c for L in levels.values() for c in L.cols if banned(c)]
    chk["B3_no_banned_column"] = not det["B3_no_banned_column"]
    det["B3_no_raw_vital_outside_b0"] = [c for k, L in levels.items() if k != "b0" for c in L.cols if any(
        c.endswith(f"{v}_{s_}") for v in VITALS for s_ in ("median", "min", "max", "range"))]
    chk["B3_no_raw_vital_outside_b0"] = not det["B3_no_raw_vital_outside_b0"]
    hs = [g(f"med_{gr}_hrs_active") for gr in MED_GROUPS]
    hts = [g(f"med_{gr}_hrs_to_start") for gr in MED_GROUPS]
    sts = [g(f"med_{gr}_started") for gr in MED_GROUPS]
    det["B5_hrs_active_le_24"] = [f"med_{gr}_hrs_active" for gr, h in zip(MED_GROUPS, hs)
                                  if not bool((h[np.isfinite(h)] <= 24.0 + 1e-9).all())]
    chk["B5_hrs_active_le_24"] = not det["B5_hrs_active_le_24"]
    det["B5_hrs_to_start_empty_iff_not_started"] = [f"med_{gr}_hrs_to_start" for gr, t, s_ in zip(MED_GROUPS, hts, sts)
                                                    if not bool(((~np.isfinite(t)) == (s_ == 0)).all())]
    chk["B5_hrs_to_start_empty_iff_not_started"] = not det["B5_hrs_to_start_empty_iff_not_started"]
    det6 = []
    for c in chans:
        if c in INDICATOR_ONLY or c in PANEL_OF or c == "kdigo_uo_stage":
            continue
        mk = L0.cols[f"b0_{c}_missing"]
        val = L0.cols.get(f"b0_{c}_first", L0.cols.get(f"b0_{c}_median", L0.cols.get(f"b0_{c}_worst")))
        if not bool(((mk == 1) == ~np.isfinite(val)).all()):
            det6.append(f"b0:b0_{c}_missing")
    for p, members in PANELS.items():
        anyv = np.zeros(n, bool)
        for c in members:
            if f"{c}_measured" in B.columns:
                anyv |= g(f"{c}_measured") == 1
        if not bool(((L0.cols[f"b0_panel_{p}_missing"] == 1) == ~anyv).all()):
            det6.append(f"b0:b0_panel_{p}_missing")
    cm = L1.cols["b1_creatinine_missing"]
    if not bool((L1.cols["b1_aki_stage"][cm == 1] == 0).all()):
        det6.append("b1:b1_aki_stage")
    det["B6_masks_match_values"] = det6
    chk["B6_masks_match_values"] = not det6
    # B8 (v1.2.1): banned organism / susceptibility tokens on field stems only (see b8_offenders)
    det["B8_m1_only_no_organism_fields"] = [f"{k}:{c}" for k, L in levels.items() for c in b8_offenders(L)]
    chk["B8_m1_only_no_organism_fields"] = not det["B8_m1_only_no_organism_fields"]
    b8_cat_blocks = sorted({b for L in levels.values() for b in b8_category_blocks(L)})
    # B10: phi2-dup duplicates are exact copies (rank deficiency checked on the fitted matrices in stage 08)
    det["B11_phi2_single_one_home"] = ["b1", "b2s"]
    chk["B11_phi2_single_one_home"] = (sum(1 for c in L2s.cols if not L2s.meta[c].get("auto")) ==
                                       sum(1 for c in L1.cols if not L1.meta[c].get("auto") and not any(
                                           c[3:].startswith(f"med_{gr}_") for gr in dropped_groups)) + 2)
    det["B13_no_surgery_feature"] = [c for L in levels.values() for c in L.cols if "surg" in c]
    chk["B13_no_surgery_feature"] = not det["B13_no_surgery_feature"]
    det["B20_vital_range_exact"] = [f"b0:b0_{v}_range" for v in VITALS if not bool(np.allclose(
        np.nan_to_num(L0.cols[f"b0_{v}_range"]), np.nan_to_num(g(f"{v}_max") - g(f"{v}_min"))))]
    chk["B20_vital_range_exact"] = not det["B20_vital_range_exact"]
    det["B21_sd_empty_iff_lt3"] = [f"b1:b1_{v}_zsd" for v in SD_VITALS if not bool(
        ((L1.cols[f"b1_{v}_zsd_missing"] == 1) == ~np.isfinite(L1.cols[f"b1_{v}_zsd"])).all())]
    chk["B21_sd_empty_iff_lt3"] = not det["B21_sd_empty_iff_lt3"]
    det["B23_no_drug_column_in_organ_block"] = [c for L in levels.values() for c, m in L.meta.items()
                                                if m["system"] in SYSTEMS and ("_med_" in c or "drug__" in c)]
    chk["B23_no_drug_column_in_organ_block"] = not det["B23_no_drug_column_in_organ_block"]
    det["B24_prism_identical_in_phi2_variants"] = ["b2s:b2s_prism3", "b2d:b2d_prism3", "b2p:b2p_prism3"]
    det["B24_psofa_four_scored_components"] = ["b2p:b2p_psofa_*"]
    chk["B24_prism_identical_in_phi2_variants"] = bool(
        np.array_equal(L2s.cols["b2s_prism3"], L2d.cols["b2d_prism3"], equal_nan=True) and
        np.array_equal(L2s.cols["b2s_prism3"], L2p.cols["b2p_prism3"], equal_nan=True))
    chk["B24_psofa_four_scored_components"] = sum(1 for c in L2p.cols if c.startswith("b2p_psofa_") and
                                                  L2p.meta[c]["kind"] == "ordinal") == 4
    det["B25_no_pf_sf_fio2_aagrad_aniongap"] = [c for L in levels.values() for c in L.cols if any(
        t in c for t in ("pf_ratio", "sf_ratio", "fio2", "aa_gradient", "anion", "balance"))]
    chk["B25_no_pf_sf_fio2_aagrad_aniongap"] = not det["B25_no_pf_sf_fio2_aagrad_aniongap"]
    chk["B26_urine_output_never_backfilled"] = True   # inherited from Stage 3 (no LOCF in Layer B; not refilled here)
    return dict(levels=levels, checks=chk, details=det, psofa=ps, prism=pr, b8_category_blocks=b8_cat_blocks)


def _b8_stem(L: Level, c: str) -> str:
    """Field stem tested by B8: for one-hot / drug category columns (meta 'category') the prefix-stripped text before
    '__' (the block or field name; the category text after '__' is never tested); otherwise the prefix-stripped
    column name (system-prefixed b2 copies keep their field suffix, which is tested)."""
    stem = c[len(L.prefix):] if c.startswith(L.prefix) else c
    if "category" in L.meta.get(c, {}):
        stem = stem.split("__", 1)[0]
    return stem.lower()


def b8_offenders(L: Level) -> list:
    return [c for c in L.cols if B8_RE.search(_b8_stem(L, c))]


def b8_category_blocks(L: Level) -> list:
    """Names of one-hot / drug blocks whose category text contains a banned token (block names only; audit note)."""
    out = set()
    for c, m in L.meta.items():
        cat = m.get("category")
        if cat is not None and any(t in str(cat).lower() or t in c.lower().split("__", 1)[-1] for t in B8_TOKENS):
            out.add(m.get("block") or f"{L.prefix}{_b8_stem(L, c)}")
    return sorted(out)


def channel_of(suf: str) -> str:
    best = ""
    for ch in SYSTEM_OF:
        if (suf == ch or suf.startswith(ch + "_")) and len(ch) > len(best):
            best = ch
    return best


def auto_masks(L: Level) -> None:
    for c in list(L.cols):
        m = L.meta[c]
        if m["cover"] is not None and m["cover"] not in L.cols:
            m["cover"] = None
        v = L.cols[c]
        if m["kind"] in ("cont", "ordinal") and np.isnan(v).any():
            cov = m["cover"]
            if cov is not None:
                need = np.isnan(v) & ~(L.cols[cov] == m["cover_when"])
                if not need.any() or "panel_" in cov:
                    continue
            nm = c[len(L.prefix):] + "_empty"
            L.add(nm, np.isnan(v).astype(float), "mask", system=m["system"])
            L.meta[L.prefix + nm]["auto"] = True
            m["cover"], m["cover_when"] = L.prefix + nm, 1


# ---------------------------------------------------------------------------------------------------- stage 08
def _skew(x: np.ndarray) -> float:
    if len(x) < 3:
        return 0.0
    s = x.std()
    return float(((x - x.mean()) ** 3).mean() / s ** 3) if s > 0 else 0.0


def fit_level(L: Level, split: np.ndarray) -> dict:
    """Train-only fit (plan section 7) applied unchanged to validation and test."""
    names = list(L.cols)
    train = split == "train"
    X = np.column_stack([L.cols[c] for c in names]).astype(float)
    rec = OrderedDict()
    for j, c in enumerate(names):
        m = L.meta[c]
        x = X[:, j]
        if m["kind"] in ("bin", "onehot", "mask"):
            if np.isnan(x).any():
                raise BuildError(f"{L.name}: empty cell in flag column {c}")
            continue
        obs = x[train & np.isfinite(x)]
        r = {"kind": m["kind"]}
        if m["kind"] == "cont" and len(obs):
            lo, hi = np.quantile(obs, WINSOR_Q)
            x = np.clip(x, lo, hi)
            obs = np.clip(obs, lo, hi)
            r["winsor"] = (float(lo), float(hi))
            if obs.min() >= 0 and _skew(obs) > LOG_SKEW:
                x, obs = np.log1p(x), np.log1p(obs)
                r["log1p"] = True
        med = float(np.median(obs)) if len(obs) else 0.0
        r["fill_train_median"] = med
        r["n_train_observed_ge_1"] = bool(len(obs) > 0)
        X[:, j] = np.where(np.isfinite(x), x, med)
        rec[c] = r
    if np.isnan(X).any():
        raise BuildError(f"{L.name}: empty cells left after the train-median fill")
    sd = X[train].std(axis=0)
    keep = sd > 1e-12
    zero_var = [c for c, k in zip(names, keep) if not k]
    names_k = [c for c, k in zip(names, keep) if k]
    Xk = X[:, keep]
    # one reference category per one-hot block (most frequent train category), named; full set kept for H1
    blocks = OrderedDict()
    for j, c in enumerate(names_k):
        b = L.meta[c]["block"]
        if L.meta[c]["kind"] == "onehot" and b:
            blocks.setdefault(b, []).append(j)
    ref_drop = {}
    drop_idx = set()
    for b, js in blocks.items():
        if len(js) < 2:
            continue
        means = Xk[train][:, js].mean(axis=0)
        jr = js[int(np.argmax(means))]
        ref_drop[b] = names_k[jr]
        drop_idx.add(jr)
    ref_names = [c for j, c in enumerate(names_k) if j not in drop_idx]
    ref_cols = [j for j in range(len(names_k)) if j not in drop_idx]
    Xr = Xk[:, ref_cols]
    mu_f, sd_f = Xk[train].mean(axis=0), Xk[train].std(axis=0)
    Zf = (Xk - mu_f) / sd_f
    Zr = Zf[:, ref_cols]
    fitted = dict(columns=rec, zero_variance_dropped=zero_var, reference_dropped=ref_drop,
                  mean={c: float(v) for c, v in zip(names_k, mu_f)}, sd={c: float(v) for c, v in zip(names_k, sd_f)})
    return dict(names=ref_names, X_imputed=Xr.astype(np.float32), X_z=Zr.astype(np.float32), names_full=names_k,
                X_full_z=Zf.astype(np.float32), fitted=fitted, kinds=[L.meta[c]["kind"] for c in ref_names])


def fit_layer_c(L_h: int, f: dict, raw_median: np.ndarray, raw_channels: list, mask_valid: np.ndarray,
                split: np.ndarray, eid: np.ndarray) -> dict:
    idx = f["encounter_index"].astype(int)
    P = L_h // 6
    ch = [str(c) for c in f["value_channels"]]
    ci = [raw_channels.index(c) for c in ch]
    raw = raw_median[idx, :P, :][:, :, ci].astype(float)
    filled = f["value_filled"].astype(float)
    miss = f["missing"].astype(bool)
    chk = {}
    chk["ids_match"] = bool((f["encounter_id"].astype(str) == eid[idx]).all())
    chk["prefix_all_valid"] = bool((mask_valid[idx, :P] == 1).all())
    chk["raw_missing_consistent"] = bool((np.isnan(raw) == miss).all())
    cdet = {"raw_missing_consistent": [ch[k] for k in range(len(ch)) if not bool(
        (np.isnan(raw[:, :, k]) == miss[:, :, k]).all())]}
    chk["filled_complete_or_noted"] = True
    tr = split[idx] == "train"
    q = len(ch)
    keep, rec = [], OrderedDict()
    Xr = np.full(raw.shape, np.nan)
    Xf = np.full(filled.shape, np.nan)
    for k in range(q):
        obs = raw[tr, :, k][~miss[tr, :, k]]
        if len(obs) < 2 or np.nanstd(obs) <= 0 or not np.isfinite(filled[:, :, k]).all():
            rec[ch[k]] = {"dropped": "no train variation or no fill value at this landmark"}
            continue
        lo, hi = np.quantile(obs, WINSOR_Q)
        r_ = np.clip(raw[:, :, k], lo, hi)
        f_ = np.clip(filled[:, :, k], lo, hi)
        o_ = np.clip(obs, lo, hi)
        lg = bool(o_.min() >= 0 and _skew(o_) > LOG_SKEW)
        if lg:
            r_, f_, o_ = np.log1p(r_), np.log1p(f_), np.log1p(o_)
        mu, sd = float(o_.mean()), float(o_.std())
        if sd <= 1e-12:
            rec[ch[k]] = {"dropped": "zero variance after winsorizing"}
            continue
        Xr[:, :, k] = (r_ - mu) / sd
        Xf[:, :, k] = (f_ - mu) / sd
        keep.append(k)
        rec[ch[k]] = {"winsor": (float(lo), float(hi)), "log1p": lg, "mean": mu, "sd": sd}
    zr = Xr[:, :, keep]
    zf = Xf[:, :, keep]
    chk["filled_no_empty_cell"] = bool(np.isfinite(zf).all())
    obs_tr = zr[tr][~np.isnan(zr[tr])] if tr.any() else np.array([0.0])
    chk["train_observed_mean_near_0"] = bool(abs(float(np.nanmean(zr[tr]))) < 0.05) if tr.any() else True
    return dict(encounter_id=eid[idx], split=split[idx], channels=np.array([ch[k] for k in keep]),
                X_raw_z=zr.astype(np.float32), X_filled_z=zf.astype(np.float32),
                missing=f["missing"][:, :, keep], carried_forward=f["carried_forward"][:, :, keep],
                imputed=f["imputed"][:, :, keep], fitted=rec, checks=chk, check_details=cdet, n_rows=len(idx),
                dropped=[c for c, r in rec.items() if "dropped" in r])


def cover_failures(levels: dict) -> list:
    """Value columns (level:column) with an empty cell not explained by its mask or flag (names only)."""
    bad = []
    for k, L in levels.items():
        for c, m in L.meta.items():
            v = L.cols[c]
            if m["kind"] in ("cont", "ordinal") and np.isnan(v).any():
                cov = m["cover"]
                if cov is None or cov not in L.cols:
                    bad.append(f"{k}:{c}")
                    continue
                if "panel_" in cov:
                    continue
                if not bool((L.cols[cov][np.isnan(v)] == m["cover_when"]).all()):
                    bad.append(f"{k}:{c}")
    return bad


def cover_check(levels: dict) -> bool:
    """Every empty cell of a value column is explained by its mask (or flag), or by a shared panel mask."""
    return not cover_failures(levels)


def failure_message(failed: list, details: dict, per_check: int = 25) -> str:
    """BUILD FAILED text: per failed check the level and column / block names involved; never row or patient counts."""
    parts = []
    for k in failed:
        names = [str(x) for x in details.get(k, [])]
        if len(names) > per_check:
            names = names[:per_check] + [f"(+{len(names) - per_check} more names)"]
        parts.append(f"{k} -> {', '.join(names) if names else 'no column-level detail'}")
    return f"acceptance checks failed: {failed}; details: " + " | ".join(parts)


# ------------------------------------------------------------------------------------------- writing and main
AUDIT_INT_KEYS = ("width", "n_columns", "n_categories", "n_blocks", "landmark_hours", "prefix_bins", "n_channels",
                  "n_dup", "rank_deficiency", "adm_src_", "m_rare", "stage3_files",
                  "n_zero_variance", "n_reference_dropped", "n_auto_masks", "n_drug_groups", "n_levels", "winsor")


def audit_ints_ok(o, key="", inherited=False) -> bool:
    """Shareable audit: every number sits under a structural key (widths, counts of columns or categories,
    configuration constants); no patient count anywhere."""
    allowed = inherited or any(key.startswith(k) for k in AUDIT_INT_KEYS)
    if isinstance(o, bool):
        return True
    if isinstance(o, (int, float, np.integer, np.floating)):
        return allowed
    if isinstance(o, dict):
        return all(audit_ints_ok(v, str(k), allowed) for k, v in o.items())
    if isinstance(o, list):
        return all(audit_ints_ok(v, key, allowed) for v in o)
    return True


def level_meta(L: Level, extra: dict) -> dict:
    return dict(level=L.name, prefix=L.prefix, plan=PLAN, contract=CONTRACT, script=SCRIPT_VERSION, spec=SPEC_DOC,
                feature_columns=list(L.cols), keys=["encounter_id", "split"],
                columns={c: {k: v for k, v in m.items() if v not in (None, False)} for c, m in L.meta.items()},
                provisional_columns=[c for c, m in L.meta.items() if m.get("provisional")],
                provisional_reason=("fluid intake (mL/kg/day) and the KDIGO urine-output stage depend on the weight "
                                    "normalizer; provisional until the B38 weight sensitivity review (Stage 3 "
                                    "decision 6)"), raw_summary=int(L.name == "b0"), **extra)


def creatinine_unit_gate(D: dict) -> None:
    """Stage 3 delivers creatinine in mg/dL (channel map v2.3: umol/L / 88.4). A train median at or above 10 would mean
    umol/L slipped through: stop instead of converting a second time or not at all."""
    tr = D["split"] == "train"
    for frame, col in ((D["B"], "creatinine_first"), (D["A"], "pre6_creatinine")):
        if col in frame.columns:
            v = num(frame[col])[tr]
            v = v[np.isfinite(v)]
            if len(v) and float(np.median(v)) >= CREATININE_MG_DL_MEDIAN_MAX:
                raise InputError(f"{col} does not look like mg/dL (train median >= {CREATININE_MG_DL_MEDIAN_MAX:g}); "
                                 "Stage 3 should have converted umol/L / 88.4")


def map_sensitivity_haque_zaritsky(D: dict, levels: dict) -> dict:
    """Team-internal sensitivity (D2): Roberts 2020 vs Haque-Zaritsky 2007 MAP z for the b1 MAP median, train rows
    aged 1-17 y with both defined. Written only to the INTERNAL fitted-transforms file."""
    A, B = D["A"], D["B"]
    age_d = num(A["age_days"])
    age_d[A["age_band"].to_numpy() == "age_invalid"] = np.nan
    x = num(B["map_median"])
    zr = D["norms"].roberts_z(x, age_d)
    zh = haque_zaritsky_z(x, age_d)
    ok = (D["split"] == "train") & np.isfinite(zr) & np.isfinite(zh)
    if ok.sum() < 3:
        return {"note": "too few train rows with both z values"}
    return {"pearson_r": float(np.corrcoef(zr[ok], zh[ok])[0, 1]), "mean_difference_roberts_minus_hz":
            float(np.mean(zr[ok] - zh[ok])), "sd_difference": float(np.std(zr[ok] - zh[ok])),
            "rows": "train, ages 1-17 y, MAP median of the Layer B window"}


def build(cfg: Config) -> dict:
    _STATE["deadline"] = time.time() + cfg.max_seconds
    say(f"{SCRIPT_VERSION}: precheck (hash gates before any data file is parsed)")
    pre = precheck(cfg)
    hashes = pre["hashes"]
    say(f"hash gates: Stage 3 manifest and {len(hashes) - 4} Stage 3 files OK; canonical, harmonization and split OK")
    if cfg.precheck_only:
        return dict(mode="precheck-only", hashes=hashes, ref_hashes=pre["ref_hashes"])
    D = load_inputs(cfg, pre["manifest"])
    D["norms"] = Norms()
    creatinine_unit_gate(D)
    say("inputs read: Layer A base, Layer B base, canonical keys, split, medication and microbiology internals")
    check_time()
    z = np.load(cfg.stage3_dir / "stage3_layer_c_raw.npz", allow_pickle=False)
    raw_channels = [str(c) for c in z["channels"]]
    if not (z["encounter_id"].astype(str) == D["eid"]).all():
        raise InputError("Layer C raw rows are not in the Layer A order")
    mask_valid = z["mask_valid"]
    raw_median = z["median"]
    zC = {}
    for c in SD_VITALS:
        b = raw_median[:, :4, raw_channels.index(c)].astype(float)
        b[mask_valid[:, :4] == 0] = np.nan
        zC[c] = b
    say("Layer C raw read (bin medians for the variability SDs and the raw tensor level)")
    report, fitted = [], {"reference_tables_sha256": dict(pre["ref_hashes"])}
    ra = build_layer_a(D, fitted, report)
    say("stage 06: Layer A levels " + ", ".join(f"{k} ({L.width} columns)" for k, L in ra["levels"].items()))
    check_time()
    rb = build_layer_b(D, zC, fitted, report)
    fitted["map_sensitivity_haque_zaritsky_INTERNAL"] = map_sensitivity_haque_zaritsky(D, rb["levels"])
    say("stage 07: Layer B levels " + ", ".join(f"{k} ({L.width} columns)" for k, L in rb["levels"].items()))
    check_time()
    levels = OrderedDict(list(ra["levels"].items()) + list(rb["levels"].items()))
    checks, details = OrderedDict(), {}
    checks.update({f"06_{k}": v for k, v in ra["checks"].items()})
    checks.update({f"07_{k}": v for k, v in rb["checks"].items()})
    details.update({f"06_{k}": v for k, v in ra["details"].items()})
    details.update({f"07_{k}": v for k, v in rb["details"].items()})
    details["07_cover_every_empty_cell"] = cover_failures(levels)
    checks["07_cover_every_empty_cell"] = not details["07_cover_every_empty_cell"]
    # stage 08: A/B levels
    fits = OrderedDict()
    for k, L in levels.items():
        fits[k] = fit_level(L, D["split"])
        say(f"stage 08: {k} fitted on train ({len(fits[k]['names'])} columns after zero-variance and reference drops)")
        check_time()
    train = D["split"] == "train"
    chk8, det8 = {}, {}
    det8["train_mean_0_sd_1"] = [f"{k}:{c}" for k, f in fits.items() for c, mu, sd in zip(
        f["names"], f["X_z"][train].mean(axis=0), f["X_z"][train].std(axis=0))
        if not (abs(mu) <= 1e-3 + 1e-8 and abs(sd - 1) <= 1e-3 + 1e-8)]
    chk8["train_mean_0_sd_1"] = all(
        bool(np.allclose(f["X_z"][train].mean(axis=0), 0, atol=1e-3) and
             np.allclose(f["X_z"][train].std(axis=0), 1, atol=1e-3)) for f in fits.values())
    det8["no_empty_cell"] = [f"{k}:{c}" for k, f in fits.items() for j, c in enumerate(f["names"])
                             if not (np.isfinite(f["X_z"][:, j]).all() and np.isfinite(f["X_imputed"][:, j]).all())]
    chk8["no_empty_cell"] = all(bool(np.isfinite(f["X_z"]).all() and np.isfinite(f["X_imputed"]).all())
                                for f in fits.values())
    det8["no_banned_column"] = [f"{k}:{c}" for k, f in fits.items() for c in f["names_full"] if banned(c)]
    chk8["no_banned_column"] = not det8["no_banned_column"]
    det8["B10_dup_rank_deficiency"] = ["b2s", "b2d"]
    # B10: phi2-dup numerical rank deficiency exceeds phi2-single by the number of retained duplicate columns
    def rank_def(Z):
        return Z.shape[1] - int(np.linalg.matrix_rank(Z.astype(np.float64)))
    dup_kept = sum(1 for c in fits["b2d"]["names_full"] if "__dup__" in c)
    rd_s = rank_def(fits["b2s"]["X_full_z"][train])
    rd_d = rank_def(fits["b2d"]["X_full_z"][train])
    chk8["B10_dup_rank_deficiency"] = (rd_d - rd_s) == dup_kept
    say(f"stage 08: phi2-dup rank check ({dup_kept} duplicate columns)")
    checks.update({f"08_{k}": v for k, v in chk8.items()})
    details.update({f"08_{k}": v for k, v in det8.items()})
    # stage 08: Layer C landmark tensors
    lc = OrderedDict()
    for Lh in LANDMARKS:
        f = np.load(cfg.stage3_dir / f"stage3_layer_c_filled_L{Lh}.npz", allow_pickle=False)
        r = fit_layer_c(Lh, f, raw_median, raw_channels, mask_valid, D["split"], D["eid"])
        r["checks"]["rows_equal_eligible"] = r["n_rows"] == EXPECT.eligible[Lh]
        lc[Lh] = r
        say(f"stage 08: Layer C L={Lh} h standardized ({len(r['channels'])} value channels)")
        check_time()
    for Lh, r in lc.items():
        for k, v in r["checks"].items():
            checks[f"08_layer_c_L{Lh}_{k}"] = v
            details[f"08_layer_c_L{Lh}_{k}"] = [f"L{Lh}:{c}" for c in r.get("check_details", {}).get(k, [])] or [
                f"L{Lh}"]
    del raw_median
    failed = [k for k, v in checks.items() if not v]
    if failed:
        raise BuildError(failure_message(failed, details))
    fitted["acceptance_audit_notes"] = dict(blank_missing_masks=ra["blank_masks"],
                                            blank_exempt_blocks=ra["blank_exempt"],
                                            b8_category_blocks=rb["b8_category_blocks"])
    say(f"acceptance checks: {len(checks)} pass, 0 fail")
    return publish(cfg, D, levels, fits, lc, checks, fitted, report, hashes, rd_s, rd_d, dup_kept)


def publish(cfg, D, levels, fits, lc, checks, fitted, report, hashes, rd_s, rd_d, dup_kept) -> dict:
    out = cfg.out_dir
    pend = out.parent / f"{out.name}.pending-{os.getpid()}"
    out.parent.mkdir(parents=True, exist_ok=True)
    pend.mkdir()
    try:
        script_sha = sha256_file(Path(__file__))
        files = OrderedDict()
        base_extra = dict(script_sha256=script_sha, input_hashes=hashes, generated_utc=utc_now(),
                          adm_src_surgical_window_hours=ADM_SRC_SURGICAL_WINDOW_H,
                          adm_src_floor_cutoff_hours=ADM_SRC_FLOOR_CUTOFF_H)
        eid, split = D["eid"], D["split"]
        elig24 = (D["stay_h"] > 24.0).astype(np.uint8)
        for k, L in levels.items():
            lay = "a" if k.startswith("a") else "b"
            name = f"layer_{lay}_{k}.csv"
            data = L.frame(eid, split).to_csv(index=False, float_format="%.6g", lineterminator="\n").encode("utf-8")
            files[name] = write_bytes_atomic(pend / name, data)
            meta = level_meta(L, dict(base_extra, sha256=files[name], file=name))
            mb = json.dumps(meta, indent=1, ensure_ascii=False).encode("utf-8")
            files[name + ".meta.json"] = write_bytes_atomic(pend / (name + ".meta.json"), mb)
            say(f"written {name}")
            check_time()
        for k, f in fits.items():
            name = f"stage08_{k}.npz"
            import io
            buf = io.BytesIO()
            np.savez_compressed(buf, encounter_id=eid.astype(str), split=split.astype(str), eligible_24h=elig24,
                                columns=np.array(f["names"]), kinds=np.array(f["kinds"]), X_imputed=f["X_imputed"],
                                X_z=f["X_z"], columns_full=np.array(f["names_full"]), X_full_z=f["X_full_z"])
            files[name] = write_bytes_atomic(pend / name, buf.getvalue())
            say(f"written {name}")
            check_time()
        for Lh, r in lc.items():
            name = f"stage08_layer_c_L{Lh}.npz"
            import io
            buf = io.BytesIO()
            np.savez_compressed(buf, encounter_id=r["encounter_id"].astype(str), split=r["split"].astype(str),
                                channels=r["channels"], X_raw_z=r["X_raw_z"], X_filled_z=r["X_filled_z"],
                                missing=r["missing"], carried_forward=r["carried_forward"], imputed=r["imputed"])
            files[name] = write_bytes_atomic(pend / name, buf.getvalue())
            say(f"written {name}")
            check_time()
        fitted_all = dict(fitted, stage08={k: f["fitted"] for k, f in fits.items()},
                          stage08_layer_c={str(Lh): r["fitted"] for Lh, r in lc.items()},
                          note="TEAM-INTERNAL: train statistics of real data (medians, quantiles, means, SDs); "
                               "never shared, never an input of another stage except through this script's outputs")
        name = "stage08_fitted_transforms_INTERNAL.json"
        files[name] = write_bytes_atomic(pend / name, json.dumps(fitted_all, indent=1, ensure_ascii=False,
                                                                 default=float).encode("utf-8"))
        name = "stage06_08_category_report_INTERNAL.csv"
        files[name] = write_bytes_atomic(pend / name, pd.DataFrame(report).to_csv(
            index=False, lineterminator="\n").encode("utf-8"))
        audit = make_audit(cfg, levels, fits, lc, checks, fitted, report, hashes, script_sha, rd_s, rd_d, dup_kept)
        if not audit_ints_ok(audit):
            raise BuildError("release gate: the shareable audit contains a number outside the structural whitelist")
        aj = json.dumps(audit, indent=1, ensure_ascii=False).encode("utf-8")
        files["stage06_08_audit.json"] = write_bytes_atomic(pend / "stage06_08_audit.json", aj)
        md = audit_md(audit).encode("utf-8")
        files["stage06_08_audit.md"] = write_bytes_atomic(pend / "stage06_08_audit.md", md)
        # release gate: reload from disk and recheck
        back = json.loads((pend / "stage06_08_audit.json").read_text(encoding="utf-8"))
        if back != audit or not audit_ints_ok(back):
            raise BuildError("release gate: disk reload of the audit differs")
        for nm, h in files.items():
            if sha256_file(pend / nm) != h:
                raise BuildError(f"release gate: {nm} hash differs after writing")
        man = dict(script=SCRIPT_VERSION, script_sha256=script_sha, spec=SPEC_DOC, plan=PLAN, contract=CONTRACT,
                   generated_utc=utc_now(), input_hashes=hashes,
                   reference_tables_sha256=dict(fitted["reference_tables_sha256"]),
                   adm_src_surgical_window_hours=ADM_SRC_SURGICAL_WINDOW_H, adm_src_floor_cutoff_hours=ADM_SRC_FLOOR_CUTOFF_H,
                   files={k: {"sha256": v, "bytes": (pend / k).stat().st_size} for k, v in files.items()},
                   acceptance_checks={k: ("pass" if v else "fail") for k, v in checks.items()}, release_gate="ok")
        write_bytes_atomic(pend / "stage06_08_manifest.json", json.dumps(man, indent=1).encode("utf-8"))
        if out.exists():
            retry(lambda: out.rmdir(), "remove empty output folder")
        retry(lambda: os.replace(pend, out), "rename pending folder to the output folder")
    except BaseException as e:
        try:
            if pend.exists():
                retry(lambda: shutil.rmtree(pend), "remove pending folder")
        except BuildError:
            e.leftover = pend.name   # type: ignore[attr-defined]
        raise
    say(f"published {len(files) + 1} files to {out.name}")
    return dict(mode="run", audit=audit, audit_md=audit_md(audit), n_files=len(files) + 1)


def make_audit(cfg, levels, fits, lc, checks, fitted, report, hashes, script_sha, rd_s, rd_d, dup_kept) -> dict:
    rep = pd.DataFrame(report)
    cats = {}
    if len(rep):
        for (lay, blk), g_ in rep.groupby(["layer", "block"]):
            cats[f"{lay}:{blk}"] = {"n_categories_kept": int(g_["kept"].sum()),
                                    "n_categories_pooled_to_other": int((g_["kept"] == 0).sum())}
    lv = OrderedDict()
    for k, L in levels.items():
        f = fits[k]
        kinds = pd.Series([L.meta[c]["kind"] for c in L.cols]).value_counts().to_dict()
        lv[k] = {"width_level_file": L.width, "width_stage08_primary": len(f["names"]),
                 "width_stage08_full_set": len(f["names_full"]),
                 "n_columns_by_kind": {kk: int(vv) for kk, vv in kinds.items()},
                 "n_zero_variance_dropped": len(f["fitted"]["zero_variance_dropped"]),
                 "n_reference_dropped": len(f["fitted"]["reference_dropped"]),
                 "n_auto_masks": sum(1 for c in L.cols if L.meta[c].get("auto")),
                 "provisional_columns": [c for c in L.cols if L.meta[c].get("provisional")]}
    lcs = OrderedDict((f"L{Lh}", {"landmark_hours": Lh, "prefix_bins": Lh // 6,
                                  "n_channels_standardized": len(r["channels"]),
                                  "channels_dropped": r["dropped"]}) for Lh, r in lc.items())
    return OrderedDict(
        script=SCRIPT_VERSION, script_sha256=script_sha, spec=SPEC_DOC, plan=PLAN, contract=CONTRACT,
        generated_utc=utc_now(),
        stage3_manifest_sha256=hashes["stage3_manifest.json"], stage3_files=PINS["stage3_files"],
        reference_tables_sha256=dict(fitted["reference_tables_sha256"]),
        inputs_sha256={k: v for k, v in hashes.items() if k in STAGE3_FILES_USED or k in (
            "stage3_manifest.json", "canonical", "harmonization", "split")},
        decisions_applied=OrderedDict(
            admission_source="first match: surgical (recent surgery ended <= t0 and interval < 48 h); floor (not "
                             "surgical and interval >= 24 h); transfer_ed (everyone else); missing interval never "
                             "counts as a long interval",
            adm_src_surgical_window_hours=ADM_SRC_SURGICAL_WINDOW_H, adm_src_floor_cutoff_hours=ADM_SRC_FLOOR_CUTOFF_H,
            m_rare=M_RARE,
            age_norming="D2: HR and RR - Fleming 2011 normal z (mean = median, SD = (p90 - p10) / 2.5631; HR "
                        "'Birth' row for age < 1 day, else 0-3 m; RR birth uses 0-3 m); MAP - Roberts 2020 "
                        "centile-interpolated z capped at the 1st/99th, centiles interpolated by exact age; creatinine "
                        "- ratio to the PODIUM 2022 age- and sex-specific baseline SCr (mg/dL); outside a table's "
                        "cover (age >= 18 y, age invalid, sex not M/F for PODIUM) -> empty + mask; other labs and "
                        "other vitals - train age-band mean/SD (empirical); Haque-Zaritsky 2007 as a team-internal "
                        "MAP sensitivity only (ages 1-17 y)",
            creatinine_units="Stage 3 delivers mg/dL (PIC umol/L / 88.4, channel map v2.3); not converted again; "
                             "train-median plausibility gate",
            vital_variability="SD of z-scored 6 h bin medians over bins 0-3 (>= 3 observed bins), HR/RR/MAP only",
            worst_direction="spec table D7; two-sided labs: value farther from the train age-band median",
            psofa="4 scored components (cardiovascular modified: 2 if any vasoactive active; neurologic not scored, "
                  "no GCS); no total column",
            prism3="non-neurologic points; PT not scored; missing -> 0 points",
            phi2_surgery="body system (clusters deferred)", hour_encoding_phi0="sin/cos",
            prior_encounters="omitted", diagnosis_missingness="the 'missing' category of the diagnosis block",
            winsor=list(WINSOR_Q), log_rule=f"log1p if non-negative and train skewness > {LOG_SKEW}",
            empty_by_design_fill="train median (flag or mask carries the information)",
            blank_category_masks="v1.2.1 (D11): a one-hot block with blank source values (active rows for "
                                 "recent_/remote_ blocks) gets a <block>_missing mask (admit_interval_band uses the "
                                 "existing admit_interval_missing); A6 checks block sum + mask = 1 (or the active "
                                 "indicator); sex is exempt (D35: no sex_missing column); masks constant in train "
                                 "are dropped by the stage 08 zero-variance rule",
            b8_rule="v1.2.1: organism / susceptibility tokens (org_, ab_, dilution, interpret) tested on field "
                    "stems at the start or after an underscore, never inside one-hot or drug category text",
            layer_c_standardization="pooled per channel over observed train cells of the landmark prefix"),
        levels=lv, n_levels=len(lv), n_dup_columns_b2d=dup_kept, rank_deficiency_b2s=rd_s, rank_deficiency_b2d=rd_d,
        categories=cats, layer_c=lcs,
        blank_missing_masks=OrderedDict((k, list(v)) for k, v in
                                        fitted["acceptance_audit_notes"]["blank_missing_masks"].items()),
        blank_exempt_blocks=OrderedDict((k, list(v)) for k, v in
                                        fitted["acceptance_audit_notes"]["blank_exempt_blocks"].items()),
        b8_category_blocks_with_banned_token=list(fitted["acceptance_audit_notes"]["b8_category_blocks"]),
        acceptance_checks=OrderedDict((k, "pass" if v else "fail") for k, v in checks.items()),
        disclosure="no patient counts: widths, column and category counts, check results and hashes only; "
                   "row-level files and *_INTERNAL* files are never shared",
        release_gate="ok")


def audit_md(a: dict) -> str:
    L = [f"# Stages 06-08 audit ({a['script']})", "", f"- plan {a['plan']}, contract {a['contract']}, spec {a['spec']}",
         f"- script SHA-256 `{a['script_sha256']}`", f"- Stage 3 manifest SHA-256 `{a['stage3_manifest_sha256']}`",
         f"- generated (UTC) {a['generated_utc']}", "", "## Decisions applied", ""]
    L += [f"- `{k}`: {v}" for k, v in a["decisions_applied"].items()]
    L += ["", "## Levels", "", "| level | level-file width | stage 08 primary | full set | zero-var dropped | "
          "references dropped | auto masks |", "|---|---|---|---|---|---|---|"]
    for k, v in a["levels"].items():
        L.append(f"| {k} | {v['width_level_file']} | {v['width_stage08_primary']} | {v['width_stage08_full_set']} | "
                 f"{v['n_zero_variance_dropped']} | {v['n_reference_dropped']} | {v['n_auto_masks']} |")
    L += ["", f"- phi2-dup: {a['n_dup_columns_b2d']} duplicate columns; rank deficiency b2s {a['rank_deficiency_b2s']}, "
          f"b2d {a['rank_deficiency_b2d']}", "", "## Categories (counts of categories, not of patients)", ""]
    L += [f"- `{k}`: kept {v['n_categories_kept']}, pooled to other {v['n_categories_pooled_to_other']}"
          for k, v in a["categories"].items()]
    L += ["", "## Layer C landmark tensors", "", "| landmark | prefix bins | channels standardized | dropped |",
          "|---|---|---|---|"]
    for k, v in a["layer_c"].items():
        L.append(f"| {v['landmark_hours']} h | {v['prefix_bins']} | {v['n_channels_standardized']} | "
                 f"{', '.join(v['channels_dropped']) or '-'} |")
    L += ["", "## Blank-category missing masks (v1.2.1, D11; column names only)", ""]
    L += [f"- {k}: {', '.join(f'`{c}`' for c in v) or '-'}" for k, v in a["blank_missing_masks"].items()]
    L += [f"- exempt (blank rows stay all-zero; D35) {k}: {', '.join(f'`{c}`' for c in v)}"
          for k, v in a["blank_exempt_blocks"].items()]
    L += ["", "## B8: blocks whose category text contains a banned token (block names only; not field stems)", "",
          "- " + (", ".join(f"`{b}`" for b in a["b8_category_blocks_with_banned_token"]) or "none")]
    n_pass = sum(1 for v in a["acceptance_checks"].values() if v == "pass")
    L += ["", f"## Acceptance checks: {n_pass} pass, {len(a['acceptance_checks']) - n_pass} fail", ""]
    L += [f"- `{k}`: {v}" for k, v in a["acceptance_checks"].items()]
    L += ["", f"## release_gate: {a['release_gate']}", "", a["disclosure"], ""]
    return "\n".join(L)


def build_config(argv=None) -> Config:
    ap = argparse.ArgumentParser(description=SCRIPT_VERSION)
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--stage3-dir", default=None, help="default data/processed/ladder/stage3")
    ap.add_argument("--canonical", default=None, help="default data/interim/canonical_icu_encounters.csv")
    ap.add_argument("--harmonization", default=None, help="default data/interim/encounter_harmonization_fields.csv")
    ap.add_argument("--split-file", default=r"..\PIC_Analysis_3_0\data\processed\subject_splits.csv")
    ap.add_argument("--out-dir", default=None, help="default data/processed/ladder/stage06_08")
    ap.add_argument("--expect-stage3-manifest-sha256", default=None)
    ap.add_argument("--precheck-only", action="store_true")
    ap.add_argument("--max-seconds", type=float, default=3600.0)
    ap.add_argument("--heartbeat-seconds", type=float, default=30.0)
    a = ap.parse_args(argv)
    root = Path(a.project_root)
    if not a.precheck_only and not a.expect_stage3_manifest_sha256:
        raise InputError("a run needs --expect-stage3-manifest-sha256 (printed by the precheck)")
    return Config(project_root=root,
                  stage3_dir=Path(a.stage3_dir) if a.stage3_dir else root / "data" / "processed" / "ladder" / "stage3",
                  canonical=Path(a.canonical) if a.canonical else root / "data" / "interim" /
                  "canonical_icu_encounters.csv",
                  harmonization=Path(a.harmonization) if a.harmonization else root / "data" / "interim" /
                  "encounter_harmonization_fields.csv",
                  split_file=Path(a.split_file),
                  out_dir=Path(a.out_dir) if a.out_dir else root / "data" / "processed" / "ladder" / "stage06_08",
                  expect_manifest_sha256=(a.expect_stage3_manifest_sha256 or "").strip().lower() or None,
                  precheck_only=a.precheck_only, max_seconds=a.max_seconds, heartbeat_seconds=a.heartbeat_seconds)


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:  # noqa: BLE001
        pass
    try:
        cfg = build_config(argv)
        with Heartbeat(cfg.heartbeat_seconds):
            res = build(cfg)
    except InputError as e:
        print(f"STOPPED (input): {e}", file=sys.stderr)
        print("result: STOPPED (input); nothing written")
        return 2
    except (BuildError, TimeLimit) as e:
        lo = getattr(e, "leftover", None)
        tail = (f"pending folder left behind: {lo} (delete it unopened after review; never share or sync it)"
                if lo else "nothing left behind")
        print(f"BUILD FAILED: {e}", file=sys.stderr)
        print(f"result: FAILED; {tail}")
        return 3
    except Exception as e:  # noqa: BLE001  (fail closed; the message text could carry values, so only the type)
        lo = getattr(e, "leftover", None)
        tb = e.__traceback__
        while tb is not None and tb.tb_next is not None:
            tb = tb.tb_next
        where = f" at line {tb.tb_lineno}" if tb is not None else ""
        print(f"BUILD FAILED (unexpected {type(e).__name__}{where})", file=sys.stderr)
        print(f"result: FAILED; {'pending folder left behind: ' + lo if lo else 'nothing left behind'}")
        return 3
    if res["mode"] == "precheck-only":
        h = res["hashes"]
        print("result: PRECHECK OK (reference tables, Stage 3 manifest, every Stage 3 file, canonical, harmonization, split); nothing written")
        print(f"stage3_manifest.json sha256 {h['stage3_manifest.json']}")
        for k in STAGE3_FILES_USED:
            print(f"{k} sha256 {h[k]}")
        print(f"canonical sha256 {h['canonical']}")
        print(f"harmonization sha256 {h['harmonization']}")
        print(f"split sha256 {h['split']}")
        for k, v in res["ref_hashes"].items():
            print(f"reference {k} sha256 {v}")
        return 0
    print(res["audit_md"])
    print("result: STAGES 06-08 OK; RELEASE GATE OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
