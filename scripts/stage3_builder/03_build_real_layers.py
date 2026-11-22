#!/usr/bin/env python3
"""
Stage 3 real-layer builder (Foundations / PIC, plan v13, contract v13): base Layers A, B and C from the raw tables.
Implements `stage3_builder_spec.md` v1.0.4 (Nov 17, 2026).
v1.0.4 (after the real Stage 2 v1.0.3 release; Vinai, Nov 17, 2026): (1) the joint (cross-release) audit runs against
the published Stage 2 v1.0.3 release (feature_sparsity_frequency_audit_v1.0.3.json, hash-pinned) with the Stage 2
v1.0.3 suppression module v1.5 (stage2_v1_0_3_suppress.py, hash-pinned; heartbeat, capped local LPs, convergence
stop); (2) links file v3 (stage3_links_v3.csv) = links v2 + every Stage 2 v1.0.3 internal equality (L5a/L5b/L5c/L6/L7)
as s2: rows; (3) a hash-pinned TEAM-INTERNAL constants file (--stage2-constants, stage3_stage2_constants_TEAM_INTERNAL.csv):
Stage 2 cells a reader already knows from outside the release (the n_ge18 structural zeros, the plan v13 stop-event
invalid-DOB count) and the cells they fix exactly ('derived'), entered as known values; values are never printed or
written, outputs record the row count only; (4) per-channel token table: 'positive+' / 'positive++' read as positive
for stool_occult_blood, urine_leukocyte_esterase, urine_nitrite and csf_protein_pandy only; '<3.00' read as negative
for resp_panel_mycoplasma only, with the Layer B flag column resp_panel_mycoplasma_below_cutoff (not a feature);
(5) labels v1.0.4 (audit provenance block release_v1_0_4; cross_release_v1_0_4). Drop list v1 unchanged.
v1.0.3 (after the real v1.0.2 run; Vinai, Nov 17, 2026): (1) ICD-O-3 morphology codes are recognized in every realistic
written form of the first code (the clinical-edition M98000/3 with its behaviour suffix was missed by v1.0.2, which
matched the whole field; all behaviour codes map to oncology); (2) ketones_urine is dropped as a feature channel
(92 feature channels; map v2.3 unchanged), listed in the audit as excluded with its reason; (3) labels v1.0.3 (audit
provenance block release_v1_0_3). The joint check with Stage 2 v1.0.2 (links v2, drops v1, suppression module v1.3)
is unchanged.
v1.0.2 (cross-release disclosure; spec §16): the aggregate audit is gated jointly with the published Stage 2 v1.0.2
release, a hash-pinned input (--stage2-release, --expect-stage2-sha256; any mismatch stops with exit 2). Its published
cells enter as known values, its '<10', 'S' and hidden cells as unknowns, linked to the Stage 3 cells by the
hash-pinned links file v2 (stage3_links_v2.csv); subtrees listed in the hash-pinned drop list (stage3_drops_v1.csv)
are withheld as a whole but stay unknowns in every equation. Cells already exposed by Stage 2 v1.0.2 alone are
reported (flagged for Stage 2 v1.0.3) and may not be narrowed further. Audit v1.1: age shows only valid / invalid;
the admission-interval missing / non-negative pair is withheld. Layers and other outputs are unchanged from v1.0.1.
v1.0.2 also folds in (unreleased, label unchanged): the growth-reference weight normalizer (spec §10) and Stage 3 review
decisions 2-9 (Vinai, Nov 17, 2026; spec §15): respiratory rate 0 = missing and below/above-gate counts; train-median
fallback for empty landmarks; unparsed qualitative text listed team-internally (strings only) and flagged; base
excess / bicarbonate de-duplicated by draw; fluid intake / KDIGO flagged provisional; drugs at t0 from t0 - 6 h with no
future ENDDATE; ICD-O morphology -> oncology; no EDREGTIME feature and no VISIT_ID linkage.
v1.0.1 (Windows file safety; outputs otherwise identical to v1.0): every create/replace/delete in the output area goes
through fs_op() (short retries on transient sharing violations from cloud sync, antivirus or the indexer); memory
maps are released and checked before their files are zipped or deleted; on any failure the frames of the failed
build are cleared (so memory maps close) before the pending folder is removed; BUILD FAILED names the operation,
errno/WinError and the path relative to the output folder's parent, never values; a pending folder that cannot be
removed is reported instead of "nothing left behind".

One common builder of the real layers (V8): Layer A base fields at t0, Layer B base summaries over
[t0, min(OUTTIME, t0 + 24 h)), and the Layer C full grid of all cohort encounters x 56 six-hour bins x the 92 feature
channels of channel map v2.3 (94 map channels; weight 1014 is a normalizer only, V12; ketones_urine excluded in v1.0.3), with mask_valid, per-channel
masks created before carry-forward, one filled copy per landmark (14), the weight normalizer, and an aggregate audit.

Safety (fail closed):
  * every input is checked against its full 64-character SHA-256 (and, where listed, its row count) before it is read;
    D_ICD_DIAGNOSES and SURGERY_VITAL_SIGNS are gated on the full hashes recorded at the precheck (V14);
  * never read: DIAGNOSES_ICD, OR_EXAM_REPORTS, the retired layer_c_features.csv, any *_raw_LOCAL* file;
  * cohort, split, deaths and the 14 landmark counts must reproduce; invalid ages may not exceed 1%;
  * every output is built in memory and written to a pending folder; acceptance checks and the aggregate release gate
    (primary <10, complementary S, residual check, LP interval audit, declared links, joint with the Stage 2 v1.0.3
    release; the hash-pinned Stage 2 v1.0.3 suppression module v1.5) run before anything is moved into place; any failure removes every file of the run;
  * progress goes to stderr with counts only (never values, identifiers or timestamps); a heartbeat line appears
    whenever 30 s pass without one (v1.0.4), so a silent console for 3 minutes means a hang.
Exit codes: 0 ok; 2 STOPPED (input): nothing written; 3 BUILD/RELEASE GATE FAILED: nothing left behind.
Run real data only as written in the spec (Copilot steps), unedited.
"""
from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import shutil
import sys
import time
import traceback
import unicodedata
import weakref
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

SCRIPT_VERSION = "03_build_real_layers v1.0.4"
SPEC_DOC = "stage3_builder_spec v1.0.4"
AUDIT_SPEC = "stage3_audit v1.1 (stage3_builder_spec v1.0.4; plan v13; contract v13)"
RELEASE_KEY = "release_v1_0_4"      # audit provenance block (the pinned module's own key is release_v1_0_3)
DROP_NOTE_KEY = "dropped_subtrees_v1_0_4"
CROSS_NOTE_KEY = "cross_release_v1_0_4"
FLAG_KEY = "flagged_for_stage2_v1_0_4"   # manifest: Stage 2 cells already exposed by Stage 2 alone
PLAN_VERSION, CONTRACT_VERSION = "v13", "v13"

# ----------------------------------------------------------------------------------------------------------------------
# Pinned inputs (contract v13 Table P; channel map v2.3). Full SHA-256 only; a changed input needs a new pin.
# ----------------------------------------------------------------------------------------------------------------------
PINS = {
    "canonical": "e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418",
    "harmonization": "ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d",
    "split": "60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31",
    "complexity_map": "5c993b88adfc32ae82c98de19d51c8295d07d34efe7ab29d83a131ca0ab9cfb5",
    "phi1_medians": "e73c8b5a8d780b51edfd0935abdb05335f57e3a0cfc0c4c52eaaa85888db0047",
    "channel_map": "7541d7ea67aa9631c984a5b2c7bb63bf17e6e38b9013a7b09ea42fe754536ec2",   # map v2.3
    "drug_map": "c557c3a12a9e1bbaf4e95c54b75affd4b3310cfd498cb6a485a13ae4e8ff37d7",
    "specimen_map": "caed429de4f2d484b40785d31ecbe3375f308db58a332d108c247b5ce331e3a5",
    "organism_map": "087e5ba3f7263b55e78b1fc08e41e6165bfa2f8282e94f10a249eca1eb6817ce",
    "neuro_map": "3cc8f4acfe9df29037b60d903ad4204c40db5bc76923de8b64234ccb2553b9d9",
    "growth_ref": "8e53ca2d4fd9f053c2fa56b153e90b4ea63d090804add22414678fd80b88b54e",  # weight-for-age ref v1 (§10)
    "v17_module": "1de2b17b2d894caff6a8362a9c72a458967691a8d2c7a533db863e87ad8dad94",  # public copy (cohort pins omitted); original private pin was ce894b20ad009369c2f21cd6db61a0379cd2db888dbf22309b21c393409de6d0
    "suppress_module": "63d903d1ceef80a61f257a08aff4acd4a8c7fd5c2361d1e36de1c0cf603dc751",   # Stage 2 v1.0.3 module v1.5
}
# v1.0.2: the published Stage 2 release is a hash-pinned input of the joint (cross-release) disclosure audit; its
# published cells are known constants, its '<10'/'S'/hidden cells unknowns. v1.0.4: Stage 2 v1.0.3 (real run Nov 17, 2026;
# data\interim\stage2_v1.0.3\feature_sparsity_frequency_audit_v1.0.3.json).
STAGE2_RELEASE_SHA256 = "6621cc3fcc79c7d5711a56fba0510790748f0d8ddbb1e8a3ad73a974fe88aef2"
STAGE2_RELEASE_LABEL = "Stage 2 v1.0.3 (feature_sparsity_frequency_audit_v1.0.3.json)"
XR = "s2:"   # prefix of Stage 2 cells in the links file and in the joint equation system
# Refused on sight (superseded maps; v2.2 and v2.1 carry the pf_ratio differencing note)
REFUSED_SHA256 = {
    "1e7c6c98f20518a824483e690002e22a8e8c5c8d2c5b81233c2182deffa8077a": "channel map v2.2 (superseded by v2.3)",
    "7d1c060a228ab0772f35ae27e43a30971bcc13c36c98fcce2e369640b91bc899": "channel map v2.1 (superseded)",
    "56e8bf72fcbcaac765a39c2e705aedef61323e94841dbdcc8b2657191af370b6": "channel map v1 (superseded)",
}
RAW_PINS = {  # table -> (sha256, rows or None); contract v13 Table P
    "PATIENTS": ("85b016fcc3050e5b5b6200103ed74ad1798720c881e26c60f1769be9f336ae0d", None),
    "ADMISSIONS": ("1451ae8ea0ff463269b348522047a4bc99efd2665c770e9c5cef0ffe2647e1ae", None),
    "ICUSTAYS": ("449119e698acaaf948921c34375db535b9fd3f6eb6e996a86069da95e90eea1b", None),
    "SURGERY_INFO": ("04450bff4592d676a7ebe2dad22d9e6c2a6f8edd07dff223b182eb0002a38d83", None),
    "CHARTEVENTS": ("fe209d47d15825f180624f2cefb9bc38f8513e1ac0b3709d57aef37ae88f328e", None),
    "LABEVENTS": ("784ecb929a358ee1a18fc9beafbdf67f30a91646b780d9983714a5e9ce0018c5", None),
    "PRESCRIPTIONS": ("61fa43dcb1ed977f4219c902b01977da2d898a16a2c33640997ad8337d2eb2c6", None),
    "MICROBIOLOGYEVENTS": ("e19cf9d41732f74911039acae6b2cd7bf5c5dd65179cbd31aa8634dd1e6bb45c", None),
    "D_ITEMS": ("a7b65cf5d61ca84f6ce5d07d3775b01116473176a980463d4fe9df1dd12e1591", None),
    "D_LABITEMS": ("d4557986fe9c5475d0db6506fdf12801ecfb0e2e7a3c01caf6ac86bbab69ed60", None),
}
PRECHECK_PINNED = ("D_ICD_DIAGNOSES", "SURGERY_VITAL_SIGNS")   # full hashes given on the command line (V14)
NEVER_READ = ("DIAGNOSES_ICD", "OR_EXAM_REPORTS", "layer_c_features")
MAP_FILES = {"channel_map": "feature_channel_map_v2_3.csv", "drug_map": "drug_group_map_v1.csv",
             "specimen_map": "specimen_map_v1.csv", "organism_map": "organism_map_v1.csv",
             "neuro_map": "neuro_admitting_dx_icd10_v1.csv", "growth_ref": "growth_weight_for_age_ref_v1.csv"}
MAP_ROWS = {"channel_map": 686, "drug_map": 343, "specimen_map": 71, "organism_map": 101, "neuro_map": 67,
            "growth_ref": 3926}


@dataclass
class Expect:
    """Reproduction targets (contract v13 §H1, §H2, §H4, §H5).
    Cohort, split, death, eligibility and surgery-count pins are omitted from this public copy.
    Restore them from a local configs/cohort_pins.json before a real PIC run."""
    encounters: int | None = None
    subjects: int | None = None
    died: int | None = None
    splits: dict = field(default_factory=dict)
    eligible: dict = field(default_factory=dict)
    layer_c_surgeries: int | None = None
    max_invalid_age_frac: float = 0.01
    channels: int = 94                 # channels defined by map v2.3 (unchanged, hash-pinned)
    feature_channels: int = 92         # v1.0.3: 94 - weight (normalizer only) - ketones_urine (excluded)
    systems: dict = field(default_factory=lambda: {
        "Cardiovascular": 8, "Respiratory": 8, "Renal": 7, "Hepatic": 7, "GI": 3, "Hematologic": 18,
        "Infectious": 23, "Neurologic": 5, "Electrolytes/homeostasis": 7, "Endocrine/metabolic": 8})
    blood_gas_items: frozenset = frozenset({"5211", "5218", "5223", "5224", "5226", "5227", "5230", "5235", "5236",
                                            "5237", "5238", "5239", "5248", "5249", "5252"})


MIN_CELL = 10
SPLITS = ("train", "validation", "test")
BIN_H, N_BINS, MAX_H, LAYER_B_H = 6.0, 56, 336.0, 24.0
LANDMARKS = (24, 48, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288, 312, 336)
CORE_LANDMARKS = (24, 48, 72, 120, 168, 240, 336)
PRE_ICU_H = 6.0                  # Layer A pre-ICU block [t0 - 6 h, t0) (decision 19); 24 h sensitivity also written
PRE_ICU_SENS_H = 24.0
LAB_RELEASE_H = 1.0              # labs available at CHARTTIME + 1 h (decision 15)
KEEP_LO_H = -24.0                # earliest availability time kept (24 h pre-ICU sensitivity; weight window)
EVER_LOOKBACK_H = 6.0            # ever_received pre-t0 lookback (V20)
VASO_CUTOFF_H, DEFAULT_CUTOFF_H = 12.0, 24.0
RECENT_LO_H, RECENT_HI_H = -24.0, 6.0   # recent surgery T in [t0 - 24 h, t0 + 6 h] (S10)
SURGERY_MAX_H = 24.0             # Q12: surgery durations above 24 h fall back to the anesthesia window
BAND_SHORT, BAND_LONG = 0.75, 1.33
UO_CAP_H = 24.0
MED_GROUPS = ("vasoactive", "antibiotic", "antifungal", "antiviral", "sedative_analgesic", "non_opioid_analgesic",
              "antiemetic", "diuretic", "osmotherapy", "electrolyte_replacement", "steroid", "nmb", "aed")
DRUG_MAP_GROUPS = (set(MED_GROUPS) - {"osmotherapy"}) | {"none", "osmotic_conditional"}
SPECIMEN_GROUPS = ("blood", "urine", "respiratory", "csf", "stool", "ureaplasma", "culture_site_unknown", "other",
                   "unmapped")
AGE_BANDS = ("neonate", "infant", "toddler", "preschool", "school", "adolescent", "ge18", "age_invalid")
# Weight normalizer (spec §10; Vinai Nov 17, 2026): charted 1014, else the sex-specific weight-for-age median (LMS M)
# of the frozen growth reference at exact age (linear in age days between table ages, no extrapolation): China 2009
# (Li et al. 2009; LMS of Zong & Li 2013 Table 2; 0-216 months) first, then WHO (2006 standards by day 0-1856 d; 2007
# reference by month 61-120); age_invalid, sex not M/F or age beyond every reference -> the train-median tiers.
WEIGHT_RULE_TEXT = ("charted weight 1014 in [t0 - 24 h, end of the Layer B window); else the sex-specific median "
                    "weight-for-age (LMS M) at exact age, linear between table ages, from the China 2009 growth "
                    "standard (Li et al., Zhonghua Er Ke Za Zhi 2009;47(7):487-492; LMS of Zong & Li, PLoS ONE "
                    "2013;8(3):e59569 Table 2; 0-18 y), else WHO (Child Growth Standards 2006, 0-5 y; Reference 2007, "
                    "5-10 y); else (age_invalid, sex not M/F, or age beyond 18 y) the train medians by age band x sex, "
                    "age band, sex, overall (n >= 10; never by band for age_invalid). No gestational-age correction.")
TRAIN_TIERS = ("train_tier_age_sex", "train_tier_age_band", "train_tier_sex_only", "train_tier_overall")
WEIGHT_SOURCES = ("charted", "cn_ref", "who_ref") + TRAIN_TIERS
GROWTH_REFS = ("cn_li2009", "who")            # primary chain order; 'who' alone is the WHO sensitivity arm
GROWTH_REF_SOURCE = {"cn_li2009": "cn_ref", "who": "who_ref"}
GROWTH_REF_COLUMNS = ["reference", "component", "sex", "age_days", "age_months", "L", "M", "S"]

# Carry-forward families (contract v13 LOCF table; V21, V28; Stage 3 Q9). Cap in bins; hours for time since last.
NO_LOCF_VALUE = {"heart_rate", "map", "respiratory_rate", "temperature", "spo2", "sbp", "dbp", "pain_score",
                 "glucose_bedside"}
LOCF_4 = {"creatinine", "bun", "hemoglobin", "platelets", "wbc", "ast", "bilirubin_total", "bilirubin_direct",
          "albumin", "ggt", "inr", "ptt", "fibrinogen", "d_dimer", "anc", "lymphocytes_abs", "cystatin_c", "amylase",
          "lipase", "ldh", "crp", "urine_output", "fluid_intake", "kdigo_uo_stage"}
LOCF_1 = {"ph", "lactate", "bicarbonate", "sodium", "potassium", "glucose", "chloride", "calcium_total",
          "calcium_ionized", "magnesium", "phosphate", "pco2", "po2", "base_excess", "sao2_bloodgas",
          "cortisol", "free_t4", "insulin", "osmolality", "csf_wbc", "csf_rbc", "csf_protein_pandy",
          "stool_occult_blood", "urine_wbc", "urine_nitrite", "urine_leukocyte_esterase", "resp_panel_rsv",
          "resp_panel_parainfluenza", "resp_panel_influenza", "resp_panel_metapneumovirus", "resp_panel_mycoplasma",
          "systemic_viral_cmv", "systemic_viral_ebv", "systemic_viral_hsv", "systemic_viral_enterovirus",
          "bacteremia_pcr"}
# Qualitative channels: value = positivity (0/1) from a fixed token table; anything else is measured, value empty
QUALITATIVE = {"stool_occult_blood", "urine_nitrite", "urine_leukocyte_esterase", "csf_protein_pandy",
               "resp_panel_rsv", "resp_panel_parainfluenza", "resp_panel_influenza", "resp_panel_metapneumovirus",
               "resp_panel_mycoplasma", "systemic_viral_cmv", "systemic_viral_ebv", "systemic_viral_hsv",
               "systemic_viral_enterovirus", "bacteremia_pcr"}
QUAL_NEG = {"阴性", "阴性(-)", "阴性（-）", "-", "(-)", "（-）", "neg", "negative", "未检出", "not detected", "nd"}
QUAL_POS = {"阳性", "阳性(+)", "阳性（+）", "+", "(+)", "（+）", "++", "+++", "++++", "pos", "positive", "检出", "detected",
            "弱阳性", "weak positive", "1+", "2+", "3+", "4+"}
# v1.0.4 (Vinai, Nov 17, 2026; from the team-internal unparsed list of the real v1.0.3 build): per-channel additions to
# the token table, read only on the named channels (same normal form: trimmed, lower case, spaces removed).
QUAL_POS_BY_CHANNEL = {ch: {"positive+", "positive++"} for ch in
                       ("stool_occult_blood", "urine_leukocyte_esterase", "urine_nitrite", "csf_protein_pandy")}
QUAL_NEG_BY_CHANNEL = {"resp_panel_mycoplasma": {"<3.00"}}
# a negative read from a below-cutoff result: flagged per encounter in Layer B (<channel>_below_cutoff; not a feature)
QUAL_BELOW_CUTOFF = {"resp_panel_mycoplasma": {"<3.00"}}
BELOW_CUTOFF_FLAG_SUFFIX = "_below_cutoff"
NORMALIZER_ONLY = {"weight"}
# v1.0.3 (Vinai, Nov 17, 2026): channels of map v2.3 dropped as features. The map file is unchanged (its pin, its 94
# channels and the Table 5.1 system counts are still checked); an excluded channel's items are not read, and the
# channel is in no layer, filled copy, train median or audit count. Listed by name with its reason in the audit.
EXCLUDED_CHANNELS = {
    "ketones_urine": ("dropped as a feature (Vinai, Nov 17, 2026): urine ketones (LABEVENTS 5356) had no train median "
                      "at any landmark in the real v1.0.2 build (too few train values), so its filled copies could not "
                      "be imputed; the map v2.3 row is kept unchanged"),
}
SPECIAL_CHART = {"urine_output": "1008", "fluid_intake": "1009", "weight": "1014"}
MAIN_ITEM_RULE = {"heart_rate": "1003"}   # D9: 1003 preferred over 1002 at the same time
# Stage 3 review decisions (Vinai, Nov 17, 2026; spec §15 items 2-9)
# 2: a respiratory rate of 0 is "not recorded" (missing, status zero_missing), not an out-of-gate value; the audit
#    splits out-of-gate rows into below / above the gate (overall scope only) for these channels.
ZERO_AS_MISSING = ("respiratory_rate",)
GATE_SIDE_SPLIT = ("respiratory_rate", "wbc")
# The split is released only if the overall out-of-gate total is >= 20: below that both sides can be < 10, and a total
# of two '<10' cells (published, or suppressed and so >= 10 to a reader) pins them; no suppression pattern protects
# it. The rule depends only on the out-of-gate total (published, or 'S' >= 10: the rule then tells a reader 10-19,
# still 10 integers).
GATE_SIDE_MIN_TOTAL = 2 * MIN_CELL
# 5: base excess and bicarbonate are de-duplicated by draw like pCO2 and pH (one value per encounter and CHARTTIME,
#    the map's main item kept: 5211 and 5248).
DEDUPE_BY_DRAW_EXTRA = ("base_excess", "bicarbonate")
# 6: still produced, flagged provisional in the feature metadata until the B38 weight sensitivity is reviewed.
PROVISIONAL_CHANNELS = ("fluid_intake", "kdigo_uo_stage")
PROVISIONAL_NOTE = ("provisional (Vinai, Nov 17, 2026): fluid intake (mL/kg/day) and the KDIGO urine-output stage "
                    "depend on the weight normalizer; still produced, to be confirmed after the B38 weight "
                    "sensitivity review")
# 8: ICD-O-3 morphology codes in the admitting-diagnosis field go to the oncology group when D_ICD_DIAGNOSES does not
#    map them; 'missing' (no code) stays its own level. v1.0.3: every realistic written form of the FIRST code is
#    recognized: 'M' + 4-digit histology 8000-9989, then optionally one extension digit (the Chinese national
#    clinical-edition form M98000) and/or a behaviour digit 0/1/2/3/6/9 written after '/' (M9800/3, M98000/3) or as a
#    6th digit (M980003); '.' after the 4 histology digits (M9800.0), an optional '-' after 'M', upper/lower case,
#    full-width characters and spaces around separators are accepted. v1.0.2 matched the whole field, so the
#    clinical-edition form with its behaviour suffix (M98000/3) failed while its first code 'M98000' was listed as
#    unmapped. All behaviour codes (0 benign, 1 uncertain, 2 in situ, 3/6/9 malignant) map to oncology (decision 8).
ICDO_MORPH_RE = re.compile(r"M-?(?P<histology>[89][0-9]{3})(?:\.?(?P<ext>[0-9])(?P<beh6>[012369])?)?"
                           r"(?:/(?P<behaviour>[012369]))?")
ICDO_CODE_SPLIT = re.compile(r"[,;|，；、\s]+")      # separators between codes (as norm_icd, but '/' is kept)
ICDO_LEVEL = "icdo_morphology"
ONCOLOGY_DOMAIN = "oncology_hematology"
# 4: unparsed qualitative text, listed for the team (strings only, no counts) in a separate TEAM-INTERNAL file
UNPARSED_MIN = 10                 # a string seen fewer times on train rows is not listed
UNPARSED_MAX_CHARS = 32           # longer strings are truncated (…)
UNPARSED_SIMPLE_NUMERIC = re.compile(r"[<>≤≥]?[0-9]{1,4}(?:\.[0-9]{1,3})?[%+]?|1:[0-9]{1,5}|[0-9]\+")
UNPARSED_FILE = "stage3_unparsed_text_TEAM_INTERNAL.md"
# Exact unit_conversion texts of channel map v2.3 -> (per-row unit factors to the gated unit or None, final factor).
# A per-row factor of None drops the row (counted). Any other text stops the run.
CONVERSIONS = {
    "mmol/L -> mg/dL x18.016": (None, 18.016),
    "g/L -> g/dL /10": (None, 0.1),
    "urea mmol/L -> BUN mg/dL x2.801": (None, 2.801),
    "mmol/L -> mg/dL x4.008": (None, 4.008),
    "umol/L -> mg/dL /88.4": (None, 1.0 / 88.4),
    "umol/L -> mg/dL /17.1": (None, 1.0 / 17.1),
    "mmol/L as is; 'ug/dl' rows dropped": ({"mmol/l": 1.0, "": 1.0, "ug/dl": None}, 1.0),
    "mmol/L -> mg/dL x3.097": (None, 3.097),
    "g/L /10 -> g/dL": (None, 0.1),
    "kPa -> mmHg x7.5 if any": ({"kpa": 7.5, "mmhg": 1.0, "": 1.0}, 1.0),
    "per-row VALUEUOM: g/dL x10 -> g/L for the 20-250 g/L gate; final unit g/dL = g/L /10":
        ({"g/l": 1.0, "g/dl": 10.0, "": 1.0}, 0.1),
    "mIU/L = uU/mL": (None, 1.0),
    "ug/dL -> nmol/L x27.59": (None, 27.59),
    "per-row VALUEUOM: ng/mL x1; pg/mL /1000 -> ng/mL": ({"ng/ml": 1.0, "": 1.0, "pg/ml": 0.001}, 1.0),
}
AGGREGATION_TEXT = {
    "1008": "mL/kg/h = volume / weight / hours since the previous entry (capped at 24 h; first entry: hours since t0); "
            "value in the bin containing its charting time; no back-fill; carried forward at most 24 h; missing days stay missing",
    "1009": "mL/kg/day: daily total of charted intake / weight (charted weight 1014, else imputed from age and gender, A2.12)",
}
DERIVED_TEXT = {"map": "MAP = (SBP + 2 x DBP) / 3",
                "kdigo_uo_stage": "stage 3 if urine output < 0.3 mL/kg/h over 24 h or anuria >= 12 h; daily average "
                                  "< 0.5 mL/kg/h as a looser stand-in for stages 1-2; else 0"}
MAP_STATUS_USED = ("assigned", "assigned_merged", "assigned_merged_conditional", "assigned_merged_indicator_only",
                   "derived")
_NULLS = {"", "nan", "none", "null", "nat"}


class InputError(Exception):
    """Exit 2: input, pin, schema or reproduction problem. Nothing is written."""


class GateError(Exception):
    """Exit 3: an acceptance check or the release gate failed. Every file of the run is removed."""


class TimeLimit(Exception):
    """--max-seconds reached (exit 3, nothing left behind)."""


# ----------------------------------------------------------------------------------------------------------------------
# Logging (stderr, counts only)
# ----------------------------------------------------------------------------------------------------------------------
_CLOCK = {"t0": time.monotonic(), "deadline": math.inf}


_SM_REF: dict = {"m": None}   # v1.0.4: the loaded suppression module (heartbeat clock and print lock)


def log(msg: str) -> None:
    m = _SM_REF["m"]
    if m is None:
        print(f"[{time.monotonic() - _CLOCK['t0']:8.1f}s] {msg}", file=sys.stderr, flush=True)
        return
    with m._PRINT_LOCK:            # one line at a time with the heartbeat thread; resets its silence clock
        now = time.monotonic()
        m._LAST_LINE[0] = now
        print(f"[{now - _CLOCK['t0']:8.1f}s] {msg}", file=sys.stderr, flush=True)
    m.beat("Stage 3 build: " + msg.split(";")[0].split(":")[0][:60], "")


def check_time() -> None:
    if time.monotonic() > _CLOCK["deadline"]:
        raise TimeLimit()


# ----------------------------------------------------------------------------------------------------------------------
# File operations in the output area (v1.0.1; Windows-safe). Windows refuses to delete or rename a file while any
# handle or memory map is open on it (Linux allows both), and cloud sync, antivirus or the search indexer may hold a new
# file for a moment. Every such operation retries briefly on a transient sharing/access error and, if it still fails,
# raises FileOpError naming the operation, errno/WinError and the path relative to the output folder's parent.
# ----------------------------------------------------------------------------------------------------------------------
FS_RETRY_DELAYS: tuple = (0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 15.0, 30.0)   # seconds; about 61 s in total
TRANSIENT_WINERRORS = {5, 32, 33}   # access denied, sharing violation, lock violation
_FS_BASE: dict = {"dir": None}       # paths are reported relative to this folder (the output folder's parent)


def rel_path(p) -> str:
    """Path relative to the output folder's parent (never an absolute path); outside it: '<outside>/' + name."""
    p = Path(p)
    base = _FS_BASE["dir"]
    if base is not None:
        for a, b in ((p, Path(base)), (p.absolute(), Path(base).absolute())):
            try:
                return str(a.relative_to(b))
            except ValueError:
                pass
    return str(Path("<outside>") / p.name)


class FileOpError(Exception):
    """Exit 3: a file operation in the output area failed (after retries for transient errors)."""

    def __init__(self, op: str, paths, exc: BaseException | None = None, attempts: int = 1, seconds: float = 0.0,
                 reason: str | None = None):
        self.op, self.paths, self.exc, self.attempts, self.seconds = op, [rel_path(x) for x in paths], exc, attempts, seconds
        self.errno = getattr(exc, "errno", None)
        self.winerror = getattr(exc, "winerror", None)
        self.reason = reason
        super().__init__(self.describe())

    def describe(self) -> str:
        kind = type(self.exc).__name__ if self.exc is not None else "error"
        bits = [f"{kind} during {self.op}: {' -> '.join(self.paths)}"]
        if self.reason:
            bits.append(self.reason)
        if self.errno is not None:
            bits.append(f"errno {self.errno}")
        if self.winerror is not None:
            bits.append(f"WinError {self.winerror}")
        bits.append(f"{self.attempts} attempt{'s' if self.attempts != 1 else ''} over {self.seconds:.1f} s")
        return "; ".join(bits)


def is_transient(e: BaseException) -> bool:
    return isinstance(e, PermissionError) or getattr(e, "winerror", None) in TRANSIENT_WINERRORS


def fs_op(op: str, fn, *paths):
    """Run fn() (a file operation on paths); retry transient sharing/access errors with backoff; else FileOpError."""
    t0 = time.monotonic()
    delays = tuple(FS_RETRY_DELAYS)
    for attempt in range(1, len(delays) + 2):
        try:
            return fn()
        except FileOpError:
            raise
        except OSError as e:
            delay = delays[attempt - 1] if attempt <= len(delays) else None
            if not is_transient(e) or delay is None or time.monotonic() + delay > _CLOCK["deadline"]:
                raise FileOpError(op, paths, e, attempt, time.monotonic() - t0) from e
            wn = f", WinError {e.winerror}" if getattr(e, "winerror", None) is not None else ""
            log(f"file operation retry: {op} {rel_path(paths[0])} ({type(e).__name__}{wn}); attempt {attempt}, "
                f"waiting {delay:g} s")
            gc.collect()
            time.sleep(delay)


def remove_tree(folder: Path, op: str = "remove folder") -> None:
    """shutil.rmtree with per-entry retries (fs_op) for transient Windows sharing violations."""
    def retry(func, path):
        if func in (os.unlink, os.remove, os.rmdir):
            fs_op(f"{op} ({func.__name__})", lambda: func(path), Path(path))
            return True
        return False

    def onexc(func, path, exc):
        if not retry(func, path):
            raise exc

    def onerror(func, path, exc_info):
        if not retry(func, path):
            raise exc_info[1]

    kw = {"onexc": onexc} if sys.version_info >= (3, 12) else {"onerror": onerror}
    try:
        shutil.rmtree(folder, **kw)
    except FileOpError:
        raise
    except OSError as e:
        raise FileOpError(op, [folder], e) from e


def release_frames(e: BaseException | None) -> None:
    """Clear the locals of the finished frames of a failed build (and its chained exceptions) so that every memory map
    they hold is closed before the pending folder is removed (Windows cannot delete a mapped file)."""
    seen = set()
    while e is not None and id(e) not in seen:
        seen.add(id(e))
        if e.__traceback__ is not None:
            traceback.clear_frames(e.__traceback__)
        e = e.__cause__ or e.__context__
    gc.collect()


_WHERE_SKIP = {"fs_op", "remove_tree", "retry", "onexc", "onerror", "<lambda>", "zip_members", "w", "_fsync_path"}


def where(e: BaseException) -> str:
    """The two innermost lines of this script in the traceback, file helpers skipped ('f line N > g line M')."""
    locs = []
    tb = e.__traceback__
    me = os.path.abspath(__file__)
    while tb is not None:
        co = tb.tb_frame.f_code
        if os.path.abspath(co.co_filename) == me and co.co_name not in _WHERE_SKIP:
            locs.append(f"{co.co_name} line {tb.tb_lineno}")
        tb = tb.tb_next
    return " > ".join(locs[-2:])


# ----------------------------------------------------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------------------------------------------------
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def is_full_hash(s: str) -> bool:
    return bool(re.fullmatch(r"[0-9a-f]{64}", str(s or "")))


def check_name(p: Path, role: str) -> None:
    low = p.name.lower()
    if "_raw_local" in low:
        raise InputError(f"{role} '{p.name}' is a _raw_LOCAL file: never used as input")
    for bad in NEVER_READ:
        if low.startswith(bad.lower()):
            raise InputError(f"{role} '{p.name}' is never read by Stage 3")


def gate_file(role: str, path: Path, expected: str) -> str:
    path = Path(path)
    check_name(path, role)
    if not path.is_file():
        raise InputError(f"{role} not found: {path}")
    if not is_full_hash(expected):
        raise InputError(f"pinned {role} hash is not a full 64-character SHA-256")
    got = sha256_file(path)
    if got in REFUSED_SHA256:
        raise InputError(f"{role}: refused input ({REFUSED_SHA256[got]}); use channel map v2.3")
    if got != expected:
        raise InputError(f"{role} SHA-256 mismatch (expected {expected}, got {got})")
    return got


def vnorm_id(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.strip("'\"").str.strip()
    s = s.where(~s.str.lower().isin(_NULLS), "")
    m = s.str.fullmatch(r"[+-]?\d+\.0*")
    if m.any():
        u = {v: str(int(v.split(".")[0])) for v in pd.unique(s[m])}
        s = s.where(~m, s.map(u))
    return s


def vclean(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
    return s.where(~s.str.lower().isin(_NULLS), "")


def fast_time(s: pd.Series) -> pd.Series:
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


def has_clock(s: pd.Series) -> pd.Series:
    return s.astype(str).str.contains(r"\d{1,2}:\d{2}", regex=True)


def norm_unit(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.lower().str.replace("μ", "u").str.replace("µ", "u").str.replace(" ", "")
    return s.where(~s.isin(_NULLS), "")


def resolve(header: list[str], candidates: Iterable[str], label: str, required: bool = True) -> str | None:
    up = {h.upper(): h for h in header}
    for c in candidates:
        if c.upper() in up:
            return up[c.upper()]
    if required:
        raise InputError(f"{label}: none of the columns {list(candidates)} found")
    return None


def hours_between(t: pd.Series | np.ndarray, t0: pd.Series | np.ndarray) -> np.ndarray:
    a = np.asarray(pd.Series(t).to_numpy(dtype="datetime64[ns]"), dtype="datetime64[ns]").astype("int64")
    b = np.asarray(pd.Series(t0).to_numpy(dtype="datetime64[ns]"), dtype="datetime64[ns]").astype("int64")
    out = (a - b) / 3.6e12
    bad = (a == np.iinfo(np.int64).min) | (b == np.iinfo(np.int64).min)
    out = out.astype(float)
    out[bad] = np.nan
    return out


def r1(x: float) -> float:
    return float(round(float(x), 1))


def read_header(path: Path) -> list[str]:
    return list(pd.read_csv(path, dtype=str, nrows=0, encoding="utf-8-sig").columns)


def read_csv_str(path: Path, usecols: list[str] | None = None) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False, usecols=usecols, encoding="utf-8-sig")


def load_pinned_module(path: Path, expected: str, tag: str):
    got = gate_file(tag, path, expected)
    name = f"_stage3_pinned_{tag}_{got[:12]}"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def age_band(days: np.ndarray) -> np.ndarray:
    days = np.asarray(days, dtype=float)
    y = days / 365.25
    band = np.select([days < 28, y < 1, y < 3, y < 6, y < 12, y < 18],
                     ["neonate", "infant", "toddler", "preschool", "school", "adolescent"], "ge18")
    return np.where(np.isfinite(days) & (days >= 0), band, "age_invalid")


def norm_icd(code: Any) -> tuple[str, str]:
    """Normalized code (no dots/spaces) and status ok / missing / invalid / multiple (first code used)."""
    s = str(code if code is not None else "").strip().upper()
    if s.lower() in _NULLS:
        return "", "missing"
    parts = [p for p in re.split(r"[,;/|，；、\s]+", s) if p]
    first = re.sub(r"[^A-Z0-9]", "", parts[0]) if parts else ""
    if not re.match(r"^[A-Z][0-9]{2}", first):
        return first, "invalid"
    return first, ("multiple" if len(parts) > 1 else "ok")


# ----------------------------------------------------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------------------------------------------------
@dataclass
class Config:
    raw_dir: Path
    canonical_path: Path
    harmonization_path: Path
    split_path: Path
    complexity_map_path: Path
    phi1_medians_path: Path
    maps_dir: Path
    v17_module_path: Path
    suppress_module_path: Path
    links_path: Path | None
    out_dir: Path
    expect_d_icd_sha256: str
    expect_svs_sha256: str
    pins: dict = field(default_factory=lambda: dict(PINS))
    raw_pins: dict = field(default_factory=lambda: dict(RAW_PINS))
    expect: Expect = field(default_factory=Expect)
    surgery_expected_missing: dict | None = None
    split_subject_col: str = "SUBJECT_ID"
    split_col: str = "split"
    harm_t0_col: str | None = None
    harm_out_col: str | None = None
    chunksize: int = 500_000
    max_seconds: float = 7200.0
    precheck_only: bool = False
    links_sha256: str | None = None
    # v1.0.2 cross-release inputs (spec §16). stage2_pin is the script constant; tests may point it at a fixture file.
    stage2_release_path: Path | None = None
    expect_stage2_sha256: str | None = None
    stage2_pin: str = STAGE2_RELEASE_SHA256
    drops_path: Path | None = None
    drops_sha256: str | None = None
    # v1.0.4: team-internal constants (Stage 2 cells known from outside the release); None = no constants (tests)
    stage2_constants_path: Path | None = None
    stage2_constants_sha256: str | None = None
    max_new_suppressions: int = 2500       # suppression module v1.5 convergence caps (Stage 2 D8 defaults)
    max_component_cells: int = 20000
    heartbeat_seconds: float = 30.0


def raw_path(cfg: Config, name: str) -> Path:
    return Path(cfg.raw_dir) / f"{name}.csv"


def map_path(cfg: Config, key: str) -> Path:
    return Path(cfg.maps_dir) / MAP_FILES[key]


def check_out_dir(cfg: Config) -> None:
    out = Path(cfg.out_dir)
    check_name(out, "output folder")
    for p in (out, out.parent):
        if "_raw_local" in p.name.lower():
            raise InputError("output folder name contains _raw_LOCAL")
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise InputError(f"output folder exists and is not empty: {out} (Stage 3 never overwrites; use a new folder)")
    for bad in ("layer_a_features.csv", "layer_b_features.csv", "layer_c_features.csv"):
        if (out / bad).exists():
            raise InputError(f"output folder holds {bad}; Stage 3 never writes next to the retired files")


def precheck_hashes(cfg: Config) -> dict:
    """Every hash gate, before any data file is read. Returns {role: sha256}."""
    h: dict[str, str] = {}
    for role, p in (("canonical", cfg.canonical_path), ("harmonization", cfg.harmonization_path),
                    ("split", cfg.split_path), ("complexity_map", cfg.complexity_map_path),
                    ("phi1_medians", cfg.phi1_medians_path)):
        h[role] = gate_file(role, Path(p), cfg.pins[role])
    for key in MAP_FILES:
        h[key] = gate_file(key, map_path(cfg, key), cfg.pins[key])
    h["v17_module"] = gate_file("v17_module", cfg.v17_module_path, cfg.pins["v17_module"])
    h["suppress_module"] = gate_file("suppress_module", cfg.suppress_module_path, cfg.pins["suppress_module"])
    for name, (sha, _rows) in cfg.raw_pins.items():
        h[f"raw:{name}"] = gate_file(f"raw {name}", raw_path(cfg, name), sha)
    for name, sha in (("D_ICD_DIAGNOSES", cfg.expect_d_icd_sha256), ("SURGERY_VITAL_SIGNS", cfg.expect_svs_sha256)):
        if not is_full_hash(str(sha or "").lower()):
            raise InputError(f"--expect-{'d-icd' if name.startswith('D_') else 'svs'}-sha256 must be the full "
                             f"64-character SHA-256 recorded at the precheck (Step 1)")
        h[f"raw:{name}"] = gate_file(f"raw {name}", raw_path(cfg, name), str(sha).lower())
    if cfg.links_path is not None:
        lp = Path(cfg.links_path)
        check_name(lp, "links file")
        if not lp.is_file():
            raise InputError(f"links file not found: {lp}")
        got = sha256_file(lp)
        if cfg.links_sha256 is not None and got != cfg.links_sha256:
            raise InputError("links file SHA-256 mismatch")
        h["links"] = got
    if cfg.drops_path is not None:
        dp = Path(cfg.drops_path)
        check_name(dp, "drops file")
        if not dp.is_file():
            raise InputError(f"drops file not found: {dp}")
        got = sha256_file(dp)
        if cfg.drops_sha256 is not None and got != cfg.drops_sha256:
            raise InputError("drops file SHA-256 mismatch")
        h["drops"] = got
    if cfg.stage2_constants_path is not None:
        cp = Path(cfg.stage2_constants_path)
        check_name(cp, "constants file")
        if not cp.is_file():
            raise InputError(f"Stage 2 constants file not found: {cp.name}")
        got = sha256_file(cp)
        if cfg.stage2_constants_sha256 is not None and got != cfg.stage2_constants_sha256:
            raise InputError("Stage 2 constants file SHA-256 mismatch (team-internal file; nothing printed)")
        h["stage2_constants"] = got      # kept out of every output (build() removes it; see read_constants)
    h["stage2_release"] = gate_stage2_release(cfg)
    log(f"hash gates: {len(h)} inputs OK (full 64-character SHA-256)")
    return h


def gate_stage2_release(cfg: Config) -> str:
    """v1.0.2: the Stage 2 release JSON (v1.0.4: v1.0.3) is required and hash-gated (fail closed, exit 2): the expected hash on
    the command line must be the full SHA-256, must equal the script's pin, and must equal the file's hash."""
    if cfg.stage2_release_path is None:
        raise InputError("the Stage 2 v1.0.3 release JSON is required (--stage2-release) for the joint disclosure audit")
    exp = str(cfg.expect_stage2_sha256 or "").strip().lower()
    if not is_full_hash(exp):
        raise InputError("--expect-stage2-sha256 must be the full 64-character SHA-256 of the Stage 2 v1.0.3 release JSON")
    if not is_full_hash(str(cfg.stage2_pin)) or exp != cfg.stage2_pin:
        raise InputError(f"--expect-stage2-sha256 is not the pinned {STAGE2_RELEASE_LABEL} hash")
    return gate_file("Stage 2 release JSON", Path(cfg.stage2_release_path), exp)


def count_rows(path: Path, chunksize: int) -> int:
    first = read_header(path)[0]
    n = 0
    for ch in pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False, usecols=[first],
                          encoding="utf-8-sig", chunksize=chunksize):
        n += len(ch)
    return n


# ----------------------------------------------------------------------------------------------------------------------
# Cohort (contract H1, H2, H4, H5)
# ----------------------------------------------------------------------------------------------------------------------
@dataclass
class Cohort:
    enc: pd.DataFrame        # canonical file order; enc_i = row position
    eh: pd.DataFrame         # key 'SUBJECT|HADM' -> enc_i, t0, OUT, W
    n: int
    split: np.ndarray
    out_h: np.ndarray
    W: np.ndarray
    nvalid: np.ndarray       # number of valid bins (mask_valid = 1 iff 6k < out_h)
    eligible: dict

    @property
    def mask_valid(self) -> np.ndarray:
        return (np.arange(N_BINS)[None, :] < self.nvalid[:, None]).astype(np.uint8)

    def train(self) -> np.ndarray:
        return self.split == "train"


def load_cohort(cfg: Config, A) -> tuple[Cohort, dict]:
    ex = cfg.expect
    ch = A.read_header(cfg.canonical_path)
    need = ["encounter_id", "SUBJECT_ID", "HADM_ID", "hadm_id_list"]
    miss = [c for c in need if c not in ch]
    if miss:
        raise InputError(f"canonical CSV missing columns: {miss}")
    died_c = resolve(ch, ("Died",), "canonical Died", required=False)
    canon = A.read_csv_str(cfg.canonical_path, usecols=need + ([died_c] if died_c else []))
    canon["encounter_id"] = canon["encounter_id"].map(lambda v: str(v).strip())
    if (canon["encounter_id"] == "").any() or canon["encounter_id"].duplicated().any():
        raise InputError("canonical encounter_id empty or not unique")
    canon["SUBJECT_ID"] = canon["SUBJECT_ID"].map(A.norm_id)
    canon["HADM_ID"] = canon["HADM_ID"].map(A.norm_id)
    canon["hadm_list"] = canon["hadm_id_list"].map(A.parse_hadm_list)
    empty = canon["hadm_list"].map(len) == 0
    canon.loc[empty, "hadm_list"] = pd.Series([[h] if h else [] for h in canon.loc[empty, "HADM_ID"]],
                                              index=canon.index[empty], dtype=object)
    if died_c:
        dv = canon[died_c].map(lambda v: str(v).strip().lower())
        if (~dv.isin({"0", "1", "0.0", "1.0", "true", "false"})).any():
            raise InputError("canonical Died has values that are not 0/1")
        canon["died"] = dv.isin({"1", "1.0", "true"})
    else:
        raise InputError("canonical file has no Died column (needed for the deaths re-assertion)")
    hh = A.read_header(cfg.harmonization_path)
    if "encounter_id" not in hh:
        raise InputError("harmonization file has no encounter_id column")
    try:
        t0c = A.resolve_column(hh, cfg.harm_t0_col, A.HARM_T0_CANDIDATES, "T0", A.HARM_T0_FALLBACK)
        outc = A.resolve_column(hh, cfg.harm_out_col, A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK)
    except A.AuditError as e:
        raise InputError(str(e)) from e
    fcu = resolve(hh, ("first_careunit", "FIRST_CAREUNIT"), "harmonization first_careunit", required=False)
    harm = A.read_csv_str(cfg.harmonization_path, usecols=list(dict.fromkeys(
        ["encounter_id", t0c, outc] + ([fcu] if fcu else []))))
    harm["encounter_id"] = harm["encounter_id"].map(lambda v: str(v).strip())
    if harm["encounter_id"].duplicated().any():
        raise InputError("harmonization encounter_id not unique")
    if set(harm["encounter_id"]) != set(canon["encounter_id"]):
        raise InputError("harmonization encounter_id set differs from the canonical set")
    hdf = pd.DataFrame({"encounter_id": harm["encounter_id"], "t0": fast_time(harm[t0c]),
                        "OUT": fast_time(harm[outc]),
                        "first_careunit_h": vclean(harm[fcu]) if fcu else ""})
    if hdf["t0"].isna().any() or hdf["OUT"].isna().any():
        raise InputError("harmonized t0 or OUT missing or unparseable for some encounters")
    if (hdf["OUT"] < hdf["t0"]).any():
        raise InputError("harmonized OUT earlier than t0 for some encounters")
    enc = canon.merge(hdf, on="encounter_id", how="left", validate="one_to_one")
    sh = A.read_header(cfg.split_path)
    for c in (cfg.split_subject_col, cfg.split_col):
        if c not in sh:
            raise InputError(f"split file missing column {c!r}")
    spl = A.read_csv_str(cfg.split_path, usecols=[cfg.split_subject_col, cfg.split_col])
    spl = pd.DataFrame({"SUBJECT_ID": spl[cfg.split_subject_col].map(A.norm_id),
                        "split": spl[cfg.split_col].map(lambda v: A.SPLIT_LABELS.get(str(v).strip().lower(), "?"))})
    if (spl["split"] == "?").any():
        raise InputError("split file has unknown split labels")
    spl = spl.drop_duplicates()
    if spl["SUBJECT_ID"].duplicated().any():
        raise InputError("a subject appears in more than one split")
    enc = enc.merge(spl, on="SUBJECT_ID", how="left")
    if enc["split"].isna().any():
        raise InputError("some cohort encounters have no split")
    n_enc, n_subj = len(enc), enc["SUBJECT_ID"].nunique()
    sc = {s: int((enc["split"] == s).sum()) for s in SPLITS}
    if n_enc != ex.encounters or n_subj != ex.subjects or sc != dict(ex.splits):
        raise InputError("cohort/split totals do not reproduce contract H1/H4 (counts withheld from the log)")
    if int(enc["died"].sum()) != ex.died:
        raise InputError("canonical Died total does not reproduce contract H1")
    enc = enc.reset_index(drop=True)
    enc["enc_i"] = np.arange(len(enc))
    out_h = hours_between(enc["OUT"], enc["t0"])
    elig = {L: int((out_h > L).sum()) for L in LANDMARKS}
    if elig != {int(k): int(v) for k, v in ex.eligible.items()}:
        raise InputError("landmark eligible counts do not reproduce contract H5 (counts withheld from the log)")
    enc["out_h"] = out_h
    enc["W"] = np.minimum(out_h, MAX_H)
    rows = [(f"{s}|{h}", i) for s, lst, i in zip(enc["SUBJECT_ID"], enc["hadm_list"], enc["enc_i"]) for h in lst]
    eh = pd.DataFrame(rows, columns=["key", "enc_i"]).drop_duplicates()
    eh = eh.merge(enc[["enc_i", "t0", "OUT", "split", "W"]], on="enc_i", how="left")
    nvalid = np.minimum(N_BINS, np.ceil(out_h / BIN_H - 1e-12)).astype(int)
    co = Cohort(enc=enc, eh=eh, n=n_enc, split=enc["split"].to_numpy(), out_h=out_h, W=enc["W"].to_numpy(),
                nvalid=nvalid, eligible=elig)
    log(f"cohort re-assertion OK: encounters, subjects, deaths, splits and the 14 landmark eligible counts reproduce")
    return co, {"harmonization_columns": {"t0": t0c, "OUT": outc, "first_careunit": fcu or ""}}


# ----------------------------------------------------------------------------------------------------------------------
# Linking (contract H3; plan I13)
# ----------------------------------------------------------------------------------------------------------------------
ROW_CLASSES = ("in_stay", "cohort_admission_outside_stay", "no_cohort_admission", "cohort_admission_time_missing")


def link_chunk(subj: pd.Series, hadm: pd.Series, t: pd.Series, eh: pd.DataFrame, offset_h: float, keep_lo: float,
               classes: dict, keep_all: bool = False) -> pd.DataFrame:
    """SUBJECT_ID and HADM_ID in hadm_id_list, then time. Classes use the availability time (t + offset) against
    [t0, OUT). Kept rows: chart time th in [keep_lo, W) (or all linked timed rows with keep_all). One encounter per
    row: the stay that contains it, else the encounter it precedes most closely, else the latest."""
    n = len(subj)
    df = pd.DataFrame({"r": np.arange(n), "key": (subj + "|" + hadm).to_numpy(), "t": t.to_numpy()})
    ok = ((subj != "") & (hadm != "")).to_numpy()
    m = df[ok].merge(eh, on="key", how="inner")
    cohort_r = pd.unique(m["r"])
    classes["no_cohort_admission"] += n - len(cohort_r)
    tm = m["t"].isna()
    tmiss_r = pd.unique(m.loc[tm, "r"])
    classes["cohort_admission_time_missing"] += len(tmiss_r)
    mv = m[~tm].copy()
    mv["th"] = hours_between(mv["t"], mv["t0"])
    ta = mv["th"].to_numpy() + offset_h
    mv["stay"] = (ta >= 0) & (ta < hours_between(mv["OUT"], mv["t0"]))
    st = mv[mv["stay"]].drop_duplicates("r")
    classes["in_stay"] += len(st)
    classes["cohort_admission_outside_stay"] += len(cohort_r) - len(tmiss_r) - len(st)
    keep = mv if keep_all else mv[(mv["th"] >= keep_lo) & (mv["th"] < mv["W"])]
    keep = keep.assign(_pre=keep["th"] < 0, _abs=keep["th"].abs())
    keep = keep.sort_values(["r", "stay", "_pre", "_abs"], ascending=[True, False, False, True]).drop_duplicates("r")
    return keep[["r", "enc_i", "th", "stay"]].reset_index(drop=True)


def read_chunks(cfg: Config, name: str, cols: list[str]):
    return pd.read_csv(raw_path(cfg, name), dtype=str, keep_default_na=False, na_filter=False, usecols=cols,
                       encoding="utf-8-sig", chunksize=cfg.chunksize)


def finish_rows(cfg: Config, name: str, classes: dict, rows: int) -> None:
    if sum(classes.values()) != rows:
        raise GateError(f"{name}: row classes do not add up to the rows read")
    exp = cfg.raw_pins.get(name, (None, None))[1]
    if exp is not None and rows != exp:
        raise InputError(f"{name}: row count differs from contract Table P (count withheld)")


# ----------------------------------------------------------------------------------------------------------------------
# Maps
# ----------------------------------------------------------------------------------------------------------------------
CHANNEL_MAP_COLUMNS = ["table", "itemid", "label", "category", "units", "train_encounters", "train_rows", "channel",
                       "table51_system", "triggered", "contract16", "unit_conversion", "plausibility_low",
                       "plausibility_high", "source_itemid", "merge_rule", "stage2_check", "status", "confidence",
                       "decision", "notes", "mask_group"]


def read_map(cfg: Config, A, key: str) -> pd.DataFrame:
    p = map_path(cfg, key)
    d = A.read_csv_str(p)
    d = d.apply(lambda s: s.map(lambda v: str(v).strip()))
    if len(d) != MAP_ROWS[key]:
        raise InputError(f"map {p.name}: row count differs from its pin")
    if key == "channel_map" and list(d.columns) != CHANNEL_MAP_COLUMNS:
        raise InputError("channel map columns differ from the v2.3 schema")
    return d


@dataclass
class Channel:
    name: str
    system: str
    items: list
    table: str                    # CHARTEVENTS / LABEVENTS / DERIVED
    triggered: bool
    indicator_only: bool
    qualitative: bool
    family: str                   # no_locf_value / locf_4 / locf_1 / triggered / indicator_only
    cap_bins: int
    cap_h: float
    main_item: str | None
    dedupe_same_time: bool
    units: str
    value_class: bool


def parse_channel_map(cm: pd.DataFrame, ex: Expect) -> dict:
    """Channels of map v2.3: rows with a used status, a channel and a Table 5.1 system (94; checked with the Table 5.1
    system counts before any exclusion). v1.0.3: the EXCLUDED_CHANNELS (ketones_urine) are then removed with their
    items, and weight 1014 is the normalizer only (92 feature channels). Unknown conversion texts, unparseable gates or
    a family gap stop the run."""
    use = cm[cm["status"].isin(MAP_STATUS_USED) & (cm["channel"] != "") & (cm["table51_system"] != "")]
    items: dict[tuple[str, str], dict] = {}
    chans: dict[str, dict] = {}
    derived: dict[str, dict] = {}
    for _, r in use.iterrows():
        ch, tab, it = r["channel"], r["table"], r["itemid"]
        lo = float(r["plausibility_low"]) if r["plausibility_low"] else None
        hi = float(r["plausibility_high"]) if r["plausibility_high"] else None
        conv = r["unit_conversion"]
        if tab == "DERIVED":
            if conv != DERIVED_TEXT.get(ch):
                raise InputError(f"derived channel {ch}: definition text differs from the spec")
            derived[ch] = {"low": lo, "high": hi}
        elif it in AGGREGATION_TEXT:
            if conv != AGGREGATION_TEXT[it]:
                raise InputError(f"item {it}: aggregation text differs from the spec")
        elif conv and conv not in CONVERSIONS:
            raise InputError(f"item {it}: unrecognized unit_conversion text (stop; spec §6.3)")
        info = {"table": tab, "itemid": it, "channel": ch, "low": lo, "high": hi, "conv": conv,
                "status": r["status"], "mask_group": r["mask_group"], "triggered": r["triggered"] == "triggered",
                "merge_rule": r["merge_rule"]}
        if tab != "DERIVED":
            items[(tab, it)] = info
        c = chans.setdefault(ch, {"items": [], "tables": set(), "system": r["table51_system"], "triggered": False,
                                  "indicator_only": False, "statuses": set(), "merge_rules": set(), "units": r["units"]})
        if c["system"] != r["table51_system"]:
            raise InputError(f"channel {ch} has two Table 5.1 systems")
        c["items"].append(it)
        c["tables"].add(tab)
        c["triggered"] |= info["triggered"]
        c["indicator_only"] |= r["status"] == "assigned_merged_indicator_only"
        c["statuses"].add(r["status"])
        if r["merge_rule"]:
            c["merge_rules"].add(r["merge_rule"])
    if len(chans) != ex.channels:
        raise InputError(f"channel map v2.3 defines {len(chans)} channels, expected {ex.channels}")
    sysc: dict[str, int] = {}
    for c in chans.values():
        sysc[c["system"]] = sysc.get(c["system"], 0) + 1
    if sysc != dict(ex.systems):
        raise InputError("channel counts per Table 5.1 system differ from plan v13 Table 5.1")
    for ch in EXCLUDED_CHANNELS:                 # v1.0.3: dropped after the map checks (fail closed if absent)
        if ch not in chans:
            raise InputError(f"excluded channel {ch} is not a channel of map v2.3")
        del chans[ch]
        items = {k: i for k, i in items.items() if i["channel"] != ch}
        derived.pop(ch, None)
    bg = {it for (t, it), i in items.items() if i["mask_group"] == "blood_gas"}
    if bg != set(ex.blood_gas_items):
        raise InputError("blood-gas mask group differs from the 15 items of map v2.3")
    out: dict[str, Channel] = {}
    for ch, c in chans.items():
        tabs = sorted(c["tables"])
        if len(tabs) != 1:
            raise InputError(f"channel {ch} mixes tables")
        main = None
        dedupe = False
        for mr in c["merge_rules"]:
            m = re.findall(r"main item (\d+)", mr)
            if m:
                if main not in (None, m[0]):
                    raise InputError(f"channel {ch}: conflicting main items")
                main = m[0]
            if "deduplicate by draw" in mr:
                dedupe = True
        if ch in MAIN_ITEM_RULE:
            main, dedupe = MAIN_ITEM_RULE[ch], True
        if ch in DEDUPE_BY_DRAW_EXTRA:             # decision 5 (Nov 17, 2026): same rule as pCO2 / pH
            if main is None:
                raise InputError(f"channel {ch}: de-duplication by draw needs a main item in the map")
            dedupe = True
        if main is not None and main not in c["items"]:
            raise InputError(f"channel {ch}: main item not among its items")
        trig, ind = c["triggered"], c["indicator_only"]
        if trig:
            fam, cap_b, cap_h = "triggered", 0, 24.0
        elif ind:
            fam, cap_b, cap_h = "indicator_only", 0, 24.0
        elif ch in NO_LOCF_VALUE:
            fam, cap_b, cap_h = "no_locf_value", 0, 24.0
        elif ch in LOCF_4:
            fam, cap_b, cap_h = "locf_4", 4, 24.0
        elif ch in LOCF_1:
            fam, cap_b, cap_h = "locf_1", 1, 6.0
        elif ch in NORMALIZER_ONLY:
            fam, cap_b, cap_h = "normalizer_only", 0, 24.0
        else:
            raise InputError(f"channel {ch} has no carry-forward family (spec §7)")
        out[ch] = Channel(name=ch, system=c["system"], items=list(c["items"]), table=tabs[0], triggered=trig,
                          indicator_only=ind, qualitative=ch in QUALITATIVE, family=fam, cap_bins=cap_b, cap_h=cap_h,
                          main_item=main, dedupe_same_time=dedupe, units=c["units"],
                          value_class=fam in ("no_locf_value", "locf_4", "locf_1"))
    for need in ("map", "kdigo_uo_stage"):
        g = derived.get(need)
        if not g or g["low"] is None or g["high"] is None:
            raise InputError(f"derived channel {need} has no gate in the map")
    feats = [ch for ch in out if ch not in NORMALIZER_ONLY]
    if len(feats) != ex.feature_channels:
        raise InputError(f"{len(feats)} feature channels, expected {ex.feature_channels} (94 map channels minus the "
                         f"weight normalizer and {len(EXCLUDED_CHANNELS)} excluded)")
    value_class = [ch for ch in feats if out[ch].value_class]
    return {"items": items, "channels": out, "derived": derived, "feature_channels": feats,
            "value_channels": value_class, "blood_gas_items": sorted(bg), "map_channels": ex.channels,
            "excluded_channels": dict(EXCLUDED_CHANNELS),
            "panels": {g: sorted({it for (t, it), i in items.items() if i["mask_group"] == g})
                       for g in ("blood_gas", "cytokines", "teg", "urine_dipstick")}}


def parse_drug_map(dm: pd.DataFrame) -> dict:
    out = {}
    for _, r in dm.iterrows():
        g = r["group"]
        if g not in DRUG_MAP_GROUPS:
            raise InputError(f"drug map: unknown group {g!r}")
        other = None
        if g == "osmotic_conditional":
            m = re.search(r"otherwise\s*->\s*([a-z_]+)", r["conditional_rule"])
            if (not m or m.group(1) not in MED_GROUPS or "neuro admitting dx -> osmotherapy" not in r["conditional_rule"]
                    or "sensitivity: all -> osmotherapy" not in r["conditional_rule"]):
                raise InputError("drug map: osmotic rule text differs from D3")
            other = m.group(1)
        key = f"{r['drug_name']} | {r['drug_name_en']}"
        if key in out:
            raise InputError("drug map: duplicate key")
        out[key] = {"group": g, "otherwise": other}
    return out


def parse_specimen_map(sm: pd.DataFrame) -> dict:
    out = {}
    for _, r in sm.iterrows():
        if r["specimen_group"] not in SPECIMEN_GROUPS:
            raise InputError(f"specimen map: unknown group {r['specimen_group']!r}")
        if r["base_code"] != r["spec_itemid"].split("+")[0]:
            raise InputError("specimen map: base_code is not the first code of spec_itemid")
        out[f"{r['spec_itemid']} | {r['spec_type_desc']}"] = {"group": r["specimen_group"], "base": r["base_code"]}
    return out


def neuro_matcher(neuro: pd.DataFrame):
    three, four = set(), set()
    for _, r in neuro.iterrows():
        if r["status"] not in ("include", "partial"):
            continue
        rule, c3 = r["match_rule"].strip(), r["icd10_3char"].strip().upper()
        if rule == "3-character prefix":
            three.add(c3)
        elif rule.startswith("4-character prefix"):
            codes = re.findall(r"([A-Z][0-9]{2})\.([0-9])", rule)
            if not codes:
                raise InputError("neuro set: cannot parse a 4-character rule")
            four |= {a + b for a, b in codes}
        else:
            raise InputError("neuro set: unknown match_rule")

    def match(code: str) -> bool:
        if not code or not re.match(r"^[A-Z][0-9]{2}", code):
            return False
        return code[:3] in three or any(code.startswith(f) for f in four)
    return match


def icdo_first_code(code: Any) -> str:
    """First code of a diagnosis field in the form the ICD-O rule reads (v1.0.3): NFKC (full-width -> ASCII), upper
    case, spaces around '/', '.', '-' and after a leading 'M' removed, then cut at the first code separator (',', ';',
    '|', '、' or a space; '/' is NOT a separator here because it introduces the ICD-O behaviour digit)."""
    s = unicodedata.normalize("NFKC", str(code if code is not None else "")).strip().upper()
    s = re.sub(r"\s*([/.\-])\s*", r"\1", s)
    s = re.sub(r"^M\s+(?=[0-9])", "M", s)
    parts = [p for p in ICDO_CODE_SPLIT.split(s) if p]
    return parts[0] if parts else ""


def icdo_parts(code: Any) -> dict | None:
    """Parsed ICD-O-3 morphology of the first code, or None: {'histology': '9800', 'extension': '0' or '',
    'behaviour': '3' or ''} (behaviour from '/b' or a 6th digit; a lone 5th digit is reported as the extension digit
    because PIC's clinical-edition form M98000 carries the behaviour after '/')."""
    m = ICDO_MORPH_RE.fullmatch(icdo_first_code(code))
    if not m:
        return None
    return {"histology": m.group("histology"), "extension": m.group("ext") or "",
            "behaviour": m.group("behaviour") or m.group("beh6") or ""}


def is_icdo_morphology(code: Any) -> bool:
    """ICD-O-3 morphology code (decision 8; v1.0.3 rule, see ICDO_MORPH_RE): the first code of the field is 'M' + 4
    digits in 8000-9989 with an optional extension digit and an optional behaviour digit (M98000, M98000/3, M9800/3,
    M9800, m 9800 / 3, M980003). Used only for codes D_ICD_DIAGNOSES does not map, so ICD-10 chapter M codes in the
    dictionary keep their ICD-10 meaning; dotted ICD-10 forms (M80.000) and histology below 8000 never match."""
    return icdo_parts(code) is not None


def icd_domain(c3: str) -> str:
    """Draft Appendix C diagnosis-only domain (plan v13; draft defaults, provisional)."""
    if not re.fullmatch(r"[A-Z][0-9]{2}", c3 or ""):
        return "other_medical"
    L, n = c3[0], int(c3[1:])
    if L == "Q" and 20 <= n <= 28:
        return "cardiac_congenital_heart"
    if L == "P":
        return "neonatal_perinatal"
    if L in ("S", "T"):
        return "surgical_trauma"
    if L in ("J", "A", "B"):
        return "respiratory_infection"
    if L == "G":
        return "neurological"
    if L == "C" or (L == "D" and 50 <= n <= 89):
        return "oncology_hematology"
    if L in ("K", "N"):
        return "gi_renal"
    return "other_medical"


def read_d_icd(cfg: Config, A) -> dict:
    """D_ICD_DIAGNOSES: ICD10_CODE_CN -> ICD10_CODE (the hash was gated at the precheck)."""
    p = raw_path(cfg, "D_ICD_DIAGNOSES")
    hdr = A.read_header(p)
    cn = resolve(hdr, ("ICD10_CODE_CN",), "D_ICD_DIAGNOSES")
    en = resolve(hdr, ("ICD10_CODE",), "D_ICD_DIAGNOSES")
    d = A.read_csv_str(p, usecols=[cn, en])
    out: dict[str, str] = {}
    for a, b in zip(d[cn], d[en]):
        k = re.sub(r"\s+", "", str(a).strip().upper())
        v = re.sub(r"[^A-Z0-9]", "", str(b).strip().upper())
        if k and v and k not in out:
            out[k] = v
    return out


def read_phi1_medians(cfg: Config, threshold: str = "10") -> dict:
    """Frozen φ1 medians (csv module: keys may hold commas inside quotes)."""
    with open(cfg.phi1_medians_path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows or set(rows[0]) != {"threshold", "key", "level", "train_n", "median_hours"}:
        raise InputError("φ1 medians file has an unexpected header")
    out = {}
    for r in rows:
        if r["threshold"].strip() != threshold:
            continue
        try:
            out[r["key"].strip()] = float(r["median_hours"])
        except ValueError:
            raise InputError("φ1 medians file has a non-numeric median")
    if "all" not in out:
        raise InputError("φ1 medians file has no 'all' median")
    return out


# ----------------------------------------------------------------------------------------------------------------------
# Static tables: ICUSTAYS, PATIENTS, ADMISSIONS
# ----------------------------------------------------------------------------------------------------------------------
def read_icustays(cfg: Config, A, co: Cohort) -> dict:
    name = "ICUSTAYS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    cols = [resolve(hdr, (c,), name) for c in ("SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "INTIME", "OUTTIME")]
    fcu = resolve(hdr, ("FIRST_CAREUNIT",), name)
    d = A.read_csv_str(p, usecols=cols + [fcu])
    finish_rows(cfg, name, {"in_stay": len(d)}, len(d))
    st = pd.DataFrame({"key": vnorm_id(d[cols[0]]) + "|" + vnorm_id(d[cols[1]]), "ICUSTAY_ID": vnorm_id(d[cols[2]]),
                       "IN": fast_time(d[cols[3]]), "OUTS": fast_time(d[cols[4]]), "unit": vclean(d[fcu])})
    m = st.merge(co.eh[["key", "enc_i", "t0", "OUT"]], on="key", how="inner")
    tol = pd.Timedelta(minutes=1)
    m = m[(m["IN"] >= m["t0"] - tol) & (m["OUTS"] <= m["OUT"] + tol)].sort_values(["enc_i", "IN"], kind="mergesort")
    m = m.drop_duplicates(["enc_i", "ICUSTAY_ID"])
    first = m.drop_duplicates("enc_i", keep="first").set_index("enc_i")["unit"]
    stays = {int(i): set(g["ICUSTAY_ID"]) for i, g in m.groupby("enc_i")}
    fu = np.array([first.get(i, "") for i in range(co.n)], dtype=object)
    h = co.enc["first_careunit_h"].to_numpy(dtype=object) if "first_careunit_h" in co.enc else np.full(co.n, "")
    fu = np.where(h != "", h, fu)
    return {"stays": stays, "first_careunit": fu, "n_without_stay": int(co.n - m["enc_i"].nunique())}


def classify_ages(t0: pd.Series, dob: pd.Series, dob_raw: pd.Series, has_row: np.ndarray):
    raw = vclean(dob_raw.fillna("")).to_numpy()
    missing = (~np.asarray(has_row, dtype=bool)) | (raw == "")
    unparse = (~missing) & dob.isna().to_numpy()
    age = hours_between(pd.to_datetime(t0), dob) / 24.0
    fin = np.isfinite(age)
    unparse = unparse | ((~missing) & (~unparse) & (~fin))
    neg = (~missing) & (~unparse) & fin & (age < 0)
    valid = (~missing) & (~unparse) & (~neg)
    status = np.select([missing, unparse, neg], ["missing", "unparseable", "negative"], "valid")
    mag = np.where(neg, -age, np.nan)
    nbin = np.select([neg & (mag < 1.0), neg & (mag <= 30.0), neg], ["lt1d", "1to30d", "gt30d"], "")
    return status, np.where(valid, age, np.nan), nbin, np.where(neg & (mag < 1.0), 0.0, np.where(valid, age, np.nan))


def read_patients(cfg: Config, A, co: Cohort) -> pd.DataFrame:
    name = "PATIENTS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, gc, dc = (resolve(hdr, (c,), name) for c in ("SUBJECT_ID", "GENDER", "DOB"))
    d = A.read_csv_str(p, usecols=[sc, gc, dc])
    finish_rows(cfg, name, {"in_stay": len(d)}, len(d))
    pt = pd.DataFrame({"SUBJECT_ID": vnorm_id(d[sc]), "sex": vclean(d[gc]).str.upper(), "DOB_raw": d[dc],
                       "DOB": fast_time(d[dc])})
    if pt["SUBJECT_ID"].duplicated().any():
        raise InputError("PATIENTS SUBJECT_ID not unique")
    m = co.enc[["enc_i", "SUBJECT_ID", "t0"]].merge(pt, on="SUBJECT_ID", how="left", indicator=True)
    m = m.sort_values("enc_i", kind="mergesort").reset_index(drop=True)
    if len(m) != co.n:
        raise InputError("PATIENTS merge changed the encounter count")
    status, age_days, nbin, age_sens = classify_ages(m["t0"], m["DOB"], m["DOB_raw"], m["_merge"].eq("both").to_numpy())
    n_inv = int((status != "valid").sum())
    if n_inv > int(cfg.expect.max_invalid_age_frac * co.n):
        raise InputError("invalid DOB/age above the 1% safety gate (D20; count withheld)")
    sx = m["sex"].fillna("").to_numpy()
    sex = np.where(np.isin(sx, ["M", "F"]), sx, "unrecognized")
    return pd.DataFrame({"sex": sex, "age_days": age_days, "age_status": status, "age_neg_bin": nbin,
                         "age_days_sens_lt1d0": age_sens})


def read_admissions(cfg: Config, A, co: Cohort) -> pd.DataFrame:
    """Earliest linked admission (min ADMITTIME over hadm_id_list): ADMITTIME, ICD10_CODE_CN, EDREGTIME."""
    name = "ADMISSIONS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc, at = (resolve(hdr, (c,), name) for c in ("SUBJECT_ID", "HADM_ID", "ADMITTIME"))
    dxc = resolve(hdr, ("ICD10_CODE_CN",), name)
    edc = resolve(hdr, ("EDREGTIME",), name, required=False)
    d = A.read_csv_str(p, usecols=[c for c in (sc, hc, at, dxc, edc) if c])
    finish_rows(cfg, name, {"in_stay": len(d)}, len(d))
    ad = pd.DataFrame({"key": vnorm_id(d[sc]) + "|" + vnorm_id(d[hc]), "ADMITTIME": fast_time(d[at]),
                       "dx": d[dxc].astype(str), "ED": fast_time(d[edc]) if edc else pd.NaT,
                       "ed_raw": vclean(d[edc]) if edc else ""})
    m = co.eh[["key", "enc_i"]].merge(ad, on="key", how="inner").sort_values(["enc_i", "ADMITTIME"],
                                                                         na_position="last", kind="mergesort")
    first = m.drop_duplicates("enc_i", keep="first").set_index("enc_i")
    out = pd.DataFrame({"enc_i": np.arange(co.n)})
    for c in ("ADMITTIME", "dx", "ED", "ed_raw"):
        out[c] = out["enc_i"].map(first[c])
    out["dx"] = out["dx"].fillna("")
    out["ed_raw"] = out["ed_raw"].fillna("")
    out["has_row"] = out["enc_i"].isin(first.index)
    out["edc_present"] = bool(edc)
    return out


# ----------------------------------------------------------------------------------------------------------------------
# Event tables: CHARTEVENTS / LABEVENTS (values), PRESCRIPTIONS, MICROBIOLOGYEVENTS
# ----------------------------------------------------------------------------------------------------------------------
KEEP_LO_EVENTS_H = -26.0     # chart times kept from t0 - 26 h (24 h pre-ICU sensitivity with the + 1 h lab release)


def per_row_units_needed(cmap: dict) -> set[tuple[str, str]]:
    out = set()
    for k, i in cmap["items"].items():
        conv = CONVERSIONS.get(i["conv"])
        if conv and conv[0] is not None:
            out.add(k)
    return out


def read_values(cfg: Config, A, co: Cohort, name: str, cmap: dict, stays: dict) -> tuple[pd.DataFrame, dict]:
    """All rows are classified (row classes); rows of mapped items with chart time in [t0 - 26 h, W) are kept.
    ICUSTAY_ID disagreements are counted on in-stay kept rows; the time-based assignment wins (I13)."""
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc = resolve(hdr, ("SUBJECT_ID",), name), resolve(hdr, ("HADM_ID",), name)
    ic, tc = resolve(hdr, ("ITEMID",), name), resolve(hdr, ("CHARTTIME",), name)
    vc = resolve(hdr, ("VALUE",), name, required=False)
    nc = resolve(hdr, ("VALUENUM",), name, required=False)
    uc = resolve(hdr, ("VALUEUOM",), name, required=False)
    icu = resolve(hdr, ("ICUSTAY_ID",), name, required=False)
    used = [c for c in (sc, hc, ic, tc, vc, nc, uc, icu) if c]
    wanted = {it for (t, it) in cmap["items"] if t == name}
    qual_items = {it for (t, it), i in cmap["items"].items() if t == name and i["channel"] in QUALITATIVE}
    unit_items = {it for (t, it) in per_row_units_needed(cmap) if t == name}
    off = LAB_RELEASE_H if name == "LABEVENTS" else 0.0
    classes = {c: 0 for c in ROW_CLASSES}
    icu_chk = {"n_rows_checked": 0, "n_icustay_blank": 0, "n_icustay_other_stay": 0, "n_icustay_agrees": 0}
    rows, parts = 0, []
    for raw in read_chunks(cfg, name, used):
        check_time()
        rows += len(raw)
        subj, hadm, t = vnorm_id(raw[sc]), vnorm_id(raw[hc]), fast_time(raw[tc])
        k = link_chunk(subj, hadm, t, co.eh, off, KEEP_LO_EVENTS_H, classes)
        if not len(k):
            continue
        r = k["r"].to_numpy()
        item = vnorm_id(raw[ic].iloc[r]).to_numpy()
        sel = np.isin(item, list(wanted))
        if not sel.any():
            continue
        r, e, th, stay, item = r[sel], k["enc_i"].to_numpy()[sel], k["th"].to_numpy()[sel], k["stay"].to_numpy()[sel], item[sel]
        if icu:
            iv = vnorm_id(raw[icu].iloc[r]).to_numpy()
            ins = np.flatnonzero(stay)
            icu_chk["n_rows_checked"] += len(ins)
            blank = iv[ins] == ""
            icu_chk["n_icustay_blank"] += int(blank.sum())
            agree = np.array([iv[j] in stays.get(int(e[j]), ()) for j in ins[~blank]], dtype=bool)
            icu_chk["n_icustay_agrees"] += int(agree.sum())
            icu_chk["n_icustay_other_stay"] += int((~agree).sum())
        num = pd.to_numeric(raw[nc].iloc[r], errors="coerce").to_numpy(dtype=float).copy() if nc else np.full(len(r), np.nan)
        txt = vclean(raw[vc].iloc[r]).to_numpy(dtype=object) if vc else np.full(len(r), "", dtype=object)
        if vc:
            miss = np.isnan(num)
            if miss.any():
                num[miss] = pd.to_numeric(pd.Series(txt[miss]), errors="coerce").to_numpy(dtype=float)
        isq = np.isin(item, list(qual_items))
        txt = np.where(isq, txt, "")
        unit = norm_unit(raw[uc].iloc[r]).to_numpy(dtype=object) if uc else np.full(len(r), "", dtype=object)
        unit = np.where(np.isin(item, list(unit_items)), unit, "")
        parts.append(pd.DataFrame({"enc_i": e.astype(np.int32), "tc": th, "t": th + off, "stay": stay, "item": item,
                                   "vnum": num, "vtext": txt, "unit": unit}))
    finish_rows(cfg, name, classes, rows)
    ev = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int32), "tc": pd.Series(dtype=float), "t": pd.Series(dtype=float),
         "stay": pd.Series(dtype=bool), "item": pd.Series(dtype=object), "vnum": pd.Series(dtype=float),
         "vtext": pd.Series(dtype=object), "unit": pd.Series(dtype=object)})
    ev["item"] = ev["item"].astype(str)
    log(f"{name}: rows read and classified; mapped rows in the kept window retained")
    return ev, {"row_classes": classes, "rows": rows, "icustay": icu_chk if icu else None}


def assign_date_only_start(th: np.ndarray, date_only: np.ndarray) -> np.ndarray:
    """Date-only STARTDATE (midnight): t0 on the t0 date, else the first bin start at or after midnight; dates
    before the t0 date keep their midnight (plan §5.3 (e))."""
    th = np.asarray(th, dtype=float).copy()
    do = np.asarray(date_only, dtype=bool)
    later = do & (th > 0)
    same = do & (th <= 0) & (th > -24)
    th[later] = np.ceil(th[later] / BIN_H - 1e-12) * BIN_H
    th[same] = 0.0
    return th


def read_prescriptions(cfg: Config, A, co: Cohort, dmap: dict, neuro_enc: np.ndarray) -> tuple[pd.DataFrame, dict]:
    name = "PRESCRIPTIONS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc = resolve(hdr, ("SUBJECT_ID",), name), resolve(hdr, ("HADM_ID",), name)
    stc = resolve(hdr, ("STARTDATE", "STARTTIME"), name)
    enc_ = resolve(hdr, ("ENDDATE", "ENDTIME"), name, required=False)
    dn = resolve(hdr, ("DRUG_NAME",), name, required=False)
    de = resolve(hdr, ("DRUG_NAME_EN",), name, required=False)
    if not dn and not de:
        raise InputError("PRESCRIPTIONS: neither DRUG_NAME nor DRUG_NAME_EN found")
    used = [c for c in (sc, hc, stc, enc_, dn, de) if c]
    classes = {c: 0 for c in ROW_CLASSES}
    rows, parts = 0, []
    t0_by = co.enc.set_index("enc_i")["t0"]
    for raw in read_chunks(cfg, name, used):
        check_time()
        rows += len(raw)
        subj, hadm = vnorm_id(raw[sc]), vnorm_id(raw[hc])
        s_raw = raw[stc]
        t = fast_time(s_raw)
        clock = has_clock(s_raw).to_numpy()
        k = link_chunk(subj, hadm, t, co.eh, 0.0, -1e9, classes, keep_all=True)
        if not len(k):
            continue
        r, e = k["r"].to_numpy(), k["enc_i"].to_numpy()
        th = assign_date_only_start(k["th"].to_numpy(), ~clock[r])
        keep = (th >= -PRE_ICU_SENS_H) & (th < co.W[e])
        if not keep.any():
            continue
        rk, ek = r[keep], e[keep]
        a = vclean(raw[dn].iloc[rk]) if dn else pd.Series("", index=raw.index[rk])
        b = vclean(raw[de].iloc[rk]) if de else pd.Series("", index=raw.index[rk])
        key = (a + " | " + b).to_numpy()
        if enc_:
            e_raw = raw[enc_].iloc[rk]
            et = fast_time(e_raw)
            eclock = has_clock(e_raw).to_numpy()
            et = et + pd.to_timedelta(np.where(eclock | et.isna().to_numpy(), 0, 24), unit="h").to_numpy()
            eh_ = hours_between(et, t0_by.reindex(ek))
        else:
            eh_ = np.full(len(rk), np.nan)
        parts.append(pd.DataFrame({"enc_i": ek.astype(np.int64), "S": th[keep], "E": eh_, "key": key}))
    finish_rows(cfg, name, classes, rows)
    docs = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int64), "S": pd.Series(dtype=float), "E": pd.Series(dtype=float),
         "key": pd.Series(dtype=object)})
    info = docs["key"].map(dmap)
    base = info.map(lambda d: d["group"] if isinstance(d, dict) else "unmapped")
    other = info.map(lambda d: d["otherwise"] if isinstance(d, dict) else None)
    osm = (base == "osmotic_conditional").to_numpy()
    neuro = neuro_enc[docs["enc_i"].to_numpy()] if len(docs) else np.zeros(0, dtype=bool)
    docs["group"] = np.where(osm, np.where(neuro, "osmotherapy", other.fillna("").to_numpy()), base.to_numpy())
    n_unmapped = int((base == "unmapped").sum())
    log(f"{name}: rows read and classified; documentation rows kept from t0 - 24 h")
    return docs, {"row_classes": classes, "rows": rows, "n_docs_unmapped_drug": n_unmapped}


def read_micro(cfg: Config, A, co: Cohort, smap: dict) -> tuple[pd.DataFrame, dict]:
    name = "MICROBIOLOGYEVENTS"
    p = raw_path(cfg, name)
    hdr = A.read_header(p)
    sc, hc = resolve(hdr, ("SUBJECT_ID",), name), resolve(hdr, ("HADM_ID",), name)
    tc = resolve(hdr, ("CHARTTIME",), name)
    si, sd = resolve(hdr, ("SPEC_ITEMID",), name, required=False), resolve(hdr, ("SPEC_TYPE_DESC",), name)
    classes = {c: 0 for c in ROW_CLASSES}
    rows, parts = 0, []
    for raw in read_chunks(cfg, name, [c for c in (sc, hc, tc, si, sd) if c]):
        check_time()
        rows += len(raw)
        subj, hadm, t = vnorm_id(raw[sc]), vnorm_id(raw[hc]), fast_time(raw[tc])
        k = link_chunk(subj, hadm, t, co.eh, 0.0, 0.0, classes)
        if not len(k):
            continue
        r = k["r"].to_numpy()
        sid = vnorm_id(raw[si].iloc[r]) if si else pd.Series("", index=raw.index[r])
        sdesc = vclean(raw[sd].iloc[r])
        parts.append(pd.DataFrame({"enc_i": k["enc_i"].to_numpy().astype(np.int64), "th": k["th"].to_numpy(),
                                   "spec_key": (sid + " | " + sdesc).to_numpy(), "spec_itemid": sid.to_numpy()}))
    finish_rows(cfg, name, classes, rows)
    ev = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int64), "th": pd.Series(dtype=float), "spec_key": pd.Series(dtype=object),
         "spec_itemid": pd.Series(dtype=object)})
    sinfo = ev["spec_key"].map(smap)
    mapped = sinfo.notna().to_numpy()
    ev["base"] = np.where(mapped, sinfo.map(lambda d: d["base"] if isinstance(d, dict) else ""),
                          ev["spec_itemid"].astype(str).str.split("+").str[0].to_numpy())
    ev["group"] = np.where(mapped, sinfo.map(lambda d: d["group"] if isinstance(d, dict) else ""), "unmapped")
    col = ev.drop_duplicates(["enc_i", "th", "base"]).reset_index(drop=True)
    log(f"{name}: rows read and classified; collections in [t0, W) retained")
    return col, {"row_classes": classes, "rows": rows, "n_rows_specimen_unmapped": int((~mapped).sum())}


# ----------------------------------------------------------------------------------------------------------------------
# Channel values: unit registry, gates, merges, derived channels (spec §6.3-§6.6)
# ----------------------------------------------------------------------------------------------------------------------
def qual_norm(txt: np.ndarray) -> pd.Series:
    """Token-table normal form of qualitative text: stripped, lower case, no spaces."""
    return pd.Series(txt, dtype=object).astype(str).str.strip().str.lower().str.replace(" ", "")


def qual_parse(txt: np.ndarray, channel: str | None = None) -> np.ndarray:
    """Positivity from the fixed token table (provisional, spec §6.5): 1 positive, 0 negative, NaN unparsed (decision
    4, Nov 17, 2026: unparsed text stays missing, flagged). v1.0.4: plus the per-channel additions
    (QUAL_POS_BY_CHANNEL, QUAL_NEG_BY_CHANNEL), read only on their own channel."""
    s = qual_norm(txt)
    neg = {t.replace(" ", "").lower() for t in QUAL_NEG | QUAL_NEG_BY_CHANNEL.get(channel, set())}
    pos = {t.replace(" ", "").lower() for t in QUAL_POS | QUAL_POS_BY_CHANNEL.get(channel, set())}
    return np.where(s.isin(pos), 1.0, np.where(s.isin(neg), 0.0, np.nan))


def qual_below_cutoff(txt: np.ndarray, channel: str | None) -> np.ndarray:
    """v1.0.4: rows read as negative from a below-cutoff result (QUAL_BELOW_CUTOFF) -> bool per row."""
    toks = {t.replace(" ", "").lower() for t in QUAL_BELOW_CUTOFF.get(channel, set())}
    if not toks:
        return np.zeros(len(txt), dtype=bool)
    return qual_norm(txt).isin(toks).to_numpy()


def classify_rows(info: dict, g: pd.DataFrame, qualitative: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Value in the final unit, the row status and the gate side (-1 below, +1 above, 0 otherwise) for one item (gate
    in the gated unit, then the final factor). Decision 2: a 0 of a ZERO_AS_MISSING channel is 'zero_missing'."""
    v = g["vnum"].to_numpy(dtype=float).copy()
    st = np.full(len(g), "measurement", dtype=object)
    side = np.zeros(len(g), dtype=np.int8)
    conv = CONVERSIONS.get(info["conv"]) if info["conv"] else None
    final = 1.0
    if conv is not None:
        per_row, final = conv
        if per_row is not None:
            u = g["unit"].astype(str).to_numpy()
            known = np.isin(u, list(per_row))
            f = np.array([per_row[x] if (x in per_row and per_row[x] is not None) else np.nan for x in u], dtype=float)
            st[~known] = "unit_unrecognized"
            st[known & np.isnan(f)] = "unit_dropped"
            v = v * f
    if qualitative:
        q = qual_parse(g["vtext"].to_numpy(), info.get("channel"))
        free = st == "measurement"
        st[free & np.isnan(q)] = "qual_unparsed"
        return q, st, side
    has = ~np.isnan(v)
    free = st == "measurement"
    if info.get("channel") in ZERO_AS_MISSING:
        st[free & has & (v == 0)] = "zero_missing"
        free = st == "measurement"
    if info["low"] is not None or info["high"] is not None:
        lo = -np.inf if info["low"] is None else info["low"]
        hi = np.inf if info["high"] is None else info["high"]
        inside = has & (v >= lo) & (v <= hi)
        st[free & ~has] = "nonnumeric"
        out = free & has & ~inside
        st[out] = "out_of_gate"
        side[out & (v < lo)] = -1
        side[out & (v > hi)] = 1
    else:
        st[free & ~has] = "nonnumeric"
    return v * final, st, side


@dataclass
class Values:
    obs: dict                 # channel -> DataFrame(enc_i, t, tc, item, v) (measurement rows; final units)
    rowstat: pd.DataFrame     # enc_i, t, channel, item, status (rows of mapped items, kept window)
    weight: np.ndarray
    weight_source: np.ndarray
    uo_info: dict
    weight_arms: "WeightArms | None" = None
    unparsed: "pd.DataFrame | None" = None   # enc_i, t, channel, text: qualitative rows outside the token table
    below_cutoff: "pd.DataFrame | None" = None   # v1.0.4: enc_i, t, channel: negatives read from a below-cutoff result


def _empty_obs() -> pd.DataFrame:
    return pd.DataFrame({"enc_i": pd.Series(dtype=np.int64), "t": pd.Series(dtype=float), "tc": pd.Series(dtype=float),
                         "item": pd.Series(dtype=object), "v": pd.Series(dtype=float)})


def parse_growth_ref(d: pd.DataFrame) -> dict:
    """Frozen weight-for-age reference (maps/growth_weight_for_age_ref_v1.csv, hash-pinned): {(reference, sex):
    (age_days, M)} with strictly increasing ages. Any schema deviation stops the run (exit 2)."""
    if list(d.columns) != GROWTH_REF_COLUMNS:
        raise InputError("growth reference columns differ from the v1 schema")
    if not (set(d["reference"]) <= set(GROWTH_REFS) and set(d["sex"]) <= {"M", "F"}):
        raise InputError("growth reference: unexpected reference or sex level")
    out = {}
    for (ref, sx), g in d.groupby(["reference", "sex"], sort=True):
        try:
            a = g["age_days"].astype(float).to_numpy()
            m = g["M"].astype(float).to_numpy()
        except ValueError:
            raise InputError("growth reference: non-numeric age or median") from None
        if len(a) < 2 or not (np.diff(a) > 0).all() or a[0] != 0.0 or not (np.isfinite(m) & (m > 0)).all():
            raise InputError("growth reference: ages must start at 0 and increase; medians must be positive")
        out[(ref, sx)] = (a, m)
    if set(out) != {(r, s) for r in GROWTH_REFS for s in ("M", "F")}:
        raise InputError("growth reference: every reference needs both sexes")
    return out


def growth_ref_median(gref: dict, ref: str, age_days: np.ndarray, sex: np.ndarray) -> np.ndarray:
    """Median weight (LMS M, kg) of reference `ref` at exact age (days) and sex, linear between table ages; NaN when
    the age is not finite, negative or beyond the table, or the sex is not M/F (no extrapolation)."""
    age_days = np.asarray(age_days, dtype=float)
    out = np.full(len(age_days), np.nan)
    for sx in ("M", "F"):
        a, m = gref[(ref, sx)]
        sel = (np.asarray(sex) == sx) & np.isfinite(age_days) & (age_days >= 0) & (age_days <= a[-1])
        out[sel] = np.interp(age_days[sel], a, m)
    return out


@dataclass
class WeightArms:
    weight: np.ndarray              # primary: charted -> cn_ref -> who_ref -> train tiers
    source: np.ndarray
    weight_train_tier: np.ndarray   # sensitivity (B38): charted -> train tiers (the v1.0.2 rule before Nov 17)
    source_train_tier: np.ndarray
    weight_who: np.ndarray          # sensitivity (B38): charted -> WHO -> train tiers
    source_who: np.ndarray


def impute_weights(co: Cohort, age_band_arr: np.ndarray, sex: np.ndarray, wt: np.ndarray,
                   age_days: np.ndarray | None = None, gref: dict | None = None) -> WeightArms:
    """Charted weight first. Primary arm: else the growth-reference median for exact age and sex (China 2009, then
    WHO), else the TRAIN-median tiers. The tiers (age band x sex, age band, sex, overall; each n >= 10 train
    encounters with charted weight) never use an age band for age_invalid (decision 20). Two sensitivity arms for
    B38: train tiers only, and WHO before the tiers."""
    tw, ts_ = _train_tier_weights(co, age_band_arr, sex, wt)
    n = len(wt)
    if gref is None or age_days is None:
        raise GateError("internal: the weight normalizer needs the growth reference and the age in days")
    valid = (np.asarray(age_band_arr) != "age_invalid")
    ad = np.where(valid, np.asarray(age_days, dtype=float), np.nan)
    med = {r: growth_ref_median(gref, r, ad, sex) for r in GROWTH_REFS}

    def arm(refs):
        w, src = wt.copy(), np.where(np.isfinite(wt), "charted", "").astype(object)
        for r in refs:
            fill = (src == "") & np.isfinite(med[r])
            w[fill], src[fill] = med[r][fill], GROWTH_REF_SOURCE[r]
        rest = src == ""
        w[rest], src[rest] = tw[rest], ts_[rest]
        return w, src
    w1, s1 = arm(GROWTH_REFS)
    w3, s3 = arm(("who",))
    for w, s in ((w1, s1), (tw, ts_), (w3, s3)):
        if len(w) != n or (s == "").any() or not (np.isfinite(w) & (w > 0)).all():
            raise GateError("weight could not be filled for some encounters (no train charted weight)")
    return WeightArms(weight=w1, source=s1, weight_train_tier=tw, source_train_tier=ts_, weight_who=w3, source_who=s3)


def _train_tier_weights(co: Cohort, age_band_arr: np.ndarray, sex: np.ndarray, wt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Charted weight, else TRAIN medians: age band x sex, age band, sex, overall (each n >= 10); age_invalid is never
    imputed from an age band (decision 20)."""
    tr = co.train() & np.isfinite(wt)
    df = pd.DataFrame({"w": wt[tr], "ab": age_band_arr[tr], "sx": sex[tr]})
    da = df[df["ab"] != "age_invalid"]
    g2 = da.groupby(["ab", "sx"])["w"].agg(["median", "size"])
    g1 = da.groupby("ab")["w"].agg(["median", "size"])
    gs = df.groupby("sx")["w"].agg(["median", "size"])
    overall = float(df["w"].median()) if len(df) else np.nan
    w = wt.copy()
    src = np.where(np.isfinite(w), "charted", "").astype(object)
    for i in np.flatnonzero(~np.isfinite(w)):
        ab, sx = age_band_arr[i], sex[i]
        if ab != "age_invalid" and (ab, sx) in g2.index and g2.loc[(ab, sx), "size"] >= MIN_CELL:
            w[i], src[i] = g2.loc[(ab, sx), "median"], "train_tier_age_sex"
        elif ab != "age_invalid" and ab in g1.index and g1.loc[ab, "size"] >= MIN_CELL:
            w[i], src[i] = g1.loc[ab, "median"], "train_tier_age_band"
        elif sx in gs.index and gs.loc[sx, "size"] >= MIN_CELL:
            w[i], src[i] = gs.loc[sx, "median"], "train_tier_sex_only"
        else:
            w[i], src[i] = overall, ("train_tier_overall" if np.isfinite(overall) else "")
    if (src == "").any():
        raise GateError("weight could not be filled for some encounters (no train charted weight)")
    return w, src


def build_values(co: Cohort, chart: pd.DataFrame, lab: pd.DataFrame, cmap: dict, first_unit: np.ndarray,
                 age_band_arr: np.ndarray, sex: np.ndarray, age_days: np.ndarray | None = None,
                 gref: dict | None = None) -> Values:
    chans: dict[str, Channel] = cmap["channels"]
    rs, meas, unp, blc = [], [], [], []
    nicu = np.array([str(u).strip().upper() == "NICU" for u in first_unit], dtype=bool)
    for table, ev in (("CHARTEVENTS", chart), ("LABEVENTS", lab)):
        for it, g in ev.groupby("item", sort=True):
            info = cmap["items"].get((table, it))
            if info is None:
                continue
            ch = chans[info["channel"]]
            v, st, side = classify_rows(info, g, ch.qualitative)
            e = g["enc_i"].to_numpy()
            if info["status"] == "assigned_merged_conditional":     # 5255: NICU only (D24)
                st = np.where((st == "measurement") & ~nicu[e], "dropped_non_nicu", st)
            if ch.indicator_only:
                v = np.full(len(v), np.nan)
            d = pd.DataFrame({"enc_i": e.astype(np.int64), "t": g["t"].to_numpy(), "tc": g["tc"].to_numpy(),
                              "item": it, "channel": ch.name, "v": v, "status": st, "side": side})
            rs.append(d)
            if ch.qualitative:
                u = st == "qual_unparsed"
                if u.any():
                    unp.append(pd.DataFrame({"enc_i": e[u].astype(np.int64), "t": g["t"].to_numpy()[u],
                                             "channel": ch.name, "text": qual_norm(g["vtext"].to_numpy()[u]).to_numpy()}))
                bc = (st == "measurement") & qual_below_cutoff(g["vtext"].to_numpy(), ch.name)
                if bc.any():
                    blc.append(pd.DataFrame({"enc_i": e[bc].astype(np.int64), "t": g["t"].to_numpy()[bc],
                                             "channel": ch.name}))
    rowstat = pd.concat(rs, ignore_index=True) if rs else pd.DataFrame(
        columns=["enc_i", "t", "tc", "item", "channel", "v", "status", "side"])
    unparsed = pd.concat(unp, ignore_index=True) if unp else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int64), "t": pd.Series(dtype=float), "channel": pd.Series(dtype=object),
         "text": pd.Series(dtype=object)})
    below_cutoff = pd.concat(blc, ignore_index=True) if blc else pd.DataFrame(
        {"enc_i": pd.Series(dtype=np.int64), "t": pd.Series(dtype=float), "channel": pd.Series(dtype=object)})
    m = rowstat[rowstat["status"] == "measurement"]
    obs: dict[str, pd.DataFrame] = {}
    dd_drop = []
    for name, g in m.groupby("channel", sort=True):
        ch = chans[name]
        if ch.dedupe_same_time and ch.main_item:
            g = g.assign(_o=(g["item"] != ch.main_item).astype(int)).sort_values(["enc_i", "tc", "_o"], kind="mergesort")
            dup = g.duplicated(["enc_i", "tc"], keep="first")
            dd_drop.append(g[dup].index.to_numpy())
            g = g[~dup].drop(columns="_o")
        obs[name] = g[["enc_i", "t", "tc", "item", "v"]].sort_values(["enc_i", "t"], kind="mergesort").reset_index(drop=True)
    if dd_drop:
        idx = np.concatenate(dd_drop)
        rowstat.loc[idx, "status"] = "dedupe_same_draw"
    # MAP (D8): same-CHARTTIME gated SBP/DBP pairs, DBP < SBP, MAP within its gate
    sb = obs.get("sbp", _empty_obs()).drop_duplicates(["enc_i", "tc"])
    db = obs.get("dbp", _empty_obs()).drop_duplicates(["enc_i", "tc"])
    pr = sb[["enc_i", "tc", "t", "v"]].merge(db[["enc_i", "tc", "v"]], on=["enc_i", "tc"], suffixes=("_s", "_d"))
    mv = (pr["v_s"].to_numpy() + 2 * pr["v_d"].to_numpy()) / 3.0
    gm = cmap["derived"]["map"]
    okm = (pr["v_d"].to_numpy() < pr["v_s"].to_numpy()) & (mv >= gm["low"]) & (mv <= gm["high"])
    obs["map"] = pr.assign(v=mv, item="derived:map")[okm][["enc_i", "t", "tc", "item", "v"]] \
        .sort_values(["enc_i", "t"], kind="mergesort").reset_index(drop=True)
    map_info = {"n_pairs": len(pr), "n_dropped": int((~okm).sum())}
    # weight (V12, decision 20): first gated 1014 value from t0 - 24 h to the end of the Layer B window
    wo = obs.pop("weight", _empty_obs())
    bend = np.minimum(co.out_h, LAYER_B_H)
    wsel = wo[(wo["tc"] >= -24.0) & (wo["tc"] < bend[wo["enc_i"].to_numpy()])].sort_values(["enc_i", "tc"], kind="mergesort")
    wf = wsel.drop_duplicates("enc_i", keep="first")
    wt = np.full(co.n, np.nan)
    wt[wf["enc_i"].to_numpy()] = wf["v"].to_numpy()
    arms = impute_weights(co, age_band_arr, sex, wt, age_days, gref)
    weight, wsrc = arms.weight, arms.source
    # urine output (D30): mL/kg/h = volume / weight / hours since the previous entry (cap 24 h; first: since t0)
    uo = obs.get("urine_output", _empty_obs())
    uo = uo[(uo["tc"] >= 0) & (uo["tc"] < co.W[uo["enc_i"].to_numpy()])].sort_values(["enc_i", "tc"], kind="mergesort")
    prev = uo.groupby("enc_i")["tc"].shift(1).fillna(0.0).to_numpy()
    iv = np.minimum(uo["tc"].to_numpy() - prev, UO_CAP_H)
    vol = uo["v"].to_numpy()
    okr = iv > 0
    rate = np.where(okr, vol / weight[uo["enc_i"].to_numpy()] / np.where(okr, iv, 1.0), np.nan)
    uo_rate = uo.assign(v=rate, vol=vol, iv=iv)[okr]
    obs["urine_output"] = uo_rate[["enc_i", "t", "tc", "item", "v"]].reset_index(drop=True)
    # KDIGO urine-output stage (provisional encoding, spec §6.6): evaluated at each entry with >= 24 h since t0
    kd = []
    for e, g in uo_rate.groupby("enc_i", sort=True):
        tcs, vols, ivs = g["tc"].to_numpy(), g["vol"].to_numpy(), g["iv"].to_numpy()
        for j in range(len(tcs)):
            te = tcs[j]
            if te < 24.0:
                continue
            win = (tcs > te - 24.0) & (tcs <= te)
            avg = vols[win].sum() / weight[e] / 24.0
            if vols[j] == 0 and ivs[j] >= 12.0:
                s = 3.0
            elif avg < 0.3:
                s = 3.0
            elif avg < 0.5:
                s = 1.0
            else:
                s = 0.0
            kd.append((e, te, te, "derived:kdigo_uo_stage", s))
    obs["kdigo_uo_stage"] = (pd.DataFrame(kd, columns=["enc_i", "t", "tc", "item", "v"]) if kd else _empty_obs())
    # fluid intake (D17, D30): mL/kg/day = charted intake over the past 24 h (from t0 on) / weight (provisional)
    fi = obs.get("fluid_intake", _empty_obs())
    fi = fi[(fi["tc"] >= 0) & (fi["tc"] < co.W[fi["enc_i"].to_numpy()])].sort_values(["enc_i", "tc"], kind="mergesort")
    fv = []
    for e, g in fi.groupby("enc_i", sort=True):
        tcs, vols = g["tc"].to_numpy(), g["v"].to_numpy()
        for j in range(len(tcs)):
            win = (tcs > tcs[j] - 24.0) & (tcs <= tcs[j])
            fv.append((e, tcs[j], tcs[j], "1009", vols[win].sum() / weight[e]))
    obs["fluid_intake"] = pd.DataFrame(fv, columns=["enc_i", "t", "tc", "item", "v"]) if fv else _empty_obs()
    for name in cmap["feature_channels"]:
        obs.setdefault(name, _empty_obs())
        obs[name]["enc_i"] = obs[name]["enc_i"].astype(np.int64)
    return Values(obs=obs, rowstat=rowstat, weight=weight, weight_source=wsrc,
                  uo_info={"map": map_info, "n_uo_zero_interval": int((~okr).sum())}, weight_arms=arms,
                  unparsed=unparsed, below_cutoff=below_cutoff)


# ----------------------------------------------------------------------------------------------------------------------
# Layer C: bins, masks, time since last, carry-forward (spec §8-§9)
# ----------------------------------------------------------------------------------------------------------------------
def in_grid(co: Cohort, e: np.ndarray, t: np.ndarray, upto: float | None = None) -> np.ndarray:
    ok = (t >= 0) & (t < co.W[e])
    if upto is not None:
        ok &= t < upto
    return ok


def bin_summaries(co: Cohort, ob: pd.DataFrame, upto: float | None = None) -> dict:
    """Per (encounter, bin): median, min, max, last value (latest availability time), its time, measured."""
    shp = (co.n, N_BINS)
    out = {k: np.full(shp, np.nan, dtype=np.float64) for k in ("median", "min", "max", "last", "last_t")}
    out["measured"] = np.zeros(shp, dtype=bool)
    if not len(ob):
        return out
    e, t = ob["enc_i"].to_numpy(), ob["t"].to_numpy()
    d = ob[in_grid(co, e, t, upto)]
    if not len(d):
        return out
    d = d.assign(k=np.floor(d["t"].to_numpy() / BIN_H).astype(int)).sort_values(["enc_i", "k", "t"], kind="mergesort")
    g = d.groupby(["enc_i", "k"], sort=False)
    ix = g.size().index
    ei, ki = ix.get_level_values(0).to_numpy(), ix.get_level_values(1).to_numpy()
    out["measured"][ei, ki] = True
    out["median"][ei, ki] = g["v"].median().to_numpy()
    out["min"][ei, ki] = g["v"].min().to_numpy()
    out["max"][ei, ki] = g["v"].max().to_numpy()
    out["last"][ei, ki] = g["v"].last().to_numpy()
    out["last_t"][ei, ki] = g["t"].max().to_numpy()
    return out


def bin_end(co: Cohort) -> np.ndarray:
    te = BIN_H * (np.arange(N_BINS)[None, :] + 1.0)
    return np.minimum(te, co.out_h[:, None])


def time_since_last(co: Cohort, last_t: np.ndarray, measured: np.ndarray, cap_h: float) -> np.ndarray:
    lt = np.where(measured, last_t, np.nan)
    ff = pd.DataFrame(lt).ffill(axis=1).to_numpy()
    tsl = np.minimum(cap_h, bin_end(co) - ff)
    tsl = np.where(np.isnan(tsl), cap_h, tsl)
    return np.where(co.mask_valid.astype(bool), tsl, np.nan)


def locf_fill(med: np.ndarray, last: np.ndarray, cap_bins: int) -> tuple[np.ndarray, np.ndarray]:
    """Carry the last observed value of the most recent observed bin within the cap (bins). Carried values never
    seed a later carry. Returns (filled values, carried flag)."""
    val = med.copy()
    carried = np.zeros(med.shape, dtype=bool)
    obs = ~np.isnan(med)
    for j in range(1, cap_bins + 1):
        src_obs = np.zeros(med.shape, dtype=bool)
        src_obs[:, j:] = obs[:, :-j]
        src_val = np.full(med.shape, np.nan)
        src_val[:, j:] = last[:, :-j]
        take = np.isnan(val) & src_obs
        val[take] = src_val[take]
        carried |= take
    return val, carried


class NpyStore:
    """Memory-mapped .npy members in a pending folder, zipped into one .npz at the end (no whole-array RAM copy)."""

    def __init__(self, folder: Path):
        self.folder = folder
        fs_op("create folder", lambda: folder.mkdir(parents=True, exist_ok=False), folder)
        self.members: dict[str, Path] = {}
        self.arrays: list = []
        self.maps: list = []   # (member path, weak reference to its memory map), checked before zipping/deleting

    def create(self, name: str, shape: tuple, dtype, fill) -> np.ndarray:
        p = self.folder / f"{name}.npy"
        a = fs_op("create memory map", lambda: np.lib.format.open_memmap(p, mode="w+", dtype=dtype, shape=shape), p)
        a[...] = fill
        self.members[name] = p
        self.arrays.append(a)
        self.maps.append((p, weakref.ref(a)))
        return a

    def put(self, name: str, arr: np.ndarray) -> None:
        p = self.folder / f"{name}.npy"
        fs_op("write member", lambda: np.save(p, arr, allow_pickle=False), p)
        self.members[name] = p

    def close(self) -> None:
        """Flush and drop the store's own memory maps, then check that none is still open (a reference held anywhere
        else keeps the file mapped; Windows then refuses to delete it). Callers drop their references first."""
        for a in self.arrays:
            a.flush()
        self.arrays.clear()
        a = None  # noqa: F841
        gc.collect()
        alive = [p for p, r in self.maps if r() is not None]
        if alive:
            raise FileOpError("release memory map", alive[:1], reason=f"{len(alive)} memory map(s) still referenced")

    def write_npz(self, path: Path) -> None:
        """Callers must drop their own references to the memmaps first (Windows cannot delete mapped files)."""
        self.close()
        tmp = path.with_name(path.name + ".tmp")

        def zip_members():
            with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1, allowZip64=True) as z:
                for name, p in self.members.items():
                    z.write(p, arcname=f"{name}.npy")

        fs_op("write zip", zip_members, tmp)
        fs_op("fsync", lambda: _fsync_path(tmp), tmp)
        fs_op("replace", lambda: os.replace(tmp, path), tmp, path)
        remove_tree(self.folder, "remove member folder")


# ----------------------------------------------------------------------------------------------------------------------
# Medication state machine (S3; plan §5.3) and microbiology / surgery events
# ----------------------------------------------------------------------------------------------------------------------
def group_cutoff(group: str) -> float:
    return VASO_CUTOFF_H if group == "vasoactive" else DEFAULT_CUTOFF_H


def build_courses(enc: np.ndarray, key: np.ndarray, S: np.ndarray, E: np.ndarray, cutoff: np.ndarray) -> pd.DataFrame:
    """Courses per (encounter, course key); inputs sorted by (enc, key, S). A documentation joins the open course if
    S - last <= cutoff and no ENDDATE of the course (after its own start) is at or before S. Stop = min(last + cutoff,
    earliest valid ENDDATE)."""
    n = len(S)
    res = []
    i = 0
    while i < n:
        ge, gk, c = enc[i], key[i], cutoff[i]
        start = last = S[i]
        endE = E[i] if (E[i] == E[i] and E[i] > S[i]) else math.inf
        j = i + 1
        while j < n and enc[j] == ge and key[j] == gk:
            sj, ej = S[j], E[j]
            if sj - last <= c and sj < endE:
                last = sj
                if ej == ej and ej > sj and ej < endE:
                    endE = ej
            else:
                res.append((ge, gk, start, last, endE, c))
                start = last = sj
                endE = ej if (ej == ej and ej > sj) else math.inf
            j += 1
        res.append((ge, gk, start, last, endE, c))
        i = j
    cr = pd.DataFrame(res, columns=["enc_i", "key", "start", "last", "endE", "cutoff"])
    cr["stop"] = np.minimum(cr["last"] + cr["cutoff"], cr["endE"])
    return cr


def courses_for(docs: pd.DataFrame, lo: float, as_of: float | None = None) -> pd.DataFrame:
    d = docs[docs["group"].isin(MED_GROUPS) & (docs["S"] >= lo)]
    if as_of is not None:
        d = d[d["S"] <= as_of].copy()
        d.loc[d["E"] > as_of, "E"] = np.nan
    d = d.sort_values(["enc_i", "key", "S"], kind="mergesort")
    cut = d["group"].map(group_cutoff).to_numpy(dtype=float)
    cr = build_courses(d["enc_i"].to_numpy(), d["key"].to_numpy(), d["S"].to_numpy(dtype=float),
                       d["E"].to_numpy(dtype=float), cut)
    gk = d.drop_duplicates(["enc_i", "key"])[["enc_i", "key", "group"]]
    if not len(cr):
        return cr.assign(group=pd.Series(dtype=object))
    return cr.merge(gk, on=["enc_i", "key"], how="left", validate="many_to_one")


def med_state_at(cr: pd.DataFrame, docs: pd.DataFrame, n: int, t: np.ndarray, lo: float, upto: float | None = None
                 ) -> dict:
    """State features per encounter at times t (array length n) for every group: active, hours_on_course,
    hours_since_stopped, ever_received (documentation in [lo, t])."""
    out = {}
    for g in MED_GROUPS:
        c = cr[cr["group"] == g]
        e = c["enc_i"].to_numpy(dtype=np.int64)
        te = t[e]
        act = (c["start"].to_numpy() <= te) & (te <= c["last"].to_numpy() + c["cutoff"].to_numpy()) & (te < c["endE"].to_numpy())
        active = np.zeros(n, dtype=bool)
        np.logical_or.at(active, e[act], True)
        hoc = np.zeros(n)
        if act.any():
            np.maximum.at(hoc, e[act], (te - c["start"].to_numpy())[act])
        stp = c["stop"].to_numpy()
        stopped = stp <= te
        last_stop = np.full(n, -np.inf)
        if stopped.any():
            np.maximum.at(last_stop, e[stopped], stp[stopped])
        hss = np.where(np.isfinite(last_stop) & ~active, t - last_stop, np.nan)
        dg = docs[(docs["group"] == g) & (docs["S"] >= lo)]
        de = dg["enc_i"].to_numpy(dtype=np.int64)
        ev = np.zeros(n, dtype=bool)
        sel = dg["S"].to_numpy() <= t[de]
        np.logical_or.at(ev, de[sel], True)
        out[g] = {"active": active, "hours_on_course": np.where(active, hoc, 0.0), "hours_since_stopped": hss,
                  "ever_received": ev}
    return out


def med_bins(co: Cohort, docs: pd.DataFrame, as_of: float | None = None) -> dict:
    """Layer C medication arrays (n, 56) per group, evaluated at each bin end (Q10: active = any time in the bin's
    valid overlap; started in the start bin; other state features at bin end)."""
    lo = -EVER_LOOKBACK_H
    cr = courses_for(docs, lo, as_of)
    valid = co.mask_valid.astype(bool)
    te = bin_end(co)
    res = {g: {k: np.zeros((co.n, N_BINS), dtype=np.float32) for k in
               ("active", "started", "hours_on_course", "hours_since_stopped", "ever_received")} for g in MED_GROUPS}
    for g in MED_GROUPS:
        c = cr[cr["group"] == g]
        e = c["enc_i"].to_numpy(dtype=np.int64)
        st, sp = c["start"].to_numpy(), c["stop"].to_numpy()
        for k in range(N_BINS):
            a = (st < te[e, k]) & (sp > BIN_H * k) & valid[e, k]
            res[g]["active"][e[a], k] = 1
            s = (st >= BIN_H * k) & (st < te[e, k]) & (st >= 0) & valid[e, k]
            res[g]["started"][e[s], k] = 1
    for k in range(N_BINS):
        stt = med_state_at(cr, docs, co.n, te[:, k], lo)
        for g in MED_GROUPS:
            res[g]["hours_on_course"][:, k] = stt[g]["hours_on_course"]
            res[g]["hours_since_stopped"][:, k] = stt[g]["hours_since_stopped"]
            res[g]["ever_received"][:, k] = stt[g]["ever_received"]
    vaso = np.zeros((co.n, N_BINS), dtype=np.float32)
    c = cr[cr["group"] == "vasoactive"]
    for k in range(N_BINS):
        vaso[:, k] = max_concurrent(c, co.n, np.full(co.n, BIN_H * k), te[:, k])
    for g in MED_GROUPS:
        for f in res[g]:
            arr = res[g][f]
            arr[~valid] = 0 if f != "hours_since_stopped" else np.nan
    vaso[~valid] = 0
    return {"groups": res, "vaso_concurrent": vaso, "courses": cr}


def future_enddate_differences(co: Cohort, docs: pd.DataFrame) -> int:
    """Decision 7 (Nov 17, 2026) as a run check: the medication state at an assessment time never uses a stop time
    (ENDDATE) after that time. Drugs at t0 (documentation in [t0 - 6 h, t0)) and the Layer B end state are rebuilt
    with every ENDDATE after the assessment time removed; any difference to the state the builder uses is counted.
    (Layer C bins are covered by the as-of rebuild check, which drops ENDDATEs after the landmark.)"""
    n = co.n
    wb = np.minimum(co.out_h, LAYER_B_H)
    d0 = docs[(docs["S"] >= -PRE_ICU_H) & (docs["S"] < 0)]
    cases = ((np.zeros(n), -PRE_ICU_H, d0, courses_for(d0, -PRE_ICU_H, as_of=0.0)),
             (wb, -EVER_LOOKBACK_H, docs, courses_for(docs, -EVER_LOOKBACK_H)))
    diff = 0
    for t, lo, d, cr_used in cases:
        dm = d.copy()
        e = dm["enc_i"].to_numpy(dtype=np.int64)
        dm.loc[dm["E"].to_numpy(dtype=float) > t[e], "E"] = np.nan
        a = med_state_at(cr_used, d, n, t, lo)
        b = med_state_at(courses_for(dm, lo), dm, n, t, lo)
        for g in MED_GROUPS:
            for f in ("active", "hours_on_course", "hours_since_stopped", "ever_received"):
                x, y = np.asarray(a[g][f], dtype=float), np.asarray(b[g][f], dtype=float)
                diff += int((~((x == y) | (np.isnan(x) & np.isnan(y)))).sum())
    return diff


def max_concurrent(cr: pd.DataFrame, n: int, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    res = np.zeros(n, dtype=int)
    if not len(cr):
        return res
    e = cr["enc_i"].to_numpy(dtype=np.int64)
    a = np.maximum(cr["start"].to_numpy(), lo[e])
    b = np.minimum(cr["stop"].to_numpy(), hi[e])
    ok = b > a
    if not ok.any():
        return res
    ev = pd.DataFrame({"e": np.r_[e[ok], e[ok]], "t": np.r_[a[ok], b[ok]],
                       "d": np.r_[np.ones(ok.sum(), dtype=int), -np.ones(ok.sum(), dtype=int)]})
    ev = ev.sort_values(["e", "t", "d"], kind="mergesort")
    ev["c"] = ev.groupby("e")["d"].cumsum()
    m = ev.groupby("e")["c"].max()
    res[m.index.to_numpy()] = m.to_numpy()
    return res


def active_measure(cr: pd.DataFrame, n: int, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Length of the union of course intervals inside [lo, hi) per encounter (hrs_active)."""
    out = np.zeros(n)
    if not len(cr):
        return out
    e = cr["enc_i"].to_numpy(dtype=np.int64)
    a = np.maximum(cr["start"].to_numpy(), lo[e])
    b = np.minimum(cr["stop"].to_numpy(), hi[e])
    ok = b > a
    d = pd.DataFrame({"e": e[ok], "a": a[ok], "b": b[ok]}).sort_values(["e", "a"], kind="mergesort")
    for ei, g in d.groupby("e", sort=False):
        tot, cs, ce = 0.0, None, None
        for x, y in zip(g["a"].to_numpy(), g["b"].to_numpy()):
            if cs is None or x > ce:
                if cs is not None:
                    tot += ce - cs
                cs, ce = x, y
            else:
                ce = max(ce, y)
        tot += ce - cs
        out[int(ei)] = tot
    return out


# ----------------------------------------------------------------------------------------------------------------------
# Surgery (v1.7 module: ROW_ID key, anaesthesia-end timing T, linkage; φ1 band from the frozen medians)
# ----------------------------------------------------------------------------------------------------------------------
def surgery_frame(cfg: Config, A, co: Cohort, med1: dict) -> tuple[pd.DataFrame, dict]:
    acfg = A.AuditConfig(canonical_path=Path(cfg.canonical_path), harmonization_path=Path(cfg.harmonization_path),
                         map_path=Path(cfg.complexity_map_path), split_path=Path(cfg.split_path),
                         surgery_info_path=raw_path(cfg, "SURGERY_INFO"), out_dir=Path(cfg.out_dir),
                         admissions_path=raw_path(cfg, "ADMISSIONS"), harm_t0_col=cfg.harm_t0_col,
                         harm_out_col=cfg.harm_out_col, split_subject_col=cfg.split_subject_col,
                         split_col=cfg.split_col, expected_surgery_rows=cfg.raw_pins["SURGERY_INFO"][1],
                         expected_encounters=cfg.expect.encounters, expected_subjects=cfg.expect.subjects,
                         expected_split_encounters=dict(cfg.expect.splits), verbose=False)
    if cfg.surgery_expected_missing is not None:
        acfg.expected_missing = dict(cfg.surgery_expected_missing)
    acfg.pinned_sha256 = {k: cfg.pins[k] for k in ("canonical", "harmonization", "complexity_map", "split")}
    try:
        inp = A.load_inputs(acfg)
        s, _meta = A.build_frame(inp)
    except A.AuditError as e:
        raise InputError(f"surgery inputs (v1.7 module): {e}") from e
    surg = inp["surg"]
    AS, _, _ = A.parse_time(surg["ANES_START_TIME"])
    s = s.sort_values("sid", kind="mergesort").reset_index(drop=True)
    s["AS"] = AS.to_numpy()[s["sid"].to_numpy()]
    vis = surg["VISIT_ID"].map(A.norm_id).to_numpy()[s["sid"].to_numpy()]
    visit = {"n_rows": len(s), "n_visit_equals_hadm": int((vis == s["HADM_ID"].to_numpy()).sum()),
             "n_visit_blank": int((vis == "").sum())}
    visit["n_visit_differs"] = visit["n_rows"] - visit["n_visit_equals_hadm"] - visit["n_visit_blank"]
    cm = inp["map"].set_index("key")
    s["standardized_name"] = s["name"].map(cm["standardized_name"]).fillna("").astype(str)
    idx = co.enc.set_index("encounter_id")["enc_i"]
    s["enc_i"] = s["encounter_id"].map(idx)
    lk = s[(s["link_class"] == "linked") & s["enc_i"].notna()].copy()
    lk["enc_i"] = lk["enc_i"].astype(np.int64)
    t0 = co.enc["t0"].to_numpy()[lk["enc_i"].to_numpy()]
    lk["T_h"] = hours_between(lk["T"], t0)
    d_be = hours_between(lk["E"], lk["B"])
    d_an = hours_between(lk["AE"], lk["AS"])
    prim = np.isfinite(d_be) & (d_be > 0) & (d_be <= SURGERY_MAX_H)
    reason = np.where(prim, "", np.where(~np.isfinite(d_be), "missing", np.where(d_be <= 0, "order", "gt24")))
    fb = (~prim) & np.isfinite(d_an) & (d_an > 0) & (d_an <= SURGERY_MAX_H)
    lk["duration_h"] = np.where(prim, d_be, np.where(fb, d_an, np.nan))
    lk["duration_source"] = np.where(prim, "surgery_window", np.where(fb, "anesthesia_window", "none"))
    lk["fallback_reason"] = reason
    lk["band_missing_flag"] = (~prim & ~fb).astype(np.uint8)
    keys = [("type:" + A.norm_name(sn), "own_type") for sn in lk["standardized_name"]]
    lvl, ref = [], []
    for (kt, _), bs, cg, cl in zip(keys, lk["body_system"], lk["cgroup"], lk["classified"]):
        cands = [(kt, "own_type")]
        if cl:
            cands.append((f"system_x_group:{bs} x {cg}", "system_x_group"))
        cands += [(f"system:{bs}", "system"), ("all", "all")]
        for k, l in cands:
            if k in med1:
                lvl.append(l)
                ref.append(med1[k])
                break
    lk["band_ref_level"] = lvl
    r = lk["duration_h"].to_numpy() / np.array(ref, dtype=float)
    band = np.where(~np.isfinite(r), "missing", np.where(r < BAND_SHORT, "short", np.where(r > BAND_LONG, "long", "expected")))
    lk["phi1_band"] = band
    lk = lk.sort_values(["enc_i", "T_h", "ROW_ID"], kind="mergesort").reset_index(drop=True)
    return lk, {"visit": visit, "n_linked": len(lk)}


# ----------------------------------------------------------------------------------------------------------------------
# Layer A base (spec §5)
# ----------------------------------------------------------------------------------------------------------------------
def interval_band(h: np.ndarray) -> np.ndarray:
    return np.where(~np.isfinite(h), "", np.where(h < 24, "lt24h", np.where(h < 72, "24to72h",
                                                                            np.where(h < 168, "72to168h", "ge168h"))))


def shift_of(hour: np.ndarray) -> np.ndarray:
    return np.where((hour >= 8) & (hour < 15), "day", np.where((hour >= 15) & (hour < 23), "evening", "night"))


def build_layer_a(co: Cohort, pt: pd.DataFrame, ad: pd.DataFrame, icu: dict, dicd: dict, neuro_match, surg: pd.DataFrame,
                  vals: Values, cmap: dict, docs: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    n = co.n
    t0 = co.enc["t0"]
    A = pd.DataFrame({"encounter_id": co.enc["encounter_id"].to_numpy(), "split": co.split,
                      "sex": pt["sex"].to_numpy(), "first_careunit": icu["first_careunit"]})
    A["age_status"] = pt["age_status"].to_numpy()
    A["age_invalid"] = (A["age_status"] != "valid").astype(np.uint8)
    A["age_neg_bin"] = pt["age_neg_bin"].to_numpy()
    A["age_days"] = pt["age_days"].to_numpy()
    A["age_months"] = A["age_days"] / (365.25 / 12.0)
    A["age_band"] = age_band(A["age_days"].to_numpy())
    A["age_days_sens_lt1d0"] = pt["age_days_sens_lt1d0"].to_numpy()
    A["age_band_sens_lt1d0"] = age_band(A["age_days_sens_lt1d0"].to_numpy())
    raw_iv = hours_between(t0, ad["ADMITTIME"])
    neg = np.isfinite(raw_iv) & (raw_iv < 0)
    A["hours_admit_to_icu_intime"] = np.where(neg, 0.0, raw_iv)                     # Q4: clipped at 0
    A["admit_interval_clipped"] = neg.astype(np.uint8)
    A["admit_interval_missing"] = (~np.isfinite(raw_iv)).astype(np.uint8)
    A["admit_interval_band"] = interval_band(A["hours_admit_to_icu_intime"].to_numpy())
    A["intime_hour"] = t0.dt.hour.to_numpy()
    A["intime_dow"] = t0.dt.dayofweek.to_numpy()
    A["intime_shift"] = shift_of(A["intime_hour"].to_numpy())
    A["intime_weekend"] = (A["intime_dow"] >= 5).astype(np.uint8)
    # diagnosis (Q8): ICD10_CODE_CN -> D_ICD_DIAGNOSES ICD10_CODE; unmapped codes -> explicit 'unmapped' levels
    st, c3, chap, mapped, raw_first = [], [], [], [], []
    for code in ad["dx"].to_numpy():
        s = str(code or "").strip().upper()
        nrm, status = norm_icd(s)
        st.append(status)
        parts = [p for p in re.split(r"[,;/|，；、\s]+", s) if p]
        first = parts[0] if parts else ""
        tgt = dicd.get(first, "")
        raw_first.append(nrm)
        if status == "missing":
            c3.append("missing"); chap.append("missing"); mapped.append(0)
        elif tgt and re.match(r"^[A-Z][0-9]{2}", tgt):
            c3.append(tgt[:3]); chap.append(tgt[0]); mapped.append(1)
        elif is_icdo_morphology(s):               # decision 8 (Nov 17, 2026): ICD-O morphology -> oncology group
            c3.append(ICDO_LEVEL); chap.append(ICDO_LEVEL); mapped.append(0)
        elif re.match(r"^[A-Z]", nrm or ""):
            c3.append(f"unmapped_{nrm[0]}"); chap.append(nrm[0]); mapped.append(0)
        else:
            c3.append("unmapped"); chap.append("unmapped"); mapped.append(0)
    A["dx_code_status"] = st
    A["dx_mapped"] = np.array(mapped, dtype=np.uint8)
    A["dx_icd10_3char"] = c3
    A["dx_chapter"] = chap
    A["dx_missing"] = (A["dx_code_status"] == "missing").astype(np.uint8)
    A["dx_icdo_morphology"] = (A["dx_icd10_3char"] == ICDO_LEVEL).astype(np.uint8)
    A["dx_domain"] = [icd_domain(c) if m else ("missing" if c == "missing" else
                                               ONCOLOGY_DOMAIN if c == ICDO_LEVEL else "other_medical")
                      for c, m in zip(c3, mapped)]
    A["dx_neuro_set"] = np.array([neuro_match(c) for c in raw_first], dtype=np.uint8)
    # surgery: recent (T in [t0 - 24 h, t0 + 6 h]) and remote (same HADM, T < t0 - 24 h)
    rec = surg[surg["timing"] == "recent_layer_a"]
    rl = rec.drop_duplicates("enc_i", keep="last").set_index("enc_i")
    rem = surg[surg["timing"] == "remote_same_hadm"]
    ml = rem.drop_duplicates("enc_i", keep="last").set_index("enc_i")
    ix = np.arange(n)
    A["recent_surgery"] = np.isin(ix, rec["enc_i"].to_numpy()).astype(np.uint8)
    A["recent_surgery_name"] = pd.Series(ix).map(rl["standardized_name"]).fillna("").to_numpy()
    A["recent_hours_T_to_t0"] = -pd.Series(ix).map(rl["T_h"]).to_numpy(dtype=float)
    A["recent_body_system"] = pd.Series(ix).map(rl["body_system"]).fillna("").to_numpy()
    A["recent_complexity_group"] = pd.Series(ix).map(rl["cgroup"]).fillna("").to_numpy()
    A["recent_phi1_band"] = pd.Series(ix).map(rl["phi1_band"]).fillna("").to_numpy()
    A["recent_band_ref_level"] = pd.Series(ix).map(rl["band_ref_level"]).fillna("").to_numpy()
    A["recent_duration_source"] = pd.Series(ix).map(rl["duration_source"]).fillna("").to_numpy()
    A["recent_band_missing_flag"] = pd.Series(ix).map(rl["band_missing_flag"]).fillna(0).to_numpy(dtype=np.uint8)
    A["remote_surgery"] = np.isin(ix, rem["enc_i"].to_numpy()).astype(np.uint8)
    A["remote_surgery_name"] = pd.Series(ix).map(ml["standardized_name"]).fillna("").to_numpy()
    A["remote_body_system"] = pd.Series(ix).map(ml["body_system"]).fillna("").to_numpy()
    pre_any = rem.groupby("enc_i")["pre_admit"].max() if len(rem) else pd.Series(dtype=bool)
    A["remote_pre_admit_flag"] = pd.Series(ix).map(pre_any).fillna(False).astype(np.uint8).to_numpy()
    # admission source components (V17; Q5: X is chosen before stage 06; components only). Decision 9 (Nov 17,
    # 2026): EDREGTIME is not an ED marker feature (no adm_src_ed_component); it stays an audit aggregate only.
    or_comp = np.isin(ix, rec.loc[rec["T_h"] <= 0, "enc_i"].to_numpy())
    A["adm_src_interval_h"] = A["hours_admit_to_icu_intime"]
    A["adm_src_or_component"] = or_comp.astype(np.uint8)
    # pre-ICU lab block [t0 - 6 h, t0) (decision 19) and the 24 h sensitivity: last released value + mask
    lab_value = [c for c in cmap["value_channels"] if cmap["channels"][c].table == "LABEVENTS"]
    X: dict[str, np.ndarray] = {}
    for w, tag in ((PRE_ICU_H, "pre6"), (PRE_ICU_SENS_H, "pre24")):
        for ch in lab_value:
            ob = vals.obs[ch]
            sel = ob[(ob["t"] >= -w) & (ob["t"] < 0)].sort_values(["enc_i", "t"], kind="mergesort")
            last = sel.drop_duplicates("enc_i", keep="last")
            v = np.full(n, np.nan)
            v[last["enc_i"].to_numpy()] = last["v"].to_numpy()
            X[f"{tag}_{ch}"] = v
            X[f"{tag}_{ch}_missing"] = np.isnan(v).astype(np.uint8)
    # drugs active at t0 (S3; documentation before t0 within the lookback)
    for w, tag in ((PRE_ICU_H, "t0"), (PRE_ICU_SENS_H, "t0_lb24")):
        d = docs[(docs["S"] >= -w) & (docs["S"] < 0)]
        cr = courses_for(d, -w, as_of=0.0)
        stt = med_state_at(cr, d, n, np.zeros(n), -w)
        for g in MED_GROUPS:
            X[f"med_{g}_active_{tag}"] = stt[g]["active"].astype(np.uint8)
    A = pd.concat([A, pd.DataFrame(X, index=A.index)], axis=1)
    info = {"n_negative_interval_clipped": int(neg.sum()), "n_interval_missing": int((~np.isfinite(raw_iv)).sum())}
    return A, info


# ----------------------------------------------------------------------------------------------------------------------
# Layer B base: [t0, min(OUTTIME, t0 + 24 h)) (spec §6)
# ----------------------------------------------------------------------------------------------------------------------
def build_layer_b(co: Cohort, vals: Values, cmap: dict, docs: pd.DataFrame, micro: pd.DataFrame) -> pd.DataFrame:
    n = co.n
    wb = np.minimum(co.out_h, LAYER_B_H)
    B = {"encounter_id": co.enc["encounter_id"].to_numpy()}
    for ch in cmap["feature_channels"]:
        c = cmap["channels"][ch]
        ob = vals.obs[ch]
        e = ob["enc_i"].to_numpy()
        d = ob[(ob["t"] >= 0) & (ob["t"] < wb[e])].sort_values(["enc_i", "t"], kind="mergesort")
        meas = np.zeros(n, dtype=np.uint8)
        meas[d["enc_i"].to_numpy()] = 1
        B[f"{ch}_measured"] = meas
        if c.indicator_only:
            continue
        g = d.groupby("enc_i")["v"]
        for stat, ser in (("first", g.first()), ("last", g.last()), ("min", g.min()), ("max", g.max()),
                          ("median", g.median())):
            v = np.full(n, np.nan)
            v[ser.index.to_numpy()] = ser.to_numpy()
            B[f"{ch}_{stat}"] = v
        if ch == "hemoglobin":
            hb = np.zeros(n, dtype=np.uint8)
            hb[d.loc[d["item"] == "5257", "enc_i"].to_numpy()] = 1
            B["hemoglobin_src_bloodgas"] = hb
        if ch == "bilirubin_total":
            bb = np.zeros(n, dtype=np.uint8)
            bb[d.loc[d["item"] == "5255", "enc_i"].to_numpy()] = 1
            B["bilirubin_total_src_5255"] = bb
    rsm = vals.rowstat[vals.rowstat["status"] == "measurement"]
    for pg, its in cmap["panels"].items():
        d = rsm[rsm["item"].isin(its)]
        e = d["enc_i"].to_numpy(dtype=np.int64)
        d = d[(d["t"] >= 0) & (d["t"] < wb[e])]
        pm = np.zeros(n, dtype=np.uint8)
        pm[d["enc_i"].to_numpy(dtype=np.int64)] = 1
        B[f"panel_{pg}_measured"] = pm
    # decision 4: unparsed qualitative text stays missing; flag (not a feature column) per qualitative channel
    up = vals.unparsed if vals.unparsed is not None else pd.DataFrame(columns=["enc_i", "t", "channel"])
    ue = up["enc_i"].to_numpy(dtype=np.int64)
    upb = up[(up["t"].to_numpy(dtype=float) >= 0) & (up["t"].to_numpy(dtype=float) < wb[ue])]
    for ch in qual_channels(cmap):
        f = np.zeros(n, dtype=np.uint8)
        f[upb.loc[upb["channel"] == ch, "enc_i"].to_numpy(dtype=np.int64)] = 1
        B[f"{ch}{UNPARSED_FLAG_SUFFIX}"] = f
    # v1.0.4: a negative read from a below-cutoff result in the Layer B window (value stays 0; not a feature column)
    bl = vals.below_cutoff if vals.below_cutoff is not None else pd.DataFrame(columns=["enc_i", "t", "channel"])
    be = bl["enc_i"].to_numpy(dtype=np.int64)
    blb = bl[(bl["t"].to_numpy(dtype=float) >= 0) & (bl["t"].to_numpy(dtype=float) < wb[be])]
    for ch in sorted(QUAL_BELOW_CUTOFF):
        if ch not in cmap["feature_channels"]:
            continue
        f = np.zeros(n, dtype=np.uint8)
        f[blb.loc[blb["channel"] == ch, "enc_i"].to_numpy(dtype=np.int64)] = 1
        B[f"{ch}{BELOW_CUTOFF_FLAG_SUFFIX}"] = f
    lo = -EVER_LOOKBACK_H
    cr = courses_for(docs, lo)
    stt = med_state_at(cr, docs, n, wb, lo)
    for g in MED_GROUPS:
        c = cr[cr["group"] == g]
        e = c["enc_i"].to_numpy(dtype=np.int64)
        s = c["start"].to_numpy()
        stv = (s >= 0) & (s < wb[e])
        h2s = np.full(n, np.inf)
        if stv.any():
            np.minimum.at(h2s, e[stv], s[stv])
        B[f"med_{g}_started"] = np.isfinite(h2s).astype(np.uint8)
        B[f"med_{g}_hrs_to_start"] = np.where(np.isfinite(h2s), h2s, np.nan)
        B[f"med_{g}_hrs_active"] = active_measure(c, n, np.zeros(n), wb)
        B[f"med_{g}_active_end"] = stt[g]["active"].astype(np.uint8)
        B[f"med_{g}_hours_on_course"] = stt[g]["hours_on_course"]
        B[f"med_{g}_hours_since_stopped"] = stt[g]["hours_since_stopped"]
        B[f"med_{g}_ever_received"] = stt[g]["ever_received"].astype(np.uint8)
    B["med_vasoactive_max_concurrent"] = max_concurrent(cr[cr["group"] == "vasoactive"], n, np.zeros(n), wb)
    me = micro["enc_i"].to_numpy(dtype=np.int64)
    mb = micro[(micro["th"] >= 0) & (micro["th"] < wb[me])]
    anyc = np.zeros(n, dtype=np.uint8)
    anyc[mb["enc_i"].to_numpy(dtype=np.int64)] = 1
    B["micro_any_culture"] = anyc
    h1 = mb.groupby("enc_i")["th"].min()
    v = np.full(n, np.nan)
    v[h1.index.to_numpy(dtype=np.int64)] = h1.to_numpy()
    B["micro_hours_to_first_culture"] = v
    for sg in SPECIMEN_GROUPS:
        f = np.zeros(n, dtype=np.uint8)
        f[mb.loc[mb["group"] == sg, "enc_i"].to_numpy(dtype=np.int64)] = 1
        B[f"micro_{sg}"] = f
    return pd.DataFrame(B)


# ----------------------------------------------------------------------------------------------------------------------
# Layer C build, filled copies per landmark, acceptance checks (spec §8-§10, §12)
# ----------------------------------------------------------------------------------------------------------------------
UNPARSED_FLAG_SUFFIX = "_unparsed_text"


def qual_channels(cmap: dict) -> list:
    return [ch for ch in cmap["feature_channels"] if cmap["channels"][ch].qualitative]


def landmark_train_medians(med: np.ndarray, train: np.ndarray, out_h: np.ndarray) -> tuple[dict, dict]:
    """Decision 3 (Nov 17, 2026). med: (n, N_BINS) observed bin medians of one channel (NaN = no value), train: (n,)
    bool, out_h: (n,) hours to the window end. The landmark train median at L is the median over observed bins of
    train encounters eligible at L (OUTTIME > t0 + L) in bins 0 .. L/6 - 1. If it is empty, the nearest earlier
    landmark with a non-empty median of its own is used; if there is none, the pooled all-landmark train median (every
    observed bin of a train encounter inside at least one eligible prefix); if that is empty too the value stays
    empty ('none'). Train rows only. Returns ({L: value}, {L: source})."""
    val, src, own_last = {}, {}, None
    pooled = None
    for L in LANDMARKS:
        nb = L // int(BIN_H)
        rows = np.flatnonzero((out_h > L) & train)
        v = med[rows, :nb]
        v = v[~np.isnan(v)]
        if len(v):
            m, s = float(np.median(v)), "landmark"
            own_last = (L, m)
        elif own_last is not None:
            m, s = own_last[1], f"earlier_landmark_L{own_last[0]:03d}"
        else:
            if pooled is None:
                reach = np.zeros(len(out_h), dtype=int)
                for L2 in LANDMARKS:
                    reach = np.where(out_h > L2, L2 // int(BIN_H), reach)
                pm = (np.arange(med.shape[1])[None, :] < reach[:, None]) & train[:, None] & ~np.isnan(med)
                pv = med[pm]
                pooled = float(np.median(pv)) if len(pv) else np.nan
            m, s = pooled, ("pooled_all_landmarks" if np.isfinite(pooled) else "none")
        val[L], src[L] = m, s
    return val, src


BANNED_CHANNELS = ("pf_ratio", "sf_ratio", "fio2", "aa_gradient", "pf_ratio_lab_crosscheck", "anion_gap", "height")


def build_layer_c(co: Cohort, vals: Values, cmap: dict, docs: pd.DataFrame, micro: pd.DataFrame, surg: pd.DataFrame,
                  pending: Path, fails: list) -> dict:
    feats = cmap["feature_channels"]
    vch = cmap["value_channels"]
    q, qv = len(feats), len(vch)
    n = co.n
    valid = co.mask_valid.astype(bool)
    raw = NpyStore(pending / "_raw_members")
    raw.put("encounter_id", co.enc["encounter_id"].to_numpy(dtype=str))
    raw.put("mask_valid", co.mask_valid)
    raw.put("channels", np.array(feats, dtype=str))
    arr = {k: raw.create(k, (n, N_BINS, q), np.float32, np.nan) for k in ("median", "min", "max", "last", "time_since_last")}
    measured = raw.create("measured", (n, N_BINS, q), np.uint8, 0)
    missing = raw.create("missing", (n, N_BINS, q), np.uint8, 1)
    fstores = {}
    for L in LANDMARKS:
        nb = L // int(BIN_H)
        rows = np.flatnonzero(co.out_h > L)
        st = NpyStore(pending / f"_filled_{L}")
        st.put("encounter_index", rows.astype(np.int64))
        st.put("encounter_id", co.enc["encounter_id"].to_numpy(dtype=str)[rows])
        st.put("value_channels", np.array(vch, dtype=str))
        fstores[L] = {"store": st, "rows": rows, "nb": nb,
                      "value": st.create("value_filled", (len(rows), nb, qv), np.float32, np.nan),
                      "carried": st.create("carried_forward", (len(rows), nb, qv), np.uint8, 0),
                      "imputed": st.create("imputed", (len(rows), nb, qv), np.uint8, 0),
                      "missing": st.create("missing", (len(rows), nb, qv), np.uint8, 1),
                      "median": np.full(qv, np.nan), "median_source": np.full(qv, "none", dtype=object)}
    tr = co.train()
    fallbacks: list[tuple[str, int, str]] = []
    checks = {"missing_mask_mismatch": 0, "locf_cap_violations": 0, "padding_violations": 0, "as_of_differences": 0,
              "lab_release_violations": 0, "filled_padding_or_gap": 0, "train_median_unavailable": 0}
    for ci, ch in enumerate(feats):
        check_time()
        c = cmap["channels"][ch]
        ob = vals.obs[ch]
        if c.table == "LABEVENTS" and len(ob):
            checks["lab_release_violations"] += int((np.abs(ob["t"].to_numpy() - ob["tc"].to_numpy() - LAB_RELEASE_H) > 1e-9).sum())
        S = bin_summaries(co, ob)
        meas = S["measured"] & valid
        measured[:, :, ci] = meas.astype(np.uint8)
        if c.indicator_only:
            missing[:, :, ci] = (~meas).astype(np.uint8)
        else:
            for k in ("median", "min", "max", "last"):
                arr[k][:, :, ci] = np.where(valid, S[k], np.nan).astype(np.float32)
            med = np.where(valid, S["median"], np.nan)
            missing[:, :, ci] = np.isnan(med).astype(np.uint8)
            checks["missing_mask_mismatch"] += int(((np.isnan(med) != (missing[:, :, ci] == 1))).sum())
        if c.value_class:
            arr["time_since_last"][:, :, ci] = time_since_last(co, S["last_t"], meas, c.cap_h).astype(np.float32)
            vi = vch.index(ch)
            med = np.where(valid, S["median"], np.nan)
            last = np.where(valid, S["last"], np.nan)
            filled, carried = locf_fill(med, last, c.cap_bins)
            # cap check: a carried cell must have an observed bin within the cap before it
            if carried.any():
                okc = np.zeros(med.shape, dtype=bool)
                for j in range(1, c.cap_bins + 1):
                    okc[:, j:] |= ~np.isnan(med[:, :-j])
                checks["locf_cap_violations"] += int((carried & ~okc).sum())
            lm_val, lm_src = landmark_train_medians(med, tr, co.out_h)
            for L, fs in fstores.items():
                rows, nb = fs["rows"], fs["nb"]
                pre = med[rows, :nb]
                m_L = lm_val[L]
                fs["median"][vi] = m_L
                fs["median_source"][vi] = lm_src[L]
                if lm_src[L] != "landmark":
                    fallbacks.append((ch, L, lm_src[L]))
                fv = filled[rows, :nb]
                cf = carried[rows, :nb]
                imp = np.isnan(fv)
                if not np.isfinite(m_L):
                    checks["train_median_unavailable"] += 1
                    imp = np.zeros_like(imp)
                fv = np.where(np.isnan(fv), m_L, fv)
                fs["value"][:, :, vi] = fv.astype(np.float32)
                fs["carried"][:, :, vi] = cf.astype(np.uint8)
                fs["imputed"][:, :, vi] = imp.astype(np.uint8)
                fs["missing"][:, :, vi] = np.isnan(pre).astype(np.uint8)
                checks["filled_padding_or_gap"] += int((~valid[rows, :nb]).sum())
        # padding rule: values empty, masks 1, flags 0
        pv = ~valid
        checks["padding_violations"] += int((measured[:, :, ci][pv] != 0).sum() + (missing[:, :, ci][pv] != 1).sum())
        if not c.indicator_only:
            checks["padding_violations"] += int((~np.isnan(arr["median"][:, :, ci][pv])).sum())
        # as-of rebuild at the core landmarks (check 13): truncating events at L leaves the prefix unchanged
        for L in CORE_LANDMARKS:
            nb = L // int(BIN_H)
            S2 = bin_summaries(co, ob, upto=float(L))
            a, b = S["median"][:, :nb], S2["median"][:, :nb]
            checks["as_of_differences"] += int((~((a == b) | (np.isnan(a) & np.isnan(b)))).sum())
            checks["as_of_differences"] += int((S["measured"][:, :nb] != S2["measured"][:, :nb]).sum())
    # panel masks, hemoglobin blood-gas own mask, source flags
    rsm = vals.rowstat[vals.rowstat["status"] == "measurement"]
    pnames = list(cmap["panels"])
    pm = np.ones((n, N_BINS, len(pnames)), dtype=np.uint8)
    for pi, pg in enumerate(pnames):
        d = rsm[rsm["item"].isin(cmap["panels"][pg])][["enc_i", "t", "v"]]
        S = bin_summaries(co, d.assign(v=0.0))
        pm[:, :, pi] = (~(S["measured"] & valid)).astype(np.uint8)
    raw.put("panel_missing", pm)
    raw.put("panels", np.array(pnames, dtype=str))
    # decision 4: unparsed qualitative text (value stays missing) flagged per bin; not a feature channel
    qch = qual_channels(cmap)
    qu = np.zeros((n, N_BINS, len(qch)), dtype=np.uint8)
    up = vals.unparsed
    for qi, ch_ in enumerate(qch):
        if up is None or not len(up):
            break
        d = up[up["channel"] == ch_]
        Su = bin_summaries(co, pd.DataFrame({"enc_i": d["enc_i"].to_numpy(dtype=np.int64),
                                             "t": d["t"].to_numpy(dtype=float), "v": 0.0}))
        qu[:, :, qi] = (Su["measured"] & valid).astype(np.uint8)
    raw.put("qual_unparsed", qu)
    raw.put("qual_channels", np.array(qch, dtype=str))
    hb = vals.obs["hemoglobin"]
    Sg = bin_summaries(co, hb[hb["item"] == "5257"])
    raw.put("hemoglobin_src_bloodgas", (Sg["measured"] & valid).astype(np.uint8))
    raw.put("hemoglobin_bloodgas_missing", (~(Sg["measured"] & valid)).astype(np.uint8))
    bi = vals.obs["bilirubin_total"]
    raw.put("bilirubin_total_src_5255", (bin_summaries(co, bi[bi["item"] == "5255"])["measured"] & valid).astype(np.uint8))
    # medications (state machine at bin end; Q10)
    mb = med_bins(co, docs)
    raw.put("med_groups", np.array(MED_GROUPS, dtype=str))
    for f in ("active", "started", "hours_on_course", "hours_since_stopped", "ever_received"):
        a = np.stack([mb["groups"][g][f] for g in MED_GROUPS], axis=2)
        raw.put(f"med_{f}", a.astype(np.float32 if f in ("hours_on_course", "hours_since_stopped") else np.uint8))
    raw.put("med_vasoactive_concurrent", mb["vaso_concurrent"].astype(np.uint8))
    act = np.stack([mb["groups"][g]["active"] for g in MED_GROUPS], axis=2)
    sta = np.stack([mb["groups"][g]["started"] for g in MED_GROUPS], axis=2)
    checks["started_not_active"] = int(((sta == 1) & (act == 0)).sum())
    checks["med_as_of_differences"] = 0
    for L in CORE_LANDMARKS:
        check_time()
        nb = L // int(BIN_H)
        mbL = med_bins(co, docs[docs["S"] <= L], as_of=float(L))
        for g in MED_GROUPS:
            for f in ("active", "started", "hours_on_course", "hours_since_stopped", "ever_received"):
                a, b = mb["groups"][g][f][:, :nb], mbL["groups"][g][f][:, :nb]
                checks["med_as_of_differences"] += int((~((a == b) | (np.isnan(a) & np.isnan(b)))).sum())
    # microbiology M1 (collection time)
    me = micro["enc_i"].to_numpy(dtype=np.int64)
    mc = np.zeros((n, N_BINS, len(SPECIMEN_GROUPS)), dtype=np.uint8)
    for si, sg in enumerate(SPECIMEN_GROUPS):
        d = micro[micro["group"] == sg]
        S = bin_summaries(co, pd.DataFrame({"enc_i": d["enc_i"].to_numpy(dtype=np.int64), "t": d["th"].to_numpy(),
                                            "v": 0.0}))
        mc[:, :, si] = (S["measured"] & valid).astype(np.uint8)
    raw.put("micro_collected", mc)
    raw.put("micro_any_collected", mc.max(axis=2))
    raw.put("specimen_groups", np.array(SPECIMEN_GROUPS, dtype=str))
    # surgery in the bin of T (Layer C surgeries: t0 + 6 h < T <= OUTTIME)
    lc = surg[surg["timing"] == "during_icu_layer_c"]
    sb = np.zeros((n, N_BINS), dtype=np.uint8)
    e, th = lc["enc_i"].to_numpy(dtype=np.int64), lc["T_h"].to_numpy()
    k = np.floor(th / BIN_H).astype(int)
    k = np.minimum(k, co.nvalid[e] - 1)
    ing = (th < MAX_H + 1e-9) & (k >= 0) & (k < N_BINS)
    np.maximum.at(sb, (e[ing], k[ing]), 1)
    raw.put("surgery_in_bin", sb)
    checks["layer_c_surgery_T_le_6h"] = int((th <= RECENT_HI_H).sum())
    checks["n_layer_c_surgeries"] = len(lc)
    checks["n_layer_c_surgeries_beyond_336h"] = int((~ing).sum())
    arr = measured = missing = None  # noqa: F841  (release the memmaps before zipping)
    raw.write_npz(pending / "stage3_layer_c_raw.npz")
    fs_src = {L: list(fs["median_source"]) for L, fs in fstores.items()}
    for L in list(fstores):
        fs = fstores.pop(L)
        st, med = fs["store"], fs["median"]
        fs.clear()
        st.put("train_median", med.astype(np.float64))
        st.put("train_median_source", np.array(fs_src.pop(L), dtype=str))
        st.write_npz(pending / f"stage3_layer_c_filled_L{L}.npz")
    return {"checks": checks, "courses": mb["courses"], "act": None, "median_fallbacks": fallbacks}


# ----------------------------------------------------------------------------------------------------------------------
# Aggregate audit document (spec §11) and release gate (spec §12; Stage 2 v1.0.3 module v1.5, hash-pinned)
# ----------------------------------------------------------------------------------------------------------------------
def _esc(k: str) -> str:
    return str(k).replace("~", "~0").replace("/", "~1")


class Agg:
    """Builds the audit JSON with its disclosure block: count fields (n_*), identities (total = sum of parts),
    hidden differences (hidden:*), ratios (pct_*), quantiles (q_*) and the frozen (already public) cells."""

    def __init__(self):
        self.doc: dict = {}
        self.count_fields: list[str] = []
        self.identities: list[tuple[str, list[str]]] = []
        self.hidden: list[str] = []
        self.ratios: list[list[str]] = []
        self.quantiles: list[list[str]] = []
        self.frozen: set[str] = set()

    def put(self, path: tuple, value: Any) -> str:
        d = self.doc
        for p in path[:-1]:
            d = d.setdefault(p, {})
        d[path[-1]] = value
        return "/" + "/".join(_esc(p) for p in path)

    def n(self, path: tuple, value: int, frozen: bool = False) -> str:
        if not str(path[-1]).startswith("n_"):
            raise GateError("internal: count key without the n_ prefix")
        ptr = self.put(path, int(value))
        self.count_fields.append(ptr)
        if frozen:
            self.frozen.add(ptr)
        return ptr

    def ident(self, total: str, parts: list[str]) -> None:
        self.identities.append((total, list(parts)))

    def hide(self, name: str) -> str:
        h = f"hidden:{name}"
        self.hidden.append(h)
        return h

    def pct(self, path: tuple, num: str, den: str, vals: dict) -> None:
        v = r1(100.0 * vals[num] / vals[den]) if vals[den] > 0 else "suppressed"
        ptr = self.put(path, v)
        self.ratios.append([ptr, num, den])

    def quant(self, path: tuple, values: np.ndarray, nptr: str) -> None:
        if len(values):
            q = np.percentile(values, [10, 25, 50, 75, 90])
            v = {k: float(round(float(x), 3)) for k, x in zip(("p10", "p25", "p50", "p75", "p90"), q)}
        else:
            v = "suppressed"
        ptr = self.put(path, v)
        self.quantiles.append([ptr, nptr])

    def values(self) -> dict:
        out = {}
        for ptr in self.count_fields:
            d: Any = self.doc
            for p in ptr.split("/")[1:]:
                d = d[p.replace("~1", "/").replace("~0", "~")]
            out[ptr] = d
        return out


def build_audit(co: Cohort, pt: pd.DataFrame, A_df: pd.DataFrame, ad: pd.DataFrame, a_info: dict, surg: pd.DataFrame,
                s_info: dict, vals: Values, cmap: dict, raw_info: dict, c_checks: dict) -> Agg:
    g = Agg()
    tot = g.n(("cohort", "n_encounters"), co.n, True)
    g.n(("cohort", "n_subjects"), co.enc["SUBJECT_ID"].nunique(), True)
    g.n(("cohort", "n_died"), int(co.enc["died"].sum()), True)
    sp = [g.n(("cohort", "split", s, "n_encounters"), int((co.split == s).sum()), True) for s in SPLITS]
    g.ident(tot, sp)
    for L in LANDMARKS:
        g.n(("landmarks", f"L{L:03d}", "n_eligible"), co.eligible[L], True)
    for name, info in raw_info.items():
        rt = g.n(("raw_rows", name, "n_rows"), info["rows"], True)
        rc = info.get("row_classes")
        if rc and name in ("CHARTEVENTS", "LABEVENTS", "PRESCRIPTIONS", "MICROBIOLOGYEVENTS"):
            a = g.n(("event_linkage", name, "n_rows_in_stay"), rc["in_stay"])
            b = g.n(("event_linkage", name, "n_rows_outside_every_window"), info["rows"] - rc["in_stay"])
            g.ident(rt, [a, b])
        icu = info.get("icustay")
        if icu:
            parts = [g.n(("icustay_id_check", name, k), v) for k, v in icu.items() if k != "n_rows_checked"]
            g.ident(g.n(("icustay_id_check", name, "n_rows_checked"), icu["n_rows_checked"]), parts)
    # age (decision 14; audit v1.1): only valid / invalid are released. The invalid detail (missing, unparseable,
    # negative, negative magnitude bins) is built so that it enters the joint audit as unknowns, and is withheld by the
    # hash-pinned drop list (/age/invalid_detail; spec §16).
    stt = pt["age_status"].to_numpy()
    if not set(np.unique(stt)) <= {"valid", "missing", "unparseable", "negative"}:
        raise GateError("internal: unexpected age status class")
    a_valid = g.n(("age", "n_valid"), int((stt == "valid").sum()))
    a_inv = g.n(("age", "n_invalid"), int((stt != "valid").sum()))
    g.ident(tot, [a_valid, a_inv])
    det = [g.n(("age", "invalid_detail", f"n_{c}"), int((stt == c).sum())) for c in ("missing", "unparseable", "negative")]
    g.ident(a_inv, det)
    nb = pt["age_neg_bin"].to_numpy()
    g.ident(det[2], [g.n(("age", "invalid_detail", "negative_magnitude", f"n_{b}"), int((nb == b).sum()))
                     for b in ("lt1d", "1to30d", "gt30d")])
    # hospital-to-ICU interval (Q4; audit v1.1): the clipped count is released; non-negative and missing are withheld
    # as a pair (/admission_interval/withheld): with the frozen cohort total either one gives the other.
    clipped = A_df["admit_interval_clipped"].to_numpy() == 1
    miss = A_df["admit_interval_missing"].to_numpy() == 1
    g.ident(tot, [g.n(("admission_interval", "n_negative_clipped_to_0"), int(clipped.sum())),
                  g.n(("admission_interval", "withheld", "n_nonnegative"), int((~clipped & ~miss).sum())),
                  g.n(("admission_interval", "withheld", "n_missing"), int(miss.sum()))])
    # EDREGTIME placeholder audit (Q6) and VISIT_ID agreement (Q7): aggregates only, not features
    edp = ad["ED"].notna().to_numpy()
    eq = edp & (ad["ED"].to_numpy() == ad["ADMITTIME"].to_numpy())
    pres = g.n(("edregtime", "n_present"), int(edp.sum()))
    e1 = g.n(("edregtime", "n_equals_admittime"), int(eq.sum()))
    e2 = g.n(("edregtime", "n_differs_from_admittime"), int((edp & ~eq).sum()))
    g.ident(pres, [e1, e2])
    g.ident(tot, [pres, g.n(("edregtime", "n_missing"), int((~edp).sum()))])
    vv = g.values()
    g.pct(("edregtime", "pct_equals_admittime_of_present"), e1, pres, vv)
    vi = s_info["visit"]
    vt = "/raw_rows/SURGERY_INFO/n_rows"
    g.ident(vt, [g.n(("visit_id", "SURGERY_INFO", "n_visit_equals_hadm"), vi["n_visit_equals_hadm"]),
                 g.n(("visit_id", "SURGERY_INFO", "n_visit_differs"), vi["n_visit_differs"]),
                 g.n(("visit_id", "SURGERY_INFO", "n_visit_blank"), vi["n_visit_blank"])])
    # diagnosis coverage (Q8)
    # diagnosis coverage (Q8; decision 8, Nov 17, 2026: ICD-O morphology codes -> oncology group, own count; 'missing'
    # = no code, separate from 'unmapped')
    c3 = A_df["dx_icd10_3char"].to_numpy()
    mapped = A_df["dx_mapped"].to_numpy() == 1
    dmiss = c3 == "missing"
    icdo = c3 == ICDO_LEVEL
    unp = c3 == "unmapped"
    unc = ~mapped & ~dmiss & ~unp & ~icdo
    p_unc = g.n(("diagnosis", "n_unmapped_with_chapter"), int(unc.sum()))
    p_unp = g.n(("diagnosis", "n_unmapped_unparsed"), int(unp.sum()))
    g.ident(tot, [g.n(("diagnosis", "n_mapped"), int(mapped.sum())), g.n(("diagnosis", "n_missing"), int(dmiss.sum())),
                  g.n(("diagnosis", "n_icdo_morphology_oncology"), int(icdo.sum())), p_unc, p_unp])
    unm = ~mapped & ~dmiss & ~icdo
    raw_codes = pd.Series([norm_icd(x)[0] for x in ad["dx"].to_numpy()])[unm]
    vc = raw_codes.value_counts()
    big = vc[vc >= MIN_CELL]
    listed = [g.n(("diagnosis", "unmapped_codes", str(k) if k else "(blank)", "n_encounters"), int(v)) for k, v in big.items()]
    pooled = g.n(("diagnosis", "n_unmapped_other_codes_pooled"), int(vc[vc < MIN_CELL].sum()))
    un_tot = g.n(("diagnosis", "n_unmapped_total"), int(unm.sum()))
    g.ident(un_tot, listed + [pooled])
    g.ident(un_tot, [p_unc, p_unp])
    # surgery (Q12): timing, duration source, fallback reasons, band reference level
    lk = g.n(("surgery", "n_linked"), len(surg))
    tm = surg["timing"].to_numpy()
    g.ident(lk, [g.n(("surgery", "timing", "n_recent_layer_a"), int((tm == "recent_layer_a").sum())),
                 g.n(("surgery", "timing", "n_remote_same_hadm"), int((tm == "remote_same_hadm").sum())),
                 g.n(("surgery", "timing", "n_layer_c"), int((tm == "during_icu_layer_c").sum()), True)])
    ds = surg["duration_source"].to_numpy()
    fbr = surg["fallback_reason"].to_numpy()
    fb = g.n(("surgery", "duration", "n_anesthesia_window_fallback"), int((ds == "anesthesia_window").sum()))
    g.ident(lk, [g.n(("surgery", "duration", "n_surgery_window"), int((ds == "surgery_window").sum())), fb,
                 g.n(("surgery", "duration", "n_band_missing_flagged"), int((ds == "none").sum()))])
    sel = ds == "anesthesia_window"
    g.ident(fb, [g.n(("surgery", "duration", "fallback_reason", f"n_{r}"), int((sel & (fbr == r)).sum()))
                 for r in ("gt24", "order", "missing")])
    rl = surg["band_ref_level"].to_numpy()
    g.ident(lk, [g.n(("surgery", "band_reference", f"n_{r}"), int((rl == r).sum()))
                 for r in ("own_type", "system_x_group", "system", "all")])
    # weight source (audit only; B38). Released: charted and imputed (imputed = cohort - charted, as public as the
    # charted count). The per-source split (cn_ref, who_ref, the four train tiers) is built for the joint audit and
    # withheld by the drop list (/weight_source/by_source): every age_invalid encounter without charted weight falls to
    # the train tiers by construction, so a released split would give the age_invalid-with-charted-weight count by
    # subtraction from /age/n_invalid whenever the valid-age tier count is known from outside (spec §10, §16.4).
    ws = vals.weight_source
    w_ch = g.n(("weight_source", "n_charted"), int((ws == "charted").sum()))
    w_im = g.n(("weight_source", "n_imputed"), int((ws != "charted").sum()))
    g.ident(tot, [w_ch, w_im])
    g.ident(w_im, [g.n(("weight_source", "by_source", f"n_{s}"), int((ws == s).sum())) for s in WEIGHT_SOURCES[1:]])
    # plausibility gates (Q13): train and overall only; validation + test as hidden differences
    rsx = vals.rowstat
    e = rsx["enc_i"].to_numpy(dtype=np.int64)
    ing = (rsx["t"].to_numpy() >= 0) & (rsx["t"].to_numpy() < co.W[e])
    rsx = rsx[ing]
    trm = co.split[rsx["enc_i"].to_numpy(dtype=np.int64)] == "train"
    for ch in cmap["feature_channels"]:
        c = cmap["channels"][ch]
        its = [i for (t, i), inf in cmap["items"].items() if inf["channel"] == ch]
        gated = any(cmap["items"][(c.table, i)]["low"] is not None or cmap["items"][(c.table, i)]["high"] is not None
                    for i in its) if c.table != "DERIVED" else False
        sub = rsx["channel"].to_numpy() == ch
        stv = rsx["status"].to_numpy()
        numeric = sub & np.isin(stv, ["measurement", "out_of_gate", "dedupe_same_draw", "dropped_non_nicu",
                                      "zero_missing"])
        if c.qualitative:
            vq = rsx["v"].to_numpy()
            base = ("plausibility", ch, "train")
            parts = [g.n(base + ("n_rows_positive",), int((sub & trm & (stv != "qual_unparsed") & (vq == 1)).sum())),
                     g.n(base + ("n_rows_negative",), int((sub & trm & (stv != "qual_unparsed") & (vq == 0)).sum())),
                     g.n(base + ("n_rows_unparsed",), int((sub & trm & (stv == "qual_unparsed")).sum()))]
            g.ident(g.n(base + ("n_rows_text",), sum(g.values()[p] for p in parts)), parts)
            continue
        if gated:
            # decision 2 (Nov 17, 2026): numeric = in gate + out of gate (+ zeros read as missing for ZERO_AS_MISSING);
            # out of gate = below + above at the overall scope for GATE_SIDE_SPLIT (train keeps the total only)
            zch = ch in ZERO_AS_MISSING
            names = ("in_gate", "out_of_gate", "numeric") + (("zero_as_missing",) if zch else ())
            side = rsx["side"].to_numpy() if "side" in rsx else np.zeros(len(rsx), dtype=np.int8)
            oog = numeric & (stv == "out_of_gate")
            cells = {}
            for scope, msk in (("train", trm), ("overall", np.ones(len(rsx), dtype=bool))):
                a = g.n(("plausibility", ch, scope, "n_rows_in_gate"),
                        int((numeric & msk & ~np.isin(stv, ["out_of_gate", "zero_missing"])).sum()))
                b = g.n(("plausibility", ch, scope, "n_rows_out_of_gate"), int((oog & msk).sum()))
                t_ = g.n(("plausibility", ch, scope, "n_rows_numeric"), int((numeric & msk).sum()))
                parts_ = [a, b]
                cl = [a, b, t_]
                if zch:
                    z = g.n(("plausibility", ch, scope, "n_rows_zero_as_missing"),
                            int((numeric & msk & (stv == "zero_missing")).sum()))
                    parts_.append(z)
                    cl.append(z)
                g.ident(t_, parts_)
                if ch in GATE_SIDE_SPLIT and scope == "overall" and int((oog & msk).sum()) >= GATE_SIDE_MIN_TOTAL:
                    # both identities are declared: out = below + above, and numeric = in + below + above (+ zeros),
                    # so a small below/above cell has publishable neighbours (in gate, numeric) to protect it with
                    lo_ = g.n(("plausibility", ch, scope, "n_rows_below_gate"), int((oog & (side < 0)).sum()))
                    hi_ = g.n(("plausibility", ch, scope, "n_rows_above_gate"), int((oog & (side > 0)).sum()))
                    g.ident(b, [lo_, hi_])
                    g.ident(t_, [a, lo_, hi_] + parts_[2:])
                cells[scope] = cl
            for j, nm in enumerate(names):
                g.ident(cells["overall"][j], [cells["train"][j], g.hide(f"{ch}_{nm}_validation_plus_test")])
        elif c.value_class or c.triggered:
            ob = vals.obs[ch]
            oe = ob["enc_i"].to_numpy(dtype=np.int64)
            sel = in_grid(co, oe, ob["t"].to_numpy()) & (co.split[oe] == "train")
            v = ob["v"].to_numpy()[sel]
            v = v[np.isfinite(v)]
            np_ = g.n(("plausibility", ch, "train_ungated", "n_values"), len(v))
            g.quant(("plausibility", ch, "train_ungated", "q_values"), v if len(v) >= MIN_CELL else np.zeros(0), np_)
        conv = [cmap["items"][(c.table, i)]["conv"] for i in its] if c.table != "DERIVED" else []
        if any(CONVERSIONS.get(x, (None, 1))[0] is not None for x in conv if x):
            g.n(("unit_rows", ch, "n_rows_unit_dropped"), int((sub & (stv == "unit_dropped")).sum()))
            g.n(("unit_rows", ch, "n_rows_unit_unrecognized"), int((sub & (stv == "unit_unrecognized")).sum()))
    return g


def render_markdown(doc: dict) -> str:
    L = [f"# Stage 3 real-layer builder: aggregate audit", "",
         f"- spec: {doc.get('spec')}", f"- release gate: {doc.get('release_gate')}",
         "- every count is >= 10, '<10' (0-9) or 'S' (complementary suppression); percentages and quantiles marked "
         "'suppressed' rest on a suppressed count", ""]

    def walk(o: Any, path: str) -> None:
        if isinstance(o, dict):
            if set(o) == {"p10", "p25", "p50", "p75", "p90"}:
                L.append(f"- `{path}`: " + ", ".join(f"{k} {o[k]}" for k in ("p10", "p25", "p50", "p75", "p90")))
                return
            for k, v in o.items():
                walk(v, f"{path}/{k}")
        else:
            L.append(f"- `{path}`: {o}")
    for top, v in doc.items():
        if top in ("disclosure", "input_hashes", "maps"):
            continue
        L.append(f"## {top}")
        L.append("")
        walk(v, top)
        L.append("")
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------------------------------------------------
# Cross-release (joint) disclosure audit with the published Stage 2 release (v1.0.2; spec §16; v1.0.4: Stage 2 v1.0.3)
# ----------------------------------------------------------------------------------------------------------------------
@dataclass
class CrossRelease:
    doc: dict                 # the Stage 2 v1.0.3 release JSON (hash-gated at the precheck)
    sha256: str
    links: list               # (kind 'eq'|'le', total, parts, note); pointers 's2:/...' are Stage 2 cells
    drops: list               # (subtree, reason) of the hash-pinned drop list
    constants: list = field(default_factory=list)   # v1.0.4: (s2 pointer, int | 'derived', row); team-internal


def json_ptr_get(doc: Any, ptr: str) -> Any:
    d = doc
    for k in ptr.split("/")[1:]:
        k = k.replace("~1", "/").replace("~0", "~")
        if not isinstance(d, dict) or k not in d:
            raise KeyError(ptr)
        d = d[k]
    return d


def read_links_v2(path: Path) -> list[tuple[str, str, tuple[str, ...], str]]:
    """Links file v2 (spec §16.2). CSV columns kind,total,parts,note; kind 'eq' (total = sum of parts) or 'le'
    (total >= sum of parts, entered with a hidden slack cell); parts ';'-separated; '#' lines are comments.
    Pointers start with '/' (Stage 3 audit count cell) or 's2:/' (Stage 2 release cell). v1.0.4: links file v3 has
    the same format (rows whose pointers are all 's2:' are Stage 2 internal equalities)."""
    df = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    need = {"kind", "total", "parts", "note"}
    if not need.issubset(df.columns):
        raise InputError(f"links file {Path(path).name} needs columns {sorted(need)} (links file v2)")
    out, seen = [], set()
    for i, row in df.iterrows():
        kind, t = row["kind"].strip(), row["total"].strip()
        parts = tuple(q.strip() for q in row["parts"].split(";") if q.strip())
        ok = lambda c: c.startswith("/") or c.startswith(XR + "/")   # noqa: E731
        if kind not in ("eq", "le") or not parts or not ok(t) or not all(ok(q) for q in parts) or t in parts:
            raise InputError(f"links file {Path(path).name}, row {i + 2}: needs kind eq|le and pointers starting "
                             f"with '/' or '{XR}/'")
        key = (kind, t, tuple(sorted(parts)))
        if key in seen:
            raise InputError(f"links file {Path(path).name}, row {i + 2}: duplicate link")
        seen.add(key)
        out.append((kind, t, parts, row["note"].strip()))
    return out


def load_cross_release(cfg: Config, SM, sha: str) -> CrossRelease:
    """Read the hash-gated Stage 2 release, the links file v2 and the drop list; every check fails closed as an
    input error (exit 2) before any data file is read."""
    path = Path(cfg.stage2_release_path)
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        raise InputError("Stage 2 release JSON could not be parsed") from None
    if not isinstance(doc, dict) or doc.get("release_gate") != "ok":
        raise InputError("Stage 2 release JSON does not carry release_gate 'ok'")
    dis = doc.get("disclosure")
    if not isinstance(dis, dict) or not all(isinstance(dis.get(k), list)
                                            for k in ("count_fields", "identities", "hidden_cells")):
        raise InputError("Stage 2 release JSON has no disclosure registry (count_fields, identities, hidden_cells)")
    links = read_links_v2(Path(cfg.links_path)) if cfg.links_path is not None else []
    if not links:
        raise InputError("links file lists no cross-release link")
    s2c = set(dis["count_fields"])
    for _k, t, parts, _n in links:
        for c in (t,) + parts:
            if c.startswith(XR):
                ptr = c[len(XR):]
                if ptr not in s2c:
                    raise InputError(f"links file Stage 2 pointer is not a count field of the Stage 2 release: {c}")
                try:
                    json_ptr_get(doc, ptr)
                except KeyError:
                    raise InputError(f"links file Stage 2 pointer not found in the Stage 2 release: {c}") from None
    drops = SM.read_drops(Path(cfg.drops_path)) if cfg.drops_path is not None else []
    consts = read_constants(Path(cfg.stage2_constants_path), doc) if cfg.stage2_constants_path is not None else []
    n_s2 = sum(1 for _k, t, parts, _n in links if all(c.startswith(XR) for c in (t,) + parts))
    log(f"cross-release inputs: {STAGE2_RELEASE_LABEL} OK ({len(dis['count_fields'])} count fields), "
        f"{len(links)} declared links ({n_s2} Stage 2 internal), {len(drops)} dropped subtrees, "
        f"{len(consts)} team-internal Stage 2 constants ({sum(1 for _p, v, _r in consts if v == 'derived')} derived)")
    return CrossRelease(doc=doc, sha256=sha, links=links, drops=drops, constants=consts)


def read_constants(path: Path, doc: dict) -> list[tuple[str, Any, int]]:
    """v1.0.4: TEAM-INTERNAL constants file (CSV pointer,value,reason; '#' lines are comments): Stage 2 cells a reader
    already knows from outside the Stage 2 release (value = a non-negative integer) or that those values fix exactly
    through the Stage 2 identities and the Stage 2 internal links ('derived'). Pointers 's2:/...' must be count cells
    of the Stage 2 release (not hidden). An integer must agree with the release: equal to a published value, 0-9 for
    a '<10' cell, >= 10 for an 'S' cell. Fail closed (exit 2); messages name rows only, never values."""
    try:
        df = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    except Exception:  # noqa: BLE001
        raise InputError(f"constants file {path.name} could not be read") from None
    if not {"pointer", "value", "reason"}.issubset(df.columns):
        raise InputError(f"constants file {path.name} needs columns ['pointer', 'value', 'reason']")
    dis = doc["disclosure"]
    s2c, s2h = set(dis["count_fields"]), set(dis["hidden_cells"])
    out, seen = [], set()
    for i, r in df.iterrows():
        row = i + 2
        ptr, val, why = r["pointer"].strip(), r["value"].strip(), r["reason"].strip()
        if not ptr.startswith(XR + "/") or not why:
            raise InputError(f"constants file row {row}: needs a pointer starting with '{XR}/' and a reason")
        if ptr in seen:
            raise InputError(f"constants file row {row}: duplicate pointer")
        seen.add(ptr)
        q = ptr[len(XR):]
        if q not in s2c or q in s2h:
            raise InputError(f"constants file row {row}: pointer is not a count cell of the Stage 2 release")
        if re.fullmatch(r"\d+", val):
            v: Any = int(val)
            lab = json_ptr_get(doc, q)
            if isinstance(lab, int) and not isinstance(lab, bool):
                if lab != v:
                    raise InputError(f"constants file row {row}: value contradicts the published Stage 2 value "
                                     f"(nothing printed)")
            elif lab == "<10":
                if v >= MIN_CELL:
                    raise InputError(f"constants file row {row}: value contradicts the Stage 2 '<10' label")
            elif lab == "S":
                if v < MIN_CELL:
                    raise InputError(f"constants file row {row}: value contradicts the Stage 2 'S' label")
            else:
                raise InputError(f"constants file row {row}: unexpected Stage 2 count value")
        elif val == "derived":
            v = "derived"
        else:
            raise InputError(f"constants file row {row}: value must be a non-negative integer or 'derived'")
        out.append((ptr, v, row))
    if not out:
        raise InputError(f"constants file {path.name} lists no constant")
    return out


def derive_constants(xr: "CrossRelease") -> dict[str, int]:
    """v1.0.4: values of the constants file: integers as given; each 'derived' cell fixed exactly by the Stage 2
    identities and the Stage 2 internal links (all pointers 's2:') over published and constant cells: exact LP
    min = max over the cell's connected block of equations, every unknown free (no bounds), as the Stage 2 v1.0.3
    module derives its public cells. A derived value must agree with the release label. Fail closed (exit 2); values
    are never printed."""
    from scipy import sparse as _sp
    from scipy.optimize import linprog as _lp
    dis = xr.doc["disclosure"]
    pub: dict[str, int] = {}
    for q in dis["count_fields"]:
        v = json_ptr_get(xr.doc, q)
        if isinstance(v, int) and not isinstance(v, bool):
            pub[XR + q] = v
    known = dict(pub)
    for ptr, v, _row in xr.constants:
        if v != "derived":
            known[ptr] = v
    eqs = [(XR + t, tuple(XR + x for x in ps)) for t, ps in dis["identities"]]
    eqs += [(t, ps) for k, t, ps, _n in xr.links if k == "eq" and all(c.startswith(XR) for c in (t,) + ps)]
    out = {ptr: v for ptr, v, _r in xr.constants if v != "derived"}
    todo = [(ptr, row) for ptr, v, row in xr.constants if v == "derived"]
    for _ in range(len(todo) + 1):
        if not todo:
            break
        by_cell: dict[str, list[int]] = defaultdict(list)
        for i, (t, ps) in enumerate(eqs):
            for c in (t,) + ps:
                if c not in known:
                    by_cell[c].append(i)
        left = []
        for ptr, row in todo:
            if ptr in known:                      # published in the release: its value is known already
                out[ptr] = known[ptr]
                continue
            rows_, cells, stack = set(), {ptr}, [ptr]
            while stack:
                for i in by_cell.get(stack.pop(), []):
                    if i in rows_:
                        continue
                    rows_.add(i)
                    t, ps = eqs[i]
                    for q in (t,) + ps:
                        if q not in known and q not in cells:
                            cells.add(q)
                            stack.append(q)
            if not rows_:
                left.append((ptr, row))
                continue
            cl = sorted(cells)
            col = {q: j for j, q in enumerate(cl)}
            coef: dict = defaultdict(float)
            b = np.zeros(len(rows_))
            for k, i in enumerate(sorted(rows_)):
                t, ps = eqs[i]
                for q, a in [(t, 1.0)] + [(q, -1.0) for q in ps]:
                    if q in known:
                        b[k] -= a * known[q]
                    else:
                        coef[(k, col[q])] += a
            kk = list(coef)
            A = _sp.csr_matrix(([coef[x] for x in kk], ([x[0] for x in kk], [x[1] for x in kk])),
                               shape=(len(rows_), len(cl)))
            e = np.zeros(len(cl))
            e[col[ptr]] = 1.0
            free = [(None, None)] * len(cl)
            lo = _lp(e, A_eq=A, b_eq=b, bounds=free, method="highs")
            hi = _lp(-e, A_eq=A, b_eq=b, bounds=free, method="highs")
            if lo.status != 0 or hi.status != 0 or abs(lo.x[col[ptr]] - hi.x[col[ptr]]) > 1e-6:
                left.append((ptr, row))
                continue
            v = int(round(lo.x[col[ptr]]))
            lab = json_ptr_get(xr.doc, ptr[len(XR):])
            if v < 0 or (lab == "<10" and v >= MIN_CELL) or (lab == "S" and v < MIN_CELL):
                raise InputError(f"constants file row {row}: derived value contradicts the Stage 2 label "
                                 f"(nothing printed)")
            known[ptr] = v
            out[ptr] = v
        if len(left) == len(todo):
            todo = left
            break
        todo = left
    if todo:
        raise InputError(f"{len(todo)} 'derived' constant(s) are not fixed by the Stage 2 identities and links over "
                         f"published and constant cells (first row {min(r for _p, r in todo)})")
    return out


class _Chain:
    """Status view: the caller's status first, then the fixed Stage 2 / slack statuses (disjoint keys)."""
    __slots__ = ("a", "b")

    def __init__(self, a, b):
        self.a, self.b = a, b

    def __getitem__(self, k):
        try:
            return self.a[k]
        except KeyError:
            return self.b[k]

    def get(self, k, d=None):
        try:
            return self[k]
        except KeyError:
            return d

    def __contains__(self, k):
        return self.get(k) is not None


def joint_parts(SM, L, xr: CrossRelease) -> tuple[list, dict, dict, list]:
    """Stage 2 cells and relations of the joint system: every declared link, plus the Stage 2 identities reachable
    from the linked Stage 2 cells through Stage 2 unknowns. Returns (relations, extra_status, extra_values,
    link_rels). Stage 2 published cells are 'const' (known, never suppressible); '<10' -> P, 'S' -> S, hidden -> H;
    'le' links get a hidden slack cell."""
    dis = xr.doc["disclosure"]
    s2_count, s2_hidden = set(dis["count_fields"]), set(dis["hidden_cells"])
    est: dict[str, str] = {}
    evl: dict[str, int] = {}
    cvals = derive_constants(xr) if xr.constants else {}    # v1.0.4: team-internal constants (never printed)

    def cell(ptr: str) -> str:
        key = XR + ptr
        if key in est:
            return est[key]
        if key in cvals:
            st, evl[key] = "const", int(cvals[key])
        elif ptr in s2_hidden:
            st = "H"
        elif ptr in s2_count:
            v = json_ptr_get(xr.doc, ptr)
            if isinstance(v, int) and not isinstance(v, bool):
                st, evl[key] = "const", v
            elif v == "<10":
                st = "P"
            elif v == "S":
                st = "S"
            else:
                raise InputError(f"Stage 2 release: unexpected count value at {ptr}")
        else:
            st = "H"     # used by a Stage 2 identity but not released: an unknown, a reader knows only >= 0
        est[key] = st
        return st
    ids = dis["identities"]
    by: dict[str, list[int]] = defaultdict(list)
    for i, (t, parts) in enumerate(ids):
        for c in [t] + list(parts):
            by[c].append(i)
    seeds = sorted({c[len(XR):] for _k, t, parts, _n in xr.links for c in (t,) + parts if c.startswith(XR)}
                   | {c[len(XR):] for c in cvals})
    seen: set[str] = set()
    take: set[int] = set()
    stack = list(seeds)
    while stack:
        c = stack.pop()
        if c in seen:
            continue
        seen.add(c)
        cell(c)
        for i in by.get(c, ()):
            if i in take:
                continue
            take.add(i)
            t, parts = ids[i]
            for x in [t] + list(parts):
                if cell(x) in SM.UNKNOWN and x not in seen:
                    stack.append(x)
    rels = [(XR + ids[i][0], tuple(XR + x for x in ids[i][1]), "stage2 identity") for i in sorted(take)]
    link_rels = []
    for j, (kind, t, parts, _note) in enumerate(xr.links):
        for c in (t,) + parts:
            if not c.startswith(XR) and L.status.get(c) not in ("pub", "P", "D"):
                raise InputError(f"links file Stage 3 pointer is not a count cell of the Stage 3 audit: {c}")
        if kind == "le":
            sl = f"{XR}hidden:link_slack_{j:03d}"
            est[sl] = "H"
            parts = parts + (sl,)
        link_rels.append((t, parts, "cross-release link"))
    return rels, est, evl, link_rels


def joint_system(SM, relations, values, kinds, params, extra: dict, tolerated: dict):
    """The pinned module's System with the fixed Stage 2 statuses overlaid, and the compromised-cell rule: a cell
    that already fails the protection rule on public information alone (Stage 2 v1.0.3 + constants + frozen public totals, every
    other Stage 3 cell unknown) is 'tolerated' with its baseline integer interval; it fails the joint audit only if
    Stage 3 narrows it further (or adds a residual failure the baseline did not have)."""

    def tighter(st, c: str, lo: float, hi: float) -> bool:
        lo0, hi0, _why0 = tolerated[c]
        a_lo, a_hi = SM.int_bounds(*SM.apriori(st[c], params))   # what the published label alone tells a reader
        lo0, hi0 = max(lo0, a_lo), min(hi0, a_hi)
        li, hi_i = SM.int_bounds(lo, hi)
        return li > lo0 or hi_i < hi0

    class JointSystem(SM.System):
        def audit(self, status, rel_ids=None, label=None):
            st = _Chain(status, extra)
            res = SM.System.audit(self, st, rel_ids, label)
            gone = []
            for c in sorted(c for c in res.unsafe if c in tolerated):
                why = res.unsafe[c]
                if "residual" in why and "residual" not in tolerated[c][2] and tolerated[c][0] != tolerated[c][1]:
                    continue        # a residual failure the baseline did not have (unless the value is public)
                rels = SM.System.component_rels(self, [c], st)
                lo, hi = SM.System.target_bounds(self, st, rels, [c])[c] if rels else SM.apriori(st[c], self.params)
                if not tighter(st, c, lo, hi):
                    gone.append(c)
            if gone:
                for c in gone:
                    del res.unsafe[c]
                res.kind_fail = {}
                for c in res.unsafe:
                    res.kind_fail[st[c]] = res.kind_fail.get(st[c], 0) + 1
                res.residual_groups = [k for k in res.residual_groups
                                       if any(m in res.unsafe for m, _ in self.rels[k])]
            return res

        def target_bounds(self, status, rel_ids, targets, lp_cap=None):
            st = _Chain(status, extra)
            out = SM.System.target_bounds(self, st, rel_ids, targets, lp_cap=lp_cap)
            for c in targets:
                if out.get(c) is None:          # v1.5: LP hit its cap -> unknown, never treated as safe
                    continue
                if c in tolerated and not tighter(st, c, *out[c]):
                    out[c] = SM.apriori(st[c], self.params)   # no narrowing beyond public information
            return out

        def local_rels(self, seeds, status, *a, **k):     # v1.5 (module v1.5 candidate checks)
            return SM.System.local_rels(self, seeds, _Chain(status, extra), *a, **k)

        def residual_fails(self, status, cell):
            st = _Chain(status, extra)
            if cell in tolerated and ("residual" in tolerated[cell][2] or tolerated[cell][0] == tolerated[cell][1]):
                return False
            return SM.System.residual_fails(self, st, cell)

        def component_rels(self, seeds, status):
            return SM.System.component_rels(self, seeds, _Chain(status, extra))

        def _solve_components(self, comps, status, label):
            """v1.3 screening proved a cell safe from the range of a few feasible points. For H/D cells the rule
            ('exactly 0' is safe) is not monotone: points that all sit at 0 do not prove the cell is 0 (its true
            range may be 0..8). Such cells get their exact maximum here (one HiGHS LP each; spec §16.6). Module v1.5
            has the same fix natively; kept as a second, independent check (v1.0.4)."""
            out = SM.System._solve_components(self, comps, status, label)
            for crels in comps:
                sus = sorted({v for u, _ in crels for v in u
                              if status[v] in ("H", "D") and v in out and out[v][1] <= SM.EPS})
                if not sus:
                    continue
                vars_ = sorted({v for u, _ in crels for v in u})
                col = {v: j for j, v in enumerate(vars_)}
                r_i, c_i, dat, rhs = [], [], [], []
                for r, (u, b) in enumerate(crels):
                    for v, a in u.items():
                        r_i.append(r)
                        c_i.append(col[v])
                        dat.append(a)
                    rhs.append(b)
                from scipy import sparse as _sp
                from scipy.optimize import linprog as _lp
                A = _sp.csr_matrix((dat, (r_i, c_i)), shape=(len(crels), len(vars_)))
                bnd = [SM.apriori(status[v], self.params) for v in vars_]
                bnd = [(lo, SM.BIG if hi == SM.INF else hi) for lo, hi in bnd]
                for v in sus:
                    e = np.zeros(len(vars_))
                    e[col[v]] = -1.0
                    SM.check_time()
                    self.lp_solves += 1
                    r = _lp(e, A_eq=A, b_eq=np.array(rhs, dtype=float), bounds=bnd, method="highs")
                    if r.status != 0:
                        raise GateError(f"joint audit: exact bound LP failed (status {r.status})")
                    hi = float(r.x[col[v]])
                    out[v] = (out[v][0], SM.INF if hi >= SM.BIG / 2 else hi)
            return out

    return JointSystem(relations, values, kinds, params)


def compromised_baseline(SM, relations, values, kinds, params, extra, L, frozen) -> dict:
    """Public-information baseline: Stage 2 v1.0.3 as released (+ team-internal constants), frozen (already public) Stage 3 cells known, every
    other Stage 3 count cell unknown (D). Cells failing the rule here cannot be protected by Stage 3; returns
    {cell: (lo, hi, reason)} with integer bounds."""
    base = {}
    for c, st in L.status.items():
        base[c] = ("pub" if c in frozen else "D") if st in ("pub", "P", "S", "D") else st
    plain = joint_system(SM, relations, values, kinds, params, extra, {})
    res = plain.audit(base, label="baseline audit (Stage 2 v1.0.3 + constants + public totals, Stage 3 withheld)")
    out = {}
    st = _Chain(base, extra)
    for c, why in sorted(res.unsafe.items()):
        rels = SM.System.component_rels(plain, [c], st)
        lo, hi = SM.System.target_bounds(plain, st, rels, [c])[c] if rels else SM.apriori(st[c], params)
        li, hi_i = SM.int_bounds(lo, hi)
        out[c] = (li, hi_i, why)
    # Stage 3 cells (>= 10, not withheld) that public information already pins: exactly -> treated as frozen (published,
    # never suppressed); to within the protection width -> tolerated with that interval if they are suppressed.
    seeds = [c for c in L.status if not c.startswith(XR) and any(
        c == t or c in parts for t, parts, src in relations if src == "cross-release link")]
    rels = SM.System.component_rels(plain, seeds, st) if seeds else set()
    cand = sorted({m for k in rels for m, _ in plain.rels[k]
                   if not m.startswith(XR) and L.status.get(m) == "pub" and m not in frozen})
    public = set()
    if cand:
        # pub cells are known in the baseline status; ask for their bounds with each one set unknown (D), batched
        st_c = _Chain(SM.Overlay(base, {c: "D" for c in cand}), extra)
        rels_c = SM.System.component_rels(plain, cand, st_c)
        bnd = SM.System.target_bounds(plain, st_c, rels_c, cand)
        for c in cand:
            li, hi_i = SM.int_bounds(*bnd[c])
            if li == hi_i:
                public.add(c)                       # exactly public already: published, never suppressed
            elif not SM.is_safe("S", *bnd[c], params):
                out[c] = (li, hi_i, "interval")     # if suppressed, it may not be narrowed beyond this interval
    return out, public


def release_gate(SM, g: Agg, links_path: Path | None, base: dict, pending: Path,
                 xr: CrossRelease | None = None, max_new: int = 2500,
                 max_component: int = 20000) -> tuple[dict, str, dict]:
    """Complementary suppression and the release gate, both from the hash-pinned Stage 2 v1.0.3 module (v1.5).
    v1.0.2: with `xr` (always in a real build) the audit is joint with the published Stage 2 release (v1.0.4: v1.0.3; spec §16)
    and the drop list is applied. Errors raised inside the module are re-raised as this script's GateError /
    InputError (fail closed); with `xr`, published cells contradicting a declared link fail the gate (exit 3)."""
    try:
        return _release_gate(SM, g, links_path, base, pending, xr, max_new, max_component)
    except (GateError, InputError, TimeLimit):
        raise
    except SM.TimeLimitError:
        raise TimeLimit() from None
    except SM.InputError as e:
        if xr is not None:
            raise GateError(f"cross-release: published Stage 2 v1.0.3 and Stage 3 cells are inconsistent under the "
                            f"declared links or identities (suppression module: {e}); check stage3_links_v3.csv and "
                            f"the constants file") from None
        raise InputError(f"suppression module: {e}") from None
    except SM.GateError as e:
        raise GateError(f"suppression module: {e}") from None


def _release_gate(SM, g: Agg, links_path: Path | None, base: dict, pending: Path,
                  xr: CrossRelease | None = None, max_new: int = 2500,
                  max_component: int = 20000) -> tuple[dict, str, dict]:
    doc = dict(base)
    for k, v in g.doc.items():
        doc[k] = v
    doc["suppression"] = {"rule": "primary <10 (zeros included), complementary S, residual check, LP interval audit, "
                                  "declared links; frozen public cells", "min_cell": MIN_CELL}
    doc["disclosure"] = {"count_fields": list(g.count_fields), "identities": [[t, p] for t, p in g.identities],
                         "hidden_cells": list(g.hidden), "ratios": g.ratios, "quantiles": g.quantiles, "dependents": []}
    doc["release_gate"] = "pending"
    bad_top = [k for k in doc if k not in SM.NON_COUNT_TOP and not isinstance(doc[k], dict)]
    if bad_top:
        raise GateError("internal: non-dict top-level key outside the non-count list")
    declared = set(g.count_fields)
    if declared != set(SM.detect_count_fields(doc)):
        raise GateError("count-field schema mismatch in the audit document")
    vals = g.values()
    status = {p: ("P" if v < MIN_CELL else "pub") for p, v in vals.items()}
    for h in g.hidden:
        status[h] = "H"
    L = SM.Loaded(doc=doc, md=None, count_fields=list(g.count_fields), hidden=list(g.hidden),
                  identities=[(t, tuple(p)) for t, p in g.identities], ratios=[tuple(r) for r in g.ratios],
                  quantiles=[tuple(q) for q in g.quantiles], dependents=[], status=dict(status), values=dict(vals),
                  primary_new=[])
    # v1.3 gate_document compares the output disclosure to L.out_disclosure (the input registry minus the dropped
    # pointers; apply_drops builds it). Without a drop list it is the input disclosure.
    params = SM.Params()
    frozen = set(g.frozen)
    for c in frozen:
        if L.status.get(c) != "pub":
            raise GateError(f"a frozen (already public) cell would be below 10: {c}")
    xinfo = None
    if xr is None:
        L.out_disclosure = L.doc.get("disclosure")
        links = SM.read_links(links_path) if links_path is not None else []
        link_rels = SM.resolve_links(L, links)
        system_rels = [(t, p, "identity") for t, p in L.identities] + list(link_rels)
        extra, evals, tolerated = {}, {}, {}
    else:
        # dropped subtrees: every cell under them becomes D (an unknown in every equation, never published); the
        # linked Stage 3 cells inside them stay in the joint equations on purpose (link_cells = empty set)
        SM.apply_drops(L, xr.drops, set(), frozen)
        if not xr.drops:
            L.out_disclosure = L.doc.get("disclosure")
        s2_rels, extra, evals, link_rels = joint_parts(SM, L, xr)
        system_rels = [(t, p, "identity") for t, p in L.identities] + link_rels + s2_rels
    values = dict(L.values)
    values.update(evals)
    kinds = dict(L.status)
    kinds.update(extra)
    if xr is not None:
        tolerated, via_s2 = compromised_baseline(SM, system_rels, values, kinds, params, extra, L, frozen)
        frozen = frozen | via_s2
        s2_bad = sorted(c[len(XR):] for c in tolerated if c.startswith(XR) and ":hidden:" not in c)
        s2_hid = sum(1 for c in tolerated if c.startswith(XR + "hidden:"))
        s3_bad = sorted(c for c in tolerated if not c.startswith(XR) and L.status.get(c) in ("P", "S"))
        s3_wd = sum(1 for c in tolerated if not c.startswith(XR) and L.status.get(c) == "D")
        xinfo = {"stage2_release": STAGE2_RELEASE_LABEL, "stage2_release_sha256": xr.sha256,
                 "links": f"{len(link_rels)} declared links (links file v3), {len(system_rels) - len(L.identities) - len(link_rels)} "
                          f"Stage 2 identities in the joint system",
                 "stage2_constants": f"{len(xr.constants)} team-internal Stage 2 constant(s) "
                                     f"({sum(1 for _p, v, _r in xr.constants if v == 'derived')} derived) as known "
                                     f"values; values and file hash not recorded",
                 "dropped_subtrees": [t for t, _ in xr.drops],
                 "compromised_stage2_cells": s2_bad,
                 "compromised_stage2_hidden_cells": f"{s2_hid} hidden Stage 2 difference(s)",
                 "compromised_stage3_cells": s3_bad,
                 "compromised_stage3_withheld_cells": f"{s3_wd} withheld (dropped) Stage 3 cell(s)",
                 "stage3_cells_public_via_stage2": f"{len(via_s2)} Stage 3 cell(s) already pinned exactly by Stage 2 "
                                                   f"v1.0.3 + constants + public totals (published, never suppressed); "
                                                   f"{sum(1 for c in tolerated if not c.startswith(XR) and L.status.get(c) == 'pub')} "
                                                   f"more pinned to within the protection width",
                 "compromised_rule": "cells already failing the protection rule on public information alone (Stage 2 "
                                     "v1.0.3 as released + team-internal constants + frozen public totals, every other "
                                     "Stage 3 cell unknown) are flagged for Stage 2 v1.0.4; the joint gate fails if "
                                     "Stage 3 narrows any of them further than that baseline interval"}
        log(f"cross-release baseline: {len(s2_bad)} Stage 2 cell(s), {s2_hid} hidden Stage 2 difference(s) and "
            f"{len(s3_bad) + s3_wd} Stage 3 cell(s) ({s3_wd} withheld) already fail on Stage 2 v1.0.3 + constants + public "
            f"totals (flagged for Stage 2 v1.0.4; Stage 3 must not narrow them further)")
    system = joint_system(SM, system_rels, values, kinds, params, extra, tolerated)
    st = dict(L.status)
    st.update(extra)
    sup = SM.Suppressor(system, st, frozen)
    sup.max_new = int(max_new)                 # v1.0.4: module v1.5 convergence caps (NotConverging -> exit 3)
    sup.max_component = int(max_component)
    res = sup.run()
    out = SM.apply_statuses(L, st)
    if xinfo is not None:
        out.setdefault("notes", {})
        out["notes"][CROSS_NOTE_KEY] = (
            f"audited jointly with the published {STAGE2_RELEASE_LABEL}, SHA-256 {xr.sha256}: its published cells "
            f"and {len(xr.constants)} team-internal constant(s) as known values, its '<10', 'S' and hidden cells as "
            f"unknowns, {len(link_rels)} declared links; no Stage 3 "
            f"or Stage 2 suppressed cell narrows below the protection rule. Cells already exposed by Stage 2 v1.0.3 "
            f"alone are flagged for Stage 2 v1.0.4 and are not narrowed further: "
            + ("; ".join(xinfo["compromised_stage2_cells"]) if xinfo["compromised_stage2_cells"] else "none")
            + (f"; plus {s2_hid} hidden Stage 2 difference(s)" if s2_hid else ""))
    prov = {"script": SCRIPT_VERSION, "suppression_module": SM.SCRIPT_VERSION,
            "new_complementary_cells": len([c for c in sup.rule_of if not c.startswith(XR)]), "rounds": sup.rounds}
    if xinfo is not None:
        prov["cross_release"] = xinfo
    out, _derived = SM.finish_document(out, st, L, res, prov, "generic")
    out = relabel_release(out, SM)
    out["release_gate"] = "ok"
    md = render_markdown(out)
    jp, mp = pending / "stage3_audit.json", pending / "stage3_audit.md"
    write_text(jp, json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    write_text(mp, md)
    re_doc = json.loads(jp.read_text(encoding="utf-8"))
    re_md = mp.read_text(encoding="utf-8")
    system2 = joint_system(SM, system_rels, values, kinds, params, extra, tolerated)
    fails, _res2 = SM.gate_document(re_doc, None, L, system2, frozen, "generic", params, True)
    if render_markdown(re_doc) != re_md:
        fails.append("markdown on disk is not the rendering of the audit JSON")
    if fails:
        raise GateError(f"release gate: {len(fails)} failure(s); first: {fails[0]}")
    log(f"release gate OK: {len(g.count_fields)} count cells ({len(L.dropped)} withheld by the drop list), "
        f"{len(g.hidden)} hidden differences, {len(system_rels)} relations (joint with Stage 2: "
        f"{'yes' if xr is not None else 'no'}), {len(frozen)} frozen cells (incl. public via Stage 2)")
    return re_doc, re_md, prov


def relabel_release(doc: dict, SM) -> dict:
    """v1.0.3/v1.0.4: the pinned suppression module (v1.5) names its provenance block and drop-list note after its
    own release (release_v1_0_3, dropped_subtrees_v1_0_3); the Stage 3 v1.0.4 audit carries them as release_v1_0_4
    and dropped_subtrees_v1_0_4 (same content, the note's release label corrected). Key order kept."""
    old_key = getattr(SM, "PROVENANCE_KEY", "release_v1_0_3")
    old_rel = old_key.replace("release_", "").replace("_", ".")            # e.g. v1.0.3
    old_note = "dropped_subtrees_" + old_key.replace("release_", "")
    out: dict = {}
    for k, v in doc.items():
        out[RELEASE_KEY if k == old_key else k] = v
    notes = out.get("notes")
    if isinstance(notes, dict) and old_note in notes:
        out["notes"] = {(DROP_NOTE_KEY if k == old_note else k):
                        (v.replace(f"withheld from the {old_rel} release", "withheld from the Stage 3 v1.0.4 release")
                         if k == old_note else v) for k, v in notes.items()}
    return out


# ----------------------------------------------------------------------------------------------------------------------
# Writing (atomic, pending folder, meta.json per file, manifest)
# ----------------------------------------------------------------------------------------------------------------------
def _fsync_path(path: Path) -> None:
    with open(path, "rb+") as f:
        os.fsync(f.fileno())


def write_text(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")

    def w():
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())

    fs_op("write", w, tmp)
    fs_op("replace", lambda: os.replace(tmp, path), tmp, path)


def write_csv(path: Path, df: pd.DataFrame) -> None:
    tmp = path.with_name(path.name + ".tmp")
    fs_op("write", lambda: df.to_csv(tmp, index=False, encoding="utf-8", lineterminator="\n", float_format="%.6g"), tmp)
    fs_op("fsync", lambda: _fsync_path(tmp), tmp)
    fs_op("replace", lambda: os.replace(tmp, path), tmp, path)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def median_fallback_note(fb: list) -> dict | str:
    """Decision 3: channel x landmark pairs whose landmark train median was empty, by name with the fallback used
    (no counts)."""
    if not fb:
        return "none (every value channel has a train median of its own at every landmark)"
    out: dict = {}
    for ch, L, src in fb:
        out.setdefault(ch, {})[f"L{L:03d}"] = src
    return out


def unparsed_text_report(vals: Values, co: Cohort, cmap: dict, doc: dict) -> str:
    """Decision 4 (Nov 17, 2026): TEAM-INTERNAL list of distinct unparsed qualitative strings per channel (token-table
    normal form), strings only. Counted on train rows in [t0, W); a string seen fewer than 10 times is not listed; a
    string with digits other than a simple numeric form (up to 4 integer digits, a decimal, %, +, a 1:n titre) is not
    listed (could be an identifier); longer strings are truncated. A channel is listed only if its train unparsed
    count is published in the gated audit, so a listed string adds nothing a reader could not see there. Nothing
    says how many strings were left out."""
    up = vals.unparsed if vals.unparsed is not None else pd.DataFrame(columns=["enc_i", "t", "channel", "text"])
    e = up["enc_i"].to_numpy(dtype=np.int64)
    t = up["t"].to_numpy(dtype=float)
    keep = (t >= 0) & (t < co.W[e]) & (co.split[e] == "train") if len(up) else np.zeros(0, dtype=bool)
    up = up[keep]
    lines = ["# Stage 3 unparsed qualitative text (TEAM-INTERNAL)", "",
             "**Team-internal: do not share outside the project team. Free text can carry names or numbers; Vinai "
             "reads this file before it is pasted anywhere.**", "",
             f"Distinct strings in the token-table normal form (trimmed, lower case, spaces removed), per qualitative "
             f"channel, from train rows in [t0, W). Listed only if seen at least {UNPARSED_MIN} times; strings with "
             f"digits other than simple numeric forms are not listed; strings longer than {UNPARSED_MAX_CHARS} "
             f"characters are truncated (…). No counts. A channel is listed only where its train unparsed count is "
             f"published in stage3_audit.json. Unparsed text stays missing (flag `qual_unparsed` in Layer C, "
             f"`<channel>{UNPARSED_FLAG_SUFFIX}` in Layer B). v1.0.4 extends the token table per channel ('positive+' / "
             f"'positive++' -> positive on stool_occult_blood, urine_leukocyte_esterase, urine_nitrite and "
             f"csf_protein_pandy; '<3.00' -> negative on resp_panel_mycoplasma, flagged "
             f"`resp_panel_mycoplasma{BELOW_CUTOFF_FLAG_SUFFIX}` in Layer B), so those strings are parsed and no "
             f"longer listed here.", ""]
    for ch in qual_channels(cmap):
        lines.append(f"## {ch}")
        try:
            pub = json_ptr_get(doc, f"/plausibility/{ch}/train/n_rows_unparsed")
        except KeyError:
            pub = None
        if not (isinstance(pub, int) and not isinstance(pub, bool)):
            lines += ["", "(not listed: the train unparsed count is not published in the audit)", ""]
            continue
        vc = up.loc[up["channel"] == ch, "text"].value_counts()
        ok = []
        for s_ in sorted(str(x) for x in vc[vc >= UNPARSED_MIN].index):
            if re.search(r"[0-9]", s_) and not UNPARSED_SIMPLE_NUMERIC.fullmatch(s_):
                continue
            ok.append(s_ if len(s_) <= UNPARSED_MAX_CHARS else s_[:UNPARSED_MAX_CHARS] + "…")
        lines.append("")
        lines += [f"- `{x.replace('`', chr(39))}`" for x in sorted(set(ok))] or ["(no string meets the listing rule)"]
        lines.append("")
    return "\n".join(lines)


A_SENSITIVITY_PREFIXES = ("pre24_", "age_days_sens", "age_band_sens")
A_AUDIT_COLUMNS = ("encounter_id", "split", "age_status", "age_neg_bin", "admit_interval_clipped", "admit_interval_missing",
                   "dx_code_status", "recent_band_ref_level", "recent_duration_source", "recent_band_missing_flag")


def file_meta(path: Path, feature_columns: list, input_hashes: dict, extra: dict | None = None) -> dict:
    m = {"file": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size, "plan": PLAN_VERSION,
         "contract": CONTRACT_VERSION, "script": SCRIPT_VERSION, "spec": SPEC_DOC,
         "script_sha256": sha256_file(Path(__file__)), "channel_map": "v2.3", "feature_columns": list(feature_columns),
         "generated_utc": utc_now(), "input_hashes": dict(input_hashes)}
    if extra:
        m.update(extra)
    return m


def build(cfg: Config) -> dict:
    """Run Stage 3. Returns a small summary (counts of files, check results). Raises InputError / GateError."""
    if cfg.raw_dir is None or cfg.split_path is None:
        raise InputError("--raw-dir and --split-file are required (or set PIC_RAW_DIR and SPLIT_FILE_PATH); no local default path is shipped")
    _CLOCK["t0"] = time.monotonic()
    _CLOCK["deadline"] = _CLOCK["t0"] + float(cfg.max_seconds)
    log(f"{SCRIPT_VERSION}: precheck (hash gates before any read)")
    check_out_dir(cfg)
    hashes = precheck_hashes(cfg)
    A = load_pinned_module(Path(cfg.v17_module_path), cfg.pins["v17_module"], "v17")
    SM = load_pinned_module(Path(cfg.suppress_module_path), cfg.pins["suppress_module"], "suppress")
    SM._T0[0] = _CLOCK["t0"]
    SM._DEADLINE[0] = _CLOCK["deadline"]
    hashes.pop("stage2_constants", None)     # team-internal file: its hash is never written to an output
    _SM_REF["m"] = SM
    try:
        with SM.Heartbeat(float(cfg.heartbeat_seconds)):
            return _build_main(cfg, A, SM, hashes)
    finally:
        _SM_REF["m"] = None


def _build_main(cfg: Config, A, SM, hashes: dict) -> dict:
    xr = load_cross_release(cfg, SM, hashes["stage2_release"])
    co, cinfo = load_cohort(cfg, A)
    maps = {k: read_map(cfg, A, k) for k in MAP_FILES}
    cmap = parse_channel_map(maps["channel_map"], cfg.expect)
    dmap = parse_drug_map(maps["drug_map"])
    smap = parse_specimen_map(maps["specimen_map"])
    neuro = neuro_matcher(maps["neuro_map"])
    gref = parse_growth_ref(maps["growth_ref"])
    med1 = read_phi1_medians(cfg)
    log(f"maps parsed: {cmap['map_channels']} channels ({len(cmap['feature_channels'])} feature channels, "
        f"{len(cmap['value_channels'])} value-class; excluded: {', '.join(cmap['excluded_channels']) or 'none'}), "
        f"blood-gas mask {len(cmap['blood_gas_items'])} items")
    if cfg.precheck_only:
        for name, (_sha, rows) in cfg.raw_pins.items():
            if rows is not None and count_rows(raw_path(cfg, name), cfg.chunksize) != rows:
                raise InputError(f"{name}: row count differs from contract Table P")
        log("precheck-only: row counts OK; nothing written")
        return {"mode": "precheck-only", "hashes": hashes}
    out = Path(cfg.out_dir)
    _FS_BASE["dir"] = out.parent
    pending = out.parent / f"{out.name}.pending-{os.getpid()}"
    if pending.exists():
        raise InputError("a pending folder from an earlier run exists; remove it after review")
    fs_op("create folder", lambda: pending.mkdir(parents=True), pending)
    try:
        res = _build_into(cfg, A, SM, co, cinfo, cmap, dmap, smap, neuro, med1, hashes, pending, xr, gref)
        if out.exists():
            fs_op("remove empty output folder", out.rmdir, out)
        fs_op("rename pending folder", lambda: os.replace(pending, out), pending, out)
        log(f"release written: {res['n_files']} files")
        return res
    except BaseException as e:
        release_frames(e)   # close the memory maps held by the failed build's frames before deleting
        try:
            if pending.exists():
                remove_tree(pending, "remove pending folder")
        except Exception:  # noqa: BLE001  (report below; the original error is the one raised)
            pass
        if pending.exists():
            e.stage3_leftover = rel_path(pending)
        raise


def _build_into(cfg, A, SM, co, cinfo, cmap, dmap, smap, neuro, med1, hashes, pending: Path,
                xr: CrossRelease | None = None, gref: dict | None = None) -> dict:
    icu = read_icustays(cfg, A, co)
    pt = read_patients(cfg, A, co)
    ad = read_admissions(cfg, A, co)
    dicd = read_d_icd(cfg, A)
    neuro_enc = np.array([neuro(norm_icd(x)[0]) for x in ad["dx"].to_numpy()], dtype=bool)
    log("static tables read (ICUSTAYS, PATIENTS, ADMISSIONS, D_ICD_DIAGNOSES)")
    docs, rx_info = read_prescriptions(cfg, A, co, dmap, neuro_enc)
    micro, mi_info = read_micro(cfg, A, co, smap)
    chart, ch_info = read_values(cfg, A, co, "CHARTEVENTS", cmap, icu["stays"])
    lab, lb_info = read_values(cfg, A, co, "LABEVENTS", cmap, icu["stays"])
    ab = age_band(pt["age_days"].to_numpy())
    vals = build_values(co, chart, lab, cmap, icu["first_careunit"], ab, pt["sex"].to_numpy(),
                        pt["age_days"].to_numpy(dtype=float), gref)
    del chart, lab
    gc.collect()
    log("channel values built (gates, units, merges, MAP, urine output, fluid intake, KDIGO); weight normalizer filled")
    surg, s_info = surgery_frame(cfg, A, co, med1)
    exp486 = cfg.expect.layer_c_surgeries
    n_c = int((surg["timing"] == "during_icu_layer_c").sum())
    if exp486 is not None and n_c != exp486:
        raise GateError("Layer C surgery count does not reproduce the pin (acceptance check 12; count withheld)")
    log("surgery linked (v1.7 module; ROW_ID; anaesthesia end T); Layer C surgery count reproduces")
    A_df, a_info = build_layer_a(co, pt, ad, icu, dicd, neuro, surg, vals, cmap, docs)
    B_df = build_layer_b(co, vals, cmap, docs, micro)
    log(f"Layer A base ({A_df.shape[1]} columns) and Layer B base ({B_df.shape[1]} columns) built")
    fails: list[str] = []
    cres = build_layer_c(co, vals, cmap, docs, micro, surg, pending, fails)
    ck = cres["checks"]
    log("Layer C raw grid and 14 filled landmark copies written to the pending folder")
    # acceptance checks (spec §10; plan §6.8)
    ck["med_future_enddate_differences"] = future_enddate_differences(co, docs)
    must_zero = ("missing_mask_mismatch", "locf_cap_violations", "padding_violations", "as_of_differences",
                 "lab_release_violations", "filled_padding_or_gap", "started_not_active", "med_as_of_differences",
                 "layer_c_surgery_T_le_6h", "med_future_enddate_differences")
    results = {}
    for k in must_zero:
        results[k] = "pass" if ck[k] == 0 else "fail"
    feats = cmap["feature_channels"]
    results["mask_valid_not_a_feature"] = "pass" if "mask_valid" not in feats else "fail"
    results["weight_not_a_feature"] = "pass" if "weight" not in feats else "fail"
    wa = vals.weight_arms
    ref_ok = all(np.allclose(wa.weight[wa.source == GROWTH_REF_SOURCE[r]],
                             growth_ref_median(gref, r, pt["age_days"].to_numpy(dtype=float), pt["sex"].to_numpy())[
                                 wa.source == GROWTH_REF_SOURCE[r]]) for r in GROWTH_REFS)
    ref_ok &= not ((wa.source != "charted") & np.isin(wa.source, TRAIN_TIERS, invert=True)
                   & (ab == "age_invalid")).any()
    results["weight_reference_rule"] = "pass" if ref_ok else "fail"
    results["no_banned_channel"] = "pass" if not set(BANNED_CHANNELS) & set(feats) else "fail"
    results["eligible_prefix_all_valid"] = "pass" if all((co.nvalid[co.out_h > L] >= L // 6).all() for L in LANDMARKS) else "fail"
    results["blood_gas_mask_15"] = "pass" if len(cmap["blood_gas_items"]) == 15 else "fail"
    results["feature_channels_92"] = "pass" if (len(feats) == cfg.expect.feature_channels
                                                and not set(EXCLUDED_CHANNELS) & set(feats)) else "fail"
    results["layer_a_rows_canonical_order"] = "pass" if list(A_df["encounter_id"]) == list(co.enc["encounter_id"]) else "fail"
    results["layer_b_rows_canonical_order"] = "pass" if list(B_df["encounter_id"]) == list(co.enc["encounter_id"]) else "fail"
    results["train_median_available_all_channels"] = "pass" if ck["train_median_unavailable"] == 0 else "noted"
    # decisions 2, 5, 8, 9 (Nov 17, 2026)
    rst = vals.rowstat
    zbad = rst[(rst["status"] == "zero_missing") & ~rst["channel"].isin(ZERO_AS_MISSING)]
    zleft = sum(int((vals.obs[c]["v"].to_numpy() == 0).sum()) for c in ZERO_AS_MISSING if c in vals.obs)
    results["respiratory_rate_zero_is_missing"] = "pass" if (len(zbad) == 0 and zleft == 0) else "fail"
    dd_ok = True
    for c in ("pco2", "ph") + DEDUPE_BY_DRAW_EXTRA:
        ob = vals.obs.get(c)
        if ob is not None and len(ob) and ob.duplicated(["enc_i", "tc"]).any():
            dd_ok = False
    results["blood_gas_one_value_per_draw"] = "pass" if dd_ok else "fail"
    icdo_ok = bool(((A_df["dx_icd10_3char"] == ICDO_LEVEL) == (A_df["dx_domain"].eq(ONCOLOGY_DOMAIN)
                                                                & A_df["dx_icdo_morphology"].eq(1))).all()
                   & (A_df.loc[A_df["dx_code_status"] == "missing", "dx_icd10_3char"] == "missing").all())
    results["icdo_morphology_oncology_missing_separate"] = "pass" if icdo_ok else "fail"
    bad_cols = [c for c in list(A_df.columns) + list(B_df.columns)
                if re.search(r"edreg|(^|_)ed_component|visit", c, flags=re.I)]
    results["no_edregtime_or_visit_id_feature"] = "pass" if not bad_cols else "fail"
    bad = [k for k, v in results.items() if v == "fail"]
    if bad:
        raise GateError(f"acceptance checks failed: {bad}")
    log(f"acceptance checks: {sum(v == 'pass' for v in results.values())} pass, 0 fail")
    # tabular outputs
    a_feats = [c for c in A_df.columns if c not in A_AUDIT_COLUMNS and not c.startswith(A_SENSITIVITY_PREFIXES)
               and not c.endswith("_t0_lb24")]
    b_flags = [c for c in B_df.columns if c.endswith((UNPARSED_FLAG_SUFFIX, BELOW_CUTOFF_FLAG_SUFFIX))]
    b_feats = [c for c in B_df.columns if c != "encounter_id" and c not in b_flags]
    b_prov = [c for c in b_feats if c.startswith(tuple(f"{x}_" for x in PROVISIONAL_CHANNELS))]
    write_csv(pending / "stage3_layer_a_base.csv", A_df)
    write_csv(pending / "stage3_layer_b_base.csv", B_df)
    wa = vals.weight_arms
    W = pd.DataFrame({"encounter_id": co.enc["encounter_id"].to_numpy(), "weight_kg": vals.weight,
                      "weight_source": vals.weight_source,
                      "weight_kg_train_tier": wa.weight_train_tier, "weight_source_train_tier": wa.source_train_tier,
                      "weight_kg_who": wa.weight_who, "weight_source_who": wa.source_who})
    write_csv(pending / "stage3_weight_normalizer.csv", W)
    sl = surg[["encounter_id", "ROW_ID", "timing", "T_h", "standardized_name", "body_system", "cgroup", "duration_h",
               "duration_source", "fallback_reason", "band_ref_level", "phi1_band", "band_missing_flag", "pre_admit"]]
    write_csv(pending / "stage3_surgery_long_INTERNAL.csv", sl)
    cr = cres["courses"]
    write_csv(pending / "stage3_medication_courses_INTERNAL.csv",
              cr.assign(encounter_id=co.enc["encounter_id"].to_numpy()[cr["enc_i"].to_numpy(dtype=np.int64)])[
                  ["encounter_id", "group", "key", "start", "last", "endE", "cutoff", "stop"]])
    write_csv(pending / "stage3_microbiology_collections_INTERNAL.csv",
              micro.assign(encounter_id=co.enc["encounter_id"].to_numpy()[micro["enc_i"].to_numpy(dtype=np.int64)])[
                  ["encounter_id", "th", "base", "group"]])
    # aggregate audit + release gate
    raw_info = {"PATIENTS": {"rows": cfg.raw_pins["PATIENTS"][1] or 0}, "ADMISSIONS": {"rows": cfg.raw_pins["ADMISSIONS"][1] or 0},
                "ICUSTAYS": {"rows": cfg.raw_pins["ICUSTAYS"][1] or 0}, "SURGERY_INFO": {"rows": s_info["visit"]["n_rows"]},
                "CHARTEVENTS": ch_info, "LABEVENTS": lb_info, "PRESCRIPTIONS": rx_info, "MICROBIOLOGYEVENTS": mi_info}
    g = build_audit(co, pt, A_df, ad, a_info, surg, s_info, vals, cmap, raw_info, ck)
    base = {"spec": AUDIT_SPEC, "input_hashes": dict(hashes),
            "maps": {"channel_map": "feature_channel_map_v2_3.csv (v2.3)", "feature_channels": len(feats),
                     "value_channels": len(cmap["value_channels"]), "map_channels": cmap["map_channels"],
                     "excluded_channels": dict(cmap["excluded_channels"])},
            "cohort_checks": results,
            "config": {"bins": N_BINS, "bin_hours": BIN_H, "lab_release_hours": LAB_RELEASE_H,
                       "pre_icu_hours": PRE_ICU_H, "landmarks": list(LANDMARKS), "core_landmarks": list(CORE_LANDMARKS)},
            "notes": {"gates": "plausibility gates marked † in plan Table 5.1 are provisional; out-of-gate counts are "
                               "for the gate lock (V13, V18)",
                      "weight_rule": WEIGHT_RULE_TEXT,
                      "scope": "aggregate only; counts below 10 shown as '<10'; complementary cells as 'S'",
                      "provisional_features": {"channels": list(PROVISIONAL_CHANNELS), "reason": PROVISIONAL_NOTE},
                      "excluded_channels": dict(cmap["excluded_channels"]),      # v1.0.3: names and reasons only
                      "respiratory_rate_zero": "a respiratory rate of 0 is read as not recorded (missing), not as "
                                               "out of gate (Vinai, Nov 17, 2026); counted as n_rows_zero_as_missing",
                      "train_median_fallback": median_fallback_note(cres["median_fallbacks"]),
                      "unparsed_text": f"distinct unparsed qualitative strings are listed for the team in {UNPARSED_FILE} "
                                       f"(team-internal; strings only, no counts)",
                      "qualitative_tokens": "v1.0.4: 'positive+' / 'positive++' are read as positive on stool_occult_blood, "
                                            "urine_leukocyte_esterase, urine_nitrite and csf_protein_pandy; '<3.00' is read "
                                            "as negative on resp_panel_mycoplasma (below the assay cutoff; flagged in "
                                            "Layer B); counted in n_rows_positive / n_rows_negative"}}
    if xr is None:
        raise GateError("internal: the joint (cross-release) audit needs the Stage 2 v1.0.3 release")
    doc, md, prov = release_gate(SM, g, Path(cfg.links_path) if cfg.links_path else None, base, pending, xr,
                                 cfg.max_new_suppressions, cfg.max_component_cells)
    write_text(pending / UNPARSED_FILE, unparsed_text_report(vals, co, cmap, doc))
    # meta files and manifest
    feat_c = {"channels": feats, "value_channels": cmap["value_channels"]}
    metas = {
        "stage3_layer_a_base.csv": file_meta(pending / "stage3_layer_a_base.csv", a_feats, hashes,
                                             {"sensitivity_columns": [c for c in A_df.columns if c.startswith(A_SENSITIVITY_PREFIXES) or c.endswith("_t0_lb24")]}),
        "stage3_layer_b_base.csv": file_meta(pending / "stage3_layer_b_base.csv", b_feats, hashes,
                                             {"provisional_columns": b_prov, "provisional_reason": PROVISIONAL_NOTE,
                                              "flag_columns": b_flags,
                                              "flag_note": "<channel>_unparsed_text: unparsed qualitative text in the "
                                                           "Layer B window (value stays missing); "
                                                           "<channel>_below_cutoff (v1.0.4): a negative read from a "
                                                           "below-cutoff result ('<3.00') in the Layer B window; flags "
                                                           "are not feature columns"}),
        "stage3_weight_normalizer.csv": file_meta(pending / "stage3_weight_normalizer.csv", [], hashes,
                                                  {"use": "normalizer only (V12); weight_source for audits and B38",
                                                   "rule": WEIGHT_RULE_TEXT,
                                                   "growth_reference": {"file": MAP_FILES["growth_ref"],
                                                                        "sha256": hashes["growth_ref"],
                                                                        "references": list(GROWTH_REFS)},
                                                   "sensitivity_columns": ["weight_kg_train_tier",
                                                                           "weight_source_train_tier",
                                                                           "weight_kg_who", "weight_source_who"]}),
        "stage3_layer_c_raw.npz": file_meta(pending / "stage3_layer_c_raw.npz", feats, hashes,
                                            {"masks": ["mask_valid", "missing", "panel_missing",
                                                       "hemoglobin_bloodgas_missing"],
                                             "flags": ["qual_unparsed (with qual_channels; not a feature)"],
                                             "provisional_channels": list(PROVISIONAL_CHANNELS),
                                             "provisional_reason": PROVISIONAL_NOTE, **feat_c}),
    }
    for L in LANDMARKS:
        nm = f"stage3_layer_c_filled_L{L}.npz"
        metas[nm] = file_meta(pending / nm, cmap["value_channels"], hashes,
                              {"landmark_hours": L, "prefix_bins": L // 6,
                               "fill": "bin median; else last observed value within the LOCF cap; else the train "
                                       "median of train eligible prefixes at L (imputed flag); if that is empty the "
                                       "nearest earlier landmark's own train median, else the pooled all-landmark "
                                       "train median (train_median_source)",
                               "train_median_fallbacks": {c: f"{sname}" for c, LL, sname in cres["median_fallbacks"]
                                                          if LL == L},
                               "provisional_channels": [c for c in PROVISIONAL_CHANNELS if c in cmap["value_channels"]],
                               "provisional_reason": PROVISIONAL_NOTE})
    for nm in ("stage3_surgery_long_INTERNAL.csv", "stage3_medication_courses_INTERNAL.csv",
               "stage3_microbiology_collections_INTERNAL.csv"):
        metas[nm] = file_meta(pending / nm, [], hashes, {"use": "internal long-form; never shared"})
    metas[UNPARSED_FILE] = file_meta(pending / UNPARSED_FILE, [], hashes,
                                     {"use": "TEAM-INTERNAL review list (decision 4): strings only, no counts; may "
                                             "hold free text; Vinai reads it before it is pasted anywhere"})
    for nm, m in metas.items():
        write_text(pending / f"{nm}.meta.json", json.dumps(m, indent=1) + "\n")
    files = {}
    for p in sorted(pending.iterdir()):
        if p.is_file():
            files[p.name] = {"sha256": sha256_file(p), "bytes": p.stat().st_size}
    manifest = {"script": SCRIPT_VERSION, "script_sha256": sha256_file(Path(__file__)), "spec": SPEC_DOC,
                "plan": PLAN_VERSION, "contract": CONTRACT_VERSION, "generated_utc": utc_now(),
                "input_hashes": dict(hashes),
                "precheck_recorded": {"D_ICD_DIAGNOSES": hashes["raw:D_ICD_DIAGNOSES"],
                                      "SURGERY_VITAL_SIGNS": hashes["raw:SURGERY_VITAL_SIGNS"]},
                "files": files, "release_gate": "ok", "acceptance_checks": results,
                "cross_release": {"stage2_release": STAGE2_RELEASE_LABEL, "stage2_release_sha256": xr.sha256,
                                  "links_sha256": hashes.get("links"), "drops_sha256": hashes.get("drops"),
                                  "dropped_subtrees": [t for t, _ in xr.drops],
                                  "stage2_constants": prov["cross_release"]["stage2_constants"],
                                  FLAG_KEY: prov["cross_release"]["compromised_stage2_cells"],
                                  "flagged_hidden_stage2_differences": prov["cross_release"]["compromised_stage2_hidden_cells"]}}
    write_text(pending / "stage3_manifest.json", json.dumps(manifest, indent=1) + "\n")
    return {"mode": "build", "n_files": len(files) + 1, "checks": results, "audit_md": md,
            "pending_name": pending.name}


# ----------------------------------------------------------------------------------------------------------------------
# Command line (real run: pins and expectations are fixed; only paths, the two precheck hashes and the time limit vary)
# ----------------------------------------------------------------------------------------------------------------------
def build_config(argv: list[str] | None = None) -> Config:
    here = Path(__file__).resolve()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--raw-dir", default=os.environ.get("PIC_RAW_DIR"))
    ap.add_argument("--split-file", default=os.environ.get("SPLIT_FILE_PATH"))
    ap.add_argument("--phi1-medians", default=None)
    ap.add_argument("--maps-dir", default=str(here.parent / "maps"))
    ap.add_argument("--links", default=str(here.parent / "stage3_links_v3.csv"))
    ap.add_argument("--stage2-constants", default=str(here.parent / "stage3_stage2_constants_TEAM_INTERNAL.csv"),
                    help="TEAM-INTERNAL Stage 2 constants file (hash-pinned; values never printed)")
    ap.add_argument("--drops", default=str(here.parent / "stage3_drops_v1.csv"))
    ap.add_argument("--stage2-release", required=True,
                    help="Stage 2 v1.0.3 release JSON (feature_sparsity_frequency_audit_v1.0.3.json)")
    ap.add_argument("--expect-stage2-sha256", required=True,
                    help="full SHA-256 of the Stage 2 v1.0.3 release JSON (must equal the script's pin)")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--expect-d-icd-sha256", required=True)
    ap.add_argument("--expect-svs-sha256", required=True)
    ap.add_argument("--precheck-only", action="store_true")
    ap.add_argument("--max-seconds", type=float, default=7200.0)
    ap.add_argument("--max-new-suppressions", type=int, default=2500,
                    help="suppression module v1.5 convergence cap (stop, exit 3, if exceeded)")
    ap.add_argument("--max-component-cells", type=int, default=20000,
                    help="suppression module v1.5 convergence cap on one equation component")
    ap.add_argument("--heartbeat-seconds", type=float, default=30.0,
                    help="a heartbeat line whenever this long passes without a progress line")
    a = ap.parse_args(argv)
    root = Path(a.project_root)
    interim = root / "data" / "interim"
    return Config(raw_dir=Path(a.raw_dir) if a.raw_dir else None, canonical_path=interim / "canonical_icu_encounters.csv",
                  harmonization_path=interim / "encounter_harmonization_fields.csv", split_path=Path(a.split_file) if a.split_file else None,
                  complexity_map_path=Path(a.maps_dir) / "surgery_complexity_map_v1.csv",
                  phi1_medians_path=Path(a.phi1_medians) if a.phi1_medians else Path(a.maps_dir) / "phi1_band_support_frozen_medians_MAP_v1.csv",
                  maps_dir=Path(a.maps_dir),
                  v17_module_path=here.parent / "audit_surgery_complexity_duration.py",
                  suppress_module_path=here.parent / "stage2_v1_0_3_suppress.py",
                  links_path=Path(a.links), out_dir=Path(a.out_dir) if a.out_dir else root / "data" / "processed" / "ladder" / "stage3",
                  expect_d_icd_sha256=a.expect_d_icd_sha256.strip().lower(),
                  expect_svs_sha256=a.expect_svs_sha256.strip().lower(), max_seconds=a.max_seconds,
                  precheck_only=a.precheck_only, links_sha256=LINKS_SHA256,
                  stage2_release_path=Path(a.stage2_release), expect_stage2_sha256=a.expect_stage2_sha256.strip().lower(),
                  drops_path=Path(a.drops), drops_sha256=DROPS_SHA256,
                  stage2_constants_path=Path(a.stage2_constants), stage2_constants_sha256=CONSTANTS_SHA256,
                  max_new_suppressions=a.max_new_suppressions, max_component_cells=a.max_component_cells,
                  heartbeat_seconds=a.heartbeat_seconds)


LINKS_SHA256: str | None = "cb4c68746156950dc61476e44868c5c30720c17255a69562b2e7f3ec18fd79c0"   # stage3_links_v3.csv (spec §16)
DROPS_SHA256: str | None = "4c491c8d51c68fe31e41cf321bdb4f26d53ed4a3eaed36cb447cddae1ea829f6"   # stage3_drops_v1.csv (spec §16)
# v1.0.4: stage3_stage2_constants_TEAM_INTERNAL.csv (team-internal; never recorded in an output)
CONSTANTS_SHA256: str | None = "1cba9684e5f3eee02927b3511bced224767a392b0a9145a41e5246a0b5cd4563"


def main(argv: list[str] | None = None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:  # noqa: BLE001
        pass
    def left(e: BaseException, clean: str) -> str:
        lo = getattr(e, "stage3_leftover", None)
        return clean if not lo else (f"pending folder left behind: {lo} (could not be removed; delete it after review, "
                                     f"never share or sync it)")

    try:
        cfg = build_config(argv)
        res = build(cfg)
    except InputError as e:
        print(f"STOPPED (input): {e}", file=sys.stderr)
        print(f"result: STOPPED (input); {left(e, 'nothing written')}")
        return 2
    except TimeLimit as e:
        print(f"BUILD/RELEASE GATE FAILED: time limit (--max-seconds) reached; {left(e, 'nothing left behind')}",
              file=sys.stderr)
        print(f"result: FAILED (time limit); {left(e, 'nothing left behind')}")
        return 3
    except GateError as e:
        print(f"BUILD/RELEASE GATE FAILED: {e}", file=sys.stderr)
        print(f"result: FAILED; {left(e, 'nothing left behind')}")
        return 3
    except FileOpError as e:
        loc = where(e)
        print(f"BUILD FAILED (file operation: {e.describe()}{'; at ' + loc if loc else ''}); "
              f"{left(e, 'nothing left behind')}", file=sys.stderr)
        print(f"result: FAILED; {left(e, 'nothing left behind')}")
        return 3
    except Exception as e:  # noqa: BLE001  (any other failure also fails closed; the message text is not printed)
        bits = []
        if isinstance(e, OSError):
            for k in ("errno", "winerror"):
                if getattr(e, k, None) is not None:
                    bits.append(f"{'errno' if k == 'errno' else 'WinError'} {getattr(e, k)}")
            if getattr(e, "filename", None):
                bits.append(f"path {rel_path(e.filename)}")
        loc = where(e)
        if loc:
            bits.append(f"at {loc}")
        detail = f": {'; '.join(bits)}" if bits else ""
        print(f"BUILD FAILED (unexpected {type(e).__name__}{detail}); {left(e, 'nothing left behind')}", file=sys.stderr)
        print(f"result: FAILED; {left(e, 'nothing left behind')}")
        return 3
    if res["mode"] == "precheck-only":
        print("result: PRECHECK OK (hashes, row counts, cohort, maps); nothing written")
        for k in ("raw:D_ICD_DIAGNOSES", "raw:SURGERY_VITAL_SIGNS"):
            print(f"{k.split(':')[1]} sha256 {res['hashes'][k]}")
        return 0
    print(res["audit_md"])
    print("result: STAGE 3 OK; RELEASE GATE OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
