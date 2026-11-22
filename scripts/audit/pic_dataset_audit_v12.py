#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pic_dataset_audit_v12.py
========================
Aggregate-only audit of the PIC database for Foundations Plan v12
(PIC Foundations; spec: revised_plan_v12.md, v12_change_list.md).

Purpose
-------
Resolve the bracketed v12 audit items before stages 06/07/07a/08x/15 are coded:

  (A) Layer A        interval and X thresholds, ED registration, surgery vs INTIME,
                     diagnosis missingness, sex missing, Chinese-ICD mapping, rare
                     3-char categories, chapter definition, split chapters (Q, D, I),
                     care unit, prior ICU admissions, inferred admission source.
  (B) Dictionaries   D_LABITEMS / D_ITEMS keyword search for the v12 audit targets and
                     every Table 5.1 item (candidate ITEMIDs, ambiguity flags).
  (C) Layers B/C     per-item and per-channel coverage 0-24 h and per 6 h bin to 168 h,
                     >=3 values in 24 h, units, plausibility, panel co-recording,
                     HR charting gaps, weight completeness, urine output / intake
                     interval structure.
  (D) Release timing extra time columns in LABEVENTS / MICROBIOLOGYEVENTS and their lag
                     from CHARTTIME (safeguard S9).
  (E) Drugs          vasoactive dose fields, NMB and antiepileptic mapping, ENDDATE
                     missingness.

Every report section carries a "resolves" list naming the v12 bracket / open item.

Privacy and data-use rules (plan §14, contract §12)
---------------------------------------------------
* Outputs are aggregates only: counts, percentages, quantiles.  No SUBJECT_ID,
  HADM_ID, ICUSTAY_ID, encounter_id, or any row-level value is written.
* Counts 1..9 are shown as "<10" (SMALL_CELL); percentages built on a suppressed
  numerator are suppressed too.  Quantiles are reported only when n >= MIN_N_QUANTILE
  and never include the minimum or maximum (single-patient values).
* Dictionary labels (D_LABITEMS / D_ITEMS), drug names, surgery names and organism
  names are vocabulary, not patient data; their counts are still suppressed.
* Note: primary suppression only.  Complementary suppression (a suppressed cell that
  can be recovered from a total) is NOT applied; review the Markdown before sharing.

How to run (Vinai, locally, credentialed environment)
-----------------------------------------------------
    python audit/pic_dataset_audit_v12.py                 # full audit
    python audit/pic_dataset_audit_v12.py --sections A,B  # subset
    python audit/pic_dataset_audit_v12.py --self-test     # synthetic tests only
    pytest audit/pic_dataset_audit_v12.py                 # same tests via pytest

Outputs: <project>/data/interim/pic_dataset_audit_v12.json and .md

This file was written without being run against PIC (no data was available to the
author).  Bracketed items [like this] are names or choices to confirm.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# =====================================================================================
# CONFIG  (edit here; everything below reads from this block)
# =====================================================================================
CONFIG: Dict[str, Any] = {
    # ---------------- paths ----------------
    # Project root; None = auto (walk up from this file / cwd looking for configs/ or data/).
    "PROJECT_ROOT": None,
    # Contract §2: the raw root is chosen in configs/data_source.json.
    "DATA_SOURCE_CONFIG": "configs/data_source.json",
    # [key names inside data_source.json: confirm]; first key that holds an existing dir wins.
    "DATA_SOURCE_KEYS": ["raw_root", "data_root", "pic_root", "root", "path", "raw_data_root"],
    # Fallback candidates, relative to the project root (contract: data/raw or sibling PICDataset/).
    "RAW_ROOT_CANDIDATES": ["../PICDataset", "data/raw", "../PICData", "../Analysis_PIC_Final/data/raw", "PICDataset", "../PIC", "data/PICDataset"],  # first four from Vinai's draft audit
    "RAW_ROOT_ENV_VAR": "PIC_ROOT",        # overrides everything if set
    "RAW_ROOT_OVERRIDE": None,             # or an absolute path string
    "RAW_SEARCH_DEPTH": 2,                 # sub-folder depth searched for ICUSTAYS.csv(.gz)
    # Canonical cohort (contract §1.1).  Column names [confirm against the file header].
    "COHORT_FILE": "data/interim/canonical_icu_encounters.csv",
    "COHORT_COLUMNS": {"encounter_id": "encounter_id", "subject_id": "SUBJECT_ID",
                       "intime": "INTIME", "outtime": "OUTTIME", "hadm_id": "HADM_ID",
                       "hadm_id_list": "hadm_id_list"},
    "SPLIT_FILE": "data/processed/subject_splits.csv",           # SUBJECT_ID + split [column names: confirm]
    "LABELS_FILE": "data/processed/encounter_labels.csv",        # fallback: encounter_id + split
    "SPLIT_COLUMN_CANDIDATES": ["split", "SPLIT", "partition"],
    # If the canonical file is missing, rebuild the cohort in memory with the contract §4.2
    # 24 h running-episode collapse (never written to disk).
    "REBUILD_COHORT_IF_MISSING": True,
    "EXPECTED_N_ENCOUNTERS": None,  # public copy: pin omitted
    "EXPECTED_SPLIT_COUNTS": None,
    "OUT_DIR": "data/interim",
    "OUT_STEM": "pic_dataset_audit_v12",

    # ---------------- reading ----------------
    "CHUNKSIZE": 500_000,          # lower if memory is tight (panel matrix is per chunk)
    "ENCODING": "utf-8",
    # Expected PIC timestamp format; rows that fail are retried with format="mixed".
    "DATETIME_FORMAT": "%Y-%m-%d %H:%M:%S",
    "SAMPLE_CAP": 100_000,       # per-channel value reservoir for quantiles
    "ITEM_SAMPLE_CAP": 20_000,   # per-ITEMID value reservoir
    "RANDOM_SEED": 6310,

    # ---------------- privacy ----------------
    "SMALL_CELL": 10,            # counts 1..SMALL_CELL-1 -> "<10"
    "SUPPRESS_ZERO": False,
    "MIN_N_QUANTILE": 10,
    "QUANTILES": [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99],
    "TOP_N_NAMES": 30,

    # ---------------- windows (plan §5-§6) ----------------
    "LAYER_B_HOURS": 24,
    "WINDOW_HOURS": 168,
    "BIN_HOURS": 6,
    "PRE_HOURS": 24,             # look-back kept for weight / pre-ICU context only
    "MIN_VALUES_VARIABILITY": 3,

    # ---------------- Layer A (plan §4) ----------------
    "X_GRID_HOURS": [6, 12, 24, 48, 72, 168],     # candidate floor thresholds X
    "INTERVAL_BAND_EDGES": [24, 72, 168],          # v12 bands <24, 24-<72, 72-<168, >=168
    "RARE_CATEGORY_MIN_TRAIN_PATIENTS": 20,
    "SURGERY_NEAR_INTIME_MINUTES": [1, 15, 60],    # timestamp-rounding check around INTIME
    # Appendix A DRAFT domain map choices (bracketed in the plan).
    "DOMAIN_CHAPTER_I": "other_medical",           # or "cardiac_congenital_heart"
    "DOMAIN_D00_D49": "other_medical",             # or "oncology_hematology"
    "DOMAIN_EXTERNAL_CAUSE_V_Y": "surgical_trauma",  # or "other_medical"
    # v12 timing categories
    "DAY_START_HOUR": 8, "EVENING_START_HOUR": 15, "NIGHT_START_HOUR": 23,
    # EDREGTIME: which linked admission row supplies it [confirm]: "earliest" or "any"
    "ED_SOURCE_ROW": "earliest",

    # ---------------- Layer B/C channels (plan Table 5.1) ----------------
    # itemids: contract §11 where known (native units, before conversion).  Channels with an
    # empty list use the Section B dictionary candidates (flagged "auto (review)").
    # range: plausibility range in NATIVE units [draft unless note == "contract"].
    # Filled below by _build_channels().
    "CHANNELS": {},
    "USE_DICTIONARY_CANDIDATES_FOR_EMPTY_CHANNELS": True,
    "MAX_AUTO_CANDIDATES_PER_CHANNEL": 8,

    # ---------------- Section E drugs ----------------
    "DRUG_GROUPS": {},           # filled by _build_drug_groups()
    "DOSE_COLUMN_CANDIDATES": ["DOSE_VAL_RX", "DOSE_UNIT_RX", "PROD_STRENGTH", "DRUG_FORM",
                               "FORM_VAL_DISP", "FORM_UNIT_DISP", "ROUTE", "RATE", "RATEUOM",
                               "DOSE", "DOSAGE", "FREQUENCY"],

    # ---------------- table file stems (case-insensitive; .csv or .csv.gz) ----------------
    "TABLES": {"PATIENTS": "PATIENTS", "ADMISSIONS": "ADMISSIONS", "ICUSTAYS": "ICUSTAYS",
               "CHARTEVENTS": "CHARTEVENTS", "LABEVENTS": "LABEVENTS",
               "PRESCRIPTIONS": "PRESCRIPTIONS", "SURGERY_INFO": "SURGERY_INFO",
               "MICROBIOLOGYEVENTS": "MICROBIOLOGYEVENTS", "D_ICD_DIAGNOSES": "D_ICD_DIAGNOSES",
               "D_LABITEMS": "D_LABITEMS", "D_ITEMS": "D_ITEMS",
               "OUTPUTEVENTS": "OUTPUTEVENTS", "INPUTEVENTS": "INPUTEVENTS"},
}


def _ch(system, source, itemids=None, en=(), cn=(), rng=None, note="draft", triggered=False,
        collect_times=False, sum24=False, fluid_hint=None):
    return {"system": system, "source": source, "itemids": list(itemids or []),
            "keywords_en": list(en), "keywords_cn": list(cn), "range": rng, "range_note": note,
            "triggered": triggered, "collect_times": collect_times, "sum24": sum24,
            "fluid_hint": fluid_hint}


def _build_channels() -> Dict[str, Dict[str, Any]]:
    """Table 5.1 channels.  Contract §11 ITEMIDs and gates where they exist; everything
    else is a keyword search target with a DRAFT native-unit range."""
    C = {}
    # ---- Cardiovascular
    C["hr"] = _ch("cardiovascular", "chart", ["1002", "1003", "SV1", "SV2"], ["heart rate", "HR", "pulse"], ["心率", "脉搏"], (30, 300), "contract", collect_times=True)
    C["sbp"] = _ch("cardiovascular", "chart", ["1016", "SV3", "SV10"], ["systolic", "SBP", "NBP s", "ABP s"], ["收缩压"], (20, 300))
    C["dbp"] = _ch("cardiovascular", "chart", ["1015", "SV4", "SV11"], ["diastolic", "DBP"], ["舒张压"], (5, 250))
    C["map"] = _ch("cardiovascular", "chart", ["SV5", "SV12"], ["mean arterial", "MAP", "mean BP"], ["平均压", "平均动脉压"], (20, 200), "contract")
    C["lactate"] = _ch("cardiovascular", "lab", ["5227"], ["lactate", "lactic acid", "LAC"], ["乳酸"], (0.1, 30), "contract")
    C["troponin_t"] = _ch("cardiovascular", "lab", [], ["troponin T", "TnT", "cTnT", "hs-TnT"], ["肌钙蛋白T"], (0, 100), triggered=True)
    C["troponin_i"] = _ch("cardiovascular", "lab", [], ["troponin I", "TnI", "cTnI"], ["肌钙蛋白I"], (0, 500), triggered=True)
    C["bnp"] = _ch("cardiovascular", "lab", [], ["BNP", "NT-proBNP", "natriuretic"], ["脑钠肽", "利钠肽"], (0, 100000), triggered=True)
    # ---- Respiratory
    C["rr"] = _ch("respiratory", "chart", ["1004", "SV13"], ["respiratory rate", "RR", "resp rate"], ["呼吸频率", "呼吸"], (4, 150), "contract")
    C["spo2"] = _ch("respiratory", "chart", ["1006", "SV17"], ["SpO2", "oxygen saturation", "pulse ox"], ["血氧饱和度", "经皮血氧"], (40, 100), "contract")
    C["ph"] = _ch("respiratory", "lab", ["5237", "5238", "5374"], ["pH"], ["酸碱度"], (6.5, 7.8), "contract")
    C["pco2"] = _ch("respiratory", "lab", [], ["pCO2", "PaCO2", "carbon dioxide partial"], ["二氧化碳分压"], (5, 200))
    C["po2"] = _ch("respiratory", "lab", [], ["pO2", "PaO2", "oxygen partial"], ["氧分压"], (5, 700))
    C["base_excess"] = _ch("respiratory", "lab", [], ["base excess", "BE", "BEecf", "BEb", "ABE", "SBE"], ["碱剩余", "碱超"], (-40, 40))
    C["sao2"] = _ch("respiratory", "lab", [], ["SaO2", "sO2", "O2 sat", "oxygen saturation"], ["氧饱和度"], (10, 100))
    C["bicarbonate"] = _ch("respiratory", "lab", ["5224", "5248"], ["bicarbonate", "HCO3", "cHCO3", "SB", "AB"], ["碳酸氢根", "碳酸氢盐"], (2, 60), "contract")
    C["fio2"] = _ch("respiratory", "both", [], ["FiO2", "FIO2", "FO2", "inspired oxygen", "oxygen concentration"], ["吸入氧浓度", "氧浓度"], (21, 100))
    # ---- Renal
    C["creatinine"] = _ch("renal", "lab", ["5032", "5041", "6954"], ["creatinine", "CREA", "Cr"], ["肌酐"], (4.42, 2210), "contract (umol/L)")
    C["bun"] = _ch("renal", "lab", ["5033", "6477"], ["urea", "BUN", "urea nitrogen"], ["尿素", "尿素氮"], (0.357, 71.4), "contract (mmol/L)")
    C["cystatin_c"] = _ch("renal", "lab", [], ["cystatin", "Cys-C", "CysC"], ["胱抑素"], (0.1, 15))
    C["urine_output"] = _ch("renal", "output", [], ["urine", "urine output", "UO"], ["尿量", "尿"], (0, 5000), collect_times=True, sum24=True)
    C["fluid_intake"] = _ch("renal", "input", [], ["intake", "input"], ["入量"], (0, 10000), collect_times=True, sum24=True)
    # ---- Hepatic
    C["alt"] = _ch("hepatic", "lab", [], ["ALT", "alanine aminotransferase", "GPT", "SGPT"], ["丙氨酸氨基转移酶", "谷丙转氨酶"], (1, 20000))
    C["ast"] = _ch("hepatic", "lab", [], ["AST", "aspartate aminotransferase", "GOT", "SGOT"], ["天门冬氨酸氨基转移酶", "天冬氨酸氨基转移酶", "谷草转氨酶"], (1, 30000))
    C["bilirubin_total"] = _ch("hepatic", "lab", [], ["total bilirubin", "TBIL", "T-BIL", "bilirubin total"], ["总胆红素"], (0.5, 1000))
    C["bilirubin_direct"] = _ch("hepatic", "lab", [], ["direct bilirubin", "DBIL", "D-BIL", "conjugated bilirubin"], ["直接胆红素", "结合胆红素"], (0, 700))
    C["albumin"] = _ch("hepatic", "lab", [], ["albumin", "ALB"], ["白蛋白"], (5, 70), fluid_hint="blood")
    C["ggt"] = _ch("hepatic", "lab", [], ["GGT", "gamma-glutamyl", "gamma glutamyl", "γ-GT"], ["谷氨酰转肽酶", "谷氨酰转移酶"], (1, 3000))
    C["ammonia"] = _ch("hepatic", "lab", [], ["ammonia", "NH3", "blood ammonia"], ["血氨"], (1, 1500), triggered=True)
    # ---- GI
    C["amylase"] = _ch("gi", "lab", [], ["amylase", "AMY", "AMS"], ["淀粉酶"], (1, 10000))
    C["lipase"] = _ch("gi", "lab", [], ["lipase", "LPS", "LIP"], ["脂肪酶"], (1, 20000))
    C["stool_occult_blood"] = _ch("gi", "lab", [], ["occult blood", "FOB", "OB", "guaiac"], ["隐血", "潜血"], None)
    # ---- Hematologic
    C["hemoglobin"] = _ch("hematologic", "lab", ["5099", "5257"], ["hemoglobin", "haemoglobin", "HGB", "Hb"], ["血红蛋白"], (10, 300), "contract (g/L)")
    C["platelets"] = _ch("hematologic", "lab", ["5129"], ["platelet", "PLT"], ["血小板"], (1, 2000), "contract")
    C["lymphocytes_abs"] = _ch("hematologic", "lab", [], ["lymphocyte absolute", "LYMPH#", "LY#", "absolute lymphocyte"], ["淋巴细胞绝对值", "淋巴细胞计数"], (0, 100))
    C["ldh"] = _ch("hematologic", "lab", [], ["LDH", "lactate dehydrogenase", "LD"], ["乳酸脱氢酶"], (10, 50000))
    C["inr"] = _ch("hematologic", "lab", [], ["INR", "international normalized ratio"], ["国际标准化比值", "国际标准化比率"], (0.5, 15))
    C["ptt"] = _ch("hematologic", "lab", [], ["APTT", "PTT", "partial thromboplastin"], ["部分凝血活酶时间"], (10, 250))
    C["fibrinogen"] = _ch("hematologic", "lab", [], ["fibrinogen", "FIB", "Fg"], ["纤维蛋白原"], (0.1, 15))
    C["d_dimer"] = _ch("hematologic", "lab", [], ["D-dimer", "D dimer", "DD", "D-D"], ["D-二聚体", "D二聚体"], (0, 200))
    C["ferritin"] = _ch("hematologic", "lab", [], ["ferritin", "FER", "SF"], ["铁蛋白"], (1, 200000), triggered=True)
    C["teg"] = _ch("hematologic", "lab", [], ["TEG", "thromboelastography", "R time"], ["血栓弹力图", "血栓弹力"], None, triggered=True)
    # ---- Infectious
    C["temperature"] = _ch("infectious", "chart", ["1001", "1031", "SV14", "SV15"], ["temperature", "temp", "T"], ["体温"], (28, 43), "contract")
    C["wbc"] = _ch("infectious", "lab", ["5141"], ["white blood cell", "WBC", "leukocyte"], ["白细胞"], (0.1, 300), "contract", fluid_hint="blood")
    C["anc"] = _ch("infectious", "lab", [], ["neutrophil absolute", "NEUT#", "NE#", "absolute neutrophil", "ANC"], ["中性粒细胞绝对值", "中性粒细胞计数"], (0, 200))
    C["crp"] = _ch("infectious", "lab", [], ["CRP", "C-reactive protein", "hs-CRP"], ["C反应蛋白", "C-反应蛋白"], (0, 500))
    C["urine_wbc"] = _ch("infectious", "lab", [], ["urine WBC", "urine white", "leukocyte esterase", "LEU"], ["尿白细胞", "白细胞(尿)"], None, fluid_hint="urine")
    C["urine_nitrite"] = _ch("infectious", "lab", [], ["nitrite", "NIT"], ["亚硝酸盐"], None, fluid_hint="urine")
    C["procalcitonin"] = _ch("infectious", "lab", [], ["procalcitonin", "PCT"], ["降钙素原"], (0, 1000), triggered=True)
    C["il6"] = _ch("infectious", "lab", [], ["IL-6", "interleukin-6", "interleukin 6"], ["白介素-6", "白细胞介素-6", "白细胞介素6"], (0, 100000), triggered=True)
    C["il10"] = _ch("infectious", "lab", [], ["IL-10", "interleukin-10", "interleukin 10"], ["白介素-10", "白细胞介素-10"], (0, 100000), triggered=True)
    C["tnf_alpha"] = _ch("infectious", "lab", [], ["TNF", "TNF-α", "TNF-alpha", "tumor necrosis factor"], ["肿瘤坏死因子"], (0, 100000), triggered=True)
    C["ifn_gamma"] = _ch("infectious", "lab", [], ["IFN-γ", "IFN-gamma", "interferon gamma", "interferon-γ"], ["干扰素", "γ-干扰素"], (0, 100000), triggered=True)
    C["beta_d_glucan"] = _ch("infectious", "lab", [], ["beta-D-glucan", "β-D-glucan", "BDG", "G test", "(1,3)-β-D"], ["葡聚糖", "G试验"], (0, 10000), triggered=True)
    # ---- Neurologic
    C["pain_score"] = _ch("neurologic", "chart", [], ["pain", "pain score", "FLACC", "NRS", "CRIES"], ["疼痛", "疼痛评分"], (0, 10))
    C["csf_wbc"] = _ch("neurologic", "lab", [], ["CSF WBC", "CSF white", "CSF cell count", "cerebrospinal fluid white"], ["脑脊液白细胞", "脑脊液细胞"], None, fluid_hint="csf")
    C["csf_protein"] = _ch("neurologic", "lab", [], ["CSF protein", "cerebrospinal fluid protein"], ["脑脊液蛋白", "脑脊液总蛋白"], None, fluid_hint="csf")
    C["gcs"] = _ch("neurologic", "chart", [], ["GCS", "Glasgow", "coma scale", "eye opening", "verbal response", "motor response"], ["格拉斯哥", "昏迷评分", "昏迷"], (3, 15))
    C["nse"] = _ch("neurologic", "lab", [], ["NSE", "neuron specific enolase", "neuron-specific enolase"], ["神经元特异性烯醇化酶"], (0, 1000), triggered=True)
    # ---- Electrolytes / homeostasis
    C["sodium"] = _ch("electrolytes", "lab", ["5062", "5230"], ["sodium", "Na", "Na+"], ["钠"], (90, 190), "contract", fluid_hint="blood")
    C["potassium"] = _ch("electrolytes", "lab", ["5226"], ["potassium", "K", "K+"], ["钾"], (1, 12), "contract", fluid_hint="blood")
    C["chloride"] = _ch("electrolytes", "lab", [], ["chloride", "Cl", "Cl-"], ["氯"], (50, 160), fluid_hint="blood")
    C["calcium_total"] = _ch("electrolytes", "lab", [], ["calcium", "Ca", "total calcium"], ["钙", "总钙"], (0.5, 5.0), fluid_hint="blood")
    C["magnesium"] = _ch("electrolytes", "lab", [], ["magnesium", "Mg"], ["镁"], (0.1, 5.0), fluid_hint="blood")
    C["phosphate"] = _ch("electrolytes", "lab", [], ["phosphate", "phosphorus", "PHOS"], ["磷", "无机磷"], (0.1, 6.0), fluid_hint="blood")
    C["anion_gap"] = _ch("electrolytes", "lab", [], ["anion gap", "AG"], ["阴离子间隙"], (-10, 60))
    # ---- Endocrine / metabolic
    C["glucose"] = _ch("endocrine_metabolic", "both", ["1011", "5047", "5223"], ["glucose", "GLU", "blood sugar"], ["葡萄糖", "血糖"], (0.5, 60), "contract (mmol/L)")
    C["cortisol"] = _ch("endocrine_metabolic", "lab", [], ["cortisol", "COR"], ["皮质醇"], None)
    C["tsh"] = _ch("endocrine_metabolic", "lab", [], ["TSH", "thyroid stimulating", "thyrotropin"], ["促甲状腺激素", "促甲状腺"], None)
    C["free_t4"] = _ch("endocrine_metabolic", "lab", [], ["free T4", "FT4", "free thyroxine"], ["游离甲状腺素", "游离T4"], None)
    C["insulin"] = _ch("endocrine_metabolic", "lab", [], ["insulin", "INS"], ["胰岛素"], None)
    C["ketones"] = _ch("endocrine_metabolic", "lab", [], ["ketone", "ketones", "beta-hydroxybutyrate", "β-hydroxybutyrate", "KET", "BHB"], ["酮体", "β-羟丁酸", "羟丁酸"], None)
    C["osmolality_serum"] = _ch("endocrine_metabolic", "lab", [], ["osmolality", "osmolarity", "OSM"], ["渗透压"], (200, 450), fluid_hint="blood")
    # ---- Supporting audit targets (not Table 5.1 features)
    C["weight"] = _ch("support", "chart", [], ["weight", "body weight", "WT"], ["体重"], (0.3, 200), collect_times=False)
    C["height"] = _ch("support", "chart", [], ["height", "length", "HT"], ["身高", "身长"], (25, 220))
    C["pupils"] = _ch("support", "chart", [], ["pupil", "pupils", "pupillary"], ["瞳孔"], None)
    C["ventilation"] = _ch("support", "both", [], ["ventilator", "ventilation", "mechanical ventilation", "PEEP", "tidal volume", "vent mode", "intubation", "CPAP", "BiPAP", "MAP aw", "PIP"], ["呼吸机", "机械通气", "通气模式", "插管", "潮气量", "呼气末正压"], None)
    C["urine_sodium"] = _ch("support", "lab", [], ["urine sodium", "UNa", "sodium urine"], ["尿钠"], None, fluid_hint="urine")
    return C


def _build_drug_groups() -> Dict[str, Dict[str, List[str]]]:
    """Drug keyword groups (English + Chinese).  'exclude' resolves substring clashes
    (e.g. 肾上腺素 is inside 去甲肾上腺素 and 异丙肾上腺素)."""
    return {
        "vasoactive": {
            "en": ["dopamine", "dobutamine", "epinephrine", "adrenaline", "norepinephrine", "noradrenaline",
                   "milrinone", "vasopressin", "isoproterenol", "isoprenaline", "phenylephrine", "levosimendan",
                   "terlipressin", "pituitrin"],
            "cn": ["多巴胺", "多巴酚丁胺", "肾上腺素", "去甲肾上腺素", "米力农", "血管加压素", "垂体后叶素",
                   "异丙肾上腺素", "去氧肾上腺素", "苯肾上腺素", "左西孟旦", "特利加压素"],
            "exclude": ["肾上腺素受体", "epinephrine dental", "肾上腺色腙"],
        },
        "nmb": {
            "en": ["vecuronium", "rocuronium", "cisatracurium", "atracurium", "succinylcholine", "suxamethonium",
                   "pancuronium", "mivacurium"],
            "cn": ["维库溴铵", "罗库溴铵", "顺阿曲库铵", "阿曲库铵", "琥珀胆碱", "琥珀酰胆碱", "泮库溴铵", "米库氯铵"],
            "exclude": [],
        },
        "antiepileptic": {
            "en": ["levetiracetam", "phenobarbital", "phenobarbitone", "phenytoin", "fosphenytoin", "valproate",
                   "valproic", "carbamazepine", "oxcarbazepine", "topiramate", "lamotrigine", "lacosamide",
                   "clonazepam", "vigabatrin", "zonisamide", "nitrazepam", "clobazam"],
            "cn": ["左乙拉西坦", "苯巴比妥", "苯妥英", "磷苯妥英", "丙戊酸", "卡马西平", "奥卡西平", "托吡酯",
                   "拉莫三嗪", "拉考沙胺", "氯硝西泮", "氨己烯酸", "唑尼沙胺", "硝西泮", "氯巴占"],
            "exclude": [],
        },
        # Dual-use: benzodiazepines used both as sedatives and as rescue antiepileptics.  Reported
        # separately so the sedative vs AED assignment is an explicit decision [decide].
        "benzodiazepine_dual_use": {
            "en": ["midazolam", "diazepam", "lorazepam"],
            "cn": ["咪达唑仑", "地西泮", "劳拉西泮"],
            "exclude": [],
        },
    }


# Table 5.1 pathogen families (ordered flag + positive flag); searched in D_LABITEMS and in
# MICROBIOLOGYEVENTS text columns.  "bacteremia_pcr" naming in PIC is unknown [bracket].
CONFIG["PATHOGEN_FAMILIES"] = {
    "respiratory_panel": {"en": ["RSV", "respiratory syncytial", "parainfluenza", "influenza", "adenovirus",
                                 "metapneumovirus", "hMPV", "mycoplasma"],
                          "cn": ["呼吸道合胞病毒", "副流感", "流感", "流感病毒", "腺病毒", "偏肺病毒", "支原体", "肺炎支原体"]},
    "systemic_viral": {"en": ["CMV", "cytomegalovirus", "EBV", "Epstein-Barr", "HSV", "herpes simplex", "enterovirus"],
                       "cn": ["巨细胞病毒", "EB病毒", "单纯疱疹", "肠道病毒"]},
    "bacteremia_pcr": {"en": ["bacterial DNA", "blood PCR", "bacteremia PCR", "pathogen PCR", "mNGS", "16S"],
                       "cn": ["细菌DNA", "病原体核酸", "宏基因组", "血培养PCR", "细菌核酸"]},
}
CONFIG["CHANNELS"] = _build_channels()
CONFIG["DRUG_GROUPS"] = _build_drug_groups()

# v12 open items each section resolves (plan §18 and v12_change_list.md conflict numbers).
RESOLVES: Dict[str, List[str]] = {
    "A0_cohort": ["cohort and split sizes omitted from this public copy (plan §2)"],
    "A1_interval": ["[handling of negative values: audit] (§4 table; conflict #22)",
                    "[X set at audit] (§4 admission source)", "check A(5) band partition"],
    "A2_ed": ["[EDREGTIME source table] (§18)", "ED rule of inferred admission source (§4)"],
    "A3_surgery": ["surgery at exactly INTIME -> Layer A only (§4; check A(13), B(13))",
                   "[whether a surgery at exactly INTIME counts as OR: decide] (conflict #15)",
                   "[surgery name / procedure-code field] (§18)"],
    "A4_diagnosis_missing": ["[diagnosis missingness indicator: decided at audit] (§4)"],
    "A5_sex": ["[audit: count of empty sex] (§4; conflict #14)"],
    "A6_icd_mapping": ["Chinese-ICD mapping audit (§4 inputs; check A(14))",
                       "[$m = 20$ confirmation] rare 3-char categories (§4; check A(7))",
                       "[chapter = first letter or WHO chapter: confirm] (Appendix A; conflict #3)",
                       "within-(chapter, domain) pooling for split chapters Q, D, I [confirm] (conflict #2)",
                       "Appendix A: chapter I placement, D00-D49 placement (conflict #3)",
                       "check A(8): diagnosis-only domain vs contract domain cross-tab"],
    "A7_careunit": ["first_careunit counts (§4 fixed field)"],
    "A8_prior_icu": ["[prior ICU admissions: include or omit] (§4 candidate)"],
    "A9_admission_source": ["[X set at audit]; admission-source rule counts (§4; check A(12))"],
    "A10_timing": ["check A(5): shifts, day type, L2 timing categories and boundary values (§4)"],
    "B_dictionary": ["[ITEMID list per channel; audit] (§5 inputs, Table 5.1)",
                     "GCS / pupils / FiO2 / ventilation / vasoactive dose / weight search (checks B(25), stage 07a)",
                     "[cortisol, TSH / free T4, insulin, ketones, serum osmolality] (Table 5.1 endocrine)",
                     "stool occult blood (Table 5.1 GI)", "mis-mapped urine sodium vs serum sodium (data-quality flags)",
                     "[troponin I folded or dropped] (Table 5.1)"],
    "C_coverage": ["0-24 h coverage recomputed on this cohort (Table 5.1 caption; conflict #9)",
                   "[GI standalone or merged: audit] (Table 5.1)", "per-bin coverage for LOCF families (§6)"],
    "C_variability": ["variability needs >=3 values (§5 φ1, §6 φ1; check B(21))"],
    "C_units_plausibility": ["unit table and native-unit ranges (§5 φ0; check B(14))",
                             "[direction of 'worst' per lab] context (§18)"],
    "C_panels": ["panel-driven coverage flag (data-quality flags)", "triggered vs routine classification (§5 Table 5.1)"],
    "C_hr_gaps": ["HR coverage 76% suggests CHARTEVENTS gaps (data-quality flags)"],
    "C_weight": ["[weight source and missingness: audit] (Table 5.1 renal; check B(26))"],
    "C_fluids": ["[urine output rate definition] (S2, §5; conflict #5)", "net fluid balance / % fluid overload (Table 5.1)"],
    "D_release": ["S9: which timestamps exist; release = CHARTTIME + turnaround if none (conflict #8)",
                  "pathogen positive flags need release timing (conflict #26)"],
    "E_drugs": ["vasoactive dose availability for pSOFA [audit] (conflict #25)",
                "[neuromuscular blockers and antiepileptics as new drug groups: decide] (conflict #7)",
                "ENDDATE missingness (§17, state machine)"],
}


# =====================================================================================
# Helpers: privacy, column matching, reading
# =====================================================================================
def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sup(n: Any) -> Any:
    """Small-cell suppression for one count."""
    if n is None or (isinstance(n, float) and math.isnan(n)):
        return None
    n = int(n)
    k = CONFIG["SMALL_CELL"]
    if n == 0:
        return f"<{k}" if CONFIG["SUPPRESS_ZERO"] else 0
    return f"<{k}" if n < k else n


def pct(num: Any, den: Any, digits: int = 1) -> Any:
    """Percentage; suppressed when the numerator is a small cell (or the denominator is)."""
    if num is None or den is None or den == 0:
        return None
    num, den = int(num), int(den)
    k = CONFIG["SMALL_CELL"]
    if (0 < num < k) or (0 < den < k) or (num == 0 and CONFIG["SUPPRESS_ZERO"]):
        return "suppressed"
    return round(100.0 * num / den, digits)


def cnt_pct(num: Any, den: Any) -> Dict[str, Any]:
    return {"n": sup(num), "pct": pct(num, den)}


def sup_counts(counter: Dict[Any, int], top: Optional[int] = None) -> Dict[str, Any]:
    """Suppress a {category: count} mapping, sorted by count desc."""
    items = sorted(((str(k), int(v)) for k, v in counter.items()), key=lambda kv: -kv[1])
    if top is not None:
        items = items[:top]
    return {k: sup(v) for k, v in items}


def quantiles(values: Iterable[float]) -> Dict[str, Any]:
    """Quantiles without min/max; suppressed below MIN_N_QUANTILE."""
    arr = np.asarray(list(values) if not isinstance(values, np.ndarray) else values, dtype=float)
    arr = arr[np.isfinite(arr)]
    n = int(arr.size)
    if n < CONFIG["MIN_N_QUANTILE"]:
        return {"n": sup(n), "quantiles": "suppressed (n too small)"}
    qs = np.quantile(arr, CONFIG["QUANTILES"])
    return {"n": sup(n), "mean": round(float(arr.mean()), 4),
            "quantiles": {f"p{int(round(q * 100)):02d}": round(float(v), 4) for q, v in zip(CONFIG["QUANTILES"], qs)}}


def tokens(name: str) -> List[str]:
    """Split a column name into lowercase word tokens (snake, kebab, camel, digits)."""
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", str(name))
    return [t for t in re.split(r"[^A-Za-z0-9]+", s.lower()) if t]


def role_match(colname: str, role_tokens: Sequence[str], exclude_tokens: Sequence[str] = ()) -> bool:
    """Word-boundary role matching.  Fixes the draft's substring misfires:
    'age' no longer matches 'stage'/'dosage', 'dt' no longer matches 'width'/'adt...',
    'los' no longer matches 'close'.  A multi-word role ('begin time') must appear as
    consecutive tokens, or as one fused token ('begintime')."""
    toks = tokens(colname)
    joined = "".join(toks)
    if any(t in toks for t in exclude_tokens):
        return False
    for role in role_tokens:
        r = tokens(role)
        if not r:
            continue
        if len(r) == 1 and r[0] in toks:
            return True
        n = len(r)
        if n > 1 and any(toks[i:i + n] == r for i in range(len(toks) - n + 1)):
            return True
        if n > 1 and "".join(r) in toks:
            return True
        if n > 1 and joined == "".join(r):
            return True
    return False


def find_col(columns: Sequence[str], exact: Sequence[str] = (), roles: Sequence[str] = (),
             exclude: Sequence[str] = ()) -> Optional[str]:
    """Case-insensitive exact name first, then word-boundary role match."""
    up = {str(c).upper(): c for c in columns}
    for e in exact:
        if e.upper() in up:
            return up[e.upper()]
    for c in columns:
        if role_match(c, roles, exclude):
            return c
    return None


def find_cols(columns: Sequence[str], roles: Sequence[str], exclude: Sequence[str] = ()) -> List[str]:
    return [c for c in columns if role_match(c, roles, exclude)]


EN_BOUNDARY = r"(?<![A-Za-z0-9]){}(?![A-Za-z0-9])"


def has_cjk(s: str) -> bool:
    return bool(re.search(r"[\u3400-\u9fff]", s))


def keyword_hits(text: Any, keywords: Sequence[str]) -> List[str]:
    """English keywords match on word boundaries (case-insensitive); keywords containing
    CJK characters match as substrings (Chinese has no word boundaries)."""
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return []
    s = str(text)
    out = []
    for k in keywords:
        if not k:
            continue
        if has_cjk(k):
            if k in s:
                out.append(k)
        elif re.search(EN_BOUNDARY.format(re.escape(k)), s, flags=re.IGNORECASE):
            out.append(k)
    return out


def norm_id(s: pd.Series) -> pd.Series:
    """Normalize ID strings: strip, drop trailing '.0', uppercase (ITEMIDs include 'SV1')."""
    return s.astype("string").str.strip().str.replace(r"\.0$", "", regex=True).str.upper()


def to_dt(s: pd.Series) -> pd.Series:
    fmt = CONFIG["DATETIME_FORMAT"]
    s = s.astype("string").str.strip()
    out = pd.to_datetime(s, format=fmt, errors="coerce") if fmt else pd.to_datetime(s, errors="coerce")
    bad = (out.isna() & s.notna() & (s != "")).fillna(False).astype(bool)
    if bad.any():
        try:
            out.loc[bad] = pd.to_datetime(s[bad], format="mixed", errors="coerce")
        except (TypeError, ValueError):
            out.loc[bad] = pd.to_datetime(s[bad], errors="coerce")
    return out


def to_num(s: pd.Series) -> pd.Series:
    """Numeric float64 with NaN (never a nullable dtype, so .to_numpy(float) is safe)."""
    o = pd.Series(s, copy=False).astype(object)
    o = o.where(pd.notna(o), np.nan)
    return pd.to_numeric(o, errors="coerce").astype(float)


def blank(s: pd.Series) -> pd.Series:
    """True where a value is missing or empty/whitespace/'nan'/'null'."""
    t = s.astype("string").str.strip()
    out = t.isna() | t.eq("") | t.str.lower().isin(["nan", "null", "none", "na", "n/a", "<na>"])
    return out.fillna(True).astype(bool)


# ---------------- dataset discovery (kept from the draft's idea) ----------------
def find_project_root() -> Path:
    if CONFIG["PROJECT_ROOT"]:
        return Path(CONFIG["PROJECT_ROOT"]).expanduser().resolve()
    starts = [Path(__file__).resolve().parent, Path.cwd().resolve()]
    for start in starts:
        for p in [start, *start.parents]:
            if (p / "configs").is_dir() or (p / "data" / "interim").is_dir() or (p / "data" / "processed").is_dir():
                return p
    return Path.cwd().resolve()


def _dir_has_pic(d: Path) -> bool:
    try:
        return any(f.name.upper() in ("ICUSTAYS.CSV", "ICUSTAYS.CSV.GZ") for f in d.iterdir() if f.is_file())
    except OSError:
        return False


def _search_pic(d: Path, depth: int) -> Optional[Path]:
    if not d.is_dir():
        return None
    if _dir_has_pic(d):
        return d
    if depth <= 0:
        return None
    try:
        subs = sorted(p for p in d.iterdir() if p.is_dir() and not p.name.startswith("."))
    except OSError:
        return None
    for s in subs:
        hit = _search_pic(s, depth - 1)
        if hit:
            return hit
    return None


def find_raw_root(project: Path) -> Tuple[Optional[Path], str]:
    env = os.environ.get(CONFIG["RAW_ROOT_ENV_VAR"] or "")
    cands: List[Tuple[Path, str]] = []
    if env:
        cands.append((Path(env), f"env {CONFIG['RAW_ROOT_ENV_VAR']}"))
    if CONFIG["RAW_ROOT_OVERRIDE"]:
        cands.append((Path(CONFIG["RAW_ROOT_OVERRIDE"]), "CONFIG override"))
    ds = project / CONFIG["DATA_SOURCE_CONFIG"]
    if ds.is_file():
        try:
            j = json.loads(ds.read_text(encoding="utf-8"))
            flat: Dict[str, Any] = {}

            def _walk(o: Any, pre: str = "") -> None:
                if isinstance(o, dict):
                    for k, v in o.items():
                        _walk(v, f"{pre}{k}.")
                elif isinstance(o, str):
                    flat[pre.rstrip(".")] = o
            _walk(j)
            for key in CONFIG["DATA_SOURCE_KEYS"]:
                for fk, v in flat.items():
                    if fk.split(".")[-1].lower() == key.lower():
                        p = Path(v)
                        cands.append(((project / p) if not p.is_absolute() else p, f"{CONFIG['DATA_SOURCE_CONFIG']}:{fk}"))
            for fk, v in flat.items():  # any other string that is a directory
                p = Path(v)
                p = (project / p) if not p.is_absolute() else p
                cands.append((p, f"{CONFIG['DATA_SOURCE_CONFIG']}:{fk} (unrecognized key)"))
        except (OSError, ValueError) as e:
            log(f"could not read {ds}: {e}")
    for c in CONFIG["RAW_ROOT_CANDIDATES"]:
        cands.append((project / c, "candidate"))
    for p, how in cands:
        hit = _search_pic(p.expanduser(), CONFIG["RAW_SEARCH_DEPTH"])
        if hit:
            return hit.resolve(), how
    return None, "not found"


def table_path(raw: Path, stem: str) -> Optional[Path]:
    """Case-insensitive lookup of STEM.csv or STEM.csv.gz in the raw root."""
    want = {f"{stem}.csv".upper(), f"{stem}.csv.gz".upper()}
    hits = [f for f in raw.iterdir() if f.is_file() and f.name.upper() in want]
    hits.sort(key=lambda f: (not f.name.lower().endswith(".csv"), f.name))  # prefer plain csv
    return hits[0] if hits else None


def read_header(path: Path) -> List[str]:
    return list(pd.read_csv(path, nrows=0, compression="infer", encoding=CONFIG["ENCODING"],
                            encoding_errors="replace").columns)


def resolve_cols(header: Sequence[str], wanted: Sequence[str]) -> Dict[str, str]:
    """Map canonical UPPERCASE names -> actual header names (case-insensitive)."""
    up = {str(c).strip().upper(): c for c in header}
    return {w: up[w.upper()] for w in wanted if w.upper() in up}


def read_small(path: Path, wanted: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Read a small table as strings; columns renamed to UPPERCASE canonical names."""
    header = read_header(path)
    if wanted is None:
        usecols = header
        ren = {c: str(c).strip().upper() for c in header}
    else:
        m = resolve_cols(header, wanted)
        usecols = list(m.values())
        ren = {v: k for k, v in m.items()}
    df = pd.read_csv(path, usecols=usecols, dtype=str, compression="infer", encoding=CONFIG["ENCODING"],
                     encoding_errors="replace", keep_default_na=True)
    return df.rename(columns=ren)


def iter_chunks(path: Path, wanted: Sequence[str]) -> Iterator[pd.DataFrame]:
    """Chunked reader with usecols; yields frames with UPPERCASE canonical column names."""
    header = read_header(path)
    m = resolve_cols(header, wanted)
    ren = {v: k for k, v in m.items()}
    reader = pd.read_csv(path, usecols=list(m.values()), dtype=str, chunksize=CONFIG["CHUNKSIZE"],
                         compression="infer", encoding=CONFIG["ENCODING"], encoding_errors="replace",
                         low_memory=False)
    for i, chunk in enumerate(reader):
        if i % 5 == 0:
            log(f"  {path.name}: chunk {i} ({(i * CONFIG['CHUNKSIZE']):,} rows read)")
        yield chunk.rename(columns=ren)


class Reservoir:
    """Approximately uniform value sample across chunks (for quantiles; never written)."""

    def __init__(self, cap: int, rng: np.random.Generator):
        self.cap, self.rng, self.seen = cap, rng, 0
        self.buf = np.empty(0, dtype=float)

    def add(self, vals: np.ndarray) -> None:
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return
        self.seen += vals.size
        keep_p = min(1.0, self.cap / max(self.seen, 1))
        if keep_p < 1.0:
            vals = vals[self.rng.random(vals.size) < keep_p]
            if self.buf.size:
                self.buf = self.buf[self.rng.random(self.buf.size) < (self.cap / (self.cap + vals.size + 1))]
        self.buf = np.concatenate([self.buf, vals])
        if self.buf.size > self.cap:
            self.buf = self.rng.choice(self.buf, self.cap, replace=False)


# =====================================================================================
# Domain helpers: ICD-10, domains, timing, admission source, bins
# =====================================================================================
WHO_CHAPTERS = [  # (chapter, first code, last code), WHO ICD-10 2019
    ("I", "A00", "B99"), ("II", "C00", "D48"), ("III", "D50", "D89"), ("IV", "E00", "E90"),
    ("V", "F00", "F99"), ("VI", "G00", "G99"), ("VII", "H00", "H59"), ("VIII", "H60", "H95"),
    ("IX", "I00", "I99"), ("X", "J00", "J99"), ("XI", "K00", "K93"), ("XII", "L00", "L99"),
    ("XIII", "M00", "M99"), ("XIV", "N00", "N99"), ("XV", "O00", "O99"), ("XVI", "P00", "P96"),
    ("XVII", "Q00", "Q99"), ("XVIII", "R00", "R99"), ("XIX", "S00", "T98"), ("XX", "V01", "Y98"),
    ("XXI", "Z00", "Z99"), ("XXII", "U00", "U99"),
]
ICD3_RE = re.compile(r"^[A-Z][0-9]{2}$")


def icd3(code: Any) -> Optional[str]:
    """Three-character ICD-10 category from a code like 'J18.900' or 'j189'; None if invalid."""
    if code is None or (isinstance(code, float) and math.isnan(code)):
        return None
    s = re.sub(r"[^A-Za-z0-9]", "", str(code)).upper()
    if len(s) < 3:
        return None
    c = s[:3]
    return c if ICD3_RE.match(c) else None


def _code_key(c3: str) -> Tuple[str, int]:
    return c3[0], int(c3[1:3])


def who_chapter(c3: Optional[str]) -> Optional[str]:
    if not c3:
        return None
    k = _code_key(c3)
    for ch, lo, hi in WHO_CHAPTERS:
        if _code_key(lo) <= k <= _code_key(hi):
            return ch
    return None


def in_range(c3: str, lo: str, hi: str) -> bool:
    return _code_key(lo) <= _code_key(c3) <= _code_key(hi)


def domain_dx_only(c3: Optional[str]) -> str:
    """Appendix A DRAFT diagnosis-only domain (care unit not used)."""
    if not c3:
        return "other_medical"
    L = c3[0]
    if in_range(c3, "Q20", "Q28"):
        return "cardiac_congenital_heart"
    if L == "I":
        return CONFIG["DOMAIN_CHAPTER_I"]
    if L == "P":
        return "neonatal_perinatal"
    if L in "ST":
        return "surgical_trauma"
    if L in "VWXY":
        return CONFIG["DOMAIN_EXTERNAL_CAUSE_V_Y"]
    if L in "JAB":
        return "respiratory_infection"
    if L == "G":
        return "neurological"
    if L == "C" or in_range(c3, "D50", "D89"):
        return "oncology_hematology"
    if L == "D":
        return CONFIG["DOMAIN_D00_D49"]
    if L in "KN":
        return "gi_renal"
    return "other_medical"


def domain_contract(c3: Optional[str], careunit: Any) -> str:
    """Contract §8.2 first-match domain (uses the care unit)."""
    cu = str(careunit).strip().upper() if careunit is not None and not (isinstance(careunit, float) and math.isnan(careunit)) else ""
    L = c3[0] if c3 else ""
    if cu == "CICU" or (c3 and in_range(c3, "Q20", "Q28")):
        return "cardiac_congenital_heart"
    if cu == "NICU" or L == "P":
        return "neonatal_perinatal"
    if cu == "SICU" or (L and L in "ST"):
        return "surgical_trauma"
    if L and L in "JAB":
        return "respiratory_infection"
    if L == "G":
        return "neurological"
    if L == "C" or (c3 and L == "D" and in_range(c3, "D50", "D89")):
        return "oncology_hematology"
    if L and L in "KN":
        return "gi_renal"
    return "other_medical"


def interval_band(h: float) -> Optional[str]:
    """v12 bands; a boundary value falls in the band whose lower end it equals."""
    if h is None or not np.isfinite(h):
        return None
    e1, e2, e3 = CONFIG["INTERVAL_BAND_EDGES"]
    if h < e1:
        return f"<{e1}h"
    if h < e2:
        return f"{e1}-<{e2}h"
    if h < e3:
        return f"{e2}-<{e3}h"
    return f">={e3}h"


def shift_of(hour: int) -> str:
    d, e, n = CONFIG["DAY_START_HOUR"], CONFIG["EVENING_START_HOUR"], CONFIG["NIGHT_START_HOUR"]
    if d <= hour < e:
        return "day_08-15"
    if e <= hour < n:
        return "evening_15-23"
    return "night_23-08"


def timing_l2(hour: int, dow: int) -> str:
    wk = "weekday" if dow < 5 else "weekend"
    dt = "daytime" if shift_of(hour) == "day_08-15" else "evening_night"
    return f"{wk}_{dt}"


def admission_source(interval_h: np.ndarray, surg_before: np.ndarray, ed_before: np.ndarray,
                     x_hours: float) -> np.ndarray:
    """v12 first-match rule: floor (interval >= X) > OR (surgery before INTIME) >
    ED (EDREGTIME before INTIME) > other.  Missing interval cannot trigger 'floor'."""
    iv = np.asarray(interval_h, dtype=float)
    out = np.full(iv.shape, "other", dtype=object)
    sb = np.asarray(surg_before, dtype=bool)
    eb = np.asarray(ed_before, dtype=bool)
    floor = np.isfinite(iv) & (iv >= x_hours)
    out[eb] = "ED"
    out[sb] = "OR"
    out[floor] = "floor"
    return out


def n_open_bins(intime: pd.Series, win_end: pd.Series) -> np.ndarray:
    """Number of 6 h bins overlapping W = [INTIME, min(OUTTIME, INTIME+168h))."""
    dur = (win_end - intime).dt.total_seconds().to_numpy() / 3600.0
    nb = np.ceil(np.clip(dur, 0, None) / CONFIG["BIN_HOURS"])
    nmax = CONFIG["WINDOW_HOURS"] // CONFIG["BIN_HOURS"]
    return np.clip(np.nan_to_num(nb, nan=0), 0, nmax).astype(int)


def bits_to_bins(mask: np.ndarray, nbins: int) -> np.ndarray:
    """(n_enc,) int64 bitmasks -> (n_enc, nbins) bool matrix."""
    m = np.asarray(mask, dtype=np.int64)
    return ((m[:, None] >> np.arange(nbins, dtype=np.int64)[None, :]) & 1).astype(bool)


def gap_hours(enc_idx: np.ndarray, hours: np.ndarray) -> np.ndarray:
    """Within-encounter gaps between consecutive observation times (sorted)."""
    if len(enc_idx) < 2:
        return np.empty(0)
    order = np.lexsort((hours, enc_idx))
    e, h = np.asarray(enc_idx)[order], np.asarray(hours, dtype=float)[order]
    same = e[1:] == e[:-1]
    return (h[1:] - h[:-1])[same]


# =====================================================================================
# Cohort
# =====================================================================================
def rebuild_cohort_from_icustays(icu: pd.DataFrame) -> pd.DataFrame:
    """Contract §4.2: 24 h running-episode collapse (in memory only)."""
    d = icu.copy()
    d["SUBJECT_ID"] = norm_id(d["SUBJECT_ID"])
    d["INTIME"], d["OUTTIME"] = to_dt(d["INTIME"]), to_dt(d["OUTTIME"])
    d = d[d["SUBJECT_ID"].notna() & d["INTIME"].notna() & d["OUTTIME"].notna() & (d["OUTTIME"] > d["INTIME"])]
    d["_icu_num"] = to_num(d.get("ICUSTAY_ID", pd.Series(index=d.index, dtype=str)))
    d = d.sort_values(["SUBJECT_ID", "INTIME", "_icu_num"])
    rows = []
    for sid, g in d.groupby("SUBJECT_ID", sort=False):
        cur = None
        for r in g.itertuples(index=False):
            if cur is None or r.INTIME > cur["OUTTIME"] + pd.Timedelta(hours=24):
                if cur is not None:
                    rows.append(cur)
                cur = {"SUBJECT_ID": sid, "INTIME": r.INTIME, "OUTTIME": r.OUTTIME,
                       "hadm": {getattr(r, "HADM_ID", None)}}
            else:
                cur["OUTTIME"] = max(cur["OUTTIME"], r.OUTTIME)
                cur["hadm"].add(getattr(r, "HADM_ID", None))
        if cur is not None:
            rows.append(cur)
    out = pd.DataFrame(rows)
    out["HADM_LIST"] = out["hadm"].apply(lambda s: sorted(str(x) for x in s if x is not None and str(x) != "nan"))
    return out.drop(columns=["hadm"])


def load_cohort(project: Path, raw: Path, notes: List[str]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    cc = CONFIG["COHORT_COLUMNS"]
    info: Dict[str, Any] = {}
    cpath = project / CONFIG["COHORT_FILE"]
    if cpath.is_file():
        df = read_small(cpath)
        cols = list(df.columns)
        sid = find_col(cols, [cc["subject_id"]], ["subject id"])
        it = find_col(cols, [cc["intime"]], ["intime", "icu intime"])
        ot = find_col(cols, [cc["outtime"]], ["outtime", "icu outtime"])
        hl = find_col(cols, [cc["hadm_id_list"]], ["hadm id list"])
        hd = find_col(cols, [cc["hadm_id"]], ["hadm id"], exclude=["list"])
        eid = find_col(cols, [cc["encounter_id"]], ["encounter id"])
        if not (sid and it and ot):
            raise SystemExit(f"Cohort file {cpath} lacks SUBJECT_ID / INTIME / OUTTIME [confirm column names in CONFIG['COHORT_COLUMNS']]")
        coh = pd.DataFrame({"SUBJECT_ID": norm_id(df[sid]), "INTIME": to_dt(df[it]), "OUTTIME": to_dt(df[ot])})
        coh["ENCOUNTER_ID"] = df[eid].astype(str) if eid else np.arange(len(df)).astype(str)
        if hl:
            def _parse(x: Any) -> List[str]:
                try:
                    v = json.loads(x) if isinstance(x, str) else []
                    return [re.sub(r"\.0$", "", str(i)).strip() for i in (v if isinstance(v, list) else [v])]
                except ValueError:
                    return [t.strip() for t in re.split(r"[;,|\s]+", str(x)) if t.strip()]
            coh["HADM_LIST"] = df[hl].apply(_parse)
        elif hd:
            coh["HADM_LIST"] = norm_id(df[hd]).apply(lambda x: [] if pd.isna(x) else [str(x)])
        else:
            coh["HADM_LIST"] = [[] for _ in range(len(coh))]
            notes.append("cohort file has no hadm_id_list/HADM_ID; admissions linked by time containment [bracket]")
        info["cohort_source"] = str(CONFIG["COHORT_FILE"])
    elif CONFIG["REBUILD_COHORT_IF_MISSING"]:
        p = table_path(raw, CONFIG["TABLES"]["ICUSTAYS"])
        icu = read_small(p, ["SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "INTIME", "OUTTIME"])
        coh = rebuild_cohort_from_icustays(icu)
        coh["ENCOUNTER_ID"] = np.arange(len(coh)).astype(str)
        info["cohort_source"] = "rebuilt in memory from ICUSTAYS (contract §4.2 24 h collapse); canonical file not found [confirm path]"
        notes.append(info["cohort_source"])
    else:
        raise SystemExit(f"Cohort file not found: {cpath}")

    # ---- split
    coh["SPLIT"] = "unknown"
    spath, lpath = project / CONFIG["SPLIT_FILE"], project / CONFIG["LABELS_FILE"]
    if spath.is_file():
        s = read_small(spath)
        scol = find_col(list(s.columns), CONFIG["SPLIT_COLUMN_CANDIDATES"], ["split"])
        sidc = find_col(list(s.columns), ["SUBJECT_ID"], ["subject id"])
        if scol and sidc:
            m = dict(zip(norm_id(s[sidc]), s[scol].astype(str).str.strip().str.lower()))
            coh["SPLIT"] = coh["SUBJECT_ID"].map(m).fillna("unknown")
            info["split_source"] = str(CONFIG["SPLIT_FILE"])
    if (coh["SPLIT"] == "unknown").all() and lpath.is_file():
        s = read_small(lpath)
        scol = find_col(list(s.columns), CONFIG["SPLIT_COLUMN_CANDIDATES"], ["split"])
        eidc = find_col(list(s.columns), ["encounter_id"], ["encounter id"])
        if scol and eidc:
            m = dict(zip(s[eidc].astype(str), s[scol].astype(str).str.strip().str.lower()))
            coh["SPLIT"] = coh["ENCOUNTER_ID"].map(m).fillna("unknown")
            info["split_source"] = str(CONFIG["LABELS_FILE"])
    coh["SPLIT"] = coh["SPLIT"].replace({"validation": "val", "valid": "val", "dev": "val"})
    if (coh["SPLIT"] == "unknown").all():
        notes.append("split file not found or unreadable: train-only statistics fall back to all encounters [bracket: split path]")
        info["split_source"] = "none"

    coh = coh.sort_values(["SUBJECT_ID", "INTIME"]).reset_index(drop=True)
    coh["ENC_IDX"] = np.arange(len(coh))
    coh["WIN_END"] = coh[["OUTTIME"]].assign(c=coh["INTIME"] + pd.Timedelta(hours=CONFIG["WINDOW_HOURS"])).min(axis=1)
    coh["N_OPEN_BINS"] = n_open_bins(coh["INTIME"], coh["WIN_END"])
    info["n_encounters"] = int(len(coh))
    info["n_subjects"] = int(coh["SUBJECT_ID"].nunique())
    exp_n = CONFIG.get("EXPECTED_N_ENCOUNTERS")
    info["matches_expected_n"] = None if exp_n is None else bool(len(coh) == exp_n)
    sc = coh["SPLIT"].value_counts().to_dict()
    info["split_counts"] = {k: int(v) for k, v in sc.items()}
    exp_s = CONFIG.get("EXPECTED_SPLIT_COUNTS") or {}
    info["split_matches_expected"] = None if not exp_s else all(sc.get(k, 0) == v for k, v in exp_s.items())
    info["invalid_times"] = int((coh["INTIME"].isna() | coh["OUTTIME"].isna()).sum())
    return coh, info


def train_mask(coh: pd.DataFrame) -> np.ndarray:
    m = (coh["SPLIT"] == "train").to_numpy()
    return m if m.any() else np.ones(len(coh), dtype=bool)


def attribute(df: pd.DataFrame, enc: pd.DataFrame, time_col: str = "CHARTTIME",
              pre_hours: Optional[float] = None) -> pd.DataFrame:
    """Attach rows to encounters by SUBJECT_ID and time: keep INTIME - pre_hours <= t < WIN_END.
    Encounters of one subject are > 24 h apart (contract §4.2), so the windows are disjoint.
    Adds ENC_IDX and H (hours from INTIME)."""
    pre = CONFIG["PRE_HOURS"] if pre_hours is None else pre_hours
    d = df.copy()
    d["SUBJECT_ID"] = norm_id(d["SUBJECT_ID"])
    d["_T"] = to_dt(d[time_col])
    d = d[d["_T"].notna()]
    m = d.merge(enc[["SUBJECT_ID", "ENC_IDX", "INTIME", "WIN_END"]], on="SUBJECT_ID", how="inner")
    m["H"] = (m["_T"] - m["INTIME"]).dt.total_seconds() / 3600.0
    m = m[(m["H"] >= -pre) & (m["_T"] < m["WIN_END"])]
    return m


# =====================================================================================
# Section A: Layer A
# =====================================================================================
def header_scan(raw: Path) -> Dict[str, List[str]]:
    out = {}
    for f in sorted(raw.iterdir()):
        if f.is_file() and re.search(r"\.csv(\.gz)?$", f.name, re.I):
            try:
                out[f.name] = read_header(f)
            except Exception as e:  # noqa: BLE001
                out[f.name] = [f"<unreadable: {type(e).__name__}>"]
    return out


def section_A(raw: Path, coh: pd.DataFrame, headers: Dict[str, List[str]], notes: List[str]) -> Dict[str, Any]:
    T = CONFIG["TABLES"]
    R: Dict[str, Any] = {}
    n = len(coh)
    tm = train_mask(coh)

    # ---------- admissions ----------
    ap = table_path(raw, T["ADMISSIONS"])
    adm_hdr = read_header(ap)
    ed_col = find_col(adm_hdr, ["EDREGTIME"], ["ed reg time", "edreg time", "ed registration time", "ed reg"])
    edout_col = find_col(adm_hdr, ["EDOUTTIME"], ["ed out time", "edout time"])
    icd_col = find_col(adm_hdr, ["ICD10_CODE_CN"], ["icd10 code cn", "icd code cn"])
    dxt_col = find_col(adm_hdr, ["DIAGNOSIS"], ["diagnosis"], exclude=["code", "icd"])
    want = ["SUBJECT_ID", "HADM_ID", "ADMITTIME", "DISCHTIME"] + [c for c in [ed_col, edout_col, icd_col, dxt_col] if c]
    adm = read_small(ap, want)
    adm["SUBJECT_ID"], adm["HADM_ID"] = norm_id(adm["SUBJECT_ID"]), norm_id(adm["HADM_ID"])
    adm["ADMITTIME"] = to_dt(adm["ADMITTIME"])
    ED = ed_col.upper() if ed_col else None
    if ED:
        adm[ED] = to_dt(adm[ED])

    # explode cohort -> linked admissions
    link = coh[["ENC_IDX", "SUBJECT_ID", "INTIME", "WIN_END", "HADM_LIST"]].explode("HADM_LIST").rename(columns={"HADM_LIST": "HADM_ID"})
    link = link[link["HADM_ID"].notna()]
    la = link.merge(adm, on=["SUBJECT_ID", "HADM_ID"], how="inner")
    if la.empty:  # fallback: admissions containing INTIME
        notes.append("no hadm_id_list match; admissions linked by ADMITTIME <= INTIME <= DISCHTIME [bracket]")
        adm["DISCHTIME"] = to_dt(adm["DISCHTIME"])
        la = coh[["ENC_IDX", "SUBJECT_ID", "INTIME", "WIN_END"]].merge(adm, on="SUBJECT_ID")
        la = la[(la["ADMITTIME"] <= la["INTIME"]) & (la["INTIME"] <= la["DISCHTIME"] + pd.Timedelta(hours=24))]
    la = la.sort_values(["ENC_IDX", "ADMITTIME"])
    first = la.groupby("ENC_IDX").head(1).set_index("ENC_IDX")
    min_admit = la.groupby("ENC_IDX")["ADMITTIME"].min()
    coh = coh.copy()
    coh["MIN_ADMIT"] = coh["ENC_IDX"].map(min_admit)
    interval = (coh["INTIME"] - coh["MIN_ADMIT"]).dt.total_seconds() / 3600.0
    coh["INTERVAL_H"] = interval

    # ---------- A0 cohort linkage ----------
    R["A0_linkage"] = {"resolves": RESOLVES["A0_cohort"],
                       "encounters_with_linked_admission": cnt_pct(coh["MIN_ADMIT"].notna().sum(), n)}

    # ---------- A1 interval ----------
    iv = interval.to_numpy(dtype=float)
    fin = np.isfinite(iv)
    neg = fin & (iv < 0)
    bands = Counter(interval_band(x) for x in iv[fin])
    R["A1_interval"] = {
        "resolves": RESOLVES["A1_interval"],
        "definition": "INTIME - min(linked ADMITTIME), hours (contract §8.1)",
        "missing": cnt_pct((~fin).sum(), n),
        "negative": cnt_pct(neg.sum(), n),
        "negative_by_magnitude": {
            "-1h to <0": cnt_pct((fin & (iv < 0) & (iv >= -1)).sum(), n),
            "-24h to <-1h": cnt_pct((fin & (iv < -1) & (iv >= -24)).sum(), n),
            "< -24h": cnt_pct((fin & (iv < -24)).sum(), n)},
        "exactly_zero": cnt_pct((fin & (iv == 0)).sum(), n),
        "distribution_hours_all": quantiles(iv[fin]),
        "distribution_hours_nonnegative": quantiles(iv[fin & (iv >= 0)]),
        "v12_bands": {k: cnt_pct(v, fin.sum()) for k, v in sorted(bands.items(), key=lambda kv: str(kv[0]))},
        "share_at_or_above_X": {f"X={x}h": cnt_pct((fin & (iv >= x)).sum(), n) for x in CONFIG["X_GRID_HOURS"]},
        "boundary_hits_exact": {f"{e}h": sup(int((fin & (iv == e)).sum())) for e in CONFIG["INTERVAL_BAND_EDGES"]},
    }

    # ---------- A2 ED registration ----------
    ed_tables = {f: [c for c in cols if role_match(c, ["edregtime", "ed reg time", "ed registration", "edreg"])]
                 for f, cols in headers.items()}
    ed_tables = {f: c for f, c in ed_tables.items() if c}
    A2: Dict[str, Any] = {"resolves": RESOLVES["A2_ed"], "tables_with_ED_registration_column": ed_tables,
                          "column_used": f"ADMISSIONS.{ed_col}" if ed_col else None}
    ed_before = np.zeros(n, dtype=bool)
    if ED:
        if CONFIG["ED_SOURCE_ROW"] == "any":
            edv = la.groupby("ENC_IDX")[ED].min()
        else:
            edv = first[ED]
        e = coh["ENC_IDX"].map(edv)
        lag = (coh["INTIME"] - e).dt.total_seconds() / 3600.0
        has = e.notna().to_numpy()
        ed_before = (has & (lag.to_numpy() > 0))
        A2.update({
            "row_rule": CONFIG["ED_SOURCE_ROW"],
            "present": cnt_pct(has.sum(), n),
            "before_INTIME": cnt_pct(ed_before.sum(), n),
            "exactly_at_INTIME": cnt_pct((has & (lag.to_numpy() == 0)).sum(), n),
            "after_INTIME": cnt_pct((has & (lag.to_numpy() < 0)).sum(), n),
            "hours_ED_reg_to_INTIME_when_before": quantiles(lag.to_numpy()[ed_before]),
            "ED_reg_before_min_ADMITTIME": cnt_pct(((e < coh["MIN_ADMIT"]) & e.notna()).sum(), n),
        })
    else:
        notes.append("no ED registration column found in ADMISSIONS [bracket: EDREGTIME source table]")
    R["A2_ed"] = A2

    # ---------- A3 surgery vs INTIME ----------
    sp = table_path(raw, T["SURGERY_INFO"])
    surg_before = np.zeros(n, dtype=bool)
    surg_at = np.zeros(n, dtype=bool)
    if sp:
        shdr = read_header(sp)
        bcol = find_col(shdr, ["SURGERY_BEGIN_TIME"], ["surgery begin time", "begin time", "oper begin time", "start time"])
        name_cols = find_cols(shdr, ["surgery name", "operation name", "oper name", "procedure", "name", "icd", "code"],
                              exclude=["id", "subject", "hadm", "visit", "doctor", "anesthesia", "anes"])
        want = ["SUBJECT_ID", "HADM_ID"] + ([bcol] if bcol else []) + name_cols[:4]
        su = read_small(sp, want)
        su["SUBJECT_ID"] = norm_id(su["SUBJECT_ID"])
        if "HADM_ID" in su:
            su["HADM_ID"] = norm_id(su["HADM_ID"])
        A3: Dict[str, Any] = {"resolves": RESOLVES["A3_surgery"], "begin_column": bcol,
                              "name_or_code_columns_found": name_cols, "rows": sup(len(su))}
        if bcol:
            B = bcol.upper()
            su[B] = to_dt(su[B])
            A3["begin_time_missing"] = cnt_pct(su[B].isna().sum(), len(su))
            if "HADM_ID" in su and link["HADM_ID"].notna().any():
                sl = link[["ENC_IDX", "SUBJECT_ID", "HADM_ID", "INTIME", "WIN_END"]].merge(su, on=["SUBJECT_ID", "HADM_ID"])
                A3["linkage"] = "SUBJECT_ID + linked HADM_ID"
            else:
                sl = coh[["ENC_IDX", "SUBJECT_ID", "INTIME", "WIN_END", "MIN_ADMIT"]].merge(su, on="SUBJECT_ID")
                sl = sl[(sl[B] >= sl["MIN_ADMIT"]) & (sl[B] < sl["WIN_END"])]
                A3["linkage"] = "SUBJECT_ID + begin time in [min ADMITTIME, window end) [bracket]"
            sl = sl[sl[B].notna()]
            dh = (sl[B] - sl["INTIME"]).dt.total_seconds() / 3600.0
            cat = np.select([dh < 0, dh == 0, (dh > 0) & (dh < 24), dh >= 24],
                            ["before_INTIME", "exactly_at_INTIME", "after_INTIME_<24h", "after_INTIME_>=24h"], "other")
            sl = sl.assign(_cat=cat, _dh=dh)
            A3["surgeries_by_timing"] = {k: sup(v) for k, v in Counter(cat).items()}
            enc_any = {k: set(sl.loc[sl["_cat"] == k, "ENC_IDX"]) for k in ["before_INTIME", "exactly_at_INTIME", "after_INTIME_<24h"]}
            A3["encounters_with"] = {k: cnt_pct(len(v), n) for k, v in enc_any.items()}
            surg_before[list(enc_any["before_INTIME"])] = True
            surg_at[list(enc_any["exactly_at_INTIME"])] = True
            near = {}
            for mins in CONFIG["SURGERY_NEAR_INTIME_MINUTES"]:
                near[f"within_±{mins}min_of_INTIME"] = sup(int((dh.abs() <= mins / 60.0).sum()))
            A3["timestamp_rounding_check"] = near
            A3["hours_before_INTIME_when_before"] = quantiles(-dh[dh < 0].to_numpy())
            for nc in name_cols[:2]:
                A3[f"top_values_{nc}"] = sup_counts(Counter(su[nc.upper()].dropna().astype(str)), CONFIG["TOP_N_NAMES"])
                A3[f"distinct_values_{nc}"] = int(su[nc.upper()].nunique())
        else:
            notes.append("SURGERY_INFO has no begin-time column [bracket]")
        R["A3_surgery"] = A3
    else:
        notes.append("SURGERY_INFO table not found [bracket]")

    # ---------- A4/A6 diagnosis ----------
    codes = coh["ENC_IDX"].map(first[icd_col.upper()]) if icd_col else pd.Series([None] * n)
    miss = blank(codes)
    A4: Dict[str, Any] = {"resolves": RESOLVES["A4_diagnosis_missing"], "column": f"ADMISSIONS.{icd_col} (earliest linked admission)",
                          "missing": cnt_pct(miss.sum(), n),
                          "missing_by_split": {s: cnt_pct((miss & (coh["SPLIT"] == s)).sum(), (coh["SPLIT"] == s).sum())
                                               for s in sorted(coh["SPLIT"].unique())},
                          "multiple_codes_in_cell": cnt_pct(codes.astype("string").str.contains(r"[;,|]", regex=True, na=False).sum(), n)}
    if dxt_col:
        dtext = coh["ENC_IDX"].map(first[dxt_col.upper()])
        A4["free_text_DIAGNOSIS_present_when_code_missing"] = cnt_pct((miss & ~blank(dtext)).sum(), max(int(miss.sum()), 1))
    R["A4_diagnosis_missing"] = A4

    dp = table_path(raw, T["D_ICD_DIAGNOSES"])
    A6: Dict[str, Any] = {"resolves": RESOLVES["A6_icd_mapping"]}
    c3 = pd.Series([None] * n, dtype=object)
    if dp is not None:
        dic = read_small(dp)
        cn = find_col(list(dic.columns), ["ICD10_CODE_CN"], ["icd10 code cn"])
        en = find_col(list(dic.columns), ["ICD10_CODE"], ["icd10 code"], exclude=["cn"])
        if cn and en:
            dk = dic[[cn, en]].copy()
            dk["_k"] = dk[cn].astype("string").str.strip().str.upper()
            dup = dk.groupby("_k")[en].nunique()
            mp = dk.drop_duplicates("_k").set_index("_k")[en]
            key = codes.astype("string").str.strip().str.upper()
            mapped = key.map(mp)
            present = ~miss
            unm = present & mapped.isna()
            c3 = mapped.apply(icd3)
            invalid = present & mapped.notna() & c3.isna()
            cn3 = key.apply(icd3)
            A6.update({
                "dictionary_rows": sup(len(dic)), "dictionary_cn_codes_with_multiple_icd10": sup(int((dup > 1).sum())),
                "codes_present": cnt_pct(present.sum(), n),
                "unmapped_codes_present_not_in_dictionary": cnt_pct(unm.sum(), n),
                "unmapped_top_code_prefixes": sup_counts(Counter(key[unm].str[:3].dropna()), CONFIG["TOP_N_NAMES"]),
                "mapped_but_invalid_icd10_format": cnt_pct(invalid.sum(), n),
                "valid_3char_category": cnt_pct(c3.notna().sum(), n),
                "cn_code_3char_equals_icd10_3char": cnt_pct(((cn3 == c3) & c3.notna()).sum(), int(c3.notna().sum())),
                "unmapped_but_cn_code_looks_like_icd10": cnt_pct((unm & cn3.notna()).sum(), n),
            })
        else:
            notes.append("D_ICD_DIAGNOSES lacks ICD10_CODE_CN / ICD10_CODE [bracket]")
    coh["C3"] = c3.values
    coh["LETTER"] = coh["C3"].str[0]
    coh["WHO"] = coh["C3"].apply(who_chapter)
    coh["DOM"] = coh["C3"].apply(domain_dx_only)

    # rare categories (train patients)
    trn = coh[tm & coh["C3"].notna()]
    pat = trn.groupby("C3")["SUBJECT_ID"].nunique()
    kmin = CONFIG["RARE_CATEGORY_MIN_TRAIN_PATIENTS"]
    rare = set(pat[pat < kmin].index)
    vocab = set(pat.index)
    nontrain = coh[~tm & coh["C3"].notna()]
    A6["rare_categories"] = {
        "threshold_train_patients": kmin,
        "train_categories": len(vocab),
        "categories_below_threshold": len(rare),
        "categories_kept": len(vocab) - len(rare),
        "encounters_in_rare_categories_all_splits": cnt_pct(coh["C3"].isin(rare).sum(), n),
        "sensitivity_thresholds": {f"<{k}": int((pat < k).sum()) for k in (10, 20, 30, 50)},
        "val_test_encounters_with_unseen_category": cnt_pct((~nontrain["C3"].isin(vocab)).sum(), len(nontrain)),
    }
    # chapter question: first letter vs WHO chapter
    xt = Counter(zip(coh["LETTER"].dropna(), coh.loc[coh["LETTER"].notna(), "WHO"].fillna("none")))
    letters_multi = defaultdict(set)
    who_multi = defaultdict(set)
    for (L, W), _ in xt.items():
        letters_multi[L].add(W)
        who_multi[W].add(L)
    A6["chapter_definition"] = {
        "letter_x_WHO_chapter_counts": {f"{L}|{W}": sup(v) for (L, W), v in sorted(xt.items())},
        "letters_spanning_WHO_chapters": {L: sorted(v) for L, v in letters_multi.items() if len(v) > 1},
        "WHO_chapters_spanning_letters": {W: sorted(v) for W, v in who_multi.items() if len(v) > 1},
        "D00_D48_vs_D50_D89": {"D00-D48": sup(int(coh["C3"].dropna().apply(lambda c: c[0] == "D" and in_range(c, "D00", "D48")).sum())),
                               "D49": sup(int((coh["C3"] == "D49").sum())),
                               "D50-D89": sup(int(coh["C3"].dropna().apply(lambda c: c[0] == "D" and in_range(c, "D50", "D89")).sum()))},
        "H00_H59_vs_H60_H95": {"H00-H59": sup(int(coh["C3"].dropna().apply(lambda c: c[0] == "H" and in_range(c, "H00", "H59")).sum())),
                               "H60-H95": sup(int(coh["C3"].dropna().apply(lambda c: c[0] == "H" and in_range(c, "H60", "H95")).sum()))},
    }
    # chapters spanning domains (Q, D, I and any other)
    cells = Counter(zip(coh["LETTER"].dropna(), coh.loc[coh["LETTER"].notna(), "DOM"]))
    span = defaultdict(dict)
    for (L, D), v in cells.items():
        span[L][D] = v
    split_letters = {L: d for L, d in span.items() if len(d) > 1}
    rare_cells = Counter((c[0], domain_dx_only(c)) for c in rare)
    A6["chapters_spanning_domains"] = {
        "domain_options": {"chapter_I": CONFIG["DOMAIN_CHAPTER_I"], "D00_D49": CONFIG["DOMAIN_D00_D49"],
                           "V_Y": CONFIG["DOMAIN_EXTERNAL_CAUSE_V_Y"]},
        "encounters_per_letter_domain_cell": {L: {D: sup(v) for D, v in d.items()} for L, d in sorted(split_letters.items())},
        "rare_categories_per_cell": {f"{L}|{D}": v for (L, D), v in sorted(rare_cells.items())},
        "letters_whose_rare_categories_span_domains": sorted({L for (L, D) in rare_cells if len([1 for (l2, _) in rare_cells if l2 == L]) > 1}),
        "chapter_I_encounters_by_careunit": "see A7 cross-tab",
    }
    R["A6_icd_mapping"] = A6

    # ---------- A5 sex ----------
    pp = table_path(raw, T["PATIENTS"])
    pat_df = read_small(pp, ["SUBJECT_ID", "GENDER", "DOB"])
    pat_df["SUBJECT_ID"] = norm_id(pat_df["SUBJECT_ID"])
    g = coh["SUBJECT_ID"].map(dict(zip(pat_df["SUBJECT_ID"], pat_df["GENDER"].astype("string").str.strip().str.upper())))
    R["A5_sex"] = {"resolves": RESOLVES["A5_sex"],
                   "values": sup_counts(Counter(g.fillna("<missing>"))),
                   "missing_or_empty": cnt_pct(blank(g).sum(), n),
                   "not_M_or_F": cnt_pct((~g.isin(["M", "F"])).sum(), n)}

    # ---------- A7 care unit ----------
    ip = table_path(raw, T["ICUSTAYS"])
    icu = read_small(ip, ["SUBJECT_ID", "ICUSTAY_ID", "INTIME", "FIRST_CAREUNIT"])
    icu["SUBJECT_ID"] = norm_id(icu["SUBJECT_ID"])
    icu["INTIME"] = to_dt(icu["INTIME"])
    icu["_n"] = to_num(icu["ICUSTAY_ID"])
    first_cu = coh[["ENC_IDX", "SUBJECT_ID", "INTIME"]].merge(icu, on=["SUBJECT_ID", "INTIME"], how="left")
    first_cu = first_cu.sort_values(["ENC_IDX", "_n"]).groupby("ENC_IDX").head(1).set_index("ENC_IDX")["FIRST_CAREUNIT"]
    cu = coh["ENC_IDX"].map(first_cu).astype("string").str.strip().str.upper()
    coh["CU"] = cu
    contract_dom = [domain_contract(a, b) for a, b in zip(coh["C3"], coh["CU"])]
    cross = Counter(zip(contract_dom, coh["DOM"]))
    R["A7_careunit"] = {"resolves": RESOLVES["A7_careunit"] + RESOLVES["A6_icd_mapping"][-1:],
                        "first_careunit_counts": sup_counts(Counter(cu.fillna("<missing>"))),
                        "contract_domain_x_dx_only_domain": {f"{a} | {b}": sup(v) for (a, b), v in sorted(cross.items())},
                        "agreement": cnt_pct(sum(v for (a, b), v in cross.items() if a == b), n),
                        "chapter_I_by_careunit": sup_counts(Counter(cu[coh["LETTER"] == "I"].fillna("<missing>"))),
                        "Q20_Q28_by_careunit": sup_counts(Counter(cu[coh["C3"].fillna("").apply(lambda c: bool(c) and in_range(c, "Q20", "Q28"))].fillna("<missing>")))}

    # ---------- A8 prior ICU admissions ----------
    prior = coh.groupby("SUBJECT_ID").cumcount().to_numpy()
    icu_valid = icu[icu["INTIME"].notna()]
    raw_prior = coh[["ENC_IDX", "SUBJECT_ID", "INTIME"]].merge(icu_valid[["SUBJECT_ID", "INTIME"]].rename(columns={"INTIME": "T2"}), on="SUBJECT_ID")
    raw_prior = raw_prior[raw_prior["T2"] < raw_prior["INTIME"]].groupby("ENC_IDX").size()
    rp = coh["ENC_IDX"].map(raw_prior).fillna(0).astype(int).to_numpy()
    R["A8_prior_icu"] = {"resolves": RESOLVES["A8_prior_icu"],
                         "prior_encounters_in_cohort": {k: cnt_pct(v, n) for k, v in
                                                        {"0": (prior == 0).sum(), "1": (prior == 1).sum(), "2": (prior == 2).sum(), ">=3": (prior >= 3).sum()}.items()},
                         "prior_raw_icu_stays": {k: cnt_pct(v, n) for k, v in
                                                 {"0": (rp == 0).sum(), "1": (rp == 1).sum(), ">=2": (rp >= 2).sum()}.items()},
                         "note": "splits are subject-grouped, so the count never crosses a split"}

    # ---------- A9 admission source ----------
    A9: Dict[str, Any] = {"resolves": RESOLVES["A9_admission_source"], "by_X": {}, "by_X_surgery_at_INTIME_counts_as_OR": {}}
    for x in CONFIG["X_GRID_HOURS"]:
        src = admission_source(iv, surg_before, ed_before, x)
        A9["by_X"][f"X={x}h"] = {k: cnt_pct(v, n) for k, v in Counter(src).items()}
        src2 = admission_source(iv, surg_before | surg_at, ed_before, x)
        A9["by_X_surgery_at_INTIME_counts_as_OR"][f"X={x}h"] = {k: cnt_pct(v, n) for k, v in Counter(src2).items()}
        A9.setdefault("conflicts", {})[f"X={x}h"] = {
            "ED_before_and_interval>=X": cnt_pct((ed_before & fin & (iv >= x)).sum(), n),
            "surgery_before_and_interval>=X": cnt_pct((surg_before & fin & (iv >= x)).sum(), n)}
    A9["missing_interval_cannot_be_floor"] = cnt_pct((~fin).sum(), n)
    R["A9_admission_source"] = A9

    # ---------- A10 timing ----------
    hr = coh["INTIME"].dt.hour
    dow = coh["INTIME"].dt.dayofweek
    ok = hr.notna()
    sh = Counter(shift_of(int(h)) for h in hr[ok])
    l2 = Counter(timing_l2(int(h), int(d)) for h, d in zip(hr[ok], dow[ok]))
    exact = coh["INTIME"].dt.minute.eq(0) & coh["INTIME"].dt.second.eq(0)
    R["A10_timing"] = {"resolves": RESOLVES["A10_timing"],
                       "L1_shift": {k: cnt_pct(v, int(ok.sum())) for k, v in sh.items()},
                       "L1_day_type": {"weekday": cnt_pct((dow < 5).sum(), n), "weekend": cnt_pct((dow >= 5).sum(), n)},
                       "L2_categories": {k: cnt_pct(v, int(ok.sum())) for k, v in l2.items()},
                       "boundary_exact_08_15_23": {f"{b:02d}:00:00": sup(int((exact & hr.eq(b)).sum())) for b in (8, 15, 23)},
                       "intime_on_the_hour_share": cnt_pct(exact.sum(), n)}
    return R, coh


# =====================================================================================
# Section B: dictionary search (D_LABITEMS, D_ITEMS)
# =====================================================================================
def _text_cols(df: pd.DataFrame) -> List[str]:
    """Label-like columns: English and Chinese/pinyin labels, abbreviations, fluid, category, units."""
    roles = ["label", "name", "abbreviation", "abbr", "fluid", "category", "unit", "unitname", "specimen",
             "loinc", "desc", "description", "cn", "en", "pinyin", "py"]
    cols = [c for c in df.columns if role_match(c, roles, exclude=["id"])]
    return cols or [c for c in df.columns if c != "ITEMID"]


def section_B(raw: Path) -> Tuple[Dict[str, Any], Dict[str, List[Tuple[str, str]]], pd.DataFrame]:
    """Returns (report, auto_candidates {channel: [(source, itemid)]}, dictionary frame)."""
    T = CONFIG["TABLES"]
    frames = []
    for tab, src in [("D_LABITEMS", "lab"), ("D_ITEMS", "chart")]:
        p = table_path(raw, T[tab])
        if p is None:
            continue
        d = read_small(p)
        if "ITEMID" not in d.columns:
            continue
        d["ITEMID"] = norm_id(d["ITEMID"])
        d["_TABLE"] = tab
        d["_SRC"] = src
        if tab == "D_ITEMS":  # LINKSTO says which event table the item lives in
            lk = find_col(list(d.columns), ["LINKSTO"], ["links to", "linksto"])
            if lk:
                lt = d[lk].astype("string").str.lower().fillna("")
                d.loc[lt.str.contains("output"), "_SRC"] = "output"
                d.loc[lt.str.contains("input"), "_SRC"] = "input"
                d.loc[lt.str.contains("lab"), "_SRC"] = "lab"
        frames.append(d)
    if not frames:
        return {"error": "no dictionary tables found"}, {}, pd.DataFrame()
    dic = pd.concat(frames, ignore_index=True, sort=False)
    tcols = sorted({c for f in frames for c in _text_cols(f)})
    dic["_TEXT"] = dic[tcols].astype("string").fillna("").agg(" | ".join, axis=1)
    fluid_col = find_col(tcols, ["FLUID"], ["fluid", "specimen", "sample type"])
    unit_col = find_col(tcols, ["UNITNAME", "UNIT"], ["unit", "unitname", "valueuom"])

    REPORT: Dict[str, Any] = {"resolves": RESOLVES["B_dictionary"], "text_columns_searched": tcols,
                              "items_in_dictionaries": {f["_TABLE"].iloc[0]: int(len(f)) for f in frames},
                              "channels": {}}
    auto: Dict[str, List[Tuple[str, str]]] = {}
    fluid_words = {"urine": ["urine", "尿"], "csf": ["csf", "cerebrospinal", "脑脊液"],
                   "blood": ["blood", "serum", "plasma", "血", "血清", "血浆"]}
    records = dic[list(dict.fromkeys(["ITEMID", "_TABLE", "_SRC", "_TEXT"] + [c for c in (fluid_col, unit_col) if c]))].to_dict("records")
    for ch, spec in CONFIG["CHANNELS"].items():
        kw = spec["keywords_en"] + spec["keywords_cn"]
        hits = []
        for r in records:
            h = keyword_hits(r["_TEXT"], kw)
            if h:
                hits.append((r, h))
        cands = []
        for r, h in hits:
            rec = {"table": r["_TABLE"], "source": r["_SRC"], "itemid": str(r["ITEMID"]),
                   "text": str(r["_TEXT"])[:160], "matched": h}
            if fluid_col:
                rec["fluid"] = str(r.get(fluid_col, ""))
            if unit_col:
                rec["unit"] = str(r.get(unit_col, ""))
            # fluid consistency flag
            if spec.get("fluid_hint"):
                txt = rec["text"].lower() + " " + rec.get("fluid", "").lower()
                want = any(w in txt for w in fluid_words[spec["fluid_hint"]])
                other = [k for k, ws in fluid_words.items() if k != spec["fluid_hint"] and k != "blood" and any(w in txt for w in ws)]
                rec["fluid_flag"] = "matches hint" if want and not other else ("other fluid: " + ",".join(other) if other else "unknown fluid")
            # strength: number of distinct keywords hit, CJK hit, and full-label hit
            rec["score"] = len(set(h)) + (1 if any(has_cjk(x) for x in h) else 0)
            cands.append(rec)
        cands.sort(key=lambda c: -c["score"])
        configured = set(spec["itemids"])
        status = ("configured (contract)" if configured else
                  ("not found" if not cands else ("single candidate (review)" if len(cands) == 1 else f"ambiguous: {len(cands)} candidates (review)")))
        if configured:
            missing_cfg = sorted(configured - set(dic["ITEMID"]))
            if missing_cfg:
                status += f"; configured ITEMIDs absent from dictionaries: {missing_cfg}"
        REPORT["channels"][ch] = {"system": spec["system"], "triggered": spec["triggered"],
                                  "configured_itemids": sorted(configured), "status": status,
                                  "n_candidates": len(cands), "candidates": cands[:25]}
        if not configured and cands:
            if spec.get("fluid_hint") == "blood":
                pref = [c for c in cands if "urine" not in c.get("fluid_flag", "") and "csf" not in c.get("fluid_flag", "")]
            elif spec.get("fluid_hint"):
                pref = [c for c in cands if c.get("fluid_flag") == "matches hint"]
            else:
                pref = cands
            auto[ch] = [(c["source"], c["itemid"]) for c in (pref or cands)[:CONFIG["MAX_AUTO_CANDIDATES_PER_CHANNEL"]]]
    # pathogen families in the lab dictionary
    REPORT["pathogen_families_in_dictionaries"] = {}
    for fam, kw in CONFIG["PATHOGEN_FAMILIES"].items():
        hits = [{"table": r["_TABLE"], "itemid": str(r["ITEMID"]), "text": str(r["_TEXT"])[:160],
                 "matched": keyword_hits(r["_TEXT"], kw["en"] + kw["cn"])}
                for r in records if keyword_hits(r["_TEXT"], kw["en"] + kw["cn"])]
        REPORT["pathogen_families_in_dictionaries"][fam] = {"n_candidates": len(hits), "candidates": hits[:40]}
    # serum vs urine sodium side by side
    na = REPORT["channels"].get("sodium", {}).get("candidates", []) + REPORT["channels"].get("urine_sodium", {}).get("candidates", [])
    REPORT["sodium_serum_vs_urine"] = [{k: c.get(k) for k in ("table", "itemid", "text", "fluid", "unit", "fluid_flag")} for c in na]
    REPORT["note"] = ("Candidates are keyword hits, not decisions.  Section C reports 0-24 h coverage for each "
                      "candidate ITEMID so the true item can be chosen; record the choice in CONFIG['CHANNELS'].")
    return REPORT, auto, dic


# =====================================================================================
# Sections C and D: event pass (CHARTEVENTS, LABEVENTS, OUTPUTEVENTS, INPUTEVENTS)
# =====================================================================================
class KeyAcc:
    """Aggregates for one key (a channel or a single ITEMID).  Encounter-level arrays live
    in memory only and are reduced to counts before anything is written."""

    def __init__(self, n: int, cap: int, rng: np.random.Generator, rng_range: Optional[Tuple[float, float]],
                 collect_times: bool = False, sum24: bool = False):
        self.mask = np.zeros(n, dtype=np.int64)
        self.n24 = np.zeros(n, dtype=np.int32)
        self.pre_any = np.zeros(n, dtype=np.uint8)
        self.n_pre = 0
        self.n_win = 0
        self.units: Counter = Counter()
        self.n_num = self.n_nonnum = self.below = self.above = 0
        self.res = Reservoir(cap, rng)
        self.range = rng_range
        self.collect_times = collect_times
        self.times: List[Tuple[np.ndarray, np.ndarray]] = []
        self.sum24 = np.zeros(n, dtype=float) if sum24 else None

    def update(self, enc: np.ndarray, h: np.ndarray, val: np.ndarray, raw_blank: np.ndarray, uom: pd.Series) -> None:
        pre = h < 0
        win = ~pre
        self.n_pre += int(pre.sum())
        if pre.any():
            self.pre_any[enc[pre]] = 1
        if not win.any():
            return
        e, hh, v = enc[win], h[win], val[win]
        self.n_win += int(win.sum())
        b = np.floor(hh / CONFIG["BIN_HOURS"]).astype(np.int64)
        b = np.clip(b, 0, CONFIG["WINDOW_HOURS"] // CONFIG["BIN_HOURS"] - 1)
        np.bitwise_or.at(self.mask, e, np.left_shift(np.int64(1), b))
        in24 = hh < CONFIG["LAYER_B_HOURS"]
        np.add.at(self.n24, e[in24], 1)
        isnum = np.isfinite(v)
        self.n_num += int(isnum.sum())
        self.n_nonnum += int((~isnum & ~raw_blank[win]).sum())
        if self.range is not None:
            lo, hi = self.range
            self.below += int((isnum & (v < lo)).sum())
            self.above += int((isnum & (v > hi)).sum())
        self.res.add(v[isnum])
        self.units.update(uom[win].fillna("<none>").astype(str).str.strip().tolist())
        if self.collect_times:
            self.times.append((e.copy(), hh.copy()))
        if self.sum24 is not None:
            ok = in24 & isnum
            np.add.at(self.sum24, e[ok], v[ok])


class EventPass:
    def __init__(self, coh: pd.DataFrame, auto: Dict[str, List[Tuple[str, str]]]):
        self.coh = coh
        self.n = len(coh)
        self.rng = np.random.default_rng(CONFIG["RANDOM_SEED"])
        self.keys: Dict[str, KeyAcc] = {}
        self.keymap: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))  # src -> itemid -> keys
        self.channel_origin: Dict[str, str] = {}
        self.track_all_items = {"output", "input"}
        self.panel: Dict[str, Dict[str, Any]] = {}
        self.release: Dict[str, Dict[str, Any]] = {}
        chans = CONFIG["CHANNELS"]
        for ch, spec in chans.items():
            ids = [(s, i) for i in spec["itemids"] for s in (["chart", "lab"] if spec["source"] == "both" else [spec["source"]])]
            origin = "configured"
            if not ids and CONFIG["USE_DICTIONARY_CANDIDATES_FOR_EMPTY_CHANNELS"]:
                ids = auto.get(ch, [])
                origin = "auto (dictionary candidates, review)" if ids else "none"
            self.channel_origin[ch] = origin
            if not ids and spec["source"] not in ("input",):
                continue
            key = f"ch:{ch}"
            self.keys[key] = KeyAcc(self.n, CONFIG["SAMPLE_CAP"], self.rng, spec["range"], spec["collect_times"], spec["sum24"])
            for s, i in ids:
                self.keymap[s][str(i).upper()].append(key)
                ik = f"it:{s}:{str(i).upper()}"
                if ik not in self.keys:
                    self.keys[ik] = KeyAcc(self.n, CONFIG["ITEM_SAMPLE_CAP"], self.rng, spec["range"])
                self.keymap[s][str(i).upper()].append(ik)
        self.channels_by_src = defaultdict(list)
        for ch, spec in chans.items():
            if f"ch:{ch}" in self.keys:
                srcs = ["chart", "lab"] if spec["source"] == "both" else [spec["source"]]
                for s in srcs:
                    self.channels_by_src[s].append(ch)

    # ---------------------------------------------------------------
    def _ensure_item(self, src: str, itemid: str) -> None:
        ik = f"it:{src}:{itemid}"
        if ik not in self.keys:
            self.keys[ik] = KeyAcc(self.n, CONFIG["ITEM_SAMPLE_CAP"], self.rng, None)
            self.keymap[src][itemid].append(ik)

    def run_table(self, path: Path, src: str, extra_time_cols: Sequence[str] = ()) -> Dict[str, Any]:
        hdr = read_header(path)
        item_c = find_col(hdr, ["ITEMID"], ["item id"])
        time_c = find_col(hdr, ["CHARTTIME"], ["chart time"])
        vnum_c = find_col(hdr, ["VALUENUM"], ["value num"])
        val_c = find_col(hdr, ["VALUE", "AMOUNT"], ["value", "amount", "volume"], exclude=["num", "uom", "unit"])
        uom_c = find_col(hdr, ["VALUEUOM", "AMOUNTUOM"], ["value uom", "amount uom", "uom", "unit"])
        if not time_c:
            return {"error": f"{path.name}: no CHARTTIME column [bracket]"}
        wanted = [c for c in ["SUBJECT_ID", item_c, time_c, vnum_c, val_c, uom_c, *extra_time_cols] if c]
        info = {"file": path.name, "columns": hdr, "itemid_column": item_c, "value_column": vnum_c or val_c,
                "unit_column": uom_c, "rows_read": 0, "rows_tracked": 0, "rows_attributed": 0}
        U = {c: c.upper() for c in wanted}
        rel = {c.upper(): {"res": Reservoir(CONFIG["SAMPLE_CAP"], self.rng), "n": 0, "missing": 0, "neg": 0, "zero": 0}
               for c in extra_time_cols}
        for chunk in iter_chunks(path, wanted):
            info["rows_read"] += len(chunk)
            if item_c:
                chunk["_ITEM"] = norm_id(chunk[U[item_c]])
                if src in self.track_all_items:
                    for it in chunk["_ITEM"].dropna().unique():
                        self._ensure_item(src, str(it))
                chunk = chunk[chunk["_ITEM"].isin(list(self.keymap[src].keys()))]
            else:
                chunk["_ITEM"] = "__ALL__"
                self._ensure_item(src, "__ALL__")
                for ch in self.channels_by_src.get(src, []):
                    if f"ch:{ch}" not in self.keymap[src]["__ALL__"]:
                        self.keymap[src]["__ALL__"].append(f"ch:{ch}")
            if chunk.empty:
                continue
            info["rows_tracked"] += len(chunk)
            m = attribute(chunk, self.coh, U[time_c])
            if m.empty:
                continue
            info["rows_attributed"] += len(m)
            raw_v = m[U[val_c]] if val_c else pd.Series([None] * len(m), index=m.index)
            v = to_num(m[U[vnum_c]]) if vnum_c else to_num(raw_v)
            if vnum_c and val_c:
                v = v.fillna(to_num(raw_v))
            m["_V"] = v.to_numpy(dtype=float)
            m["_BLANK"] = blank(raw_v).to_numpy() if val_c else ~np.isfinite(m["_V"].to_numpy())
            m["_UOM"] = m[U[uom_c]] if uom_c else None
            # release / store timing (Section D), cohort rows only
            for c in extra_time_cols:
                C = c.upper()
                t2 = to_dt(m[C])
                lag = (t2 - m["_T"]).dt.total_seconds().to_numpy() / 3600.0
                r = rel[C]
                r["n"] += len(m)
                r["missing"] += int(np.isnan(lag).sum())
                r["neg"] += int((lag < 0).sum())
                r["zero"] += int((lag == 0).sum())
                r["res"].add(lag)
            # per key updates
            for itemid, g in m.groupby("_ITEM", sort=False):
                keys = self.keymap[src].get(str(itemid), [])
                if not keys:
                    continue
                enc = g["ENC_IDX"].to_numpy()
                h = g["H"].to_numpy(dtype=float)
                val = g["_V"].to_numpy(dtype=float)
                bl = g["_BLANK"].to_numpy(dtype=bool)
                uom = g["_UOM"] if uom_c else pd.Series([None] * len(g))
                for k in dict.fromkeys(keys):
                    self.keys[k].update(enc, h, val, bl, uom.reset_index(drop=True) if not uom_c else uom)
            self._panel_update(m, src)
        info["release_columns"] = {c: {"rows": sup(r["n"]), "missing": cnt_pct(r["missing"], r["n"]),
                                       "negative_lag": cnt_pct(r["neg"], r["n"]), "zero_lag": cnt_pct(r["zero"], r["n"]),
                                       "lag_hours": quantiles(r["res"].buf)} for c, r in rel.items()}
        return info

    def _panel_update(self, m: pd.DataFrame, src: str) -> None:
        """Co-recording of channels at one (subject, CHARTTIME) group, window rows only."""
        if src not in ("lab", "chart"):
            return
        chs = self.channels_by_src.get(src, [])
        if not chs:
            return
        idx = {c: i for i, c in enumerate(chs)}
        rows = []
        w = m[m["H"] >= 0]
        for itemid, g in w.groupby("_ITEM", sort=False):
            for k in self.keymap[src].get(str(itemid), []):
                if k.startswith("ch:") and k[3:] in idx:
                    rows.append(pd.DataFrame({"g": g["SUBJECT_ID"].astype(str) + "|" + g["_T"].astype(str), "c": idx[k[3:]]}))
        if not rows:
            return
        d = pd.concat(rows, ignore_index=True).drop_duplicates()
        gid, _ = pd.factorize(d["g"])
        G, Cn = int(gid.max()) + 1, len(chs)
        M = np.zeros((G, Cn), dtype=np.float32)
        M[gid, d["c"].to_numpy()] = 1.0
        co = (M.T @ M).astype(np.int64)
        size = M.sum(axis=1)
        st = self.panel.setdefault(src, {"channels": chs, "co": np.zeros((Cn, Cn), dtype=np.int64),
                                         "in_group_ge2": np.zeros(Cn, dtype=np.int64),
                                         "in_group_ge5": np.zeros(Cn, dtype=np.int64),
                                         "groups": 0, "group_size": Counter()})
        st["co"] += co
        st["in_group_ge2"] += (M[size >= 2]).sum(axis=0).astype(np.int64)
        st["in_group_ge5"] += (M[size >= 5]).sum(axis=0).astype(np.int64)
        st["groups"] += G
        st["group_size"].update(Counter(np.minimum(size, 20).astype(int).tolist()))

    # ---------------------------------------------------------------
    def summarize_key(self, key: str, full: bool) -> Dict[str, Any]:
        a = self.keys[key]
        coh = self.coh
        n = self.n
        tm = train_mask(coh)
        nb = CONFIG["WINDOW_HOURS"] // CONFIG["BIN_HOURS"]
        has24 = a.n24 > 0
        long24 = ((coh["WIN_END"] - coh["INTIME"]).dt.total_seconds() / 3600.0 >= CONFIG["LAYER_B_HOURS"]).to_numpy()
        ge3 = a.n24 >= CONFIG["MIN_VALUES_VARIABILITY"]
        out: Dict[str, Any] = {
            "records_in_window_0_168h": sup(a.n_win), "records_pre_INTIME_24h": sup(a.n_pre),
            "coverage_0_24h_all": cnt_pct(has24.sum(), n),
            "coverage_0_24h_train": cnt_pct((has24 & tm).sum(), tm.sum()),
            "coverage_0_24h_if_window_ge_24h": cnt_pct((has24 & long24).sum(), long24.sum()),
            "ge3_values_0_24h_all": cnt_pct(ge3.sum(), n),
            "ge3_values_0_24h_among_covered": cnt_pct(ge3.sum(), has24.sum()),
            "values_per_covered_encounter_0_24h": quantiles(a.n24[has24]),
            "units": sup_counts(a.units, 8),
            "numeric_records": sup(a.n_num), "non_numeric_nonblank_records": sup(a.n_nonnum),
        }
        if a.range is not None:
            out["plausibility"] = {"range_native": list(a.range), "below": cnt_pct(a.below, a.n_num), "above": cnt_pct(a.above, a.n_num)}
        out["value_distribution"] = quantiles(a.res.buf)
        if full:
            B = bits_to_bins(a.mask, nb)
            openb = np.arange(nb)[None, :] < coh["N_OPEN_BINS"].to_numpy()[:, None]
            out["coverage_per_6h_bin"] = [pct(int((B[:, j] & openb[:, j]).sum()), int(openb[:, j].sum())) for j in range(nb)]
            out["coverage_per_6h_bin_denominators"] = [sup(int(openb[:, j].sum())) for j in range(nb)]
        return out

    def summarize(self) -> Dict[str, Any]:
        R: Dict[str, Any] = {"resolves": RESOLVES["C_coverage"] + RESOLVES["C_variability"] + RESOLVES["C_units_plausibility"],
                             "channels": {}, "items": {}}
        for ch, spec in CONFIG["CHANNELS"].items():
            key = f"ch:{ch}"
            entry = {"system": spec["system"], "triggered": spec["triggered"], "origin": self.channel_origin.get(ch),
                     "range_note": spec["range_note"]}
            if key in self.keys:
                entry.update(self.summarize_key(key, full=True))
                entry["itemids_used"] = sorted({k.split(":", 2)[2] for s in self.keymap for i, ks in self.keymap[s].items()
                                                for k in ks if key in ks and k.startswith("it:")})
            else:
                entry["status"] = "no ITEMIDs configured or found"
            R["channels"][ch] = entry
        for key in self.keys:
            if key.startswith("it:"):
                s = self.summarize_key(key, full=False)
                if self.keys[key].n_win or self.keys[key].n_pre:
                    R["items"][key[3:]] = s
        # system-level: share of encounters with any channel of the system in 0-24 h
        sysenc: Dict[str, np.ndarray] = {}
        for ch, spec in CONFIG["CHANNELS"].items():
            if spec["system"] == "support" or f"ch:{ch}" not in self.keys:
                continue
            a = self.keys[f"ch:{ch}"]
            sysenc[spec["system"]] = sysenc.get(spec["system"], np.zeros(self.n, bool)) | (a.n24 > 0)
        R["systems_any_channel_0_24h"] = {s: cnt_pct(v.sum(), self.n) for s, v in sysenc.items()}
        return R

    def panels(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {"resolves": RESOLVES["C_panels"],
                               "note": "groups = distinct (subject, CHARTTIME) within 0-168 h; groups split across read chunks are counted twice (minor)"}
        for src, st in self.panel.items():
            chs = st["channels"]
            diag = np.diag(st["co"])
            per = {}
            for i, c in enumerate(chs):
                if diag[i] == 0:
                    continue
                partners = sorted(((chs[j], st["co"][i, j]) for j in range(len(chs)) if j != i and st["co"][i, j] > 0),
                                  key=lambda kv: -kv[1])[:6]
                per[c] = {"group_records": sup(int(diag[i])),
                          "share_in_groups_ge2_channels": pct(int(st["in_group_ge2"][i]), int(diag[i])),
                          "share_in_groups_ge5_channels": pct(int(st["in_group_ge5"][i]), int(diag[i])),
                          "top_co_recorded": {p: pct(int(v), int(diag[i])) for p, v in partners}}
            out[src] = {"groups": sup(st["groups"]), "group_size_distribution": sup_counts(st["group_size"]), "channels": per}
        return out

    def hr_gaps(self) -> Dict[str, Any]:
        a = self.keys.get("ch:hr")
        if a is None or not a.times:
            return {"resolves": RESOLVES["C_hr_gaps"], "status": "no HR records"}
        e = np.concatenate([t[0] for t in a.times])
        h = np.concatenate([t[1] for t in a.times])
        gaps = gap_hours(e, h)
        out = {"resolves": RESOLVES["C_hr_gaps"], "gap_hours_0_168h": quantiles(gaps)}
        out["gap_share_over"] = {f">{x}h": cnt_pct(int((gaps > x).sum()), gaps.size) for x in (1, 2, 4, 6, 12, 24)}
        in24 = h < CONFIG["LAYER_B_HOURS"]
        g24 = gap_hours(e[in24], h[in24])
        out["gap_hours_0_24h"] = quantiles(g24)
        # largest within-24h gap per encounter including window edges
        df = pd.DataFrame({"e": e[in24], "h": h[in24]})
        if not df.empty:
            st = df.groupby("e")["h"].agg(["min", "max", "count"])
            out["first_HR_hours_after_INTIME"] = quantiles(st["min"].to_numpy())
            out["encounters_with_HR_0_24h"] = cnt_pct(len(st), self.n)
        out["records_per_hour_of_window"] = round(len(h) / max(1.0, float(((self.coh["WIN_END"] - self.coh["INTIME"]).dt.total_seconds() / 3600).sum())), 4)
        return out

    def weight_and_fluids(self) -> Dict[str, Any]:
        n = self.n
        W: Dict[str, Any] = {"resolves": RESOLVES["C_weight"]}
        wa = self.keys.get("ch:weight")
        uo = self.keys.get("ch:urine_output")
        fi = self.keys.get("ch:fluid_intake")
        if wa is not None:
            w24 = (wa.n24 > 0) | (wa.pre_any > 0)
            wany = (wa.mask > 0) | (wa.pre_any > 0)
            W.update({"weight_in_-24h_to_24h": cnt_pct(w24.sum(), n), "weight_any_-24h_to_168h": cnt_pct(wany.sum(), n),
                      "weight_values": quantiles(wa.res.buf), "weight_units": sup_counts(wa.units, 6)})
            if uo is not None:
                u2 = uo.n24 >= 2
                W["weight_and_ge2_urine_records_0_24h (mL/kg/h feasible)"] = cnt_pct((w24 & u2).sum(), n)
                if fi is not None:
                    W["weight_and_intake_and_output_0_24h (%FO feasible)"] = cnt_pct((w24 & (uo.n24 > 0) & (fi.n24 > 0)).sum(), n)
        else:
            W["status"] = "no weight item found (see Section B)"
        F: Dict[str, Any] = {"resolves": RESOLVES["C_fluids"]}
        for name, a in (("urine_output", uo), ("fluid_intake", fi)):
            if a is None:
                F[name] = {"status": "no items found"}
                continue
            d: Dict[str, Any] = {"encounters_with_records_0_24h": cnt_pct((a.n24 > 0).sum(), n),
                                 "encounters_with_ge2_records_0_24h": cnt_pct((a.n24 >= 2).sum(), n),
                                 "records_per_encounter_0_24h": quantiles(a.n24[a.n24 > 0]),
                                 "record_values": quantiles(a.res.buf), "units": sup_counts(a.units, 6)}
            if a.times:
                e = np.concatenate([t[0] for t in a.times])
                h = np.concatenate([t[1] for t in a.times])
                g = gap_hours(e[h < 24], h[h < 24])
                d["gap_between_records_hours_0_24h"] = quantiles(g)
                d["record_minute_of_hour_on_the_hour"] = pct(int(np.isclose(np.mod(h * 60, 60), 0, atol=1e-6).sum()), len(h))
            if a.sum24 is not None:
                tot = a.sum24[a.n24 > 0]
                d["total_0_24h_per_encounter"] = quantiles(tot)
            F[name] = d
        return {"C_weight": W, "C_fluids": F}


# =====================================================================================
# Section D (microbiology part) and Section E (drugs)
# =====================================================================================
def _kw_regex(en: Sequence[str], cn: Sequence[str]) -> Optional[re.Pattern]:
    parts = [EN_BOUNDARY.format(re.escape(k)) for k in en if k] + [re.escape(k) for k in cn if k]
    return re.compile("|".join(parts), flags=re.IGNORECASE) if parts else None


def extra_time_columns(hdr: Sequence[str], base: str = "CHARTTIME") -> List[str]:
    """Time-like columns other than CHARTTIME (result / report / store / verify / receive)."""
    roles = ["time", "date", "dt", "datetime", "timestamp", "storetime", "reporttime", "resulttime",
             "verifiedtime", "receivetime", "releasetime", "signtime", "approvetime"]
    return [c for c in hdr if c.upper() != base.upper() and role_match(c, roles, exclude=["id", "uom", "unit"])]


def section_D_micro(raw: Path, coh: pd.DataFrame) -> Dict[str, Any]:
    p = table_path(raw, CONFIG["TABLES"]["MICROBIOLOGYEVENTS"])
    if p is None:
        return {"status": "MICROBIOLOGYEVENTS not found"}
    hdr = read_header(p)
    time_c = find_col(hdr, ["CHARTTIME"], ["chart time"])
    if not time_c:
        return {"status": "MICROBIOLOGYEVENTS has no CHARTTIME column [bracket]", "columns": hdr}
    extra = extra_time_columns(hdr)
    text_cols = [c for c in hdr if role_match(c, ["spec type desc", "org name", "test name", "name", "desc", "description",
                                                   "item", "label"], exclude=["id", "itemid", "ab"])]
    org_c = find_col(hdr, ["ORG_NAME"], ["org name", "organism"])
    wanted = list(dict.fromkeys(["SUBJECT_ID", time_c, *extra, *text_cols] + ([org_c] if org_c else [])))
    fam_re = {f: _kw_regex(k["en"], k["cn"]) for f, k in CONFIG["PATHOGEN_FAMILIES"].items()}
    rng = np.random.default_rng(CONFIG["RANDOM_SEED"] + 1)
    rel = {c.upper(): {"res": Reservoir(CONFIG["SAMPLE_CAP"], rng), "n": 0, "missing": 0, "neg": 0, "zero": 0} for c in extra}
    n = len(coh)
    coll24 = np.zeros(n, bool)
    fam_enc24 = {f: np.zeros(n, bool) for f in fam_re}
    fam_enc = {f: np.zeros(n, bool) for f in fam_re}
    fam_rows = Counter()
    fam_vals: Dict[str, Counter] = {f: Counter() for f in fam_re}
    rows = org_rows = 0
    for chunk in iter_chunks(p, [c for c in wanted if c]):
        m = attribute(chunk, coh, time_c.upper(), pre_hours=0)
        if m.empty:
            continue
        rows += len(m)
        e = m["ENC_IDX"].to_numpy()
        h = m["H"].to_numpy()
        coll24[e[h < 24]] = True
        if org_c:
            org_rows += int((~blank(m[org_c.upper()])).sum())
        for c in extra:
            C = c.upper()
            lag = (to_dt(m[C]) - m["_T"]).dt.total_seconds().to_numpy() / 3600.0
            r = rel[C]
            r["n"] += len(m)
            r["missing"] += int(np.isnan(lag).sum())
            r["neg"] += int((lag < 0).sum())
            r["zero"] += int((lag == 0).sum())
            r["res"].add(lag)
        txt = m[[c.upper() for c in text_cols]].astype("string").fillna("").agg(" | ".join, axis=1) if text_cols else pd.Series([""] * len(m))
        for f, rx in fam_re.items():
            if rx is None:
                continue
            hit = txt.str.contains(rx, na=False).to_numpy()
            if hit.any():
                fam_rows[f] += int(hit.sum())
                fam_enc[f][e[hit]] = True
                fam_enc24[f][e[hit & (h < 24)]] = True
                fam_vals[f].update(txt[hit].str.slice(0, 80).tolist())
    return {"resolves": RESOLVES["D_release"], "file": p.name, "columns": hdr,
            "time_columns_beyond_CHARTTIME": extra,
            "release_columns": {c: {"rows": sup(r["n"]), "missing": cnt_pct(r["missing"], r["n"]),
                                    "negative_lag": cnt_pct(r["neg"], r["n"]), "zero_lag": cnt_pct(r["zero"], r["n"]),
                                    "lag_hours": quantiles(r["res"].buf)} for c, r in rel.items()},
            "rows_in_cohort_window": sup(rows),
            "rows_with_organism_name": cnt_pct(org_rows, rows),
            "encounters_with_collection_0_24h": cnt_pct(coll24.sum(), n),
            "pathogen_families": {f: {"rows": sup(fam_rows[f]), "encounters_0_24h": cnt_pct(fam_enc24[f].sum(), n),
                                      "encounters_0_168h": cnt_pct(fam_enc[f].sum(), n),
                                      "top_matched_text": sup_counts(fam_vals[f], 15)} for f in fam_re},
            "note": "no extra time column means release = CHARTTIME + assumed turnaround (S9 / M3 grid)" if not extra else ""}


def section_E(raw: Path, coh: pd.DataFrame, headers: Dict[str, List[str]]) -> Dict[str, Any]:
    p = table_path(raw, CONFIG["TABLES"]["PRESCRIPTIONS"])
    if p is None:
        return {"status": "PRESCRIPTIONS not found"}
    hdr = read_header(p)
    st_c = find_col(hdr, ["STARTDATE", "STARTTIME"], ["start date", "start time"])
    en_c = find_col(hdr, ["ENDDATE", "ENDTIME"], ["end date", "end time"])
    name_cols = [c for c in hdr if role_match(c, ["drug name", "drug name en", "drug name generic", "drug", "generic name"],
                                              exclude=["id", "type", "form", "strength", "unit"])]
    dose_cols = [c for c in hdr if c.upper() in {d.upper() for d in CONFIG["DOSE_COLUMN_CANDIDATES"]}
                 or role_match(c, ["dose", "dosage", "rate", "strength", "speed"], exclude=["id"])]
    wanted = list(dict.fromkeys(["SUBJECT_ID", st_c, en_c, *name_cols, *dose_cols]))
    groups = CONFIG["DRUG_GROUPS"]
    grx = {g: _kw_regex(v["en"], v["cn"]) for g, v in groups.items()}
    exrx = {g: (re.compile("|".join(re.escape(x) for x in v["exclude"]), re.I) if v["exclude"] else None) for g, v in groups.items()}
    n = len(coh)
    tot_rows = 0
    end_missing_all = 0
    start_date_only = 0
    G = {g: {"rows_all": 0, "rows_window": 0, "end_missing": 0, "enc24": np.zeros(n, bool), "enc": np.zeros(n, bool),
             "names": Counter(), "agents": Counter(), "dose_nonblank": Counter(), "dose_units": Counter()} for g in groups}
    overlap = Counter()
    for chunk in iter_chunks(p, [c for c in wanted if c]):
        tot_rows += len(chunk)
        if en_c:
            end_missing_all += int(blank(chunk[en_c.upper()]).sum())
        if st_c:
            sv = chunk[st_c.upper()].astype("string").str.strip()
            start_date_only += int((sv.str.len() <= 10).sum() + sv.str.endswith("00:00:00").sum())
        txt = chunk[[c.upper() for c in name_cols]].astype("string").fillna("").agg(" | ".join, axis=1) if name_cols else pd.Series([""] * len(chunk), index=chunk.index)
        hits = {}
        for g, rx in grx.items():
            h = txt.str.contains(rx, na=False) if rx is not None else pd.Series(False, index=chunk.index)
            if exrx[g] is not None:
                h &= ~txt.str.contains(exrx[g], na=False)
            hits[g] = h
        combo = pd.DataFrame(hits)
        for k, v in combo[combo.sum(axis=1) > 1].astype(int).astype(str).agg("".join, axis=1).value_counts().items():
            overlap["+".join(g for g, b in zip(combo.columns, k) if b == "1")] += int(v)
        for g, h in hits.items():
            if not h.any():
                continue
            sub = chunk[h.to_numpy()]
            st = G[g]
            st["rows_all"] += len(sub)
            if en_c:
                st["end_missing"] += int(blank(sub[en_c.upper()]).sum())
            for nc in name_cols[:2]:
                st["names"].update(sub[nc.upper()].dropna().astype(str).str.strip().tolist())
            for t in txt[h].tolist():
                for kw in keyword_hits(t, groups[g]["en"] + groups[g]["cn"]):
                    st["agents"][kw.lower()] += 1
            if st_c:
                m = attribute(sub, coh, st_c.upper())
                m = m[m["H"] >= 0]
                if not m.empty:
                    st["rows_window"] += len(m)
                    e, hh = m["ENC_IDX"].to_numpy(), m["H"].to_numpy()
                    st["enc"][e] = True
                    st["enc24"][e[hh < 24]] = True
                    for dc in dose_cols:
                        D = dc.upper()
                        st["dose_nonblank"][dc] += int((~blank(m[D])).sum())
                        if role_match(dc, ["unit", "uom"]):
                            st["dose_units"].update(m[D].dropna().astype(str).str.strip().tolist())
    dose_like_elsewhere = {f: [c for c in cols if role_match(c, ["dose", "dosage", "rate", "speed", "amount", "ug kg min"], exclude=["id"])]
                           for f, cols in headers.items() if not f.upper().startswith("PRESCRIPTIONS")}
    out = {"resolves": RESOLVES["E_drugs"], "file": p.name, "name_columns": name_cols, "dose_columns": dose_cols,
           "start_column": st_c, "end_column": en_c, "rows": sup(tot_rows),
           "ENDDATE_missing_all_rows": cnt_pct(end_missing_all, tot_rows),
           "STARTDATE_date_only_or_midnight": cnt_pct(start_date_only, tot_rows),
           "group_overlaps_rows": sup_counts(overlap),
           "dose_or_rate_columns_in_other_tables": {k: v for k, v in dose_like_elsewhere.items() if v},
           "groups": {}}
    for g, st in G.items():
        out["groups"][g] = {
            "rows_all": sup(st["rows_all"]), "rows_in_cohort_window": sup(st["rows_window"]),
            "ENDDATE_missing": cnt_pct(st["end_missing"], st["rows_all"]),
            "encounters_0_24h": cnt_pct(st["enc24"].sum(), n), "encounters_0_168h": cnt_pct(st["enc"].sum(), n),
            "agents_matched": sup_counts(st["agents"]),
            "top_drug_names": sup_counts(st["names"], CONFIG["TOP_N_NAMES"]),
            "dose_field_nonblank_in_window": {c: pct(v, st["rows_window"]) for c, v in st["dose_nonblank"].items()},
            "dose_units_in_window": sup_counts(st["dose_units"], 10),
        }
    out["note"] = ("pSOFA / vasoactive-inotrope dosing needs a rate (e.g. mcg/kg/min).  If only DOSE_VAL_RX / "
                   "PROD_STRENGTH exist, the vasoactive component cannot be scored as published -> 'modified pSOFA'.")
    return out


# =====================================================================================
# Report writing
# =====================================================================================
def _json_default(o: Any) -> Any:
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if isinstance(o, (pd.Timestamp,)):
        return str(o)
    if isinstance(o, set):
        return sorted(o)
    return str(o)


def _md_value(v: Any) -> str:
    if isinstance(v, dict) and set(v.keys()) == {"n", "pct"}:
        return f"{v['n']} ({v['pct']}%)" if isinstance(v["pct"], (int, float)) else f"{v['n']} ({v['pct']})"
    if isinstance(v, (list, tuple)) and len(v) > 12:
        return ", ".join(str(x) for x in v[:12]) + f", ... ({len(v)} total)"
    if isinstance(v, dict) and len(v) <= 8 and all(not isinstance(x, (dict, list)) for x in v.values()):
        return "; ".join(f"{k}: {x}" for k, x in v.items())
    return json.dumps(v, ensure_ascii=False, default=_json_default) if isinstance(v, (dict, list)) else str(v)


def _is_cp(v: Any) -> bool:
    return isinstance(v, dict) and set(v.keys()) == {"n", "pct"}


def _md_dict(d: Dict[str, Any], lines: List[str], depth: int = 0, max_depth: int = 2) -> None:
    pad = "  " * depth
    for k, v in d.items():
        if k in ("resolves", "columns", "candidates", "coverage_per_6h_bin_denominators"):
            continue
        nested = isinstance(v, dict) and not _is_cp(v) and (
            len(v) > 8 or any(isinstance(x, dict) and not _is_cp(x) for x in v.values()))
        if nested and depth < max_depth:
            lines.append(f"{pad}- **{k}**")
            _md_dict(v, lines, depth + 1, max_depth)
        else:
            lines.append(f"{pad}- {k}: {_md_value(v)}")


def write_markdown(report: Dict[str, Any], path: Path) -> None:
    L: List[str] = []
    meta = report["meta"]
    L += ["# PIC dataset audit for Foundations Plan v12 (aggregates only)", "",
          f"Run: {meta['run_local']} (local) / {meta['run_utc']} UTC. Script {meta['script']} v{meta['script_version']}.",
          f"Small cells: counts below {CONFIG['SMALL_CELL']} shown as '<{CONFIG['SMALL_CELL']}'; quantiles omit min/max and need n >= {CONFIG['MIN_N_QUANTILE']}.",
          "Review before sharing: complementary suppression is not applied.", ""]
    if report.get("notes"):
        L += ["## Notes and brackets raised during the run", ""] + [f"- {x}" for x in report["notes"]] + [""]
    L += ["## Cohort", ""]
    _md_dict(report.get("cohort", {}), L)
    L.append("")
    for sec, title in [("A", "Section A: Layer A"), ("B", "Section B: dictionary search"), ("C", "Section C: Layers B/C"),
                       ("D", "Section D: release timing"), ("E", "Section E: drugs")]:
        if sec not in report:
            continue
        L += [f"## {title}", ""]
        S = report[sec]
        if sec == "B":
            L += ["Resolves: " + "; ".join(RESOLVES["B_dictionary"]), "",
                  "| Channel | System | Trig. | Status | Configured | Top candidates (ITEMID: label) |", "|---|---|---|---|---|---|"]
            for ch, e in S.get("channels", {}).items():
                top = "; ".join(f"{c['itemid']}: {c['text'][:40]}" for c in e.get("candidates", [])[:3]).replace("|", "/")
                L.append(f"| {ch} | {e['system']} | {'Y' if e['triggered'] else ''} | {e['status']} | {', '.join(e['configured_itemids'])} | {top} |")
            L.append("")
            L += ["### Serum vs urine sodium candidates", ""]
            for c in S.get("sodium_serum_vs_urine", []):
                L.append(f"- {c.get('table')} {c.get('itemid')}: {str(c.get('text'))[:80]} (fluid: {c.get('fluid')}; {c.get('fluid_flag')})")
            L += ["", "### Pathogen families in dictionaries", ""]
            for f, e in S.get("pathogen_families_in_dictionaries", {}).items():
                L.append(f"- {f}: {e['n_candidates']} candidates; " + "; ".join(f"{c['itemid']}: {c['text'][:40]}" for c in e["candidates"][:5]))
            L.append("")
            continue
        if sec == "C":
            ch = S.get("C_coverage", {}).get("channels", {})
            L += ["Resolves: " + "; ".join(RESOLVES["C_coverage"] + RESOLVES["C_variability"]), "",
                  "| Channel | System | Origin | Cov. 0-24h | Cov. 0-24h (train) | >=3 values 24h | Bin 0 | Bin 4 | Bin 12 | Main unit | Below / above range |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
            for name, e in ch.items():
                if "coverage_0_24h_all" not in e:
                    L.append(f"| {name} | {e['system']} | {e.get('origin')} | not found | | | | | | | |")
                    continue
                bins = e.get("coverage_per_6h_bin", [None] * 28)
                unit = next(iter(e.get("units", {}) or {"": ""}))
                pl = e.get("plausibility", {})
                L.append(f"| {name} | {e['system']} | {e.get('origin')} | {_md_value(e['coverage_0_24h_all'])} | "
                         f"{_md_value(e['coverage_0_24h_train'])} | {_md_value(e['ge3_values_0_24h_all'])} | {bins[0]} | {bins[4]} | {bins[12]} | "
                         f"{unit} | {pl.get('below', {}).get('pct', '')} / {pl.get('above', {}).get('pct', '')} |")
            L.append("")
            for sub in ("systems_any_channel_0_24h",):
                L += [f"### {sub}", ""] + [f"- {k}: {_md_value(v)}" for k, v in S.get("C_coverage", {}).get(sub, {}).items()] + [""]
            for sub in ("C_panels", "C_hr_gaps", "C_weight", "C_fluids", "D_lab_chart_release"):
                if sub in S:
                    L += [f"### {sub}", "", "Resolves: " + "; ".join(S[sub].get("resolves", [])), ""]
                    _md_dict(S[sub], L)
                    L.append("")
            continue
        for sub, body in S.items():
            if isinstance(body, dict):
                L += [f"### {sub}", ""]
                if body.get("resolves"):
                    L += ["Resolves: " + "; ".join(body["resolves"]), ""]
                _md_dict(body, L)
                L.append("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def atomic_write(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


# =====================================================================================
# Main
# =====================================================================================
SCRIPT_VERSION = "1.0 (v12 audit, 2026-11-09)"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Aggregate-only PIC audit for Foundations Plan v12")
    ap.add_argument("--sections", default="A,B,C,D,E", help="comma list of A,B,C,D,E")
    ap.add_argument("--project-root", default=None)
    ap.add_argument("--raw-root", default=None)
    ap.add_argument("--chunksize", type=int, default=None)
    ap.add_argument("--self-test", action="store_true", help="run synthetic self-tests only (no data read)")
    args = ap.parse_args(argv)
    if args.self_test:
        return run_self_tests()
    if args.project_root:
        CONFIG["PROJECT_ROOT"] = args.project_root
    if args.raw_root:
        CONFIG["RAW_ROOT_OVERRIDE"] = args.raw_root
    if args.chunksize:
        CONFIG["CHUNKSIZE"] = args.chunksize
    secs = {x.strip().upper() for x in args.sections.split(",") if x.strip()}

    project = find_project_root()
    raw, how = find_raw_root(project)
    if raw is None:
        log("PIC raw root not found; set PIC_ROOT, --raw-root, or CONFIG['RAW_ROOT_OVERRIDE']")
        return 2
    log(f"project root: {project}")
    log(f"raw root: {raw} ({how})")
    notes: List[str] = []
    headers = header_scan(raw)
    report: Dict[str, Any] = {
        "meta": {"script": Path(__file__).name, "script_version": SCRIPT_VERSION,
                 "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                 "run_local": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
                 "raw_root_found_via": how, "sections": sorted(secs),
                 "spec": "revised_plan_v12.md; v12_change_list.md",
                 "privacy": {"small_cell": CONFIG["SMALL_CELL"], "min_n_quantile": CONFIG["MIN_N_QUANTILE"],
                             "patient_level_rows_written": False},
                 "config_snapshot": {k: v for k, v in CONFIG.items() if k not in ("CHANNELS", "DRUG_GROUPS")},
                 "channel_config": {k: {kk: vv for kk, vv in v.items() if kk in ("system", "source", "itemids", "range", "range_note", "triggered")}
                                    for k, v in CONFIG["CHANNELS"].items()}},
        "table_headers": headers,
        "notes": notes,
    }
    missing = [t for t, stem in CONFIG["TABLES"].items() if table_path(raw, stem) is None]
    if missing:
        notes.append(f"tables not found in raw root: {missing} [bracket: file names]")

    coh, cinfo = load_cohort(project, raw, notes)
    report["cohort"] = {"resolves": RESOLVES["A0_cohort"], **cinfo}
    log(f"cohort: {len(coh):,} encounters")

    if "A" in secs:
        log("Section A ...")
        A, coh = section_A(raw, coh, headers, notes)
        report["A"] = A
    auto: Dict[str, List[Tuple[str, str]]] = {}
    if "B" in secs or "C" in secs:
        log("Section B ...")
        B, auto, _ = section_B(raw)
        if "B" in secs:
            report["B"] = B
    if "C" in secs or "D" in secs:
        log("Sections C/D event pass ...")
        ev = EventPass(coh, auto)
        C: Dict[str, Any] = {}
        D_lab: Dict[str, Any] = {"resolves": RESOLVES["D_release"]}
        T = CONFIG["TABLES"]
        for tab, src in [("CHARTEVENTS", "chart"), ("LABEVENTS", "lab"), ("OUTPUTEVENTS", "output"), ("INPUTEVENTS", "input")]:
            p = table_path(raw, T[tab])
            if p is None:
                continue
            extra = extra_time_columns(read_header(p))
            log(f"  {tab}: extra time columns {extra}")
            info = ev.run_table(p, src, extra)
            D_lab[tab] = {"columns": info.get("columns"), "time_columns_beyond_CHARTTIME": extra,
                          "release_columns": info.get("release_columns", {}),
                          "rows_read": sup(info.get("rows_read", 0)), "rows_tracked": sup(info.get("rows_tracked", 0)),
                          "rows_attributed_to_cohort": sup(info.get("rows_attributed", 0))}
        if "C" in secs:
            C["C_coverage"] = ev.summarize()
            C["C_panels"] = ev.panels()
            C["C_hr_gaps"] = ev.hr_gaps()
            C.update(ev.weight_and_fluids())
            report["C"] = C
        if "D" in secs:
            log("Section D microbiology ...")
            report["D"] = {"D_lab_chart_release": D_lab, "D_micro": section_D_micro(raw, coh)}
    elif "D" in secs:
        report["D"] = {"D_micro": section_D_micro(raw, coh)}
    if "E" in secs:
        log("Section E ...")
        report["E"] = {"E_drugs": section_E(raw, coh, headers)}

    out_dir = project / CONFIG["OUT_DIR"]
    out_dir.mkdir(parents=True, exist_ok=True)
    jpath = out_dir / f"{CONFIG['OUT_STEM']}.json"
    mpath = out_dir / f"{CONFIG['OUT_STEM']}.md"
    atomic_write(jpath, json.dumps(report, ensure_ascii=False, indent=2, default=_json_default))
    tmp = mpath.with_suffix(".md.tmp")
    write_markdown(report, tmp)
    os.replace(tmp, mpath)
    log(f"wrote {jpath}")
    log(f"wrote {mpath}")
    return 0


# =====================================================================================
# Self-tests (synthetic data only; run with --self-test or pytest)
# =====================================================================================
def _toy_cohort() -> pd.DataFrame:
    coh = pd.DataFrame({
        "SUBJECT_ID": ["1", "1", "2"],
        "INTIME": pd.to_datetime(["2100-01-01 10:00:00", "2100-02-01 08:00:00", "2100-01-05 23:30:00"]),
        "OUTTIME": pd.to_datetime(["2100-01-02 22:00:00", "2100-02-20 08:00:00", "2100-01-06 01:00:00"]),
        "SPLIT": ["train", "train", "test"], "HADM_LIST": [["10"], ["11"], ["20"]],
    })
    coh["ENC_IDX"] = np.arange(3)
    coh["WIN_END"] = pd.concat([coh["OUTTIME"], coh["INTIME"] + pd.Timedelta(hours=168)], axis=1).min(axis=1)
    coh["N_OPEN_BINS"] = n_open_bins(coh["INTIME"], coh["WIN_END"])
    return coh


def test_role_match_word_boundaries():
    assert role_match("AGE", ["age"]) and role_match("age_years_at_intime", ["age"])
    assert not role_match("STAGE", ["age"]) and not role_match("DOSAGE", ["age"]) and not role_match("aki_stage", ["age"])
    assert role_match("event_dt", ["dt"]) and not role_match("WIDTH", ["dt"]) and not role_match("ADT_FLAG", ["dt"])
    assert role_match("icu_los_hours", ["los"]) and not role_match("CLOSE_TIME", ["los"]) and not role_match("closed_bin", ["los"])
    assert role_match("SURGERY_BEGIN_TIME", ["begin time"]) and role_match("surgeryBeginTime", ["begin time"])
    assert role_match("EDREGTIME", ["edregtime"]) and not role_match("EDOUTTIME", ["edregtime"])
    assert not role_match("ICUSTAY_ID", ["time"], exclude=["id"])


def test_find_col_case_insensitive():
    cols = ["subject_id", "ItemID", "charttime", "valuenum"]
    assert find_col(cols, ["ITEMID"]) == "ItemID"
    assert find_col(cols, ["CHARTTIME"]) == "charttime"
    assert find_col(["DOSAGE", "AGE_YRS"], [], ["age"]) == "AGE_YRS"
    assert resolve_cols(["Subject_Id", "VALUE"], ["SUBJECT_ID", "VALUE", "X"]) == {"SUBJECT_ID": "Subject_Id", "VALUE": "VALUE"}


def test_suppression_and_quantiles():
    assert sup(0) == 0 and sup(3) == "<10" and sup(9) == "<10" and sup(10) == 10
    assert pct(5, 100) == "suppressed" and pct(50, 100) == 50.0 and pct(0, 100) == 0.0
    q = quantiles(np.arange(5))
    assert q["quantiles"] == "suppressed (n too small)"
    q = quantiles(np.arange(100, dtype=float))
    assert "p50" in q["quantiles"] and "min" not in q["quantiles"] and "max" not in q["quantiles"]


def test_keyword_hits_boundaries_and_cjk():
    assert keyword_hits("Phosphate (serum)", ["pH"]) == []
    assert keyword_hits("Blood gas pH", ["pH"]) == ["pH"]
    assert keyword_hits("Na+ (serum)", ["Na"]) == ["Na"]
    assert keyword_hits("钠(血清)", ["钠"]) == ["钠"]
    assert keyword_hits("norepinephrine 4mg", ["epinephrine"]) == []
    assert keyword_hits("Epinephrine 1mg", ["epinephrine"]) == ["epinephrine"]
    assert keyword_hits(float("nan"), ["x"]) == []


def test_norm_id_and_reservoir():
    s = norm_id(pd.Series(["5227.0", " sv1 ", "1002"]))
    assert list(s) == ["5227", "SV1", "1002"]
    r = Reservoir(50, np.random.default_rng(0))
    for _ in range(10):
        r.add(np.arange(100, dtype=float))
    assert r.buf.size <= 50 and r.seen == 1000


def test_icd_chapters_and_domains():
    assert icd3("j18.900") == "J18" and icd3("12") is None and icd3("1AB") is None
    assert who_chapter("D48") == "II" and who_chapter("D50") == "III" and who_chapter("D49") is None
    assert who_chapter("H60") == "VIII" and who_chapter("T07") == "XIX"
    assert domain_dx_only("Q21") == "cardiac_congenital_heart" and domain_dx_only("Q65") == "other_medical"
    assert domain_dx_only("D61") == "oncology_hematology" and domain_dx_only("D10") == CONFIG["DOMAIN_D00_D49"]
    assert domain_dx_only("I50") == CONFIG["DOMAIN_CHAPTER_I"] and domain_dx_only(None) == "other_medical"
    assert domain_contract("J18", "CICU") == "cardiac_congenital_heart" and domain_contract("I50", "PICU") == "other_medical"


def test_bands_timing_and_admission_source():
    assert interval_band(23.9) == "<24h" and interval_band(24) == "24-<72h" and interval_band(168) == ">=168h"
    assert interval_band(float("nan")) is None
    assert shift_of(8) == "day_08-15" and shift_of(15) == "evening_15-23" and shift_of(23) == "night_23-08" and shift_of(7) == "night_23-08"
    assert timing_l2(10, 5) == "weekend_daytime" and timing_l2(16, 0) == "weekday_evening_night"
    iv = np.array([30.0, 2.0, 2.0, 2.0, np.nan])
    sb = np.array([True, True, False, False, False])
    eb = np.array([True, True, True, False, True])
    assert list(admission_source(iv, sb, eb, 24)) == ["floor", "OR", "ED", "other", "ED"]


def test_attribution_bins_and_gaps():
    coh = _toy_cohort()
    assert list(coh["N_OPEN_BINS"]) == [6, 28, 1]
    ev = pd.DataFrame({"SUBJECT_ID": ["1", "1", "1", "1", "2", "3"],
                       "CHARTTIME": ["2100-01-01 09:00:00", "2100-01-01 10:00:00", "2100-01-02 11:00:00",
                                     "2100-01-03 10:00:00", "2100-01-05 23:45:00", "2100-01-01 10:00:00"]})
    m = attribute(ev, coh)
    # pre-INTIME row kept (H=-1), two window rows for enc 0, one past OUTTIME dropped, enc 2 row kept, subject 3 dropped
    assert sorted(m["ENC_IDX"].tolist()) == [0, 0, 0, 2]
    acc = KeyAcc(len(coh), 100, np.random.default_rng(0), (0, 10))
    acc.update(m["ENC_IDX"].to_numpy(), m["H"].to_numpy(), np.array([1.0, 2.0, 50.0, np.nan]),
               np.array([False, False, False, True]), pd.Series(["u", "u", "u", None]))
    B = bits_to_bins(acc.mask, 28)
    assert B[0, 0] and B[0, 4] and not B[0, 1] and B[2, 0]
    assert acc.n_pre == 1 and acc.n24[0] == 1 and acc.above == 1
    assert np.allclose(gap_hours(np.array([0, 0, 0, 1]), np.array([5.0, 1.0, 2.0, 3.0])), [1.0, 3.0])


def test_cohort_rebuild_24h_rule():
    icu = pd.DataFrame({"SUBJECT_ID": ["1", "1", "1"], "HADM_ID": ["10", "10", "11"], "ICUSTAY_ID": ["1", "2", "3"],
                        "INTIME": ["2100-01-01 00:00:00", "2100-01-03 00:00:00", "2100-01-10 00:00:00"],
                        "OUTTIME": ["2100-01-02 06:00:00", "2100-01-04 00:00:00", "2100-01-11 00:00:00"]})
    c = rebuild_cohort_from_icustays(icu)
    assert len(c) == 2 and c.iloc[0]["HADM_LIST"] == ["10"]


def run_self_tests() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"FAIL {t.__name__}: {type(e).__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} self-tests passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
