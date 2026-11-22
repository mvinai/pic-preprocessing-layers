#!/usr/bin/env python3
"""Encounter-harmonization consistency check: READ-ONLY, AGGREGATE-ONLY.

Compares canonical_icu_encounters.csv, encounter_harmonization_fields.csv, the subject split file and (optionally)
raw ICUSTAYS.csv against the documented definitions (contract v12 H1-H7, v12.1 amendment A1.2/A1.4, preprocessing
log Steps 1.x / 4.6, specs). It prints counts only: no IDs, no timestamps, no row-level data. Every count 1-9 prints
as "<10". A zero prints as 0 (a check that passes). Only disagreement counts and documented totals are printed, never
two overlapping totals whose difference would expose a small number.

Usage (from the project root, Windows paths are fine):
    python scripts/check_encounter_harmonization_consistency.py
    python scripts/check_encounter_harmonization_consistency.py --icustays "C:/.../PICDataset/ICUSTAYS.csv" --json-out data/interim/harmonization_consistency_check.json
    add --surgery-vitals "C:/.../PICDataset/SURGERY_VITAL_SIGNS.csv" for section 13 (approved vitals -> ROW_ID assignment, Nov 11)
Section 12 (needs SURGERY_INFO) applies the approved surgery timing rules (Vinai, Nov 11, 2026): timing time T =
ANES_END_TIME (fallback SURGERY_END_TIME, counted); recent_layer_a if t0-24h <= T <= t0+6h (Layer A only), during_icu_layer_c
if t0+6h < T <= OUTTIME (merge gaps included; in_gap is a descriptor), remote_same_hadm if T < t0-24h (later encounter),
one surgery -> one encounter. Section 13 (needs SURGERY_VITAL_SIGNS) assigns each vitals row to the ROW_ID in its
(SUBJECT_ID, HADM_ID, OPER_ID) group whose [ANES_START, ANES_END] contains MONITOR_TIME (+-30 min sensitivity).
Exit codes: 0 = ran (disagreements are reported, not raised); 2 = input or hash failure.
"""
from __future__ import annotations

import os

import argparse
import ast
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CANONICAL_SHA256 = "e5feeb22b5457ab89c4c73235e6a4722d538aec5f88595b6f049004171286418"
HARMONIZATION_SHA256 = "ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d"
DEFAULT_SPLIT_FILE = os.environ.get("SPLIT_FILE_PATH")
DEFAULT_RAW_DIR = os.environ.get("PIC_RAW_DIR")
TOL_H = 1.0 / 60.0           # 1 minute, in hours
LANDMARKS = (24, 48, 72)
EXPECTED = {}  # public copy: cohort/split pins omitted; see configs/cohort_pins.example.json
# Documented name -> real header (naming-only mismatches are reported, not fixed)
DOC_NAMES = {"icu_span_hours": "span_icu_hours", "sum_icu_stay_hours": "icu_time_hours", "n_gaps": "gap_count_between_stays",
             "total_gap_hours": "gap_hours_total", "age_at_t0": "age_years_at_intime", "age_at_t0_days": "age_years_at_intime",
             "prior_encounters_n": "prior_encounters", "death_attributed": "death_attribution_flag", "OUT": "OUTTIME"}
SPLIT_LABELS = {"train": "train", "training": "train", "val": "validation", "valid": "validation", "validation": "validation",
                "dev": "validation", "test": "test", "testing": "test"}
TRUE = {"1", "1.0", "true", "t", "y", "yes"}


class InputError(Exception):
    pass


def sup(n) -> int | str:
    """Primary suppression for printed counts: 1-9 -> '<10'. Zero stays 0 (a passing check)."""
    n = int(n)
    return "<10" if 0 < n < 10 else n


def partition(parts: dict) -> dict:
    """Suppress a partition whose total is published elsewhere: 1-9 -> '<10', and if any part is suppressed, the
    smallest remaining published part is also withheld ('S'), so no small part can be recovered by subtraction."""
    raw = {k: int(v) for k, v in parts.items()}
    out = {k: sup(v) for k, v in raw.items()}
    if any(o == "<10" for o in out.values()):
        pub = sorted((v, k) for k, v in raw.items() if out[k] != "<10" and v > 0)
        if pub:
            out[pub[0][1]] = "S"
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_id(x) -> str:
    s = str(x).strip().strip("'\"").strip()
    if s.lower() in {"", "nan", "none", "null", "nat"}:
        return ""
    m = re.fullmatch(r"([+-]?\d+)\.0*", s)
    return str(int(m.group(1))) if m else s


def parse_list(v) -> tuple:
    s = str(v).strip()
    if s.lower() in {"", "nan", "none", "null", "[]", "()"}:
        return ()
    items = None
    if s[0] in "[(":
        for loader in (json.loads, ast.literal_eval):
            try:
                obj = loader(s)
            except Exception:  # noqa: BLE001
                continue
            items = list(obj) if isinstance(obj, (list, tuple, set)) else [obj]
            break
    if items is None:
        items = [t for t in re.split(r"[,;|\s]+", s.strip("[]()")) if t]
    out = []
    for it in items:
        n = norm_id(it)
        if n and n not in out:
            out.append(n)
    return tuple(sorted(out))


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False, encoding="utf-8-sig")


def ts(s: pd.Series) -> pd.Series:
    s = s.map(lambda v: str(v).strip()).replace({"": None, "nan": None, "NaT": None})
    try:
        return pd.to_datetime(s, errors="coerce", format="mixed")
    except (TypeError, ValueError):
        return pd.to_datetime(s, errors="coerce")


def num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.map(lambda v: str(v).strip()), errors="coerce")


def flag(s: pd.Series) -> pd.Series:
    return s.map(lambda v: str(v).strip().lower() in TRUE)


def hours(a: pd.Series, b: pd.Series) -> pd.Series:
    return (a - b).dt.total_seconds() / 3600.0


class Report:
    def __init__(self) -> None:
        self.sections: dict[str, dict] = {}

    def add(self, section: str, key: str, value) -> None:
        self.sections.setdefault(section, {})[key] = value

    def count(self, section: str, key: str, n) -> None:
        self.add(section, key, sup(n))

    def vs_expected(self, section: str, key: str, n, expected) -> None:
        """Documented totals: print the total (>=10 always here) and whether it equals the documented value.
        If expected is None (public copy, pins omitted), record the value only."""
        if expected is None:
            self.add(section, key, {"value": sup(n), "documented": None, "matches": None})
            return
        self.add(section, key, {"value": sup(n), "documented": expected, "matches": int(n) == int(expected)})

    def text(self) -> str:
        lines = ["=== Encounter harmonization consistency check (aggregate only; counts 1-9 shown as <10) ==="]
        for sec, d in self.sections.items():
            lines.append(f"\n[{sec}]")
            for k, v in d.items():
                lines.append(f"  {k}: {v}")
        return "\n".join(lines)


def check(canonical: Path, harmonization: Path, split_file: Path | None, icustays: Path | None,
          expect_canonical: str | None = CANONICAL_SHA256, expect_harm: str | None = HARMONIZATION_SHA256,
          expected: dict | None = None, surgery_info: Path | None = None, surgery_vitals: Path | None = None) -> Report:
    exp = dict(EXPECTED, **(expected or {}))
    R = Report()
    # ---------------------------------------------------------------- hashes before reading
    for label, p, e in (("canonical", canonical, expect_canonical), ("harmonization", harmonization, expect_harm)):
        if not Path(p).exists():
            raise InputError(f"{label} file not found")
        got = sha256(p)
        if e is not None and got != e:
            raise InputError(f"{label} SHA-256 mismatch (expected {e[:12]}..., got {got[:12]}...)")
        R.add("0_inputs", f"{label}_sha256", got)
    c = read(canonical)
    h = read(harmonization)
    R.add("0_inputs", "canonical_columns", list(c.columns))
    R.add("0_inputs", "harmonization_columns", list(h.columns))
    present = set(h.columns) | set(c.columns)
    R.add("1_naming", "documented_name -> real header",
          {d: (r if r in present else "ABSENT") + (" (documented name itself present)" if d in present else "")
           for d, r in DOC_NAMES.items()})
    for col in ("encounter_id", "SUBJECT_ID", "t0", "INTIME", "OUTTIME"):
        if col not in h.columns:
            raise InputError(f"harmonization file has no {col!r} column")
    for col in ("encounter_id", "SUBJECT_ID", "INTIME", "OUTTIME"):
        if col not in c.columns:
            raise InputError(f"canonical file has no {col!r} column")
    if c["encounter_id"].duplicated().any() or h["encounter_id"].duplicated().any():
        raise InputError("encounter_id not unique")
    R.vs_expected("2_rows", "canonical_encounters", len(c), exp.get("encounters"))
    R.count("2_rows", "harmonization_rows_minus_canonical_rows_abs", abs(len(h) - len(c)))
    R.count("2_rows", "encounter_ids_only_in_one_file", len(set(c["encounter_id"]) ^ set(h["encounter_id"])))
    m = c.merge(h, on="encounter_id", how="inner", suffixes=("_c", "_h"))
    R.vs_expected("2_rows", "subjects", m["SUBJECT_ID_c"].map(norm_id).nunique(), exp.get("subjects"))

    # ---------------------------------------------------------------- shared columns: canonical vs harmonization
    shared = [col for col in c.columns if col in h.columns and col != "encounter_id"]
    dis = {}
    for col in shared:
        a, b = m[f"{col}_c"], m[f"{col}_h"]
        if col in ("INTIME", "OUTTIME", "hospital_admit_time", "hospital_disch_time"):
            ta, tb = ts(a), ts(b)
            d = ~((ta == tb) | (ta.isna() & tb.isna()))
        elif col == "hadm_id_list":
            d = a.map(parse_list) != b.map(parse_list)
        elif col in ("SUBJECT_ID", "HADM_ID", "ICUSTAY_ID"):
            d = a.map(norm_id) != b.map(norm_id)
        elif col == "split":
            d = a.map(lambda v: SPLIT_LABELS.get(v.strip().lower(), v)) != b.map(lambda v: SPLIT_LABELS.get(v.strip().lower(), v))
        else:
            na, nb = num(a), num(b)
            both = na.notna() & nb.notna()
            d = (both & ((na - nb).abs() > 1e-6)) | (~both & (a.str.strip() != b.str.strip()))
        dis[col] = sup(int(d.sum()))
    R.add("3_canonical_vs_harmonization", "disagreeing_rows_per_shared_column", dis)

    # ---------------------------------------------------------------- time definitions (harmonization file)
    t0, hin, hout = ts(m["t0"]), ts(m["INTIME_h"]), ts(m["OUTTIME_h"])
    cin = ts(m["INTIME_c"])
    R.count("4_times", "t0_missing_or_unparseable", t0.isna().sum())
    R.count("4_times", "harm_OUTTIME_missing_or_unparseable", hout.isna().sum())
    R.count("4_times", "t0_ne_harm_INTIME", (t0 != hin).sum())
    R.count("4_times", "t0_ne_canonical_INTIME", (t0 != cin).sum())
    R.count("4_times", "OUTTIME_le_t0", (hout <= t0).sum())
    span_obs = hours(hout, t0)
    for name, col in (("span_icu_hours", "span_icu_hours"), ("icu_time_hours", "icu_time_hours"),
                      ("icu_los_hours(harm)", "icu_los_hours_h" if "icu_los_hours_h" in m else "icu_los_hours"),
                      ("icu_los_hours(canonical)", "icu_los_hours_c" if "icu_los_hours_c" in m else "icu_los_hours")):
        if col in m.columns:
            v = num(m[col])
            R.count("4_times", f"abs(harm OUTTIME - t0 - {name}) > 1 min", ((span_obs - v).abs() > TOL_H).sum())
            R.count("4_times", f"{name}_missing", v.isna().sum())
        else:
            R.add("4_times", f"{name}", "column absent")
    if {"span_icu_hours", "icu_time_hours"} <= set(m.columns):
        sp, it = num(m["span_icu_hours"]), num(m["icu_time_hours"])
        R.count("4_times", "icu_time_hours > span_icu_hours + 1 min", (it > sp + TOL_H).sum())
        R.count("4_times", "icu_time_hours < span_icu_hours - 1 min (encounters with gap time)", (it < sp - TOL_H).sum())
        if "gap_hours_total" in m.columns:
            g = num(m["gap_hours_total"]).fillna(0)
            R.count("4_times", "abs(span - icu_time - gap_hours_total) > 1 min", ((sp - it - g).abs() > TOL_H).sum())
    for col in ("icu_los_hours_c", "icu_los_hours_h", "icu_los_hours"):
        if col in m.columns and {"span_icu_hours", "icu_time_hours"} <= set(m.columns):
            v = num(m[col])
            is_span = ((v - num(m["span_icu_hours"])).abs() <= TOL_H)
            is_time = ((v - num(m["icu_time_hours"])).abs() <= TOL_H)
            R.add("4_times", f"{col}_matches", {"span_icu_hours": sup(is_span.sum()), "icu_time_hours": sup(is_time.sum()),
                                                "neither": sup((~is_span & ~is_time).sum())})  # overlapping sets, not a partition
            break
    R.vs_expected("4_times", "span_gt_1440h", (span_obs > 1440).sum(), exp.get("span_gt_1440h"))
    if "gap_count_between_stays" in m.columns and "encounter_n_stays_h" in m.columns:
        gc, ns = num(m["gap_count_between_stays"]).fillna(0), num(m["encounter_n_stays_h"]).fillna(1)
        R.count("5_merge", "gap_count > encounter_n_stays - 1", (gc > ns - 1).sum())
        R.count("5_merge", "encounters_with_gap_count > 0", (gc > 0).sum())
    if "encounter_n_stays_h" in m.columns:
        ns = num(m["encounter_n_stays_h"])
        R.add("5_merge", "encounter_n_stays_distribution", partition({"1": (ns == 1).sum(), "2": (ns == 2).sum(), "3+": (ns >= 3).sum()}))
    hl = m["hadm_id_list_h"].map(parse_list) if "hadm_id_list_h" in m.columns else m["hadm_id_list"].map(parse_list)
    R.vs_expected("5_merge", "multi_hadm_encounters", (hl.map(len) > 1).sum(), exp.get("multi_hadm_encounters"))
    hcol = "HADM_ID_h" if "HADM_ID_h" in m.columns else "HADM_ID"
    R.count("5_merge", "representative_HADM_ID_not_in_hadm_id_list", sum(norm_id(a) not in b for a, b in zip(m[hcol], hl)))
    nla = "n_linked_admissions_h" if "n_linked_admissions_h" in m.columns else ("n_linked_admissions" if "n_linked_admissions" in m.columns else None)
    if nla:
        R.count("5_merge", "n_linked_admissions_ne_len(hadm_id_list)", (num(m[nla]) != hl.map(len)).sum())
    fmts = m["hadm_id_list_h"] if "hadm_id_list_h" in m.columns else m["hadm_id_list"]
    R.add("5_merge", "hadm_id_list_text_formats", partition({
        "bracketed": fmts.str.strip().str.startswith("[").sum(),
        "delimited": (~fmts.str.strip().str.startswith("[") & fmts.str.contains(r"[,;|]")).sum(),
        "plain_single": (~fmts.str.strip().str.startswith("[") & ~fmts.str.contains(r"[,;|]")).sum()}))

    # ---------------------------------------------------------------- landmark eligibility indicators
    for L in LANDMARKS:
        col = f"in_icu_at_{L}h"
        if col not in m.columns:
            R.add("6_landmarks", col, "column absent")
            continue
        fl = flag(m[col])
        R.vs_expected("6_landmarks", f"{col}_total", fl.sum(), exp.get(col))
        R.count("6_landmarks", f"{col} != (OUTTIME > t0 + {L}h strictly)", (fl != (hout > t0 + pd.Timedelta(hours=L))).sum())
        R.count("6_landmarks", f"{col} != (OUTTIME >= t0 + {L}h)", (fl != (hout >= t0 + pd.Timedelta(hours=L))).sum())

    # ---------------------------------------------------------------- deaths and attribution
    subj = m["SUBJECT_ID_c"].map(norm_id)
    died_col = "Died_h" if "Died_h" in m.columns else "Died"
    died = flag(m[died_col])
    R.vs_expected("7_deaths", "Died_total", died.sum(), exp.get("Died"))
    if "death_attribution_flag" in m.columns:
        daf = flag(m["death_attribution_flag"])
        R.vs_expected("7_deaths", "death_attribution_flag_total", daf.sum(), exp.get("death_attribution_flag"))
        R.count("7_deaths", "flag_without_Died", (daf & ~died).sum())
        order = pd.DataFrame({"s": subj, "t0": t0, "hl": hl}).reset_index()
        order = order.sort_values(["s", "t0"])
        last_of_subject = ~order.duplicated("s", keep="last")
        # last of its admission: no later encounter of the same subject shares any HADM
        later_shares = []
        for _, g in order.groupby("s", sort=False):
            sets = [set(x) for x in g["hl"]]
            later_shares += [any(sets[i] & sets[j] for j in range(i + 1, len(sets))) for i in range(len(sets))]
        order["not_last_subject"] = ~last_of_subject.to_numpy()
        order["not_last_admission"] = later_shares
        order = order.set_index("index").sort_index()
        R.count("7_deaths", "Died_and_not_last_encounter_of_subject", (died & order["not_last_subject"]).sum())
        R.count("7_deaths", "Died_and_not_last_encounter_of_its_admission", (died & order["not_last_admission"]).sum())
        R.count("7_deaths", "flag != (Died & not last of subject)", (daf != (died & order["not_last_subject"])).sum())
        R.count("7_deaths", "flag != (Died & not last of admission)", (daf != (died & order["not_last_admission"])).sum())

    # ---------------------------------------------------------------- next encounter / readmission rules
    o = pd.DataFrame({"s": subj, "t0": t0, "out": hout, "hl": hl}).sort_values(["s", "t0"])
    nxt = o.groupby("s").shift(-1)
    has_next = nxt["t0"].notna()
    gap = hours(nxt["t0"], o["out"])
    R.count("8_next_encounter", "next_encounter_starts_within_24h_of_OUTTIME (must be 0)", (has_next & (gap <= 24)).sum())
    same_h = [bool(set(a) & set(b)) if isinstance(b, tuple) else False for a, b in zip(o["hl"], nxt["hl"])]
    R.count("8_next_encounter", "next_encounter_in_same_HADM (same-HADM ICU readmissions, all splits, all encounters)", int(np.sum(same_h)))
    if "prior_encounters" in m.columns:
        rank = o.groupby("s").cumcount()
        pe = num(m["prior_encounters"]).reindex(o.index)
        R.count("8_next_encounter", "prior_encounters != number of earlier encounters of the subject", (pe != rank).sum())

    # ---------------------------------------------------------------- splits
    if "split" in c.columns and "split" in h.columns:
        sc = c.set_index("encounter_id")["split"].map(lambda v: SPLIT_LABELS.get(v.strip().lower(), "?"))
        R.count("9_split", "canonical_split_unknown_label", (sc == "?").sum())
        for sp in ("train", "validation", "test"):
            R.vs_expected("9_split", f"canonical_{sp}_encounters", (sc == sp).sum(), exp.get(sp))
    if split_file is not None and Path(split_file).exists():
        sf = read(split_file)
        sf = pd.DataFrame({"s": sf["SUBJECT_ID"].map(norm_id), "sp": sf["split"].map(lambda v: SPLIT_LABELS.get(v.strip().lower(), "?"))})
        R.count("9_split", "split_file_subjects_in_more_than_one_split", (sf.drop_duplicates().groupby("s").size() > 1).sum())
        mp = sf.drop_duplicates("s").set_index("s")["sp"]
        enc_sp = subj.map(mp)
        R.count("9_split", "cohort_encounters_missing_from_split_file", enc_sp.isna().sum())
        for side, col in (("canonical", "split_c"), ("harmonization", "split_h")):
            if col in m.columns:
                lab = m[col].map(lambda v: SPLIT_LABELS.get(v.strip().lower(), "?"))
                R.count("9_split", f"{side}_split != split_file", (lab != enc_sp).sum())
        R.add("9_split", "split_file_sha256", sha256(split_file))
    else:
        R.add("9_split", "split_file", "not provided / not found")

    # ---------------------------------------------------------------- age, admit-time
    if "age_years_at_intime_h" in m.columns or "age_years_at_intime" in m.columns:
        a = num(m["age_years_at_intime_h"] if "age_years_at_intime_h" in m.columns else m["age_years_at_intime"])
        R.count("10_other", "age_years_missing", a.isna().sum())
        R.count("10_other", "age_years_negative", (a < 0).sum())
    hat = ts(m["hospital_admit_time_h"] if "hospital_admit_time_h" in m.columns else m["hospital_admit_time"]) if (
        "hospital_admit_time_h" in m.columns or "hospital_admit_time" in m.columns) else None
    if hat is not None:
        d = hours(hat, t0)
        R.count("10_other", "hospital_admit_time_after_t0", (d > 0).sum())
        R.count("10_other", "hospital_admit_time_after_t0_by_over_1h", (d > 1).sum())

    # ---------------------------------------------------------------- optional: rebuild from raw ICUSTAYS
    if icustays is not None and Path(icustays).exists():
        raw = read(icustays)
        st = pd.DataFrame({"s": raw["SUBJECT_ID"].map(norm_id), "hadm": raw["HADM_ID"].map(norm_id),
                           "id": raw["ICUSTAY_ID"].map(norm_id), "in": ts(raw["INTIME"]), "out": ts(raw["OUTTIME"])})
        valid = st["out"] > st["in"]
        R.vs_expected("11_rebuild_from_ICUSTAYS", "valid_raw_stays (OUTTIME > INTIME)", valid.sum(), exp.get("valid_raw_stays"))
        st = st[valid].sort_values(["s", "in", "out"]).reset_index(drop=True)
        enc_no, run_end, cur_s, n = [], None, None, -1
        nested = overlaps = 0
        for s_, i_, o_ in zip(st["s"], st["in"], st["out"]):
            if s_ != cur_s or i_ > run_end + pd.Timedelta(hours=24):
                n += 1
                cur_s, run_end = s_, o_
            else:
                if i_ < run_end:
                    overlaps += 1
                    nested += int(o_ <= run_end)
                run_end = max(run_end, o_)
            enc_no.append(n)
        st["e"] = enc_no
        st["dur"] = hours(st["out"], st["in"])
        g = st.groupby("e")
        rb = pd.DataFrame({"s": g["s"].first(), "t0": g["in"].min(), "max_out": g["out"].max(),
                           "last_stay_out": g["out"].last(), "n_stays": g.size(),
                           "first_id": g["id"].first(), "hadms": g["hadm"].agg(lambda x: tuple(sorted(set(v for v in x if v)))),
                           "sum_stay_h": g["dur"].sum()})
        R.vs_expected("11_rebuild_from_ICUSTAYS", "rebuilt_encounters (running max OUTTIME, <=24 h)", len(rb), exp.get("encounters"))
        R.count("11_rebuild_from_ICUSTAYS", "overlapping_stays", overlaps)
        R.count("11_rebuild_from_ICUSTAYS", "nested_stays", nested)
        R.count("11_rebuild_from_ICUSTAYS", "encounters_where_max_OUTTIME != last_stay_OUTTIME (nested last stay)",
                (rb["max_out"] != rb["last_stay_out"]).sum())
        key = pd.DataFrame({"s": subj, "t0": t0, "out": hout, "hl": hl,
                            "ns": num(m["encounter_n_stays_h"]) if "encounter_n_stays_h" in m.columns else np.nan,
                            "icu_time": num(m["icu_time_hours"]) if "icu_time_hours" in m.columns else np.nan})
        j = key.merge(rb, on=["s", "t0"], how="outer", indicator=True)
        R.count("11_rebuild_from_ICUSTAYS", "cohort_encounters_without_rebuilt_match_on_(SUBJECT_ID, t0)", (j["_merge"] == "left_only").sum())
        R.count("11_rebuild_from_ICUSTAYS", "rebuilt_encounters_not_in_cohort", (j["_merge"] == "right_only").sum())
        b = j[j["_merge"] == "both"]
        R.count("11_rebuild_from_ICUSTAYS", "harm_OUTTIME != rebuilt max OUTTIME", (b["out"] != b["max_out"]).sum())
        R.count("11_rebuild_from_ICUSTAYS", "harm_OUTTIME != OUTTIME of the last stay by INTIME", (b["out"] != b["last_stay_out"]).sum())
        R.count("11_rebuild_from_ICUSTAYS", "encounter_n_stays != rebuilt", (b["ns"] != b["n_stays"]).sum())
        R.count("11_rebuild_from_ICUSTAYS", "hadm_id_list != rebuilt stay HADM set", (b["hl"] != b["hadms"]).sum())
        R.count("11_rebuild_from_ICUSTAYS", "abs(icu_time_hours - sum of stay durations) > 1 min", ((b["icu_time"] - b["sum_stay_h"]).abs() > TOL_H).sum())
    else:
        R.add("11_rebuild_from_ICUSTAYS", "status", "skipped (pass --icustays to run)")
        st = None

    # ---------------------------------------------------------------- optional: surgery timing (v1.7 rules, Vinai Nov 11)
    if surgery_info is not None and Path(surgery_info).exists():
        S12 = "12_surgery_vs_encounters"
        si = read(surgery_info)
        sg = pd.DataFrame({"ROW_ID": si["ROW_ID"].map(norm_id), "s": si["SUBJECT_ID"].map(norm_id), "h": si["HADM_ID"].map(norm_id),
                           "o": si["OPER_ID"].map(norm_id), "B": ts(si["SURGERY_BEGIN_TIME"]), "E": ts(si["SURGERY_END_TIME"]),
                           "AS": ts(si["ANES_START_TIME"]), "AE": ts(si["ANES_END_TIME"])})
        R.count(S12, "SURGERY_INFO_ROW_ID_blank", (sg["ROW_ID"] == "").sum())
        R.count(S12, "SURGERY_INFO_ROW_ID_repeated", sg["ROW_ID"].duplicated().sum())
        g = sg.groupby(["s", "h", "o"]).size()
        R.count(S12, "shared_(SUBJECT_ID,HADM_ID,OPER_ID)_groups", (g >= 2).sum())
        R.count(S12, "rows_in_shared_groups", g[g >= 2].sum())
        # timing time T = anaesthesia end, fallback surgery end (counted)
        sg["T"] = sg["AE"].where(sg["AE"].notna(), sg["E"])
        R.count(S12, "timing_fallback_to_SURGERY_END_TIME", (sg["AE"].isna() & sg["E"].notna()).sum())
        R.count(S12, "no_timing_time (anaesthesia end and surgery end both missing)", sg["T"].isna().sum())
        enc_k = pd.DataFrame({"e": np.arange(len(m)), "s": subj, "t0": t0, "out": hout, "hl": hl})
        ek = enc_k.explode("hl").rename(columns={"hl": "h"})
        j = sg[sg["T"].notna()].reset_index().merge(ek, on=["s", "h"], how="inner")
        # rule 1: t0 - 24 h <= T <= OUT (recent if T <= t0 + 6 h, else Layer C); rule 2: next encounter with t0 > T (remote)
        w = j[(j["T"] >= j["t0"] - pd.Timedelta(hours=24)) & (j["T"] <= j["out"])].sort_values(["index", "t0"]).drop_duplicates("index")
        r2 = j[~j["index"].isin(w["index"]) & (j["t0"] > j["T"])].sort_values(["index", "t0"]).drop_duplicates("index")
        a = pd.concat([w.assign(rule=1), r2.assign(rule=2)])
        after = hours(a["T"], a["t0"])
        a["cls"] = np.where(a["rule"] == 2, "remote_same_hadm", np.where(after <= 6, "recent_layer_a", "during_icu_layer_c"))
        n_cand = j.groupby("index").size()
        post = set(n_cand.index) - set(a["index"])
        R.add(S12, "one_surgery_one_encounter_classes", partition({
            "recent_layer_a (end in [t0-24h, t0+6h]; Layer A only)": (a["cls"] == "recent_layer_a").sum(),
            "remote_same_hadm (end < t0-24h; Layer A)": (a["cls"] == "remote_same_hadm").sum(),
            "during_icu_layer_c (end in (t0+6h, OUT])": (a["cls"] == "during_icu_layer_c").sum(),
            "post_icu_same_hadm": len(post)}))
        R.count(S12, "recent_ending_after_t0 (end in (t0, t0+6h]; would have been Layer C by begin time)",
                ((a["cls"] == "recent_layer_a") & (after > 0)).sum())
        R.count(S12, "surgeries_with_2+_candidate_encounters", (n_cand >= 2).sum())
        # explicit reuse rule: a Layer C / recent surgery of one encounter is Layer A history of a later same-HADM encounter
        later = a[a["cls"] != "remote_same_hadm"].merge(ek[["e", "s", "h", "t0"]].rename(columns={"e": "e2", "t0": "t0_2"}), on=["s", "h"])
        R.count(S12, "assigned_surgeries_that_are_also_Layer_A_history_of_a_later_same_HADM_encounter",
                later[later["t0_2"] > later["t0"]]["index"].nunique())
        if st is not None:
            lc = a[a["cls"] == "during_icu_layer_c"][["index", "s", "T"]].merge(st[["s", "in", "out"]], on="s", how="left")
            in_stay = set(lc[(lc["T"] >= lc["in"]) & (lc["T"] <= lc["out"])]["index"])
            n_lc = int((a["cls"] == "during_icu_layer_c").sum())
            R.count(S12, "in_gap (Layer C surgeries ending in a merge gap; descriptor only, kept in Layer C)", n_lc - len(in_stay))
        else:
            R.add(S12, "in_gap", "skipped (needs --icustays)")
        # ------------------------------------------------------------ optional: SURGERY_VITAL_SIGNS -> ROW_ID assignment
        if surgery_vitals is not None and Path(surgery_vitals).exists():
            hdr = pd.read_csv(surgery_vitals, nrows=0, encoding="utf-8-sig").columns.tolist()
            mt = [c for c in ("MONITOR_TIME", "MONITORTIME") if c in hdr]
            if len(mt) != 1:
                raise InputError("SURGERY_VITAL_SIGNS monitor-time column must resolve to exactly one of MONITOR_TIME / MONITORTIME")
            sv = pd.read_csv(surgery_vitals, dtype=str, keep_default_na=False, na_filter=False, encoding="utf-8-sig",
                             usecols=["SUBJECT_ID", "HADM_ID", "OPER_ID", mt[0]])
            v = pd.DataFrame({"s": sv["SUBJECT_ID"].map(norm_id), "h": sv["HADM_ID"].map(norm_id), "o": sv["OPER_ID"].map(norm_id),
                              "t": ts(sv[mt[0]])}).reset_index()
            V = "13_vitals_to_ROW_ID"
            keys = g.rename("gsize").reset_index()
            vk = v.merge(keys, on=["s", "h", "o"], how="left")
            R.count(V, "vitals_rows_total", len(v))
            R.count(V, "vitals_rows_with_no_SURGERY_INFO_key", vk["gsize"].isna().sum())
            vk = vk[vk["gsize"].notna()]
            R.count(V, "vitals_rows_in_shared_OPER_ID_groups", (vk["gsize"] >= 2).sum())
            # approved rule (Nov 11): each row -> the ROW_ID in its key group whose [ANES_START, ANES_END] contains MONITOR_TIME.
            # A missing anaesthesia bound falls back to the surgery begin/end time (same fallback as the timing rule; counted).
            sg["wS"] = sg["AS"].where(sg["AS"].notna(), sg["B"])
            sg["wE"] = sg["AE"].where(sg["AE"].notna(), sg["E"])
            R.count(V, "surgeries_whose_window_uses_a_surgery_time_fallback", ((sg["AS"].isna() & sg["B"].notna()) | (sg["AE"].isna() & sg["E"].notna())).sum())
            R.count(V, "surgeries_with_no_usable_window", (sg["wS"].isna() | sg["wE"].isna()).sum())
            for label, tol in (("anesthesia_window", 0.0), ("anesthesia_window_pm30min", 0.5)):
                c = vk[["index", "s", "h", "o", "t", "gsize"]].merge(sg[["s", "h", "o", "wS", "wE"]], on=["s", "h", "o"], how="left")
                hit = (c["t"] >= c["wS"] - pd.Timedelta(hours=tol)) & (c["t"] <= c["wE"] + pd.Timedelta(hours=tol))
                nh = hit.groupby(c["index"]).sum().reindex(vk["index"], fill_value=0).to_numpy()
                shared = (vk["gsize"] >= 2).to_numpy()
                # disjoint scopes (no overlapping totals that could be differenced)
                for scope, msk in (("unique_key_rows", ~shared), ("shared_group_rows", shared)):
                    R.add(V, f"{scope}_by_{label}", partition({"assigned_to_exactly_one_ROW_ID": (nh[msk] == 1).sum(),
                                                               "ambiguous_2+_windows": (nh[msk] >= 2).sum(),
                                                               "unassigned_0_windows": (nh[msk] == 0).sum()}))
        elif surgery_vitals is not None:
            R.add("13_vitals_to_ROW_ID", "status", "SURGERY_VITAL_SIGNS not found")
    else:
        R.add("12_surgery_vs_encounters", "status", "skipped (pass --surgery-info to run)")
    return R


def main(argv=None) -> int:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    here = Path(__file__).resolve()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", default=str(here.parents[1]))
    ap.add_argument("--canonical")
    ap.add_argument("--harmonization")
    ap.add_argument("--split-file", default=DEFAULT_SPLIT_FILE)
    ap.add_argument("--icustays", default=(str(Path(DEFAULT_RAW_DIR) / "ICUSTAYS.csv") if DEFAULT_RAW_DIR else None))
    ap.add_argument("--surgery-info", default=(str(Path(DEFAULT_RAW_DIR) / "SURGERY_INFO.csv") if DEFAULT_RAW_DIR else None))
    ap.add_argument("--surgery-vitals", default=None, help="optional SURGERY_VITAL_SIGNS.csv (large) for section 13")
    ap.add_argument("--json-out")
    ap.add_argument("--expect-canonical-sha", default=CANONICAL_SHA256)
    ap.add_argument("--expect-harm-sha", default=HARMONIZATION_SHA256)
    a = ap.parse_args(argv)
    root = Path(a.project_root)
    can = Path(a.canonical) if a.canonical else root / "data" / "interim" / "canonical_icu_encounters.csv"
    harm = Path(a.harmonization) if a.harmonization else root / "data" / "interim" / "encounter_harmonization_fields.csv"
    try:
        R = check(can, harm, Path(a.split_file) if a.split_file else None, Path(a.icustays) if a.icustays else None,
                  a.expect_canonical_sha or None, a.expect_harm_sha or None,
                  surgery_info=Path(a.surgery_info) if a.surgery_info else None,
                  surgery_vitals=Path(a.surgery_vitals) if a.surgery_vitals else None)
    except InputError as e:
        print(f"CHECK STOPPED: {e}", file=sys.stderr)
        return 2
    print(R.text())
    if a.json_out:
        Path(a.json_out).write_text(json.dumps(R.sections, indent=1, default=str), encoding="utf-8")
        print(f"\naggregate JSON written: {a.json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
