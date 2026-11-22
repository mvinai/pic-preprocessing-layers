"""Tests for feature_sparsity_frequency_audit.py (Stage 2). SYNTHETIC DATA ONLY - never real data.

Run:  python -m pytest -q test_feature_sparsity_frequency_audit.py
The end-to-end tests build a fake project with synthetic_stage2_data.py (pins and expected counts are overridden
for the synthetic files); the unit tests exercise every decided rule on small hand-made inputs."""
from __future__ import annotations

import json
import math
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import feature_sparsity_frequency_audit as M  # noqa: E402
import synthetic_stage2_data as S  # noqa: E402


def _find_v17() -> Path:
    for p in (os.environ.get("V17_MODULE_PATH"), HERE / M.V17_MODULE_FILE, HERE.parent / M.V17_MODULE_FILE,
              HERE.parent / "duration_audit_v2" / M.V17_MODULE_FILE):
        if p and Path(p).exists():
            return Path(p)
    pytest.skip("v1.7 module audit_surgery_complexity_duration.py not found next to the script")


V17 = _find_v17()
A = M.load_v17(V17)


# ----------------------------------------------------------------------------- helpers
def mk_cohort(out_h, split=None, gaps=None):
    n = len(out_h)
    out_h = np.asarray(out_h, dtype=float)
    t0 = pd.Timestamp("2016-01-01 06:00:00")
    enc = pd.DataFrame({"enc_i": np.arange(n), "encounter_id": [f"e{i}" for i in range(n)],
                        "SUBJECT_ID": [str(100 + i) for i in range(n)], "hadm_list": [[str(500 + i)] for i in range(n)],
                        "t0": t0, "OUT": [t0 + pd.Timedelta(hours=float(h)) for h in out_h],
                        "split": split or ["train"] * n, "out_h": out_h, "W": np.minimum(out_h, M.MAX_H)})
    return M.Cohort.build(enc, gaps or {})


def ev_df(rows):
    """rows: (enc_i, th, item, vnum, unit[, vtext])"""
    out = []
    for r in rows:
        e, th, it, v, u = r[:5]
        out.append({"enc_i": e, "th": float(th), "item": it, "vnum": np.nan if v is None else float(v),
                    "vtext": r[5] if len(r) > 5 else ("" if v is None else str(v)), "unit": u})
    d = pd.DataFrame(out, columns=["enc_i", "th", "item", "vnum", "vtext", "unit"])
    return d.astype({"enc_i": np.int32})


@pytest.fixture(scope="module")
def maps(tmp_path_factory):
    d = tmp_path_factory.mktemp("maps")
    S.write_maps(d)
    rd = lambda f: pd.read_csv(d / f, dtype=str, keep_default_na=False)  # noqa: E731
    return {"channel": M.parse_channel_map(rd("feature_channel_map_v1.csv")),
            "drug": M.parse_drug_map(rd("drug_group_map_v1.csv")),
            "specimen": M.parse_specimen_map(rd("specimen_map_v1.csv")),
            "organism": M.parse_organism_map(rd("organism_map_v1.csv")),
            "neuro": rd("neuro_admitting_dx_icd10_v1.csv"), "dir": d}


# ============================================================================= pure rules
def test_psofa_cut_points_pf_d19_and_sf_d15():
    pf = M.psofa_score(np.array([400, 399.9, 300, 299.9, 200, 199.9, 100, 99.9, 5]), M.PSOFA_PF)
    assert pf.tolist() == [0, 1, 1, 2, 2, 3, 3, 4, 4]
    sf = M.psofa_score(np.array([292, 291.9, 264, 263.9, 221, 220.9, 148, 147.9]), M.PSOFA_SF)
    assert sf.tolist() == [0, 1, 1, 2, 2, 3, 3, 4]


def test_outcome_groups_edges_and_died():
    out_h = np.array([24 + 24, 24 + 24.01, 24 + 120, 24 + 120.01, 24 + 336, 24 + 336.5, 100.0])
    died = np.array([0, 0, 0, 0, 0, 0, 1], dtype=bool)
    g = M.outcome_group(died, out_h, 24)
    # R=24 -> level 0; R=24.01 -> 1; R=120 -> 3; R=120.01 -> 4; R=336 -> 6; R>336 -> 7; died -> 8
    assert g.tolist() == ["0-2", "0-2", "3-5", "3-5", "6-7", "6-7", "8"]


def test_age_bands_boundaries():
    d = np.array([0, 27.9, 28, 364, 365.25, 3 * 365.25, 6 * 365.25, 12 * 365.25, 18 * 365.25])
    assert M.age_band(d).tolist() == ["neonate", "neonate", "infant", "infant", "toddler", "preschool", "school",
                                      "adolescent", "ge18"]


def test_age_band_invalid_for_nan_and_negative():
    d = np.array([0.0, np.nan, -0.5, np.inf, -40.0])
    assert M.age_band(d).tolist() == ["neonate", "age_invalid", "age_invalid", "age_invalid", "age_invalid"]
    assert "age_invalid" in M.AGE_BANDS


def test_classify_ages_missing_unparseable_negative_bins_and_valid():
    t0 = pd.Series(pd.to_datetime(["2016-01-10 00:00:00"] * 8))
    raw = pd.Series(["", "NaT", "not-a-date", "2016-01-10 12:00:00", "2016-01-20 00:00:00",
                     "2016-02-20 00:00:00", "2008-01-10 00:00:00", "2012-01-10 00:00:00"])
    parsed = M.fast_time(raw)
    has = np.array([True, True, True, True, True, True, True, False])
    status, age, nbin = M.classify_ages(t0, parsed, raw, has)
    assert status.tolist() == ["missing", "missing", "unparseable", "negative", "negative", "negative", "valid",
                               "missing"]
    assert nbin.tolist() == ["", "", "", "lt1d", "1to30d", "gt30d", "", ""]
    assert np.isnan(age[:6]).all() and np.isnan(age[7])
    assert np.isfinite(age[6]) and age[6] == pytest.approx(8 * 365.25, rel=0, abs=2)
    assert M.age_band(age).tolist() == ["age_invalid"] * 6 + ["school"] + ["age_invalid"]


def test_invalid_age_rate_stops_above_one_percent():
    M.check_invalid_age_rate(4, 400)
    M.check_invalid_age_rate(100, 10000)
    with pytest.raises(M.AuditError, match="safety gate"):
        M.check_invalid_age_rate(5, 400)
    with pytest.raises(M.AuditError, match="D20"):
        M.check_invalid_age_rate(101, 10000)
    assert M.invalid_age_limit(10000) == 100


def test_weight_fallback_skips_age_invalid_band():
    # 12 train neonates M charted at 3 kg, 12 train infants F at 8 kg; two uncharted age_invalid rows.
    n = 26
    out_h = [48.0] * n
    split = ["train"] * n
    co = mk_cohort(out_h, split=split)
    ab = np.array(["neonate"] * 12 + ["infant"] * 12 + ["age_invalid", "age_invalid"])
    sx = np.array(["M"] * 12 + ["F"] * 12 + ["M", "unrecognized"])
    la = pd.DataFrame({"enc_i": np.arange(n), "age_band": ab, "sex": sx})
    wt = np.array([3.0] * 12 + [8.0] * 12 + [np.nan, np.nan])
    w, src = M.impute_weights(co, la, wt)
    assert (src[:24] == "charted").all()
    assert src[24] == "imputed_sex_only" and w[24] == pytest.approx(3.0)
    assert src[25] == "imputed_overall" and w[25] == pytest.approx(np.median(wt[:24]))
    # a valid-age uncharted neonate M uses the age-band x sex median, not a later fallback
    wt2 = wt.copy(); wt2[0] = np.nan
    w2, src2 = M.impute_weights(co, la, wt2)
    assert src2[0] == "imputed_age_sex" and w2[0] == pytest.approx(3.0)


def test_icd_normalisation_statuses():
    assert M.norm_icd("q28.200") == ("Q28200", "ok")
    assert M.norm_icd("") == ("", "missing")
    assert M.norm_icd("XYZ")[1] == "invalid"
    assert M.norm_icd("J18.9,G40.1") == ("J189", "multiple")


def test_neuro_q28_four_character_rule(maps):
    m = M.neuro_matcher(maps["neuro"])
    yes = ["S06001", "Q282", "Q2820", "Q2830", "G00", "A39000"]
    no = ["Q280", "Q28", "Q2800", "G401", "G91", "J189", "", "123"]
    assert all(m(c) for c in yes)
    assert not any(m(c) for c in no)


def test_neuro_unknown_rule_stops(maps):
    bad = maps["neuro"].copy()
    bad.loc[bad["icd10_3char"] == "Q28", "match_rule"] = "some other rule"
    with pytest.raises(M.AuditError):
        M.neuro_matcher(bad)


def test_icd_domain_draft_defaults():
    assert M.icd_domain("Q21") == "cardiac_congenital_heart"
    assert M.icd_domain("Q05") == "other_medical"       # Q00-Q07 default (open point)
    assert M.icd_domain("I50") == "other_medical"       # chapter I default
    assert M.icd_domain("D30") == "other_medical" and M.icd_domain("D61") == "oncology_hematology"
    assert M.icd_domain("") == "other_medical"


def test_locf_families_and_caps():
    assert M.locf_family("creatinine", False) == ("4_bins_24h", 24.0)
    assert M.locf_family("sodium", False) == ("1_bin_6h", 6.0)
    assert M.locf_family("pf_ratio", False) == ("1_bin_6h", 6.0)
    assert M.locf_family("troponin_t", True) == ("triggered_no_carry", 24.0)
    assert M.locf_family("heart_rate", False) == ("no_carry", 24.0)


def test_date_only_start_assignment():
    # t0 06:00 -> the t0 date midnight is -6 h; later midnights 18, 42 h; earlier midnight -30 h
    th = np.array([-6.0, 18.0, 42.0, -30.0, 3.5])
    do = np.array([True, True, True, True, False])
    got = M.assign_date_only_start(th, do)
    assert got.tolist() == [0.0, 18.0, 42.0, -30.0, 3.5]
    # a midnight that is not a bin start rounds up to the next bin start
    assert M.assign_date_only_start(np.array([20.0]), np.array([True])).tolist() == [24.0]


# ============================================================================= map parsing
def test_channel_map_parse_merge_mains_and_gates(maps):
    cm = maps["channel"]
    assert cm["merge_main"]["creatinine"] == "5041" and cm["merge_main"]["po2"] == "5239"
    assert cm["merge_main"]["lipase"] == "6796"
    assert cm["derived_gates"]["map"] == {"low": 15.0, "high": 200.0}
    assert ("CHARTEVENTS", "1010") not in cm["items"]  # dropped rows unused
    assert cm["channels"]["heart_rate"]["items"] == ["1002", "1003"]


def test_unknown_conversion_text_stops(maps):
    df = pd.read_csv(maps["dir"] / "feature_channel_map_v1.csv", dtype=str, keep_default_na=False)
    df.loc[df["itemid"] == "5099", "unit_conversion"] = "g/L -> g/dL /10 (edited)"
    with pytest.raises(M.AuditError, match="unrecognized unit_conversion"):
        M.parse_channel_map(df)


def test_drug_map_osmotic_rule_parse(maps):
    dm = maps["drug"]
    assert dm["20%甘露醇注射液 | 20% Mannitol Injection"]["otherwise"] == "diuretic"
    assert dm["10%氯化钠注射液 | 10% Sodium Chloride Injection"]["otherwise"] == "electrolyte_replacement"
    df = pd.read_csv(maps["dir"] / "drug_group_map_v1.csv", dtype=str, keep_default_na=False)
    df.loc[df["group"] == "osmotic_conditional", "conditional_rule"] = "otherwise -> diuretic"
    with pytest.raises(M.AuditError):
        M.parse_drug_map(df)


def test_specimen_base_code_must_be_first_token(maps):
    df = pd.read_csv(maps["dir"] / "specimen_map_v1.csv", dtype=str, keep_default_na=False)
    df.loc[1, "base_code"] = "LIS0204"
    with pytest.raises(M.AuditError):
        M.parse_specimen_map(df)


def test_organism_result_types():
    assert M.organism_result_types("Gram-positive coccus") == ("gram_class", "identification")
    assert M.organism_result_types("smear finding (not an isolate)") == ("gram_class",)
    assert M.organism_result_types("normal flora / mixed") == ("negative_final",)
    assert M.organism_result_types("blank (no organism recorded)") == ()
    with pytest.raises(M.AuditError):
        M.organism_result_types("something new")


def test_real_v1_maps_parse_if_available():
    """Parses the REAL v1 maps (decision files, no patient rows) when present; skipped otherwise."""
    d = Path(os.environ.get("MAPS_DIR", HERE.parents[1] / "data" / "interim" / "maps"))
    if not (d / M.MAP_FILES["channel"]).exists():
        pytest.skip("real v1 maps not present")
    for k, fn in M.MAP_FILES.items():
        assert M.sha256_file(d / fn) == M.MAP_SHA256[k], fn
    rd = lambda f: pd.read_csv(d / f, dtype=str, keep_default_na=False).apply(lambda s: s.str.strip())  # noqa: E731
    cm = M.parse_channel_map(rd(M.MAP_FILES["channel"]))
    assert cm["merge_main"]["hemoglobin"] == "5099"
    M.parse_drug_map(rd(M.MAP_FILES["drug"]))
    M.parse_specimen_map(rd(M.MAP_FILES["specimen"]))
    M.parse_organism_map(rd(M.MAP_FILES["organism"]))
    m = M.neuro_matcher(rd(M.MAP_FILES["neuro"]))
    assert m("Q282") and not m("Q280")


# ============================================================================= values
def test_unit_conversions_per_row(maps):
    cm = maps["channel"]
    d = ev_df([(0, 1, "5257", 12.0, "g/dl"), (0, 2, "5257", 120.0, "g/l"), (0, 3, "5257", 12.0, "g/l"),
               (0, 4, "5257", 5.0, "mg/dl")])
    v, st = M.classify_item_rows(cm["items"][("LABEVENTS", "5257")], d)
    assert v[:2].tolist() == [120.0, 120.0]
    assert st.tolist() == ["measurement", "measurement", "out_of_gate", "unit_unrecognized"]
    d = ev_df([(0, 1, "6525", 50.0, "pg/ml"), (0, 1, "6525", 0.05, "ng/ml")])
    v, st = M.classify_item_rows(cm["items"][("LABEVENTS", "6525")], d)
    assert np.allclose(v, [0.05, 0.05]) and (st == "measurement").all()


def test_magnesium_ug_dl_rows_held_out(maps):
    d = ev_df([(0, 1, "5059", 0.8, "mmol/l"), (0, 2, "5059", 2.0, "ug/dl")])
    _, st = M.classify_item_rows(maps["channel"]["items"][("LABEVENTS", "5059")], d)
    assert st.tolist() == ["measurement", "unit_pending"]


def test_ungated_text_item_counts_text(maps):
    d = ev_df([(0, 1, "6473", None, "", "阴性"), (0, 2, "6473", None, "", "")])
    _, st = M.classify_item_rows(maps["channel"]["items"][("LABEVENTS", "6473")], d)
    assert st.tolist() == ["measurement", "blank"]


def test_hr_merge_keeps_1003_at_same_time(maps):
    chart = ev_df([(0, 1, "1003", 120, "bpm"), (0, 1, "1002", 118, "bpm"), (0, 2, "1002", 110, "bpm"),
                   (1, 1, "1003", 100, "bpm")])
    vb = M.build_values(chart, ev_df([]), maps["channel"])
    hr = vb.events["heart_rate"]
    assert sorted(zip(hr["enc_i"], hr["tc"], hr["item"])) == [(0, 1.0, "1003"), (0, 2.0, "1002"), (1, 1.0, "1003")]
    assert len(vb.hr["pairs"]) == 1 and vb.hr["pairs"]["v"].iloc[0] == 118


def test_map_derivation_and_rejections(maps):
    chart = ev_df([(0, 1, "1016", 90, "mmHg"), (0, 1, "1015", 60, "mmHg"),      # valid: (90+120)/3 = 70
                   (0, 2, "1016", 70, "mmHg"), (0, 2, "1015", 70, "mmHg"),      # dbp >= sbp
                   (0, 3, "1016", 260, "mmHg"), (0, 3, "1015", 60, "mmHg"),     # sbp outside its gate
                   (0, 4, "1016", 90, "mmHg")])                                  # unpaired
    vb = M.build_values(chart, ev_df([]), maps["channel"])
    pr = vb.map_pairs["pairs"].sort_values("tc")
    assert pr["reason"].tolist() == ["valid", "dbp_ge_sbp", "outside_gate"]
    assert vb.events["map"]["v"].tolist() == [70.0]


def test_pf_same_draw_and_po2_preference(maps):
    lab = ev_df([(0, 1, "5222", 50, "%"), (0, 1, "5239", 100, "mmhg"), (0, 1, "5244", 80, "mmhg"),
                 (0, 1, "5243", 199, "mmhg"),
                 (0, 2, "5222", 50, "%"), (0, 2.5, "5239", 100, "mmhg")])   # not the same draw -> no P/F
    vb = M.build_values(ev_df([]), lab, maps["channel"])
    pf = vb.pf["pairs"]
    assert len(pf) == 1 and pf["v"].iloc[0] == 200.0 and pf["v_lab"].iloc[0] == 199
    assert vb.events["pf_ratio"]["t"].tolist() == [2.0]   # lab release +1 h


def test_sf_pairing_5min_spo2_le_97_and_no_pf(maps):
    m = 1 / 60
    lab = ev_df([(0, 1, "5222", 40, "%"),                                   # SpO2 4 min later, 92 -> S/F 230
                 (0, 10, "5222", 40, "%"),                                  # nearest SpO2 is 5.5 min away -> none
                 (0, 20, "5222", 40, "%"),                                  # SpO2 99 -> rejected (> 97)
                 (0, 30, "5222", 50, "%"), (0, 30, "5239", 100, "mmhg"),    # valid P/F at the draw -> no S/F
                 (0, 40, "5222", 40, "%")])                                 # SpO2 exactly 97 -> kept
    chart = ev_df([(0, 1 + 4 * m, "1006", 92, "%"), (0, 10 + 5.5 * m, "1006", 90, "%"), (0, 20 - 2 * m, "1006", 99, "%"),
                   (0, 30, "1006", 90, "%"), (0, 40 + 5 * m, "1006", 97, "%")])
    vb = M.build_values(chart, lab, maps["channel"])
    r = vb.sf["rows"].sort_values("tc")
    assert r["reason"].tolist() == ["sf", "no_spo2_within_5min", "spo2_gt_97", "sf"]
    sf = vb.events["sf_ratio"].sort_values("tc")
    assert np.allclose(sf["v"], [230.0, 242.5])
    # availability = later of the FiO2 release (+1 h) and the SpO2 chart time
    assert np.allclose(sf["t"], [2.0, 41.0])


# ============================================================================= event features
def test_bin_and_landmark_boundaries():
    co = mk_cohort([24.0, 24.5, 100.0])
    e = np.array([1, 1, 2, 2, 2])
    t = np.array([24.0, 23.99, 0.0, 5.999, 6.0])
    F = M.event_features(co, e, t)
    assert not co.elig[24][0] and co.elig[24][1]                 # OUT = t0 + 24 h is not eligible at 24 h
    assert F["avail"][1.0][("full", 24)][1]                      # 23.99 counts at 24 h
    occ = F["occ"]
    assert occ[1, 4] == 1 and occ[1, 3] == 1                     # t = 24 h opens bin 4
    assert occ[2, 0] == 2 and occ[2, 1] == 1                     # bins are [6k, 6k + 6)


def test_lab_release_offsets_change_availability():
    co = mk_cohort([100.0])
    tc = np.array([23.0])
    F = M.event_features(co, np.array([0]), tc + 1.0, alt={0.0: tc + 0.0, 2.0: tc + 2.0})
    assert not F["avail"][1.0][("full", 24)][0]                  # released at exactly 24 h -> not available at 24 h
    assert F["avail"][0.0][("full", 24)][0]
    assert not F["avail"][2.0][("full", 24)][0]


def test_pre_icu_windows_half_open():
    co = mk_cohort([50.0, 50.0, 50.0])
    F = M.event_features(co, np.array([0, 1, 2]), np.array([-6.0, -6.01, 0.0]))
    assert F["pre6"].tolist() == [True, False, False]
    assert F["pre24"].tolist() == [True, True, False]


def test_last24_window_and_last_age():
    co = mk_cohort([100.0, 100.0])
    F = M.event_features(co, np.array([0, 1]), np.array([47.0, 23.99]))
    assert F["avail"][1.0][("last24", 48)].tolist() == [True, False]
    assert F["last_age"][48][0] == pytest.approx(1.0)


def test_gap_bins_and_interval_crossing():
    co = mk_cohort([100.0], gaps={0: [(20.0, 37.0)]})
    assert co.gap_bins[0].tolist()[:7] == [False, False, False, False, True, True, False]  # bins 24-30, 30-36 in gap
    F = M.event_features(co, np.array([0, 0, 0]), np.array([10.0, 19.0, 40.0]))
    assert F["int_cross"].tolist() == [False, True]


def test_window_cut_at_336h():
    co = mk_cohort([400.0])
    F = M.event_features(co, np.array([0, 0]), np.array([335.9, 336.0]))
    assert F["cnt"][0] == 1


def test_duplicate_same_time_rows_counted_as_zero_intervals():
    co = mk_cohort([100.0])
    F = M.event_features(co, np.array([0, 0, 0]), np.array([5.0, 5.0, 7.0]))
    assert sorted(F["int_d"].tolist()) == [0.0, 2.0]
    assert F["int_med"][0] == 2.0  # zero intervals excluded from the per-encounter median


# ============================================================================= medications (S3)
def courses(rows, cutoff=24.0):
    d = pd.DataFrame(rows, columns=["enc_i", "key", "S", "E"]).sort_values(["enc_i", "key", "S"], kind="mergesort")
    return M.build_courses(d["enc_i"].to_numpy(), d["key"].to_numpy(), d["S"].to_numpy(float), d["E"].to_numpy(float),
                           np.full(len(d), cutoff))


def test_single_dose_course_and_activity():
    cr = courses([(0, "a", 10.0, np.nan)])
    assert len(cr) == 1 and cr["stop"].iloc[0] == 34.0
    assert M.active_at(cr, 34.0).tolist() == [True] and M.active_at(cr, 34.01).tolist() == [False]
    assert M.active_at(cr, 9.99).tolist() == [False]


def test_gap_within_and_beyond_cutoff():
    cr = courses([(0, "a", 0.0, np.nan), (0, "a", 24.0, np.nan), (0, "a", 48.01, np.nan)])
    assert cr["n_docs"].tolist() == [2, 1]


def test_enddate_before_next_doc_ends_course_and_after_is_ignored():
    cr = courses([(0, "a", 0.0, 5.0), (0, "a", 10.0, np.nan)])
    assert cr["n_docs"].tolist() == [1, 1] and cr["stop"].iloc[0] == 5.0 and cr["end_type"].iloc[0] == "enddate"
    cr = courses([(0, "a", 0.0, 100.0), (0, "a", 10.0, np.nan)])
    assert cr["n_docs"].tolist() == [2] and cr["stop"].iloc[0] == 34.0 and cr["end_type"].iloc[0] == "cutoff"


def test_enddate_at_or_before_start_ignored():
    cr = courses([(0, "a", 10.0, 10.0), (0, "b", 10.0, 9.0)])
    assert np.isinf(cr["endE"]).all() and (cr["stop"] == 34.0).all()


def test_as_of_48h_identity():
    rng = np.random.default_rng(3)
    n = 60
    rows = []
    for e in range(n):
        for _ in range(int(rng.integers(1, 6))):
            s = float(rng.uniform(-20, 120))
            E = float(s + rng.uniform(-5, 60)) if rng.random() < 0.4 else np.nan
            rows.append((e, f"k{int(rng.integers(0, 3))}", s, E, rng.choice(["vasoactive", "antibiotic"])))
    docs = pd.DataFrame(rows, columns=["enc_i", "key", "S", "E", "group"]).assign(subgroup="")
    assert M.as_of_identity_check(docs, n) == 0


def test_osmotic_resolution_by_admitting_dx(maps):
    docs = pd.DataFrame({"enc_i": [0, 1, 2, 2], "key": ["20%甘露醇注射液 | 20% Mannitol Injection"] * 2 +
                         ["10%氯化钠注射液 | 10% Sodium Chloride Injection", "头孢曲松 | Ceftriaxone"]})
    nm = M.neuro_matcher(maps["neuro"])
    neuro = np.array([nm(M.norm_icd(c)[0]) for c in ("S06.001", "Q28.000", "Q28.201")])
    r = M.resolve_drug_groups(docs, maps["drug"], neuro)
    assert r["group"].tolist() == ["osmotherapy", "diuretic", "osmotherapy", "antibiotic"]
    assert r["group_sens"].tolist() == ["osmotherapy", "osmotherapy", "osmotherapy", "antibiotic"]


def test_vasoactive_concurrent_count():
    cr = courses([(0, "dopa", 0.0, np.nan), (0, "epi", 5.0, np.nan), (0, "mil", 30.0, np.nan), (1, "dopa", 0.0, np.nan)],
                 cutoff=12.0)
    mc = M.max_concurrent(cr, 3, np.zeros(3), np.full(3, 24.0))
    assert mc.tolist() == [2, 1, 0]


# ============================================================================= microbiology
def test_base_code_collections_and_susceptibility_only_with_addons(maps):
    ev = pd.DataFrame({"enc_i": [0, 0, 1, 2], "th": [5.0, 5.0, 5.0, 5.0],
                       "spec_key": ["LIS0001 | 血培养", "LIS0001+LIS0204 | 血培养,常规药敏定量试验(MIC)", "LIS0001 | 血培养",
                                    "LIS0002 | 尿培养"],
                       "spec_itemid": ["LIS0001", "LIS0001+LIS0204", "LIS0001", "LIS0002"],
                       "org_key": ["MIC1001 | 金黄色葡萄球菌", "MIC1001 | 金黄色葡萄球菌", "MIC1001 | 金黄色葡萄球菌", "MIC1002 | 无细菌生长"]})
    col, info = M.micro_collections(ev, maps["specimen"], maps["organism"])
    assert len(col) == 3                                       # the add-on row joins its base collection
    c0 = col[col["enc_i"] == 0].iloc[0]
    assert c0["base"] == "LIS0001" and c0["susc"] and c0["ident"] and c0["gram"]
    c1 = col[col["enc_i"] == 1].iloc[0]
    assert not c1["susc"] and c1["ident"]                      # no add-on key -> no susceptibility result
    c2 = col[col["enc_i"] == 2].iloc[0]
    assert c2["neg"] and not c2["gram"] and c2["group"] == "urine"
    assert set(col["base"]) == {"LIS0001", "LIS0002"}          # never the '+' add-on code


# ============================================================================= disclosure
def test_qdict_requires_ten_values():
    assert M.qdict(np.arange(9)) is None
    assert M.qdict(np.arange(10))["p50"] == 4.5


def test_dob_breakdown_suppression_on_small_invalid_counts():
    """Small invalid-status cells are <10; the large valid complement is complementary-S; total stays published."""
    out = M.Out(A)
    tot = out.cell(("denominators", "n_encounters"), 100)
    out.frozen.add(tot)
    sc = out.partition(("layer_a", "dob_breakdown", "status"),
                       {"missing": 3, "unparseable": 2, "negative": 1, "valid": 94}, total=tot)
    n_inv = out.cell(("layer_a", "dob_breakdown", "n_invalid"), 6)
    out.ident(n_inv, [sc["missing"], sc["unparseable"], sc["negative"]])
    st, supp = M.apply_suppression_frozen(A, out.reg, out.frozen)
    assert st[tot] == "pub" and st[sc["valid"]] == "S"
    for k in ("missing", "unparseable", "negative"):
        assert st[sc[k]] == "P"
    assert st[n_inv] == "P" and supp["recoverable"] == 0


def test_frozen_totals_never_complementary_suppressed():
    out = M.Out(A)
    tot = out.cell(("x", "n_total"), 40)
    out.frozen.add(tot)
    a = out.cell(("x", "n_a"), 3)
    b = out.cell(("x", "n_b"), 37)
    out.ident(tot, [a, b])
    st, supp = M.apply_suppression_frozen(A, out.reg, out.frozen)
    assert st[a] == "P" and st[b] == "S" and st[tot] == "pub" and supp["recoverable"] == 0


def test_frozen_total_under_ten_is_an_error():
    out = M.Out(A)
    tot = out.cell(("x", "n_total"), 4)
    out.frozen.add(tot)
    with pytest.raises(M.ReleaseGateError):
        M.apply_suppression_frozen(A, out.reg, out.frozen)


def test_count_key_prefix_enforced():
    out = M.Out(A)
    with pytest.raises(M.AuditError):
        out.cell(("x", "count"), 5)
    with pytest.raises(M.AuditError):
        out.val(("x", "n_label"), "text")


def test_cli_outcome_flag_default_off():
    # Path flags are required (no local filesystem default is shipped).
    base = ["--project-root", "/tmp/x", "--raw-dir", "/tmp/x/raw", "--split-file", "/tmp/x/splits.csv"]
    assert M.build_config(base).outcome_descriptors is False
    assert M.build_config(base + ["--outcome-descriptors"]).outcome_descriptors is True


# ============================================================================= end-to-end (synthetic project)
@pytest.fixture(scope="module")
def project(tmp_path_factory):
    root = tmp_path_factory.mktemp("s2")
    return S.make_project(root)


def _run(project, out_dir, **kw):
    cfg = S.make_config(M, project, out_dir, V17, **kw)
    M.run_audit(cfg, sentinels=[S.SENTINEL_SUBJECT])
    return cfg, json.loads(cfg.shared_json.read_text(encoding="utf-8")), cfg.shared_md.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def run_off(project, tmp_path_factory):
    return _run(project, tmp_path_factory.mktemp("out_off"))


@pytest.fixture(scope="module")
def run_on(project, tmp_path_factory):
    return _run(project, tmp_path_factory.mktemp("out_on"), outcome_descriptors=True)


def _published_counts(doc):
    res = []

    def walk(x, p):
        if isinstance(x, dict):
            for k, v in x.items():
                if p and p[0] in M.NON_COUNT_TOP:
                    continue
                walk(v, p + (k,))
        elif p and M.COUNT_KEY_RE.match(str(p[-1])):
            res.append((p, x))
    walk(doc, ())
    return res


def test_e2e_release_gate_and_files(run_off):
    cfg, doc, md = run_off
    assert doc["release_gate"] == "ok" and doc["outcome_descriptors_flag"] is False
    assert cfg.local_json.exists() and cfg.local_json.name.endswith("_raw_LOCAL.json")
    shared = sorted(p.name for p in cfg.out_dir.iterdir())
    assert shared == sorted([cfg.shared_json.name, cfg.shared_md.name, cfg.local_json.name])
    assert M.release_gate(A, cfg.shared_json, cfg.shared_md, "ok", [S.SENTINEL_SUBJECT]) == []
    assert "outcome_descriptors" not in doc and "outcome_denominators" not in doc


def test_e2e_no_count_under_ten_published(run_off):
    _, doc, _ = run_off
    for p, v in _published_counts(doc):
        assert v in ("<10", "S") or (isinstance(v, int) and v >= 10), p


def test_e2e_no_identifier_or_timestamp_in_shared(run_off, project):
    cfg, doc, md = run_off
    txt = cfg.shared_json.read_text(encoding="utf-8") + md
    assert S.SENTINEL_SUBJECT not in txt
    enc = project["enc"]
    for s in enc["encounter_id"].head(50):
        assert s not in txt
    assert not M.TIMESTAMP_RE.search(txt)


def test_e2e_split_totals_published(run_off, project):
    _, doc, _ = run_off
    for s, n in project["split_counts"].items():
        assert doc["denominators"][s]["n_encounters"] == n


def test_e2e_sodium_ever_matches_slow_reference(run_off, project):
    """Independent slow reference: train encounters with a gated sodium value released (+1 h) in [t0, W)."""
    _, doc, _ = run_off
    enc = project["enc"]
    lab = pd.read_csv(project["raw"] / "LABEVENTS.csv", dtype=str, keep_default_na=False)
    lab = lab[lab["ITEMID"] == "5230"]
    by = {(r.SUBJECT_ID, r.HADM_ID): r for r in enc.itertuples()}
    hit = set()
    for _, row in lab.iterrows():
        e = by.get((row["SUBJECT_ID"], row["HADM_ID"]))
        if e is None or e.split != "train":
            continue
        v = float(row["VALUENUM"]) if row["VALUENUM"] else math.nan
        if not (100 <= v <= 190):
            continue
        t = (pd.Timestamp(row["CHARTTIME"]) - e.t0).total_seconds() / 3600 + 1.0
        W = min((e.OUT - e.t0).total_seconds() / 3600, 336.0)
        if 0 <= t < W:
            hit.add(e.encounter_id)
    got = doc["channels"]["sodium"]["train"]["n_ever_336h"]
    assert got == (len(hit) if len(hit) >= 10 else "<10") or got == "S"


def test_e2e_reconciliation_classes_sum(run_off, project):
    _, doc, _ = run_off
    rec = doc["reconciliation"]["CHARTEVENTS"]
    assert rec["n_total"] == project["raw_rows"]["CHARTEVENTS"]


def test_e2e_dob_breakdown_present_and_aggregate(run_off):
    _, doc, md = run_off
    db = doc["layer_a"]["dob_breakdown"]
    assert "status" in db and "negative_bins" in db and "by_split" in db and "by_first_careunit" in db
    assert "by_split_and_first_careunit" in db and "weight_imputation" in db
    assert "n_invalid" in db and "n_age_days" in db and "q_age_days" in db
    assert "age_invalid" in db["weight_imputation"]
    assert "safety_gate" in db and "1%" in db["safety_gate"]
    assert "sex" in db["weight_imputation_rule"] and "age_invalid" in db["weight_imputation_rule"]
    assert "## Layer A DOB breakdown (D20)" in md
    # default synthetic patients have valid DOBs; invalid cells are small and must not appear as raw 0-9
    for k in ("n_missing", "n_unparseable", "n_negative", "n_valid"):
        v = db["status"][k]
        assert v in ("<10", "S") or (isinstance(v, int) and v >= 10)


def test_e2e_map_checks_present(run_off):
    _, doc, _ = run_off
    mc = doc["map_checks"]
    for k in ("row_status", "units", "duplicates", "magnesium_ug_dl_rows", "heart_rate", "map", "pf_ratio", "sf_ratio",
              "oxygenation_coverage", "charting_interval", "prescription_start_check", "micro_fill_rates"):
        assert k in mc, k
    assert set(mc["units"]["LABEVENTS_5257"]) >= {"n_g/dl", "n_g/l"}
    assert mc["heart_rate"]["rule"].startswith("1003 kept")
    assert "ORG_NAME" in mc["micro_fill_rates"]["columns"]
    assert "n_in_stay_date_only" in mc["prescription_start_check"]
    assert "n_sf_available_without_pf" in mc["oxygenation_coverage"]["train"]["72h"]


def test_e2e_charting_interval_never_imputes_zero_days(run_off):
    _, doc, _ = run_off
    ci = doc["map_checks"]["charting_interval"]["urine_output"]
    assert "n_complete_days_without_entry" in ci
    # daily totals are computed only over days with at least one entry
    n_days = ci["n_complete_days_train"]
    n_tot = ci.get("n_daily_total_ml_per_kg_per_day") or ci.get("n_daily_total_ml")
    if isinstance(n_days, int) and isinstance(n_tot, int) and isinstance(ci["n_complete_days_without_entry"], int):
        assert n_tot == n_days - ci["n_complete_days_without_entry"]


def test_e2e_osmotic_and_vasoactive_sections(run_off):
    _, doc, _ = run_off
    assert set(doc["osmotic_resolution"]["train"]) == {"primary", "sensitivity_all_osmotherapy"}
    assert "subgroup_sensitivity" in doc["vasoactive_count"]["train"]["24h"]
    assert doc["notes"]["as_of_48h_identity"] == "ok"


def test_e2e_micro_m3_grid(run_off):
    _, doc, _ = run_off
    m3 = doc["microbiology"]["train"]["m3"]["72h"]
    assert set(m3) == {"grid_x0.5", "grid_x1", "grid_x1.5"}
    assert set(m3["grid_x1"]) >= {"gram_class", "identification", "susceptibility", "negative_final"}


def test_e2e_outcome_descriptors_only_with_flag(run_on):
    cfg, doc, md = run_on
    assert doc["outcome_descriptors_flag"] is True and doc["release_gate"] == "ok"
    assert "outcome_denominators" in doc and "outcome_descriptors" in doc
    assert "outcome1" in doc["outcome_descriptors"]["sodium"]
    assert M.release_gate(A, cfg.shared_json, cfg.shared_md, "ok", [S.SENTINEL_SUBJECT]) == []


# ----------------------------------------------------------------------------- release gate tampering
def _tamper(run_off, tmp_path, fn):
    cfg, doc, md = run_off
    j, m = tmp_path / "a.json", tmp_path / "a.md"
    d = json.loads(json.dumps(doc))
    md2 = fn(d, md)
    j.write_text(json.dumps(d, ensure_ascii=False), encoding="utf-8")
    m.write_text(md2 if md2 is not None else M.render_markdown(d), encoding="utf-8")
    return M.release_gate(A, j, m, "ok", [S.SENTINEL_SUBJECT])


def test_gate_catches_small_count(run_off, tmp_path):
    def f(d, md):
        d["channels"]["sodium"]["train"]["n_ever_336h"] = 5
    assert _tamper(run_off, tmp_path, f)


def test_gate_catches_unsuppressed_cell(run_off, tmp_path):
    def f(d, md):
        p = next(p for p, v in _published_counts(d) if v == "<10")
        x = d
        for k in p[:-1]:
            x = x[k]
        x[p[-1]] = 12
    assert _tamper(run_off, tmp_path, f)


def test_gate_catches_timestamp_and_sentinel(run_off, tmp_path):
    assert _tamper(run_off, tmp_path, lambda d, md: d["notes"].__setitem__("x", "2016-01-01 06:00:00"))
    assert _tamper(run_off, tmp_path, lambda d, md: d["notes"].__setitem__("x", S.SENTINEL_SUBJECT))


def test_gate_catches_markdown_mismatch(run_off, tmp_path):
    assert _tamper(run_off, tmp_path, lambda d, md: md + "\nextra")


def test_gate_catches_outcome_keys_with_flag_off(run_off, tmp_path):
    assert _tamper(run_off, tmp_path, lambda d, md: d.__setitem__("outcome_denominators", {}))


def test_gate_catches_quantile_for_small_group(run_off, tmp_path):
    def f(d, md):
        tr = d["channels"]["sodium"]["train"]
        tr["n_first_hours"] = "<10"
    assert _tamper(run_off, tmp_path, f)


# ----------------------------------------------------------------------------- hash and count gates
def _copy_project(project, tmp_path):
    root = tmp_path / "p"
    shutil.copytree(project["root"], root)
    p = dict(project)
    p.update(root=root, raw=root / "raw", interim=root / "data" / "interim")
    return p


@pytest.mark.parametrize("target", ["canonical", "map", "raw", "complexity"])
def test_hash_mismatch_stops(project, tmp_path, target):
    p = _copy_project(project, tmp_path)
    f = {"canonical": p["interim"] / "canonical_icu_encounters.csv", "map": p["interim"] / "maps" / "drug_group_map_v1.csv",
         "raw": p["raw"] / "CHARTEVENTS.csv", "complexity": p["root"] / "config" / "surgery_complexity_map_v1.csv"}[target]
    with open(f, "a", encoding="utf-8") as fh:
        fh.write("\n")
    cfg = S.make_config(M, p, tmp_path / "out", V17)
    with pytest.raises(M.HashMismatchError):
        M.run_audit(cfg)
    assert not cfg.shared_json.exists()


def test_v17_module_hash_mismatch_stops(tmp_path):
    bad = tmp_path / M.V17_MODULE_FILE
    bad.write_text(V17.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(M.HashMismatchError):
        M.load_v17(bad)


def test_record_only_table_row_count_gated(project, tmp_path):
    p = _copy_project(project, tmp_path)
    cfg = S.make_config(M, p, tmp_path / "out", V17)
    cfg.expected_raw_rows["ADMISSIONS"] += 1
    with pytest.raises(M.AuditError):
        M.run_audit(cfg)
    assert not cfg.shared_json.exists()


def test_main_exit_code_2_on_bad_input(tmp_path):
    rc = M.main(["--project-root", str(tmp_path), "--raw-dir", str(tmp_path / "raw"), "--split-file", str(tmp_path / "s.csv"),
                 "--v17-module", str(V17)])
    assert rc == 2


# ============================================================================= added coverage (spec 13)
def test_link_requires_admission_in_hadm_list_for_pre_icu():
    co = mk_cohort([50.0, 50.0])
    t0 = co.enc["t0"].iloc[0]
    subj = pd.Series(["100", "100", "100", "101"])
    hadm = pd.Series(["500", "999", "500", "501"])            # row 1: same subject, other admission -> not linked
    t = pd.Series([t0 - pd.Timedelta(hours=3), t0 - pd.Timedelta(hours=3), pd.NaT, t0 + pd.Timedelta(hours=60)])
    classes = {c: 0 for c in M.ROW_CLASSES}
    k = M.link_chunk(subj, hadm, t, co.eh, 0.0, M.KEEP_LO_H, classes)
    assert k["r"].tolist() == [0]
    assert classes["no_cohort_admission"] == 1 and classes["cohort_admission_time_missing"] == 1
    assert classes["cohort_admission_outside_stay"] == 2 and sum(classes.values()) == 4


def test_lab_release_offset_moves_row_class():
    co = mk_cohort([10.0])
    t0 = co.enc["t0"].iloc[0]
    for off, cls in ((0.0, "in_stay_train"), (1.0, "cohort_admission_outside_stay")):
        classes = {c: 0 for c in M.ROW_CLASSES}
        M.link_chunk(pd.Series(["100"]), pd.Series(["500"]), pd.Series([t0 + pd.Timedelta(hours=9.5)]), co.eh, off,
                     M.KEEP_LO_H, classes)
        assert classes[cls] == 1


def test_merge_group_keeps_source_itemid(maps):
    lab = ev_df([(0, 1, "5041", 50, "μmol/l"), (0, 2, "5032", 60, "μmol/l")])
    vb = M.build_values(ev_df([]), lab, maps["channel"])
    assert vb.events["creatinine"]["item"].tolist() == ["5041", "5032"]


def test_empty_bins_are_never_zero_values():
    """An empty urine-output bin stays empty (no row, no value): occupancy 0 and no last value, never a 0 mL entry."""
    co = mk_cohort([48.0, 48.0])
    F = M.event_features(co, np.array([0]), np.array([3.0]))
    assert F["occ"][1].sum() == 0 and not F["ever_w"][1] and np.isnan(F["last_age"][24][1])
    assert F["occ"][0, 0] == 1 and F["occ"][0, 1:8].sum() == 0


def test_last_obs_age_against_locf_cap():
    co = mk_cohort([100.0, 100.0])
    F = M.event_features(co, np.array([0, 1]), np.array([66.0, 47.9]))
    ages = F["last_age"][72]
    fam, cap = M.locf_family("sodium", False)
    assert ages.tolist() == pytest.approx([6.0, 24.1]) and ((ages <= cap).tolist() == [True, False])
    fam, cap = M.locf_family("creatinine", False)
    assert (ages <= cap).tolist() == [True, False]


def test_pair_block_co_measurement():
    out = M.Out(A)
    X = np.zeros((40, 3), dtype=bool)
    X[:30, 0] = True
    X[10:35, 1] = True
    X[:2, 2] = True
    M._pair_block(out, ("cm",), X, ["a", "b", "c"], "draws")
    r = out.reg.raw
    assert r[A.P("cm", "n_draws")] == 40
    assert r[A.P("cm", "pairs", "a__b", "n_draws_both")] == 20
    dep = {p: v for p, _, v in out.kss}
    assert dep[A.P("cm", "pairs", "a__b", "x_jaccard")] == round(20 / 35, 3)
    rat = {p: (n, d) for p, n, d, _ in out.ratios}
    assert rat[A.P("cm", "pairs", "a__b", "pct_b_given_a")] == (A.P("cm", "pairs", "a__b", "n_draws_both"),
                                                                  A.P("cm", "channels", "a", "n_draws_with_channel"))


def test_e2e_comeasurement_sections(run_off):
    _, doc, _ = run_off
    cmx = doc["comeasurement"]["train"]
    assert "lab_draws" in cmx and "vitals_same_charttime" in cmx and "bin_cooccupancy_bins_0_11" in cmx
    v = cmx["vitals_same_charttime"]
    assert "heart_rate__respiratory_rate" in v["pairs"]


def test_e2e_metric10_strata(run_off):
    _, doc, _ = run_off
    st = doc["drug_groups"]["antibiotic"]["train"]["strata"]
    assert set(st) == {"age_band", "surgical_status", "first_careunit"}


def test_e2e_surgery_link_classes(run_off, project):
    _, doc, _ = run_off
    lc = doc["surgery"]["link_classes"]
    assert any(k.startswith("n_") for k in json.dumps(lc).split('"'))
    assert "layer_c" in doc["surgery"] and "layer_a" in doc["surgery"]


def test_gate_failure_messages_are_specific(run_off, tmp_path):
    def f(d, md):
        d["channels"]["sodium"]["train"]["n_ever_336h"] = 5
    fails = _tamper(run_off, tmp_path, f)
    assert any("below 10" in x or "identity" in x or "md" in x for x in fails), fails
