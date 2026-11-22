"""pytest suite for audit_surgery_complexity_duration.py (spec v1.7 + fix list incl. Amendments A, B and C).

SYNTHETIC DATA ONLY (synthetic_data.py). No real PIC file is read. Synthetic inputs get their own SHA-256 values
through AuditConfig.pinned_sha256 (test override), so the hash check is still exercised on every run.

Run:  python -m pytest -q test_audit_surgery_complexity_duration.py
"""
from __future__ import annotations

import copy
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit_surgery_complexity_duration as A  # noqa: E402
import synthetic_data as S  # noqa: E402


# =============================================================================
# Fixtures
# =============================================================================
@pytest.fixture(scope="session")
def project(tmp_path_factory):
    root = tmp_path_factory.mktemp("syn")
    return S.make_project(root, n_subjects=600, seed=7)


@pytest.fixture(scope="session")
def e2e(project, tmp_path_factory):
    """One full end-to-end synthetic run (hash checks -> linkage -> stats -> suppression -> release gate)."""
    out = tmp_path_factory.mktemp("out")
    cfg = S.make_config(project, out, bootstrap_reps=30, verbose=False)
    doc = A.run_audit(cfg, sentinels=[S.SENTINEL_SUBJECT])
    return {"cfg": cfg, "doc": doc, "out": out}


@pytest.fixture(scope="session")
def frame(project, tmp_path_factory):
    cfg = S.make_config(project, tmp_path_factory.mktemp("f"), bootstrap_reps=0, verbose=False)
    inp = A.load_inputs(cfg)
    f, meta = A.build_frame(inp)
    return f, meta, inp, cfg


def _analysis(proj, tmp, reps=0):
    cfg = S.make_config(proj, tmp, bootstrap_reps=reps, verbose=False)
    inp = A.load_inputs(cfg)
    f, meta = A.build_frame(inp)
    an = A.Analysis(cfg, f, meta, inp["info"])
    an.run()
    return an, f


# =============================================================================
# Parsing (known defects: double parsing, json.loads on repr, IDs not str)
# =============================================================================
@pytest.mark.parametrize("v,exp", [
    ("[1, 2]", ["1", "2"]), ("['1', '2']", ["1", "2"]), ("[1.0, 2.0]", ["1", "2"]), ("1;2", ["1", "2"]),
    ("1,2", ["1", "2"]), ("1|2", ["1", "2"]), ("123", ["123"]), ("123.0", ["123"]), ("", []), ("nan", []),
    ([1, 2], ["1", "2"]), (("3", 3.0), ["3"]), ("[]", []), ("(5, 6)", ["5", "6"]), ("[ '7' , \"8\" ]", ["7", "8"]),
])
def test_parse_hadm_list(v, exp):
    assert A.parse_hadm_list(v) == exp


def test_parse_hadm_list_is_idempotent_on_parsed_lists():
    once = A.parse_hadm_list("['10', '11']")
    assert A.parse_hadm_list(once) == once  # a parsed list passed again is not crashed on (pd.isna on list defect)


@pytest.mark.parametrize("v,exp", [("123", "123"), ("123.0", "123"), (123, "123"), (123.0, "123"), (" 123 ", "123"),
                                   ("", ""), (float("nan"), ""), (None, ""), ("'45'", "45"), ("0", "0")])
def test_norm_id(v, exp):
    assert A.norm_id(v) == exp


def test_csv_read_as_str(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text("SUBJECT_ID,HADM_ID\n123.0,0045\n,NA\n", encoding="utf-8")
    df = A.read_csv_str(p)
    assert df["SUBJECT_ID"].tolist() == ["123.0", ""] and df["HADM_ID"].tolist() == ["0045", "NA"]


def test_name_normalisation():
    assert A.norm_name("  ASD   Closure ") == "asd closure"


# =============================================================================
# Hashes (four pinned, checked before reading)
# =============================================================================
def test_pinned_hashes_are_full_and_four():
    assert set(A.PINNED_SHA256) == {"canonical", "harmonization", "complexity_map", "split"}
    assert all(len(h) == 64 and all(c in "0123456789abcdef" for c in h) for h in A.PINNED_SHA256.values())


@pytest.mark.parametrize("which", ["canonical", "harmonization", "complexity_map", "split"])
def test_hash_mismatch_fails(project, tmp_path, which):
    cfg = S.make_config(project, tmp_path, bootstrap_reps=0, verbose=False)
    bad = dict(cfg.pinned_sha256)
    bad[which] = "0" * 64
    cfg.pinned_sha256 = bad
    with pytest.raises(A.HashMismatchError):
        A.load_inputs(cfg)


def test_hash_checked_before_reading(project, tmp_path):
    """A modified (even unparseable) input fails on the hash, never on parsing."""
    proj = copy.deepcopy(project)
    dst = tmp_path / "canon.csv"
    dst.write_bytes(b"\x00\x01 not a csv at all")
    proj["paths"]["canonical"] = dst
    cfg = S.make_config(proj, tmp_path, bootstrap_reps=0, verbose=False)
    with pytest.raises(A.HashMismatchError):
        A.load_inputs(cfg)


def test_short_pinned_hash_rejected(project, tmp_path):
    cfg = S.make_config(project, tmp_path, bootstrap_reps=0, verbose=False)
    cfg.pinned_sha256 = dict(cfg.pinned_sha256, split=cfg.pinned_sha256["split"][:12])
    with pytest.raises(A.HashMismatchError):
        A.load_inputs(cfg)


def test_production_hashes_are_the_default():
    # Path flags are required (no local filesystem default is shipped).
    cfg = A.build_config([
        "--project-root", "/nonexistent",
        "--raw-dir", "/nonexistent/raw",
        "--split-file", "/nonexistent/splits.csv",
    ])
    assert cfg.pinned_sha256 == A.PINNED_SHA256
    assert cfg.bootstrap_reps == 2000 and cfg.expected_surgery_rows is None


# =============================================================================
# Input validation
# =============================================================================
def test_missing_harmonization_column_fails(project, tmp_path):
    proj = copy.deepcopy(project)
    h = pd.read_csv(proj["paths"]["harmonization"], dtype=str).rename(columns={"t0": "something_else"})
    p = tmp_path / "harm.csv"
    h.to_csv(p, index=False)
    proj["paths"]["harmonization"] = p
    proj["hashes"]["harmonization"] = A.sha256_file(p)
    cfg = S.make_config(proj, tmp_path, bootstrap_reps=0, verbose=False)
    with pytest.raises(A.AuditError, match="harmonization T0"):
        A.load_inputs(cfg)
    cfg.harm_t0_col = "something_else"  # explicit override works
    assert A.load_inputs(cfg)["info"]["checks"]["harmonization_columns"]["t0"] == "something_else"


def test_wrong_expected_rows_fails(project, tmp_path):
    cfg = S.make_config(project, tmp_path, bootstrap_reps=0, verbose=False)
    cfg.expected_surgery_rows += 1
    with pytest.raises(A.AuditError, match="rows"):
        A.load_inputs(cfg)


def test_wrong_split_totals_fail(project, tmp_path):
    cfg = S.make_config(project, tmp_path, bootstrap_reps=0, verbose=False)
    cfg.expected_split_encounters = dict(cfg.expected_split_encounters, train=cfg.expected_split_encounters["train"] + 1)
    with pytest.raises(A.AuditError, match="cohort/split"):
        A.load_inputs(cfg)


def test_criteria_config_validated(monkeypatch):
    A.validate_criteria()
    monkeypatch.setattr(A, "TAU_B_SIGNED_MAX", 0.3)
    with pytest.raises(A.AuditError):
        A.validate_criteria()


# =============================================================================
# Linkage (one surgery, one encounter)
# =============================================================================
def _enc(rows):
    df = pd.DataFrame(rows, columns=["encounter_id", "SUBJECT_ID", "hadm_list", "t0", "OUT", "split"])
    df["t0"] = pd.to_datetime(df["t0"]).astype("datetime64[ns]")
    df["OUT"] = pd.to_datetime(df["OUT"]).astype("datetime64[ns]")
    return df


def _surg(rows):
    """T = timing time (anaesthesia end; v1.7)."""
    df = pd.DataFrame(rows, columns=["sid", "SUBJECT_ID", "HADM_ID", "T"])
    df["T"] = pd.to_datetime(df["T"]).astype("datetime64[ns]")
    return df


T = pd.Timestamp("2020-01-10 00:00:00")


def _two_encounters_same_hadm():
    return _enc([["e1", "1", ["10"], T, T + pd.Timedelta(hours=48), "train"],
                 ["e2", "1", ["10"], T + pd.Timedelta(hours=100), T + pd.Timedelta(hours=150), "train"]])


@pytest.mark.parametrize("offset_h,cls,timing,enc", [
    # offsets are anaesthesia END relative to e1 t0 (e1: t0..t0+48 h; e2: t0+100..t0+150 h; same HADM)
    (-24.0, "linked", "recent_layer_a", "e1"),          # exactly t0 - 24 h -> recent (Layer A)
    (-24.01, "linked", "remote_same_hadm", "e1"),       # just over 24 h before t0 -> remote
    (0.0, "linked", "recent_layer_a", "e1"),            # ends at t0 -> recent (old B = t0 edge is moot)
    (6.0, "linked", "recent_layer_a", "e1"),            # ends exactly t0 + 6 h -> still Layer A only
    (6.01, "linked", "during_icu_layer_c", "e1"),       # just after t0 + 6 h -> Layer C
    (48.0, "linked", "during_icu_layer_c", "e1"),       # ends at OUT -> Layer C
    (80.0, "linked", "recent_layer_a", "e2"),           # between encounters, 20 h before e2 -> later encounter, recent
    (60.0, "linked", "remote_same_hadm", "e2"),         # between encounters, 40 h before e2 -> later encounter, remote
    (49.0, "linked", "remote_same_hadm", "e2"),         # 1 h after e1 OUT -> later encounter (not e1)
    (103.0, "linked", "recent_layer_a", "e2"),          # 3 h after e2 t0 -> recent for e2
    (120.0, "linked", "during_icu_layer_c", "e2"),      # Layer C of the second encounter
    (151.0, "post_icu_same_hadm", "", ""),              # after the last OUT -> rule 3
])
def test_linkage_rules(offset_h, cls, timing, enc):
    s = _surg([[0, "1", "10", T + pd.Timedelta(hours=offset_h)]])
    r = A.link_surgeries(s, _two_encounters_same_hadm()).iloc[0]
    assert r["link_class"] == cls and r["timing"] == timing and r["encounter_id"] == enc


def test_linkage_one_surgery_one_encounter_and_classes():
    enc = _two_encounters_same_hadm()
    enc = pd.concat([enc, _enc([["e3", "2", ["20", "21"], T, T + pd.Timedelta(hours=10), "test"],
                                ["e4", "3", ["30"], pd.NaT, pd.NaT, "train"]])], ignore_index=True)
    s = _surg([[0, "1", "10", T + pd.Timedelta(hours=1)],      # during e1 (both e1/e2 candidates)
               [1, "1", "99", T],                              # other HADM (earlier/later split by first t0)
               [2, "77", "10", T],                             # not in cohort
               [3, "1", "10", pd.NaT],                         # invalid time
               [4, "2", "21", T + pd.Timedelta(hours=2)],      # second HADM in a multi-HADM list
               [5, "3", "30", T]])                             # encounter has no harmonized times
    r = A.link_surgeries(s, enc).set_index("sid")
    assert len(r) == 6 and r.index.is_unique
    assert r.loc[0, "encounter_id"] == "e1" and r.loc[0, "n_candidates"] == 2
    assert r.loc[1, "link_class"] == "other_hadm" and r.loc[1, "other_hadm_side"] == "later_hadm"
    assert r.loc[2, "link_class"] == "not_in_cohort"
    assert r.loc[3, "link_class"] == "invalid_time"
    assert r.loc[4, "encounter_id"] == "e3" and r.loc[4, "split"] == "test"
    assert r.loc[5, "link_class"] == "not_linkable_missing_encounter_time"


def test_recent_is_layer_a_only_never_layer_c(frame):
    """v1.7: surgeries ending in [t0 - 24 h, t0 + 6 h] are recent_layer_a and never also during_icu_layer_c."""
    f, _, _, _ = frame
    lk = f[f["link_class"] == "linked"]
    assert lk["ROW_ID"].is_unique
    after = (lk["T"] - lk["t0"]).dt.total_seconds() / 3600
    rec = lk["timing"] == "recent_layer_a"
    assert rec.any() and ((after >= -24) & (after <= 6))[rec].all()
    assert (lk.loc[lk["timing"] == "during_icu_layer_c"].pipe(lambda d: (d["T"] - d["t0"]).dt.total_seconds() / 3600) > 6).all()
    assert (lk.loc[lk["timing"] == "remote_same_hadm"].pipe(lambda d: (d["T"] - d["t0"]).dt.total_seconds() / 3600) < -24).all()
    assert lk["recent_after_t0"].sum() >= 10 and (lk["recent_after_t0"] <= rec).all()
    fam = A.FAMILIES["d_during_icu"][1](lk)
    assert not (fam & rec.to_numpy()).any()


def test_timing_uses_anesthesia_end_with_surgery_end_fallback(frame, project):
    f, _, inp, _ = frame
    raw = inp["surg"]
    ae_blank = raw["ANES_END_TIME"].str.strip() == ""
    e_blank = raw["SURGERY_END_TIME"].str.strip() == ""
    assert (ae_blank & ~e_blank).sum() >= 1 and (ae_blank & e_blank).sum() >= 1
    fb = f.loc[(ae_blank & ~e_blank).to_numpy()]
    assert fb["timing_fallback"].all() and (fb["T"] == fb["E"]).all()
    assert (f.loc[(ae_blank & e_blank).to_numpy(), "link_class"] == "invalid_time").all()
    ok = f.loc[(~ae_blank).to_numpy()]
    assert (ok["T"] == ok["AE"]).all()
    # a missing BEGIN time no longer blocks linkage (timing is by the end); the duration is still invalid
    b_blank = (raw["SURGERY_BEGIN_TIME"].str.strip() == "").to_numpy()
    assert (f.loc[b_blank & (~ae_blank).to_numpy(), "link_class"] != "invalid_time").all()


def test_timing_fallback_count_published(e2e):
    p = e2e["doc"]["population"]["pooled"]
    for k in ("timing_fallback_surgery_end", "recent_ending_after_t0", "pre_admit_time", "straddling_t0"):
        assert p[k] in ("<10", "S") or (isinstance(p[k], int) and p[k] >= 10)
    assert "begin_equals_t0" not in p
    assert set(p["by_timing_class"]) == {"recent_layer_a", "remote_same_hadm", "during_icu_layer_c"}


def test_not_in_cohort_counted_once(frame):
    f, _, _, _ = frame
    assert (f["link_class"] == "not_in_cohort").sum() == 25
    assert f["sid"].is_unique


def test_reconciliation_sums_exactly(e2e, project, frame):
    f, _, _, _ = frame
    counts = f["link_class"].value_counts()
    assert set(counts.index) <= set(A.LINK_CLASSES)
    assert int(counts.sum()) == project["n_surgeries"]
    raw = json.loads(e2e["cfg"].local_json.read_text(encoding="utf-8"))
    assert sum(raw["link_class_counts"].values()) == project["n_surgeries"]
    assert raw["raw_counts"]["/reconciliation/surgery_info_rows"] == project["n_surgeries"]


def test_reconciliation_identity_failure_raises():
    reg = A.Registry()
    t = reg.count("/a", 10)
    with pytest.raises(A.ReconciliationError):
        reg.identity(t, [reg.count("/b", 4), reg.count("/c", 5)])


def test_duplicate_row_id_fails_in_build_frame(frame):
    _, _, inp, _ = frame
    inp2 = dict(inp)
    surg = inp["surg"].copy()
    surg.loc[1, "ROW_ID"] = surg.loc[0, "ROW_ID"]
    inp2["surg"] = surg
    with pytest.raises(A.ReconciliationError):
        A.build_frame(inp2)


def test_blank_row_id_fails_in_build_frame(frame):
    _, _, inp, _ = frame
    inp2 = dict(inp)
    surg = inp["surg"].copy()
    surg.loc[2, "ROW_ID"] = "  "
    inp2["surg"] = surg
    with pytest.raises(A.ReconciliationError):
        A.build_frame(inp2)


def _rewrite_surgery_info(project, tmp_path, edit):
    proj = copy.deepcopy(project)
    si = pd.read_csv(project["paths"]["surgery_info"], dtype=str, keep_default_na=False)
    edit(si)
    dst = tmp_path / "SURGERY_INFO.csv"
    si.to_csv(dst, index=False)
    proj["paths"]["surgery_info"] = dst
    return proj


def test_duplicate_row_id_stops_load_inputs(project, tmp_path):
    def edit(si):
        si.loc[5, "ROW_ID"] = si.loc[4, "ROW_ID"]
    proj = _rewrite_surgery_info(project, tmp_path, edit)
    with pytest.raises(A.AuditError, match="ROW_ID is not unique"):
        A.load_inputs(S.make_config(proj, tmp_path / "o", bootstrap_reps=0, verbose=False))


def test_float_like_row_id_duplicate_detected(project, tmp_path):
    def edit(si):
        si.loc[5, "ROW_ID"] = si.loc[4, "ROW_ID"] + ".0"
    proj = _rewrite_surgery_info(project, tmp_path, edit)
    with pytest.raises(A.AuditError, match="ROW_ID"):
        A.load_inputs(S.make_config(proj, tmp_path / "o", bootstrap_reps=0, verbose=False))


def test_blank_row_id_stops_load_inputs(project, tmp_path):
    def edit(si):
        si.loc[3, "ROW_ID"] = ""
    proj = _rewrite_surgery_info(project, tmp_path, edit)
    proj["missing"]["ROW_ID"] = 1  # get past the missingness check: the ROW_ID check itself must stop the run
    with pytest.raises(A.AuditError, match="ROW_ID is blank"):
        A.load_inputs(S.make_config(proj, tmp_path / "o", bootstrap_reps=0, verbose=False))


def _independent_key_truth(project):
    si = pd.read_csv(project["paths"]["surgery_info"], dtype=str, keep_default_na=False)
    nid = lambda v: v.strip()[:-2] if v.strip().endswith(".0") else v.strip()  # noqa: E731
    k = pd.DataFrame({"s": si["SUBJECT_ID"].map(nid), "h": si["HADM_ID"].map(nid), "o": si["OPER_ID"].map(nid),
                      "b": si["SURGERY_BEGIN_TIME"].str.strip(),
                      "n": si["SURGERY_NAME"].map(lambda v: " ".join(v.lower().split()))})
    g = k.groupby(["s", "h", "o"]).size()
    d = k[k["b"] != ""].groupby(["s", "h", "b", "n"]).size()
    return {"shared_oper_id_groups": int((g >= 2).sum()), "rows_in_shared_oper_id_groups": int(g[g >= 2].sum()),
            "possible_true_duplicate_groups": int((d >= 2).sum()), "rows_in_possible_true_duplicate_groups": int(d[d >= 2].sum())}


def test_shared_oper_id_rows_all_kept(frame, project):
    f, meta, inp, _ = frame
    assert len(f) == project["n_surgeries"] == len(inp["surg"])
    assert f["ROW_ID"].is_unique and (f["ROW_ID"] != "").all()
    truth = _independent_key_truth(project)
    assert truth["shared_oper_id_groups"] >= 10 and truth["possible_true_duplicate_groups"] >= 1  # generator plants both
    for k, v in truth.items():
        assert meta["key_diag"][k] == v, k
    assert sum(meta["key_diag"]["shared_oper_id_group_sizes"].values()) == truth["shared_oper_id_groups"]
    # distinct surgeries that share an OPER_ID are all linked/analysed independently (nothing collapsed)
    shared = f[f.duplicated(["SUBJECT_ID", "HADM_ID", "OPER_ID"], keep=False)]
    assert (shared["link_class"] == "linked").sum() >= 10
    lk = shared[shared["link_class"] == "linked"]
    assert lk.duplicated(["SUBJECT_ID", "HADM_ID", "OPER_ID"], keep=False).any()


def test_key_diagnostics_published_suppressed(e2e, project):
    kd = e2e["doc"]["reconciliation"]["surgery_key_ROW_ID"]
    truth = _independent_key_truth(project)
    assert set(kd) == set(truth)
    for k, v in kd.items():
        assert v in ("<10", "S") or (isinstance(v, int) and v >= 10 and v == truth[k])
    local = json.loads(e2e["cfg"].local_json.read_text(encoding="utf-8"))
    assert local["extras"]["surgery_key_diagnostics"]["shared_oper_id_groups"] == truth["shared_oper_id_groups"]
    assert local["checks"]["surgery_key"] == "ROW_ID"


def test_no_shared_oper_ids_gives_zero_diagnostics(tmp_path):
    proj = S.make_project(tmp_path / "p", n_subjects=120, seed=11, shared_oper_every=0, true_dup_every=0)
    cfg = S.make_config(proj, tmp_path / "o", bootstrap_reps=0, verbose=False)
    f, meta = A.build_frame(A.load_inputs(cfg))
    assert meta["key_diag"]["shared_oper_id_groups"] == 0 and meta["key_diag"]["possible_true_duplicate_groups"] == 0


# ---------------------------------------------------------------- harmonization column defaults (v1.6)
def test_resolve_out_uses_outtime_fallback():
    hdr = ["encounter_id", "INTIME", "OUTTIME", "t0", "span_icu_hours"]
    assert A.resolve_column(hdr, None, A.HARM_T0_CANDIDATES, "T0", A.HARM_T0_FALLBACK) == "t0"
    assert A.resolve_column(hdr, None, A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK) == "OUTTIME"


def test_resolve_out_prefers_primary_over_fallback_and_explicit_wins():
    hdr = ["encounter_id", "OUT", "OUTTIME", "t0"]
    assert A.resolve_column(hdr, None, A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK) == "OUT"
    assert A.resolve_column(hdr, "OUTTIME", A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK) == "OUTTIME"
    with pytest.raises(A.AuditError):
        A.resolve_column(["encounter_id", "t0"], None, A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK)
    with pytest.raises(A.AuditError):
        A.resolve_column(hdr, "nope", A.HARM_OUT_CANDIDATES, "OUT", A.HARM_OUT_FALLBACK)


def test_real_harmonization_layout_runs_without_flags(tmp_path):
    proj = S.make_project(tmp_path / "p", n_subjects=150, seed=5, harm_layout="real")
    hdr = pd.read_csv(proj["paths"]["harmonization"], nrows=0).columns.tolist()
    assert "OUT" not in hdr and "OUTTIME" in hdr and "t0" in hdr
    cfg = S.make_config(proj, tmp_path / "o", bootstrap_reps=0, verbose=False)
    inp = A.load_inputs(cfg)
    assert inp["info"]["checks"]["harmonization_columns"] == {"t0": "t0", "OUT": "OUTTIME"}
    assert inp["info"]["checks"]["OUT_minus_t0_ne_span_icu_hours_gt_1min"] == 0


def test_out_vs_span_disagreement_is_counted(tmp_path):
    proj = S.make_project(tmp_path / "p", n_subjects=150, seed=5, harm_layout="real")
    h = pd.read_csv(proj["paths"]["harmonization"], dtype=str, keep_default_na=False)
    h.loc[:2, "span_icu_hours"] = "0.5"
    h.to_csv(proj["paths"]["harmonization"], index=False)
    proj["hashes"]["harmonization"] = __import__("hashlib").sha256(proj["paths"]["harmonization"].read_bytes()).hexdigest()
    inp = A.load_inputs(S.make_config(proj, tmp_path / "o", bootstrap_reps=0, verbose=False))
    assert inp["info"]["checks"]["OUT_minus_t0_ne_span_icu_hours_gt_1min"] == 3


def test_float_like_ids_link(frame):
    f, _, inp, _ = frame
    raw = inp["surg"]
    fl = raw["SUBJECT_ID"].str.endswith(".0")
    assert fl.any()
    assert (f.loc[fl.to_numpy(), "link_class"] == "linked").any()


# =============================================================================
# Durations, bands, groups
# =============================================================================
def test_durations_and_fallback():
    t = pd.Timestamp("2020-01-01 10:00")
    B = pd.Series([t, t, t, pd.NaT, pd.NaT, t])
    E = pd.Series([t + pd.Timedelta(hours=2), t, t - pd.Timedelta(hours=1), pd.NaT, pd.NaT, pd.NaT])
    AS = pd.Series([t - pd.Timedelta(hours=1)] * 4 + [pd.NaT, t])
    AE = pd.Series([t + pd.Timedelta(hours=9)] * 4 + [pd.NaT, t])
    d = A.compute_durations(B, E, AS, AE)
    assert d["D"].iloc[0] == 2.0 and not d["fallback"].iloc[0]          # fallback NOT used when B, E present
    assert d["invalid_reason"].tolist()[1:3] == ["zero_duration", "begin_after_end"]
    assert d["D"].iloc[3] == 10.0 and d["fallback"].iloc[3]
    assert d["invalid_reason"].iloc[4] == "missing_window"
    assert d["invalid_reason"].iloc[5] == "fallback_nonpositive"


@pytest.mark.parametrize("h,b", [(1.99, "lt2h"), (2.0, "2to4h"), (4.0, "2to4h"), (4.01, "gt4h")])
def test_bands(h, b):
    assert A.band_of(h) == b


def test_group_mapping():
    assert A.primary_group("Cardiac", "2", "major open", "G1") == "G1"
    assert A.primary_group("Cardiac", "NC", "major open", "NC") == "NC"
    assert A.primary_group("GI/abdominal", "", "catheter/interventional", "I1") == "I1"
    with pytest.raises(A.AuditError):
        A.primary_group("Cardiac", "4", "major open", "G2")  # map disagrees -> stop


def test_map_na_literal_kept(frame):
    _, _, inp, _ = frame
    assert (inp["map"]["intensity"] == "NA").any()  # 'NA' is a string, not a missing value


# =============================================================================
# Verdict truth table (signed D2)
# =============================================================================
@pytest.mark.parametrize("tau,d2", [(0.5, "not_met"), (0.2, "not_met"), (0.19, "met"), (0.0, "met"), (-0.5, "met"),
                                    (float("nan"), "not_assessable"), (None, "not_assessable")])
def test_d2_signed(tau, d2):
    v = A.evaluate_criteria([1.0, 2.0, 3.0], tau, 0.1, 0.1, 0.1)
    assert v["D2"] == d2
    assert v["big_discrepancy"] == ("Y" if d2 == "met" else "N")
    assert v["label"] == "approved_criteria_2026-11-10_D2signed"


def test_d2_negative_is_not_abs():
    assert A.evaluate_criteria([1.0, 2.0], -0.9, 0.0, 0.0, 0.0)["D2"] == "met"   # |tau| would say not_met


@pytest.mark.parametrize("meds,d1", [([1, 2, 3], "not_met"), ([1, 2, 2], "met"), ([3, 2, 1], "met"), ([1, 3, 2], "met"),
                                     ([1], "not_assessable"), (None, "not_assessable"), ([1, float("nan")], "not_assessable")])
def test_d1(meds, d1):
    assert A.evaluate_criteria(meds, 0.5, 0.1, 0.1, 0.1)["D1"] == d1


@pytest.mark.parametrize("share,d3", [(0.25, "not_met"), (0.2501, "met"), (0.0, "not_met"), (None, "not_assessable")])
def test_d3(share, d3):
    assert A.evaluate_criteria([1, 2], 0.5, share, 0.0, 0.0)["D3"] == d3


def test_d4_d5_never_change_big():
    v = A.evaluate_criteria([1, 2, 3], 0.5, 0.1, 0.99, 0.99)
    assert v["D4"] == "met" and v["D5"] == "met" and v["big_discrepancy"] == "N"


def test_all_not_assessable():
    v = A.evaluate_criteria(None, float("nan"), None, None, None)
    assert v["big_discrepancy"] == "not_assessable"


def test_mixed_not_assessable_and_not_met_is_N():
    v = A.evaluate_criteria(None, 0.6, None, None, None)
    assert v["D1"] == "not_assessable" and v["big_discrepancy"] == "N"


# =============================================================================
# Tau at surgery level
# =============================================================================
def test_tau_is_surgery_level_not_group_medians():
    rng = np.random.default_rng(0)
    base = rng.normal(3, 1.5, 40)
    x = np.repeat([0, 1, 2], 40).astype(float)
    y = np.concatenate([base, base + 0.01, base + 0.02])  # medians strictly increasing, heavy overlap
    t = A.kendall_tau_b_surgery_level(x, y)
    assert math.isclose(t, sps.kendalltau(x, y, variant="b")[0])
    med_tau = sps.kendalltau([0, 1, 2], [np.median(y[x == k]) for k in range(3)])[0]
    assert med_tau == 1.0 and t < 0.5


def test_tau_nan_cases():
    assert math.isnan(A.kendall_tau_b_surgery_level(np.zeros(20), np.arange(20.0)))
    assert math.isnan(A.kendall_tau_b_surgery_level(np.array([0.0, 1.0]), np.array([1.0, 2.0])))


def test_published_tau_matches_surgery_level_train(e2e, frame):
    f, meta, _, _ = frame
    d = f[f["analysed"] & f["valid"] & (f["body_system"] == "Cardiac") & (f["split"] == "train")]
    d = d[d["cgroup"].isin(["G1", "G2", "G3"])]
    x = d["cgroup"].map({"G1": 0, "G2": 1, "G3": 2}).to_numpy(float)
    exp = sps.kendalltau(x, d["D"].to_numpy(float), variant="b")[0]
    v = e2e["doc"]["systems"]["Cardiac"]["verdict"]
    assert v["tau_b"] == round(exp, 3)
    assert isinstance(v["tau_b_ci"], list) and v["tau_b_ci"][0] <= v["tau_b"] <= v["tau_b_ci"][1]


def test_inverse_system_flagged(e2e):
    v = e2e["doc"]["systems"]["Hepatobiliary"]["verdict"]
    assert v["tau_b"] < 0 and v["D2"] == "met" and v["D1"] == "met" and v["big_discrepancy"] == "Y"
    assert e2e["doc"]["systems"]["Cardiac"]["verdict"]["big_discrepancy"] == "N"


def test_cluster_bootstrap_resamples_subjects():
    rng = np.random.default_rng(1)
    inv = np.array([0, 0, 0, 1, 2, 2])
    w = A.Bootstrapper.cluster_weights(inv, 3, rng)
    assert w[0] == w[1] == w[2] and w[4] == w[5]  # all surgeries of a subject move together


def test_bootstrap_deterministic():
    subj = np.repeat(np.arange(30), 2).astype(str)
    x = np.tile([0.0, 1.0], 30)
    y = x + np.random.default_rng(2).normal(0, 1, 60)
    b1, b2 = A.Bootstrapper(200, 5), A.Bootstrapper(200, 5)
    f = lambda idx: A.kendall_tau_b_surgery_level(x[idx], y[idx])  # noqa: E731
    assert b1.ci(subj, f, "k")["ci"] == b2.ci(subj, f, "k")["ci"]


def test_references_are_train_only(project, tmp_path):
    """Shifting validation durations must not change any train threshold/verdict input."""
    shifted = S.make_project(tmp_path / "shift", n_subjects=600, seed=7, val_shift_hours=5.0)
    a1, _ = _analysis(project, tmp_path / "o1")
    a2, _ = _analysis(shifted, tmp_path / "o2")
    for bs in ("Cardiac", "GI/abdominal", "Hepatobiliary"):
        v1, v2 = a1.stats[A.P("systems", bs, "verdict")], a2.stats[A.P("systems", bs, "verdict")]
        for k in ("d1_medians_train", "tau_b", "D1", "D2"):
            assert v1[k] == v2[k]
    assert a1.stats[A.P("systems", "Cardiac", "groups", "G1", "validation", "median")] != \
        a2.stats[A.P("systems", "Cardiac", "groups", "G1", "validation", "median")]


# =============================================================================
# Groups < 10 and NC are reported, not dropped
# =============================================================================
def test_small_groups_and_nc_reported(e2e):
    doc = e2e["doc"]
    g = doc["systems"]["Cardiac"]["groups"]
    assert set(g) == {"G1", "G2", "G3", "NC"} and set(g["NC"]) == set(A.ALL_SPLITS)
    neuro = doc["systems"]["Neurosurgical"]
    assert "I2" in neuro["groups"] and neuro["groups"]["I2"]["train"]["n"] == "<10"
    assert neuro["groups"]["I2"]["train"]["median"] == "<10"
    assert neuro["verdict"]["big_discrepancy"] == "not_assessable" and neuro["not_assessable_reason"] == "groups_below_10"
    assert doc["systems"]["Vascular access~1line placement".replace("~1", "/")]["not_assessable_reason"] == "single_group"


def test_sensitivities_present(e2e):
    c = e2e["doc"]["systems"]["Cardiac"]
    for k in ("e_four_level_or_stat5", "a_single_component", "b_no_cross_boundary", "c_no_gt24_no_fallback", "d_during_icu"):
        assert k in c["sensitivity"]
    for k in ("review_flag_Y", "review_flag_N", "STAT_group_crosses_Y", "STAT_group_crosses_N"):
        assert k in c["flagged_vs_unflagged"]
    assert set(c["sensitivity"]["e_four_level_or_stat5"]["groups"]) == {"1", "2", "3", "4", "5", "NC"}
    assert "stratified_tau_b" in e2e["doc"]["pooled"] and "across_system" in e2e["doc"]
    assert "major open" in e2e["doc"]["across_system"]


def test_discordance_only_names_ge_10(e2e, frame):
    f, _, _, _ = frame
    d = f[f["analysed"] & f["valid"] & (f["body_system"] == "Neurosurgical")]
    listed = set(e2e["doc"]["systems"]["Neurosurgical"]["discordant_names"])
    assert listed == set(d.groupby("name").size().loc[lambda s: s >= 10].index)
    assert "endoscopic third ventriculostomy" not in listed


# =============================================================================
# Disclosure: suppression only on counts, complementary non-recoverable, gate
# =============================================================================
def _walk(o, parts=()):
    """Yields (JSON pointer, value) for every leaf."""
    if isinstance(o, dict):
        for k, v in o.items():
            yield from _walk(v, parts + (k,))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from _walk(v, parts + (i,))
    else:
        yield A.P(*parts), o


def test_suppression_only_on_count_fields(e2e):
    doc = e2e["doc"]
    counts = set(doc["disclosure"]["count_fields"])
    body = {k: v for k, v in doc.items() if k != "disclosure"}
    for path, v in _walk(body):
        if v == "S":
            assert path in counts, path
        if isinstance(v, int) and not isinstance(v, bool) and path in counts:
            assert v >= 10
    # medians of groups whose n is complementary-suppressed stay published
    g2 = doc["systems"]["Cardiac"]["groups"]["G2"]["train"]
    assert isinstance(g2["median"], float)
    # criteria, thresholds and taus never carry a suppression marker
    assert doc["criteria"]["tau_b_signed_max"] == 0.2
    for bs in ("Cardiac", "GI/abdominal", "Hepatobiliary"):
        assert isinstance(doc["systems"][bs]["verdict"]["tau_b"], float)
        assert all(isinstance(m, float) for m in doc["systems"][bs]["verdict"]["d1_medians_train"])


def test_no_small_counts_and_no_min_max(e2e):
    doc = e2e["doc"]
    for ptr in doc["disclosure"]["count_fields"]:
        v = A.get_ptr(doc, ptr)
        assert v in ("<10", "S") or (isinstance(v, int) and v >= 10)
    assert not [p for p, _ in _walk(doc) if p.rsplit("/", 1)[-1] in ("min", "max")]


def test_zero_policy_zero_is_marked(e2e):
    raw = json.loads(e2e["cfg"].local_json.read_text(encoding="utf-8"))["raw_counts"]
    zeros = [c for c, v in raw.items() if v == 0 and not c.startswith("hidden:")]
    assert zeros
    assert all(A.get_ptr(e2e["doc"], c) in ("<10", "S") for c in zeros)


def test_counters_are_real(e2e):
    doc = e2e["doc"]
    vals = [A.get_ptr(doc, p) for p in doc["disclosure"]["count_fields"]]
    s = doc["suppression"]
    assert s["primary_cells"] == vals.count("<10") and s["complementary_cells"] == vals.count("S")
    assert s["complementary_cells"] > 0 and s["recoverable"] == 0


def test_complementary_suppression_simple_table():
    reg = A.Registry()
    t = reg.count("/t", 25)
    reg.identity(t, [reg.count("/a", 20), reg.count("/b", 5)])
    known = {"/t": 25, "/a": 20}
    assert A.recoverability_check(reg.identities, known, {"/b": "P"})[1] == ["/b"]  # recoverable by subtraction
    status, supp = A.apply_suppression(reg)
    assert status["/b"] == "P" and "S" in status.values() and supp["recoverable"] == 0
    known = {c: reg.raw[c] for c, s in status.items() if s == "pub"}
    unknown = {c: s for c, s in status.items() if s != "pub"}
    assert A.recoverability_check(reg.identities, known, unknown)[1] == []


def test_complementary_suppression_two_way_table():
    """2x3 table with margins; a single small cell must not be recoverable from rows, columns or the grand total."""
    reg = A.Registry()
    vals = [[30, 4, 22], [15, 40, 18]]
    cells = [[reg.count(f"/c{i}{j}", vals[i][j]) for j in range(3)] for i in range(2)]
    rows = [reg.count(f"/r{i}", sum(vals[i])) for i in range(2)]
    cols = [reg.count(f"/k{j}", vals[0][j] + vals[1][j]) for j in range(3)]
    g = reg.count("/g", sum(map(sum, vals)))
    for i in range(2):
        reg.identity(rows[i], cells[i])
    for j in range(3):
        reg.identity(cols[j], [cells[0][j], cells[1][j]])
    reg.identity(g, rows)
    reg.identity(g, cols)
    status, _ = A.apply_suppression(reg)
    assert status["/c01"] == "P"
    known = {c: reg.raw[c] for c, s in status.items() if s == "pub"}
    unknown = {c: s for c, s in status.items() if s != "pub"}
    _, unsafe, bounds = A.recoverability_check(reg.identities, known, unknown)
    assert unsafe == []
    lo, hi = bounds["/c01"]
    assert hi - lo >= 1  # interval, not a point


def test_hidden_difference_protected():
    """Primary n=50 published, sensitivity n=45 published -> hidden difference 5 would be disclosed."""
    reg = A.Registry()
    p = reg.count("/p", 50)
    reg.subset(reg.count("/s", 45), p, "x")
    _, unsafe, _ = A.recoverability_check(reg.identities, {"/p": 50, "/s": 45}, {h: "H" for h in reg.hidden})
    assert unsafe
    status, _ = A.apply_suppression(reg)
    assert "S" in status.values()


def test_shared_output_is_not_recoverable_from_disk(e2e):
    assert A.release_gate(e2e["cfg"].shared_json, e2e["cfg"].shared_md, "ok", [S.SENTINEL_SUBJECT]) == []


def _copy_shared(e2e, tmp_path):
    j = tmp_path / e2e["cfg"].shared_json.name
    m = tmp_path / e2e["cfg"].shared_md.name
    shutil.copy(e2e["cfg"].shared_json, j)
    shutil.copy(e2e["cfg"].shared_md, m)
    return j, m


def test_gate_fails_on_injected_small_cell(e2e, tmp_path):
    j, m = _copy_shared(e2e, tmp_path)
    doc = json.loads(j.read_text(encoding="utf-8"))
    target = next(p for p in doc["disclosure"]["count_fields"] if A.get_ptr(doc, p) == "<10")
    A.set_ptr(doc, target, 5)
    doc["suppression"]["primary_cells"] -= 1
    j.write_text(json.dumps(doc), encoding="utf-8")
    m.write_text(A.render_markdown(doc), encoding="utf-8")
    fails = A.release_gate(j, m, "ok")
    assert any("below 10" in f for f in fails)


def test_gate_fails_on_undeclared_count_field(e2e, tmp_path):
    j, m = _copy_shared(e2e, tmp_path)
    doc = json.loads(j.read_text(encoding="utf-8"))
    doc["population"]["pooled"]["linked_extra"] = {"n": 7}
    j.write_text(json.dumps(doc), encoding="utf-8")
    m.write_text(A.render_markdown(doc), encoding="utf-8")
    assert any("schema" in f for f in A.release_gate(j, m, "ok"))


def test_gate_fails_on_sentinel_and_edited_md(e2e, tmp_path):
    j, m = _copy_shared(e2e, tmp_path)
    assert A.release_gate(j, m, "ok", [S.SENTINEL_SUBJECT]) == []
    m.write_text(m.read_text(encoding="utf-8") + "\nsubject 987654321\n", encoding="utf-8")
    fails = A.release_gate(j, m, "ok", [S.SENTINEL_SUBJECT])
    assert any("sentinel" in f for f in fails) and any("markdown" in f for f in fails)


def test_gate_fails_on_recoverable_cell(e2e, tmp_path):
    """Un-suppress one complementary cell: the solver must detect the now-recoverable cells."""
    j, m = _copy_shared(e2e, tmp_path)
    doc = json.loads(j.read_text(encoding="utf-8"))
    raw = json.loads(e2e["cfg"].local_json.read_text(encoding="utf-8"))["raw_counts"]
    s_cells = [p for p in doc["disclosure"]["count_fields"] if A.get_ptr(doc, p) == "S"]
    for target in s_cells:
        trial = copy.deepcopy(doc)
        A.set_ptr(trial, target, raw[target])
        trial["suppression"]["complementary_cells"] -= 1
        j.write_text(json.dumps(trial), encoding="utf-8")
        m.write_text(A.render_markdown(trial), encoding="utf-8")
        fails = A.release_gate(j, m, "ok")
        if fails:
            assert any("recoverab" in f for f in fails)
            return
    pytest.fail("no complementary cell was necessary")


def test_publish_failure_leaves_no_shared_file(e2e, tmp_path):
    cfg = copy.deepcopy(e2e["cfg"])
    cfg.out_dir = tmp_path
    doc = json.loads(e2e["cfg"].shared_json.read_text(encoding="utf-8"))
    doc["release_gate"] = "pending"
    target = next(p for p in doc["disclosure"]["count_fields"] if A.get_ptr(doc, p) == "<10")
    A.set_ptr(doc, target, 3)
    cfg.shared_json.write_text("stale file from an earlier run", encoding="utf-8")
    with pytest.raises(A.ReleaseGateError):
        A.publish(cfg, doc)
    assert sorted(p.name for p in tmp_path.iterdir()) == []


def test_main_exit_codes(monkeypatch):
    def boom_gate(cfg, sentinels=()):
        raise A.ReleaseGateError("x")

    def boom_hash(cfg, sentinels=()):
        raise A.HashMismatchError("x")
    argv = [
        "--project-root", "/nonexistent",
        "--raw-dir", "/nonexistent/raw",
        "--split-file", "/nonexistent/splits.csv",
    ]
    monkeypatch.setattr(A, "run_audit", boom_gate)
    assert A.main(argv) == 3
    monkeypatch.setattr(A, "run_audit", boom_hash)
    assert A.main(argv) == 2


def test_outputs_and_local_separation(e2e):
    out = e2e["out"]
    names = sorted(p.name for p in out.iterdir())
    assert names == ["surgery_complexity_duration_audit_MAP_v1.json", "surgery_complexity_duration_audit_MAP_v1.md",
                     "surgery_complexity_duration_audit_MAP_v1_raw_LOCAL.json"]
    doc = e2e["doc"]
    assert doc["release_gate"] == "ok" and doc["spec"] == "surgery_complexity_duration_audit v1.7"
    assert doc["criteria"]["d2_rule"].startswith("met if tau_b < ") and doc["criteria"]["signed"] is True
    assert doc["bootstrap"]["replicates"] == 30 and doc["bootstrap"]["cluster"] == "SUBJECT_ID"
    assert set(doc["input_hashes"]) == {"canonical", "harmonization", "complexity_map", "split"}
    text = e2e["cfg"].shared_json.read_text(encoding="utf-8") + e2e["cfg"].shared_md.read_text(encoding="utf-8")
    assert S.SENTINEL_SUBJECT not in text and "LOCAL ONLY" not in text
    assert "LOCAL ONLY" in e2e["cfg"].local_json.read_text(encoding="utf-8")


def test_markdown_is_rendering_of_json(e2e):
    assert e2e["cfg"].shared_md.read_text(encoding="utf-8") == A.render_markdown(e2e["doc"])


def test_canonical_split_column_must_agree(project, tmp_path):
    proj = copy.deepcopy(project)
    c = pd.read_csv(proj["paths"]["canonical"], dtype=str, keep_default_na=False)
    sp = pd.read_csv(proj["paths"]["split"], dtype=str).set_index("SUBJECT_ID")["split"]
    c["split"] = c["SUBJECT_ID"].map(sp)
    p = tmp_path / "canon.csv"
    c.to_csv(p, index=False)
    proj["paths"]["canonical"] = p
    proj["hashes"]["canonical"] = A.sha256_file(p)
    cfg = S.make_config(proj, tmp_path, bootstrap_reps=0, verbose=False)
    assert A.load_inputs(cfg)["info"]["checks"]["canonical_split_disagreements"] == 0
    c.loc[0, "split"] = "test" if c.loc[0, "split"] != "test" else "train"
    c.to_csv(p, index=False)
    proj["hashes"]["canonical"] = A.sha256_file(p)
    cfg = S.make_config(proj, tmp_path, bootstrap_reps=0, verbose=False)
    with pytest.raises(A.AuditError, match="disagrees"):
        A.load_inputs(cfg)
