"""pytest suite for phi1_band_support_audit.py. SYNTHETIC DATA ONLY.

Synthetic inputs come from ../duration_audit_v2/synthetic_data.py (no real PIC file is read). They get their own
SHA-256 values through BandConfig.pinned_sha256 (test override), so the hash gates still run on every test.
The duration audit v1.7 module is loaded from ../duration_audit_v2 after its pinned SHA-256 is verified.

Run:  python -m pytest -q test_phi1_band_support_audit.py
"""
from __future__ import annotations

import contextlib
import io
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
V17_DIR = HERE.parent / "duration_audit_v2"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(V17_DIR))
import phi1_band_support_audit as B  # noqa: E402
import synthetic_data as S  # noqa: E402

V17 = V17_DIR / B.V17_MODULE_FILE


def make_cfg(proj: dict, out: Path, **kw) -> B.BandConfig:
    p = proj["paths"]
    base = dict(canonical_path=p["canonical"], harmonization_path=p["harmonization"], map_path=p["complexity_map"],
                split_path=p["split"], surgery_info_path=p["surgery_info"], admissions_path=p["admissions"],
                out_dir=Path(out), v17_module_path=V17, pinned_sha256=dict(proj["hashes"]),
                expected_surgery_rows=proj["n_surgeries"], expected_encounters=proj["n_encounters"],
                expected_subjects=proj["n_subjects"], expected_split_encounters=proj["split_encounters"],
                expected_map_rows=proj["map_rows"], expected_missing=proj["missing"], verbose=False)
    base.update(kw)
    return B.BandConfig(**base)


@pytest.fixture(scope="session")
def A():
    return B.load_v17(V17)


@pytest.fixture(scope="session")
def project(tmp_path_factory):
    return S.make_project(tmp_path_factory.mktemp("syn"), n_subjects=900, seed=13, harm_layout="real")


@pytest.fixture(scope="session")
def e2e(project, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    cfg = make_cfg(project, out, verbose=True)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        doc = B.run_audit(cfg, sentinels=[S.SENTINEL_SUBJECT])
    return {"cfg": cfg, "doc": doc, "out": out, "stdout": buf.getvalue()}


@pytest.fixture(scope="session")
def prepared(project, A, tmp_path_factory):
    cfg = make_cfg(project, tmp_path_factory.mktemp("p"))
    inp = A.load_inputs(B.v17_config(cfg, A))
    frame, _ = A.build_frame(inp)
    return B.prepare(frame, inp["map"], A.UNCLASSIFIED)


def _train(rows):
    """rows: (type, system, group, D) -> train frame."""
    return pd.DataFrame(rows, columns=["type", "system", "group", "D"]).assign(split="train")


def _toy():
    rows = [("A", "S1", "G1", 1.0 + 0.1 * i) for i in range(12)]       # own type at 10, not at 20
    rows += [("B", "S1", "G1", 2.0), ("B", "S1", "G1", 2.2), ("B", "S1", "G1", 2.4)]  # pool S1xG1 = 15
    rows += [("C", "S1", "G2", 5.0), ("C", "S1", "G2", 5.5)]          # pool S1xG2 = 2; system S1 = 17
    rows += [("D", "S2", "I1", 9.0), ("D", "S2", "I1", 9.5), ("D", "S2", "I1", 10.0)]  # system S2 = 3 -> all
    return _train(rows)


# =============================================================================
# Bands and cut points
# =============================================================================
@pytest.mark.parametrize("r,band", [(0.0001, "short"), (0.7499, "short"), (0.749999999, "short"), (0.75, "expected"),
                                    (1.0, "expected"), (1.33, "expected"), (1.3300001, "long"), (1.3301, "long"), (7.0, "long")])
def test_band_boundaries(r, band):
    assert B.band_of_ratio(r) == band


@pytest.mark.parametrize("m", [2.21, 1.37, 0.51, 3.3333333, 0.1, 6.84])
def test_cut_points_inclusive_after_ratio(m):
    assert B.band_of_ratio(B.ratio(0.75 * m, m)) == "expected"
    assert B.band_of_ratio(B.ratio(1.33 * m, m)) == "expected"
    assert B.band_of_ratio(B.ratio(0.7499 * m, m)) == "short"
    assert B.band_of_ratio(B.ratio(1.3301 * m, m)) == "long"


def test_config_values_are_the_approved_ones():
    assert (B.SHORT_BELOW, B.LONG_ABOVE, B.PRIMARY_MIN_N, B.SENSITIVITY_MIN_N) == (0.75, 1.33, 10, 20)
    assert B.LEVELS == ("own_type", "system_x_group", "system", "all")


# =============================================================================
# Fallback order and thresholds
# =============================================================================
def test_fallback_order_threshold_10():
    tr = _toy()
    fr = B.freeze_medians(tr, 10)
    k, m = fr.reference("A", "S1", "G1")
    assert k == ("own_type", "A") and m == pytest.approx(np.median(tr[tr.type == "A"].D))
    k, m = fr.reference("B", "S1", "G1")
    assert k == ("system_x_group", "S1", "G1") and m == pytest.approx(np.median(tr[(tr.system == "S1") & (tr.group == "G1")].D))
    k, m = fr.reference("C", "S1", "G2")
    assert k == ("system", "S1") and m == pytest.approx(np.median(tr[tr.system == "S1"].D))
    k, m = fr.reference("D", "S2", "I1")
    assert k == ("all",) and m == pytest.approx(np.median(tr.D))


def test_fallback_order_threshold_20():
    tr = _toy()
    fr = B.freeze_medians(tr, 20)
    assert fr.own_types == set()
    for t, bs, g in (("A", "S1", "G1"), ("B", "S1", "G1"), ("C", "S1", "G2"), ("D", "S2", "I1")):
        assert fr.reference(t, bs, g)[0] == ("all",)  # S1xG1 = 15, S1 = 17, S2 = 3, all = 20


def test_threshold_boundary_exactly_min_n():
    rows = [("A", "S1", "G1", 1.0)] * 10 + [("B", "S1", "G1", 3.0)] * 9 + [("C", "S2", "I1", 2.0)] * 20
    fr = B.freeze_medians(_train(rows), 10)
    assert fr.reference("A", "S1", "G1")[0] == ("own_type", "A")          # n = 10 -> own
    assert fr.reference("B", "S1", "G1")[0] == ("system_x_group", "S1", "G1")  # n = 9 -> pool of 19
    fr20 = B.freeze_medians(_train(rows), 20)
    assert fr20.reference("C", "S2", "I1")[0] == ("own_type", "C")        # n = 20 -> own at 20
    assert fr20.reference("A", "S1", "G1")[0] == ("all",)                 # S1 pool = 19 < 20


def test_unseen_type_goes_to_fallback():
    fr = B.freeze_medians(_toy(), 10)
    assert fr.reference("NEW", "S1", "G1")[0] == ("system_x_group", "S1", "G1")
    assert fr.reference("NEW", "S9", "X")[0] == ("all",)


def test_all_level_needs_min_n_and_train_only():
    with pytest.raises(B.AuditError):
        B.freeze_medians(_train([("A", "S1", "G1", 1.0)] * 9), 10)
    bad = _toy()
    bad.loc[0, "split"] = "test"
    with pytest.raises(B.AuditError):
        B.freeze_medians(bad, 10)


def test_assign_bands_uses_reference_median():
    tr = _toy()
    fr = B.freeze_medians(tr, 10)
    m = fr.type_stats["A"][1]
    df = pd.DataFrame({"type": ["A", "A", "A"], "system": ["S1"] * 3, "group": ["G1"] * 3,
                       "D": [0.75 * m, 1.33 * m, 1.34 * m]})
    out = B.assign_bands(df, fr)
    assert list(out.band) == ["expected", "expected", "long"] and set(out.level) == {"own_type"}


# =============================================================================
# Train-only medians (no validation/test leakage)
# =============================================================================
def test_medians_ignore_validation_and_test(prepared):
    bd = prepared[prepared.banded]
    tr = bd[bd.split == "train"][["type", "system", "group", "D", "split"]]
    f1 = B.freeze_medians(tr, 10)
    bd2 = bd.copy()
    bd2.loc[bd2.split != "train", "D"] = bd2.loc[bd2.split != "train", "D"] * 7 + 3
    tr2 = bd2[bd2.split == "train"][["type", "system", "group", "D", "split"]]
    f2 = B.freeze_medians(tr2, 10)
    assert f1.type_stats == f2.type_stats and f1.sg_stats == f2.sg_stats and f1.all_stats == f2.all_stats


def test_e2e_frozen_medians_equal_train_medians(e2e, prepared):
    local = json.loads(e2e["cfg"].local_json.read_text(encoding="utf-8"))
    tr = prepared[prepared.banded & (prepared.split == "train")]
    for row in local["extras"]["frozen_medians_min_n_10"]:
        if row["level"] == "own_type":
            t = row["key"].removeprefix("type:")
            assert row["median_hours"] == pytest.approx(float(np.median(tr[tr.type == t].D)))
            assert row["train_n"] == int((tr.type == t).sum())
    assert local["extras"]["checks_train_rows_used_for_medians"] == len(tr)


def test_tertiles_train_only_and_own_types_only():
    tr = _toy()
    cuts = B.freeze_tertiles(tr, {"A"})
    assert set(cuts) == {"A"}
    q1, q2 = np.quantile(tr[tr.type == "A"].D, [1 / 3, 2 / 3])
    assert cuts["A"] == (pytest.approx(q1), pytest.approx(q2))
    assert B.tertile_band(q1 - 1e-9, *cuts["A"]) == "short"
    assert B.tertile_band(q1, *cuts["A"]) == "expected" and B.tertile_band(q2, *cuts["A"]) == "expected"
    assert B.tertile_band(q2 + 1e-9, *cuts["A"]) == "long"


# =============================================================================
# Durations: surgery window, anaesthesia fallback, invalid and > 24 h
# =============================================================================
def test_duration_surgery_window_then_anesthesia_fallback(A):
    t = lambda s: pd.Series(pd.to_datetime(s))  # noqa: E731
    Bt = t(["2020-01-01 10:00", None, "2020-01-01 10:00", "2020-01-01 10:00", None])
    Et = t(["2020-01-01 12:00", "2020-01-01 12:00", "2020-01-01 09:00", "2020-01-01 10:00", None])
    AS = t(["2020-01-01 09:00", "2020-01-01 09:30", "2020-01-01 08:00", "2020-01-01 08:00", None])
    AE = t(["2020-01-01 13:00", "2020-01-01 12:30", "2020-01-01 11:00", "2020-01-01 11:00", None])
    d = A.compute_durations(Bt, Et, AS, AE)
    assert d.D[0] == pytest.approx(2.0) and not d.fallback[0]          # surgery window wins
    assert d.D[1] == pytest.approx(3.0) and d.fallback[1]              # B missing -> AE - AS
    assert np.isnan(d.D[2]) and d.invalid_reason[2] == "begin_after_end"
    assert np.isnan(d.D[3]) and d.invalid_reason[3] == "zero_duration"
    assert np.isnan(d.D[4]) and d.invalid_reason[4] == "missing_window"


def _fake_frame():
    return pd.DataFrame({
        "ROW_ID": ["1", "2", "3", "4", "5"], "SUBJECT_ID": ["9"] * 5, "name": ["x", "x", "x", "x", "x"],
        "link_class": ["linked", "linked", "linked", "linked", "post_icu_same_hadm"],
        "split": ["train", "train", "train", "train", ""], "mapped": [True] * 5,
        "D": [2.0, 3.0, np.nan, 30.0, 2.0], "valid": [True, True, False, True, True],
        "fallback": [False, True, False, False, False], "gt24": [False, False, False, True, False],
        "invalid_reason": ["", "", "zero_duration", "", ""], "body_system": ["Cardiac"] * 5, "cgroup": ["G1"] * 5})


def test_prepare_population_flags():
    cmap = pd.DataFrame({"key": ["x"], "standardized_name": ["type x"]})
    df = B.prepare(_fake_frame(), cmap, "Unclassified")
    assert list(df.banded) == [True, True, False, True, False]   # fallback kept; invalid excluded; >24 h kept; unlinked out
    assert set(df.type) == {"type x"}


def test_prepare_unclassified_group_and_reserved_names():
    f = _fake_frame()
    f["body_system"] = "Unclassified"; f["cgroup"] = ""
    df = B.prepare(f, pd.DataFrame({"key": ["x"], "standardized_name": ["y"]}), "Unclassified")
    assert set(df.group) == {B.UNCLASSIFIED_GROUP}
    with pytest.raises(B.AuditError):
        B.prepare(_fake_frame(), pd.DataFrame({"key": ["x"], "standardized_name": ["long"]}), "Unclassified")


def test_e2e_invalid_excluded_and_gt24_counted(e2e, prepared):
    raw = json.loads(e2e["cfg"].local_json.read_text(encoding="utf-8"))["raw_counts"]
    for sp in B.SPLITS:
        d = prepared[(prepared.link_class == "linked") & (prepared.split == sp) & prepared.mapped]
        assert raw[f"/population/{sp}/valid"] == int(d.valid.sum())
        assert raw[f"/population/{sp}/invalid_duration"] == int((~d.valid).sum())
        assert raw[f"/population/{sp}/d_gt_24h"] == int(d.gt24.sum())
        assert raw[f"/population/{sp}/duration_fallback_anes"] == int((d.valid & d.fallback).sum())
        assert raw[f"/bands/min_n_10/{sp}/overall/n"] == int(d.valid.sum())  # > 24 h stays in the bands
    assert sum(raw[f"/population/{sp}/duration_fallback_anes"] for sp in B.SPLITS) > 0


# =============================================================================
# End to end: outputs, suppression, release gate
# =============================================================================
def test_e2e_outputs_and_release_gate(e2e, A):
    cfg, doc = e2e["cfg"], e2e["doc"]
    assert doc["release_gate"] == "ok" and doc["spec"] == B.SPEC_VERSION
    for p in (cfg.shared_json, cfg.shared_md, cfg.shared_csv, cfg.local_json, cfg.local_csv):
        assert p.exists()
    assert sorted(x.name for x in cfg.out_dir.iterdir()) == sorted(
        x.name for x in (cfg.shared_json, cfg.shared_md, cfg.shared_csv, cfg.local_json, cfg.local_csv))
    assert B.release_gate(A, cfg.shared_json, cfg.shared_md, cfg.shared_csv, "ok", [S.SENTINEL_SUBJECT]) == []
    assert doc["suppression"]["recoverable"] == 0 and doc["suppression"]["primary_cells"] > 0
    assert "raw_LOCAL" not in cfg.shared_md.read_text(encoding="utf-8") + cfg.shared_csv.read_text(encoding="utf-8")


def test_e2e_no_published_count_below_10(e2e):
    doc = e2e["doc"]
    for ptr in B.detect_count_fields(doc):
        cur = doc
        for part in ptr.split("/")[1:]:
            cur = cur[part.replace("~1", "/").replace("~0", "~")]
        assert cur in ("<10", "S") or (isinstance(cur, int) and cur >= 10)


def test_e2e_reconciliation_and_structure(e2e, project):
    doc = e2e["doc"]
    raw = json.loads(e2e["cfg"].local_json.read_text(encoding="utf-8"))["raw_counts"]
    assert raw["/reconciliation/surgery_info_rows"] == project["n_surgeries"]
    for k in (10, 20):
        for sp in B.SPLITS:
            assert sum(raw[f"/median_source/min_n_{k}/{sp}/{lv}/surgeries"] for lv in B.LEVELS) == raw[f"/population/{sp}/valid"]
            assert sum(raw[f"/bands/min_n_{k}/{sp}/overall/{b}"] for b in B.BANDS) == raw[f"/population/{sp}/valid"]
    assert set(doc["bands"]["min_n_10"]["train"]) == {"overall", "by_system", "by_group"}
    assert set(doc["unseen_in_train"]) == {"validation", "test"}
    assert set(doc["tertile_sensitivity"]) == set(B.SPLITS)
    for sp in B.SPLITS:
        for t, o in doc["own_types"][sp].items():
            assert o["any_band_below_10"] in ("Y", "N")


def test_e2e_frozen_csv(e2e):
    rows = list(pd.read_csv(e2e["cfg"].shared_csv, dtype=str, keep_default_na=False).itertuples(index=False))
    assert rows and all(r.train_n != "<10" for r in rows)
    assert {r.threshold for r in rows} == {"10", "20"}
    assert any(r.key == "all" for r in rows)


def test_console_summary_is_short(e2e):
    lines = [x for x in e2e["stdout"].splitlines() if x.strip()]
    assert len(lines) <= 14 and lines[-1] == "release gate ok"
    assert S.SENTINEL_SUBJECT not in e2e["stdout"]


def _copy_shared(e2e, tmp):
    cfg = e2e["cfg"]
    paths = []
    for p in (cfg.shared_json, cfg.shared_md, cfg.shared_csv):
        q = Path(tmp) / p.name
        shutil.copy(p, q)
        paths.append(q)
    return paths


def _first_published(doc):
    for ptr in sorted(B.detect_count_fields(doc)):
        parts = [p.replace("~1", "/").replace("~0", "~") for p in ptr.split("/")[1:]]
        cur = doc
        for p in parts[:-1]:
            cur = cur[p]
        if isinstance(cur[parts[-1]], int):
            return cur, parts[-1]
    raise AssertionError


@pytest.mark.parametrize("tamper", ["small_count", "pct_without_count", "min_key", "timestamp", "md_edit", "csv_edit", "csv_small_n", "sentinel"])
def test_release_gate_catches_tampering(e2e, A, tmp_path, tamper):
    j, m, c = _copy_shared(e2e, tmp_path)
    doc = json.loads(j.read_text(encoding="utf-8"))
    if tamper == "small_count":
        cur, k = _first_published(doc); cur[k] = 5
    elif tamper == "pct_without_count":
        o = doc["bands"]["min_n_10"]["train"]["overall"]
        o["short"] = "S"; o["pct_short"] = 12.3
    elif tamper == "min_key":
        doc["config"]["min"] = 0.1
    elif tamper == "timestamp":
        doc["config"]["note"] = "2020-01-01 10:00"
    elif tamper == "sentinel":
        doc["config"]["note"] = S.SENTINEL_SUBJECT
    if tamper in ("small_count", "pct_without_count", "min_key", "timestamp", "sentinel"):
        j.write_text(json.dumps(doc), encoding="utf-8")
        m.write_text(B.render_markdown(doc), encoding="utf-8")
        c.write_text(B.render_csv(doc), encoding="utf-8")
    elif tamper == "md_edit":
        m.write_text(m.read_text(encoding="utf-8") + "\nextra\n", encoding="utf-8")
    elif tamper == "csv_edit":
        c.write_text(c.read_text(encoding="utf-8").replace(",all,", ",ALL,"), encoding="utf-8")
    elif tamper == "csv_small_n":
        c.write_text(c.read_text(encoding="utf-8") + "10,type:x,own_type,<10,1.0\n", encoding="utf-8")
    assert B.release_gate(A, j, m, c, "ok", [S.SENTINEL_SUBJECT]) != []


def test_complementary_suppression_protects_differencing(A):
    reg = A.Registry()
    t = reg.count("/x/total", 40)
    a = reg.count("/x/a", 5)
    b = reg.count("/x/b", 35)
    reg.identity(t, [a, b])
    status, supp = A.apply_suppression(reg)
    assert status[a] == "P" and (status[t] == "S" or status[b] == "S") and supp["recoverable"] == 0


def test_superset_identity_hides_rest(A):
    an = B.BandAnalysis(A, pd.DataFrame(), {})
    tot = an.c(30, "t", "n")
    m1 = an.c(25, "t", "m1")
    an.superset(tot, [m1], "x")
    status, _ = A.apply_suppression(an.reg)
    assert status[tot] == "S" or status[m1] == "S"   # rest = 5 would otherwise be recoverable


# =============================================================================
# Hash gates stop the run
# =============================================================================
def test_input_hash_mismatch_stops_run(project, tmp_path):
    bad = dict(project["hashes"]); bad["complexity_map"] = "0" * 64
    cfg = make_cfg(project, tmp_path / "o", pinned_sha256=bad)
    with pytest.raises(B.HashMismatchError):
        B.run_audit(cfg)
    assert not (tmp_path / "o").exists() or not any((tmp_path / "o").iterdir())


def test_v17_module_hash_mismatch_stops_run(project, tmp_path):
    edited = tmp_path / B.V17_MODULE_FILE
    edited.write_bytes(V17.read_bytes() + b"\n# edited\n")
    cfg = make_cfg(project, tmp_path / "o", v17_module_path=edited)
    with pytest.raises(B.HashMismatchError):
        B.run_audit(cfg)
    assert not (tmp_path / "o").exists()


def test_v17_module_missing_stops_run(project, tmp_path):
    with pytest.raises(B.AuditError):
        B.run_audit(make_cfg(project, tmp_path / "o", v17_module_path=tmp_path / "nope.py"))


def test_cli_with_default_pins_on_synthetic_inputs_exits_2(project, tmp_path, capsys):
    p = project["paths"]
    rc = B.main(["--canonical", str(p["canonical"]), "--harmonization", str(p["harmonization"]), "--map", str(p["complexity_map"]),
                 "--split-file", str(p["split"]), "--raw-dir", str(p["surgery_info"].parent), "--out-dir", str(tmp_path / "o"),
                 "--v17-module", str(V17)])
    assert rc == 2 and "AUDIT STOPPED" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()
