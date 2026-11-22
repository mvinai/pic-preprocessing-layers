"""Synthetic-data tests for 06_08_build_levels_preprocess.py (stages 06-08, plan v13).

The Stage 3 v1.0.4 fixture (PIC-format, synthetic) is built once per session with the Stage 3 test module in
../stage3_builder/, then this script runs on its 46 outputs. No real patient data is read. The pins for the real
canonical file and split file are switched off for the fixture only (the real run keeps them)."""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
S3 = HERE.parent / "stage3_builder"
spec = importlib.util.spec_from_file_location("s0608", HERE / "06_08_build_levels_preprocess.py")
M = importlib.util.module_from_spec(spec)
sys.modules["s0608"] = M
spec.loader.exec_module(M)

REAL_PINS = dict(M.PINS)


def _build_stage3_fixture(work: Path) -> Path:
    """Stage 3 v1.0.4 outputs on the synthetic PIC-format fixture (same recipe as the Stage 3 tests)."""
    sys.path.insert(0, str(S3))
    import test_03_build_real_layers as T  # noqa: E402
    B, FX = T.B, T.FX
    fx = FX.make_project(work / "proj")
    ph = work / "ph.json"
    fx = dict(fx, s2=ph, s2_sha=FX.make_stage2_release(ph))
    cap = {}
    real = B.release_gate

    class Captured(Exception):
        pass

    def grab(SM, g, *a, **k):
        cap.update(g.values())
        raise Captured()

    B.release_gate = grab
    try:
        try:
            T.run_build(T.make_cfg(fx, work / "cap"))
        except Captured:
            pass
    finally:
        B.release_gate = real
    s2 = work / "feature_sparsity_frequency_audit_v1.0.3.json"
    sha = FX.make_stage2_release(s2, cap, fx["expect"]["splits"], protect=True)
    cp = work / "consts.csv"
    csha = FX.make_stage2_constants(cp, stop_event=False)
    fx = dict(fx, s2=s2, s2_sha=sha)
    with contextlib.redirect_stdout(io.StringIO()):
        T.run_build(T.make_cfg(fx, work / "out", stage2_constants_path=cp, stage2_constants_sha256=csha))
    return work


@pytest.fixture(scope="session")
def fx(tmp_path_factory):
    work = _build_stage3_fixture(tmp_path_factory.mktemp("s3"))
    out = work / "out"
    A = pd.read_csv(out / "stage3_layer_a_base.csv", dtype=str, keep_default_na=False)
    M.PINS.update(canonical=None, split=None)          # fixture only: the real files are pinned in the script
    M.EXPECT.encounters = len(A)
    M.EXPECT.splits = A["split"].value_counts().to_dict()
    M.EXPECT.eligible = {L: len(np.load(out / f"stage3_layer_c_filled_L{L}.npz")["encounter_id"])
                         for L in M.LANDMARKS}
    man = M.sha256_file(out / "stage3_manifest.json")
    return dict(work=work, s3=out, man=man,
                base=["--project-root", str(work / "proj"), "--stage3-dir", str(out),
                      "--split-file", str(work / "proj" / "subject_splits.csv")])


def run(fx, out, *extra):
    buf, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
        rc = M.main(fx["base"] + ["--out-dir", str(out)] + list(extra))
    return rc, buf.getvalue() + err.getvalue()


@pytest.fixture(scope="session")
def built(fx, tmp_path_factory):
    out = tmp_path_factory.mktemp("run") / "stage06_08"
    rc, txt = run(fx, out, "--expect-stage3-manifest-sha256", fx["man"])
    assert rc == 0, txt
    return dict(out=out, txt=txt)


# --------------------------------------------------------------------------------------------------- pins
def test_real_pins_full_sha256():
    for k in ("stage3_script_sha256", "channel_map", "phi1_medians", "canonical", "split"):
        v = REAL_PINS[k]
        assert isinstance(v, str) and len(v) == 64 and all(c in "0123456789abcdef" for c in v), k
    assert REAL_PINS["stage3_script_version"] == "03_build_real_layers v1.0.4"
    assert REAL_PINS["stage3_files"] == 46


def test_real_expectations():
    e = M.Expect()
    assert e.encounters is None and e.splits == {} and e.eligible == {}


# --------------------------------------------------------------------------------------------------- gates
def test_precheck_writes_nothing(fx, tmp_path):
    out = tmp_path / "o"
    rc, txt = run(fx, out, "--precheck-only")
    assert rc == 0, txt
    assert "PRECHECK OK" in txt and f"stage3_manifest.json sha256 {fx['man']}" in txt
    assert not out.exists() and not list(tmp_path.glob("*.pending-*"))


def test_run_needs_manifest_hash_and_x_option_is_gone(fx, tmp_path):
    rc, txt = run(fx, tmp_path / "o")
    assert rc == 2 and "expect-stage3-manifest-sha256" in txt
    with pytest.raises(SystemExit):     # the single-X rule was replaced by decision D1 (fixed, recorded settings)
        M.build_config(["--adm-src-x-hours", "48", "--expect-stage3-manifest-sha256", fx["man"]])


def test_wrong_manifest_hash_stops(fx, tmp_path):
    rc, txt = run(fx, tmp_path / "o", "--expect-stage3-manifest-sha256", "0" * 64)
    assert rc == 2 and "differs" in txt and not (tmp_path / "o").exists()


def test_tampered_stage3_file_stops(fx, tmp_path):
    s3 = tmp_path / "s3"
    shutil.copytree(fx["s3"], s3)
    p = s3 / "stage3_layer_b_base.csv"
    p.write_bytes(p.read_bytes() + b"\n")
    base = [a if a != str(fx["s3"]) else str(s3) for a in fx["base"]]
    with contextlib.redirect_stdout(io.StringIO()) as buf, contextlib.redirect_stderr(io.StringIO()) as err:
        rc = M.main(base + ["--out-dir", str(tmp_path / "o"),
                            "--expect-stage3-manifest-sha256", fx["man"]])
    assert rc == 2, buf.getvalue() + err.getvalue()
    assert not (tmp_path / "o").exists()


def test_extra_file_in_stage3_folder_stops(fx, tmp_path):
    s3 = tmp_path / "s3"
    shutil.copytree(fx["s3"], s3)
    (s3 / "stray.csv").write_text("x\n")
    base = [a if a != str(fx["s3"]) else str(s3) for a in fx["base"]]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        rc = M.main(base + ["--out-dir", str(tmp_path / "o"), "--precheck-only"])
    assert rc == 2


def test_nonempty_out_dir_refused(fx, tmp_path):
    out = tmp_path / "o"
    out.mkdir()
    (out / "keep.txt").write_text("x")
    rc, txt = run(fx, out, "--expect-stage3-manifest-sha256", fx["man"])
    assert rc == 2 and (out / "keep.txt").read_text() == "x"


# --------------------------------------------------------------------------------------------------- full run
def test_full_run_files_and_checks(built):
    out = built["out"]
    names = {p.name for p in out.iterdir()}
    for lv in ("a0", "a1", "a2"):
        assert f"layer_a_{lv}.csv" in names and f"stage08_{lv}.npz" in names
    for lv in ("b0", "b1", "b2s", "b2d", "b2p"):
        assert f"layer_b_{lv}.csv" in names and f"stage08_{lv}.npz" in names
    for L in M.LANDMARKS:
        assert f"stage08_layer_c_L{L}.npz" in names
    assert "result: STAGES 06-08 OK; RELEASE GATE OK" in built["txt"]
    a = json.loads((out / "stage06_08_audit.json").read_text(encoding="utf-8"))
    assert a["release_gate"] == "ok"
    ac = a["acceptance_checks"]
    assert len(ac) > 100 and all(v == "pass" for v in ac.values()), [k for k, v in ac.items() if v != "pass"]
    man = json.loads((out / "stage06_08_manifest.json").read_text(encoding="utf-8"))
    for f, d in man["files"].items():
        assert M.sha256_file(out / f) == d["sha256"], f
    assert set(man["files"]) | {"stage06_08_manifest.json"} == names


def test_audit_has_no_patient_counts(built):
    a = json.loads((built["out"] / "stage06_08_audit.json").read_text(encoding="utf-8"))
    assert M.audit_ints_ok(a)
    txt = json.dumps(a)
    for bad in ("n_train_subjects", "n_patients", "prevalence", "n_rows", "eligible_n"):
        assert bad not in txt


def test_no_banned_or_id_columns(built):
    out = built["out"]
    for p in list(out.glob("layer_*.csv")):
        cols = pd.read_csv(p, nrows=0).columns.tolist()
        feats = [c for c in cols if c not in ("encounter_id", "split")]
        assert not any(M.banned(c) for c in feats), p.name
        assert not any(c.lower() in ("subject_id", "hadm_id", "died", "label") for c in feats)
        assert not any("weight" in c for c in feats)


def test_stage08_train_mean0_sd1_and_no_empty(built):
    for p in built["out"].glob("stage08_[ab]*.npz"):
        z = np.load(p, allow_pickle=False)
        tr = z["split"] == "train"
        for key in ("X_z", "X_full_z"):
            X = z[key]
            assert np.isfinite(X).all(), p.name
            assert np.allclose(X[tr].astype(np.float64).mean(0), 0, atol=1e-5), p.name  # stored as float32
            sd = X[tr].astype(np.float64).std(0, ddof=0)
            assert np.allclose(sd, 1, atol=1e-4), p.name
        assert np.isfinite(z["X_imputed"]).all()


def test_layer_c_standardized_on_observed_train(built):
    for L in (24, 48, 336):
        z = np.load(built["out"] / f"stage08_layer_c_L{L}.npz", allow_pickle=False)
        assert z["X_raw_z"].shape[1] == L // 6
        assert np.isfinite(z["X_filled_z"]).all()
        miss = z["missing"].astype(bool)
        assert np.isnan(z["X_raw_z"][miss]).all() and np.isfinite(z["X_raw_z"][~miss]).all()
        tr = z["split"] == "train"
        X = z["X_raw_z"][tr]
        mu = np.nanmean(X.reshape(-1, X.shape[-1]), axis=0)
        assert np.all(np.abs(mu[np.isfinite(mu)]) < 1e-4)


def test_b2d_rank_deficiency_equals_kept_duplicates(built):
    a = json.loads((built["out"] / "stage06_08_audit.json").read_text(encoding="utf-8"))
    assert a["acceptance_checks"]["08_B10_dup_rank_deficiency"] == "pass"


def test_rows_and_order_identical_across_levels(built):
    ids = None
    for p in built["out"].glob("layer_*.csv"):
        e = pd.read_csv(p, usecols=["encounter_id"], dtype=str)["encounter_id"].tolist()
        ids = ids or e
        assert e == ids, p.name
    for p in built["out"].glob("stage08_[ab]*.npz"):
        assert list(np.load(p, allow_pickle=False)["encounter_id"]) == ids


# --------------------------------------------------------------------------------------------------- leakage
def test_test_rows_do_not_change_fitted_transforms(fx, tmp_path, monkeypatch):
    """Perturb every numeric value of validation and test rows after the inputs are read: every fitted transform
    (vocabularies, age-band references, winsor limits, medians, means, SDs, Layer C stats) must be unchanged."""
    real_load = M.load_inputs

    def fitted_of(out):
        f = json.loads((out / "stage08_fitted_transforms_INTERNAL.json").read_text(encoding="utf-8"))
        f.pop("generated_utc", None)
        return f

    def run_with(perturb, out):
        def load(cfg, man):
            D = real_load(cfg, man)
            if perturb:
                held = D["split"] != "train"
                meas = ("_first", "_last", "_min", "_max", "_median")
                for key in ("A", "B"):
                    T = D[key]
                    for c in T.columns:
                        if not ((key == "B" and c.endswith(meas)) or
                                (key == "A" and c.startswith("pre6_") and not c.endswith("_missing"))):
                            continue
                        v = pd.to_numeric(T[c], errors="coerce")
                        if v.nunique() <= 5:          # qualitative / ordinal results keep their codes
                            continue
                        ok = v.notna().to_numpy() & held
                        if ok.any():
                            T.loc[ok, c] = (v[ok] * 1.5 + 0.25).map(lambda x: f"{x:.6g}")
            return D
        monkeypatch.setattr(M, "load_inputs", load)
        rc, txt = run(fx, out, "--expect-stage3-manifest-sha256", fx["man"])
        monkeypatch.setattr(M, "load_inputs", real_load)
        return rc, txt

    rc1, t1 = run_with(False, tmp_path / "a")
    rc2, t2 = run_with(True, tmp_path / "b")
    assert rc1 == 0, t1
    assert rc2 == 0, t2[-3000:]
    f1, f2 = fitted_of(tmp_path / "a"), fitted_of(tmp_path / "b")
    assert json.dumps(f1, sort_keys=True) == json.dumps(f2, sort_keys=True)
    za = np.load(tmp_path / "a" / "stage08_b1.npz")
    zb = np.load(tmp_path / "b" / "stage08_b1.npz")
    held = za["split"] != "train"
    assert np.array_equal(za["X_z"][~held], zb["X_z"][~held])     # train rows identical
    assert not np.array_equal(za["X_z"][held], zb["X_z"][held])   # the perturbation reached held-out rows


# --------------------------------------------------------------------------------------------------- units
def test_admission_source_d1_rule():
    assert M.ADM_SRC_SURGICAL_WINDOW_H == 48.0 and M.ADM_SRC_FLOOR_CUTOFF_H == 24.0
    cases = [  # (interval h, recent surgery ended <= t0, expected)
        (30.0, 1, "surgical"), (30.0, 0, "floor"), (10.0, 0, "transfer_ed"), (60.0, 1, "floor"),
        (47.999, 1, "surgical"), (48.0, 1, "floor"),           # surgical window: interval < 48 h
        (24.0, 0, "floor"), (23.999, 0, "transfer_ed"),        # floor cutoff: interval >= 24 h
        (24.0, 1, "surgical"), (10.0, 1, "surgical"), (0.0, 0, "transfer_ed"), (48.0, 0, "floor"),
        (np.nan, 1, "surgical"), (np.nan, 0, "transfer_ed"),   # missing interval never counts as long
    ]
    iv = np.array([c[0] for c in cases])
    sc = np.array([c[1] for c in cases])
    assert M.adm_source(iv, sc).tolist() == [c[2] for c in cases]


def test_admission_source_block_in_a1(built):
    out = built["out"]
    cols = pd.read_csv(out / "layer_a_a1.csv", nrows=0).columns.tolist()
    blk = [c for c in cols if c.startswith("a1_admission_source__")]
    assert sorted(blk) == sorted(f"a1_admission_source__{v}" for v in ("surgical", "floor", "transfer_ed"))
    a = json.loads((out / "stage06_08_audit.json").read_text(encoding="utf-8"))
    ac = a["acceptance_checks"]
    assert ac["06_A12_admission_source_partition"] == "pass"
    assert ac["06_A13_admission_source_not_identical_to_interval_band"] == "pass"
    d = a["decisions_applied"]
    assert d["adm_src_surgical_window_hours"] == 48.0 and d["adm_src_floor_cutoff_hours"] == 24.0
    assert "adm_src_x_hours" not in json.dumps(a)
    df = pd.read_csv(out / "layer_a_a1.csv", usecols=blk + [c for c in cols if c.startswith("a1_admit_interval_band__")])
    assert (df[blk].sum(axis=1) == 1).all()
    for c in blk:
        for b in [x for x in df.columns if x.startswith("a1_admit_interval_band__")]:
            assert not df[c].equals(df[b]) and not df[c].equals(1 - df[b])


def test_a13_duplicate_check_detects_a_copy_of_the_interval_band():
    iv = np.array([1.0, 30.0, 50.0, 200.0, 10.0, 30.0])
    bands = [(iv < 24).astype(float), ((iv >= 24) & (iv < 72)).astype(float), (iv >= 72).astype(float)]

    def block(sc):
        src = M.adm_source(iv, np.asarray(sc))
        return np.column_stack([(src == v).astype(float) for v in M.ADM_SRC_LEVELS])

    assert not M.adm_src_not_duplicate(block([0] * 6), bands, iv)     # no surgery: floor == interval >= 24 h
    assert M.adm_src_not_duplicate(block([0, 1, 0, 0, 1, 0]), bands, iv)


def test_psofa_boundaries():
    age = np.array([0.5, 6, 18, 30, 100, 180, 300, np.nan])
    one = np.ones(len(age))
    r = M.psofa(age, one * 150, one * 1.19, one * 200, one * 0, one * 0)
    assert (r["coagulation"] == 0).all() and (r["hepatic"] == 0).all()
    assert np.isnan(r["renal"][-1]) and np.isnan(r["cardiovascular"][-1])
    r = M.psofa(age[:1], [149], [1.2], [45.9], [0.8], [0])
    assert r["coagulation"][0] == 1 and r["hepatic"][0] == 1 and r["cardiovascular"][0] == 1 and r["renal"][0] == 1
    r = M.psofa([6], [19.9], [12.0], [np.nan], [1.2], [1])
    assert r["coagulation"][0] == 4 and r["hepatic"][0] == 4 and r["cardiovascular"][0] == 2 and r["renal"][0] == 4
    r = M.psofa([100, 100, 100, 100], [np.nan] * 4, [np.nan] * 4, [np.nan] * 4, [0.69, 1.09, 1.79, 2.6], [0] * 4)
    assert (r["coagulation"] == 0).all() and r["renal"].tolist() == [0, 1, 2, 4]   # 60-143 mo: 0.7/1.1/1.8/2.6
    r = M.psofa([100], [50], [6.0], [64.9], [1.8], [0])
    assert r["coagulation"][0] == 2 and r["hepatic"][0] == 3 and r["cardiovascular"][0] == 1 and r["renal"][0] == 3


def test_prism3_boundaries():
    nan = [np.nan]

    def p(band, **kw):
        args = dict(sbp_min=nan, hr_max=nan, t_min=nan, t_max=nan, ph_min=nan, ph_max=nan, hco3_min=nan,
                    hco3_max=nan, po2_min=nan, pco2_max=nan, glu_max=nan, k_max=nan, cr_max=nan, bun_max=nan,
                    wbc_min=nan, plt_min=nan, ptt_max=nan)
        args.update({k: [v] for k, v in kw.items()})
        return M.prism3([band], **args)[0]

    assert p("school") == 0
    assert p("age_invalid") != p("age_invalid")      # NaN
    assert p("neonate", sbp_min=39.9) == 7 and p("neonate", sbp_min=40) == 3 and p("neonate", sbp_min=55.1) == 0
    assert p("adolescent", hr_max=155.1) == 4 and p("adolescent", hr_max=145) == 3 and p("adolescent", hr_max=144) == 0
    assert p("infant", t_min=32.9) == 3 and p("infant", t_max=40.1) == 3 and p("infant", t_max=40) == 0
    assert p("school", ph_min=6.99) == 6 and p("school", ph_min=7.28) == 2 and p("school", hco3_min=4.9) == 6
    assert p("school", hco3_max=34.1) == 4 and p("school", ph_max=7.56) == 3 and p("school", ph_max=7.48) == 2
    assert p("school", po2_min=41.9) == 6 and p("school", po2_min=42) == 3 and p("school", pco2_max=75.1) == 3
    assert p("school", pco2_max=50) == 1 and p("school", glu_max=200.1) == 2 and p("school", k_max=6.95) == 3
    assert p("neonate", cr_max=0.86) == 2 and p("adolescent", cr_max=1.3) == 0 and p("neonate", bun_max=12) == 3
    assert p("school", bun_max=14.9) == 0 and p("school", wbc_min=2.9) == 4
    assert p("school", plt_min=49) == 5 and p("school", plt_min=99) == 4 and p("school", plt_min=200) == 2
    assert p("neonate", ptt_max=85.1) == 3 and p("infant", ptt_max=57.1) == 3 and p("infant", ptt_max=57) == 0


def test_rare_categories_pooled_train_only():
    vals = np.array(["a"] * 25 + ["b"] * 25 + ["c"] * 5)
    train = np.array([True] * 25 + [False] * 25 + [True] * 5)
    subj = np.array([f"s{i}" for i in range(55)])
    voc = M.fit_vocab(vals, train, subj)
    assert voc == ["a"]                              # b only in held-out rows, c below 20 subjects
    assert M.pooled(vals, voc).tolist().count("other") == 30


def test_stage3_id_and_split_normalization():
    assert [M.norm_id(x) for x in ("123", " 123 ", "123.0", "'123'", "nan", "")] == ["123"] * 4 + ["", ""]
    assert M.SPLIT_LABELS["val"] == "validation" and M.SPLIT_LABELS["dev"] == "validation"


def test_tampered_harmonization_stops(fx, tmp_path):
    proj = tmp_path / "proj"
    shutil.copytree(fx["work"] / "proj", proj)
    h = proj / "data" / "interim" / "encounter_harmonization_fields.csv"
    h.write_bytes(h.read_bytes() + b"\n")
    base = [a if a != str(fx["work"] / "proj") else str(proj) for a in fx["base"]]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
        rc = M.main(base + ["--out-dir", str(tmp_path / "o"), "--precheck-only"])
    assert rc == 2 and "harmonization" in err.getvalue()


# --------------------------------------------------------------------------------------------- D2 age norms (v1.2.0)
Z90 = 1.2815515655446004
Z99 = 2.3263478740408408
N = M.Norms()
DPM = M.DAYS_PER_MONTH


def test_reference_tables_pinned_and_present():
    assert set(M.REF_FILES) == {"fleming2011_hr_centiles.csv", "fleming2011_rr_centiles.csv",
                                "roberts2020_map_percentiles.csv", "podium2022_baseline_scr.csv"}
    for name, pin in M.REF_FILES.items():
        assert len(pin) == 64 and M.sha256_file(HERE / name) == pin, name


def test_tampered_reference_table_stops(fx, tmp_path, monkeypatch):
    ref = tmp_path / "ref"
    ref.mkdir()
    for name in M.REF_FILES:
        shutil.copy(HERE / name, ref / name)
    p = ref / "podium2022_baseline_scr.csv"
    p.write_bytes(p.read_bytes().replace(b"0.57", b"0.58"))
    monkeypatch.setattr(M, "REF_DIR", ref)
    rc, txt = run(fx, tmp_path / "o", "--precheck-only")
    assert rc == 2 and "podium2022_baseline_scr.csv" in txt


def test_fleming_z_median_and_p90():
    age = np.array([2.5 * 365.25])                      # HR 2-3 y: p10 92, p50 110, p90 128 (symmetric)
    assert N.fleming_z("heart_rate", [110.0], age)[0] == pytest.approx(0.0)
    assert N.fleming_z("heart_rate", [128.0], age)[0] == pytest.approx(Z90)
    assert N.fleming_z("heart_rate", [92.0], age)[0] == pytest.approx(-Z90)
    mu, sd = N.fleming_params("respiratory_rate", np.array([10 * 365.25]))   # RR 8-12 y: p10 16, p50 19, p90 22
    assert mu[0] == 19 and sd[0] == pytest.approx((22 - 16) / (2 * Z90))


def test_fleming_birth_rule_and_cover():
    ages = np.array([0.0, 0.5, 1.0, 45.0, 216 * DPM - 1, 216 * DPM, np.nan, -1.0])
    mu, _ = N.fleming_params("heart_rate", ages)
    assert mu[:5].tolist() == [127, 127, 143, 143, 73]   # 'Birth' only for age < 1 day; then 0-3 m; 15-18 y
    assert np.isnan(mu[5:]).all()                         # >= 18 y, missing, negative -> empty
    mu_r, _ = N.fleming_params("respiratory_rate", ages)
    assert mu_r[:4].tolist() == [43, 43, 43, 43]           # RR at birth uses the 0-3 m row
    assert np.isnan(mu_r[5:]).all()


def test_roberts_z_median_p90_and_cap():
    a = np.full(4, 60.0)                                  # an evaluation age: 1-3 months row
    z = N.roberts_z([59.0, 76.0, 20.0, 150.0], a)        # p50 59, p90 76, below p01, above p99
    assert z[0] == pytest.approx(0.0) and z[1] == pytest.approx(Z90)
    assert z[2] == pytest.approx(-Z99) and z[3] == pytest.approx(Z99)


def test_roberts_interpolation_by_age_and_cover():
    cen = N.roberts_centiles(np.array([97.5, 0.0, 10.0, 6400.0, 6575.0, 6576.0, np.nan]))
    j50 = N.rob_cols.index("p50")
    assert cen[0, j50] == pytest.approx((59 + 62) / 2)    # midway between the 60 d (59) and 135 d (62) rows
    assert cen[1, j50] == 56 and cen[2, j50] == 56        # before the first evaluation age (15 d): first row
    assert cen[3, j50] == 76 and cen[4, j50] == 76        # after the last evaluation age up to 6575 d: last row
    assert np.isnan(cen[5]).all() and np.isnan(cen[6]).all()
    assert N.roberts_z([60.5], [97.5])[0] == pytest.approx(0.0)
    assert np.isnan(N.roberts_z([60.0], [np.nan])[0])


def test_podium_sex_and_age_lookup():
    am = np.array([0.5, 0.5, 13.0, 150.0, 150.0, 227.9, 228.0, 50.0, np.nan])
    sx = ["M", "F", "m", "M", "F", "F", "M", "U", "M"]
    b = N.podium_baseline(am, sx)
    assert b[:6].tolist() == [0.57, 0.62, 0.32, 0.65, 0.58, 0.58]
    assert np.isnan(b[6:]).all()                          # >= 19 y, unknown sex, missing age -> empty
    assert N.creatinine_ratio([1.3], [150.0], ["M"])[0] == pytest.approx(2.0)


def test_creatinine_unit_conversion_and_gate():
    assert M.creatinine_mg_dl([88.4, 44.2], "umol/L").tolist() == pytest.approx([1.0, 0.5])
    assert M.creatinine_mg_dl([88.4], "\u03bcmol/L")[0] == pytest.approx(1.0)
    assert M.creatinine_mg_dl([0.4], "mg/dL")[0] == pytest.approx(0.4)
    with pytest.raises(M.InputError):
        M.creatinine_mg_dl([1.0], "mmol/L")
    assert M.STAGE3_CREATININE_UNIT == "mg/dL"
    D = dict(split=np.array(["train"] * 3), A=pd.DataFrame({"pre6_creatinine": ["0.3", "0.4", ""]}),
             B=pd.DataFrame({"creatinine_first": ["0.3", "0.5", "0.4"]}))
    M.creatinine_unit_gate(D)                             # mg/dL passes
    D["B"] = pd.DataFrame({"creatinine_first": ["35", "40", "30"]})
    with pytest.raises(M.InputError, match="mg/dL"):
        M.creatinine_unit_gate(D)                         # umol/L-looking values stop the run


def test_haque_zaritsky_sensitivity_only():
    z = M.haque_zaritsky_z([1.5 * 5 + 55, 1.5 * 5 + 40, 60.0], [5 * 365.25, 5 * 365.25, 0.5 * 365.25])
    assert z[0] == pytest.approx(0.0) and z[1] == pytest.approx(-1.6448536, abs=1e-6) and np.isnan(z[2])


def test_d2_columns_in_levels(built):
    out = built["out"]
    a1 = pd.read_csv(out / "layer_a_a1.csv", nrows=0).columns.tolist()
    a2 = pd.read_csv(out / "layer_a_a2.csv", nrows=0).columns.tolist()
    b1 = pd.read_csv(out / "layer_b_b1.csv", nrows=0).columns.tolist()
    assert "a1_pre6_creatinine_ratio_podium" in a1 and "a1_pre6_creatinine_z" not in a1
    assert "a2_pre6_renal__creatinine_ratio_podium" in a2
    assert {"b1_creatinine_first_ratio_podium", "b1_creatinine_worst_ratio_podium"} <= set(b1)
    assert not any(c.startswith("b1_creatinine_") and c.endswith("_z") for c in b1)
    for v in ("heart_rate", "respiratory_rate", "map"):
        assert f"b1_{v}_agez_missing" in b1 and f"b1_{v}_median_z" in b1
    a = json.loads((out / "stage06_08_audit.json").read_text(encoding="utf-8"))
    assert a["reference_tables_sha256"] == M.REF_FILES
    assert "pearson_r" not in json.dumps(a)              # the Haque-Zaritsky sensitivity stays team-internal
    f = json.loads((out / "stage08_fitted_transforms_INTERNAL.json").read_text(encoding="utf-8"))
    assert "map_sensitivity_haque_zaritsky_INTERNAL" in f


def test_d2_values_match_reference_on_fixture(built, fx):
    """Recompute the published z / ratio for the fixture rows from the Stage 3 base and compare with b1."""
    A = pd.read_csv(fx["s3"] / "stage3_layer_a_base.csv", dtype=str, keep_default_na=False)
    B = pd.read_csv(fx["s3"] / "stage3_layer_b_base.csv", dtype=str, keep_default_na=False)
    b1 = pd.read_csv(built["out"] / "layer_b_b1.csv")
    inv = A["age_band"].to_numpy() == "age_invalid"
    ad = pd.to_numeric(A["age_days"], errors="coerce").to_numpy().copy()
    ad[inv] = np.nan
    am = pd.to_numeric(A["age_months"], errors="coerce").to_numpy().copy()
    am[inv] = np.nan
    hr = N.fleming_z("heart_rate", pd.to_numeric(B["heart_rate_median"], errors="coerce").to_numpy(), ad)
    mp = N.roberts_z(pd.to_numeric(B["map_min"], errors="coerce").to_numpy(), ad)
    cr = N.creatinine_ratio(pd.to_numeric(B["creatinine_first"], errors="coerce").to_numpy(), am, A["sex"])
    for exp, col in ((hr, "b1_heart_rate_median_z"), (mp, "b1_map_min_z"), (cr, "b1_creatinine_first_ratio_podium")):
        got = b1[col].to_numpy()
        assert np.array_equal(np.isnan(exp), np.isnan(got)), col
        assert np.allclose(exp[~np.isnan(exp)], got[~np.isnan(got)], rtol=1e-5, atol=1e-5), col
    assert ((b1["b1_map_agez_missing"] == 1) == b1["b1_map_median_z"].isna()).all()


# --------------------------------------------------------------------------------------------------- v1.2.1
def _run_patched(fx, out, monkeypatch, edit):
    """Full run with the Stage 3 inputs edited in memory after they are read (synthetic fixture only)."""
    real_load = M.load_inputs

    def load(cfg, man):
        D = real_load(cfg, man)
        edit(D)
        return D
    monkeypatch.setattr(M, "load_inputs", load)
    try:
        return run(fx, out, "--expect-stage3-manifest-sha256", fx["man"])
    finally:
        monkeypatch.setattr(M, "load_inputs", real_load)


def _blank_categories(D):
    A = D["A"]
    n = len(A)
    idx = np.arange(n)
    every = lambda k, off: (idx % k) == off  # noqa: E731  (spread over train, validation and test rows)
    A.loc[every(7, 1), "admit_interval_band"] = ""
    A.loc[every(11, 2), "sex"] = ""
    A.loc[every(13, 3), "first_careunit"] = ""
    rs = pd.to_numeric(A["recent_surgery"], errors="coerce").fillna(0).to_numpy() == 1
    assert rs.sum() >= 4
    A.loc[rs & every(3, 0), "recent_body_system"] = ""
    A.loc[~rs & every(5, 4), "recent_body_system"] = ""      # inactive rows: blank but not a missing value


@pytest.fixture(scope="session")
def blank_run(fx, tmp_path_factory):
    mp = pytest.MonkeyPatch()
    out = tmp_path_factory.mktemp("blank") / "stage06_08"
    try:
        rc, txt = _run_patched(fx, out, mp, _blank_categories)
    finally:
        mp.undo()
    return dict(rc=rc, txt=txt, out=out)


def test_v121_blank_categories_get_missing_masks_and_a6_passes(blank_run):
    assert blank_run["rc"] == 0, blank_run["txt"][-3000:]
    a = json.loads((blank_run["out"] / "stage06_08_audit.json").read_text(encoding="utf-8"))
    ac = a["acceptance_checks"]
    assert ac["06_A6_blocks_sum_to_one"] == "pass" and ac["06_A4_fixed_fields_identical"] == "pass"
    assert all(v == "pass" for v in ac.values())
    bm = a["blank_missing_masks"]
    for k in ("a0", "a1", "a2"):
        assert f"{k}_first_careunit_missing" in bm[k]
        assert f"{k}_sex_missing" not in bm[k]                      # D35: no sex_missing column
        assert f"{k}_sex" in a["blank_exempt_blocks"][k]
    for k in ("a1", "a2"):
        assert f"{k}_admit_interval_missing" in bm[k] and f"{k}_recent_body_system_missing" in bm[k]
    cols = {k: pd.read_csv(blank_run["out"] / f"layer_a_{k}.csv", nrows=0).columns.tolist() for k in ("a0", "a1", "a2")}
    for k in ("a1", "a2"):
        assert cols[k].count(f"{k}_admit_interval_missing") == 1      # existing mask reused, never duplicated
        assert f"{k}_admit_interval_band_missing" not in cols[k]
    assert M.audit_ints_ok(a)
    md = (blank_run["out"] / "stage06_08_audit.md").read_text(encoding="utf-8")
    assert "a1_recent_body_system_missing" in md


def test_v121_blank_masks_identical_and_cover_active_rows(blank_run):
    out = blank_run["out"]
    fr = {k: pd.read_csv(out / f"layer_a_{k}.csv") for k in ("a0", "a1", "a2")}   # synthetic fixture only
    for blk in ("first_careunit", "sex"):
        ref = fr["a0"][[c for c in fr["a0"].columns if c.startswith(f"a0_{blk}")]].to_numpy()
        for k in ("a1", "a2"):
            assert np.array_equal(ref, fr[k][[c for c in fr[k].columns if c.startswith(f"{k}_{blk}")]].to_numpy())
    a1 = fr["a1"]
    blk = a1[[c for c in a1.columns if c.startswith("a1_recent_body_system__")]].sum(axis=1).to_numpy()
    mk = a1["a1_recent_body_system_missing"].to_numpy()
    assert np.array_equal(blk + mk, a1["a1_early_surgery_flag"].to_numpy())
    assert (mk[a1["a1_early_surgery_flag"].to_numpy() == 0] == 0).all()


def test_v121_onehot_missing_mask_unit():
    L = M.Level("a9", "a9_")
    L.onehot("unit", np.array(["x", "", "y", "x"], dtype=object), ["x", "y"])
    assert L.block_missing == {"a9_unit": "a9_unit_missing"}
    assert L.cols["a9_unit_missing"].tolist() == [0, 1, 0, 0]
    L.onehot("recent_sys", np.array(["", "", "y", "x"], dtype=object), ["x", "y"], active=np.array([0, 1, 1, 1.0]))
    assert L.cols["a9_recent_sys_missing"].tolist() == [0, 1, 0, 0]
    L.onehot("full", np.array(["x", "y"], dtype=object), ["x", "y"])
    assert "a9_full_missing" not in L.cols and "a9_full" not in L.block_missing
    L.onehot("sex", np.array(["F", ""], dtype=object), ["F"])
    assert "a9_sex_missing" not in L.cols and L.block_exempt == ["a9_sex"]


def _b0(cols):
    L = M.Level("b0", "b0_")
    for nm, cat in cols:
        L.add(nm, np.zeros(3), "bin" if cat is not None else "cont", system="drug" if cat else None)
        if cat is not None:
            L.meta["b0_" + nm]["category"] = cat
    return L


def test_v121_b8_drug_names_with_tab_and_mab_pass():
    L = _b0([("drug__amoxicillin_tab_500mg_1a2b3c4d", "amoxicillin tab 500mg"),
             ("drug__rituximab_inj_5e6f7a8b", "rituximab inj"), ("drug__other", "other"),
             ("heart_rate_first", None), ("micro_culture_count", None)])
    assert M.b8_offenders(L) == []
    assert M.b8_category_blocks(L) == ["b0_drug"]


def test_v121_b8_real_banned_stems_still_fail():
    L = _b0([("org_name", None), ("micro_ab_susceptibility", None), ("mic_dilution", None),
             ("interpretation_flag", None), ("lab_value", None), ("drug__org_x_1a2b3c4d", "org x")])
    assert M.b8_offenders(L) == ["b0_org_name", "b0_micro_ab_susceptibility", "b0_mic_dilution",
                                 "b0_interpretation_flag"]
    L2 = M.Level("b2s", "b2s_")
    L2.add("infectious__org_name", np.zeros(3), "cont", system="infectious")
    assert M.b8_offenders(L2) == ["b2s_infectious__org_name"]
    L3 = M.Level("b9", "b9_")
    L3.onehot("org_species", np.array(["a", "b"], dtype=object), ["a", "b"])
    assert M.b8_offenders(L3) == ["b9_org_species__a", "b9_org_species__b"]


def test_v121_b8_full_run_drug_names_pass_and_org_field_fails(fx, tmp_path, monkeypatch):
    def rename(D):
        k = D["meds"]["key"].value_counts().index[:2].tolist()
        D["meds"]["key"] = D["meds"]["key"].replace({k[0]: "Paracetamol tab_ 500 mg", k[1]: "Rituximab_ inj"})
    rc, txt = _run_patched(fx, tmp_path / "ok", monkeypatch, rename)
    assert rc == 0, txt[-3000:]
    a = json.loads((tmp_path / "ok" / "stage06_08_audit.json").read_text(encoding="utf-8"))
    assert a["acceptance_checks"]["07_B8_m1_only_no_organism_fields"] == "pass"
    assert "b0_drug" in a["b8_category_blocks_with_banned_token"]
    real_b = M.build_layer_b

    def with_org(*args, **kw):
        r = real_b(*args, **kw)
        L0 = r["levels"]["b0"]
        L0.add("org_name", np.zeros(len(next(iter(L0.cols.values())))), "bin", system="infectious")
        r["details"]["B8_m1_only_no_organism_fields"] = [f"{k}:{c}" for k, L in r["levels"].items()
                                                        for c in M.b8_offenders(L)]
        r["checks"]["B8_m1_only_no_organism_fields"] = not r["details"]["B8_m1_only_no_organism_fields"]
        return r
    monkeypatch.setattr(M, "build_layer_b", with_org)
    rc, txt = run(fx, tmp_path / "bad", "--expect-stage3-manifest-sha256", fx["man"])
    assert rc == 3 and "07_B8_m1_only_no_organism_fields -> b0:b0_org_name" in txt
    assert not (tmp_path / "bad").exists()


def test_v121_failure_message_names_level_and_block(fx, tmp_path, monkeypatch):
    def bad_dow(D):
        A = D["A"]
        A.loc[np.arange(len(A)) % 9 == 4, "intime_dow"] = "9"      # not blank, not a weekday: an uncovered row
    rc, txt = _run_patched(fx, tmp_path / "o", monkeypatch, bad_dow)
    assert rc == 3, txt[-2000:]
    line = [ln for ln in txt.splitlines() if ln.startswith("BUILD FAILED")][0]
    assert "06_A6_blocks_sum_to_one -> a0:a0_intime_dow" in line
    assert "result: FAILED; nothing left behind" in txt and not (tmp_path / "o").exists()
    msg = line.split("details:", 1)[1]
    assert not any(tok.isdigit() for tok in msg.replace(",", " ").replace("|", " ").split())   # names only


def test_v121_failure_message_format():
    m = M.failure_message(["06_A6_blocks_sum_to_one", "08_x"], {"06_A6_blocks_sum_to_one": ["a1:a1_sex"]})
    assert "06_A6_blocks_sum_to_one -> a1:a1_sex" in m and "08_x -> no column-level detail" in m
