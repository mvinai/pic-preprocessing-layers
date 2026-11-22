"""pytest suite for 03_build_real_layers.py (Stage 3 builder v1.0.4).

Every test runs on the small PIC-format fixture written by stage3_fixture_pic.py (synthetic values only; Q3).
No test reads PIC data, the team synthetic dataset, the real Stage 2 release, the real constants file or any
_raw_LOCAL file; the Stage 2 v1.0.3-format release used by the cross-release tests is synthetic
(stage3_fixture_pic.make_stage2_release; protected by the pinned module v1.5 as Stage 2 v1.0.3 was gated) and so is
the constants file (stage3_fixture_pic.make_stage2_constants).
Run from this folder:  python -m pytest -q test_03_build_real_layers.py
"""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stage3_fixture_pic as FX  # noqa: E402
import stage3_joint_check as JC  # noqa: E402

_spec = importlib.util.spec_from_file_location("stage3_builder_03", HERE / "03_build_real_layers.py")
B = importlib.util.module_from_spec(_spec)
sys.modules["stage3_builder_03"] = B
_spec.loader.exec_module(B)

MAP_V23 = "7541d7ea67aa9631c984a5b2c7bb63bf17e6e38b9013a7b09ea42fe754536ec2"
MAP_V22 = "1e7c6c98"   # prefix only, for the refused-list check


def make_cfg(fx: dict, out: Path, **over) -> "B.Config":
    pins = dict(B.PINS)
    pins.update(fx["pins"])
    e = fx["expect"]
    ex = B.Expect(encounters=e["encounters"], subjects=e["subjects"], died=e["died"], splits=e["splits"],
                  eligible=e["eligible"], layer_c_surgeries=None)
    root, interim = Path(fx["root"]), Path(fx["interim"])
    kw = dict(raw_dir=Path(fx["raw"]), canonical_path=interim / "canonical_icu_encounters.csv",
              harmonization_path=interim / "encounter_harmonization_fields.csv", split_path=root / "subject_splits.csv",
              complexity_map_path=root / "config" / "surgery_complexity_map_v1.csv",
              phi1_medians_path=interim / "phi1_band_support_frozen_medians_MAP_v1.csv", maps_dir=interim / "maps",
              v17_module_path=HERE / "audit_surgery_complexity_duration.py",
              suppress_module_path=HERE / "stage2_v1_0_3_suppress.py", links_path=HERE / "stage3_links_v3.csv",
              out_dir=out, expect_d_icd_sha256=fx["d_icd"], expect_svs_sha256=fx["svs"], pins=pins,
              raw_pins=dict(fx["raw_pins"]), expect=ex, links_sha256=B.LINKS_SHA256,
              drops_path=HERE / "stage3_drops_v1.csv", drops_sha256=B.DROPS_SHA256,
              stage2_release_path=fx.get("s2"), expect_stage2_sha256=fx.get("s2_sha"),
              stage2_pin=fx.get("s2_sha") or B.STAGE2_RELEASE_SHA256,
              stage2_constants_path=fx.get("consts"), stage2_constants_sha256=fx.get("consts_sha"))
    kw.update(over)
    return B.Config(**kw)


class _Captured(Exception):
    pass


def add_stage2_release(fx: dict, folder: Path, stop_event: bool = False, **kw) -> dict:
    """Synthetic Stage 2 v1.0.3-format release consistent with the fixture's Stage 3 audit values: (1) a placeholder
    release (same structure) passes the precheck while the Stage 3 audit values are captured just before the release
    gate; (2) the consistent release is written (default: protected by the pinned module v1.5, n_ge18 known) and its
    hash becomes the test pin (stands in for the real pin); (3) a synthetic constants file (n_ge18 zeros; with
    stop_event also the invalid-DOB total and its two derived cells)."""
    kw.setdefault("protect", True)
    folder.mkdir(parents=True, exist_ok=True)
    ph = folder / "stage2_placeholder.json"
    fx = dict(fx, s2=ph, s2_sha=FX.make_stage2_release(ph))
    cap: dict = {}
    real = B.release_gate

    def grab(SM, g, *a, **k):
        cap.update(g.values())
        raise _Captured()
    B.release_gate = grab
    try:
        with pytest.raises(_Captured):
            run_build(make_cfg(fx, folder / "capture_out"))
    finally:
        B.release_gate = real
    s2 = folder / "feature_sparsity_frequency_audit_v1.0.3.json"
    sha = FX.make_stage2_release(s2, cap, fx["expect"]["splits"], **kw)
    cp = folder / "constants_synthetic.csv"
    return dict(fx, s2=s2, s3=cap, s2_sha=sha, consts=cp, consts_sha=FX.make_stage2_constants(cp, stop_event=stop_event))


def run_build(cfg) -> tuple[dict, str]:
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        res = B.build(cfg)
    return res, err.getvalue()


@pytest.fixture(scope="session")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("fx")
    fx = add_stage2_release(FX.make_project(root / "proj"), root / "s2")
    out = root / "out"
    res, err = run_build(make_cfg(fx, out))
    z = np.load(out / "stage3_layer_c_raw.npz", allow_pickle=False)
    A = pd.read_csv(out / "stage3_layer_a_base.csv", dtype=str, keep_default_na=False)
    Bb = pd.read_csv(out / "stage3_layer_b_base.csv", dtype=str, keep_default_na=False)
    eids = list(A["encounter_id"])
    return dict(fx=fx, out=out, res=res, err=err, z=z, A=A, B=Bb, eids=eids,
                idx={k: eids.index(v) for k, v in fx["planted"].items()})


@pytest.fixture(scope="session")
def built_padding(tmp_path_factory):
    root = tmp_path_factory.mktemp("fxpad")
    fx = add_stage2_release(FX.make_project(root / "proj", padding_extremes=True), root / "s2")
    out = root / "out"
    run_build(make_cfg(fx, out))
    return out


def ch(z, name):
    return list(z["channels"]).index(name)


def grp(z, name):
    return list(z["med_groups"]).index(name)


# ------------------------------------------------------------------------------------------- pins and gates
def test_pins_full_sha256_and_map_v23():
    assert B.PINS["channel_map"] == MAP_V23
    for k, v in B.PINS.items():
        assert B.is_full_hash(v), k
    for name, (sha, _rows) in B.RAW_PINS.items():
        assert B.is_full_hash(sha), name
    assert B.is_full_hash(B.LINKS_SHA256) and B.is_full_hash(B.DROPS_SHA256)
    assert B.STAGE2_RELEASE_SHA256 == "6621cc3fcc79c7d5711a56fba0510790748f0d8ddbb1e8a3ad73a974fe88aef2"   # v1.0.3
    assert B.PINS["suppress_module"] == "63d903d1ceef80a61f257a08aff4acd4a8c7fd5c2361d1e36de1c0cf603dc751"   # v1.5
    assert B.is_full_hash(B.CONSTANTS_SHA256)
    assert any(h.startswith(MAP_V22) for h in B.REFUSED_SHA256)


def test_links_file_pinned():
    assert B.sha256_file(HERE / "stage3_links_v3.csv") == B.LINKS_SHA256
    assert B.sha256_file(HERE / "stage3_drops_v1.csv") == B.DROPS_SHA256
    links = B.read_links_v2(HERE / "stage3_links_v3.csv")
    assert len(links) == 186 and {k for k, *_ in links} == {"eq", "le"}
    s2_only = [r for r in links if all(c.startswith(B.XR) for c in (r[1],) + r[2])]
    assert len(s2_only) == 76 and all(k == "eq" for k, *_ in s2_only)


def test_links_v3_is_v2_plus_stage2_internal_equalities():
    """v1.0.4: links v3 = links v2 (byte-identical rows) + the Stage 2 v1.0.3 internal equalities (L5a/L5b/L5c/L6/L7)
    not already in v2; the six v2 L5 rows are kept once, so all 76 are present."""
    body = lambda p: [ln for ln in p.read_text(encoding="utf-8").splitlines() if not ln.startswith("#")]  # noqa: E731
    v2, v3 = body(HERE / "stage3_links_v2.csv"), body(HERE / "stage3_links_v3.csv")
    assert v3[:len(v2)] == v2 and len(v3) == len(v2) + 70
    l3 = B.read_links_v2(HERE / "stage3_links_v3.csv")
    new = l3[116:]
    assert all("Stage 2 links v1.0.3 L" in n for *_, n in new)
    tags = [n.split("Stage 2 links v1.0.3 ")[1].split(":")[0] for *_, n in new]
    assert {t: tags.count(t) for t in set(tags)} == {"L5a": 18, "L5c": 15, "L6": 34, "L7": 3}
    assert sum(1 for *_, n in l3[:116] if n.startswith(("Stage 2 internal",))) == 6      # 3 L5a + 3 L5b in v2


def test_cli_defaults_links_v3_and_drops():
    cfg = B.build_config(["--expect-d-icd-sha256", "0" * 64, "--expect-svs-sha256", "0" * 64,
                          "--stage2-release", "s2.json", "--expect-stage2-sha256", "A" * 64])
    assert cfg.links_path.name == "stage3_links_v3.csv" and cfg.links_sha256 == B.LINKS_SHA256
    assert cfg.suppress_module_path.name == "stage2_v1_0_3_suppress.py"
    assert cfg.stage2_constants_path.name == "stage3_stage2_constants_TEAM_INTERNAL.csv"
    assert cfg.stage2_constants_sha256 == B.CONSTANTS_SHA256
    assert (cfg.max_new_suppressions, cfg.max_component_cells, cfg.heartbeat_seconds) == (2500, 20000, 30.0)
    assert cfg.drops_path.name == "stage3_drops_v1.csv" and cfg.drops_sha256 == B.DROPS_SHA256
    assert cfg.stage2_pin == B.STAGE2_RELEASE_SHA256 and cfg.expect_stage2_sha256 == "a" * 64


def test_hash_mismatch_fails_closed(built, tmp_path):
    fx = built["fx"]
    maps = tmp_path / "maps"
    shutil.copytree(Path(fx["interim"]) / "maps", maps)
    p = maps / B.MAP_FILES["channel_map"]
    p.write_bytes(p.read_bytes() + b"\n")
    out = tmp_path / "out"
    with pytest.raises(B.InputError, match="SHA-256 mismatch"):
        run_build(make_cfg(fx, out, maps_dir=maps))
    assert not out.exists()
    assert not list(tmp_path.glob("out.pending-*"))


def test_refused_map_hash(built, tmp_path, monkeypatch):
    fx = built["fx"]
    maps = tmp_path / "maps"
    shutil.copytree(Path(fx["interim"]) / "maps", maps)
    real = B.sha256_file(maps / B.MAP_FILES["channel_map"])
    monkeypatch.setitem(B.REFUSED_SHA256, real, "test: stands in for channel map v2.2")
    out = tmp_path / "out"
    with pytest.raises(B.InputError, match="refused input"):
        run_build(make_cfg(fx, out, maps_dir=maps))
    assert not out.exists()


def test_short_hash_rejected(built, tmp_path):
    with pytest.raises(B.InputError, match="full 64-character"):
        run_build(make_cfg(built["fx"], tmp_path / "out", expect_d_icd_sha256=built["fx"]["d_icd"][:12]))


def test_raw_local_and_never_read_names(tmp_path):
    p = tmp_path / "encounters_raw_LOCAL.csv"
    p.write_text("x\n")
    with pytest.raises(B.InputError, match="_raw_LOCAL"):
        B.gate_file("canonical", p, "0" * 64)
    for bad in B.NEVER_READ:
        q = tmp_path / f"{bad}.csv"
        q.write_text("x\n")
        with pytest.raises(B.InputError):
            B.gate_file("x", q, "0" * 64)


def test_nonempty_out_dir_refused(built, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "keep.txt").write_text("x")
    with pytest.raises(B.InputError, match="not empty"):
        run_build(make_cfg(built["fx"], out))
    assert (out / "keep.txt").read_text() == "x"


def test_precheck_only_writes_nothing(built, tmp_path):
    out = tmp_path / "out"
    res, _ = run_build(make_cfg(built["fx"], out, precheck_only=True))
    assert res["mode"] == "precheck-only"
    assert not out.exists()


def test_cli_requires_precheck_hashes():
    with pytest.raises(SystemExit) as e:
        B.build_config([])
    assert e.value.code == 2


@pytest.mark.parametrize("drop", ["--stage2-release", "--expect-stage2-sha256"])
def test_cli_requires_stage2_flags(drop):
    argv = {"--expect-d-icd-sha256": "0" * 64, "--expect-svs-sha256": "0" * 64, "--stage2-release": "s2.json",
            "--expect-stage2-sha256": "0" * 64}
    argv.pop(drop)
    with pytest.raises(SystemExit) as e:
        B.build_config([x for kv in argv.items() for x in kv])
    assert e.value.code == 2


def test_cli_real_pins_on_fixture_exit_2(built, tmp_path, capsys):
    fx = built["fx"]
    out = tmp_path / "out"
    rc = B.main(["--project-root", str(fx["root"]), "--raw-dir", str(fx["raw"]),
                 "--split-file", str(Path(fx["root"]) / "subject_splits.csv"), "--out-dir", str(out),
                 "--expect-d-icd-sha256", fx["d_icd"], "--expect-svs-sha256", fx["svs"],
                 "--stage2-release", str(fx["s2"]), "--expect-stage2-sha256", B.STAGE2_RELEASE_SHA256])
    cap = capsys.readouterr()
    assert rc == 2
    assert "nothing written" in cap.out
    assert not out.exists()


# ------------------------------------------------------------------------------------------- build result
def test_build_checks_pass(built):
    checks = built["res"]["checks"]
    assert checks and all(v in ("pass", "noted") for v in checks.values()), checks
    assert checks["feature_channels_92"] == "pass" and "feature_channels_93" not in checks
    assert checks["blood_gas_mask_15"] == "pass"
    assert len(checks) == 24                                             # 24 acceptance checks (v1.0.3)


def test_release_written_and_manifest_hashes(built):
    out = built["out"]
    man = json.loads((out / "stage3_manifest.json").read_text(encoding="utf-8"))
    assert man["release_gate"] == "ok"
    files = man["files"]
    assert len([n for n in files if n.startswith("stage3_layer_c_filled_L") and n.endswith(".npz")]) == 14
    for name, f in files.items():
        assert B.sha256_file(out / name) == f["sha256"], name
        mp = out / (name + ".meta.json")
        if name.endswith(".meta.json"):
            continue
        if name.startswith("stage3_layer") or name.endswith("_INTERNAL.csv") or name == "stage3_weight_normalizer.csv":
            meta = json.loads(mp.read_text(encoding="utf-8"))
            assert meta["sha256"] == f["sha256"] and meta["channel_map"] == "v2.3"
    assert man["input_hashes"]["channel_map"] == MAP_V23
    assert not list(out.parent.glob("out.pending-*"))


def test_stderr_progress_has_no_identifiers(built):
    err = built["err"]
    assert "E1000" not in err and "SUBJECT_ID=" not in err
    for line in err.splitlines():
        assert line.startswith("[") or line.startswith("  "), line


def test_92_feature_channels_and_blood_gas(built):
    z = built["z"]
    meta = json.loads((built["out"] / "stage3_layer_c_raw.npz.meta.json").read_text(encoding="utf-8"))
    assert len(z["channels"]) == 92 and "weight" not in list(z["channels"])   # 94 mapped; weight = normalizer only
    assert "ketones_urine" not in list(z["channels"])                          # v1.0.3: excluded
    feat = set(meta["feature_columns"])
    assert "weight" not in feat and "mask_valid" not in feat
    assert "blood_gas" in list(z["panels"])
    for f in ("stage3_layer_a_base.csv", "stage3_layer_b_base.csv"):
        m = json.loads((built["out"] / (f + ".meta.json")).read_text(encoding="utf-8"))
        assert not any(c.startswith("weight") for c in m["feature_columns"])


# ------------------------------------------------------------------------------------------- planted values
def test_heart_rate_bin_median_and_gate(built):
    z, i = built["z"], built["idx"]["P0"]
    k = ch(z, "heart_rate")
    assert z["median"][i, 0, k] == pytest.approx(110.0)          # 1003 100 & 120 (1002 90 deduplicated)
    assert np.isnan(z["median"][i, 1, k]) and z["missing"][i, 1, k] == 1   # 500 bpm out of gate -> missing


def test_creatinine_unit_and_bins(built):
    z, i = built["z"], built["idx"]["P0"]
    k = ch(z, "creatinine")
    assert z["median"][i, 3, k] == pytest.approx(1.0)            # 88.4 umol/L at 22.5 h, released 23.5 h
    assert z["median"][i, 4, k] == pytest.approx(2.0)            # 176.8 umol/L at 23.5 h, released 24.5 h


def test_locf_slow_lab_cap_and_impute(built):
    out, i = built["out"], built["idx"]["P0"]
    f = np.load(out / "stage3_layer_c_filled_L72.npz")
    r = list(f["encounter_index"]).index(i)
    k = list(f["value_channels"]).index("creatinine")
    v, cf, im = f["value_filled"][r, :, k], f["carried_forward"][r, :, k], f["imputed"][r, :, k]
    assert v.shape[0] == 12
    assert list(cf[5:9]) == [1, 1, 1, 1] and cf[9] == 0         # 4-bin slow-lab cap
    assert np.allclose(v[5:9], 2.0)
    assert list(im[9:]) == [1, 1, 1] and np.allclose(v[9:], f["train_median"][k])
    assert list(im[:3]) == [1, 1, 1]                             # before the first value: train median


def test_lactate_one_bin_cap(built):
    out, i = built["out"], built["idx"]["P0"]
    f = np.load(out / "stage3_layer_c_filled_L24.npz")
    r = list(f["encounter_index"]).index(i)
    k = list(f["value_channels"]).index("lactate")
    assert f["value_filled"][r, 1, k] == pytest.approx(2.0) and f["carried_forward"][r, 1, k] == 1
    assert f["imputed"][r, 2, k] == 1


def test_no_carry_forward_q9(built):
    out, i = built["out"], built["idx"]["P0"]
    f = np.load(out / "stage3_layer_c_filled_L24.npz")
    vc = list(f["value_channels"])
    for c in ("glucose_bedside", "sbp", "dbp", "pain_score"):
        if c in vc:
            assert int(f["carried_forward"][:, :, vc.index(c)].sum()) == 0, c


def test_landmark_specific_medians(built):
    out = built["out"]
    a = np.load(out / "stage3_layer_c_filled_L24.npz")["train_median"]
    b = np.load(out / "stage3_layer_c_filled_L336.npz")["train_median"]
    assert a.shape == b.shape and not np.allclose(a, b, equal_nan=True)


def test_map_from_pair(built):
    z, i = built["z"], built["idx"]["P0"]
    assert z["median"][i, 0, ch(z, "map")] == pytest.approx(70.0)


def test_pre_icu_block_not_seeded(built):
    A, i, z = built["A"], built["idx"]["P0"], built["z"]
    assert float(A.loc[i, "pre6_sodium"]) == pytest.approx(140.0)
    assert float(A.loc[i, "pre24_sodium"]) == pytest.approx(140.0)
    assert np.isnan(z["median"][i, 0, ch(z, "sodium")])


def test_weight_and_urine_rate(built):
    out, i, z = built["out"], built["idx"]["P0"], built["z"]
    W = pd.read_csv(out / "stage3_weight_normalizer.csv")
    row = W[W["encounter_id"] == built["eids"][i]].iloc[0]
    assert row["weight_kg"] == pytest.approx(10.0) and row["weight_source"] == "charted"
    assert z["median"][i, 1, ch(z, "urine_output")] == pytest.approx(2.0)
    assert z["median"][i, 5, ch(z, "kdigo_uo_stage")] == pytest.approx(3.0)


def test_after_out_not_used(built):
    z, i = built["z"], built["idx"]["P0"]
    assert int(z["missing"][i, :, ch(z, "temperature")].min()) == 1


def test_medication_state_machine(built):
    z, i, A = built["z"], built["idx"]["P0"], built["A"]
    g = grp(z, "vasoactive")
    assert list(z["med_active"][i, :4, g]) == [1, 1, 1, 0]
    assert z["med_started"][i, 0, g] == 1
    assert not np.any((z["med_started"] == 1) & (z["med_active"] == 0))
    assert A.loc[i, "med_antibiotic_active_t0"] == "1" and A.loc[i, "med_vasoactive_active_t0"] == "0"


def test_osmotic_resolution_and_date_only(built):
    M = pd.read_csv(built["out"] / "stage3_medication_courses_INTERNAL.csv")
    e2 = built["eids"][built["idx"]["P2"]]
    m = M[(M["encounter_id"] == e2) & M["key"].str.contains("Mannitol")]
    assert set(m["group"]) == {"diuretic"}                       # G40 is not a neuro admitting diagnosis
    th = B.assign_date_only_start(np.array([-30.0, -5.0, 0.0, 13.0, 18.0]), np.array([1, 1, 1, 1, 0], bool))
    assert list(th) == [-30.0, 0.0, 0.0, 18.0, 18.0]


def test_surgery_timing_and_layer_c(built):
    z, i, A = built["z"], built["idx"]["P0"], built["A"]
    assert A.loc[i, "recent_surgery"] == "1"
    assert list(z["surgery_in_bin"][i, :6]) == [0, 0, 0, 0, 0, 1]


def test_q12_duration_fallbacks(built):
    S = pd.read_csv(built["out"] / "stage3_surgery_long_INTERNAL.csv")
    e1 = built["eids"][built["idx"]["P1"]]
    e2 = built["eids"][built["idx"]["P2"]]
    r1 = S[S["encounter_id"] == e1].iloc[0]
    r2 = S[S["encounter_id"] == e2].iloc[0]
    assert r1["duration_source"] == "anesthesia_window" and r1["fallback_reason"] == "gt24"
    assert r2["phi1_band"] == "missing" and int(r2["band_missing_flag"]) == 1


def test_q4_clip_and_q8_unmapped(built):
    A, idx = built["A"], built["idx"]
    i = idx["P0"]
    assert A.loc[i, "admit_interval_clipped"] == "1" and float(A.loc[i, "hours_admit_to_icu_intime"]) == 0.0
    assert A.loc[i, "dx_icd10_3char"].startswith("unmapped")
    assert A.loc[idx["P1"], "dx_icd10_3char"].startswith("unmapped")
    assert A.loc[idx["P2"], "dx_icd10_3char"] == "G40"


def test_bilirubin_5255_nicu_only(built):
    Bb, idx = built["B"], built["idx"]
    assert Bb.loc[idx["P4_NICU"], "bilirubin_total_src_5255"] == "1"
    assert Bb.loc[idx["P5"], "bilirubin_total_src_5255"] == "0"


def test_missing_dob_age_invalid(built):
    A, i = built["A"], built["idx"]["P5"]
    assert A.loc[i, "age_band"] == "age_invalid"


def test_padding_extremes_change_nothing(built, built_padding):
    a, b = built["out"], built_padding
    for f in ("stage3_layer_a_base.csv", "stage3_layer_b_base.csv"):
        assert (a / f).read_bytes() == (b / f).read_bytes(), f
    for f in ["stage3_layer_c_raw.npz"] + [f"stage3_layer_c_filled_L{L}.npz" for L in range(24, 337, 24)]:
        x, y = np.load(a / f), np.load(b / f)
        assert x.files == y.files
        for m in x.files:
            assert np.array_equal(x[m], y[m], equal_nan=x[m].dtype.kind == "f"), (f, m)


# ------------------------------------------------------------------------------------------- disclosure
def _walk(x, p=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _walk(v, f"{p}/{k}")
    elif isinstance(x, list):
        for j, v in enumerate(x):
            yield from _walk(v, f"{p}[{j}]")
    else:
        yield p, x


def test_audit_has_no_small_cells(built):
    doc = json.loads((built["out"] / "stage3_audit.json").read_text(encoding="utf-8"))
    assert doc["release_gate"] == "ok"
    for p, v in _walk(doc):
        if "/n_" in p and isinstance(v, int) and not isinstance(v, bool):
            assert v >= 10, p
    md = (built["out"] / "stage3_audit.md").read_text(encoding="utf-8")
    assert B.render_markdown(doc) == md


def _sm():
    return B.load_pinned_module(HERE / "stage2_v1_0_3_suppress.py", B.PINS["suppress_module"], "suppress_module")


def test_complementary_suppression_blocks_differencing(tmp_path):
    SM = _sm()
    g = B.Agg()
    t = g.n(("x", "n_total"), 100)
    a = g.n(("x", "a", "n_rows"), 95)
    b = g.n(("x", "b", "n_rows"), 5)
    g.ident(t, [a, b])
    doc, md, _ = B.release_gate(SM, g, None, {"spec": "test"}, tmp_path)
    assert doc["x"]["b"]["n_rows"] == "<10"
    assert doc["x"]["a"]["n_rows"] == "S" or doc["x"]["n_total"] == "S"


def test_frozen_cells_cannot_protect_fails(tmp_path):
    SM = _sm()
    g = B.Agg()
    t = g.n(("x", "n_total"), 100, frozen=True)
    a = g.n(("x", "a", "n_rows"), 95, frozen=True)
    b = g.n(("x", "b", "n_rows"), 5)
    g.ident(t, [a, b])
    with pytest.raises(B.GateError, match="cannot protect"):
        B.release_gate(SM, g, None, {"spec": "test"}, tmp_path)


def test_cannot_protect_pinned_identity_v12(tmp_path):
    """v1.3 (same as v1.2 without a drop list; as in synthetic_stall_v101.py): a '<10' part of a small identity whose
    total and other part are frozen, inside a larger linked component with many publishable cells. v1.1 kept suppressing
    unrelated cells; v1.3 must stop with 'cannot protect' and write nothing."""
    SM = _sm()
    g = B.Agg()
    top = g.n(("big", "n_total"), 5000)
    parts = [g.n(("big", f"g{j:02d}", "n_rows"), 100 + j) for j in range(40)]
    rest = g.n(("big", "rest", "n_rows"), 5000 - sum(100 + j for j in range(40)))
    g.ident(top, parts + [rest])
    t = g.n(("big", "g00", "split", "n_total"), 100)
    a = g.n(("big", "g00", "split", "a", "n_rows"), 95, frozen=True)
    b = g.n(("big", "g00", "split", "b", "n_rows"), 5)
    g.frozen.add(t)
    g.ident(t, [a, b])
    g.ident(parts[0], [t])
    with pytest.raises(B.GateError, match="cannot protect"):
        B.release_gate(SM, g, None, {"spec": "test"}, tmp_path)
    assert not (tmp_path / "stage3_audit.json").exists()


@pytest.mark.parametrize("msg", ["cannot protect 4 cell(s): test", "suppression not converging: test"])
def test_v12_gate_errors_exit_3_nothing_left(built, tmp_path, capsys, monkeypatch, msg):
    """The v1.3 exit-3 messages raised inside the pinned module fail the whole Stage 3 run closed (exit 3, no
    output folder, no pending folder)."""
    real_loader = B.load_pinned_module

    def loader(path, expected, tag):
        mod = real_loader(path, expected, tag)
        if tag == "suppress":
            def run(self):
                raise mod.GateError(msg)
            monkeypatch.setattr(mod.Suppressor, "run", run)
        return mod

    monkeypatch.setattr(B, "load_pinned_module", loader)
    cfg = make_cfg(built["fx"], tmp_path / "out")
    monkeypatch.setattr(B, "build_config", lambda argv=None: cfg)
    rc = B.main([])
    cap = capsys.readouterr()
    assert rc == 3
    assert "BUILD/RELEASE GATE FAILED: suppression module: " + msg.split(":")[0] in cap.err
    assert "nothing left behind" in cap.out
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob("out.pending-*"))


def test_frozen_small_cell_fails(tmp_path):
    SM = _sm()
    g = B.Agg()
    g.n(("x", "n_rows"), 4, frozen=True)
    with pytest.raises(B.GateError, match="frozen"):
        B.release_gate(SM, g, None, {"spec": "test"}, tmp_path)


# ------------------------------------------------------------------------------------------- helpers
def test_locf_fill_does_not_chain():
    med = np.array([[1.0, np.nan, np.nan, np.nan, 5.0, np.nan]])
    v, c = B.locf_fill(med, med, 2)
    assert np.allclose(v[0, :3], 1.0) and np.isnan(v[0, 3]) and v[0, 5] == 5.0
    assert list(c[0]) == [False, True, True, False, False, True]


def test_age_band_and_icd():
    assert list(B.age_band(np.array([10, 100, 800, -1, np.nan]))) == ["neonate", "infant", "toddler", "age_invalid", "age_invalid"]
    assert B.norm_icd("J18.900") == ("J18900", "ok")
    assert B.norm_icd("")[1] == "missing"
    assert B.norm_icd("XYZ")[1] == "invalid"
    assert B.norm_icd("J18.9, K35.8")[1] == "multiple"


def test_time_limit_fails_closed(built, tmp_path):
    out = tmp_path / "out"
    with pytest.raises(B.TimeLimit):
        run_build(make_cfg(built["fx"], out, max_seconds=0.0))
    assert not out.exists()
    assert not list(tmp_path.glob("out.pending-*"))


# ------------------------------------------------------------------------------------------- Windows file semantics (v1.0.1)
# Windows refuses to delete or rename a file while this process (or cloud sync, antivirus, the indexer) holds a handle or
# a memory map on it; Linux allows both, so a plain Linux run cannot catch such bugs. WinFS patches os.replace/rename/
# unlink/remove/rmdir to raise PermissionError (errno 13, WinError 32) when the target, or anything under a target
# folder, is open or memory-mapped by this process (read from /proc), and can inject transient or persistent
# sharing violations to stand in for cloud-sync/antivirus locks. On Windows itself the OS enforces these rules natively,
# so WinFS only injects the locks there.
def _semantics() -> str:
    if Path("/proc/self/maps").exists() and Path("/proc/self/fd").is_dir():
        return "simulated"
    if os.name == "nt":
        return "native"
    pytest.skip("needs /proc (or Windows) for Windows sharing semantics")


def _perm(path) -> PermissionError:
    e = PermissionError(13, "simulated Windows sharing violation", str(path))
    e.winerror = 32
    return e


class WinFS:
    def __init__(self, monkeypatch, fail=None):
        self.mode = _semantics()
        self.fail = fail or (lambda op, path: False)
        self.blocked: list = []
        real = {k: getattr(os, k) for k in ("replace", "rename", "unlink", "remove", "rmdir")}

        def absp(path, dir_fd=None):
            path = os.fspath(path)
            if dir_fd is not None and not os.path.isabs(path):
                return os.path.join(os.readlink(f"/proc/self/fd/{dir_fd}"), path)
            return os.path.join(os.path.realpath(os.path.dirname(os.path.abspath(path))), os.path.basename(path))

        def busy() -> set:
            out = set()
            for fd in os.listdir("/proc/self/fd"):
                try:
                    out.add(os.readlink(f"/proc/self/fd/{fd}"))
                except OSError:
                    pass
            with open("/proc/self/maps") as f:
                for line in f:
                    parts = line.split(maxsplit=5)
                    if len(parts) == 6 and parts[5].startswith("/"):
                        out.add(parts[5].strip().removesuffix(" (deleted)"))
            return out

        def check(op, *targets):
            b = None
            for t in targets:
                if self.fail(op, t):
                    self.blocked.append((op, t))
                    raise _perm(t)
                if self.mode == "native":
                    continue
                b = busy() if b is None else b
                isdir = os.path.isdir(t)
                if any((x == t and not isdir) or x.startswith(t + os.sep) for x in b):
                    self.blocked.append((op, t))
                    raise _perm(t)

        def replace(src, dst, *, src_dir_fd=None, dst_dir_fd=None):
            check("replace", absp(src, src_dir_fd), absp(dst, dst_dir_fd))
            return real["replace"](src, dst, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd)

        def rename(src, dst, *, src_dir_fd=None, dst_dir_fd=None):
            check("replace", absp(src, src_dir_fd), absp(dst, dst_dir_fd))
            return real["rename"](src, dst, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd)

        def unlink(path, *, dir_fd=None):
            check("unlink", absp(path, dir_fd))
            return real["unlink"](path, dir_fd=dir_fd)

        def remove(path, *, dir_fd=None):
            check("unlink", absp(path, dir_fd))
            return real["remove"](path, dir_fd=dir_fd)

        def rmdir(path, *, dir_fd=None):
            check("rmdir", absp(path, dir_fd))
            return real["rmdir"](path, dir_fd=dir_fd)

        for k, f in (("replace", replace), ("rename", rename), ("unlink", unlink), ("remove", remove), ("rmdir", rmdir)):
            monkeypatch.setattr(os, k, f)


class Transient:
    """Fail the first n matching operations (op, path predicate), then let them through."""

    def __init__(self, rules):
        self.rules = [[op, pred, n] for op, pred, n in rules]
        self.hits = 0

    def __call__(self, op, path):
        for r in self.rules:
            if r[0] == op and r[1](path) and r[2] > 0:
                r[2] -= 1
                self.hits += 1
                return True
        return False


def test_fs_op_retries_transient_then_succeeds(monkeypatch, capsys):
    monkeypatch.setitem(B._CLOCK, "deadline", float("inf"))   # an earlier test may have set a past deadline
    monkeypatch.setattr(B, "FS_RETRY_DELAYS", (0.0, 0.0, 0.0))
    n = {"k": 0}

    def fn():
        n["k"] += 1
        if n["k"] < 3:
            raise _perm("x.npz.tmp")
        return 7

    assert B.fs_op("replace", fn, Path("x.npz.tmp")) == 7
    assert n["k"] == 3
    assert capsys.readouterr().err.count("file operation retry: replace") == 2


def test_fs_op_persistent_names_operation_errno_and_relative_path(tmp_path, monkeypatch):
    monkeypatch.setitem(B._CLOCK, "deadline", float("inf"))   # an earlier test may have set a past deadline
    monkeypatch.setattr(B, "FS_RETRY_DELAYS", (0.0, 0.0, 0.0))
    monkeypatch.setitem(B._FS_BASE, "dir", tmp_path)
    src = tmp_path / "out.pending-1" / "stage3_layer_c_raw.npz.tmp"

    def fn():
        raise _perm(src)

    with pytest.raises(B.FileOpError) as ei:
        B.fs_op("replace", fn, src, src.with_suffix(""))
    msg = str(ei.value)
    assert "PermissionError during replace: " + str(Path("out.pending-1") / "stage3_layer_c_raw.npz.tmp") in msg
    assert "errno 13" in msg and "WinError 32" in msg and "4 attempts" in msg
    assert str(tmp_path) not in msg


def test_fs_op_non_transient_not_retried(tmp_path, monkeypatch):
    monkeypatch.setitem(B._CLOCK, "deadline", float("inf"))   # an earlier test may have set a past deadline
    monkeypatch.setattr(B, "FS_RETRY_DELAYS", (0.0, 0.0, 0.0))
    n = {"k": 0}

    def fn():
        n["k"] += 1
        raise FileNotFoundError(2, "missing", "y")

    with pytest.raises(B.FileOpError, match="1 attempt over"):
        B.fs_op("replace", fn, tmp_path / "y")
    assert n["k"] == 1


def test_remove_tree_retries_locked_member(tmp_path, monkeypatch):
    monkeypatch.setitem(B._CLOCK, "deadline", float("inf"))   # an earlier test may have set a past deadline
    monkeypatch.setattr(B, "FS_RETRY_DELAYS", (0.0, 0.0, 0.0))
    d = tmp_path / "_raw_members"
    d.mkdir()
    for nm in ("median.npy", "missing.npy"):
        (d / nm).write_bytes(b"x" * 64)
    tr = Transient([("unlink", lambda p: p.endswith("median.npy"), 2), ("rmdir", lambda p: p.endswith("_raw_members"), 1)])
    WinFS(monkeypatch, tr)
    B.remove_tree(d, "remove member folder")
    assert not d.exists() and tr.hits == 3


def test_npystore_will_not_delete_a_mapped_member(tmp_path, monkeypatch):
    monkeypatch.setitem(B._CLOCK, "deadline", float("inf"))   # an earlier test may have set a past deadline
    """A memory map still referenced when its store is zipped is caught by name (old code: Windows PermissionError on
    the delete, no path); with every reference dropped the zip and the delete succeed under Windows semantics."""
    monkeypatch.setattr(B, "FS_RETRY_DELAYS", (0.0,))
    monkeypatch.setitem(B._FS_BASE, "dir", tmp_path)
    WinFS(monkeypatch)
    st = B.NpyStore(tmp_path / "pend" / "_filled_24")
    a = st.create("value_filled", (3, 4, 2), np.float32, np.nan)
    view = a[0]
    a = None
    with pytest.raises(B.FileOpError, match="release memory map") as ei:
        st.write_npz(tmp_path / "pend" / "stage3_layer_c_filled_L24.npz")
    assert str(Path("pend") / "_filled_24" / "value_filled.npy") in str(ei.value)
    del view
    st2 = B.NpyStore(tmp_path / "pend" / "_filled_48")
    b = st2.create("value_filled", (3, 4, 2), np.float32, np.nan)
    b[1, 2, 1] = 5.0
    b = None
    st2.put("train_median", np.array([1.0, 2.0]))
    st2.write_npz(tmp_path / "pend" / "stage3_layer_c_filled_L48.npz")
    assert not (tmp_path / "pend" / "_filled_48").exists()
    with np.load(tmp_path / "pend" / "stage3_layer_c_filled_L48.npz") as z:
        assert z["value_filled"][1, 2, 1] == 5.0 and np.isnan(z["value_filled"][0, 0, 0])


def test_failure_in_layer_c_leaves_nothing_under_windows_semantics(built, tmp_path, capsys, monkeypatch):
    """A persistent sharing violation while Layer C is written: exit 3, the message names the operation, errno and the
    relative path (no absolute path), and the pending folder is removed although Layer C memory maps were open when the
    error was raised (v1.0 left them mapped during cleanup, so Windows kept the pending folder)."""
    monkeypatch.setattr(B, "FS_RETRY_DELAYS", (0.0, 0.0))
    cfg = make_cfg(built["fx"], tmp_path / "out")
    monkeypatch.setattr(B, "build_config", lambda argv=None: cfg)
    WinFS(monkeypatch, lambda op, p: op == "replace" and p.endswith("stage3_layer_c_raw.npz.tmp"))
    rc = B.main([])
    cap = capsys.readouterr()
    assert rc == 3
    assert "BUILD FAILED (file operation: PermissionError during replace: out.pending-" in cap.err
    assert "stage3_layer_c_raw.npz.tmp" in cap.err and "errno 13" in cap.err and "WinError 32" in cap.err
    assert "3 attempts" in cap.err and "build_layer_c line" in cap.err
    assert str(tmp_path) not in cap.err and str(tmp_path) not in cap.out
    assert "result: FAILED; nothing left behind" in cap.out
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob("out.pending-*"))


def test_windows_semantics_build_with_transient_locks(built, tmp_path, monkeypatch):
    """Whole build under Windows semantics with cloud-sync-like transient locks on a zip replace, member deletes and
    the final folder rename: it retries, succeeds, and writes the same release as the plain build."""
    monkeypatch.setattr(B, "FS_RETRY_DELAYS", (0.0, 0.0, 0.0))
    tr = Transient([("replace", lambda p: p.endswith("stage3_layer_c_raw.npz.tmp"), 2),
                    ("unlink", lambda p: os.sep + "_raw_members" + os.sep in p, 2),
                    ("unlink", lambda p: os.sep + "_filled_336" + os.sep in p, 1),
                    ("replace", lambda p: ".pending-" in os.path.basename(p), 1),
                    ("replace", lambda p: p.endswith("stage3_layer_a_base.csv.tmp"), 1)])
    WinFS(monkeypatch, tr)
    out = tmp_path / "out"
    res, err = run_build(make_cfg(built["fx"], out))
    assert res["mode"] == "build" and tr.hits == 7
    # 7 injected locks; the first failed delete inside rmtree goes straight to the retry handler unlogged, so 5 logs
    # (on Windows itself antivirus or the indexer may add real, harmless retries)
    assert err.count("file operation retry:") >= 5 and str(tmp_path) not in err
    assert not list(tmp_path.glob("out.pending-*"))
    ref = built["out"]
    man = json.loads((out / "stage3_manifest.json").read_text(encoding="utf-8"))
    man0 = json.loads((ref / "stage3_manifest.json").read_text(encoding="utf-8"))
    assert set(man["files"]) == set(man0["files"])
    for nm in ("stage3_layer_a_base.csv", "stage3_layer_b_base.csv", "stage3_weight_normalizer.csv"):
        assert man["files"][nm]["sha256"] == man0["files"][nm]["sha256"], nm
    with np.load(out / "stage3_layer_c_raw.npz") as z1, np.load(ref / "stage3_layer_c_raw.npz") as z0:
        assert set(z1.files) == set(z0.files)
        for k in z0.files:
            assert np.array_equal(z1[k], z0[k], equal_nan=z0[k].dtype.kind == "f"), k


def test_unexpected_oserror_message_has_errno_and_relative_path(tmp_path, capsys, monkeypatch):
    monkeypatch.setitem(B._FS_BASE, "dir", tmp_path)
    monkeypatch.setattr(B, "build_config", lambda argv=None: None)

    def boom(cfg):
        raise _perm(tmp_path / "out.pending-9" / "stage3_layer_a_base.csv")

    monkeypatch.setattr(B, "build", boom)
    rc = B.main([])
    cap = capsys.readouterr()
    assert rc == 3
    assert "BUILD FAILED (unexpected PermissionError: errno 13" in cap.err
    assert "path " + str(Path("out.pending-9") / "stage3_layer_a_base.csv") in cap.err and str(tmp_path) not in cap.err


# ------------------------------------------------------------------------------------------- cross-release (v1.0.2)
# The Stage 2 release used here is synthetic (make_stage2_release): consistent with the fixture's Stage 3 audit under
# stage3_links_v2.csv, primary and simple complementary suppression applied; its hash stands in for the real pin.
DROPS = [r[0] for r in (ln.split(",", 1) for ln in (HERE / "stage3_drops_v1.csv").read_text(encoding="utf-8").splitlines()
                        if ln and not ln.startswith("#")) if r[0] != "subtree"]


def _fails_input(fx, tmp_path, match, exc=None, **over):
    out = tmp_path / "out"
    with pytest.raises(exc or B.InputError, match=match):
        run_build(make_cfg(fx, out, **over))
    assert not out.exists()
    assert not list(tmp_path.glob("out.pending-*"))


def test_stage2_hash_gate_modified_file(built, tmp_path):
    fx = built["fx"]
    bad = tmp_path / "feature_sparsity_frequency_audit_v1.0.3.json"
    bad.write_bytes(Path(fx["s2"]).read_bytes() + b"\n")
    _fails_input(fx, tmp_path, "Stage 2 release JSON SHA-256 mismatch", stage2_release_path=bad)


@pytest.mark.parametrize("case,match", [("short", "full 64-character"), ("not_pin", "not the pinned"),
                                        ("missing", "required")])
def test_stage2_hash_gate_expectations(built, tmp_path, case, match):
    fx = built["fx"]
    over = {"short": dict(expect_stage2_sha256=fx["s2_sha"][:16]),
            "not_pin": dict(stage2_pin=B.STAGE2_RELEASE_SHA256),    # the fixture file is not the pinned release
            "missing": dict(stage2_release_path=None)}[case]
    _fails_input(fx, tmp_path, match, **over)


def test_stage2_hash_gate_cli_exit_2(built, tmp_path, capsys, monkeypatch):
    cfg = make_cfg(built["fx"], tmp_path / "out", expect_stage2_sha256="0" * 64)
    monkeypatch.setattr(B, "build_config", lambda argv=None: cfg)
    rc = B.main([])
    cap = capsys.readouterr()
    assert rc == 2 and "STOPPED (input)" in cap.err and "nothing written" in cap.out
    assert not (tmp_path / "out").exists()


def test_stage2_release_must_be_gated_and_linked(built, tmp_path):
    fx = built["fx"]
    doc = json.loads(Path(fx["s2"]).read_text(encoding="utf-8"))
    doc["release_gate"] = "pending"
    p = tmp_path / "s2_pending.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    sha = B.sha256_file(p)
    _fails_input(fx, tmp_path, "release_gate 'ok'", stage2_release_path=p, expect_stage2_sha256=sha, stage2_pin=sha)
    doc["release_gate"] = "ok"
    doc["disclosure"]["count_fields"].remove("/surgery/link_classes/n_linked")
    q = tmp_path / "s2_unlinked.json"
    q.write_text(json.dumps(doc), encoding="utf-8")
    sha = B.sha256_file(q)
    _fails_input(fx, tmp_path, "not a count field of the Stage 2 release", stage2_release_path=q,
                 expect_stage2_sha256=sha, stage2_pin=sha)


def test_dropped_subtrees_absent(built):
    doc = json.loads((built["out"] / "stage3_audit.json").read_text(encoding="utf-8"))
    md = (built["out"] / "stage3_audit.md").read_text(encoding="utf-8")
    assert len(DROPS) == 20
    dis = doc["disclosure"]
    for t in DROPS:
        with pytest.raises(KeyError):
            B.json_ptr_get(doc, t)
        assert f"`{t[1:]}/" not in md, t
        assert not any(c.startswith(t + "/") for c in dis["count_fields"]), t
        assert not any(c.startswith(t + "/") for tt, ps in dis["identities"] for c in [tt] + ps), t
    assert set(doc["age"]) == {"n_valid", "n_invalid"}
    assert "n_negative" not in json.dumps(doc["age"]) and "negative_magnitude" not in md
    assert set(doc["admission_interval"]) == {"n_negative_clipped_to_0"}
    assert set(doc["weight_source"]) == {"n_charted", "n_imputed"}           # per-source split withheld (§10)
    assert set(doc["event_linkage"]) == {"LABEVENTS", "PRESCRIPTIONS"}
    for ch in ("heart_rate", "temperature", "pain_score"):
        assert set(doc["plausibility"][ch]) == {"overall"}     # overall kept for the gate lock
    assert "train" in doc["plausibility"]["platelets"]
    assert "dropped_subtrees_v1_0_4" in doc["notes"] and "cross_release_v1_0_4" in doc["notes"]
    assert not any(k.endswith(("v1_0_2", "v1_0_3")) for k in doc["notes"])
    assert "Stage 3 v1.0.4 release" in doc["notes"]["dropped_subtrees_v1_0_4"]
    assert "differencing_fix_v1_0_3" not in doc["notes"]                    # Stage 2 profile only


@pytest.fixture(scope="session")
def built_planted(built, tmp_path_factory):
    """Planted cross-release leak: in the synthetic Stage 2 release, surgery link class 'other admission' is 4 ('<10')
    and the linked count is 'S', so Stage 2 alone is safe; the Stage 3 linked-surgery count (published in a plain
    build) would pin the 4 through the declared link."""
    root = tmp_path_factory.mktemp("planted")
    fx, s3 = built["fx"], built["fx"]["s3"]
    s2 = root / "feature_sparsity_frequency_audit_v1.0.3.json"
    sha = FX.make_stage2_release(s2, s3, fx["expect"]["splits"],
                                 override={"/surgery/link_classes/n_other_hadm": 4,
                                           "/surgery/link_classes/n_total": s3["/surgery/n_linked"] + 4},
                                 statuses={"/surgery/link_classes/n_linked": "S"})
    fxp = dict(fx, s2=s2, s2_sha=sha)
    out = root / "out"
    run_build(make_cfg(fxp, out))
    return dict(fx=fxp, out=out)


def test_cross_release_leak_caught(built, built_planted, tmp_path, capsys):
    s3 = built["fx"]["s3"]
    plain = json.loads((built["out"] / "stage3_audit.json").read_text(encoding="utf-8"))
    doc = json.loads((built_planted["out"] / "stage3_audit.json").read_text(encoding="utf-8"))
    assert s3["/surgery/n_linked"] >= 10 and plain["surgery"]["n_linked"] == s3["/surgery/n_linked"]
    assert doc["surgery"]["n_linked"] == "S"                       # the joint gate suppressed it
    man = json.loads((built_planted["out"] / "stage3_manifest.json").read_text(encoding="utf-8"))
    assert "/surgery/link_classes/n_other_hadm" not in man["cross_release"]["flagged_for_stage2_v1_0_4"]
    links = HERE / "stage3_links_v3.csv"
    ap = built_planted["out"] / "stage3_audit.json"
    cs = str(built["fx"]["consts"])
    assert JC.main(str(ap), str(built_planted["fx"]["s2"]), str(links), cs) == 0
    # control: the values a build without the joint audit publishes are caught by the independent check
    for c in doc["disclosure"]["count_fields"]:
        if B.json_ptr_get(doc, c) == "S" and isinstance(B.json_ptr_get(plain, c), int):
            parts = c.split("/")[1:]
            d = doc
            for k in parts[:-1]:
                d = d[k]
            d[parts[-1]] = B.json_ptr_get(plain, c)
    ctl = tmp_path / "control_audit.json"
    ctl.write_text(json.dumps(doc), encoding="utf-8")
    capsys.readouterr()
    assert JC.main(str(ctl), str(built_planted["fx"]["s2"]), str(links), cs) == 1
    assert "NARROWED P s2:/surgery/link_classes/n_other_hadm: joint (4, 4)" in capsys.readouterr().out


def test_cross_release_inconsistent_fails_closed(built, tmp_path, capsys, monkeypatch):
    fx, s3 = built["fx"], built["fx"]["s3"]
    s2 = tmp_path / "s2_inconsistent.json"
    sha = FX.make_stage2_release(s2, s3, fx["expect"]["splits"],
                                 override={"/surgery/link_classes/n_linked": s3["/surgery/n_linked"] + 20,
                                           "/surgery/link_classes/n_total": s3["/surgery/n_linked"] + 45})
    cfg = make_cfg(dict(fx, s2=s2, s2_sha=sha), tmp_path / "out")
    monkeypatch.setattr(B, "build_config", lambda argv=None: cfg)
    rc = B.main([])
    cap = capsys.readouterr()
    assert rc == 3 and "cross-release: published Stage 2 v1.0.3 and Stage 3 cells are inconsistent" in cap.err
    assert "nothing left behind" in cap.out
    assert not (tmp_path / "out").exists() and not list(tmp_path.glob("out.pending-*"))


@pytest.fixture(scope="session")
def built_legacy(tmp_path_factory):
    """A Stage 2 release with only the simple v1.0.2-style complementary suppression (protect=False): it leaks inside
    Stage 2 alone (the invalid-age cells follow from published valid ages and the public cohort size)."""
    root = tmp_path_factory.mktemp("legacy")
    fx = add_stage2_release(FX.make_project(root / "proj"), root / "s2", protect=False)
    out = root / "out"
    run_build(make_cfg(fx, out))
    return dict(fx=fx, out=out)


def test_compromised_stage2_cells_tolerated_and_flagged(built_legacy):
    """Stage 3 cannot repair a leak inside Stage 2 alone: the build passes, the cells are flagged for Stage 2 v1.0.4,
    and the independent check finds nothing narrowed beyond what Stage 2 already shows."""
    out, fx = built_legacy["out"], built_legacy["fx"]
    man = json.loads((out / "stage3_manifest.json").read_text(encoding="utf-8"))
    doc = json.loads((out / "stage3_audit.json").read_text(encoding="utf-8"))
    xr = man["cross_release"]
    assert man["release_gate"] == "ok" and xr["stage2_release_sha256"] == fx["s2_sha"]
    assert "/layer_a/dob_breakdown/n_invalid" in xr["flagged_for_stage2_v1_0_4"]
    assert doc["release_v1_0_4"]["cross_release"]["compromised_stage2_cells"] == xr["flagged_for_stage2_v1_0_4"]
    assert "release_v1_0_3" not in doc and doc["release_v1_0_4"]["script"] == "03_build_real_layers v1.0.4"
    assert doc["release_v1_0_4"]["suppression_module"] == "stage2_v1_0_3_suppress v1.5"
    assert "/layer_a/dob_breakdown/n_invalid" in doc["notes"]["cross_release_v1_0_4"]
    h = man["input_hashes"]
    assert h["stage2_release"] == fx["s2_sha"] and h["links"] == B.LINKS_SHA256 and h["drops"] == B.DROPS_SHA256
    assert JC.main(str(out / "stage3_audit.json"), str(fx["s2"]), str(HERE / "stage3_links_v3.csv"),
                   str(fx["consts"])) == 0


def test_joint_with_v1_0_3_style_release_is_clean(built, capsys):
    """v1.0.4 acceptance (synthetic): with a Stage 2 v1.0.3-style release (protected by module v1.5 under its own
    identities + internal links, n_ge18 known, table sizes public), links v3 and the constants file, the joint gate
    passes with nothing compromised and no hidden difference, and the independent check finds 0 / 0."""
    man = json.loads((built["out"] / "stage3_manifest.json").read_text(encoding="utf-8"))
    xr = man["cross_release"]
    assert man["release_gate"] == "ok"
    assert xr["flagged_for_stage2_v1_0_4"] == []
    assert xr["flagged_hidden_stage2_differences"].startswith("0 ")
    assert xr["stage2_constants"].startswith("3 team-internal Stage 2 constant(s) (0 derived)")
    capsys.readouterr()
    assert JC.main(str(built["out"] / "stage3_audit.json"), str(built["fx"]["s2"]),
                   str(HERE / "stage3_links_v3.csv"), str(built["fx"]["consts"])) == 0
    line = capsys.readouterr().out.splitlines()[0]
    assert "; 0 cell(s) narrow only as far as public information alone" in line and line.endswith(
        "; 0 cell(s) narrowed further by Stage 3")
    assert "release gate OK" in built["err"]


# ------------------------------------------------------------------- v1.0.2 weight normalizer: growth reference (§10)
GREF = HERE / "maps" / "growth_weight_for_age_ref_v1.csv"
DPM = 365.25 / 12.0


def _gref():
    return B.parse_growth_ref(pd.read_csv(GREF, dtype=str, keep_default_na=False))


def test_growth_ref_pinned_and_parsed():
    assert B.PINS["growth_ref"] == B.sha256_file(GREF) and B.is_full_hash(B.PINS["growth_ref"])
    assert B.MAP_FILES["growth_ref"] == GREF.name and B.MAP_ROWS["growth_ref"] == 3926
    g = _gref()
    assert set(g) == {("cn_li2009", "M"), ("cn_li2009", "F"), ("who", "M"), ("who", "F")}
    assert len(g[("cn_li2009", "M")][0]) == 46 and g[("cn_li2009", "M")][0][-1] == pytest.approx(216 * DPM)
    assert len(g[("who", "F")][0]) == 1857 + 60 and g[("who", "F")][0][-1] == pytest.approx(120 * DPM)


def test_growth_ref_known_table_values():
    g = _gref()
    f = lambda ref, sex, days: float(B.growth_ref_median(g, ref, np.array([days]), np.array([sex]))[0])  # noqa: E731
    # China 2009 (Zong & Li 2013 Table 2, column M = median)
    for sex, mo, m in (("M", 0, 3.32), ("F", 0, 3.21), ("M", 12, 10.05), ("F", 12, 9.40), ("M", 60, 18.98),
                       ("F", 120, 31.76), ("M", 216, 61.40), ("F", 216, 51.41)):
        assert f("cn_li2009", sex, mo * DPM) == pytest.approx(m, abs=1e-9), (sex, mo)
    # WHO 2006 standards (by day) and WHO 2007 reference (by month), column M
    assert f("who", "M", 0) == pytest.approx(3.3464) and f("who", "F", 0) == pytest.approx(3.2322)
    assert f("who", "M", 61 * DPM) == pytest.approx(18.5057) and f("who", "M", 120 * DPM) == pytest.approx(31.1586)
    assert f("who", "F", 120 * DPM) == pytest.approx(31.8578)


def test_growth_ref_interpolation():
    g = _gref()
    f = lambda ref, sex, days: float(B.growth_ref_median(g, ref, np.array([days]), np.array([sex]))[0])  # noqa: E731
    assert f("cn_li2009", "M", 7 * DPM) == pytest.approx((8.41 + 9.05) / 2)            # between 6 and 8 months
    assert f("cn_li2009", "M", 13.5 * DPM) == pytest.approx(10.05 + 0.5 * (10.68 - 10.05))
    assert f("cn_li2009", "F", 10.0) == pytest.approx(3.21 + (10.0 / DPM) * (4.20 - 3.21))
    assert f("cn_li2009", "F", 207 * DPM) == pytest.approx((51.20 + 51.31) / 2)
    # across the WHO 2006 (day 1856) / 2007 (61 months = 1856.67 days) join
    lo, hi = 18.4968, 18.5057
    assert lo <= f("who", "M", 1856.3) <= hi


def test_growth_ref_no_extrapolation_or_invalid_input():
    g = _gref()
    days = np.array([216 * DPM + 1, -1.0, np.nan, 400.0, 120 * DPM + 1])
    sex = np.array(["M", "F", "M", "unrecognized", "F"])
    cn = B.growth_ref_median(g, "cn_li2009", days, sex)
    who = B.growth_ref_median(g, "who", days, sex)
    assert np.isnan(cn[:4]).all() and np.isfinite(cn[4])
    assert np.isnan(who).all()


class _Co:
    def __init__(self, train):
        self._t = np.asarray(train, dtype=bool)

    def train(self):
        return self._t


def _chain_case():
    """30 train encounters with charted weight (10 per sex x band toddler and 10 F school), then five targets."""
    rng = np.random.default_rng(7)
    ab = ["toddler"] * 20 + ["school"] * 10
    sx = ["M"] * 10 + ["F"] * 10 + ["F"] * 10
    days = [800.0] * 20 + [3000.0] * 10
    wt = list(rng.uniform(10, 14, 20)) + list(rng.uniform(25, 35, 10))
    # targets: valid M toddler; 19.2-year-old F; age_invalid F; unrecognized-sex toddler; charted 9 kg
    ab += ["toddler", "ge18", "age_invalid", "toddler", "toddler"]
    sx += ["M", "F", "F", "unrecognized", "M"]
    days += [400.0, 7013.0, np.nan, 400.0, 400.0]
    wt += [np.nan, np.nan, np.nan, np.nan, 9.0]
    train = [True] * 30 + [False] * 5
    return _Co(train), np.array(ab, dtype=object), np.array(sx, dtype=object), np.array(wt), np.array(days)


def test_weight_chain_flags_and_arms():
    co, ab, sx, wt, days = _chain_case()
    g = _gref()
    arms = B.impute_weights(co, ab, sx, wt, days, g)
    s, w = arms.source[30:], arms.weight[30:]
    assert list(s) == ["cn_ref", "train_tier_sex_only", "train_tier_sex_only", "train_tier_age_band", "charted"]
    assert w[0] == pytest.approx(10.05 + (400.0 / DPM - 12) / 3 * (10.68 - 10.05))   # 13.14 months, China boys
    fem = np.median(wt[10:30])                         # train F charted weights (toddler and school)
    assert w[1] == pytest.approx(fem) and w[2] == pytest.approx(fem)
    assert w[3] == pytest.approx(np.median(wt[:20])) and w[4] == 9.0
    # sensitivity arms: train tiers only (the rule before Nov 17) and WHO before the tiers
    assert list(arms.source_train_tier[30:]) == ["train_tier_age_sex", "train_tier_sex_only", "train_tier_sex_only",
                                                 "train_tier_age_band", "charted"]
    assert arms.weight_train_tier[30] == pytest.approx(np.median(wt[:10]))
    assert list(arms.source_who[30:]) == ["who_ref", "train_tier_sex_only", "train_tier_sex_only",
                                          "train_tier_age_band", "charted"]
    assert arms.weight_who[30] == pytest.approx(float(np.interp(400.0, *g[("who", "M")])))
    assert set(arms.source) <= set(B.WEIGHT_SOURCES) and (arms.weight > 0).all()
    assert (arms.weight[:30] == wt[:30]).all() and (arms.source[:30] == "charted").all()


def test_weight_chain_who_fallback_when_china_table_ends():
    """who_ref is reached only where the China table does not cover the age (with the full table: never, since China
    covers 0-216 months and WHO 0-120); shown with the China table cut at 24 months."""
    co, ab, sx, wt, days = _chain_case()
    g = _gref()
    a, m = g[("cn_li2009", "M")]
    g2 = dict(g)
    g2[("cn_li2009", "M")] = (a[a <= 24 * DPM + 1e-9], m[a <= 24 * DPM + 1e-9])
    days = days.copy()
    days[30] = 5 * 365.25                                   # 5-year-old boy, beyond the cut China table
    arms = B.impute_weights(co, ab, sx, wt, days, g2)
    assert arms.source[30] == "who_ref"
    assert arms.weight[30] == pytest.approx(float(np.interp(5 * 365.25, *g[("who", "M")])))


def test_growth_ref_schema_refused():
    d = pd.read_csv(GREF, dtype=str, keep_default_na=False)
    with pytest.raises(B.InputError, match="columns"):
        B.parse_growth_ref(d.drop(columns=["S"]))
    bad = d.copy()
    bad.loc[3, "age_days"] = "0.0000"                       # China girls: second age back to 0 -> not increasing
    with pytest.raises(B.InputError, match="increase"):
        B.parse_growth_ref(bad)
    bad = d.copy()
    bad.loc[0, "sex"] = "U"
    with pytest.raises(B.InputError, match="sex"):
        B.parse_growth_ref(bad)


def test_growth_ref_hash_gate(built, tmp_path):
    fx = built["fx"]
    maps = tmp_path / "maps"
    shutil.copytree(Path(fx["interim"]) / "maps", maps)
    p = maps / B.MAP_FILES["growth_ref"]
    p.write_text(p.read_text(encoding="utf-8").replace("61.4000", "61.5000"), encoding="utf-8")
    out = tmp_path / "out"
    with pytest.raises(B.InputError, match="SHA-256 mismatch"):
        run_build(make_cfg(fx, out, maps_dir=maps))
    assert not out.exists()


def test_weight_normalizer_on_build(built):
    out, A = built["out"], built["A"]
    W = pd.read_csv(out / "stage3_weight_normalizer.csv", keep_default_na=False)
    assert list(W.columns) == ["encounter_id", "weight_kg", "weight_source", "weight_kg_train_tier",
                               "weight_source_train_tier", "weight_kg_who", "weight_source_who"]
    assert list(W["encounter_id"]) == list(A["encounter_id"])
    g = _gref()
    days = pd.to_numeric(A["age_days"], errors="coerce").to_numpy()
    sex = A["sex"].to_numpy()
    band = A["age_band"].to_numpy()
    src = W["weight_source"].to_numpy()
    assert set(src) <= set(B.WEIGHT_SOURCES) and {"charted", "cn_ref"} <= set(src)
    cn = src == "cn_ref"
    assert np.allclose(W["weight_kg"].to_numpy()[cn], B.growth_ref_median(g, "cn_li2009", days, sex)[cn])
    nonch = src != "charted"
    ok_ref = (band != "age_invalid") & np.isin(sex, ["M", "F"]) & (days <= 216 * DPM)
    assert (src[nonch & ok_ref] == "cn_ref").all()
    assert np.isin(src[nonch & ~ok_ref], B.TRAIN_TIERS).all()
    p5 = built["idx"]["P5"]                                  # missing DOB -> age_invalid -> train tier
    assert band[p5] == "age_invalid" and src[p5] in B.TRAIN_TIERS[2:]
    st = W["weight_source_train_tier"].to_numpy()
    assert (st[~nonch] == "charted").all() and np.isin(st[nonch], B.TRAIN_TIERS).all()
    sw = W["weight_source_who"].to_numpy()
    assert (sw[~nonch] == "charted").all() and set(sw) <= {"charted", "who_ref", *B.TRAIN_TIERS}
    assert (sw[nonch & ok_ref & (days <= 120 * DPM)] == "who_ref").all()
    for c in ("weight_kg", "weight_kg_train_tier", "weight_kg_who"):
        assert (W[c].to_numpy(dtype=float) > 0).all()
        assert np.allclose(W[c].to_numpy(dtype=float)[~nonch], W["weight_kg"].to_numpy(dtype=float)[~nonch])
    assert built["res"]["checks"]["weight_reference_rule"] == "pass"
    meta = json.loads((out / "stage3_weight_normalizer.csv.meta.json").read_text(encoding="utf-8"))
    assert meta["growth_reference"]["sha256"] == B.PINS["growth_ref"] and meta["feature_columns"] == []
    doc = json.loads((out / "stage3_audit.json").read_text(encoding="utf-8"))
    assert "China 2009" in doc["notes"]["weight_rule"]


# ----------------------------------------------------- Stage 3 review decisions 2-9 (Vinai, Nov 17, 2026; spec §15)
def _doc(built):
    return json.loads((built["out"] / "stage3_audit.json").read_text(encoding="utf-8"))


def _meta(built, name):
    return json.loads((built["out"] / (name + ".meta.json")).read_text(encoding="utf-8"))


def test_classify_rows_zero_missing_and_gate_side():
    g = pd.DataFrame({"vnum": [0.0, 2.0, 150.0, 30.0, np.nan], "vtext": [""] * 5, "unit": [""] * 5})
    info = {"channel": "respiratory_rate", "low": 5.0, "high": 120.0, "conv": ""}
    _v, st, side = B.classify_rows(info, g, False)
    assert list(st) == ["zero_missing", "out_of_gate", "out_of_gate", "measurement", "nonnumeric"]
    assert list(side) == [0, -1, 1, 0, 0]
    _v, st2, side2 = B.classify_rows(dict(info, channel="wbc", low=0.1, high=200.0), g, False)
    assert st2[0] == "out_of_gate" and side2[0] == -1               # a WBC of 0 stays out of gate (below)


def test_respiratory_rate_zero_missing_and_sides_on_build(built):
    z, i, Bb, s3 = built["z"], built["idx"]["P0"], built["B"], built["fx"]["s3"]
    k = ch(z, "respiratory_rate")
    assert z["median"][i, 0, k] == pytest.approx(30.0) and z["min"][i, 0, k] == pytest.approx(30.0)  # 0, 2, 150 dropped
    assert float(Bb.loc[i, "respiratory_rate_min"]) == pytest.approx(30.0)
    b = "/plausibility/respiratory_rate"
    for sc in ("train", "overall"):
        assert s3[f"{b}/{sc}/n_rows_zero_as_missing"] >= 1
        assert s3[f"{b}/{sc}/n_rows_numeric"] == sum(s3[f"{b}/{sc}/n_rows_{x}"] for x in ("in_gate", "out_of_gate",
                                                                                         "zero_as_missing"))
    assert s3[f"{b}/overall/n_rows_below_gate"] >= 1
    assert s3[f"{b}/overall/n_rows_below_gate"] + s3[f"{b}/overall/n_rows_above_gate"] == s3[f"{b}/overall/n_rows_out_of_gate"]
    assert f"{b}/train/n_rows_below_gate" not in s3                    # split at the overall scope only
    doc = _doc(built)
    assert set(doc["plausibility"]["respiratory_rate"]) == {"overall"}  # train stays withheld (drop list)
    assert {"n_rows_zero_as_missing", "n_rows_below_gate", "n_rows_above_gate"} <= set(doc["plausibility"]["respiratory_rate"]["overall"])
    ids = {(t, tuple(p)) for t, p in doc["disclosure"]["identities"]}
    o = f"{b}/overall/n_rows_"
    assert (o + "out_of_gate", (o + "below_gate", o + "above_gate")) in ids
    assert (o + "numeric", (o + "in_gate", o + "below_gate", o + "above_gate", o + "zero_as_missing")) in ids
    links = (HERE / "stage3_links_v2.csv").read_text(encoding="utf-8")
    assert ("eq,s2:/map_checks/row_status/respiratory_rate/n_out_of_gate,/plausibility/respiratory_rate/train/"
            "n_rows_out_of_gate;/plausibility/respiratory_rate/train/n_rows_zero_as_missing,") in links


def test_wbc_gate_sides_only_from_20_rows(built):
    z, i, s3 = built["z"], built["idx"]["P0"], built["fx"]["s3"]
    assert z["median"][i, 0, ch(z, "wbc")] == pytest.approx(12.0)      # 0.05, 0 (below) and 250 (above) are missing
    b = "/plausibility/wbc/overall"
    out = s3[f"{b}/n_rows_out_of_gate"]
    assert (f"{b}/n_rows_below_gate" in s3) == (out >= B.GATE_SIDE_MIN_TOTAL)
    assert "n_rows_zero_as_missing" not in _doc(built)["plausibility"]["wbc"]["overall"]
    assert not any(k.startswith("/plausibility/wbc/train/n_rows_below") for k in s3)


def test_base_excess_bicarbonate_one_value_per_draw(built):
    z, i = built["z"], built["idx"]["P0"]
    assert z["median"][i, 1, ch(z, "base_excess")] == pytest.approx(-2.0)   # 5211 kept over 5249 (same CHARTTIME)
    assert z["median"][i, 1, ch(z, "bicarbonate")] == pytest.approx(20.0)   # 5248 kept over 5224
    cm = pd.read_csv(HERE / "maps" / "feature_channel_map_v2_3.csv", dtype=str, keep_default_na=False,
                     encoding="utf-8-sig").apply(lambda s: s.map(lambda v: str(v).strip()))
    e = built["fx"]["expect"]
    ex = B.Expect(encounters=e["encounters"], subjects=e["subjects"], died=e["died"], splits=e["splits"],
                  eligible=e["eligible"], layer_c_surgeries=None)
    chans = B.parse_channel_map(cm, ex)["channels"]
    for c, m in (("ph", "5237"), ("pco2", "5235"), ("base_excess", "5211"), ("bicarbonate", "5248")):
        assert chans[c].dedupe_same_time and chans[c].main_item == m, c
    assert built["res"]["checks"]["blood_gas_one_value_per_draw"] == "pass"


def test_landmark_train_median_fallback_rules():
    n = 4
    out_h = np.array([30.0, 60.0, 400.0, 400.0])
    train = np.array([True, True, True, False])
    med = np.full((n, B.N_BINS), np.nan)
    med[2, 5] = 7.0                     # first train value at 30-36 h: the L24 prefix (bins 0-3) is empty
    med[3, 0] = 99.0                    # validation value: never used
    val, src = B.landmark_train_medians(med, train, out_h)
    assert src[24] == "pooled_all_landmarks" and val[24] == 7.0
    assert src[48] == "landmark" and val[48] == 7.0 and src[336] == "landmark"
    med2 = np.full((n, B.N_BINS), np.nan)
    med2[0, 0] = 5.0                    # eligible at L24 only
    med2[3, 1] = 50.0
    val2, src2 = B.landmark_train_medians(med2, train, out_h)
    assert src2[24] == "landmark" and val2[24] == 5.0
    assert src2[48] == "earlier_landmark_L024" and val2[48] == 5.0 and src2[336] == "earlier_landmark_L024"
    med3 = np.full((n, B.N_BINS), np.nan)
    med3[3, 0] = 1.0                    # validation only
    val3, src3 = B.landmark_train_medians(med3, train, out_h)
    assert set(src3.values()) == {"none"} and all(np.isnan(v) for v in val3.values())


def test_train_median_fallback_recorded_on_build(built):
    out, doc = built["out"], _doc(built)
    fb = doc["notes"]["train_median_fallback"]
    allowed = {"landmark", "pooled_all_landmarks", "none"} | {f"earlier_landmark_L{x:03d}" for x in B.LANDMARKS}
    for L in B.LANDMARKS:
        f = np.load(out / f"stage3_layer_c_filled_L{L}.npz")
        src = [str(x) for x in f["train_median_source"]]
        assert len(src) == len(f["value_channels"]) and set(src) <= allowed
        for c, s in zip(f["value_channels"], src):
            if s != "landmark":
                assert fb[str(c)][f"L{L:03d}"] == s
            if s != "none":
                assert np.isfinite(f["train_median"][list(f["value_channels"]).index(c)])
    assert fb["kdigo_uo_stage"]["L024"] == "pooled_all_landmarks"    # KDIGO needs 24 h of history: L24 is empty
    assert all(isinstance(v, str) for d in fb.values() for v in d.values())       # names and types only, no counts


def test_unparsed_text_flag_and_team_internal_list(built):
    out, z, i, Bb = built["out"], built["z"], built["idx"]["P0"], built["B"]
    qc = list(z["qual_channels"])
    assert len(qc) == 14 and z["qual_unparsed"].shape == (len(built["eids"]), B.N_BINS, 14)
    k = qc.index("stool_occult_blood")
    assert z["qual_unparsed"][i, 5, k] == 1                            # P0 strings at 30-35 h (+1 h release)
    assert z["measured"][i, 5, ch(z, "stool_occult_blood")] == 0 and z["missing"][i, 5, ch(z, "stool_occult_blood")] == 1
    col = "stool_occult_blood" + B.UNPARSED_FLAG_SUFFIX
    assert col in Bb.columns and set(Bb[col]) <= {"0", "1"} and "1" in set(Bb[col])
    mb = _meta(built, "stage3_layer_b_base.csv")
    assert col in mb["flag_columns"] and col not in mb["feature_columns"]
    md = (out / B.UNPARSED_FILE).read_text(encoding="utf-8")
    assert md.startswith("# Stage 3 unparsed qualitative text (TEAM-INTERNAL)") and "do not share outside" in md
    sec = md.split("## stool_occult_blood")[1].split("##")[0]
    assert "- `1:100`" in sec and "- `weak?`" in sec and "- `seeattachedreportforfulldetailsp…`" in sec
    for bad in ("床号12", "12345", "rare", "seeattachedreportforfulldetailsplease"):
        assert bad not in md, bad
    body = [ln for ln in md.splitlines() if ln and not ln.startswith(("#", "**", "Distinct", "(no string", "(not listed"))]
    assert all(ln.startswith("- `") for ln in body)                    # strings only: no counts, no hidden tally
    man = json.loads((out / "stage3_manifest.json").read_text(encoding="utf-8"))
    assert B.UNPARSED_FILE in man["files"]


def test_unparsed_text_listing_rules_unit():
    from types import SimpleNamespace as NS
    rows = ([("urine_nitrite", 0, "abc")] * 10 + [("urine_nitrite", 1, "zzz")] * 10 + [("urine_nitrite", 0, "x9y")] * 10
            + [("urine_nitrite", 0, "<0.5")] * 10 + [("urine_nitrite", 0, "123456")] * 10
            + [("urine_nitrite", 0, "few")] * 9 + [("stool_occult_blood", 0, "abc")] * 12)
    up = pd.DataFrame({"enc_i": [r[1] for r in rows], "t": 5.0, "channel": [r[0] for r in rows], "text": [r[2] for r in rows]})
    co = NS(W=np.array([100.0, 100.0]), split=np.array(["train", "validation"]))
    cmap = {"feature_channels": ["stool_occult_blood", "urine_nitrite"],
            "channels": {"stool_occult_blood": NS(qualitative=True), "urine_nitrite": NS(qualitative=True)}}
    doc = {"plausibility": {"stool_occult_blood": {"train": {"n_rows_unparsed": "<10"}},
                            "urine_nitrite": {"train": {"n_rows_unparsed": 59}}}}
    md = B.unparsed_text_report(NS(unparsed=up), co, cmap, doc)
    un = md.split("## urine_nitrite")[1]
    assert "- `abc`" in un and "- `<0.5`" in un
    for bad in ("zzz", "x9y", "123456", "few"):                          # validation, digits, long number, < 10
        assert bad not in md
    st = md.split("## stool_occult_blood")[1].split("##")[0]
    assert "abc" not in st and "not published" in st                   # suppressed count: channel not listed


def test_fluid_intake_kdigo_provisional_metadata(built):
    mb = _meta(built, "stage3_layer_b_base.csv")
    assert {"fluid_intake_median", "kdigo_uo_stage_max"} <= set(mb["provisional_columns"])
    assert all(c.startswith(("fluid_intake_", "kdigo_uo_stage_")) for c in mb["provisional_columns"])
    assert set(mb["provisional_columns"]) <= set(mb["feature_columns"]) and "B38" in mb["provisional_reason"]
    for nm in ("stage3_layer_c_raw.npz", "stage3_layer_c_filled_L24.npz"):
        assert _meta(built, nm)["provisional_channels"] == ["fluid_intake", "kdigo_uo_stage"], nm
    assert _doc(built)["notes"]["provisional_features"]["channels"] == ["fluid_intake", "kdigo_uo_stage"]


def test_drugs_at_t0_lookback_and_future_enddate(built):
    A, i = built["A"], built["idx"]["P1"]
    assert A.loc[i, "med_sedative_analgesic_active_t0"] == "1"         # t0 - 5 h; its ENDDATE t0 + 30 h is not used
    assert A.loc[i, "med_diuretic_active_t0"] == "0"                   # ENDDATE t0 - 1 h is known at t0: stopped
    assert A.loc[i, "med_antibiotic_active_t0"] == "0"                 # t0 - 10 h: outside the 6 h look-back
    assert A.loc[i, "med_antibiotic_active_t0_lb24"] == "1"            # inside the 24 h sensitivity
    assert built["res"]["checks"]["med_future_enddate_differences"] == "pass"


def test_icdo_morphology_to_oncology_missing_separate(built):
    A, eids, s3 = built["A"], built["eids"], built["fx"]["s3"]
    for i in (built["idx"]["P3"], eids.index("E100007_1")):            # M98000, M9800/3
        assert A.loc[i, "dx_icd10_3char"] == "icdo_morphology" and A.loc[i, "dx_domain"] == "oncology_hematology"
        assert A.loc[i, "dx_icdo_morphology"] == "1" and A.loc[i, "dx_mapped"] == "0"
    i7 = eids.index("E100008_1")                                       # M80000 is in D_ICD_DIAGNOSES: ICD-10 wins
    assert A.loc[i7, "dx_icd10_3char"] == "M80" and A.loc[i7, "dx_mapped"] == "1"
    miss = A["dx_code_status"] == "missing"
    assert miss.any() and (A.loc[miss, "dx_icd10_3char"] == "missing").all() and (A.loc[miss, "dx_domain"] == "missing").all()
    assert s3["/diagnosis/n_icdo_morphology_oncology"] == 2 + len(FX.DX_ICDO_PLANTED)
    assert not any(k.startswith("/diagnosis/unmapped_codes/M") for k in s3)
    for c, ok in (("M98000", True), ("M9800/3", True), ("m 8000 3", True), ("M80.000", False), ("M54", False),
                  ("M7000", False), ("C91.0", False), ("", False)):
        assert B.is_icdo_morphology(c) is ok, c


# ------------------------------------------------------------------- v1.0.3: ICD-O real formats, ketones excluded
REAL_UNMAPPED_V102 = ("M98000", "M95000", "M80001", "M90800", "M89700", "M95900", "M91200", "M89600", "M91700",
                      "M89000", "M98010")       # exactly as listed in diagnosis/unmapped_codes of the real v1.0.2 audit


@pytest.mark.parametrize("code", REAL_UNMAPPED_V102 + tuple(f"{c}/3" for c in REAL_UNMAPPED_V102) + (
    "M9800/3", "M9800/0", "M9800", "m98000/3", "  m 98000 / 3 ", "M-9800/3", "M9800.0", "M980003", "M98000/3,C92.0",
    "M98000 C92.0", "M98000;C92.0", "M98000、C92.0", "\uff2d\uff19\uff18\uff10\uff10\uff10\uff0f\uff13", "M91200/0",
    "M80001/1", "M80002/2", "M98000/6", "M98000/9"))
def test_icdo_real_formats_recognized(code):
    assert B.is_icdo_morphology(code), code
    p = B.icdo_parts(code)
    assert p["histology"][0] in "89" and len(p["histology"]) == 4


@pytest.mark.parametrize("code", ("M80.000", "M80.0", "M54", "M54.500", "M7000", "M79990/3", "M98000/5", "M98000/3x",
                                  "M9800//3", "C92.0,M98000/3", "C91.000", "XYZ", "", "M", "M98", "8000/3"))
def test_icdo_non_morphology_not_recognized(code):
    assert not B.is_icdo_morphology(code), code


def test_icdo_behaviour_parts():
    assert B.icdo_parts("M98000/3") == {"histology": "9800", "extension": "0", "behaviour": "3"}
    assert B.icdo_parts("M9800/3") == {"histology": "9800", "extension": "", "behaviour": "3"}
    assert B.icdo_parts("M98000") == {"histology": "9800", "extension": "0", "behaviour": ""}
    assert B.icdo_parts("M980003") == {"histology": "9800", "extension": "0", "behaviour": "3"}
    assert B.icdo_parts("M91200/0")["behaviour"] == "0"          # benign: still oncology (decision 8, all behaviours)


def test_v102_rule_missed_suffix_form():
    """Root cause of the v1.0.2 miss: the old rule matched the whole field, so 'M98000/3' (clinical-edition code with
    its behaviour suffix) failed while norm_icd's first code 'M98000' was listed as unmapped."""
    import re
    old = re.compile(r"M[89][0-9]{3}(?:[/.]?[0-9])?")
    assert not old.fullmatch("M98000/3") and B.norm_icd("M98000/3")[0] == "M98000"
    assert B.is_icdo_morphology("M98000/3")


def test_icdo_real_format_fixture_rows_go_to_oncology(built):
    A, eids, s3 = built["A"], built["eids"], built["fx"]["s3"]
    for si, code in FX.DX_ICDO_PLANTED.items():
        i = eids.index(f"E{100001 + si}_1")
        assert A.loc[i, "dx_icd10_3char"] == "icdo_morphology", code
        assert A.loc[i, "dx_domain"] == "oncology_hematology" and A.loc[i, "dx_icdo_morphology"] == "1", code
        assert A.loc[i, "dx_mapped"] == "0" and A.loc[i, "dx_missing"] == "0", code
    for c in REAL_UNMAPPED_V102:                                  # none of the real codes is listed as unmapped
        assert f"/diagnosis/unmapped_codes/{c}/n_encounters" not in s3
    assert not any(k.startswith("/diagnosis/unmapped_codes/M") for k in s3)
    # 'missing' stays its own level, separate from unmapped and from ICD-O
    tot = sum(s3[f"/diagnosis/{k}"] for k in ("n_mapped", "n_missing", "n_icdo_morphology_oncology",
                                              "n_unmapped_with_chapter", "n_unmapped_unparsed"))
    assert tot == len(A) and s3["/diagnosis/n_missing"] == int((A["dx_code_status"] == "missing").sum())
    assert built["res"]["checks"]["icdo_morphology_oncology_missing_separate"] == "pass"


def test_ketones_urine_excluded_everywhere(built):
    out, z, A, Bb = built["out"], built["z"], built["A"], built["B"]
    assert "ketones_urine" in B.EXCLUDED_CHANNELS and B.Expect().feature_channels == 92
    assert "ketones_urine" not in list(z["channels"])
    assert not [c for c in list(A.columns) + list(Bb.columns) if "ketone" in c.lower()]
    for L in B.LANDMARKS:
        f = np.load(out / f"stage3_layer_c_filled_L{L}.npz")
        assert "ketones_urine" not in [str(x) for x in f["value_channels"]]
    doc = _doc(built)
    assert doc["maps"]["feature_channels"] == 92 and doc["maps"]["map_channels"] == 94
    assert set(doc["maps"]["excluded_channels"]) == {"ketones_urine"}
    assert "dropped as a feature" in doc["maps"]["excluded_channels"]["ketones_urine"]
    assert not any(ch.isdigit() for ch in doc["maps"]["excluded_channels"]["ketones_urine"].replace("5356", "")
                   .replace("v2.3", "").replace("v1.0.2", "").replace("Nov 17, 2026", ""))   # no counts
    assert doc["notes"]["excluded_channels"] == doc["maps"]["excluded_channels"]
    rest = dict(doc, notes={k: v for k, v in doc["notes"].items() if k != "excluded_channels"})
    assert "ketones_urine" not in json.dumps({k: v for k, v in rest.items() if k != "maps"})
    fb = doc["notes"]["train_median_fallback"]
    assert not isinstance(fb, dict) or "ketones_urine" not in fb
    md = (out / "stage3_audit.md").read_text(encoding="utf-8")
    assert [ln for ln in md.splitlines() if "ketones_urine" in ln] == [
        f"- `notes/excluded_channels/ketones_urine`: {B.EXCLUDED_CHANNELS['ketones_urine']}"]   # only as excluded
    assert "maps parsed: 94 channels (92 feature channels" in built["err"] and "excluded: ketones_urine" in built["err"]


def test_excluded_channel_must_exist_in_map():
    saved = dict(B.EXCLUDED_CHANNELS)
    try:
        B.EXCLUDED_CHANNELS["not_a_channel"] = "x"
        cm = pd.read_csv(HERE / "maps" / "feature_channel_map_v2_3.csv", dtype=str, keep_default_na=False,
                         encoding="utf-8-sig").apply(lambda s: s.map(lambda v: str(v).strip()))
        with pytest.raises(B.InputError, match="excluded channel"):
            B.parse_channel_map(cm, B.Expect(feature_channels=91))
    finally:
        B.EXCLUDED_CHANNELS.clear()
        B.EXCLUDED_CHANNELS.update(saved)


def test_train_median_check_passes_with_pooled_kdigo(built):
    """kdigo_uo_stage at L24 falls back to the pooled train median (finite): that is a fallback, not 'unavailable';
    the check is 'noted' only when a landmark has no median at all (ketones_urine in the real v1.0.2 build)."""
    doc = _doc(built)
    assert doc["notes"]["train_median_fallback"]["kdigo_uo_stage"]["L024"] == "pooled_all_landmarks"
    assert built["res"]["checks"]["train_median_available_all_channels"] == "pass"
    assert sum(v == "pass" for v in built["res"]["checks"].values()) == 24


def test_no_edregtime_feature_and_hadm_only_linkage(built):
    A, Bb, i = built["A"], built["B"], built["idx"]["P5"]
    bad = [c for c in list(A.columns) + list(Bb.columns) if "ed_component" in c or "edreg" in c.lower() or "visit" in c.lower()]
    assert not bad
    assert A.loc[i, "remote_surgery"] == "1"          # VISIT_ID differs, HADM_ID matches: linked (by HADM_ID)
    assert A.loc[i, "recent_surgery"] == "0"          # VISIT_ID equals this HADM, HADM_ID does not: never linked
    S = pd.read_csv(built["out"] / "stage3_surgery_long_INTERNAL.csv")
    assert len(S[S["encounter_id"] == built["eids"][i]]) == 1
    doc = _doc(built)
    assert "n_present" in doc["edregtime"] and "SURGERY_INFO" in doc["visit_id"]   # aggregates only
    assert built["res"]["checks"]["no_edregtime_or_visit_id_feature"] == "pass"


# ------------------------------------------------------------------- v1.0.4: per-channel tokens, below-cutoff flag
def test_qual_tokens_per_channel_unit():
    t = np.array(["positive+", "Positive ++", "POSITIVE+", "<3.00", " < 3.00 ", "<3.0", "阳性", "neg"], dtype=object)
    for chn in ("stool_occult_blood", "urine_leukocyte_esterase", "urine_nitrite", "csf_protein_pandy"):
        q = B.qual_parse(t, chn)
        assert list(q[:3]) == [1.0, 1.0, 1.0] and np.isnan(q[3]) and np.isnan(q[5]), chn
    q = B.qual_parse(t, "resp_panel_mycoplasma")
    assert np.isnan(q[0]) and np.isnan(q[1]) and q[3] == 0.0 and q[4] == 0.0 and np.isnan(q[5])   # '<3.0' not listed
    q = B.qual_parse(t, "resp_panel_rsv")                     # a qualitative channel without additions
    assert np.isnan(q[:6]).all() and q[6] == 1.0 and q[7] == 0.0
    q = B.qual_parse(t)                                       # no channel: the global table only
    assert np.isnan(q[:6]).all()
    assert list(B.qual_below_cutoff(t, "resp_panel_mycoplasma")) == [False, False, False, True, True, False, False, False]
    assert not B.qual_below_cutoff(t, "stool_occult_blood").any()
    assert set(B.QUAL_POS_BY_CHANNEL) == {"stool_occult_blood", "urine_leukocyte_esterase", "urine_nitrite",
                                          "csf_protein_pandy"}
    assert B.QUAL_NEG_BY_CHANNEL == {"resp_panel_mycoplasma": {"<3.00"}} == B.QUAL_BELOW_CUTOFF


def test_qual_tokens_in_build_and_below_cutoff_flag(built):
    z, i, Bb, s3 = built["z"], built["idx"]["P2"], built["B"], built["fx"]["s3"]
    rows = {(it, h): t for it, h, t in FX.QUAL_V104_PLANTED}
    assert rows[(5326, 2.0)] == "positive+"
    assert z["median"][i, 0, ch(z, "stool_occult_blood")] == 1.0                 # positive+ at 2 h (+1 h release)
    assert z["median"][i, 0, ch(z, "urine_nitrite")] == 1.0                      # 'Positive ++' at 3 h
    assert z["median"][i, 0, ch(z, "csf_protein_pandy")] == 1.0
    assert z["median"][i, 1, ch(z, "resp_panel_mycoplasma")] == 0.0              # '<3.00' at 5 h -> bin 1 (6 h)
    assert z["median"][i, 6, ch(z, "urine_leukocyte_esterase")] == 1.0          # 40 h (+1 h) -> bin 6
    k = list(z["qual_channels"]).index("resp_panel_mycoplasma")
    assert z["qual_unparsed"][i, 5, k] == 1                                      # 'positive+' stays unparsed here
    assert z["qual_unparsed"][i, 1, k] == 0
    col = "resp_panel_mycoplasma" + B.BELOW_CUTOFF_FLAG_SUFFIX
    assert col in Bb.columns and Bb.loc[i, col] == "1" and set(Bb[col]) <= {"0", "1"}
    assert Bb.loc[i, "resp_panel_mycoplasma" + B.UNPARSED_FLAG_SUFFIX] == "0"      # 30 h is outside the B window
    mb = _meta(built, "stage3_layer_b_base.csv")
    assert col in mb["flag_columns"] and col not in mb["feature_columns"] and "below_cutoff" in mb["flag_note"]
    flags = [c for c in Bb.columns if c.endswith((B.UNPARSED_FLAG_SUFFIX, B.BELOW_CUTOFF_FLAG_SUFFIX))]
    assert len(flags) == 15 and sorted(mb["flag_columns"]) == sorted(flags)        # 14 unparsed + 1 below-cutoff
    assert len(mb["feature_columns"]) == Bb.shape[1] - 1 - 15
    doc = _doc(built)
    assert "positive+" in doc["notes"]["qualitative_tokens"] and "<3.00" in doc["notes"]["qualitative_tokens"]
    md = (built["out"] / B.UNPARSED_FILE).read_text(encoding="utf-8")
    assert "positive+" not in md.split("Distinct")[0] and "v1.0.4 extends the token table" in md
    for sec in md.split("## ")[1:]:
        assert "- `positive+`" not in sec and "- `<3.00`" not in sec


def test_layer_b_column_count_v104(built):
    """Layer B gains one flag column (resp_panel_mycoplasma_below_cutoff); Layer A is unchanged."""
    a_meta, b_meta = _meta(built, "stage3_layer_a_base.csv"), _meta(built, "stage3_layer_b_base.csv")
    assert built["B"].shape[1] == len(b_meta["feature_columns"]) + len(b_meta["flag_columns"]) + 1
    assert "resp_panel_mycoplasma_below_cutoff" in list(built["B"].columns)


# ------------------------------------------------------------------- v1.0.4: team-internal Stage 2 constants
def _xr(built, consts_rows, tmp_path):
    p = tmp_path / "c.csv"
    with open(p, "w", encoding="utf-8") as f:
        f.write("pointer,value,reason\n" + "".join(f"{a},{b},{c}\n" for a, b, c in consts_rows))
    doc = json.loads(Path(built["fx"]["s2"]).read_text(encoding="utf-8"))
    links = B.read_links_v2(HERE / "stage3_links_v3.csv")
    return p, doc, links


def test_constants_read_and_derive(built, tmp_path):
    rel = json.loads(Path(built["fx"]["s2"]).read_text(encoding="utf-8"))
    inv_label = B.json_ptr_get(rel, "/layer_a/dob_breakdown/n_invalid")
    s3 = built["fx"]["s3"]
    true_inv = s3["/age/n_invalid"]
    rows = [(f"s2:{p}", "0", "zero") for p in FX.S2_GE18] + [
        ("s2:/layer_a/dob_breakdown/n_invalid", str(true_inv), "stop event"),
        ("s2:/layer_a/dob_breakdown/status/n_valid", "derived", "cohort minus invalid"),
        ("s2:/layer_a/dob_breakdown/n_age_days", "derived", "valid")]
    p, doc, links = _xr(built, rows, tmp_path)
    c = B.read_constants(p, doc)
    assert [v for _p, v, _r in c][-2:] == ["derived", "derived"] and len(c) == 6
    xr = B.CrossRelease(doc=doc, sha256="x", links=links, drops=[], constants=c)
    out = B.derive_constants(xr)
    assert out["s2:/layer_a/dob_breakdown/status/n_valid"] == s3["/age/n_valid"]
    assert out["s2:/layer_a/dob_breakdown/n_age_days"] == s3["/age/n_valid"]
    assert out["s2:/layer_a/train/age_band/n_ge18"] == 0
    assert inv_label in ("<10", "S") or isinstance(inv_label, int)


@pytest.mark.parametrize("rows,match", [
    ([("/layer_a/train/age_band/n_ge18", "0", "no prefix")], "pointer starting with 's2:/'"),
    ([("s2:/layer_a/train/age_band/n_ge18", "0", "")], "and a reason"),
    ([("s2:/layer_a/train/age_band/n_ge18", "0", "a"), ("s2:/layer_a/train/age_band/n_ge18", "0", "b")], "duplicate"),
    ([("s2:/no/such/cell", "0", "x")], "not a count cell"),
    ([("s2:/layer_a/train/age_band/n_ge18", "12", "x")], "contradicts the Stage 2 '<10' label"),
    ([("s2:/layer_a/train/age_band/n_ge18", "-1", "x")], "non-negative integer or 'derived'"),
    ([("s2:/denominators/n_encounters", "1", "x")], "contradicts the published Stage 2 value"),
    ([("s2:/layer_a/dob_breakdown/status/n_missing", "derived", "x")], "not fixed by the Stage 2 identities"),
])
def test_constants_fail_closed(built, tmp_path, rows, match):
    p, doc, links = _xr(built, rows, tmp_path)
    with pytest.raises(B.InputError, match=match):
        c = B.read_constants(p, doc)
        B.derive_constants(B.CrossRelease(doc=doc, sha256="x", links=links, drops=[], constants=c))


def test_constants_messages_never_print_values(built, tmp_path):
    p, doc, links = _xr(built, [("s2:/denominators/n_encounters", "987654", "x")], tmp_path)
    with pytest.raises(B.InputError) as e:
        B.read_constants(p, doc)
    assert "987654" not in str(e.value) and str(B.json_ptr_get(doc, "/denominators/n_encounters")) not in str(e.value)


def test_constants_file_hash_gate_and_not_recorded(built, tmp_path):
    fx = built["fx"]
    _fails_input(fx, tmp_path, "constants file SHA-256 mismatch", stage2_constants_sha256="0" * 64)
    _fails_input(fx, tmp_path, "constants file not found", stage2_constants_path=tmp_path / "nope.csv")
    sha = fx["consts_sha"]
    for p in built["out"].iterdir():
        if p.is_file() and p.suffix in (".json", ".md"):
            assert sha not in p.read_text(encoding="utf-8"), p.name
    man = json.loads((built["out"] / "stage3_manifest.json").read_text(encoding="utf-8"))
    assert "stage2_constants" not in man["input_hashes"]
    assert "values and file hash not recorded" in man["cross_release"]["stage2_constants"]


@pytest.fixture(scope="session")
def built_stop(tmp_path_factory):
    """The real v1.0.4 constants shape on the fixture: n_ge18 zeros + the stop-event invalid-DOB total + 2 derived
    cells. The fixture's synthetic Stage 2 release could not be gated with the invalid total public (it is < 10 here,
    as in the Stage 2 mirror), so Stage 3 knows more than Stage 2 assumed: cells it exposes are flagged, the gate
    still passes and nothing is narrowed beyond public information (incl. the constants)."""
    root = tmp_path_factory.mktemp("stop")
    fx = add_stage2_release(FX.make_project(root / "proj"), root / "s2", stop_event=True)
    out = root / "out"
    res, err = run_build(make_cfg(fx, out))
    return dict(fx=fx, out=out, err=err)


def test_constants_stop_event_build(built_stop):
    out, fx = built_stop["out"], built_stop["fx"]
    man = json.loads((out / "stage3_manifest.json").read_text(encoding="utf-8"))
    xr = man["cross_release"]
    assert man["release_gate"] == "ok"
    assert xr["stage2_constants"].startswith("6 team-internal Stage 2 constant(s) (2 derived)")
    assert "/layer_a/dob_breakdown/n_invalid" not in xr["flagged_for_stage2_v1_0_4"]   # a constant is known, not flagged
    assert "6 team-internal Stage 2 constants (2 derived)" in built_stop["err"]
    assert JC.main(str(out / "stage3_audit.json"), str(fx["s2"]), str(HERE / "stage3_links_v3.csv"),
                   str(fx["consts"])) == 0


# ------------------------------------------------------------------- v1.0.4: module v1.5 integration
def test_joint_system_v15_api(built, tmp_path):
    SM = _sm()
    assert SM.SCRIPT_VERSION == "stage2_v1_0_3_suppress v1.5" and SM.PROVENANCE_KEY == "release_v1_0_3"
    rels = [("s2:t", ("a", "s2:b"), "x")]
    extra = {"s2:t": "const", "s2:b": "P"}
    sysj = B.joint_system(SM, rels, {"s2:t": 12, "a": 5}, {"a": "P", **extra}, SM.Params(), extra, {})
    st = {"a": "P"}
    assert sysj.local_rels(["a"], st) == {0}
    tb = sysj.target_bounds(st, {0}, ["a"], lp_cap=5.0)
    assert tb["a"] is not None and tb["a"][0] >= 3 - 1e-6 and tb["a"][1] <= 9 + 1e-6


def test_relabel_release_v104():
    SM = _sm()
    doc = {"x": 1, "release_v1_0_3": {"a": 1},
           "notes": {"dropped_subtrees_v1_0_3": "2 subtree(s) withheld from the v1.0.3 release as a whole", "k": "v"}}
    out = B.relabel_release(doc, SM)
    assert list(out) == ["x", "release_v1_0_4", "notes"]
    assert out["notes"] == {"dropped_subtrees_v1_0_4": "2 subtree(s) withheld from the Stage 3 v1.0.4 release as a whole",
                            "k": "v"}


def test_heartbeat_lines_during_build(built):
    """Every progress line goes through the module's print lock and resets its silence clock; the heartbeat thread
    is started for the build (a heartbeat line appears only after 30 s without a line)."""
    assert "03_build_real_layers v1.0.4: precheck" in built["err"]
    assert B._SM_REF["m"] is None                     # cleared after the build


def test_heartbeat_prints_when_silent(built, tmp_path, monkeypatch):
    fx = built["fx"]
    real = B.load_cohort

    def slow(cfg, A):
        time_ = __import__("time")
        time_.sleep(2.5)
        return real(cfg, A)
    monkeypatch.setattr(B, "load_cohort", slow)
    cfg = make_cfg(fx, tmp_path / "out", heartbeat_seconds=0.6, precheck_only=True)
    _res, err = run_build(cfg)
    assert "heartbeat: still working - Stage 3 build: cross-release inputs" in err
    assert "LP solves so far" in err and "time left" in err


def test_convergence_caps_reach_the_module(built, tmp_path, capsys, monkeypatch):
    cfg = make_cfg(built["fx"], tmp_path / "out", max_new_suppressions=1)
    monkeypatch.setattr(B, "build_config", lambda argv=None: cfg)
    rc = B.main([])
    cap = capsys.readouterr()
    assert rc == 3 and "exceed the limit of 1" in cap.err and "nothing left behind" in cap.out
    assert not (tmp_path / "out").exists() and not list(tmp_path.glob("out.pending-*"))
