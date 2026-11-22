"""Tests for stage2_followup_checks.py. SYNTHETIC DATA ONLY - never real data.

Run:  python -m pytest -q test_stage2_followup_checks.py
Builds the synthetic Stage 2 project (../feature_sparsity_frequency_audit/synthetic_stage2_data.py), appends planted
bilirubin 5075/5255 and FiO2 rows, copies the decided channel map v2 (a decision file, no patient rows) and overrides
the pins for the synthetic files only."""
from __future__ import annotations

import dataclasses
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "feature_sparsity_frequency_audit"))
import stage2_followup_checks as F  # noqa: E402
import synthetic_stage2_data as S  # noqa: E402

M = F.load_stage2(HERE / F.STAGE2_FILE)
MAP_V2 = Path(__file__).resolve().parents[1] / "stage3_builder" / "maps" / "feature_channel_map_v2_3.csv"

def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plant(proj: dict, seed: int = 11) -> dict:
    rng = np.random.default_rng(seed)
    raw = proj["raw"]
    lb = pd.read_csv(raw / "LABEVENTS.csv", dtype=str, keep_default_na=False)
    enc = proj["enc"]
    add = []
    for r in enc.itertuples():
        sp = (r.OUT - r.t0).total_seconds() / 3600
        for k in rng.uniform(0, min(sp, 200), size=int(rng.integers(0, 4))):
            tt = S.ts((r.t0 + pd.Timedelta(hours=float(k))).floor("s"))
            add.append([r.SUBJECT_ID, r.HADM_ID, "5075", tt, f"{rng.uniform(5, 300):.1f}", "μmol/L"])
            if r.unit == "NICU" or rng.random() < 0.2:
                add.append([r.SUBJECT_ID, r.HADM_ID, "5255", tt, f"{rng.uniform(20, 400):.1f}", "μmol/L"])
        if r.unit == "NICU" and rng.random() < 0.5:
            tt = S.ts((r.t0 + pd.Timedelta(hours=2)).floor("s"))
            add.append([r.SUBJECT_ID, r.HADM_ID, "5255", tt, f"{rng.uniform(50, 500):.1f}", "μmol/L"])
        for k in rng.uniform(0, min(sp, 100), size=int(rng.integers(0, 3))):
            tt = S.ts((r.t0 + pd.Timedelta(hours=float(k))).floor("s"))
            fi = rng.choice(["0.4", "0.6", "1", "0", "5", "45", "250", "abc"])
            add.append([r.SUBJECT_ID, r.HADM_ID, "5222", tt, fi, "%"])
            add.append([r.SUBJECT_ID, r.HADM_ID, "5239", tt, f"{rng.uniform(40, 200):.0f}", "mmHg"])
    ad = pd.DataFrame(add, columns=["SUBJECT_ID", "HADM_ID", "ITEMID", "CHARTTIME", "VALUE", "VALUEUOM"])
    ad["VALUENUM"] = pd.to_numeric(ad["VALUE"], errors="coerce").map(lambda v: "" if pd.isna(v) else f"{v:g}")
    lb = pd.concat([lb, ad[lb.columns]], ignore_index=True)
    lb.to_csv(raw / "LABEVENTS.csv", index=False)
    proj["raw_rows"]["LABEVENTS"] = len(lb)
    proj["raw_sha"]["LABEVENTS"] = S.sha(raw / "LABEVENTS.csv")
    shutil.copy(MAP_V2, proj["interim"] / "maps" / F.CHANNEL_MAP_V2_FILE)
    return proj


@pytest.fixture(scope="module")
def proj(tmp_path_factory):
    return _plant(S.make_project(tmp_path_factory.mktemp("proj"), n_subjects=600, seed=7))


def make_cfg(proj, out_dir: Path, **kw) -> F.FUConfig:
    s2 = S.make_config(M, proj, out_dir, HERE / M.V17_MODULE_FILE)
    ov = {f.name: getattr(s2, f.name) for f in dataclasses.fields(s2)}
    for k in ("raw_dir", "canonical_path", "harmonization_path", "split_path", "complexity_map_path", "maps_dir",
              "out_dir", "v17_module_path", "harm_t0_col", "harm_out_col", "split_col", "chunksize", "verbose", "raw_sha256"):
        ov.pop(k)
    raw = proj["raw"]
    pins = {"LABEVENTS": proj["raw_sha"]["LABEVENTS"],
            **{n: S.sha(raw / f"{n}.csv") for n in ("PATIENTS", "ADMISSIONS", "ICUSTAYS", "SURGERY_INFO")}}
    # Public copy ships the v2.3 map with coverage counts redacted, so the real-run v2 hash pin is overridden.
    base = dict(raw_dir=s2.raw_dir, canonical_path=s2.canonical_path, harmonization_path=s2.harmonization_path,
                split_path=s2.split_path, complexity_map_path=s2.complexity_map_path, maps_dir=s2.maps_dir,
                out_dir=out_dir, stage2_path=HERE / F.STAGE2_FILE, v17_module_path=HERE / M.V17_MODULE_FILE,
                harm_t0_col="t0", harm_out_col="OUTTIME", verbose=False, raw_sha256=pins, s2_overrides=ov,
                channel_map_v2_sha256=_sha(MAP_V2) if MAP_V2.is_file() else F.CHANNEL_MAP_V2_SHA256)
    base.update(kw)
    return F.FUConfig(**base)


@pytest.fixture(scope="module")
def run(proj, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    cfg = make_cfg(proj, out)
    return cfg, F.run_followup(cfg)


# ----------------------------------------------------------------------------- pure rules
def test_fio2_bands_edges():
    b = F.fio2_band(np.array([0, 1, 1.0001, 20.9, 21, 100, 100.01]))
    assert b.tolist() == ["le_1", "le_1", "gt_1_lt_21", "gt_1_lt_21", "21_to_100", "21_to_100", "gt_100"]


def test_pf_category_edges():
    c = F.pf_category(np.array([10, 9.99, 700, 700.01, np.inf]), np.array([50, 50, 50, 50, 0]), 10, 700)
    assert c.tolist() == ["in_gate", "below_gate", "in_gate", "above_gate", "fio2_le_0"]


# ----------------------------------------------------------------------------- end to end
@pytest.mark.skip(reason="needs the original feature_channel_map_v2.csv (not shipped; public copy has v2.3)")
def test_run_publishes_gated_outputs(run):
    cfg, doc = run
    assert doc["release_gate"] == "ok" and doc["spec"] == F.SCRIPT_VERSION
    assert cfg.shared_md.read_text(encoding="utf-8") == F.render_markdown(doc)
    assert not list(cfg.out_dir.glob("*raw_LOCAL*")) and not list(cfg.out_dir.glob("*.pending.*"))
    assert doc["suppression"]["recoverable"] == 0
    assert doc["maps"]["channel"]["sha256"] == cfg.channel_map_v2_sha256
    assert {k for k in doc["raw_tables"]} == {"LABEVENTS", "PATIENTS", "ADMISSIONS", "ICUSTAYS", "SURGERY_INFO"}


def _walk_counts(o):
    if isinstance(o, dict):
        for k, v in o.items():
            if k.startswith("n_") and not isinstance(v, (dict, list)):
                yield v
            elif k not in ("disclosure",):
                yield from _walk_counts(v)


@pytest.mark.skip(reason="needs the original feature_channel_map_v2.csv (not shipped; public copy has v2.3)")
def test_no_published_small_counts(run):
    _, doc = run
    body = {k: v for k, v in doc.items() if k in ("bilirubin", "fio2_format")}
    vals = list(_walk_counts(body))
    assert vals and all(v in ("<10", "S") or (isinstance(v, int) and v >= 10) for v in vals)


@pytest.mark.skip(reason="needs the original feature_channel_map_v2.csv (not shipped; public copy has v2.3)")
def test_identities_and_independent_counts(proj, run):
    cfg, doc = run
    b = doc["bilirubin"]
    a, n, nn = (b["all_train"]["n_encounters"], b["by_first_careunit"]["NICU"]["n_encounters"],
                b["by_first_careunit"]["non_NICU"]["n_encounters"])
    if all(isinstance(x, int) for x in (a, n, nn)):
        assert a == n + nn
    assert a == int((proj["enc"]["split"] == "train").sum())
    p = doc["fio2_format"]["pairs_percent_rule"]
    if all(isinstance(p.get(k), int) for k in ("n_pairs", "n_in_gate", "n_out_of_gate")):
        assert p["n_pairs"] == p["n_in_gate"] + p["n_out_of_gate"]
    # planted FiO2 values 0.4 / 0.6 / 1 / 0 fall in the <= 1 band, 250 above 100
    bands = doc["fio2_format"]["rows_5222"]["numeric_by_band"]
    assert bands["n_le_1"] not in (0,) and bands["n_gt_100"] != 0


@pytest.mark.skip(reason="needs the original feature_channel_map_v2.csv (not shipped; public copy has v2.3)")
def test_release_gate_catches_tamper(run, tmp_path):
    cfg, doc = run
    j, m = tmp_path / "x.json", tmp_path / "x.md"
    d = json.loads(json.dumps(doc))
    d["bilirubin"]["all_train"]["n_encounters"] = 5
    j.write_text(json.dumps(d), encoding="utf-8")
    m.write_text(F.render_markdown(d), encoding="utf-8")
    A = M.load_v17(HERE / M.V17_MODULE_FILE)
    fails = F.release_gate(M, A, j, m, "ok")
    assert any("below 10" in f for f in fails)


@pytest.mark.parametrize("what", ["stage2", "map", "patients"])
def test_hash_gates_stop(proj, tmp_path, what):
    kw = {}
    if what == "stage2":
        kw["stage2_sha256"] = "0" * 64
    elif what == "map":
        kw["channel_map_v2_sha256"] = "0" * 64
    cfg = make_cfg(proj, tmp_path, **kw)
    if what == "patients":
        cfg.raw_sha256 = dict(cfg.raw_sha256, PATIENTS="0" * 64)
    with pytest.raises(F.FollowupError):
        F.run_followup(cfg)
    assert not cfg.shared_json.exists()
