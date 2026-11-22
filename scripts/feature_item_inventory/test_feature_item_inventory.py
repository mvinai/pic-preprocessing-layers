"""pytest suite for feature_item_inventory.py. SYNTHETIC DATA ONLY.

Inputs come from synthetic_inventory_data.py (no real PIC file is read). Synthetic inputs get their own SHA-256 values
and expected counts through InvConfig test overrides, so every gate still runs. The duration audit v1.7 module is
loaded from ../duration_audit_v2 after its pinned SHA-256 is verified.

Run:  python -m pytest -q test_feature_item_inventory.py
"""
from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import json
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
V17_DIR = HERE.parent / "duration_audit_v2"
sys.path.insert(0, str(HERE))
import feature_item_inventory as F  # noqa: E402
import synthetic_inventory_data as S  # noqa: E402

V17 = V17_DIR / F.V17_MODULE_FILE


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def make_cfg(proj: dict, out: Path, **kw) -> F.InvConfig:
    base = dict(raw_dir=proj["raw"], canonical_path=proj["canonical"], harmonization_path=proj["harmonization"],
                split_path=proj["split"], out_dir=Path(out), v17_module_path=V17,
                pinned_sha256={"canonical": sha(proj["canonical"]), "harmonization": sha(proj["harmonization"]),
                               "split": sha(proj["split"])},
                harm_t0_col="t0", harm_out_col="OUTTIME",
                expected_encounters=proj["expected_encounters"], expected_subjects=proj["expected_subjects"],
                expected_split_encounters=proj["split_counts"], expected_eligible=proj["eligible"],
                expected_raw_rows=proj["raw_rows"], chunksize=700, verbose=False)
    base.update(kw)
    return F.InvConfig(**base)


@pytest.fixture(scope="session")
def A():
    return F.load_v17(V17)


@pytest.fixture(scope="session")
def project(tmp_path_factory):
    return S.make_project(tmp_path_factory.mktemp("syn"))


@pytest.fixture(scope="session")
def e2e(project, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    cfg = make_cfg(project, out, verbose=True)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        doc = F.run_inventory(cfg, sentinels=[S.SENTINEL_SUBJECT])
    local = json.loads(cfg.local_json.read_text(encoding="utf-8"))
    return {"cfg": cfg, "doc": doc, "out": out, "stdout": buf.getvalue(), "raw": local["raw_counts"]}


def rawc(e2e, A, *parts):
    return e2e["raw"][A.P(*parts)]


def key_itemid(r):
    return F.vnorm_id(pd.Series([r["ITEMID"]]))[0] or F.BLANK


def key_drug(r):
    k = f"{r['DRUG_NAME'].strip()} | {r['DRUG_NAME_EN'].strip()}"
    return k if k.strip(" |") else F.BLANK


def key_micro(idc, nc):
    def f(r):
        k = f"{r[idc].strip()} | {r[nc].strip()}"
        return k if k.strip(" |") else F.BLANK
    return f


REF_SPECS = {
    "CHARTEVENTS": ("CHARTEVENTS", key_itemid, "CHARTTIME", 0.0),
    "LABEVENTS": ("LABEVENTS", key_itemid, "CHARTTIME", 1.0),
    "OUTPUTEVENTS": ("OUTPUTEVENTS", key_itemid, "CHARTTIME", 0.0),
    "INPUTEVENTS": ("INPUTEVENTS", lambda r: f"unit:{r['AMOUNTUOM'].strip()}", "CHARTTIME", 0.0),
    "PRESCRIPTIONS": ("PRESCRIPTIONS", key_drug, "STARTDATE", 0.0),
    "MICROBIOLOGYEVENTS_specimen": ("MICROBIOLOGYEVENTS", key_micro("SPEC_ITEMID", "SPEC_TYPE_DESC"), "CHARTTIME", 0.0),
    "MICROBIOLOGYEVENTS_organism": ("MICROBIOLOGYEVENTS", key_micro("ORG_ITEMID", "ORG_NAME"), "CHARTTIME", 0.0),
}


@pytest.fixture(scope="session")
def reference(project):
    return {inv: S.reference_counts(project, t, kf, tc, off) for inv, (t, kf, tc, off) in REF_SPECS.items()}


# ----------------------------------------------------------------------------- pure helpers
def test_vnorm_id_matches_v17_norm_id(A):
    vals = ["123", " 123.0 ", "'7'", "007", "007.0", "+5.0", "", "nan", "None", "12.00", "abc", "-5.0", "1.5"]
    got = F.vnorm_id(pd.Series(vals)).tolist()
    assert got == [A.norm_id(v) for v in vals]


def test_fast_time_formats():
    t = F.fast_time(pd.Series(["2016-01-01 10:00:00", "2016-01-01T10:30", "", "garbage", "nan"]))
    assert t[0] == pd.Timestamp("2016-01-01 10:00:00") and t[1] == pd.Timestamp("2016-01-01 10:30:00")
    assert t[2:].isna().all()


def test_resolve_is_case_insensitive_and_required():
    assert F.resolve(["subject_id", "X"], ("SUBJECT_ID",), "t") == "subject_id"
    assert F.resolve(["X"], ("Y",), "t", required=False) is None
    with pytest.raises(F.AuditError):
        F.resolve(["X"], ("Y",), "t")


def test_item_sort_key_orders_numeric_then_text():
    ks = ["b | x", "100", "20", F.BLANK]
    assert sorted(ks, key=F.item_sort_key) == ["20", "100", F.BLANK, "b | x"]


# ----------------------------------------------------------------------------- end to end
def test_end_to_end_release_ok(e2e):
    cfg, doc = e2e["cfg"], e2e["doc"]
    assert doc["release_gate"] == "ok" and doc["spec"] == F.SPEC_VERSION
    for p in [cfg.shared_json, cfg.shared_md, cfg.local_json, *(cfg.template(k) for k in F.TEMPLATE_FILES)]:
        assert p.exists()
    assert not list(Path(e2e["out"]).glob("*.pending.*")) and not list(Path(e2e["out"]).glob("*.tmp"))
    shared = [p.name for p in Path(e2e["out"]).iterdir() if p != cfg.local_json]
    assert all("_raw_LOCAL" not in n for n in shared)


def test_all_inventories_present(e2e):
    assert set(e2e["doc"]["inventories"]) == set(F.INVENTORIES)
    assert set(e2e["doc"]["tables"]) == set(F.TABLES)


@pytest.mark.parametrize("inv", list(REF_SPECS))
def test_item_counts_match_independent_reference(e2e, reference, A, inv):
    ref = reference[inv]["items"]
    got_items = set(e2e["doc"]["inventories"][inv]["items"])
    assert got_items == set(ref)
    for it, r in ref.items():
        for k in ("rows", "encounters", "encounters_layer_b", "encounters_pre_icu_6h"):
            assert rawc(e2e, A, "inventories", inv, "items", it, "train", k) == r[k], (inv, it, k)


@pytest.mark.parametrize("inv", ["CHARTEVENTS", "LABEVENTS", "PRESCRIPTIONS", "MICROBIOLOGYEVENTS_specimen"])
def test_row_classes_match_reference_and_reconcile(e2e, reference, A, inv):
    t = F.INV_TABLE[inv]
    cls = reference[inv]["classes"]
    for c, v in cls.items():
        assert rawc(e2e, A, "tables", t, "row_classes", c) == v, (t, c)
    assert rawc(e2e, A, "tables", t, "row_classes", "rows_total") == sum(cls.values())


def test_noise_rows_classified(e2e, A):
    c = lambda k: rawc(e2e, A, "tables", "CHARTEVENTS", "row_classes", k)  # noqa: E731
    assert c("no_cohort_admission") >= 45          # non-cohort subject, blank HADM, other HADM
    assert c("cohort_admission_time_missing") >= 15
    assert rawc(e2e, A, "tables", "LABEVENTS", "row_classes", "cohort_admission_time_missing") >= 15


def test_lab_release_offset_and_windows(e2e, A, project):
    n = len(project["planted"])
    assert n >= 10
    b = ("inventories", "LABEVENTS", "items", S.PLANT_LAB, "train")
    assert rawc(e2e, A, *b, "rows") == 2 * n                 # released at t0 and at t0+24h; the third is released after OUT
    assert rawc(e2e, A, *b, "encounters") == n
    assert rawc(e2e, A, *b, "encounters_layer_b") == n       # the release at t0+24h is outside [t0, t0+24h)
    assert rawc(e2e, A, *b, "encounters_pre_icu_6h") == n    # release at t0-6h in; t0-6h-1s out


def test_train_only_counts(e2e, A, project):
    n_train = project["split_counts"]["train"]
    assert rawc(e2e, A, "denominators", "train_encounters") == n_train
    for inv in F.INVENTORIES:
        t = F.INV_TABLE[inv]
        assert rawc(e2e, A, "tables", t, "train", "encounters") <= n_train
        tot = sum(rawc(e2e, A, "inventories", inv, "items", it, "train", "rows") for it in e2e["doc"]["inventories"][inv]["items"])
        assert tot == rawc(e2e, A, "tables", t, "row_classes", "in_stay_train")


def test_not_in_dictionary_item(e2e, A):
    it = e2e["doc"]["inventories"]["CHARTEVENTS"]["items"].get(S.CHART_NOT_IN_DICT)
    if it is not None:
        assert it["in_dictionary"] is False and it["label"] == ""
        assert rawc(e2e, A, "inventories", "CHARTEVENTS", "rows_not_in_dictionary") == \
            rawc(e2e, A, "inventories", "CHARTEVENTS", "items", S.CHART_NOT_IN_DICT, "train", "rows")


def test_inputevents_without_itemid_keyed_by_unit(e2e):
    items = set(e2e["doc"]["inventories"]["INPUTEVENTS"]["items"])
    assert items and all(k.startswith("unit:") for k in items)


def test_drug_and_micro_fields(e2e):
    rx = e2e["doc"]["inventories"]["PRESCRIPTIONS"]["items"]
    assert rx["维库溴铵 | Vecuronium"]["drug_name_en"] == "Vecuronium"
    org = e2e["doc"]["inventories"]["MICROBIOLOGYEVENTS_organism"]["items"]
    assert F.BLANK in org and "80023 | STAPHYLOCOCCUS AUREUS" in org
    assert org["80023 | STAPHYLOCOCCUS AUREUS"]["org_itemid"] == "80023"


def test_raw_table_hashes_recorded(e2e, project):
    rt = e2e["doc"]["raw_tables"]
    for t in list(F.TABLES) + ["D_ITEMS", "D_LABITEMS"]:
        assert rt[t]["sha256"] == sha(project["raw"] / f"{t}.csv")
    md = e2e["cfg"].shared_md.read_text(encoding="utf-8")
    assert rt["LABEVENTS"]["sha256"] in md


def test_input_hashes_recorded(e2e, project):
    h = e2e["doc"]["input_hashes"]
    assert h == {"canonical": sha(project["canonical"]), "harmonization": sha(project["harmonization"]), "split": sha(project["split"])}


# ----------------------------------------------------------------------------- disclosure
def test_published_counts_never_below_10(e2e):
    doc = e2e["doc"]
    for ptr in doc["disclosure"]["count_fields"]:
        v = e2e["doc"]
        for part in ptr.strip("/").split("/"):
            v = v[part.replace("~1", "/").replace("~0", "~")]
        assert v in ("<10", "S") or (isinstance(v, int) and v >= 10)


def test_suppression_present_and_not_recoverable(e2e):
    s = e2e["doc"]["suppression"]
    assert s["primary_cells"] > 0 and s["recoverable"] == 0
    assert "<10" in e2e["cfg"].shared_md.read_text(encoding="utf-8")


def test_raw_values_of_suppressed_cells_are_0_to_9_or_complement(e2e):
    """v1.7 rule: counts 0-9 are primary-suppressed as <10 (zero included)."""
    doc, raw = e2e["doc"], e2e["raw"]
    for ptr in doc["disclosure"]["count_fields"]:
        v = doc
        for part in ptr.strip("/").split("/"):
            v = v[part.replace("~1", "/").replace("~0", "~")]
        if v == "<10":
            assert 0 <= raw[ptr] <= 9
        elif isinstance(v, int):
            assert v == raw[ptr]


def test_units_hidden_unless_rows_published(e2e):
    for d in e2e["doc"]["inventories"].values():
        for o in d["items"].values():
            if not isinstance(o["train"]["rows"], int):
                assert o["units"] == []


def test_templates_list_only_items_with_10_or_more_encounters(e2e, A):
    cfg, raw = e2e["cfg"], e2e["raw"]
    rows = list(csv.DictReader(io.StringIO(cfg.template("channel").read_text(encoding="utf-8"))))
    ids = {(r["table"], r["itemid"]) for r in rows}
    assert ("CHARTEVENTS", "220045") in ids
    assert ("CHARTEVENTS", "777001") not in ids
    for r in rows:
        assert r["train_encounters"] != "<10"
        inv = r["table"]
        assert raw[A.P("inventories", inv, "items", r["itemid"], "train", "encounters")] >= 10
        assert r["channel"] == "" and r["notes"] == ""


def test_template_headers_and_drug_groups(e2e):
    cfg = e2e["cfg"]
    for k, cols in F.TEMPLATE_COLUMNS.items():
        assert cfg.template(k).read_text(encoding="utf-8").splitlines()[0] == ",".join(cols)
    drows = list(csv.DictReader(io.StringIO(cfg.template("drug").read_text(encoding="utf-8"))))
    names = {r["drug_name_en"] for r in drows}
    assert "Rare drug/x" not in names and "Dopamine" in names
    assert all(r["group"] == "" for r in drows)
    groups = e2e["doc"]["config"]["drug_groups"]
    assert "nmb" in groups and "aed" in groups


def test_search_hits(e2e):
    h = e2e["doc"]["search_hits"]
    assert "CHARTEVENTS:224639" in h["weight"]
    assert "CHARTEVENTS:226559" in h["urine_output"]
    assert "CHARTEVENTS:220739" in h["gcs"]
    assert "PRESCRIPTIONS:维库溴铵 | Vecuronium" in h["nmb"]
    assert "PRESCRIPTIONS:左乙拉西坦 | Levetiracetam" in h["aed"]
    assert "LABEVENTS:51001" in h["endocrine"]


def test_sentinel_absent_and_no_timestamps(e2e):
    text = "".join(p.read_text(encoding="utf-8") for p in e2e["out"].iterdir() if "_raw_LOCAL" not in p.name)
    assert S.SENTINEL_SUBJECT not in text
    assert not F.TIMESTAMP_RE.search(text)


def test_items_listed_in_key_order(e2e):
    for d in e2e["doc"]["inventories"].values():
        ks = list(d["items"])
        assert ks == sorted(ks, key=F.item_sort_key)


def test_console_summary_short(e2e):
    out = e2e["stdout"]
    assert "release_gate=ok" in out and "recoverable=0" in out
    summary = out[out.index("==="):]
    assert len(summary.splitlines()) <= 15


def test_deterministic_rerun(project, e2e, tmp_path):
    cfg = make_cfg(project, tmp_path)
    F.run_inventory(cfg, sentinels=[S.SENTINEL_SUBJECT])
    assert cfg.shared_json.read_bytes() == e2e["cfg"].shared_json.read_bytes()
    assert cfg.shared_md.read_bytes() == e2e["cfg"].shared_md.read_bytes()


# ----------------------------------------------------------------------------- release gate tampering
def _copy_shared(e2e, dst: Path):
    cfg = e2e["cfg"]
    j, m = dst / cfg.shared_json.name, dst / cfg.shared_md.name
    shutil.copy(cfg.shared_json, j)
    shutil.copy(cfg.shared_md, m)
    t = {}
    for k in F.TEMPLATE_FILES:
        t[k] = dst / cfg.template(k).name
        shutil.copy(cfg.template(k), t[k])
    return j, m, t


def _first(doc, pred):
    for ptr in doc["disclosure"]["count_fields"]:
        v = doc
        for part in ptr.strip("/").split("/"):
            v = v[part.replace("~1", "/").replace("~0", "~")]
        if pred(v):
            return ptr
    raise AssertionError("no matching cell")


def test_gate_passes_on_untouched_copy(e2e, A, tmp_path):
    j, m, t = _copy_shared(e2e, tmp_path)
    assert F.release_gate(A, j, m, t, "ok", [S.SENTINEL_SUBJECT]) == []


def _tamper_json(e2e, A, tmp_path, fn, rerender=True):
    j, m, t = _copy_shared(e2e, tmp_path)
    doc = json.loads(j.read_text(encoding="utf-8"))
    fn(doc)
    j.write_text(json.dumps(doc, ensure_ascii=False), encoding="utf-8")
    if rerender:
        m.write_text(F.render_markdown(doc), encoding="utf-8")
        for k, text in F.render_templates(doc).items():
            t[k].write_text(text, encoding="utf-8")
    return F.release_gate(A, j, m, t, "ok", [S.SENTINEL_SUBJECT])


def test_gate_rejects_small_count(e2e, A, tmp_path):
    def f(doc):
        A.set_ptr(doc, _first(doc, lambda v: isinstance(v, int)), 5)
    assert any("below 10" in x for x in _tamper_json(e2e, A, tmp_path, f))


def test_gate_rejects_unsuppressed_primary(e2e, A, tmp_path):
    def f(doc):
        A.set_ptr(doc, _first(doc, lambda v: v == "<10"), 12)
    assert _tamper_json(e2e, A, tmp_path, f)


def test_gate_rejects_counter_change(e2e, A, tmp_path):
    def f(doc):
        doc["suppression"]["primary_cells"] += 1
    assert any("counters" in x for x in _tamper_json(e2e, A, tmp_path, f))


def test_gate_rejects_undeclared_count_field(e2e, A, tmp_path):
    def f(doc):
        doc["tables"]["CHARTEVENTS"]["train"]["rows"] = 50
    assert any("schema" in x for x in _tamper_json(e2e, A, tmp_path, f))


def test_gate_rejects_min_max_and_timestamps(e2e, A, tmp_path):
    def f(doc):
        doc["inventories"]["CHARTEVENTS"]["items"]["220045"]["max"] = 1
        doc["inventories"]["CHARTEVENTS"]["items"]["220045"]["label"] = "2016-01-01 10:00"
    fails = _tamper_json(e2e, A, tmp_path, f)
    assert any("min/max" in x for x in fails) and any("timestamp" in x for x in fails)


def test_gate_rejects_units_on_suppressed_item(e2e, A, tmp_path):
    def f(doc):
        for o in doc["inventories"]["CHARTEVENTS"]["items"].values():
            if not isinstance(o["train"]["rows"], int):
                o["units"] = ["u"]
                return
        raise AssertionError
    assert any("units" in x for x in _tamper_json(e2e, A, tmp_path, f))


def test_gate_rejects_sentinel(e2e, A, tmp_path):
    def f(doc):
        doc["config"]["note"] = S.SENTINEL_SUBJECT
    assert any("sentinel" in x for x in _tamper_json(e2e, A, tmp_path, f))


def test_gate_rejects_edited_markdown(e2e, A, tmp_path):
    j, m, t = _copy_shared(e2e, tmp_path)
    m.write_text(m.read_text(encoding="utf-8") + "\nextra\n", encoding="utf-8")
    assert any("markdown" in x for x in F.release_gate(A, j, m, t, "ok"))


def test_gate_rejects_template_with_rare_item(e2e, A, tmp_path):
    j, m, t = _copy_shared(e2e, tmp_path)
    t["channel"].write_text(t["channel"].read_text(encoding="utf-8") + "CHARTEVENTS,777001,Rare,,,<10,<10,,,,,,,,\n", encoding="utf-8")
    fails = F.release_gate(A, j, m, t, "ok")
    assert any("template channel" in x for x in fails) and any("fewer than 10" in x for x in fails)


def test_gate_rejects_wrong_status_and_raw_local_name(e2e, A, tmp_path):
    j, m, t = _copy_shared(e2e, tmp_path)
    assert any("header" in x for x in F.release_gate(A, j, m, t, "pending"))
    j2 = tmp_path / "feature_item_inventory_raw_LOCAL.json"
    shutil.copy(j, j2)
    assert any("_raw_LOCAL" in x for x in F.release_gate(A, j2, m, t, "ok"))


def test_publish_gate_failure_leaves_no_shared_file(project, tmp_path, monkeypatch):
    cfg = make_cfg(project, tmp_path)
    monkeypatch.setattr(F, "release_gate", lambda *a, **k: ["forced failure"])
    with pytest.raises(F.ReleaseGateError):
        F.run_inventory(cfg)
    left = [p.name for p in tmp_path.iterdir()]
    assert left == ["feature_item_inventory_raw_LOCAL.json"]


# ----------------------------------------------------------------------------- input gates and exit codes
def _main_with(cfg, monkeypatch, capsys):
    monkeypatch.setattr(F, "build_config", lambda argv=None: cfg)
    rc = F.main([])
    return rc, capsys.readouterr().err


@pytest.mark.parametrize("which", ["canonical", "harmonization", "split"])
def test_hash_gate_exit_2(project, tmp_path, monkeypatch, capsys, which):
    cfg = make_cfg(project, tmp_path)
    cfg.pinned_sha256[which] = "0" * 64
    rc, err = _main_with(cfg, monkeypatch, capsys)
    assert rc == 2 and "SHA-256" in err
    assert not list(tmp_path.iterdir())


def test_v17_module_hash_exit_2(project, tmp_path, monkeypatch, capsys):
    cfg = make_cfg(project, tmp_path, v17_module_sha256="1" * 64)
    rc, _ = _main_with(cfg, monkeypatch, capsys)
    assert rc == 2


def test_pinned_constants():
    assert F.PINNED_SHA256["harmonization"] == "ffd1b8ae6a13e491be2f7a6631eb4efc1f323656795b80b91e105a0644209d6d"
    assert F.PINNED_SHA256["split"] == "60196792457af718155a878f25b3d502ae12dbc71ccabb4c2709518c62e40e31"
    assert F.EXPECTED_ELIGIBLE == {}
    assert F.InvConfig.__dataclass_fields__["expected_encounters"].default is None


@pytest.mark.parametrize("kw", [{"expected_encounters": 1}, {"expected_subjects": 1},
                                {"expected_split_encounters": {"train": 1, "validation": 1, "test": 1}}])
def test_cohort_mismatch_exit_2(project, tmp_path, monkeypatch, capsys, kw):
    rc, err = _main_with(make_cfg(project, tmp_path, **kw), monkeypatch, capsys)
    assert rc == 2 and "cohort" in err


def test_eligibility_mismatch_exit_2(project, tmp_path, monkeypatch, capsys):
    el = dict(project["eligible"])
    el[336] += 1
    rc, err = _main_with(make_cfg(project, tmp_path, expected_eligible=el), monkeypatch, capsys)
    assert rc == 2 and "eligible" in err


def test_raw_row_count_mismatch_exit_2(project, tmp_path, monkeypatch, capsys):
    rr = dict(project["raw_rows"])
    rr["LABEVENTS"] += 1
    rc, err = _main_with(make_cfg(project, tmp_path, expected_raw_rows=rr), monkeypatch, capsys)
    assert rc == 2 and "LABEVENTS" in err
    assert not (tmp_path / "feature_item_inventory.json").exists()


def test_missing_raw_table_exit_2(project, tmp_path, monkeypatch, capsys):
    raw2 = tmp_path / "raw"
    shutil.copytree(project["raw"], raw2)
    (raw2 / "OUTPUTEVENTS.csv").unlink()
    rc, err = _main_with(make_cfg(project, tmp_path / "o", raw_dir=raw2), monkeypatch, capsys)
    assert rc == 2 and "OUTPUTEVENTS" in err


def test_release_gate_failure_exit_3(project, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(F, "release_gate", lambda *a, **k: ["forced failure"])
    rc, err = _main_with(make_cfg(project, tmp_path), monkeypatch, capsys)
    assert rc == 3 and "RELEASE GATE FAILED" in err


def test_cli_defaults_come_from_env_and_flags(monkeypatch):
    monkeypatch.delenv("PIC_RAW_DIR", raising=False)
    monkeypatch.delenv("SPLIT_FILE_PATH", raising=False)
    with pytest.raises(F.AuditError, match="raw-dir"):
        F.build_config([])
    monkeypatch.setenv("PIC_RAW_DIR", "EnvRaw")
    monkeypatch.setenv("SPLIT_FILE_PATH", "EnvSplit.csv")
    cfg = F.build_config([])
    assert str(cfg.raw_dir) == "EnvRaw" and str(cfg.split_path) == "EnvSplit.csv"
    cfg2 = F.build_config(["--raw-dir", "X", "--split-file", "Y.csv", "--harm-t0-col", "t0", "--harm-out-col", "OUTTIME"])
    assert str(cfg2.raw_dir) == "X" and cfg2.harm_t0_col == "t0" and cfg2.harm_out_col == "OUTTIME"
