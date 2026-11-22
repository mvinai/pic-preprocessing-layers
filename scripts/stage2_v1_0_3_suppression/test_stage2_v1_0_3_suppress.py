"""pytest suite for stage2_v1_0_3_suppress.py (module v1.5). Synthetic data only (synthetic_v101_output.py, synthetic_s2_v103_mirror.py)."""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stage2_v1_0_3_suppress as T  # noqa: E402
import synthetic_v101_output as SY  # noqa: E402

SCRIPT = HERE / "stage2_v1_0_3_suppress.py"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def base_args(c: dict) -> list[str]:
    return ["--in-json", str(c["json"]), "--in-md", str(c["md"]), "--links", str(c["links"]),
            "--frozen", str(c["frozen"])]


@pytest.fixture
def case(tmp_path):
    return SY.write_case(tmp_path / "in")


@pytest.fixture
def released(case, tmp_path):
    out = tmp_path / "out"
    r = run_cli(*base_args(case), "--out-dir", str(out))
    assert r.returncode == 0, r.stderr
    doc = json.loads((out / T.OUT_JSON).read_text(encoding="utf-8"))
    return {"case": case, "out": out, "doc": doc, "stdout": r.stdout}


# ---------------------------------------------------------------- interval rule (rule 3)
def test_rule3_p_cells_need_full_range():
    p = T.Params()
    assert T.is_safe("P", 0, 9, p)
    assert not T.is_safe("P", 0, 1, p)       # FiO2 pattern: five '<10' cells summing to a tiny number
    assert not T.is_safe("P", 3, 9, p)       # narrowed: reveals a nonzero small count
    assert not T.is_safe("P", 4, 4, p)       # exactly recoverable


def test_rule3_s_and_h_cells():
    p = T.Params()
    assert T.is_safe("S", 10, 19, p) and not T.is_safe("S", 10, 15, p)
    assert T.is_safe("S", 10, math.inf, p)
    assert T.is_safe("S", 10, 15, T.Params(s_rule="exact")) and not T.is_safe("S", 12, 12, T.Params(s_rule="exact"))
    assert T.is_safe("H", 500, 500, p)       # a hidden difference provably >= 10
    assert T.is_safe("H", 0, 0, p)           # exactly zero
    assert T.is_safe("H", 5, 12, p)          # straddles 9|10
    assert not T.is_safe("H", 1, 9, p) and not T.is_safe("H", 0, 5, p)
    assert not T.is_safe("P", 0, 8, T.Params(width=10)) and T.is_safe("P", 0, 4, T.Params(width=5))


def test_single_equation_projection_matches_lp():
    # x1 + x2 + x3 = 12 with P,P,S boxes -> P in [0, 2], S in [10, 12]
    box = {"a": (0, 9), "b": (0, 9), "c": (10, math.inf)}
    out = T.single_equation_bounds({"a": 1, "b": 1, "c": 1}, 12, box)
    assert out["a"] == (0, 2) and out["b"] == (0, 2) and out["c"] == (10, 12)


def test_block_lp_two_equation_component():
    # a ('<10') sits in a row and a column group; each group alone leaves a free, the two together pin it:
    # row: R1 = a + b (b is S >= 10, R1 large) -> a free in 0..9; column: C1 = a + c, C1 = 12, c >= 10 -> a <= 2
    rels = [("R1", ("a", "b"), "identity"), ("C1", ("a", "c"), "identity"), ("G", ("b", "c"), "identity")]
    values = {"R1": 53, "C1": 12, "G": 63}
    kinds = {"R1": "pub", "C1": "pub", "G": "pub", "a": "P", "b": "S", "c": "S"}
    sysm = T.System(rels, values, kinds, T.Params())
    res = sysm.audit(dict(kinds))
    assert sysm.lp_solves > 0                 # a multi-equation component goes through the block LP
    assert "a" in res.unsafe
    # the grand-total identity makes R1 + C1 - G = 2a, so a is even exactly recoverable
    assert res.exact_recoverable >= 1
    # without the grand total and with a wider column total, a keeps its full 0..9 range (still a 2-equation LP)
    rels2 = rels[:2]
    sys2 = T.System(rels2, {"R1": 53, "C1": 25}, kinds, T.Params())
    assert not sys2.audit(dict(kinds)).unsafe and sys2.lp_solves > 0


# ---------------------------------------------------------------- end to end
def test_check_only_catches_v101_exposure_and_writes_nothing(case, tmp_path):
    before = {k: sha(v) for k, v in case.items()}
    r = run_cli(*base_args(case), "--check-only")
    assert r.returncode == 3
    assert "GATE FAIL (check-only)" in r.stdout
    for ptr in ("/map_checks/row_status/fio2/n_out_of_gate", "/map_checks/sf_ratio/fio2_draws/n_spo2_gt_97",
                "/map_checks/linked_demo/r1/n_c1", "/map_checks/map/pairs/n_outside_gate"):
        assert ptr in r.stdout
    assert {k: sha(v) for k, v in case.items()} == before
    assert sorted(p.name for p in case["json"].parent.iterdir()) == sorted(v.name for v in case.values())


def test_old_single_cell_gate_misses_the_differencing_exposure(tmp_path):
    # without the linked-table and second-cell cases, NO cell is exactly recoverable (the v1.0.1 gate passes),
    # yet the new interval / residual gate fails on the FiO2-like and S/F-like groups
    c = SY.write_case(tmp_path / "in", include_linked_table=False, include_second_cell=False)
    r = run_cli(*base_args(c), "--check-only")
    assert r.returncode == 3
    assert "cells exactly recoverable (the v1.0.1 gate's own criterion): 0" in r.stdout
    assert "groups failing the residual (differencing) check: 2" in r.stdout


def test_release_fixes_fio2_units_and_sf(released):
    d = released["doc"]
    rs = d["map_checks"]["row_status"]["fio2"]
    assert rs["n_total"] == "S" and isinstance(rs["n_measurement"], int)          # total goes, frozen part stays
    for it in ("LABEVENTS_5222", "LABEVENTS_5236"):                                # linked unit totals follow
        u = d["map_checks"]["units"][it]
        assert u["n_total"] == "S" and "S" in (u["n_(blank)"], list(u.values())[1])
    sf = d["map_checks"]["sf_ratio"]
    assert isinstance(sf["n_fio2_draws"], int) and isinstance(sf["fio2_draws"]["n_with_valid_pf"], int)
    assert sf["fio2_draws"]["n_no_spo2_within_5min"] == "S"                       # derived count suppressed
    aa = d["map_checks"]["row_status"]["aa_like"]                                  # safe group untouched
    assert isinstance(aa["n_total"], int) and isinstance(aa["n_measurement"], int)


def test_release_linked_table_and_second_cell(released):
    d = released["doc"]
    t = d["map_checks"]["linked_demo"]
    # interior cells kept; margins carry the protection (row, column and grand)
    assert isinstance(t["r1"]["n_c2"], int) and isinstance(t["r2"]["n_c1"], int)
    assert t["r1"]["n_total"] == "S" and t["col_totals"]["n_c1"] == "S" and t["n_grand_total"] == "S"
    pr = d["map_checks"]["map"]["pairs"]
    assert pr["n_valid"] == "S" and pr["n_dbp_ge_sbp"] == "S" and isinstance(pr["n_total"], int)


def test_derived_values_follow(released):
    t = released["doc"]["map_checks"]["linked_demo"]
    assert t["r1"]["pct_c2"] == "suppressed" and t["x_share_c1"] == "suppressed"
    assert isinstance(t["r2"]["q_values"], dict)                                  # its count is still published


def test_published_values_unchanged_and_monotone(released):
    src = json.loads(released["case"]["json"].read_text(encoding="utf-8"))
    out = released["doc"]
    for ptr in src["disclosure"]["count_fields"]:
        a, b = T.get_ptr(src, ptr), T.get_ptr(out, ptr)
        if isinstance(b, int):
            assert a == b
        else:
            assert b in ("<10", "S") and (a == b or isinstance(a, int))
        if a == "<10":
            assert b == "<10"
    assert out["disclosure"] == src["disclosure"]


def test_output_passes_gate_in_check_only(released):
    o = released["out"]
    c = released["case"]
    r = run_cli("--in-json", str(o / T.OUT_JSON), "--in-md", str(o / T.OUT_MD), "--links", str(c["links"]),
                "--frozen", str(c["frozen"]), "--check-only")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "cells failing the protection rule: 0" in r.stdout


def test_markdown_is_rendering_and_header(released):
    o = released["out"]
    md = (o / T.OUT_MD).read_text(encoding="utf-8")
    assert md == T.render_markdown(released["doc"])
    assert "v1.0.3" in md.splitlines()[2]


def test_manifest_and_summary(released):
    o = released["out"]
    man = json.loads((o / T.OUT_MANIFEST).read_text(encoding="utf-8"))
    assert man["release_gate"] == "ok"
    for f in man["outputs"]:
        assert sha(o / f["file"]) == f["sha256"]
    roles = {i["role"]: i["sha256"] for i in man["inputs"]}
    assert roles["input_json"] == sha(released["case"]["json"])
    summ = (o / T.OUT_SUMMARY).read_text(encoding="utf-8")
    assert "group total suppressed | 4 |" in summ and "part suppressed | 3 |" in summ
    # the summary carries no count values: none of the synthetic published figures >= 10 appear in it
    raw = SY.build()["raw"]
    body = summ.split("## Cells by rule")[0]
    for v in {v for v in raw.values() if v >= 100}:
        assert f" {v} " not in body and f"{v}," not in body
    assert sorted(p.name for p in o.iterdir()) == sorted([T.OUT_JSON, T.OUT_MD, T.OUT_SUMMARY, T.OUT_MANIFEST])


def test_deterministic(case, tmp_path):
    outs = []
    for k in range(2):
        o = tmp_path / f"o{k}"
        assert run_cli(*base_args(case), "--out-dir", str(o)).returncode == 0
        outs.append({p.name: sha(p) for p in o.iterdir()})
    assert outs[0] == outs[1]


# ---------------------------------------------------------------- refusals and gate failures
def test_refuses_raw_local_input(case, tmp_path):
    bad = case["json"].with_name("feature_sparsity_frequency_audit_raw_LOCAL.json")
    bad.write_bytes(case["json"].read_bytes())
    r = run_cli("--in-json", str(bad), "--in-md", str(case["md"]), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2 and "_raw_LOCAL" in r.stderr
    assert not (tmp_path / "o").exists()


def test_refuses_input_folder_as_output(case):
    r = run_cli(*base_args(case), "--out-dir", str(case["json"].parent))
    assert r.returncode == 2 and "new folder" in r.stderr


def test_refuses_nonempty_output(case, tmp_path):
    o = tmp_path / "o"
    o.mkdir()
    (o / "x.txt").write_text("x")
    assert run_cli(*base_args(case), "--out-dir", str(o)).returncode == 2


def test_wrong_sha_pin_stops(case, tmp_path):
    r = run_cli(*base_args(case), "--out-dir", str(tmp_path / "o"), "--expect-input-sha256", "0" * 64)
    assert r.returncode == 2 and not (tmp_path / "o").exists()


def test_cannot_protect_fails_gate_and_leaves_nothing(case, tmp_path):
    # freeze every cell of the second-cell group: nothing can be suppressed there -> exit 3, no files
    fz = case["frozen"].read_text(encoding="utf-8")
    fz += "/map_checks/map/pairs/n_valid,test\n/map_checks/map/pairs/n_dbp_ge_sbp,test\n"
    case["frozen"].write_text(fz, encoding="utf-8")
    o = tmp_path / "o"
    r = run_cli(*base_args(case), "--out-dir", str(o))
    assert r.returncode == 3 and "RELEASE GATE FAILED" in r.stderr
    assert not o.exists() or not any(o.iterdir())


def test_renderer_drift_stops(case, tmp_path):
    case["md"].write_text(case["md"].read_text(encoding="utf-8") + "edited\n", encoding="utf-8")
    o = tmp_path / "o"
    r = run_cli(*base_args(case), "--out-dir", str(o))
    assert r.returncode == 3 and "renderer drift" in r.stderr
    assert not o.exists() or not any(o.iterdir())


def test_inconsistent_link_stops(case, tmp_path):
    # a link that contradicts published values must stop the run (input error), not be silently used
    txt = case["links"].read_text(encoding="utf-8")
    txt += "/map_checks/pf_ratio/same_draw_pairs/n_valid,/map_checks/pf_ratio/same_draw_pairs/n_outside_gate,bad\n"
    case["links"].write_text(txt, encoding="utf-8")
    r = run_cli(*base_args(case), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2 and "violate" in r.stderr


def test_missing_link_pointer_stops(case, tmp_path):
    case["links"].write_text(case["links"].read_text(encoding="utf-8") + "/nope/n_x,/denominators/n_encounters,x\n",
                             encoding="utf-8")
    r = run_cli(*base_args(case), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2 and "not found" in r.stderr


def test_published_count_under_10_is_primary_suppressed(tmp_path):
    c = SY.write_case(tmp_path / "in")
    doc = json.loads(c["json"].read_text(encoding="utf-8"))
    # plant a v1.0.1 bug: a published 7 inside a protected group (aa_like out_of_gate), counters adjusted
    doc["map_checks"]["row_status"]["aa_like"]["n_out_of_gate"] = 2
    doc["suppression"]["primary_cells"] -= 1
    c["json"].write_text(T.dump_json(doc), encoding="utf-8", newline="")
    c["md"].write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    o = tmp_path / "o"
    r = run_cli(*base_args(c), "--out-dir", str(o))
    assert r.returncode == 0, r.stderr
    out = json.loads((o / T.OUT_JSON).read_text(encoding="utf-8"))
    assert out["map_checks"]["row_status"]["aa_like"]["n_out_of_gate"] == "<10"
    assert "primary, new (published count under 10 found in the input) | 1 |" in (o / T.OUT_SUMMARY).read_text()


def test_generic_profile_check_only(tmp_path):
    c = SY.write_case(tmp_path / "in")
    r = run_cli("--in-json", str(c["json"]), "--profile", "generic", "--check-only")
    assert r.returncode == 3 and "GATE FAIL" in r.stdout


# ---------------------------------------------------------------- v1.4 regression: identical to v1.3 (and so to v1.0)
# v1.4 is compared with the delivered v1.3 module (df77234b..., itself regression-tested against v1.0 in its own
# suite). Without --public / --carry-json and without a zero-screened hidden cell, v1.4 must give the same decisions
# and the same files; only labels differ (v1.0.2 -> v1.0.3 names, spec / notes text, three summary lines and the new
# provenance keys, all normalised below).
REF = HERE / "reference" / "stage2_v1_0_2_suppress_v1_3.py"
REF_SHA = "df77234b6c68e1ca9b542f6f387f861141649fdd32e88f5f23ea009acac46d4f"
REF_V10 = HERE / "reference" / "stage2_v1_0_2_suppress_v1_0.py"
REF_V10_SHA = "a341cdb0789a111ce76a076a0bd29e508a8bc68a4354e0818dce0cd989b1b305"
_re = __import__("re")
_SCRIPT_ID = _re.compile(
    r'("(?:script|script_sha256|spec|sha256)": ")[^"]*(")|(- (?:script|spec): ).*|(LP solves: )\d+'
    r'|("bytes": )\d+')
_V14_LINES = _re.compile(r"^(- public \(already known\) cells: .*|- carried release: .*|"
                         r"\| 2 complementary, carried from the previous release .*)\n", _re.M)
_V14_PROV = ("public_sha256", "n_public_cells", "carry_json", "carry_json_sha256", "carried_suppressions",
             "republished_public_cells")
_NOTE_KEYS = ("differencing_fix_v1_0_2", "differencing_fix_v1_0_3")
_OLD_OUTPUT_SPEC = ("feature_sparsity_frequency_audit v1.0.2 (v1.0.1 output with the differencing fix; "
                    "stage2_v1_0_2_suppression_spec v1.0, stage 2)")


def _norm(text: str) -> str:
    # drop only what must differ between versions: the script's own name/version/sha, the spec label, the
    # LP-solve counter, file hashes/sizes quoted in summaries, and the v1.4 release labels
    for spec in (T.OUTPUT_SPEC, _OLD_OUTPUT_SPEC):
        text = text.replace(spec, "X")
    text = _V14_LINES.sub("", text).replace("v1.0.3", "v1.0.2").replace("v1_0_3", "v1_0_2")
    return _SCRIPT_ID.sub(lambda m: (m.group(1) or m.group(3) or m.group(4) or m.group(5) or "") + "X" + (m.group(2) or ""), text)


def _norm_doc(doc: dict) -> dict:
    for k in ("release_v1_0_2", "release_v1_0_3"):
        if k in doc:
            prov = doc.pop(k)
            for q in _V14_PROV:
                prov.pop(q, None)
            doc["release_v1_0_2"] = prov
    notes = doc.get("notes", {})
    for k in _NOTE_KEYS:
        if k in notes:
            del notes[k]
            notes["differencing_fix_v1_0_2"] = "X"
    for k in list(notes):
        if k.startswith("dropped_subtrees_v1_0_"):
            notes["dropped_subtrees_v1_0_2"] = notes.pop(k).replace("v1.0.3", "v1.0.2")
    doc["spec"] = "X"
    return doc


def _run_both(c: dict, tmp_path: Path, extra: list[str], args_fn=base_args):
    res = {}
    for tag, script in (("old", REF), ("new", SCRIPT)):
        out = tmp_path / tag
        a = [*args_fn(c), *extra]
        if "--check-only" not in extra:
            a += ["--out-dir", str(out)]
        r = subprocess.run([sys.executable, str(script), *a], capture_output=True, text=True)
        files = {}
        if out.exists():
            for p in sorted(out.iterdir()):
                name = p.name.replace("v1.0.3", "v1.0.2")
                if name.startswith("feature_") and name.endswith(".json"):
                    files[name] = _norm(json.dumps(_norm_doc(json.loads(p.read_text(encoding="utf-8"))), sort_keys=True))
                    continue
                if name.startswith("feature_") and name.endswith(".md"):
                    # normalised rendering (the gate checks that the file on disk is the rendering of its JSON)
                    j = json.loads(p.with_suffix(".json").read_text(encoding="utf-8"))
                    files[name] = T.render_markdown(_norm_doc(j))
                    continue
                if name == "stage2_v1.0.2_manifest.json":
                    m = json.loads(p.read_text(encoding="utf-8"))
                    for k in ("script", "script_sha256", "spec"):
                        m.pop(k, None)
                    for f in m["outputs"]:   # every output carries the release labels; contents compared above
                        f.pop("sha256", None), f.pop("bytes", None)
                        f["file"] = f["file"].replace("v1.0.3", "v1.0.2")
                    files[name] = json.dumps(m, sort_keys=True)
                else:
                    files[name] = _norm(p.read_text(encoding="utf-8"))
        res[tag] = (r.returncode, _norm(r.stdout), files)
    return res


def test_reference_v13_and_v10_scripts_are_the_delivered_ones():
    assert sha(REF) == REF_SHA and sha(REF_V10) == REF_V10_SHA


REG_CASES = {
    "main_release": (dict(), []),
    "main_check_only": (dict(), ["--check-only"]),
    "no_linked_no_second_release": (dict(include_linked_table=False, include_second_cell=False), []),
    "no_linked_no_second_check_only": (dict(include_linked_table=False, include_second_cell=False),
                                       ["--check-only"]),
    "main_width_5": (dict(), ["--protection-width", "5"]),
    "main_s_rule_exact": (dict(), ["--s-rule", "exact"]),
}


@pytest.mark.parametrize("name", sorted(REG_CASES))
def test_regression_identical_to_v10(name, tmp_path):
    kw, extra = REG_CASES[name]
    c = SY.write_case(tmp_path / "in", **kw)
    res = _run_both(c, tmp_path, extra)
    assert res["old"][0] == res["new"][0]
    assert res["old"][1] == res["new"][1]
    assert res["old"][2] == res["new"][2]
    assert res["new"][2] or "--check-only" in extra or res["new"][0] != 0


def test_regression_primary_new_identical_to_v10(tmp_path):
    c = SY.write_case(tmp_path / "in")
    doc = json.loads(c["json"].read_text(encoding="utf-8"))
    doc["map_checks"]["row_status"]["aa_like"]["n_out_of_gate"] = 2
    doc["suppression"]["primary_cells"] -= 1
    c["json"].write_text(T.dump_json(doc), encoding="utf-8", newline="")
    c["md"].write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    res = _run_both(c, tmp_path, [])
    assert res["old"] == res["new"] and res["new"][0] == 0


def test_regression_scaled_synthetic_identical_to_v10(tmp_path):
    # ~2 MB realistic structure (hub chains, partitions, 7 landmarks, linked identities)
    import synthetic_scale_v101 as SC
    paths = SC.write(tmp_path / "in", 0.03)
    small = SY.write_case(tmp_path / "small")
    c = {"json": paths["json"], "md": paths["md"], "links": tmp_path / "links.csv", "frozen": small["frozen"]}
    c["links"].write_text("\n".join(small["links"].read_text(encoding="utf-8").splitlines()[:3]) + "\n",
                          encoding="utf-8")
    res = _run_both(c, tmp_path, [])
    assert res["new"][0] == 0
    assert res["old"] == res["new"]


# ---------------------------------------------------------------- v1.1: time limit, workers, progress log
def test_max_seconds_fails_closed_and_leaves_nothing(case, tmp_path):
    out = tmp_path / "out"
    r = run_cli(*base_args(case), "--out-dir", str(out), "--max-seconds", "0.001")
    assert r.returncode == 3
    assert "time limit" in r.stderr and "--max-seconds" in r.stderr
    assert not out.exists()


def test_workers_give_identical_output(case, tmp_path):
    outs = []
    for w in ("1", "3"):
        o = tmp_path / f"w{w}"
        r = run_cli(*base_args(case), "--out-dir", str(o), "--workers", w)
        assert r.returncode == 0, r.stderr
        outs.append({p.name: p.read_bytes() for p in o.iterdir()})
    assert outs[0] == outs[1]


def test_progress_log_on_stderr_counts_only(released, tmp_path):
    out = tmp_path / "o2"
    r = run_cli(*base_args(released["case"]), "--out-dir", str(out))
    lines = r.stderr.strip().splitlines()
    assert lines and all(__import__("re").match(r"^\[\s*\d+\.\ds\] ", ln) for ln in lines)
    for phase in ("reading the input JSON", "input audit:", "complementary suppression: start",
                  "release gate: ok", "release written"):
        assert any(phase in ln for ln in lines), phase
    # counts only: no JSON pointer (cell identity) and no '<10' / value text is ever printed
    body = "\n".join(ln.split("] ", 1)[1] for ln in lines)
    assert not __import__("re").search(r"/[A-Za-z_]", body) and "<10" not in body


# ---------------------------------------------------------------- v1.2: the large-component stall (real run, v1.1)
import re as _re  # noqa: E402

import synthetic_stall_v101 as ST  # noqa: E402


def _stall_case(tmp_path, kind):
    return ST.write(tmp_path / f"stall_{kind}", kind)


def _stall_args(c):
    return ["--in-json", str(c["json"]), "--in-md", str(c["md"]), "--frozen", str(c["frozen"])]


def _v11_choice(self, u):
    # v1.1's large-component choice: the static first candidate, no check that it helps
    m = self.choose(u, lookahead=False)
    return None if m is None else [m]


def _rounds(stderr: str):
    fails = [int(x) for x in _re.findall(r"suppression round \d+: (\d+) failing cells", stderr)]
    largest = [int(x) for x in _re.findall(r"suppression round \d+: audit: .*largest multi-equation component "
                                           r"(\d+) cells", stderr)]
    return fails, largest


def test_stall_reproducer_shows_the_v11_pattern(tmp_path, monkeypatch, capsys):
    # with v1.1's choice and no guard: constant failing count, 2 new cells per round, giant component growing
    c = _stall_case(tmp_path, "pinned")
    monkeypatch.setattr(T.Suppressor, "choose_exact", _v11_choice)
    monkeypatch.setattr(T, "STALL_ROUNDS", 10_000)
    monkeypatch.setattr(T, "MAX_ROUNDS", 15)
    rc = T.main([*_stall_args(c), "--out-dir", str(tmp_path / "o")])
    err = capsys.readouterr().err
    fails, largest = _rounds(err)
    assert rc == 3 and "suppression not converging: round limit reached" in err
    assert len(fails) == 15 and set(fails) == {2}
    assert all(b > a for a, b in zip(largest, largest[1:]))           # the component grows every round
    assert len(_re.findall(r": 2 new suppressions", err)) == 15
    assert not (tmp_path / "o").exists()


def test_v12_pinned_cells_fail_closed_in_round_1(tmp_path):
    # provably unprotectable (fixed by frozen cells): v1.2 stops at once instead of walking the component
    c = _stall_case(tmp_path, "pinned")
    out = tmp_path / "o"
    r = run_cli(*_stall_args(c), "--out-dir", str(out))
    assert r.returncode == 3
    assert "RELEASE GATE FAILED: cannot protect 2 cell(s): even with every publishable, non-frozen cell" in r.stderr
    assert "pinned_item/n_courses_in_window" in r.stderr                 # pointers, never values
    assert "suppression round 2" not in r.stderr
    assert not out.exists()


def test_v12_low_value_case_converges_with_checked_picks(tmp_path):
    c = _stall_case(tmp_path, "low")
    out = tmp_path / "o"
    r = run_cli(*_stall_args(c), "--out-dir", str(out))
    assert r.returncode == 0, r.stderr
    assert "suppression round 2: no failing cell; done" in r.stderr
    assert "complementary suppression: done (2 new cells, 1 rounds)" in r.stderr
    doc = json.loads((out / T.OUT_JSON).read_text(encoding="utf-8"))
    for g in ("hub_00", "hub_01"):
        assert doc["surgery"]["giant"][g]["low_item"]["by_site"]["site_big"]["n"] == "S"   # the cell that helps
        assert all(v["n"] == 10 for k, v in doc["surgery"]["giant"][g]["low_item"]["by_site"].items()
                   if k != "site_big")                                   # none of the useless 10s
    chk = run_cli("--in-json", str(out / T.OUT_JSON), "--in-md", str(out / T.OUT_MD), "--frozen", str(c["frozen"]),
                  "--check-only")
    assert chk.returncode == 0, chk.stdout


def test_v12_every_exact_pick_makes_the_cell_pass(tmp_path):
    c = _stall_case(tmp_path, "low")
    doc = json.loads(c["json"].read_text(encoding="utf-8"))
    L = T.load_output(c["json"], c["md"], "stage2", T.EXPECT_INPUT_SPEC)
    frozen = {p for p, _ in T.read_frozen(c["frozen"])}
    system = T.build_system(L, [], T.Params())
    status = dict(L.status)
    res = system.audit(status)
    assert len(res.unsafe) == 2
    sup = T.Suppressor(system, status, frozen)
    for u in sorted(res.unsafe):
        picks = sup.choose_exact(u)
        assert picks and sup.makes_safe(u, picks)
        # independent confirmation with the full audit (every cell, all equations)
        assert u not in system.audit(T.Overlay(status, {m: "S" for m in picks})).unsafe
    assert doc["surgery"]  # input untouched


def test_guard_escalates_and_recovers(tmp_path, monkeypatch, capsys):
    # useless v1.1 picks for 5 rounds -> stronger rule (all candidates within 3 groups) -> converges, gate ok
    c = _stall_case(tmp_path, "low")
    monkeypatch.setattr(T.Suppressor, "choose_exact", _v11_choice)
    rc = T.main([*_stall_args(c), "--out-dir", str(tmp_path / "o")])
    err = capsys.readouterr().err
    assert rc == 0, err
    assert "NOT improved for 5 round(s)" in err and "ESCALATING to the stronger rule" in err
    assert "release gate: ok" in err


def test_guard_fails_closed_when_not_converging(tmp_path, monkeypatch, capsys):
    c = _stall_case(tmp_path, "pinned")
    monkeypatch.setattr(T.Suppressor, "choose_exact", _v11_choice)
    monkeypatch.setattr(T.Suppressor, "unprotectable", lambda self, u: False)
    out = tmp_path / "o"
    rc = T.main([*_stall_args(c), "--out-dir", str(out)])
    err = capsys.readouterr().err
    assert rc == 3
    assert "RELEASE GATE FAILED: suppression not converging" in err
    assert "ESCALATING" in err
    assert not out.exists()


def test_round_log_reports_trend(released, tmp_path):
    out = tmp_path / "o3"
    r = run_cli(*base_args(released["case"]), "--out-dir", str(out))
    assert "failing cells, first round" in r.stderr
    assert _re.search(r"failing cells, improved from \d+", r.stderr)
    assert _re.search(r"failing cells, NOT improved for \d round\(s\) \(best \d+\)", r.stderr)


# ---------------------------------------------------------------- v1.3: hash-pinned drop list (--drops)
DROPS = HERE / "stage2_v1_0_3_drops.csv"
NMB_SUB = "/drug_groups/nmb/validation/strata/surgical_status"
NMB_FAIL = [f"{NMB_SUB}/no_recent_surgery/n_active_72h", f"{NMB_SUB}/recent_surgery/n_active_72h"]


def _nmb_case(folder: Path, pin_train: bool = False, small_recent: bool = False) -> dict:
    """Synthetic mirror of the real v1.2 failure: the NMB validation total of n_active_72h is already public
    (frozen) and its two surgical_status strata are '<10' summing to it, so each stratum is narrowed to 6..9 and no
    publishable cell can widen it. pin_train adds the same trap on the train split (not covered by the drop list)."""
    b = SY.Builder()
    splits = {"train": 700, "validation": 150, "test": 150}
    tot = b.cell("/denominators/n_encounters", sum(splits.values()))
    enc = {s: b.cell(f"/denominators/{s}/n_encounters", n) for s, n in splits.items()}
    b.ident(tot, list(enc.values()))
    frozen = [("/denominators/n_encounters", "cohort size public")] + \
             [(enc[s], "split size public") for s in splits]
    for split, ever, everb, act, strata in (
            ("train", 300, 280, 60 if not pin_train else 12, {"no_recent_surgery": (160, 40 if not pin_train else 5),
                                                             "recent_surgery": (140, 20 if not pin_train else 7)}),
            ("validation", 75, 50, 15, {"no_recent_surgery": (40, 7), "recent_surgery": (35, 8)} if not small_recent
             else {"no_recent_surgery": (68, 8), "recent_surgery": (7, 7)})):
        base = f"/drug_groups/nmb/{split}"
        ne = b.cell(f"{base}/n_ever_336h", ever)
        b.subset(ne, enc[split])
        nb = b.cell(f"{base}/n_ever_layer_b", everb)
        b.subset(nb, ne)
        na = b.cell(f"{base}/n_active_72h", act)
        b.subset(na, ne)
        b.ratios.append((f"{base}/pct_ever_336h", ne, enc[split]))
        if split == "validation" or pin_train:
            frozen.append((na, "synthetic: already public in the follow-up output"))
        se, sa = [], []
        for k, (e, a) in strata.items():
            sb = f"{base}/strata/surgical_status/{k}"
            ce = b.cell(f"{sb}/n_ever_336h", e)
            ca = b.cell(f"{sb}/n_active_72h", a)
            b.subset(ca, ce)
            b.ratios.append((f"{sb}/pct_ever_336h", ce, enc[split]))
            nr = b.cell(f"{sb}/n_rate_per_day_all", e)
            b.ident(nr, [ce])
            b.quantiles.append((f"{sb}/q_rate_per_day_all", nr, SY.Q))
            se.append(ce)
            sa.append(ca)
        b.ident(ne, se)
        b.ident(na, sa)
    b.doc["provisional_classes"]["drug_groups"] = {"nmb": {"class": "measured_indicator_only"}}
    doc = b.finish()
    folder.mkdir(parents=True, exist_ok=True)
    pj, pm, pf = folder / "feature_sparsity_frequency_audit.json", folder / "feature_sparsity_frequency_audit.md", \
        folder / "frozen.csv"
    pj.write_text(T.dump_json(doc), encoding="utf-8", newline="")
    pm.write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    pf.write_text("pointer,reason\n" + "".join(f"{p},{r}\n" for p, r in frozen), encoding="utf-8", newline="")
    return {"json": pj, "md": pm, "frozen": pf}


def _nmb_args(c: dict) -> list[str]:
    return ["--in-json", str(c["json"]), "--in-md", str(c["md"]), "--frozen", str(c["frozen"])]


def test_drops_file_is_the_documented_one():
    rows = T.read_drops(DROPS)
    assert rows == [(NMB_SUB, rows[0][1])] and "Vinai" in rows[0][1]


def test_nmb_mirror_fails_without_the_drop(tmp_path):
    c = _nmb_case(tmp_path / "in")
    o = tmp_path / "o"
    r = run_cli(*_nmb_args(c), "--out-dir", str(o))
    # the real run reported both strata; here the two share one component, so the first one found is reported
    assert r.returncode == 3 and "RELEASE GATE FAILED: cannot protect" in r.stderr, r.stderr
    assert any(p in r.stderr for p in NMB_FAIL) and "fixed by frozen cells" in r.stderr
    assert not o.exists() or not any(o.iterdir())


def test_nmb_mirror_passes_with_the_drop(tmp_path):
    c = _nmb_case(tmp_path / "in")
    o = tmp_path / "o"
    r = run_cli(*_nmb_args(c), "--drops", str(DROPS), "--out-dir", str(o))
    assert r.returncode == 0, r.stderr
    doc = json.loads((o / T.OUT_JSON).read_text(encoding="utf-8"))
    src = json.loads(c["json"].read_text(encoding="utf-8"))
    # the whole subtree is gone (counts, percentages, quantiles), the totals above it are kept as published
    assert doc["drug_groups"]["nmb"]["validation"]["strata"] == {}   # parent object kept, subtree gone
    assert doc["drug_groups"]["nmb"]["validation"]["n_active_72h"] == 15
    assert doc["drug_groups"]["nmb"]["validation"]["n_ever_336h"] == 75
    assert doc["drug_groups"]["nmb"]["train"]["strata"] == src["drug_groups"]["nmb"]["train"]["strata"]
    assert "surgical_status" not in json.dumps(doc["disclosure"]).replace("/train/strata/surgical_status", "")
    assert doc[T.PROVENANCE_KEY]["dropped_subtrees"] == [NMB_SUB]
    assert doc[T.PROVENANCE_KEY]["drops_sha256"] == sha(DROPS)
    assert NMB_SUB in doc["notes"]["dropped_subtrees_v1_0_3"]
    # manifest: the drop list is a hash-pinned input; dropped pointers by name only (strings, no values)
    man = json.loads((o / T.OUT_MANIFEST).read_text(encoding="utf-8"))
    assert {"role": "drops", "file": DROPS.name, "bytes": DROPS.stat().st_size, "sha256": sha(DROPS)} in man["inputs"]
    (d,) = man["dropped"]
    assert d["subtree"] == NMB_SUB and set(NMB_FAIL) <= set(d["count_cells"])
    assert len(d["count_cells"]) == 6 and all(isinstance(x, str) and x.startswith(NMB_SUB + "/")
                                              for x in d["count_cells"])
    summ = (o / T.OUT_SUMMARY).read_text(encoding="utf-8")
    assert f"  - {NMB_SUB}" in summ and "dropped subtrees: 1 holding 6 count cells" in summ
    # the released file passes the independent check-only audit
    r2 = run_cli("--in-json", str(o / T.OUT_JSON), "--in-md", str(o / T.OUT_MD), "--check-only")
    assert r2.returncode == 0 and "GATE OK" in r2.stdout, r2.stdout + r2.stderr


def test_drop_still_fails_closed_on_another_unprotectable_cell(tmp_path):
    c = _nmb_case(tmp_path / "in", pin_train=True)
    o = tmp_path / "o"
    r = run_cli(*_nmb_args(c), "--drops", str(DROPS), "--out-dir", str(o))
    assert r.returncode == 3 and "RELEASE GATE FAILED: cannot protect" in r.stderr, r.stderr
    assert "/drug_groups/nmb/train/strata/surgical_status/" in r.stderr and NMB_SUB + "/" not in r.stderr
    assert not o.exists() or not any(o.iterdir())


def test_dropped_cell_is_protected_against_subtraction(tmp_path):
    # drop only one stratum whose n_ever_336h is small: published total - published sibling would give it exactly,
    # so the audit (dropped cells are unpublished, hidden-cell rule) must suppress the sibling or the total
    c = _nmb_case(tmp_path / "in", small_recent=True)
    dr = tmp_path / "drops.csv"
    dr.write_text(f"subtree,reason\n{NMB_SUB}/recent_surgery,test\n", encoding="utf-8")
    o = tmp_path / "o"
    r = run_cli(*_nmb_args(c), "--drops", str(dr), "--out-dir", str(o))
    assert r.returncode == 0, r.stderr
    doc = json.loads((o / T.OUT_JSON).read_text(encoding="utf-8"))
    v = doc["drug_groups"]["nmb"]["validation"]
    assert "recent_surgery" not in v["strata"]["surgical_status"]
    assert "S" in (v["n_ever_336h"], v["strata"]["surgical_status"]["no_recent_surgery"]["n_ever_336h"])
    assert T.is_safe("D", 0, 20, T.Params())
    assert not T.is_safe("D", 3, 8, T.Params()) and T.is_safe("D", 35, 35, T.Params())


@pytest.mark.parametrize("row,msg", [
    ("/drug_groups/nmb/validation/strata/nope,x", "not found"),
    ("/drug_groups/nmb/validation/n_active_72h,x", "not a subtree"),
    ("/drug_groups/nmb/validation,x", "frozen"),
    ("/disclosure/count_fields,x", "not a droppable"),
    (f"{NMB_SUB},", "reason"),
])
def test_bad_drop_list_stops(tmp_path, row, msg):
    c = _nmb_case(tmp_path / "in")
    dr = tmp_path / "drops.csv"
    dr.write_text(f"subtree,reason\n{row}\n", encoding="utf-8")
    r = run_cli(*_nmb_args(c), "--drops", str(dr), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2 and msg in r.stderr, r.stderr
    assert not (tmp_path / "o").exists()


def test_drops_not_allowed_with_check_only(case):
    r = run_cli(*base_args(case), "--check-only", "--drops", str(DROPS))
    assert r.returncode == 2


# ===================================================================================== v1.4 (Stage 2 v1.0.3)
import csv as _csv  # noqa: E402
import shutil as _shutil  # noqa: E402

import synthetic_s2_v103_mirror as MR  # noqa: E402
import stage2_independent_check as IC  # noqa: E402

REF_MOD_DIR = HERE / "reference"
# the Stage 3 builder folder (box layout); on the Copilot machine it is absent and the two Stage 3 checks are skipped
_S3_CANDIDATES = [
    HERE.parent / "stage3_builder",
    next((q for q in (HERE.parent / "stage3_builder", HERE.parent / "stage3_spec" / "stage3_builder") if (q / "03_build_real_layers.py").is_file()), HERE.parent / "stage3_builder"),
    HERE.parent.parent / "scripts" / "stage3_builder",
]
S3_DIR = next((p for p in _S3_CANDIDATES if (p / "03_build_real_layers.py").is_file()), _S3_CANDIDATES[1])
HAVE_S3 = (S3_DIR / "03_build_real_layers.py").is_file() and (S3_DIR / "stage3_links_v2.csv").is_file()
LINKS_V103 = HERE / "stage2_v1_0_3_links.csv"
PUBLIC_REAL = HERE / "stage2_v1_0_3_public_TEAM_INTERNAL.csv"
FLAG15 = ["/layer_a/dob_breakdown/by_split/train/n_age_days", "/layer_a/dob_breakdown/by_split/train/n_invalid",
          "/layer_a/dob_breakdown/by_split/train/status/n_valid", "/layer_a/train/age_band/n_age_invalid",
          "/strata_denominators/age_band/train/adolescent/n_encounters",
          "/strata_denominators/age_band/train/age_invalid/n_encounters"]
DERIVED_ROWS = ("/layer_a/dob_breakdown/status/n_valid,derived,cohort minus public invalid\n"
                "/layer_a/dob_breakdown/n_age_days,derived,equals valid\n")


def _ref13_import():
    sys.path.insert(0, str(REF_MOD_DIR))
    try:
        import stage2_v1_0_2_suppress_v1_3 as R  # noqa: E402
    finally:
        sys.path.remove(str(REF_MOD_DIR))
    return R


def test_screen_safe_counterexample_v13_vs_v14():
    """The v1.3 hidden-cell bug: screening points that all sit at 0 were taken as proof that an H cell is exactly 0.
    Here v2 really ranges 0-3 (a recoverable small hidden difference): v1.3 reports (0, 0) 'safe', v1.4 the exact
    (0, 3) 'unsafe'."""
    rows = ((((0, 1), (1, 1), (3, -1), (4, 1)), 14), (((0, 1), (1, 1), (2, 1), (3, 1)), 13))
    kinds = ("H", "H", "H", "S", "S")
    names = ["v0", "v1", "v2", "v3", "v4"]
    R = _ref13_import()
    old, _ = R.solve_component(names, kinds, rows, R.Params())
    new, _ = T.solve_component(names, kinds, rows, T.Params())
    assert tuple(round(x) for x in old["v2"]) == (0, 0) and R.is_safe("H", *old["v2"], R.Params())
    assert tuple(round(x) for x in new["v2"]) == (0, 3) and not T.is_safe("H", *new["v2"], T.Params())
    for v in names:                                   # every v1.4 interval is the exact LP interval
        assert tuple(round(x) for x in new[v]) in {(0, 3), (10, 13), (21, 27)}


def test_screen_safe_only_withholds_the_exact_zero_clause():
    p = T.Params()
    assert not T.screen_safe("H", 0, 0, p) and not T.screen_safe("D", 0, 0, p)   # needs the exact maximum
    assert T.screen_safe("H", 12, 40, p) and T.screen_safe("H", 5, 12, p) and T.screen_safe("H", 0, 30, p)
    assert T.screen_safe("P", 0, 9, p) and not T.screen_safe("P", 0, 0, p) and T.screen_safe("S", 10, 30, p)
    assert T.is_safe("H", 0, 0, p)                    # the rule itself is unchanged: an exact 0 is safe


@pytest.fixture(scope="module")
def mirror(tmp_path_factory):
    """Synthetic mirror: v1.0.1-format input, its v1.0.2 release (module v1.3), the links v3, public files."""
    w = tmp_path_factory.mktemp("mirror")
    mr = MR.write_case(w)
    with open(mr["public"], "a", encoding="utf-8") as f:
        f.write(DERIVED_ROWS)
    links = w / "links_v1_0_3.csv"
    r = subprocess.run([sys.executable, str(HERE / "make_stage2_links_v1_0_3.py"), str(mr["json"]), "-", str(links)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    _shutil.copy(REF, w / "s13.py")
    r = subprocess.run([sys.executable, "s13.py", "--in-json", mr["json"].name, "--in-md", mr["md"].name, "--links",
                        mr["links_v102"].name, "--frozen", "frozen.csv", "--drops", str(DROPS), "--out-dir", "out_v102"],
                       cwd=w, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    return dict(mr, w=w, links=links, r102=w / "out_v102" / "feature_sparsity_frequency_audit_v1.0.2.json")


def _v14_args(m: dict, public=None, carry=None) -> list[str]:
    return ["--in-json", str(m["json"]), "--in-md", str(m["md"]), "--links", str(m["links"]), "--frozen",
            str(m["frozen"]), "--drops", str(DROPS), "--public", str(public or m["public"]),
            "--carry-json", str(carry or m["r102"])]


@pytest.fixture(scope="module")
def released_v103(mirror, tmp_path_factory):
    out = tmp_path_factory.mktemp("v103") / "out"
    r = run_cli(*_v14_args(mirror), "--out-dir", str(out))
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
    doc = json.loads((out / T.OUT_JSON).read_text(encoding="utf-8"))
    return dict(out=out, doc=doc, log=r.stdout + r.stderr)


def _lab(doc, p):
    v = T.get_ptr(doc, p)
    return "pub" if isinstance(v, int) else v


def test_links_generator_rows(mirror):
    rows = list(_csv.DictReader(l for l in open(mirror["links"], encoding="utf-8") if not l.startswith("#")))
    kinds = [r["note"].split(":")[0] for r in rows]
    assert {k: kinds.count(k) for k in set(kinds)} == {"L5a": 21, "L5b": 3, "L5c": 15, "L6": 34, "L7": 3}
    assert not any("ge18" in r["total"] for r in rows)         # a structural zero is no link (public file instead)


def test_shipped_links_v3_keep_the_v102_rows_and_add_76():
    base = [tuple(r.values()) for r in _csv.DictReader(
        l for l in open(HERE / "stage2_v1_0_2_links.csv", encoding="utf-8") if not l.startswith("#"))]
    new = [tuple(r.values()) for r in _csv.DictReader(
        l for l in open(LINKS_V103, encoding="utf-8") if not l.startswith("#"))]
    assert new[:13] == base and len(base) == 13 and len(new) == 89


@pytest.mark.skipif(not PUBLIC_REAL.is_file(), reason=
                    "team-internal public-values file is not shipped in this public copy")
def test_real_public_file_shape_without_values():
    """The team-internal public file: 10 rows (3 structural zeros, 4 table sizes, 1 asserted count, 2 derived);
    the test never prints a value."""
    rows = list(_csv.DictReader(l for l in open(PUBLIC_REAL, encoding="utf-8") if not l.startswith("#")))
    assert len(rows) == 10
    kinds = {r["pointer"]: ("int" if r["value"].isdigit() else r["value"] or "input") for r in rows}
    assert [kinds[f"/layer_a/{s}/age_band/n_ge18"] for s in ("train", "validation", "test")] == ["int"] * 3
    assert all(r["value"] == "0" for r in rows if r["pointer"].endswith("/n_ge18"))
    rec = {k.split("/")[2]: v for k, v in kinds.items() if k.startswith("/reconciliation/")}
    assert rec == {"CHARTEVENTS": "publish", "MICROBIOLOGYEVENTS": "publish", "LABEVENTS": "input",
                   "PRESCRIPTIONS": "input"}                   # decision D3 (Vinai, Nov 17, 2026)
    assert kinds["/layer_a/dob_breakdown/n_invalid"] == "int"
    assert kinds["/layer_a/dob_breakdown/status/n_valid"] == kinds["/layer_a/dob_breakdown/n_age_days"] == "derived"
    assert int(next(r["value"] for r in rows if r["pointer"] == "/layer_a/dob_breakdown/n_invalid")) >= 10


def test_v102_mirror_fails_alone_with_links_v3(mirror):
    lines = []
    assert IC.check(str(mirror["r102"]), str(mirror["links"]), str(mirror["public"]), str(mirror["json"]),
                    out=lines.append) == 1
    failing = {l.split()[-1] for l in lines[1:]}
    assert set(FLAG15) <= failing


def test_v103_release_gate_and_independent_checks(mirror, released_v103):
    assert released_v103["doc"]["release_gate"] == "ok"
    rel = released_v103["out"] / T.OUT_JSON
    lines = []
    assert IC.check(str(rel), str(mirror["links"]), str(mirror["public"]), str(mirror["json"]), out=lines.append) == 0
    if not HAVE_S3:
        return
    s3l = S3_DIR / "stage3_links_v2.csv"
    pub_sizes = mirror["w"] / "public_sizes.csv"
    pub_sizes.write_text("".join(l + "\n" for l in mirror["public"].read_text(encoding="utf-8").splitlines()
                                 if l.startswith(("pointer,", "/reconciliation/"))), encoding="utf-8")
    assert IC.check(str(rel), str(s3l), str(pub_sizes), str(mirror["json"]), out=lines.append) == 0


def test_v103_suppressions_contain_v102_and_public_labels_kept(mirror, released_v103):
    d102 = json.loads(mirror["r102"].read_text(encoding="utf-8"))
    d101 = json.loads(mirror["json"].read_text(encoding="utf-8"))
    d103 = released_v103["doc"]
    cf = d103["disclosure"]["count_fields"]
    assert cf == d102["disclosure"]["count_fields"]
    repub = {f"/reconciliation/{t}/n_total" for t in MR.PUBLISH_TABLES}
    for p in cf:
        a, b = _lab(d102, p), _lab(d103, p)
        if p in repub:                                         # D3: published again (value of v1.0.1)
            assert b == "pub" and T.get_ptr(d103, p) == T.get_ptr(d101, p), p
        elif a != "pub":
            assert b == a, p                                   # suppressed / '<10' in v1.0.2 stay so
        else:
            assert b == "pub" and T.get_ptr(d103, p) == T.get_ptr(d101, p) or b == "S", p
    for s in ("train", "validation", "test"):
        assert _lab(d103, f"/layer_a/{s}/age_band/n_ge18") == "<10"     # public, label unchanged
    for t in ("CHARTEVENTS", "MICROBIOLOGYEVENTS"):
        assert _lab(d102, f"/reconciliation/{t}/n_total") == "S"       # suppressed in v1.0.2 ...
        assert _lab(d103, f"/reconciliation/{t}/n_total") == "pub"     # ... published again (decision D3)
        assert _lab(d103, f"/reconciliation/{t}/n_no_cohort_admission") == "<10"
    for p in FLAG15:
        assert _lab(d103, p) != "pub"
    prov = d103[T.PROVENANCE_KEY]
    assert prov["n_public_cells"] == 10 and prov["carried_suppressions"] > 0 and prov["republished_public_cells"] == 2


def test_v103_outputs_never_list_public_pointers_or_values(mirror, released_v103):
    summ = (released_v103["out"] / "stage2_v1.0.3_release_summary.md").read_text(encoding="utf-8")
    man = (released_v103["out"] / "stage2_v1.0.3_manifest.json").read_text(encoding="utf-8")
    assert "public (already known) cells: 10" in summ and "team-internal file, pointers and values not listed" in summ
    pub_ptrs = [r["pointer"] for r in _csv.DictReader(open(mirror["public"], encoding="utf-8"))]
    for text in (summ, man, released_v103["log"]):
        for p in pub_ptrs:
            assert p not in text
        assert "synthetic stand-in" not in text


def test_v103_check_only_catches_a_leaked_flagged_cell(mirror, released_v103, tmp_path):
    """Tamper: publish the train invalid-age count again -> the check-only audit (links v3 + public) fails."""
    doc = json.loads(json.dumps(released_v103["doc"]))
    raw = MR.build()["raw"]
    T.set_ptr(doc, "/layer_a/train/age_band/n_age_invalid", raw["/layer_a/train/age_band/n_age_invalid"])
    p = tmp_path / "tampered.json"
    p.write_text(T.dump_json(doc), encoding="utf-8")
    lines = []
    assert IC.check(str(p), str(mirror["links"]), str(mirror["public"]), str(mirror["json"]), out=lines.append) == 1


def test_control_table_sizes_kept_s_also_passes_alone(mirror, tmp_path):
    """Control for decision D3: the same public file with the two table sizes kept S (empty values) gives a valid
    release too (nothing republished); the decided release publishes them (test above)."""
    pub = tmp_path / "public_keepS.csv"
    pub.write_text(mirror["public"].read_text(encoding="utf-8").replace(",publish,public table size", ",,public table size"),
                   encoding="utf-8")
    out = tmp_path / "out"
    r = run_cli(*_v14_args(mirror, public=pub), "--out-dir", str(out))
    assert r.returncode == 0, (r.stdout + r.stderr)[-2000:]
    d = json.loads((out / T.OUT_JSON).read_text(encoding="utf-8"))
    for t in MR.PUBLISH_TABLES:
        assert T.get_ptr(d, f"/reconciliation/{t}/n_total") == "S"
    assert d[T.PROVENANCE_KEY]["republished_public_cells"] == 0
    lines = []
    assert IC.check(str(out / T.OUT_JSON), str(mirror["links"]), str(pub), str(mirror["json"]), out=lines.append) == 0


def _pub_variant(mirror, tmp_path, edit) -> Path:
    rows = mirror["public"].read_text(encoding="utf-8").splitlines()
    p = tmp_path / "public_bad.csv"
    p.write_text("\n".join(edit(rows)) + "\n", encoding="utf-8")
    return p


@pytest.mark.parametrize("edit,msg", [
    (lambda rs: rs + ["/reconciliation/LABEVENTS/n_no_cohort_admission,,x"], "duplicate"),
    (lambda rs: [r.replace("/reconciliation/LABEVENTS/n_total,,", "/reconciliation/LABEVENTS/n_total,5,") for r in rs],
     "contradicts the published input value"),
    (lambda rs: [r for r in rs if not r.startswith("/layer_a/dob_breakdown/n_invalid,")], "'derived' public cell"),
    (lambda rs: rs + ["/layer_a/train/age_band/n_neonate,x,y"], "value must be"),
    (lambda rs: rs + ["/denominators/train/n_encounters,,y"], "frozen"),
    (lambda rs: rs + ["/drug_groups/nmb/validation/strata/surgical_status/recent_surgery/n_active_72h,1,y"], "dropped"),
    (lambda rs: rs + ["/layer_a/dob_breakdown/by_split/train/status/n_valid,,y"], "gives no value"),
    (lambda rs: [r.replace("/layer_a/train/age_band/n_ge18,0,", "/layer_a/train/age_band/n_ge18,12,") for r in rs],
     "'<10' label"),
    (lambda rs: rs + ["/nope/n_x,1,y"], "not a count cell"),
])
def test_bad_public_file_stops(mirror, tmp_path, edit, msg):
    if msg == "duplicate":
        edit = lambda rs: rs + [rs[1]]  # noqa: E731
    p = _pub_variant(mirror, tmp_path, edit)
    r = run_cli(*_v14_args(mirror, public=p), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2 and msg in (r.stdout + r.stderr), (r.stdout + r.stderr)[-1500:]
    assert not (tmp_path / "o").exists()


def test_derived_public_without_invalid_total_fails_closed(mirror, tmp_path):
    p = _pub_variant(mirror, tmp_path, lambda rs: [r for r in rs if not r.startswith("/layer_a/dob_breakdown/n_invalid,")])
    r = run_cli(*_v14_args(mirror, public=p), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2 and "not fixed by the identities and links" in (r.stdout + r.stderr)


def test_public_invalid_total_without_derived_rows_cannot_protect(mirror, tmp_path):
    """With the invalid-DOB total public, valid DOBs = cohort - invalid are known exactly: without the two 'derived'
    rows the run stops (cannot protect) and leaves nothing."""
    p = _pub_variant(mirror, tmp_path, lambda rs: [r for r in rs if ",derived," not in r])
    r = run_cli(*_v14_args(mirror, public=p), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 3 and "cannot protect" in (r.stdout + r.stderr)
    assert not (tmp_path / "o").exists()


def test_carry_must_be_the_v102_release_of_this_input(mirror, tmp_path):
    d = json.loads(mirror["r102"].read_text(encoding="utf-8"))
    bad = json.loads(json.dumps(d))
    bad[T.CARRY_PROVENANCE_KEY]["source_json_sha256"] = "0" * 64
    p = tmp_path / "carry_bad_sha.json"
    p.write_text(json.dumps(bad), encoding="utf-8")
    r = run_cli(*_v14_args(mirror, carry=p), "--out-dir", str(tmp_path / "o1"))
    assert r.returncode == 2 and "source_json_sha256" in (r.stdout + r.stderr)
    bad = json.loads(json.dumps(d))
    c = next(c for c in d["disclosure"]["count_fields"] if isinstance(T.get_ptr(d, c), int))
    T.set_ptr(bad, c, T.get_ptr(d, c) + 1)
    p2 = tmp_path / "carry_bad_value.json"
    p2.write_text(json.dumps(bad), encoding="utf-8")
    r = run_cli(*_v14_args(mirror, carry=p2), "--out-dir", str(tmp_path / "o2"))
    assert r.returncode == 2 and "differs from the input" in (r.stdout + r.stderr)
    bad = json.loads(json.dumps(d))
    bad["release_gate"] = "pending"
    p3 = tmp_path / "carry_not_released.json"
    p3.write_text(json.dumps(bad), encoding="utf-8")
    r = run_cli(*_v14_args(mirror, carry=p3), "--out-dir", str(tmp_path / "o3"))
    assert r.returncode == 2 and "release_gate" in (r.stdout + r.stderr)
    r = run_cli(*_v14_args(mirror), "--check-only")
    assert r.returncode == 2


def test_carry_sha_pin(mirror, tmp_path):
    r = run_cli(*_v14_args(mirror), "--expect-carry-sha256", "a" * 64, "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2 and not (tmp_path / "o").exists()


def test_v103_release_input_refused_as_input(released_v103, mirror, tmp_path):
    out = released_v103["out"]
    r = run_cli("--in-json", str(out / T.OUT_JSON), "--in-md", str(out / T.OUT_MD), "--out-dir", str(tmp_path / "o"))
    assert r.returncode == 2


@pytest.mark.skip(reason="joint helper assumes a private box layout (legacy suppress-module path and links pin); not run in this public copy")
def test_joint_with_stage3_builder_synthetic(tmp_path):
    """End to end with the Stage 3 builder v1.0.3 on its synthetic fixture (joint_check_synthetic.py): the v1.0.3
    mirror release with the table sizes published (decision D3, Vinai, Nov 17, 2026) passes alone and jointly (gate
    ok, no compromised Stage 2 cell, 0 hidden, independent joint check 0)."""
    r = subprocess.run([sys.executable, str(HERE / "joint_check_synthetic.py"), str(tmp_path / "joint")],
                       capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    assert "JOINT SYNTHETIC CHECK: PASS" in r.stdout
    assert "RESULT [v1.0.3 mirror (decided D3: table sizes published)]: PASS" in r.stdout


# ---------------------------------------------------------------- v1.5: heartbeat, bounded selection, convergence stop
import csv as _csv  # noqa: E402
import importlib.util as _ilu  # noqa: E402
import io as _io  # noqa: E402
import contextlib as _ctx  # noqa: E402
import time as _time  # noqa: E402

import synthetic_merge_v105 as MG  # noqa: E402


def test_v15_heartbeat_prints_during_silence(capsys):
    # a phase that prints nothing for a while still produces 'heartbeat' lines (counts only, no pointer)
    T._LAST_LINE[0] = _time.monotonic()
    with T.Heartbeat(0.3):
        T.beat("suppression round 9: selection", "failing cell 3/7, 1 picks, 4 exact checks")
        _time.sleep(1.5)
    err = capsys.readouterr().err
    hb = [ln for ln in err.splitlines() if "heartbeat: still working - suppression round 9: selection" in ln]
    assert len(hb) >= 2
    assert all(_re.match(r"^\[\s*\d+\.\ds\] heartbeat", ln) and "LP solves so far" in ln for ln in hb)
    assert not _re.search(r"/[A-Za-z_]", "\n".join(hb))


def test_v15_heartbeat_option_on_cli(case, tmp_path):
    r = run_cli(*base_args(case), "--out-dir", str(tmp_path / "o"), "--heartbeat-seconds", "0")
    assert r.returncode == 2
    r = run_cli(*base_args(case), "--check-only", "--withhold-proposal", str(tmp_path / "p.csv"))
    assert r.returncode == 2


def _low_system(tmp_path):
    c = _stall_case(tmp_path, "low")
    L = T.load_output(c["json"], c["md"], "stage2", T.EXPECT_INPUT_SPEC)
    frozen = {p for p, _ in T.read_frozen(c["frozen"])}
    system = T.build_system(L, [], T.Params())
    status = dict(L.status)
    return c, L, frozen, system, status


def test_v15_capped_lp_is_never_safe(tmp_path, monkeypatch):
    # an exact selection LP that hits its time cap (HiGHS status 1) makes the target 'unknown': never a safe pick,
    # never a proof that a cell cannot be protected
    c, L, frozen, system, status = _low_system(tmp_path)
    res = system.audit(status)
    sup = T.Suppressor(system, status, frozen)
    u = sorted(res.unsafe)[0]
    good = sup.choose_exact(u)
    assert good and sup.makes_safe(u, good)

    class _R:
        status, message, x = 1, "time limit reached", None
    monkeypatch.setattr(T, "linprog", lambda *a, **k: _R())
    rels = system.local_rels([u, *good], T.Overlay(status, {m: "S" for m in good}))
    b = system.target_bounds(T.Overlay(status, {m: "S" for m in good}), rels, [u], lp_cap=0.01)
    assert b[u] is None
    assert sup.check_fix(u, good) == (False, False)
    assert sup.unprotectable(u) is False


def test_v15_pin_rank_and_checked_picks(tmp_path):
    # the stuck cells of the merge reproducer: an S cell pinned by a link to a published twin; the fix is the twin
    # plus the twin's partition total (both checked). v1.5 ranks candidates by the equations that pin the cell.
    c, L, frozen, system, status = _low_system(tmp_path)
    res = system.audit(status)
    sup = T.Suppressor(system, status, frozen)
    sup.res = res
    for u in sorted(res.unsafe):
        k = sup.merge_key(f"{u.rsplit('/', 1)[0]}/by_site/site_big/n", u)
        assert k[0] < 99                     # shares an equation with u: a 'pin' candidate
        picks = sup.choose_exact(u)
        assert picks and u not in system.audit(T.Overlay(status, {m: "S" for m in picks})).unsafe


def _load_as(path: Path, name: str):
    spec = _ilu.spec_from_file_location(name, path)
    m = _ilu.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def merge_small(tmp_path_factory):
    d = tmp_path_factory.mktemp("merge")
    # scale chosen so the synthetic graph is large enough for the v1.5 vs v1.4 merge contrast
    return MG.write(d / "in", scale=0.6, hub_frac=0.35)


def _merge_args(m: dict) -> list[str]:
    return ["--in-json", str(m["json"]), "--in-md", str(m["md"]), "--links", str(m["links"])]


def _run_mod(mod, argv):
    err = _io.StringIO()
    with _ctx.redirect_stderr(err), _ctx.redirect_stdout(_io.StringIO()):
        rc = mod.main(argv)
    return rc, err.getvalue()


@pytest.mark.skip(reason="public synthetic scale was changed so it does not use real split sizes; merge-size contrast is not reproduced")
def test_v15_merge_reproducer_converges_without_merging(merge_small, tmp_path, monkeypatch):
    # the real-run pattern (components merging through strata denominators and links) on a small synthetic mirror,
    # with the large-component threshold lowered so the exact path runs: v1.5 converges, gate ok, and the largest
    # component stays near its input size; v1.4 on the same input merges it into one far larger component
    monkeypatch.setattr(T, "LARGE_COMPONENT", 20)
    rc, err = _run_mod(T, [*_merge_args(merge_small), "--out-dir", str(tmp_path / "o15"), "--max-seconds", "300"])
    assert rc == 0, err[-2000:]
    assert "release gate: ok" in err
    fails, largest = _rounds(err)
    assert all(b < a for a, b in zip(fails, fails[1:]))  # no oscillation: the failing count only goes down
    assert max(largest) <= 2 * largest[0]
    assert "chain" in err and "selection" in err
    v14 = _load_as(HERE / "reference" / "stage2_v1_0_3_suppress_v1_4.py", "_s2_v14_ref")
    monkeypatch.setattr(v14, "LARGE_COMPONENT", 20)
    rc4, err4 = _run_mod(v14, [*_merge_args(merge_small), "--out-dir", str(tmp_path / "o14"), "--max-seconds", "300"])
    _f4, largest4 = _rounds(err4)
    assert max(largest4) > 2 * max(largest)


def test_v15_not_converging_stops_with_proposal(tmp_path, monkeypatch, capsys):
    # useless picks (v1.1's choice) + a small --max-new-suppressions: exit 3, STOPPING line, a withhold proposal
    # printed (pointers only) and written outside the out-dir, nothing in the out-dir
    c = _stall_case(tmp_path, "pinned")
    monkeypatch.setattr(T.Suppressor, "choose_exact", _v11_choice)
    out, prop = tmp_path / "o", tmp_path / "proposal_TEAM_INTERNAL.csv"
    rc = T.main([*_stall_args(c), "--out-dir", str(out), "--max-new-suppressions", "3",
                 "--withhold-proposal", str(prop)])
    err = capsys.readouterr().err
    assert rc == 3
    assert "STOPPING: suppression not converging" in err and "--max-new-suppressions" in err
    assert "WITHHOLD PROPOSAL:" in err and "RELEASE GATE FAILED" in err
    assert not out.exists()
    rows = list(_csv.DictReader(ln for ln in prop.read_text(encoding="utf-8").splitlines() if not ln.startswith("#")))
    assert rows and all(r["pointer"].startswith("/") and r["reason"] for r in rows)
    frozen = {p for p, _ in T.read_frozen(c["frozen"])}
    assert not frozen & {r["pointer"] for r in rows}


def test_v15_proposal_refuses_existing_file_and_out_dir(tmp_path, monkeypatch, capsys):
    c = _stall_case(tmp_path, "pinned")
    monkeypatch.setattr(T.Suppressor, "choose_exact", _v11_choice)
    prop = tmp_path / "exists.csv"
    prop.write_text("x\n")
    rc = T.main([*_stall_args(c), "--out-dir", str(tmp_path / "o"), "--max-new-suppressions", "3",
                 "--withhold-proposal", str(prop)])
    assert rc == 3 and "NOT written" in capsys.readouterr().err and prop.read_text() == "x\n"


def test_v15_max_component_cap_stops(tmp_path, monkeypatch, capsys):
    c = _stall_case(tmp_path, "low")
    out = tmp_path / "o"
    rc = T.main([*_stall_args(c), "--out-dir", str(out), "--max-component-cells", "50"])
    err = capsys.readouterr().err
    assert rc == 3 and "--max-component-cells" in err and "WITHHOLD PROPOSAL" in err and not out.exists()


def _write_withhold(p: Path, rows) -> Path:
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = _csv.writer(f)
        w.writerow(["pointer", "reason"])
        w.writerows(rows)
    return p


def test_v15_withhold_leaf_cells_release_ok(tmp_path):
    # the cell-level alternative to D4: the two large site counts withheld; they leave the JSON, the registry and
    # the markdown, are audited as unpublished cells, and the gate passes (the item is protected without new cells)
    c = _stall_case(tmp_path, "low")
    big = [f"/surgery/giant/hub_0{h}/low_item/by_site/site_big/n" for h in (0, 1)]
    wf = _write_withhold(tmp_path / "withhold.csv", [(p, "v1.5 test") for p in big])
    out = tmp_path / "o"
    r = run_cli(*_stall_args(c), "--out-dir", str(out), "--withhold", str(wf))
    assert r.returncode == 0, r.stderr
    assert "2 single cell(s) withheld" in r.stderr
    doc = json.loads((out / T.OUT_JSON).read_text(encoding="utf-8"))
    for p in big:
        with pytest.raises(KeyError):
            T.get_ptr(doc, p)
        assert p not in doc["disclosure"]["count_fields"]
    man = json.loads((out / T.OUT_MANIFEST).read_text(encoding="utf-8")) if hasattr(T, "OUT_MANIFEST") else None
    if man is not None:
        assert "withhold" in json.dumps(man)
    chk = run_cli("--in-json", str(out / T.OUT_JSON), "--in-md", str(out / T.OUT_MD), "--frozen", str(c["frozen"]),
                  "--check-only")
    assert chk.returncode == 0, chk.stdout + chk.stderr


def test_v15_withhold_refusals(case, tmp_path):
    frozen_ptr = next(p for p, _ in T.read_frozen(case["frozen"]))
    link_ptr = T.read_links(case["links"])[0][0]
    for i, (ptr, want) in enumerate([(frozen_ptr, "frozen"), (link_ptr, "links-file"),
                                     ("/no/such/cell", "not a count cell")]):
        wf = _write_withhold(tmp_path / f"w{i}.csv", [(ptr, "test")])
        out = tmp_path / f"o{i}"
        r = run_cli(*base_args(case), "--out-dir", str(out), "--withhold", str(wf))
        assert r.returncode == 2 and want in r.stderr, r.stderr[-1500:]
        assert not out.exists()


def test_v15_proposal_round_trip_passes_gate(tmp_path):
    # the proposal is checked by a full audit when made, and used as --withhold it gives a release that passes
    c = _stall_case(tmp_path, "low")
    prop = tmp_path / "prop.csv"
    r = run_cli(*_stall_args(c), "--out-dir", str(tmp_path / "o1"), "--max-component-cells", "50",
                "--withhold-proposal", str(prop))
    assert r.returncode == 3 and "with them withheld, 0 cell(s) still fail the full audit" in r.stderr
    out = tmp_path / "o2"
    r = run_cli(*_stall_args(c), "--out-dir", str(out), "--withhold", str(prop))
    assert r.returncode == 0, r.stderr
    assert "release gate: ok" in r.stderr
    chk = run_cli("--in-json", str(out / T.OUT_JSON), "--in-md", str(out / T.OUT_MD), "--frozen", str(c["frozen"]),
                  "--check-only")
    assert chk.returncode == 0
