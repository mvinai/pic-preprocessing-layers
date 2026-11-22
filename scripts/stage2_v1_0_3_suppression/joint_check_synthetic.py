#!/usr/bin/env python3
"""joint_check_synthetic.py: end-to-end synthetic check of Stage 2 v1.0.3 against the Stage 3 builder v1.0.3
(SYNTHETIC ONLY: the Stage 3 PIC-format fixture and a Stage 2 mirror built from it; never real data).

  1. Stage 3 fixture project (stage3_fixture_pic.make_project); its Stage 3 audit values are captured just before the
     builder's release gate (as test_03_build_real_layers.add_stage2_release does).
  2. A v1.0.1-format Stage 2 mirror consistent with those values under stage3_links_v2.csv
     (synthetic_s2_v103_mirror.build(values=...)).
  3. v1.0.2 mirror release: the pinned module v1.3 (reference copy), v1.0.2 configuration (no links, frozen, drops).
  4. v1.0.3 mirror release: module v1.4 with --carry-json (step 3), the generated links v3, --public, frozen, drops.
  5. Each release: (a) Stage 2 alone, independent check (stage2_independent_check.py) with the links v3 and public
     values; (b) jointly: the Stage 3 builder v1.0.3 with that release (pin overridden as in its tests) -> RELEASE
     GATE, the compromised Stage 2 list (flagged_for_stage2_v1_0_3), and the independent joint check
     stage3_joint_check.py (links v2).
Prints one line per result; exit 0 when the decided v1.0.3 (D3: table sizes published) passes everything (alone: 0 failing; joint: gate ok, nothing
compromised, joint check 0)."""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
S3 = next((q for q in (HERE.parent / "stage3_builder", HERE.parent / "stage3_spec" / "stage3_builder") if (q / "03_build_real_layers.py").is_file()), HERE.parent / "stage3_builder")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(S3))
import stage3_fixture_pic as FX  # noqa: E402
import stage3_joint_check as JC  # noqa: E402
import synthetic_s2_v103_mirror as MR  # noqa: E402
import stage2_independent_check as IC  # noqa: E402

_spec = importlib.util.spec_from_file_location("stage3_builder_03", S3 / "03_build_real_layers.py")
B = importlib.util.module_from_spec(_spec)
sys.modules["stage3_builder_03"] = B
_spec.loader.exec_module(B)
PY = sys.executable


def make_cfg(fx, out, **over):
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
              v17_module_path=S3 / "audit_surgery_complexity_duration.py",
              suppress_module_path=S3 / "stage2_v1_0_3_suppress.py", links_path=S3 / "stage3_links_v2.csv",
              out_dir=out, expect_d_icd_sha256=fx["d_icd"], expect_svs_sha256=fx["svs"], pins=pins,
              raw_pins=dict(fx["raw_pins"]), expect=ex, links_sha256=B.LINKS_SHA256,
              drops_path=S3 / "stage3_drops_v1.csv", drops_sha256=B.DROPS_SHA256,
              stage2_release_path=fx.get("s2"), expect_stage2_sha256=fx.get("s2_sha"),
              stage2_pin=fx.get("s2_sha") or B.STAGE2_RELEASE_SHA256)
    kw.update(over)
    return B.Config(**kw)


class _Captured(Exception):
    pass


def build3(cfg):
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        res = B.build(cfg)
    return res, err.getvalue()


def capture_s3(fx, folder):
    folder.mkdir(parents=True, exist_ok=True)
    ph = folder / "stage2_placeholder.json"
    fx = dict(fx, s2=ph, s2_sha=FX.make_stage2_release(ph))
    cap = {}
    real = B.release_gate

    def grab(SM, g, links_path, base, pending, xr=None):
        cap.update(g.values())
        raise _Captured()
    B.release_gate = grab
    try:
        try:
            build3(make_cfg(fx, folder / "capture_out"))
        except _Captured:
            pass
    finally:
        B.release_gate = real
    if not cap:
        raise SystemExit("could not capture the Stage 3 audit values")
    return cap


def sh(cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def joint(fx, s2, out, label):
    sha = B.sha256_file(s2)
    cfg = make_cfg(dict(fx, s2=s2, s2_sha=sha), out)
    res, err = build3(cfg)
    man = json.loads((out / "stage3_manifest.json").read_text(encoding="utf-8"))
    xr = man["cross_release"]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        jc = JC.main(str(out / "stage3_audit.json"), str(s2), str(S3 / "stage3_links_v2.csv"))
    print(f"[{label}] joint with Stage 3 builder: RELEASE GATE {man['release_gate'].upper()}; compromised Stage 2 "
          f"cells {len(xr['flagged_for_stage2_v1_0_3'])} ({xr['flagged_hidden_stage2_differences']}); "
          f"independent joint check exit {jc}: {buf.getvalue().splitlines()[0]}")
    return man, xr, jc, buf.getvalue()


def main(work: Path) -> int:
    work = work.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    fx = FX.make_project(work / "proj")
    s3 = capture_s3(fx, work / "cap")
    print(f"Stage 3 fixture built; {len(s3)} Stage 3 audit values captured (synthetic)")
    case = work / "s2"
    mr = MR.write_case(case, values=s3, splits=fx["expect"]["splits"])
    rc, log = sh([PY, str(HERE / "make_stage2_links_v1_0_3.py"), str(mr["json"]), "-", str(case / "links_v1_0_3.csv")],
                 HERE)
    assert rc == 0, log
    shutil.copy(HERE / "reference" / "stage2_v1_0_2_suppress_v1_3.py", case / "s13.py")
    drops = HERE / "stage2_v1_0_3_drops.csv"
    rc, log = sh([PY, "s13.py", "--in-json", mr["json"].name, "--in-md", mr["md"].name, "--links",
                  mr["links_v102"].name, "--frozen", "frozen.csv", "--drops", str(drops), "--out-dir", "out_v102"], case)
    print(f"[v1.0.2 mirror] module v1.3 run exit {rc}")
    assert rc == 0, log[-3000:]
    r102 = case / "out_v102" / "feature_sparsity_frequency_audit_v1.0.2.json"
    n_inv = mr["raw"]["/layer_a/dob_breakdown/n_invalid"]
    if n_inv >= 10:
        with open(mr["public"], "a", encoding="utf-8") as f:
            f.write("/layer_a/dob_breakdown/status/n_valid,derived,cohort size (frozen) minus the public invalid-DOB count\n")
            f.write("/layer_a/dob_breakdown/n_age_days,derived,equals the valid-DOB count\n")
    else:   # the fixture's invalid-DOB total is a '<10' cell, not a count quoted elsewhere: not public here
        rows = [l for l in mr["public"].read_text(encoding="utf-8").splitlines()
                if not l.startswith("/layer_a/dob_breakdown/n_invalid,")]
        mr["public"].write_text("\n".join(rows) + "\n", encoding="utf-8")
        print("[v1.0.3 mirror] fixture invalid-DOB total < 10: not treated as public (the mirror case covers that path)")
    pub_keep = case / "public_keepS_control.csv"      # control only: D3 decided 'publish' (Vinai, Nov 17, 2026)
    pub_keep.write_text(mr["public"].read_text(encoding="utf-8").replace(",publish,public table size",
                                                                         ",,public table size"), encoding="utf-8")
    rels = [("v1.0.2 mirror", r102, None)]
    for label, pub, od in (("v1.0.3 mirror (control: table sizes kept S)", pub_keep, "out_v103_keepS"),
                           ("v1.0.3 mirror (decided D3: table sizes published)", mr["public"], "out_v103")):
        rc, log = sh([PY, str(HERE / "stage2_v1_0_3_suppress.py"), "--in-json", mr["json"].name, "--in-md",
                      mr["md"].name, "--links", "links_v1_0_3.csv", "--frozen", "frozen.csv", "--drops", str(drops),
                      "--public", pub.name, "--carry-json", str(r102), "--out-dir", od], case)
        gate = [l for l in log.splitlines() if "RELEASE GATE" in l or "release gate:" in l]
        print(f"[{label}] module v1.4 run exit {rc}; {gate[-1].strip() if gate else ''}")
        assert rc == 0, log[-3000:]
        rels.append((label, case / od / "feature_sparsity_frequency_audit_v1.0.3.json", pub))
    ok = {}
    for k, (label, rel, pub) in enumerate(rels):
        lines = []
        rc_a = IC.check(str(rel), str(case / "links_v1_0_3.csv"), str(pub or mr["public"]), str(mr["json"]),
                        out=lines.append)
        print(f"[{label}] Stage 2 alone (links v3 + public): exit {rc_a}: {lines[0]}")
        lines2 = []
        rc_b = IC.check(str(rel), str(S3 / "stage3_links_v2.csv"), None, str(mr["json"]), out=lines2.append)
        print(f"[{label}] Stage 2 alone (Stage 3 links v2 internal rows): exit {rc_b}: {lines2[0]}")
        man, xr, jc, _ = joint(fx, rel, work / f"s3_out_{k}", label)
        if xr["flagged_for_stage2_v1_0_3"]:
            print(f"  [{label}] compromised (pointers): " + "; ".join(xr["flagged_for_stage2_v1_0_3"]))
        ok[label] = (rc_a == 0 and rc_b == 0 and man["release_gate"] == "ok" and not xr["flagged_for_stage2_v1_0_3"]
                     and xr["flagged_hidden_stage2_differences"].startswith("0 ") and jc == 0)
    for label, v in ok.items():
        print(f"RESULT [{label}]: {'PASS' if v else 'NOT CLEAN'}")
    ok = ok[rels[-1][0]]
    print("JOINT SYNTHETIC CHECK: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="s2v103_joint_"))))
