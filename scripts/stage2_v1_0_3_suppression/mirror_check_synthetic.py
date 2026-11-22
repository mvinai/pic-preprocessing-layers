#!/usr/bin/env python3
"""mirror_check_synthetic.py: the v1.0.2 -> v1.0.3 pipeline on the default synthetic mirror (SYNTHETIC ONLY), whose
v1.0.1 label pattern reproduces the real cells flagged by Stage 3 v1.0.3 (train invalid-age / adolescent / >= 18
cells; CHARTEVENTS and MICROBIOLOGYEVENTS reconciliation).

  1. mirror v1.0.1-format output (synthetic_s2_v103_mirror.write_case) and links v3 (make_stage2_links_v1_0_3.py)
  2. v1.0.2 mirror release with module v1.3 (reference copy), v1.0.2 configuration
  3. Stage 2 alone, independent check (stage2_independent_check.py): v1.0.2 mirror with links v3 + public values
     (expect failures, including the flagged families), and with the Stage 3 links v2 internal rows + public table sizes
  4. v1.0.3 mirror release with module v1.4 (--carry-json, links v3, --public incl. the invalid-DOB total and the two
     'derived' cells), with the decided D3 (CHARTEVENTS / MICROBIOLOGYEVENTS sizes published again; the mirror's
     public file) and, as a control, the same file with those sizes kept S
  5. the same independent checks on each v1.0.3 mirror release (expect 0 failing)
Exit 0 when both v1.0.3 releases pass and the v1.0.2 release fails as expected."""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
S3L = next((q for q in (HERE.parent / "stage3_builder", HERE.parent / "stage3_spec" / "stage3_builder") if (q / "03_build_real_layers.py").is_file()), HERE.parent / "stage3_builder") / "stage3_links_v2.csv"
sys.path.insert(0, str(HERE))
import synthetic_s2_v103_mirror as MR  # noqa: E402
import stage2_independent_check as IC  # noqa: E402

PY = sys.executable
FLAGGED = re.compile(r"/(layer_a/dob_breakdown/by_split/train/(n_age_days|n_invalid|status/n_valid)|layer_a/train/age_band/"
                     r"(n_age_invalid|n_ge18)|reconciliation/(CHARTEVENTS|MICROBIOLOGYEVENTS)/(n_cohort_admission_time_"
                     r"missing|n_no_cohort_admission|n_total)|strata_denominators/age_band/train/(adolescent|age_invalid)/"
                     r"(n_eligible_72h|n_encounters))$")


def sh(cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def check(label, rel, links, public, layout):
    lines = []
    rc = IC.check(str(rel), str(links), str(public) if public else None, str(layout), out=lines.append)
    fam = sorted(l.split()[-1] for l in lines[1:] if FLAGGED.search(l.split()[-1]))
    print(f"[{label}] exit {rc}: {lines[0]}" + (f"; of the 15 flagged families failing: {len(fam)}" if rc else ""))
    for f in fam:
        print(f"    {f}")
    return rc


def main(work: Path) -> int:
    work = work.resolve()
    if work.exists():
        shutil.rmtree(work)
    mr = MR.write_case(work)
    with open(mr["public"], "a", encoding="utf-8") as f:
        f.write("/layer_a/dob_breakdown/status/n_valid,derived,cohort size (frozen) minus the public invalid-DOB count\n")
        f.write("/layer_a/dob_breakdown/n_age_days,derived,equals the valid-DOB count\n")
    pub_keep = work / "public_keepS_control.csv"      # control only: D3 decided 'publish' (Vinai, Nov 17, 2026)
    pub_keep.write_text(mr["public"].read_text(encoding="utf-8").replace(",publish,public table size",
                                                                         ",,public table size"), encoding="utf-8")
    pub_sizes = work / "public_table_sizes.csv"
    pub_sizes.write_text("".join(l + "\n" for l in mr["public"].read_text(encoding="utf-8").splitlines()
                                 if l.startswith(("pointer,", "/reconciliation/"))), encoding="utf-8")
    rc, log = sh([PY, str(HERE / "make_stage2_links_v1_0_3.py"), str(mr["json"]), "-", str(work / "links_v1_0_3.csv")],
                 HERE)
    print(log.strip())
    shutil.copy(HERE / "reference" / "stage2_v1_0_2_suppress_v1_3.py", work / "s13.py")
    drops = HERE / "stage2_v1_0_3_drops.csv"
    rc, log = sh([PY, "s13.py", "--in-json", mr["json"].name, "--in-md", mr["md"].name, "--links", mr["links_v102"].name,
                  "--frozen", "frozen.csv", "--drops", str(drops), "--out-dir", "out_v102"], work)
    print(f"[v1.0.2 mirror] module v1.3 run exit {rc}")
    r102 = work / "out_v102" / "feature_sparsity_frequency_audit_v1.0.2.json"
    a = check("v1.0.2 mirror, alone, links v3 + public", r102, work / "links_v1_0_3.csv", mr["public"], mr["json"])
    b = check("v1.0.2 mirror, alone, Stage 3 links v2 internal rows + public table sizes", r102, S3L, pub_sizes,
              mr["json"])
    ok = a == 1 and b == 1
    for label, pub, od in (("control: table sizes kept S", pub_keep, "out_v103_keepS"),
                           ("decided D3: table sizes published", mr["public"], "out_v103")):
        rc, log = sh([PY, str(HERE / "stage2_v1_0_3_suppress.py"), "--in-json", mr["json"].name, "--in-md",
                      mr["md"].name, "--links", "links_v1_0_3.csv", "--frozen", "frozen.csv", "--drops", str(drops),
                      "--public", pub.name, "--carry-json", str(r102), "--out-dir", od], work)
        gate = [l.strip() for l in log.splitlines() if "release gate:" in l or "RELEASE GATE" in l]
        print(f"[v1.0.3 mirror, {label}] module v1.4 run exit {rc}; {gate[-1] if gate else log[-500:]}")
        if rc:
            ok = False
            continue
        r = work / od / "feature_sparsity_frequency_audit_v1.0.3.json"
        a = check(f"v1.0.3 mirror ({label}), alone, links v3 + public", r, work / "links_v1_0_3.csv", pub, mr["json"])
        b = check(f"v1.0.3 mirror ({label}), alone, Stage 3 links v2 internal rows + public table sizes", r, S3L,
                  pub_sizes, mr["json"])
        ok = ok and a == 0 and b == 0
    print("MIRROR SYNTHETIC CHECK: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="s2v103_mirror_"))))
