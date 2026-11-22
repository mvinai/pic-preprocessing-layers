"""pytest for check_encounter_harmonization_consistency.py. SYNTHETIC DATA ONLY (built below; no PIC file is read).

The generator builds raw ICU stays (single stays, overlaps, nested last stays, gaps <= 24 h, same-HADM readmissions
> 24 h, multi-HADM merges), derives the canonical and harmonization files with the documented rules, then plants
specific defects and checks that the aggregate counters detect them."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_encounter_harmonization_consistency as C  # noqa: E402

FMT = "%Y-%m-%d %H:%M:%S"
SENTINEL = "987654321"


def make_stays(n_subjects=400, seed=3):
    rng = np.random.default_rng(seed)
    rows, hadm = [], 700000
    base = pd.Timestamp("2016-01-01 06:00:00")
    for i in range(n_subjects):
        s = SENTINEL if i == 0 else str(200000 + i)
        t = base + pd.Timedelta(days=int(rng.integers(0, 1500)), hours=int(rng.integers(0, 24)))
        hadm += 1
        h = str(hadm)
        kind = rng.choice(["single", "overlap", "nested", "gap", "readmit", "multi_hadm"], p=[0.6, 0.06, 0.06, 0.08, 0.12, 0.08])
        if i == 1:
            kind = "exact24"
        if kind == "exact24":
            rows.append((s, h, t, t + pd.Timedelta(hours=24)))
            continue
        d1 = float(rng.uniform(5, 300))
        rows.append((s, h, t, t + pd.Timedelta(hours=d1)))
        o1 = t + pd.Timedelta(hours=d1)
        if kind == "overlap":
            rows.append((s, h, o1 - pd.Timedelta(hours=2), o1 + pd.Timedelta(hours=float(rng.uniform(5, 50)))))
        elif kind == "nested":
            rows.append((s, h, t + pd.Timedelta(hours=1), t + pd.Timedelta(hours=min(3.0, d1 - 1))))
        elif kind == "gap":
            g = float(rng.uniform(0.5, 23.5))
            rows.append((s, h, o1 + pd.Timedelta(hours=g), o1 + pd.Timedelta(hours=g + float(rng.uniform(5, 80)))))
        elif kind == "readmit":
            g = float(rng.uniform(30, 300))
            rows.append((s, h, o1 + pd.Timedelta(hours=g), o1 + pd.Timedelta(hours=g + float(rng.uniform(5, 80)))))
        elif kind == "multi_hadm":
            hadm += 1
            g = float(rng.uniform(0.5, 20))
            rows.append((s, str(hadm), o1 + pd.Timedelta(hours=g), o1 + pd.Timedelta(hours=g + 30)))
    st = pd.DataFrame(rows, columns=["SUBJECT_ID", "HADM_ID", "INTIME", "OUTTIME"])
    st["INTIME"], st["OUTTIME"] = st["INTIME"].dt.floor("s"), st["OUTTIME"].dt.floor("s")
    st.insert(0, "ICUSTAY_ID", [str(300000 + k) for k in range(len(st))])
    return st


def build_encounters(st):
    st = st.sort_values(["SUBJECT_ID", "INTIME", "OUTTIME"]).reset_index(drop=True)
    enc, cur, end = [], None, None
    n = -1
    for s, i, o in zip(st["SUBJECT_ID"], st["INTIME"], st["OUTTIME"]):
        if s != cur or i > end + pd.Timedelta(hours=24):
            n += 1; cur, end = s, o
        else:
            end = max(end, o)
        enc.append(n)
    st["e"] = enc
    st["dur"] = (st["OUTTIME"] - st["INTIME"]).dt.total_seconds() / 3600
    out = []
    for e, d in st.groupby("e"):
        t0, last = d["INTIME"].min(), d["OUTTIME"].max()
        # union of intervals -> gap hours
        iv = sorted(zip(d["INTIME"], d["OUTTIME"]))
        gaps, gh, run = 0, 0.0, iv[0][1]
        for a, b in iv[1:]:
            if a > run:
                gaps += 1; gh += (a - run).total_seconds() / 3600
            run = max(run, b)
        hl = sorted(set(d["HADM_ID"]))
        out.append({"SUBJECT_ID": d["SUBJECT_ID"].iloc[0], "HADM_ID": d["HADM_ID"].iloc[0], "hadm_list": hl,
                    "ICUSTAY_ID": d["ICUSTAY_ID"].iloc[0], "encounter_n_stays": len(d), "n_linked_admissions": len(hl),
                    "INTIME": t0, "OUTTIME": last, "span": (last - t0).total_seconds() / 3600, "gaps": gaps, "gh": gh})
    e = pd.DataFrame(out).sort_values(["SUBJECT_ID", "INTIME"]).reset_index(drop=True)
    e["seq"] = e.groupby("SUBJECT_ID").cumcount()
    e["encounter_id"] = e["SUBJECT_ID"] + "_" + (e["seq"] + 1).astype(str)
    return e


def write_project(root: Path, mutate=None, st=None):
    st = make_stays() if st is None else st
    e = build_encounters(st)
    rng = np.random.default_rng(5)
    subjects = sorted(e["SUBJECT_ID"].unique())
    split = {s: ("train" if rng.random() < 0.7 else ("val" if rng.random() < 0.5 else "test")) for s in subjects}
    sp_long = {"train": "train", "val": "validation", "test": "test"}
    e["split"] = e["SUBJECT_ID"].map(split).map(sp_long)
    # deaths: on some last encounters, plus 3 deaths on a non-last encounter of the same HADM (attribution case)
    e["Died"] = 0
    last = ~e.duplicated("SUBJECT_ID", keep="last")
    e.loc[last & (rng.random(len(e)) < 0.08), "Died"] = 1
    nxt_same = [bool(set(a) & set(b)) for a, b in zip(e["hadm_list"], e.groupby("SUBJECT_ID")["hadm_list"].shift(-1).map(lambda v: v if isinstance(v, list) else []))]
    cand = [i for i, x in enumerate(nxt_same) if x][:3]
    e.loc[cand, "Died"] = 1
    e["death_attribution_flag"] = ((e["Died"] == 1) & ~last).astype(int)
    tfmt = lambda x: x.strftime(FMT)  # noqa: E731
    fmts = ["json", "repr", "plain"]
    def hl_text(lst, k):
        f = fmts[k % 3] if len(lst) > 1 else "plain"
        return {"json": "[" + ", ".join(lst) + "]", "repr": str(lst), "plain": ",".join(lst)}[f]
    can = pd.DataFrame({
        "encounter_id": e["encounter_id"], "SUBJECT_ID": e["SUBJECT_ID"], "HADM_ID": e["HADM_ID"],
        "hadm_id_list": [hl_text(l, k) for k, l in enumerate(e["hadm_list"])], "ICUSTAY_ID": e["ICUSTAY_ID"],
        "encounter_n_stays": e["encounter_n_stays"].astype(str), "n_linked_admissions": e["n_linked_admissions"].astype(str),
        "INTIME": e["INTIME"].map(tfmt), "OUTTIME": e["OUTTIME"].map(tfmt), "icu_los_hours": e["span"].round(4).astype(str),
        "hospital_admit_time": (e["INTIME"] - pd.Timedelta(hours=5)).map(tfmt),
        "hospital_disch_time": (e["OUTTIME"] + pd.Timedelta(days=4)).map(tfmt), "hospital_los_hours": "100",
        "age_years_at_intime": "3.5", "HOSPITAL_EXPIRE_FLAG_any": e["Died"].astype(str), "Died": e["Died"].astype(str),
        "split": e["split"]})
    harm = can.drop(columns=["split"]).copy()
    t0, out = e["INTIME"], e["OUTTIME"]
    harm["t0"] = can["INTIME"]
    for L in (24, 48, 72):
        harm[f"in_icu_at_{L}h"] = (out > t0 + pd.Timedelta(hours=L)).astype(int).astype(str)
    harm["gap_count_between_stays"] = e["gaps"].astype(str)
    harm["gap_hours_total"] = e["gh"].round(4).astype(str)
    harm["first_careunit"] = "PICU"
    harm["prior_encounters"] = e["seq"].astype(str)
    harm["death_attribution_flag"] = e["death_attribution_flag"].astype(str)
    harm["icu_time_hours"] = (e["span"] - e["gh"]).round(4).astype(str)
    harm["span_icu_hours"] = e["span"].round(4).astype(str)
    harm["split"] = can["split"]
    sf = pd.DataFrame({"SUBJECT_ID": subjects, "split": [split[s] for s in subjects]})
    raw = st.assign(INTIME=st["INTIME"].map(tfmt), OUTTIME=st["OUTTIME"].map(tfmt))
    if mutate:
        mutate(can, harm, sf, raw, e)
    root.mkdir(parents=True, exist_ok=True)
    paths = {k: root / f"{k}.csv" for k in ("canonical", "harm", "split", "icustays")}
    can.to_csv(paths["canonical"], index=False)
    harm.to_csv(paths["harm"], index=False)
    sf.to_csv(paths["split"], index=False)
    raw.to_csv(paths["icustays"], index=False)
    sha = {k: hashlib.sha256(paths[k].read_bytes()).hexdigest() for k in ("canonical", "harm")}
    exp = {"encounters": len(can), "subjects": len(subjects), "Died": int(e["Died"].sum()),
           "death_attribution_flag": int(e["death_attribution_flag"].sum()),
           "multi_hadm_encounters": int((e["hadm_list"].map(len) > 1).sum()),
           "span_gt_1440h": int((e["span"] > 1440).sum()), "valid_raw_stays": len(st),
           **{f"in_icu_at_{L}h": int((out > t0 + pd.Timedelta(hours=L)).sum()) for L in (24, 48, 72)},
           **{k: int((can["split"] == k).sum()) for k in ("train", "validation", "test")}}
    return paths, sha, exp, e


def run(tmp_path, mutate=None):
    paths, sha, exp, e = write_project(tmp_path, mutate)
    R = C.check(paths["canonical"], paths["harm"], paths["split"], paths["icustays"], sha["canonical"], sha["harm"], exp)
    return R, e


def zero_disagreements(sec: dict):
    return {k: v for k, v in sec.items() if isinstance(v, (int, str)) and v not in (0,) and "distribution" not in k}


def test_clean_synthetic_passes(tmp_path):
    R, e = run(tmp_path)
    s = R.sections
    assert all(v == 0 for v in s["3_canonical_vs_harmonization"]["disagreeing_rows_per_shared_column"].values())
    t = s["4_times"]
    assert t["t0_ne_harm_INTIME"] == 0 and t["abs(harm OUTTIME - t0 - span_icu_hours) > 1 min"] == 0
    assert t["icu_time_hours > span_icu_hours + 1 min"] == 0 and t["abs(span - icu_time - gap_hours_total) > 1 min"] == 0
    for L in (24, 48, 72):
        assert s["6_landmarks"][f"in_icu_at_{L}h != (OUTTIME > t0 + {L}h strictly)"] == 0
        assert s["6_landmarks"][f"in_icu_at_{L}h_total"]["matches"]
    assert s["6_landmarks"]["in_icu_at_24h != (OUTTIME >= t0 + 24h)"] == "<10"  # the planted exact-24 h encounter
    d = s["7_deaths"]
    assert d["Died_total"]["matches"] and d["death_attribution_flag_total"]["matches"]
    assert d["flag != (Died & not last of subject)"] == 0 and d["flag != (Died & not last of admission)"] == 0
    assert s["8_next_encounter"]["next_encounter_starts_within_24h_of_OUTTIME (must be 0)"] == 0
    assert s["8_next_encounter"]["prior_encounters != number of earlier encounters of the subject"] == 0
    assert s["9_split"]["canonical_split != split_file"] == 0 and s["9_split"]["harmonization_split != split_file"] == 0
    r = s["11_rebuild_from_ICUSTAYS"]
    assert r["rebuilt_encounters (running max OUTTIME, <=24 h)"]["matches"]
    assert r["harm_OUTTIME != rebuilt max OUTTIME"] == 0
    assert r["encounters_where_max_OUTTIME != last_stay_OUTTIME (nested last stay)"] != 0  # nested stays exist
    assert r["harm_OUTTIME != OUTTIME of the last stay by INTIME"] == r["encounters_where_max_OUTTIME != last_stay_OUTTIME (nested last stay)"]
    assert s["1_naming"]["documented_name -> real header"]["icu_span_hours"] == "span_icu_hours"


def test_detects_last_stay_outtime(tmp_path):
    """Plant the 3.0-style 'last stay OUTTIME' (shorter for nested last stays) in both files."""
    def mut(can, harm, sf, raw, e):
        st = raw.assign(i=pd.to_datetime(raw["INTIME"]), o=pd.to_datetime(raw["OUTTIME"])).sort_values(["SUBJECT_ID", "i"])
        for k, r in e.iterrows():
            d = st[(st["SUBJECT_ID"] == r["SUBJECT_ID"]) & (st["i"] >= r["INTIME"]) & (st["i"] <= r["OUTTIME"])]
            lo = d["o"].iloc[-1]
            if lo != r["OUTTIME"]:
                can.loc[k, "OUTTIME"] = harm.loc[k, "OUTTIME"] = lo.strftime(FMT)
    R, _ = run(tmp_path, mut)
    t, r = R.sections["4_times"], R.sections["11_rebuild_from_ICUSTAYS"]
    assert t["abs(harm OUTTIME - t0 - span_icu_hours) > 1 min"] != 0
    assert r["harm_OUTTIME != rebuilt max OUTTIME"] != 0 and r["harm_OUTTIME != OUTTIME of the last stay by INTIME"] == 0


def test_detects_split_and_flag_and_prior_defects(tmp_path):
    def mut(can, harm, sf, raw, e):
        harm.loc[5, "split"] = "test" if harm.loc[5, "split"] != "test" else "train"
        harm.loc[7, "in_icu_at_48h"] = "1" if harm.loc[7, "in_icu_at_48h"] == "0" else "0"
        harm.loc[9, "prior_encounters"] = "4"
        harm.loc[11, "OUTTIME"] = "2030-01-01 00:00:00"
    R, _ = run(tmp_path, mut)
    s = R.sections
    assert s["9_split"]["harmonization_split != split_file"] == "<10"
    assert s["3_canonical_vs_harmonization"]["disagreeing_rows_per_shared_column"]["OUTTIME"] == "<10"
    assert s["6_landmarks"]["in_icu_at_48h != (OUTTIME > t0 + 48h strictly)"] == "<10"
    assert s["8_next_encounter"]["prior_encounters != number of earlier encounters of the subject"] == "<10"


def test_patient_vs_admission_last_definition(tmp_path):
    """A death on a subject's non-last encounter whose later encounter is in a DIFFERENT HADM differs between the
    log's 'not last encounter of the patient' and the contract's 'not last encounter of the admission'."""
    def mut(can, harm, sf, raw, e):
        s_multi = e.groupby("SUBJECT_ID").filter(lambda d: len(d) >= 2)
        for k, r in s_multi.iterrows():
            nxt = e.loc[k + 1] if k + 1 in e.index and e.loc[k + 1, "SUBJECT_ID"] == r["SUBJECT_ID"] else None
            if nxt is not None and not set(r["hadm_list"]) & set(nxt["hadm_list"]):
                return  # none in the generator: readmissions share HADM
        # build one: give a readmission encounter a different HADM
        k = next(i for i, x in enumerate(e.duplicated("SUBJECT_ID", keep="last")) if x)
        can.loc[k + 1, "hadm_id_list"] = harm.loc[k + 1, "hadm_id_list"] = "999999"
        can.loc[k + 1, "HADM_ID"] = harm.loc[k + 1, "HADM_ID"] = "999999"
        can.loc[k, "Died"] = harm.loc[k, "Died"] = "1"
        harm.loc[k, "death_attribution_flag"] = "1"
    R, _ = run(tmp_path, mut)
    d = R.sections["7_deaths"]
    assert d["flag != (Died & not last of subject)"] == 0
    assert d["flag != (Died & not last of admission)"] == "<10"


def test_hash_mismatch_stops(tmp_path):
    paths, sha, exp, _ = write_project(tmp_path)
    with pytest.raises(C.InputError):
        C.check(paths["canonical"], paths["harm"], paths["split"], None, "0" * 64, sha["harm"], exp)


def test_output_has_no_ids_or_timestamps(tmp_path, capsys):
    paths, sha, exp, e = write_project(tmp_path)
    out = tmp_path / "o.json"
    rc = C.main(["--canonical", str(paths["canonical"]), "--harmonization", str(paths["harm"]), "--split-file", str(paths["split"]),
                 "--icustays", str(paths["icustays"]), "--expect-canonical-sha", sha["canonical"], "--expect-harm-sha", sha["harm"],
                 "--json-out", str(out)])
    text = capsys.readouterr().out + out.read_text(encoding="utf-8")
    assert rc == 0
    assert SENTINEL not in text and "200001" not in text and "700001" not in text
    assert not __import__("re").search(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}", text)
    for v in json.loads(out.read_text(encoding="utf-8")).values():
        for x in (v.values() if isinstance(v, dict) else []):
            assert not (isinstance(x, int) and not isinstance(x, bool) and 0 < x < 10)


def test_suppression_helper():
    assert C.sup(0) == 0 and C.sup(1) == "<10" and C.sup(9) == "<10" and C.sup(10) == 10


# ---------------------------------------------------------------- sections 12-13 (surgery timing, vitals -> ROW_ID)
def write_surgery(root: Path, st: pd.DataFrame, e: pd.DataFrame):
    """Planted truth for the v1.7 timing rules (anaesthesia end T; fallback surgery end):
    - every encounter: one surgery; on the first 12 encounters it ends at t0 + 3 h (recent_layer_a, ends after t0),
      otherwise at t0 + 7 h when OUT allows (during_icu_layer_c);
    - each multi-stay encounter with a gap starting after t0 + 7 h: one surgery ending mid-gap (Layer C, in_gap);
    - 15 first encounters: one surgery ending at t0 - 48 h (remote_same_hadm);
    - 12 rows with no anaesthesia end (fallback to surgery end), 3 rows with neither (no timing time);
    - 15 shared-OPER_ID pairs inside the encounter (Layer C); 3 pairs overlap (ambiguous vitals rows), 2 extra
      vitals rows outside every window (unassigned); unique-key vitals: one inside each window + 3 outside."""
    rows, vit, rid = [], [], 0
    truth = {"recent": 0, "recent_after_t0": 0, "layer_c": 0, "in_gap": 0, "remote": 0, "fallback": 0, "no_time": 0,
             "shared": 0, "reuse": 0, "sh_assigned": 0, "sh_amb": 0, "sh_unassigned": 0, "u_assigned": 0, "u_unassigned": 0}

    def add(s, h, o, T, dur=1.0, ae=True, e_ok=True):
        nonlocal rid
        rid += 1
        B = T - pd.Timedelta(hours=dur, minutes=20)
        rows.append({"ROW_ID": str(rid), "SUBJECT_ID": s, "HADM_ID": h, "VISIT_ID": "1", "OPER_ID": o,
                     "ANES_START_TIME": (B - pd.Timedelta(minutes=20)).strftime(FMT), "ANES_END_TIME": T.strftime(FMT) if ae else "",
                     "SURGERY_BEGIN_TIME": B.strftime(FMT),
                     "SURGERY_END_TIME": (T - pd.Timedelta(minutes=20)).strftime(FMT) if e_ok else "",
                     "SURGERY_NAME": "x", "ANES_METHOD": "g", "SURGERY_POSITION": "s"})
        return B - pd.Timedelta(minutes=20), T
    later_same = [bool(set(a) & set(b)) for a, b in zip(e["hadm_list"], e.groupby("SUBJECT_ID")["hadm_list"].shift(-1).map(lambda v: v if isinstance(v, list) else []))]
    stays = st.sort_values(["SUBJECT_ID", "INTIME"])
    n_fb = 0
    for k, r in e.iterrows():
        s, h, t0, out = r["SUBJECT_ID"], r["HADM_ID"], r["INTIME"], r["OUTTIME"]
        if k < 12:
            AS, AE = add(s, h, f"1{k}", t0 + pd.Timedelta(hours=3))
            truth["recent"] += 1; truth["recent_after_t0"] += 1; truth["reuse"] += int(later_same[k])
        elif out >= t0 + pd.Timedelta(hours=7):
            fb = n_fb < 12 and k % 5 == 0
            n_fb += int(fb)
            AS, AE = add(s, h, f"1{k}", t0 + pd.Timedelta(hours=7), ae=not fb)
            truth["layer_c"] += 1; truth["fallback"] += int(fb); truth["reuse"] += int(later_same[k])
            T_eff = t0 + pd.Timedelta(hours=7) - (pd.Timedelta(minutes=20) if fb else pd.Timedelta(0))
            ss = stays[stays["SUBJECT_ID"] == s]
            truth["in_gap"] += int(not ((ss["INTIME"] <= T_eff) & (T_eff <= ss["OUTTIME"])).any())  # short first stay
            # with the fallback, T = surgery end = t0 + 6 h 40 min: still Layer C
        else:
            continue
        vit.append((s, h, f"1{k}", AS + pd.Timedelta(minutes=10))); truth["u_assigned"] += 1
        if k in (20, 21, 22):
            vit.append((s, h, f"1{k}", AS - pd.Timedelta(hours=5))); truth["u_unassigned"] += 1
        if r["gaps"] > 0:
            d = stays[(stays["SUBJECT_ID"] == s) & (stays["INTIME"] >= t0) & (stays["INTIME"] <= out)]
            run = d["OUTTIME"].iloc[0]
            for a_, b_ in zip(d["INTIME"].iloc[1:], d["OUTTIME"].iloc[1:]):
                if a_ > run + pd.Timedelta(minutes=10):
                    mid = run + (a_ - run) / 2
                    if mid > t0 + pd.Timedelta(hours=7):
                        add(s, d["HADM_ID"].iloc[0], f"9{k}", mid, dur=0.02)
                        truth["layer_c"] += 1; truth["in_gap"] += 1; truth["reuse"] += int(later_same[k])
                    break
                run = max(run, b_)
    firsts = e[e["seq"] == 0].iloc[30:45]
    for k, r in firsts.iterrows():
        add(r["SUBJECT_ID"], r["HADM_ID"], f"7{k}", r["INTIME"] - pd.Timedelta(hours=48)); truth["remote"] += 1
    for k, r in e[e["seq"] == 0].iloc[50:53].iterrows():
        add(r["SUBJECT_ID"], r["HADM_ID"], f"8{k}", r["INTIME"] + pd.Timedelta(hours=1), ae=False, e_ok=False); truth["no_time"] += 1
    singles = e[(e["encounter_n_stays"] == 1) & (e["span"] > 14) & (e.index >= 60)].head(15)
    for j, (k, r) in enumerate(singles.iterrows()):
        s, h = r["SUBJECT_ID"], r["HADM_ID"]
        T1 = r["INTIME"] + pd.Timedelta(hours=9)
        T2 = T1 + (pd.Timedelta(minutes=30) if j < 3 else pd.Timedelta(hours=3))
        AS1, _ = add(s, h, f"5{k}", T1); AS2, _ = add(s, h, f"5{k}", T2)
        truth["shared"] += 1; truth["layer_c"] += 2; truth["reuse"] += 2 * int(later_same[k])
        for AS_ in (AS1, AS2):
            vit.append((s, h, f"5{k}", AS_ + pd.Timedelta(minutes=45)))
        if j < 3:
            truth["sh_amb"] += 2
        else:
            truth["sh_assigned"] += 2
        if j in (5, 6):
            vit.append((s, h, f"5{k}", T1 - pd.Timedelta(days=1))); truth["sh_unassigned"] += 1
    si = pd.DataFrame(rows)
    sv = pd.DataFrame(vit, columns=["SUBJECT_ID", "HADM_ID", "OPER_ID", "MONITOR_TIME"])
    sv["MONITOR_TIME"] = sv["MONITOR_TIME"].map(lambda x: x.strftime(FMT))
    sv.insert(0, "ROW_ID", [str(i + 1) for i in range(len(sv))])
    sv["ITEMID"], sv["VALUE"] = "1", "80"
    si.to_csv(root / "SURGERY_INFO.csv", index=False)
    sv.to_csv(root / "SURGERY_VITAL_SIGNS.csv", index=False)
    return truth


def test_surgery_sections(tmp_path):
    st = make_stays()
    paths, sha, exp, e = write_project(tmp_path, st=st)
    truth = write_surgery(tmp_path, st, e)
    R = C.check(paths["canonical"], paths["harm"], paths["split"], paths["icustays"], sha["canonical"], sha["harm"], exp,
                surgery_info=tmp_path / "SURGERY_INFO.csv", surgery_vitals=tmp_path / "SURGERY_VITAL_SIGNS.csv")
    s12, s13 = R.sections["12_surgery_vs_encounters"], R.sections["13_vitals_to_ROW_ID"]
    assert s12["SURGERY_INFO_ROW_ID_repeated"] == 0
    assert s12["shared_(SUBJECT_ID,HADM_ID,OPER_ID)_groups"] == C.sup(truth["shared"])
    assert truth["in_gap"] >= 10 and truth["layer_c"] >= 10 and truth["recent"] >= 10 and truth["remote"] >= 10
    assert s12["one_surgery_one_encounter_classes"] == C.partition({
        "recent_layer_a (end in [t0-24h, t0+6h]; Layer A only)": truth["recent"],
        "remote_same_hadm (end < t0-24h; Layer A)": truth["remote"],
        "during_icu_layer_c (end in (t0+6h, OUT])": truth["layer_c"], "post_icu_same_hadm": 0})
    assert s12["recent_ending_after_t0 (end in (t0, t0+6h]; would have been Layer C by begin time)"] == C.sup(truth["recent_after_t0"])
    assert s12["in_gap (Layer C surgeries ending in a merge gap; descriptor only, kept in Layer C)"] == C.sup(truth["in_gap"])
    assert s12["timing_fallback_to_SURGERY_END_TIME"] == C.sup(truth["fallback"])
    assert s12["no_timing_time (anaesthesia end and surgery end both missing)"] == C.sup(truth["no_time"])
    assert s12["assigned_surgeries_that_are_also_Layer_A_history_of_a_later_same_HADM_encounter"] == C.sup(truth["reuse"])
    sh = s13["shared_group_rows_by_anesthesia_window"]
    assert sh == C.partition({"assigned_to_exactly_one_ROW_ID": truth["sh_assigned"], "ambiguous_2+_windows": truth["sh_amb"],
                              "unassigned_0_windows": truth["sh_unassigned"]})
    assert sh["ambiguous_2+_windows"] == "<10" and sh["assigned_to_exactly_one_ROW_ID"] == "S"
    u = s13["unique_key_rows_by_anesthesia_window"]
    assert u == C.partition({"assigned_to_exactly_one_ROW_ID": truth["u_assigned"], "ambiguous_2+_windows": 0,
                             "unassigned_0_windows": truth["u_unassigned"]})
    assert s13["shared_group_rows_by_anesthesia_window_pm30min"]["unassigned_0_windows"] == "<10"
    assert s13["surgeries_whose_window_uses_a_surgery_time_fallback"] == C.sup(truth["fallback"])
    assert s13["surgeries_with_no_usable_window"] == C.sup(truth["no_time"])


def test_partition_complementary():
    assert C.partition({"a": 100, "b": 50, "c": 3}) == {"a": 100, "b": "S", "c": "<10"}
    assert C.partition({"a": 100, "b": 50, "c": 0}) == {"a": 100, "b": 50, "c": 0}
    assert C.partition({"a": 100, "b": 5}) == {"a": "S", "b": "<10"}
