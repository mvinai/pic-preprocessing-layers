"""Synthetic Stage 2 v1.0.1-format output that mirrors the cross-release leak structures of the real release
(SYNTHETIC VALUES ONLY; tests and the box demo; never real data).

It reuses the Builder of synthetic_v101_output.py and writes the structures behind the Stage 3 v1.0.3 findings:

  L5  layer_a/<split>/age_band (8 bands; no encounter >= 18 years, so 'ge18' is no strata level) beside
      strata_denominators/age_band/<split>/<level> (n_encounters, n_eligible_72h, n_open_window_ge_6h) and
      strata/age_band/<split>/<level>/<channel>/n_ever_layer_b (heart_rate seen in every encounter, so its hidden
      complement is 0 and the published child is a lower bound equal to the level size); first care unit and
      recent-surgery strata beside their layer_a / surgery counts
  L4/L6 layer_a/dob_breakdown overall, by split, by first care unit, by split x care unit; weight imputation (+ age_invalid)
  L3  reconciliation/<table> (in-stay per split, outside stay, two small not-in-stay classes, total = public table size)
  plus every Stage 2 cell that stage3_links_v2.csv names (row status of 40 channels, admission-interval bands and
      negatives, diagnosis code status, surgery link classes, Layer C surgeries), and an NMB-like validation subtree
      whose two '<10' strata are fixed by a frozen total (the v1.0.2 drop).

`values` (dict pointer -> int) from a Stage 3 audit build (`s3`) makes the release consistent with that build under
stage3_links_v2.csv (joint test); without it the default 'mirror' values reproduce the real v1.0.1 label pattern of the
15 cells Stage 3 v1.0.3 flagged (train invalid-age / adolescent / >= 18; CHARTEVENTS and MICROBIOLOGYEVENTS
reconciliation). v1.0.1 labels: primary '<10' (zeros included) and the v1.0.1 complementary rule (no single suppressed
cell or small hidden difference recoverable from one identity), plus a few forced 'S' labels as in the real file."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import stage2_v1_0_3_suppress as T
from synthetic_v101_output import Builder

SPLITS = ("train", "validation", "test")
BANDS = ("neonate", "infant", "toddler", "preschool", "school", "adolescent", "ge18", "age_invalid")
UNITS = ("CICU", "General ICU", "NICU", "PICU", "SICU")
TABLES = ("CHARTEVENTS", "LABEVENTS", "PRESCRIPTIONS", "MICROBIOLOGYEVENTS")
L1 = ("temperature heart_rate respiratory_rate spo2 urine_output fluid_intake glucose_bedside pain_score dbp sbp "
      "albumin alt ast bun bilirubin_direct ggt ldh magnesium phosphate bilirubin_total anc hemoglobin "
      "lymphocytes_abs platelets wbc ptt d_dimer fibrinogen inr base_excess chloride bicarbonate potassium lactate "
      "sodium pco2 ph sao2_bloodgas crp cystatin_c").split()
STRATA_CH = ("heart_rate", "temperature", "lactate")
NMB = "/drug_groups/nmb/validation"
PUBLISH_TABLES = ("CHARTEVENTS", "MICROBIOLOGYEVENTS")   # decision D3 (Vinai, Nov 17, 2026)


def _dist(v: int, n: int) -> list[int]:
    q = v // n
    return [v - q * (n - 1)] + [q] * (n - 1)


# default mirror values (invented; shaped like the real pattern, never real figures)
MIRROR = {
    "splits": {"train": 700, "validation": 150, "test": 150},
    "bands": {"train": [150, 200, 120, 90, 80, 40, 0, 20], "validation": [35, 40, 30, 20, 15, 6, 0, 4],
              "test": [34, 42, 28, 20, 16, 7, 0, 3]},
    # missing, unparseable, negative (lt1d, 1to30d, gt30d) per split; invalid = age_invalid band
    "dob": {"train": (3, 2, (12, 2, 1)), "validation": (1, 0, (3, 0, 0)), "test": (0, 1, (2, 0, 0))},
    "recon": {"CHARTEVENTS": ((500000, 100000, 100000), 20000, 3, 2),
              "LABEVENTS": ((300000, 60000, 60000), 9000, 40, 12),
              "PRESCRIPTIONS": ((80000, 16000, 15000), 3000, 25, 11),
              "MICROBIOLOGYEVENTS": ((3000, 700, 650), 400, 4, 1)},
}


def build(values: dict | None = None, splits: dict | None = None, mirror: dict | None = None) -> dict:
    """values: Stage 3 audit pointer -> true value (joint test); None = default mirror values."""
    V = values or {}
    M = mirror or MIRROR
    sp_n = dict(splits or M["splits"])
    b = Builder()
    b.doc["notes"] = {"synthetic": "synthetic mirror of the v1.0.1 structure (test data only)"}
    force: dict[str, str] = {}
    g = lambda p, d: int(V.get(p, d))  # noqa: E731
    tot = b.cell("/denominators/n_encounters", sum(sp_n.values()))
    enc = {s: b.cell(f"/denominators/{s}/n_encounters", sp_n[s]) for s in SPLITS}
    b.ident(tot, list(enc.values()))
    b.doc["cohort_checks"] = {"encounters": sum(sp_n.values())}
    # ---- age bands per split (joint: spread from the Stage 3 invalid total; mirror: fixed pattern)
    if values is not None:
        inv_tot = g("/age/n_invalid", 3)
        inv = dict(zip(SPLITS, _dist(inv_tot, 3)))
        bands = {}
        for s in SPLITS:
            rest = sp_n[s] - inv[s]
            share = [0.2, 0.3, 0.2, 0.12, 0.12]
            bb = [int(rest * x) for x in share]
            bb.append(rest - sum(bb))                      # adolescent
            bands[s] = bb + [0, inv[s]]
        miss_t = g("/age/invalid_detail/n_missing", 1)
        unp_t = g("/age/invalid_detail/n_unparseable", 0)
        nb_t = [g(f"/age/invalid_detail/negative_magnitude/n_{k}", 0) for k in ("lt1d", "1to30d", "gt30d")]
        assert miss_t + unp_t + sum(nb_t) == inv_tot, "Stage 3 invalid-age detail does not add up"
        # spread the classes over splits in order (deterministic)
        pool = ["missing"] * miss_t + ["unparseable"] * unp_t + ["lt1d"] * nb_t[0] + ["1to30d"] * nb_t[1] + \
               ["gt30d"] * nb_t[2]
        dob = {}
        k = 0
        for s in SPLITS:
            take = pool[k:k + inv[s]]
            k += inv[s]
            dob[s] = (take.count("missing"), take.count("unparseable"),
                      (take.count("lt1d"), take.count("1to30d"), take.count("gt30d")))
    else:
        bands = {s: list(M["bands"][s]) for s in SPLITS}
        dob = M["dob"]
    for s in SPLITS:
        assert sum(bands[s]) == sp_n[s], s
        m_, u_, nb = dob[s]
        assert m_ + u_ + sum(nb) == bands[s][7], s
    levels = [bd for i, bd in enumerate(BANDS) if sum(bands[s][i] for s in SPLITS) > 0]
    # ---- denominators per split: eligible 72 h and open window >= 6 h (with hidden complements)
    e72, w6 = {}, {}
    for s in SPLITS:
        e72[s] = b.cell(f"/denominators/{s}/eligible/72h/n_eligible", 0)
        w6[s] = b.cell(f"/denominators/{s}/n_open_window_ge_6h", 0)
    # ---- layer_a age bands + strata age band (+ channel children)
    for s in SPLITS:
        la = {bd: b.cell(f"/layer_a/{s}/age_band/n_{bd}", bands[s][i]) for i, bd in enumerate(BANDS)}
        b.ident(enc[s], list(la.values()))
        st_enc, st_e72, st_w6 = [], [], []
        for lv in levels:
            n = bands[s][BANDS.index(lv)]
            base = f"/strata_denominators/age_band/{s}/{lv}"
            ce = b.cell(f"{base}/n_encounters", n)
            el = b.cell(f"{base}/n_eligible_72h", max(0, n - max(1, n // 6)))
            ow = b.cell(f"{base}/n_open_window_ge_6h", max(0, n - max(1, n // 10)))
            b.subset(el, ce)
            b.subset(ow, ce)
            st_enc.append(ce), st_e72.append(el), st_w6.append(ow)
            for ch, frac in zip(STRATA_CH, (1.0, 0.9, 0.4)):
                nv = int(round(n * frac))
                if ch == "temperature":
                    nv = max(0, n - max(12, n // 10)) if n >= 24 else n
                c = b.cell(f"/strata/age_band/{s}/{lv}/{ch}/n_ever_layer_b", nv)
                b.subset(c, ce)
        b.raw[e72[s]] = sum(b.raw[c] for c in st_e72)
        b.raw[w6[s]] = sum(b.raw[c] for c in st_w6)
        b.ident(e72[s], st_e72)
        b.ident(w6[s], st_w6)
        b.subset(e72[s], enc[s])
        b.subset(w6[s], enc[s])
        b.ident(enc[s], st_enc)
    # ---- first care unit (layer_a + strata) and recent surgery (surgery/layer_a + strata)
    cu = {}
    for s in SPLITS:
        n = sp_n[s]
        parts = [int(n * x) for x in (0.2, 0.1, 0.35, 0.25)]
        parts.append(n - sum(parts))
        cu[s] = dict(zip(UNITS, parts))
        la = [b.cell(f"/layer_a/{s}/first_careunit/n_{u}", cu[s][u]) for u in UNITS]
        b.ident(enc[s], la)
        stc = []
        for u in UNITS:
            ce = b.cell(f"/strata_denominators/first_careunit/{s}/{u}/n_encounters", cu[s][u])
            stc.append(ce)
        b.ident(enc[s], stc)
        rec = max(12, n // 5)
        r = b.cell(f"/surgery/layer_a/{s}/n_encounters_recent", rec)
        b.subset(r, enc[s])
        sr = b.cell(f"/strata_denominators/surgical_status/{s}/recent_surgery/n_encounters", rec)
        sn = b.cell(f"/strata_denominators/surgical_status/{s}/no_recent_surgery/n_encounters", n - rec)
        b.ident(enc[s], [sn, sr])
    # ---- DOB breakdown (overall, by split, by care unit, by split x care unit) + weight imputation
    d = "/layer_a/dob_breakdown"

    def dob_block(base, counts, total):
        m_, u_, nb = counts
        st = {"missing": m_, "unparseable": u_, "negative": sum(nb), "valid": b.raw[total] - m_ - u_ - sum(nb)}
        cells = {k: b.cell(f"{base}/status/n_{k}", v) for k, v in st.items()}
        b.ident(total, list(cells.values()))
        ni = b.cell(f"{base}/n_invalid", m_ + u_ + sum(nb))
        b.ident(ni, [cells["missing"], cells["unparseable"], cells["negative"]])
        na = b.cell(f"{base}/n_age_days", st["valid"])
        b.ident(na, [cells["valid"]])
        bins = [b.cell(f"{base}/negative_bins/n_{k}", v) for k, v in zip(("lt1d", "1to30d", "gt30d"), nb)]
        b.ident(cells["negative"], bins)
        return ni
    tot_dob = (sum(dob[s][0] for s in SPLITS), sum(dob[s][1] for s in SPLITS),
               tuple(sum(dob[s][2][i] for s in SPLITS) for i in range(3)))
    ni_all = dob_block(d, tot_dob, tot)
    ni_s = {s: dob_block(f"{d}/by_split/{s}", dob[s], enc[s]) for s in SPLITS}

    def unit_split(counts):   # invalid encounters by care unit: negatives to NICU, missing / unparseable to PICU
        m_, u_, nb = counts
        out = {u: (0, 0, 0) for u in UNITS}
        out["NICU"] = (0, 0, sum(nb))
        out["PICU"] = (m_, u_, 0)
        return out

    def unit_block(base, c):
        cells = {k: b.cell(f"{base}/status/n_{k}", v) for k, v in zip(("missing", "unparseable", "negative"), c)}
        t = b.cell(f"{base}/status/n_total", sum(c))
        b.ident(t, list(cells.values()))
        ni = b.cell(f"{base}/n_invalid", sum(c))
        b.ident(ni, list(cells.values()))
        return ni
    uall = {u: [0, 0, 0] for u in UNITS}
    for s in SPLITS:
        us = unit_split(dob[s])
        parts = []
        for u in UNITS:
            parts.append(unit_block(f"{d}/by_split_and_first_careunit/{s}/{u}", us[u]))
            for i in range(3):
                uall[u][i] += us[u][i]
        b.ident(ni_s[s], parts)
    b.ident(ni_all, [unit_block(f"{d}/by_first_careunit/{u}", uall[u]) for u in UNITS])
    inv_all = b.raw[ni_all]
    wsrc = ("charted", "imputed_age_sex", "imputed_age_band", "imputed_sex_only", "imputed_overall")
    n_all = b.raw[tot]
    if values is not None:      # joint: every invalid-age encounter has a charted weight (Stage 3 'le' links hold)
        w_inv = [inv_all, 0, 0, 0, 0]
    else:
        w_inv = [max(0, inv_all - 7), 0, 0, min(6, inv_all), 1 if inv_all >= 7 else 0]
        w_inv[0] = inv_all - sum(w_inv[1:])
    w_all = [int(n_all * 0.8), int(n_all * 0.12), 25, 15 + w_inv[3], w_inv[4] + 0]
    w_all[0] = n_all - sum(w_all[1:])
    b.partition(f"{d}/weight_imputation", dict(zip(wsrc, w_all)))
    wi = [b.cell(f"{d}/weight_imputation/age_invalid/n_{k}", v) for k, v in zip(wsrc, w_inv)]
    b.ident(ni_all, wi)
    # ---- reconciliation (public table sizes)
    for t in TABLES:
        base = f"/reconciliation/{t}"
        if values is not None:
            ins = _dist(g(f"/event_linkage/{t}/n_rows_in_stay", 3000), 3)
            outs = _dist(g(f"/event_linkage/{t}/n_rows_outside_every_window", 30), 3)
            out_stay, no_adm, t_miss = outs
        else:
            ins, out_stay, no_adm, t_miss = M["recon"][t]
        parts = {**{f"in_stay_{s}": x for s, x in zip(SPLITS, ins)}, "cohort_admission_outside_stay": out_stay,
                 "no_cohort_admission": no_adm, "cohort_admission_time_missing": t_miss}
        b.partition(base, parts)
    # ---- the Stage 2 cells named by stage3_links_v2.csv (joint test; also present in the mirror)
    for ch in L1:
        base = f"/map_checks/row_status/{ch}"
        meas = g(f"/plausibility/{ch}/train/n_rows_in_gate", 1000)
        oog = g(f"/plausibility/{ch}/train/n_rows_out_of_gate", 20) + int(V.get(f"/plausibility/{ch}/train/n_rows_zero_as_missing", 0))
        b.partition(base, {"measurement": meas, "out_of_gate": oog, "nonnumeric": 12, "blank": 3, "unit_pending": 0,
                           "unit_unrecognized": 0})
    miss = _dist(g("/admission_interval/withheld/n_missing", 0), 3)
    neg = _dist(g("/admission_interval/n_negative_clipped_to_0", 0), 3)
    dxm = _dist(g("/diagnosis/n_missing", 0), 3)
    lc = _dist(g("/surgery/timing/n_layer_c", 60), 3)
    for j, s in enumerate(SPLITS):
        h = f"/layer_a/{s}/hours_admit_to_icu_band"
        hm = b.cell(f"{h}/n_missing", miss[j])
        hl = b.cell(f"{h}/n_lt24h", sp_n[s] - miss[j])
        b.ident(enc[s], [hm, hl])
        ng = b.cell(f"/layer_a/{s}/n_hours_admit_to_icu_negative", neg[j])
        b.subset(ng, enc[s])
        c = f"/layer_a/{s}/diagnosis_code_status"
        b.ident(enc[s], [b.cell(f"{c}/n_ok", sp_n[s] - dxm[j]), b.cell(f"{c}/n_missing", dxm[j])])
        b.cell(f"/surgery/layer_c_extra/{s}/n_surgeries", lc[j])
    lk = g("/surgery/n_linked", 40)
    b.partition("/surgery/link_classes", {"linked": lk, "other_hadm": 25})
    # ---- NMB-like validation subtree (two '<10' strata fixed by a frozen total: the v1.0.2 drop)
    nt = b.cell(f"{NMB}/n_active_72h", 15)
    b.ident(nt, [b.cell(f"{NMB}/strata/surgical_status/no_recent_surgery/n_active_72h", 8),
                 b.cell(f"{NMB}/strata/surgical_status/recent_surgery/n_active_72h", 7)])
    # ---- v1.0.1 labels: the real file's forced S labels, then the v1.0.1 complementary rule
    if values is None:
        for p in (f"{d}/by_split/train/status/n_valid", f"{d}/status/n_valid",
                  "/strata_denominators/age_band/train/adolescent/n_encounters",
                  "/strata_denominators/age_band/train/adolescent/n_eligible_72h"):
            force[p] = "S"
    b.forced.update(force)
    _v101_complementary(b, frozen={tot, *enc.values(), nt})
    return {"doc": b.finish(), "raw": b.raw}


def _v101_complementary(b: Builder, frozen: set[str]) -> None:
    """The v1.0.1 gate's rule: no suppressed cell (and no small hidden difference) is recoverable from ONE identity.
    While an identity has exactly one unknown that is a suppressed count (or a hidden difference of 1-9), suppress its
    smallest published, non-frozen member."""
    for _ in range(200):
        changed = False
        for t, parts in b.identities:
            cells = [t] + list(parts)
            unk = [c for c in cells if c.startswith("hidden:") or b.status(c) != "pub"]
            if len(unk) != 1:
                continue
            u = unk[0]
            if u.startswith("hidden:") and not (1 <= b.raw[u] <= 9):
                continue
            pick = sorted((b.raw[c], c) for c in cells if not c.startswith("hidden:") and b.status(c) == "pub"
                          and c not in frozen)
            if pick:
                b.forced[pick[0][1]] = "S"
                changed = True
        if not changed:
            return
    raise RuntimeError("v1.0.1-style complementary rule did not settle")


FROZEN = [("/denominators/n_encounters", "cohort size public"),
          ("/denominators/train/n_encounters", "split size public"),
          ("/denominators/validation/n_encounters", "split size public"),
          ("/denominators/test/n_encounters", "split size public"),
          (f"{NMB}/n_active_72h", "synthetic NMB-like total: already public (fixes its two '<10' strata)")]


def write_case(folder: Path, values: dict | None = None, splits: dict | None = None) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    r = build(values, splits)
    doc = r["doc"]
    pj = folder / "feature_sparsity_frequency_audit.json"
    pm = folder / "feature_sparsity_frequency_audit.md"
    pj.write_text(T.dump_json(doc), encoding="utf-8", newline="")
    pm.write_text(T.render_markdown(doc), encoding="utf-8", newline="")
    pf = folder / "frozen.csv"
    with open(pf, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["pointer", "reason"])
        w.writerows(FROZEN)
    pl = folder / "links_v1_0_2_mirror.csv"           # the v1.0.2 run of the mirror: no cross-table links
    pl.write_text("total,parts,note\n", encoding="utf-8")
    raw = r["raw"]
    pp = folder / "public_TEAM_INTERNAL.csv"          # what v1.0.3 treats as already known (mirror values)
    with open(pp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["pointer", "value", "reason"])
        for s in SPLITS:
            w.writerow([f"/layer_a/{s}/age_band/n_ge18", 0, "structural zero (no strata level)"])
        for t in TABLES:     # decision D3: CHARTEVENTS and MICROBIOLOGYEVENTS sizes are published again
            w.writerow([f"/reconciliation/{t}/n_total", "publish" if t in PUBLISH_TABLES else "", "public table size"])
        w.writerow(["/layer_a/dob_breakdown/n_invalid", raw["/layer_a/dob_breakdown/n_invalid"],
                    "synthetic stand-in for the v1.0 stop-event count"])
    return {"json": pj, "md": pm, "frozen": pf, "links_v102": pl, "public": pp, "raw": raw}


if __name__ == "__main__":
    import sys
    out = write_case(Path(sys.argv[1] if len(sys.argv) > 1 else "mirror_case"))
    print(json.dumps({k: str(v) for k, v in out.items() if k != "raw"}, indent=1))
