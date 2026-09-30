# -*- coding: utf-8 -*-
"""Supplemental_Note_15_check_numbers.py: independent number check for Supplemental Note 15 (post hoc revision analysis).

Recomputes every number printed in
    results/final_notes/Supplemental_Note_15_matched_background_genuine_chemistry.md
from the stored simulation outputs in results/B_overadjust_sim/ (per-table replicate files where they exist,
summary files otherwise), writes one row per number to
    results/final_notes/Supplemental_Note_15_number_check.csv
(id, full-precision value, text as printed, source file, selection), and then extracts every numeric token from
the Note and reports any token that is not the printed text of a row of that table or a listed identifier or
design constant. Exit code 1 if an unexplained token is found.

Read-only with respect to results/B_overadjust_sim/ and the manuscript folder (whose Supplemental Data 4 and 5 are
used only as a cross-check of the observed values). Run from this folder: python Supplemental_Note_15_check_numbers.py
(set PYTHONDONTWRITEBYTECODE=1 to avoid a __pycache__). Written 2026-09-30.
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
W = HERE.parents[1] if HERE.name == "final_notes" else Path(
    r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
B = W / "results" / "B_overadjust_sim"
OUT = W / "results" / "final_notes"
NOTE = W.parent / "supplemental" / "Supplemental_Note_15_matched_background_genuine_chemistry.md"  # repository copy
if not NOTE.exists():
    NOTE = OUT / "Supplemental_Note_15_matched_background_genuine_chemistry.md"
CSV_OUT = OUT / "Supplemental_Note_15_number_check.csv"

ROWS: list[dict] = []
OBS_SFE006 = None


def fmt(v, kind):
    v = float(v)
    if kind == "f2":
        return "%.2f" % v
    if kind == "f3":
        return "%.3f" % v
    if kind == "f4":
        return "%.4f" % v
    if kind == "f1":
        return "%.1f" % v
    if kind == "p0":
        return "%.0f" % (100 * v)
    if kind == "p1":
        return "%.1f" % (100 * v)
    if kind == "i":
        return "{:,}".format(int(round(v)))
    if kind == "r10":
        return "{:,}".format(int(round(v / 10.0) * 10))
    if kind == "r100":
        return "{:,}".format(int(round(v / 100.0) * 100))
    raise ValueError(kind)


def add(key, value, kind, source, selection):
    ROWS.append({"id": key, "value": float(value), "text": fmt(value, kind), "format": kind,
                 "source": source, "selection": selection})


def reps(name, cols=None):
    return pd.read_csv(B / name, usecols=cols, low_memory=False)


def band(r):
    r = np.asarray(r, dtype=float)
    return np.median(r), np.percentile(r, 2.5), np.percentile(r, 97.5)


def add_band(key, r, source, selection):
    m, lo, hi = band(r)
    add(key + "_median", m, "f2", source, selection)
    add(key + "_p2.5", lo, "f2", source, selection)
    add(key + "_p97.5", hi, "f2", source, selection)
    return m, lo, hi


# ------------------------------------------------------------------ universe and VIS10 (recomputed)
prep = json.loads((B / "prep_summary.json").read_text(encoding="utf-8"))
ver = prep["verification"]
add("candidate_peptides", ver["candidate_rows"], "i", "prep_summary.json", "verification.candidate_rows")
add("candidate_match_pct", ver["rows_sequence_match"] / ver["candidate_rows"], "p1", "prep_summary.json",
    "rows_sequence_match / candidate_rows")
add("proteins_dropped", ver["proteins_absent_from_fasta"] + ver["proteins_with_any_mismatch"], "i",
    "prep_summary.json", "absent + any mismatch")

u = pd.read_csv(B / "universe_cysteines.csv.gz", low_memory=False)
add("universe_cys", len(u), "i", "universe_cysteines.csv.gz", "rows")
add("universe_prot", u["protein"].nunique(), "i", "universe_cysteines.csv.gz", "unique proteins")
add("universe_detected", u["detected"].sum(), "i", "universe_cysteines.csv.gz", "detected == 1")
add("universe_detected_pct", u["detected"].mean(), "p1", "universe_cysteines.csv.gz", "mean(detected)")
assert len(u) == prep["n_universe_cysteines"] and int(u["detected"].sum()) == prep["n_detected_cysteines"]

VIS10 = prep["vis10_names"]
X = StandardScaler().fit_transform(u[VIS10].to_numpy(dtype=float))
groups = u["protein"].to_numpy()
rng = np.random.default_rng(20260930)  # fold seed of 01_prep_universe.py (provenance.json seeds.cv_fold_seed_prep)
uprot = np.unique(groups)
fold_of = dict(zip(uprot, rng.permutation(len(uprot)) % 5))
fold = np.asarray([fold_of[g] for g in groups])


def cv_auc(t):
    oof = np.zeros(len(t))
    for k in range(5):
        tr, te = fold != k, fold == k
        oof[te] = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs").fit(X[tr], t[tr]).predict_proba(X[te])[:, 1]
    return float(roc_auc_score(t, oof))


def cv_r2(t):
    oof = np.zeros(len(t))
    for k in range(5):
        tr, te = fold != k, fold == k
        oof[te] = LinearRegression().fit(X[tr], t[tr]).predict(X[te])
    return float(r2_score(t, oof))


auc_det = cv_auc(u["detected"].to_numpy())
add("vis10_auc_detection", auc_det, "f3", "universe_cysteines.csv.gz (recomputed; prep_summary.json agrees)",
    "logistic regression on standardized VIS10, 5-fold CV, folds by protein")
assert abs(auc_det - prep["diagnostics"]["vis10_predicts_detection_auc_cv5"]) < 1e-6
ATTR = {"SFE006": "a1_SFE006_KR", "SNO021": "a1b_SNO021_K", "SNO016": "a2_SNO016_KRH",
        "PERS009": "a2b_PERS009_AKRV", "SNO006": "a3_SNO006_DE3", "SFE002": "a3b_SFE002_E"}
aucs = {}
for k, a in ATTR.items():
    aucs[k] = cv_auc(u[a].to_numpy().astype(int))
    stored = [r for r in prep["attributes"] if r["attribute"] == a][0]["vis10_predicts_attribute_auc_cv5"]
    assert abs(aucs[k] - stored) < 1e-6, (k, aucs[k], stored)
    add("vis10_auc_" + k, aucs[k], "f2", "universe_cysteines.csv.gz (recomputed; prep_summary.json agrees)",
        "same procedure, target = attribute")
basic = [aucs[k] for k in ("SFE006", "SNO021", "SNO016", "PERS009")]
acid = [aucs[k] for k in ("SNO006", "SFE002")]
add("vis10_auc_basic_min", min(basic), "f2", "as above", "min over the four K/R attributes")
add("vis10_auc_basic_max", max(basic), "f2", "as above", "max over the four K/R attributes")
add("vis10_auc_acidic_min", min(acid), "f2", "as above", "min over the two acidic attributes")
add("vis10_auc_acidic_max", max(acid), "f2", "as above", "max over the two acidic attributes")

kr = pd.read_csv(B / "universe_krcount.csv.gz", low_memory=False)
assert len(kr) == len(u)
kr = kr.set_index(["protein", "position"]).loc[list(zip(u["protein"], u["position"]))].reset_index() \
    if "position" in kr.columns and "position" in u.columns else kr
kcol = "kr_count20" if "kr_count20" in kr.columns else [c for c in kr.columns if "count20" in c.lower()][0]
kr20 = kr[kcol].to_numpy(dtype=float)
add("kr20_sd", kr20.std(), "f1", "universe_krcount.csv.gz", "population SD of the K/R count within +/-20")
r2 = cv_r2(kr20)
add("vis10_r2_kr20", r2, "f2", "universe_krcount.csv.gz (recomputed; krcount_universe_summary.json agrees)",
    "linear regression on standardized VIS10, same folds")
ks = json.loads((B / "krcount_universe_summary.json").read_text(encoding="utf-8"))
assert abs(r2 - [c for c in ks["columns"] if c["column"] == "kr_count20"][0]["vis10_predicts_r2_cv5"]) < 1e-6
corr = {k: float(np.corrcoef(kr20, u[a].to_numpy(dtype=float))[0, 1]) for k, a in ATTR.items()}
add("r_kr20_basic_min", min(corr[k] for k in ("SFE006", "SNO021", "SNO016", "PERS009")), "f2",
    "universe_krcount.csv.gz + universe_cysteines.csv.gz", "Pearson r, min over K/R attributes")
add("r_kr20_basic_max", max(corr[k] for k in ("SFE006", "SNO021", "SNO016", "PERS009")), "f2",
    "as above", "Pearson r, max over K/R attributes")
add("r_kr20_acidic_absmax", max(abs(corr[k]) for k in ("SNO006", "SFE002")), "f2", "as above",
    "max |r| over acidic attributes")

# ------------------------------------------------------------------ tables without chemistry, observed detection
s0 = reps("reps_S0.csv.gz")
s0a = s0[s0.attribute == "a1_SFE006_KR"]
add("table_npos", s0a.n_pos.mean(), "r10", "reps_S0.csv.gz", "mean positives per table, 500 tables")
add("table_nsites", s0a.n_obs.mean(), "r100", "reps_S0.csv.gz", "mean sites per table")
add("table_nprot", s0a.n_proteins.mean(), "r10", "reps_S0.csv.gz", "mean proteins per table")
add("auc_nochem_observed", s0a.prop_auc_before.mean(), "f2", "reps_S0.csv.gz", "mean propensity AUC")
st = reps("reps_steep.csv.gz")
add("auc_nochem_steeper", st[(st.scenario == "S0steep") & (st.attribute == "a1_SFE006_KR")].prop_auc_before.mean(),
    "f2", "reps_steep.csv.gz", "S0steep mean propensity AUC")

# ------------------------------------------------------------------ loss beyond the sequence-independent attribute
s1 = reps("reps_S1.csv.gz")
perm0 = reps("reps_perm.csv.gz")
perm1 = reps("reps_permgrid.csv.gz")
cells = {
    "observed": (s0, s1[s1.scenario == "S1"]),
    "steeper": (st[st.scenario == "S0steep"], st[st.scenario == "S1steep"]),
    "permuted": (perm0[perm0.scenario == "S0perm"], perm1[perm1.scenario == "S1permgrid"]),
}
LAT = "a4_LATENT30"
EFFS = (0.5, 1.0, 1.5)
loss = {}
ntab = {}
for model, (n0, n1) in cells.items():
    def part(att, eff):
        z0 = n0[n0.attribute == att]
        z1 = n1[(n1.attribute == att) & (n1.planted_attribute == att) & (np.isclose(n1.effect_log2_or, eff))] \
            if "planted_attribute" in n1.columns else n1[(n1.attribute == att) & np.isclose(n1.effect_log2_or, eff)]
        ntab[(model, att, eff)] = (len(z0), len(z1))
        a_ = (z1.matched.mean() - z0.matched.mean()) / eff
        b_ = (z1.matched.mean() - z0.matched.mean()) / (z1.baseline.mean() - z0.baseline.mean())
        return a_, b_
    for k, att in ATTR.items():
        vals = []
        for eff in EFFS:
            a_, b_ = part(att, eff)
            la, lb = part(LAT, eff)
            vals += [1 - a_ / la, 1 - b_ / lb]
        loss[(model, k)] = (min(vals), max(vals))
        add(f"loss_{model}_{k}_min", min(vals), "p0", "reps_S0/S1, reps_steep, reps_perm/permgrid",
            "1 - (genuine-part retention / latent's), normalizations A and B, planted 0.5/1.0/1.5")
        add(f"loss_{model}_{k}_max", max(vals), "p0", "as above", "as above")
    for cls, ks_ in (("offsets", ("SFE006", "SNO021")), ("composition", ("SNO016", "PERS009")),
                     ("acidic", ("SNO006", "SFE002"))):
        add(f"loss_{model}_{cls}_min", min(loss[(model, k)][0] for k in ks_), "p0", "as above", "class min")
        add(f"loss_{model}_{cls}_max", max(loss[(model, k)][1] for k in ks_), "p0", "as above", "class max")
# cross-check against the stored summary
ov = pd.read_csv(B / "overadjustment_summary.csv")
for model, (dm, sc) in {"observed": ("real", "S1"), "steeper": ("steep", "S1steep"),
                        "permuted": ("permuted_within_protein", "S1permgrid")}.items():
    for k, att in ATTR.items():
        z = ov[(ov.detection_model == dm) & (ov.planted_scenario == sc) & (ov.attribute == att)
               & (ov.effect_log2_or.isin(EFFS))]
        v = np.r_[1 - z.matched_genuine_part_retention_relative_to_latent,
                  1 - z.matched_increment_normalised_retention_relative_to_latent]
        assert abs(v.min() - loss[(model, k)][0]) < 1e-6 and abs(v.max() - loss[(model, k)][1]) < 1e-6, (model, k)
# Monte Carlo standard errors of the individual loss estimates (stored delta-method SEs)
for model, (dm, sc) in {"observed": ("real", "S1"), "steeper": ("steep", "S1steep"),
                        "permuted": ("permuted_within_protein", "S1permgrid")}.items():
    z = ov[(ov.detection_model == dm) & (ov.planted_scenario == sc) & (ov.attribute.isin(ATTR.values()))
           & (ov.effect_log2_or.isin(EFFS))]
    se = np.r_[z.matched_genuine_part_retention_relative_to_latent_mcse,
               z.matched_increment_normalised_retention_relative_to_latent_mcse]
    add(f"mcse_{model}_min", se.min(), "p1", "overadjustment_summary.csv", "min stored MC SE, points")
    add(f"mcse_{model}_max", se.max(), "p1", "overadjustment_summary.csv", "max stored MC SE, points")
# table counts per cell
add("ntab_planted_observed", ntab[("observed", ATTR["SFE006"], 1.0)][1], "i", "reps_S1.csv.gz", "tables per cell")
add("ntab_planted_steeper", ntab[("steeper", ATTR["SFE006"], 1.0)][1], "i", "reps_steep.csv.gz", "tables per cell")
add("ntab_planted_permuted", ntab[("permuted", ATTR["SFE006"], 1.0)][1], "i", "reps_permgrid.csv.gz", "per cell")
add("ntab_null_observed", ntab[("observed", ATTR["SFE006"], 1.0)][0], "i", "reps_S0.csv.gz", "tables")
add("ntab_null_steeper", ntab[("steeper", ATTR["SFE006"], 1.0)][0], "i", "reps_steep.csv.gz", "tables")
add("ntab_null_permuted", ntab[("permuted", ATTR["SFE006"], 1.0)][0], "i", "reps_perm.csv.gz", "tables")
# the sequence-independent attribute's own retention of a planted effect (matched / planted), observed detection
lat = s1[(s1.scenario == "S1") & (s1.attribute == LAT)]
lr = [lat[np.isclose(lat.effect_log2_or, e)].matched.mean() / e for e in EFFS]
add("latent_matched_over_truth_min", min(lr), "f2", "reps_S1.csv.gz", "mean matched / planted, 0.5/1.0/1.5")
add("latent_matched_over_truth_max", max(lr), "f2", "reps_S1.csv.gz", "as above")

# ------------------------------------------------------------------ observed re-tests
obs = pd.read_csv(B / "observed_claims.csv")
OBS = {"SFE006": 'SFE-006 (Yang 2014, own data)', "SNO021": "SNO-021", "PERS009": "PERS-009",
       "SNO006": "SNO-006", "SFE002": "SFE-002"}
for k, c in OBS.items():
    r = obs[obs.claim == c].iloc[0]
    add(f"obs_retention_{k}", r.retention_matched_over_baseline, "f2", "observed_claims.csv", c)
    add(f"obs_baseline_{k}", r.baseline, "f2", "observed_claims.csv", c)
    add(f"obs_npos_{k}", r.n_pos, "i", "observed_claims.csv", c)
r = obs[obs.claim == OBS["SFE006"]].iloc[0]
add("obs_baseline_SFE006_3dp", r.baseline, "f3", "observed_claims.csv", "SFE-006 baseline")
aucs_obs = [obs[obs.claim == c].iloc[0].propensity_auc for c in OBS.values()]
add("obs_auc_min", min(aucs_obs), "f2", "observed_claims.csv", "propensity AUC, the five re-tests with a retention")
add("obs_auc_max", max(aucs_obs), "f2", "observed_claims.csv", "as above")
add("obs_npos_SNO014", obs[obs.claim == "SNO-014"].iloc[0].n_pos, "i", "observed_claims.csv", "SNO-014 positives")
# cross-check the observed values against the manuscript's Supplemental Data (read-only)
SUPP = W.parent / "supplemental"  # the repository's Supplemental Data
if SUPP.exists():
    for fname, claims in (("Supplemental_Data_4_retest_round_b.csv", ("SNO-021", "PERS-009")),
                          ("Supplemental_Data_5_retest_round_d.csv", ("SNO-006", "SFE-002"))):
        sd = pd.read_csv(SUPP / fname)
        for c in claims:
            rows = sd[sd["claim_id"].astype(str).str.strip() == c]
            want = obs[obs.claim == c].iloc[0]
            ok = [(abs(float(r.baseline_log2_or) - float(want.baseline)) < 5e-5
                   and abs(float(r.matched_log2_or) - float(want.matched)) < 5e-5) for r in rows.itertuples()]
            print("supplemental cross-check", c, fname, "rows", len(rows), "baseline and matched agree in", sum(ok))
            assert any(ok), c

# ------------------------------------------------------------------ chemistry on the attribute, calibrated
dis = reps("reps_discrim.csv.gz", ["scenario", "attribute", "effect_log2_or", "rate", "n_pos", "baseline", "matched"])
for sc, key in (("DM1real", "observed"), ("DM1steep", "steeper"), ("DM1perm", "permuted"), ("DM1r16", "rate16")):
    z = dis[(dis.scenario == sc) & (dis.attribute == "a1_SFE006_KR")]
    add_band(f"attr_calib_SFE006_{key}", z.matched / z.baseline, "reps_discrim.csv.gz",
             f"{sc}, {len(z)} tables at the interpolated strength")
ovs = pd.read_csv(B / "observed_vs_simulated.csv")
s1g = reps("reps_S1grid.csv.gz", ["scenario", "attribute", "planted_attribute", "effect_log2_or", "baseline", "matched"])
for k in ("SNO021", "PERS009"):
    row = ovs[(ovs.claim == OBS[k]) & (ovs.detection_model == "real")].iloc[0]
    g = row.nearest_grid_effect
    z = s1g[(s1g.attribute == ATTR[k]) & (s1g.planted_attribute == ATTR[k]) & np.isclose(s1g.effect_log2_or, g)]
    m, lo, hi = add_band(f"attr_calib_{k}_observed", z.matched / z.baseline, "reps_S1grid.csv.gz",
                         f"S1grid, nearest grid effect {g}, {len(z)} tables")
    assert abs(m - row.sim_retention_p50_at_nearest) < 1e-9
below = ovs[ovs.claim.isin([OBS["SFE006"], OBS["SNO021"], OBS["PERS009"]])].observed_retention_below_sim_p2_5.all()
assert below, "observed retention not below all attribute-chemistry bands"

# ------------------------------------------------------------------ chemistry on the K/R count within +/-20
pg = reps("reps_proxygrid.csv.gz", ["scenario", "attribute", "effect_log2_or", "baseline", "matched"])
obs_b = float(obs[obs.claim == OBS["SFE006"]].iloc[0].baseline)
for sc, key in (("P1grid", "observed"), ("P1steepgrid", "steeper"), ("P1permgrid", "permuted")):
    z = pg[(pg.scenario == sc) & (pg.attribute == "a1_SFE006_KR")].groupby("effect_log2_or").baseline.mean()
    bstar = float(np.interp(obs_b, z.values, z.index.values)) if np.all(np.diff(z.values) > 0) else np.nan
    add(f"count_b_{key}", bstar, "f2", "reps_proxygrid.csv.gz", f"{sc}: b with mean SFE-006 baseline = observed")
for sc, key in (("DM2real", "observed"), ("DM2steep", "steeper"), ("DM2perm", "permuted"), ("DM2r16", "rate16")):
    z = dis[(dis.scenario == sc) & (dis.attribute == "a1_SFE006_KR")]
    add_band(f"count_calib_SFE006_{key}", z.matched / z.baseline, "reps_discrim.csv.gz",
             f"{sc}, {len(z)} tables at the interpolated strength")
    add(f"count_calib_SFE006_{key}_n", len(z), "i", "reps_discrim.csv.gz", sc)
    if key == "rate16":
        add("rate16_npos", z.n_pos.mean(), "r100", "reps_discrim.csv.gz", "mean positives per table, DM2r16")
retn = [np.median(dis[(dis.scenario == s) & (dis.attribute == "a1_SFE006_KR")].eval("matched / baseline"))
        for s in ("DM2real", "DM2steep", "DM2perm")]
add("count_calib_SFE006_min", min(retn), "f2", "reps_discrim.csv.gz", "min of the three medians")
add("count_calib_SFE006_max", max(retn), "f2", "reps_discrim.csv.gz", "max of the three medians")
# SNO-021: maximum mean baseline on the grids and the band at the strongest strength
for sc, key in (("P1grid", "observed"), ("P1steepgrid", "steeper"), ("P1permgrid", "permuted")):
    z = pg[(pg.scenario == sc) & (pg.attribute == "a1b_SNO021_K")]
    add(f"sno021_count_maxbase_{key}", z.groupby("effect_log2_or").baseline.mean().max(), "f2",
        "reps_proxygrid.csv.gz", f"{sc}: max over b of the mean SNO-021 baseline")
    zz = z[np.isclose(z.effect_log2_or, 6.0)]
    add_band(f"sno021_count_b6_{key}", zz.matched / zz.baseline, "reps_proxygrid.csv.gz", f"{sc}, b = 6, {len(zz)} tables")
    add(f"sno021_count_b6_{key}_n", len(zz), "i", "reps_proxygrid.csv.gz", sc)
# PERS-009
for sc, key in (("P1grid", "observed"), ("P1steepgrid", "steeper"), ("P1permgrid", "permuted")):
    z = pg[(pg.scenario == sc) & (pg.attribute == "a2b_PERS009_AKRV")]
    add(f"pers009_count_maxbase_{key}", z.groupby("effect_log2_or").baseline.mean().max(), "f2",
        "reps_proxygrid.csv.gz", f"{sc}: max over b of the mean PERS-009 baseline")
z = pg[(pg.scenario == "P1permgrid") & (pg.attribute == "a2b_PERS009_AKRV")]
zm = z.groupby("effect_log2_or").baseline.mean()
pobs = float(obs[obs.claim == "PERS-009"].iloc[0].baseline)
i = int(np.argmax(zm.values >= pobs))
bst = float(np.interp(pobs, zm.values[i - 1:i + 1], zm.index.values[i - 1:i + 1]))
add("pers009_count_b_permuted", bst, "f2", "reps_proxygrid.csv.gz", "interpolated b, permuted detection")
zz = z[np.isclose(z.effect_log2_or, 3.0)]
add_band("pers009_count_b3_permuted", zz.matched / zz.baseline, "reps_proxygrid.csv.gz", f"P1permgrid, b = 3, {len(zz)}")
add("pers009_count_b3_permuted_n", len(zz), "i", "reps_proxygrid.csv.gz", "tables")
# acidic: largest rise of the baseline over b = 0 in any model
rise = []
for sc in ("P1grid", "P1steepgrid", "P1permgrid"):
    for att in (ATTR["SNO006"], ATTR["SFE002"]):
        zm = pg[(pg.scenario == sc) & (pg.attribute == att)].groupby("effect_log2_or").baseline.mean()
        rise.append(zm.max() - zm.loc[0.0])
add("acidic_count_rise_max", max(rise), "f2", "reps_proxygrid.csv.gz", "max over models, attributes, b of mean baseline minus b=0")
add("count_grid_bmax", pg.effect_log2_or.max(), "i", "reps_proxygrid.csv.gz", "strongest b simulated")

# ------------------------------------------------------------------ other counts (observed detection)
d2 = reps("reps_discrim2.csv.gz", ["scenario", "attribute", "effect_log2_or", "baseline", "matched"])
z = d2[(d2.scenario == "DM2akrv10") & (d2.attribute == "a1_SFE006_KR")]
add("akrv10_b", z.effect_log2_or.iloc[0], "f2", "reps_discrim2.csv.gz", "DM2akrv10 planted b")
add_band("akrv10_SFE006", z.matched / z.baseline, "reps_discrim2.csv.gz", f"DM2akrv10, {len(z)} tables")
k5 = reps("reps_kr5cal.csv.gz", ["scenario", "attribute", "effect_log2_or", "baseline", "matched"])
z = k5[(k5.scenario == "DM2kr5") & (k5.attribute == "a1_SFE006_KR")]
add("kr5_b", z.effect_log2_or.iloc[0], "f2", "reps_kr5cal.csv.gz", "DM2kr5 planted b")
add_band("kr5_SFE006", z.matched / z.baseline, "reps_kr5cal.csv.gz", f"DM2kr5, {len(z)} tables")

# PERS-009 under the other two counts (observed detection): interpolated b and band at the nearest grid strength
psens = pd.concat([reps("reps_proxysens.csv.gz", ["scenario", "attribute", "effect_log2_or", "baseline", "matched"]),
                   reps("reps_proxysensext.csv.gz", ["scenario", "attribute", "effect_log2_or", "baseline", "matched"])])
pobs = float(obs[obs.claim == "PERS-009"].iloc[0].baseline)
for sc, key in (("P1kr5grid", "kr5"), ("P1akrv10grid", "akrv10")):
    z = psens[(psens.scenario == sc) & (psens.attribute == "a2b_PERS009_AKRV")]
    zm = z.groupby("effect_log2_or").baseline.mean().sort_index()
    i = int(np.argmax(zm.values >= pobs))
    bst = float(np.interp(pobs, zm.values[i - 1:i + 1], zm.index.values[i - 1:i + 1])) if i > 0 else float(zm.index[0])
    grid = np.asarray(zm.index.values, dtype=float)
    near = grid[np.argmin(np.abs(grid - bst))]
    zz = z[np.isclose(z.effect_log2_or, near)]
    add(f"pers009_{key}_b", bst, "f2", "reps_proxysens(ext).csv.gz", f"{sc}: interpolated b for the PERS-009 baseline")
    m, lo, hi = add_band(f"pers009_{key}", zz.matched / zz.baseline, "reps_proxysens(ext).csv.gz",
                         f"{sc}, nearest grid b {near}, {len(zz)} tables")
    assert lo > float(obs[obs.claim == "PERS-009"].iloc[0].retention_matched_over_baseline), key

# ------------------------------------------------------------------ length selection, no chemistry
ls = reps("reps_lenselcal.csv.gz", ["scenario", "attribute", "baseline", "matched", "rate"])
z = ls[(ls.scenario == "L0cal") & (ls.attribute == "a1_SFE006_KR")]
add("lensel_n", len(z), "i", "reps_lenselcal.csv.gz", "L0cal tables")
add("lensel_baseline", z.baseline.mean(), "f3", "reps_lenselcal.csv.gz", "L0cal mean SFE-006 baseline")
add_band("lensel_retention", z.matched / z.baseline, "reps_lenselcal.csv.gz", "L0cal per-table retention")

# ------------------------------------------------------------------ statements the Note makes about these numbers
V = {r["id"]: r["value"] for r in ROWS}
o_sfe, o_sno21, o_pers = V["obs_retention_SFE006"], V["obs_retention_SNO021"], V["obs_retention_PERS009"]
for key in ("observed", "steeper", "permuted"):  # observed SFE-006 retention below every attribute-chemistry range
    assert o_sfe < V[f"attr_calib_SFE006_{key}_p2.5"], key
    assert V[f"count_calib_SFE006_{key}_p2.5"] <= o_sfe <= V[f"count_calib_SFE006_{key}_p97.5"], key  # inside
    assert V[f"sno021_count_b6_{key}_p2.5"] <= o_sno21 <= V[f"sno021_count_b6_{key}_p97.5"], key
    assert V[f"sno021_count_maxbase_{key}"] < V["obs_baseline_SNO021"], key  # SNO-021 baseline not reached
assert min(V[f"sno021_count_b6_{k}_median"] for k in ("observed", "steeper", "permuted")) <= o_sno21 \
    <= max(V[f"sno021_count_b6_{k}_median"] for k in ("observed", "steeper", "permuted"))
assert V["pers009_count_maxbase_observed"] < V["obs_baseline_PERS009"] > V["pers009_count_maxbase_steeper"]
assert V["pers009_count_maxbase_permuted"] >= V["obs_baseline_PERS009"]
assert V["pers009_count_b3_permuted_p2.5"] <= o_pers <= V["pers009_count_b3_permuted_p97.5"]
assert V["kr5_SFE006_median"] < o_sfe < V["akrv10_SFE006_median"]  # "on either side of the observed 0.43"
assert V["kr5_SFE006_p2.5"] <= o_sfe <= V["kr5_SFE006_p97.5"] and V["akrv10_SFE006_p2.5"] <= o_sfe <= V["akrv10_SFE006_p97.5"]
assert V["lensel_retention_p97.5"] < o_sfe  # length selection removed more than the re-test
assert V["obs_auc_min"] > V["auc_nochem_steeper"] > V["auc_nochem_observed"]
print("all statement checks passed")

# ------------------------------------------------------------------ design constants (documented, not results)
DESIGN = [
    ("design_count_window", 20, "i", "phase2_claim_cohorts DIG25 column pep_cleavage_sites_within_20", "+/-20 residues"),
    ("design_positive_rate", 0.20, "p0", "provenance.json settings.positive_rate_among_detected_primary", "20%"),
    ("design_latent_prevalence", 0.30, "p0", "provenance.json settings.latent_prevalence", "30%"),
    ("design_caliper_sd", 0.2, "f1", "provenance.json settings.caliper_sd", "0.2 SD of the logit"),
    ("design_planted_0.5", 0.5, "f1", "overadjustment_summary.csv effect_log2_or", "planted effect"),
    ("design_planted_1.0", 1.0, "f1", "overadjustment_summary.csv effect_log2_or", "planted effect"),
    ("design_planted_1.5", 1.5, "f1", "overadjustment_summary.csv effect_log2_or", "planted effect"),
    ("design_attr_grid_lo", 0, "i", "reps_S1grid.csv.gz effect_log2_or", "grid start"),
    ("design_attr_grid_hi", 3.0, "f1", "reps_S1grid.csv.gz effect_log2_or", "grid end"),
    ("design_attr_grid_step", 0.25, "f2", "reps_S1grid.csv.gz effect_log2_or", "grid step"),
    ("design_count_grid_4", 4, "i", "provenance.json revision_round2.proxy_grid_b", "0 to 4 in steps of 0.5"),
    ("design_count_grid_5", 5, "i", "provenance.json revision_round2.proxy_grid_b", "then 5"),
    ("design_count_grid_step", 0.5, "f1", "provenance.json revision_round2.proxy_grid_b", "step"),
    ("design_pers009_grid_b", 3.0, "f1", "reps_proxygrid.csv.gz", "grid strength nearest to 2.97"),
]
g1 = sorted(s1g.effect_log2_or.unique())
assert g1[0] == 0 and abs(max(e for e in g1 if e <= 3.0) - 3.0) < 1e-9 and 0.25 in g1
for d in DESIGN:
    add(*d)

out = pd.DataFrame(ROWS)
OUT.mkdir(parents=True, exist_ok=True)
out.to_csv(CSV_OUT, index=False, encoding="utf-8")
print("wrote", CSV_OUT, len(out), "rows")
with pd.option_context("display.width", 250, "display.max_rows", 500, "display.max_colwidth", 60):
    print(out[["id", "value", "text", "source"]].to_string())

# ------------------------------------------------------------------ token check of the Note
if NOTE.exists():
    t = NOTE.read_text(encoding="utf-8")
    t = re.sub(r"`[^`]*`", " ", t)                                     # code spans: claim ids, paths, files
    for s in ("K or R at −10, −8 to −6, −4, −2 or +4 to +8", "K at −10, +1, +5, +6 or +10",
              "at least 3 of the ten positions within ±5 are K, R or H",
              "at least 3 of the ten positions within ±5 are A, K, R or V", "D or E at +3",
              "E at −4, −3, +1, +3, +4 or +5"):
        assert s in t, s
        t = t.replace(s, " ")
    for pat in (r"Supplemental Note \d+", r"Supplemental Data \d+ and \d+", r"Table S15\.1", r"UP000000589",
                r"PXD063463", r"2026_03", r"VIS10", r"DIG25", r"log2", r"sha256 [0-9a-f]{8}…", r"\b[0-9a-f]{8}…",
                r"±\d+", r"2\.5th", r"97\.5th", r"1:1", r"R²", r"sha256", r"25-feature", r"(?m)^\d+\. "):
        t = re.sub(pat, " ", t)
    allowed = set(out["text"].astype(str).str.replace("−", "-"))
    bad = []
    for m in re.finditer(r"[−+-]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?", t):
        tok = m.group(0).replace("−", "-").lstrip("+")
        if tok not in allowed:
            bad.append((tok, t[max(0, m.start() - 40):m.end() + 20].replace("\n", " ")))
    print("numeric tokens checked; unexplained:", len(bad))
    for b in bad:
        print("  ", b)
    sys.exit(1 if bad else 0)
