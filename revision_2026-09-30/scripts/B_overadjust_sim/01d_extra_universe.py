# -*- coding: utf-8 -*-
"""B_overadjust_sim step 1d (POST HOC, added in the revision after adversarial verification).

Extra per-cysteine columns for the universe of step 1, written to a SEPARATE file so that the stored
universe (and every stored replicate) stays unchanged:

  a1c_KRH_p5p8     K, R or H at any of +5..+8 (repository function phase2_claim_cohorts.flank_residue_flag,
                   imported verbatim). Exploratory Artifact-1 analogue: the manuscript describes the Artifact-1
                   enrichment as K/R/H residues enriched at offsets of roughly +5 to +8 (Results, Artifact 1),
                   but its exact binary feature could not be recovered from the repository.
  m3nc_best_logit  per-cysteine logit of the best peptide-level detection probability from the step-1
                   peptide model refitted WITHOUT the two composition terms (frac_basic, frac_acidic); used by
                   the 'steep_nocomp' selection model (positives selected on a detectability score that
                   contains no peptide composition term).
  vis_det_logit    logit of the VIS10 detection model (detected ~ VIS10 over the universe; the repository's
                   rt.propensity_score, identical to sim_engine's det_logit); used by the 'steep_vis'
                   selection model, in which positives are selected on a score that VIS10 represents exactly.

The candidate table, protein filter and model settings are those of 01_prep_universe.py (same code path);
the refitted model's coefficients are written to results/B_overadjust_sim/extra_universe_summary.json.
Writes results/B_overadjust_sim/universe_extra.csv.gz (same row order as universe_cysteines.csv.gz).
"""
from __future__ import annotations

import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import phase2_claim_cohorts as p2c  # noqa: E402  (repo, read-only)
import phase2b_claim_cohorts as p2b  # noqa: E402
import run_phase2_claims_under_detectability_control as rt  # noqa: E402

OUT = os.path.join(cb.RESULTS, "universe_extra.csv.gz")
OUT_JSON = os.path.join(cb.RESULTS, "extra_universe_summary.json")
KRH_P5P8_OFFSETS = [5, 6, 7, 8]
M3_FEATS_FULL = ["length", "missed_cleavages", "mass", "gravy", "frac_basic", "frac_acidic",
                 "dist_n_term", "dist_c_term", "n_cys"]
M3_FEATS_NOCOMP = [f for f in M3_FEATS_FULL if f not in ("frac_basic", "frac_acidic")]
VIS10_NAMES = ["pep_len", "pep_mass", "pep_gravy", "pep_detectable_length", "pep_detectable_mass",
               "pep_detectable_both", "pep_log_len", "pep_mc1_detectable", "pep_mc2_detectable",
               "pep_detectable_any_missed_cleavage"]


def best_logit_per_cys(kept_rows, p_cand, keys):
    best_p = {}
    for acc, st, pep, pp in zip(kept_rows["protein"], kept_rows["start"], kept_rows["sequence"], p_cand):
        for i, ch in enumerate(pep):
            if ch == "C":
                k = (acc, int(st) + i)
                if pp > best_p.get(k, -1.0):
                    best_p[k] = float(pp)
    bp = np.asarray([best_p.get(k, np.nan) for k in keys], dtype=float)
    bp = np.clip(bp, 1e-6, 1 - 1e-6)
    return np.log(bp / (1 - bp))


def main():
    u = pd.read_csv(cb.UNIVERSE, dtype={"protein": str})
    keys = list(zip(u["protein"], u["position"].astype(int)))
    fasta = p2b.read_fasta_gz(cb.FASTA)
    seq_of = {a: fasta[a] for a in u["protein"].unique()}

    # 1. exploratory Artifact-1 analogue attribute (repository function, verbatim)
    a1c = p2c.flank_residue_flag(seq_of, keys, KRH_P5P8_OFFSETS, {"K", "R", "H"})

    # 2. composition-free peptide detection model (same candidate filter as step 1)
    cand = pd.read_csv(cb.CANDIDATES, dtype={"protein": str, "sequence": str})
    sv = cand["protein"].map(fasta).to_numpy()
    starts, ends, peps = cand["start"].to_numpy(), cand["end"].to_numpy(), cand["sequence"].to_numpy()
    ok = np.zeros(len(cand), dtype=bool)
    for i in range(len(cand)):
        s = sv[i]
        if isinstance(s, str):
            ok[i] = s[starts[i] - 1:ends[i]] == peps[i]
    cand["fasta_match"] = ok
    per_prot = cand.groupby("protein")["fasta_match"].agg(["size", "sum"])
    keep_prot = {p for p, r in per_prot.iterrows() if p in fasta and r["sum"] == r["size"]}
    kept = cand[cand["protein"].isin(keep_prot)]
    y = kept["empirical_detected"].to_numpy()
    fits = {}
    logits = {}
    for name, feats in (("full", M3_FEATS_FULL), ("nocomp", M3_FEATS_NOCOMP)):
        X = StandardScaler().fit_transform(kept[feats].to_numpy(dtype=float))
        m = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs").fit(X, y)
        p = m.predict_proba(X)[:, 1]
        fits[name] = {"features": feats, "candidate_level_auc_in_sample": float(roc_auc_score(y, p)),
                      "coefficients_standardised": dict(zip(feats, [float(c) for c in m.coef_[0]]))}
        logits[name] = best_logit_per_cys(kept, p, keys)
    # the refitted full model must reproduce the stored m3_best_logit (determinism check)
    full_diff = float(np.nanmax(np.abs(logits["full"] - u["m3_best_logit"].to_numpy(dtype=float))))

    # 3. VIS10 detection-model logit (as sim_engine computes det_logit)
    vis = u[VIS10_NAMES].to_numpy(dtype=float)
    det = u["detected"].to_numpy().astype(int)
    vis_logit = rt.propensity_score(vis, det)[1]

    out = pd.DataFrame({"protein": u["protein"], "position": u["position"], "a1c_KRH_p5p8": a1c,
                        "m3nc_best_logit": logits["nocomp"], "vis_det_logit": vis_logit})
    out.to_csv(OUT, index=False, compression="gzip")

    dmask = det == 1
    def corr_det(v):
        return float(np.corrcoef(v[dmask], u["m3_best_logit"].to_numpy(dtype=float)[dmask])[0, 1])
    summ = {
        "item": cb.ITEM, "label": "POST HOC (revision after adversarial verification); not registered",
        "a1c_KRH_p5p8": {"definition": "K, R or H at any of +5..+8 (phase2_claim_cohorts.flank_residue_flag)",
                         "prevalence_all": float(a1c.mean()), "prevalence_detected": float(a1c[dmask].mean()),
                         "prevalence_undetected": float(a1c[~dmask].mean())},
        "peptide_models": fits,
        "full_model_refit_max_abs_diff_vs_stored_m3_best_logit": full_diff,
        "corr_among_detected_with_m3_best_logit": {"m3nc_best_logit": corr_det(logits["nocomp"]),
                                                    "vis_det_logit": corr_det(vis_logit)},
        "n_rows": int(len(out)),
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
