"""G_plmsnosite_rescore, step 4 (system python): build a PROPOSED replacement for
Source_Data_text_detectability_baseline_plmsnosite.csv from this item's outputs.

POST HOC revision analysis (2026-09-30); not registered. The MCP package is not modified; the proposal
is written to results/G_plmsnosite_rescore/proposed_Source_Data_text_detectability_baseline_plmsnosite.csv.

Revised after adversarial verification round 1: long format (one row per number cited in the proposed
Results / Methods / Note 5 text), so that the within-protein and lysine/arginine-stratified numbers the
revised text cites have a source row; the tabulated pLMSNOSite value is described as *matching* the
AUROC of thresholded calls (not "shown to be" it). Registered intervals (seed 20260924) are kept for the
registered VIS10 / DIG25 / label-permuted AUROCs; every other interval comes from the paired
protein-clustered bootstrap (5,000 resamples, seed 20260930) unless stated.

Revised after adversarial verification round 2: the within-protein rows now use the manuscript's estimand
(unweighted mean of per-protein AUROCs over test proteins carrying both labels; step 33) as the primary
figure; the pooled-pair values of round 1 are kept in a block labelled as a sensitivity, next to the
Supplemental Note 12 resampling scheme, a second seed and a descriptive breakdown by protein size.
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

REG = "registered interval: protein-clustered bootstrap, 5,000, seed 20260924 (one generator, DIG25->VIS10->NC)"
PAIR = "paired protein-clustered bootstrap, 5,000, seed 20260930 (resamples shared by all scores)"
S_REG = "registered analysis (2026-09-24), reproduced 2026-09-30 by a re-implementation"
S_PH = "post hoc (revision 2026-09-30)"


def main():
    m = pd.read_csv(C.OUT / "metrics_main.csv").set_index("score")
    rec = pd.read_csv(C.OUT / "recovery_ratios.csv")
    dif = pd.read_csv(C.OUT / "paired_differences.csv")
    thr = pd.read_csv(C.OUT / "plmsnosite_threshold_metrics.csv").set_index("quantity")
    rob = pd.read_csv(C.OUT / "robustness_other.csv").set_index("quantity")
    ter = pd.read_csv(C.OUT / "detectability_learning_tertiles.csv")
    sp = pd.read_csv(C.OUT / "detectability_learning_spearman.csv")
    lr = pd.read_csv(C.OUT / "detectability_learning_logistic.csv")
    kr = pd.read_csv(C.OUT / "kr_stratified_summary.csv").set_index("statistic")
    rep = json.loads((C.OUT / "detectability_reproduction.json").read_text(encoding="utf-8"))["results"]
    diag = json.loads((C.OUT / "reproduction_diagnostics.json").read_text(encoding="utf-8"))
    npc = json.loads((C.OUT / "numpy_forward_pass_check.json").read_text(encoding="utf-8"))
    shared = diag["H1_shared_generator"]["DIG25->VIS10->NC"]
    n_all, n_pos, n_prot = int(m.loc["pLMSNOSite", "n_sites"]), int(m.loc["pLMSNOSite", "n_pos"]), int(m.loc["pLMSNOSite", "n_proteins"])

    rows = []

    def add(block, quantity, model, value, lo=None, hi=None, method="", n=None, npos=None, nprot=None, status=S_PH,
            source=""):
        rows.append({"block": block, "quantity": quantity, "model": model, "value": value, "ci_low": lo, "ci_high": hi,
                     "interval_method": method, "n_sites": n, "n_positives": npos, "n_proteins": nprot,
                     "status": status, "source_file": source})

    # --- registered detectability baselines
    for det in ("VIS10", "DIG25"):
        add("registered baseline", "test AUROC", det, rep[det]["auroc_sklearn"], shared[det]["ci"][0], shared[det]["ci"][1],
            REG, n_all, n_pos, n_prot, S_REG, "detectability_reproduction.json; reproduction_diagnostics.json")
        add("registered baseline", "test AUPRC", det, rep[det]["auprc_sklearn"], m.loc[det, "auprc_ci_low"],
            m.loc[det, "auprc_ci_high"], PAIR + " (the registered analysis gave no AUPRC interval)", n_all, n_pos, n_prot,
            S_REG, "detectability_reproduction.json; metrics_main.csv")
        add("registered baseline", "test AUROC, interval from the shared post hoc resamples", det, m.loc[det, "auroc"],
            m.loc[det, "auroc_ci_low"], m.loc[det, "auroc_ci_high"], PAIR, n_all, n_pos, n_prot, S_PH, "metrics_main.csv")
    add("registered baseline", "test AUROC, label-permuted negative control", "DIG25 label-permuted",
        rep["DIG25_label_permuted"]["auroc"], shared["NC"]["ci"][0], shared["NC"]["ci"][1], REG, n_all, n_pos, n_prot,
        S_REG, "detectability_reproduction.json; reproduction_diagnostics.json")
    add("registered baseline", "test prevalence", "", m.loc["pLMSNOSite", "auprc_chance"], n=n_all, npos=n_pos,
        nprot=n_prot, status="data", source="metrics_main.csv")
    for i, v in enumerate(C.MANUSCRIPT["PC_folds"], 1):
        add("registered baseline", f"positive control, grouped 5-fold CV on training, fold {i} (as registered)", "DIG25", v,
            status="registered value as reported in Note 5; fold assignment not recorded; not reproduced",
            source="Supplemental Note 5")
    add("registered baseline", "positive control, mean of 5 folds (as registered)", "DIG25", C.MANUSCRIPT["PC_mean"],
        status="registered value as reported in Note 5; not reproduced", source="Supplemental Note 5")
    gk = rep["DIG25_groupkfold5_train"]
    for i, v in enumerate(gk["fold_auroc"], 1):
        add("registered baseline", f"positive control re-run, unshuffled protein-grouped 5-fold CV, fold {i}", "DIG25", v,
            status=S_PH, source="detectability_reproduction.json")
    add("registered baseline", "positive control re-run, mean of 5 folds", "DIG25", gk["mean"], status=S_PH,
        source="detectability_reproduction.json")

    # --- re-scored predictor
    for key, label in (("pLMSNOSite", "pLMSNOSite released models, probabilities"),
                       ("pLMSNOSite_call_0.5", "pLMSNOSite released models, calls thresholded at 0.5"),
                       ("ProtT5_arm", "pLMSNOSite ProtT5 arm alone"), ("Embedding_arm", "pLMSNOSite embedding arm alone"),
                       ("rankavg_pLMSNOSite_DIG25", "unfitted rank average of pLMSNOSite and DIG25")):
        for met in ("auroc", "auprc"):
            add("re-scored predictor", f"test {met.upper()}", label, m.loc[key, met], m.loc[key, f"{met}_ci_low"],
                m.loc[key, f"{met}_ci_high"], PAIR, n_all, n_pos, n_prot, S_PH, "metrics_main.csv")
    for q, lab in (("sn", "sensitivity at 0.5"), ("sp", "specificity at 0.5"),
                   ("ba", "(Sn+Sp)/2 at 0.5 = AUROC of thresholded calls = 'AUROC' printed by evaluate_model.py"),
                   ("gmean", "geometric mean of Sn and Sp at 0.5"), ("mcc", "MCC at 0.5")):
        add("re-scored predictor", lab, "pLMSNOSite released models", thr.loc[q, "value"], thr.loc[q, "ci_low"],
            thr.loc[q, "ci_high"], PAIR, n_all, n_pos, n_prot, S_PH, "plmsnosite_threshold_metrics.csv")
    for q in ("tp", "fn", "tn", "fp"):
        add("re-scored predictor", f"confusion matrix count {q.upper()} at 0.5", "pLMSNOSite released models",
            thr.loc[q, "value"], status=S_PH, source="plmsnosite_threshold_metrics.csv; evaluate_model_verbatim_stdout.txt")
    add("re-scored predictor", "max |numpy - TensorFlow| probability (TensorFlow-free re-computation)",
        "pLMSNOSite released models", npc["max_abs_diff_vs_tensorflow"]["stacked_probability"], status=S_PH,
        source="numpy_forward_pass_check.json")

    # --- paired differences
    for comp in ("pLMSNOSite - VIS10", "pLMSNOSite - DIG25", "rankavg_pLMSNOSite_DIG25 - pLMSNOSite", "DIG25 - VIS10"):
        for met in ("auroc", "auprc"):
            r = dif[(dif.comparison == comp) & (dif.metric == met)].iloc[0]
            add("paired difference", f"{met.upper()} difference", comp.replace("rankavg_pLMSNOSite_DIG25", "rank average"),
                r.difference, r.ci_low, r.ci_high, PAIR, n_all, n_pos, n_prot, S_PH, "paired_differences.csv")

    # --- recovery ratios
    for det in ("VIS10", "DIG25"):
        for prefix, lab, st in (
                ("pLMSNOSite re-scored score-based AUROC", "recovered share of above-chance AUROC, vs re-scored pLMSNOSite", S_PH),
                ("pLMSNOSite re-scored AUPRC", "recovered share of AUPRC above prevalence, vs re-scored pLMSNOSite", S_PH),
                ("pLMSNOSite tabulated", "recovered share of above-chance AUROC, vs tabulated pLMSNOSite 0.754 (registered denominator)",
                 "registered statistic, recomputed at full precision"),
                ("PreSNO tabulated", "recovered share of above-chance AUROC, vs tabulated PreSNO 0.756 (quoted; superseded)",
                 "first-version statistic, recomputed; not used"),
                ("DeepNitro tabulated", "recovered share of above-chance AUROC, vs tabulated DeepNitro 0.731 (quoted; superseded)",
                 "first-version statistic, recomputed; not used")):
            r = rec[(rec.detectability_model == det) & (rec.denominator.str.startswith(prefix))].iloc[0]
            meth = PAIR if "re-scored" in prefix else "numerator-only bootstrap (denominator a fixed quoted constant), seed 20260930"
            add("recovery ratio", lab, det, r.recovery, r.ci_low, r.ci_high, meth, n_all, n_pos, n_prot, st, "recovery_ratios.csv")

    # --- within-protein AUROC (revision round 2): the manuscript's estimand, the unweighted mean of per-protein
    #     AUROCs over test proteins carrying both labels (repo protein_rows -> within_protein_auc_mean), is primary;
    #     the pooled-pair value of round 1 is kept only as a labelled sensitivity.
    wpm = pd.read_csv(C.OUT / "withinprot_mean_summary.csv")
    wjs = json.loads((C.OUT / "withinprot_mean.json").read_text(encoding="utf-8"))
    dg = wjs["diagnostics"]
    nb, sb, pb = dg["n_proteins_with_both_labels"], dg["n_sites_in_them"], dg["n_positives_in_them"]
    wp_mean = wpm.estimand.str.startswith("mean of per-protein AUROCs")
    wp_pool = wpm.estimand.str.startswith("pooled within-protein pairs")
    sch_main = wpm.scheme.str.startswith("all 278 test proteins resampled, seed 20260930")
    sch_n12 = wpm.scheme.str.startswith("only proteins with both labels resampled")
    sch_alt = wpm.scheme.str.startswith("all 278 test proteins resampled, seed 20260924")
    PAIR_WP = (PAIR + "; per replicate, multiplicity-weighted mean of the per-protein AUROCs of the resampled proteins "
               "that carry both labels")
    for sel, block, meth, st, tag in (
            (wp_mean & sch_main, "within-protein (manuscript estimand, primary)", PAIR_WP, S_PH,
             "mean of per-protein AUROCs"),
            (wp_mean & sch_n12, "within-protein, sensitivity: Supplemental Note 12 resampling scheme",
             "paired bootstrap over the 220 proteins with both labels only, 5,000, seed 20260930", S_PH + "; sensitivity",
             "mean of per-protein AUROCs"),
            (wp_mean & sch_alt, "within-protein, sensitivity: second bootstrap seed",
             PAIR_WP.replace("20260930", "20260924"), S_PH + "; Monte Carlo check", "mean of per-protein AUROCs"),
            (wp_pool & sch_main, "within-protein, sensitivity: pooled within-protein pairs (round-1 estimand)", PAIR,
             S_PH + "; sensitivity (weights proteins by their number of positive-negative pairs)", "pooled pairs")):
        for _, r in wpm[sel].iterrows():
            if r.quantity.endswith("(secondary score)") or r.quantity.endswith("(secondary)"):
                continue
            add(block, f"{r.quantity} [{tag}]", r.model, r.value, r.ci_low, r.ci_high, meth, sb, pb, nb, st,
                "withinprot_mean_summary.csv")
    ppool = wjs["point_pooled_within_protein_auroc"]
    emb_note = (f"; descriptive, not in the pre-listed comparisons; with pooled pairs the embedding arm gives "
                f"{ppool['Embedding_arm']:.3f} against DIG25 {ppool['DIG25']:.3f}")
    for _, r in wpm[wp_mean & sch_main & wpm.model.isin(["Embedding_arm", "Embedding_arm - DIG25",
                                                          "Embedding_arm - pLMSNOSite"])].iterrows():
        add("within-protein, secondary (descriptive)", f"{r.quantity} [mean of per-protein AUROCs]",
            r.model.replace("Embedding_arm", "pLMSNOSite embedding arm alone"), r.value, r.ci_low, r.ci_high, PAIR_WP,
            sb, pb, nb, S_PH + emb_note, "withinprot_mean_summary.csv")
    size = pd.read_csv(C.OUT / "withinprot_mean_by_size.csv")
    for _, r in size.iterrows():
        g = f"pair-count tertile {int(r.pair_count_tertile)} ({r.pair_count_range} within-protein pairs per protein)"
        add("within-protein, by protein size (descriptive)", f"share of all within-protein pairs, {g}", "",
            r.share_of_all_within_protein_pairs, n=int(r.n_sites), npos=int(r.n_pos), nprot=int(r.n_proteins),
            status=S_PH + "; descriptive", source="withinprot_mean_by_size.csv")
        for k, lab in (("pLMSNOSite", "pLMSNOSite released models"), ("VIS10", "VIS10"), ("DIG25", "DIG25")):
            add("within-protein, by protein size (descriptive)", f"mean of per-protein AUROCs, {g}", lab,
                r[f"mean_per_protein_auroc_{k}"], n=int(r.n_sites), npos=int(r.n_pos), nprot=int(r.n_proteins),
                status=S_PH + "; descriptive", source="withinprot_mean_by_size.csv")
    # --- robustness: all flagged sites excluded (cited in the proposed Note; added in revision round 2, the round-1
    #     file lacked these rows although the round-1 Note cited them)
    rbs = pd.read_csv(C.OUT / "robustness_subsets.csv").set_index("subset")
    wps = pd.read_csv(C.OUT / "withinprot_mean_subsets.csv").set_index("subset")
    ex = rbs.loc["excl. all of the above"]
    exw = wps.loc["excl. all of the above"]
    blk = "robustness: all flagged test sites excluded"
    nx, px, kx = int(ex.n_sites), int(ex.n_pos), int(ex.n_clusters)
    add(blk, "test AUROC", "pLMSNOSite released models", ex.auroc_pLMSNOSite, ex.auroc_pLMSNOSite_ci_low,
        ex.auroc_pLMSNOSite_ci_high, PAIR.replace("(resamples shared by all scores)", "over the retained proteins"),
        nx, px, kx, S_PH, "robustness_subsets.csv")
    for k in ("VIS10", "DIG25"):
        add(blk, "AUROC difference", f"pLMSNOSite - {k}", ex[f"auroc_diff_pLMSNOSite_minus_{k}"],
            ex[f"auroc_diff_{k}_ci_low"], ex[f"auroc_diff_{k}_ci_high"], PAIR, nx, px, kx, S_PH, "robustness_subsets.csv")
        add(blk, "recovered share of above-chance AUROC, vs re-scored pLMSNOSite", k, ex[f"recovery_{k}"],
            ex[f"recovery_{k}_ci_low"], ex[f"recovery_{k}_ci_high"], PAIR, nx, px, kx, S_PH, "robustness_subsets.csv")
        add(blk, "within-protein AUROC difference [mean of per-protein AUROCs]", f"pLMSNOSite - {k}",
            exw[f"diff_pLMSNOSite_minus_{k}"], exw[f"diff_{k}_ci_low"], exw[f"diff_{k}_ci_high"], PAIR_WP,
            int(exw.n_sites_in_them), int(exw.n_pos_in_them), int(exw.n_proteins_both_labels), S_PH,
            "withinprot_mean_subsets.csv")
        add(blk, "recovered share of pLMSNOSite's above-chance within-protein AUROC [mean of per-protein AUROCs]", k,
            exw[f"recovery_{k}"], exw[f"recovery_{k}_ci_low"], exw[f"recovery_{k}_ci_high"], PAIR_WP,
            int(exw.n_sites_in_them), int(exw.n_pos_in_them), int(exw.n_proteins_both_labels), S_PH,
            "withinprot_mean_subsets.csv")

    # --- base rates of the detectability indicators (cited in the proposed Note)
    trk = pd.read_csv(C.OUT / "detectability_tracking.csv")
    for _, r in trk[trk.analysis == "base rate"].iterrows():
        lab = ("share on a peptide detectable with <=2 missed cleavages" if r.indicator.endswith("any_missed_cleavage")
               else "share on a detectable fully cleaved peptide (7-30 residues, 700-3,500 Da)")
        add("base rates", f"{lab} ({int(r.n_indicator_1)} of {int(r.n)})", r["class"], r.value, n=int(r.n),
            npos=int(r.n) if r["class"] == "positives" else 0, status="data", source="detectability_tracking.csv")

    # the pooled values re-derived in step 33 equal step 30's (asserted there); keep step 30's file as a cross-reference
    assert abs(rob.loc["pooled within-protein AUROC, pLMSNOSite", "value"]
               - wpm[wp_pool & sch_main & (wpm.model == "pLMSNOSite")].value.iloc[0]) < 1e-9

    # --- negatives only: score vs K/R count and detectability
    nneg = n_all - n_pos
    krmap = [
        ("N1|spearman_score_kr|negatives|pLMSNOSite", "Spearman with the number of K and R within 20 residues, negatives", "pLMSNOSite"),
        ("N1|spearman_score_kr|negatives|Embedding_arm", "Spearman with the number of K and R within 20 residues, negatives", "embedding arm"),
        ("N1|spearman_score_kr|negatives|ProtT5_arm", "Spearman with the number of K and R within 20 residues, negatives", "ProtT5 arm"),
        ("N1|spearman_indicator_kr|negatives|le2_missed_cleavages",
         "Spearman of 'peptide detectable with <=2 missed cleavages' with the K/R count, negatives", "indicator"),
        ("N3|spearman_pLMSNOSite_DIG25_crude|negatives", "Spearman of pLMSNOSite with DIG25 score, negatives", "pLMSNOSite ~ DIG25"),
        ("N3|spearman_pLMSNOSite_DIG25_withinKR|negatives", "same, pooled within K/R-count strata", "pLMSNOSite ~ DIG25"),
        ("N3|spearman_pLMSNOSite_VIS10_crude|negatives", "Spearman of pLMSNOSite with VIS10 score, negatives", "pLMSNOSite ~ VIS10"),
        ("N3|spearman_pLMSNOSite_VIS10_withinKR|negatives", "same, pooled within K/R-count strata", "pLMSNOSite ~ VIS10"),
        ("N2|auroc_indicator_crude|negatives|le2_missed_cleavages|pLMSNOSite",
         "AUROC for 'on a peptide detectable with <=2 missed cleavages', negatives", "pLMSNOSite"),
        ("N2|auroc_indicator_KRstrat|negatives|le2_missed_cleavages|pLMSNOSite", "same, K/R-count-stratified", "pLMSNOSite"),
        ("N2|auroc_indicator_crude|negatives|le2_missed_cleavages|Embedding_arm",
         "AUROC for 'on a peptide detectable with <=2 missed cleavages', negatives", "embedding arm"),
        ("N2|auroc_indicator_KRstrat|negatives|le2_missed_cleavages|Embedding_arm", "same, K/R-count-stratified", "embedding arm"),
        ("N2|auroc_indicator_crude|negatives|le2_missed_cleavages|ProtT5_arm",
         "AUROC for 'on a peptide detectable with <=2 missed cleavages', negatives", "ProtT5 arm"),
        ("N2|auroc_indicator_KRstrat|negatives|le2_missed_cleavages|ProtT5_arm", "same, K/R-count-stratified", "ProtT5 arm"),
        ("N2|auroc_indicator_crude|negatives|fully_cleaved|pLMSNOSite",
         "AUROC for 'on a detectable fully cleaved peptide', negatives", "pLMSNOSite"),
        ("N2|auroc_indicator_KRstrat|negatives|fully_cleaved|pLMSNOSite", "same, K/R-count-stratified", "pLMSNOSite"),
        ("N4|spearman_score_kr|negatives_indicator1|pLMSNOSite",
         "Spearman with the K/R count, negatives on a peptide detectable with <=2 missed cleavages", "pLMSNOSite"),
        ("N4|spearman_score_kr|negatives_indicator0|pLMSNOSite",
         "Spearman with the K/R count, negatives not on such a peptide", "pLMSNOSite"),
    ]
    for key, lab, model in krmap:
        r = kr.loc[key]
        add("negatives: detectability and K/R count", lab, model, r.value, r.ci_low, r.ci_high, PAIR, nneg, 0, None, S_PH,
            "kr_stratified_summary.csv")
    for key, lab, model in (
            ("L1|label_auroc_KRstrat|all|pLMSNOSite", "label AUROC, K/R-count-stratified", "pLMSNOSite released models"),
            ("L1|label_auroc_KRstrat|all|DIG25", "label AUROC, K/R-count-stratified", "DIG25"),
            ("L1|label_auroc_KRstrat|all|VIS10", "label AUROC, K/R-count-stratified", "VIS10"),
            ("L1|recovery_VIS10_vs_pLMSNOSite_KRstrat|all", "recovered share of above-chance AUROC, K/R-count-stratified", "VIS10"),
            ("L1|recovery_DIG25_vs_pLMSNOSite_KRstrat|all", "recovered share of above-chance AUROC, K/R-count-stratified", "DIG25"),
            ("L1|label_auroc_diff_pLMSNOSite_minus_VIS10_KRstrat|all", "label AUROC difference, K/R-count-stratified", "pLMSNOSite - VIS10"),
            ("L1|label_auroc_diff_pLMSNOSite_minus_DIG25_KRstrat|all", "label AUROC difference, K/R-count-stratified", "pLMSNOSite - DIG25")):
        r = kr.loc[key]
        add("K/R-count-stratified label AUROC", lab, model, r.value, r.ci_low, r.ci_high, PAIR, n_all, n_pos, n_prot, S_PH,
            "kr_stratified_summary.csv")
    for _, r in sp.iterrows():
        if r["class"] != "positives":  # negatives are reported above (identical values, same resamples)
            continue
        add("positives: score agreement", "Spearman, positives only", r["pair"], r.spearman_rho, r.ci_low, r.ci_high,
            PAIR, int(r.n), int(r.n), None, S_PH, "detectability_learning_spearman.csv")

    # --- tertiles and logistic (descriptive)
    for _, r in ter[(ter.stratifier == "DIG25_tertile") & (ter.score == "pLMSNOSite")].iterrows():
        add("DIG25-score tertiles", f"AUROC within DIG25 tertile {int(r.tertile)} (cut points {r.cut_points})",
            "pLMSNOSite released models", r.auroc, r.ci_low, r.ci_high, PAIR, int(r.n_sites), int(r.n_pos),
            int(r.n_proteins), S_PH, "detectability_learning_tertiles.csv")
    for _, r in lr.iterrows():
        if r.model in ("M1: pLMSNOSite", "M3: pLMSNOSite + DIG25"):
            add("logistic regression (test set, descriptive)", f"odds ratio per SD of {r.term}, {r.model}",
                r.term.replace("logit_", ""), r.or_per_sd, r.or_per_sd_ci_low, r.or_per_sd_ci_high,
                "protein-cluster-robust (statsmodels), not bootstrap", n_all, n_pos, n_prot, S_PH,
                "detectability_learning_logistic.csv")

    # --- tabulated values from the pLMSNOSite paper
    add("tabulated (quoted)", "test AUROC tabulated in the pLMSNOSite paper", "pLMSNOSite", C.MANUSCRIPT["published_auroc_pLMSNOSite"],
        status="quoted; matches the 'AUROC' that the released evaluate_model.py computes from calls thresholded at 0.5 "
               "(0.7537 on re-running); the paper's specificity 0.773 differs from the re-run's 0.772 (2,447/3,168; "
               "0.773 needs 1 to 3 more true negatives)",
        source="Pratyush et al. 2023; evaluate_model_verbatim_stdout.txt")
    for name in ("PreSNO", "DeepNitro"):
        add("tabulated (quoted)", "test AUROC tabulated in the pLMSNOSite paper", name, C.MANUSCRIPT[f"published_auroc_{name}"],
            status="quoted; not recomputable (no released predictions), computation unverified; not used",
            source="Pratyush et al. 2023")

    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / "proposed_Source_Data_text_detectability_baseline_plmsnosite.csv", index=False, float_format="%.10g")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    print(out[["block", "quantity", "model", "value", "ci_low", "ci_high"]].round(4).to_string())
    print(len(out), "rows")


if __name__ == "__main__":
    main()
