"""G_plmsnosite_rescore, step 3c (system python): within-protein AUROC under the manuscript's estimand.

POST HOC revision analysis (2026-09-30, revision round 2); not registered, not pre-specified.

Why this step exists. Verification round 2 (problem 1, major) found that round 1 reported a POOLED
within-protein AUROC (all within-protein positive-negative pairs pooled over proteins, i.e. proteins
weighted by their pair count), whereas the manuscript's "within-protein AUC" (Results, ranker refit;
Supplemental Notes 11 and 12; repository scripts/common.py protein_rows -> within_protein_auc_mean)
is the UNWEIGHTED MEAN of per-protein AUROCs over the proteins that carry both labels. This step makes
the manuscript's estimand the primary within-protein figure and keeps the pooled-pair value only as a
labelled sensitivity.

Per-protein AUROCs are computed by the repository's own `protein_rows` (repo/scripts/common.py,
imported by path with bytecode writing disabled; nothing is written to the repository), and checked
against an independent tie-corrected Mann-Whitney count.

Bootstrap, primary. The SAME 5,000 protein-clustered resamples as step 30 (all 278 test proteins,
sorted UniProt accessions, numpy default_rng(20260930), one rng.integers(0, 278, 278) per replicate).
In each replicate the estimand is the multiplicity-weighted mean of the per-protein AUROCs of the
resampled proteins that carry both labels (a protein drawn twice counts twice, exactly as two distinct
proteins would). Models are not refit. Identity with step 30's resamples is asserted against the
stored replicate file (pooled within-protein and across-protein AUROC replicates reproduced).

Sensitivities (labelled as such in every output):
  (a) the Supplemental Note 12 scheme: resample only the proteins that carry both labels
      (rng.integers(0, 220, 220), seed 20260930), as in repo scripts/analyse_public_refit_2026-09-21.py;
  (b) the pooled-pair estimand of round 1 (pair-weighted mean of the same per-protein AUROCs, identical
      to step 30's values);
  (c) a second bootstrap seed (20260924) for the primary scheme;
  (d) the data-quality subsets and identical-sequence clustering of step 30.
Descriptive: pair-count distribution, within-protein tie shares, per-protein wins, and the mean
per-protein AUROC by protein size (tertiles of the number of within-protein pairs).
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

import sys

sys.dont_write_bytecode = True  # never write __pycache__ into the read-only repository

import importlib.util  # noqa: E402
import json  # noqa: E402
import pathlib  # noqa: E402
import platform  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy  # noqa: E402
import sklearn  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

SEED = 20260930
SEED_ALT = 20260924
REPS = 5000
REPO_COMMON = pathlib.Path(r"C:/Users/admin/Desktop/小论文/_cys_repo_work/repo/scripts/common.py")

PRIMARY = ("pLMSNOSite", "VIS10", "DIG25")
SECONDARY = ("ProtT5_arm", "Embedding_arm", "rankavg_pLMSNOSite_DIG25", "rankavg_pLMSNOSite_VIS10")
ALL_SCORES = PRIMARY + SECONDARY
EST_MEAN = "mean of per-protein AUROCs (manuscript estimand; primary)"
EST_POOL = "pooled within-protein pairs (round-1 estimand; sensitivity)"
SCH_MAIN = "all 278 test proteins resampled, seed 20260930, resamples shared with step 30 (primary)"
SCH_N12 = "only proteins with both labels resampled, seed 20260930 (Supplemental Note 12 scheme; sensitivity)"
SCH_ALT = "all 278 test proteins resampled, seed 20260924 (Monte Carlo check)"


def load_repo_protein_rows():
    spec = importlib.util.spec_from_file_location("repo_common_protein_rows", REPO_COMMON)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.protein_rows


def mann_whitney_auc(pos, neg):
    d = pos[:, None] - neg[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size), int((d == 0).sum()), int(d.size)


def per_protein_table(df, scores, protein_rows):
    """Per-protein AUROC of every score over proteins carrying both labels (repo protein_rows),
    with pair counts, tie counts and an independent Mann-Whitney check."""
    y = df["Target"].to_numpy().astype(int)
    prot = df["UniProt"].to_numpy().astype(str)
    base = None
    cols = {}
    for k in scores:
        s = df[k].to_numpy(dtype=float)
        rows = protein_rows(y, s, prot, prot)
        acc = [r["accession"] for r in rows]
        if base is None:
            base = acc
        assert acc == base, "protein_rows returned a different protein set for another score"
        cols[f"auc_{k}"] = [r["auc"] for r in rows]
    tab = pd.DataFrame({"UniProt": base})
    extra = {"n_sites": [], "n_pos": [], "n_neg": [], "n_pairs": []}
    ties = {k: [] for k in scores}
    max_dev = 0.0
    for a in base:
        idx = np.flatnonzero(prot == a)
        yy = y[idx].astype(bool)
        extra["n_sites"].append(len(idx))
        extra["n_pos"].append(int(yy.sum()))
        extra["n_neg"].append(int((~yy).sum()))
        extra["n_pairs"].append(int(yy.sum() * (~yy).sum()))
        for k in scores:
            s = df[k].to_numpy(dtype=float)[idx]
            auc, nt, npair = mann_whitney_auc(s[yy], s[~yy])
            ties[k].append(nt)
            max_dev = max(max_dev, abs(auc - cols[f"auc_{k}"][len(extra["n_sites"]) - 1]))
    for c, v in extra.items():
        tab[c] = v
    for k in scores:
        tab[f"auc_{k}"] = cols[f"auc_{k}"]
        tab[f"tied_pairs_{k}"] = ties[k]
    assert max_dev < 1e-12, max_dev
    return tab, max_dev


def boot_primary(tab, codes_of_protein, n_clusters, seed, reps=REPS, scores=ALL_SCORES, stored=None):
    """Weighted mean of per-protein AUROCs over the both-label proteins of each replicate, on the
    protein-clustered resamples of all clusters (C.multiplicities, per_replicate mode)."""
    A = tab[[f"auc_{k}" for k in scores]].to_numpy(float)
    pairs = tab["n_pairs"].to_numpy(float)
    mean_b = np.empty((reps, len(scores)))
    pool_b = np.empty((reps, len(scores)))
    n_both = np.empty(reps)
    for b, mult in enumerate(C.multiplicities(n_clusters, reps, seed, "per_replicate")):
        m = mult[codes_of_protein].astype(float)
        n_both[b] = m.sum()
        mean_b[b] = (m @ A) / m.sum() if m.sum() > 0 else np.nan
        pool_b[b] = (m * pairs) @ A / (m * pairs).sum() if (m * pairs).sum() > 0 else np.nan
    return ({k: mean_b[:, i] for i, k in enumerate(scores)}, {k: pool_b[:, i] for i, k in enumerate(scores)}, n_both)


def boot_note12(tab, seed, reps=REPS, scores=ALL_SCORES):
    """Supplemental Note 12 scheme: resample only the proteins that carry both labels."""
    A = tab[[f"auc_{k}" for k in scores]].to_numpy(float)
    kb = len(A)
    rng = np.random.default_rng(seed)
    out = np.empty((reps, len(scores)))
    for b in range(reps):
        i = rng.integers(0, kb, kb)
        out[b] = A[i].mean(0)
    return {k: out[:, j] for j, k in enumerate(scores)}


def summarize(point, reps, estimand, scheme, n_prot, n_sites, n_pos, scores=PRIMARY):
    rows = []
    for k in scores:
        lo, hi = C.pctl(reps[k])
        rows.append({"estimand": estimand, "scheme": scheme, "quantity": "within-protein AUROC", "model": k,
                     "value": point[k], "ci_low": lo, "ci_high": hi, "share_replicates_le_0": None,
                     "n_invalid_replicates": int((~np.isfinite(reps[k])).sum()),
                     "n_proteins_both_labels": n_prot, "n_sites": n_sites, "n_positives": n_pos})
    for a, b in (("pLMSNOSite", "VIS10"), ("pLMSNOSite", "DIG25"), ("DIG25", "VIS10")):
        if a in scores and b in scores:
            d = reps[a] - reps[b]
            lo, hi = C.pctl(d)
            rows.append({"estimand": estimand, "scheme": scheme, "quantity": "within-protein AUROC difference",
                         "model": f"{a} - {b}", "value": point[a] - point[b], "ci_low": lo, "ci_high": hi,
                         "share_replicates_le_0": C.boot_p_le(d), "n_invalid_replicates": int((~np.isfinite(d)).sum()),
                         "n_proteins_both_labels": n_prot, "n_sites": n_sites, "n_positives": n_pos})
    for k in ("VIS10", "DIG25"):
        if k in scores:
            den = reps["pLMSNOSite"] - 0.5
            r = np.where(den > 0, (reps[k] - 0.5) / np.where(den > 0, den, 1.0), np.nan)
            lo, hi = C.pctl(r)
            rows.append({"estimand": estimand, "scheme": scheme,
                         "quantity": "recovered share of pLMSNOSite's above-chance within-protein AUROC",
                         "model": k, "value": (point[k] - 0.5) / (point["pLMSNOSite"] - 0.5), "ci_low": lo,
                         "ci_high": hi, "share_replicates_le_0": None,
                         "n_invalid_replicates": int((~np.isfinite(r)).sum()),
                         "n_proteins_both_labels": n_prot, "n_sites": n_sites, "n_positives": n_pos})
    return rows


def main():
    t0 = time.time()
    inputs = {}
    det = pd.read_csv(C.OUT / "detectability_test_scores.csv")
    plm = pd.read_csv(C.OUT / "plmsnosite_test_scores.csv")
    flg = pd.read_csv(C.OUT / "site_flags_and_combined_scores.csv")
    for f in ("detectability_test_scores.csv", "plmsnosite_test_scores.csv", "site_flags_and_combined_scores.csv",
              "robustness_other.csv", "bootstrap_replicates_seed20260930.csv.gz"):
        inputs[f"results/G_plmsnosite_rescore/{f}"] = C.sha256_file(C.OUT / f)
    inputs["repo/scripts/common.py (imported: protein_rows)"] = C.sha256_file(REPO_COMMON)
    for other in (plm, flg):
        assert (det.UniProt.values == other.UniProt.values).all()
        assert (det.Position.values == other.Position.values).all()
        assert (det.Target.values == other.Target.values).all()
    df = det[["row", "UniProt", "Position", "Target"]].copy()
    df["pLMSNOSite"] = plm["plmsnosite_prob"].to_numpy()
    df["VIS10"] = det["vis10_score"].to_numpy()
    df["DIG25"] = det["dig25_score"].to_numpy()
    df["ProtT5_arm"] = plm["prott5_base_prob"].to_numpy()
    df["Embedding_arm"] = plm["embedding_base_prob"].to_numpy()
    df["rankavg_pLMSNOSite_DIG25"] = flg["rankavg_pLMSNOSite_DIG25"].to_numpy()
    df["rankavg_pLMSNOSite_VIS10"] = flg["rankavg_pLMSNOSite_VIS10"].to_numpy()
    for c in ("seq_cluster", "flag_shared_id", "flag_seq_in_train", "flag_window_in_train", "flag_prott5_cropped",
              "flag_duplicate_site"):
        df[c] = flg[c].to_numpy()

    protein_rows = load_repo_protein_rows()
    checks = {}

    # ---------------------------------------------------------------- per-protein table (all sites)
    tab, max_dev = per_protein_table(df, ALL_SCORES, protein_rows)
    checks["protein_rows_vs_mann_whitney_max_abs_diff"] = max_dev
    codes, labels = C.protein_codes(df.UniProt.tolist())
    K = len(labels)
    pos_of = {a: i for i, a in enumerate(labels)}
    both_codes = np.array([pos_of[a] for a in tab.UniProt], dtype=np.int64)
    n_prot_both, n_sites_both, n_pos_both = len(tab), int(tab.n_sites.sum()), int(tab.n_pos.sum())
    point_mean = {k: float(tab[f"auc_{k}"].mean()) for k in ALL_SCORES}
    point_pool = {k: float((tab[f"auc_{k}"] * tab.n_pairs).sum() / tab.n_pairs.sum()) for k in ALL_SCORES}

    # the pooled estimand re-derived here equals step 30's values
    rob = pd.read_csv(C.OUT / "robustness_other.csv").set_index("quantity")
    dev = max(abs(point_pool[k] - rob.loc[f"pooled within-protein AUROC, {k}", "value"]) for k in PRIMARY)
    assert dev < 1e-9, dev
    checks["pooled_point_vs_step30_max_abs_diff"] = dev

    # ---------------------------------------------------------------- primary bootstrap (shared resamples)
    mean_b, pool_b, n_both_b = boot_primary(tab, both_codes, K, SEED)
    stored = pd.read_csv(C.OUT / "bootstrap_replicates_seed20260930.csv.gz")
    dev_pool = max(float(np.max(np.abs(pool_b[k] - stored[f"withinprot|{k}"].to_numpy()))) for k in PRIMARY)
    assert dev_pool < 1e-9, dev_pool
    ranked = C.RankedScore(df["pLMSNOSite"].to_numpy(float), df.Target.to_numpy())
    au = np.array([ranked.auroc(m[codes].astype(float)) for m in C.multiplicities(K, REPS, SEED, "per_replicate")])
    dev_au = float(np.max(np.abs(au - stored["auroc|pLMSNOSite"].to_numpy())))
    assert dev_au < 1e-9, dev_au
    checks["replicates_identical_to_step30"] = {
        "pooled_within_protein_max_abs_diff_vs_stored": dev_pool,
        "across_protein_auroc_pLMSNOSite_max_abs_diff_vs_stored": dev_au,
        "stored_file": "bootstrap_replicates_seed20260930.csv.gz (written by 30_compare.py)"}

    # one replicate expanded explicitly: a protein drawn twice enters as two distinct proteins
    mult1 = next(C.multiplicities(K, 1, 1))
    parts = []
    for code in np.flatnonzero(mult1):
        sub = df.loc[codes == code]
        for j in range(int(mult1[code])):
            s = sub.copy()
            s["UniProt"] = f"{labels[code]}#{j}"
            parts.append(s)
    exp = pd.concat(parts, ignore_index=True)
    y_e = exp.Target.to_numpy().astype(int)
    p_e = exp.UniProt.to_numpy().astype(str)
    dev_e = 0.0
    for k in PRIMARY:
        rows_e = protein_rows(y_e, exp[k].to_numpy(float), p_e, p_e)
        m_exp = float(np.mean([r["auc"] for r in rows_e]))
        m = mult1[both_codes].astype(float)
        m_w = float(m @ tab[f"auc_{k}"].to_numpy(float) / m.sum())
        dev_e = max(dev_e, abs(m_exp - m_w))
    assert dev_e < 1e-12, dev_e
    checks["weighted_mean_vs_explicitly_expanded_replicate_seed1_max_abs_diff"] = dev_e

    # ---------------------------------------------------------------- sensitivities
    n12_b = boot_note12(tab, SEED)
    alt_mean_b, _, _ = boot_primary(tab, both_codes, K, SEED_ALT)
    rows = []
    rows += summarize(point_mean, mean_b, EST_MEAN, SCH_MAIN, n_prot_both, n_sites_both, n_pos_both)
    rows += summarize(point_mean, n12_b, EST_MEAN, SCH_N12, n_prot_both, n_sites_both, n_pos_both)
    rows += summarize(point_mean, alt_mean_b, EST_MEAN, SCH_ALT, n_prot_both, n_sites_both, n_pos_both)
    rows += summarize(point_pool, pool_b, EST_POOL, SCH_MAIN, n_prot_both, n_sites_both, n_pos_both)
    # secondary scores (primary scheme only)
    for k in SECONDARY:
        lo, hi = C.pctl(mean_b[k])
        rows.append({"estimand": EST_MEAN, "scheme": SCH_MAIN, "quantity": "within-protein AUROC (secondary score)",
                     "model": k, "value": point_mean[k], "ci_low": lo, "ci_high": hi, "share_replicates_le_0": None,
                     "n_invalid_replicates": 0, "n_proteins_both_labels": n_prot_both, "n_sites": n_sites_both,
                     "n_positives": n_pos_both})
    for a, b in (("rankavg_pLMSNOSite_DIG25", "pLMSNOSite"), ("rankavg_pLMSNOSite_VIS10", "pLMSNOSite"),
                 ("Embedding_arm", "DIG25"), ("ProtT5_arm", "DIG25"), ("Embedding_arm", "pLMSNOSite")):
        d = mean_b[a] - mean_b[b]
        lo, hi = C.pctl(d)
        rows.append({"estimand": EST_MEAN, "scheme": SCH_MAIN, "quantity": "within-protein AUROC difference (secondary)",
                     "model": f"{a} - {b}", "value": point_mean[a] - point_mean[b], "ci_low": lo, "ci_high": hi,
                     "share_replicates_le_0": C.boot_p_le(d), "n_invalid_replicates": 0,
                     "n_proteins_both_labels": n_prot_both, "n_sites": n_sites_both, "n_positives": n_pos_both})
    summ = pd.DataFrame(rows)
    summ.to_csv(C.OUT / "withinprot_mean_summary.csv", index=False, float_format="%.10g")

    # ---------------------------------------------------------------- (d) subsets / clustering (primary scheme)
    subsets = {
        "all sites": np.ones(len(df), bool),
        "excl. 4 proteins shared by ID with training": ~df.flag_shared_id.to_numpy(bool),
        "excl. proteins whose sequence occurs in training": ~df.flag_seq_in_train.to_numpy(bool),
        "excl. sites whose 37-residue window occurs in training": ~df.flag_window_in_train.to_numpy(bool),
        "excl. sites with cropped ProtT5 features (New_Position != Position)": ~df.flag_prott5_cropped.to_numpy(bool),
        "excl. duplicated test entry": ~df.flag_duplicate_site.to_numpy(bool),
        "excl. all of the above": ~(df.flag_shared_id.to_numpy(bool) | df.flag_seq_in_train.to_numpy(bool)
                                    | df.flag_window_in_train.to_numpy(bool) | df.flag_prott5_cropped.to_numpy(bool)
                                    | df.flag_duplicate_site.to_numpy(bool)),
    }
    sub_rows = []
    for sname, mask in list(subsets.items()) + [("all sites, clusters = identical sequences", None)]:
        if mask is None:
            sub, ccol = df, "seq_cluster"
        else:
            sub, ccol = df.loc[mask].reset_index(drop=True), "UniProt"
        t_s, _ = per_protein_table(sub, PRIMARY, protein_rows)
        ccodes, clabels = C.protein_codes(sub[ccol].tolist())
        cluster_of = dict(zip(sub.UniProt.astype(str), sub[ccol].astype(str)))
        cpos = {c: i for i, c in enumerate(clabels)}
        cb = np.array([cpos[cluster_of[a]] for a in t_s.UniProt], dtype=np.int64)
        mb, _, _ = boot_primary(t_s, cb, len(clabels), SEED, scores=PRIMARY)
        pm = {k: float(t_s[f"auc_{k}"].mean()) for k in PRIMARY}
        row = {"subset": sname, "cluster_unit": ccol, "n_sites": len(sub), "n_pos": int(sub.Target.sum()),
               "n_clusters": len(clabels), "n_proteins_both_labels": len(t_s),
               "n_sites_in_them": int(t_s.n_sites.sum()), "n_pos_in_them": int(t_s.n_pos.sum())}
        for k in PRIMARY:
            lo, hi = C.pctl(mb[k])
            row.update({f"wp_mean_{k}": pm[k], f"wp_mean_{k}_ci_low": lo, f"wp_mean_{k}_ci_high": hi})
        for k in ("VIS10", "DIG25"):
            d = mb["pLMSNOSite"] - mb[k]
            lo, hi = C.pctl(d)
            row.update({f"diff_pLMSNOSite_minus_{k}": pm["pLMSNOSite"] - pm[k], f"diff_{k}_ci_low": lo,
                        f"diff_{k}_ci_high": hi})
            r = (mb[k] - 0.5) / (mb["pLMSNOSite"] - 0.5)
            lo, hi = C.pctl(r)
            row.update({f"recovery_{k}": (pm[k] - 0.5) / (pm["pLMSNOSite"] - 0.5), f"recovery_{k}_ci_low": lo,
                        f"recovery_{k}_ci_high": hi})
        sub_rows.append(row)
    subs = pd.DataFrame(sub_rows)
    assert abs(subs.iloc[0]["wp_mean_pLMSNOSite"] - point_mean["pLMSNOSite"]) < 1e-15
    subs.to_csv(C.OUT / "withinprot_mean_subsets.csv", index=False, float_format="%.10g")

    # ---------------------------------------------------------------- descriptive: sizes, ties, wins
    pairs = tab.n_pairs.to_numpy(float)
    order = np.sort(pairs)[::-1]
    q1, q2 = np.quantile(pairs, [1 / 3, 2 / 3])
    grp = np.where(pairs <= q1, 1, np.where(pairs <= q2, 2, 3))
    size_rows = []
    for g in (1, 2, 3):
        sel = grp == g
        r = {"pair_count_tertile": g, "pair_count_range": f"{int(pairs[sel].min())}-{int(pairs[sel].max())}",
             "n_proteins": int(sel.sum()), "n_sites": int(tab.n_sites[sel].sum()), "n_pos": int(tab.n_pos[sel].sum()),
             "share_of_all_within_protein_pairs": float(pairs[sel].sum() / pairs.sum())}
        for k in PRIMARY:
            r[f"mean_per_protein_auroc_{k}"] = float(tab.loc[sel, f"auc_{k}"].mean())
            r[f"pooled_auroc_{k}"] = float((tab.loc[sel, f"auc_{k}"] * pairs[sel]).sum() / pairs[sel].sum())
        size_rows.append(r)
    size_df = pd.DataFrame(size_rows)
    size_df.to_csv(C.OUT / "withinprot_mean_by_size.csv", index=False, float_format="%.10g")

    wins = {}
    for a, b in (("pLMSNOSite", "VIS10"), ("pLMSNOSite", "DIG25"), ("DIG25", "VIS10")):
        d = tab[f"auc_{a}"] - tab[f"auc_{b}"]
        wins[f"{a} vs {b}"] = {"higher": int((d > 1e-12).sum()), "equal": int((d.abs() <= 1e-12).sum()),
                               "lower": int((d < -1e-12).sum())}
    tie = {k: {"pooled_share_of_within_protein_pairs_tied": float(tab[f"tied_pairs_{k}"].sum() / pairs.sum()),
               "mean_per_protein_share_tied": float((tab[f"tied_pairs_{k}"] / pairs).mean()),
               "proteins_with_all_pairs_tied": int((tab[f"tied_pairs_{k}"] == pairs).sum())} for k in PRIMARY}
    excluded = df.loc[~df.UniProt.isin(set(tab.UniProt))]
    ex_by = excluded.groupby("UniProt").Target.agg(["size", "sum"])
    diag = {
        "n_test_proteins": K,
        "n_proteins_with_both_labels": n_prot_both, "n_sites_in_them": n_sites_both, "n_positives_in_them": n_pos_both,
        "n_proteins_excluded": int(len(ex_by)), "n_sites_excluded": int(len(excluded)),
        "n_positives_excluded": int(excluded.Target.sum()),
        "excluded_proteins_all_negative": int((ex_by["sum"] == 0).sum()),
        "excluded_proteins_all_positive": int((ex_by["sum"] == ex_by["size"]).sum()),
        "pairs_per_protein": {"min": float(pairs.min()), "q25": float(np.quantile(pairs, 0.25)),
                              "median": float(np.median(pairs)), "q75": float(np.quantile(pairs, 0.75)),
                              "max": float(pairs.max()), "total": float(pairs.sum())},
        "share_of_pairs_in_10_largest_proteins": float(order[:10].sum() / order.sum()),
        "proteins_with_one_positive": int((tab.n_pos == 1).sum()),
        "proteins_with_at_most_4_pairs": int((pairs <= 4).sum()),
        "pair_count_tertile_cut_points": [float(q1), float(q2)],
        "replicate_count_of_both_label_proteins": {"min": float(n_both_b.min()), "median": float(np.median(n_both_b)),
                                                    "max": float(n_both_b.max())},
        "within_protein_ties": tie,
        "per_protein_wins": wins,
    }
    tab_out = tab.copy()
    tab_out.insert(1, "seq_cluster", df.drop_duplicates("UniProt").set_index("UniProt").loc[tab.UniProt, "seq_cluster"].values)
    tab_out["pair_count_tertile"] = grp
    tab_out.to_csv(C.OUT / "withinprot_mean_per_protein.csv", index=False, float_format="%.17g")

    summary = {
        "label": ("POST HOC revision analysis 2026-09-30, revision round 2 (G_plmsnosite_rescore); not registered, "
                  "not pre-specified"),
        "estimand_primary": ("unweighted mean of per-protein AUROCs over the test proteins that carry both labels "
                             "(repo scripts/common.py protein_rows -> within_protein_auc_mean; the manuscript's "
                             "within-protein AUC)"),
        "estimand_sensitivity": ("pooled within-protein AUROC: concordant within-protein (positive, negative) pairs over "
                                 "all within-protein pairs = pair-weighted mean of the per-protein AUROCs (round 1)"),
        "bootstrap_primary": {"replicates": REPS, "seed": SEED, "unit": "protein (UniProt accession), sorted labels",
                              "draws": "numpy default_rng(seed); rng.integers(0, 278, 278) per replicate (all test proteins)",
                              "replicate_estimand": "multiplicity-weighted mean of per-protein AUROCs of resampled "
                                                    "proteins carrying both labels",
                              "shared_with": "30_compare.py (identity asserted) and 32_kr_stratified.py",
                              "interval": "95% percentile", "models_refit": False},
        "bootstrap_sensitivity_note12": {"replicates": REPS, "seed": SEED,
                                         "draws": "rng.integers(0, 220, 220) per replicate over the proteins with both labels"},
        "point_mean_of_per_protein_auroc": point_mean,
        "point_pooled_within_protein_auroc": point_pool,
        "diagnostics": diag,
        "checks": checks,
        "inputs_sha256": inputs,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                        "scipy": scipy.__version__, "scikit-learn": sklearn.__version__, "OMP_NUM_THREADS": "4"},
        "runtime_s": round(time.time() - t0, 1),
    }
    (C.OUT / "withinprot_mean.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")

    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    print(summ[["estimand", "scheme", "quantity", "model", "value", "ci_low", "ci_high", "share_replicates_le_0"]]
          .assign(estimand=lambda x: x.estimand.str[:22], scheme=lambda x: x.scheme.str[:28],
                  quantity=lambda x: x.quantity.str[:34]).round(4).to_string())
    print(subs.round(4).T.to_string())
    print(size_df.round(4).T.to_string())
    print(json.dumps(diag, indent=1))
    print(json.dumps(checks, indent=1))
    print("runtime", summary["runtime_s"], "s")


if __name__ == "__main__":
    main()
