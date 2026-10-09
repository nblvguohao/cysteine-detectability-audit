"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific: revision after verification (round 2).

Problem raised (major): c03 split each arm's CAM sites by whether their cysteine was identified in the TRYPTIC arm
without hydroxylamine (HydN) but kept the UNSPLIT background. Identification in a tryptic digest raises the share of
cysteines with a K/R within 6-12 residues in every class of cysteine (sites, observed non-sites, unobserved cysteines),
so a subset selected on it looks 'tryptic' against an unselected background whatever its site status. Here the
background is restricted to the same identification stratum as the positives.

(a) Split sets, stratum-matched background (primary, as specified by the verifier). For each arm (HydP, 2026_03-filtered
    tables, --expand-background), rows are split by identification in the tryptic HydN arm (position key; for the
    expanded, unlisted cysteines the key is looked up in the HydN table directly). Within each stratum the positives are
    the CAM sites of that stratum and the background is every label-0 row of the SAME stratum (proteome: every other
    cysteine of the identified proteins; observed: the identified cysteines without CAM). Own rule and, for the
    non-tryptic arms, the trypsin rule; both bands; both backgrounds. Sensitivity: the same split by identification in
    ANY HydN arm (pooled).
(b) Candidate hydroxylamine-dependent sets within their own identification stratum: strict_<src> (identified without CAM
    in the HydN source) and nonspecific_<src> (CAM in both) against the label-0 rows identified in that HydN source
    (matched = same protease; pooled = any arm); CAM sites of the stratum outside the set are excluded (c03's 'exclude').
    Lenient sets span both strata (strict plus not identified) and are not re-audited here.
(c) Descriptive: share of rows flagged by each rule and band, by class (site / observed non-site / unobserved cysteine)
    and tryptic-HydN stratum (the mechanism behind the problem).
(d) Descriptive reconciliation for the stored CAM sets: Mantel-Haenszel log2 odds ratio (Cys-Audit stats.mh_log2_or,
    as used by the tool's claim re-test) stratified by tryptic-HydN identification, by abundance decile (c05 edges,
    Global intensity, + missing class), and by both; 'none' is the crude (uncorrected) odds ratio from the same counts.
    Protein-clustered bootstrap (Cys-Audit stats.multiplicities), 5,000 replicates, 97.5% intervals (as the cleavage
    bands), seed 20260930 + 4000 + running offset (listed per row). Stratifying on identification conditions on a
    variable that depends on cleavage geometry and on abundance; the abundance strata are the check on that.

Every cleavage statistic in (a) and (b) is Cys-Audit v0.2.2 checks.cleavage.run, unchanged, with the seed of every other
cleavage statistic of this item (20260922; bands seed + 101 + j), 5,000 replicates, 97.5% intervals; fewer than 20
positives or 20 background rows is UNDECIDABLE by the tool's fixed minimum.
Outputs: t11_stratum_matched_cleavage.csv, t11_flag_shares_by_stratum.csv, t11_mh_stratified.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402
import c05_depth_confound as D5  # noqa: E402  (decile edges, identical to c05)

from cys_audit import stats  # noqa: E402
from cys_audit.checks import cleavage  # noqa: E402
from cys_audit.constants import CLEAVAGE_BANDS, CLEAVAGE_CI_LEVEL  # noqa: E402
from cys_audit.io import background_mask, read_sites  # noqa: E402

MH_SEED0 = C.NEW_SEED + 4000


def run_set(ds, rule, bg):
    """Cys-Audit cleavage.run, unchanged (same call as c03.run_set), one row per band."""
    res = cleavage.run(ds, {"protease": rule, "background": bg, "seed": C.TOOL_SEED, "reps": C.REPS})
    flags, has_seq = cleavage.band_flags(ds, rule)
    bgm = background_mask(ds, bg) & has_seq
    pos = (ds.label == 1) & has_seq
    rows = []
    bands = (res.get("details") or {}).get("bands") or {}
    for b in CLEAVAGE_BANDS:
        t = bands.get(b)
        rows.append({"band": b, "n_positive": int(pos.sum()), "n_background": int(bgm.sum()),
                     "n_clusters": res.get("n_clusters"),
                     "estimate": None if t is None else t["estimate"],
                     "ci_low": None if t is None else t["ci"][0], "ci_high": None if t is None else t["ci"][1],
                     "status": res["status"] if t is None else t["status"],
                     "share_flagged_positive": float(flags[b][pos].mean()) if pos.any() else np.nan,
                     "share_flagged_background": float(flags[b][bgm].mean()) if bgm.any() else np.nan,
                     "below_minimum": t is None,
                     "reason": res["reason"] if t is None else t["reason"]})
    return rows


def keyset(path):
    d = pd.read_csv(path, sep="\t", dtype={"protein": str})
    return set(zip(d.protein, d.position.astype(int)))


def mh_boot(codes, K, s_index, S, ispos, flag, reps, seed, level):
    """Mantel-Haenszel log2 OR over strata with a protein-clustered bootstrap, following Cys-Audit verdict.py:
    strata are informative when they hold both classes in the full sample; each replicate uses the informative strata
    that still hold both classes."""
    npos_s = np.bincount(s_index, weights=ispos.astype(float), minlength=S)
    nbg_s = np.bincount(s_index, weights=(~ispos).astype(float), minlength=S)
    informative = (npos_s > 0) & (nbg_s > 0)
    orphaned_pos = float(npos_s[~informative].sum() / npos_s.sum()) if npos_s.sum() else np.nan
    orphaned_bg = float(nbg_s[~informative].sum() / nbg_s.sum()) if nbg_s.sum() else np.nan
    inf_rows = informative[s_index]
    cell = np.zeros((K, S * 4), dtype=np.int64)
    for j, sel in enumerate((ispos & flag, ispos & ~flag, ~ispos & flag, ~ispos & ~flag)):
        m = sel & inf_rows
        np.add.at(cell, (codes[m], s_index[m] * 4 + j), 1)
    tot = cell.sum(axis=0).reshape(S, 4)
    point = stats.mh_log2_or(tot[informative])
    vals = []
    for w in stats.multiplicities(K, reps, seed):
        tab = (w @ cell).reshape(w.shape[0], S, 4)
        for r in range(w.shape[0]):
            t = tab[r][informative]
            t = t[t.sum(axis=1) > 0]
            has_both = (t[:, 0] + t[:, 1] > 0) & (t[:, 2] + t[:, 3] > 0)
            vals.append(stats.mh_log2_or(t[has_both]) if has_both.any() else np.nan)
    vals = np.asarray(vals, dtype=float)
    lo, hi = stats.percentile_interval(vals, level)
    return {"estimate": point, "ci_low": lo, "ci_high": hi, "n_strata": S, "n_informative_strata": int(informative.sum()),
            "orphaned_positive_share": orphaned_pos, "orphaned_background_share": orphaned_bg,
            "n_undefined_reps": int((~np.isfinite(vals)).sum())}


def main():
    C.ensure_fasta()
    hydn_paths = {a: C.filtered_path(a, "HydN") for a in C.ARMS}
    C.record_inputs("c11_stratum_matched", [C.filtered_path(a, "HydP") for a in C.ARMS] + list(hydn_paths.values())
                    + [f"{C.RES}/t2_site_table_{a}.tsv" for a in C.ARMS] + [C.FASTA])
    keysN = {a: keyset(hydn_paths[a]) for a in C.ARMS}
    keys_pooled = set().union(*keysN.values())
    # c05 decile edges (log10 Global intensity over the union of HydP proteins of the filtered base)
    st5 = D5.load()
    e10, _ = D5.edges_from(st5, 10)

    out, fl, mh = [], [], []
    mh_seed = MH_SEED0
    for arm in C.ARMS:
        ds = read_sites(C.filtered_path(arm, "HydP"), fasta=C.FASTA, expand_background=True)
        k = list(zip(ds.protein, ds.position.tolist()))
        inT = np.array([x in keysN["Trypsin"] for x in k])
        inP = np.array([x in keys_pooled for x in k])
        inM = np.array([x in keysN[arm] for x in k])
        cam = ds.label == 1
        obs = ds.detected == 1
        # consistency with the c02 site table (listed rows) and with c03's candidate sets
        st = pd.read_csv(f"{C.RES}/t2_site_table_{arm}.tsv", sep="\t", dtype={"protein": str})
        idx = {(p, int(q)): i for i, (p, q) in enumerate(zip(st.protein, st.position))}
        rowmap = np.array([idx.get(x, -1) for x in k])
        listed = rowmap >= 0
        assert listed.sum() == len(st) and np.array_equal(listed, obs), "listed rows must be the identified rows"
        for col, m in (("in_hydn_trypsin_arm", inT), ("in_hydn_pooled", inP), ("in_hydn_matched", inM)):
            assert np.array_equal(m[listed], st[col].values[rowmap[listed]] == 1), col
        spec = {}
        for s in ("strict_matched", "nonspecific_matched", "strict_pooled", "nonspecific_pooled"):
            v = np.zeros(len(ds), dtype=bool)
            v[listed] = st[s].values[rowmap[listed]] == 1
            spec[s] = v
        assert not (spec["strict_matched"] & ~inM).any() and not (spec["strict_pooled"] & ~inP).any()
        rules = [C.OWN_RULE[arm]] + ([] if arm == "Trypsin" else ["trypsin"])

        # (a) split sets with stratum-matched backgrounds
        jobs = []
        for strat, S in (("trypsin_hydn", inT), ("pooled_hydn", inP)):
            for side, m in (("in", S), ("not_in", ~S)):
                jobs.append((f"{side}_{strat}", "split_stratum_matched", strat, m, None))
        # (b) candidate sets within their identification stratum (exclude handling inside the stratum)
        for s, S, strat in (("strict_matched", inM, "matched_hydn"), ("nonspecific_matched", inM, "matched_hydn"),
                            ("strict_pooled", inP, "pooled_hydn"), ("nonspecific_pooled", inP, "pooled_hydn")):
            m = S & ~(cam & ~spec[s])
            jobs.append((s, "candidate_in_stratum", strat, m, spec[s]))
        for name, design, strat, m, lab in jobs:
            d2 = ds.subset(m)
            if lab is not None:
                d2.label = lab[m].astype(ds.label.dtype)
            for rule in rules:
                for bg in ("proteome", "observed"):
                    for r in run_set(d2, rule, bg):
                        r.update({"arm": arm, "positive_set": name, "design": design, "stratum_variable": strat,
                                  "rule": rule, "rule_type": "own" if rule == C.OWN_RULE[arm] else "trypsin",
                                  "background": bg, "n_rows_in_subset": int(m.sum()),
                                  "n_sites_in_set": int((d2.label == 1).sum())})
                        out.append(r)

        # (c) flag shares by class and tryptic-HydN stratum
        cls = {"site": cam, "observed_nonsite": (~cam) & obs, "unobserved": (~cam) & (~obs)}
        for rule in rules:
            flags, has_seq = cleavage.band_flags(ds, rule)
            for b in CLEAVAGE_BANDS:
                for cn, cm in cls.items():
                    for sn, sm in (("in_trypsin_hydn", inT), ("not_in_trypsin_hydn", ~inT), ("all", np.ones(len(ds), bool))):
                        mm = cm & sm & has_seq
                        fl.append({"arm": arm, "rule": rule, "rule_type": "own" if rule == C.OWN_RULE[arm] else "trypsin",
                                   "band": b, "class": cn, "stratum": sn, "n": int(mm.sum()),
                                   "n_flagged": int(flags[b][mm].sum()),
                                   "share_flagged": float(flags[b][mm].mean()) if mm.any() else np.nan})

        # (d) Mantel-Haenszel reconciliation for the stored CAM set
        la = np.log10(ds.abundance)
        dec = D5.assign_bins(pd.DataFrame({"log10_abundance": la}), e10)          # 0..9, 10 = missing
        strata_defs = {"none": np.zeros(len(ds), dtype=np.int64),
                       "trypsin_hydn": inT.astype(np.int64),
                       "abundance_decile": dec.astype(np.int64),
                       "trypsin_hydn_x_abundance_decile": (inT.astype(np.int64) * 11 + dec).astype(np.int64)}
        for rule in rules:
            flags, has_seq = cleavage.band_flags(ds, rule)
            for bg in ("proteome", "observed"):
                bgm = background_mask(ds, bg) & has_seq
                pos = cam & has_seq
                keep = pos | bgm
                codes, labels = stats.cluster_index([c for c, kk in zip(ds.cluster, keep) if kk])
                kidx = np.flatnonzero(keep)
                for sname, sv in strata_defs.items():
                    lv = np.unique(sv[kidx])
                    s_index = np.searchsorted(lv, sv[kidx])
                    for j, b in enumerate(CLEAVAGE_BANDS):
                        res = mh_boot(codes, len(labels), s_index, len(lv), pos[kidx], flags[b][kidx], C.REPS, mh_seed,
                                      CLEAVAGE_CI_LEVEL)
                        res.update({"arm": arm, "set": "cam", "rule": rule,
                                    "rule_type": "own" if rule == C.OWN_RULE[arm] else "trypsin", "background": bg,
                                    "band": b, "strata": sname, "n_positive": int(pos.sum()), "n_background": int(bgm.sum()),
                                    "n_clusters": len(labels), "seed": mh_seed,
                                    "positives_in_trypsin_hydn_share": float(inT[pos].mean()),
                                    "background_in_trypsin_hydn_share": float(inT[bgm].mean())})
                        mh.append(res)
                        mh_seed += 1
        print("done", arm, flush=True)

    cols = ["arm", "positive_set", "design", "stratum_variable", "rule", "rule_type", "background", "band",
            "n_rows_in_subset", "n_sites_in_set", "n_positive", "n_background", "n_clusters", "estimate", "ci_low",
            "ci_high", "status", "below_minimum", "share_flagged_positive", "share_flagged_background", "reason"]
    df = pd.DataFrame(out)[cols]
    df.to_csv(f"{C.RES}/t11_stratum_matched_cleavage.csv", index=False)
    fdf = pd.DataFrame(fl)
    fdf.to_csv(f"{C.RES}/t11_flag_shares_by_stratum.csv", index=False)
    mcols = ["arm", "set", "rule", "rule_type", "background", "band", "strata", "n_positive", "n_background", "n_clusters",
             "positives_in_trypsin_hydn_share", "background_in_trypsin_hydn_share", "n_strata", "n_informative_strata",
             "orphaned_positive_share", "orphaned_background_share", "estimate", "ci_low", "ci_high", "n_undefined_reps",
             "seed"]
    mdf = pd.DataFrame(mh)[mcols]
    mdf.to_csv(f"{C.RES}/t11_mh_stratified.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.max_rows", 1000)
    print(df[["arm", "positive_set", "rule_type", "background", "band", "n_positive", "n_background", "estimate",
              "ci_low", "ci_high", "status"]].to_string(index=False))
    print(fdf[fdf.band == "distal_6_12"].to_string(index=False))
    print(mdf[["arm", "rule_type", "background", "band", "strata", "estimate", "ci_low", "ci_high",
               "orphaned_positive_share"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
