"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific, task 2.

Hydroxylamine-specific site definitions for the four-protease ABE deposit PXD063463, from the arms processed without
hydroxylamine (HydN). A HydP 'site' is a cysteine carrying carbamidomethyl (CAM, localisation >= 0.75) in the HydP arm
(the stored definition). New definitions, written after the stored results had been seen:

  strict   CAM in HydP AND the same cysteine identified in the HydN arm WITHOUT CAM
  lenient  CAM in HydP AND NOT CAM in the HydN arm (not identified in HydN, or identified without CAM)
  nonspecific (control) CAM in HydP AND CAM in HydN
  untested  CAM in HydP and not identified in HydN (lenient minus strict)
  cam_protein_absent (descriptive) CAM in HydP on a protein never identified in the HydN arm(s): a presence/absence
           stand-in for the deposit's own protein-level enrichment contrast, which needs intensities these tables lack

HydN source
  matched  (primary, as specified) the HydN arm digested with the same protease
  pooled   (sensitivity) the union of the four HydN arms. Rationale, ASSUMED rather than verified: Supplemental
           Note 4 describes one enrichment split across four proteases; if each enrichment (with or without
           hydroxylamine) was labelled (NEM block, +/- hydroxylamine, capture, iodoacetamide) before the split, whether
           a cysteine carries CAM without hydroxylamine should not depend on the protease that later digested it.
           The deposit's protocol could not be checked offline and the tables do not show it.
Matching key
  position (primary) (protein accession, residue number)
  window   (sensitivity) the 21-residue sequence window centred on the cysteine in the 2026_03 FASTA, which also
           matches a cysteine reported under a different razor protein with identical local sequence

Revision after verification (round 1): the site tables also carry in_hydn_trypsin_arm / cam_hydn_trypsin_arm
(identification / CAM in the TRYPTIC HydN arm, position key) for every arm, because the pooled HydN identifications are
95.5% tryptic and selection through that arm has to be examined separately (t2_selection_summary.csv; c03).

Naming note: 'strict' and 'lenient' are presence/absence definitions of CANDIDATE hydroxylamine-dependent sites. They
do not measure hydroxylamine dependence, which needs intensities the converted tables do not carry.

Outputs: t2_site_table_<arm>.tsv (one row per HydP-identified cysteine), t2_specific_counts.csv, t2_hydn_summary.csv,
t2_selection_summary.csv, in results/<ID>/ for the FASTA-filtered tables (default) or results/<ID>/unfiltered_base/ with
C_BASE=unfiltered (the deposit-converted tables as stored, the base of the manuscript's current counts).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

from cys_audit import stats  # noqa: E402
from cys_audit.io import read_fasta  # noqa: E402


def keyset(df, key):
    return set(df[key])


def clustered_share(proteins, success, seed):
    """Protein-clustered percentile bootstrap 95% interval of a share (Cys-Audit stats.boot_proportion)."""
    success = np.asarray(success, dtype=bool)
    if len(success) == 0:
        return np.nan, np.nan, np.nan
    codes, labels = stats.cluster_index(list(proteins))
    pt, (lo, hi), _ = stats.boot_proportion(codes, len(labels), success, C.REPS, seed, 0.95)
    return pt, lo, hi


def main():
    seqs = read_fasta(C.ensure_fasta())
    os.makedirs(C.OUT, exist_ok=True)
    C.record_inputs(f"c02_specific_sites[{C.BASE}]", [C.arm_path(a, h) for a in C.ARMS for h in ("HydP", "HydN")]
                    + [C.FASTA])
    tabs = {}
    for a in C.ARMS:
        for h in ("HydP", "HydN"):
            d = C.read_arm(a, h)
            d["window"] = [C.window(seqs.get(p), int(pos)) for p, pos in zip(d.protein, d.position)]
            d["key"] = list(zip(d.protein, d.position))
            # a row whose protein is absent from the FASTA has no window: fall back to its position key, so that
            # absent-window rows never match each other
            d["wkey"] = [w if w is not None else ("pos", p, int(q)) for w, p, q in zip(d.window, d.protein, d.position)]
            tabs[(a, h)] = d
    # pooled HydN
    pooled = pd.concat([tabs[(a, "HydN")] for a in C.ARMS], ignore_index=True)
    sources = {}
    for a in C.ARMS:
        N = tabs[(a, "HydN")]
        sources[(a, "matched")] = N
        sources[(a, "pooled")] = pooled
    count_rows, hyd_rows, sel_rows = [], [], []
    seed = C.NEW_SEED
    for a in C.ARMS:
        P = tabs[(a, "HydP")].copy()
        out = P[["protein", "position", "label", "detected", "n_cys_peptide", "abundance", "peptide_mass", "window"]].copy()
        out = out.rename(columns={"label": "cam_hydp"})
        for src in ("matched", "pooled"):
            N = sources[(a, src)]
            for key in ("key", "wkey"):
                idN = keyset(N, key)
                camN = keyset(N[N.label == 1], key)
                tag = f"{src}" if key == "key" else f"{src}_window"
                in_n = np.array([k in idN for k in P[key]])
                cam_n = np.array([k in camN for k in P[key]])
                cam_p = P.label.values == 1
                out[f"in_hydn_{tag}"] = in_n.astype(int)
                out[f"cam_hydn_{tag}"] = cam_n.astype(int)
                out[f"strict_{tag}"] = (cam_p & in_n & ~cam_n).astype(int)
                out[f"lenient_{tag}"] = (cam_p & ~cam_n).astype(int)
                out[f"nonspecific_{tag}"] = (cam_p & cam_n).astype(int)
                n_cam = int(cam_p.sum())
                n_testable = int((cam_p & in_n).sum())
                n_nonspec = int((cam_p & cam_n).sum())
                # share of HydP sites also CAM in HydN (denominator: all HydP sites)
                sub = P.loc[cam_p]
                f_all = clustered_share(sub.protein, cam_n[cam_p], seed); seed += 1
                # share among HydP sites whose cysteine was identified in HydN at all (testable)
                sel = cam_p & in_n
                f_test = clustered_share(P.protein[sel], cam_n[sel], seed) if sel.any() else (np.nan,) * 3; seed += 1
                wl_all = stats.wilson(n_nonspec, n_cam)
                wl_test = stats.wilson(n_nonspec, n_testable) if n_testable else (np.nan, np.nan)
                count_rows.append({
                    "arm": a, "hydn_source": src, "match_key": "position" if key == "key" else "window21",
                    "hydp_identified": len(P), "hydp_cam_sites": n_cam,
                    "hydn_identified": len(set(N[key])), "hydn_cam": len(camN),
                    "hydp_identified_also_identified_in_hydn": int(in_n.sum()),
                    "hydp_sites_identified_in_hydn": n_testable,
                    "nonspecific_cam_in_both": n_nonspec,
                    "strict": int((cam_p & in_n & ~cam_n).sum()),
                    "lenient": int((cam_p & ~cam_n).sum()),
                    "untested_not_identified_in_hydn": int((cam_p & ~in_n).sum()),
                    "hydn_cam_not_cam_in_hydp": len(camN - set(P.loc[cam_p, key])),
                    "frac_sites_also_cam_in_hydn": n_nonspec / n_cam,
                    "frac_sites_also_cam_in_hydn_ci_low": f_all[1], "frac_sites_also_cam_in_hydn_ci_high": f_all[2],
                    "frac_sites_also_cam_in_hydn_wilson_low": wl_all[0], "frac_sites_also_cam_in_hydn_wilson_high": wl_all[1],
                    "frac_testable_sites_cam_in_hydn": n_nonspec / n_testable if n_testable else np.nan,
                    "frac_testable_ci_low": f_test[1], "frac_testable_ci_high": f_test[2],
                    "frac_testable_wilson_low": wl_test[0], "frac_testable_wilson_high": wl_test[1],
                    "frac_hydp_identified_testable": in_n.mean(),
                })
        # protein-level stand-in for the deposit's own enrichment contrast: the protein was never identified in the
        # HydN arm(s) at all (presence/absence only; the tables carry no HydP/HydN intensities)
        for src in ("matched", "pooled"):
            protN = set(sources[(a, src)].protein)
            out[f"protein_in_hydn_{src}"] = P.protein.isin(protN).astype(int).values
            out[f"cam_protein_absent_{src}"] = ((P.label.values == 1) & ~P.protein.isin(protN).values).astype(int)
        # Revision after verification (round 1): identification in the TRYPTIC HydN arm specifically. The pooled HydN
        # identifications are 95.5% tryptic, so selection through that arm acts on every pooled definition and on the
        # strict trypsin set; c03 splits each arm's sites by this column (position key).
        NT = tabs[("Trypsin", "HydN")]
        idT, camT = set(NT.key), set(NT[NT.label == 1].key)
        in_t = np.array([k in idT for k in P.key])
        cam_t = np.array([k in camT for k in P.key])
        out["in_hydn_trypsin_arm"] = in_t.astype(int)
        out["cam_hydn_trypsin_arm"] = cam_t.astype(int)
        cam_p = P.label.values == 1
        in_m = out["in_hydn_matched"].values == 1
        in_pool = out["in_hydn_pooled"].values == 1
        strict_pool = out["strict_pooled"].values == 1
        sel_rows.append({"arm": a, "hydp_sites": int(cam_p.sum()),
                         "sites_identified_in_own_hydn_arm": int((cam_p & in_m).sum()),
                         "sites_identified_in_tryptic_hydn_arm": int((cam_p & in_t).sum()),
                         "sites_identified_in_pooled_hydn": int((cam_p & in_pool).sum()),
                         "sites_identified_in_pooled_hydn_only_via_tryptic_arm":
                             int((cam_p & in_pool & in_t & ~in_m).sum()) if a != "Trypsin" else 0,
                         "sites_cam_in_tryptic_hydn_arm": int((cam_p & cam_t).sum()),
                         "strict_pooled": int(strict_pool.sum()),
                         "strict_pooled_identified_in_own_hydn_arm": int((strict_pool & in_m).sum()),
                         "strict_pooled_identified_in_tryptic_hydn_arm": int((strict_pool & in_t).sum())})
        out.to_csv(f"{C.OUT}/t2_site_table_{a}.tsv", sep="\t", index=False, lineterminator="\n")
        # HydN arm on its own: CAM share of identified (the specificity-control readout)
        N = tabs[(a, "HydN")]
        ci_n = clustered_share(N.protein, N.label.values == 1, seed); seed += 1
        ci_p = clustered_share(P.protein, P.label.values == 1, seed); seed += 1
        prot_abs = ~P.protein.isin(set(pooled.protein))
        hyd_rows.append({"arm": a, "hydp_identified": len(P), "hydp_cam": int((P.label == 1).sum()),
                         "hydp_proteins_absent_from_all_hydn": P.protein[prot_abs].nunique(),
                         "hydp_cys_on_proteins_absent_from_all_hydn": int(prot_abs.sum()),
                         "hydp_cam_on_proteins_absent_from_all_hydn": int((P.label[prot_abs] == 1).sum()),
                         "hydp_share": ci_p[0], "hydp_share_ci_low": ci_p[1], "hydp_share_ci_high": ci_p[2],
                         "hydn_identified": len(N), "hydn_cam": int((N.label == 1).sum()),
                         "hydn_share": ci_n[0], "hydn_share_ci_low": ci_n[1], "hydn_share_ci_high": ci_n[2],
                         "hydn_proteins": N.protein.nunique(), "hydp_proteins": P.protein.nunique(),
                         "ratio_identified_hydp_over_hydn": len(P) / len(N)})
    cnt = pd.DataFrame(count_rows)
    cnt.to_csv(f"{C.OUT}/t2_specific_counts.csv", index=False)
    hyd = pd.DataFrame(hyd_rows)
    hyd.to_csv(f"{C.OUT}/t2_hydn_summary.csv", index=False)
    sel = pd.DataFrame(sel_rows)
    sel.to_csv(f"{C.OUT}/t2_selection_summary.csv", index=False)
    print(sel.to_string(index=False))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(hyd.to_string(index=False))
    print(cnt[["arm", "hydn_source", "match_key", "hydp_cam_sites", "hydp_sites_identified_in_hydn",
               "nonspecific_cam_in_both", "strict", "lenient", "untested_not_identified_in_hydn",
               "frac_sites_also_cam_in_hydn", "frac_sites_also_cam_in_hydn_ci_low", "frac_sites_also_cam_in_hydn_ci_high",
               "frac_testable_sites_cam_in_hydn", "frac_testable_ci_low", "frac_testable_ci_high",
               "frac_hydp_identified_testable"]].to_string(index=False))


if __name__ == "__main__":
    main()
