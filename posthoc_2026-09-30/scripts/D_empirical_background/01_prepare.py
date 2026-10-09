"""Step 1 (POST HOC revision analysis, item D_empirical_background): build the cysteine-level table of the
PXD063463 trypsin acyl-biotin-exchange arm and verify every sequence against UniProt 2026_03.

Inputs (read-only)
  inputs/phase4_inputs/pxd063463_Trypsin_HydP.tsv   +hydroxylamine arm (label = carbamidomethyl site, loc >= 0.75)
  inputs/phase4_inputs/pxd063463_Trypsin_HydN.tsv   -hydroxylamine arm (specificity check)
  inputs/repo_results/phase3_candidate_peptides_2026-09-22.csv   Cys-containing Trypsin/P candidates (2025 FASTA
                                                    coordinates) with the global-arm identification flag
  external/UP000000589_10090.fasta.gz               UniProt 2026_03 mouse reference proteome

Outputs (results/D_empirical_background/)
  cys_table_trypsin_arm.csv      one row per cysteine of every identified protein that has a 2026_03 sequence
  cand_cys_map.csv               candidate-row -> contained cysteine (analysis proteins only, 2026_03 coordinates)
  prepare_checks.json            sequence verification and counts

Definitions
  positive (primary)      label = 1 in HydP (carbamidomethyl with localisation >= 0.75), exactly as in Supplemental Note 4
  positive (HA-specific)  label = 1 in HydP AND identified in HydN without carbamidomethyl (detected = 1, label = 0)
  proteome background     every other cysteine of the identified proteins (Cys-Audit --expand-background)
  observed background     identified in HydP, label 0
  kr_prox / kr_dist       trypsin-competent residue (K/R not before P) within 1-3 / 6-12 residues, either side
                          (Cys-Audit proteases.band_flag with rule 'trypsin', the code behind Fig. 2)
  de_prox / de_dist       D or E within 1-3 / 6-12 residues, either side (negative-control feature; rule 'gluc_de')
  theo_trypsin            in >= 1 fully tryptic peptide (K/R, not before P), <= 2 missed cleavages, 7-30 residues
                          (the manuscript's Detectability model; Cys-Audit proteases.digest)
  theo_trypsinP           the same with Trypsin/P (no proline rule; sensitivity)
  cand_n / cand_theo      number of containing phase-3 candidates / any containing candidate with length 7-30
  emp_obs                 in >= 1 phase-3 candidate identified in the global (unenriched) arm
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

from common import (IN_CAND, IN_FASTA_GZ, IN_HYDN, IN_HYDP, RESULTS, dump_json, read_fasta_gz)
from cys_audit import proteases  # noqa: E402  (path set in common)

WATER = 18.010565
MONO = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768,
        "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858,
        "K": 128.09496, "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111,
        "Y": 163.06333, "W": 186.07931}


def phase3_candidates(seq):
    """Re-implementation of the phase-3 candidate generator (Trypsin/P, 0-2 missed cleavages, 7-80 residues,
    <= 4600 Da, must contain C), used only to check that the stored candidate table matches the 2026_03 sequence."""
    bounds = sorted(set([0] + [i + 1 for i, ch in enumerate(seq) if ch in "KR"] + [len(seq)]))
    nb = len(bounds) - 1
    out = set()
    for i in range(nb):
        for mc in (0, 1, 2):
            j = i + mc
            if j >= nb:
                break
            s, e = bounds[i], bounds[j + 1]
            if e - s < 7 or e - s > 80:
                continue
            pep = seq[s:e]
            if "C" not in pep:
                continue
            if WATER + sum(MONO.get(r, 0.0) for r in pep) > 4600.0:
                continue
            out.add((s + 1, e, pep))
    return out


def covered(n, peps):
    """Boolean array (1-based index 1..n) marking residues covered by any (start1, end1) peptide."""
    cov = np.zeros(n + 2, dtype=np.int64)
    for s, e in peps:
        cov[s] += 1
        cov[e + 1] -= 1
    return np.cumsum(cov)[: n + 1] > 0


def main():
    seqs = read_fasta_gz(IN_FASTA_GZ)
    P = pd.read_csv(IN_HYDP, sep="\t", dtype={"protein": str})
    N = pd.read_csv(IN_HYDN, sep="\t", dtype={"protein": str})
    C = pd.read_csv(IN_CAND, dtype={"protein": str})
    C["cand_row"] = np.arange(len(C))
    checks = {"fasta_entries_2026_03": len(seqs), "hydp_rows": int(len(P)), "hydp_positives": int(P.label.sum()),
              "hydp_proteins": int(P.protein.nunique()), "hydn_rows": int(len(N)), "hydn_positives": int(N.label.sum())}

    # --- verify site tables against the 2026_03 sequences
    for name, T in (("hydp", P), ("hydn", N)):
        absent = ~T.protein.isin(seqs.keys())
        wrong = sum(1 for p, pos in zip(T.protein, T.position)
                    if p in seqs and (pos > len(seqs[p]) or seqs[p][pos - 1] != "C"))
        checks[f"{name}_rows_protein_absent_2026_03"] = int(absent.sum())
        checks[f"{name}_positives_protein_absent_2026_03"] = int(T.label[absent].sum())
        checks[f"{name}_proteins_absent_2026_03"] = sorted(T.protein[absent].unique().tolist())
        checks[f"{name}_rows_not_pointing_at_C"] = int(wrong)

    prots = sorted(p for p in P.protein.unique() if p in seqs)
    pset = set(prots)
    checks["analysis_proteins"] = len(prots)

    # --- verify the candidate table against the 2026_03 sequences (analysis proteins), remap changed proteins
    Ca = C[C.protein.isin(pset)].copy()
    ok = np.array([seqs[p][s - 1:e] == q for p, s, e, q in zip(Ca.protein, Ca.start, Ca.end, Ca.sequence)])
    Ca["seq_ok_2026"] = ok
    changed = sorted(Ca.protein[~ok].unique().tolist())
    new_start = Ca.start.to_numpy().copy()
    keep = ok.copy()
    remap_log = []
    for k in np.flatnonzero(~ok):
        p, q = Ca.protein.iat[k], Ca.sequence.iat[k]
        hits = [i for i in range(len(seqs[p]) - len(q) + 1) if seqs[p].startswith(q, i)]
        if len(hits) == 1:
            new_start[k] = hits[0] + 1
            keep[k] = True
        remap_log.append({"protein": p, "start_2025": int(Ca.start.iat[k]), "sequence": q,
                          "remapped_start_2026": (hits[0] + 1) if len(hits) == 1 else None})
    Ca["start26"] = new_start
    Ca["end26"] = Ca.start26 + Ca.sequence.str.len() - 1
    Ca = Ca[keep].copy()
    checks["candidates_analysis_proteins"] = int(ok.size)
    checks["candidates_mismatching_2026_03"] = int((~ok).sum())
    checks["candidates_remapped_uniquely"] = int(sum(1 for r in remap_log if r["remapped_start_2026"] is not None))
    checks["candidates_dropped_unmappable"] = int(sum(1 for r in remap_log if r["remapped_start_2026"] is None))
    checks["proteins_with_changed_sequence"] = changed

    # regenerate candidates for unchanged proteins: the stored table must equal the generator's output
    n_same = n_diff = 0
    for p in prots:
        if p in changed:
            continue
        want = set(zip(C.start[C.protein == p], C.end[C.protein == p], C.sequence[C.protein == p]))
        got = phase3_candidates(seqs[p])
        if want == got:
            n_same += 1
        else:
            n_diff += 1
    checks["regenerated_candidate_sets_identical"] = n_same
    checks["regenerated_candidate_sets_different"] = n_diff

    # --- candidate -> contained cysteine map (2026_03 coordinates)
    rows = []
    for r, p, s, q in zip(Ca.cand_row, Ca.protein, Ca.start26, Ca.sequence):
        for i, ch in enumerate(q):
            if ch == "C":
                rows.append((r, p, s + i))
    M = pd.DataFrame(rows, columns=["cand_row", "protein", "position"])
    M = M.merge(C[["cand_row", "empirical_detected", "theoretical_detectable", "length", "missed_cleavages"]],
                on="cand_row", how="left")
    M.to_csv(f"{RESULTS}/cand_cys_map.csv", index=False)

    agg = M.groupby(["protein", "position"]).agg(cand_n=("cand_row", "size"),
                                                  cand_theo=("theoretical_detectable", "max"),
                                                  emp_obs=("empirical_detected", "max"),
                                                  emp_n=("empirical_detected", "sum")).reset_index()

    # --- cysteine table
    Pk = P.set_index(["protein", "position"])
    Nk = N.set_index(["protein", "position"])
    recs = []
    for p in prots:
        s = seqs[p]
        n = len(s)
        m_kr = proteases.competent_mask(s, "trypsin")
        m_de = proteases.competent_mask(s, "gluc_de")
        theo = covered(n, proteases.digest(s, "trypsin", missed=2, min_len=7, max_len=30))
        theoP = covered(n, proteases.digest(s, "trypsin/p", missed=2, min_len=7, max_len=30))
        for i, ch in enumerate(s, 1):
            if ch != "C":
                continue
            key = (p, i)
            inP = key in Pk.index
            inN = key in Nk.index
            recs.append({
                "protein": p, "position": i, "protein_length": n,
                "label": int(Pk.at[key, "label"]) if inP else 0,
                "detected": int(inP),
                "n_cys_peptide": float(Pk.at[key, "n_cys_peptide"]) if inP else np.nan,
                "hydn_detected": int(inN),
                "hydn_label": int(Nk.at[key, "label"]) if inN else 0,
                "kr_prox": int(proteases.band_flag(s, i, (1, 3), "trypsin", m_kr)),
                "kr_dist": int(proteases.band_flag(s, i, (6, 12), "trypsin", m_kr)),
                "de_prox": int(proteases.band_flag(s, i, (1, 3), "gluc_de", m_de)),
                "de_dist": int(proteases.band_flag(s, i, (6, 12), "gluc_de", m_de)),
                "theo_trypsin": int(theo[i]), "theo_trypsinP": int(theoP[i]),
                "seq_changed_since_2025": int(p in changed),
            })
    T = pd.DataFrame(recs)
    T = T.merge(agg, on=["protein", "position"], how="left")
    for c in ("cand_n", "cand_theo", "emp_obs", "emp_n"):
        T[c] = T[c].fillna(0).astype(int)
    T["ha_specific"] = ((T.label == 1) & (T.hydn_detected == 1) & (T.hydn_label == 0)).astype(int)
    T.to_csv(f"{RESULTS}/cys_table_trypsin_arm.csv", index=False)

    # every HydP row whose protein has a 2026_03 sequence must land in the table as detected
    landed = set(zip(T.protein[T.detected == 1], T.position[T.detected == 1]))
    want = set((p, pos) for p, pos in zip(P.protein, P.position) if p in pset)
    checks["hydp_rows_with_sequence"] = len(want)
    checks["hydp_rows_landed"] = len(want & landed)
    pos_all = T.label == 1
    checks.update({
        "cys_rows": int(len(T)), "positives": int(pos_all.sum()), "identified": int(T.detected.sum()),
        "observed_background": int(((T.label == 0) & (T.detected == 1)).sum()),
        "proteome_background": int((T.label == 0).sum()),
        "ha_specific_positives": int(T.ha_specific.sum()),
        "hydp_positive_also_positive_in_hydn": int(((T.label == 1) & (T.hydn_label == 1)).sum()),
        "hydp_positive_not_identified_in_hydn": int(((T.label == 1) & (T.hydn_detected == 0)).sum()),
        "share_positive_theo_trypsin": float(T.theo_trypsin[pos_all].mean()),
        "share_positive_theo_trypsinP": float(T.theo_trypsinP[pos_all].mean()),
        "share_positive_in_any_candidate": float((T.cand_n[pos_all] > 0).mean()),
        "share_positive_emp_obs": float(T.emp_obs[pos_all].mean()),
        "share_proteome_bg_theo_trypsin": float(T.theo_trypsin[T.label == 0].mean()),
        "share_proteome_bg_emp_obs": float(T.emp_obs[T.label == 0].mean()),
        "share_observed_bg_emp_obs": float(T.emp_obs[(T.label == 0) & (T.detected == 1)].mean()),
        "proteins_with_any_emp_obs_cys": int(T.groupby("protein").emp_obs.max().sum()),
        "remap_log": remap_log,
    })
    dump_json(checks, f"{RESULTS}/prepare_checks.json")
    for k, v in checks.items():
        if k != "remap_log":
            print(k, v)


if __name__ == "__main__":
    main()
