"""Phase 3: theoretical vs empirical Cys detectability on Dataset B's Global (unenriched) arm.

Every parameter is read from the pre-registration; none is redefined here. See
scripts/preregister_phase3_empirical_detectability_2026-09-22.py for the full specification: digestion
rule, candidate filters, protein universe, target, features, model, comparisons and pass rule.

Gate A0 refuses to run if the protocol JSON's sha256 differs from what the pre-registration's own audit
recorded, or if any declared input file's sha256 has changed since registration.

OUTPUTS (new names)
  results/phase3_candidate_peptides_2026-09-22.csv          one row per candidate Cys-containing peptide
  results/phase3_empirical_detectability_summary_2026-09-22.json
  results/phase3_empirical_detectability_2026-09-22_audit.json
Interpreter: project venv (numpy, scikit-learn).
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression

csv.field_size_limit(2 ** 28)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTOCOL = os.path.join(ROOT, "protocols", "phase3_empirical_detectability_preregistration_2026-09-22.json")
PREREG_AUDIT = os.path.join(ROOT, "results", "phase3_empirical_detectability_prereg_2026-09-22_audit.json")
DATA_DIR = os.path.join(ROOT, "external", "phase1_public_20260922", "PXD063463")
OUT_CAND = os.path.join(ROOT, "results", "phase3_candidate_peptides_2026-09-22.csv")
OUT_SUM = os.path.join(ROOT, "results", "phase3_empirical_detectability_summary_2026-09-22.json")
OUT_AUD = os.path.join(ROOT, "results", "phase3_empirical_detectability_2026-09-22_audit.json")

# standard monoisotopic residue masses (Da) and Kyte-Doolittle hydrophobicity -- textbook constants
MONO = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768,
        "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858,
        "K": 128.09496, "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111,
        "Y": 163.06333, "W": 186.07931}
WATER = 18.010565
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2,
      "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def read_fasta(path):
    seqs, acc, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                m = re.match(r"^>(sp|tr)\|([^|]+)\|", line)
                acc = m.group(2) if m else None
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def dread(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def acc_set(rows, col):
    s = set()
    for r in rows:
        if r.get("Reverse") == "+" or r.get("Potential contaminant") == "+":
            continue
        v = r[col]
        for tok in v.split(";"):
            m = re.match(r"^(sp|tr)\|([^|]+)\|", tok.strip())
            if m:
                s.add(m.group(2))
    return s


def cleave_sites(seq):
    """Trypsin/P: cut C-terminal to every K or R (site index = residue AFTER the cut, 0-based)."""
    return [i + 1 for i, ch in enumerate(seq) if ch in "KR"]


def candidates_for_protein(acc, seq, mc_set, min_len, max_len_cap, max_mass, must_residue):
    bounds = [0] + cleave_sites(seq) + [len(seq)]
    bounds = sorted(set(bounds))
    nb = len(bounds) - 1
    out = []
    for i in range(nb):
        for mc in mc_set:
            j = i + mc
            if j >= nb:
                break
            start, end = bounds[i], bounds[j + 1]
            L = end - start
            if L < min_len or L > max_len_cap:
                continue
            pep = seq[start:end]
            if must_residue not in pep:
                continue
            mass = WATER + sum(MONO.get(r, 0.0) for r in pep)
            if mass > max_mass:
                continue
            out.append({"protein": acc, "start": start + 1, "end": end, "seq": pep, "length": L,
                       "missed_cleavages": mc, "mass": mass})
    return out


def main():
    for p in (OUT_CAND, OUT_SUM, OUT_AUD):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    with open(PREREG_AUDIT, encoding="utf-8") as fh:
        pa = json.load(fh)
    if pa["protocol_sha256"] != sha(PROTOCOL):
        sys.exit("REFUSE: pre-registration protocol changed since it was registered")
    for rel, want in pa["declared_inputs_sha256"].items():
        fp = os.path.join(DATA_DIR, rel)
        if sha(fp) != want:
            sys.exit(f"REFUSE: declared input changed since registration: {rel}")
    with open(PROTOCOL, encoding="utf-8") as fh:
        P = json.load(fh)

    dg = P["digestion"]
    pep_rows = dread(os.path.join(DATA_DIR, P["dataset"]["files"]["peptides"]))
    pg_own = dread(os.path.join(DATA_DIR, P["dataset"]["files"]["proteinGroups_own_arm"]))
    pg_cross = dread(os.path.join(DATA_DIR, P["dataset"]["files"]["proteinGroups_cross_arm"]))
    fasta = read_fasta(os.path.join(DATA_DIR, P["dataset"]["files"]["fasta"]))

    universe = acc_set(pep_rows, "Leading razor protein") | acc_set(pg_cross, "Protein IDs")
    universe = sorted(a for a in universe if a in fasta)

    # ground truth: exact (protein, start, sequence) triples from the Global peptide table
    truth = set()
    for r in pep_rows:
        if r.get("Reverse") == "+" or r.get("Potential contaminant") == "+":
            continue
        m = re.match(r"^(sp|tr)\|([^|]+)\|", r["Leading razor protein"])
        if not m:
            continue
        try:
            truth.add((m.group(2), int(r["Start position"]), r["Sequence"]))
        except (KeyError, ValueError):
            continue

    # protein abundance from the Global arm's own quantification
    abund = {}
    for r in pg_own:
        if r.get("Reverse") == "+" or r.get("Potential contaminant") == "+":
            continue
        try:
            inten = float(r.get("Intensity") or 0)
        except ValueError:
            inten = 0.0
        for tok in r["Majority protein IDs"].split(";"):
            mm = re.match(r"^(sp|tr)\|([^|]+)\|", tok.strip())
            if mm and mm.group(2) not in abund:
                abund[mm.group(2)] = inten

    abe_ids = acc_set(pg_cross, "Protein IDs")

    mc_set = set(dg["missed_cleavages_generated"])
    rows = []
    for acc in universe:
        seq = fasta[acc]
        cands = candidates_for_protein(acc, seq, mc_set, dg["min_length_generated"],
                                       dg["max_length_generation_cap"], dg["max_monoisotopic_mass_da"], dg["must_contain_residue"])
        for c in cands:
            pep = c["seq"]
            n = len(pep)
            gravy = sum(KD.get(r, 0.0) for r in pep) / n
            frac_basic = sum(pep.count(r) for r in "KRH") / n
            frac_acidic = sum(pep.count(r) for r in "DE") / n
            n_cys = pep.count("C")
            dist_n = c["start"] - 1
            dist_c = len(seq) - c["end"]
            a = abund.get(acc)
            log_abund = math.log10(1 + a) if a is not None else 0.0
            abund_missing = int(a is None)
            abe_id = int(acc in abe_ids)
            detected = int((acc, c["start"], pep) in truth)
            rows.append({"protein": acc, "start": c["start"], "end": c["end"], "sequence": pep,
                        "length": n, "missed_cleavages": c["missed_cleavages"], "mass": round(c["mass"], 4),
                        "gravy": round(gravy, 5), "frac_basic": round(frac_basic, 5), "frac_acidic": round(frac_acidic, 5),
                        "dist_n_term": dist_n, "dist_c_term": dist_c, "n_cys": n_cys,
                        "log_abundance": round(log_abund, 5), "abundance_missing": abund_missing,
                        "abe_identified": abe_id, "theoretical_detectable": int(7 <= n <= 30),
                        "empirical_detected": detected})

    with open(OUT_CAND, "w", newline="", encoding="utf-8") as fh:
        fields = list(rows[0].keys())
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # ---------------- analysis ----------------
    feat_names = P["features"]
    geom_names = P["cleavage_geometry_features"]
    X_full = np.array([[r[f] for f in feat_names] for r in rows], dtype=float)
    y = np.array([r["empirical_detected"] for r in rows], dtype=int)
    theo = np.array([r["theoretical_detectable"] for r in rows], dtype=int)
    proteins = np.array([r["protein"] for r in rows])
    length = np.array([r["length"] for r in rows])
    mc = np.array([r["missed_cleavages"] for r in rows])

    def auc_rank(scores, labels):
        order = np.argsort(scores, kind="mergesort")
        ranks = np.empty(len(scores))
        s_sorted = scores[order]
        i = 0
        rnk = 1
        while i < len(s_sorted):
            j = i
            while j < len(s_sorted) and s_sorted[j] == s_sorted[i]:
                j += 1
            avg = (rnk + (rnk + (j - i) - 1)) / 2
            ranks[order[i:j]] = avg
            rnk += (j - i)
            i = j
        n_pos, n_neg = int(labels.sum()), int((1 - labels).sum())
        u = ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2
        return float(u / (n_pos * n_neg))

    def cv_oof_auc(cols, seed_off):
        idx = [feat_names.index(c) for c in cols]
        Xc = X_full[:, idx]
        mu, sd = Xc.mean(0), Xc.std(0)
        sd[sd == 0] = 1.0
        Xs = (Xc - mu) / sd
        uniq_prot = sorted(set(proteins.tolist()))
        rng = np.random.default_rng(P["master_seed"] + seed_off)
        perm = rng.permutation(len(uniq_prot))
        folds_of_protein = {uniq_prot[perm[i]]: i % P["model"]["cv_folds"] for i in range(len(uniq_prot))}
        fold_id = np.array([folds_of_protein[p] for p in proteins])
        oof = np.zeros(len(y), dtype=float)
        per_fold_auc = []
        for k in range(P["model"]["cv_folds"]):
            te = fold_id == k
            tr = ~te
            if y[tr].sum() == 0 or (1 - y[tr]).sum() == 0:
                continue
            clf = LogisticRegression(max_iter=2000)
            clf.fit(Xs[tr], y[tr])
            p = clf.predict_proba(Xs[te])[:, 1]
            oof[te] = p
            if y[te].sum() > 0 and (1 - y[te]).sum() > 0:
                per_fold_auc.append(auc_rank(p, y[te]))
        return auc_rank(oof, y), per_fold_auc, oof, fold_id

    m1_auc = auc_rank(theo.astype(float), y)
    m2_auc, m2_fold, m2_oof, fold_id = cv_oof_auc(geom_names, seed_off=1)
    m3_auc, m3_fold, m3_oof, _ = cv_oof_auc(feat_names, seed_off=1)

    def cluster_boot_diff(vals_a, vals_b, prots, reps, seed):
        uniq = sorted(set(prots.tolist()))
        idx_by = {c: [] for c in uniq}
        for i, c in enumerate(prots):
            idx_by[c].append(i)
        rng = np.random.default_rng(seed)
        n = len(uniq)
        draws = rng.integers(0, n, size=(reps, n))
        cl = np.array(uniq, dtype=object)
        out = np.empty(reps)
        for r in range(reps):
            rows_idx = []
            for c in cl[draws[r]]:
                rows_idx.extend(idx_by[c])
            rows_idx = np.array(rows_idx)
            out[r] = vals_a(rows_idx) - vals_b(rows_idx)
        return out

    # D1: length-window direction
    def rate(mask_fn):
        def f(idx):
            m = mask_fn(idx)
            sub_y = y[idx][m]
            return sub_y.mean() if len(sub_y) else 0.0
        return f
    in_win = rate(lambda idx: (length[idx] >= 7) & (length[idx] <= 30))
    out_win = rate(lambda idx: length[idx] > 30)
    d1_point = in_win(np.arange(len(y))) - out_win(np.arange(len(y)))
    d1_draws = cluster_boot_diff(in_win, out_win, proteins, P["bootstrap_reps"], P["master_seed"] + 2)
    d1_lo, d1_hi = float(np.quantile(d1_draws, 0.025)), float(np.quantile(d1_draws, 0.975))

    # D2: missed-cleavage direction
    mc_rates = [float(y[mc == v].mean()) for v in (0, 1, 2)]
    def spearman3(a, b):
        def rk(v):
            order = np.argsort(v, kind="mergesort")
            r = np.empty(len(v))
            r[order] = np.arange(1, len(v) + 1)
            return r
        ra, rb = rk(np.array(a, float)), rk(np.array(b, float))
        return float(np.corrcoef(ra, rb)[0, 1])
    d2_rho = spearman3([0, 1, 2], mc_rates)

    # D3: informativeness, AUC(M3) - AUC(M2), per fold then averaged, bootstrapped over proteins on OOF preds
    def auc_from_oof(oof_vec):
        def f(idx):
            yy = y[idx]
            if yy.sum() == 0 or (1 - yy).sum() == 0:
                return 0.5
            return auc_rank(oof_vec[idx], yy)
        return f
    m2_f = auc_from_oof(m2_oof)
    m3_f = auc_from_oof(m3_oof)
    d3_point = m3_auc - m2_auc
    d3_draws = cluster_boot_diff(m3_f, m2_f, proteins, P["bootstrap_reps"], P["master_seed"] + 3)
    d3_lo, d3_hi = float(np.quantile(d3_draws, 0.025)), float(np.quantile(d3_draws, 0.975))

    # D4: theoretical-flag coefficient sign in a refit including it
    idx_all = [feat_names.index(c) for c in feat_names]
    X_d4 = np.hstack([X_full, theo.reshape(-1, 1)])
    mu, sd = X_d4.mean(0), X_d4.std(0)
    sd[sd == 0] = 1.0
    Xs_d4 = (X_d4 - mu) / sd
    clf_d4 = LogisticRegression(max_iter=2000)
    clf_d4.fit(Xs_d4, y)
    d4_coefs = dict(zip(feat_names + ["theoretical_detectable_flag"], clf_d4.coef_[0].tolist()))
    d4_theo_coef = d4_coefs["theoretical_detectable_flag"]

    pass_rule = (d1_lo > 0 and d2_rho <= -0.5 and d3_lo > 0 and d3_point >= 0.03 and d4_theo_coef > 0)

    summary = {
        "n_proteins_universe": len(universe), "n_candidates": len(rows), "n_positive": int(y.sum()),
        "n_theoretical_detectable": int(theo.sum()),
        "M1_theoretical_flag_auc": round(m1_auc, 4),
        "M2_geometry_only_auc": round(m2_auc, 4), "M2_per_fold_auc": [round(a, 4) for a in m2_fold],
        "M3_full_model_auc": round(m3_auc, 4), "M3_per_fold_auc": [round(a, 4) for a in m3_fold],
        "D1_length_direction": {"point": round(d1_point, 4), "ci95": [round(d1_lo, 4), round(d1_hi, 4)], "pass": bool(d1_lo > 0)},
        "D2_missed_cleavage_direction": {"rates_by_mc_0_1_2": [round(r, 4) for r in mc_rates], "spearman": round(d2_rho, 4), "pass": bool(d2_rho <= -0.5)},
        "D3_informativeness": {"point": round(d3_point, 4), "ci95": [round(d3_lo, 4), round(d3_hi, 4)], "pass": bool(d3_lo > 0 and d3_point >= 0.03)},
        "D4_theoretical_flag_coefficient": {"value": round(d4_theo_coef, 4), "pass": bool(d4_theo_coef > 0)},
        "M3_full_model_coefficients": {k: round(v, 4) for k, v in dict(zip(feat_names, clf_d4.coef_[0][:len(feat_names)].tolist())).items()},
        "pass": bool(pass_rule)
    }
    with open(OUT_SUM, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    audit = {"script": "scripts/run_phase3_empirical_detectability_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "protocol_sha256": sha(PROTOCOL), "prereg_audit_sha256": sha(PREREG_AUDIT),
             "summary": summary,
             "outputs": {os.path.relpath(OUT_CAND, ROOT): sha(OUT_CAND), os.path.relpath(OUT_SUM, ROOT): sha(OUT_SUM)}}
    with open(OUT_AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
