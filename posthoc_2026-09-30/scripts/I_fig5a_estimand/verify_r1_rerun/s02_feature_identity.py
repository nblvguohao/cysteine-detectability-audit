"""Step 2 - what the DIG25 complement (CPL15) encodes about peptide visibility.

POST HOC revision analysis (2026-09-30), item I_fig5a_estimand. Not registered, not pre-specified.

The feature code is imported unchanged from the original repository
(scripts/run_cross_protease_detectability_probe.py; its LF-normalised sha256 is checked against
the value recorded in results/ptm_detectability_census_audit.json). The checks are properties
of the feature DEFINITIONS, so any protein sequences exercise them. Two sequence sources are
used, neither of which is a Figure 5a cohort:
  (1) synthetic sequences, seed 20260930, residue frequencies uniform over 20 amino acids;
  (2) the first 3,000 entries of the mouse reference proteome in external/ (read-only),
      used only to give realistic agreement rates for the near-identities.
Identities are asserted exactly; near-identities are reported as agreement fractions.

Run:  python -B s02_feature_identity.py
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import sys

import numpy as np

import common as C

sys.dont_write_bytecode = True  # never write __pycache__ into the read-only repository


def load_probe():
    spec = importlib.util.spec_from_file_location("probe_readonly", C.PROBE_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_fasta_gz(path, limit):
    seqs, name, buf = [], None, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                    if len(seqs) >= limit:
                        return seqs
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and len(seqs) < limit:
        seqs.append((name, "".join(buf)))
    return seqs


def synthetic(n, seed):
    rng = np.random.default_rng(seed)
    aa = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    out = []
    for i in range(n):
        L = int(rng.integers(60, 1500))
        s = "".join(rng.choice(aa, size=L))
        out.append((f"synthetic_{i}", s))
    return out


def feature_matrix(probe, seqs):
    rule = probe.PROTEASES["Trypsin"]
    rows = []
    for _, s in seqs:
        if "C" not in s:
            continue
        b = probe.boundaries_for(s, rule)
        for i, ch in enumerate(s):
            if ch == "C":
                rows.append(probe.features(s, i, b, rule))
    return np.asarray(rows, dtype=np.float64)


def checks(probe, X):
    f = {n: X[:, k] for k, n in enumerate(probe.FEATURE_NAMES)}
    Lrec = f["pep_dist_to_peptide_n_term"] + f["pep_dist_to_peptide_c_term"] + 1.0
    inwin = lambda L: (L >= probe.DETECTABLE_LENGTH[0]) & (L <= probe.DETECTABLE_LENGTH[1])
    mc1_len_rule = (inwin(f["pep_mc1_left_len"]) | inwin(f["pep_mc1_right_len"])).astype(float)
    mc2_len_rule = (inwin(f["pep_mc2_left_len"]) | inwin(f["pep_mc2_right_len"])).astype(float)
    any_len_rule = np.maximum.reduce([inwin(Lrec).astype(float), mc1_len_rule, mc2_len_rule])
    from scipy.stats import spearmanr, pearsonr
    out = {
        "n_cysteines": int(X.shape[0]),
        "exact_identities": {
            "pep_offset_in_peptide == pep_dist_to_peptide_n_term (duplicate column)":
                bool(np.array_equal(f["pep_offset_in_peptide"], f["pep_dist_to_peptide_n_term"])),
            "pep_len == pep_dist_to_peptide_n_term + pep_dist_to_peptide_c_term + 1":
                bool(np.array_equal(f["pep_len"], Lrec)),
            "pep_log_len == log1p(reconstructed length)": bool(np.allclose(f["pep_log_len"], np.log1p(Lrec), rtol=0, atol=1e-12)),
            "pep_detectable_length == 1[7 <= reconstructed length <= 30]":
                bool(np.array_equal(f["pep_detectable_length"], inwin(Lrec).astype(float))),
        },
        "near_identities_agreement_fraction": {
            "pep_detectable_both vs length-window flag from CPL15": float(np.mean(f["pep_detectable_both"] == inwin(Lrec))),
            "pep_detectable_mass vs length-window flag from CPL15": float(np.mean(f["pep_detectable_mass"] == inwin(Lrec))),
            "pep_mc1_detectable vs length-only rule on CPL15 mc1 lengths": float(np.mean(f["pep_mc1_detectable"] == mc1_len_rule)),
            "pep_mc2_detectable vs length-only rule on CPL15 mc2 lengths": float(np.mean(f["pep_mc2_detectable"] == mc2_len_rule)),
            "pep_detectable_any_missed_cleavage vs length-only rule on CPL15": float(np.mean(f["pep_detectable_any_missed_cleavage"] == any_len_rule)),
        },
        "correlations": {
            "spearman(pep_mass, reconstructed length)": float(spearmanr(f["pep_mass"], Lrec).statistic),
            "pearson(pep_mass, reconstructed length)": float(pearsonr(f["pep_mass"], Lrec).statistic),
        },
        "vis10_columns_exactly_determined_by_cpl15": ["pep_len", "pep_log_len", "pep_detectable_length"],
        "vis10_columns_not_determined_by_cpl15": ["pep_mass (near-determined by length)", "pep_gravy"],
    }
    return out


def main():
    raw_sha, lf_sha = C.sha256(C.PROBE_SCRIPT), C.sha256_lf(C.PROBE_SCRIPT)
    recorded = json.loads(C.CENSUS_AUDIT.read_text(encoding="utf-8"))["imported_feature_definitions"]["sha256"]
    probe = load_probe()
    vis = sorted(probe.VISIBILITY_ONLY)
    cpl = [n for n in probe.FEATURE_NAMES if n not in probe.VISIBILITY_ONLY]
    res = {
        "label": "POST HOC revision analysis 2026-09-30, item I_fig5a_estimand; feature-definition check",
        "probe_script": str(C.PROBE_SCRIPT), "probe_sha256_raw": raw_sha, "probe_sha256_lf_normalised": lf_sha,
        "probe_sha256_recorded_in_census_audit": recorded, "probe_matches_recorded_after_lf_normalisation": lf_sha == recorded,
        "n_features_dig25": len(probe.FEATURE_NAMES), "vis10": vis, "cpl15": cpl,
        "window_for_kr_count": "residues index-20..index+20 (41-residue window centred on the cysteine), K or R counted, proline blocking ignored",
        "seed_synthetic": C.SEED,
    }
    syn = synthetic(2000, C.SEED)
    res["synthetic"] = checks(probe, feature_matrix(probe, syn))
    mouse = read_fasta_gz(C.MOUSE_FASTA, 3000)
    res["mouse_first_3000_entries"] = checks(probe, feature_matrix(probe, mouse))
    res["mouse_first_entry"] = mouse[0][0]
    ok = all(res["synthetic"]["exact_identities"].values()) and all(res["mouse_first_3000_entries"]["exact_identities"].values())
    res["all_exact_identities_hold"] = bool(ok)
    (C.OUT / "feature_identity_check.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    C.record_inputs("s02_feature_identity", [C.PROBE_SCRIPT, C.CENSUS_AUDIT, C.MOUSE_FASTA],
                    {"probe_sha256_lf_normalised": lf_sha, "seed": C.SEED})
    print(json.dumps(res, indent=1)[:4000])


if __name__ == "__main__":
    main()
