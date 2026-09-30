"""Non-null validation of the within-protein normal approximation used for binary features.

The manuscript's binary-feature statistic compares the observed count of feature-positive
residues in the positive set with its within-protein expectation, using a normal
approximation (Experimental Procedures, "Site-selection tests"). The validation reported in
the manuscript used one class of one dataset in which there was no effect
(permutation P = 0.7435). This script asks whether the same approximation holds where there
IS an effect, at the parameters of the two cohorts the binary statistics are reported for.

It is a simulation at those parameters, NOT a re-analysis of the cohorts: the per-protein
distribution of cysteines is not released with the manuscript, so it is modelled as a
negative binomial matched to the cohort's mean, at three dispersions.

Design.  Residues are generated protein by protein. Residue j of protein i carries the
binary feature with probability p_i ~ 0.1 + 0.9 Beta(3,1); it enters the positive set with
probability q (no bias) or with probability q*bias if it carries the feature (bias). Given
the realised margins c_i (cysteines), k_i (feature-positive) and n_i (positives), the exact
within-protein permutation distribution of the statistic is
    X_i ~ Hypergeometric(c_i, k_i, n_i),
because permutation of the feature labels holds the positives fixed.

Statistics compared.  Analytic: z = (sum X_i - E) / sqrt(V) with
E = sum n_i p_i and V = sum n_i p_i (1 - p_i), p_i = k_i / c_i, and a two-sided normal P.
Permutation: the two-sided P of the same z over 400 within-protein permutations, using the
finite-population-corrected variance of the hypergeometric.
"""
import csv
import numpy as np
from math import erfc, sqrt

RNG = np.random.default_rng(20260927)
N_PERM = 400
N_REP = 25
BIASES = [1.0, 1.15, 1.35]          # 1.0 = null; >1 = positives favour the feature
DISPERSIONS = [0.5, 1.0, 2.0]
COHORTS = {'rice PXD072089': (5004, 737, 640), 'human PXD044043': (18567, 3206, 1519)}


def simulate(rng, n_prot, mean_cys, size, q, bias):
    """One dataset: returns c_i, k_i, n_i, X_i per protein."""
    p_nb = size / (size + mean_cys)
    c = np.maximum(rng.negative_binomial(size, p_nb, n_prot) + 1, 2)
    idx = np.repeat(np.arange(n_prot), c)
    p_feat = np.clip(0.1 + 0.9 * rng.beta(3, 1, size=n_prot), 0.02, 0.98)
    feat = rng.random(idx.size) < p_feat[idx]
    q_i = np.where(feat, q * bias, q)
    pos = rng.random(idx.size) < np.minimum(q_i, 1.0)
    c_i = np.bincount(idx, minlength=n_prot).astype(np.int64)
    k_i = np.bincount(idx[feat], minlength=n_prot).astype(np.int64)
    n_i = np.bincount(idx[pos], minlength=n_prot).astype(np.int64)
    X_i = np.bincount(idx[pos & feat], minlength=n_prot).astype(np.int64)
    return c_i, k_i, n_i, X_i


def stats(c, k, n, O):
    p = np.divide(k, c, out=np.zeros_like(k, dtype=float), where=c > 0)
    E = float((n * p).sum())
    V = float((n * p * (1 - p)).sum())
    return ((O - E) / sqrt(V), E, V) if V > 0 else (0.0, E, V)


def permutation_p(rng, c, k, n, z_obs, B):
    draws = rng.hypergeometric(k[None, :], (c - k)[None, :], n[None, :], size=(B, c.size))
    O = draws.sum(axis=1)
    p = np.divide(k, c, out=np.zeros_like(k, dtype=float), where=c > 0)
    fpc = np.divide(c - n, c - 1, out=np.zeros_like(k, dtype=float), where=c > 1)
    E = float((n * p).sum())
    V = float((n * p * (1 - p) * fpc).sum())
    if V <= 0:
        return 1.0, 0.0
    z = (O - E) / sqrt(V)
    return float((np.abs(z) >= abs(z_obs)).mean()), float(z.mean())


rows = []
print(f'{N_REP} datasets x {N_PERM} within-protein permutations per cell\n')
print(f"{'cohort':<20s} {'disp':>5s} {'bias':>5s} {'z(analytic)':>12s} {'z(perm)':>9s} "
      f"{'mean|dP|':>9s} {'P>0.01':>7s} {'P>0.05':>7s}")
for cname, (res, pos, prot) in COHORTS.items():
    q0, mean_cys = pos / res, res / prot
    for size in DISPERSIONS:
        for bias in BIASES:
            dP, za_, zp_ = [], [], []
            for _ in range(N_REP):
                c, k, n, X = simulate(RNG, prot, mean_cys, size, q0, bias)
                za, _, _ = stats(c, k, n, X.sum())
                pp, zp = permutation_p(RNG, c, k, n, za, N_PERM)
                pa = float(erfc(abs(za) / sqrt(2)))
                dP.append(abs(pa - pp)); za_.append(za); zp_.append(zp)
            dP = np.array(dP)
            r = dict(cohort=cname, dispersion=size, bias=bias,
                     mean_z_analytic=float(np.mean(za_)), mean_z_permutation=float(np.mean(zp_)),
                     mean_abs_P_diff=float(dP.mean()),
                     share_diff_gt_0p01=float((dP > 0.01).mean()),
                     share_diff_gt_0p05=float((dP > 0.05).mean()))
            rows.append(r)
            print(f"{cname:<20s} {size:>5} {bias:>5} {r['mean_z_analytic']:>12.2f} "
                  f"{r['mean_z_permutation']:>9.2f} {r['mean_abs_P_diff']:>9.4f} "
                  f"{r['share_diff_gt_0p01']:>7.2f} {r['share_diff_gt_0p05']:>7.2f}")

with open('analysis/nonnull_permutation_check.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print('\nwrote analysis/nonnull_permutation_check.csv')
