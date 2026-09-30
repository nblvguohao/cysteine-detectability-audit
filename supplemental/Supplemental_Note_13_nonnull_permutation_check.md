# Supplemental Note 13. Is the normal approximation adequate where there is an effect?

The binary-feature statistic compares the observed number of feature-positive residues in the positive set
with its within-protein expectation, using a normal approximation (Experimental Procedures, "Site-selection
tests"). The validation reported in the main text was run on the polar STNQ class of the S-nitrosylation
cohort, where the permutation P of 0.7435 agreed with the analytic value — that is, in a class with **no
effect**. Agreement at P ≈ 0.74 does not show that the approximation holds at the effect sizes the statistic
is used to report (z = −2.61, +3.34, −8.78, +2.67 in Artifact 2; +0.67 to +15.72 in Artifact 1).

This note reports a check in the non-null regime.

**What was done.** Residues were generated protein by protein at the parameters of the two cohorts whose
binary statistics are reported (PXD072089: 5,004 cysteines in 581 proteins, 737 positives; PXD044043: 18,567
in 1,519, 3,206 positives). Residues per protein were drawn from a negative binomial matched to the cohort
mean at three dispersions (0.5, 1.0, 2.0); each residue carried the binary feature with a per-protein
probability drawn from 0.1 + 0.9 Beta(3,1); it entered the positive set with probability q, or with
probability q times a bias factor when it carried the feature. Both the null (bias 1.0) and two biased
regimes (1.15 and 1.35) were run, 25 datasets per cell and 400 within-protein permutations per dataset.

Given the realised margins c_i (cysteines), k_i (feature-positive) and n_i (positives) per protein, the exact
within-protein permutation distribution of the statistic is X_i ~ Hypergeometric(c_i, k_i, n_i), because
permuting feature labels holds the positives fixed. The analytic statistic uses E = Σ n_i p_i and
V = Σ n_i p_i (1 − p_i) with p_i = k_i / c_i; the permutation statistic uses the finite-population-corrected
variance. The two two-sided P values were compared dataset by dataset.

**Result.** In the regime of the reported statistics the two agree closely, and the agreement improves with
the size of the effect:

| mean |z| | cells | mean absolute difference in P | differed by >0.01 | differed by >0.05 |
|---|---|---|---|---|
| 3.0 – 3.5 | 6 | ≤ 0.0053 | ≤ 12% of datasets | ≤ 4% |
| 6.6 – 6.9 | 3 | 0.0000 | 0% | 0% |
| null (bias 1.0) | 6 | ≤ 0.0203 | — | — |

Under the null the two differ by up to 0.02, but both values are above 0.9 there and no reading changes.
The residual disagreement is therefore concentrated where it cannot affect a conclusion.

**Limits.** This is a simulation at the cohorts' parameters, not a re-analysis of them: the per-protein
distribution of cysteines is not released with the manuscript and is modeled, at three dispersions, so the
result is a check on the approximation rather than a measurement of these cohorts. The reported statistics
were computed on the authors' analytical inputs, which this note does not reproduce. Effect sizes were
generated through a selection bias applied at the residue level, which is one of several ways an effect of
that magnitude could arise.

**Reproduction.** `analysis/nonnull_permutation_check.py` (fixed seed 20260927, 400 permutations, 25 datasets
per cell); stored output `Source_Data_text_nonnull_permutation_check.csv`.
