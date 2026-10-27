import os
import tempfile
import unittest

import numpy as np

import _helpers
from cys_audit import synthetic
from cys_audit.constants import CLAIM_PASS, CLAIM_UNDECIDABLE
from cys_audit.io import read_sites
from cys_audit.propensity import fit_logit, prepare, strata_from_logit
from cys_audit.verdict import claim_retest


class Fit(unittest.TestCase):
    def test_recovers_sign_and_scale(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=(4000, 2))
        eta = 1.5 * x[:, 0] - 1.0 * x[:, 1] - 0.5
        y = (rng.random(4000) < 1 / (1 + np.exp(-eta))).astype(float)
        w, fitted = fit_logit(x, y)
        self.assertAlmostEqual(w[1], 1.5, delta=0.15)
        self.assertAlmostEqual(w[2], -1.0, delta=0.15)
        self.assertAlmostEqual(w[0], -0.5, delta=0.15)
        self.assertEqual(fitted.shape, (4000,))

    def test_prepare_imputes_and_drops(self):
        X, kept, notes = prepare({"a": np.array([1.0, np.nan, 3.0]), "b": np.array([2.0, 2.0, 2.0])})
        self.assertEqual(kept, ["a"])
        self.assertEqual(X.shape, (3, 1))
        self.assertTrue(any("median" in n for n in notes) and any("constant" in n for n in notes))

    def test_quintiles(self):
        s = strata_from_logit(np.arange(100.0))
        self.assertEqual(sorted(set(s.tolist())), [0, 1, 2, 3, 4])
        self.assertEqual(int((s == 0).sum()), 20)


def _confounded(genuine, seed=21):
    """z drives both the label and the feature; the genuine variant adds a direct label->feature link."""
    seqs, cols = synthetic.generate("none", 0.0, seed, n_proteins=500)
    rng = np.random.default_rng(seed + 1)
    n = len(cols["protein"])
    z = rng.normal(size=n)
    lab = (rng.random(n) < 1 / (1 + np.exp(-(2.0 * z - 1.8)))).astype(int)
    feat_eta = 1.6 * z - 0.2 + (1.2 * (lab - 0.5) if genuine else 0.0)
    feat = (rng.random(n) < 1 / (1 + np.exp(-feat_eta))).astype(int)
    cols = {"protein": cols["protein"], "position": cols["position"], "label": lab, "detected": np.ones(n, dtype=int),
            "feature": feat, "z": z}
    d = tempfile.mkdtemp()
    synthetic.write(seqs, cols, d)
    return read_sites(os.path.join(d, "sites.tsv"), fasta=os.path.join(d, "proteome.fasta"))


class ClaimControl(unittest.TestCase):
    def test_confounded_claim_passes_without_but_not_with_covariate(self):
        ds = _confounded(genuine=False)
        weak = claim_retest(ds, _helpers.cfg(claim_feature="feature"))
        strong = claim_retest(ds, _helpers.cfg(claim_feature="feature", claim_covariates=["z"]))
        self.assertEqual(weak["status"], CLAIM_PASS)        # the band flags cannot see z
        self.assertNotEqual(strong["status"], CLAIM_PASS)   # the propensity control removes it
        self.assertIn("propensity quintile", strong["details"]["strata_from"][0])

    def test_genuine_claim_survives_with_covariate(self):
        ds = _confounded(genuine=True)
        r = claim_retest(ds, _helpers.cfg(claim_feature="feature", claim_covariates=["z"]))
        self.assertEqual(r["status"], CLAIM_PASS, r)

    def test_bad_covariates(self):
        ds = _confounded(genuine=False)
        self.assertEqual(claim_retest(ds, _helpers.cfg(claim_feature="feature", claim_covariates=["nope"]))["status"],
                         CLAIM_UNDECIDABLE)
        self.assertEqual(claim_retest(ds, _helpers.cfg(claim_feature="feature", claim_covariates=["feature"]))["status"],
                         CLAIM_UNDECIDABLE)


if __name__ == "__main__":
    unittest.main()
