"""The random control must not degenerate into the baseline.

Phase 5b found that under stratification, where no row is discarded, the "size-matched random control" was a
random subset of everything, i.e. everything, so it reproduced the baseline to every digit and an underpowered
result was read as an effect that disappeared under control. These tests are the two positive controls of
protocols/cys_audit_random_control_preregistration_2026-09-22.json: the new branch must fire when nothing is
dropped, and it must produce something that is NOT the baseline.
"""
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import _helpers  # noqa: E402
from cys_audit import proteases, synthetic  # noqa: E402
from cys_audit.constants import RANDOM_CONTROL_RETAINED_LIMIT  # noqa: E402
from cys_audit.io import read_sites  # noqa: E402
from cys_audit.verdict import claim_retest  # noqa: E402


def _dataset(artefact, strength, feature_fn, seed=5):
    seqs, cols = synthetic.generate(artefact, strength, seed, n_proteins=500)
    rng = np.random.default_rng(seed + 1)
    cols["feature"] = feature_fn(seqs, cols, rng)
    d = tempfile.mkdtemp()
    synthetic.write(seqs, cols, d)
    return read_sites(os.path.join(d, "sites.tsv"), fasta=os.path.join(d, "proteome.fasta"))


def _flag(seqs, cols, rng):
    masks = {p: proteases.competent_mask(seqs[p], "trypsin") for p in seqs}
    return np.array([int(proteases.band_flag(seqs[p], int(q), (1, 3), "trypsin", masks[p]))
                     for p, q in zip(cols["protein"], cols["position"])])


def _independent_feature(seqs, cols, rng):
    """A feature unrelated to the covariates, so every stratum holds both classes, nothing is discarded and
    the control retains the whole dataset - precisely the situation in which the old size-matched control
    became the baseline."""
    return rng.integers(0, 2, size=len(cols["protein"]))


class RandomControlDoesNotDegenerate(unittest.TestCase):
    def setUp(self):
        ds = _dataset("none", 0.0, _independent_feature)
        self.result = claim_retest(ds, _helpers.cfg(claim_feature="feature"))
        self.details = self.result.get("details", {})

    def test_random_control_is_reported_with_its_kind(self):
        rc = self.details.get("random_control")
        self.assertIsNotNone(rc, "a claim retest must say which random control it used")
        self.assertIn(rc["kind"], ("size_matched_random", "permuted_stratum"))
        self.assertIn("retained_fraction", rc)

    def test_branch_follows_the_declared_threshold(self):
        rc = self.details["random_control"]
        expected = ("size_matched_random" if rc["retained_fraction"] <= RANDOM_CONTROL_RETAINED_LIMIT
                    else "permuted_stratum")
        self.assertEqual(rc["kind"], expected)

    def test_control_is_not_the_baseline_when_nothing_was_dropped(self):
        """PC1: the failure this change exists to remove."""
        rc = self.details["random_control"]
        self.assertEqual(rc["kind"], "permuted_stratum", "fixture no longer exercises the new branch")
        self.assertEqual(rc["retained_fraction"], 1.0, "fixture must retain every row")
        baseline = self.details["baseline"]["estimate"]
        self.assertNotAlmostEqual(
            rc["estimate"], baseline, places=6,
            msg="the random control reproduced the baseline: it is not a control")

    def test_permuted_control_covers_zero(self):
        """Permuting the strata destroys the covariate-label relation, so the control's interval should
        include zero on a dataset with no planted artefact."""
        rc = self.details["random_control"]
        lo, hi = rc["ci"]
        self.assertLessEqual(lo, 0.0)
        self.assertGreaterEqual(hi, 0.0)

    def test_size_matched_branch_still_used_when_rows_are_dropped(self):
        """PC2: the unchanged branch. A planted cleavage artefact leaves orphaned strata, rows are dropped,
        and the size-matched control remains the right null."""
        ds = _dataset("cleavage", 1.0, _flag)
        rc = claim_retest(ds, _helpers.cfg(claim_feature="feature"))["details"]["random_control"]
        self.assertEqual(rc["kind"], "size_matched_random")
        self.assertLessEqual(rc["retained_fraction"], RANDOM_CONTROL_RETAINED_LIMIT)


if __name__ == "__main__":
    unittest.main()
