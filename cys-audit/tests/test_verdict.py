import os
import tempfile
import unittest

import numpy as np

import _helpers
from cys_audit import proteases, synthetic
from cys_audit.constants import CLAIM_FAIL, CLAIM_NOT_EVALUATED, CLAIM_PASS, CLAIM_UNDECIDABLE
from cys_audit.io import read_sites
from cys_audit.verdict import VERDICT_TO_CLAIM, claim_retest, claim_verdict

# All 13 return paths, copied from VERDICT_BRANCHES in the analysis tree's scripts/ptm_detectability_diagnostics.py
BRANCHES = [
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.4, 1.2), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.8, "precheck_blocking": True}, "undecidable"),
    ({"baseline": None, "baseline_ci": None, "controlled_ci": None, "random_control_ci": None,
      "controlled_point": None}, "out_of_instrument_scope"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.4, 1.2), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.8, "claim_direction": "null_no_preference"}, "undecidable"),
    ({"baseline": 0.1, "baseline_ci": (-0.3, 0.5), "controlled_ci": (0.2, 0.9), "random_control_ci": (-0.2, 0.4),
      "controlled_point": 0.55, "claim_direction": "null_no_preference"}, "reverses"),
    ({"baseline": 0.1, "baseline_ci": (-0.3, 0.5), "controlled_ci": (-0.2, 0.4), "random_control_ci": (-0.2, 0.4),
      "controlled_point": 0.1, "claim_direction": "null_no_preference"}, "survives"),
    ({"baseline": -1.0, "baseline_ci": (-1.5, -0.5), "controlled_ci": (-1.2, -0.4), "random_control_ci": (-1.5, -0.5),
      "controlled_point": -0.8}, "baseline_contradicts_claim"),
    ({"baseline": 0.2, "baseline_ci": (-0.3, 0.7), "controlled_ci": (-0.4, 0.6), "random_control_ci": (-0.4, 0.6),
      "controlled_point": 0.1, "baseline_reproduced": False}, "undecidable"),
    ({"baseline": 0.2, "baseline_ci": (-0.3, 0.7), "controlled_ci": (-0.4, 0.6), "random_control_ci": (-0.4, 0.6),
      "controlled_point": 0.1, "baseline_reproduced": True}, "undecidable"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (-1.2, -0.4), "random_control_ci": (0.5, 1.5),
      "controlled_point": -0.8}, "reverses"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.1, 0.4), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.3}, "attenuated"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.6, 1.3), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.9}, "survives"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (-0.1, 0.5), "random_control_ci": (0.6, 1.4),
      "controlled_point": 0.2}, "vanishes"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (-0.1, 0.5), "random_control_ci": (-0.2, 1.4),
      "controlled_point": 0.2}, "undecidable"),
]


class VerdictRule(unittest.TestCase):
    def test_all_13_paths(self):
        reasons = set()
        for kw, exp in BRANCHES:
            v = claim_verdict(**kw)
            self.assertEqual(v["verdict"], exp, kw)
            reasons.add(v["reason"])
        self.assertEqual(len(reasons), 13)

    def test_mapping_total(self):
        for _, exp in BRANCHES:
            self.assertIn(exp, VERDICT_TO_CLAIM)


def _claim_dataset(artefact, strength, feature_fn, seed=5):
    seqs, cols = synthetic.generate(artefact, strength, seed, n_proteins=500)
    rng = np.random.default_rng(seed + 1)
    cols["feature"] = feature_fn(seqs, cols, rng)
    d = tempfile.mkdtemp()
    synthetic.write(seqs, cols, d)
    return read_sites(os.path.join(d, "sites.tsv"), fasta=os.path.join(d, "proteome.fasta"))


class ClaimRetest(unittest.TestCase):
    def test_not_evaluated_without_feature(self):
        ds, _ = _helpers.synth_dataset("none", 0.0)
        self.assertEqual(claim_retest(ds, _helpers.cfg())["status"], CLAIM_NOT_EVALUATED)

    def test_artefactual_claim_does_not_pass(self):
        # the "feature" IS the cleavage flag, and labels were planted on it: the control must remove it
        def flag(seqs, cols, rng):
            masks = {p: proteases.competent_mask(seqs[p], "trypsin") for p in seqs}
            return np.array([int(proteases.band_flag(seqs[p], int(q), (1, 3), "trypsin", masks[p]))
                             for p, q in zip(cols["protein"], cols["position"])])
        ds = _claim_dataset("cleavage", 1.0, flag)
        r = claim_retest(ds, _helpers.cfg(claim_feature="feature"))
        self.assertIn(r["status"], (CLAIM_FAIL, CLAIM_UNDECIDABLE))

    def test_genuine_claim_passes(self):
        # feature enriched among positives independently of every detectability covariate
        def genuine(seqs, cols, rng):
            lab = np.asarray(cols["label"])
            return (rng.random(len(lab)) < np.where(lab == 1, 0.7, 0.3)).astype(int)
        ds = _claim_dataset("none", 0.0, genuine)
        r = claim_retest(ds, _helpers.cfg(claim_feature="feature"))
        self.assertEqual(r["status"], CLAIM_PASS, r)


if __name__ == "__main__":
    unittest.main()
