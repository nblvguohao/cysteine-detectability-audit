import unittest

import numpy as np

import _helpers
from cys_audit import stats
from cys_audit.checks import abundance
from cys_audit.constants import FAIL, UNDECIDABLE, WARNING


class AUC(unittest.TestCase):
    def test_known_values(self):
        self.assertAlmostEqual(stats.auc([1, 2, 3, 4], [0, 0, 1, 1]), 1.0)
        self.assertAlmostEqual(stats.auc([1, 2, 3, 4], [1, 1, 0, 0]), 0.0)
        self.assertAlmostEqual(stats.auc([1, 1, 1, 1], [1, 0, 1, 0]), 0.5)
        self.assertAlmostEqual(stats.auc([1, 2, 2, 3], [0, 1, 0, 1]), 0.875)


class Detector(unittest.TestCase):
    def test_planted(self):
        ds, _ = _helpers.synth_dataset("abundance", 1.0)
        r = abundance.run(ds, _helpers.cfg())
        self.assertIn(r["status"], (FAIL, WARNING))
        self.assertGreater(r["estimate"], 0.5)

    def test_null(self):
        ds, _ = _helpers.synth_dataset("abundance", 0.0)
        self.assertNotEqual(abundance.run(ds, _helpers.cfg())["status"], FAIL)

    def test_low_coverage_undecidable(self):
        ds, _ = _helpers.synth_dataset("abundance", 1.0)
        ds.abundance = ds.abundance.copy()
        ds.abundance[: int(0.5 * len(ds))] = np.nan
        self.assertEqual(abundance.run(ds, _helpers.cfg())["status"], UNDECIDABLE)


if __name__ == "__main__":
    unittest.main()
