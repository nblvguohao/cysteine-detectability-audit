import unittest

import _helpers  # noqa: F401
from cys_audit.constants import FAIL, PASS, UNDECIDABLE, WARNING
from cys_audit.status import status_from_interval


class StatusRule(unittest.TestCase):
    def test_grid(self):
        d = 0.5
        cases = [((0.05, -0.2, 0.3), PASS), ((0.1, -0.6, 0.7), UNDECIDABLE), ((1.2, 0.8, 1.6), FAIL),
                 ((-1.2, -1.6, -0.6), FAIL), ((0.6, 0.2, 1.0), WARNING), ((0.3, 0.1, 0.4), WARNING),
                 ((-0.3, -0.45, -0.01), WARNING), ((0.5, 0.5, 0.9), FAIL), ((0.0, -0.5, 0.2), UNDECIDABLE)]
        for (pt, lo, hi), exp in cases:
            self.assertEqual(status_from_interval(pt, lo, hi, 0.0, d)[0], exp, (pt, lo, hi))

    def test_nonfinite(self):
        self.assertEqual(status_from_interval(float("nan"), 0, 1, 0, 0.5)[0], UNDECIDABLE)
        self.assertEqual(status_from_interval(None, 0, 1, 0, 0.5)[0], UNDECIDABLE)

    def test_auc_null(self):
        self.assertEqual(status_from_interval(0.51, 0.47, 0.54, 0.5, 0.05)[0], PASS)
        self.assertEqual(status_from_interval(0.70, 0.66, 0.74, 0.5, 0.05)[0], FAIL)


if __name__ == "__main__":
    unittest.main()
