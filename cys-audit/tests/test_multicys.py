import unittest

import _helpers
from cys_audit.checks import multicys
from cys_audit.constants import FAIL, UNDECIDABLE, WARNING


class Detector(unittest.TestCase):
    def test_planted(self):
        ds, _ = _helpers.synth_dataset("multicys", 1.0)
        r = multicys.run(ds, _helpers.cfg(background="observed"))
        self.assertIn(r["status"], (FAIL, WARNING))
        self.assertGreater(r["details"]["share_multi_only_positive"], r["details"]["share_multi_only_background"])

    def test_null(self):
        ds, _ = _helpers.synth_dataset("multicys", 0.0)
        self.assertNotEqual(multicys.run(ds, _helpers.cfg(background="observed"))["status"], FAIL)

    def test_requires_columns(self):
        ds, _ = _helpers.synth_dataset("cleavage", 0.0)   # no n_cys_peptide column
        self.assertEqual(multicys.run(ds, _helpers.cfg())["status"], UNDECIDABLE)


if __name__ == "__main__":
    unittest.main()
