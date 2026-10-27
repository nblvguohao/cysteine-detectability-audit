import os
import tempfile
import unittest

import _helpers
from cys_audit.checks import overlap
from cys_audit.constants import FAIL, PASS, UNDECIDABLE
from cys_audit.io import read_sites


class Detector(unittest.TestCase):
    def test_full_depletion_fails(self):
        ds, _ = _helpers.synth_dataset("overlap", 1.0)
        self.assertEqual(overlap.run(ds, _helpers.cfg())["status"], FAIL)

    def test_no_depletion_passes(self):
        ds, _ = _helpers.synth_dataset("overlap", 0.0)
        self.assertEqual(overlap.run(ds, _helpers.cfg())["status"], PASS)

    def test_exact_counts(self):
        d = tempfile.mkdtemp()
        rows = [["P1", 1, 1, 1], ["P1", 5, 1, 1], ["P2", 3, 0, 1], ["P2", 9, 0, 0], ["P3", 2, 1, 1]]
        _helpers.write_tsv(os.path.join(d, "s.tsv"), ["protein", "position", "label", "detected"], rows)
        ds = read_sites(os.path.join(d, "s.tsv"))
        r = overlap.run(ds, _helpers.cfg())
        self.assertEqual(r["details"]["detected"], 4)
        self.assertEqual(r["details"]["positive_and_detected"], 3)
        self.assertAlmostEqual(r["estimate"], 0.75)

    def test_no_detected_column(self):
        ds, _ = _helpers.synth_dataset("cleavage", 0.0)
        ds.detected = None
        self.assertEqual(overlap.run(ds, _helpers.cfg())["status"], UNDECIDABLE)


if __name__ == "__main__":
    unittest.main()
