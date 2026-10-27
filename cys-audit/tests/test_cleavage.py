import unittest

import _helpers
from cys_audit import proteases
from cys_audit.checks import cleavage
from cys_audit.constants import FAIL, UNDECIDABLE, WARNING


class Rules(unittest.TestCase):
    def test_trypsin_proline(self):
        m = proteases.competent_mask("AKPAKA", "trypsin")
        self.assertEqual(m, [False, False, False, False, True, False])
        m = proteases.competent_mask("AKPAKA", "trypsin/p")
        self.assertEqual(m, [False, True, False, False, True, False])

    def test_band(self):
        s = "AAKCAAAAAAAA"          # K at 3, C at 4 -> distance 1
        self.assertTrue(proteases.band_flag(s, 4, (1, 3), "trypsin"))
        self.assertFalse(proteases.band_flag(s, 4, (6, 12), "trypsin"))
        self.assertTrue(proteases.band_flag("CAAAAAAD", 1, (6, 12), "aspn"))

    def test_digest(self):
        peps = proteases.digest("AAAAAAAKBBBBBBBBRPCCCCCCCK", "trypsin", missed=0, min_len=1, max_len=50)
        self.assertEqual(peps, [(1, 8), (9, 26)])


class Detector(unittest.TestCase):
    def test_planted_artefact_detected(self):
        ds, _ = _helpers.synth_dataset("cleavage", 1.0)
        r = cleavage.run(ds, _helpers.cfg())
        self.assertIn(r["status"], (FAIL, WARNING))
        self.assertGreater(r["details"]["bands"]["proximal_1_3"]["estimate"], 0.5)

    def test_null_not_flagged(self):
        ds, _ = _helpers.synth_dataset("cleavage", 0.0)
        r = cleavage.run(ds, _helpers.cfg())
        self.assertNotIn(r["status"], (FAIL,))

    def test_needs_fasta(self):
        ds, _ = _helpers.synth_dataset("cleavage", 1.0)
        ds.fasta = None
        self.assertEqual(cleavage.run(ds, _helpers.cfg())["status"], UNDECIDABLE)

    def test_wrong_protease_rule_weaker(self):
        ds, _ = _helpers.synth_dataset("cleavage", 1.0)
        t = cleavage.run(ds, _helpers.cfg())["details"]["bands"]["proximal_1_3"]["estimate"]
        g = cleavage.run(ds, _helpers.cfg(protease="gluc"))["details"]["bands"]["proximal_1_3"]["estimate"]
        self.assertGreater(abs(t), abs(g))


if __name__ == "__main__":
    unittest.main()
