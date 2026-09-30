import unittest

import numpy as np

import _helpers
from cys_audit.checks import searchspace
from cys_audit.constants import FAIL, PASS, UNDECIDABLE, WARNING


def run(**kw):
    return searchspace.run(None, _helpers.cfg(**kw))


class SearchSpace(unittest.TestCase):
    def test_delta_and_resolving_mass(self):
        d = abs(searchspace.DELTAS["persulfidation"] - searchspace.DELTAS["dioxidation"])
        self.assertAlmostEqual(d, 0.017758, places=6)
        self.assertAlmostEqual(searchspace.resolving_mass(d, 4.5), 3946.24, places=2)  # = 0.017758 Da / 4.5 ppm (tree CLAUDE.md s.9.71.4)

    def test_rules(self):
        self.assertEqual(run(search_mods=["Sulfide"])["status"], FAIL)                       # absent, 4.5 ppm, up to 4600 Da
        self.assertEqual(run(search_mods=["Sulfide"], precursor_ppm=2.0)["status"], WARNING)   # absent, M* 8879 > range
        self.assertEqual(run(search_mods=["Sulfide", "Dioxidation"])["status"], WARNING)      # present, not resolvable above M*
        self.assertEqual(run(search_mods=["Sulfide", "Dioxidation"], precursor_ppm=2.0)["status"], PASS)
        self.assertEqual(run(identity_readout="chemistry_inferred")["status"], UNDECIDABLE)
        self.assertEqual(run(precursor_ppm=None)["status"], UNDECIDABLE)
        self.assertEqual(run(modification="S-nitrosylation")["status"], UNDECIDABLE)

    def test_observed_masses(self):
        ds, _ = _helpers.synth_dataset("overlap", 0.0, n_proteins=50)
        ds.peptide_mass = np.full(len(ds), 1500.0)
        self.assertEqual(searchspace.run(ds, _helpers.cfg())["status"], WARNING)   # all below M*: protected by tolerance only
        ds.peptide_mass[ds.label == 1] = 4200.0
        self.assertEqual(searchspace.run(ds, _helpers.cfg())["status"], FAIL)

    def test_aliases(self):
        self.assertEqual(searchspace.canonical("Sulfide (C)"), "persulfidation")
        self.assertEqual(searchspace.canonical("Dioxidation (C)"), "dioxidation")
        self.assertEqual(searchspace.canonical("sulfinylation"), "dioxidation")


if __name__ == "__main__":
    unittest.main()
