import os
import tempfile
import unittest

import _helpers
from cys_audit.io import InputError, read_sites


class Validation(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        with open(os.path.join(self.d, "p.fasta"), "w") as fh:
            fh.write(">sp|P1|X_HUMAN\nMACKCAAAC\n>P2\nMCC\n")

    def path(self, rows, header=("protein", "position", "label")):
        p = os.path.join(self.d, "s.tsv")
        _helpers.write_tsv(p, header, rows)
        return p

    def test_missing_column(self):
        with self.assertRaises(InputError):
            read_sites(self.path([["P1", 3]], header=("protein", "position")))

    def test_duplicate(self):
        with self.assertRaises(InputError):
            read_sites(self.path([["P1", 3, 1], ["P1", 3, 0]]))

    def test_bad_label(self):
        with self.assertRaises(InputError):
            read_sites(self.path([["P1", 3, 2]]))

    def test_fasta_residue_check(self):
        with self.assertRaises(InputError):
            read_sites(self.path([["P1", 2, 1], ["P1", 4, 0]]), fasta=os.path.join(self.d, "p.fasta"))
        ds = read_sites(self.path([["P1", 3, 1], ["P1", 5, 0]]), fasta=os.path.join(self.d, "p.fasta"))
        self.assertEqual(len(ds), 2)

    def test_expand_background(self):
        ds = read_sites(self.path([["P1", 3, 1], ["P2", 2, 1]]), fasta=os.path.join(self.d, "p.fasta"),
                        expand_background=True)
        self.assertEqual(sorted(zip(ds.protein, ds.position.tolist())),
                         [("P1", 3), ("P1", 5), ("P1", 9), ("P2", 2), ("P2", 3)])
        self.assertEqual(int(ds.label.sum()), 2)

    def test_abundance_must_be_protein_level(self):
        with self.assertRaises(InputError):
            read_sites(self.path([["P1", 3, 1, 5.0], ["P1", 5, 0, 6.0]],
                                 header=("protein", "position", "label", "abundance")))


if __name__ == "__main__":
    unittest.main()
