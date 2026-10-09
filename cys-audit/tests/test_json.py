import hashlib
import json
import os
import tempfile
import unittest

import _helpers
from cys_audit import cli, synthetic
from cys_audit.constants import STATUSES

REQUIRED_TOP = {"tool", "version", "dataset", "modification", "protease", "background", "seed", "input",
                "decision_constants", "tests", "dataset_summary", "claim", "overall_claim_status", "environment"}
REQUIRED_TEST = {"test", "status", "reason", "estimate", "ci", "n_positive", "n_background", "details"}


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def load(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


class EndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        synthetic.simulate("none", 0.0, 3, cls.d, n_proteins=200)
        cls.args = ["audit", "--input", os.path.join(cls.d, "sites.tsv"), "--fasta", os.path.join(cls.d, "proteome.fasta"),
                    "--protease", "trypsin", "--modification", "persulfidation", "--background", "observed",
                    "--precursor-ppm", "4.5", "--search-mods", "Sulfide", "--reps", "300", "--dataset", "unit"]

    def run_to(self, out):
        self.assertEqual(cli.main(self.args + ["--output", out]), 0)
        return load(os.path.join(out, "audit.json"))

    def test_schema_and_files(self):
        out = os.path.join(self.d, "r1")
        rec = self.run_to(out)
        self.assertTrue(REQUIRED_TOP <= set(rec))
        self.assertEqual(set(rec["tests"]), {"cleavage_geometry", "multi_cysteine", "protein_abundance",
                                             "positive_detected_overlap", "search_space"})
        for t in rec["tests"].values():
            self.assertTrue(REQUIRED_TEST <= set(t), t["test"])
            self.assertIn(t["status"], STATUSES)
        self.assertIn(rec["overall_claim_status"], {"PASS", "ATTENUATED", "UNDECIDABLE", "FAIL", "NOT_EVALUATED"})
        for f in ("audit.json", "audit.html", "summary.csv", "manifest.json"):
            self.assertTrue(os.path.exists(os.path.join(out, f)), f)
        man = load(os.path.join(out, "manifest.json"))
        for f, h in man["outputs"].items():
            self.assertEqual(sha(os.path.join(out, f)), h, f)
        # strict JSON: no NaN tokens
        with open(os.path.join(out, "audit.json"), encoding="utf-8") as fh:
            self.assertNotIn("NaN", fh.read())

    def test_deterministic(self):
        a = os.path.join(self.d, "ra")
        b = os.path.join(self.d, "rb")
        self.run_to(a)
        self.run_to(b)
        for f in ("audit.json", "summary.csv", "audit.html"):
            self.assertEqual(sha(os.path.join(a, f)), sha(os.path.join(b, f)), f)

    def test_bad_input_exit_code(self):
        p = os.path.join(self.d, "bad.tsv")
        _helpers.write_tsv(p, ["protein", "position"], [["P", 1]])
        self.assertEqual(cli.main(["audit", "--input", p, "--protease", "trypsin", "--output",
                                   os.path.join(self.d, "rbad")]), 2)


if __name__ == "__main__":
    unittest.main()
