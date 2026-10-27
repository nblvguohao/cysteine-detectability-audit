import importlib.util
import os
import unittest

import _helpers

EX = os.path.join(os.path.dirname(_helpers.HERE), "examples", "convert_maxquant_abe.py")
spec = importlib.util.spec_from_file_location("conv_abe", EX)
conv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(conv)


class ProbabilityString(unittest.TestCase):
    def test_two_cys(self):
        self.assertEqual(conv.cys_probs("AAANFFSASC(0.114)VPC(0.886)ADQSSFPK", "AAANFFSASCVPCADQSSFPK"),
                         {9: 0.114, 12: 0.886})

    def test_one_annotated(self):
        self.assertEqual(conv.cys_probs("AAKALDKRQAHLC(1)VLASNCDEPMYVKLVE", "AAKALDKRQAHLCVLASNCDEPMYVKLVE"), {12: 1.0})

    def test_accession(self):
        self.assertEqual(conv.acc("sp|P63323|RS12_MOUSE"), "P63323")
        self.assertTrue(conv.bad_protein("CON__Q0IIK2"))


if __name__ == "__main__":
    unittest.main()
