import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "build_smultixcan_covariance",
    ROOT / "build_smultixcan_covariance.py",
)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class MatchingTests(unittest.TestCase):
    def test_multiallelic_exact_alt_and_orientation(self):
        model = {
            "rs1": {
                "specs": {("A", "G")},
                "tissues": {"T1", "T2"},
                "genes": {"g"},
                "dbs": set(),
            }
        }
        pvar = {
            "rs1": [
                {"chromosome": "1", "position": 10, "ref": "A", "alts": ("C", "G")}
            ]
        }
        resolved, report = M.resolve_targets(model, pvar)
        self.assertEqual(resolved["rs1"]["pvar_alt"], "G")
        self.assertEqual(resolved["rs1"]["selected_alt_index"], 2)
        self.assertEqual(resolved["rs1"]["reference_effect_allele"], "G")
        self.assertEqual(report[0]["multiallelic_reference"], "true")

    def test_split_multiallelic_rows_select_exact_model_alt(self):
        model = {
            "rs2": {
                "specs": {("A", "AC")},
                "tissues": {"T1"},
                "genes": {"g"},
                "dbs": set(),
            }
        }
        pvar = {
            "rs2": [
                {"chromosome": "1", "position": 20, "ref": "A", "alts": ("AC",)},
                {"chromosome": "1", "position": 20, "ref": "A", "alts": ("ACAC",)},
            ]
        }
        resolved, report = M.resolve_targets(model, pvar)
        self.assertEqual(resolved["rs2"]["pvar_alt"], "AC")
        self.assertEqual(report[0]["multiallelic_reference"], "true")

    def test_indel_requires_exact_sequence(self):
        self.assertEqual(M.match_model_alleles("A", "AT", "A", "AT"), ("exact", "AT"))
        self.assertIsNone(M.match_model_alleles("T", "AT", "A", "AT"))

    def test_palindrome_removed(self):
        self.assertEqual(M.match_model_alleles("A", "T", "A", "T"), ("palindrome_removed", None))

    def test_nonpalindromic_complement(self):
        self.assertEqual(M.match_model_alleles("T", "C", "A", "G"), ("complement", "G"))

    def test_effect_dosage_covariance_ddof_one(self):
        retained = {
            "rs1": np.array([0.0, 1.0, 2.0]),
            "rs2": np.array([2.0, 1.0, 0.0]),
        }
        rsids, matrix, eigenvalues = M.covariance_matrix(retained)
        self.assertEqual(rsids, ["rs1", "rs2"])
        np.testing.assert_allclose(matrix, [[1.0, -1.0], [-1.0, 1.0]])
        self.assertGreaterEqual(eigenvalues.min(), -1e-12)

    def test_vcf_multiallelic_counts_selected_effect_allele(self):
        with tempfile.TemporaryDirectory() as directory:
            vcf = Path(directory) / "x.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
                "1\t10\trs1\tA\tC,G\t.\t.\t.\tGT\t0/2\t1/2\t2/2\n"
            )
            resolved = {
                "rs1": {
                    "chromosome": "1",
                    "position": 10,
                    "pvar_ref": "A",
                    "pvar_alt": "G",
                    "reference_effect_allele": "G",
                }
            }
            samples, dosages = M.parse_effect_dosages(vcf, resolved)
            self.assertEqual(samples, ["s1", "s2", "s3"])
            np.testing.assert_array_equal(dosages["rs1"], [1.0, 1.0, 2.0])


if __name__ == "__main__":
    unittest.main()
