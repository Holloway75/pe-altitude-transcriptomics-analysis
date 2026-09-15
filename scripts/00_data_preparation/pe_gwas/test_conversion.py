#!/usr/bin/env python3
"""Unit tests for allele matching; no external genomic files are required."""
import importlib.util
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "convert", HERE / "convert_gwas_hg38_to_hg19.py")
convert = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(convert)


class AlleleMatchingTest(unittest.TestCase):
    def test_same_strand_direct_and_swapped(self):
        self.assertEqual(convert.try_match("C", "A", "A", "C"),
                         ("snp_same", "C", "A"))
        self.assertEqual(convert.try_match("A", "C", "A", "C"),
                         ("snp_same", "A", "C"))

    def test_negative_strand_preserves_physical_effect(self):
        self.assertEqual(convert.try_match("T", "G", "A", "C"),
                         ("snp_neg", "A", "C"))
        self.assertEqual(convert.try_match("G", "T", "A", "C"),
                         ("snp_neg", "C", "A"))

    def test_mismatch(self):
        self.assertIsNone(convert.try_match("A", "C", "A", "G"))

    def test_indel_representation_is_trimmed(self):
        # 1KG REF/ALT=AT/A and GWAS other/effect=T/empty after trimming.
        self.assertEqual(convert.try_match("A", "AT", "AT", "A"),
                         ("indel", "A", "AT"))
        self.assertEqual(convert.try_match("AT", "A", "AT", "A"),
                         ("indel", "AT", "A"))

    def test_symbolic_alleles_invalid(self):
        self.assertFalse(convert.valid_allele("<DEL>"))
        self.assertFalse(convert.valid_allele("."))
        self.assertTrue(convert.valid_allele("AC"))


if __name__ == "__main__":
    unittest.main()
