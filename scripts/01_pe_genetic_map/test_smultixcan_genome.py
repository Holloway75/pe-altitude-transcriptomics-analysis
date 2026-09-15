#!/usr/bin/env python3
"""Fast unit tests for the genome-wide S-MultiXcan skill helpers."""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from smultixcan_genome_core import bh, classify_smultixcan  # noqa: E402
from run_all_gene_smultixcan import compare_gene_universes  # noqa: E402


def valid_row(**updates):
    row = {
        "status": "0", "n": "4", "n_indep": "3", "pvalue": "0.01",
        "eigen_max": "2", "eigen_min_kept": "0.1",
    }
    row.update(updates)
    return row


class GenomeCoreTests(unittest.TestCase):
    def test_bh_preserves_input_order_and_monotonicity(self):
        observed = bh([0.04, 0.001, 0.02])
        self.assertEqual(observed, [0.04, 0.003, 0.03])

    def test_success_gate(self):
        self.assertEqual(classify_smultixcan(valid_row()), "success")

    def test_degenerate_cross_tissue_is_not_success(self):
        self.assertEqual(
            classify_smultixcan(valid_row(n="1", n_indep="1")),
            "degenerate_not_cross_tissue",
        )

    def test_nonzero_status_mapping(self):
        self.assertEqual(
            classify_smultixcan(valid_row(status="-2")),
            "no_spredixcan_results",
        )

    def test_nonfinite_diagnostics_fail(self):
        self.assertEqual(
            classify_smultixcan(valid_row(eigen_min_kept=str(math.nan))),
            "invalid_matrix_diagnostics",
        )

    def test_outcome_specific_missing_genes_are_audited_not_unexpected(self):
        missing, unexpected = compare_gene_universes(
            {"ENSG1", "ENSG2", "ENSG3"}, {"ENSG1", "ENSG3"}
        )
        self.assertEqual(missing, ["ENSG2"])
        self.assertEqual(unexpected, [])

    def test_returned_gene_outside_covariance_set_is_unexpected(self):
        missing, unexpected = compare_gene_universes(
            {"ENSG1", "ENSG2"}, {"ENSG1", "ENSG9"}
        )
        self.assertEqual(missing, ["ENSG2"])
        self.assertEqual(unexpected, ["ENSG9"])


if __name__ == "__main__":
    unittest.main()
