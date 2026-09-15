#!/usr/bin/env python3
"""Shared helpers for the guarded genome-wide S-MultiXcan skill scripts."""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
from pathlib import Path


OUTPUT_COLUMNS = [
    "gene", "gene_name", "pvalue", "n", "n_indep", "p_i_best",
    "t_i_best", "p_i_worst", "t_i_worst", "eigen_max", "eigen_min",
    "eigen_min_kept", "z_min", "z_max", "z_mean", "z_sd", "tmi",
    "status",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else path.open(mode, newline="")


def write_tsv(path: Path, fieldnames, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bh(pvalues):
    """Benjamini-Hochberg values in the original input order."""
    if not pvalues:
        return []
    indexed = sorted(enumerate(map(float, pvalues)), key=lambda x: x[1])
    result = [math.nan] * len(indexed)
    running = 1.0
    for rank in range(len(indexed), 0, -1):
        index, value = indexed[rank - 1]
        running = min(running, value * len(indexed) / rank, 1.0)
        result[index] = running
    return result


def classify_smultixcan(row):
    status = str(row.get("status", ""))
    if status not in {"0", "0.0"}:
        return {
            "-1": "no_data",
            "-2": "no_spredixcan_results",
            "-3": "no_model_product",
            "-4": "pvalue_numerical_underflow",
            "-5": "singular_covariance",
            "-6": "inverse_error",
            "-7": "complex_covariance",
            "-8": "inadequate_inverse",
        }.get(status, "model_failed")
    try:
        n = int(float(row["n"]))
        n_indep = int(float(row["n_indep"]))
        eigen_max = float(row["eigen_max"])
        eigen_min_kept = float(row["eigen_min_kept"])
        pvalue = float(row["pvalue"])
    except (KeyError, TypeError, ValueError):
        return "invalid_matrix_diagnostics"
    if n < 2 or n_indep < 2:
        return "degenerate_not_cross_tissue"
    if not (
        math.isfinite(pvalue)
        and 0.0 < pvalue <= 1.0
        and math.isfinite(eigen_max)
        and math.isfinite(eigen_min_kept)
        and eigen_max > 0.0
        and eigen_min_kept > 0.0
        and math.isfinite(eigen_max / eigen_min_kept)
    ):
        return "invalid_matrix_diagnostics"
    return "success"


def normalized_gene(value: str) -> str:
    return str(value).split(".", 1)[0]

