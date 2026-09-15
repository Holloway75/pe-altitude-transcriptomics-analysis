#!/usr/bin/env python3
"""Threaded, result-equivalent runner for the archived Module 3 ALT kernel.

The archived implementation is loaded without modification.  Only ``donor_z``
is replaced by a row-partitioned implementation: every gene is evaluated by
the same formula and the input random-direction array is unchanged.  This is
useful for the expensive 999,999-draw reproduction run.  Equivalence must be
checked against the archived runner with identical seed, draw count, pathway
fraction, and chunk size before formal use.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy import special


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ARCHIVED_RUNNER = (
    PROJECT_ROOT
    / "archive/03_bidirectional_convergence/scripts/run_alt_continuous.py"
)


def load_archived_runner():
    spec = importlib.util.spec_from_file_location("module3_alt_archived", ARCHIVED_RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load archived runner: {ARCHIVED_RUNNER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def donor_z_block(contrast: np.ndarray, directions: np.ndarray) -> np.ndarray:
    n = contrast.shape[1]
    projections = contrast @ directions
    total_ss = np.sum(contrast * contrast, axis=1)[:, None]
    residual_ss = np.maximum(
        total_ss - projections * projections,
        np.finfo(np.float64).tiny,
    )
    t_value = projections / np.sqrt(residual_ss / (n - 1))
    probability = np.clip(
        special.stdtr(n - 1, t_value),
        np.finfo(np.float64).eps,
        1 - np.finfo(np.float64).eps,
    )
    return special.ndtri(probability)


def build_threaded_donor_z(workers: int):
    pool = ThreadPoolExecutor(max_workers=workers)

    def donor_z(contrast: np.ndarray, directions: np.ndarray) -> np.ndarray:
        if directions.shape[1] == 1 or workers == 1:
            return donor_z_block(contrast, directions)
        block_size = math.ceil(contrast.shape[0] / workers)
        blocks = [
            contrast[start : start + block_size]
            for start in range(0, contrast.shape[0], block_size)
        ]
        return np.vstack(list(pool.map(lambda block: donor_z_block(block, directions), blocks)))

    return donor_z


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--registry", type=Path, required=True)
    result.add_argument("--contrast", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    result.add_argument("--draws", type=int, required=True)
    result.add_argument("--seed", type=int, required=True)
    result.add_argument("--pathway-fraction", type=float, required=True)
    result.add_argument("--chunk-size", type=int, default=250)
    result.add_argument("--result-class", default="CALIBRATION_BENCHMARK_NOT_FORMAL")
    result.add_argument("--workers", type=int, default=8)
    return result


if __name__ == "__main__":
    args = parser().parse_args()
    archived = load_archived_runner()
    archived.donor_z = build_threaded_donor_z(args.workers)
    delattr(args, "workers")
    archived.main(args)
