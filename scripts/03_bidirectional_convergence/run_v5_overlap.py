#!/usr/bin/env python3
"""Run the frozen Module 3 v5 directional-overlap analysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats


QUADRANTS = (
    ("alt_up_pe_positive", "up", "positive"),
    ("alt_up_pe_negative", "up", "negative"),
    ("alt_down_pe_positive", "down", "positive"),
    ("alt_down_pe_negative", "down", "negative"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = columns or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def holm(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    adjusted = [1.0] * n
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (n - rank) * values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def fisher_result(a: int, alt_size: int, pe_size: int, universe_size: int) -> tuple[float, float, float, float]:
    b = alt_size - a
    c = pe_size - a
    d = universe_size - a - b - c
    table = np.asarray([[a, b], [c, d]], dtype=np.int64)
    p = float(stats.fisher_exact(table, alternative="greater").pvalue)
    estimate = stats.contingency.odds_ratio(table, kind="conditional")
    interval = estimate.confidence_interval(confidence_level=0.95, alternative="two-sided")
    return float(estimate.statistic), float(interval.low), float(interval.high), p


def fixed_tail_indices(values: np.ndarray, n_up: int, n_down: int) -> tuple[np.ndarray, np.ndarray]:
    if n_up <= 0 or n_down <= 0 or n_up + n_down > values.size:
        raise ValueError("invalid fixed tail sizes")
    up = np.argpartition(values, values.size - n_up)[-n_up:]
    down = np.argpartition(values, n_down - 1)[:n_down]
    return up, down


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--v4-registry-dir", type=Path, required=True)
    parser.add_argument("--module2-release", type=Path, required=True)
    parser.add_argument("--limma-prior", type=Path, required=True)
    parser.add_argument("--locked-dir", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()

    project = args.project_root.resolve()
    config_path = args.config.resolve()
    v4_registry = args.v4_registry_dir.resolve()
    module2_release = args.module2_release.resolve()
    limma_prior_path = args.limma_prior.resolve()
    locked_dir = args.locked_dir.resolve()
    result_dir = args.result_dir.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))

    for directory in (locked_dir, result_dir):
        if directory.exists() and any(directory.iterdir()):
            raise SystemExit(f"refusing to overwrite non-empty release: {directory}")
        directory.mkdir(parents=True, exist_ok=True)

    universe_path = v4_registry / "formal_gene_universe.tsv"
    mapping_audit_path = v4_registry / "gene_mapping_audit.tsv"
    formal_stats_path = module2_release / "gene_statistics/ALT_GSE103927_ALT16_TOTAL.tsv.gz"
    donor_path = module2_release / "donor_contrasts/ALT_GSE103927_ALT16_TOTAL.tsv.gz"
    context_stats_path = module2_release / "gene_statistics/ALT_GSE333506_HAN_WEEK4_TOTAL.tsv.gz"
    pe_path = project / "results/01_pe_genetic_map/CURRENT/pe_whole_blood.tsv"
    required = [config_path, universe_path, mapping_audit_path, formal_stats_path, donor_path, context_stats_path, pe_path, limma_prior_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("missing input(s): " + ", ".join(missing))

    universe = read_tsv(universe_path)
    if len(universe) != 8101 or len({row["hgnc_id"] for row in universe}) != len(universe):
        raise ValueError("v5 requires the reconstructed unique 8,101-gene v4 universe")
    formal_stats = {row["hgnc_id"]: row for row in read_tsv(formal_stats_path)}
    limma_prior = {row["hgnc_id"]: row for row in read_tsv(limma_prior_path)}
    pe_stats = {row["gene_id"].split(".")[0]: row for row in read_tsv(pe_path)}
    context_by_symbol = {row["hgnc_symbol"]: row for row in read_tsv(context_stats_path)}

    donor_rows = read_tsv(donor_path)
    donor_by_id = {row["hgnc_id"]: row for row in donor_rows}
    donor_columns = [name for name in donor_rows[0] if name not in {"hgnc_id", "hgnc_symbol"}]
    if len(donor_columns) != 21:
        raise ValueError(f"expected 21 donors, found {len(donor_columns)}")

    genes: list[dict[str, object]] = []
    matrix = np.empty((len(universe), len(donor_columns)), dtype=np.float64)
    df_residual = np.empty(len(universe), dtype=np.float64)
    df_prior = np.empty(len(universe), dtype=np.float64)
    s2_prior = np.empty(len(universe), dtype=np.float64)
    stdev_unscaled = np.empty(len(universe), dtype=np.float64)
    for index, row in enumerate(universe):
        hgnc_id = row["hgnc_id"]
        alt = formal_stats[hgnc_id]
        pe = pe_stats[row["jti_ensembl_gene_id"].split(".")[0]]
        donor = donor_by_id[hgnc_id]
        prior = limma_prior[hgnc_id]
        matrix[index, :] = [float(donor[column]) for column in donor_columns]
        df_residual[index] = float(prior["df_residual"])
        df_prior[index] = float(prior["df_prior"])
        s2_prior[index] = float(prior["s2_prior"])
        stdev_unscaled[index] = float(prior["stdev_unscaled"])
        context = context_by_symbol.get(row["hgnc_symbol"])
        genes.append(
            {
                **row,
                "alt_module2_log2fc": float(alt["module2_moderated_logFC"]),
                "alt_module2_t": float(alt["module2_moderated_t"]),
                "alt_module2_p": float(alt["module2_p"]),
                "alt_module2_fdr": float(alt["module2_fdr"]),
                "pe_q_gene_tissue": float(pe["q_gene_tissue"]),
                "context_gse333506_log2fc": float(context["module2_moderated_logFC"]) if context else math.nan,
                "context_gse333506_t": float(context["module2_moderated_t"]) if context else math.nan,
                "context_gse333506_fdr": float(context["module2_fdr"]) if context else math.nan,
            }
        )

    if not np.isfinite(matrix).all():
        raise ValueError("donor matrix has non-finite values")
    if not all(np.isfinite(x).all() for x in (df_residual, df_prior, s2_prior, stdev_unscaled)):
        raise ValueError("limma prior table has non-finite values")
    observed_means = matrix.mean(axis=1)
    module2_means = np.asarray([float(row["alt_module2_log2fc"]) for row in genes])
    max_mean_delta = float(np.max(np.abs(observed_means - module2_means)))
    if max_mean_delta > 1e-8:
        raise ValueError(f"donor means do not reproduce module2 log2FC: {max_mean_delta}")
    observed_variance = matrix.var(axis=1, ddof=1)
    observed_s2_post = (df_residual * observed_variance + df_prior * s2_prior) / (df_residual + df_prior)
    reconstructed_moderated_t = observed_means / (np.sqrt(observed_s2_post) * stdev_unscaled)
    module2_moderated_t = np.asarray([float(row["alt_module2_t"]) for row in genes])
    max_moderated_t_delta = float(np.max(np.abs(reconstructed_moderated_t - module2_moderated_t)))
    if max_moderated_t_delta > 1e-8:
        raise ValueError(f"frozen limma priors do not reproduce moderated t: {max_moderated_t_delta}")

    alt_fc = module2_means
    alt_fdr = np.asarray([float(row["alt_module2_fdr"]) for row in genes])
    alt_t = np.asarray([float(row["alt_module2_t"]) for row in genes])
    pe_z = np.asarray([float(row["pe_z"]) for row in genes])
    pe_q = np.asarray([float(row["pe_q_gene_tissue"]) for row in genes])
    strict_alt_up = (alt_fdr < config["strict_alt_fdr_max"]) & (alt_fc >= config["strict_alt_abs_log2fc_min"])
    strict_alt_down = (alt_fdr < config["strict_alt_fdr_max"]) & (alt_fc <= -config["strict_alt_abs_log2fc_min"])
    fdr_alt_up = (alt_fdr < config["strict_alt_fdr_max"]) & (alt_fc > 0)
    fdr_alt_down = (alt_fdr < config["strict_alt_fdr_max"]) & (alt_fc < 0)
    strict_pe_positive = (pe_q < config["strict_pe_q_max"]) & (pe_z > 0)
    strict_pe_negative = (pe_q < config["strict_pe_q_max"]) & (pe_z < 0)
    tail_n = int(math.floor(config["tail_fraction"] * len(genes)))
    alt_top_idx, alt_bottom_idx = fixed_tail_indices(alt_t, tail_n, tail_n)
    pe_top_idx, pe_bottom_idx = fixed_tail_indices(pe_z, tail_n, tail_n)
    tail_alt_up = np.zeros(len(genes), dtype=bool); tail_alt_up[alt_top_idx] = True
    tail_alt_down = np.zeros(len(genes), dtype=bool); tail_alt_down[alt_bottom_idx] = True
    tail_pe_positive = np.zeros(len(genes), dtype=bool); tail_pe_positive[pe_top_idx] = True
    tail_pe_negative = np.zeros(len(genes), dtype=bool); tail_pe_negative[pe_bottom_idx] = True

    scenarios = {
        "PRIMARY_STRICT": (strict_alt_up, strict_alt_down, strict_pe_positive, strict_pe_negative),
        "SENSITIVITY_ALT_FDR_ONLY": (fdr_alt_up, fdr_alt_down, strict_pe_positive, strict_pe_negative),
        "SENSITIVITY_BOTH_AXES_5PCT_TAILS": (tail_alt_up, tail_alt_down, tail_pe_positive, tail_pe_negative),
    }

    rng = np.random.default_rng(int(config["sign_flip_seed"]))
    draws = int(config["sign_flip_draws"])
    signs = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=(draws, len(donor_columns)), replace=True)
    sign_path = locked_dir / "donor_sign_flip_matrix.npz"
    np.savez_compressed(sign_path, signs=signs, donor_ids=np.asarray(donor_columns), seed=np.asarray([config["sign_flip_seed"]]))

    squared_sums = np.square(matrix).sum(axis=1)[:, None]
    n = matrix.shape[1]
    null_counts = {scenario: np.zeros((draws, 4), dtype=np.int32) for scenario in scenarios}
    set_sizes = {
        scenario: (int(up.sum()), int(down.sum()), int(pos.sum()), int(neg.sum()))
        for scenario, (up, down, pos, neg) in scenarios.items()
    }
    pe_masks = {
        scenario: (pos.astype(np.int8), neg.astype(np.int8))
        for scenario, (_, _, pos, neg) in scenarios.items()
    }
    chunk_size = int(config["monte_carlo_chunk_size"])
    for start in range(0, draws, chunk_size):
        stop = min(draws, start + chunk_size)
        signed = signs[start:stop, :].T.astype(np.float64, copy=False)
        means = matrix @ signed / n
        variances = np.maximum((squared_sums - n * np.square(means)) / (n - 1), 0.0)
        moderated_variances = (
            df_residual[:, None] * variances + df_prior[:, None] * s2_prior[:, None]
        ) / (df_residual[:, None] + df_prior[:, None])
        standard_errors = np.sqrt(moderated_variances) * stdev_unscaled[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            t_values = np.divide(means, standard_errors, out=np.zeros_like(means), where=standard_errors > 0)
        for scenario, sizes in set_sizes.items():
            n_up, n_down, _, _ = sizes
            pe_pos, pe_neg = pe_masks[scenario]
            up_idx = np.argpartition(t_values, t_values.shape[0] - n_up, axis=0)[-n_up:, :]
            down_idx = np.argpartition(t_values, n_down - 1, axis=0)[:n_down, :]
            counts = null_counts[scenario]
            counts[start:stop, 0] = pe_pos[up_idx].sum(axis=0)
            counts[start:stop, 1] = pe_neg[up_idx].sum(axis=0)
            counts[start:stop, 2] = pe_pos[down_idx].sum(axis=0)
            counts[start:stop, 3] = pe_neg[down_idx].sum(axis=0)

    result_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    boundary_rows: list[dict[str, object]] = []
    gene_memberships: dict[str, list[str]] = {str(row["hgnc_id"]): [] for row in genes}
    for scenario, masks in scenarios.items():
        alt_up, alt_down, pe_pos, pe_neg = masks
        mask_lookup = {"up": alt_up, "down": alt_down, "positive": pe_pos, "negative": pe_neg}
        raw_empirical: list[float] = []
        raw_fisher: list[float] = []
        scenario_rows: list[dict[str, object]] = []
        for quadrant_index, (quadrant, altitude_direction, pe_direction) in enumerate(QUADRANTS):
            alt_mask = mask_lookup[altitude_direction]
            pe_mask = mask_lookup[pe_direction]
            overlap_mask = alt_mask & pe_mask
            observed = int(overlap_mask.sum())
            null = null_counts[scenario][:, quadrant_index]
            exceedances = int(np.count_nonzero(null >= observed))
            empirical_p = (exceedances + 1) / (draws + 1)
            odds_ratio, ci_low, ci_high, fisher_p = fisher_result(observed, int(alt_mask.sum()), int(pe_mask.sum()), len(genes))
            null_mean = float(null.mean())
            row = {
                "analysis_id": config["analysis_id"],
                "scenario": scenario,
                "inference_role": "primary" if scenario == "PRIMARY_STRICT" else "sensitivity",
                "quadrant": quadrant,
                "altitude_direction": altitude_direction,
                "pe_direction": pe_direction,
                "universe_n": len(genes),
                "altitude_set_n": int(alt_mask.sum()),
                "pe_set_n": int(pe_mask.sum()),
                "observed_overlap_n": observed,
                "fisher_odds_ratio": odds_ratio,
                "fisher_ci95_low": ci_low,
                "fisher_ci95_high": ci_high,
                "fisher_p_one_sided": fisher_p,
                "empirical_draws": draws,
                "empirical_exceedances": exceedances,
                "empirical_p_one_sided": empirical_p,
                "empirical_mcse": math.sqrt(empirical_p * (1 - empirical_p) / (draws + 1)),
                "null_mean_overlap": null_mean,
                "null_q95": float(np.quantile(null, 0.95, method="higher")),
                "null_q99": float(np.quantile(null, 0.99, method="higher")),
                "empirical_enrichment": observed / null_mean if null_mean > 0 else math.inf,
            }
            scenario_rows.append(row)
            raw_empirical.append(empirical_p)
            raw_fisher.append(fisher_p)
            membership = np.flatnonzero(overlap_mask)
            for gene_index in membership:
                gene = genes[gene_index]
                gene_memberships[str(gene["hgnc_id"])].append(f"{scenario}:{quadrant}")
                gene_rows.append(
                    {
                        "analysis_id": config["analysis_id"],
                        "scenario": scenario,
                        "quadrant": quadrant,
                        "hgnc_id": gene["hgnc_id"],
                        "hgnc_symbol": gene["hgnc_symbol"],
                        "alt_log2fc": gene["alt_module2_log2fc"],
                        "alt_t": gene["alt_module2_t"],
                        "alt_fdr": gene["alt_module2_fdr"],
                        "pe_z": gene["pe_z"],
                        "pe_q": gene["pe_q_gene_tissue"],
                        "context_gse333506_log2fc": gene["context_gse333506_log2fc"],
                        "context_gse333506_t": gene["context_gse333506_t"],
                        "context_gse333506_fdr": gene["context_gse333506_fdr"],
                    }
                )
        empirical_holm = holm(raw_empirical)
        fisher_holm = holm(raw_fisher)
        for index, row in enumerate(scenario_rows):
            row["empirical_holm_p"] = empirical_holm[index]
            row["fisher_holm_p"] = fisher_holm[index]
            null = null_counts[scenario][:, index]
            threshold = None
            for overlap in range(0, min(int(row["altitude_set_n"]), int(row["pe_set_n"])) + 1):
                raw_p = (1 + int(np.count_nonzero(null >= overlap))) / (draws + 1)
                if min(1.0, 4 * raw_p) < 0.05:
                    threshold = overlap
                    break
            boundary_rows.append(
                {
                    "analysis_id": config["analysis_id"],
                    "scenario": scenario,
                    "quadrant": row["quadrant"],
                    "universe_n": len(genes),
                    "altitude_set_n": row["altitude_set_n"],
                    "pe_set_n": row["pe_set_n"],
                    "null_mean_overlap": row["null_mean_overlap"],
                    "observed_overlap_n": row["observed_overlap_n"],
                    "minimum_overlap_for_bonferroni_upper_bound_lt_0_05": threshold if threshold is not None else "NA",
                    "minimum_enrichment_over_null_mean": threshold / float(row["null_mean_overlap"]) if threshold is not None and float(row["null_mean_overlap"]) > 0 else "NA",
                    "boundary_note": "Conservative sufficient boundary using 4*p_empirical; exact Holm significance depends on all four observed P values",
                }
            )
        result_rows.extend(scenario_rows)

    candidate_trigger = any(float(row["empirical_holm_p"]) < 0.05 for row in result_rows if row["scenario"] == "PRIMARY_STRICT")
    for index, row in enumerate(genes):
        row["strict_alt_up"] = bool(strict_alt_up[index])
        row["strict_alt_down"] = bool(strict_alt_down[index])
        row["strict_pe_positive"] = bool(strict_pe_positive[index])
        row["strict_pe_negative"] = bool(strict_pe_negative[index])
        row["overlap_memberships"] = "|".join(gene_memberships[str(row["hgnc_id"])]) or "NA"

    universe_output = locked_dir / "gene_universe_and_mapping.tsv"
    write_tsv(universe_output, genes)
    write_tsv(result_dir / "directional_overlap.tsv", result_rows)
    if gene_rows:
        write_tsv(result_dir / "directional_overlap_genes.tsv", gene_rows)
    else:
        write_tsv(result_dir / "directional_overlap_genes.tsv", [{"analysis_id": config["analysis_id"], "scenario": "NA", "quadrant": "NA", "hgnc_id": "NA", "hgnc_symbol": "NA", "note": "no observed overlaps"}])
    write_tsv(result_dir / "detectability_boundary.tsv", boundary_rows)

    input_sources = {
        "config": config_path,
        "v4_reconstructed_universe": universe_path,
        "v4_mapping_audit": mapping_audit_path,
        "formal_alt_statistics": formal_stats_path,
        "formal_alt_donor_contrasts": donor_path,
        "context_alt_statistics": context_stats_path,
        "pe_whole_blood": pe_path,
        "gse103927_limma_empirical_bayes_prior": limma_prior_path,
    }
    registry = {
        **config,
        "execution_status": "OVERLAP_COMPLETE",
        "universe_n": len(genes),
        "donor_n": len(donor_columns),
        "donor_ids": donor_columns,
        "historical_v4_universe_n": 8101,
        "historical_v4_universe_preserved": True,
        "hgnc_recovery_note": "The missing 2026-07-27 file was replaced by the frozen official current HGNC file; the archived v4 freeze code reproduced the identical 8,101-gene universe.",
        "input_sources": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in input_sources.items()},
        "locked_outputs": {
            "gene_universe_and_mapping": {"path": str(universe_output), "sha256": sha256(universe_output)},
            "donor_sign_flip_matrix": {"path": str(sign_path), "sha256": sha256(sign_path)},
        },
        "result_outputs": {
            "directional_overlap": str(result_dir / "directional_overlap.tsv"),
            "directional_overlap_genes": str(result_dir / "directional_overlap_genes.tsv"),
            "detectability_boundary": str(result_dir / "detectability_boundary.tsv"),
        },
        "candidate_trigger_passed": candidate_trigger,
        "candidate_genes_file_created": False,
        "validation_metrics": {
            "max_donor_mean_vs_module2_log2fc_delta": max_mean_delta,
            "max_reconstructed_vs_module2_moderated_t_delta": max_moderated_t_delta
        },
    }
    write_json(locked_dir / "module3_v5_analysis_registry.json", registry)
    print(json.dumps({"status": "PASS", "universe_n": len(genes), "draws": draws, "candidate_trigger_passed": candidate_trigger}))


if __name__ == "__main__":
    main()
