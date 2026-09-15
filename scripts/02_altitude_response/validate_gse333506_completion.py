#!/usr/bin/env python3

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out = np.empty(n)
    out[order] = adjusted
    return out


def verify_manifest(directory: Path):
    manifest = pd.read_csv(directory / "artifact_manifest.tsv", sep="\t")
    failures = []
    for row in manifest.itertuples(index=False):
        path = directory / row.file
        if not path.exists():
            failures.append(f"missing:{row.file}")
        elif path.stat().st_size != row.bytes:
            failures.append(f"size:{row.file}")
        elif sha256(path) != row.sha256:
            failures.append(f"sha256:{row.file}")
    return failures


def verify_source_manifest(path: Path):
    manifest = pd.read_csv(path, sep="\t")
    failures = []
    for row in manifest.itertuples(index=False):
        source = Path(row.path)
        if not source.exists():
            failures.append(f"missing:{source}")
        elif source.stat().st_size != row.bytes:
            failures.append(f"size:{source}")
        elif sha256(source) != row.sha256:
            failures.append(f"sha256:{source}")
    return failures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--frozen-dir", required=True, type=Path)
    parser.add_argument("--completion-dir", required=True, type=Path)
    parser.add_argument("--active-results-root", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    checks = {}
    errors = []

    sample = pd.read_csv(args.result_dir / "sample_metadata.tsv", sep="\t")
    primary = sample[sample["include_primary"]]
    checks["sample_contract"] = {
        "n_total": int(len(sample)),
        "n_primary": int(len(primary)),
        "n_primary_donors": int(primary["donor_id"].nunique()),
        "timepoint_counts": primary["timepoint"].value_counts().to_dict(),
        "population": sorted(primary["population"].unique().tolist()),
    }
    if len(sample) != 42 or len(primary) != 14:
        errors.append("sample_count_contract")
    if primary["donor_id"].nunique() != 7:
        errors.append("donor_count_contract")
    if set(primary["timepoint"]) != {"baseline", "week 4"}:
        errors.append("timepoint_contract")
    if set(primary["population"]) != {"Han Chinese"}:
        errors.append("population_contract")

    pairs = pd.read_csv(args.result_dir / "donor_altitude_mapping.tsv", sep="\t")
    if len(pairs) != 7 or pairs["donor_id"].nunique() != 7:
        errors.append("pair_mapping_contract")

    de = pd.read_csv(args.result_dir / "differential_expression.tsv", sep="\t")
    recalculated = bh_adjust(de["P.Value"].to_numpy())
    max_fdr_error = float(np.max(np.abs(recalculated - de["adj.P.Val"])))
    checks["differential_expression"] = {
        "n_genes": int(len(de)),
        "n_unique_symbols": int(de["gene_symbol"].nunique()),
        "fdr_0_05": int((de["adj.P.Val"] < 0.05).sum()),
        "up_fdr_0_05": int(
            ((de["adj.P.Val"] < 0.05) & (de["logFC"] > 0)).sum()
        ),
        "down_fdr_0_05": int(
            ((de["adj.P.Val"] < 0.05) & (de["logFC"] < 0)).sum()
        ),
        "max_bh_recalculation_error": max_fdr_error,
    }
    if len(de) != 12101 or de["gene_symbol"].nunique() != len(de):
        errors.append("gene_universe_contract")
    if max_fdr_error > 1e-12:
        errors.append("bh_recalculation")
    if not np.array_equal(
        de["direction"].to_numpy(),
        np.where(de["logFC"] > 0, "up", "down"),
    ):
        errors.append("direction_encoding")

    summary = pd.read_csv(args.result_dir / "contrast_summary.tsv", sep="\t").iloc[0]
    for name, value in {
        "n_genes_tested": len(de),
        "fdr_0.05": (de["adj.P.Val"] < 0.05).sum(),
        "up_fdr_0.05": ((de["adj.P.Val"] < 0.05) & (de["logFC"] > 0)).sum(),
        "down_fdr_0.05": ((de["adj.P.Val"] < 0.05) & (de["logFC"] < 0)).sum(),
    }.items():
        if int(summary[name]) != int(value):
            errors.append(f"contrast_summary:{name}")

    qc = pd.read_csv(args.result_dir / "sample_qc.tsv", sep="\t")
    checks["sample_qc"] = {
        "n_samples": int(len(qc)),
        "all_good_grid": bool((qc["is_good_grid"] == 1).all()),
        "pca_outliers": int(qc["pca_outlier_first5"].sum()),
        "minimum_within_pair_correlation":
            float(qc["within_pair_correlation"].min()),
        "minimum_median_sample_correlation":
            float(qc["median_correlation_to_other_samples"].min()),
    }
    if len(qc) != 14 or not (qc["is_good_grid"] == 1).all():
        errors.append("raw_array_qc")

    loo = pd.read_csv(args.result_dir / "leave_one_donor_out.tsv", sep="\t")
    sig_loo = loo[loo["primary_fdr"] < 0.05]
    checks["leave_one_donor_out"] = {
        "n_primary_fdr": int(len(sig_loo)),
        "all_sign_concordant": bool((sig_loo["loo_sign_concordance"] == 1).all()),
        "all_max_p_lt_0_05": bool((sig_loo["loo_max_p"] < 0.05).all()),
    }

    for filename in (
        "sample_metadata.tsv",
        "donor_altitude_mapping.tsv",
        "probe_mapping_audit.tsv",
        "source_manifest.tsv",
    ):
        if sha256(args.result_dir / filename) != sha256(args.frozen_dir / filename):
            errors.append(f"frozen_result_mismatch:{filename}")

    result_manifest_failures = verify_manifest(args.result_dir)
    completion_manifest_failures = verify_manifest(args.completion_dir)
    source_failures = verify_source_manifest(args.result_dir / "source_manifest.tsv")
    checks["manifest"] = {
        "result_failures": result_manifest_failures,
        "completion_failures": completion_manifest_failures,
        "source_failures": source_failures,
    }
    errors.extend(f"result_manifest:{x}" for x in result_manifest_failures)
    errors.extend(f"completion_manifest:{x}" for x in completion_manifest_failures)
    errors.extend(f"source_manifest:{x}" for x in source_failures)

    comparison = pd.read_csv(
        args.completion_dir / "adapted_layer_gene_comparison.tsv", sep="\t"
    )
    agreement = pd.read_csv(
        args.completion_dir / "adapted_layer_agreement.tsv", sep="\t"
    )
    recomputed = []
    for estimand, frame in comparison.groupby("g196_estimand"):
        recomputed.append(
            {
                "g196_estimand": estimand,
                "n_common_genes": len(frame),
                "pearson_log2fc":
                    frame["g196_log2fc"].corr(frame["g333_log2fc"],
                                              method="pearson"),
                "spearman_log2fc":
                    frame["g196_log2fc"].rank().corr(
                        frame["g333_log2fc"].rank(), method="pearson"
                    ),
                "sign_concordance": frame["same_direction"].mean(),
            }
        )
    recomputed = pd.DataFrame(recomputed)
    merged = agreement.merge(recomputed, on="g196_estimand",
                             suffixes=("_reported", "_recomputed"))
    agreement_error = max(
        float(np.max(np.abs(
            merged[f"{metric}_reported"] - merged[f"{metric}_recomputed"]
        )))
        for metric in ("pearson_log2fc", "spearman_log2fc",
                       "sign_concordance")
    )
    checks["adapted_layer"] = {
        "agreement_max_recalculation_error": agreement_error,
        "reported": agreement.to_dict(orient="records"),
    }
    # R and pandas use slightly different floating-point tie handling for
    # Spearman ranks; 1e-6 is far below the reported precision.
    if agreement_error > 1e-6:
        errors.append("adapted_layer_agreement_recalculation")

    retired_named_paths = [
        str(path) for path in args.active_results_root.rglob("*GSE75665*")
    ]
    checks["retired_dataset_paths"] = retired_named_paths
    if retired_named_paths:
        errors.append("retired_GSE75665_named_path_in_active_results")

    report = {
        "analysis_id": "m2_gse333506_han_week4_paired_v1_20260727",
        "completion_id": "m2_adapted_layer_completion_v1_20260727",
        "status": "PASS" if not errors else "FAIL",
        "checks": checks,
        "errors": errors,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": report["status"], "errors": errors},
                     ensure_ascii=False))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
