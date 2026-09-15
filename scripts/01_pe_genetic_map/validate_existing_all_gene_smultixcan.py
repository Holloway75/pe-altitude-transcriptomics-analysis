#!/usr/bin/env python3
"""Revalidate existing genome-wide S-MultiXcan outputs under outcome-subset BH rules."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from smultixcan_genome_core import bh, classify_smultixcan, normalized_gene, sha256  # noqa: E402
from run_all_gene_smultixcan import OUTPUT_COLUMNS  # noqa: E402


def covariance_genes(path: Path) -> set[str]:
    with gzip.open(path, "rt", newline="") as handle:
        genes = {
            normalized_gene(row["gene"])
            for row in csv.DictReader(handle, delimiter="\t")
            if row["build_status"] == "success"
        }
    if not genes:
        raise RuntimeError("no valid covariance genes")
    return genes


def rows_by_gene(path: Path, expected_prefix: list[str] | None = None):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if expected_prefix and (reader.fieldnames or [])[: len(expected_prefix)] != expected_prefix:
            raise RuntimeError(f"unexpected schema: {path}")
        rows = list(reader)
    genes = [normalized_gene(row["gene"]) for row in rows]
    if len(genes) != len(set(genes)):
        raise RuntimeError(f"duplicate gene rows: {path}")
    return dict(zip(genes, rows))


def close(a: float, b: float, tolerance: float = 5e-15) -> bool:
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def validate_dataset(source: Path, valid_genes: set[str], output: Path):
    raw_path = source / "smultixcan.raw.tsv"
    enriched_path = source / "smultixcan_all_genes.tsv"
    method_path = source / "method.json"
    for path in (raw_path, enriched_path, method_path):
        if not path.is_file():
            raise RuntimeError(f"required artifact absent: {path}")

    raw = rows_by_gene(raw_path, OUTPUT_COLUMNS)
    enriched = rows_by_gene(enriched_path, ["dataset_id"] + OUTPUT_COLUMNS)
    returned = set(raw)
    missing = sorted(valid_genes - returned)
    unexpected = sorted(returned - valid_genes)
    successful = sorted(gene for gene, row in raw.items() if classify_smultixcan(row) == "success")
    expected_q = dict(zip(successful, bh([float(raw[gene]["pvalue"]) for gene in successful])))

    classification_errors = []
    q_errors = []
    raw_enriched_errors = []
    assigned = 0
    for gene, raw_row in raw.items():
        out = enriched.get(gene)
        if out is None:
            raw_enriched_errors.append(gene)
            continue
        for field in OUTPUT_COLUMNS:
            if out[field] != raw_row[field]:
                raw_enriched_errors.append(f"{gene}:{field}")
                break
        observed_class = out.get("classification", "")
        expected_class = classify_smultixcan(raw_row)
        if observed_class != expected_class:
            classification_errors.append(gene)
        observed_q = out.get("qvalue_global_bh", "NA")
        if gene in expected_q:
            assigned += observed_q != "NA"
            if observed_q == "NA" or not close(float(observed_q), expected_q[gene]):
                q_errors.append(gene)
        elif observed_q != "NA":
            q_errors.append(gene)

    output.mkdir(parents=True, exist_ok=False)
    audit_path = output / "gene-execution-audit.tsv"
    with audit_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "raw_row_present", "classification", "bh_eligible"])
        for gene in sorted(valid_genes | returned):
            if gene not in raw:
                writer.writerow([gene, "false", "outcome_unavailable", "false"])
            else:
                classification = classify_smultixcan(raw[gene])
                writer.writerow([gene, "true", classification, str(classification == "success").lower()])

    checks = {
        "covariance_resource_has_valid_genes": bool(valid_genes),
        "raw_exact_schema": True,
        "raw_gene_rows_unique": True,
        "enriched_gene_rows_unique": True,
        "raw_enriched_gene_sets_equal": set(raw) == set(enriched),
        "raw_enriched_statistics_equal": not raw_enriched_errors,
        "returned_genes_subset_of_valid_covariance_genes": not unexpected,
        "all_valid_covariance_genes_audited": sum(1 for _ in open(audit_path)) - 1 == len(valid_genes | returned),
        "classifications_recomputed_equal": not classification_errors,
        "successful_cross_tissue_family_nonempty": bool(successful),
        "bh_family_qvalue_closed": assigned == len(successful),
        "bh_qvalues_recomputed_equal": not q_errors,
    }
    report = {
        "schema_version": "m1.5-outcome-subset-validation-v2",
        "dataset_id": source.name,
        "status": "pass" if all(checks.values()) else "fail",
        "checks": checks,
        "counts": {
            "valid_covariance_genes": len(valid_genes),
            "returned_gene_rows": len(raw),
            "outcome_unavailable_genes": len(missing),
            "unexpected_returned_genes": len(unexpected),
            "global_bh_family_size": len(successful),
            "global_bh_qvalues_assigned": assigned,
            "global_bh_significant_q_lt_0.05": sum(expected_q[x] < 0.05 for x in successful),
        },
        "examples": {
            "outcome_unavailable": missing[:20],
            "unexpected_returned": unexpected[:20],
            "classification_errors": classification_errors[:20],
            "qvalue_errors": q_errors[:20],
            "raw_enriched_errors": raw_enriched_errors[:20],
        },
        "source_artifacts": {
            "raw": {"path": str(raw_path.resolve()), "sha256": sha256(raw_path)},
            "enriched": {"path": str(enriched_path.resolve()), "sha256": sha256(enriched_path)},
            "method": {"path": str(method_path.resolve()), "sha256": sha256(method_path)},
        },
        "audit": {"path": str(audit_path.resolve()), "sha256": sha256(audit_path)},
    }
    report_path = output / "validation-report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report, report_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--covariance-resource", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    covariance_resource = args.covariance_resource.resolve()
    output = args.output.resolve()
    if output.exists():
        raise RuntimeError(f"output must not exist: {output}")
    resource_validation = covariance_resource / "validation-report.json"
    gene_qc = covariance_resource / "gene_qc.tsv.gz"
    if json.loads(resource_validation.read_text()).get("status") != "pass":
        raise RuntimeError("covariance resource validation is not pass")
    valid_genes = covariance_genes(gene_qc)
    output.mkdir(parents=True)
    reports = []
    report_paths = []
    for dataset in sorted(path for path in source_root.iterdir() if path.is_dir()):
        report, path = validate_dataset(dataset, valid_genes, output / dataset.name)
        reports.append(report)
        report_paths.append(path)
    overall = {
        "schema_version": "m1.5-outcome-subset-validation-v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass" if reports and all(x["status"] == "pass" for x in reports) else "fail",
        "rule": "BH separately per GWAS across returned rows classified success; outcome-unavailable valid-covariance genes excluded",
        "covariance_resource": {
            "path": str(covariance_resource),
            "validation_sha256": sha256(resource_validation),
            "gene_qc_sha256": sha256(gene_qc),
            "valid_genes": len(valid_genes),
        },
        "source_root": str(source_root),
        "datasets": [
            {"dataset_id": x["dataset_id"], "status": x["status"], **x["counts"]}
            for x in reports
        ],
        "dataset_validation_sha256": {path.parent.name: sha256(path) for path in report_paths},
        "validator": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
    }
    (output / "validation-report.json").write_text(json.dumps(overall, indent=2, sort_keys=True) + "\n")
    print(json.dumps(overall, sort_keys=True))
    if overall["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
