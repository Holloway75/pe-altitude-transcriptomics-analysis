#!/usr/bin/env python3
"""Run guarded genome-wide S-MultiXcan sequentially for explicit GWAS inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from smultixcan_genome_core import (  # noqa: E402
    OUTPUT_COLUMNS,
    bh,
    classify_smultixcan,
    normalized_gene,
    sha256,
)


AUDIT_FIELDS = ["gene", "raw_row_present", "classification", "status"]
EXTRA_FIELDS = [
    "dataset_id", "condition_number", "classification", "qvalue_global_bh",
    "global_bh_significant", "direction_available", "signed_statistic", "direction",
]


def read_manifest(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ["dataset_id", "normalized_gwas", "spredixcan_folder"]
        if reader.fieldnames != required:
            raise RuntimeError(f"dataset manifest schema must be exactly {required}")
        rows = list(reader)
    if not rows or len({x["dataset_id"] for x in rows}) != len(rows):
        raise RuntimeError("dataset manifest must contain unique nonempty datasets")
    return rows


def read_covariance_genes(path: Path):
    genes = set()
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["build_status"] == "success":
                genes.add(normalized_gene(row["gene"]))
    if not genes:
        raise RuntimeError("covariance resource has no successful genes")
    return genes


def read_retained_rsids(path: Path):
    with gzip.open(path, "rt", newline="") as handle:
        return {row["rsid"] for row in csv.DictReader(handle, delimiter="\t")}


def write_cleared(normalized_gwas: Path, retained, output: Path):
    opener = gzip.open if str(normalized_gwas).endswith(".gz") else open
    found = set()
    with opener(normalized_gwas, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "rsid" not in (reader.fieldnames or []):
            raise RuntimeError(f"normalized GWAS lacks rsid: {normalized_gwas}")
        for row in reader:
            rsid = row.get("rsid", "")
            if rsid in retained:
                found.add(rsid)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["rsid"])
        writer.writerows((x,) for x in sorted(found))
    if not found:
        raise RuntimeError(f"no retained model SNPs intersect GWAS: {normalized_gwas}")
    return len(found)


def validate_spredixcan_folder(path: Path, expected_tissues):
    files = sorted(path.glob("*.csv"))
    if not files:
        raise RuntimeError(f"S-PrediXcan CSV folder empty: {path}")
    tissues = [x.stem for x in files]
    if len(tissues) != len(set(tissues)):
        raise RuntimeError(f"duplicate S-PrediXcan tissue names: {path}")
    if set(tissues) != set(expected_tissues):
        raise RuntimeError(
            f"S-PrediXcan/model tissue mismatch: missing={sorted(set(expected_tissues)-set(tissues))[:5]} "
            f"unexpected={sorted(set(tissues)-set(expected_tissues))[:5]}"
        )
    return files


def command(args, spredixcan_folder, covariance, cleared, raw_output):
    # Genome-wide mode intentionally omits --throw: the vendored interface then
    # isolates per-gene Python exceptions. Omitted valid-covariance genes are
    # audited as outcome-specific unavailable genes and excluded from BH.
    return [
        str(args.python_metaxcan), str(args.interface),
        "--models_folder", str(args.models_folder),
        "--models_name_pattern", "JTI_(.*).db",
        "--metaxcan_folder", str(spredixcan_folder),
        "--metaxcan_file_name_parse_pattern", "(.*).csv",
        "--snp_covariance", str(covariance),
        "--cleared_snps", str(cleared),
        "--cutoff_eigen_ratio", "1e-6",
        "--trimmed_ensemble_id",
        "--output", str(raw_output),
    ]


def parse_raw(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise RuntimeError(f"unexpected S-MultiXcan schema: {reader.fieldnames}")
        rows = list(reader)
    genes = [normalized_gene(x["gene"]) for x in rows]
    if len(genes) != len(set(genes)):
        raise RuntimeError("duplicate genes in S-MultiXcan output")
    return rows


def compare_gene_universes(covariance_genes, returned_genes):
    """Return outcome-unavailable and invalid unexpected genes.

    Outcome-specific omissions are allowed; returned genes must remain a subset
    of the valid common-reference covariance genes.
    """
    return (
        sorted(set(covariance_genes) - set(returned_genes)),
        sorted(set(returned_genes) - set(covariance_genes)),
    )


def enrich_rows(dataset_id, rows):
    enriched = []
    for row in rows:
        classification = classify_smultixcan(row)
        condition = "NA"
        try:
            emax = float(row["eigen_max"])
            emin = float(row["eigen_min_kept"])
            if math.isfinite(emax) and math.isfinite(emin) and emin > 0:
                condition = format(emax / emin, ".17g")
        except (TypeError, ValueError):
            pass
        enriched.append(
            {
                "dataset_id": dataset_id,
                **row,
                "condition_number": condition,
                "classification": classification,
                "qvalue_global_bh": "NA",
                "global_bh_significant": "false",
                "direction_available": "false",
                "signed_statistic": "NA",
                "direction": "NA",
            }
        )
    successful = [x for x in enriched if x["classification"] == "success"]
    for row, qvalue in zip(successful, bh([float(x["pvalue"]) for x in successful])):
        row["qvalue_global_bh"] = format(qvalue, ".17g")
        row["global_bh_significant"] = str(qvalue < 0.05).lower()
    enriched.sort(
        key=lambda x: (
            0 if x["classification"] == "success" else 1,
            float(x["pvalue"]) if x["classification"] == "success" else math.inf,
            x["gene"],
        )
    )
    return enriched


def write_enriched(path: Path, rows):
    fields = ["dataset_id"] + OUTPUT_COLUMNS + EXTRA_FIELDS[1:]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def run_dataset(args, entry, covariance, retained, covariance_genes, covariance_method_hash, model_tissues):
    dataset_id = entry["dataset_id"]
    normalized_gwas = Path(entry["normalized_gwas"]).resolve()
    spredixcan = Path(entry["spredixcan_folder"]).resolve()
    if not normalized_gwas.is_file() or not spredixcan.is_dir():
        raise RuntimeError(f"dataset inputs absent: {dataset_id}")
    tissue_files = validate_spredixcan_folder(spredixcan, model_tissues)
    output = args.output_root / dataset_id
    output.mkdir(parents=True, exist_ok=False)
    cleared = output / "cleared_snps.tsv"
    n_cleared = write_cleared(normalized_gwas, retained, cleared)
    raw = output / "smultixcan.raw.tsv"
    cmd = command(args, spredixcan, covariance, cleared, raw)
    completed = subprocess.run(cmd, cwd=args.project_root, capture_output=True, text=True)
    (output / "stdout.txt").write_text(completed.stdout)
    (output / "stderr.txt").write_text(completed.stderr)
    if completed.returncode or not raw.is_file():
        raise RuntimeError(f"S-MultiXcan failed for {dataset_id}: exit={completed.returncode}; {completed.stderr[-2000:]}")
    rows = parse_raw(raw)
    enriched = enrich_rows(dataset_id, rows)
    result = output / "smultixcan_all_genes.tsv"
    write_enriched(result, enriched)
    returned = {normalized_gene(x["gene"]): x for x in enriched}
    audit = []
    for gene in sorted(covariance_genes | set(returned)):
        row = returned.get(gene)
        audit.append(
            {
                "gene": gene,
                "raw_row_present": str(row is not None).lower(),
                "classification": row["classification"] if row else "outcome_unavailable",
                "status": row["status"] if row else "NA",
            }
        )
    audit_path = output / "gene_execution_audit.tsv"
    with audit_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(audit)
    successful = [x for x in enriched if x["classification"] == "success"]
    missing, unexpected = compare_gene_universes(covariance_genes, returned)
    assigned_q = [x for x in enriched if x["qvalue_global_bh"] != "NA"]
    checks = {
        "interface_exit_zero": completed.returncode == 0,
        "exact_output_schema": True,
        "unique_output_gene_rows": len(returned) == len(rows),
        "returned_genes_subset_of_valid_covariance_genes": not unexpected,
        "valid_covariance_gene_universe_audited": len(audit) == len(covariance_genes | set(returned)),
        "successful_cross_tissue_family_nonempty": bool(successful),
        "all_success_rows_status_zero": all(x["status"] in {"0", "0.0"} for x in successful),
        "all_success_rows_n_ge_2": all(int(float(x["n"])) >= 2 for x in successful),
        "all_success_rows_n_indep_ge_2": all(int(float(x["n_indep"])) >= 2 for x in successful),
        "all_success_rows_have_global_bh": all(x["qvalue_global_bh"] != "NA" for x in successful),
        "bh_family_qvalue_closed": len(successful) == len(assigned_q),
        "all_rows_unsigned": all(x["direction_available"] == "false" and x["direction"] == "NA" for x in enriched),
    }
    validation = {
        "dataset_id": dataset_id,
        "stage": "M1.5_genome_wide_smultixcan",
        "status": "pass" if all(checks.values()) else "fail",
        "checks": checks,
        "counts": {
            "covariance_genes": len(covariance_genes), "returned_rows": len(rows),
            "successful_cross_tissue_genes": len(successful),
            "global_bh_family_size": len(successful),
            "global_bh_qvalues_assigned": len(assigned_q),
            "global_bh_significant_q_lt_0.05": sum(float(x["qvalue_global_bh"]) < 0.05 for x in successful),
            "cleared_rsids": n_cleared, "missing_gene_rows": len(missing),
            "unexpected_gene_rows": len(unexpected),
        },
        "classification_counts": {
            value: sum(x["classification"] == value for x in enriched)
            for value in sorted({x["classification"] for x in enriched})
        },
        "missing_gene_examples": missing[:20], "unexpected_gene_examples": unexpected[:20],
    }
    validation_path = output / "validation-report.json"
    validation_path.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    method = {
        "dataset_id": dataset_id, "stage": "M1.5_genome_wide_smultixcan",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "SMulTiXcan-core/eigen_ratio-1e-6/genome-wide-guarded-v1",
        "multiple_testing": "BH across every successful nondegenerate S-MultiXcan gene for this GWAS",
        "direction_available": False,
        "exception_policy": "vendored per-gene catch; omitted valid-covariance genes are outcome_unavailable and excluded from BH",
        "command": cmd,
        "inputs": {
            "normalized_gwas": {"path": str(normalized_gwas), "sha256": sha256(normalized_gwas)},
            "spredixcan_files": [{"path": str(x), "sha256": sha256(x)} for x in tissue_files],
            "covariance_method_sha256": covariance_method_hash,
            "interface": {"path": str(args.interface), "sha256": sha256(args.interface)},
            "genome_runner": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
            "cleared_snps": {"path": str(cleared), "sha256": sha256(cleared)},
        },
        "outputs": {
            "raw": {"path": str(raw), "sha256": sha256(raw)},
            "all_genes": {"path": str(result), "sha256": sha256(result)},
            "gene_execution_audit": {"path": str(audit_path), "sha256": sha256(audit_path)},
            "validation_report": {"path": str(validation_path), "sha256": sha256(validation_path)},
        },
    }
    (output / "method.json").write_text(json.dumps(method, indent=2, sort_keys=True) + "\n")
    return validation


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--project-root", type=Path, required=True)
    result.add_argument("--dataset-manifest", type=Path, required=True)
    result.add_argument("--covariance-resource", type=Path, required=True)
    result.add_argument("--models-folder", type=Path, required=True)
    result.add_argument("--python-metaxcan", type=Path, required=True)
    result.add_argument("--interface", type=Path, required=True)
    result.add_argument("--output-root", type=Path, required=True)
    result.add_argument("--continue-on-dataset-failure", action="store_true")
    result.add_argument("--allow-test-subset", action="store_true", help="Permit a covariance resource explicitly labelled test_subset")
    return result


def main():
    args = parser().parse_args()
    for key in ("project_root", "dataset_manifest", "covariance_resource", "models_folder", "python_metaxcan", "interface", "output_root"):
        setattr(args, key, getattr(args, key).resolve())
    if args.output_root.exists() and any(args.output_root.iterdir()):
        raise RuntimeError(f"output root must be new or empty: {args.output_root}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    covariance_method = args.covariance_resource / "method.json"
    covariance_validation = args.covariance_resource / "validation-report.json"
    covariance = args.covariance_resource / "snp_covariance.txt.gz"
    retained_path = args.covariance_resource / "retained_rsids.tsv.gz"
    gene_qc = args.covariance_resource / "gene_qc.tsv.gz"
    for path in (args.dataset_manifest, covariance_method, covariance_validation, covariance, retained_path, gene_qc, args.python_metaxcan, args.interface):
        if not path.exists():
            raise RuntimeError(f"required input absent: {path}")
    if json.loads(covariance_validation.read_text()).get("status") != "pass":
        raise RuntimeError("covariance resource validation is not pass")
    method = json.loads(covariance_method.read_text())
    if method.get("scope") != "genome_wide" and not args.allow_test_subset:
        raise RuntimeError("formal runner requires a full chromosomes 1-22 genome_wide covariance resource")
    for key, path in (
        ("covariance", covariance), ("retained_rsids", retained_path),
        ("gene_qc", gene_qc), ("validation_report", covariance_validation),
    ):
        if sha256(path) != method["outputs"][key]["sha256"]:
            raise RuntimeError(f"covariance resource hash mismatch: {key}")
    model_tissues = set()
    for item in method["inputs"]["model_databases"]:
        current = args.models_folder / Path(item["path"]).name
        if not current.is_file() or sha256(current) != item["sha256"]:
            raise RuntimeError(f"JTI model differs from covariance resource: {current}")
        model_tissues.add(item["tissue"])
    retained = read_retained_rsids(retained_path)
    covariance_genes = read_covariance_genes(gene_qc)
    entries = read_manifest(args.dataset_manifest)
    validations = []
    failures = []
    for entry in entries:  # Deliberately serial; do not parallelize GWAS outcomes.
        try:
            validations.append(
                run_dataset(args, entry, covariance, retained, covariance_genes, sha256(covariance_method), model_tissues)
            )
        except Exception as error:
            failures.append({"dataset_id": entry["dataset_id"], "error": str(error)})
            if not args.continue_on_dataset_failure:
                raise
    overall = {
        "status": "pass" if not failures and validations and all(x["status"] == "pass" for x in validations) else "fail",
        "dataset_order": [x["dataset_id"] for x in entries],
        "datasets_completed": [x["dataset_id"] for x in validations],
        "failures": failures,
        "sequential_execution": True,
    }
    (args.output_root / "validation-report.json").write_text(json.dumps(overall, indent=2, sort_keys=True) + "\n")
    if overall["status"] != "pass":
        raise RuntimeError("one or more genome-wide S-MultiXcan datasets failed validation")
    print(json.dumps(overall, sort_keys=True))


if __name__ == "__main__":
    main()
