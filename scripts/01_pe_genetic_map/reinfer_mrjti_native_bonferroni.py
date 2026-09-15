#!/usr/bin/env python3
"""Recompute native MR-JTI Bonferroni CIs from a validated frozen bootstrap run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
from collections import Counter, defaultdict
from pathlib import Path

NA = "NA"


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row.get(key, NA) for key in fields} for row in rows)
    os.replace(tmp, path)


def write_json(path: Path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quantile_type7(values, probability):
    ordered = sorted(values)
    index = (len(ordered) - 1) * probability
    low = math.floor(index)
    high = math.ceil(index)
    if low == high:
        return ordered[low]
    return ordered[low] + (index - low) * (ordered[high] - ordered[low])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    if output.exists():
        raise SystemExit(f"output already exists: {output}")
    source_validation = json.loads((source / "independent-validation-report.json").read_text())
    if source_validation.get("status") != "pass":
        raise SystemExit("source run lacks passing independent validation")
    shutil.copytree(source, output, copy_function=os.link)

    manifest = read_tsv(output / "candidate-manifest.tsv")
    results = read_tsv(output / "mrjti-global-bh.tsv")
    family_sizes = Counter(row["dataset_id"] for row in manifest)
    expected = {(row["dataset_id"], row["tissue"], row["gene_id"]) for row in manifest}
    if expected != {(row["dataset_id"], row["tissue"], row["gene_id"]) for row in results}:
        raise SystemExit("source candidate/result closure failed")

    statuses = Counter()
    for row in results:
        dataset = row["dataset_id"]
        family_size = family_sizes[dataset]
        family_alpha = 0.05 / family_size
        task_dir = output / "tasks" / dataset / row["tissue"] / row["gene_id"]
        row["task_dir"] = str(task_dir)
        row["mrjti_family_size"] = family_size
        row["mrjti_family_alpha"] = format(family_alpha, ".17g")
        row["pvalue"] = row["q_mrjti"] = NA
        if row["status"] == "success_inference_descriptive":
            draws = [float(x["expression_beta"]) for x in read_tsv(task_dir / "bootstrap.tsv")]
            result_path = task_dir / "mrjti-result.tsv"
            result = read_tsv(result_path)[0]
            original = float(result["original_lasso_coefficient"])
            upper_quantile = quantile_type7(draws, 1 - family_alpha / 2)
            lower_quantile = quantile_type7(draws, family_alpha / 2)
            ci_low = 2 * original - upper_quantile
            ci_high = 2 * original - lower_quantile
            significance = "sig" if ci_low * ci_high > 0 else "nonsig"
            result.update({
                "ci_low": format(ci_low, ".17g"), "ci_high": format(ci_high, ".17g"),
                "ci_method": "vendored_basic_residual_bootstrap_bonferroni",
                "ci_significance": significance, "n_genes": str(family_size),
                "family_alpha": format(family_alpha, ".17g"),
                "inference_status": "native_bonferroni_ci_selection_dependent",
            })
            write_tsv(result_path, list(result), [result])
            row.update({
                "status": "success_inference_bonferroni_selection_dependent",
                "ci_low": result["ci_low"], "ci_high": result["ci_high"],
                "ci_method": result["ci_method"], "ci_significance": significance,
                "mrjti_supported": str(significance == "sig").lower(),
                "inference_status": "native_bonferroni_ci_selection_dependent",
            })
        else:
            row["ci_significance"] = row["mrjti_supported"] = NA
        task_result = task_dir / "task-result.json"
        write_json(task_result, row)
        statuses[row["status"]] += 1

    fields = list(results[0])
    anchor = fields.index("mrjti_supported") + 1
    for field in reversed(("ci_significance", "mrjti_family_size", "mrjti_family_alpha")):
        if field not in fields:
            fields.insert(anchor, field)
    write_tsv(output / "mrjti-global-bh.tsv", fields, results)

    method = json.loads((output / "method.json").read_text())
    project_root = Path(__file__).resolve().parents[2]
    spredixcan_runs = {
        "FIGSHARE_22680904_v2": "figshare_22680904_v2_spredixcan_hg19_v1_20260718",
        "GCST90269904": "gcst90269904_spredixcan_hg19_v1_20260718",
        "GCST90301704": "gcst90301704_spredixcan_hg19_v1_20260718",
        "GCST90454233": "gcst90454233_spredixcan_hg19_v1_20260718",
    }
    for evidence in method["source_evidence"]:
        dataset = evidence["dataset_id"]
        parent = project_root / "results/01_pe_genetic_map/analyses" / spredixcan_runs[dataset] / "M1.3_spredixcan" / dataset
        local_paths = {
            "table": parent / "spredixcan_all_tissues.tsv",
            "method": parent / "method.json",
            "validation": parent / "validation-report.json",
            "normalized_gwas": project_root / "data/processed/01_pe_genetic_map/M1.1_prepare_gwas" / dataset / "metaxcan.tsv.gz",
        }
        for key, path in local_paths.items():
            if not path.is_file():
                raise SystemExit(f"missing rebound input: {path}")
            evidence[key] = str(path)
            evidence[key + "_sha256"] = sha256(path)
    protocol = project_root / "docs/protocol/SNP_MATCHING_RULES.md"
    method["matching_protocol"]["path"] = str(protocol)
    method["matching_protocol"]["sha256"] = sha256(protocol)
    method.update({
        "schema_version": "m1.6-mrjti-native-bonferroni-v2",
        "source_scientific_bundle_sha256": source_validation["scientific_bundle_sha256"],
        "selection_dependency": "same outcome used for candidate selection and MR-JTI; Bonferroni significance is conditional on the selected family and is not independent causal confirmation",
    })
    method["candidate_manifest"] = str(output / "candidate-manifest.tsv")
    method["mrjti"] = {
        "folds": 5, "bootstrap": 500, "seed": 20260714, "minimum_snps": 20,
        "family_definition": "planned selected gene-by-tissue tasks within each GWAS, including unavailable tasks",
        "family_sizes": dict(sorted(family_sizes.items())),
        "ci": "vendored basic residual bootstrap with native Bonferroni alpha=0.05/n_genes",
        "pvalue": None, "qvalue": None, "support": "CI_significance sig/nonsig; selection-dependent",
    }
    method["formal_publication_modified"] = False
    write_json(output / "method.json", method)
    write_json(output / "validation-report.json", {
        "status": "pending_independent_validation", "errors": [], "candidate_tasks": len(manifest),
        "result_tasks": len(results), "status_counts": dict(statuses),
        "native_bonferroni_ci_available": True, "calibrated_p_or_q_available": False,
    })
    (output / "independent-validation-report.json").unlink(missing_ok=True)
    print(json.dumps({"family_sizes": dict(family_sizes), "status_counts": dict(statuses)}, sort_keys=True))


if __name__ == "__main__":
    main()
