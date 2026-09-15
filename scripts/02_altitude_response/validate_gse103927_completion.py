#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_manifest(path: Path, base: Path, checks: list):
    for row in rows(path):
        target = Path(row.get("path", "")) if "path" in row else base / row["file"]
        if not target.is_absolute():
            target = (base / target).resolve()
        checks.append((target.exists(), f"manifest target exists: {target}"))
        if target.exists():
            checks.append(
                (target.stat().st_size == int(row["bytes"]),
                 f"manifest byte count: {target}")
            )
            checks.append(
                (sha256(target) == row["sha256"],
                 f"manifest SHA-256: {target}")
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse103927-dir", required=True, type=Path)
    parser.add_argument("--comparison-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    gdir = args.gse103927_dir.resolve()
    cdir = args.comparison_dir.resolve()
    checks = []

    required_g = {
        "ANALYSIS_SUMMARY.md", "analysis_registry.json",
        "artifact_manifest.tsv", "cell_composition_proxy_changes.tsv",
        "cell_composition_proxy_summary.tsv",
        "composition_proxy_sensitivity.tsv", "contrast_summary.tsv",
        "differential_expression.tsv", "donor_altitude_mapping.tsv",
        "leave_one_donor_out.tsv", "model_objects.rds", "pca.png",
        "probe_mapping_audit.tsv", "sample_metadata.tsv", "sample_qc.tsv",
        "source_manifest.tsv", "top_genes.tsv", "volcano.png",
    }
    required_c = {
        "ANALYSIS_SUMMARY.md", "analysis_registry.json",
        "artifact_manifest.tsv", "composition_proxy_direction_agreement.tsv",
        "composition_proxy_gene_direction_comparison.tsv",
        "gene_direction_agreement.tsv", "gene_direction_comparison.tsv",
        "gene_effect_scatter.png", "pathway_direction_agreement.tsv",
        "pathway_direction_comparison.tsv", "pathway_effect_scatter.png",
        "significant_gene_overlap.tsv", "source_manifest.tsv",
    }
    checks.append((required_g <= {p.name for p in gdir.iterdir()},
                   "GSE103927 required artifacts"))
    checks.append((required_c <= {p.name for p in cdir.iterdir()},
                   "comparison required artifacts"))

    de = rows(gdir / "differential_expression.tsv")
    summary = rows(gdir / "contrast_summary.tsv")[0]
    checks.append((len(de) == 18224, "GSE103927 tested-gene count"))
    checks.append((len({r["gene_symbol"] for r in de}) == len(de),
                   "GSE103927 unique gene symbols"))
    checks.append((len({r["hgnc_id"] for r in de}) == len(de),
                   "GSE103927 unique HGNC IDs"))
    pvalues = [float(r["P.Value"]) for r in de]
    fdr = [float(r["adj.P.Val"]) for r in de]
    checks.append((all(0 <= x <= 1 for x in pvalues + fdr),
                   "GSE103927 valid P/FDR range"))
    checks.append((pvalues == sorted(pvalues),
                   "GSE103927 results sorted by raw P"))
    checks.append((sum(x < 0.05 for x in fdr) == int(summary["fdr_0.05"]),
                   "GSE103927 FDR count matches summary"))

    metadata = rows(gdir / "sample_metadata.tsv")
    pairs = rows(gdir / "donor_altitude_mapping.tsv")
    checks.append((sum(r["include_primary"] == "TRUE" for r in metadata) == 42,
                   "42 primary samples"))
    checks.append((len(pairs) == 21, "21 complete donor pairs"))
    checks.append((len({r["donor_id"] for r in pairs}) == 21,
                   "21 unique paired donors"))

    qc = rows(gdir / "sample_qc.tsv")
    checks.append((len(qc) == 42, "42 sample-QC rows"))
    checks.append((not any(r["pca_outlier_first5"] == "TRUE" for r in qc),
                   "no PCA outlier flags"))
    checks.append((min(float(r["within_pair_correlation"]) for r in qc) > 0.95,
                   "all within-pair correlations > 0.95"))

    loo = rows(gdir / "leave_one_donor_out.tsv")
    checks.append((len(loo) == len(de), "one LOO summary per tested gene"))
    de_fdr = {r["gene_symbol"] for r in de if float(r["adj.P.Val"]) < 0.05}
    loo_fdr = [r for r in loo if r["gene_symbol"] in de_fdr]
    checks.append(
        (all(math.isclose(float(r["loo_sign_concordance"]), 1.0)
             for r in loo_fdr),
         "all primary FDR genes retain direction in every LOO model")
    )

    comparison = rows(cdir / "gene_direction_comparison.tsv")
    agreement = rows(cdir / "gene_direction_agreement.tsv")
    all_agreement = next(
        r for r in agreement if r["stratum"] == "all_common_current_HGNC"
    )
    checks.append((len(comparison) == 11111, "common HGNC universe size"))
    checks.append((len({r["hgnc_id"] for r in comparison}) == len(comparison),
                   "comparison HGNC IDs unique"))
    observed_sign = sum(r["same_direction"] == "TRUE" for r in comparison)
    checks.append(
        (observed_sign == round(float(all_agreement["sign_concordance"]) *
                                len(comparison)),
         "gene sign-concordance summary matches rows")
    )
    checks.append(
        (int(all_agreement["both_fdr"]) ==
         sum(r["both_fdr"] == "TRUE" for r in comparison),
         "both-FDR count matches comparison rows")
    )

    conditional = rows(cdir / "composition_proxy_gene_direction_comparison.tsv")
    checks.append((len(conditional) == len(comparison),
                   "conditional comparison uses same universe"))
    pathways = rows(cdir / "pathway_direction_comparison.tsv")
    pathway_agreement = rows(cdir / "pathway_direction_agreement.tsv")[0]
    checks.append((len(pathways) == int(pathway_agreement["n_pathways"]),
                   "pathway count matches agreement"))

    check_manifest(gdir / "artifact_manifest.tsv", gdir, checks)
    check_manifest(gdir / "source_manifest.tsv", gdir, checks)
    check_manifest(cdir / "artifact_manifest.tsv", cdir, checks)
    check_manifest(cdir / "source_manifest.tsv", cdir, checks)

    failures = [name for passed, name in checks if not passed]
    report = {
        "status": "PASS" if not failures else "FAIL",
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failures),
        "failures": failures,
        "gse103927_dir": str(gdir),
        "comparison_dir": str(cdir),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
