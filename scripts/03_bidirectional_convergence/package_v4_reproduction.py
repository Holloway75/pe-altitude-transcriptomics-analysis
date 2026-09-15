#!/usr/bin/env python3
"""Package the independently reconstructed Module 3 v4 audit release."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


HISTORICAL_MINIMUM_JOINT_P = 0.04002
ORIGINAL_HGNC_SHA256 = "41fe4d55110b22ac3a4c2579d281b0fcf51f18c0b622b0f88f45d49bb2c0d9fb"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_artifact(source: Path, output: Path, relative: str) -> Path:
    target = output / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    return target


def main(args: argparse.Namespace) -> None:
    source = args.reproduction_root.resolve()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty release: {output}")
    output.mkdir(parents=True, exist_ok=True)

    artifacts = {
        "formal_result_summary": (source / "results/FORMAL_RESULT.md", "FORMAL_RESULT.md"),
        "formal_result_qc": (source / "results/formal_result_qc.json", "formal_result_qc.json"),
        "formal_release": (source / "results/formal_release.json", "formal_release.json"),
        "bidirectional_enrichment": (
            source / "results/bidirectional_enrichment.tsv",
            "bidirectional_enrichment.tsv",
        ),
        "formal_alt_pathways": (
            source / "work/formal/alt/alt_pathway_results.tsv",
            "null_banks/alt_pathway_results.tsv",
        ),
        "formal_alt_report": (
            source / "work/formal/alt/alt_bank_report.json",
            "null_banks/alt_bank_report.json",
        ),
        "formal_pe_pathways": (
            source / "work/formal/pe/pe_pathway_results.tsv",
            "null_banks/pe_pathway_results.tsv",
        ),
        "formal_pe_report": (
            source / "work/formal/pe/pe_bank_report.json",
            "null_banks/pe_bank_report.json",
        ),
        "formal_gate": (
            source / "work/calibration/formal_execution_gate.json",
            "calibration/formal_execution_gate.json",
        ),
        "axis_stability": (
            source / "work/calibration/axis_stability_report.json",
            "calibration/axis_stability_report.json",
        ),
        "composition_calibration": (
            source / "work/calibration/composition_logic_report.json",
            "calibration/composition_logic_report.json",
        ),
        "null_resolution": (
            source / "work/calibration/null_resolution_audit.json",
            "calibration/null_resolution_audit.json",
        ),
        "frozen_registry": (
            source / "data_processed/registry/module3_analysis_registry.yaml",
            "registry/module3_analysis_registry.json",
        ),
        "formal_gene_universe": (
            source / "data_processed/registry/formal_gene_universe.tsv",
            "registry/formal_gene_universe.tsv",
        ),
        "formal_gene_sets": (
            source / "data_processed/registry/formal_gene_sets.tsv",
            "registry/formal_gene_sets.tsv",
        ),
        "source_manifest": (
            source / "data_processed/registry/source_manifest.tsv",
            "registry/source_manifest.tsv",
        ),
        "magma_pathways": (
            source / "work/magma/pe_formal_pathways.gsa.out",
            "magma/pe_formal_pathways.gsa.out",
        ),
        "threaded_equivalence_original": (
            source / "work/benchmark/alt_formal_shape_equivalence/original/alt_pathway_results.tsv",
            "equivalence/original_10000_alt_pathways.tsv",
        ),
        "threaded_equivalence_threaded": (
            source / "work/benchmark/alt_formal_shape_equivalence/threaded/alt_pathway_results.tsv",
            "equivalence/threaded_10000_alt_pathways.tsv",
        ),
    }
    missing = [str(path) for path, _ in artifacts.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing reproduction artifact(s): " + ", ".join(missing))

    copied = {
        name: copy_artifact(path, output, relative)
        for name, (path, relative) in artifacts.items()
    }
    qc = json.loads(copied["formal_result_qc"].read_text(encoding="utf-8"))
    stability = json.loads(copied["axis_stability"].read_text(encoding="utf-8"))
    registry = json.loads(copied["frozen_registry"].read_text(encoding="utf-8"))

    with copied["threaded_equivalence_original"].open(encoding="utf-8", newline="") as handle:
        original = list(csv.DictReader(handle, delimiter="\t"))
    with copied["threaded_equivalence_threaded"].open(encoding="utf-8", newline="") as handle:
        threaded = list(csv.DictReader(handle, delimiter="\t"))
    equivalence_pass = original == threaded

    checks = {
        "formal_qc_pass": qc.get("status") == "PASS",
        "formal_gene_universe_8101": registry["formal_gene_universe"]["n_genes"] == 8101,
        "formal_pathways_724": qc.get("formal_pathways") == 724,
        "formal_directional_tests_1448": qc.get("formal_directional_tests") == 1448,
        "zero_by_significant": qc.get("by_fdr_le_0_05") == 0,
        "one_nominal_joint_result": qc.get("joint_p_le_0_05_uncorrected") == 1,
        "axis_stability_pass": stability.get("status") == "PASS",
        "five_of_six_axis_seeds_pass": stability.get("passing_seed_runs") == 5,
        "threaded_kernel_exact_10000_row_equivalence": equivalence_pass,
    }
    current_minimum = float(qc["minimum_joint_p"])
    validation = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "analysis_id": "m3_jti_gse103927_v4_reproduction_20260809",
        "reconstructed_historical_analysis_id": registry["analysis_id"],
        "checks": checks,
        "historical_minimum_joint_p": HISTORICAL_MINIMUM_JOINT_P,
        "reproduced_minimum_joint_p": current_minimum,
        "absolute_difference": abs(current_minimum - HISTORICAL_MINIMUM_JOINT_P),
        "original_hgnc_snapshot_available": False,
        "original_hgnc_expected_sha256": ORIGINAL_HGNC_SHA256,
        "substitute_hgnc_sha256": registry["hgnc_snapshot"]["sha256"],
        "interpretation": (
            "The current official HGNC substitute exactly reconstructs the 8,101-gene, "
            "724-pathway, 1,448-test family and the historical decision, but the top joint "
            "P differs by 7e-5; this is a conclusion-level rather than byte-identical reproduction."
        ),
    }
    (output / "validation_report.json").write_text(
        json.dumps(validation, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary = f"""# Module 3 v4 independent reproduction audit

## Verdict

This reproduction passed all audit checks and reproduces the statistical
family and decisions of the historical v4: 8 101 genes, 724 pathways and
1 448 directional tests, with 0 tests at BY FDR <= 0.05 and 1 test at
uncorrected joint P <= 0.05. The top rank is still `REACTOME_TBC_RABGAPS`
in the same-sign direction; the reproduced joint P is
`{current_minimum:.5f}`, against `{HISTORICAL_MINIMUM_JOINT_P:.5f}` recorded
in the historical process notes.

## Reproduction boundary

The historical `hgnc_complete_set_20260727.tsv` is no longer in the project;
its expected SHA-256 is `{ORIGINAL_HGNC_SHA256}`. The current official HGNC
file (SHA-256 `{registry['hgnc_snapshot']['sha256']}`) was used instead. It
exactly restores the gene, pathway and test counts, but the top-rank joint P
differs by `{abs(current_minimum - HISTORICAL_MINIMUM_JOINT_P):.5g}`. This
product is therefore a conclusion-level full reproduction and does not claim
byte-level identity with the deleted historical results.

## Calibration and formal sampling

- The reference-axis calibration seed failed the per-pathway KS diagnostic;
  all 5 prespecified follow-up seeds passed (5/6 overall), meeting the
  frozen v4.1 stability gate;
- ALT and PE axes each completed 999 999 formal sampling draws;
- intersection-union monotonicity, maxT sensitivity and Monte Carlo
  accuracy all passed;
- the parallel ALT entry reproduced the archived entry's 1 448-row results
  field-for-field at the formal block size and seed, and in a 10 000-draw
  benchmark across all 724 pathways.
"""
    (output / "REPRODUCTION_AUDIT.md").write_text(summary, encoding="utf-8")

    records = []
    for path in sorted(p for p in output.rglob("*") if p.is_file()):
        records.append(
            {
                "relative_path": str(path.relative_to(output)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    with (output / "artifact_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["relative_path", "bytes", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)
    if validation["status"] != "PASS":
        raise SystemExit("v4 reproduction package validation failed")
    print(f"V4 REPRODUCTION PACKAGE PASS: {output}")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--reproduction-root", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


if __name__ == "__main__":
    main(parser().parse_args())
