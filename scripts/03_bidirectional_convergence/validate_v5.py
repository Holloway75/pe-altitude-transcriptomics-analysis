#!/usr/bin/env python3
"""Validate, summarize, and manifest the Module 3 v5 release."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--locked-dir", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()
    locked = args.locked_dir.resolve(); results = args.result_dir.resolve()
    registry_path = locked / "module3_v5_analysis_registry.json"
    universe_path = locked / "gene_universe_and_mapping.tsv"
    signs_path = locked / "donor_sign_flip_matrix.npz"
    overlap_path = results / "directional_overlap.tsv"
    overlap_genes_path = results / "directional_overlap_genes.tsv"
    boundary_path = results / "detectability_boundary.tsv"
    signatures_path = results / "hypoxia_altitude_signatures.tsv"
    core_path = results / "signature_core_genes.tsv"
    signature_registry_path = locked / "compact_signature_registry.tsv"
    required = [registry_path, universe_path, signs_path, overlap_path, overlap_genes_path, boundary_path, signatures_path, core_path, signature_registry_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing: raise SystemExit("missing output(s): " + ", ".join(missing))

    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    universe = read_tsv(universe_path); overlap = read_tsv(overlap_path)
    boundaries = read_tsv(boundary_path); signatures = read_tsv(signatures_path)
    checks: list[dict[str, object]] = []
    failures: list[str] = []
    def check(check_id: str, passed: bool, detail: str) -> None:
        checks.append({"check_id": check_id, "status": "PASS" if passed else "FAIL", "detail": detail})
        if not passed: failures.append(f"{check_id}: {detail}")

    check("universe_count", len(universe) == 8101, f"observed={len(universe)} expected=8101")
    check("universe_unique_hgnc", len({row['hgnc_id'] for row in universe}) == 8101, "HGNC IDs must be unique")
    check("overlap_rows", len(overlap) == 12, f"observed={len(overlap)} expected=12")
    check("boundary_rows", len(boundaries) == 12, f"observed={len(boundaries)} expected=12")
    check("signature_rows", len(signatures) == len(registry["compact_signatures"]), f"observed={len(signatures)} expected={len(registry['compact_signatures'])}")
    check("signature_ids", {row["signature_id"] for row in signatures} == set(registry["compact_signatures"]), "signature family equals frozen config")
    signs = np.load(signs_path, allow_pickle=False)
    matrix = signs["signs"]
    check("sign_matrix_shape", matrix.shape == (registry["sign_flip_draws"], 21), f"shape={matrix.shape}")
    check("sign_matrix_values", set(np.unique(matrix).tolist()) == {-1, 1}, f"values={np.unique(matrix).tolist()}")
    check("same_draws_all_quadrants", all(int(row["empirical_draws"]) == registry["sign_flip_draws"] for row in overlap), "all overlap rows use frozen draw count")
    p_formula_ok = all(abs(float(row["empirical_p_one_sided"]) - ((int(row["empirical_exceedances"]) + 1) / (int(row["empirical_draws"]) + 1))) < 1e-14 for row in overlap)
    check("empirical_p_plus_one", p_formula_ok, "P=(exceedances+1)/(draws+1)")
    primary = [row for row in overlap if row["scenario"] == "PRIMARY_STRICT"]
    check("primary_four_quadrants", len(primary) == 4, f"observed={len(primary)}")
    trigger = any(float(row["empirical_holm_p"]) < 0.05 for row in primary)
    check("candidate_trigger_registry", trigger == bool(registry["candidate_trigger_passed"]), f"computed={trigger} registry={registry['candidate_trigger_passed']}")
    candidate_path = results / "candidate_genes.tsv"
    check("candidate_file_stop_rule", candidate_path.exists() == trigger, f"trigger={trigger} candidate_file_exists={candidate_path.exists()}")
    check("finite_overlap_statistics", all(math.isfinite(float(row[field])) for row in overlap for field in ("fisher_p_one_sided", "empirical_p_one_sided", "empirical_holm_p", "fisher_holm_p", "null_mean_overlap")), "all required statistics finite")

    strict = {row["quadrant"]: row for row in primary}
    wars_rows = [row for row in read_tsv(overlap_genes_path) if row.get("scenario") == "PRIMARY_STRICT" and row.get("hgnc_symbol") == "WARS1"]
    check("wars1_strict_overlap_present", len(wars_rows) == 1 and wars_rows[0]["quadrant"] == "alt_up_pe_negative", f"rows={len(wars_rows)}")

    summary_path = results / "MODULE3_V5_EXPLORATORY_SUMMARY.md"
    strict_lines = []
    for quadrant in ("alt_up_pe_positive", "alt_up_pe_negative", "alt_down_pe_positive", "alt_down_pe_negative"):
        row = strict[quadrant]
        strict_lines.append(f"| `{quadrant}` | {row['altitude_set_n']} | {row['pe_set_n']} | {row['observed_overlap_n']} | {float(row['fisher_p_one_sided']):.6g} | {float(row['empirical_p_one_sided']):.6g} | {float(row['empirical_holm_p']):.6g} |")
    sig_lines = [f"| `{row['signature_id']}` | {row['formal_direction']} | {float(row['formal_p']):.6g} | {float(row['formal_bh_fdr']):.6g} | {row['context_status']} |" for row in signatures]
    wars = strict["alt_up_pe_negative"]
    conclusion = "At least one exploratory quadrant passed the frozen empirical Holm gate." if trigger else "No exploratory quadrant passed the frozen empirical Holm gate; no candidate list was generated."
    summary_path.write_text(
        "# Module 3 v5 exploratory directional-overlap summary\n\n"
        f"- Analysis ID: `{registry['analysis_id']}`\n"
        "- Status: post hoc exploratory analysis after the negative v4 confirmatory test\n"
        f"- Universe: {len(universe):,} one-to-one HGNC genes\n"
        f"- Sign-flip draws: {registry['sign_flip_draws']:,}; the full gene vector of each donor shared one sign per draw\n"
        "- Multiplicity: empirical Holm across four primary quadrants\n\n"
        "## Primary strict overlap\n\n"
        "| Quadrant | ALT set | PE set | Overlap | Fisher P | Empirical P | Empirical Holm P |\n"
        "|---|---:|---:|---:|---:|---:|---:|\n" + "\n".join(strict_lines) + "\n\n"
        f"The ALT-up/PE-negative strict overlap contained {wars['observed_overlap_n']} gene(s); WARS1 was present. {conclusion}\n\n"
        "## Frozen compact oxygen-response signatures\n\n"
        "| Signature | Direction | CAMERA P | BH FDR | GSE333506 context |\n"
        "|---|---|---:|---:|---|\n" + "\n".join(sig_lines) + "\n\n"
        "These signatures are biological positive controls for the GSE103927 altitude axis and are not additional PE-overlap tests.\n",
        encoding="utf-8"
    )

    validation_path = results / "validation_report.json"
    payload = {"analysis_id": registry["analysis_id"], "status": "PASS" if not failures else "FAIL", "checks_total": len(checks), "checks_passed": sum(row["status"] == "PASS" for row in checks), "failures": failures, "checks": checks}
    validation_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if failures: raise SystemExit("validation failed: " + "; ".join(failures))

    registry["execution_status"] = "VALIDATED_COMPLETE"
    registry["validation_report"] = {"path": str(validation_path), "sha256": sha256(validation_path)}
    registry["summary"] = {"path": str(summary_path), "sha256": sha256(summary_path)}
    registry_path.write_text(json.dumps(registry, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    manifest_rows = []
    for root in (locked, results):
        for path in sorted(root.iterdir()):
            if path.is_file() and path.name != "MANIFEST.sha256.tsv":
                manifest_rows.append({"release_area": "locked_inputs" if root == locked else "results", "file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_tsv(locked / "MANIFEST.sha256.tsv", manifest_rows)
    print(json.dumps({"status": "PASS", "checks": len(checks), "candidate_trigger_passed": trigger}))


if __name__ == "__main__":
    main()

