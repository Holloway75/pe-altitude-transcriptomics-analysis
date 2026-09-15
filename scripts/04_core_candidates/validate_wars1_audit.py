#!/usr/bin/env python3
"""Validator for the frozen Module 4 WARS1 post-selection audit release.

Frozen release:
    results/04_core_candidates/analyses/m4_wars1_post_selection_audit_v1_20260802/

Checks
------
1. manifest_consistent            every artifact_manifest.tsv row exists with matching
                                  bytes and SHA-256; no unexpected extra files.
2. release_metadata_valid         formal_release.json parses; analysis_id matches the
                                  release directory; required decision fields present.
3. coloc_recomputed               coloc.abf is rerun from the frozen harmonized variants
                                  for all three priors and reproduces the stored PP.H4
                                  values (requires Rscript with the coloc package).
4. region_qc_consistent           region QC metrics recompute from the harmonized file.
5. nearby_wars1_matches_standard  the WARS1 row of the nearby-gene table equals the
                                  standard-prior row of the sensitivity table.
6. decision_recorded              the decision table still records do_not_upgrade.
7. qc_stress_frozen               the frozen QC-stress supplement (audit A-3,
                                  2026-09-10) exists with hash closure, the
                                  relaxed+tightened 48-cell grid, a uniformly
                                  negative overall direction, and the 6-allocation
                                  exact permutation table.

Reruns of the audit scripts themselves write to work/ and never into the release;
this validator only reads the release plus one temporary R script and writes the
validation_report.json.

Usage:
    python3 validate_wars1_audit.py [--release-dir DIR] [--qc-stress-dir DIR]
                                    [--rscript PATH] [--write-report]
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ANALYSIS_ID = "m4_wars1_post_selection_audit_v1_20260802"
EXPECTED_FILES_EXTRA = {"artifact_manifest.tsv", "validation_report.json"}

R_RECOMPUTE = r"""
suppressPackageStartupMessages({library(data.table); library(coloc)})
args <- commandArgs(trailingOnly = TRUE)
rel <- args[[1]]; out <- args[[2]]
x <- fread(file.path(rel, "colocalization", "wars1_coloc_harmonized.tsv"))
setorder(x, BP, SNP)
d1 <- list(beta = x$beta_exp, varbeta = x$se_exp^2, snp = x$SNP,
           position = x$BP, type = "quant", N = 670, MAF = x$maf)
d2 <- list(beta = x$beta_out_aligned, varbeta = x$se_out^2, snp = x$SNP,
           position = x$BP, type = "cc", N = 611484, s = 16349/611484, MAF = x$maf)
rows <- list()
for (p12 in c(1e-6, 1e-5, 1e-4)) {
  fit <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12)
  sm <- fit$summary
  rows[[length(rows) + 1]] <- data.table(
    p12 = p12, nsnps = as.integer(sm[["nsnps"]]),
    pp_h0 = as.numeric(sm[["PP.H0.abf"]]), pp_h1 = as.numeric(sm[["PP.H1.abf"]]),
    pp_h2 = as.numeric(sm[["PP.H2.abf"]]), pp_h3 = as.numeric(sm[["PP.H3.abf"]]),
    pp_h4 = as.numeric(sm[["PP.H4.abf"]]))
}
fwrite(rbindlist(rows), out, sep = "\t")
cat("coloc version:", as.character(packageVersion("coloc")), "\n")
"""


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class CheckList:
    def __init__(self):
        self.checks = {}
        self.failures = []

    def record(self, name: str, ok: bool, detail: str = ""):
        self.checks[name] = bool(ok)
        if not ok:
            self.failures.append(f"{name}: {detail}" if detail else name)
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail and not ok else ""))


def check_manifest(release: Path, cl: CheckList):
    manifest = read_tsv(release / "artifact_manifest.tsv")
    problems = []
    listed = set()
    for row in manifest:
        rel = row["relative_path"]
        listed.add(rel)
        target = release / rel
        if not target.is_file():
            problems.append(f"missing {rel}")
            continue
        if target.stat().st_size != int(row["bytes"]):
            problems.append(f"size mismatch {rel}")
        if sha256_of(target) != row["sha256"]:
            problems.append(f"sha256 mismatch {rel}")
    actual = {
        str(p.relative_to(release))
        for p in release.rglob("*")
        if p.is_file()
    }
    unexpected = actual - listed - EXPECTED_FILES_EXTRA
    for extra in sorted(unexpected):
        problems.append(f"unexpected file {extra}")
    cl.record("manifest_consistent", not problems, "; ".join(problems))
    return not problems


def check_metadata(release: Path, cl: CheckList):
    try:
        meta = json.loads((release / "formal_release.json").read_text())
        ok = (
            meta.get("analysis_id") == ANALYSIS_ID
            and meta.get("status") == "FROZEN_POST_SELECTION_AUDIT"
            and meta.get("overall_decision", {}).get("wars1_upgraded_to_core_candidate") is False
            and meta.get("overall_decision", {}).get("candidate_genes_tsv_generated") is False
            and release.name == ANALYSIS_ID
        )
        cl.record("release_metadata_valid", ok, "unexpected metadata content")
    except Exception as exc:  # noqa: BLE001
        cl.record("release_metadata_valid", False, str(exc))


def check_coloc_recompute(release: Path, cl: CheckList, rscript: str):
    stored = {
        float(r["p12"]): r
        for r in read_tsv(release / "colocalization" / "wars1_coloc_abf_prior_sensitivity.tsv")
    }
    with tempfile.TemporaryDirectory() as tmp:
        r_path = Path(tmp) / "recompute.R"
        out_path = Path(tmp) / "recomputed.tsv"
        r_path.write_text(R_RECOMPUTE)
        proc = subprocess.run(
            [rscript, str(r_path), str(release), str(out_path)],
            capture_output=True, text=True,
        )
        if proc.returncode != 0 or not out_path.is_file():
            cl.record("coloc_recomputed", False,
                      (proc.stderr or proc.stdout or "no output")[-400:])
            return
        recomputed = read_tsv(out_path)
    problems = []
    if len(recomputed) != len(stored):
        problems.append(f"row count {len(recomputed)} != {len(stored)}")
    for row in recomputed:
        ref = stored.get(float(row["p12"]))
        if ref is None:
            problems.append(f"unexpected p12 {row['p12']}")
            continue
        if int(row["nsnps"]) != int(ref["nsnps"]):
            problems.append(f"nsnps mismatch at p12={row['p12']}")
        for key in ("pp_h0", "pp_h1", "pp_h2", "pp_h3", "pp_h4"):
            if abs(float(row[key]) - float(ref[key])) > 1e-12:
                problems.append(f"{key} mismatch at p12={row['p12']}")
    cl.record("coloc_recomputed", not problems, "; ".join(problems))


def check_region_qc(release: Path, cl: CheckList):
    qc = {r["metric"]: r["value"] for r in read_tsv(release / "colocalization" / "wars1_coloc_region_qc.tsv")}
    harm = read_tsv(release / "colocalization" / "wars1_coloc_harmonized.tsv")
    problems = []
    if str(len(harm)) != str(qc.get("harmonized_maf_gt_0.01")):
        problems.append("harmonized row count mismatch")
    lead_eqtl = min(harm, key=lambda r: float(r["p_exp"]))
    lead_gwas = min(harm, key=lambda r: float(r["p_out"]))
    if abs(float(lead_eqtl["p_exp"]) - float(qc["min_eqtl_p"])) > 1e-9 * float(qc["min_eqtl_p"]):
        problems.append("min_eqtl_p mismatch")
    if lead_eqtl["SNP"] != qc["lead_eqtl_snp"]:
        problems.append("lead_eqtl_snp mismatch")
    if abs(float(lead_gwas["p_out"]) - float(qc["min_gwas_p"])) > 1e-9 * float(qc["min_gwas_p"]):
        problems.append("min_gwas_p mismatch")
    if lead_gwas["SNP"] != qc["lead_gwas_snp"]:
        problems.append("lead_gwas_snp mismatch")
    if qc["region_start"] != "100300125" or qc["region_end"] != "101300125" or qc["region_chr"] != "14":
        problems.append("region definition mismatch")
    cl.record("region_qc_consistent", not problems, "; ".join(problems))


def check_nearby(release: Path, cl: CheckList):
    nearby = read_tsv(release / "colocalization" / "nearby_gene_coloc_standard_prior.tsv")
    prior = read_tsv(release / "colocalization" / "wars1_coloc_abf_prior_sensitivity.tsv")
    wars_nearby = next((r for r in nearby if r["gene"] == "WARS1"), None)
    standard = next((r for r in prior if r["prior_id"] == "standard"), None)
    ok = (
        wars_nearby is not None
        and standard is not None
        and wars_nearby["nsnps"] == standard["nsnps"]
        and abs(float(wars_nearby["pp_h4"]) - float(standard["pp_h4"])) < 1e-12
        and abs(float(wars_nearby["pp_h3"]) - float(standard["pp_h3"])) < 1e-12
        and len(nearby) == 4
    )
    cl.record("nearby_wars1_matches_standard", ok, "nearby-gene table inconsistent with standard prior")


def check_decision(release: Path, cl: CheckList):
    rows = read_tsv(release / "three_low_cost_validation_summary.tsv")
    overall = next((r for r in rows if r["validation"] == "overall_candidate_decision"), None)
    coloc_row = next((r for r in rows if r["validation"] == "Whole_Blood_eQTL_PE_colocalization"), None)
    ok = (
        overall is not None
        and overall["result"] == "do_not_upgrade"
        and coloc_row is not None
        and coloc_row["result"] == "inconclusive_negative"
    )
    cl.record("decision_recorded", ok, "decision table no longer records do_not_upgrade")


QC_STRESS_ID = "m4_wars1_qc_stress_v2_20260910"


def check_qc_stress(qc_dir: Path, cl: CheckList):
    """Validate the frozen QC-stress supplement (audit A-3, 2026-09-10).

    Verifies the manifest hash closure, the presence of the baseline and
    relaxation grid levels, the sign of the overall direction in every grid
    cell, and the 6-allocation exact permutation table.
    """
    problems = []
    if not qc_dir.is_dir():
        cl.record("qc_stress_frozen", False, f"missing directory {qc_dir}")
        return
    manifest_path = qc_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        cl.record("qc_stress_frozen", False, f"manifest unreadable: {exc}")
        return
    if manifest.get("analysis_id") != qc_dir.name or manifest.get("analysis_id") != QC_STRESS_ID:
        problems.append("analysis_id mismatch")
    recorded = manifest.get("outputs", {})
    present = {p.name for p in qc_dir.iterdir() if p.is_file()} - {"manifest.json", "validation_report.json"}
    if set(recorded) != present:
        problems.append(f"file closure mismatch: {sorted(set(recorded) ^ present)}")
    for name, meta in recorded.items():
        p = qc_dir / name
        if not p.is_file() or p.stat().st_size != meta["size_bytes"] or sha256_of(p) != meta["sha256"]:
            problems.append(f"digest mismatch: {name}")
    grid = read_tsv(qc_dir / "qc_threshold_grid.tsv") if (qc_dir / "qc_threshold_grid.tsv").is_file() else []
    if len(grid) != 48:
        problems.append(f"grid has {len(grid)} rows, expected 48")
    genes = {int(r["min_genes"]) for r in grid}
    libs = {int(r["min_lib"]) for r in grid}
    mts = {round(float(r["max_mt"]), 2) for r in grid}
    if not {100, 200, 500, 1000} <= genes or not {250, 500, 1000, 2500} <= libs or not {0.10, 0.20, 0.30} <= mts:
        problems.append("grid lacks relaxed or tightened levels")
    baseline = [r for r in grid if r["min_genes"] == "200" and r["min_lib"] == "500" and round(float(r["max_mt"]), 2) == 0.20]
    if len(baseline) != 1 or abs(float(baseline[0]["log2_pe_vs_control"]) - (-1.17)) >= 0.005:
        problems.append("baseline row does not reproduce frozen section 3.5 value -1.17")
    finite = [float(r["log2_pe_vs_control"]) for r in grid if r["log2_pe_vs_control"] not in ("", "NA", None)]
    if len(finite) != len(grid) or not all(v < 0 for v in finite):
        problems.append("non-finite or non-negative grid directions present")
    perm = read_tsv(qc_dir / "exact_donor_permutation.tsv") if (qc_dir / "exact_donor_permutation.tsv").is_file() else []
    if len(perm) != 6 or sum(1 for r in perm if r["matches_observed"] == "TRUE") != 1:
        problems.append("exact permutation table is not the 6 allocations with one observed match")
    cl.record("qc_stress_frozen", not problems, "; ".join(problems))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    default_release = Path(__file__).resolve().parents[2] / (
        "results/04_core_candidates/analyses/" + ANALYSIS_ID
    )
    ap.add_argument("--release-dir", default=str(default_release))
    ap.add_argument("--qc-stress-dir", default=str(
        Path(__file__).resolve().parents[2] / "results/04_core_candidates/analyses" / QC_STRESS_ID
    ), help="frozen QC-stress supplement directory (audit A-3)")
    ap.add_argument("--rscript", default=shutil.which("Rscript") or "Rscript")
    ap.add_argument("--write-report", action="store_true",
                    help="write validation_report.json into the release directory")
    args = ap.parse_args()

    release = Path(args.release_dir).resolve()
    if not release.is_dir():
        print(f"ERROR: release directory not found: {release}", file=sys.stderr)
        return 1

    print(f"Validating {release}")
    cl = CheckList()
    check_manifest(release, cl)
    check_metadata(release, cl)
    check_coloc_recompute(release, cl, args.rscript)
    check_region_qc(release, cl)
    check_nearby(release, cl)
    check_decision(release, cl)
    check_qc_stress(Path(args.qc_stress_dir).resolve(), cl)

    report = {
        "analysis_id": ANALYSIS_ID,
        "validated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "release_dir": str(release),
        "checks": cl.checks,
        "checks_passed": sum(cl.checks.values()),
        "checks_total": len(cl.checks),
        "failures": cl.failures,
        "status": "PASS" if not cl.failures else "FAIL",
    }
    if args.write_report:
        out = release / "validation_report.json"
        out.write_text(json.dumps(report, indent=1) + "\n")
        print(f"report written: {out}")
    else:
        print(json.dumps(report, indent=1))
    return 0 if not cl.failures else 1


if __name__ == "__main__":
    sys.exit(main())
