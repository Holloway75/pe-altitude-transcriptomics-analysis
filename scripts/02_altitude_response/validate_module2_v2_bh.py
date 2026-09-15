#!/usr/bin/env python3
"""Independent value-level validator for the CURRENT Module 2 deliverables (v2).

Written 2026-09-10 for pre-submission audit B-10. The legacy
validate_module2.py/seal_module2.py pair targets a retired payload schema whose
whitelist contains the four datasets now excluded from the study
(GSE196728/GSE75665/GSE103940/GSE46480); it cannot validate the v2 axes and is
kept only as a historical record. This script fills the gap for the current
frozen deliverables:

1. manifest_closure     every file in the v2 input-freeze source_manifest.tsv
                        exists with matching bytes and sha256.
2. bh_recomputed        BH FDR is independently recomputed from P.Value over the
                        full gene family of each axis DE table and matches the
                        stored adj.P.Val exactly (tolerance 1e-12 relative).
3. counts_match_audit   FDR<0.05 gene counts equal the audited v2 numbers
                        (GSE103927 3,183; GSE333506 2,103).
4. direction_consistent direction column agrees with sign(logFC), and the
                        significant_fdr_0.05 flag agrees with FDR<0.05.

Read-only with respect to all frozen directories; the report is written to
work/02_altitude_response/validate_module2_v2_bh_<date>/report.json.

Usage:
    python3 scripts/02_altitude_response/validate_module2_v2_bh.py
        [--project-root DIR] [--freeze-dir DIR] [--report-dir DIR]
"""

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

EXPECTED_FDR_COUNTS = {
    "GSE103927": 3183,
    "GSE333506": 2103,
}
AXIS_REL = {
    "GSE103927": "results/02_altitude_response/analyses/m2_gse103927_alt16_paired_v2_20260910/GSE103927/differential_expression.tsv",
    "GSE333506": "results/02_altitude_response/analyses/m2_gse333506_han_week4_paired_v2_20260910/GSE333506/differential_expression.tsv",
}
FREEZE_REL = "data/processed/02_altitude_response/releases/m2_module3_input_freeze_v2_20260910"
REL_TOL = 1e-12


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bh_from_pvals(pvals):
    """BH adjusted p-values (limma/step-up convention), computed independently."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        val = min(prev, pvals[i] * n / rank)
        q[i] = val
        prev = val
    return q


def check_manifest(root: Path, freeze: Path, checks, failures):
    manifest_path = freeze / "source_manifest.tsv"
    if not manifest_path.is_file():
        failures.append("source_manifest.tsv missing")
        checks["manifest_closure"] = False
        return
    problems = []
    for row in read_tsv(manifest_path):
        p = Path(row["path"])
        if not p.is_absolute():
            p = root / p
        if not p.is_file():
            problems.append(f"missing {row['path']}")
            continue
        if p.stat().st_size != int(row["bytes"]):
            problems.append(f"size mismatch {row['path']}")
        if sha256_of(p) != row["sha256"]:
            problems.append(f"sha256 mismatch {row['path']}")
    checks["manifest_closure"] = not problems
    if problems:
        failures.extend(problems)


def check_axis(root: Path, dataset: str, checks, failures):
    path = root / AXIS_REL[dataset]
    if not path.is_file():
        failures.append(f"{dataset}: DE table missing: {path}")
        return
    rows = read_tsv(path)
    pvals = [float(r["P.Value"]) for r in rows]
    stored = [float(r["adj.P.Val"]) for r in rows]
    recomputed = bh_from_pvals(pvals)

    problems = []
    max_rel = 0.0
    for i, (s, r) in enumerate(zip(stored, recomputed)):
        denom = max(1.0, abs(s))
        rel = abs(s - r) / denom
        max_rel = max(max_rel, rel)
        if rel > REL_TOL:
            problems.append(f"{dataset}: BH mismatch row {i} (stored {s}, recomputed {r})")
            if len(problems) > 5:
                break
    checks[f"bh_recomputed_{dataset}"] = not problems
    if problems:
        failures.extend(problems)

    n_sig = sum(1 for r in rows if float(r["adj.P.Val"]) < 0.05)
    ok = n_sig == EXPECTED_FDR_COUNTS[dataset]
    checks[f"counts_match_audit_{dataset}"] = ok
    if not ok:
        failures.append(f"{dataset}: FDR<0.05 count {n_sig} != audited {EXPECTED_FDR_COUNTS[dataset]}")

    dir_problems = []
    for i, r in enumerate(rows):
        expect_up = float(r["logFC"]) > 0
        if (r["direction"] == "up") != expect_up:
            dir_problems.append(f"{dataset}: direction/logFC sign disagree row {i}")
            break
        flag = str(r["significant_fdr_0.05"]).strip().lower() in {"true", "1", "yes"}
        if flag != (float(r["adj.P.Val"]) < 0.05):
            dir_problems.append(f"{dataset}: significant_fdr_0.05 flag disagrees with adj.P.Val row {i}")
            break
    checks[f"direction_consistent_{dataset}"] = not dir_problems
    failures.extend(dir_problems)
    print(f"  {dataset}: {len(rows)} genes, BH max rel diff {max_rel:.2e}, FDR<0.05 = {n_sig}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    ap.add_argument("--project-root", default=str(root))
    ap.add_argument("--freeze-dir", default=None)
    ap.add_argument("--report-dir", default=None)
    args = ap.parse_args()

    root = Path(args.project_root).resolve()
    freeze = Path(args.freeze_dir).resolve() if args.freeze_dir else root / FREEZE_REL

    checks = {}
    failures = []
    print(f"Validating module2 v2 deliverables against freeze {freeze}")
    check_manifest(root, freeze, checks, failures)
    for dataset in AXIS_REL:
        check_axis(root, dataset, checks, failures)

    report = {
        "validator": "validate_module2_v2_bh.py",
        "purpose": "audit B-10: independent BH recompute + freeze closure for module2 v2 axes",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "freeze_dir": str(freeze),
        "checks": checks,
        "checks_passed": sum(checks.values()),
        "checks_total": len(checks),
        "failures": failures,
        "status": "PASS" if not failures else "FAIL",
    }
    report_dir = Path(args.report_dir).resolve() if args.report_dir else (
        root / "work/02_altitude_response/validate_module2_v2_bh_20260910"
    )
    report_dir.mkdir(parents=True, exist_ok=True)
    out = report_dir / "report.json"
    out.write_text(json.dumps(report, indent=1) + "\n")
    print(f"{report['status']}: {report['checks_passed']}/{report['checks_total']} checks; report: {out}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
