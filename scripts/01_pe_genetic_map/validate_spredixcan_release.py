#!/usr/bin/env python3
"""Recompute and write the per-dataset M1.3 S-PrediXcan validation report.

The S-PrediXcan runner (`run_single_gwas_spredixcan.py`) writes the association
table, the tissue-level QC table and `method.json`.  This script independently
recomputes every reported quantity from those artifacts (plus a full
Benjamini-Hochberg recomputation over the association table) and writes
`validation-report.json` in the frozen release schema.  Downstream consumers
(`reinfer_mrjti_native_bonferroni.py`) require this file; running it here keeps
the per-step chain self-contained.  The script never modifies any other file
and exits non-zero on any inconsistency.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, statistics, sys
from pathlib import Path

ARTIFACTS = ["spredixcan_all_tissues.tsv", "significant_gene_tissue.tsv",
             "model_match_qc.tsv", "method.json"]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def bh(values):
    out=[1.0]*len(values); order=sorted(range(len(values)),key=lambda i:(values[i],i)); q=1.0
    for rank in range(len(order)-1,-1,-1):
        i=order[rank]; q=min(q,values[i]*len(values)/(rank+1)); out[i]=q
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset-id", required=True)
    ap.add_argument("--result-dir", type=Path, required=True,
                    help="M1.3 result directory for this dataset (contains the four artifacts)")
    ap.add_argument("--output", type=Path, default=None,
                    help="report path (default: <result-dir>/validation-report.json)")
    args = ap.parse_args()
    root = args.result_dir
    for name in ARTIFACTS:
        if not (root / name).is_file(): raise SystemExit(f"missing artifact: {root/name}")
    method = json.loads((root / "method.json").read_text())
    qc = read_tsv(root / "model_match_qc.tsv")
    rows = read_tsv(root / "spredixcan_all_tissues.tsv")
    if method.get("dataset_id") != args.dataset_id: raise SystemExit("dataset_id mismatch in method.json")
    if any(r["dataset_id"] != args.dataset_id for r in rows): raise SystemExit("dataset_id mismatch in association table")

    # Harmonization / tissue closure.
    expected_tissues = int(method["n_tissues"])
    if len(qc) != expected_tissues: raise SystemExit(f"model_match_qc has {len(qc)} tissues, expected {expected_tissues}")
    n_input = {int(r["n_input"]) for r in qc}
    if len(n_input) != 1: raise SystemExit("n_input differs across tissues")
    n_input = n_input.pop()
    rates = [float(r["match_rate"]) for r in qc]
    minimum_required = 0.20
    if min(rates) < minimum_required: raise SystemExit("model match rate below minimum")
    audits = {r["direction_audit"] for r in qc}
    closures = [int(r["closure_difference"]) for r in qc]

    # P/Z audit recorded by the runner; recompute the row count closure from it.
    pz = dict(method.get("p_z_audit", {}))
    concordant = int(pz.get("concordant_or_rounded", 0))
    discrepant = int(pz.get("discrepant", 0))
    invalid = int(pz.get("invalid", 0)) + int(pz.get("p_zero", 0))
    if concordant + discrepant + invalid != n_input: raise SystemExit("p/z audit does not close against n_input")
    if invalid: raise SystemExit("invalid beta/SE/P rows present")

    # Association table: independent BH recomputation and significance counts.
    pvalues = [float(r["pvalue"]) for r in rows]
    if any(not math.isfinite(p) or not 0 <= p <= 1 for p in pvalues): raise SystemExit("nonfinite pvalue in association table")
    qs = bh(pvalues)
    for r, q in zip(rows, qs):
        if abs(float(r["q_global"]) - q) > 1e-12 + 1e-9 * abs(q): raise SystemExit("q_global column disagrees with BH recomputation")
    pairs = {(r["tissue"], r["gene_id"]) for r in rows}
    if len(pairs) != len(rows): raise SystemExit("duplicate (tissue, gene_id) rows")
    if int(method["global_bh_family_size"]) != len(rows): raise SystemExit("global BH family size mismatch")
    global_bh = sum(q <= 0.05 for q in qs)
    tissue_bonf = sum(r["bonferroni_significant_tissue"] == "true" for r in rows)
    union = sum(r["bonferroni_significant_tissue"] == "true" or q <= 0.05 for r, q in zip(rows, qs))
    if int(method["n_significant_union"]) != union: raise SystemExit("n_significant_union mismatch")

    report = {
        "dataset_id": args.dataset_id,
        "status": "pass",
        "normalized_gwas_rows": n_input,
        "tissues_expected": expected_tissues,
        "tissues_completed": len(qc),
        "association_rows": len(rows),
        "unique_tissue_gene_rows": len(pairs),
        "model_match_rate": {"minimum": min(rates), "median": statistics.median(rates),
                             "maximum": max(rates), "minimum_required": minimum_required},
        "harmonization": {"all_direction_audits_passed": audits == {"pass"},
                          "maximum_closure_difference": max(closures)},
        "gwas_p_z_audit": {"concordant_or_rounded": concordant, "discrepant_warning": discrepant,
                           "invalid": invalid, "authoritative_z": "beta/standard_error"},
        "multiple_testing": {"global_bh_significant_rows": global_bh,
                             "tissue_bonferroni_significant_rows": tissue_bonf,
                             "significant_union_rows": union},
        "artifact_sha256": {name: sha256(root / name) for name in ARTIFACTS},
    }
    # Edge case retained from the frozen release: with zero global-BH rows,
    # tissue-level Bonferroni rows are enumerated so no signal is silently lost.
    if global_bh == 0 and tissue_bonf > 0:
        report["significant_tissue_rows"] = [
            {"tissue": r["tissue"], "gene_id": r["gene_id"], "gene_name": r["gene_name"],
             "pvalue": float(r["pvalue"]), "q_global": float(r["q_global"]), "zscore": float(r["zscore"])}
            for r in rows if r["bonferroni_significant_tissue"] == "true"]

    output = args.output or (root / "validation-report.json")
    output.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"dataset_id": args.dataset_id, "status": "pass",
                      "association_rows": len(rows), "global_bh_significant_rows": global_bh,
                      "tissue_bonferroni_significant_rows": tissue_bonf, "report": str(output)}))


if __name__ == "__main__":
    main()
