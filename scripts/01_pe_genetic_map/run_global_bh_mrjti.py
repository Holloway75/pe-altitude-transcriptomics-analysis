#!/usr/bin/env python3
"""Run MR-JTI for every outcome-specific S-PrediXcan global-BH hit.

The execution unit is GWAS x gene x tissue.  Full cis-eQTL data are cached by
gene x tissue and the fixed hg19 EUR LD background is cached by gene.  GWAS
tables are scanned once per outcome after all required eQTL rsIDs are known.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MINIMAL_PATH = ROOT / "scripts/01_pe_genetic_map/mrjti_minimal.py"
SPEC = importlib.util.spec_from_file_location("mrjti_minimal_batch_helpers", MINIMAL_PATH)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)
PREPARE_PATH = ROOT / "scripts/01_pe_genetic_map/prepare_gwas.py"
PREPARE_SPEC = importlib.util.spec_from_file_location("prepare_gwas_batch_helpers", PREPARE_PATH)
P = importlib.util.module_from_spec(PREPARE_SPEC)
assert PREPARE_SPEC.loader is not None
PREPARE_SPEC.loader.exec_module(P)

FIELDS = [
    "dataset_id", "tissue", "gene_id", "gene_name", "source_z", "source_pvalue",
    "source_q_global", "source_family_size", "selection_rule", "selection_relation",
    "status", "n_snps", "standardized_expression_coefficient",
    "original_lasso_coefficient", "beta_scale", "ci_low", "ci_high", "ci_method",
    "lambda", "cv_mean_error", "cv_error_se", "bootstrap_zero_fraction",
    "bootstrap_positive_fraction", "bootstrap_negative_fraction", "pvalue",
    "q_mrjti", "mrjti_supported", "ci_significance", "mrjti_family_size",
    "mrjti_family_alpha", "direction_concordant", "inference_status", "task_dir",
]


def parse_source(value: str):
    if "=" not in value:
        raise argparse.ArgumentTypeError("source must be DATASET_ID=/absolute/path/to/spredixcan_all_tissues.tsv")
    dataset, raw_path = value.split("=", 1)
    path = Path(raw_path).resolve()
    if not dataset or not path.is_file():
        raise argparse.ArgumentTypeError(f"invalid source: {value}")
    return dataset, path


def file_sha(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def clean_id(value: str):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def freeze_candidates(sources, output: Path):
    candidates = []
    evidence = []
    for dataset_id, table in sources:
        parent = table.parent
        method_path = parent / "method.json"
        validation_path = parent / "validation-report.json"
        if not method_path.is_file() or not validation_path.is_file():
            raise RuntimeError(f"{dataset_id}: missing S-PrediXcan method or validation report")
        method = json.loads(method_path.read_text())
        validation = json.loads(validation_path.read_text())
        if validation.get("status") != "pass" or method.get("dataset_id") != dataset_id:
            raise RuntimeError(f"{dataset_id}: upstream validation/method identity failed")
        rows = list(M.read_tsv(table))
        family = []
        seen = set()
        for row in rows:
            key = (row.get("gene_id", "").split(".")[0], row.get("tissue", ""))
            if row.get("run_status") != "success":
                continue
            if key in seen:
                raise RuntimeError(f"{dataset_id}: duplicate successful gene-tissue row: {key}")
            seen.add(key)
            try:
                p = float(row["pvalue"])
                q = float(row["q_global"])
            except (KeyError, TypeError, ValueError):
                raise RuntimeError(f"{dataset_id}: successful row lacks finite P/q: {key}")
            if not (math.isfinite(p) and 0 <= p <= 1 and math.isfinite(q) and 0 <= q <= 1):
                raise RuntimeError(f"{dataset_id}: successful row has invalid P/q: {key}")
            family.append(row)
            if q <= 0.05:
                candidates.append({
                    "dataset_id": dataset_id, "tissue": row["tissue"], "gene_id": key[0],
                    "gene_name": row.get("gene_name") or "NA", "source_z": row["zscore"],
                    "source_pvalue": row["pvalue"], "source_q_global": row["q_global"],
                    "source_family_size": len(rows),
                })
        expected = int(method.get("n_associations", -1))
        if len(family) != expected or len(rows) != expected:
            raise RuntimeError(f"{dataset_id}: global-BH family closure failed: rows={len(rows)}, family={len(family)}, method={expected}")
        # Path relocation (2026-09-09): pre-migration method.json records absolute
        # paths under the former project root. Remap onto the current ROOT by the
        # data/processed/ suffix and enforce the recorded sha256 before use.
        gwas_recorded = method["normalized_gwas"]
        gwas_path = Path(gwas_recorded)
        gwas_remapped = False
        if not gwas_path.is_file():
            marker = "data/processed/"
            idx = gwas_recorded.find(marker)
            if idx == -1:
                raise RuntimeError(f"{dataset_id}: recorded normalized GWAS not found and not relocatable: {gwas_recorded}")
            gwas_path = ROOT / gwas_recorded[idx:]
            if not gwas_path.is_file():
                raise RuntimeError(f"{dataset_id}: recorded normalized GWAS missing after relocation: {gwas_path}")
            gwas_remapped = True
        if file_sha(gwas_path) != method["normalized_gwas_sha256"]:
            raise RuntimeError(f"{dataset_id}: normalized GWAS sha256 mismatch after relocation: {gwas_path}")
        evidence.append({
            "dataset_id": dataset_id, "table": str(table), "table_sha256": file_sha(table),
            "method": str(method_path.resolve()), "method_sha256": file_sha(method_path),
            "validation": str(validation_path.resolve()), "validation_sha256": file_sha(validation_path),
            "global_bh_family_size": len(family),
            "global_bh_significant_tasks": sum(float(x["q_global"]) <= 0.05 for x in family),
            "normalized_gwas": str(gwas_path),
            "normalized_gwas_sha256": method["normalized_gwas_sha256"],
            "normalized_gwas_recorded_path": gwas_recorded,
            "normalized_gwas_path_remapped": gwas_remapped,
        })
    candidates.sort(key=lambda x: (dict((d, i) for i, (d, _) in enumerate(sources))[x["dataset_id"]], x["tissue"], x["gene_id"]))
    for i, row in enumerate(candidates, 1):
        row["task_id"] = f"T{i:05d}"
    manifest = output / "candidate-manifest.tsv"
    M.write_tsv(manifest, ["task_id", "dataset_id", "tissue", "gene_id", "gene_name", "source_z", "source_pvalue", "source_q_global", "source_family_size"], candidates)
    return candidates, evidence


def read_epi(path: Path):
    out = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            chrom, gene, gd, bp, name, strand = line.rstrip("\n").split("\t")[:6]
            out[gene.split(".")[0]].append((chrom.removeprefix("chr"), int(bp), name, strand))
    return out


def resolve_loci(candidates, env):
    by_tissue = defaultdict(set)
    for row in candidates:
        by_tissue[row["tissue"]].add(row["gene_id"])
    loci_by_pair = {}
    unavailable = {}
    tissue_epi = {}
    for tissue, genes in sorted(by_tissue.items()):
        epi_path = Path(env["GTEX_V8_BESD_DIR"]) / tissue / f"{tissue}.epi"
        index = read_epi(epi_path)
        tissue_epi[tissue] = {"path": str(epi_path.resolve()), "sha256": file_sha(epi_path)}
        for gene in genes:
            hits = list(dict.fromkeys(index.get(gene, [])))
            if len(hits) != 1:
                unavailable[(gene, tissue)] = f"besd_probe_count_{len(hits)}"
                continue
            loci_by_pair[(gene, tissue)] = hits[0]
    by_gene = defaultdict(set)
    for (gene, tissue), (chrom, bp, name, strand) in loci_by_pair.items():
        by_gene[gene].add((chrom, bp))
    conflict = {g: x for g, x in by_gene.items() if len(x) != 1}
    for gene in conflict:
        for pair in list(loci_by_pair):
            if pair[0] == gene:
                unavailable[pair] = "besd_probe_locus_conflict_across_tissues"
                del loci_by_pair[pair]
    usable_by_gene = defaultdict(set)
    for (gene, tissue), (chrom, bp, name, strand) in loci_by_pair.items():
        usable_by_gene[gene].add((chrom, bp))
    return loci_by_pair, {g: next(iter(x)) for g, x in usable_by_gene.items()}, tissue_epi, unavailable


def extract_one_eqtl(row, env, output):
    resource = output / "resources/eqtl" / clean_id(row["tissue"]) / row["gene_id"]
    raw, kept, counts, path = M.extract_eqtl(env, resource, row["gene_id"], row["tissue"], row["gene_id"])
    if not raw:
        raise RuntimeError(f"{row['tissue']}/{row['gene_id']}: empty BESD extraction")
    return (row["gene_id"], row["tissue"]), {"raw": raw, "kept": kept, "counts": dict(counts), "path": path}


def extract_eqtls(candidates, env, output, workers):
    unique = {}
    for row in candidates:
        unique.setdefault((row["gene_id"], row["tissue"]), row)
    result = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(extract_one_eqtl, row, env, output) for row in unique.values()]
        for i, future in enumerate(as_completed(futures), 1):
            key, value = future.result()
            result[key] = value
            if i % 50 == 0 or i == len(futures):
                print(f"eQTL {i}/{len(futures)}", flush=True)
    return result


def load_gwas_subset(path: Path, wanted):
    found = defaultdict(list)
    for row in M.read_tsv(path):
        rid = row.get("rsid")
        if rid in wanted:
            found[rid].append(row)
    return found


def harmonize_from_map(eqtl, refrows, gwas_map):
    refs = {r["rsid"]: r for r in refrows}
    counts = Counter(eqtl_eligible=len(eqtl), gwas_rsid_overlap=sum(rid in gwas_map for rid in eqtl))
    rows = []
    for rid, e in eqtl.items():
        if rid not in gwas_map:
            counts["excluded_missing_gwas"] += 1; continue
        if len(gwas_map[rid]) != 1:
            counts["excluded_duplicate_gwas_rsid"] += 1; continue
        if rid not in refs:
            counts["excluded_missing_ld_reference"] += 1; continue
        g, rr = gwas_map[rid][0], refs[rid]
        if str(g["chromosome"]).removeprefix("chr") != e["Chr"] or int(float(g["position"])) != int(e["BP"]) or rr["position"] != int(e["BP"]):
            counts["excluded_position_conflict"] += 1; continue
        a1, a2 = e["A1"].upper(), e["A2"].upper()
        if M.palindromic(a1, a2):
            counts["excluded_palindromic_snp"] += 1; continue
        refsign, refstatus = M.pair_match(a1, a2, rr["ref"], rr["alt"], True)
        if not refsign:
            counts["excluded_reference_allele_conflict"] += 1; continue
        if not M.finite(g["beta"], g["standard_error"], g["pvalue"]) or float(g["standard_error"]) <= 0:
            counts["excluded_invalid_gwas_stat"] += 1; continue
        sign, status = M.pair_match(a1, a2, g["effect_allele"], g["non_effect_allele"], True)
        if not sign:
            counts["excluded_gwas_allele_conflict"] += 1; continue
        rows.append({
            "rsid": rid, "internal_id": rr["internal_id"], "genotype_id": rr["genotype_id"],
            "effect_allele": a1, "other_allele": a2, "chromosome": e["Chr"], "position": e["BP"],
            "variant_type": rr["variant_type"], "harmonization": status, "reference_match": refstatus,
            "ldscore": rr["ldscore"], "eqtl_beta": float(e["b"]), "eqtl_se": float(e["SE"]),
            "eqtl_p": float(e["p"]), "gwas_beta": sign * float(g["beta"]),
            "gwas_se": float(g["standard_error"]), "gwas_p": float(g["pvalue"]),
            "gwas_source_effect_allele": g["effect_allele"], "gwas_source_other_allele": g["non_effect_allele"],
            "gwas_source_beta": float(g["beta"]), "gwas_source_se": float(g["standard_error"]),
            "gwas_source_p": float(g["pvalue"]), "gwas_alignment_sign": sign,
        })
    rows.sort(key=lambda r: (int(r["chromosome"]), int(r["position"]), r["internal_id"]))
    counts["harmonized"] = len(rows)
    return rows, counts


def build_gene_reference(gene, locus, env, output):
    chrom, bp = locus
    resource = output / "resources/ld" / gene
    map_path = resource / "ld_reference/ld-reference-variant-map.tsv"
    method_path = resource / "resource-method.json"
    prefix = resource / "ld_reference" / f"{gene}_EUR_hg19"
    score = resource / "ld_reference" / f"{gene}_EUR_hg19_ldscore.score.ld"
    if map_path.is_file() and method_path.is_file() and prefix.with_suffix(".bed").is_file() and score.is_file():
        rows = list(M.read_tsv(map_path))
        for row in rows:
            row["position"] = int(row["position"])
            for field in ("maf", "missing_rate", "ldscore"):
                row[field] = float(row[field])
        return prefix, rows, score, json.loads(method_path.read_text())
    pvar = Path(env["LD_REF_HG19_DIR"]) / f"eur_chr{chrom}.pvar"
    candidates, counts = M.reference_candidates(pvar, int(chrom), max(1, bp - 2_000_000), bp + 2_000_000)
    fasta = P.FastaReference(env["HG19_FASTA"])
    checked = []
    try:
        for row in candidates:
            if fasta.fetch(str(chrom), int(row["position"]), len(row["ref"])).upper() != row["ref"]:
                counts["reference_excluded_fasta_ref_mismatch"] += 1
                continue
            if M.palindromic(row["ref"], row["alt"]):
                counts["reference_excluded_palindromic_snp"] += 1
                continue
            pos, ref, alt = P.left_align(int(row["position"]), row["ref"], row["alt"],
                                          lambda q, n: fasta.fetch(str(chrom), q, n))
            if (pos, ref, alt) != (int(row["position"]), row["ref"], row["alt"]):
                counts["reference_excluded_not_normalized"] += 1
                continue
            checked.append(row)
    finally:
        fasta.close()
    counts["reference_pre_genotype_qc_after_fasta_and_protocol"] = len(checked)
    prefix, rows, score = M.build_reference(env, resource, checked, int(chrom), gene)
    method = {"gene_id": gene, "chromosome": chrom, "probe_bp": bp, "pvar": str(pvar.resolve()),
              "pvar_sha256": file_sha(pvar), "counts": dict(counts), "reference_post_genotype_qc": len(rows)}
    atomic_json(method_path, method)
    return prefix, rows, score, method


def prune_task(env, task_dir: Path, refprefix: Path, rows):
    task_dir.mkdir(parents=True, exist_ok=True)
    M.write_tsv(task_dir / "harmonized.tsv", list(rows[0]) if rows else ["rsid"], rows)
    M.write_tsv(task_dir / "target-genotype-ids.txt", ["genotype_id"], rows)
    pre = task_dir / "prune"
    M.run([env["PLINK2"], "--bfile", refprefix, "--extract", task_dir / "target-genotype-ids.txt",
           "--indep-pairwise", "1000kb", "1", "0.1", "--out", pre], task_dir / "plink-prune.log")
    kept = {x.strip() for x in Path(str(pre) + ".prune.in").read_text().splitlines() if x.strip()}
    pruned = [r for r in rows if r["genotype_id"] in kept]
    check = task_dir / "prune-validation"
    M.run([env["PLINK2"], "--bfile", refprefix, "--extract", Path(str(pre) + ".prune.in"),
           "--r2-unphased", "cols=id", "--ld-window-kb", "1000", "--ld-window-r2", "0.1000000001",
           "--out", check], task_dir / "plink-prune-validation.log")
    vcor = check.with_suffix(".vcor")
    if vcor.is_file() and sum(1 for _ in vcor.open()) > 1:
        raise RuntimeError("LD pruning threshold validation failed")
    fields = ["rsid", "internal_id", "genotype_id", "chromosome", "position", "variant_type",
              "effect_allele", "other_allele", "ldscore", "eqtl_beta", "eqtl_se", "eqtl_p",
              "gwas_beta", "gwas_se", "gwas_p", "gwas_source_effect_allele",
              "gwas_source_other_allele", "gwas_source_beta", "gwas_source_se", "gwas_source_p",
              "gwas_alignment_sign"]
    M.write_tsv(task_dir / "mrjti-input.tsv", fields, pruned)
    M.write_tsv(task_dir / "pruned-variants.tsv", list(pruned[0]) if pruned else fields, pruned)
    return pruned


def blank_result(row, task_dir):
    return {
        "dataset_id": row["dataset_id"], "tissue": row["tissue"], "gene_id": row["gene_id"],
        "gene_name": row["gene_name"], "source_z": row["source_z"], "source_pvalue": row["source_pvalue"],
        "source_q_global": row["source_q_global"], "source_family_size": row["source_family_size"],
        "selection_rule": "outcome-specific global BH q<=0.05 among all successful finite gene-tissue tests",
        "selection_relation": "same_outcome_selection_and_analysis", "status": "model_failed", "n_snps": 0,
        "standardized_expression_coefficient": M.NA, "original_lasso_coefficient": M.NA,
        "beta_scale": "standardized_across_retained_variants", "ci_low": M.NA, "ci_high": M.NA,
        "ci_method": M.NA, "lambda": M.NA, "cv_mean_error": M.NA, "cv_error_se": M.NA,
        "bootstrap_zero_fraction": M.NA, "bootstrap_positive_fraction": M.NA,
        "bootstrap_negative_fraction": M.NA, "pvalue": M.NA, "q_mrjti": M.NA,
        "mrjti_supported": M.NA, "ci_significance": M.NA, "mrjti_family_size": M.NA,
        "mrjti_family_alpha": M.NA, "direction_concordant": M.NA, "inference_status": "not_run",
        "task_dir": str(task_dir.resolve()),
    }


def unavailable_result(row, output, reason):
    task_dir = output / "tasks" / row["dataset_id"] / clean_id(row["tissue"]) / row["gene_id"]
    task_dir.mkdir(parents=True, exist_ok=True)
    result = blank_result(row, task_dir)
    result["status"] = "harmonization_failed"
    result["inference_status"] = f"not_run:{reason}"
    atomic_json(task_dir / "task-result.json", result)
    return result


def run_task(row, eqtl_info, ref, gwas_map, env, output, family_size):
    task_dir = output / "tasks" / row["dataset_id"] / clean_id(row["tissue"]) / row["gene_id"]
    done = task_dir / "task-result.json"
    if done.is_file():
        return json.loads(done.read_text())
    result = blank_result(row, task_dir)
    try:
        refprefix, refrows, score, refmethod = ref
        harmonized, hc = harmonize_from_map(eqtl_info["kept"], refrows, gwas_map)
        pruned = prune_task(env, task_dir, refprefix, harmonized)
        result["n_snps"] = len(pruned)
        flow = dict(eqtl_info["counts"])
        flow.update(hc)
        flow.update({"reference_post_genotype_qc": len(refrows), "pruned": len(pruned)})
        M.write_tsv(task_dir / "snp-flow.tsv", sorted(flow), [flow])
        if len(pruned) < 20:
            result["status"] = "ineligible_lt20"
            result["inference_status"] = "not_run_below_execution_threshold"
        else:
            M.run([env["RSCRIPT"], ROOT / "scripts/01_pe_genetic_map/mrjti_minimal_runner.R",
                   "--df-path", task_dir / "mrjti-input.tsv", "--vendor-script", ROOT / "tools/MR-JTI/mr/MR-JTI.r",
                   "--result-path", task_dir / "mrjti-result.tsv", "--bootstrap-path", task_dir / "bootstrap.tsv",
                   "--n-folds", "5", "--n-bootstrap", "500", "--n-genes", str(family_size),
                   "--min-snps", "20", "--seed", str(M.SEED)],
                  task_dir / "mrjti.log")
            rr = next(M.read_tsv(task_dir / "mrjti-result.tsv"))
            beta = float(rr["standardized_expression_coefficient"])
            result.update({
                "status": "success_inference_bonferroni_selection_dependent", "standardized_expression_coefficient": beta,
                "original_lasso_coefficient": rr["original_lasso_coefficient"], "ci_low": rr["ci_low"],
                "ci_high": rr["ci_high"], "ci_method": rr["ci_method"], "lambda": rr["lambda"],
                "cv_mean_error": rr["cv_mean_error"], "cv_error_se": rr["cv_error_se"],
                "bootstrap_zero_fraction": rr["bootstrap_zero_fraction"],
                "bootstrap_positive_fraction": rr["bootstrap_positive_fraction"],
                "bootstrap_negative_fraction": rr["bootstrap_negative_fraction"],
                "ci_significance": rr["ci_significance"], "mrjti_supported": str(rr["ci_significance"] == "sig").lower(),
                "mrjti_family_size": rr["n_genes"], "mrjti_family_alpha": rr["family_alpha"],
                "direction_concordant": str(beta * float(row["source_z"]) > 0).lower(),
                "inference_status": "native_bonferroni_ci_selection_dependent",
            })
    except Exception as exc:
        result["status"] = "model_failed"
        result["inference_status"] = f"failed:{type(exc).__name__}:{exc}"
    atomic_json(done, result)
    return result


def validate(output, candidates, results):
    errors = []
    expected = {(x["dataset_id"], x["gene_id"], x["tissue"]) for x in candidates}
    observed = [(x["dataset_id"], x["gene_id"], x["tissue"]) for x in results]
    if len(observed) != len(set(observed)) or set(observed) != expected:
        errors.append("candidate/result task closure failed")
    counts = Counter(x["status"] for x in results)
    for row in results:
        if any(str(row[k]) != M.NA for k in ("pvalue", "q_mrjti")):
            errors.append(f"{row['dataset_id']}/{row['tissue']}/{row['gene_id']}: unsupported P/q field populated")
        if row["status"] == "success_inference_bonferroni_selection_dependent":
            task = Path(row["task_dir"])
            boots = list(M.read_tsv(task / "bootstrap.tsv")) if (task / "bootstrap.tsv").is_file() else []
            if int(row["n_snps"]) < 20 or len(boots) != 500:
                errors.append(f"{row['dataset_id']}/{row['tissue']}/{row['gene_id']}: success artifact closure failed")
        elif row["status"] == "ineligible_lt20" and int(row["n_snps"]) >= 20:
            errors.append(f"{row['dataset_id']}/{row['tissue']}/{row['gene_id']}: ineligible count mismatch")
        elif row["status"] == "harmonization_failed" and not str(row["inference_status"]).startswith("not_run:besd_probe_"):
            errors.append(f"{row['dataset_id']}/{row['tissue']}/{row['gene_id']}: invalid harmonization-failure reason")
        elif row["status"] not in {"ineligible_lt20", "success_inference_bonferroni_selection_dependent", "harmonization_failed"}:
            errors.append(f"{row['dataset_id']}/{row['tissue']}/{row['gene_id']}: task failed ({row['inference_status']})")
    report = {"status": "pass" if not errors else "fail", "errors": errors,
              "candidate_tasks": len(candidates), "result_tasks": len(results), "status_counts": dict(counts),
              "native_bonferroni_ci_available": True, "calibrated_p_or_q_available": False, "repeatability_proven": False}
    atomic_json(output / "validation-report.json", report)
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", action="append", type=parse_source, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-tasks", type=int, help="Interface test only; output is labelled test_subset")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()) and not args.resume:
        raise SystemExit("output exists and is nonempty; use a new directory or --resume")
    output.mkdir(parents=True, exist_ok=True)
    env = M.load_env()
    candidates, evidence = freeze_candidates(args.source, output)
    full_count = len(candidates)
    if args.max_tasks is not None:
        candidates = candidates[:args.max_tasks]
        M.write_tsv(output / "candidate-manifest.tsv", ["task_id", "dataset_id", "tissue", "gene_id", "gene_name", "source_z", "source_pvalue", "source_q_global", "source_family_size"], candidates)
    loci_by_pair, loci_by_gene, epi_evidence, unavailable = resolve_loci(candidates, env)
    runnable = [x for x in candidates if (x["gene_id"], x["tissue"]) not in unavailable]
    eqtls = extract_eqtls(runnable, env, output, max(1, args.workers))
    by_dataset = defaultdict(list)
    for row in candidates:
        by_dataset[row["dataset_id"]].append(row)
    evidence_by_dataset = {x["dataset_id"]: x for x in evidence}
    references = {}
    genes = sorted({x["gene_id"] for x in runnable})
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(build_gene_reference, gene, loci_by_gene[gene], env, output): gene for gene in genes}
        for i, future in enumerate(as_completed(futures), 1):
            gene = futures[future]
            references[gene] = future.result()
            if i % 10 == 0 or i == len(futures):
                print(f"LD reference {i}/{len(futures)}", flush=True)
    results = []
    source_order = [x[0] for x in args.source]
    for dataset_id in source_order:
        tasks = by_dataset.get(dataset_id, [])
        if not tasks:
            continue
        runnable_tasks = [x for x in tasks if (x["gene_id"], x["tissue"]) not in unavailable]
        wanted = set()
        for row in runnable_tasks:
            wanted.update(eqtls[(row["gene_id"], row["tissue"])]["kept"])
        gwas_path = Path(evidence_by_dataset[dataset_id]["normalized_gwas"])
        if file_sha(gwas_path) != evidence_by_dataset[dataset_id]["normalized_gwas_sha256"]:
            raise RuntimeError(f"{dataset_id}: normalized GWAS hash changed")
        print(f"scanning GWAS {dataset_id} for {len(wanted)} eQTL rsIDs", flush=True)
        gwas_map = load_gwas_subset(gwas_path, wanted)
        dataset_results = [None] * len(tasks)
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futures = {}
            for index, row in enumerate(tasks):
                gene = row["gene_id"]
                pair = (gene, row["tissue"])
                if pair in unavailable:
                    future = pool.submit(unavailable_result, row, output, unavailable[pair])
                else:
                    future = pool.submit(run_task, row, eqtls[pair], references[gene], gwas_map, env, output, len(tasks))
                futures[future] = index
            for i, future in enumerate(as_completed(futures), 1):
                result = future.result()
                dataset_results[futures[future]] = result
                if i % 10 == 0 or i == len(futures):
                    print(f"{dataset_id} task {i}/{len(tasks)} status={result['status']}", flush=True)
        results.extend(dataset_results)
        del gwas_map
    M.write_tsv(output / "mrjti-global-bh.tsv", FIELDS, results)
    protocol = ROOT / "docs/protocol/SNP_MATCHING_RULES.md"
    method = {
        "schema_version": "m1.6-mrjti-native-bonferroni-v2", "scope": "test_subset" if args.max_tasks is not None else "all_global_bh_hits",
        "full_candidate_task_count": full_count, "executed_candidate_task_count": len(candidates),
        "unique_gene_tissue_count": len({(x['gene_id'], x['tissue']) for x in candidates}),
        "unique_gene_count": len({x['gene_id'] for x in candidates}), "source_evidence": evidence,
        "unavailable_besd_probe_tasks": len([x for x in candidates if (x['gene_id'], x['tissue']) in unavailable]),
        "unavailable_besd_probe_pairs": [{"gene_id": g, "tissue": t, "reason": reason} for (g, t), reason in sorted(unavailable.items())],
        "besd_epi_evidence": epi_evidence, "candidate_manifest": str((output / 'candidate-manifest.tsv').resolve()),
        "candidate_manifest_sha256": file_sha(output / "candidate-manifest.tsv"),
        "selection_rule": "within each GWAS, q_global<=0.05 over all unique successful finite gene-tissue S-PrediXcan tests",
        "selection_dependency": "same outcome used for candidate selection and MR-JTI; Bonferroni significance is conditional on the selected family and is not independent causal confirmation",
        "matching_protocol": {"path": str(protocol.resolve()), "declared_version": "v1.0", "sha256": file_sha(protocol),
            "target_build": "GRCh37/hg19", "reference_fasta": env["HG19_FASTA"], "coordinate_conversion": "none",
            "indel_normalization": True, "multiallelic_mode": "remove_all", "palindromic_policy": "remove",
            "strand_policy": "allow_nonpalindromic_complement", "effect_flip_policy": "beta=-beta, EAF=1-EAF, OR=1/OR",
            "rsid_policy": "required_by_downstream"},
        "variant_universes": {"eqtl": "all finite GTEx BESD cis-eQTL; no significance filter",
            "ld_background": "1000G EUR hg19 TSS +/-2Mb, MAF>=0.01, missingness<=0.02",
            "mr_target": "outcome/tissue eQTL-GWAS-reference intersection after 1000kb r2=0.1 pruning"},
        "ld_score": {"software": "GCTA", "ld_wind_kb": 1000, "r2_cutoff": 0.01, "computed_before_pruning": True},
        "mrjti": {"folds": 5, "bootstrap": 500, "seed": M.SEED, "minimum_snps": 20,
            "family_definition": "planned selected gene-by-tissue tasks within each GWAS, including unavailable tasks",
            "family_sizes": {dataset: len(rows) for dataset, rows in by_dataset.items()},
            "ci": "vendored basic residual bootstrap with native Bonferroni alpha=0.05/n_genes",
            "pvalue": None, "qvalue": None, "support": "CI_significance sig/nonsig; selection-dependent"},
        "formal_publication_modified": False,
    }
    atomic_json(output / "method.json", method)
    report = validate(output, candidates, results)
    print(json.dumps(report, sort_keys=True))
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
