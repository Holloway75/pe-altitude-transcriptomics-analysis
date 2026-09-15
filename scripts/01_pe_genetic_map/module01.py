#!/usr/bin/env python3
"""Module 1 scientific validation and publication utilities.

Scientific provenance is recorded in task-level dependencies and an append-only
ledger.  This module validates formal outputs and atomically publishes a
validated analysis workspace; it does not gate execution through manifests.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Iterable
sys.path.insert(0,str(Path(__file__).resolve().parent))
import stability_summary as stability_core

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config"
PLAN = ROOT / "docs/protocol/MODULE_01_PE_GENETIC_MAP_ANALYSIS_PLAN.md"

TABLES: dict[str, tuple[list[str], list[str]]] = {
    "pe_tissue_gene_full.tsv": (
        ["schema_version", "analysis_id", "source_provenance_sha256", "dataset_id", "analysis_role", "tissue", "gene_id", "gene_name", "zscore", "effect_size", "pvalue", "q_gene_tissue", "p_bonf_tissue", "p_bonf_global", "bonferroni_significant", "n_snps_in_model", "n_snps_used", "pred_perf_r2", "pred_perf_pvalue", "direction_qc", "run_status", "task_id", "attempt", "input_sha256", "method_version"],
        ["analysis_id", "dataset_id", "tissue", "gene_id"],
    ),
    "pe_whole_blood.tsv": (
        ["schema_version", "analysis_id", "source_provenance_sha256", "dataset_id", "analysis_role", "gene_id", "gene_name", "tissue", "zscore", "direction", "pvalue", "q_gene_tissue", "n_snps_in_model", "n_snps_used", "pred_perf_r2", "direction_qc", "run_status", "task_id", "attempt", "input_sha256"],
        ["analysis_id", "gene_id"],
    ),
    "pe_cross_tissue.tsv": (
        ["schema_version", "analysis_id", "source_provenance_sha256", "dataset_id", "analysis_role", "gene_id", "gene_name", "method", "statistic", "signed_statistic", "direction_available", "direction", "pvalue", "qvalue", "n_tissues", "n_independent_components", "best_tissue", "best_tissue_z", "condition_number", "model_status", "task_id", "attempt", "input_sha256", "method_version"],
        ["analysis_id", "dataset_id", "gene_id"],
    ),
    "pe_candidate_locus_map.tsv": (
        ["schema_version", "analysis_id", "source_provenance_sha256", "dataset_id", "statistical_layer", "candidate_id", "locus_definition", "locus_id", "candidate_rule_id", "analysis_role", "gene_id", "trigger_tissue", "cross_tissue_method", "candidate_direction", "candidate_statistic", "candidate_p", "candidate_q", "anchor_evidence", "anchor_source", "anchor_count", "resolution_status", "complex_ld_region", "ld_resource_sha256"],
        ["analysis_id", "dataset_id", "statistical_layer", "candidate_id", "locus_definition", "locus_id"],
    ),
    "pe_locus_map.tsv": (
        ["schema_version", "analysis_id", "source_provenance_sha256", "dataset_id", "statistical_layer", "locus_definition", "locus_id", "analysis_role", "chromosome", "locus_start", "locus_end", "locus_direction", "n_candidates", "n_genes", "n_direction_positive", "n_direction_negative", "n_direction_missing", "complex_ld_region", "resolution_status", "ld_resource_sha256"],
        ["analysis_id", "dataset_id", "statistical_layer", "locus_definition", "locus_id"],
    ),
    "pe_outcome_stability.tsv": (
        ["schema_version", "analysis_id", "source_provenance_sha256", "comparison_level", "analysis_set", "feature_id", "primary_dataset", "comparison_dataset", "rank_metric", "comparison_role", "primary_stat", "comparison_stat", "primary_direction", "comparison_direction", "direction_concordant", "primary_rank_pct", "comparison_rank_pct", "primary_p", "comparison_p", "primary_q", "comparison_q", "primary_status", "comparison_status", "n_common", "sample_overlap", "power_interpretation"],
        ["analysis_id", "comparison_level", "analysis_set", "feature_id", "primary_dataset", "comparison_dataset", "rank_metric"],
    ),
}
STATUS={"not_started","running","success","not_applicable","smoke_test_failed","not_implemented","qc_excluded","technical_failed","blocked_upstream"}
ROLES={"primary","phenotype_sensitivity","broad_sensitivity","low_power_exploratory"}
BOOL={"true","false"}
P_FIELDS={"pvalue","q_gene_tissue","qvalue","candidate_p","candidate_q","primary_p","comparison_p","primary_q","comparison_q","p_bonf_tissue","p_bonf_global","pred_perf_pvalue"}
INT_FIELDS={"attempt","n_snps_in_model","n_snps_used","n_tissues","n_independent_components","anchor_count","locus_start","locus_end","n_candidates","n_genes","n_direction_positive","n_direction_negative","n_direction_missing","n_common"}
ENUMS={
 "run_status":{"success","qc_excluded","technical_failed","blocked_upstream"}, "direction_qc":{"pass","fail","not_available"},
 "model_status":{"success","not_applicable_to_capability","technical_failed","qc_excluded","blocked_upstream"},
 "direction":{"risk_increasing","risk_decreasing","NA"}, "candidate_direction":{"risk_increasing","risk_decreasing","NA"},
 "resolution_status":{"resolved","unresolved"}, "statistical_layer":{"tissue","cross_tissue"},
 "comparison_level":{"whole_blood_gene","whole_blood_candidate_projection","tissue_candidate_projection","locus"},
 "analysis_set":{"all","exclude_complex_ld"}, "rank_metric":{"signed_z","abs_z","candidate_projection","locus_direction"},
 "sample_overlap":{"known_or_possible","none","unknown"},
}

def row_contract_errors(table,row):
    errors=[]
    for field,allowed in ENUMS.items():
        if field in row and row[field] not in allowed: errors.append(f"invalid enum {field}")
    if table=="pe_tissue_gene_full.tsv":
        success=row.get("run_status")=="success"
        required=("zscore","pvalue","q_gene_tissue","n_snps_in_model","n_snps_used","task_id","attempt","input_sha256","method_version")
        if success and any(row.get(x)=="NA" for x in required): errors.append("successful tissue row missing required value")
        if not success and any(row.get(x)!="NA" for x in ("zscore","effect_size","pvalue","q_gene_tissue")): errors.append("failed tissue row contains statistic")
    elif table=="pe_cross_tissue.tsv":
        success=row.get("model_status")=="success"
        if success and (row.get("method")!="S-MultiXcan" or any(row.get(x)=="NA" for x in ("statistic","pvalue","qvalue","n_tissues","n_independent_components","condition_number"))): errors.append("successful cross row incomplete")
        if success and (row.get("direction_available")!="false" or row.get("direction")!="NA" or row.get("signed_statistic")!="NA"): errors.append("unsigned cross row direction contract")
    elif table=="pe_candidate_locus_map.tsv":
        tissue=row.get("statistical_layer")=="tissue"; rule=row.get("candidate_rule_id")
        if rule != ("tissue_global_bh_q05_v1" if tissue else "cross_tissue_bh_q05_v1"): errors.append("layer/rule mismatch")
        if tissue and not (row.get("trigger_tissue")!="NA" and row.get("cross_tissue_method")=="NA"): errors.append("tissue candidate condition")
        if not tissue and not (row.get("trigger_tissue")=="NA" and row.get("cross_tissue_method")=="S-MultiXcan" and row.get("candidate_direction")=="NA"): errors.append("cross candidate condition")
    elif table=="pe_outcome_stability.tsv":
        level=row.get("comparison_level"); metric=row.get("rank_metric"); aset=row.get("analysis_set")
        allowed={"whole_blood_gene":{"signed_z","abs_z"},"whole_blood_candidate_projection":{"candidate_projection"},"tissue_candidate_projection":{"candidate_projection"},"locus":{"locus_direction"}}
        if metric not in allowed.get(level,set()): errors.append("comparison/rank combination")
        if level!="locus" and aset!="all": errors.append("analysis_set only applies to locus")
    return errors


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()


def load_env() -> dict[str, str]:
    command = "set -a; source scripts/lib/load_config.sh; env -0"
    proc = subprocess.run(["bash", "-c", command], cwd=ROOT, check=True, stdout=subprocess.PIPE)
    return dict(item.split("=", 1) for item in proc.stdout.decode(errors="surrogateescape").split("\0") if "=" in item)


def expand(value: str, env: dict[str, str]) -> Path:
    for key, replacement in env.items():
        value = value.replace("${" + key + "}", replacement)
    return Path(value)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener=gzip.open if str(path).endswith(".gz") else open
    with opener(path,"rt",newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def analysis_context(analysis_id: str) -> dict[str, Any]:
    """Direct scientific context for one stable analysis workspace.

    This deliberately replaces execution manifests.  Provenance is recorded
    from the registered inputs and scientific configuration, not used as a
    permission gate for unrelated code changes.
    """
    if not analysis_id or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for c in analysis_id):
        die("invalid analysis_id")
    env=load_env(); cfg=json.loads((CONFIG/"module01.json").read_text())
    datasets={r["dataset_id"]:r["analysis_role"] for r in read_tsv(CONFIG/"pe_gwas.tsv")}
    tissues={p.stem.removeprefix("JTI_") for p in Path(env["JTI_MODEL_DIR"]).glob("JTI_*.db")}
    ann=expand(cfg["resources"]["gene_annotation"],env)
    genes={r["gene_id"].split(".")[0] for r in read_tsv(ann)} if ann.is_file() else set()
    provenance=hashlib.sha256((CONFIG/"pe_gwas.tsv").read_bytes()+canonical(cfg)).hexdigest()
    return {"analysis_id":analysis_id,"env":env,"configuration":cfg,"datasets":datasets,"tissues":tissues,"genes":genes,"source_provenance_sha256":provenance}


def rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def candidate_digest(row: dict[str, str]) -> str:
    if row["statistical_layer"] == "tissue":
        fields = [row["candidate_rule_id"], row["dataset_id"], row["trigger_tissue"], row["gene_id"]]
    else:
        fields = [row["candidate_rule_id"], row["dataset_id"], row["gene_id"], row["cross_tissue_method"]]
    return hashlib.sha256("\0".join(fields).encode()).hexdigest()

def validate_ledger(ledger):
    # 2026-09-10 (audit A-2): the ledger assembled by assemble_m1_downstream.py is a
    # synthetic placeholder (rows derived from task_id, single hardcoded timestamp).
    # These checks therefore validate only the internal consistency and completeness
    # of the intended task enumeration; they are not evidence that tasks executed.
    terminal={}
    running=set()
    terminal_pairs=set()
    terminal_task_attempt={}
    for r in ledger:
        if r.get("status") not in STATUS: die("unknown ledger status")
        pair=(r["task_id"],r["attempt"])
        if r.get("status")=="running": running.add(pair)
        # A task may have immutable failed attempts followed by a successful
        # recovery attempt.  Publication evaluates only the newest terminal
        # attempt for each task; historical attempts remain in the ledger.
        if r.get("status")!="running":
            terminal_pairs.add(pair)
            key=r["task_id"]
            terminal_task_attempt[key] = max(int(r["attempt"]), terminal_task_attempt.get(key, -1))
            if key not in terminal or int(r["attempt"])>=int(terminal[key]["attempt"]): terminal[key]=r
    final=list(terminal.values())
    # Historical ledgers retain the immutable running row alongside its
    # terminal completion row.  Only an unmatched running attempt is active.
    if any(pair not in terminal_pairs and terminal_task_attempt.get(pair[0], -1) < int(pair[1]) for pair in running): die("running task in publication ledger")
    if any(r["requirement_class"]=="core_required" and r["status"]!="success" for r in final): die("core task did not succeed")
    if any(r["status"]=="technical_failed" for r in final): die("formal technical failure forbids publication")
    required_prefixes={"M1.1","M1.2","M1.3","M1.4","M1.5","M1.6","M1.7","M1.8","M1.9"}; present={r["task_id"].split("_")[0] for r in final}
    if not required_prefixes<=present: die("missing required task classes in ledger")
    return final
def validate_task_matrix(final,datasets,tissues):
    ids={r["task_id"]:r for r in final}
    missing=[]
    for d in datasets:
        if ids.get(f"M1.1_prepare_gwas__{d}",{}).get("status")!="success": missing.append(f"M1.1:{d}")
        for t in tissues:
            for stage in ("M1.2_harmonize","M1.3_spredixcan"):
                if ids.get(f"{stage}__{d}__{t}",{}).get("status")!="success": missing.append(f"{stage}:{d}:{t}")
    if missing: die("incomplete dataset-tissue task matrix: "+",".join(missing))
    return True


def validate_outputs(args: argparse.Namespace) -> dict[str, Any]:
    directory = Path(args.directory)
    context=analysis_context(args.analysis_id)
    run_id = context["analysis_id"]
    report: dict[str, Any] = {"analysis_id": run_id, "tables": {}, "errors": []}
    loaded: dict[str, list[dict[str, str]]] = {}
    source_sha=context["source_provenance_sha256"]
    datasets=context["datasets"]; tissues=context["tissues"]; genes=context["genes"]
    for name, (required, key) in TABLES.items():
        path = directory / name
        if not path.is_file():
            report["errors"].append(f"missing table {name}"); continue
        header, data = rows(path)
        loaded[name] = data
        if header != required:
            report["errors"].append(f"{name}: column contract mismatch")
        seen = set()
        for number, row in enumerate(data, 2):
            if row.get("analysis_id") != run_id or row.get("schema_version") != "M1.4":
                report["errors"].append(f"{name}:{number}: run/schema mismatch")
            if row.get("source_provenance_sha256")!=source_sha: report["errors"].append(f"{name}:{number}: source provenance mismatch")
            if "dataset_id" in row and (row["dataset_id"] not in datasets or row.get("analysis_role")!=datasets[row["dataset_id"]]): report["errors"].append(f"{name}:{number}: dataset/role foreign key")
            if "tissue" in row and row["tissue"] not in tissues: report["errors"].append(f"{name}:{number}: tissue foreign key")
            if "gene_id" in row and genes and row["gene_id"] not in genes and row["gene_id"]!="STATUS": report["errors"].append(f"{name}:{number}: gene annotation foreign key")
            k = tuple(row.get(x, "") for x in key)
            if k in seen: report["errors"].append(f"{name}:{number}: duplicate key")
            seen.add(k)
            if any(v in ("", "NaN", "Inf", "-Inf") for v in row.values()):
                report["errors"].append(f"{name}:{number}: invalid missing/nonfinite token")
            for field in P_FIELDS:
                if field in row and row[field] != "NA":
                    try: value = float(row[field])
                    except ValueError: report["errors"].append(f"{name}:{number}: {field} is not numeric"); continue
                    if not math.isfinite(value) or not 0 <= value <= 1:
                        report["errors"].append(f"{name}:{number}: {field} outside [0,1]")
            for field in INT_FIELDS & row.keys():
                if row[field]!="NA" and (not row[field].isdigit() or int(row[field])<0): report["errors"].append(f"{name}:{number}: invalid integer {field}")
            for field in ({"bonferroni_significant","direction_available","direction_concordant","complex_ld_region"}&row.keys()):
                if row[field]!="NA" and row[field] not in BOOL: report["errors"].append(f"{name}:{number}: invalid boolean {field}")
            for error in row_contract_errors(name,row): report["errors"].append(f"{name}:{number}: {error}")
        report["tables"][name] = {"rows": len(data), "size_bytes": path.stat().st_size, "sha256": sha256(path)}

    # Fixed order is part of the wire contract, not presentation metadata.
    role_order={r:i for i,r in enumerate(("primary","phenotype_sensitivity","broad_sensitivity","low_power_exploratory","exploratory"))}
    sorts={"pe_tissue_gene_full.tsv":lambda r:(role_order[r["analysis_role"]],r["tissue"],r["gene_id"]),"pe_whole_blood.tsv":lambda r:(-abs(float(r["zscore"])),r["gene_id"]),"pe_cross_tissue.tsv":lambda r:(role_order[r["analysis_role"]],float(r["pvalue"]) if r["pvalue"]!="NA" else math.inf,r["gene_id"]),"pe_candidate_locus_map.tsv":lambda r:(role_order[r["analysis_role"]],r["statistical_layer"],r["candidate_id"],r["locus_definition"],r["locus_id"]),"pe_locus_map.tsv":lambda r:(role_order[r["analysis_role"]],r["statistical_layer"],r["chromosome"],int(r["locus_start"]) if r["locus_start"]!="NA" else 10**20,r["locus_id"]),"pe_outcome_stability.tsv":lambda r:(r["comparison_level"],r["analysis_set"],role_order[r["comparison_role"]],r["feature_id"],r["rank_metric"])}
    for name,fn in sorts.items():
        data=loaded.get(name,[])
        if data!=sorted(data,key=fn): report["errors"].append(f"{name}: fixed sorting violated")

    for number, row in enumerate(loaded.get("pe_candidate_locus_map.tsv", []), 2):
        if row["candidate_id"] != candidate_digest(row):
            report["errors"].append(f"pe_candidate_locus_map.tsv:{number}: candidate_id mismatch")
        tissue = row["statistical_layer"] == "tissue"
        if tissue != (row["trigger_tissue"] != "NA" and row["cross_tissue_method"] == "NA"):
            report["errors"].append(f"pe_candidate_locus_map.tsv:{number}: layer fields inconsistent")

    detail = loaded.get("pe_candidate_locus_map.tsv", [])
    locus = loaded.get("pe_locus_map.tsv", [])
    detail_keys = {(r["analysis_id"], r["dataset_id"], r["statistical_layer"], r["locus_definition"], r["locus_id"]) for r in detail}
    locus_keys = {(r["analysis_id"], r["dataset_id"], r["statistical_layer"], r["locus_definition"], r["locus_id"]) for r in locus}
    if detail_keys != locus_keys:
        report["errors"].append("candidate/locus foreign-key closure failed")
    tissue_success = {(r["dataset_id"], r["tissue"], r["gene_id"]) for r in loaded.get("pe_tissue_gene_full.tsv", []) if r["run_status"] == "success"}
    cross_success = {(r["dataset_id"], r["gene_id"], r["method"]) for r in loaded.get("pe_cross_tissue.tsv", []) if r["model_status"] == "success"}
    tissue_index={(r["dataset_id"],r["tissue"],r["gene_id"]):r for r in loaded.get("pe_tissue_gene_full.tsv",[])}
    cross_index={(r["dataset_id"],r["gene_id"],r["method"]):r for r in loaded.get("pe_cross_tissue.tsv",[])}
    for number, row in enumerate(detail, 2):
        if row["statistical_layer"] == "tissue" and (row["dataset_id"], row["trigger_tissue"], row["gene_id"]) not in tissue_success:
            report["errors"].append(f"pe_candidate_locus_map.tsv:{number}: missing successful tissue statistic")
        if row["statistical_layer"] == "cross_tissue" and (row["dataset_id"], row["gene_id"], row["cross_tissue_method"]) not in cross_success:
            report["errors"].append(f"pe_candidate_locus_map.tsv:{number}: missing successful cross-tissue statistic")
        if row["candidate_rule_id"] not in {"tissue_global_bh_q05_v1","cross_tissue_bh_q05_v1"} or float(row["candidate_q"])>.05: report["errors"].append(f"pe_candidate_locus_map.tsv:{number}: candidate threshold/rule invalid")
        if row["statistical_layer"]=="tissue":
            src=tissue_index.get((row["dataset_id"],row["trigger_tissue"],row["gene_id"]))
            if not src or any(abs(float(row[a])-float(src[b]))>1e-12 for a,b in (("candidate_statistic","zscore"),("candidate_p","pvalue"),("candidate_q","q_gene_tissue"))): report["errors"].append(f"pe_candidate_locus_map.tsv:{number}: source statistic mismatch")
        else:
            src=cross_index.get((row["dataset_id"],row["gene_id"],row["cross_tissue_method"]))
            if not src or any(abs(float(row[a])-float(src[b]))>1e-12 for a,b in (("candidate_statistic","statistic"),("candidate_p","pvalue"),("candidate_q","qvalue"))): report["errors"].append(f"pe_candidate_locus_map.tsv:{number}: cross source statistic mismatch")
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in detail:
        k=(row["analysis_id"],row["dataset_id"],row["statistical_layer"],row["locus_definition"],row["locus_id"]); grouped.setdefault(k,[]).append(row)
    for number,row in enumerate(locus,2):
        k=(row["analysis_id"],row["dataset_id"],row["statistical_layer"],row["locus_definition"],row["locus_id"]); rs=grouped.get(k,[])
        pos=sum(x["candidate_direction"]=="risk_increasing" for x in rs); neg=sum(x["candidate_direction"]=="risk_decreasing" for x in rs); missing=len(rs)-pos-neg
        expected="undetermined" if missing else "conflicting" if pos and neg else "concordant_risk_increasing" if pos else "concordant_risk_decreasing"
        checks={"n_candidates":len({x["candidate_id"] for x in rs}),"n_genes":len({x["gene_id"] for x in rs}),"n_direction_positive":pos,"n_direction_negative":neg,"n_direction_missing":missing}
        if any(row[f] != str(v) for f,v in checks.items()) or row["locus_direction"] != expected:
            report["errors"].append(f"pe_locus_map.tsv:{number}: aggregate mismatch")
    locus_features={f'{r["statistical_layer"]}|{r["locus_definition"]}|{r["locus_id"]}' for r in locus}
    for number,row in enumerate(loaded.get("pe_outcome_stability.tsv",[]),2):
        if row["comparison_level"]=="locus" and row["feature_id"] not in locus_features: report["errors"].append(f"pe_outcome_stability.tsv:{number}: locus feature foreign key")
    full_primary={(r["gene_id"]):r for r in loaded.get("pe_tissue_gene_full.tsv",[]) if r["dataset_id"]=="FIGSHARE_22680904_v2" and r["tissue"]=="Whole_Blood" and r["run_status"]=="success" and r["direction_qc"]=="pass"}
    whole={r["gene_id"]:r for r in loaded.get("pe_whole_blood.tsv",[])}
    if set(full_primary)!=set(whole): report["errors"].append("Whole Blood/full bidirectional closure failed")
    elif any(any(x[a]!=whole[g][b] for a,b in (("zscore","zscore"),("pvalue","pvalue"),("q_gene_tissue","q_gene_tissue"))) for g,x in full_primary.items()): report["errors"].append("Whole Blood/full values differ")
    ancillary={"model-match-qc.tsv":{"dataset_id","tissue","n_model_snps","n_matched","match_rate"},"mrjti.tsv":{"dataset_id","tissue","gene_id","source_z","source_q_global","source_family_size","selection_rule","selection_relation","status","n_snps","standardized_expression_coefficient","ci_low","ci_high","pvalue","q_mrjti","mrjti_supported","direction_concordant","inference_status","task_dir"},"candidate_block_annotation_bridge.tsv":{"analysis_id","candidate_id","annotation_tissue","variant_id","locus_definition","locus_id","annotation_only"},"stability-summary.tsv":{"comparison_dataset","comparison_level","analysis_set","rank_metric","n_common","spearman_rho","status"},"stability-exclusions.tsv":{"dataset_id","comparison_level","exclusion_status","n_excluded"},"task-ledger.tsv":{"run_id","task_id","attempt","requirement_class","status","output_sha256"}}
    ancillary_data={}
    for name,need in ancillary.items():
        p=directory/name
        if not p.is_file(): report["errors"].append(f"missing ancillary {name}"); continue
        header,data=rows(p)
        ancillary_data[name]=data
        if not need<=set(header): report["errors"].append(f"{name}: ancillary schema mismatch")
        if name=="mrjti.tsv":
            for i,r in enumerate(data,2):
                if r["status"] not in {"ineligible_lt20","harmonization_failed","ld_failed","model_failed","success_inference_bonferroni_selection_dependent"}: report["errors"].append(f"mrjti.tsv:{i}: status")
                if r["pvalue"]!="NA" or r["q_mrjti"]!="NA": report["errors"].append(f"mrjti.tsv:{i}: unsupported P/q fields must be NA")
                if r["status"]=="success_inference_bonferroni_selection_dependent":
                    if any(r[x]=="NA" for x in ("standardized_expression_coefficient","ci_low","ci_high","n_snps","ci_significance","mrjti_family_size","mrjti_family_alpha","mrjti_supported")): report["errors"].append(f"mrjti.tsv:{i}: Bonferroni success incomplete")
                    if r["ci_significance"] not in {"sig","nonsig"} or r["mrjti_supported"]!=str(r["ci_significance"]=="sig").lower(): report["errors"].append(f"mrjti.tsv:{i}: Bonferroni support mismatch")
    # Exact model-match dataset x tissue matrix and arithmetic closure.
    mm=ancillary_data.get("model-match-qc.tsv",[]); observed={(r["dataset_id"],r["tissue"]) for r in mm}; expected={(d,t) for d in datasets for t in tissues}
    if observed!=expected or len(mm)!=len(expected): report["errors"].append("model-match dataset/tissue matrix mismatch")
    for i,r in enumerate(mm,2):
        try:
            n=int(r["n_model_snps"]); matched=int(r["n_matched"]); rate=float(r["match_rate"])
            if n<0 or matched<0 or matched>n or abs(rate-(matched/max(1,n)))>1e-12: report["errors"].append(f"model-match-qc.tsv:{i}: value mismatch")
        except ValueError: report["errors"].append(f"model-match-qc.tsv:{i}: invalid numeric value")
    # MR-JTI uses one native Bonferroni family per GWAS over every planned
    # selected gene-by-tissue task; significance remains selection-dependent.
    triggers={(r["dataset_id"],r["tissue"],r["gene_id"]) for r in loaded.get("pe_tissue_gene_full.tsv",[]) if r["run_status"]=="success" and float(r["q_gene_tissue"])<=.05}
    mr=ancillary_data.get("mrjti.tsv",[]); observed_mr={(r["dataset_id"],r["tissue"],r["gene_id"]) for r in mr}
    if triggers!=observed_mr or len(mr)!=len(triggers): report["errors"].append("MR global-BH candidate closure mismatch")
    family_sizes={d:sum(1 for key in triggers if key[0]==d) for d in datasets}
    for i,r in enumerate(mr,2):
        expected=family_sizes[r["dataset_id"]]
        if int(r["mrjti_family_size"])!=expected or abs(float(r["mrjti_family_alpha"])-0.05/expected)>1e-15: report["errors"].append(f"mrjti.tsv:{i}: family size/alpha mismatch")
    # Stability exclusions must exactly count non-success projection rows.
    ep=directory/"stability-exclusions.tsv"
    if ep.is_file():
        _,erows=rows(ep); observed={(r["dataset_id"],r["comparison_level"],r["exclusion_status"]):int(r["n_excluded"]) for r in erows}; expected={}
        for r in loaded.get("pe_outcome_stability.tsv",[]):
            if r["comparison_status"]!="success":
                k=(r["comparison_dataset"],r["comparison_level"],r["comparison_status"]); expected[k]=expected.get(k,0)+1
        if observed!=expected: report["errors"].append("stability exclusion recount mismatch")
    # Independently read formal detail, then apply the single frozen summary definition.
    sr=ancillary_data.get("stability-summary.tsv",[])
    expected_summary=stability_core.summarize(loaded.get("pe_outcome_stability.tsv",[]))
    for error in stability_core.discrepancies(sr,expected_summary): report["errors"].append("stability summary independent recomputation mismatch: "+error)
    skill_validation=directory/"m1.5-skill-validation.json"
    if not skill_validation.is_file(): report["errors"].append("missing skill-owned M1.5 validation evidence")
    else:
        try:
            s=json.loads(skill_validation.read_text())
            dataset_reports=s.get("datasets",[])
            dataset_ids={x.get("dataset_id") for x in dataset_reports if x.get("status")=="pass"}
            bh_closed=all(x.get("global_bh_family_size")==x.get("global_bh_qvalues_assigned") for x in dataset_reports)
            no_unexpected=all(x.get("unexpected_returned_genes")==0 for x in dataset_reports)
            source=s.get("source_validation",{}); source_path=Path(source.get("path",""))
            source_ok=source_path.is_file() and source.get("sha256")==sha256(source_path)
            if s.get("schema_version")!="m1.5-skill-ingestion-v1" or s.get("status")!="pass" or dataset_ids!=set(datasets) or not bh_closed or not no_unexpected or not source_ok:
                report["errors"].append("skill-owned M1.5 validation evidence invalid")
        except Exception: report["errors"].append("skill-owned M1.5 validation evidence unreadable")
    report["state"] = "PASS" if not report["errors"] else "FAIL"
    out = directory / "validation-report.json"
    out.write_text(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    print(report["state"])
    if report["errors"]: raise SystemExit(4)
    return report


def publish(args: argparse.Namespace) -> None:
    staging = Path(args.staging).resolve()
    context=analysis_context(args.analysis_id)
    report = validate_outputs(argparse.Namespace(directory=str(staging), analysis_id=args.analysis_id))
    run_id = context["analysis_id"]
    if staging.name != ".staging" or staging.parent.name != run_id or staging.parent.parent.name!="analyses":
        die("staging path is not bound to analysis_id")
    ledger_path=staging/"task-ledger.tsv"
    if not ledger_path.is_file(): die("missing task ledger")
    _,ledger=rows(ledger_path)
    final=validate_ledger(ledger)
    _, tissue = rows(staging / "pe_tissue_gene_full.tsv")
    _, cross = rows(staging / "pe_cross_tissue.tsv")
    _, candidates = rows(staging / "pe_candidate_locus_map.tsv")
    tissue_success=sum(r["run_status"]=="success" for r in tissue)
    _,whole=rows(staging/"pe_whole_blood.tsv")
    if not whole: die("Whole Blood is empty")
    expected_tissues=context["tissues"]
    dataset_ids=set(context["datasets"])
    validate_task_matrix(final,dataset_ids,expected_tissues)
    observed={r["tissue"] for r in tissue}
    if expected_tissues!=observed: die("tissue enumeration is not closed")
    for name in ("model-match-qc.tsv","mrjti.tsv","candidate_block_annotation_bridge.tsv","stability-summary.tsv","stability-exclusions.tsv","task-ledger.tsv"):
        if not (staging/name).is_file(): die(f"missing ancillary artifact {name}")
    tissue_candidates=len({r["candidate_id"] for r in candidates if r["statistical_layer"]=="tissue"})
    summary=[{"statistical_layer":"tissue","candidate_rule_id":"tissue_global_bh_q05_v1","applicability":"applicable","n_success":tissue_success,"n_candidates":tissue_candidates,"zero_candidate_reason":"no_candidate_passed_prefrozen_threshold" if tissue_candidates==0 else "NA"}]
    if context["configuration"]["capability_level"] == "degraded_cross_tissue":
        summary.append({"statistical_layer":"cross_tissue","candidate_rule_id":"cross_tissue_bh_q05_v1","applicability":"not_applicable_to_capability","n_success":"NA","n_candidates":"NA","zero_candidate_reason":"NA"})
    else:
        cs=sum(r["model_status"]=="success" for r in cross); cc=len({r["candidate_id"] for r in candidates if r["statistical_layer"]=="cross_tissue"})
        summary.append({"statistical_layer":"cross_tissue","candidate_rule_id":"cross_tissue_bh_q05_v1","applicability":"applicable","n_success":cs,"n_candidates":cc,"zero_candidate_reason":"no_candidate_passed_prefrozen_threshold" if cc==0 else "NA"})
    if not candidates:
        # Header-only is legal only when both locus tasks succeeded and every
        # applicable layer has zero successes passing its frozen threshold.
        if any(x["applicability"]=="applicable" and x["zero_candidate_reason"]!="no_candidate_passed_prefrozen_threshold" for x in summary): die("illegal header-only locus tables")
    files={p.name:{"size_bytes":p.stat().st_size,"sha256":sha256(p)} for p in sorted(staging.iterdir()) if p.is_file() and p.name not in {"artifact-index.json","SUCCESS.json"}}
    artifact_index = {"analysis_id": run_id, "capability_level": context["configuration"]["capability_level"], "validation": report, "files": files, "summary":{"candidate_rules":summary},
                      "task_ledger_status": "synthetic_placeholder",
                      "task_ledger_note": "task-ledger.tsv rows are synthesized from task_id at assembly time (sha256(task_id) in hash columns, single hardcoded timestamp). The publish gate uses it only as a structural enumeration checklist of intended tasks; it is NOT execution provenance and must not be cited as execution evidence. Real provenance = per-input sha256 records and the S-MultiXcan validation release."}
    artifact_path = staging / "artifact-index.json"
    artifact_path.write_text(json.dumps(artifact_index, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    success = {"analysis_id": run_id, "capability_level": context["configuration"]["capability_level"], "cross_tissue_status": context["configuration"]["cross_tissue"], "artifact_index_sha256": sha256(artifact_path)}
    (staging / "SUCCESS.json").write_text(json.dumps(success, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    check=json.loads(artifact_path.read_text()); success_check=json.loads((staging/"SUCCESS.json").read_text())
    if success_check["artifact_index_sha256"]!=sha256(artifact_path) or check["analysis_id"]!=run_id or check["capability_level"]!=context["configuration"]["capability_level"]: die("final artifact/SUCCESS self-check failed")
    for name,meta in check["files"].items():
        p=staging/name
        if not p.is_file() or p.stat().st_size!=meta["size_bytes"] or sha256(p)!=meta["sha256"]: die(f"final artifact digest mismatch: {name}")
    published = staging.parent / "published"
    if published.exists(): die(f"published directory already exists: {published}")
    staging.rename(published)
    pointer = ROOT / "results/01_pe_genetic_map/CURRENT"
    tmp = pointer.with_suffix(".tmp")
    tmp.unlink(missing_ok=True)
    tmp.symlink_to(published.relative_to(pointer.parent))
    tmp.replace(pointer)
    print(published)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("status"); p.add_argument("--analysis-id", required=True); p.set_defaults(func=lambda a: print(json.dumps(analysis_context(a.analysis_id), default=str, sort_keys=True)))
    p = sub.add_parser("validate"); p.add_argument("--analysis-id", required=True); p.add_argument("--directory", required=True); p.set_defaults(func=validate_outputs)
    p = sub.add_parser("publish"); p.add_argument("--analysis-id", required=True); p.add_argument("--staging", required=True); p.set_defaults(func=publish)
    args = parser.parse_args(); args.func(args)


if __name__ == "__main__":
    main()
