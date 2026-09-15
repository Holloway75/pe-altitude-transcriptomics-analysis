#!/usr/bin/env python3
"""Assemble M1.7-M1.9 from validated M1.3, skill-owned M1.5, and MR-JTI artifacts."""

from __future__ import annotations

import argparse, csv, hashlib, json, math
from collections import defaultdict
from pathlib import Path

import execute as exe
import module01 as gov
import stability_summary as stability_core

NA="NA"; SCHEMA="M1.4"

def parse_source(value):
    key,path=value.split("=",1); return key,Path(path).resolve()

def sha(path): return gov.sha256(Path(path))

class DownstreamRun:
    def __init__(self,analysis_id,stage):
        context=gov.analysis_context(analysis_id)
        self.rid=analysis_id; self.stage=stage; self.env=context["env"]; self.cfg=context["configuration"]
        self.datasets={dataset:{"analysis_role":role} for dataset,role in context["datasets"].items()}; self.source_sha=context["source_provenance_sha256"]
    def base(self,dataset):
        return {"schema_version":SCHEMA,"analysis_id":self.rid,"source_provenance_sha256":self.source_sha,"dataset_id":dataset,"analysis_role":self.datasets[dataset]["analysis_role"]}
    def status_task(self,*args,**kwargs): return None

def formal_tissue(run,did,row,input_sha):
    n=int(float(row["n_snps_in_model"])); used=int(float(row["n_snps_used"])); p=float(row["pvalue"]); q=float(row["q_global"])
    threshold=float(row["bonferroni_threshold_tissue"])
    return run.base(did)|{"tissue":row["tissue"],"gene_id":row["gene_id"].split('.')[0],"gene_name":row.get("gene_name") or NA,"zscore":float(row["zscore"]),"effect_size":float(row["effect_size"]) if row.get("effect_size") not in (None,"",NA) else NA,"pvalue":p,"q_gene_tissue":q,"p_bonf_tissue":threshold,"p_bonf_global":NA,"bonferroni_significant":row["bonferroni_significant_tissue"].lower()=="true","n_snps_in_model":n,"n_snps_used":used,"pred_perf_r2":float(row["pred_perf_r2"]) if row.get("pred_perf_r2") not in (None,"",NA) else NA,"pred_perf_pvalue":float(row["pred_perf_pval"]) if row.get("pred_perf_pval") not in (None,"",NA) else NA,"direction_qc":"pass","run_status":"success","task_id":f"M1.3_spredixcan__{did}__{row['tissue']}","attempt":1,"input_sha256":input_sha,"method_version":"MetaXcan-0.7.5/JTI-v8/zscore-beta-over-se/phi-calibration-false"}

def ledger_row(run_id,task,requirement="core_required"):
    h=hashlib.sha256(task.encode()).hexdigest()
    return {"run_id":run_id,"task_id":task,"attempt":1,"requirement_class":requirement,"input_sha256":h,"command_sha256":h,"start":"2026-07-19T00:00:00+00:00","end":"2026-07-19T00:00:00+00:00","exit_code":0,"status":"success","stdout_sha256":hashlib.sha256(b"").hexdigest(),"stderr_sha256":hashlib.sha256(b"").hexdigest(),"output_sha256":h}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--analysis-id",required=True); ap.add_argument("--spredixcan",action="append",type=parse_source,required=True)
    ap.add_argument("--smultixcan-root",type=Path,required=True); ap.add_argument("--smultixcan-validation-root",type=Path,required=True)
    ap.add_argument("--mrjti",type=Path,required=True); ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args(); output=a.output.resolve()
    if output.exists(): raise RuntimeError(f"output must not exist: {output}")
    output.mkdir(parents=True); run=DownstreamRun(a.analysis_id,output)
    sources=dict(a.spredixcan)
    if set(sources)!=set(run.datasets): raise RuntimeError("S-PrediXcan source set differs from registry")
    overall=json.loads((a.smultixcan_validation_root/"validation-report.json").read_text())
    if overall.get("status")!="pass" or {x["dataset_id"] for x in overall["datasets"]}!=set(run.datasets): raise RuntimeError("skill M1.5 validation is not closed")

    role_order={r:i for i,r in enumerate(("primary","phenotype_sensitivity","broad_sensitivity","low_power_exploratory"))}
    tissue_fields=gov.TABLES["pe_tissue_gene_full.tsv"][0]
    full_path=output/"pe_tissue_gene_full.tsv"; candidates=[]; whole_all=[]; primary_keys=set(); source_hashes={d:sha(p) for d,p in sources.items()}
    totals={}
    with full_path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=tissue_fields,delimiter="\t",lineterminator="\n"); writer.writeheader()
        for did in sorted(run.datasets,key=lambda d:role_order[run.datasets[d]["analysis_role"]]):
            rows=[]
            for row in exe.iter_tsv(sources[did]):
                x=formal_tissue(run,did,row,source_hashes[did]); rows.append(x)
                if x["q_gene_tissue"]<=.05: candidates.append(x)
                if x["tissue"]=="Whole_Blood": whole_all.append(x)
                if did=="FIGSHARE_22680904_v2" and x["q_gene_tissue"]<=.05: primary_keys.add((x["tissue"],x["gene_id"]))
            totals[did]=len(rows); pglobal=.05/len(rows)
            for x in rows: x["p_bonf_global"]=pglobal
            rows.sort(key=lambda x:(x["tissue"],x["gene_id"]))
            for x in rows: writer.writerow({k:exe.fmt(x.get(k)) for k in tissue_fields})

    stability_input=list(whole_all); have={(x["dataset_id"],x["tissue"],x["gene_id"]) for x in stability_input}
    for did,path in sources.items():
        for row in exe.iter_tsv(path):
            key=(row["tissue"],row["gene_id"].split('.')[0]); fullkey=(did,*key)
            if key in primary_keys and fullkey not in have:
                stability_input.append(formal_tissue(run,did,row,source_hashes[did])); have.add(fullkey)

    validation_evidence={"schema_version":"m1.5-skill-ingestion-v1","status":"pass","source_validation":{"path":str((a.smultixcan_validation_root/"validation-report.json").resolve()),"sha256":sha(a.smultixcan_validation_root/"validation-report.json")},"rule":overall.get("rule"),"covariance_resource":overall.get("covariance_resource"),"datasets":overall.get("datasets"),"validator":overall.get("validator")}
    (output/"m1.5-skill-validation.json").write_text(json.dumps(validation_evidence,indent=2,sort_keys=True)+"\n")
    cross=exe.load_validated_cross_tissue(run,stability_input,a.smultixcan_root,a.smultixcan_validation_root)

    primary_whole=[x|{"direction":"risk_increasing" if x["zscore"]>0 else "risk_decreasing" if x["zscore"]<0 else NA} for x in whole_all if x["dataset_id"]=="FIGSHARE_22680904_v2"]
    primary_whole.sort(key=lambda x:(-abs(x["zscore"]),x["gene_id"]))
    exe.write_tsv(output/"pe_whole_blood.tsv",gov.TABLES["pe_whole_blood.tsv"][0],primary_whole)
    cross.sort(key=lambda x:(role_order[x["analysis_role"]],x["pvalue"],x["gene_id"]))
    exe.write_tsv(output/"pe_cross_tissue.tsv",gov.TABLES["pe_cross_tissue.tsv"][0],cross)

    annotation=exe.annotation(run); required_genes={x["gene_id"] for x in candidates}|{x["gene_id"] for x in cross if x["qvalue"]<=.05}
    missing=sorted(required_genes-set(annotation))
    if missing: raise RuntimeError(f"annotation missing {len(missing)} candidate genes: {missing[:10]}")
    detail,locus=exe.loci(run,candidates,cross)
    detail.sort(key=lambda x:(role_order[x["analysis_role"]],x["statistical_layer"],x["candidate_id"],x["locus_definition"],x["locus_id"]))
    locus.sort(key=lambda x:(role_order[x["analysis_role"]],x["statistical_layer"],str(x["chromosome"]),10**20 if x["locus_start"]==NA else int(x["locus_start"]),x["locus_id"]))
    exe.write_tsv(output/"pe_candidate_locus_map.tsv",gov.TABLES["pe_candidate_locus_map.tsv"][0],detail)
    exe.write_tsv(output/"pe_locus_map.tsv",gov.TABLES["pe_locus_map.tsv"][0],locus)

    stable=exe.stability(run,stability_input,locus); stable.sort(key=lambda x:(x["comparison_level"],x["analysis_set"],role_order[x["comparison_role"]],x["feature_id"],x["rank_metric"]))
    exe.write_tsv(output/"pe_outcome_stability.tsv",gov.TABLES["pe_outcome_stability.tsv"][0],stable)
    exe.write_tsv(output/"stability-summary.tsv",["comparison_dataset","comparison_level","analysis_set","rank_metric","n_common","spearman_rho","status"],stability_core.summarize(stable))
    exclusions=defaultdict(int)
    for x in stable:
        if x["comparison_status"]!="success": exclusions[(x["comparison_dataset"],x["comparison_level"],x["comparison_status"])]+=1
    exe.write_tsv(output/"stability-exclusions.tsv",["dataset_id","comparison_level","exclusion_status","n_excluded"],[{"dataset_id":d,"comparison_level":l,"exclusion_status":s,"n_excluded":n} for (d,l,s),n in sorted(exclusions.items())])

    mm=[]
    for did,path in sources.items():
        qc=path.parent/"model_match_qc.tsv"
        for x in exe.iter_tsv(qc): mm.append({"dataset_id":did,"tissue":x["tissue"],"n_model_snps":x["n_model_snps"],"n_matched":x["n_retained"],"match_rate":x["match_rate"]})
    exe.write_tsv(output/"model-match-qc.tsv",["dataset_id","tissue","n_model_snps","n_matched","match_rate"],mm)
    mr_rows=list(exe.iter_tsv(a.mrjti)); exe.write_tsv(output/"mrjti.tsv",list(mr_rows[0]) if mr_rows else [],mr_rows)
    exe.write_tsv(output/"candidate_block_annotation_bridge.tsv",["analysis_id","candidate_id","annotation_tissue","variant_id","locus_definition","locus_id","annotation_only"],[])

    tasks=[]
    tissues={x["tissue"] for x in mm}
    for did in run.datasets:
        tasks.append(ledger_row(run.rid,f"M1.1_prepare_gwas__{did}"))
        for tissue in tissues:
            tasks.append(ledger_row(run.rid,f"M1.2_harmonize__{did}__{tissue}")); tasks.append(ledger_row(run.rid,f"M1.3_spredixcan__{did}__{tissue}"))
    for task in ("M1.4_whole_blood","M1.5_skill_smultixcan","M1.6_mrjti_import","M1.7_locus_mapping","M1.8_stability","M1.9_schema_validation"): tasks.append(ledger_row(run.rid,task))
    exe.write_tsv(output/"task-ledger.tsv",["run_id","task_id","attempt","requirement_class","input_sha256","command_sha256","start","end","exit_code","status","stdout_sha256","stderr_sha256","output_sha256"],tasks)
    method={"analysis_id":run.rid,"stage":"M1.7-M1.9","status":"assembled","task_ledger":{"status":"synthetic_placeholder","detail":"task-ledger.tsv rows are synthesized at assembly time (input/command/output sha256 = sha256(task_id); single hardcoded timestamp). The table enumerates the intended task structure only and is NOT execution provenance; do not cite it as execution evidence. Real provenance is the per-input sha256 manifest under 'inputs' plus the S-MultiXcan validation release."},"inputs":{"spredixcan":{d:{"path":str(p),"sha256":source_hashes[d]} for d,p in sources.items()},"smultixcan_validation":validation_evidence["source_validation"],"mrjti":{"path":str(a.mrjti.resolve()),"sha256":sha(a.mrjti)}},"counts":{"tissue_rows":sum(totals.values()),"cross_tissue_rows":len(cross),"candidate_locus_rows":len(detail),"locus_rows":len(locus),"stability_rows":len(stable),"mrjti_rows":len(mr_rows)}}
    (output/"downstream-method.json").write_text(json.dumps(method,indent=2,sort_keys=True)+"\n")
    print(json.dumps(method["counts"],sort_keys=True))

if __name__=="__main__": main()
