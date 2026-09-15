#!/usr/bin/env python3
"""Resume M1.7-M1.9 from already assembled formal M1.3/M1.5 tables."""
from __future__ import annotations
import argparse, json
from collections import defaultdict
from pathlib import Path
import assemble_m1_downstream as asm
import execute as exe
import module01 as gov
import stability_summary as stability_core

NA="NA"
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--analysis-id",required=True); ap.add_argument("--staging",type=Path,required=True); ap.add_argument("--spredixcan",action="append",type=asm.parse_source,required=True); ap.add_argument("--mrjti",type=Path,required=True); a=ap.parse_args()
    output=a.staging.resolve(); run=asm.DownstreamRun(a.analysis_id,output); sources=dict(a.spredixcan)
    role_order={r:i for i,r in enumerate(("primary","phenotype_sensitivity","broad_sensitivity","low_power_exploratory"))}
    candidates=[]; whole=[]; primary_keys=set()
    for x in exe.iter_tsv(output/"pe_tissue_gene_full.tsv"):
        row={k:(float(v) if k in {"zscore","effect_size","pvalue","q_gene_tissue","p_bonf_tissue","p_bonf_global","pred_perf_r2","pred_perf_pvalue"} and v!=NA else int(v) if k in {"n_snps_in_model","n_snps_used","attempt"} and v!=NA else v=="true" if k=="bonferroni_significant" else v) for k,v in x.items()}
        if row["q_gene_tissue"]<=.05: candidates.append(row)
        if row["tissue"]=="Whole_Blood": whole.append(row)
        if row["dataset_id"]=="FIGSHARE_22680904_v2" and row["q_gene_tissue"]<=.05: primary_keys.add((row["tissue"],row["gene_id"]))
    stability_input=list(whole); have={(x["dataset_id"],x["tissue"],x["gene_id"]) for x in stability_input}
    for x in exe.iter_tsv(output/"pe_tissue_gene_full.tsv"):
        key=(x["tissue"],x["gene_id"]); fullkey=(x["dataset_id"],*key)
        if key in primary_keys and fullkey not in have:
            row={k:(float(v) if k in {"zscore","effect_size","pvalue","q_gene_tissue","p_bonf_tissue","p_bonf_global","pred_perf_r2","pred_perf_pvalue"} and v!=NA else int(v) if k in {"n_snps_in_model","n_snps_used","attempt"} and v!=NA else v=="true" if k=="bonferroni_significant" else v) for k,v in x.items()}; stability_input.append(row); have.add(fullkey)
    cross=[]
    for x in exe.iter_tsv(output/"pe_cross_tissue.tsv"):
        cross.append({k:(float(v) if k in {"statistic","pvalue","qvalue","best_tissue_z","condition_number"} and v!=NA else int(v) if k in {"n_tissues","n_independent_components","attempt"} and v!=NA else False if k=="direction_available" else v) for k,v in x.items()})
    annotation=exe.annotation(run); required={x["gene_id"] for x in candidates}|{x["gene_id"] for x in cross if x["qvalue"]<=.05}; missing=sorted(required-set(annotation))
    detail,locus=exe.loci(run,candidates,cross)
    detail.sort(key=lambda x:(role_order[x["analysis_role"]],x["statistical_layer"],x["candidate_id"],x["locus_definition"],x["locus_id"])); locus.sort(key=lambda x:(role_order[x["analysis_role"]],x["statistical_layer"],str(x["chromosome"]),10**20 if x["locus_start"]==NA else int(x["locus_start"]),x["locus_id"]))
    exe.write_tsv(output/"pe_candidate_locus_map.tsv",gov.TABLES["pe_candidate_locus_map.tsv"][0],detail); exe.write_tsv(output/"pe_locus_map.tsv",gov.TABLES["pe_locus_map.tsv"][0],locus)
    stable=exe.stability(run,stability_input,locus); stable.sort(key=lambda x:(x["comparison_level"],x["analysis_set"],role_order[x["comparison_role"]],x["feature_id"],x["rank_metric"])); exe.write_tsv(output/"pe_outcome_stability.tsv",gov.TABLES["pe_outcome_stability.tsv"][0],stable)
    exe.write_tsv(output/"stability-summary.tsv",["comparison_dataset","comparison_level","analysis_set","rank_metric","n_common","spearman_rho","status"],stability_core.summarize(stable))
    exclusions=defaultdict(int)
    for x in stable:
        if x["comparison_status"]!="success": exclusions[(x["comparison_dataset"],x["comparison_level"],x["comparison_status"])]+=1
    exe.write_tsv(output/"stability-exclusions.tsv",["dataset_id","comparison_level","exclusion_status","n_excluded"],[{"dataset_id":d,"comparison_level":l,"exclusion_status":s,"n_excluded":n} for (d,l,s),n in sorted(exclusions.items())])
    mm=[]
    for did,path in sources.items():
        for x in exe.iter_tsv(path.parent/"model_match_qc.tsv"): mm.append({"dataset_id":did,"tissue":x["tissue"],"n_model_snps":x["n_model_snps"],"n_matched":x["n_retained"],"match_rate":x["match_rate"]})
    exe.write_tsv(output/"model-match-qc.tsv",["dataset_id","tissue","n_model_snps","n_matched","match_rate"],mm)
    mr=list(exe.iter_tsv(a.mrjti)); exe.write_tsv(output/"mrjti.tsv",list(mr[0]),mr); exe.write_tsv(output/"candidate_block_annotation_bridge.tsv",["analysis_id","candidate_id","annotation_tissue","variant_id","locus_definition","locus_id","annotation_only"],[])
    tasks=[]; tissues={x["tissue"] for x in mm}
    for did in run.datasets:
        tasks.append(asm.ledger_row(run.rid,f"M1.1_prepare_gwas__{did}"))
        for tissue in tissues: tasks.extend((asm.ledger_row(run.rid,f"M1.2_harmonize__{did}__{tissue}"),asm.ledger_row(run.rid,f"M1.3_spredixcan__{did}__{tissue}")))
    for task in ("M1.4_whole_blood","M1.5_skill_smultixcan","M1.6_mrjti_import","M1.7_locus_mapping","M1.8_stability","M1.9_schema_validation"): tasks.append(asm.ledger_row(run.rid,task))
    exe.write_tsv(output/"task-ledger.tsv",["run_id","task_id","attempt","requirement_class","input_sha256","command_sha256","start","end","exit_code","status","stdout_sha256","stderr_sha256","output_sha256"],tasks)
    method={"analysis_id":run.rid,"stage":"M1.7-M1.9","status":"assembled","task_ledger":{"status":"synthetic_placeholder","detail":"task-ledger.tsv rows are synthesized at assembly time (sha256 columns = sha256(task_id); hardcoded timestamp). Enumerates intended task structure only; NOT execution provenance and must not be cited as execution evidence."},"counts":{"tissue_candidates":len(candidates),"cross_tissue_rows":len(cross),"candidate_locus_rows":len(detail),"locus_rows":len(locus),"stability_rows":len(stable),"mrjti_rows":len(mr),"candidate_genes_without_frozen_annotation":len(missing)},"candidate_genes_without_frozen_annotation_examples":missing[:50]}; (output/"downstream-method.json").write_text(json.dumps(method,indent=2,sort_keys=True)+"\n"); print(json.dumps(method["counts"],sort_keys=True))
if __name__=="__main__": main()
