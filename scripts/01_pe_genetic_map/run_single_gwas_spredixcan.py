#!/usr/bin/env python3
"""One-pass tissue harmonization and S-PrediXcan for one normalized GWAS."""
from __future__ import annotations

import argparse, csv, gzip, hashlib, heapq, json, math, sqlite3, subprocess, sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).parent))
import execute as ex
import prepare_gwas as prep

RAW_COLUMNS = ["gene","gene_name","zscore","effect_size","pvalue","var_g","pred_perf_r2",
               "pred_perf_pval","pred_perf_qval","n_snps_used","n_snps_in_cov",
               "n_snps_in_model","best_gwas_p","largest_weight"]
RESULT_COLUMNS = ["dataset_id","tissue","gene_id","gene_name","zscore","effect_size","pvalue",
                  "q_global","bonferroni_threshold_tissue","bonferroni_significant_tissue",
                  "var_g","pred_perf_r2","pred_perf_pval","pred_perf_qval","n_snps_used",
                  "n_snps_in_cov","n_snps_in_model","best_gwas_p","largest_weight","run_status"]


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()


def bh(values):
    out=[1.0]*len(values); order=sorted(range(len(values)),key=lambda i:(values[i],i)); q=1.0
    for rank in range(len(order)-1,-1,-1):
        i=order[rank]; q=min(q,values[i]*len(values)/(rank+1)); out[i]=q
    return out


def write_tsv(path, fields, rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    with tmp.open("w",newline="") as f:
        w=csv.DictWriter(f,fields,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    tmp.replace(path)


def harmonize_once(gwas, models, work, seed, audit_per_status):
    maps={}; reverse=defaultdict(list); model_counts={}
    for index,(tissue,db) in enumerate(models):
        alleles=ex.model_alleles(db); maps[tissue]=alleles; model_counts[tissue]=len(alleles)
        for rsid,pair in alleles.items(): reverse[rsid].append((index,pair))
    handles=[]; writers=[]; counts=[Counter() for _ in models]; samples=defaultdict(list)
    for tissue,_ in models:
        path=work/"M1.2_metaxcan_input"/tissue/"gwas.tsv.gz"; path.parent.mkdir(parents=True,exist_ok=True)
        handle=gzip.open(path,"wt",newline=""); handles.append(handle)
        writer=csv.DictWriter(handle,prep.OUT,delimiter="\t",lineterminator="\n",extrasaction="ignore"); writer.writeheader(); writers.append(writer)
    n_input=0; pz=Counter()
    try:
        for row in ex.iter_tsv(gwas):
            n_input+=1; audit=ex.audit_p_z(row["beta"],row["standard_error"],row["pvalue"]); pz[audit]+=1
            if audit=="invalid": raise RuntimeError("invalid beta/SE/P in normalized GWAS")
            for index,pair in reverse.get(row["rsid"],()):
                tissue=models[index][0]; counts[index]["model_overlap"]+=1
                got,status=prep.harmonize(row["effect_allele"],row["non_effect_allele"],row["beta"],row.get("effect_allele_frequency"),pair[0],pair[1])
                if got is None: counts[index][status]+=1; continue
                x=dict(row); x.update(source_ea=row["effect_allele"],source_nea=row["non_effect_allele"],source_beta=row["beta"],
                    source_eaf=row.get("effect_allele_frequency","NA"),effect_allele=pair[0],non_effect_allele=pair[1],
                    beta=got["beta"],zscore=got["beta"]/float(row["standard_error"]),
                    effect_allele_frequency=got["eaf"] if got["eaf"] is not None else "NA",harmonization_status=status)
                writers[index].writerow({k:ex.fmt(x.get(k)) for k in prep.OUT}); counts[index][status]+=1; counts[index]["n_retained"]+=1
                rank=int(hashlib.sha256(f"{seed}\0{tissue}\0{status}\0{row['rsid']}\0{row['variant_id']}".encode()).hexdigest(),16)
                heap=samples[(tissue,status)]; item=(-rank,row["rsid"],x)
                if len(heap)<audit_per_status: heapq.heappush(heap,item)
                elif item[0]>heap[0][0]: heapq.heapreplace(heap,item)
    finally:
        for handle in handles: handle.close()
    qc_rows=[]
    for index,(tissue,db) in enumerate(models):
        c=counts[index]; c["n_input"]=n_input; c["n_model_snps"]=model_counts[tissue]; c["model_snp_absent"]=n_input-c["model_overlap"]
        retained_states=sum(c[x] for x in ("direct","swapped","strand_direct","strand_swapped"))
        exclusions=sum(v for k,v in c.items() if k not in {"n_input","n_model_snps","model_overlap","n_retained","direct","swapped","strand_direct","strand_swapped"})
        c["closure_difference"]=n_input-retained_states-exclusions
        if c["closure_difference"] or retained_states!=c["n_retained"]: raise RuntimeError("harmonization closure failed: "+tissue)
        audit_rows=[x for (tis,_),heap in samples.items() if tis==tissue for _,_,x in heap]
        prep.direction_audit(audit_rows,maps[tissue],seed,audit_per_status)
        rate=c["n_retained"]/max(1,c["n_model_snps"])
        if rate<0.20: raise RuntimeError(f"model match rate below 0.20: {tissue}={rate}")
        qpath=work/"M1.2_metaxcan_input"/tissue/"qc.tsv"
        write_tsv(qpath,["metric","value"],[{"metric":k,"value":v} for k,v in sorted(c.items())])
        qc_rows.append({"tissue":tissue,"n_input":n_input,"n_model_snps":c["n_model_snps"],"n_model_overlap":c["model_overlap"],
                        "n_retained":c["n_retained"],"match_rate":f"{rate:.17g}","direction_audit":"pass","closure_difference":0})
    return qc_rows,dict(pz)


def run_one(model, args, work, result_raw):
    tissue,db=model; covariance=db.with_suffix(".txt.gz"); gwas=work/"M1.2_metaxcan_input"/tissue/"gwas.tsv.gz"
    output=result_raw/f"{tissue}.csv"; log=result_raw/"logs"; log.mkdir(parents=True,exist_ok=True)
    cmd=[args.python_metaxcan,str(ROOT/"tools/MetaXcan/software/SPrediXcan.py"),"--model_db_path",str(db),"--covariance",str(covariance),
         "--gwas_file",str(gwas),"--snp_column","rsid","--effect_allele_column","effect_allele","--non_effect_allele_column","non_effect_allele",
         "--zscore_column","zscore","--beta_column","beta","--se_column","standard_error","--pvalue_column","pvalue","--output_file",str(output),
         "--remove_ens_version","--additional_output","--throw","--overwrite"]
    p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
    (log/f"{tissue}.stdout.log").write_text(p.stdout); (log/f"{tissue}.stderr.log").write_text(p.stderr)
    if p.returncode or not output.is_file(): raise RuntimeError(f"S-PrediXcan failed: {tissue}: {p.stderr[-1000:]}")
    with output.open(newline="") as f:
        rd=csv.DictReader(f)
        if rd.fieldnames!=RAW_COLUMNS: raise RuntimeError(f"unexpected S-PrediXcan schema: {tissue}: {rd.fieldnames}")
        rows=list(rd)
    return tissue,db,covariance,output,cmd,rows


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--dataset-id",required=True); ap.add_argument("--normalized-gwas",type=Path,required=True)
    ap.add_argument("--model-dir",type=Path,required=True); ap.add_argument("--python-metaxcan",required=True); ap.add_argument("--work-dir",type=Path,required=True)
    ap.add_argument("--result-dir",type=Path,required=True); ap.add_argument("--workers",type=int,default=6); ap.add_argument("--seed",type=int,default=20260714)
    ap.add_argument("--audit-per-status",type=int,default=100); args=ap.parse_args()
    models=[]
    for db in sorted(args.model_dir.glob("JTI_*.db")):
        covariance=db.with_suffix(".txt.gz")
        if not covariance.is_file(): raise SystemExit("missing covariance: "+str(covariance))
        models.append((ex.tissue(db),db))
    if len(models)!=49: raise SystemExit(f"expected 49 paired models, found {len(models)}")
    args.work_dir.mkdir(parents=True,exist_ok=True); args.result_dir.mkdir(parents=True,exist_ok=True); raw=args.result_dir/"tissues"; raw.mkdir(parents=True,exist_ok=True)
    match_qc,pz=harmonize_once(args.normalized_gwas,models,args.work_dir,args.seed,args.audit_per_status)
    completed=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        future={pool.submit(run_one,m,args,args.work_dir,raw):m[0] for m in models}
        for f in as_completed(future): completed.append(f.result())
    completed.sort(key=lambda x:x[0]); all_rows=[]; artifacts=[]; seen=set()
    for tissue,db,covariance,output,cmd,rows in completed:
        threshold=.05/max(1,len(rows))
        for row in rows:
            gid=row["gene"].split(".")[0]; key=(tissue,gid)
            if key in seen: raise RuntimeError("standardized gene collision: "+tissue+":"+gid)
            seen.add(key); z=float(row["zscore"]); p=float(row["pvalue"])
            if not math.isfinite(z) or not math.isfinite(p) or not 0<=p<=1: raise RuntimeError("nonfinite association")
            record={"dataset_id":args.dataset_id,"tissue":tissue,"gene_id":gid,"gene_name":row["gene_name"],"zscore":row["zscore"],
                    "effect_size":row["effect_size"],"pvalue":row["pvalue"],"bonferroni_threshold_tissue":f"{threshold:.17g}",
                    "bonferroni_significant_tissue":"true" if p<=threshold else "false","run_status":"success"}
            record.update({k:row[k] for k in RAW_COLUMNS if k not in {"gene","gene_name","zscore","effect_size","pvalue"}}); all_rows.append(record)
        artifacts.append({"tissue":tissue,"model_db":str(db.resolve()),"model_db_sha256":sha256(db),"covariance":str(covariance.resolve()),
                          "covariance_sha256":sha256(covariance),"output":str(output.resolve()),"output_sha256":sha256(output),"n_rows":len(rows),"command":cmd})
    qs=bh([float(r["pvalue"]) for r in all_rows])
    for row,q in zip(all_rows,qs): row["q_global"]=f"{q:.17g}"
    global_bh_family_size=len(all_rows)
    global_bh_q_values_assigned=sum(
        r["run_status"]=="success" and math.isfinite(float(r["q_global"])) for r in all_rows
    )
    if global_bh_q_values_assigned!=global_bh_family_size:
        raise RuntimeError("global BH family/q-value closure failed")
    all_rows.sort(key=lambda r:(r["tissue"],float(r["pvalue"]),r["gene_id"]))
    write_tsv(args.result_dir/"spredixcan_all_tissues.tsv",RESULT_COLUMNS,all_rows)
    write_tsv(args.result_dir/"significant_gene_tissue.tsv",RESULT_COLUMNS,[r for r in all_rows if r["bonferroni_significant_tissue"]=="true" or float(r["q_global"])<=.05])
    write_tsv(args.result_dir/"model_match_qc.tsv",["tissue","n_input","n_model_snps","n_model_overlap","n_retained","match_rate","direction_audit","closure_difference"],sorted(match_qc,key=lambda r:r["tissue"]))
    method={"schema_version":"m1.3-spredixcan-single-gwas-v1","dataset_id":args.dataset_id,"normalized_gwas":str(args.normalized_gwas.resolve()),
            "normalized_gwas_sha256":sha256(args.normalized_gwas),"method":"MetaXcan-0.7.5/JTI-v8/zscore-beta-over-se",
            "phi_calibration":False,"global_multiple_testing":"BH across all successful gene-tissue rows",
            "global_bh_family_definition":"unique (gene_id,tissue) rows with run_status=success and finite pvalue in [0,1], within this GWAS only",
            "global_bh_family_size":global_bh_family_size,"global_bh_q_values_assigned":global_bh_q_values_assigned,
            "tissue_trigger":"Bonferroni 0.05/n successful genes per tissue","seed":args.seed,"p_z_audit":pz,
            "n_tissues":len(models),"n_associations":len(all_rows),"n_significant_union":sum(r["bonferroni_significant_tissue"]=="true" or float(r["q_global"])<=.05 for r in all_rows),
            "artifacts":artifacts}
    (args.result_dir/"method.json").write_text(json.dumps(method,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:method[k] for k in ("dataset_id","n_tissues","n_associations","global_bh_family_size","global_bh_q_values_assigned","n_significant_union","p_z_audit")},sort_keys=True))


if __name__=="__main__": main()
