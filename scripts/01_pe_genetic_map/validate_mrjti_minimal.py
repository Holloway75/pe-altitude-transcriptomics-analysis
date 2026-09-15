#!/usr/bin/env python3
"""Strict closure validator with path-independent isolated-run comparison."""
import argparse, csv, gzip, hashlib, json, math, re
from collections import Counter
from pathlib import Path

RSID=re.compile(r"^rs[0-9]+$"); ALLELE=re.compile(r"^[ACGT]+$")
SPREDIXCAN_METHOD="MetaXcan-0.7.5/JTI-v8/zscore-beta-over-se/phi-calibration-false"

def iter_rows(path):
    op=gzip.open if str(path).endswith(".gz") else open
    with op(path,"rt",newline="") as f: yield from csv.DictReader(f,delimiter="\t")
def rows(path): return list(iter_rows(path))

def file_sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def close(a,b,tol=1e-12):
    try: return math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol)
    except (TypeError,ValueError): return False

def normalized_json(value):
    if isinstance(value,dict): return {k:normalized_json(v) for k,v in sorted(value.items()) if k not in {"path","created_at","started_at","finished_at"}}
    if isinstance(value,list): return [normalized_json(x) for x in value]
    return value

def scientific_artifacts(d):
    patterns=["eqtl/SH2B3.txt","ld_reference/ld-reference-variant-map.tsv","ld_reference/SH2B3_EUR_hg19_ldscore.score.ld","mrjti-minimal-snp-flow.tsv","mrjti-minimal.tsv"]
    for od in sorted((d/"outcomes").glob("*")):
        if od.is_dir():
            patterns += [str((od/x).relative_to(d)) for x in ("harmonized.tsv","pruned-variants.tsv","mrjti-input.tsv","prune.prune.in","prune-validation.vcor","mrjti-result.tsv","bootstrap.tsv") if (od/x).is_file()]
    out={rel:file_sha(d/rel) for rel in patterns if (d/rel).is_file()}
    method=normalized_json(json.loads((d/"mrjti-minimal-method.json").read_text()))
    out["mrjti-minimal-method.normalized.json"]=hashlib.sha256(json.dumps(method,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    return out

def artifact_digest(items):
    return hashlib.sha256("\n".join(f"{k}\0{v}" for k,v in sorted(items.items())).encode()).hexdigest()

def audit_artifacts(d):
    files=sorted((d/"logs").glob("*"))+sorted((d/"outcomes").glob("*/plink*.log"))+sorted((d/"outcomes").glob("*/mrjti.log"))
    return {str(p.relative_to(d)):file_sha(p) for p in files if p.is_file()}

def same_variant(a,b):
    return (a["rsid"]==b["rsid"] and str(a["chromosome"]).removeprefix("chr")==str(b["chromosome"]).removeprefix("chr")
            and int(float(a["position"]))==int(float(b["position"]))
            and {a["effect_allele"].upper(),a["other_allele"].upper()}=={b["ref"].upper(),b["alt"].upper()})

def eqtl_closes(x,e):
    return (str(x["chromosome"]).removeprefix("chr")==str(e["Chr"]).removeprefix("chr") and int(float(x["position"]))==int(float(e["BP"]))
            and (x["effect_allele"].upper(),x["other_allele"].upper())==(e["A1"].upper(),e["A2"].upper())
            and all(close(x[a],e[b]) for a,b in (("eqtl_beta","b"),("eqtl_se","SE"),("eqtl_p","p"))))

def gwas_closes(x,g):
    return (str(g["chromosome"]).removeprefix("chr")==str(x["chromosome"]).removeprefix("chr") and int(float(g["position"]))==int(float(x["position"]))
            and (g["effect_allele"].upper(),g["non_effect_allele"].upper())==(x["gwas_source_effect_allele"].upper(),x["gwas_source_other_allele"].upper())
            and all(close(x[a],g[b]) for a,b in (("gwas_source_beta","beta"),("gwas_source_se","standard_error"),("gwas_source_p","pvalue")))
            and close(x["gwas_beta"],float(x["gwas_alignment_sign"])*float(x["gwas_source_beta"])) and close(x["gwas_se"],x["gwas_source_se"]) and close(x["gwas_p"],x["gwas_source_p"]))

def validate(directory):
    d=directory.resolve(); errors=[]; required=[]
    method=json.loads((d/"mrjti-minimal-method.json").read_text()); inputs=method.get("inputs",{})
    eq_path=d/"eqtl/SH2B3.txt"; eq=rows(eq_path)
    if len(eq)!=4273: errors.append(f"expected 4273 v1-provider raw cis-eQTL rows, found {len(eq)}")
    def valid_eqtl(r):
        try: return math.isfinite(float(r["b"])) and math.isfinite(float(r["SE"])) and float(r["SE"])>0 and math.isfinite(float(r["p"]))
        except (KeyError,TypeError,ValueError): return False
    canonical=[r for r in eq if RSID.fullmatch(r["SNP"]) and valid_eqtl(r)]
    multiplicity=Counter(r["SNP"] for r in canonical); mult={rid for rid,n in multiplicity.items() if n>1}; eligible={r["SNP"]:r for r in canonical if r["SNP"] not in mult}
    ref=rows(d/"ld_reference/ld-reference-variant-map.tsv"); byrs={r["rsid"]:r for r in ref}
    score_path=d/"ld_reference/SH2B3_EUR_hg19_ldscore.score.ld"
    with score_path.open() as h:
        h.readline(); score_rows={x[0]:x[-1] for line in h if len(x:=line.split())>=2}
    if len(byrs)!=len(ref): errors.append("reference rsID is not unique")
    if len({(r["chromosome"],r["position"]) for r in ref})!=len(ref): errors.append("reference position is not unique")
    for r in ref:
        if not RSID.fullmatch(r["rsid"]): errors.append(f"noncanonical reference rsID: {r['rsid']}")
        if not ALLELE.fullmatch(r["ref"]) or not ALLELE.fullmatch(r["alt"]): errors.append(f"nonsequence reference allele: {r['rsid']}")
        if r["internal_id"]!=f"{r['chromosome']}:{r['position']}:{r['ref']}:{r['alt']}": errors.append(f"internal ID mismatch: {r['rsid']}")
        try:
            maf,miss,score=float(r["maf"]),float(r["missing_rate"]),float(r["ldscore"])
            if not (.01<=maf<=.5) or not (0<=miss<=.02) or not math.isfinite(score): raise ValueError
        except ValueError: errors.append(f"reference QC value invalid: {r['rsid']}")
        if r["genotype_id"] not in score_rows or not close(r["ldscore"],score_rows.get(r["genotype_id"])): errors.append(f"raw LD-score closure failure: {r['rsid']}")
    if len(score_rows)!=len(ref): errors.append("raw LD-score/reference row count mismatch")
    for key in ("eqtl","pvar","ld_score","reference_bed","reference_bim","reference_fam"):
        evidence=inputs.get(key,{})
        if not evidence.get("path") or not Path(evidence["path"]).is_file() or file_sha(evidence["path"])!=evidence.get("sha256"): errors.append(f"method input evidence mismatch: {key}")
    result=rows(d/"mrjti-minimal.tsv")
    if len(result)!=4: errors.append("expected four outcome results")
    selection=method.get("selection",{}); source_by={}
    if selection.get("mode")=="explicit_v2_spredixcan":
        sp=Path(selection.get("path",""))
        if not sp.is_file() or file_sha(sp)!=selection.get("sha256"): errors.append("explicit S-PrediXcan selection evidence mismatch")
        else:
            selected=[x for x in iter_rows(sp) if x.get("tissue")==method.get("tissue") and x.get("gene_id","").split(".")[0]==method.get("gene_id")]
            counts=Counter(x.get("dataset_id") for x in selected); source_by={x["dataset_id"]:x for x in selected if counts[x.get("dataset_id")]==1}
            if selection.get("method_version")!=SPREDIXCAN_METHOD or any(x.get("method_version")!=SPREDIXCAN_METHOD or x.get("run_status")!="success" for x in selected): errors.append("explicit S-PrediXcan method/run-status gate mismatch")
            if selection.get("family_size") in (None,0) or source_by.get(selection.get("selection_dataset"),{}).get("bonferroni_significant")!="true": errors.append("explicit S-PrediXcan candidate-rule gate mismatch")
    elif selection.get("mode")!="historical_candidate_identity_only": errors.append("unknown candidate-source mode")
    for r in result:
        did=r["dataset_id"]; od=d/"outcomes"/did; harmonized=rows(od/"harmonized.tsv"); pruned=rows(od/"pruned-variants.tsv"); inp=rows(od/"mrjti-input.tsv")
        if selection.get("mode")=="historical_candidate_identity_only":
            if r.get("source_z")!="NA" or r.get("direction_concordant")!="NA" or r.get("source_z_mode")!="historical_candidate_identity_only": errors.append(f"{did}: historical mode leaked source Z/direction")
        else:
            sr=source_by.get(did)
            if sr is None or not close(r.get("source_z"),sr.get("zscore")) or r.get("source_z_mode")!="explicit_v2_spredixcan": errors.append(f"{did}: explicit source Z closure failure")
        ginfo=inputs.get("normalized_gwas",{}).get(did,{}); gpath=Path(ginfo.get("path",""))
        if not gpath.is_file() or file_sha(gpath)!=ginfo.get("sha256"): errors.append(f"{did}: normalized GWAS evidence mismatch"); gby={}
        else:
            wanted={x["rsid"] for x in inp}; grows=[x for x in iter_rows(gpath) if x.get("rsid") in wanted]; gcount=Counter(x["rsid"] for x in grows)
            if any(n!=1 for n in gcount.values()): errors.append(f"{did}: normalized GWAS target rsID is not unique")
            gby={x["rsid"]:x for x in grows if gcount[x["rsid"]]==1}
        prune_ids={x.strip() for x in (od/"prune.prune.in").read_text().splitlines() if x.strip()}
        closure={"rsid","internal_id","genotype_id","chromosome","position","effect_allele","other_allele","ldscore","gwas_source_effect_allele","gwas_source_other_allele","gwas_source_beta","gwas_source_se","gwas_source_p","gwas_alignment_sign"}
        if not inp or not closure<=set(inp[0]): errors.append(f"{did}: MR input lacks closure columns: {sorted(closure-(set(inp[0]) if inp else set()))}"); input_ids=set()
        else: input_ids={x["genotype_id"] for x in inp}
        if prune_ids!={x["genotype_id"] for x in pruned} or prune_ids!=input_ids: errors.append(f"{did}: prune/input set mismatch")
        if len(inp)!=int(r["n_snps"]): errors.append(f"{did}: result/input count mismatch")
        hby={x["rsid"]:x for x in harmonized}
        for x in inp:
            if not closure<=set(x): continue
            rid=x["rsid"]; rr=byrs.get(rid); hr=hby.get(rid); er=eligible.get(rid); gr=gby.get(rid)
            if None in (rr,hr,er,gr): errors.append(f"{did}:{rid}: upstream foreign-key closure failure"); continue
            if not same_variant(x,rr): errors.append(f"{did}:{rid}: reference coordinate/allele closure failure")
            if not eqtl_closes(x,er): errors.append(f"{did}:{rid}: raw eQTL field closure failure")
            if not gwas_closes(x,gr): errors.append(f"{did}:{rid}: normalized GWAS/direction-transform closure failure")
            for key in ("internal_id","genotype_id","ldscore"):
                if str(x[key])!=str(rr[key]): errors.append(f"{did}:{rid}: {key} closure failure")
            if any(str(x[k])!=str(hr[k]) for k in x if k in hr): errors.append(f"{did}:{rid}: harmonized/input closure failure")
        vcor=od/"prune-validation.vcor"; vlog=od/"plink-prune-validation.log"
        if not vcor.is_file() or not vlog.is_file(): errors.append(f"{did}: required LD recomputation artifact missing")
        else:
            vlines=vcor.read_text().splitlines()
            if not vlines or vlines[0].split("\t")!=["#ID_A","ID_B","UNPHASED_R2"]: errors.append(f"{did}: invalid LD recomputation header")
            for line in vlines[1:]:
                try:
                    if float(line.split("\t")[2])>0.1: errors.append(f"{did}: r2 threshold violation")
                except (IndexError,ValueError): errors.append(f"{did}: malformed LD recomputation row")
            log=vlog.read_text(errors="replace")
            for token in ("--r2-unphased","--ld-window-kb 1000","--ld-window-r2 0.1000000001","SH2B3_EUR_hg19"):
                if token not in log: errors.append(f"{did}: LD validation log lacks {token}")
        if r["status"]=="success_inference_bonferroni_selection_dependent":
            rr=rows(od/"mrjti-result.tsv"); boot=rows(od/"bootstrap.tsv")
            if len(rr)!=1 or len(boot)!=int(rr[0].get("n_bootstrap",-1)): errors.append(f"{did}: result/bootstrap row closure failure")
            else:
                one=rr[0]; draws=[float(x["expression_beta"]) for x in boot]
                mapping=("standardized_expression_coefficient","original_lasso_coefficient","ci_low","ci_high","lambda","cv_mean_error","cv_error_se","bootstrap_zero_fraction","bootstrap_positive_fraction","bootstrap_negative_fraction")
                if any(not close(r[k],one[k]) for k in mapping): errors.append(f"{did}: final/per-outcome result mismatch")
                derived=(sum(draws)/len(draws),sum(x==0 for x in draws)/len(draws),sum(x>0 for x in draws)/len(draws),sum(x<0 for x in draws)/len(draws))
                if not all(close(a,b) for a,b in zip(derived,(one["standardized_expression_coefficient"],one["bootstrap_zero_fraction"],one["bootstrap_positive_fraction"],one["bootstrap_negative_fraction"]))): errors.append(f"{did}: bootstrap summary mismatch")
            if len(inp)<20: errors.append(f"{did}: Bonferroni success below 20-SNP execution threshold")
            if r.get("pvalue")!="NA" or r.get("q_mrjti")!="NA": errors.append(f"{did}: unsupported P/q field is populated")
            expected_sig="sig" if float(r["ci_low"])*float(r["ci_high"])>0 else "nonsig"
            if r.get("ci_significance")!=expected_sig or r.get("mrjti_supported")!=str(expected_sig=="sig").lower(): errors.append(f"{did}: native Bonferroni support mismatch")
            if selection.get("mode")=="explicit_v2_spredixcan" and r.get("direction_concordant")!=str(float(r["standardized_expression_coefficient"])*float(r["source_z"])>0).lower(): errors.append(f"{did}: direction concordance mismatch")
            required += [od/"mrjti-result.tsv",od/"bootstrap.tsv"]
        required += [od/x for x in ("harmonized.tsv","pruned-variants.tsv","mrjti-input.tsv","prune.prune.in","prune-validation.vcor","plink-prune-validation.log")]
    required += [d/"ld_reference/ld-reference-variant-map.tsv",d/"ld_reference/SH2B3_EUR_hg19_ldscore.score.ld",d/"mrjti-minimal-method.json",d/"mrjti-minimal-snp-flow.tsv",d/"mrjti-minimal.tsv"]
    errors += [f"missing scientific artifact: {p}" for p in required if not p.is_file()]
    science=scientific_artifacts(d); audit=audit_artifacts(d)
    return errors,{"raw_cis_eqtl":len(eq),"ld_reference_variants":len(ref),"outcomes":len(result),"scientific_artifacts":science,"scientific_bundle_sha256":artifact_digest(science),"audit_artifacts":audit,"audit_bundle_sha256":artifact_digest(audit)}

def compare_artifacts(left,right):
    keys=sorted(set(left)|set(right)); return [{"category":"method" if k.endswith(".json") else "scientific_table_or_array","artifact":k,"left":left.get(k),"right":right.get(k),"status":"match" if left.get(k)==right.get(k) else "different"} for k in keys]

def ensure_distinct(left,right):
    if Path(left).resolve()==Path(right).resolve(): raise ValueError("self-compare is not a repeatability test")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("directory",type=Path); ap.add_argument("--compare",type=Path); ap.add_argument("--report",type=Path); a=ap.parse_args()
    if a.compare:
        try: ensure_distinct(a.directory,a.compare)
        except ValueError as e: raise SystemExit(str(e))
    errors,stats=validate(a.directory)
    if a.compare:
        other_errors,other=validate(a.compare); errors += [f"comparison:{x}" for x in other_errors]
        differences=compare_artifacts(stats["scientific_artifacts"],other["scientific_artifacts"]); stats["scientific_comparison"]=differences
        if any(x["status"]!="match" for x in differences): errors.append("isolated rerun normalized scientific artifacts differ")
        stats["comparison_directory"]=str(a.compare.resolve()); stats["repeatability_proven"]=not errors
    else: stats["repeatability_proven"]=False
    report={"status":"pass" if not errors else "fail","errors":errors,**stats}
    if a.report: a.report.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True))
    if errors: raise SystemExit(1)
if __name__=="__main__": main()
