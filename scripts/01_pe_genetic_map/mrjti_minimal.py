#!/usr/bin/env python3
"""Auditable minimal MR-JTI run for Whole_Blood / SH2B3 across four outcomes."""
from __future__ import annotations

import argparse, csv, gzip, hashlib, json, math, os, re, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENE = "ENSG00000111252"
GENE_NAME = "SH2B3"
TISSUE = "Whole_Blood"
PRIMARY = "FIGSHARE_22680904_v2"
SEED = 20260714
NA = "NA"
RSID = re.compile(r"^rs[0-9]+$")
SEQUENCE_ALLELE = re.compile(r"^[ACGT]+$")
COMP = str.maketrans("ACGT", "TGCA")

def load_env():
    p=subprocess.run(["bash","-c","set -a; source scripts/lib/load_config.sh; env -0"],cwd=ROOT,check=True,stdout=subprocess.PIPE)
    return dict(x.split("=",1) for x in p.stdout.decode(errors="surrogateescape").split("\0") if "=" in x)

def run(cmd, log, cwd=ROOT):
    p=subprocess.run(list(map(str,cmd)),cwd=cwd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    log.parent.mkdir(parents=True,exist_ok=True); log.write_text(p.stdout)
    if p.returncode: raise RuntimeError(f"command failed ({p.returncode}): {' '.join(map(str,cmd))}; see {log}")

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()

def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def read_tsv(path):
    op=gzip.open if str(path).endswith(".gz") else open
    with op(path,"rt",newline="") as f: yield from csv.DictReader(f,delimiter="\t")

def finite(*xs):
    try: return all(math.isfinite(float(x)) for x in xs)
    except (TypeError,ValueError): return False

SPREDIXCAN_METHOD = "MetaXcan-0.7.5/JTI-v8/zscore-beta-over-se/phi-calibration-false"

def p_z_qc(beta,se,pvalue,relative_tolerance=.05,absolute_tolerance=1e-12):
    """Audit reported P against beta/SE without making P authoritative."""
    beta,se,pvalue=map(float,(beta,se,pvalue))
    if not finite(beta,se,pvalue) or se<=0 or not 0<=pvalue<=1: return "invalid"
    if pvalue==0: return "p_zero"
    implied=math.erfc(abs(beta/se)/math.sqrt(2))
    if abs(pvalue-implied)<=absolute_tolerance or abs(pvalue-implied)/max(pvalue,implied,1e-300)<=relative_tolerance: return "concordant_or_rounded"
    return "discrepant"

def resolve_candidate_source(path,dataset_ids):
    """Return validated v2 source rows, or explicit direction-free history."""
    if path is None:
        return ({did:None for did in dataset_ids},{"mode":"historical_candidate_identity_only","upstream_v2_complete":False,"selection_dataset":PRIMARY,"candidate_rule":"SH2B3/Whole_Blood frozen from v1 legacy upstream result; not revalidated"})
    path=Path(path); rows=list(read_tsv(path)); selected=[]
    for r in rows:
        if r.get("tissue")==TISSUE and r.get("gene_id","").split(".")[0]==GENE: selected.append(r)
    by=defaultdict(list)
    for r in selected: by[r.get("dataset_id")].append(r)
    if set(by)!=set(dataset_ids) or any(len(v)!=1 for v in by.values()): raise RuntimeError("v2 S-PrediXcan candidate rows do not close one-to-one across outcomes")
    out={}
    for did in dataset_ids:
        r=by[did][0]
        if r.get("method_version")!=SPREDIXCAN_METHOD: raise RuntimeError(f"non-v2 S-PrediXcan method for {did}")
        if r.get("run_status")!="success" or not finite(r.get("zscore")): raise RuntimeError(f"invalid v2 S-PrediXcan row for {did}")
        out[did]=r
    primary=out[PRIMARY]
    if primary.get("bonferroni_significant")!="true": raise RuntimeError("v2 primary candidate does not satisfy frozen Bonferroni rule")
    family=primary.get("bonferroni_family_size") or primary.get("family_size")
    if family in (None,"",NA): raise RuntimeError("v2 candidate source lacks Bonferroni family size")
    return out,{"mode":"explicit_v2_spredixcan","upstream_v2_complete":True,"selection_dataset":PRIMARY,"candidate_rule":"tissue Bonferroni p<=0.05/family_size","family_size":int(float(family)),"path":str(path.resolve()),"sha256":sha(path),"method_version":SPREDIXCAN_METHOD}

def is_snp(a,b): return len(a)==len(b)==1 and a in "ACGT" and b in "ACGT"
def palindromic(a,b): return is_snp(a,b) and {a,b} in ({"A","T"},{"C","G"})
def complement(a): return a.translate(COMP)[::-1]

def pair_match(a1,a2,b1,b2,allow_complement=True):
    a1,a2,b1,b2=(x.upper() for x in (a1,a2,b1,b2))
    if (b1,b2)==(a1,a2): return 1,"direct"
    if (b1,b2)==(a2,a1): return -1,"swapped"
    if allow_complement and is_snp(a1,a2) and not palindromic(a1,a2):
        if (complement(b1),complement(b2))==(a1,a2): return 1,"complement"
        if (complement(b1),complement(b2))==(a2,a1): return -1,"complement_swapped"
    return 0,"allele_conflict"

def filter_eqtl_rows(raw):
    candidates=[]; counts=Counter(raw=len(raw))
    for r in raw:
        rid=r["SNP"]
        if not RSID.fullmatch(rid): counts["excluded_noncanonical_rsid"]+=1; continue
        if not finite(r["b"],r["SE"],r["p"]) or float(r["SE"])<=0: counts["excluded_invalid_stat"]+=1; continue
        candidates.append(r)
    multiplicity=Counter(r["SNP"] for r in candidates)
    duplicate_ids={rid for rid,n in multiplicity.items() if n != 1}
    counts["excluded_duplicate_eqtl_rows"]=sum(multiplicity[rid] for rid in duplicate_ids)
    counts["excluded_duplicate_eqtl_groups"]=len(duplicate_ids)
    kept={r["SNP"]:r for r in candidates if multiplicity[r["SNP"]] == 1}
    counts["eligible_eqtl"]=len(kept)
    return kept,counts

def extract_eqtl(env,out,gene_id=GENE,tissue=TISSUE,artifact_name=GENE_NAME):
    """Extract one complete BESD cis-eQTL probe without significance filtering."""
    eqtl_dir=out/"eqtl" if (gene_id,tissue,artifact_name)==(GENE,TISSUE,GENE_NAME) else out/"eqtl"/tissue
    prefix=eqtl_dir/artifact_name; prefix.parent.mkdir(parents=True,exist_ok=True); txt=prefix.with_suffix(".txt")
    if not txt.is_file():
        run([ROOT/"tools/smr/smr","--beqtl-summary",Path(env["GTEX_V8_BESD_DIR"])/tissue/tissue,
             "--query","1","--probe",gene_id,"--out",prefix],out/"logs"/f"smr-{tissue}-{gene_id}.log")
    raw=list(read_tsv(txt)); kept,counts=filter_eqtl_rows(raw)
    return raw,kept,counts,txt

def reference_candidates(pvar,chrom,start,end):
    region=[]; ids=set(); positions=set(); counts=Counter()
    with pvar.open() as f:
        for line in f:
            if line.startswith("#"): continue
            x=line.rstrip().split("\t"); pos=int(x[1])
            if x[0]!=str(chrom) or not start<=pos<=end: continue
            counts["reference_region_raw"]+=1
            rid,ref,alt=x[2],x[3].upper(),x[4].upper()
            if not RSID.fullmatch(rid): counts["reference_excluded_noncanonical_rsid"]+=1; continue
            if "," in alt: counts["reference_excluded_multiallelic"]+=1; continue
            if not SEQUENCE_ALLELE.fullmatch(ref) or not SEQUENCE_ALLELE.fullmatch(alt):
                counts["reference_excluded_nonsequence_allele"]+=1; continue
            region.append({"chromosome":x[0],"position":pos,"rsid":rid,"ref":ref,"alt":alt})
            ids.add(rid); positions.add((x[0],pos))
    idn=Counter(); posn=Counter()
    with pvar.open() as f:
        for line in f:
            if line.startswith("#"): continue
            x=line.rstrip().split("\t"); key=(x[0],int(x[1]))
            if x[2] in ids: idn[x[2]]+=1
            if key in positions: posn[key]+=1
    keep=[]
    for r in region:
        if idn[r["rsid"]]!=1: counts["reference_excluded_nonunique_rsid"]+=1; continue
        if posn[(r["chromosome"],r["position"])]!=1: counts["reference_excluded_nonunique_position"]+=1; continue
        r["internal_id"]=f"{r['chromosome']}:{r['position']}:{r['ref']}:{r['alt']}"; keep.append(r)
    counts["reference_pre_genotype_qc"]=len(keep)
    return keep,counts

def build_reference(env,out,rows,chrom,artifact_name=GENE_NAME):
    refdir=out/"ld_reference"; refdir.mkdir(parents=True,exist_ok=True)
    extract=refdir/"reference-rsids.txt"; write_tsv(extract,["rsid"],rows)
    prefix=refdir/f"{artifact_name}_EUR_hg19"
    if not prefix.with_suffix(".bed").is_file():
        run([env["PLINK2"],"--pfile",Path(env["LD_REF_HG19_DIR"])/f"eur_chr{chrom}",
             "--extract",extract,"--min-alleles","2","--max-alleles","2","--maf","0.01","--geno","0.02",
             "--make-bed","--freq","--missing","variant-only","--out",prefix],out/"logs/plink-reference.log")
    original={(r["chromosome"],r["position"],r["ref"],r["alt"]):r for r in rows}
    original.update({(r["chromosome"],r["position"],r["alt"],r["ref"]):r for r in rows})
    freq={r["ID"]:r for r in read_tsv(prefix.with_suffix(".afreq"))}
    miss={r["ID"]:r for r in read_tsv(prefix.with_suffix(".vmiss"))}
    mapped=[]
    with prefix.with_suffix(".bim").open() as f:
        for line in f:
            c,i,cm,bp,a1,a2=line.rstrip().split("\t"); src=original.get((c,int(bp),a1,a2))
            if not src: raise RuntimeError(f"cannot close reference allele mapping for {i}")
            af=freq[i]; vm=miss[i]
            if is_snp(src["ref"],src["alt"]): variant_type="snp"
            elif len(src["ref"]) != len(src["alt"]): variant_type="indel"
            else: variant_type="mnv_or_complex_substitution"
            mapped.append(src|{"genotype_id":i,"maf":min(float(af["ALT_FREQS"]),1-float(af["ALT_FREQS"])),"missing_rate":float(vm["F_MISS"]),"variant_type":variant_type})
    score_prefix=refdir/f"{artifact_name}_EUR_hg19_ldscore"
    score_file=Path(str(score_prefix)+".score.ld")
    if not score_file.is_file():
        run([env["GCTA"],"--bfile",prefix,"--ld-score","--ld-wind","1000","--ld-rsq-cutoff","0.01","--out",score_prefix],out/"logs/gcta-ldscore.log")
    with score_file.open() as f:
        header=f.readline().split(); scores={x[0]:float(x[-1]) for line in f if (x:=line.split())}
    if len(scores)!=len(mapped): raise RuntimeError("LD-score/reference row closure failed")
    for r in mapped: r["ldscore"]=scores[r["genotype_id"]]
    write_tsv(refdir/"ld-reference-variant-map.tsv",["rsid","internal_id","genotype_id","chromosome","position","ref","alt","variant_type","maf","missing_rate","ldscore"],mapped)
    return prefix,mapped,score_file

def gwas_matches(path,eqtl):
    found=defaultdict(list)
    for r in read_tsv(path):
        if r.get("rsid") in eqtl: found[r["rsid"]].append(r)
    return found

def harmonize(dataset,gwas_path,eqtl,refrows):
    refs={r["rsid"]:r for r in refrows}; gw=gwas_matches(gwas_path,eqtl); counts=Counter(); rows=[]
    counts["eqtl_eligible"]=len(eqtl); counts["gwas_rsid_overlap"]=len(gw)
    for rid,e in eqtl.items():
        if rid not in gw: counts["excluded_missing_gwas"]+=1; continue
        if len(gw[rid])!=1: counts["excluded_duplicate_gwas_rsid"]+=1; continue
        if rid not in refs: counts["excluded_missing_ld_reference"]+=1; continue
        g=gw[rid][0]; rr=refs[rid]
        if str(g["chromosome"]).removeprefix("chr")!=e["Chr"] or int(float(g["position"]))!=int(e["BP"]) or rr["position"]!=int(e["BP"]): counts["excluded_position_conflict"]+=1; continue
        a1,a2=e["A1"].upper(),e["A2"].upper()
        if palindromic(a1,a2): counts["excluded_palindromic_snp"]+=1; continue
        refsign,refstatus=pair_match(a1,a2,rr["ref"],rr["alt"],True)
        if not refsign: counts["excluded_reference_allele_conflict"]+=1; continue
        if not finite(g["beta"],g["standard_error"],g["pvalue"]) or float(g["standard_error"])<=0: counts["excluded_invalid_gwas_stat"]+=1; continue
        sign,status=pair_match(a1,a2,g["effect_allele"],g["non_effect_allele"],True)
        if not sign: counts["excluded_gwas_allele_conflict"]+=1; continue
        rows.append({"rsid":rid,"internal_id":rr["internal_id"],"genotype_id":rr["genotype_id"],"effect_allele":a1,"other_allele":a2,"chromosome":e["Chr"],"position":e["BP"],"variant_type":rr["variant_type"],"harmonization":status,"reference_match":refstatus,"ldscore":rr["ldscore"],"eqtl_beta":float(e["b"]),"eqtl_se":float(e["SE"]),"eqtl_p":float(e["p"]),"gwas_beta":sign*float(g["beta"]),"gwas_se":float(g["standard_error"]),"gwas_p":float(g["pvalue"]),"gwas_source_effect_allele":g["effect_allele"],"gwas_source_other_allele":g["non_effect_allele"],"gwas_source_beta":float(g["beta"]),"gwas_source_se":float(g["standard_error"]),"gwas_source_p":float(g["pvalue"]),"gwas_alignment_sign":sign})
    rows.sort(key=lambda r:(int(r["chromosome"]),int(r["position"]),r["internal_id"])); counts["harmonized"]=len(rows)
    return rows,counts

def prune(env,out,dataset,refprefix,rows):
    d=out/"outcomes"/dataset; d.mkdir(parents=True,exist_ok=True)
    write_tsv(d/"harmonized.tsv",list(rows[0]) if rows else ["rsid"],rows)
    write_tsv(d/"target-genotype-ids.txt",["genotype_id"],rows)
    pre=d/"prune"
    run([env["PLINK2"],"--bfile",refprefix,"--extract",d/"target-genotype-ids.txt","--indep-pairwise","1000kb","1","0.1","--out",pre],d/"plink-prune.log")
    kept={x.strip() for x in Path(str(pre)+".prune.in").read_text().splitlines() if x.strip()}
    pruned=[r for r in rows if r["genotype_id"] in kept]
    check=d/"prune-validation"
    run([env["PLINK2"],"--bfile",refprefix,"--extract",Path(str(pre)+".prune.in"),
         "--r2-unphased","cols=id","--ld-window-kb","1000","--ld-window-r2","0.1000000001","--out",check],d/"plink-prune-validation.log")
    vcor=check.with_suffix(".vcor")
    if vcor.is_file() and sum(1 for _ in vcor.open())>1: raise RuntimeError(f"LD pruning threshold validation failed for {dataset}")
    fields=["rsid","internal_id","genotype_id","chromosome","position","variant_type","effect_allele","other_allele","ldscore","eqtl_beta","eqtl_se","eqtl_p","gwas_beta","gwas_se","gwas_p","gwas_source_effect_allele","gwas_source_other_allele","gwas_source_beta","gwas_source_se","gwas_source_p","gwas_alignment_sign"]
    write_tsv(d/"mrjti-input.tsv",fields,pruned)
    write_tsv(d/"pruned-variants.tsv",list(pruned[0]) if pruned else fields,pruned)
    return d,pruned

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--analysis-id",required=True); ap.add_argument("--output",type=Path); ap.add_argument("--spredixcan-table",type=Path); ap.add_argument("--n-genes",type=int,default=1); a=ap.parse_args()
    if a.n_genes < 1: raise SystemExit("--n-genes must be a positive planned per-GWAS family size")
    env=load_env(); aid=a.analysis_id
    out=a.output or ROOT/"results/01_pe_genetic_map/minimal-tests"/f"{aid}-SH2B3-Whole_Blood-v2"
    out.mkdir(parents=True,exist_ok=True)
    raw,eqtl,eqc,eqfile=extract_eqtl(env,out); probe=int(raw[0]["Probe_bp"]); chrom=int(raw[0]["Probe_Chr"])
    pvar=Path(env["LD_REF_HG19_DIR"])/f"eur_chr{chrom}.pvar"
    candidates,refc=reference_candidates(pvar,chrom,max(1,probe-2_000_000),probe+2_000_000)
    refprefix,refrows,scorefile=build_reference(env,out,candidates,chrom)
    datasets=[]
    with (ROOT/"config/pe_gwas.tsv").open() as f: datasets=list(csv.DictReader(f,delimiter="\t"))
    full=[]; flow=[]; source_rows,source_provenance=resolve_candidate_source(a.spredixcan_table,[x["dataset_id"] for x in datasets]); gwas_inputs={}
    for ds in datasets:
        did=ds["dataset_id"]; gwas=ROOT/"work/01_pe_genetic_map"/aid/"M1.1_prepare_gwas"/did/"metaxcan.tsv.gz"
        gwas_inputs[did]={"path":str(gwas.resolve()),"sha256":sha(gwas)}
        rows,hc=harmonize(did,gwas,eqtl,refrows); d,pruned=prune(env,out,did,refprefix,rows)
        status="ineligible_lt20" if len(pruned)<20 else "success_inference_bonferroni_selection_dependent"; source=source_rows[did]; source_z=source["zscore"] if source else NA
        relation="same_outcome_selection_and_analysis" if did==PRIMARY else "selected_in_primary_analyzed_in_other_outcome"
        result={"dataset_id":did,"tissue":TISSUE,"gene_id":GENE,"gene_name":GENE_NAME,"status":status,"n_snps":len(pruned),"standardized_expression_coefficient":NA,"original_lasso_coefficient":NA,"beta_scale":"standardized_across_retained_variants","ci_low":NA,"ci_high":NA,"ci_method":NA,"ci_significance":NA,"mrjti_family_size":a.n_genes,"mrjti_family_alpha":0.05/a.n_genes,"lambda":NA,"cv_mean_error":NA,"cv_error_se":NA,"bootstrap_zero_fraction":NA,"bootstrap_positive_fraction":NA,"bootstrap_negative_fraction":NA,"pvalue":NA,"q_mrjti":NA,"source_z":source_z,"source_z_mode":source_provenance["mode"],"direction_concordant":NA,"inference_status":"not_run","selection_dataset":PRIMARY,"selection_analysis_relation":relation,"mrjti_supported":NA,"test_scope":"minimal-test"}
        if len(pruned)>=20:
            run([env["RSCRIPT"],ROOT/"scripts/01_pe_genetic_map/mrjti_minimal_runner.R","--df-path",d/"mrjti-input.tsv","--vendor-script",ROOT/"tools/MR-JTI/mr/MR-JTI.r","--result-path",d/"mrjti-result.tsv","--bootstrap-path",d/"bootstrap.tsv","--n-folds","5","--n-bootstrap","500","--n-genes",str(a.n_genes),"--min-snps","20","--seed",str(SEED)],d/"mrjti.log")
            rr=next(read_tsv(d/"mrjti-result.tsv")); beta=float(rr["standardized_expression_coefficient"]); lo=float(rr["ci_low"]); hi=float(rr["ci_high"])
            concordant=NA if source is None else str(beta*float(source_z)>0).lower()
            result.update({"standardized_expression_coefficient":beta,"original_lasso_coefficient":rr["original_lasso_coefficient"],"ci_low":lo,"ci_high":hi,"ci_method":rr["ci_method"],"ci_significance":rr["ci_significance"],"mrjti_supported":str(rr["ci_significance"]=="sig").lower(),"lambda":rr["lambda"],"cv_mean_error":rr["cv_mean_error"],"cv_error_se":rr["cv_error_se"],"bootstrap_zero_fraction":rr["bootstrap_zero_fraction"],"bootstrap_positive_fraction":rr["bootstrap_positive_fraction"],"bootstrap_negative_fraction":rr["bootstrap_negative_fraction"],"direction_concordant":concordant,"inference_status":"native_bonferroni_ci_selection_dependent","status":"success_inference_bonferroni_selection_dependent"})
        full.append(result)
        merged=eqc+refc+hc; merged["reference_post_genotype_qc"]=len(refrows); merged["pruned"]=len(pruned)
        flow.append({"dataset_id":did,**merged})
    flow_fields=["dataset_id"]+sorted(set().union(*(x.keys() for x in flow))-{"dataset_id"})
    flow=[{k:(r.get(k,0) if k!="dataset_id" else r[k]) for k in flow_fields} for r in flow]
    write_tsv(out/"mrjti-minimal-snp-flow.tsv",flow_fields,flow)
    write_tsv(out/"mrjti-minimal.tsv",list(full[0]),full)
    method={"analysis_id":aid,"scope":"minimal-test","gene_id":GENE,"gene_name":GENE_NAME,"tissue":TISSUE,"probe_bp":probe,"cis_definition":"GTEx BESD TSS +/- 1 Mb","ld_reference_window":"TSS +/- 2 Mb (candidate-specific; never reused for another gene)","ld_reference":"1000 Genomes EUR hg19","reference_qc":{"canonical_rsid":True,"chromosome_wide_unique_rsid":True,"unique_position":True,"alleles":2,"sequence_alleles_only":True,"maf_min":.01,"missingness_max":.02,"biallelic_sequence_indels_and_mnvs":True},"ld_score":{"software":"GCTA","window_kb":1000,"r2_cutoff":.01,"computed_before_pruning":True},"pruning":{"window_kb":1000,"step":1,"r2":.1},"mrjti":{"folds":5,"bootstrap":500,"seed":SEED,"n_genes":a.n_genes,"coefficient_scale":"standardized_across_retained_variants","ci":"vendored basic residual-bootstrap; native Bonferroni alpha=0.05/n_genes","pvalue":None,"qvalue":None,"support_call":"sig/nonsig conditional on selected family","minimum_snps_is_execution_threshold_only":True},"selection":source_provenance|{"post_selection_inference":"Bonferroni-corrected within the selected family; no independent causal-confirmation claim"},"inputs":{"eqtl":{"path":str(eqfile.resolve()),"sha256":sha(eqfile)},"pvar":{"path":str(pvar.resolve()),"sha256":sha(pvar)},"ld_score":{"path":str(scorefile.resolve()),"sha256":sha(scorefile)},"reference_bed":{"path":str(refprefix.with_suffix('.bed').resolve()),"sha256":sha(refprefix.with_suffix('.bed'))},"reference_bim":{"path":str(refprefix.with_suffix('.bim').resolve()),"sha256":sha(refprefix.with_suffix('.bim'))},"reference_fam":{"path":str(refprefix.with_suffix('.fam').resolve()),"sha256":sha(refprefix.with_suffix('.fam'))},"normalized_gwas":gwas_inputs},"formal_publication_modified":False}
    (out/"mrjti-minimal-method.json").write_text(json.dumps(method,indent=2,sort_keys=True)+"\n")
    print(out)

if __name__=="__main__": main()
