#!/usr/bin/env python3
"""Read-only preflight for the project MR-JTI workflow."""
import argparse, csv, gzip, json, os, sqlite3, subprocess
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--project-root",type=Path,default=Path(__file__).resolve().parents[2]); a=ap.parse_args()
    root=a.project_root.resolve(); errors=[]
    loader=root/"scripts/lib/load_config.sh"
    if not loader.is_file(): errors.append(f"missing config loader: {loader}"); env={}
    else:
        p=subprocess.run(["bash","-c","set -a; source scripts/lib/load_config.sh; env -0"],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if p.returncode: errors.append(p.stderr.decode(errors="replace")); env={}
        else: env=dict(x.split("=",1) for x in p.stdout.decode(errors="surrogateescape").split("\0") if "=" in x)
    # The S-MultiXcan scripts below lived in the run-mrjti-candidate-analysis
    # skill and were moved into this project directory on 2026-09-02. The former
    # requirement on scripts/01_pe_genetic_map/cross_tissue_core.py was dropped
    # at the same time: that module belongs to the disabled legacy
    # tissue-covariance-union path, is absent project-wide, and nothing imports it.
    m1=root/"scripts/01_pe_genetic_map"
    required_files=[root/"config/pe_gwas.tsv",root/"tools/smr/smr",root/"tools/MR-JTI/mr/MR-JTI.r",root/"tools/MetaXcan/software/SPrediXcan.py",root/"tools/MetaXcan/software/SMulTiXcan.py",m1/"prepare_gwas.py",m1/"mrjti_minimal.py",m1/"mrjti_minimal_runner.R",m1/"validate_mrjti_minimal.py",m1/"build_smultixcan_covariance.py",m1/"build_all_gene_smultixcan_covariance.py",m1/"run_all_gene_smultixcan.py",m1/"smultixcan_genome_core.py",m1/"validate_existing_all_gene_smultixcan.py",m1/"test_smultixcan_genome.py",m1/"test_smultixcan_reference_covariance.py"]
    for x in required_files:
        if not x.is_file(): errors.append(f"missing file: {x}")
    required_env=["GTEX_V8_BESD_DIR","JTI_MODEL_DIR","LD_REF_HG19_DIR","HG19_FASTA","DBSNP_HG38_VCF","HG19_VCF_DIR","PYTHON_METAXCAN","PLINK2","GCTA","RSCRIPT","PUBLIC_LDSC_DIR"]
    resolved={}
    for key in required_env:
        value=env.get(key); resolved[key]=value
        if not value or not Path(value).exists(): errors.append(f"missing configured path {key}: {value}")
    versions={}
    for key,args in {"PLINK2":["--version"],"GCTA":["--version"],"RSCRIPT":["--version"],"PYTHON_METAXCAN":["--version"]}.items():
        value=env.get(key)
        if not value or not os.access(value,os.X_OK): errors.append(f"configured executable is not executable {key}: {value}"); continue
        p=subprocess.run([value,*args],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        lines=[x.strip() for x in p.stdout.splitlines() if x.strip() and set(x.strip())!={"*"}]
        versions[key]=next((x for x in lines if "version" in x.lower()),lines[0] if lines else f"exit={p.returncode}")
        # GCTA 1.94.1 prints its version and then exits 1 because --out is
        # absent; the banner itself is a successful executable/version probe.
        if p.returncode and not (key=="GCTA" and any("version v" in x.lower() for x in lines)): errors.append(f"cannot query {key} version: exit {p.returncode}")
    ld_dir=Path(env["LD_REF_HG19_DIR"]) if env.get("LD_REF_HG19_DIR") else None
    if ld_dir:
        for chrom in range(1,23):
            for suffix in ("pgen","pvar","psam"):
                p=ld_dir/f"eur_chr{chrom}.{suffix}"
                if not p.is_file(): errors.append(f"missing LD reference member: {p}")
    fasta=Path(env["HG19_FASTA"]) if env.get("HG19_FASTA") else None
    if fasta and not Path(str(fasta)+".fai").is_file(): errors.append(f"missing FASTA index: {fasta}.fai")
    dbsnp=Path(env["DBSNP_HG38_VCF"]) if env.get("DBSNP_HG38_VCF") else None
    if dbsnp and not any(Path(str(dbsnp)+x).is_file() for x in (".tbi",".csi")): errors.append(f"missing VCF index: {dbsnp}")
    bridge_dir=Path(env["HG19_VCF_DIR"]) if env.get("HG19_VCF_DIR") else None
    if bridge_dir:
        for chrom in range(1,23):
            bridges=sorted(bridge_dir.glob(f"ALL.chr{chrom}.*.vcf.gz"))
            if len(bridges)!=1: errors.append(f"expected one hg19 bridge VCF for chr{chrom}, found {len(bridges)}"); continue
            if not any(Path(str(bridges[0])+x).is_file() for x in (".tbi",".csi")): errors.append(f"missing bridge VCF index: {bridges[0]}")
    current=root/"results/01_pe_genetic_map/CURRENT"
    if not current.exists(): errors.append(f"missing CURRENT: {current}")
    model_dir=Path(env["JTI_MODEL_DIR"]) if env.get("JTI_MODEL_DIR") else None
    dbs=sorted(model_dir.glob("JTI_*.db")) if model_dir and model_dir.is_dir() else []
    covs=sorted(model_dir.glob("JTI_*.txt.gz")) if model_dir and model_dir.is_dir() else []
    db_stems={p.stem.removeprefix("JTI_") for p in dbs}; cov_stems={p.name.removeprefix("JTI_").removesuffix(".txt.gz") for p in covs}
    model_counts={"db":len(dbs),"covariance":len(covs),"paired_tissues":len(db_stems&cov_stems)}
    if not dbs or db_stems!=cov_stems: errors.append(f"unpaired JTI tissue stems: db_only={sorted(db_stems-cov_stems)}, covariance_only={sorted(cov_stems-db_stems)}")
    required_schema={"weights":{"rsid","gene","weight","ref_allele","eff_allele"},"extra":{"gene","genename","pred.perf.R2","n.snps.in.model","pred.perf.pval","pred.perf.qval"}}
    for db in dbs:
        try:
            with sqlite3.connect(f"file:{db}?mode=ro",uri=True) as con:
                for table,need in required_schema.items():
                    got={r[1] for r in con.execute(f"pragma table_info('{table}')")}
                    if not need<=got: errors.append(f"model schema missing {db.name}:{table}:{sorted(need-got)}")
        except sqlite3.Error as e: errors.append(f"model sqlite unreadable {db}: {e}")
    for cov in covs:
        try:
            with gzip.open(cov,"rt",newline="") as h: header=next(csv.reader(h,delimiter="\t"),[])
            if header != ["GENE","RSID1","RSID2","VALUE"]: errors.append(f"covariance header mismatch {cov.name}: {header}")
        except (OSError,EOFError) as e: errors.append(f"covariance unreadable {cov}: {e}")
    rscript=env.get("RSCRIPT")
    if rscript:
        p=subprocess.run([rscript,"-e",'quit(status=ifelse(all(vapply(c("optparse","glmnet","HDCI"),requireNamespace,logical(1),quietly=TRUE)),0,1))'])
        if p.returncode: errors.append("missing required R package among optparse, glmnet, HDCI")
    report={"status":"pass" if not errors else "fail","validation_level":"path_executable_pairing_and_basic_schema_only","scientific_validation":"not_performed","project_root":str(root),"current":os.path.realpath(current) if current.exists() else None,"configured_paths":resolved,"software_versions":versions,"jti_model_counts":model_counts,"errors":errors}
    print(json.dumps(report,indent=2,sort_keys=True))
    raise SystemExit(0 if not errors else 1)
if __name__=="__main__": main()
