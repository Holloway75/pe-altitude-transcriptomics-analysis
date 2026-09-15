#!/usr/bin/env python3
"""Governed M1.1--M1.8 executor for the frozen Module 1 protocol.

The executor is intentionally resumable only by a new attempt.  It never
publishes; module01.py remains the sole validator/publisher.
"""
from __future__ import annotations
import argparse, csv, datetime as dt, gzip, hashlib, heapq, json, math, os, signal, sqlite3, subprocess, sys, tempfile
import numpy as np
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent)); import module01 as gov
import prepare_gwas as prepmod
import stability_summary as stability_core
NA="NA"; SCHEMA="M1.4"

def digest_bytes(b): return hashlib.sha256(b).hexdigest()
def digest_file(p): return gov.sha256(Path(p))
def utc(): return dt.datetime.now(dt.timezone.utc).isoformat()
def bh(values):
    out=[None]*len(values); order=sorted(range(len(values)),key=lambda i:(values[i],i)); q=1.0
    for rank_i in range(len(order)-1,-1,-1):
        i=order[rank_i]; q=min(q,values[i]*len(values)/(rank_i+1)); out[i]=q
    return out
def fmt(x):
    if x is None: return NA
    if isinstance(x,bool): return str(x).lower()
    return f"{x:.17g}" if isinstance(x,float) else str(x)
def audit_p_z(beta,se,pvalue,relative_tolerance=.05,absolute_tolerance=1e-12):
    """P is non-authoritative; audit reasonable rounding against beta/SE Z."""
    try: beta,se,pvalue=map(float,(beta,se,pvalue))
    except (TypeError,ValueError): return "invalid"
    if not all(math.isfinite(x) for x in (beta,se,pvalue)) or se<=0 or not 0<=pvalue<=1: return "invalid"
    if pvalue==0: return "p_zero"
    implied=math.erfc(abs(beta/se)/math.sqrt(2))
    if abs(pvalue-implied)<=absolute_tolerance or abs(pvalue-implied)/max(pvalue,implied,1e-300)<=relative_tolerance: return "concordant_or_rounded"
    return "discrepant"
def write_tsv(path, fields, rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    with tmp.open("w",newline="") as h:
        w=csv.DictWriter(h,fields,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader()
        for r in rows: w.writerow({k:fmt(r.get(k)) for k in fields})
    tmp.replace(path)
def iter_tsv(path):
    """Yield TSV/CSV records without materialising an entire association file."""
    op=gzip.open if str(path).endswith(".gz") else open
    with op(path,"rt",newline="") as h:
        sample=h.read(8192); h.seek(0)
        lines=sample.splitlines()
        delim="," if lines and lines[0].count(",")>lines[0].count("\t") else "\t"
        yield from csv.DictReader(h,delimiter=delim)
def read_tsv(path):
    """Compatibility helper for the small metadata/result files only."""
    return list(iter_tsv(path))

class Run:
    def __init__(self, analysis_id):
        if not analysis_id or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for c in analysis_id): gov.die("invalid analysis_id")
        self.rid=analysis_id
        self.work=ROOT/"work/01_pe_genetic_map"/self.rid; self.stage=ROOT/"results/01_pe_genetic_map/analyses"/self.rid/".staging"
        self.work.mkdir(parents=True,exist_ok=True); self.stage.mkdir(parents=True,exist_ok=True)
        self.ledger=self.work/"task-ledger.tsv"; self.env=gov.load_env(); self.cfg=json.loads((ROOT/"config/module01.json").read_text())
        if not self.ledger.exists():
            header=["analysis_id","task_id","attempt","requirement_class","input_sha256","command_sha256","start","end","exit_code","status","stdout_sha256","stderr_sha256","output_sha256"]
            self.ledger.write_text("\t".join(header)+"\n")
        self.source_sha=digest_bytes((ROOT/"config/pe_gwas.tsv").read_bytes()+gov.canonical(self.cfg))
        self.datasets={r["dataset_id"]:r for r in gov.read_tsv(ROOT/"config/pe_gwas.tsv")}
        self.models=sorted(Path(self.env["JTI_MODEL_DIR"]).glob("JTI_*.db"))
        if not self.models: gov.die("no registered JTI models")
        # This is a cache of verified *individual input digests*, not an
        # execution manifest.  It prevents each resumed task from rereading
        # immutable multi-GB VCF resources.  Any change in filesystem identity,
        # size, or nanosecond mtime forces a fresh SHA-256 calculation.
        self.digest_cache_path=self.work/"input-digest-cache.json"
        try: self.digest_cache=json.loads(self.digest_cache_path.read_text())
        except (FileNotFoundError,json.JSONDecodeError): self.digest_cache={}
    def input_digest(self, p):
        p=Path(p); st=p.stat(); key=str(p.resolve())
        fingerprint=[st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns]
        prior=self.digest_cache.get(key)
        if prior and prior.get("fingerprint")==fingerprint: return prior["sha256"]
        value=digest_file(p)
        self.digest_cache[key]={"fingerprint":fingerprint,"sha256":value}
        tmp=self.digest_cache_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.digest_cache,sort_keys=True,separators=(",",":"))+"\n")
        tmp.replace(self.digest_cache_path)
        return value
    def task(self, task_id, requirement, inputs, command, fn, resume_outputs=(), compatible_command_hashes=()):
        base=self.work/task_id; base.mkdir(parents=True,exist_ok=True)
        attempts=[int(p.name.split("-")[-1]) for p in base.glob("attempt-*") if p.name.split("-")[-1].isdigit()]
        n=max(attempts,default=0)+1; wd=base/f"attempt-{n:02d}"; wd.mkdir()
        ih=digest_bytes("\n".join(self.input_digest(p) for p in inputs if Path(p).is_file()).encode()); ch=digest_bytes(command.encode())
        # A restart may reuse only a ledger-confirmed, byte-present success with
        # identical frozen inputs and command.  This is not a silent skip: the
        # historical successful attempt remains the provenance for downstream
        # rows, while interrupted attempts remain append-only evidence.
        resume_outputs=[Path(p) for p in resume_outputs]
        if resume_outputs and all(p.is_file() and p.stat().st_size>0 for p in resume_outputs):
            # An interrupted later attempt must not hide an earlier verified
            # success.  Find the most recent success, rather than the last
            # terminal row of any kind, so resumption is monotonic.
            latest=None
            for row in gov.read_tsv(self.ledger):
                if row["task_id"]==task_id and row["status"]=="success": latest=row
            # ``compatible_command_hashes`` is a narrowly scoped migration
            # bridge for an old ledger label.  It is permitted only where the
            # complete scientific inputs (including the implementation) are
            # byte-identical and the prior outputs are present.  New work must
            # always match the current command fingerprint exactly.
            compatible=set(compatible_command_hashes)
            if latest and latest["status"]=="success" and latest["input_sha256"]==ih and latest["command_sha256"] in compatible|{ch}:
                # If a later interrupted/failed attempt exists, reusing the
                # older success without recording recovery would leave that
                # failure as the task's apparent terminal state at publish
                # time.  Append a new, digest-verified cache recovery attempt
                # so the ledger truthfully records which state resume chose.
                terminal_attempts=[int(row["attempt"]) for row in gov.read_tsv(self.ledger)
                                   if row["task_id"]==task_id and row["status"]!="running"]
                if terminal_attempts and max(terminal_attempts)>int(latest["attempt"]):
                    now=utc(); empty=digest_bytes(b"")
                    (wd/"stdout.txt").write_bytes(b""); (wd/"stderr.txt").write_bytes(b"")
                    oh=digest_bytes("\n".join(f"{p}\0{digest_file(p)}" for p in resume_outputs).encode())
                    row=[self.rid,task_id,n,requirement,ih,ch,now,now,0,"success",empty,empty,oh]
                    with self.ledger.open("a") as h:
                        h.write("\t".join(map(str,row))+"\n"); h.flush(); os.fsync(h.fileno())
                    return n,ih,resume_outputs
                wd.rmdir()
                return int(latest["attempt"]),ih,resume_outputs
        start=utc(); status="running"; code="NA"; outputs=[]
        def append(end,status,code,oh):
            stdout=wd/"stdout.txt"; stderr=wd/"stderr.txt"
            row=[self.rid,task_id,n,requirement,ih,ch,start,end,code,status,digest_file(stdout),digest_file(stderr),oh]
            with self.ledger.open("a") as h: h.write("\t".join(map(str,row))+"\n"); h.flush(); os.fsync(h.fileno())
        (wd/"stdout.txt").write_bytes(b""); (wd/"stderr.txt").write_bytes(b"")
        append("NA","running","NA",digest_bytes(b""))
        old={s:signal.getsignal(s) for s in (signal.SIGINT,signal.SIGTERM)}
        def interrupted(sig,_): raise KeyboardInterrupt(f"signal {sig}")
        for s in old: signal.signal(s,interrupted)
        try: outputs=list(fn(wd) or []); status="success"; code=0
        except BaseException as e:
            code=128+getattr(e,"signal",1) if isinstance(e,KeyboardInterrupt) else 1; status="technical_failed"
            with (wd/"stderr.txt").open("ab") as h: h.write((type(e).__name__+": "+str(e)).encode())
        finally:
            for s,h in old.items(): signal.signal(s,h)
        oh=digest_bytes("\n".join(f"{Path(p)}\0{digest_file(p)}" for p in outputs if Path(p).is_file()).encode())
        append(utc(),status,code,oh)
        if status!="success": raise RuntimeError(f"{task_id} failed: {(wd/'stderr.txt').read_text(errors='replace')}")
        return n,ih,outputs
    def base(self,dataset): return {"schema_version":SCHEMA,"analysis_id":self.rid,"source_provenance_sha256":self.source_sha,"dataset_id":dataset,"analysis_role":self.datasets[dataset]["analysis_role"]}
    def status_task(self,task_id,requirement,status,payload):
        if status not in gov.STATUS: raise ValueError(status)
        base=self.work/task_id; base.mkdir(parents=True,exist_ok=True); n=max([int(p.name.split("-")[-1]) for p in base.glob("attempt-*")],default=0)+1; wd=base/f"attempt-{n:02d}"; wd.mkdir(); out=wd/"report.json"; out.write_text(json.dumps(payload,sort_keys=True,indent=2)+"\n"); now=utc(); h=digest_file(out)
        prefix=[self.rid,task_id,n,requirement,h,digest_bytes(task_id.encode()),now]
        running=prefix+[NA,NA,"running",digest_bytes(b""),digest_bytes(b""),digest_bytes(b"")]
        terminal=prefix+[now,0,status,digest_bytes(b""),digest_bytes(b""),h]
        with self.ledger.open("a") as f:f.write("\t".join(map(str,running))+"\n"+"\t".join(map(str,terminal))+"\n"); f.flush(); os.fsync(f.fileno())
        return n,h

def model_meta(db):
    with sqlite3.connect(db) as c:
        entities=[r[0] for r in c.execute("select distinct gene from weights")]; assert_gene_entities([{"gene_id":g.split('.')[0],"model_gene_id":g} for g in entities])
        extra={r[0].split('.')[0]:r for r in c.execute('select gene,genename,"pred.perf.R2","n.snps.in.model","pred.perf.pval" from extra')}
        snps=defaultdict(set)
        for rs,g in c.execute("select rsid,gene from weights"): snps[g.split('.')[0]].add(rs)
    return extra,snps
def model_alleles(db):
    out={}; malformed_rsids=set()
    with sqlite3.connect(db) as c:
        for rs,ref,eff in c.execute("select distinct rsid,ref_allele,eff_allele from weights"):
            pair=(str(eff).upper(),str(ref).upper())
            # Indels are valid here (e.g. ``AC``/``A``); only a token that is
            # not composed entirely of DNA bases is malformed.  If an rsID has
            # even one malformed representation, omit it for every gene rather
            # than silently assigning another gene's alleles to that weight.
            if any(not allele or set(allele)-set("ACGT") for allele in pair):
                malformed_rsids.add(rs); continue
            if rs in out and out[rs]!=pair: raise RuntimeError("model_allele_collision:"+rs)
            out[rs]=pair
    for rs in malformed_rsids: out.pop(rs,None)
    return out

def indexed_model_variants(gwas_index, model_db, model_gene_id):
    """Load anchors only for a selected candidate, never for every gene row."""
    with sqlite3.connect(model_db) as model, sqlite3.connect(gwas_index) as gwas:
        rsids=[r[0] for r in model.execute("select distinct rsid from weights where gene=?",(model_gene_id,))]
        out=[]
        for offset in range(0,len(rsids),900):
            chunk=rsids[offset:offset+900]
            out.extend(json.loads(payload) for _,payload in gwas.execute("select rsid,payload from gwas where rsid in ({})".format(",".join("?" for _ in chunk)),chunk))
    return out

def harmonize_for_model(source,out,qc,index,alleles,seed,per_status=100):
    """Stream a model-specific GWAS and retain a disk-backed rsID lookup.

    The old implementation retained the full input three times (row list, set,
    and dictionary), which is unsafe for the frozen GWAS files.  The SQLite
    index is intentionally part of the task output: it is the bounded-memory
    replacement for the former in-process ``gwas_by_rs`` dictionary.
    """
    counts=defaultdict(int); samples=defaultdict(list)
    out=Path(out); out.parent.mkdir(parents=True,exist_ok=True)
    index=Path(index); index.parent.mkdir(parents=True,exist_ok=True)
    if index.exists(): index.unlink()
    tmp=out.with_suffix(out.suffix+".tmp")
    opener=gzip.open if str(out).endswith(".gz") else open
    with sqlite3.connect(index) as db, opener(tmp,"wt",newline="") as handle:
        db.execute("create table gwas(rsid text primary key, payload text not null)")
        writer=csv.DictWriter(handle,prepmod.OUT,delimiter="\t",lineterminator="\n",extrasaction="ignore"); writer.writeheader()
        for r in iter_tsv(source):
            counts["n_input"]+=1
            pz=audit_p_z(r["beta"],r["standard_error"],r["pvalue"])
            counts[f"qc_pz_{pz}"]+=1
            if pz=="invalid": raise RuntimeError("invalid beta/SE/P after GWAS normalization")
            m=alleles.get(r["rsid"])
            if not m: counts["model_snp_absent"]+=1; continue
            got,status=prepmod.harmonize(r["effect_allele"],r["non_effect_allele"],r["beta"],r.get("effect_allele_frequency"),m[0],m[1])
            if got is None: counts[status]+=1; continue
            x=dict(r); x.update(source_ea=r["effect_allele"],source_nea=r["non_effect_allele"],source_beta=r["beta"],source_eaf=r.get("effect_allele_frequency",NA),effect_allele=m[0],non_effect_allele=m[1],beta=got["beta"],zscore=got["beta"]/float(r["standard_error"]),effect_allele_frequency=got["eaf"] if got["eaf"] is not None else NA,harmonization_status=status)
            writer.writerow({k:fmt(x.get(k)) for k in prepmod.OUT})
            db.execute("insert or replace into gwas values(?,?)",(x["rsid"],json.dumps(x,sort_keys=True,separators=(",",":"))))
            counts[status]+=1; counts["n_retained"]+=1
            # Fixed-seed bounded reservoir: retain only audit candidates, not
            # millions of otherwise identical full records.
            rank=int(hashlib.sha256(f"{seed}\0{status}\0{x['rsid']}\0{x['variant_id']}".encode()).hexdigest(),16)
            heap=samples[status]; item=(-rank,x)
            if len(heap)<per_status: heapq.heappush(heap,item)
            elif item[0]>heap[0][0]: heapq.heapreplace(heap,item)
        db.execute("create index gwas_rsid on gwas(rsid)")
        db.commit()
    excluded=sum(v for k,v in counts.items() if k not in {"n_input","n_retained","direct","swapped","strand_direct","strand_swapped"} and not k.startswith("qc_")); counts["closure_difference"]=counts["n_input"]-counts["n_retained"]-excluded
    if counts["closure_difference"]: raise RuntimeError("harmonization closure failed")
    prepmod.direction_audit([x for heap in samples.values() for _,x in heap],alleles,seed,per_status)
    tmp.replace(out)
    write_tsv(qc,["metric","value"],[{"metric":k,"value":v} for k,v in sorted(counts.items())])
    return counts
def tissue(db): return db.stem.removeprefix("JTI_")

def assert_gene_entities(rows):
    """Stop on standardized-ID collisions; do not select by P."""
    seen=defaultdict(set)
    for r in rows: seen[r["gene_id"].split(".")[0]].add(r.get("model_gene_id",r["gene_id"]))
    bad=sorted(g for g,v in seen.items() if len(v)>1)
    if bad: raise RuntimeError("gene_id_collision:"+",".join(bad))

def validate_match_rates(rows, primary, prepared, absolute=.20, relative=.80):
    by=defaultdict(dict)
    for r in rows: by[r["tissue"]][r["dataset_id"]]=float(r["match_rate"])
    for tis,vals in by.items():
        if any(v<absolute for v in vals.values()): raise RuntimeError("match_rate_alert:absolute:"+tis)
        controls=[vals[x] for x in prepared if x!=primary and x in vals]
        if len(controls)<2: raise RuntimeError("match_rate_alert:missing_prepared_controls:"+tis)
        controls.sort(); median=(controls[(len(controls)-1)//2]+controls[len(controls)//2])/2
        if vals.get(primary,-1)<relative*median: raise RuntimeError("match_rate_alert:relative:"+tis)

def mrjti_record(task_id,dataset,tissue,gene,source_z,seed,n_snps,status,**stats):
    raise RuntimeError("retired MR-JTI P/q/support record contract; use isolated descriptive minimal runner")
def execute_mrjti(run,full):
    raise RuntimeError(
        "M1.6 full MR-JTI is intentionally disabled: the former n_snps_used "
        "placeholder was scientifically invalid. Run mrjti_minimal.py first; "
        "formal publication requires a separately authorized full rerun."
    )

def bridge_cache_fingerprint(path):
    stat=Path(path).stat()
    return {"path":str(Path(path).resolve()),"device":stat.st_dev,"inode":stat.st_ino,"size":stat.st_size,"mtime_ns":stat.st_mtime_ns}

def ensure_hg19_bridge_cache(run,dataset,source,bridges):
    """Build/reuse a Stage-00-style bridge map outside result directories.

    It contains only 1kG candidate records. Final M1.1 FASTA and allele
    validation remains mandatory, so cache reuse never promotes an unchecked
    variant to a MetaXcan input.
    """
    builder=ROOT/"scripts/01_pe_genetic_map/build_hg19_bridge_index.py"
    cache_dir=Path(run.env["PE_GENETIC_INTERIM_DIR"])/"bridge-cache"/dataset
    index=cache_dir/"hg19-rsid-bridge.pkl"; report=cache_dir/"hg19-rsid-bridge.report.json"; manifest=cache_dir/"manifest.json"
    expected={"schema_version":"m1-hg19-rsid-bridge-cache-v1","dataset_id":dataset,"source_sha256":run.input_digest(source),"builder_sha256":digest_file(builder),"bridges":[bridge_cache_fingerprint(path) for path in bridges]}
    try:
        prior=json.loads(manifest.read_text())
    except (FileNotFoundError,json.JSONDecodeError): prior=None
    if index.is_file() and report.is_file() and prior==expected: return index,report,manifest
    cache_dir.mkdir(parents=True,exist_ok=True)
    command=[sys.executable,str(builder),"--gwas",str(source),"--bridge",*map(str,bridges),"--bcftools",run.env["BCFTOOLS"],"--out",str(index),"--report",str(report)]
    process=subprocess.run(command,cwd=ROOT,text=True,capture_output=True)
    (cache_dir/"build.stdout.log").write_text(process.stdout)
    (cache_dir/"build.stderr.log").write_text(process.stderr)
    if process.returncode or not index.is_file() or not report.is_file(): raise RuntimeError(f"bridge-cache build failed for {dataset}: "+process.stderr[-1000:])
    manifest_tmp=manifest.with_suffix(".tmp"); manifest_tmp.write_text(json.dumps(expected,sort_keys=True,indent=2)+"\n"); manifest_tmp.replace(manifest)
    return index,report,manifest

def execute_spredixcan(run:Run):
    full=[]; match=[]; prepared={}
    prep_script=ROOT/"scripts/01_pe_genetic_map/prepare_gwas.py"
    fasta=Path(run.env["HG19_FASTA"])
    if not fasta.is_file(): raise RuntimeError("registered hg19 FASTA provider incomplete")
    for did,d in run.datasets.items():
        source=Path(run.env["PE_GWAS_ROOT"])/d["hg19_relpath"]
        out=run.work/"M1.1_prepare_gwas"/did/"metaxcan.tsv.gz"; qc=run.work/"M1.1_prepare_gwas"/did/"qc.tsv"
        def prep(wd,source=source,out=out,qc=qc,source_build="hg19"):
            cmd=[sys.executable,str(prep_script),"--gwas",str(source),"--out",str(out),"--qc",str(qc),"--fasta",str(fasta),"--source-build",source_build]
            p=subprocess.run(cmd,cwd=ROOT,capture_output=True); (wd/"stdout.txt").write_bytes(p.stdout); (wd/"stderr.txt").write_bytes(p.stderr)
            if p.returncode: raise RuntimeError(f"prepare_gwas exit {p.returncode}")
            return [out,qc]
        run.task(f"M1.1_prepare_gwas__{did}","core_required",[source,prep_script,fasta],"prepare_gwas_v5_direct_hg19_fasta",prep,resume_outputs=[out,qc]); prepared[did]=out
    for did,base_gwas in prepared.items():
        outcome_rows=[]
        for db in run.models:
            covariance=db.with_suffix(".txt.gz")
            if not covariance.is_file(): raise RuntimeError(f"missing model covariance: {covariance}")
            tis=tissue(db); gwas=run.work/"M1.2_metaxcan_input"/did/f"{tis}.tsv.gz"; hqc=run.work/"M1.2_metaxcan_input"/did/f"{tis}.qc.tsv"; gwas_index=run.work/"M1.2_metaxcan_input"/did/f"{tis}.sqlite"; alleles=model_alleles(db)
            def harmonized(wd,base_gwas=base_gwas,gwas=gwas,hqc=hqc,gwas_index=gwas_index,alleles=alleles):
                harmonize_for_model(base_gwas,gwas,hqc,gwas_index,alleles,run.cfg["random_seed"],run.cfg["qc"]["direction_audit_per_status"])
                return [gwas,hqc,gwas_index]
            run.task(f"M1.2_harmonize__{did}__{tis}","core_required",[base_gwas,db],"harmonize_model_alleles",harmonized,resume_outputs=[gwas,hqc,gwas_index])
            raw=run.work/"M1.3_spredixcan"/did/f"{tis}.csv"
            cmd=[run.env["PYTHON_METAXCAN"],str(ROOT/"tools/MetaXcan/software/SPrediXcan.py"),"--model_db_path",str(db),"--covariance",str(covariance),"--gwas_file",str(gwas),"--snp_column","rsid","--effect_allele_column","effect_allele","--non_effect_allele_column","non_effect_allele","--zscore_column","zscore","--beta_column","beta","--se_column","standard_error","--pvalue_column","pvalue","--output_file",str(raw),"--remove_ens_version","--additional_output","--throw","--overwrite"]
            def one(wd,cmd=cmd,raw=raw):
                raw.parent.mkdir(parents=True,exist_ok=True); p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
                (wd/"stdout.txt").write_text(p.stdout); (wd/"stderr.txt").write_text(p.stderr)
                if p.returncode or not raw.is_file(): raise RuntimeError(f"SPrediXcan exit {p.returncode}")
                return [raw]
            attempt,ih,_=run.task(f"M1.3_spredixcan__{did}__{tis}","core_required",[gwas,db,covariance]," ".join(cmd),one,resume_outputs=[raw])
            meta,snps=model_meta(db); model_union=set().union(*snps.values()) if snps else set()
            model_sha=digest_file(db)
            with sqlite3.connect(gwas_index) as gdb:
                n_matched=gdb.execute("select count(*) from gwas").fetchone()[0]
                rate=n_matched/max(1,len(model_union))
                match.append({"dataset_id":did,"tissue":tis,"n_model_snps":len(model_union),"n_matched":n_matched,"match_rate":rate})
                entity_rows=[]
                for x in iter_tsv(raw):
                    original=x.get("gene") or x.get("gene_id") or ""; gid=original.split('.')[0]; z=float(x.get("zscore") or x.get("z_score")); p=float(x.get("pvalue")); e=meta.get(gid); entity_rows.append({"gene_id":gid,"model_gene_id":original})
                    used=int(float(x.get("n_snps_used") or x.get("n_snps_in_model") or (e[3] if e else 0)))
                    r=run.base(did)|{"tissue":tis,"gene_id":gid,"gene_name":e[1] if e and e[1] else NA,"zscore":z,"effect_size":float(x["effect_size"]) if x.get("effect_size") not in (None,"",NA) else NA,"pvalue":p,"n_snps_in_model":int(e[3]) if e else used,"n_snps_used":used,"pred_perf_r2":float(e[2]) if e and e[2] is not None else NA,"pred_perf_pvalue":float(e[4]) if e and e[4] is not None else NA,"direction_qc":"pass","run_status":"success","task_id":f"M1.3_spredixcan__{did}__{tis}","attempt":attempt,"input_sha256":ih,"method_version":"MetaXcan-0.7.5/JTI-v8/zscore-beta-over-se/phi-calibration-false","_used_variants":[],"_model_variants":[],"_model_sha":model_sha,"_model_db":str(db),"_gwas_index":str(gwas_index),"_model_gene_id":original}
                    outcome_rows.append(r)
            assert_gene_entities(entity_rows)
        qs=bh([r["pvalue"] for r in outcome_rows]); ng=len(outcome_rows)
        by_t=defaultdict(list)
        for i,r in enumerate(outcome_rows): r["q_gene_tissue"]=qs[i]; by_t[r["tissue"]].append(r)
        for rs in by_t.values():
            n=len(rs)
            for r in rs: r["p_bonf_tissue"]=.05/n; r["p_bonf_global"]=.05/ng; r["bonferroni_significant"]=r["pvalue"]<=.05/n
        full.extend(outcome_rows)
    validate_match_rates(match,"FIGSHARE_22680904_v2",list(prepared),run.cfg["qc"]["minimum_model_match_rate"],run.cfg["qc"]["primary_relative_match_rate"])
    return full,match

def load_validated_cross_tissue(run,full,result_root,validation_root):
    """Consume, but never recompute, skill-owned genome-wide M1.5 outputs."""
    result_root=Path(result_root).resolve(); validation_root=Path(validation_root).resolve()
    overall=validation_root/"validation-report.json"
    if not overall.is_file() or json.loads(overall.read_text()).get("status")!="pass":
        raise RuntimeError("skill M1.5 overall validation is absent or not pass")
    overall_report=json.loads(overall.read_text())
    validated={x["dataset_id"]:x for x in overall_report.get("datasets",[]) if x.get("status")=="pass"}
    if set(validated)!=set(run.datasets):
        raise RuntimeError("skill M1.5 validated dataset set differs from registry")
    skill_evidence={
        "schema_version":"m1.5-skill-ingestion-v1",
        "status":"pass",
        "source_validation":{"path":str(overall),"sha256":digest_file(overall)},
        "rule":overall_report.get("rule"),
        "covariance_resource":overall_report.get("covariance_resource"),
        "datasets":overall_report.get("datasets"),
        "validator":overall_report.get("validator"),
    }
    (run.stage/"m1.5-skill-validation.json").write_text(json.dumps(skill_evidence,sort_keys=True,indent=2)+"\n")
    allrows=[]
    tissue_index=defaultdict(list)
    for row in full:
        tissue_index[(row["dataset_id"],row["gene_id"])].append(row)
    for did in run.datasets:
        result=result_root/did/"smultixcan_all_genes.tsv"
        validation=validation_root/did/"validation-report.json"
        if not result.is_file() or not validation.is_file():
            raise RuntimeError(f"skill M1.5 artifacts absent for {did}")
        report=json.loads(validation.read_text())
        if report.get("status")!="pass" or report.get("dataset_id")!=did:
            raise RuntimeError(f"skill M1.5 dataset validation failed for {did}")
        recorded=report.get("source_artifacts",{}).get("enriched",{}).get("sha256")
        if recorded and recorded!=digest_file(result):
            raise RuntimeError(f"skill M1.5 enriched output hash mismatch for {did}")
        input_sha=digest_bytes((digest_file(result)+digest_file(validation)+digest_file(overall)).encode())
        rs=[]
        for x in iter_tsv(result):
            if x.get("classification")!="success": continue
            gid=x["gene"].split('.')[0]; tissue_rows=tissue_index.get((did,gid),[]); best=max(tissue_rows,key=lambda r:abs(r["zscore"])) if tissue_rows else None
            rs.append(run.base(did)|{"gene_id":gid,"gene_name":x.get("gene_name",NA),"method":"S-MultiXcan","statistic":-math.log10(float(x["pvalue"])),"signed_statistic":NA,"direction_available":False,"direction":NA,"pvalue":float(x["pvalue"]),"qvalue":float(x["qvalue_global_bh"]),"n_tissues":int(float(x["n"])),"n_independent_components":int(float(x["n_indep"])),"best_tissue":best["tissue"] if best else NA,"best_tissue_z":best["zscore"] if best else NA,"condition_number":float(x["condition_number"]),"model_status":"success","task_id":f"M1.5_skill_smultixcan__{did}","attempt":1,"input_sha256":input_sha,"method_version":"skill-SMulTiXcan-core/eigen_ratio-1e-6/outcome-subset-BH-v2"})
        if len(rs)!=int(report["counts"]["global_bh_family_size"]):
            raise RuntimeError(f"skill M1.5 BH family size mismatch for {did}")
        allrows.extend(rs)
        run.status_task(f"M1.5_skill_smultixcan__{did}","cross_tissue_conditional","success",{"rows":len(rs),"result":str(result),"validation":str(validation),"input_sha256":input_sha})
    return allrows

def annotation(run):
    rows=read_tsv(gov.expand(run.cfg["resources"]["gene_annotation"],run.env)); return {r["gene_id"].split('.')[0]:r for r in rows}
def interval_hits(chrom,pos,blocks):
    return [b for b in blocks if str(b["chromosome"]).removeprefix("chr")==str(chrom).removeprefix("chr") and int(b["start"])<=int(pos)<=int(b["end"])]

def candidate_anchors(candidate, used_snps, model_snps, ann, cross=False):
    """Frozen anchor priority. Cross-tissue candidates never use a fallback."""
    if used_snps: return used_snps,"used_snps",candidate["task_id"]
    if cross: return [],"anchor_unavailable",candidate["task_id"]
    if model_snps: return model_snps,"model_snps_only",candidate.get("_model_sha",candidate["input_sha256"])
    a=ann.get(candidate["gene_id"])
    if not a:return [],"anchor_unavailable",NA
    if a.get("tss") not in (None,"",NA): return [{"chromosome":a["chromosome"],"position":a["tss"],"variant_id":NA}],"tss_fallback",candidate["input_sha256"]
    if all(a.get(x) not in (None,"",NA) for x in ("chromosome","start","end")): return [{"chromosome":a["chromosome"],"position":a["start"],"end":a["end"],"variant_id":NA}],"gene_interval_fallback",candidate["input_sha256"]
    return [],"anchor_unavailable",candidate["input_sha256"]

def sensitivity_loci(candidates,ann,max_distance=1_000_000):
    points=sorted((ann[c["gene_id"]]["chromosome"],int(ann[c["gene_id"]]["tss"]),c["gene_id"],c["candidate_id"],c) for c in candidates if c["gene_id"] in ann and ann[c["gene_id"]].get("tss") not in (None,"",NA))
    out={}; cluster=[]; last=None
    def flush():
        if cluster:
            lid=f"TSS_1MB:{cluster[0][0]}:{cluster[0][1]}-{cluster[-1][1]}" 
            for *_,c in cluster: out[c["candidate_id"]]=lid
    for item in points:
        if last and (item[0]!=last[0] or item[1]-last[1]>max_distance): flush(); cluster=[]
        cluster.append(item); last=item
    flush(); return out

def spearman(xs,ys):
    if len(xs)<10:return None,"insufficient_n"
    def ranks(v):
        order=sorted(range(len(v)),key=lambda i:v[i]); r=[0.]*len(v); i=0
        while i<len(order):
            j=i+1
            while j<len(order) and v[order[j]]==v[order[i]]:j+=1
            mid=((i+1)+j)/2
            for k in order[i:j]:r[k]=mid
            i=j
        return r
    a,b=ranks(xs),ranks(ys); ma=sum(a)/len(a); mb=sum(b)/len(b); num=sum((x-ma)*(y-mb) for x,y in zip(a,b)); den=math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
    return (num/den if den else None),("success" if den else "constant_rank")
def z_direction(z): return "positive" if float(z)>0 else "negative" if float(z)<0 else NA
def loci(run,full,cross=None):
    ann=annotation(run); blocks=read_tsv(gov.expand(run.cfg["resources"]["european_ld_blocks"],run.env)); byc=defaultdict(list)
    complex_regions=read_tsv(gov.expand(run.cfg["resources"]["complex_regions"],run.env))
    for b in blocks: byc[b["chromosome"]].append(b)
    ldsha=digest_file(gov.expand(run.cfg["resources"]["european_ld_blocks"],run.env)); detail=[]
    for r in full:
        if r["run_status"]!="success" or r["q_gene_tissue"]>.05: continue
        a=ann.get(r["gene_id"]); cid=digest_bytes("\0".join(["tissue_global_bh_q05_v1",r["dataset_id"],r["tissue"],r["gene_id"]]).encode()); candidate=dict(r,candidate_id=cid)
        model_variants=r.get("_model_variants",[])
        if not model_variants and r.get("_gwas_index") and r.get("_model_db") and r.get("_model_gene_id"):
            model_variants=indexed_model_variants(r["_gwas_index"],r["_model_db"],r["_model_gene_id"])
        anchors,evidence,anchor_source=candidate_anchors(candidate,r.get("_used_variants",[]),model_variants,ann)
        hits=[]
        for v in anchors:
            hits.extend(interval_hits(v["chromosome"],v["position"],blocks))
        hits=list({b["locus_id"]:b for b in hits}.values())
        if not hits: hits=[{"locus_id":"UNRESOLVED:"+cid,"chromosome":a["chromosome"] if a else NA,"start":NA,"end":NA}]
        for b in hits:
            detail.append(run.base(r["dataset_id"])|{"statistical_layer":"tissue","candidate_id":cid,"locus_definition":"BerisaPickrell2016_EUR_hg19","locus_id":b["locus_id"],"candidate_rule_id":"tissue_global_bh_q05_v1","gene_id":r["gene_id"],"trigger_tissue":r["tissue"],"cross_tissue_method":NA,"candidate_direction":"risk_increasing" if r["zscore"]>0 else "risk_decreasing" if r["zscore"]<0 else NA,"candidate_statistic":r["zscore"],"candidate_p":r["pvalue"],"candidate_q":r["q_gene_tissue"],"anchor_evidence":evidence,"anchor_source":anchor_source,"anchor_count":len(anchors),"resolution_status":"resolved" if not b["locus_id"].startswith("UNRESOLVED") else "unresolved","complex_ld_region":False,"ld_resource_sha256":ldsha,"_chromosome":b.get("chromosome",NA),"_start":b.get("start",NA),"_end":b.get("end",NA)})
    for r in cross or []:
        if r["model_status"]!="success" or r["qvalue"]>.05: continue
        cid=digest_bytes("\0".join(["cross_tissue_bh_q05_v1",r["dataset_id"],r["gene_id"],r["method"]]).encode()); lid="UNRESOLVED:"+cid
        detail.append(run.base(r["dataset_id"])|{"statistical_layer":"cross_tissue","candidate_id":cid,"locus_definition":"BerisaPickrell2016_EUR_hg19","locus_id":lid,"candidate_rule_id":"cross_tissue_bh_q05_v1","gene_id":r["gene_id"],"trigger_tissue":NA,"cross_tissue_method":r["method"],"candidate_direction":NA,"candidate_statistic":r["statistic"],"candidate_p":r["pvalue"],"candidate_q":r["qvalue"],"anchor_evidence":"anchor_unavailable","anchor_source":r["input_sha256"],"anchor_count":0,"resolution_status":"unresolved","complex_ld_region":False,"ld_resource_sha256":ldsha,"_chromosome":NA,"_start":NA,"_end":NA})
    # Independent ±1 Mb TSS sensitivity definition; it never replaces LD blocks.
    representatives={r["candidate_id"]:r for r in detail}
    by_scope=defaultdict(list)
    for row in representatives.values(): by_scope[(row["dataset_id"],row["statistical_layer"])].append(row)
    for scoped in by_scope.values():
        sens=sensitivity_loci(scoped,ann)
        for cid,lid in sens.items():
            x=dict(representatives[cid]); a=ann[x["gene_id"]]; x.update(locus_definition="TSS_1MB_CHAIN_hg19",locus_id=lid,_chromosome=a["chromosome"],_start=max(0,int(a["tss"])-1_000_000),_end=int(a["tss"])+1_000_000); detail.append(x)
    grouped=defaultdict(list)
    for r in detail: grouped[(r["dataset_id"],r["statistical_layer"],r["locus_definition"],r["locus_id"])].append(r)
    summary=[]
    for (did,layer,definition,lid),rs in grouped.items():
        pos=sum(r["candidate_direction"]=="risk_increasing" for r in rs); neg=sum(r["candidate_direction"]=="risk_decreasing" for r in rs); miss=len(rs)-pos-neg
        direction="undetermined" if miss else "conflicting" if pos and neg else "concordant_risk_increasing" if pos else "concordant_risk_decreasing"
        chrom,start,end=rs[0]["_chromosome"],rs[0]["_start"],rs[0]["_end"]
        is_complex=False if NA in (chrom,start,end) else any(str(c["chromosome"]).removeprefix("chr")==str(chrom).removeprefix("chr") and int(c["start"])<=int(end) and int(c["end"])>=int(start) for c in complex_regions)
        for x in rs:x["complex_ld_region"]=is_complex
        summary.append(run.base(did)|{"statistical_layer":layer,"locus_definition":definition,"locus_id":lid,"chromosome":chrom,"locus_start":start,"locus_end":end,"locus_direction":direction,"n_candidates":len({r["candidate_id"] for r in rs}),"n_genes":len({r["gene_id"] for r in rs}),"n_direction_positive":pos,"n_direction_negative":neg,"n_direction_missing":miss,"complex_ld_region":is_complex,"resolution_status":rs[0]["resolution_status"],"ld_resource_sha256":ldsha})
    return detail,summary

def stability(run,full,locus_rows=None):
    by=defaultdict(dict)
    for r in full:
        if r["tissue"]=="Whole_Blood" and r["run_status"]=="success": by[r["dataset_id"]][r["gene_id"]]=r
    primary="FIGSHARE_22680904_v2"; rows=[]
    for comp in [x for x in run.datasets if x!=primary]:
        common=sorted(g for g in set(by[primary])&set(by[comp]) if by[primary][g].get("direction_qc")=="pass" and by[comp][g].get("direction_qc")=="pass"); n=len(common)
        def pct(vals, reverse=False):
            order=sorted(vals,key=lambda g:((-vals[g]) if reverse else vals[g],g)); ranks={}; i=0
            while i<len(order):
                j=i+1
                while j<len(order) and vals[order[j]]==vals[order[i]]: j+=1
                mid=((i+1)+j)/2; ranks.update({g:(mid-1)/(n-1) if n>1 else NA for g in order[i:j]}); i=j
            return ranks
        signed_a=pct({g:by[primary][g]["zscore"] for g in common}); signed_b=pct({g:by[comp][g]["zscore"] for g in common})
        abs_a=pct({g:abs(by[primary][g]["zscore"]) for g in common},True); abs_b=pct({g:abs(by[comp][g]["zscore"]) for g in common},True)
        for gid in common:
            a,b=by[primary][gid],by[comp][gid]
            for metric in ("signed_z","abs_z"):
                rows.append({"schema_version":SCHEMA,"analysis_id":run.rid,"source_provenance_sha256":run.source_sha,"comparison_level":"whole_blood_gene","analysis_set":"all","feature_id":gid,"primary_dataset":primary,"comparison_dataset":comp,"rank_metric":metric,"comparison_role":run.datasets[comp]["analysis_role"],"primary_stat":a["zscore"],"comparison_stat":b["zscore"],"primary_direction":z_direction(a["zscore"]),"comparison_direction":z_direction(b["zscore"]),"direction_concordant":NA if 0 in (a["zscore"],b["zscore"]) else a["zscore"]*b["zscore"]>0,"primary_rank_pct":signed_a[gid] if metric=="signed_z" else abs_a[gid],"comparison_rank_pct":signed_b[gid] if metric=="signed_z" else abs_b[gid],"primary_p":a["pvalue"],"comparison_p":b["pvalue"],"primary_q":a["q_gene_tissue"],"comparison_q":b["q_gene_tissue"],"primary_status":"success","comparison_status":"success","n_common":n,"sample_overlap":"known_or_possible","power_interpretation":"low_power_uncertain" if comp=="GCST90301704" else "descriptive_only"})
        # Project every primary candidate, retaining absent/non-success features.
        primary_candidates=sorted(g for g,r in by[primary].items() if r["q_gene_tissue"]<=.05)
        for gid in primary_candidates:
            a=by[primary][gid]; b=by[comp].get(gid)
            if b is None:
                cid=digest_bytes("\0".join(["tissue_global_bh_q05_v1",primary,"Whole_Blood",gid]).encode())
                rows.append({"schema_version":SCHEMA,"analysis_id":run.rid,"source_provenance_sha256":run.source_sha,"comparison_level":"whole_blood_candidate_projection","analysis_set":"all","feature_id":cid,"primary_dataset":primary,"comparison_dataset":comp,"rank_metric":"candidate_projection","comparison_role":run.datasets[comp]["analysis_role"],"primary_stat":a["zscore"],"comparison_stat":NA,"primary_direction":z_direction(a["zscore"]),"comparison_direction":NA,"direction_concordant":NA,"primary_rank_pct":NA,"comparison_rank_pct":NA,"primary_p":a["pvalue"],"comparison_p":NA,"primary_q":a["q_gene_tissue"],"comparison_q":NA,"primary_status":"success","comparison_status":"not_success","n_common":n,"sample_overlap":"known_or_possible","power_interpretation":"low_power_uncertain" if comp=="GCST90301704" else "descriptive_only"})
    # Locus features are compared only within the same layer/definition.
    loci=defaultdict(dict)
    for r in locus_rows or []: loci[r["dataset_id"]][f'{r["statistical_layer"]}|{r["locus_definition"]}|{r["locus_id"]}']=r
    for comp in [x for x in run.datasets if x!=primary]:
        for analysis_set,keep in (("all",lambda r:True),("exclude_complex_ld",lambda r:not r["complex_ld_region"])):
            common=sorted(k for k in set(loci[primary])&set(loci[comp]) if keep(loci[primary][k]) and keep(loci[comp][k])); n=len(common)
            for fid in common:
                a,b=loci[primary][fid],loci[comp][fid]
                rows.append({"schema_version":SCHEMA,"analysis_id":run.rid,"source_provenance_sha256":run.source_sha,"comparison_level":"locus","analysis_set":analysis_set,"feature_id":fid,"primary_dataset":primary,"comparison_dataset":comp,"rank_metric":"locus_direction","comparison_role":run.datasets[comp]["analysis_role"],"primary_stat":NA,"comparison_stat":NA,"primary_direction":a["locus_direction"],"comparison_direction":b["locus_direction"],"direction_concordant":a["locus_direction"]==b["locus_direction"],"primary_rank_pct":NA,"comparison_rank_pct":NA,"primary_p":NA,"comparison_p":NA,"primary_q":NA,"comparison_q":NA,"primary_status":"success","comparison_status":"success","n_common":n,"sample_overlap":"known_or_possible","power_interpretation":"low_power_uncertain" if comp=="GCST90301704" else "descriptive_only"})
    stats={(r["dataset_id"],r["tissue"],r["gene_id"]):r for r in full}; primary_candidates=[r for r in full if r["dataset_id"]==primary and r["q_gene_tissue"]<=.05]
    for comp in [x for x in run.datasets if x!=primary]:
        for a in primary_candidates:
            b=stats.get((comp,a["tissue"],a["gene_id"])); cid=digest_bytes("\0".join(["tissue_global_bh_q05_v1",primary,a["tissue"],a["gene_id"]]).encode())
            rows.append({"schema_version":SCHEMA,"analysis_id":run.rid,"source_provenance_sha256":run.source_sha,"comparison_level":"tissue_candidate_projection","analysis_set":"all","feature_id":cid,"primary_dataset":primary,"comparison_dataset":comp,"rank_metric":"candidate_projection","comparison_role":run.datasets[comp]["analysis_role"],"primary_stat":a["zscore"],"comparison_stat":b["zscore"] if b else NA,"primary_direction":z_direction(a["zscore"]),"comparison_direction":z_direction(b["zscore"]) if b else NA,"direction_concordant":NA if not b or 0 in (a["zscore"],b["zscore"]) else a["zscore"]*b["zscore"]>0,"primary_rank_pct":NA,"comparison_rank_pct":NA,"primary_p":a["pvalue"],"comparison_p":b["pvalue"] if b else NA,"primary_q":a["q_gene_tissue"],"comparison_q":b["q_gene_tissue"] if b else NA,"primary_status":"success","comparison_status":b["run_status"] if b else "not_success","n_common":NA,"sample_overlap":"known_or_possible","power_interpretation":"low_power_uncertain" if comp=="GCST90301704" else "descriptive_only"})
    return rows

def stability_summaries(rows):
    return stability_core.summarize(rows)

def main():
    # 2026-09-10 (audit B-1): this integrated `run` entry is RETIRED. It always
    # aborts at M1.6 (execute_mrjti raises unconditionally), so advertising it as
    # the documented entry was misleading. The real chain that produced the frozen
    # results is the per-step entry points listed in README.md ("Execution order"):
    # preflight.py -> prepare_gwas.py -> run_single_gwas_spredixcan.py ->
    # build_all_gene_smultixcan_covariance.py / run_all_gene_smultixcan.py ->
    # validate_existing_all_gene_smultixcan.py -> assemble_m1_downstream.py
    # (+ continue_m1_downstream.py). Fail fast with that pointer instead of
    # dying mid-run.
    raise RuntimeError(
        "execute.py integrated run is retired (audit B-1): it cannot complete "
        "because full MR-JTI (M1.6) is intentionally disabled. Use the real "
        "per-step chain documented in scripts/01_pe_genetic_map/README.md "
        "(execution order); mrjti_minimal.py / mrjti_minimal_runner.R cover the "
        "descriptive MR-JTI layer."
    )
    ap=argparse.ArgumentParser(); ap.add_argument("--analysis-id",default="m1-pe-genetic-map"); ap.add_argument("--smultixcan-root",type=Path,required=True,help="Output root created by the skill genome-wide runner"); ap.add_argument("--smultixcan-validation-root",type=Path,required=True,help="Independent pass validation created by the skill validator"); a=ap.parse_args(); run=Run(a.analysis_id)
    if run.cfg["capability_level"]!="complete": raise RuntimeError("skill-owned M1.5 ingestion requires capability_level=complete")
    full,match=execute_spredixcan(run)
    whole=sorted([r for r in full if r["dataset_id"]=="FIGSHARE_22680904_v2" and r["tissue"]=="Whole_Blood"],key=lambda r:(-abs(r["zscore"]),r["gene_id"]))
    for r in whole: r["direction"]="risk_increasing" if r["zscore"]>0 else "risk_decreasing" if r["zscore"]<0 else NA
    cross=load_validated_cross_tissue(run,full,a.smultixcan_root,a.smultixcan_validation_root); run.status_task("M1.4_whole_blood","core_required","success",{"rows":len(whole)})
    detail,locus=loci(run,full,cross); run.status_task("M1.7_locus_mapping","core_required","success",{"detail":len(detail),"loci":len(locus)})
    stable=stability(run,full,locus); run.status_task("M1.8_stability","core_required","success",{"rows":len(stable)})
    stable_summary=stability_summaries(stable)
    role={"primary":0,"phenotype_sensitivity":1,"broad_sensitivity":2,"low_power_exploratory":3,"exploratory":3}
    full.sort(key=lambda r:(role[r["analysis_role"]],r["tissue"],r["gene_id"])); cross.sort(key=lambda r:(role[r["analysis_role"]],math.inf if r["pvalue"]==NA else r["pvalue"],r["gene_id"])); detail.sort(key=lambda r:(role[r["analysis_role"]],r["statistical_layer"],r["candidate_id"],r["locus_definition"],r["locus_id"])); locus.sort(key=lambda r:(role[r["analysis_role"]],r["statistical_layer"],str(r["chromosome"]),10**20 if r["locus_start"]==NA else int(r["locus_start"]),r["locus_id"])); stable.sort(key=lambda r:(r["comparison_level"],r["analysis_set"],role[r["comparison_role"]],r["feature_id"],r["rank_metric"]))
    outputs={"pe_tissue_gene_full.tsv":full,"pe_whole_blood.tsv":whole,"pe_cross_tissue.tsv":cross,"pe_candidate_locus_map.tsv":detail,"pe_locus_map.tsv":locus,"pe_outcome_stability.tsv":stable}
    for name,rows in outputs.items(): write_tsv(run.stage/name,gov.TABLES[name][0],rows)
    write_tsv(run.stage/"model-match-qc.tsv",["dataset_id","tissue","n_model_snps","n_matched","match_rate"],match)
    mr=execute_mrjti(run,full)
    write_tsv(run.stage/"mrjti.tsv",["analysis_id","task_id","attempt","trigger_dataset","dataset_id","tissue","gene_id","execution_status","scientific_status","status","seed","n_snps","n_folds","n_bootstrap","beta","ci_low","ci_high","pvalue","q_mrjti","mrjti_supported","input_sha256","harmonization_summary","pruning_status","pruning_summary","ld_status","ld_coverage"],mr)
    write_tsv(run.stage/"candidate_block_annotation_bridge.tsv",["analysis_id","candidate_id","annotation_tissue","variant_id","locus_definition","locus_id","annotation_only"],[])
    write_tsv(run.stage/"stability-summary.tsv",["comparison_dataset","comparison_level","analysis_set","rank_metric","n_common","spearman_rho","status"],stable_summary)
    exc=defaultdict(int)
    for r in stable:
        if r["comparison_status"]!="success": exc[(r["comparison_dataset"],r["comparison_level"],r["comparison_status"])]+=1
    exclusions=[{"dataset_id":d,"comparison_level":level,"exclusion_status":status,"n_excluded":n} for (d,level,status),n in sorted(exc.items())]
    write_tsv(run.stage/"stability-exclusions.tsv",["dataset_id","comparison_level","exclusion_status","n_excluded"],exclusions)
    run.status_task("M1.6_mrjti_summary","core_required","success",{"tasks":len(mr),"supported":sum(x["mrjti_supported"] for x in mr)})
    run.status_task("M1.9_artifact_generation","core_required","success",{"tables":list(outputs)})
    (run.stage/"task-ledger.tsv").write_bytes(run.ledger.read_bytes())
    print(run.stage)
if __name__=="__main__": main()
