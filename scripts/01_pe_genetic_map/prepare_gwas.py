#!/usr/bin/env python3
"""M1.1 hg19 normalization and optional hg38-to-hg19 rsID bridging.

Already-hg19 inputs are validated directly against their declared coordinates,
source REF/ALT and the frozen hg19 FASTA.  Only explicitly declared hg38 inputs
use the deterministic dbSNP/1KG bridge.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, math, pickle, random, re, sqlite3, subprocess, tempfile
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

REQUIRED=["chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","p_value","rsid","variant_id"]
OUT=["chromosome","position","effect_allele","non_effect_allele","beta","standard_error","zscore","pvalue","effect_allele_frequency","rsid","pre_normalization_key","variant_id","harmonization_status"]
BASES=set("ACGT"); COMP=str.maketrans("ACGT","TGCA")
def op(p,m): return gzip.open(p,m+"t",newline="") if p.suffix==".gz" else p.open(m,newline="")
def canonical_header(name):
    """Normalize the conventional VCF-style leading ``#`` on column names.

    Some frozen GWAS exporters label the first field ``#chromosome``.  This is
    a header marker, not part of the schema name; retain strict validation for
    every other field after removing it.
    """
    name = (name or "").lstrip("\ufeff")
    return name[1:] if name.startswith("#") else name
def gwas_reader(handle):
    reader=csv.DictReader(handle,delimiter="\t")
    if reader.fieldnames:
        reader.fieldnames=[canonical_header(x) for x in reader.fieldnames]
    missing=set(REQUIRED)-set(reader.fieldnames or [])
    if missing: raise SystemExit("missing columns: "+",".join(sorted(missing)))
    return reader
def norm_chr(x):
    x=x.strip().removeprefix("chr").removeprefix("CHR"); return x
def complement(a): return a.translate(COMP) if len(a)==1 and a in BASES else None
def is_palindrome(a,b): return {a,b} in ({"A","T"},{"C","G"})
def minimal_alleles(pos,ref,alt):
    """VCF minimal representation; repeat left-alignment is importer/FASTA-bound."""
    pos=int(pos); ref,alt=ref.upper(),alt.upper()
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]: ref,alt=ref[:-1],alt[:-1]
    while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]: ref,alt,pos=ref[1:],alt[1:],pos+1
    return pos,ref,alt
def left_align(pos,ref,alt,fetch):
    """Repeat-aware VCF left alignment against the frozen FASTA."""
    pos,ref,alt=minimal_alleles(pos,ref,alt)
    if len(ref)==len(alt): return pos,ref,alt
    while pos>1 and ref[-1]==alt[-1]:
        prev=fetch(pos-1,1).upper()
        if len(prev)!=1: raise ValueError("reference fetch failed")
        ref,alt,pos=prev+ref[:-1],prev+alt[:-1],pos-1
        pos,ref,alt=minimal_alleles(pos,ref,alt)
    return pos,ref,alt

def source_variant_identity(row):
    """Parse a declared source CHR:POS:REF:ALT identity when available."""
    value=str(row.get("variant_id","")).strip()
    for separator in (":","_"):
        parts=value.split(separator)
        if len(parts)==4:
            chrom,pos,ref,alt=parts
            chrom=norm_chr(chrom); ref=ref.upper(); alt=alt.upper()
            if chrom in {str(i) for i in range(1,23)} and pos.isdigit() and ref and alt and not (set(ref+alt)-BASES):
                return chrom,int(pos),ref,alt
    return None

def normalize_hg19(row,reference):
    """Validate an already-hg19 row without requiring 1KG membership."""
    chrom=norm_chr(str(row.get("chromosome",""))); pos=str(row.get("base_pair_location","")).strip()
    if chrom not in {str(i) for i in range(1,23)} or not pos.isdigit(): return None,"invalid_coordinate"
    ea=str(row.get("effect_allele","")).upper().strip(); nea=str(row.get("other_allele","")).upper().strip()
    if not ea or not nea or ea==nea or set(ea+nea)-BASES: return None,"invalid_allele"
    identity=source_variant_identity(row)
    if identity is not None:
        source_chrom,source_pos,_,_=identity
        if source_chrom!=chrom or source_pos!=int(pos): return None,"source_identity_conflict"
    # Registered converters historically encoded variant_id as
    # CHR_POS_other_effect, not CHR_POS_REF_ALT.  Never infer REF/ALT ordering
    # from that field. Resolve it from the row alleles and hg19 FASTA.
    candidates=[]
    for ref,alt in ((ea,nea),(nea,ea)):
        if reference(chrom,int(pos),len(ref)).upper()==ref:
            candidates.append((ref,alt,ea,nea))
    if len(ea)==len(nea)==1:
        cea,cnea=complement(ea),complement(nea)
        for ref,alt in ((cea,cnea),(cnea,cea)):
            if reference(chrom,int(pos),len(ref)).upper()==ref:
                candidates.append((ref,alt,cea,cnea))
    candidates=list(dict.fromkeys(candidates))
    if not candidates: return None,"reference_mismatch"
    identities={(ref,alt) for ref,alt,_,_ in candidates}
    if len(identities)!=1: return None,"source_identity_ambiguous"
    ref,alt=next(iter(identities))
    aligned=[(ae,an) for r,a,ae,an in candidates if (r,a)==(ref,alt)]
    if len(set(aligned))!=1: return None,"source_allele_conflict"
    aligned_ea,aligned_nea=aligned[0]
    original=f"{chrom}:{pos}:{ref}:{alt}"
    norm_pos,norm_ref,norm_alt=left_align(int(pos),ref,alt,lambda p,n:reference(chrom,p,n))
    if reference(chrom,norm_pos,len(norm_ref)).upper()!=norm_ref: return None,"reference_mismatch"
    rec=dict(row,chromosome=chrom,base_pair_location=str(norm_pos),effect_allele=aligned_ea,
             other_allele=aligned_nea,pre_normalization_key=original)
    rec["variant_id"]=f"{chrom}:{norm_pos}:{norm_ref}:{norm_alt}"
    return rec,None

def harmonize(ea,nea,beta,eaf,effect,other,reference_ok=True):
    """Implement the frozen seven-level decision table exactly."""
    vals=[str(x).upper().strip() for x in (ea,nea,effect,other)]; ea,nea,effect,other=vals
    if not reference_ok: return None,"reference_mismatch"
    if any(not x or any(b not in BASES for b in x) for x in vals) or ea==nea or effect==other: return None,"invalid_allele"
    if (len(ea)==len(nea)==1 and is_palindrome(ea,nea)) or (len(effect)==len(other)==1 and is_palindrome(effect,other)): return None,"palindrome_excluded"
    branches=[]
    if (ea,nea)==(effect,other): branches.append(("direct",False))
    if (ea,nea)==(other,effect): branches.append(("swapped",True))
    if len(ea)==len(nea)==len(effect)==len(other)==1 and (complement(ea),complement(nea))==(effect,other): branches.append(("strand_direct",False))
    if len(ea)==len(nea)==len(effect)==len(other)==1 and (complement(ea),complement(nea))==(other,effect): branches.append(("strand_swapped",True))
    if len(branches)!=1: return None,"classification_conflict" if len(branches)>1 else "allele_mismatch"
    status,flip=branches[0]; f=None if eaf in (None,"","NA") else float(eaf)
    return {"beta":-float(beta) if flip else float(beta),"eaf":None if f is None else 1-f if flip else f},status

def normalize_and_bridge(row, dbsnp_matches=None, bridge_matches=None, reference=None):
    """Return one canonical hg19 record or a closed exclusion reason.

    Match providers are callables returning candidate dicts.  Existing rsIDs
    skip dbSNP; missing rsIDs require one exact hg38 dbSNP match.  Every record
    then requires one equivalent hg19 1KG bridge match.
    """
    chrom=norm_chr(row.get("chromosome","")); pos=str(row.get("base_pair_location","")).strip()
    if chrom not in {str(i) for i in range(1,23)} or not pos.isdigit(): return None,"invalid_coordinate"
    ea=str(row.get("effect_allele","")).upper(); nea=str(row.get("other_allele","")).upper()
    if not ea or not nea or any(b not in BASES for b in ea+nea): return None,"invalid_allele"
    rsid=str(row.get("rsid","")).strip()
    if not (rsid.startswith("rs") and rsid[2:].isdigit()):
        hits=list(dbsnp_matches(row) if dbsnp_matches else [])
        if len(hits)!=1: return None,"dbsnp_zero_match" if not hits else "dbsnp_multiple_match"
        rsid=hits[0]["rsid"]
    if bridge_matches is None: raise RuntimeError("missing frozen hg19 1KG bridge provider")
    hits=list(bridge_matches(rsid,row))
    if len(hits)!=1: return None,"bridge_zero_match" if not hits else "bridge_multiple_match"
    h=hits[0]
    if norm_chr(str(h["chromosome"])) not in {str(i) for i in range(1,23)}: return None,"invalid_coordinate"
    ref,alt=str(h["ref"]).upper(),str(h["alt"]).upper()
    if reference is None: raise RuntimeError("missing frozen FASTA provider")
    if reference(norm_chr(str(h["chromosome"])),int(h["position"]),len(ref)).upper()!=ref:return None,"reference_mismatch"
    observed={ea,nea}; equivalent=observed=={ref,alt}; aligned_ea,aligned_nea=ea,nea
    if not equivalent and len(ea)==len(nea)==len(ref)==len(alt)==1:
        comp_ea,comp_nea=complement(ea),complement(nea)
        equivalent={comp_ea,comp_nea}=={ref,alt}
        if equivalent: aligned_ea,aligned_nea=comp_ea,comp_nea
    if not equivalent:return None,"reference_mismatch"
    original=f'{norm_chr(str(h["chromosome"]))}:{h["position"]}:{ref}:{alt}'
    pos,ref,alt=left_align(h["position"],ref,alt,lambda p,n:reference(norm_chr(str(h["chromosome"])),p,n)); rec=dict(row,chromosome=norm_chr(str(h["chromosome"])),base_pair_location=str(pos),effect_allele=aligned_ea,other_allele=aligned_nea,rsid=rsid,pre_normalization_key=original)
    rec["variant_id"]=f'{rec["chromosome"]}:{pos}:{ref}:{alt}'
    return rec,None

class VariantProvider:
    def __init__(self,path,bcftools=None): self.path=Path(path); self.bcftools=bcftools; self.rows=None
    def _tsv(self):
        if self.rows is None:
            self.rows=read_provider(self.path)
        return self.rows
    def dbsnp(self,row):
        chrom=norm_chr(row["chromosome"]); pos=str(row["base_pair_location"]); ea=row["effect_allele"].upper(); nea=row["other_allele"].upper()
        return [r for r in self._query(chrom,pos) if ({r["ref"],r["alt"]}=={ea,nea} or (len(ea)==len(nea)==1 and {r["ref"],r["alt"]}=={complement(ea),complement(nea)}))]
    def bridge(self,rsid,row): return [r for r in self._query(None,None,rsid) if r["rsid"]==rsid]
    def _query(self,chrom,pos,rsid=None):
        if self.path.suffix in {".tsv",".txt"}: return self._tsv()
        if not self.bcftools: raise RuntimeError("VCF provider requires --bcftools")
        region=f"{chrom}:{pos}-{pos}" if chrom else None; cmd=[self.bcftools,"query","-f","%CHROM\t%POS\t%REF\t%ALT\t%ID\n"]
        if region: cmd += ["-r",region]
        if rsid: cmd += ["-i",f'ID="{rsid}"']
        cmd.append(str(self.path)); p=subprocess.run(cmd,text=True,capture_output=True)
        if p.returncode: raise RuntimeError("provider query failed")
        out=[]
        for line in p.stdout.splitlines():
            c,q,ref,alts,ids=line.split("\t");
            for alt in alts.split(","):
                for rid in ids.split(";"): out.append({"chromosome":norm_chr(c),"position":q,"ref":ref.upper(),"alt":alt.upper(),"rsid":rid})
        return out
class Providers:
    def __init__(self,providers): self.providers=providers
    def dbsnp(self,row): return [x for p in self.providers for x in p.dbsnp(row)]
    def bridge(self,rsid,row): return [x for p in self.providers for x in p.bridge(rsid,row)]

class IndexedBridge:
    """Read the Stage 00 rsID bridge index without rescanning 1kG VCFs.

    The index is a candidate map, not a statement that an input variant is
    valid.  ``normalize_and_bridge`` still requires exactly one candidate,
    FASTA REF agreement, allele closure, and indel normalization.
    """
    def __init__(self,path):
        with Path(path).open("rb") as handle: self.index=pickle.load(handle)
    def bridge(self,rsid,row):
        if not (rsid.startswith("rs") and rsid[2:].isdigit()): return []
        chrom=norm_chr(str(row.get("chromosome", "")))
        candidates=self.index.get(int(rsid[2:]),())
        return [{"chromosome":str(c),"position":str(p),"ref":str(ref).upper(),"alt":str(alt).upper(),"rsid":rsid}
                for c,p,ref,alt in candidates if norm_chr(str(c))==chrom]

def bridge_paths_by_chromosome(paths):
    """Bind each frozen 1KG VCF to the chromosome encoded in its filename."""
    out={}
    for path in paths:
        match=re.search(r"(?:^|[._-])chr([0-9]+)(?:[._-]|$)",Path(path).name,flags=re.I)
        if not match: raise RuntimeError(f"cannot determine bridge chromosome: {path}")
        chrom=norm_chr(match.group(1))
        if chrom in out: raise RuntimeError(f"duplicate bridge chromosome: {chrom}")
        out[chrom]=Path(path)
    return out

def build_batch_bridge_cache(gwas, bridge_paths, bcftools, db):
    """Populate an exact rsID bridge cache with one query per chromosome.

    The previous record-at-a-time provider invoked ``bcftools`` once for every
    SNP and every chromosome.  Here we first collect rsIDs by the frozen hg19
    source chromosome, then query each corresponding indexed 1KG VCF once.
    The cache deliberately preserves multiple candidate rows so the existing
    zero/one/multiple-match decision rule remains unchanged.
    """
    if not bcftools: raise RuntimeError("batch bridge cache requires --bcftools")
    by_chrom=bridge_paths_by_chromosome(bridge_paths)
    db.execute("create table bridge_cache(chromosome text, rsid text, position text, ref text, alt text, unique(chromosome,rsid,position,ref,alt))")
    with tempfile.TemporaryDirectory(prefix="m1-bridge-") as td:
        td=Path(td); handles={chrom:(td/f"chr{chrom}.rsids").open("w") for chrom in by_chrom}
        try:
            with op(Path(gwas),"r") as handle:
                for row in gwas_reader(handle):
                    chrom=norm_chr(str(row.get("chromosome", "")))
                    rsid=str(row.get("rsid", "")).strip()
                    if chrom in handles and rsid.startswith("rs") and rsid[2:].isdigit(): handles[chrom].write(rsid+"\n")
        finally:
            for handle in handles.values(): handle.close()
        for chrom,path in sorted(by_chrom.items(),key=lambda x:int(x[0])):
            ids=td/f"chr{chrom}.rsids"
            if not ids.stat().st_size: continue
            cmd=[bcftools,"query","-i",f"ID=@{ids}","-f","%CHROM\t%POS\t%REF\t%ALT\t%ID\n",str(path)]
            proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
            assert proc.stdout is not None
            for line in proc.stdout:
                c,pos,ref,alts,ids_out=line.rstrip("\n").split("\t")
                if norm_chr(c)!=chrom: raise RuntimeError("bridge VCF chromosome mismatch")
                for alt in alts.split(","):
                    for rsid in ids_out.split(";"):
                        if rsid.startswith("rs") and rsid[2:].isdigit(): db.execute("insert or ignore into bridge_cache values(?,?,?,?,?)",(chrom,rsid,pos,ref.upper(),alt.upper()))
            stderr=proc.stderr.read() if proc.stderr else ""
            if proc.wait()!=0: raise RuntimeError("batch bridge query failed: "+stderr.strip())
    db.execute("create index bridge_cache_lookup on bridge_cache(chromosome,rsid)")
    db.commit()

def cached_bridge(db):
    def lookup(rsid,row):
        chrom=norm_chr(str(row.get("chromosome", "")))
        return [{"chromosome":c,"position":p,"ref":ref,"alt":alt,"rsid":rid} for c,rid,p,ref,alt in db.execute("select chromosome,rsid,position,ref,alt from bridge_cache where chromosome=? and rsid=?",(chrom,rsid))]
    return lookup
class FastaReference:
    def __init__(self,path):
        self.path=Path(path); self.index={}; self.handle=None
        fai=Path(str(path)+".fai")
        if fai.is_file():
            for line in fai.read_text().splitlines():
                n,l,o,bw,lw=line.split("\t")[:5]; self.index[n.removeprefix("chr")]=(int(l),int(o),int(bw),int(lw))
        else:
            seq=""; name=None
            for line in self.path.read_text().splitlines():
                if line.startswith(">"):
                    if name:self.index[name]=seq
                    name=line[1:].split()[0].removeprefix("chr"); seq=""
                else: seq+=line.strip()
            if name:self.index[name]=seq
        if self.index and not isinstance(next(iter(self.index.values())),str): self.handle=self.path.open("rb")
    def fetch(self,chrom,pos,n):
        item=self.index[chrom.removeprefix("chr")]
        if isinstance(item,str): return item[pos-1:pos-1+n]
        _,offset,bw,lw=item; start=pos-1; out=b""
        h=self.handle
        if h is None: h=self.path.open("rb")
        try:
            while len(out)<n:
                q=start+len(out); h.seek(offset+(q//bw)*lw+(q%bw)); out+=h.read(min(n-len(out),bw-(q%bw)))
        finally:
            if self.handle is None: h.close()
        return out.decode()
    def close(self):
        if self.handle is not None:
            self.handle.close(); self.handle=None
def read_provider(path):
    with op(path,"r") as h:
        rows=list(csv.DictReader(h,delimiter="\t"))
    need={"chromosome","position","ref","alt","rsid"}
    if not rows or not need<=set(rows[0]): raise RuntimeError("provider schema invalid")
    return rows

def parse(row,dbsnp,bridge,reference,source_build,require_canonical_rsid=False):
    ea,oa=row["effect_allele"].upper().strip(),row["other_allele"].upper().strip()
    if not ea or not oa or any(b not in BASES for b in ea+oa) or ea==oa:return None,"invalid_allele"
    if len(ea)==len(oa)==1 and is_palindrome(ea,oa): return None,"palindrome_excluded"
    try: beta,se,p=float(row["beta"]),float(row["standard_error"]),float(row["p_value"])
    except (ValueError,TypeError): return None,"invalid_statistic"
    if not all(map(math.isfinite,(beta,se,p))) or se<=0 or not 0<=p<=1:return None,"invalid_statistic"
    if source_build not in {"hg19","hg38"}: return None,"invalid_source_build"
    if require_canonical_rsid:
        rsid=str(row.get("rsid","")).strip()
        if not (rsid.startswith("rs") and rsid[2:].isdigit()): return None,"missing_rsid"
    if source_build=="hg19":
        bridged,reason=normalize_hg19(row,reference.fetch)
    else:
        if dbsnp is None or bridge is None: return None,"missing_bridge_provider"
        bridged,reason=normalize_and_bridge(row,dbsnp.dbsnp,bridge.bridge,reference.fetch)
    if reason:return None,reason
    eaf=row.get("effect_allele_frequency","NA").strip() or "NA"
    rec=dict(chromosome=bridged["chromosome"],position=bridged["base_pair_location"],effect_allele=bridged["effect_allele"],non_effect_allele=bridged["other_allele"],beta=f"{beta:.17g}",standard_error=f"{se:.17g}",zscore=f"{beta/se:.17g}",pvalue=f"{p:.17g}",effect_allele_frequency=eaf,rsid=bridged["rsid"],pre_normalization_key=bridged["pre_normalization_key"],variant_id=bridged["variant_id"],harmonization_status="unmatched_model")
    return rec,None

def direction_audit(rows, model_alleles, seed, per_status=100, tolerance=1e-10):
    """Fixed-seed independent direction audit; raises on the first mismatch."""
    groups=defaultdict(list)
    for r in rows: groups[r["harmonization_status"]].append(r)
    rng=random.Random(seed)
    for status in sorted(groups):
        sample=groups[status] if len(groups[status])<=per_status else rng.sample(groups[status],per_status)
        for r in sample:
            m=model_alleles[r["rsid"]]
            got,s=harmonize(r["source_ea"],r["source_nea"],r["source_beta"],r.get("source_eaf"),m[0],m[1])
            if s!=status or got is None or abs(float(r["beta"])-got["beta"])>tolerance or abs(float(r["zscore"])-float(r["beta"])/float(r["standard_error"]))>tolerance: raise ValueError("direction_audit_failed")
    return True

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--gwas",type=Path,required=True); ap.add_argument("--out",type=Path,required=True); ap.add_argument("--qc",type=Path,required=True); ap.add_argument("--dbsnp",type=Path); ap.add_argument("--bridge",type=Path,nargs="*"); ap.add_argument("--bridge-index",type=Path); ap.add_argument("--fasta",type=Path,required=True); ap.add_argument("--source-build",required=True,choices=("hg19","hg38")); ap.add_argument("--bcftools"); a=ap.parse_args(); a.out.parent.mkdir(parents=True,exist_ok=True); a.qc.parent.mkdir(parents=True,exist_ok=True)
    if a.source_build=="hg19" and (a.bridge_index or a.bridge or a.dbsnp): raise SystemExit("hg19 input must not use dbSNP/1KG bridge options")
    if a.source_build=="hg38" and (not a.dbsnp or bool(a.bridge_index)==bool(a.bridge)): raise SystemExit("hg38 input requires --dbsnp and exactly one of --bridge-index or --bridge")
    if not a.fasta.is_file(): raise SystemExit("missing frozen FASTA provider")
    dbsnp=VariantProvider(a.dbsnp,a.bcftools) if a.dbsnp else None; reference=FastaReference(a.fasta)
    bridge=IndexedBridge(a.bridge_index) if a.bridge_index else Providers([VariantProvider(x,a.bcftools) for x in a.bridge]) if a.bridge else None
    counts=Counter(n_input=0,n_retained=0)
    with tempfile.NamedTemporaryFile(suffix=".sqlite") as tmp, sqlite3.connect(tmp.name) as db:
        db.execute("create table records(vkey text,sitekey text,payload text,source_line integer)")
        with op(a.gwas,"r") as h:
            rd=gwas_reader(h)
            for line,row in enumerate(rd,2):
                counts["n_input"]+=1; rec,reason=parse(row,dbsnp,bridge,reference,a.source_build,require_canonical_rsid=True)
                if reason: counts[reason]+=1; continue
                if rec["pre_normalization_key"]!=rec["variant_id"]: counts["n_key_changed_by_normalization"]+=1
                key=rec["variant_id"]; parts=key.split(":"); sitekey=":".join(parts[:2]); db.execute("insert into records values(?,?,?,?)",(key,sitekey,"\t".join(rec[x] for x in OUT),line))
        db.commit(); n_valid=db.execute("select count(*) from records").fetchone()[0]
        multiallelic={r[0] for r in db.execute("select sitekey from records group by sitekey having count(distinct vkey)>1")}
        counts["multiallelic_removed"]=sum(r[0] in multiallelic for r in db.execute("select sitekey from records"))
        conflict={r[0] for r in db.execute("select vkey from records where sitekey not in (select sitekey from records group by sitekey having count(distinct vkey)>1) group by vkey having count(distinct payload)>1")}; counts["duplicate_conflict"]=sum(r[0] in conflict and r[1] not in multiallelic for r in db.execute("select vkey,sitekey from records"))
        with op(a.out,"w") as h:
            wr=csv.writer(h,delimiter="\t",lineterminator="\n"); wr.writerow(OUT)
            for key,sitekey,payload,_ in db.execute("select vkey,sitekey,payload,min(source_line) from records group by vkey,sitekey,payload order by vkey"):
                if sitekey in multiallelic or key in conflict: continue
                wr.writerow(payload.split("\t")); counts["n_retained"]+=1
    reference.close()
    counts["duplicate_identical_collapsed"]=n_valid-counts["multiallelic_removed"]-counts["duplicate_conflict"]-counts["n_retained"]
    counts["closure_difference"]=counts["n_input"]-counts["n_retained"]-sum(v for k,v in counts.items() if k not in {"n_input","n_retained","closure_difference","n_key_changed_by_normalization"})
    with a.qc.open("w",newline="") as h: w=csv.writer(h,delimiter="\t",lineterminator="\n"); w.writerow(["metric","value"]); w.writerows(sorted(counts.items()))
    if counts["closure_difference"]: raise SystemExit("QC closure failed")
if __name__=="__main__": main()
