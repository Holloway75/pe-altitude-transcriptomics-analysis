#!/usr/bin/env python3
"""Prepare the three small, versioned mapping resources required by M1.3."""
from __future__ import annotations
import argparse, gzip, hashlib, json, os, re, urllib.request
from pathlib import Path

GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
LDBLOCK_URL = "https://bitbucket.org/nygcresearch/ldetect-data/raw/master/EUR/fourier_ls-all.bed"
HIGH_LD_SOURCE = "https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)"
HIGH_LD = """1 48000000 52000000
2 86000000 100500000
2 134500000 138000000
2 183000000 190000000
3 47500000 50000000
3 83500000 87000000
3 89000000 97500000
5 44500000 50500000
5 98000000 100500000
5 129000000 132000000
5 135500000 138500000
6 25000000 35000000
6 57000000 64000000
6 140000000 142500000
7 55000000 66000000
8 7000000 13000000
8 43000000 50000000
8 112000000 115000000
10 37000000 43000000
11 46000000 57000000
11 87500000 90500000
12 33000000 40000000
12 109500000 112000000
20 32000000 34500000
"""

def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20), b""): h.update(b)
    return h.hexdigest()

def download(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists(): return
    tmp=path.with_suffix(path.suffix+".tmp")
    with urllib.request.urlopen(url) as source, tmp.open("wb") as out:
        while b:=source.read(1<<20): out.write(b)
    tmp.replace(path)

def attrs(text):
    return dict(re.findall(r'(\w+) "([^"]+)"', text))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument(
        "--external-root",
        type=Path,
        default=Path(os.environ["EXTERNAL_DATA_ROOT"]) if os.environ.get("EXTERNAL_DATA_ROOT") else None,
        help="external data root (default: EXTERNAL_DATA_ROOT)",
    )
    a=ap.parse_args()
    if a.external_root is None:
        ap.error("--external-root is required when EXTERNAL_DATA_ROOT is unset")
    raw=a.external_root/"annotations/source/gencode.v19.annotation.gtf.gz"
    gene=a.external_root/"annotations/gencode_hg19_gene_annotation.tsv.gz"
    blocks=a.external_root/"1kGenome/hg19_eur_ld_blocks/blocks.tsv"
    complex_bed=a.external_root/"annotations/hg19_complex_ld_regions.bed"
    download(GENCODE_URL, raw)
    source_blocks=blocks.with_name("fourier_ls-all.bed.source")
    download(LDBLOCK_URL, source_blocks)
    gene.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(raw,"rt") as src, gzip.open(gene.with_suffix(".tmp"),"wt") as out:
        out.write("gene_id\tgene_name\tchromosome\tstart\tend\tstrand\tgene_type\ttss\n")
        for line in src:
            if line.startswith("#"): continue
            f=line.rstrip().split("\t")
            chrom=f[0].removeprefix("chr")
            if len(f)!=9 or f[2]!="gene" or chrom not in {str(i) for i in range(1,23)}: continue
            x=attrs(f[8]); gid=x.get("gene_id","").split(".")[0]
            if not gid: continue
            start,end=int(f[3]),int(f[4]); strand=f[6]; tss=start if strand=="+" else end
            out.write(f"{gid}\t{x.get('gene_name','NA')}\t{chrom}\t{start}\t{end}\t{strand}\t{x.get('gene_type',x.get('gene_biotype','NA'))}\t{tss}\n")
    gene.with_suffix(".tmp").replace(gene)
    blocks.parent.mkdir(parents=True, exist_ok=True)
    with source_blocks.open() as src, blocks.open("w") as out:
        out.write("chromosome\tstart\tend\tld_version\tlocus_id\n")
        for line in src:
            f=line.split()
            if len(f)<3: continue
            chrom=f[0].removeprefix("chr")
            if chrom not in {str(i) for i in range(1,23)}: continue
            start,end=int(float(f[1])),int(float(f[2])); lid=f"BP2016_EUR:{chrom}:{start}-{end}"
            out.write(f"{chrom}\t{start}\t{end}\tBerisaPickrell2016_EUR_hg19\t{lid}\n")
    complex_bed.parent.mkdir(parents=True, exist_ok=True)
    with complex_bed.open("w") as out:
        out.write("chromosome\tstart\tend\tregion_id\tsource\n")
        for i,line in enumerate(HIGH_LD.splitlines(),1):
            c,s,e=line.split(); out.write(f"{c}\t{s}\t{e}\thigh_ld_{i:02d}\t{HIGH_LD_SOURCE}\n")
    meta={"gencode":{"url":GENCODE_URL,"source_sha256":digest(raw),"artifact_sha256":digest(gene)},"ld_blocks":{"url":LDBLOCK_URL,"source_sha256":digest(source_blocks),"artifact_sha256":digest(blocks)},"complex_regions":{"url":HIGH_LD_SOURCE,"artifact_sha256":digest(complex_bed)}}
    (a.external_root/"annotations/module01-resources.json").write_text(json.dumps(meta,sort_keys=True,indent=2)+"\n")
    print(json.dumps(meta,sort_keys=True))
if __name__=="__main__": main()
