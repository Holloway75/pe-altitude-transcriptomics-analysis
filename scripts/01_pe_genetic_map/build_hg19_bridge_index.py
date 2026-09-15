#!/usr/bin/env python3
"""Build a reusable canonical-rsID to hg19 1kG candidate map for M1.1.

The map is only a lookup accelerator. ``prepare_gwas.py`` still performs the
scientific checks: unique-candidate closure, FASTA REF validation, indel
normalization, palindrome removal, and final key deduplication.
"""
from __future__ import annotations

import argparse, json, pickle, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from prepare_gwas import bridge_paths_by_chromosome, gwas_reader, norm_chr, op


def canonical_rsid(value):
    value=str(value).strip()
    return value if value.startswith("rs") and value[2:].isdigit() else None


def collect_requested_rsids(gwas, chromosomes):
    requested={chrom:set() for chrom in chromosomes}
    with op(Path(gwas),"r") as handle:
        for row in gwas_reader(handle):
            chrom=norm_chr(str(row.get("chromosome", "")))
            rsid=canonical_rsid(row.get("rsid", ""))
            if chrom in requested and rsid: requested[chrom].add(rsid)
    return requested


def build_index(gwas, bridges, bcftools):
    by_chrom=bridge_paths_by_chromosome(bridges)
    requested=collect_requested_rsids(gwas,by_chrom)
    index=defaultdict(set)
    with tempfile.TemporaryDirectory(prefix="m1-bridge-index-") as td:
        td=Path(td)
        for chrom,path in sorted(by_chrom.items(),key=lambda item:int(item[0])):
            wanted=requested[chrom]
            if not wanted: continue
            ids=td/f"chr{chrom}.rsids"
            ids.write_text("\n".join(sorted(wanted))+"\n")
            command=[bcftools,"query","-i",f"ID=@{ids}","-f","%CHROM\t%POS\t%REF\t%ALT\t%ID\n",str(path)]
            process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
            assert process.stdout is not None
            for line in process.stdout:
                c,pos,ref,alts,ids_out=line.rstrip("\n").split("\t")
                if norm_chr(c)!=chrom: raise RuntimeError("bridge VCF chromosome mismatch")
                for rsid in ids_out.split(";"):
                    if rsid not in wanted: continue
                    for alt in alts.upper().split(","):
                        if alt and alt not in {".","*"} and not alt.startswith("<"):
                            index[int(rsid[2:])].add((int(chrom),int(pos),ref.upper(),alt))
            stderr=process.stderr.read() if process.stderr else ""
            if process.wait()!=0: raise RuntimeError("bridge query failed: "+stderr.strip())
    serial={rsid:sorted(candidates) for rsid,candidates in index.items()}
    metrics={"requested_canonical_rsids":sum(map(len,requested.values())),"mapped_rsids":len(serial),"candidate_records":sum(map(len,serial.values()))}
    return serial,metrics


def flatten_bridge_groups(groups):
    return [path for group in groups for path in group]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--gwas",type=Path,required=True)
    parser.add_argument("--bridge",type=Path,nargs="+",action="append",required=True,
                        help="hg19 VCFs; repeated --bridge groups are flattened")
    parser.add_argument("--bcftools",required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--report",type=Path,required=True)
    args=parser.parse_args()
    bridges=flatten_bridge_groups(args.bridge)
    index,metrics=build_index(args.gwas,bridges,args.bcftools)
    args.out.parent.mkdir(parents=True,exist_ok=True); args.report.parent.mkdir(parents=True,exist_ok=True)
    tmp=args.out.with_suffix(args.out.suffix+".tmp")
    with tmp.open("wb") as handle: pickle.dump(index,handle,protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(args.out)
    payload={"schema_version":"m1-hg19-rsid-bridge-v1","gwas":str(args.gwas.resolve()),"metrics":metrics}
    report_tmp=args.report.with_suffix(args.report.suffix+".tmp")
    report_tmp.write_text(json.dumps(payload,sort_keys=True,indent=2)+"\n")
    report_tmp.replace(args.report)


if __name__=="__main__": main()
