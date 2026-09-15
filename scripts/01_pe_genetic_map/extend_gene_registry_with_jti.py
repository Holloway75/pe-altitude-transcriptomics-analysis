#!/usr/bin/env python3
"""Extend the M1 gene registry with coordinate-missing genes from frozen JTI DBs."""
import argparse,csv,gzip,hashlib,json,sqlite3
from pathlib import Path
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--annotation",type=Path,required=True); ap.add_argument("--models",type=Path,required=True); ap.add_argument("--report",type=Path,required=True); a=ap.parse_args()
    before=sha(a.annotation)
    with gzip.open(a.annotation,"rt",newline="") as h: rows=list(csv.DictReader(h,delimiter="\t"))
    by={x["gene_id"].split('.')[0]:x for x in rows}; model={}
    for db in sorted(a.models.glob("JTI_*.db")):
        with sqlite3.connect(db) as c:
            for gene,name,kind in c.execute("select gene,genename,geneannot from extra"):
                gid=gene.split('.')[0]; prior=model.get(gid)
                current=(name or "NA",kind or "NA")
                if prior and prior!=current and prior[0]!=current[0]: raise RuntimeError(f"JTI gene metadata conflict: {gid}")
                model[gid]=current
    added=[]
    for gid in sorted(set(model)-set(by)):
        name,kind=model[gid]; added.append({"gene_id":gid,"gene_name":name,"chromosome":"NA","start":"NA","end":"NA","strand":"NA","gene_type":"JTI_model_no_frozen_coordinate:"+kind,"tss":"NA"})
    merged=rows+added
    tmp=a.annotation.with_suffix(a.annotation.suffix+".tmp")
    with gzip.open(tmp,"wt",newline="") as h:
        w=csv.DictWriter(h,fieldnames=["gene_id","gene_name","chromosome","start","end","strand","gene_type","tss"],delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(merged)
    tmp.replace(a.annotation)
    report={"status":"pass","annotation":str(a.annotation.resolve()),"before_sha256":before,"after_sha256":sha(a.annotation),"coordinate_genes":len(rows),"jti_genes":len(model),"added_coordinate_missing_jti_genes":len(added),"total_registry_genes":len(merged)}
    a.report.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n"); print(json.dumps(report,sort_keys=True))
if __name__=="__main__": main()
