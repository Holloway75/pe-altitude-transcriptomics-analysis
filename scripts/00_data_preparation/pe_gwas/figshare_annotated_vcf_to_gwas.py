#!/usr/bin/env python3
"""Convert a dbSNP-annotated Figshare VCF to the common harmonized GWAS TSV."""
import argparse, gzip, json, os, subprocess

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--qc-json")
    args = ap.parse_args()
    fmt = r"%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/EA_ALT\t%INFO/BETA\t%INFO/SE\t%INFO/EAF\t%INFO/PVALUE\n"
    proc = subprocess.Popen([args.bcftools, "query", "-f", fmt, args.vcf], stdout=subprocess.PIPE, text=True)
    qc = {"n_input_hg38_normalized": 0, "dbsnp_matched": 0,
          "dbsnp_unmatched_or_nonunique": 0}
    with gzip.open(args.out + ".tmp", "wt") as out:
        out.write("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trsid\tvariant_id\n")
        for line in proc.stdout:
            qc["n_input_hg38_normalized"] += 1
            chrom, pos, ref, alt, rid, ea_alt, beta, se, eaf, pval = line.rstrip().split("\t")
            if rid == "." or not rid.startswith("rs") or ";" in rid:
                qc["dbsnp_unmatched_or_nonunique"] += 1
                continue
            qc["dbsnp_matched"] += 1
            chrom = chrom.removeprefix("chr")
            effect, other = (alt, ref) if ea_alt == "1" else (ref, alt)
            out.write("\t".join((chrom, pos, effect, other, beta, se, eaf, pval,
                                  rid, f"{chrom}_{pos}_{other}_{effect}")) + "\n")
    if proc.wait() != 0:
        raise SystemExit("bcftools query failed")
    os.replace(args.out + ".tmp", args.out)
    qc["accounting_closed"] = (qc["n_input_hg38_normalized"] ==
        qc["dbsnp_matched"] + qc["dbsnp_unmatched_or_nonunique"])
    if not qc["accounting_closed"]:
        raise SystemExit("dbSNP bridge accounting failed")
    if args.qc_json:
        with open(args.qc_json + ".tmp", "w") as out:
            json.dump(qc, out, indent=2, sort_keys=True); out.write("\n")
        os.replace(args.qc_json + ".tmp", args.qc_json)

if __name__ == "__main__": main()
