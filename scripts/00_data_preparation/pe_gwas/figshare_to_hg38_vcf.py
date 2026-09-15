#!/usr/bin/env python3
"""Validate Figshare METAL rows against hg38 FASTA and emit biallelic VCF."""
import argparse
import gzip
import json
import math
import os

COMP = str.maketrans("ACGT", "TGCA")
PAL = {frozenset(("A", "T")), frozenset(("C", "G"))}


class IndexedFasta:
    def __init__(self, fasta, fai):
        self.fh = open(fasta, "rb")
        self.index = {}
        with open(fai) as handle:
            for line in handle:
                name, length, offset, bases, width = line.rstrip().split("\t")[:5]
                self.index[name] = tuple(map(int, (length, offset, bases, width)))

    def sequence(self, chrom, pos, size=1):
        length, offset, bases, width = self.index[chrom]
        if pos < 1 or pos + size - 1 > length:
            raise IndexError
        result = bytearray()
        for zero in range(pos - 1, pos - 1 + size):
            self.fh.seek(offset + (zero // bases) * width + zero % bases)
            result.extend(self.fh.read(1))
        return result.decode().upper()


def finite(value):
    try:
        x = float(value)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metal", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--qc-json", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_vcf)), exist_ok=True)
    fa = IndexedFasta(args.fasta, args.fai)
    qc = {k: 0 for k in ("n_input", "invalid_columns", "nonautosomal",
          "invalid_coordinate", "invalid_allele", "invalid_statistic",
          "palindrome_excluded", "reference_mismatch", "snp_same",
          "snp_negative_strand", "indel", "n_output")}
    tmp = args.out_vcf + ".tmp"
    with gzip.open(args.metal, "rt") as src, open(tmp, "w") as out:
        header = src.readline().rstrip("\r\n").split("\t")
        required = ["Chromosome", "Position", "Allele1", "Allele2", "Freq1",
                    "Effect", "StdErr", "P-value"]
        missing = [x for x in required if x not in header]
        if missing:
            raise SystemExit(f"missing METAL columns: {missing}")
        c = {x: header.index(x) for x in required}
        out.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        for chrom in range(1, 23):
            length = fa.index[f"chr{chrom}"][0]
            out.write(f"##contig=<ID=chr{chrom},length={length}>\n")
        for key, typ, desc in (
            ("EA_ALT", "Integer", "1 when METAL Allele1 is ALT; 0 when REF"),
            ("BETA", "Float", "METAL effect for Allele1"),
            ("SE", "Float", "Standard error"), ("EAF", "Float", "Allele1 frequency"),
            ("PVALUE", "Float", "Association P value"),
            ("SOURCE_LINE", "Integer", "One-based METAL data row")):
            out.write(f'##INFO=<ID={key},Number=1,Type={typ},Description="{desc}">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for source_line, line in enumerate(src, 1):
            qc["n_input"] += 1
            f = line.rstrip("\r\n").split("\t")
            if len(f) != len(header):
                qc["invalid_columns"] += 1; continue
            try:
                chrom = int(f[c["Chromosome"]]); pos = int(f[c["Position"]])
            except ValueError:
                qc["invalid_coordinate"] += 1; continue
            if not 1 <= chrom <= 22:
                qc["nonautosomal"] += 1; continue
            a1, a2 = f[c["Allele1"]].upper(), f[c["Allele2"]].upper()
            if not a1 or not a2 or a1 == a2 or set(a1 + a2) - set("ACGT"):
                qc["invalid_allele"] += 1; continue
            beta, se = finite(f[c["Effect"]]), finite(f[c["StdErr"]])
            eaf, pval = finite(f[c["Freq1"]]), finite(f[c["P-value"]])
            if (beta is None or se is None or se <= 0 or eaf is None or
                    not 0 <= eaf <= 1 or pval is None or not 0 <= pval <= 1):
                qc["invalid_statistic"] += 1; continue
            if len(a1) == len(a2) == 1 and frozenset((a1, a2)) in PAL:
                qc["palindrome_excluded"] += 1; continue
            try:
                refbase = fa.sequence(f"chr{chrom}", pos)
            except (KeyError, IndexError):
                qc["invalid_coordinate"] += 1; continue
            status = None
            a1_is_ref = a1 == fa.sequence(f"chr{chrom}", pos, len(a1))
            a2_is_ref = a2 == fa.sequence(f"chr{chrom}", pos, len(a2))
            if a1_is_ref and not a2_is_ref:
                ref, alt, ea_alt = a1, a2, 0
            elif a2_is_ref and not a1_is_ref:
                ref, alt, ea_alt = a2, a1, 1
            elif len(a1) == len(a2) == 1:
                ca1, ca2 = a1.translate(COMP), a2.translate(COMP)
                if ca1 == refbase:
                    ref, alt, ea_alt = ca1, ca2, 0
                elif ca2 == refbase:
                    ref, alt, ea_alt = ca2, ca1, 1
                else:
                    qc["reference_mismatch"] += 1; continue
                status = "snp_negative_strand"
            else:
                qc["reference_mismatch"] += 1; continue
            if status is None:
                status = "indel" if len(a1) != 1 or len(a2) != 1 else "snp_same"
            qc[status] += 1
            info = f"EA_ALT={ea_alt};BETA={beta:.17g};SE={se:.17g};EAF={eaf:.17g};PVALUE={pval:.17g};SOURCE_LINE={source_line}"
            out.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")
            qc["n_output"] += 1
    excluded = sum(qc[k] for k in ("invalid_columns", "nonautosomal",
        "invalid_coordinate", "invalid_allele", "invalid_statistic",
        "palindrome_excluded", "reference_mismatch"))
    qc["accounting_closed"] = qc["n_input"] == qc["n_output"] + excluded
    qc["retained_breakdown_closed"] = qc["n_output"] == qc["snp_same"] + qc["snp_negative_strand"] + qc["indel"]
    if not qc["accounting_closed"] or not qc["retained_breakdown_closed"]:
        raise SystemExit("QC accounting failed")
    os.replace(tmp, args.out_vcf)
    with open(args.qc_json + ".tmp", "w") as out:
        json.dump(qc, out, indent=2, sort_keys=True); out.write("\n")
    os.replace(args.qc_json + ".tmp", args.qc_json)


if __name__ == "__main__":
    main()
