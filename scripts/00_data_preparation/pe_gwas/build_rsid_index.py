#!/usr/bin/env python3
"""
Build a 1kG hg19 rsID -> [(chr, pos, ref, alt), ...] lookup, filtered to the
rsIDs that actually appear in the given GWAS summary file(s).

Why filter: the full 1kG panel has ~84M records; only ~15-20M share an rsID
with these GWAS. Filtering up front (via awk membership test, C-level fast)
keeps the in-memory dict and pickle small and the downstream join cheap.

Output: a pickled dict {rsid_int: [(chr:int, pos:int, ref:str, alt:str), ...]}.
Multi-candidate rsIDs (same rsID at >1 1kG site) are preserved verbatim so the
converter can later disambiguate by ref/alt + chromosome.

Usage:
  python build_rsid_index.py \
      --gwas  /path/A.h.tsv.gz /path/B.h.tsv.gz \
      --vcf-dir  /path/to/1kGenome/hg19_vcf \
      --vcf-pattern 'ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz' \
      --out  data/processed/rsid_hg19_index.pkl \
      [--bcftools /usr/bin/bcftools] [--threads 4]
"""
import argparse, gzip, os, pickle, subprocess, sys, tempfile

RS_RE = None  # compiled lazily
MISSING_CHROMOSOMES = []  # populated by stream_filtered_1kg; checked in main

def find_rsid_col(header_fields):
    """Return 0-based index of the rsid column in a GWAS header line."""
    for i, f in enumerate(header_fields):
        if f.strip().lower() == "rsid":
            return i
    raise SystemExit(f"no 'rsid' column in header: {header_fields}")

def collect_gwas_rsids(gwas_paths):
    """Stream each GWAS, collect the set of rsID integer parts (rs### -> int)."""
    import re
    pat = re.compile(r"^rs(\d+)$")
    s = set()
    for p in gwas_paths:
        n = 0
        with gzip.open(p, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = find_rsid_col(header)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= col:
                    continue
                m = pat.match(f[col])
                if m:
                    s.add(int(m.group(1)))
                    n += 1
        print(f"[gwas] {os.path.basename(p)}: {n} rs-form rows; "
              f"running unique set size = {len(s)}", flush=True)
    return s

def stream_filtered_1kg(vcf_dir, pattern, set_file, bcftools, threads):
    """
    For each chr 1..22, run:
      bcftools query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n' VCF
        | awk 'NR==FNR{seen[$1]=1;next} ($1 in seen)' SET_FILE -
    Yields lines: rsid chr pos ref alt  (only rsIDs present in the GWAS set).
    The awk reads SET_FILE first (NR==FNR true), building the membership hash,
    then filters the piped bcftools stream. Re-reading the ~20M-line set file
    per chromosome is cheap (C fread, ~1s/chr).
    """
    # set_file holds bare rsID integers (rs### -> int, "rs" stripped). The bcftools
    # stream emits IDs WITH the "rs" prefix, so sub() it off before membership test.
    awk = 'NR==FNR{seen[$1]=1;next} {id=$1; sub(/^rs/,"",id); if(id in seen) print}'
    for c in range(1, 23):
        vcf = os.path.join(vcf_dir, pattern.format(chr=c))
        if not os.path.exists(vcf):
            # A missing chromosome would silently truncate the index (all its
            # rsIDs dropped from the downstream conversion). Record and let the
            # caller fail loudly at the end of streaming.
            print(f"[1kg] chr{c}: MISSING {vcf}", file=sys.stderr)
            MISSING_CHROMOSOMES.append(c)
            continue
        qcmd = [bcftools, "query", "-f", "%ID\t%CHROM\t%POS\t%REF\t%ALT\n", "-r", str(c), vcf]
        # bcftools query -> awk filter; both C-level, no Python in the hot path
        awkcmd = ["awk", awk, set_file, "-"]
        with subprocess.Popen(qcmd, stdout=subprocess.PIPE) as bcft, \
             subprocess.Popen(awkcmd, stdin=bcft.stdout, stdout=subprocess.PIPE,
                              text=True, bufsize=1 << 20) as aw:
            bcft.stdout.close()  # awk now owns the pipe; parent must not hold it
            n = 0
            for line in aw.stdout:
                yield line
                n += 1
            rc_awk = aw.wait()
            rc_bcft = bcft.wait()
        # 2026-09-10 (audit B-8): check BOTH exit codes. A non-zero bcftools
        # (e.g. truncated VCF) with a zero awk status previously passed
        # unnoticed, silently shrinking the index.
        if rc_awk != 0:
            raise SystemExit(f"[1kg] chr{c}: awk exited {rc_awk}")
        if rc_bcft != 0:
            raise SystemExit(f"[1kg] chr{c}: bcftools exited {rc_bcft}; "
                             f"refusing to write a silently truncated index")
        print(f"[1kg] chr{c}: kept {n} rs-matched records", flush=True)

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwas", nargs="+", required=True, help="harmonised hg38 GWAS .tsv.gz file(s)")
    ap.add_argument("--rsid-set", help="persistent bare-integer rsID set; create if absent, reuse if present")
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--vcf-pattern", required=True,
                    help="python format string with {chr} placeholder")
    ap.add_argument("--out", required=True, help="output .pkl path")
    ap.add_argument("--bcftools", default="/usr/bin/bcftools")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    set_file = args.rsid_set
    remove_set_file = False
    if set_file and os.path.isfile(set_file) and os.path.getsize(set_file) > 0:
        with open(set_file) as fh:
            n_ids = sum(1 for _ in fh)
        print(f"[1/3] reusing persistent rsID set: {set_file} ({n_ids} IDs)", flush=True)
    else:
        print("[1/3] collecting rsID set from GWAS ...", flush=True)
        s = collect_gwas_rsids(args.gwas)
        print(f"[1/3] unique rsIDs in GWAS union: {len(s)}", flush=True)
        if set_file:
            os.makedirs(os.path.dirname(os.path.abspath(set_file)), exist_ok=True)
            tmp_set = set_file + ".tmp"
            with open(tmp_set, "w") as tf:
                for i in sorted(s):
                    tf.write(f"{i}\n")
            os.replace(tmp_set, set_file)
        else:
            with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, dir="/tmp") as tf:
                for i in s:
                    tf.write(f"{i}\n")
                set_file = tf.name
            remove_set_file = True
    print(f"[2/3] rsID set written to {set_file} for awk membership filter", flush=True)

    index = {}
    n_recs = 0
    for line in stream_filtered_1kg(args.vcf_dir, args.vcf_pattern, set_file,
                                    args.bcftools, args.threads):
        f = line.split("\t")
        # f = [ID, CHROM, POS, REF, ALT]
        try:
            rid = int(f[0][2:])
        except (ValueError, IndexError):
            continue
        chrom = int(f[1])
        pos = int(f[2])
        ref = f[3].upper()
        # A multiallelic VCF row represents one candidate per ALT.  Keeping the
        # comma-joined ALT string would make every biallelic GWAS comparison
        # fail, including otherwise valid indels.
        for alt in f[4].rstrip("\n").upper().split(","):
            if alt in {"", ".", "*"} or alt.startswith("<"):
                continue
            index.setdefault(rid, []).append((chrom, pos, ref, alt))
            n_recs += 1

    multi = sum(1 for v in index.values() if len(v) > 1)
    print(f"[3/3] loaded {n_recs} 1kG records -> {len(index)} unique rsIDs "
          f"({multi} have >1 candidate, kept as lists for ref/alt disambiguation)",
          flush=True)

    if MISSING_CHROMOSOMES:
        raise SystemExit(
            f"[1kg] chromosome VCFs missing: {MISSING_CHROMOSOMES}; "
            f"refusing to write a silently truncated index")

    tmp_out = args.out + ".tmp"
    with open(tmp_out, "wb") as fh:
        pickle.dump(index, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_out, args.out)
    sz = os.path.getsize(args.out) / (1 << 20)
    print(f"[done] wrote {args.out} ({sz:.1f} MiB)", flush=True)

    if remove_set_file:
        os.unlink(set_file)

if __name__ == "__main__":
    main()
