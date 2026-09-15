#!/usr/bin/env python3
"""
Convert one harmonised hg38 GWAS summary file to hg19 by rsID-bridging against
the 1kG hg19 VCF (no liftover). Uses the index built by build_rsid_index.py.

Per-row logic (E=effect_allele, O=other_allele from GWAS; R=ref, A=alt from 1kG;
all upper-case):
  1. rsid must match ^rs\\d+$  -> else drop (cannot bridge; liftover forbidden).
  2. SNP palindrome (len(E)==len(O)==1 and E==comp(O), i.e. {A,T} or {C,G})
     -> drop (strand-ambiguous, standard MR practice).
  3. rsID lookup -> candidates; keep candidates on the SAME chromosome as the GWAS
     row (chr number preserved across hg38<->hg19; kills cross-chr collisions).
  4. Allele match per candidate:
       SNP (all length 1):
           {E,O}=={A,R}           -> same strand, keep E,O
           {comp(E),comp(O)}=={A,R}-> negative strand, emit comp(E),comp(O)
                                       (beta UNCHANGED: physical allele same,
                                        only the letter label moves to + strand)
       indel / mixed (any allele len>1):
           normalize (common prefix+suffix trim) both (R,A) and (O,E); if allele
           sets equal -> match, keep GWAS E,O as-is (indels have no strand).
  5. Exactly one match -> emit; >1 -> ambiguous, drop; 0 -> drop.

Output columns = input columns minus {hm_coordinate_conversion, hm_code}; with
chromosome/base_pair_location replaced by hg19, effect/other_allele possibly
complemented, and variant_id recomputed as chr_pos_other_effect.

Usage:
  python convert_gwas_hg38_to_hg19.py \
      --gwas  /path/GWAS.h.tsv.gz \
      --index data/processed/rsid_hg19_index.pkl \
      --out   /path/out.hg19.tsv.gz \
      [--bgzip /path/bgzip] [--tabix /path/tabix]
"""
import argparse, gzip, json, math, os, pickle, re, shlex, subprocess, sys

_RS = re.compile(r"^rs(\d+)$")
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

def complement(allele):
    """Return complement (uppercase) or None if any base not ACGTN."""
    out = []
    for ch in allele:
        if ch not in _COMP:
            return None
        out.append(_COMP[ch])
    return "".join(out)

def normalize_pair(a, b):
    """Trim common prefix and suffix from allele pair (a,b); return (a',b')."""
    p = 0
    while p < min(len(a), len(b)) and a[p] == b[p]:
        p += 1
    end_a, end_b = len(a), len(b)
    while end_a > p and end_b > p and a[end_a - 1] == b[end_b - 1]:
        end_a -= 1
        end_b -= 1
    return a[p:end_a], b[p:end_b]

def try_match(E, O, R, A):
    """Return (status, Eout, Oout) or None.
    status in {'snp_same','snp_neg','indel'}."""
    g_snp = (len(E) == 1 and len(O) == 1)
    k_snp = (len(R) == 1 and len(A) == 1)
    if g_snp and k_snp:
        if {E, O} == {A, R}:
            return ("snp_same", E, O)
        Ec, Oc = complement(E), complement(O)
        if Ec is not None and Oc is not None and {Ec, Oc} == {A, R}:
            return ("snp_neg", Ec, Oc)
        return None
    # indel / mixed: compare normalized allele sets
    Rt, At = normalize_pair(R, A)
    Ot, Et = normalize_pair(O, E)
    if Et == At and Ot == Rt:
        return ("indel", A, R)
    if Et == Rt and Ot == At:
        return ("indel", R, A)
    return None

def valid_allele(a):
    """Only literal DNA alleles are eligible; symbolic/missing alleles are not."""
    return bool(a) and a != "." and set(a) <= set("ACGT")

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwas", required=True, help="harmonised hg38 .tsv.gz")
    ap.add_argument("--index", required=True, help="rsid_hg19_index.pkl")
    ap.add_argument("--out", required=True, help="output .tsv.gz path")
    ap.add_argument("--qc-json", help="write machine-readable closed QC accounting")
    ap.add_argument("--bgzip", default=os.environ.get("BGZIP", "bgzip"))
    ap.add_argument("--tabix", default=os.environ.get("TABIX", "tabix"))
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    print(f"[load] reading index {args.index} ...", flush=True)
    with open(args.index, "rb") as fh:
        index = pickle.load(fh)
    print(f"[load] index: {len(index)} rsIDs", flush=True)

    # column map
    with gzip.open(args.gwas, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    col = {name: i for i, name in enumerate(header)}
    for need in ("chromosome", "base_pair_location", "effect_allele",
                 "other_allele", "beta", "rsid", "variant_id"):
        if need not in col:
            raise SystemExit(f"[err] input missing required column '{need}'; got {header}")
    DROP = {"hm_coordinate_conversion", "hm_code"}
    out_header = [h for h in header if h not in DROP]

    qc = dict(n_in=0, invalid_record=0, nonrs=0, palindrome=0, nomatch=0, ambiguous=0,
              snp_same=0, snp_neg=0, indel=0, n_out=0)

    # Write to an UNCOMPRESSED temp file with a '#'-prefixed header (tabix skips
    # lines beginning with '#'). The output positions are hg19 but rows are in
    # GWAS-input (hg38) order, so they are NOT sorted by (chr,pos). tabix requires
    # sorted input, so we sort -k1,1 -k2,2n before bgzip below.
    tmp = args.out + ".unsorted.tmp"
    tmp_fh = open(tmp, "w")
    tmp_fh.write("#" + "\t".join(out_header) + "\n")

    def emit(fields_list):
        tmp_fh.write("\t".join(fields_list) + "\n")

    with gzip.open(args.gwas, "rt") as fh:
        fh.readline()  # header already read above via separate open; re-skip
        for line in fh:
            qc["n_in"] += 1
            f = line.rstrip("\n").split("\t")
            if len(f) != len(header):
                qc["invalid_record"] += 1
                continue
            rsid_field = f[col["rsid"]]
            m = _RS.match(rsid_field)
            if not m:
                qc["nonrs"] += 1
                continue
            rid = int(m.group(1))
            E = f[col["effect_allele"]].upper()
            O = f[col["other_allele"]].upper()
            try:
                beta = float(f[col["beta"]])
            except ValueError:
                beta = math.nan
            if not valid_allele(E) or not valid_allele(O) or E == O or not math.isfinite(beta):
                qc["invalid_record"] += 1
                continue
            # 2. palindrome (SNP only)
            if len(E) == 1 and len(O) == 1 and _COMP.get(O) == E:
                qc["palindrome"] += 1
                continue
            cands = index.get(rid)
            if not cands:
                qc["nomatch"] += 1
                continue
            try:
                gchr = int(f[col["chromosome"]])
            except ValueError:
                qc["invalid_record"] += 1
                continue
            if not 1 <= gchr <= 22:
                qc["invalid_record"] += 1
                continue
            matches = []
            for (cchr, cpos, R, A) in cands:
                if cchr != gchr:
                    continue
                r = try_match(E, O, R, A)
                if r is not None:
                    matches.append((cchr, cpos, r))
            if not matches:
                qc["nomatch"] += 1
                continue
            # Duplicate identical records in a reference do not create a real
            # ambiguity; distinct coordinate/allele matches do.
            matches = list(dict.fromkeys(matches))
            if len(matches) > 1:
                qc["ambiguous"] += 1
                continue
            cchr, cpos, (status, Eout, Oout) = matches[0]
            qc[status] += 1
            # build output row preserving order of out_header
            row_out = []
            for name in out_header:
                if name == "chromosome":
                    row_out.append(str(cchr))
                elif name == "base_pair_location":
                    row_out.append(str(cpos))
                elif name == "effect_allele":
                    row_out.append(Eout)
                elif name == "other_allele":
                    row_out.append(Oout)
                elif name == "variant_id":
                    row_out.append(f"{cchr}_{cpos}_{Oout}_{Eout}")
                else:
                    row_out.append(f[col[name]])
            emit(row_out)
            qc["n_out"] += 1

    tmp_fh.close()

    # sort by (chr, pos) then bgzip. Whole-file sort keeps the '#'-prefixed header
    # first ('#' < '1' in C locale). sort -k1,1 -k2,2n groups each chr contiguously
    # and orders positions ascending within each chr -- which is what tabix needs.
    print(f"[sort+bgzip] sorting temp -> {args.out} ...", flush=True)
    sort_pipe = (f"LC_ALL=C sort -k1,1 -k2,2n --parallel=4 -S 2G "
                 f"{shlex.quote(tmp)} | {shlex.quote(args.bgzip)} -c > {shlex.quote(args.out)}")
    subprocess.run(["bash", "-c", sort_pipe], check=True)
    os.unlink(tmp)

    print(f"[index] indexing output with tabix ...", flush=True)
    # Output is a TSV with col1=chromosome, col2=base_pair_location (1-based).
    # Use explicit region columns (not -p vcf, which would mis-parse the TSV).
    subprocess.run([args.tabix, "-s", "1", "-b", "2", "-e", "2", args.out], check=True)

    kept = (qc["n_in"] - qc["invalid_record"] - qc["nonrs"] - qc["palindrome"]
            - qc["nomatch"] - qc["ambiguous"])
    print("[qc] " + " | ".join(f"{k}={v}" for k, v in qc.items()), flush=True)
    print(f"[qc] account check: in({qc['n_in']}) - invalid({qc['invalid_record']}) "
          f"- nonrs({qc['nonrs']}) "
          f"- pal({qc['palindrome']}) - nomatch({qc['nomatch']}) "
          f"- ambig({qc['ambiguous']}) = {kept}; out={qc['n_out']} "
          f"({'OK' if kept == qc['n_out'] else 'MISMATCH'})", flush=True)
    print(f"[qc] matched breakdown: snp_same={qc['snp_same']} "
          f"snp_neg={qc['snp_neg']} indel={qc['indel']}", flush=True)
    if kept != qc["n_out"] or qc["n_out"] != qc["snp_same"] + qc["snp_neg"] + qc["indel"]:
        raise SystemExit("[err] QC accounting did not close")
    if args.qc_json:
        payload = {**qc, "accounting_closed": True,
                   "matched_breakdown_closed": True}
        tmp_json = args.qc_json + ".tmp"
        with open(tmp_json, "w") as fh:
            json.dump(payload, fh, sort_keys=True, indent=2)
            fh.write("\n")
        os.replace(tmp_json, args.qc_json)

if __name__ == "__main__":
    main()
