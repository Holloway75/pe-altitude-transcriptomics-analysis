#!/usr/bin/env bash
# verify_conversion.sh OUT SRC VCF_PATTERN BCFTOOLS
# End-to-end checks for one converted file:
#   1. row-count accounting (recomputed from files; full drop accounting is in
#      the converter's .qc.json)
#   2. palindrome-zero (SNPs with {A,T} or {C,G} effect/other must be 0)
#   3. coordinate bounds (base_pair_location <= b37 chr length, from VCF ##contig)
#   4. spot-check 10 random rows via tabix range query on the 1kG VCF
#   5. beta conservation (sample 5000 rows: output beta == source beta, joined by rsid;
#      rows where effect/other were complemented for -strand are skipped)
#
# 2026-09-10 (audit B-6): checks 2-5 are now asserted and the script exits
# non-zero when any of them fails (previously it always exited 0, so a FAIL
# line in the log had no gate effect). Sampling is seeded deterministically
# ("PE20260910") so repeated verifications of the same files are reproducible.
set -uo pipefail

OUT="$1"; SRC="$2"; VCF_PATTERN="$3"; BCFTOOLS="$4"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/verify_pe_gwas.XXXXXX")
trap 'rm -rf "${WORK}"' EXIT

FAILURES=0
fail_check() { FAILURES=$((FAILURES + 1)); }

echo "########## verify: $(basename "$OUT") ##########"

# 1. row counts (informational; assertion is only that output has data rows
#    and did not grow relative to source)
N_IN=$(zcat "$SRC" 2>/dev/null | wc -l)            # incl. header
N_OUT=$(zcat "$OUT" 2>/dev/null | wc -l)
echo "[1] rows: src(incl header)=$N_IN  out(incl header)=$N_OUT"
if [[ "$N_OUT" -lt 2 || "$N_OUT" -gt "$N_IN" ]]; then
  echo "  [FAIL] implausible row counts (expect 2 <= out <= src)"
  fail_check
fi

# 2. palindrome zero
PAL=$(zcat "$OUT" 2>/dev/null | awk -F'\t' 'NR>1{
  e=$3;o=$4; if(length(e)==1&&length(o)==1&&
  ((e=="A"&&o=="T")||(e=="T"&&o=="A")||(e=="C"&&o=="G")||(e=="G"&&o=="C"))) c++
} END{print c+0}')
echo "[2] palindromic SNP rows in output: $PAL  (expect 0)"
if [[ "$PAL" -ne 0 ]]; then
  echo "  [FAIL] palindromic rows present"
  fail_check
fi

# 3. coordinate bounds (b37 chr lengths, standard GRCh37)
# NOTE: use n=split(...) not NF (NF is 0 in BEGIN before any record is read).
OVER=$(zcat "$OUT" 2>/dev/null | awk -F'\t' '
  BEGIN{
    n=split("1 249250621 2 243199373 3 198022430 4 191154276 5 180915260 6 171115067 7 159138663 8 146364022 9 141213431 10 135534747 11 135006516 12 133851895 13 115169878 14 107349540 15 102531392 16 90354753 17 81195210 18 78077248 19 59128983 20 63025520 21 48129895 22 51304566", a, " ");
    for(i=1;i<=n;i+=2) L[a[i]]=a[i+1];
  } NR>1 { if(($2+0) > (L[$1]+0)) c++ } END{print c+0}')
echo "[3] rows with pos > b37 chr length: $OVER  (expect 0)"
if [[ "$OVER" -ne 0 ]]; then
  echo "  [FAIL] coordinate-bound violations present"
  fail_check
fi

# 4. spot-check 10 random rows via tabix range query on 1kG VCF
echo "[4] spot-check 10 rows (tabix range query on 1kG, fixed seed):"
# pick 10 random data rows
SAMPLE=$(zcat "$OUT" 2>/dev/null | awk -F'\t' 'NR>1{print}' \
  | shuf -n 10 --random-source=<(yes PE20260910))
N_SAMPLED=$(printf '%s\n' "$SAMPLE" | grep -c . || true)
FAIL=0
if [[ "$N_SAMPLED" -ne 10 ]]; then
  echo "  [FAIL] sampled $N_SAMPLED rows, need 10 for the spot-check"
  FAIL=1
fi
while IFS=$'\t' read -r chr pos effa otha beta se eaf p rsid rest; do
  [[ -z "$chr" ]] && continue   # inert when nothing was sampled
  VCF=$(echo "$VCF_PATTERN" | sed "s/{chr}/$chr/")
  # range query can return multiple overlapping records (SNP + SV at same locus);
  # pick the one whose ID (col3) matches our rsid, skipping SV symbolic alleles.
  rec=$("$BCFTOOLS" view -H "$VCF" -r "${chr}:${pos}-${pos}" 2>/dev/null \
        | awk -F'\t' -v rid="$rsid" '$3==rid && $5 !~ /^</ {print; exit}')
  if [[ -z "$rec" ]]; then
    echo "  [FAIL] rsid=$rsid chr=$chr pos=$pos : no 1kG record with this ID at this position"; FAIL=1; continue
  fi
  rid=$(echo "$rec" | cut -f3)
  ref=$(echo "$rec" | cut -f4 | tr 'a-z' 'A-Z')
  alt=$(echo "$rec" | cut -f5 | tr 'a-z' 'A-Z')
  ok_id=0; ok_alleles=0
  [[ "$rid" == "$rsid" ]] && ok_id=1
  # alleles: {eff,oth} must equal {ref,alt} (allowing orientation)
  if [[ ( "$effa" == "$alt" && "$otha" == "$ref" ) || ( "$effa" == "$ref" && "$otha" == "$alt" ) ]]; then
    ok_alleles=1
  fi
  printf "  rsid=%s chr=%s pos=%s ref=%s alt=%s  id_ok=%s alleles_ok=%s\n" \
    "$rsid" "$chr" "$pos" "$ref" "$alt" "$ok_id" "$ok_alleles"
  [[ $ok_id -eq 1 && $ok_alleles -eq 1 ]] || FAIL=1
done <<< "$SAMPLE"
echo "  spot-check overall: $([ $FAIL -eq 0 ] && echo PASS || echo FAIL)"
[[ $FAIL -eq 0 ]] || fail_check

# 5. beta conservation (sample 5000 output rows, join source by rsid)
echo "[5] beta conservation (5000-sample, joined by rsid, fixed seed):"
zcat "$OUT" 2>/dev/null | awk -F'\t' 'NR>1{print $9"\t"$3"\t"$4"\t"$5}' \
  | shuf -n 5000 --random-source=<(yes PE20260910) > "${WORK}/beta_sample.tsv"
# build source rsid->(eff,oth,beta) hash for sampled rsids
zcat "$SRC" 2>/dev/null | awk -F'\t' '
  FNR==NR { want[$1]=1; next }
  FNR==1  { for(i=1;i<=NF;i++) if($i~/^rsid$/) c=i; else if($i~/^effect_allele$/) e=i; else if($i~/^other_allele$/) o=i; else if($i~/^beta$/) b=i; next }
  { if(($c) in want) print $c"\t"$e"\t"$o"\t"$b }
' "${WORK}/beta_sample.tsv" - > "${WORK}/src_beta.tsv" 2>/dev/null
# compare: for sampled output rows where alleles unchanged from source, beta must match
BETA_STATS=$(awk -F'\t' '
  FNR==NR { src_beta[$1"_"$2"_"$3]=$4; next }
  { key=$1"_"$2"_"$3;
    if(key in src_beta){
      checked++;
      if($4 != src_beta[key]) mismatch++;
    }
  }
  END{ print checked+0, mismatch+0 }
' "${WORK}/src_beta.tsv" "${WORK}/beta_sample.tsv")
read -r N_CHECKED N_MISMATCH <<< "$BETA_STATS"
echo "  checked(allele-unchanged rows)=$N_CHECKED beta_mismatch=$N_MISMATCH"
if [[ "$N_MISMATCH" -ne 0 || "$N_CHECKED" -eq 0 ]]; then
  echo "  [FAIL] beta conservation violated (or no comparable rows)"
  fail_check
fi

echo "########## end verify: $(basename "$OUT") ##########"
if [[ "$FAILURES" -ne 0 ]]; then
  echo "[verify] FAIL: ${FAILURES} check(s) failed for $(basename "$OUT")"
  exit 1
fi
echo "[verify] PASS: all checks passed for $(basename "$OUT")"
exit 0
