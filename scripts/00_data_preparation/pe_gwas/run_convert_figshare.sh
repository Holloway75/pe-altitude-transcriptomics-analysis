#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../lib/load_config.sh"

DATASET=FIGSHARE_22680904_v2
SRC="${PE_GWAS_ROOT}/${DATASET}/raw/metal_preec_European_allBiobanks_omitNone_1.txt.gz"
OUT="${PE_GWAS_ROOT}/${DATASET}/hg19/${DATASET}.hg19.tsv.gz"
WORK="${PE_GWAS_ROOT}/${DATASET}/work/rsid_bridge"
HG38_VCF="${WORK}/${DATASET}.hg38.raw.vcf"
SORTED_VCF="${WORK}/${DATASET}.hg38.sorted.vcf.gz"
NORM_VCF="${WORK}/${DATASET}.hg38.normalized.vcf.gz"
ANNOTATED_VCF="${WORK}/${DATASET}.hg38.dbsnp.vcf.gz"
HG38_TSV="${WORK}/${DATASET}.hg38.rsid.tsv.gz"
INDEX="${PROCESSED_DIR}/${DATASET}.rsid_hg19_index.pkl"
THREADS="${THREADS:-4}"

require_executable "Python" "${PYTHON_DEFAULT}"
require_executable "bcftools" "${BCFTOOLS}"
require_executable "bgzip" "${BGZIP}"
require_executable "tabix" "${TABIX}"
for p in "${SRC}" "${HG19_FASTA}" "${HG19_FASTA_FAI}" "${HG38_FASTA}" "${HG38_FASTA_FAI}" "${DBSNP_HG38_VCF}" "${DBSNP_HG38_VCF}.tbi"; do
  require_path "Figshare preprocessing input" "${p}"
done
mkdir -p "${WORK}" "$(dirname "${OUT}")"

if [[ ! -s "${HG38_VCF}" ]]; then
  "${PYTHON_DEFAULT}" "${SCRIPT_DIR}/figshare_to_hg38_vcf.py" \
    --metal "${SRC}" --fasta "${HG38_FASTA}" --fai "${HG38_FASTA_FAI}" \
    --out-vcf "${HG38_VCF}" --qc-json "${WORK}/figshare_hg38_input.qc.json"
fi
if [[ ! -s "${SORTED_VCF}" ]]; then
  "${BCFTOOLS}" sort --temp-dir "${WORK}/bcftools-sort.XXXXXX" -Oz -o "${SORTED_VCF}.tmp" "${HG38_VCF}"
  mv "${SORTED_VCF}.tmp" "${SORTED_VCF}"
  "${TABIX}" -f -p vcf "${SORTED_VCF}"
fi
if [[ ! -s "${NORM_VCF}" ]]; then
  "${BCFTOOLS}" norm -f "${HG38_FASTA}" -c e -Oz -o "${NORM_VCF}.tmp" "${SORTED_VCF}"
  mv "${NORM_VCF}.tmp" "${NORM_VCF}"
  "${TABIX}" -f -p vcf "${NORM_VCF}"
fi
if [[ ! -s "${ANNOTATED_VCF}" ]]; then
  "${BCFTOOLS}" annotate --threads "${THREADS}" -a "${DBSNP_HG38_VCF}" -c ID \
    -Oz -o "${ANNOTATED_VCF}.tmp" "${NORM_VCF}"
  mv "${ANNOTATED_VCF}.tmp" "${ANNOTATED_VCF}"
  "${TABIX}" -f -p vcf "${ANNOTATED_VCF}"
fi
if [[ ! -s "${HG38_TSV}" ]]; then
  "${PYTHON_DEFAULT}" "${SCRIPT_DIR}/figshare_annotated_vcf_to_gwas.py" \
    --vcf "${ANNOTATED_VCF}" --out "${HG38_TSV}" --bcftools "${BCFTOOLS}" \
    --qc-json "${WORK}/figshare_dbsnp_bridge.qc.json"
fi
# 2026-09-10 (audit B-8): the rsID index is valid only for the exact
# HG38_TSV bytes; a stale index would silently drop refreshed variants.
# Same input-signature gate as run_convert_pe.sh.
INDEX_MANIFEST="${INDEX}.inputs.sha256"
RSID_SET="${INDEX}.requested_rsids.txt"
INDEX_SIGNATURE=$(sha256sum "${HG38_TSV}")
CURRENT_SIGNATURE=$([[ -r "${INDEX_MANIFEST}" ]] && cat "${INDEX_MANIFEST}" || true)
if [[ ! -s "${INDEX}" || "${CURRENT_SIGNATURE}" != "${INDEX_SIGNATURE}" || "${FORCE:-0}" == 1 ]]; then
  rsid_set_signature=$([[ -r "${RSID_SET}.inputs.sha256" ]] && cat "${RSID_SET}.inputs.sha256" || true)
  if [[ "${rsid_set_signature}" != "${INDEX_SIGNATURE}" ]]; then
    rm -f "${RSID_SET}" "${RSID_SET}.inputs.sha256"
  fi
  "${PYTHON_DEFAULT}" "${SCRIPT_DIR}/build_rsid_index.py" --gwas "${HG38_TSV}" \
    --rsid-set "${RSID_SET}" --vcf-dir "${HG19_VCF_DIR}" \
    --vcf-pattern 'ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz' \
    --out "${INDEX}" --bcftools "${BCFTOOLS}" --threads "${THREADS}"
  printf '%s\n' "${INDEX_SIGNATURE}" > "${RSID_SET}.inputs.sha256"
  printf '%s\n' "${INDEX_SIGNATURE}" > "${INDEX_MANIFEST}"
fi
"${PYTHON_DEFAULT}" "${SCRIPT_DIR}/convert_gwas_hg38_to_hg19.py" \
  --gwas "${HG38_TSV}" --index "${INDEX}" --out "${OUT}" \
  --qc-json "${OUT}.qc.json" --bgzip "${BGZIP}" --tabix "${TABIX}"
bash "${SCRIPT_DIR}/verify_conversion.sh" "${OUT}" "${HG38_TSV}" \
  "${HG19_VCF_DIR}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" "${BCFTOOLS}" \
  | tee "$(dirname "${OUT}")/conversion_qc.log"
