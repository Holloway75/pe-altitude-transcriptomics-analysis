#!/usr/bin/env bash
# Convert GWAS Catalog harmonized PE files from hg38 to hg19 using the shared manifest.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../lib/load_config.sh
source "${SCRIPT_DIR}/../../lib/load_config.sh"

MANIFEST="${PE_GWAS_MANIFEST:-${CONFIG_DIR}/pe_gwas.tsv}"
INDEX="${RSID_HG19_INDEX:-${PROCESSED_DIR}/rsid_hg19_index.pkl}"
RSID_SET="${RSID_HG19_SET:-${INDEX}.requested_rsids.txt}"
VCF_PATTERN="${HG19_VCF_PATTERN:-}"
if [[ -z "${VCF_PATTERN}" ]]; then
  # Keep the Python format braces outside a ${var:-default} expansion: the
  # first `}` in {chr} would otherwise terminate the shell expansion early.
  VCF_PATTERN='ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
fi
THREADS="${THREADS:-4}"
FORCE="${FORCE:-0}"

require_executable "Python" "${PYTHON_DEFAULT}"
require_executable "bcftools" "${BCFTOOLS}"
require_executable "bgzip" "${BGZIP}"
require_executable "tabix" "${TABIX}"
require_path "PE GWAS manifest" "${MANIFEST}"
require_path "hg19 VCF directory" "${HG19_VCF_DIR}"

mapfile -t rows < <(awk -F '\t' 'NR > 1 && $7 == "gwas_catalog_harmonized" {print $1 "\t" $8 "\t" $9}' "${MANIFEST}")
if [[ ${#rows[@]} -eq 0 ]]; then
  echo "ERROR: no gwas_catalog_harmonized rows in ${MANIFEST}" >&2
  exit 1
fi

inputs=()
dataset_ids=()
for row in "${rows[@]}"; do
  IFS=$'\t' read -r dataset input_rel output_rel <<< "${row}"
  input="${PE_GWAS_ROOT}/${input_rel}"
  require_path "${dataset} hg38 input" "${input}"
  inputs+=("${input}")
  dataset_ids+=("${dataset}")
done

index_manifest="${INDEX}.inputs.sha256"
# An rsID index is valid only for the exact source bytes, not merely the same
# dataset labels.  Otherwise a refreshed GWAS can silently reuse a stale index.
expected_signature=$(sha256sum "${inputs[@]}")
current_signature=$([[ -r "${index_manifest}" ]] && cat "${index_manifest}" || true)
if [[ ! -s "${INDEX}" || "${current_signature}" != "${expected_signature}" || "${FORCE}" == 1 ]]; then
  rsid_set_manifest="${RSID_SET}.inputs.sha256"
  rsid_set_signature=$([[ -r "${rsid_set_manifest}" ]] && cat "${rsid_set_manifest}" || true)
  if [[ "${rsid_set_signature}" != "${expected_signature}" ]]; then
    rm -f "${RSID_SET}" "${rsid_set_manifest}"
  fi
  "${PYTHON_DEFAULT}" "${SCRIPT_DIR}/build_rsid_index.py" \
    --gwas "${inputs[@]}" \
    --rsid-set "${RSID_SET}" \
    --vcf-dir "${HG19_VCF_DIR}" \
    --vcf-pattern "${VCF_PATTERN}" \
    --out "${INDEX}" \
    --bcftools "${BCFTOOLS}" --threads "${THREADS}"
  printf '%s\n' "${expected_signature}" > "${rsid_set_manifest}"
  printf '%s\n' "${expected_signature}" > "${index_manifest}"
fi

for row in "${rows[@]}"; do
  IFS=$'\t' read -r dataset input_rel output_rel <<< "${row}"
  input="${PE_GWAS_ROOT}/${input_rel}"
  output="${PE_GWAS_ROOT}/${output_rel}"
  qc="$(dirname "${output}")/conversion_qc.log"
  mkdir -p "$(dirname "${output}")"
  if [[ -s "${output}" && -s "${output}.tbi" && "${FORCE}" != 1 ]]; then
    echo "[skip] ${dataset}: hg19 output already exists"
    continue
  fi
  : > "${qc}"
  "${PYTHON_DEFAULT}" "${SCRIPT_DIR}/convert_gwas_hg38_to_hg19.py" \
    --gwas "${input}" --index "${INDEX}" --out "${output}" \
    --qc-json "${output}.qc.json" \
    --bgzip "${BGZIP}" --tabix "${TABIX}" 2>&1 | tee -a "${qc}"
  bash "${SCRIPT_DIR}/verify_conversion.sh" \
    "${output}" "${input}" "${HG19_VCF_DIR}/${VCF_PATTERN}" "${BCFTOOLS}" 2>&1 | tee -a "${qc}"
done

echo "NOTE: FIGSHARE_22680904_v2 uses preprocess_mode=figshare_metal and requires its dedicated normalization step."
