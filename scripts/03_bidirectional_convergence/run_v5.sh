#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/scripts/lib/load_config.sh"

ANALYSIS_ID="m3_exploratory_shared_programs_v5_20260729"
MODE="${1:-all}"
CONFIG="${PROJECT_ROOT}/config/module03_v5.json"
MODULE2_RELEASE="${PROCESSED_DIR}/02_altitude_response/releases/m2_module3_input_freeze_v1_20260728"
V4_REGISTRY="${PROJECT_ROOT}/work/03_bidirectional_convergence/v4_registry_recovery_current_hgnc/data_processed/registry"
LOCKED_DIR="${PROCESSED_DIR}/03_bidirectional_convergence/releases/${ANALYSIS_ID}"
RESULT_DIR="${RESULTS_DIR}/03_bidirectional_convergence/analyses/${ANALYSIS_ID}"
PYTHON_DEG="${CONDA_ROOT}/envs/deg/bin/python"
RSCRIPT_DEG="${CONDA_ROOT}/envs/deg/bin/Rscript"
LIMMA_PRIOR="${PROJECT_ROOT}/work/03_bidirectional_convergence/input_recovery_20260809/GSE103927_limma_empirical_bayes_prior.tsv"

run_overlap() {
  "${PYTHON_DEG}" "${SCRIPT_DIR}/run_v5_overlap.py" \
    --project-root "${PROJECT_ROOT}" \
    --config "${CONFIG}" \
    --v4-registry-dir "${V4_REGISTRY}" \
    --module2-release "${MODULE2_RELEASE}" \
    --limma-prior "${LIMMA_PRIOR}" \
    --locked-dir "${LOCKED_DIR}" \
    --result-dir "${RESULT_DIR}"
}

run_signatures() {
  "${RSCRIPT_DEG}" "${SCRIPT_DIR}/run_compact_signatures.R" \
    --config "${CONFIG}" \
    --locked-dir "${LOCKED_DIR}" \
    --result-dir "${RESULT_DIR}" \
    --module2-release "${MODULE2_RELEASE}" \
    --hallmark-gmt "${MSIGDB_ROOT}/h.all.v2024.1.Hs.symbols.gmt" \
    --reactome-gmt "${MSIGDB_ROOT}/c2.cp.reactome.v2024.1.Hs.symbols.gmt"
}

run_validate() {
  "${PYTHON_DEG}" "${SCRIPT_DIR}/validate_v5.py" \
    --locked-dir "${LOCKED_DIR}" \
    --result-dir "${RESULT_DIR}"
}

case "${MODE}" in
  overlap) run_overlap ;;
  signatures) run_signatures ;;
  validate) run_validate ;;
  all)
    run_overlap
    run_signatures
    run_validate
    ;;
  *)
    echo "Usage: $0 {overlap|signatures|validate|all}" >&2
    exit 2
    ;;
esac
