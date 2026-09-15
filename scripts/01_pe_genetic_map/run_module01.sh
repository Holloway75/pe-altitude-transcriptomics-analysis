#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

usage() {
  cat <<'EOF'
Usage: run_module01.sh COMMAND [arguments]
  status --analysis-id ID     show registered scientific context
  validate --analysis-id ID DIRECTORY
                              validate the six formal tables
  publish --analysis-id ID STAGING
                              validate and atomically publish

RETIRED (audit B-1, 2026-09-10): the former `run`/`resume` entry delegated to
execute.py, which can never complete (full MR-JTI M1.6 is intentionally
disabled there). The real executable chain that produced the frozen results is
the per-step entry points documented in README.md ("Execution order"):
  preflight.py -> prepare_gwas.py -> run_single_gwas_spredixcan.py ->
  build_all_gene_smultixcan_covariance.py / run_all_gene_smultixcan.py ->
  validate_existing_all_gene_smultixcan.py -> assemble_m1_downstream.py
  (+ continue_m1_downstream.py); MR-JTI via mrjti_minimal.py.

task-ledger.tsv is a synthetic placeholder (see TASK_LEDGER_PLACEHOLDER_NOTE.md
in the published release); it is not execution evidence.
EOF
}
cmd="${1:-}"; shift || true
case "$cmd" in
  run|resume)
    usage >&2
    exit 3
    ;;
  status) [[ "${1:-}" == "--analysis-id" && $# == 2 ]] || { usage >&2; exit 2; }; exec python3 "$SCRIPT_DIR/module01.py" status --analysis-id "$2" ;;
  validate) [[ "${1:-}" == "--analysis-id" && $# == 3 ]] || { usage >&2; exit 2; }; exec python3 "$SCRIPT_DIR/module01.py" validate --analysis-id "$2" --directory "$3" ;;
  publish) [[ "${1:-}" == "--analysis-id" && $# == 3 ]] || { usage >&2; exit 2; }; exec python3 "$SCRIPT_DIR/module01.py" publish --analysis-id "$2" --staging "$3" ;;
  -h|--help) usage ;;
  *) usage >&2; exit 2 ;;
esac
