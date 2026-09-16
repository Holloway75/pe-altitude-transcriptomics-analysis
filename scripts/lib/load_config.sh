#!/usr/bin/env bash
# shellcheck shell=bash

if [[ -n "${PROJECT_CONFIG_LOADED:-}" ]]; then
  return 0 2>/dev/null || exit 0
fi

LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${LIB_DIR}/../.." && pwd)}"
CONFIG_DIR="${PROJECT_CONFIG_DIR:-${PROJECT_ROOT}/config}"

had_allexport=0
[[ $- == *a* ]] && had_allexport=1
set -a
for config_file in "${CONFIG_DIR}/paths.env" "${CONFIG_DIR}/analysis.env"; do
  if [[ ! -r "${config_file}" ]]; then
    echo "ERROR: configuration file is missing: ${config_file}" >&2
    return 2 2>/dev/null || exit 2
  fi
  # shellcheck disable=SC1090
  source "${config_file}"
done
((had_allexport)) || set +a

export PROJECT_ROOT CONFIG_DIR
PROJECT_CONFIG_LOADED=1
export -n PROJECT_CONFIG_LOADED 2>/dev/null || true

require_path() {
  local label=$1 path=$2
  if [[ ! -e "${path}" ]]; then
    echo "ERROR: ${label} does not exist: ${path}" >&2
    return 1
  fi
}

require_executable() {
  local label=$1 executable=$2
  if [[ "${executable}" == */* ]]; then
    [[ -x "${executable}" ]] || { echo "ERROR: ${label} is not executable: ${executable}" >&2; return 1; }
  else
    command -v "${executable}" >/dev/null 2>&1 || { echo "ERROR: ${label} is unavailable: ${executable}" >&2; return 1; }
  fi
}
