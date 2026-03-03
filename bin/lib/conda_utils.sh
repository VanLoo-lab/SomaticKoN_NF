#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------
# ensure_conda_env
#
# Usage:
#   ensure_conda_env <env_name> <env_yml> [env_root] [pipeline_root]
#
# Example:
#   ensure_conda_env mutect2 \
#     ${projectDir}/envs/mutect2.yml \
#     /rsrch6/home/genetics/zzhang18/conda-envs \
#     ${projectDir}
#
# Behavior:
#   - installs env only if missing
#   - installs into env_root/env_name
#   - activates the env
#   - safe for repeated calls
# ---------------------------------------------------------
ensure_conda_env() {

  if [[ $# -lt 2 ]]; then
    echo "Usage: ensure_conda_env <env_name> <env_yml> [env_root] [pipeline_root]" >&2
    return 2
  fi

  local ENV_NAME="$1"
  local ENV_YML="$2"
  local ENV_ROOT="${3:-}"
  local PIPELINE_ROOT="${4:-}"
  local ENV_PATH=""
  local -a CANDIDATE_ROOTS=()

  if [[ ! -f "${ENV_YML}" ]]; then
    echo "[ERROR] Conda env yml not found: ${ENV_YML}" >&2
    return 2
  fi

  # Explicit env root wins. Otherwise probe pipeline-local first, then home.
  if [[ -n "${ENV_ROOT}" ]]; then
    CANDIDATE_ROOTS+=("${ENV_ROOT}")
  else
    if [[ -n "${PIPELINE_ROOT}" ]]; then
      CANDIDATE_ROOTS+=("${PIPELINE_ROOT}/conda_envs")
    fi
    CANDIDATE_ROOTS+=("${HOME}/.conda/envs" "${HOME}/conda-envs")
  fi

  if [[ ${#CANDIDATE_ROOTS[@]} -eq 0 ]]; then
    echo "[ERROR] No conda env root candidates available" >&2
    return 2
  fi

  for root in "${CANDIDATE_ROOTS[@]}"; do
    if [[ -d "${root}/${ENV_NAME}" ]]; then
      ENV_ROOT="${root}"
      ENV_PATH="${root}/${ENV_NAME}"
      break
    fi
  done

  if [[ -z "${ENV_PATH}" ]]; then
    ENV_ROOT="${CANDIDATE_ROOTS[0]}"
    ENV_PATH="${ENV_ROOT}/${ENV_NAME}"
  fi

  echo "[INFO] Conda env root candidates: ${CANDIDATE_ROOTS[*]}"
  echo "[INFO] Selected conda env root: ${ENV_ROOT}"

  # init conda
  eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

  mkdir -p "${ENV_ROOT}"


  # ---- create env if missing ----
  if [[ ! -d "${ENV_PATH}" ]]; then
    echo "[INFO] Creating conda env '${ENV_NAME}' at ${ENV_PATH}"
    
      ## ---- lock to avoid race conditions on HPC ----
    local LOCK="${ENV_PATH}.lock"
    exec 9>"${LOCK}"
    flock -n 9 || {
    echo "[INFO] Waiting for conda env lock: ${ENV_NAME}"
    flock 9
    }

    conda env create \
      --prefix "${ENV_PATH}" \
      --file "${ENV_YML}"
  else
    echo "[INFO] Conda env '${ENV_NAME}' already exists at ${ENV_PATH}"
  fi

  # ---- activate ----
  conda activate "${ENV_PATH}"

  # ---- sanity check ----
  if [[ "${CONDA_PREFIX:-}" != "${ENV_PATH}" ]]; then
    echo "[ERROR] Failed to activate conda env: ${ENV_PATH}" >&2
    return 2
  fi

  echo "[INFO] Activated conda env: ${ENV_NAME}"
}
