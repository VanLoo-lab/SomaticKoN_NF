#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------
# ensure_conda_env
#
# Usage:
#   ensure_conda_env <env_name> <env_yml> [conda_base] [env_root]
#
# Example:
#   ensure_conda_env mutect2 \
#     ${projectDir}/envs/mutect2.yml \
#     /risapps/rhel8/miniforge3/24.5.0-0 \
#     /rsrch6/home/genetics/zzhang18/conda-envs
#
# Behavior:
#   - installs env only if missing
#   - installs into env_root/env_name
#   - activates the env
#   - safe for repeated calls
# ---------------------------------------------------------
ensure_conda_env() {

  if [[ $# -lt 2 ]]; then
    echo "Usage: ensure_conda_env <env_name> <env_yml> [conda_base] [env_root]" >&2
    return 2
  fi

  local ENV_NAME="$1"
  local ENV_YML="$2"
  local ENV_ROOT="${3:-$HOME/conda-envs}"
  local ENV_PATH="${ENV_ROOT}/${ENV_NAME}"

  if [[ ! -f "${ENV_YML}" ]]; then
    echo "[ERROR] Conda env yml not found: ${ENV_YML}" >&2
    return 2
  fi

  # init conda
  #source "${CONDA_BASE}/etc/profile.d/conda.sh"
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
