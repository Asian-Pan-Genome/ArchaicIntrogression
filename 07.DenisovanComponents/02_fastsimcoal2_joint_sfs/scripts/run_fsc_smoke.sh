#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || {
  echo "Usage: FSC_BIN=/path/to/fsc28 $0 MODEL_PREFIX" >&2
  exit 2
}

MODEL=$1
: "${FSC_BIN:?Set FSC_BIN to the fastsimcoal2 v2.8 executable}"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
MODEL_DIR="${PROJECT_DIR}/models"
OBS_SFS="${PROJECT_DIR}/data/observed_sfs/YRI10_CHB20_DEN_1SNPper100kb_MSFS.obs"
RUN_DIR="${PROJECT_DIR}/work/${MODEL}"

[[ -f "${MODEL_DIR}/${MODEL}.tpl" ]] || { echo "Missing ${MODEL}.tpl" >&2; exit 3; }
[[ -f "${MODEL_DIR}/${MODEL}.est" ]] || { echo "Missing ${MODEL}.est" >&2; exit 3; }
[[ ! -e "${RUN_DIR}" ]] || { echo "Refusing to overwrite ${RUN_DIR}" >&2; exit 4; }

mkdir -p "${RUN_DIR}"
cp "${MODEL_DIR}/${MODEL}.tpl" "${RUN_DIR}/"
cp "${MODEL_DIR}/${MODEL}.est" "${RUN_DIR}/"
cp "${OBS_SFS}" "${RUN_DIR}/${MODEL}_MSFS.obs"

cd "${RUN_DIR}"
"${FSC_BIN}" \
  -t "${MODEL}.tpl" \
  -e "${MODEL}.est" \
  -n10000 \
  -L10 \
  -m --multiSFS -M -C1 -c1 -0
