#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-python3}

"${PYTHON_BIN}" "${SCRIPT_DIR}/plot_models.py" \
  --layout paired \
  --output-dir "${PROJECT_DIR}/figures"
