#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
export AMIR_RUNTIME_DIR="${AMIR_RUNTIME_DIR:-/workspace/.runtime}"
export PATH="$AMIR_RUNTIME_DIR/bin:$PATH"
export TF_CPP_MIN_LOG_LEVEL="${TF_CPP_MIN_LOG_LEVEL:-2}"
export TF_NUM_INTRAOP_THREADS="${TF_NUM_INTRAOP_THREADS:-2}"
export TF_NUM_INTEROP_THREADS="${TF_NUM_INTEROP_THREADS:-1}"
exec .venv/bin/python -m streamlit run app.py \
  --server.address 0.0.0.0 --server.port "${AMIR_PORT:-8501}" \
  --server.headless true --browser.gatherUsageStats false
