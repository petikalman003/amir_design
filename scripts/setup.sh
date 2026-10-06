#!/usr/bin/env bash
# Reproduce the Python, model, sequence-screening, and BLAST dependencies.
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
export UV_CACHE_DIR="${UV_CACHE_DIR:-/workspace/.cache/uv}"
export UV_PYTHON_INSTALL_DIR="${UV_PYTHON_INSTALL_DIR:-/workspace/.runtime/uv-python}"
export UV_PYTHON_BIN_DIR="${UV_PYTHON_BIN_DIR:-/workspace/.runtime/bin}"
export AMIR_RUNTIME_DIR="${AMIR_RUNTIME_DIR:-/workspace/.runtime}"
export PATH="$AMIR_RUNTIME_DIR/bin:$PATH"
if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--references" ) ]]; then
  printf 'Usage: bash scripts/setup.sh [--references]\n' >&2
  exit 2
fi
command -v uv >/dev/null || { printf 'Install uv first: https://docs.astral.sh/uv/\n' >&2; exit 1; }
uv python install 3.11.12
uv sync --frozen --extra dev --python 3.11.12
bash scripts/install_r.sh

clone_dependency() {
  local remote="$1" destination="$2" commit="$3"
  if [[ ! -e "$destination" ]]; then
    git clone --no-checkout "$remote" "$destination"
    git -C "$destination" checkout --detach "$commit"
  fi
  if [[ "$(git -C "$destination" rev-parse HEAD)" != "$commit" ]]; then
    printf 'Existing dependency %s is not at the required commit; preserved without changes.\n' "$destination" >&2
    exit 1
  fi
  git -C "$destination" diff --quiet HEAD -- || {
    printf 'Existing dependency %s has tracked edits; preserved without changes.\n' "$destination" >&2
    exit 1
  }
}
mkdir -p tools
clone_dependency https://github.com/ParkSJ-91/shRNAI.git tools/shrnai 69f8576e027b8036eaeb2d039c40a2ccdbabf373
clone_dependency https://github.com/tacazares/SeedMatchR.git tools/seedmatchr 2fc0c4fba31e7306377488ac8b4c81d6445df5e2

.venv/bin/python - <<'PY'
from pathlib import Path
from amir_design.scoring import ShRNAIScorer
from amir_design.offtargets import SeedMatchRunner
ShRNAIScorer(Path('tools/shrnai'))
SeedMatchRunner(Path('tools/seedmatchr'))._check_runtime()
print('Verified shRNAI models, SeedMatchR source, and native runtime.')
PY
if [[ "${1:-}" == "--references" ]]; then
  .venv/bin/python -m amir_design.reference_cli --species both --root data/references
fi
printf 'Dependencies ready. The published miR-E insert template is included; use --references to prepare species data, then bash scripts/start.sh.\n'
