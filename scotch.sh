#!/usr/bin/env bash
set -euo pipefail
SCOTCH_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "${SCOTCH_ROOT}/src/run_pipeline.py" "$@"
