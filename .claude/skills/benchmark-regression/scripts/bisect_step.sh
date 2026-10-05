#!/usr/bin/env bash
# `git bisect run` step: benchmark the bisect worktree's HEAD and classify it against the
# reference results as good (exit 0) or bad (exit 1). Revisions that cannot be benchmarked
# (build or run failures) are skipped (exit 125).
#
# Usage, from the bisect worktree:
#   git bisect run bisect_step.sh <reference.json[,...]> <cutoff> [all|default|matrix] [cpu|gpu|amdgpu]
set -uo pipefail

REFERENCE=$1
CUTOFF=$2
GROUP=${3:-all}
ARCH=${4:-cpu}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORKDIR=${SPEEDY_BENCH_WORKDIR:-${TMPDIR:-/tmp}/speedyweather-benchmark-regression}
TREE=$(git rev-parse --show-toplevel)
SHA=$(git rev-parse HEAD)
OUTPUT="$WORKDIR/results/bisect-${SHA:0:10}.json"

# reuse a measurement of this revision from an earlier bisect step or run
if [[ ! -f $OUTPUT ]]; then
    "$SCRIPT_DIR/run_debug_benchmark.sh" "$TREE" "$OUTPUT" "$ARCH" || exit 125
fi

julia --startup-file=no "$SCRIPT_DIR/regression.jl" verdict --reference "$REFERENCE" --candidate "$OUTPUT" \
    --group "$GROUP" --cutoff "$CUTOFF"
status=$?
[[ $status -le 1 ]] && exit $status
exit 125
