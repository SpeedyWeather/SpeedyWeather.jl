#!/usr/bin/env bash
# Benchmark any revision of SpeedyWeather.jl with the debug mode of manual_benchmarking.jl,
# always using the benchmark harness of *this* checkout so that every revision is measured
# identically. The benchmarked tree is never modified: the harness is copied into a separate
# environment directory whose Project.toml points at the tree via [sources] (Julia ≥ 1.11).
#
# Usage: run_debug_benchmark.sh <revision|tree> <output.json> [cpu|gpu|amdgpu]
#   revision  any git revision; benchmarked in a detached worktree under $WORKDIR/trees
#   tree      an existing checkout, benchmarked as is (used by bisect_step.sh)
#
# WORKDIR is $SPEEDY_BENCH_WORKDIR, default ${TMPDIR:-/tmp}/speedyweather-benchmark-regression.
# The julia log is written next to the output as <output>.log.
set -euo pipefail

if [[ $# -lt 2 ]]; then
    sed -n '2,13p' "$0"
    exit 2
fi

TARGET=$1
OUTPUT=$(mkdir -p "$(dirname "$2")" && cd "$(dirname "$2")" && pwd)/$(basename "$2")
ARCH=${3:-cpu}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
HARNESS="$REPO/SpeedyWeather/benchmark"
WORKDIR=${SPEEDY_BENCH_WORKDIR:-${TMPDIR:-/tmp}/speedyweather-benchmark-regression}

julia_minor=$(julia --startup-file=no -e 'print(VERSION.minor)')
if [[ $(julia --startup-file=no -e 'print(VERSION.major)') -eq 1 && $julia_minor -lt 11 ]]; then
    echo "Julia ≥ 1.11 is required ([sources] is ignored before), found 1.$julia_minor" >&2
    exit 1
fi

if [[ -d $TARGET ]]; then
    TREE=$(cd "$TARGET" && pwd)
else
    SHA=$(git -C "$REPO" rev-parse --verify "$TARGET^{commit}")
    TREE="$WORKDIR/trees/${SHA:0:10}"
    if [[ ! -d $TREE ]]; then
        git -C "$REPO" worktree add --detach --quiet "$TREE" "$SHA"
    fi
fi
SHA=$(git -C "$TREE" rev-parse HEAD)

# slim environment with only what the debug mode loads, SpeedyWeather from the tree
ENV_DIR="$WORKDIR/envs/${SHA:0:10}-$ARCH"
mkdir -p "$ENV_DIR"
cp "$HARNESS"/{manual_benchmarking.jl,benchmark_suite.jl,define_benchmarks.jl} "$ENV_DIR/"
{
    echo "[deps]"
    grep -E '^(BenchmarkTools|Dates|JSON3|Printf|SpeedyWeather) *= *"' "$HARNESS/Project.toml"
    case $ARCH in
        gpu) echo 'CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"' ;;
        amdgpu) echo 'AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"' ;;
    esac
    echo
    echo "[sources]"
    echo "SpeedyWeather = {path = \"$TREE/SpeedyWeather\"}"
} > "$ENV_DIR/Project.toml"

echo "Benchmarking ${SHA:0:10} ($(git -C "$TREE" log -1 --format=%s)) on $ARCH → $OUTPUT"
{
    julia --project="$ENV_DIR" --startup-file=no -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
    julia --project="$ENV_DIR" --startup-file=no "$ENV_DIR/manual_benchmarking.jl" "$ARCH" --debug --output="$OUTPUT"
} > "$OUTPUT.log" 2>&1 || {
    echo "Benchmark of ${SHA:0:10} failed, last lines of $OUTPUT.log:" >&2
    tail -n 30 "$OUTPUT.log" >&2
    exit 1
}

# check that the tree's packages were benchmarked (not e.g. a registered release) and record the revision
julia --startup-file=no "$SCRIPT_DIR/regression.jl" stamp "$OUTPUT" --tree "$TREE" --revision "$SHA"
grep -E '^\| PrimitiveWetModel' "$OUTPUT.log" || true
