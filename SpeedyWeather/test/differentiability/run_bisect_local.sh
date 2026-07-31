#!/usr/bin/env bash
# Sequential bisection sweep on an x86 CPU node (interactive / no scheduler).
# Each variant runs in its OWN julia process so an Enzyme codegen abort (exit 134,
# "unhandled accumulate with partial sizes") only kills that one variant — the
# sweep continues and prints a PASS/ABORT table at the end.
#
# Run from the repo root. Set up the env ONCE first:
#   julia --project=SpeedyWeather/test/differentiability \
#         SpeedyWeather/test/differentiability/setup_diff_env.jl
# Then:
#   SpeedyWeather/test/differentiability/run_bisect_local.sh              # all variants
#   SpeedyWeather/test/differentiability/run_bisect_local.sh A2_param C2_paramtend   # subset
#
# Config via env: JULIA (default: julia). Must be Julia 1.12 on x86 to see the abort.
set -u
JULIA="${JULIA:-julia}"
PROJ="SpeedyWeather/test/differentiability"
DRIVER="$PROJ/bisect_pwet_param.jl"
LOGDIR="$PROJ/_logs"
mkdir -p "$LOGDIR"

ALL=(A1_state A2_param B1_dynonly C1_dyntend C2_paramtend C3_implicit \
     C4_diffusion C5_updateprog C6_transform C7_ocean C8_land)
VARIANTS=("${@:-}")
[ -z "${VARIANTS[*]}" ] && VARIANTS=("${ALL[@]}")

declare -A RESULT
for v in "${VARIANTS[@]}"; do
    log="$LOGDIR/bisect_${v}.log"
    echo ">>> $v  (log: $log)"
    "$JULIA" --project="$PROJ" "$DRIVER" "$v" >"$log" 2>&1
    code=$?
    if [ $code -eq 0 ] && grep -q "RESULT: SUCCESS" "$log"; then
        RESULT[$v]="PASS"
    else
        RESULT[$v]="ABORT/ERR (exit $code)"
    fi
    echo "    -> ${RESULT[$v]}"
done

echo
echo "==================== BISECTION SUMMARY ===================="
printf "%-16s %s\n" "VARIANT" "OUTCOME"
for v in "${VARIANTS[@]}"; do printf "%-16s %s\n" "$v" "${RESULT[$v]}"; done
echo "=========================================================="
echo "A RED/ABORT variant = the culprit component. See BISECTION_x86.md."
