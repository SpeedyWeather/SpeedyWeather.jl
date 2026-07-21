#!/usr/bin/env bash
# SLURM array job — the parameter-AD bisection sweep on an x86 cluster (Julia 1.12).
# One array task per variant, each in its own job so an Enzyme codegen abort
# (exit 134, "unhandled accumulate with partial sizes") isolates that component.
#
# ── BEFORE SUBMITTING ─────────────────────────────────────────────────────────
#   1. EDIT the #SBATCH --partition / --account lines below for your cluster
#      (or override at submit: sbatch --partition=... --account=... ...).
#   2. Make sure JULIA points to a 1.12 x86 build (edit the `module load` /
#      JULIA= line in the body). Confirm with `julia --version`.
#   3. Set up the env ONCE (NOT inside the array — avoids Manifest races):
#        julia --project=SpeedyWeather/test/differentiability \
#              SpeedyWeather/test/differentiability/setup_diff_env.jl
#   4. Submit FROM THE REPO ROOT:
#        sbatch SpeedyWeather/test/differentiability/bisect_pwet_param.slurm
# ──────────────────────────────────────────────────────────────────────────────
#
#SBATCH --job-name=enz112bisect
#SBATCH --array=0-10                 # 11 variants (see VARIANTS below)
#SBATCH --time=00:40:00              # each variant aborts fast or compiles in <~30min
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=SpeedyWeather/test/differentiability/_logs/bisect_%A_%a.log
### EDIT for your cluster:
###SBATCH --partition=CHANGEME
###SBATCH --account=CHANGEME

set -u
# --- Julia 1.12 x86 access: EDIT for your cluster -----------------------------
# module load julia/1.12        # if your cluster uses environment modules
JULIA="${JULIA:-julia}"
# ------------------------------------------------------------------------------

PROJ="SpeedyWeather/test/differentiability"
DRIVER="$PROJ/bisect_pwet_param.jl"
mkdir -p "$PROJ/_logs"

VARIANTS=(A1_state A2_param B1_dynonly C1_dyntend C2_paramtend C3_implicit \
          C4_diffusion C5_updateprog C6_transform C7_ocean C8_land)
V="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

echo "SLURM task $SLURM_ARRAY_TASK_ID -> variant $V on $(hostname)"
"$JULIA" --version
srun "$JULIA" --project="$PROJ" "$DRIVER" "$V"
code=$?
echo "variant $V finished with exit code $code (134 = Enzyme abort = culprit)"
exit $code
