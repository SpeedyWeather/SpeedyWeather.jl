#!/usr/bin/env bash
# SLURM array job — the parameter-AD bisection sweep on an x86 cluster (Julia 1.12).
# One array task per variant, each in its own job so an Enzyme codegen abort
# (exit 134, "unhandled accumulate with partial sizes") isolates that component.
#
# ── BEFORE SUBMITTING ─────────────────────────────────────────────────────────
#   1. FIRST, INSPECT YOUR CURRENT LOCAL SLURM SCRIPTS to get the up-to-date
#      cluster config — the --qos / --account / --partition / --mail-user /
#      `module load julia/...` values are cluster- and project-specific AND CHANGE
#      OVER TIME. Do NOT trust the example values below; verify against a recent
#      script you actually use now. Places to look (may be stale — check dates):
#          ~/Nextcloud/*/submit*.sh   (e.g. NeuralQG3/, mcbb-climate-project/ —
#          seen historically: qos=gpushort/priority, account=flai/brasnet, PIK;
#          those loaded OLD Julia 1.9/1.3 — you need a 1.12 x86 build.)
#      Copy the CURRENT values into the EDIT-marked lines below, and pick a
#      Julia 1.12 x86 module (or juliaup/local 1.12); confirm `julia --version`.
#   2. Set up the env ONCE (NOT inside the array — avoids Manifest races):
#        julia --project=SpeedyWeather/test/differentiability \
#              SpeedyWeather/test/differentiability/setup_diff_env.jl
#   3. Submit FROM THE REPO ROOT:
#        sbatch SpeedyWeather/test/differentiability/bisect_pwet_param.slurm
# ──────────────────────────────────────────────────────────────────────────────
#
#SBATCH --job-name=enz112bisect
#SBATCH --array=0-10                 # 11 variants (see VARIANTS below)
#SBATCH --time=00:40:00              # each variant aborts fast or compiles in <~30min
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=SpeedyWeather/test/differentiability/_logs/bisect_%A_%a.log
### EDIT for your cluster (copy from a sibling submit script — see note above):
###SBATCH --qos=short
###SBATCH --account=CHANGEME
###SBATCH --partition=CHANGEME
###SBATCH --mail-type=END
###SBATCH --mail-user=CHANGEME

set -u
# --- Julia 1.12 x86 access: EDIT for your cluster (see step 1) -----------------
# module purge
# module load julia/1.12        # <-- must be a 1.12 x86 build (NOT the 1.9/1.3 in the examples)
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
