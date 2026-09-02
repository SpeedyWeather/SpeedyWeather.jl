#!/bin/bash
# Submit the quadrature stability matrix as a SLURM array, one GPU job per case.
#
#   ./healpix_quadrature/submit_matrix.sh <grid> <truncation> <years> <diffusion_hours...>
#
# e.g. the decisive comparison at the configuration that is known to blow up:
#   ./healpix_quadrature/submit_matrix.sh HEALPixGrid 128 10 4
# and the Step 2 diffusion sweep:
#   ./healpix_quadrature/submit_matrix.sh HEALPixGrid 128 10 1 4 24 96
#
# Cases A…E are described in run_case.jl.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GRID="${1:-HEALPixGrid}"
TRUNC="${2:-128}"
YEARS="${3:-10}"
shift 3 || true
DIFFUSIONS=("$@")
[ ${#DIFFUSIONS[@]} -eq 0 ] && DIFFUSIONS=(4)

CASES=(A B C D E)
OUTPUT="$REPO/healpix_quadrature/runs"
mkdir -p "$OUTPUT" "$REPO/healpix_quadrature/logs"

# one line per job, read back by index inside the array task
JOBLIST="$OUTPUT/joblist_${GRID}_T${TRUNC}_${YEARS}y.txt"
: > "$JOBLIST"
for diffusion in "${DIFFUSIONS[@]}"; do
    for case in "${CASES[@]}"; do
        echo "case=$case grid=$GRID trunc=$TRUNC years=$YEARS diffusion_hours=$diffusion" >> "$JOBLIST"
    done
done
N=$(wc -l < "$JOBLIST")
echo "submitting $N jobs from $JOBLIST"

sbatch --array="1-${N}%5" <<EOF
#!/bin/bash
#SBATCH --qos=gpushort
#SBATCH --job-name=quadmatrix
#SBATCH --account=flai
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=$REPO/healpix_quadrature/logs/quad-%A_%a.log
#SBATCH --error=$REPO/healpix_quadrature/logs/quad-%A_%a.log

module load julia
ARGUMENTS=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
echo "task \${SLURM_ARRAY_TASK_ID}: \$ARGUMENTS"
julia --project=$REPO/healpix_quadrature $REPO/healpix_quadrature/run_case.jl \$ARGUMENTS output=$OUTPUT
EOF
