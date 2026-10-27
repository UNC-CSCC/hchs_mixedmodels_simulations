#!/bin/bash
#SBATCH -p general
#SBATCH --account=rc_cai_pi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.err

# Measure one replicate with srun, set BLOCK so each job is ~1 hour, then
# sbatch --array=1-<ceil(N_REPS/BLOCK)> run_stata_chunked.slurm <do_name> <out_template> [args...]
set -u
BLOCK=120
N_REPS=1000

DO_NAME="$1"; OUT_TEMPLATE="$2"; shift 2
REPO="/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK + 1 ))
LAST=$(( SLURM_ARRAY_TASK_ID * BLOCK )); [ "$LAST" -gt "$N_REPS" ] && LAST="$N_REPS"

module purge; module add stata/19
export STATATMP="${REPO}/.statatmp"
WORK="${STATATMP}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORK"; cd "$WORK"

for SIM in $(seq "$FIRST" "$LAST"); do
    OUT="${REPO}/${OUT_TEMPLATE//\{SIM\}/$SIM}"
    [ -s "$OUT" ] && continue
    timeout 120m stata-se -b do "${REPO}/scripts/stata/${DO_NAME}.do" "$SIM" "$@"
    [ -s "$OUT" ] && echo "rep $SIM done" || echo "rep $SIM FAILED" >&2
done

cd "$REPO"; rm -rf "$WORK"
