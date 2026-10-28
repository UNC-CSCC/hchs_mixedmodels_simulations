#!/bin/bash
#SBATCH -p general
#SBATCH --account=rc_cai_pi
#SBATCH --job-name=scen5_mi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.err

set -u
BLOCK=9
N_REPS=1000
VARIANT="$1"; POP="$2"

OUT_DIR="${SLURM_SUBMIT_DIR}/data/estimates/${POP}/scen5_mi/v${VARIANT}"
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK + 1 ))
LAST=$(( SLURM_ARRAY_TASK_ID * BLOCK )); [ "$LAST" -gt "$N_REPS" ] && LAST="$N_REPS"

module purge; module add stata/19
export STATATMP="${SLURM_SUBMIT_DIR}/.statatmp"
WORK="${STATATMP}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORK" "$OUT_DIR"; cd "$WORK"

for SIM in $(seq "$FIRST" "$LAST"); do
    OUT="${OUT_DIR}/scen5_mi_coef_${SIM}.csv"
    [ -s "$OUT" ] && continue
    START=$(date +%s)
    timeout 120m stata-se -b do "${SLURM_SUBMIT_DIR}/scripts/stata/scen5_mi.do" "$SIM" "$VARIANT" "$POP"
    [ -s "$OUT" ] && echo "rep $SIM done in $(( $(date +%s) - START )) s" || echo "rep $SIM FAILED after $(( $(date +%s) - START )) s" >&2
done

cd "$SLURM_SUBMIT_DIR"; rm -rf "$WORK"
