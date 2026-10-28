#!/bin/bash
#SBATCH -p general
#SBATCH --account=rc_cai_pi
#SBATCH --job-name=population_truth
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%j.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%j.err

set -u
POP="$1"
MODEL="${2:-hhsub}"

module purge; module add stata/19
export STATATMP="${SLURM_SUBMIT_DIR}/.statatmp"
WORK="${STATATMP}/${SLURM_JOB_ID}"
mkdir -p "$WORK"; cd "$WORK"

stata-se -b do "${SLURM_SUBMIT_DIR}/scripts/stata/population_truth_${MODEL}.do" "$POP"
cp "population_truth_${MODEL}.log" "${SLURM_SUBMIT_DIR}/slurm_logs/population_truth_${MODEL}_${POP}_${SLURM_JOB_ID}.log"

cd "$SLURM_SUBMIT_DIR"; rm -rf "$WORK"
