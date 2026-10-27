#!/bin/bash
#SBATCH -p general
#SBATCH --account=rc_cai_pi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.err

# One runner for every scenario. Each array task runs a block of consecutive
# replicates of one do-file; the block is sized from the measured seconds per
# replicate so each job lasts about an hour (Research Computing guidance; one
# replicate per job got this account's CPU limit throttled).
#
# Usage:
#   sbatch --array=1-<NTASKS> --job-name=<name> run_stata_chunked.slurm \
#          <do_name> <avg_sec> <out_template> [extra do-file args...]
#
#   <avg_sec>       measured seconds per replicate (from a smoke test); the block
#                   size is ceil(3600 / avg_sec), i.e. ~1 hour of work per job
#   <out_template>  output CSV path relative to the repo, with {SIM} for the replicate
#   extra args      passed to the do-file after the replicate id (variant, population)
#
# Pick the array size: NTASKS = ceil(N_REPS * avg_sec / 3600).
#
# scen9 example, dec2024, ~30 s/rep, 1000 reps -> block 120, NTASKS=9:
#   sbatch --array=1-9 --job-name=scen9_dec2024 run_stata_chunked.slurm \
#          scen9 30 'data/estimates/dec2024/scen9/scen9_coef_{SIM}.csv' dec2024
#
# scen15 example, oct2026 v1, ~400 s/rep, 1000 reps -> block 9, NTASKS=112:
#   sbatch --array=1-112 --job-name=scen15_v1_oct2026 run_stata_chunked.slurm \
#          scen15_mi 400 'data/estimates/oct2026/scen15_mi/v1/scen15_mi_coef_{SIM}.csv' 1 oct2026
#
# N_REPS (default 1000) and REP_TIMEOUT (default 120m) may be set in the environment.

set -u
DO_NAME="$1"; AVG_SEC="$2"; OUT_TEMPLATE="$3"; shift 3
EXTRA_ARGS=("$@")
N_REPS="${N_REPS:-1000}"
REP_TIMEOUT="${REP_TIMEOUT:-120m}"

REPO="/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
DO_FILE="${REPO}/scripts/stata/${DO_NAME}.do"

BLOCK=$(( (3600 + AVG_SEC - 1) / AVG_SEC ))
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK + 1 ))
LAST=$(( SLURM_ARRAY_TASK_ID * BLOCK ))
[ "${LAST}" -gt "${N_REPS}" ] && LAST="${N_REPS}"

module purge
module add stata/19

export STATATMP="${REPO}/.statatmp"
TASK_DIR="${STATATMP}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TASK_DIR}" "${REPO}/slurm_logs"
cd "${TASK_DIR}"

echo "Task ${SLURM_ARRAY_TASK_ID}: ${DO_NAME} replicates ${FIRST}-${LAST} (block ${BLOCK}), args: ${EXTRA_ARGS[*]:-none}"
N_DONE=0; N_SKIP=0; N_FAIL=0
for SIM in $(seq "${FIRST}" "${LAST}"); do
    OUT_CSV="${REPO}/${OUT_TEMPLATE//\{SIM\}/${SIM}}"
    if [ -s "${OUT_CSV}" ]; then
        N_SKIP=$(( N_SKIP + 1 )); continue
    fi
    START=$(date +%s)
    timeout "${REP_TIMEOUT}" stata-se -b do "${DO_FILE}" "${SIM}" "${EXTRA_ARGS[@]}"
    RC=$?
    cp "${DO_NAME}.log" "${REPO}/slurm_logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_rep${SIM}.log" 2>/dev/null
    if [ -s "${OUT_CSV}" ]; then
        N_DONE=$(( N_DONE + 1 ))
        echo "  replicate ${SIM}: done in $(( $(date +%s) - START )) s"
    else
        N_FAIL=$(( N_FAIL + 1 ))
        echo "  replicate ${SIM}: FAILED (exit ${RC}; 124 means it hit the ${REP_TIMEOUT} limit) after $(( $(date +%s) - START )) s" >&2
    fi
done

cd "${REPO}"
rm -rf "${TASK_DIR}"
echo "Task ${SLURM_ARRAY_TASK_ID} finished: ${N_DONE} done, ${N_SKIP} skipped (already existed), ${N_FAIL} failed"
