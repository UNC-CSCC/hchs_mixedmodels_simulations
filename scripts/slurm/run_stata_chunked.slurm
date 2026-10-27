#!/bin/bash
#SBATCH -p general
#SBATCH --account=rc_cai_pi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.err

# Bundled runner: each array task runs a CHUNK of consecutive replicates of one
# Stata scenario, so each job lasts about an hour or more (Research Computing
# guidance; one replicate per job got this account's CPU limit throttled).
#
# Usage:
#   sbatch --array=1-<ceil(R/CHUNK)> --job-name=<name> scripts/slurm/run_stata_chunked.slurm \
#          <do_name> <chunk> <out_template> [extra do-file args...]
#
#   <do_name>       do-file in scripts/stata without .do, e.g. scen15_mi
#   <chunk>         replicates per task; task t runs replicates (t-1)*chunk+1 .. t*chunk
#   <out_template>  output CSV path relative to the repo, with {SIM} for the replicate,
#                   used to skip replicates that already exist
#   extra args      passed to the do-file after the replicate id (e.g. the variant, population)
#
# scen9 example (dec2024, 1000 replicates, 60 per task -> 17 jobs):
#   sbatch --array=1-17 --job-name=scen9_dec2024 scripts/slurm/run_stata_chunked.slurm \
#          scen9 60 'data/estimates/dec2024/scen9/scen9_coef_{SIM}.csv' dec2024
#
# scen15 example (v1, oct2026, 1000 replicates, 9 per task -> 112 jobs):
#   sbatch --array=1-112 --job-name=scen15_v1_oct2026 scripts/slurm/run_stata_chunked.slurm \
#          scen15_mi 9 'data/estimates/oct2026/scen15_mi/v1/scen15_mi_coef_{SIM}.csv' 1 oct2026
#
# Optional environment variables:
#   N_REPS       total replicates (default 1000); replicates above it are not run
#   REP_TIMEOUT  per-replicate time limit (default 120m), so one non-converging
#                replicate cannot use up the whole job and block the rest of its chunk

set -u
DO_NAME="$1"; CHUNK="$2"; OUT_TEMPLATE="$3"; shift 3
EXTRA_ARGS=("$@")
N_REPS="${N_REPS:-1000}"
REP_TIMEOUT="${REP_TIMEOUT:-120m}"

REPO="/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
DO_FILE="${REPO}/scripts/stata/${DO_NAME}.do"
[ -f "${DO_FILE}" ] || { echo "No such do-file: ${DO_FILE}" >&2; exit 2; }

FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
LAST=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
[ "${LAST}" -gt "${N_REPS}" ] && LAST="${N_REPS}"

module purge
module add stata/19

export STATATMP="${REPO}/.statatmp"
TASK_DIR="${STATATMP}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TASK_DIR}" "${REPO}/slurm_logs"
cd "${TASK_DIR}"

echo "Task ${SLURM_ARRAY_TASK_ID}: ${DO_NAME} replicates ${FIRST}-${LAST}, args: ${EXTRA_ARGS[*]:-none}"
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
