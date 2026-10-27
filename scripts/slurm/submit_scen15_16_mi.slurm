#!/bin/bash
#SBATCH -p general
#SBATCH --account=rc_cai_pi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.err

# All six MI arms (scen15/scen16 x variants 1-3) on one population, in one array.
# Each task runs CHUNK consecutive replicates of one arm. With C = ceil(N_REPS/CHUNK) tasks per arm,
# task t is arm (t-1)/C (0-2 = scen15 v1-3, 3-5 = scen16 v1-3), chunk (t-1)%C.
# Usage: sbatch --array=1-<6*C> --job-name=scen15_16_mi_<pop> submit_scen15_16_mi.slurm <population> <chunk>
#   population: dec2024 or oct2026 (the requested arms use oct2026)
# Smoke test (replicate 1 of each arm, reused by the full run): MAX_REP=1 sbatch --array=1-6 ... oct2026 1000
# Optional environment: N_REPS (default 1000), MAX_REP (run only replicates up to this), REP_TIMEOUT (default 120m)
set -u
POPULATION="$1"; CHUNK="$2"
N_REPS="${N_REPS:-1000}"
REP_TIMEOUT="${REP_TIMEOUT:-120m}"
REPO="/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"

C=$(( (N_REPS + CHUNK - 1) / CHUNK ))
T=$(( SLURM_ARRAY_TASK_ID - 1 ))
ARM=$(( T / C ))
CH=$(( T % C ))
SCEN=$(( ARM < 3 ? 15 : 16 ))
VARIANT=$(( ARM % 3 + 1 ))
FIRST=$(( CH * CHUNK + 1 ))
LAST=$(( (CH + 1) * CHUNK ))
[ "${LAST}" -gt "${N_REPS}" ] && LAST="${N_REPS}"
[ "${LAST}" -gt "${MAX_REP:-${N_REPS}}" ] && LAST="${MAX_REP}"

module purge
module add stata/19

export STATATMP="${REPO}/.statatmp"
TASK_DIR="${STATATMP}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
OUT_DIR="${REPO}/data/estimates/${POPULATION}/scen${SCEN}_mi/v${VARIANT}"
mkdir -p "${TASK_DIR}" "${OUT_DIR}" "${REPO}/slurm_logs"
cd "${TASK_DIR}"

echo "Task ${SLURM_ARRAY_TASK_ID}: scen${SCEN} v${VARIANT} (${POPULATION}) replicates ${FIRST}-${LAST}"
N_DONE=0; N_SKIP=0; N_FAIL=0
for SIM in $(seq "${FIRST}" "${LAST}"); do
    OUT_CSV="${OUT_DIR}/scen${SCEN}_mi_coef_${SIM}.csv"
    if [ -s "${OUT_CSV}" ]; then
        N_SKIP=$(( N_SKIP + 1 )); continue
    fi
    START=$(date +%s)
    timeout "${REP_TIMEOUT}" stata-se -b do "${REPO}/scripts/stata/scen${SCEN}_mi.do" "${SIM}" "${VARIANT}" "${POPULATION}"
    RC=$?
    cp "scen${SCEN}_mi.log" "${REPO}/slurm_logs/scen${SCEN}_mi_v${VARIANT}_${POPULATION}_${SLURM_ARRAY_JOB_ID}_rep${SIM}.log" 2>/dev/null
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
echo "Task ${SLURM_ARRAY_TASK_ID} finished: ${N_DONE} done, ${N_SKIP} skipped, ${N_FAIL} failed"
