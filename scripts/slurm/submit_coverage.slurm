#!/bin/bash
#SBATCH --job-name=coverage_eval
#SBATCH --account=rc_cai_pi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/coverage_%j.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/coverage_%j.err

module purge
module load r/4.4.0 2>/dev/null || module load r 2>/dev/null || true

mkdir -p "${SLURM_SUBMIT_DIR}/data/summaries" "${SLURM_SUBMIT_DIR}/slurm_logs"

# One summary per scenario/population/variant that has estimates on disk.
for PRESET in 1_dec2024 2_dec2024 3_dec2024 3_oct2026 \
              4_v1_oct2026 4_v2_oct2026 4_v3_oct2026 \
              5_v1_oct2026 5_v2_oct2026 5_v3_oct2026; do
  Rscript "${SLURM_SUBMIT_DIR}/scripts/evaluation/coverage_assessment.R" "${PRESET}"
done
