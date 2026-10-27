#!/bin/bash
#SBATCH --job-name=coverage_eval
#SBATCH --account=rc_cai_pi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/coverage_%j.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/coverage_%j.err

REPO="/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"

module purge
module load r/4.4.0 2>/dev/null || module load r 2>/dev/null || true

mkdir -p "${REPO}/data/summaries" "${REPO}/slurm_logs"

# One summary per scenario/population/variant that has estimates on disk.
for PRESET in 9_dec2024 9_oct2026 \
              15_v1_oct2026 15_v2_oct2026 15_v3_oct2026 \
              16_v1_oct2026 16_v2_oct2026 16_v3_oct2026; do
  Rscript "${REPO}/scripts/evaluation/coverage_assessment.R" "${PRESET}"
done
