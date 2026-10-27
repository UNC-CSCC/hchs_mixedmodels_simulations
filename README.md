# hchs_mixedmodels_simulations

Self-contained simulation study of weighted three-level mixed models (visits in
subjects in households) for the HCHS design. Inputs are read from the shared
beibo sample/population data; all outputs live inside this repo.

## Scenarios

| Scenario | Method | Populations |
|----------|--------|-------------|
| `scen9`  | Two-step pseudo-likelihood (globally scaled weights for fixed effects & residual variance; effective-size weights for the household/subject variances) | `dec2024`, `oct2026` |
| `scen15` | Multiple imputation with **effective**-size scaled weights (`s1`) as imputation-model covariates, variants v1–v3 | `oct2026` |
| `scen16` | Multiple imputation with **cluster**-size scaled weights (`s2`) as imputation-model covariates, variants v1–v3 | `oct2026` |

MI variants: v1 = `log(w_bg) + log(hh_s)`, v2 = `+ log(sub_s)`, v3 = `log(w_bg) + log(hhsub_s)`.

## Input data (read-only, external)

- Dec2024 samples: `/work/users/b/e/beibo/v3data/Dec2024/sample_miss00pct/samplemiss_<sim>.csv`
- Oct2026 samples: `/work/users/b/e/beibo/v3data/Oct2026/sample_miss00pct/samplemiss_<sim>.csv`
- Populations (for truth): `population_3visits_Dec2024.csv`, `population_3visits_Oct2026.csv`

`scen9` clusters SEs on block groups for both populations. The MI arms cluster on
block groups for `dec2024` and on households for `oct2026`.

## Layout

```
scripts/stata/        scen9.do, scen15_mi.do, scen16_mi.do, population_truth_hhsub.do
scripts/slurm/        run_stata_chunked.slurm, submit_scen15_16_mi.slurm, submit_coverage.slurm
scripts/evaluation/   coverage_assessment.R
data/estimates/       per-replicate coefficient CSVs (gitignored)   data/estimates/<pop>/<scen>/...
data/summaries/       coverage summaries (gitignored)               data/summaries/<pop>/...
data/population_coefficients/  tracked truth CSVs (population_coef_<pop>_hhsub.csv)
slurm_logs/           job logs (gitignored)
```

## Running

All jobs bundle replicates into ~1-hour chunks (see the repo CLAUDE.md rules).
**Smoke-test one replicate of each new code path before any full array.**

Population truth (one short job per population):

```
srun --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:30:00 \
  stata-se -b do scripts/stata/population_truth_hhsub.do dec2024
```

scen9, e.g. dec2024 (1000 reps, 60 per task → 17 jobs):

```
sbatch --array=1-17 --job-name=scen9_dec2024 scripts/slurm/run_stata_chunked.slurm \
  scen9 60 'data/estimates/dec2024/scen9/scen9_coef_{SIM}.csv' dec2024
```

scen15/scen16 all six arms on oct2026 (chunk 9 → C=112 tasks/arm → array 1-672):

```
# smoke test first: replicate 1 of each arm
MAX_REP=1 sbatch --array=1-6 --job-name=scen15_16_smoke scripts/slurm/submit_scen15_16_mi.slurm oct2026 1000
# then the full run
sbatch --array=1-672 --job-name=scen15_16_mi_oct2026 scripts/slurm/submit_scen15_16_mi.slurm oct2026 9
```

Coverage summaries (after estimates exist):

```
sbatch scripts/slurm/submit_coverage.slurm
```

Preset names for `coverage_assessment.R`: `9_dec2024`, `9_oct2026`,
`15_v{1,2,3}_oct2026`, `16_v{1,2,3}_oct2026`.
