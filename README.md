# hchs_mixedmodels_simulations

Simulation study of survey-weighted three-level linear mixed models (visits in
subjects in households) under the HCHS/SOL sampling design.

## Model

```
y_gfr ~ age_strat_new + x12 + x17 + x18 + x6 + y_bmi + (1 | household) + (1 | subject)
```

Targets: the fixed effects, the household and subject variances (`var_hhid`,
`var_subid`) and the residual variance (`s2`). The truth is the same model fitted,
unweighted, to the full population.

## Methods

**scen9: two-step pseudo-likelihood** (complete data)
1. Fit with globally scaled subject weights to estimate the fixed effects and `s2`.
2. Hold those fixed and refit with the household weight and the `s1`-scaled subject
   weight to estimate `var_hhid` and `var_subid`.

SEs are clustered on block groups.

**scen15 / scen16: multiple imputation + scen9**
Visit 2–3 outcomes are set missing at random and imputed by chained regression
(M = 5), with log design weights as imputation covariates. Each completed dataset
is fitted with scen9, and the results are pooled with Rubin's rules. scen15 uses
`s1`-scaled weights as covariates and scen16 uses `s2`-scaled weights:

| Variant | Imputation weight covariates |
|---|---|
| v1 | `log(w_bg) + log(hh_s)` |
| v2 | `log(w_bg) + log(hh_s) + log(sub_s)` |
| v3 | `log(w_bg) + log(hhsub_s)` |

**Weight scaling** (Pfeffermann et al., 1998), within each cluster *j*:

- `s1`, effective sample size: `w_ij × Σw_ij / Σw_ij²`. The scaled weights sum to the effective cluster size.
- `s2`, cluster sample size: `w_ij × n_j / Σw_ij`. The scaled weights sum to the actual cluster size.

## Runs

| Scenario | Populations |
|---|---|
| scen9 | dec2024, oct2026 |
| scen15 v1–v3 | oct2026 |
| scen16 v1–v3 | oct2026 |

## Layout

```
scripts/
  stata/        scen9.do, scen15_mi.do, scen16_mi.do, population_truth_hhsub.do
  slurm/        scen9.slurm, scen15_mi.slurm, scen16_mi.slurm, submit_coverage.slurm
  evaluation/   coverage_assessment.R
data/
  population_coefficients/   truth, one CSV per population (tracked)
  estimates/<pop>/<scen>/    one CSV per replicate (not tracked)
  summaries/<pop>/           coverage tables (tracked)
```

## Running

Run everything from the repo root.

**1. Population truth**, one per population:

```
module load stata/19
srun --mem=4G --time=00:30:00 stata-se -b do scripts/stata/population_truth_hhsub.do dec2024
```

**2. Time one replicate**, for each new scenario:

```
srun --mem=2G --time=00:30:00 stata-se -b do scripts/stata/scen9.do 1 dec2024
```

**3. Submit.** Set `BLOCK` in the scenario's `.slurm` file so each job takes about
an hour, then submit with `--array=1-ceil(1000/BLOCK)`. Replicates that already
have output are skipped.

| Run | Command |
|---|---|
| scen9 | `sbatch --array=1-9 scripts/slurm/scen9.slurm dec2024` |
| scen15 | `sbatch --array=1-112 scripts/slurm/scen15_mi.slurm 1 oct2026` |
| scen16 | `sbatch --array=1-112 scripts/slurm/scen16_mi.slurm 1 oct2026` |

For the MI scenarios, the first argument is the variant (1–3).

**4. Coverage**, once the estimates exist:

```
sbatch scripts/slurm/submit_coverage.slurm
```

This reports bias, empirical vs. model SE, model-based coverage, and coverage
from the empirical SE.
