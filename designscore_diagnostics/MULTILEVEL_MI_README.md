# Multilevel MI extension

Completed: Oct2026 00pct samples 1–100, with MAR visit masking and 20 shared
imputations per sample for original scenario 3 and design-score estimation.
The [comparison report](summaries/oct2026/mlmi_joint6_m20_first100_1_100/COMPARISON_WITH_CURRENT_MI.md)
contains the results and lists the output files. Design-score subject-variance
bias fell from +2.30 to +0.44, and household-variance bias from −2.65 to −1.06.
Residual-variance bias remained near −1.06. The five-imputation sensitivity gave
similar results.

All 100 samples passed the fixed chain checks and completed both analyses.
Sixty-five used the primary chain length; 35 needed longer chains with the same
model and priors. Independent R recalculation verified the pooled results and
summaries.

All files and outputs for this extension stay inside `designscore_diagnostics`. The
existing MI code, saved imputations, estimates, and summaries remain unchanged.

## Model

One row represents a subject. The outcome vector is ordered as
`GFR1, BMI1, GFR2, BMI2, GFR3, BMI3`. The joint Gaussian model is

```text
Y_i = X_i B + V_psu[i] + U_household[i] + E_i
V_g ~ N6(0, G)    U_h ~ N6(0, H)    E_i ~ N6(0, R)
```

`G`, `H`, and `R` are unrestricted 6-by-6 covariance matrices. A PSU effect is
shared across its households; a household effect is shared across its subjects.
The six entries of a subject's remaining trajectory are correlated through `R`.
Conditional on predictors, covariance is `G` for different households in a PSU,
`G + H` for different subjects in a household, and `G + H + R` for one subject.

`R` contains stable subject differences and visit-level variation. It is **not**
the analysis-model `s2`; its diagonal is not `var_subid`. This flexible wide
representation avoids requiring equal visit variances or a common cross-visit
covariance. Original scenario 3 and design-score estimation still determine the
three analysis variance components. Extra PSU structure in imputation does not
change those analysis models or their stored population targets.

Each outcome/visit has its own regression coefficients. Predictors are baseline
age, age >=45, `x12`, `x14`, `x17`, `x18`, visit-2/3 time, full stratum indicators,
the existing three log-weight predictors, `x17 * age_group`, and the three
`x12 * collapsed_stratum` interactions used in the population risk specification.
Baseline GFR/BMI are observed members of the joint outcome vector, not duplicated
regression predictors. Constant/dependent X columns are removed deterministically.

Scope is Oct2026 `sample_miss00pct`: baseline covariates must be complete. Both
GFR and BMI are masked when `miss_ind_mar==1`, before standardization or any
sampler calculation. Visit 1 must be observed. The method does not yet impute
missing baseline covariates or binary outcomes. No observation is dropped for
having only one sampled household member or one observed visit.

## Posterior draws

The sampler combines Gibbs/data augmentation with Metropolis updates in Mata:

1. Draw the regression coefficient matrix from its matrix-normal conditional.
2. Draw each shared household vector, then each shared PSU vector, using all
   current member trajectories and the relevant covariance matrices.
3. Draw household, PSU, and person-trajectory covariance matrices.
4. Apply noncentered household/PSU covariance updates and a joint household/person
   covariance update that integrates out household effects, immediately redrawing
   those effects afterward. These additional transitions target the same model
   and priors; they address slow movement between variance components. The joint
   proposal is adapted during warmup only and then fixed.
5. Draw each person's missing outcome block from its multivariate normal
   conditional, given observed outcomes and the current parameters/effects.

Missing cells alone are written back; observed values are kept exactly. Parameters
and shared effects are drawn with uncertainty rather than replaced by fitted
means. Each retained dataset is one completed draw. Four independently seeded,
dispersedly initialized chains contribute equal numbers of completed datasets.

Outcomes are standardized using observed means/SDs, and X columns are standardized
before computation. For each covariance C the prior is
`C ~ inverse-Wishart(df=8, scale=s_C * I6)`. Default `s_H=0.05`, `s_G=0.01`, and
`s_R=0.5`; these are prior means on the standardized scale, not population truths.
The coefficient prior is `B | R ~ matrix-normal(0, I_p/ridge, R)`, with
`ridge=0.0001`. Priors are proper and exposed as run settings; small variance
components especially require sensitivity checks. The residual covariance
conditional includes both the `ridge * B'B` contribution and the additional p
degrees of freedom implied by this coefficient prior.

No population outcome values or fitted population variance targets enter the
sampler. Sampling/stratum predictors condition the imputation model; the sampler
does not treat survey weights as frequency counts or claim automatic design
consistency. The analysis weights and exact design-score inclusion probabilities
remain unchanged. The Gaussian joint model may still need refinement for Oct2026.

## Separate run stages

The new entry point is `code/execution/run_designscore_mlmi.py`. It uses the
local Stata executable and Python with NumPy/Pandas. Pooling uses the existing
R `qt()` helper. The default command prints settings, launches nothing, and
creates no run outputs:

```powershell
python code/execution/run_designscore_mlmi.py --stage plan
```

The stages can be invoked separately with explicit settings and a fresh label:

```powershell
# Generates imputations and checks them; does not run analysis models.
python code/execution/run_designscore_mlmi.py --stage impute --first 1 --last 1 --run-label mlmi_example --burn 2000 --between 1000

# Uses those validated saved imputations for both analysis methods.
python code/execution/run_designscore_mlmi.py --stage fit --first 1 --last 1 --run-label mlmi_example --burn 2000 --between 1000

# Reads saved estimates only, without fitting or imputing.
python code/execution/run_designscore_mlmi.py --stage summarize --first 1 --last 1 --run-label mlmi_example --burn 2000 --between 1000
```

`--stage all` performs imputation, validation, paired fits, and summaries.
The completed primary run used 20 imputations, four chains, 2,000 warmup iterations,
1,000 iterations between retained imputations, and a trace every five iterations.
The 35 longer runs used spacing 2,000; none needed the 4,000 fallback. All retained
the same model, priors, M, and thresholds. The original short CLI defaults remain in the frozen
source for reproducibility; they are not validated production settings.
`--workers` defaults to one. Subsequent stages must specify the same settings.
Use a new `--run-label mlmi_...` when changing the model,
priors, sampler settings, checks, or source code; existing runs are protected.

The convergence screen checks all regression coefficients, unique covariance
entries, and the augmented-data log posterior. It uses rank-normalized folded/
split R-hat <=1.01 and bulk/tail ESS >=400 by default. Burn-in is excluded. Failed
chains retain candidate data and traces but are not passed to either analysis.
Diagnostics cannot prove convergence. A failure calls for investigation and a
new, appropriately configured run; thresholds should not simply be relaxed to
obtain results.

## Identical imputations and checks

The wrapper imports completed draws with `mi import flong`, converts to standard
Stata wide MI storage, and saves a candidate file. Python checks trace completeness,
chain diagnostics, all observed cells, masking in m=0, and completeness of every
draw. A passing candidate is promoted to the single shared `mi_<sim>.dta`.

The unchanged `code/stata/scen3_designscore_mi.do` is then called **only after that
file exists and passes validation**. Its existing reuse branch reads it and skips
the old chained-equation imputer. It performs the common first-stage fit, then
original scenario 3 and design-score fits, with identical input signatures checked
for every imputation. The new driver also hashes the saved MI file before/after
fitting and seals the paired output hashes. The new summary module rejects
unvalidated or changed inputs and calls the existing Rubin pooling formulas,
while labeling the imputation method and seeds correctly.

Saved-data covariance diagnostics compare within-person residual mean squares,
three visit-pair contrasts, same-person cross-visit covariance, and household
cross-person residual products. They reuse existing complete-data coefficients
when available and do not refit anything. These descriptive checks are unweighted;
they are not new variance-component estimates. Hidden complete outcomes are
read only after imputation for these checks and are never used to select or
calibrate imputations. Agreement with simulation truth is not an acceptance gate.

## Files

| File | Purpose |
|---|---|
| `code/stata/multilevel_mi.mata` | Joint sampler and deterministic matrix helpers |
| `code/stata/multilevel_mi_collapsed.mata` | Joint household/person covariance transition |
| `code/stata/impute_multilevel_mi.do` | Masking, predictors, chains, saved MI datasets |
| `code/execution/run_designscore_mlmi.py` | Plan, generation, validation, paired fits, resume guards |
| `code/evaluation/check_multilevel_mi.py` | Chain, integrity, and covariance checks |
| `code/evaluation/summarize_designscore_mlmi.py` | Existing pooling formulas with multilevel metadata |
| `code/execution/complete_designscore_mlmi.py` | Longer-chain retries and final verification |
| `code/evaluation/collect_designscore_mlmi.py` | One validated paired result per requested sample |
| `code/evaluation/compare_multilevel_with_current_mi.py` | Matched comparison and prespecified M=5 sensitivity |
| `code/evaluation/verify_designscore_mlmi.R` | Independent recalculation of pooled inference and summaries |
| `code/execution/finalize_designscore_mlmi.py` | Input/output checks and final report generation |
| `code/stata/check_multilevel_mi_code.do` | Definition/algebra/fixed-storage checks; no models |
| `code/evaluation/test_multilevel_mi_code.py` | Diagnostic and orchestration tests; no models |

Per-run outputs use `data/estimates/oct2026/<mlmi_label>/`, with
`imputations/`, `sampler/`, `per_imputation/`, `audits/`, and `pooled/` subfolders.
Logs use `logs/oct2026/<mlmi_label>/`; summaries use
`summaries/oct2026/<mlmi_label>_<first>_<last>/`; source hashes/settings/seeds use
`provenance/<mlmi_label>_inputs.json`. Covariance traces are standardized; the
`sampler/scaling_<sim>.csv` rows give observed means/SDs in the outcome order above.

The completed first-100 comparison is saved under
`summaries/oct2026/mlmi_joint6_m20_first100_1_100/`.
Its source-run table records any longer chains; failed/interrupted candidates
remain separate. The M=5 sensitivity uses imputation indices 1, 6, 11, 16, and 2,
chosen before examining analysis results. It only repools saved fits. The existing
wide-MI results are reused unchanged. The analyses retain seven-point integration.

References: [Stata MI import](https://www.stata.com/manuals/mimiimportflong.pdf),
[longitudinal and higher-level MI](https://pmc.ncbi.nlm.nih.gov/articles/PMC11755704/),
and [MCMC diagnostics](https://mc-stan.org/docs/reference-manual/analysis.html).
