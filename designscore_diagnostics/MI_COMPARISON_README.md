# Paired MAR multiple-imputation comparison

This extension compares original scenario 3 with design-score estimation on the
**same five completed imputations per sample**. The initial run uses Oct2026
replicates 1–100 from `sample_miss00pct`. Here, `00pct` means no baseline
covariate missingness; visit missingness is still introduced by `miss_ind_mar`.
Existing source scripts, the complete-data diagnostics, and the method report
are unchanged by this extension.

## Imputation and analysis

1. Set both GFR and BMI to missing wherever `miss_ind_mar == 1`. The hidden
   simulated outcomes are excluded from imputation and estimation.
2. Reshape to one row per person. Use the existing scenario 4 **s1, variant 2**
   chained linear-regression imputation model: visit 2/3 GFR and BMI are imputed
   jointly through chained equations. Predictors are the other imputed visit
   variables, baseline GFR/BMI, age group, `x12`, `x14`, `x17`, `x18`, visit
   times, three stratum indicators, and `log(w_bg)`, `log(hh_s1)`,
   `log(sub_s1)`. Use five imputations, the default 10 burn-in iterations,
   and seed `100000 + replicate`.
3. Save the completed MI set once. For each imputation, restore all three
   visits and run the shared first-stage fit for fixed effects and residual
   variance. Original scenario 3 then uses its constrained `meglm` fit with
   `sub_s1`; design-score uses the existing Mata solver and exact household,
   person, and pair inclusion probabilities. No response-probability weighting
   is added: MI supplies the missing visits, and the sampling weights remain
   unchanged.
4. Check that analysis-data signatures match before and after both fits.
   Retain each method's own standard errors: household-cluster errors for the
   original method; the stacked, stratified PSU sandwich for design-score.

The wide imputation model models associations among visits within a person.
It **does not explicitly model household or PSU random effects**. This is a
controlled comparison using the repository's existing imputation model, not
a claim that this imputation model fully preserves all multilevel covariance.
Remaining variance-component bias can reflect imputation as well as estimation.

## Pooling and summaries

For each parameter, pool all five estimates with Rubin's rules:
the total variance is the mean within-imputation variance plus
`(1 + 1/M)` times the between-imputation variance. Fixed effects use their
original scale. Positive variance components use `log(estimate)` and the
delta-method variance `(SE/estimate)^2`; transform the pooled estimate and
interval back afterward. Thus the reported variance point estimate is the
geometric mean across imputations. The arithmetic mean is also retained.

Both methods use the same finite-imputation t critical value, with Rubin
degrees of freedom `(M - 1) * (1 + within / ((1 + 1/M) * between))^2`.
The complete-data inference is treated as asymptotic; there is no additional
small-sample degrees-of-freedom correction. This pooling deliberately differs
from the older MI scripts' mixture of arithmetic point estimates and log-scale
normal intervals.

The coverage files retain the diagnostic summary fields: mean, empirical
bias, empirical/model SE, and coverage against the stored Oct2026 population
coefficients. They also report Monte Carlo SEs, counts of usable intervals,
and the requested denominator. The paired file gives design-score minus
original differences and their Monte Carlo SEs. Missing or failed imputations
are never silently dropped or assigned zero SE. With only 100 replicates,
coverage estimates are exploratory; five imputations also leave appreciable
imputation Monte Carlo noise.

## Files and reproducibility

All paths below are relative to `designscore_diagnostics`:

| Location | Contents |
|---|---|
| `code/stata/scen3_designscore_mi.do` | Shared imputations and paired fits |
| `code/execution/run_designscore_mi.py` | Local parallel runner and resume checks |
| `code/execution/finalize_designscore_mi.py` | Combined summary, independent verification, and results note |
| `code/evaluation/summarize_designscore_mi.py` | Rubin pooling and summaries |
| `data/estimates/oct2026/mi_s1_v2_m5/` | Saved MI datasets, per-imputation estimates, pooled estimates, audits |
| `logs/oct2026/mi_s1_v2_m5/` | Stata model and validation logs |
| `summaries/oct2026/mi_s1_v2_m5_1_100/` | Final summaries, pooled results, verification, and provenance |
| `provenance/mi_s1_v2_m5_inputs.json` | Source sample checksums, seeds, and design/code checksums |

Run from this folder on Windows with Python 3.11+, R, and Stata 19 or later:

```powershell
python code/execution/run_designscore_mi.py --first 1 --last 100 --workers 16
python code/execution/finalize_designscore_mi.py --first 1 --last 100
```

The runner defaults to the local Windows Stata path. Pass `--stata` to change
it; set `R_SCRIPT` for a different Rscript executable. Completed paired fits
are reused. Saved MI sets are reused if fitting needs a retry. Run labels
separate experiments; do not reuse a label with changed inputs or code.

`verify_mi_inputs.do` checks saved imputations against the original sample
values, without fitting models. The independent R verifier recalculates all
pooled results, coverage summaries, and paired comparisons. Small numerical
unit checks cover Rubin pooling and invalid inputs.

On this machine, Stata's `rename.ado` is missing. `mi_rename_compat.do` supplies
a session-only single-variable rename wrapper around the installed `_rename`
primitive, and only when `rename` is unavailable. It does not change the
Stata installation. Large local datasets and logs remain ignored by Git,
following the existing diagnostic-folder policy. No commits or pushes are
performed by this workflow.
