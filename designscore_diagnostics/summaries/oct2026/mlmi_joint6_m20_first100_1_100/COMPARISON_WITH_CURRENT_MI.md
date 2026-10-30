# Multilevel MI versus current MI

With multilevel MI (M=20), design-score subject-variance bias changes from +2.30 to +0.44, and residual-variance bias from -1.09 to -1.06. Original scenario 3 subject-variance bias changes from +7.60 to +5.48.

Matched Oct2026 00pct samples 1-100. Both GFR and BMI were masked using miss_ind_mar. Each pair of analysis methods used identical imputations within its imputation method.

The multilevel run uses M=20; existing wide chained MI uses M=5 and was not rerun. The prespecified M=5 sensitivity uses multilevel imputation indices [1, 6, 11, 16, 2]. This subset changes pooling only; it does not refit models.

The new imputer jointly models the six GFR/BMI measurements with shared PSU and household effects and an unrestricted subject-trajectory covariance. It also adds continuous baseline age, full stratum indicators, and the specified age/stratum interactions. This comparison evaluates that complete imputation specification; it does not isolate the contribution of random effects alone. The analysis estimators and population targets are unchanged.

| Analysis | Component | Target | Current MI mean (bias) | Multilevel MI mean (bias) | Coverage: current / multilevel |
|---|---|---:|---:|---:|---:|
| original | var_hhid | 22.74 | 18.67 (-4.07) | 20.38 (-2.36) | 100% / 99% |
| original | var_subid | 175.17 | 182.77 (+7.60) | 180.65 (+5.48) | 86% / 92% |
| original | s2 | 95.09 | 94.00 (-1.09) | 94.03 (-1.06) | 88% / 88% |
| designscore | var_hhid | 22.74 | 20.09 (-2.65) | 21.68 (-1.06) | 100% / 98% |
| designscore | var_subid | 175.17 | 177.47 (+2.30) | 175.62 (+0.44) | 95% / 95% |
| designscore | s2 | 95.09 | 94.00 (-1.09) | 94.03 (-1.06) | 88% / 87% |

For the prespecified five-imputation subset, design-score biases are var_hhid: -1.05 (coverage 97%), var_subid: +0.51 (coverage 94%), s2: -1.05 (coverage 87%). The CSV comparison includes both methods and all parameters.

## Residual covariance check

These unweighted residual moments use the same saved complete-data coefficients and source covariates. They help check covariance preservation but are not fitted variance components. Each sample receives equal weight.

| Moment | Complete data | Current MI | Multilevel MI |
|---|---:|---:|---:|
| Household cross-person product | 20.85 | 18.67 | 20.27 |
| Same-subject cross-visit covariance | 203.38 | 203.04 | 202.69 |
| Within-person mean square | 95.23 | 93.83 | 93.90 |

## Checks and interpretation

All requested samples have complete paired results and passed the recorded input/chain checks. Largest R-hat: 1.0100; smallest bulk ESS: 443.8; smallest tail ESS: 581.6.

Chain lengths and source runs are recorded per replicate in `chain_diagnostics_summary.csv`. Any longer run retains the same model, priors, imputation count, and convergence thresholds; failed candidate imputations are excluded.

Coverage uses the same Rubin formulas for both analyses, with log-scale variance pooling and finite-imputation t intervals. The five-imputation sensitivity and paired Monte Carlo standard errors are in the accompanying CSV files. With 100 replicates, coverage differences retain appreciable Monte Carlo uncertainty.

Covariance-preservation diagnostics are unweighted descriptive calculations using fixed, previously saved complete-data coefficients. They are not additional fitted variance components and were not used to tune imputations to hidden simulation outcomes. Passing chain diagnostics does not establish that the imputation model is correctly specified.

## Files

- `comparison_with_current_mi.csv`: all parameters, both methods, M=20 and M=5 comparisons.
- `paired_changes_from_current_mi.csv`: paired changes and their Monte Carlo standard errors.
- `pooled_multilevel_m5_sensitivity.csv`: replicate-level five-imputation sensitivity.
- `chain_diagnostics_summary.csv`: sampler checks for all requested replicates.
- `covariance_preservation.csv`: complete versus imputed residual moments.
