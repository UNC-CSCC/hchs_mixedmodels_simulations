# Oct2026 paired MI comparison: replicates 1–100

All **100 replicates** completed with five shared imputations each.
Both methods used the same 500 completed datasets. The source folder
is `sample_miss00pct`; the MAR flag still masks visits (mean visit-2 missingness
30.0%, visit-3 missingness 45.4%).

| Parameter | Population target | Original MI bias | Original coverage | Design-score MI bias | Design-score coverage |
|---|---:|---:|---:|---:|---:|
| var_hhid | 22.738 | -4.072 | 100% | -2.649 | 100% |
| var_subid | 175.174 | +7.596 | 86% | +2.298 | 95% |
| s2 | 95.090 | -1.089 | 88% | -1.089 | 88% |

Coverage is for 95% Rubin intervals, using log-scale pooling for positive
variance components. Full fixed-effect and variance summaries include empirical
SE, mean model SE, and Monte Carlo SEs. The paired comparison isolates differences
between the analysis methods given these shared imputations.

The imputation model is the existing scenario 4 s1/v2 wide chained-regression
model. It does not explicitly include household or PSU random effects. Therefore
this run does not establish that the imputation model preserves all multilevel
covariance. With 100 replicates, coverage has substantial Monte Carlo uncertainty.

Independent R calculations reproduce the pooled estimates and summary statistics.
Checks of all saved MI datasets confirm that observed values are preserved,
MAR-flagged values are masked before imputation, and all completed imputations
contain the required visits. See `validation_mi.json` for verification records.

- [Original scenario 3 MI summary](coverage_scen3_original_mi.csv)
- [Design-score MI summary](coverage_scen3_designscore_mi.csv)
- [Paired differences](paired_comparison.csv)
- [Pooled estimates](pooled_estimates.csv)
- [Replicate audit](audit_mi.csv)
- [Methods and code guide](../../../MI_COMPARISON_README.md)
