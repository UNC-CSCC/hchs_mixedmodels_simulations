# Design score estimation for the Oct2026 three level mixed model

**Method report — 29 October 2026**

The design-score variant estimates household and subject variance components by constructing sampling-weighted versions of the Gaussian likelihood scores that would be available from the complete finite population. It retains the original scenario 3 first-stage estimates of the fixed effects and residual variance. Its second stage uses individual and pair inclusion probabilities, together with the full eligible household size, to estimate the two remaining variance components.

For the 1,000 saved Oct2026 complete-data samples, the estimated subject variance averaged **175.153937**, compared with the census target **175.174303**. The empirical bias was **−0.020366**, or **−0.0116%**, and coverage of the nominal 95% interval was **93.8%**. All 1,000 fits produced usable estimates and intervals. These results support the proposed correction to the variance estimating equations in this setting; the interval calculation remains an approximation to the sampling variance.

## Purpose and comparison with the original scenario

The original `scen3.do` uses two weighted mixed-model fits. The first uses globally scaled subject weights and supplies the fixed effects and residual variance. The second fixes those quantities and estimates the household and subject variances using subject weights scaled within each sampled household, `sub_s1`.

The archived results that motivated the change were:

| Population and method | Replicates | Subject variance target | Mean estimate | Empirical bias | Relative bias | Coverage |
|---|---:|---:|---:|---:|---:|---:|
| Dec2024 original scenario 3 | 1,000 | 117.310795 | 117.175998 | −0.134797 | −0.115% | 94.6% |
| Oct2026 original scenario 3 | 1,000 | 175.174303 | 181.174332 | +6.000029 | +3.425% | 83.0% |
| Oct2026 design score | 1,000 | 175.174303 | 175.153937 | −0.020366 | −0.0116% | 93.8% |

Sources: the archived [Dec2024 summary](../data/summaries/dec2024/coverage_scen3.csv), archived [Oct2026 summary](../data/summaries/oct2026/coverage_scen3.csv), and [design-score summary](summaries/oct2026/coverage_scen3_designscore_exact.csv). The comparison uses the saved summaries; it is not a replicate-level paired analysis with verified common input-file provenance.

Within a sampled household, the original subject weight has the form

```math
w^{(s1)}_{hi}
=
w_{hi}\,
\frac{\sum_{j\in s_h}w_{hj}}
     {\sum_{j\in s_h}w_{hj}^{2}},
\qquad
\sum_{i\in s_h} w^{(s1)}_{hi}
=
\frac{\left(\sum_{i\in s_h}w_{hi}\right)^2}
     {\sum_{i\in s_h}w_{hi}^{2}}.
```

The weights sum to an effective sample size, rather than an estimate of the number of eligible people in the household. In a household with one sampled person, the scaled subject weight is exactly one regardless of that person's original selection probability. If all sampled subject weights are equal, they also all become one.

This normalization changes how residual squares, cross-products, and cluster size enter the weighted likelihood. There is no general identity making the resulting variance scores unbiased estimates of the census scores. The design-score construction supplies that identity directly. The singleton example explains a limitation of the original weighting mechanism; it does not establish that singletons alone caused the entire Oct2026 bias.

## Model and estimation target

Let household `h` contain `M_h` eligible subjects, indexed by `i`. Each subject has exactly three complete visits, indexed by `t`. The working model is

```math
Y_{hit}=x_{hit}^{\mathsf T}\beta+u_h+v_{hi}+\varepsilon_{hit},
\qquad t=1,2,3,
```

with mutually independent Gaussian random terms

```math
u_h\sim N(0,a),\qquad
v_{hi}\sim N(0,b),\qquad
\varepsilon_{hit}\sim N(0,e).
```

The parameter names used in the code are:

| Mathematical quantity | Code or output name | Interpretation |
|---|---|---|
| `β` | Intercept and six slopes | Fixed effects |
| `a` | `var_hhid` | Household random-intercept variance |
| `b` | `var_subid` | Subject random-intercept variance |
| `e` | `s2` | Visit-level residual variance |
| `M_h` | `m_h` | Full eligible household size before subject sampling |

The fixed-effect design contains an intercept, `age_strat_new`, `x12`, `x17`, `x18`, `x6`, and `y_bmi`; the response is `y_gfr`.

The simulation evaluates estimates against the unweighted maximum-likelihood fit of this model to the saved finite population. Thus “truth” in the summary means the **census fitted parameter**, which need not equal the parameter originally used to generate the population. The Oct2026 targets are stored in [the population coefficient file](data/population_coefficients/population_coef_oct2026_hhsub.csv).

The design expectation derived below concerns recovery of this chosen census score. It is an algebraic property of sampling-weighted finite-population totals and does not require the Gaussian working model to be the exact data-generating distribution. Interpreting the fitted components as underlying model variances still depends on the working model.

## The retained first stage

The first stage fits the same complete-data model as the original scenario:

```stata
meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
    gl_hhid: , pweight(w_bghh) || ///
    gl_subid: , pweight(w_sub_glob) vce(cluster gl_hhid)
```

The visit weight is one. The household weight is the original nominal `w_bghh`. Raw subject weights are multiplied by a common normalizer:

```math
w_{hi}^{\mathrm{glob}}=c\,w_{hi},
\qquad
c=\frac{n_{\mathrm{obs}}}
        {\sum_{\text{observed visits}}w_{hi}}
 =\frac{n_{\mathrm{sub}}}
        {\sum_{\text{sampled subjects}}w_{hi}}.
```

The two expressions for `c` agree because every sampled subject contributes three visits. Stata supports weights at the separate levels of a multilevel model; its weighted-likelihood interface is described in the [official meglm manual](https://www.stata.com/manuals/memeglm.pdf).

The first stage supplies `β̂` and `ê`, which are held fixed during the second-stage point estimation. Its household and subject variances, denoted `â₁` and `b̂₁`, provide starting values and remain nuisance parameters in the standard-error calculation. The final household and subject variance estimates, `â₂` and `b̂₂`, come from the design-score stage.

The final standard errors for the fixed effects and residual variance are recomputed from the stacked equations described below. They are not simply copied from the first Stata fit.

## Sampling information required by the second stage

### Household selection and full household size

Let `π_h` be the probability that household `h` is selected through the block-group and household sampling stages, before subject selection. The outer household weight is

```math
W_h=\frac{1}{\pi_h}.
```

The household frame supplies `M_h` and `π_h`. It incorporates the actual rounded numbers of selected block groups within strata and selected households within block group and surname group. It includes households that may contribute no person to the final sample.

The second stage therefore uses exact stage 1–2 household probabilities from the frame. The first stage retains the original nominal `w_bghh` weights. This distinction is part of the implemented method.

The frame uses identifiers and design variables from the full population. Population outcome values are not used to fit a sample's design-score model.

### Subject and pair selection

Conditional on the selected stage-2 households, let `N_g` be the total number of eligible people in age group `g` and `n_g` the number selected without replacement. The two age groups are younger than 45 and age 45 or older.

For a person in group `g` and two distinct people in the same group,

```math
p_{hi}=\frac{n_g}{N_g},
\qquad
p_{hij}=\frac{n_g(n_g-1)}{N_g(N_g-1)}.
```

For a pair in different age groups, the separate age-group draws give

```math
p_{hij}
=
\frac{n_{\mathrm{younger}}}{N_{\mathrm{younger}}}
\frac{n_{\mathrm{older}}}{N_{\mathrm{older}}}.
```

These probabilities are conditional on the selected household set. The eligible counts cover **all** selected stage-2 households, including those absent from the final person-level file. Computing the denominators only from households that appear in that file would omit eligible people and give incorrect probabilities.

For the local run, [the recovery script](code/preparation/recover_designscore_counts.R) replayed the sampling draws from saved seeds and verified the sampled subject IDs and file checksums. It recovered selection metadata without replacing the sample files or generating new outcomes. The combined records are in `data/design_frames/oct2026_subject_counts.csv`.

All reported results use the explicit `exact` mode. The optional `bernoulli` mode instead uses nominal person probabilities 0.3575 and 0.55 and products for pairs. That approximation was not used for the 1,000-sample results.

## Derivation of the census variance scores

### Reduce the balanced visit records to subject means

For fixed `β`, define the average residual of subject `i`:

```math
r_{hi}
=
\frac{1}{3}\sum_{t=1}^{3}
\left(Y_{hit}-x_{hit}^{\mathsf T}\beta\right).
```

These residuals subtract the fixed-effect prediction only. They do not subtract predicted household or subject random effects.

Define

```math
d=b+\frac{e}{3},
\qquad
D_h=d+aM_h,
\qquad
r_h=(r_{h1},\ldots,r_{hM_h})^{\mathsf T}.
```

Under the working model, the covariance matrix of the subject means is

```math
V_h=dI_{M_h}+a\,\mathbf{1}_{M_h}\mathbf{1}_{M_h}^{\mathsf T}.
```

The within-subject visit contrasts depend on `e` and `β` but not on `a` or `b`. With `β` and `e` fixed, the subject means therefore contain the part of the Gaussian likelihood needed for the two variance scores.

### Express the likelihood through two household totals

Let

```math
B_h=\sum_{i=1}^{M_h}r_{hi}^{2},
\qquad
Q_h=\left(\sum_{i=1}^{M_h}r_{hi}\right)^2
   =B_h+2\sum_{i<j}r_{hi}r_{hj}.
```

The covariance identities are

```math
|V_h|=d^{M_h-1}D_h,
\qquad
V_h^{-1}
=
\frac{1}{d}I_{M_h}
-\frac{a}{dD_h}\mathbf{1}_{M_h}\mathbf{1}_{M_h}^{\mathsf T}.
```

Up to terms independent of `a` and `b`, the census household log likelihood is

```math
\ell_h(a,b\mid\beta,e)
=
-\frac12\left[
M_h\log d
+\log\left(\frac{D_h}{d}\right)
+\frac{B_h}{d}
-\frac{aQ_h}{dD_h}
\right].
```

Differentiation gives

```math
U_{a,h}
=
\frac{\partial\ell_h}{\partial a}
=
\frac12\left[
\frac{Q_h}{D_h^2}
-\frac{M_h}{D_h}
\right],
```

and

```math
U_{b,h}
=
\frac{\partial\ell_h}{\partial b}
=
\frac12\left[
\frac{B_h}{d^2}
-\frac{a(2d+aM_h)Q_h}{d^2D_h^2}
-\frac{M_h}{d}
+\frac{aM_h}{dD_h}
\right].
```

For any fixed candidate parameter values, these scores are linear in the household totals `B_h`, `Q_h`, and the count appearing in the trace terms. That structure permits direct inverse-probability estimation.

## Constructing the design weighted scores

For the sampled subjects `s_h` in household `h`, construct

```math
\widehat N_h
=
\sum_{i\in s_h}\frac{1}{p_{hi}},
\qquad
\widehat B_h
=
\sum_{i\in s_h}\frac{r_{hi}^2}{p_{hi}},
```

and

```math
\widehat Q_h
=
\widehat B_h
+
2\sum_{\substack{i<j\\i,j\in s_h}}
\frac{r_{hi}r_{hj}}{p_{hij}}.
```

The individual squared terms require first-order inclusion probabilities; the distinct-person cross-products require joint inclusion probabilities. This follows the general design-weighted estimating-equation principle discussed by [Rao, Verret, and Hidiroglou, section 3](https://www150.statcan.gc.ca/n1/pub/12-001-x/2013002/article/11887/section3-eng.htm). The particular equations here are derived from the full Gaussian census likelihood above; they are not a direct implementation of that paper's composite likelihood.

The pair-weighted quantity `Q̂_h` is not the square of an estimated total. For example, under an illustrative independent subject draw with probability one-half, a residual square receives weight two, while a distinct-person product receives weight four. Squaring an inverse-probability estimate of the residual total would incorrectly give the squared terms weight four as well. In a particular sample, `Q̂_h` can also be negative even though its census target `Q_h` is nonnegative.

The estimating equations used by the solver are

```math
\widehat U_a
=
\frac12\sum_{h\in s_{12}}W_h
\left[
\frac{\widehat Q_h}{D_h^2}
-\frac{\widehat N_h}{D_h}
\right]
=0,
```

and

```math
\widehat U_b
=
\frac12\sum_{h\in s_{12}}W_h
\left[
\frac{\widehat B_h}{d^2}
-\frac{a(2d+aM_h)\widehat Q_h}{d^2D_h^2}
-\frac{\widehat N_h}{d}
+\frac{a\widehat N_h}{dD_h}
\right]
=0.
```

Here `s₁₂` denotes the households selected at the first two sampling stages. Households with no selected person have `N̂_h = B̂_h = Q̂_h = 0`, so they need not appear as rows in the model-fitting dataset.

**Full household size remains in the covariance denominators.** The solver uses `D_h = d + aM_h` from the frame, not `d + aN̂_h`, the observed sample size, or an effective sample size. Replacing the denominator by an estimated count would introduce a nonlinear ratio and break the expectation identity below.

### Why the count terms also use an estimated total

Although the full `M_h` is known from the frame, the fitting dataset contains only households with at least one sampled person. Attaching a fixed `M_h` trace contribution only to those nonempty households would condition that contribution on appearing in the final sample.

Using `N̂_h` makes the entire contribution zero for an empty sampled household while retaining the correct expectation. This allows the equations to be evaluated from observed subjects without explicitly inserting all empty stage-2 households.

### The design expectation identity

Conditional on the selected stage-2 household set, the inverse-probability identities give

```math
E_3(\widehat N_h\mid s_{12})=M_h,\qquad
E_3(\widehat B_h\mid s_{12})=B_h,\qquad
E_3(\widehat Q_h\mid s_{12})=Q_h.
```

For fixed `a`, `b`, `β`, and `e`, the covariance denominators are fixed and the scores are linear in these estimated totals. Applying expectation over subject selection and then household selection yields

```math
E_{\mathrm{design}}\left[
\begin{pmatrix}\widehat U_a\\\widehat U_b\end{pmatrix}
\right]
=
\sum_{h\in\mathcal H}
\begin{pmatrix}U_{a,h}\\U_{b,h}\end{pmatrix},
```

where `ℋ` is the full population of households. The result requires the stated inclusion probabilities, positive probabilities for the required individual and pair contributions, and correct frame sizes.

This is a **design-unbiased estimate of the census score at fixed parameter values**. It does not imply that the nonlinear solution has zero finite-sample bias. In the implemented procedure, `β̂` and `ê` are estimated from the same sample, so the fixed-parameter identity also does not make the plug-in estimating functions exactly unbiased. The simulation evaluates these remaining effects.

## Numerical solution in Mata

The two equations are gradients of the integrated objective

```math
\widehat L(a,b\mid\beta,e)
=
-\frac12\sum_hW_h\left[
\widehat N_h\log d
+\frac{\widehat N_h}{M_h}\log\left(\frac{D_h}{d}\right)
+\frac{\widehat B_h}{d}
-\frac{a\widehat Q_h}{dD_h}
\right].
```

At fixed parameter values, the same design-expectation argument recovers the census conditional log likelihood, up to the omitted constants. This objective is used to solve the estimating equations; it is not an ordinary likelihood for the observed sample.

The implementation in [variance_designscore.mata](code/stata/variance_designscore.mata) precomputes quadratic summaries of the outcomes and covariates. Residual squares and cross-products can then be evaluated at different fixed-effect values without reimporting the data. This is useful both for optimization and for numerical differentiation of the stacked equations.

The numerical procedure is:

1. Use the first-stage variances as starting values, with additional starts obtained by multiplying them by `(0.25, 2)` and `(4, 0.5)`.
2. Optimize the interior region using `log(a)` and `log(b)`.
3. Optimize each one-variance-zero boundary separately, and evaluate the corner `a = b = 0`.
4. Reject candidates that fail convergence, numerical guard, or score checks.
5. Select the accepted candidate with the highest objective, allowing boundary selection when objectives are effectively tied.

There are nine optimization attempts across the three starts and three regions, plus the direct corner check. Each optimization permits up to 300 iterations and uses analytic first derivatives with BFGS through Mata's [optimize interface](https://www.stata.com/manuals/m-5optimize.pdf).

The implementation checks variance-scale scores rather than relying only on gradients with respect to log variances. If `v_j` denotes either variance and `S = Σ_h W_h` over the observed households, its convergence diagnostic is

```math
K
=
\max_j\frac{v_j+e/3}{S}
\begin{cases}
|\widehat U_j|,&v_j>0,\\
\max(0,\widehat U_j),&v_j=0.
\end{cases}
```

Accepted solutions require `K ≤ 10⁻⁶`. At a zero boundary, the maximizing score must be nonpositive in the feasible direction. This check prevents a very small variance from creating a misleadingly small log-scale gradient.

The multiple starts and boundary checks improve numerical robustness, but do not prove a global maximum for every possible dataset.

## Standard errors and confidence intervals

### Propagating uncertainty from both stages

Holding `β̂` and `ê` fixed during point estimation does not make them known without error. The standard-error calculation therefore stacks the first- and second-stage estimating equations.

The parameter vector has 13 entries:

```math
\theta=
\left(
\beta^{\mathsf T},
\log a_1,\log b_1,\log e,\log c,
\log a_2,\log b_2
\right)^{\mathsf T}.
```

There are seven fixed effects, three first-stage covariance parameters, the global subject-weight normalizer, and two second-stage variances. The household estimating vector contains the corresponding first-stage scores, the normalizer equation

```math
g_{c,h}=n_h-c\sum_{i\in s_h}w_{hi},
```

and the two design scores. Here `n_h` is the number of sampled subjects in the household. The normalizer equation is not multiplied by the household weight, matching the global normalization used in the Stata fit.

Let `ψ_h(θ)` be this stacked household estimating vector. The Jacobian is

```math
A
=
-\frac{\partial}{\partial\theta^{\mathsf T}}
\sum_h\psi_h(\theta)
\bigg|_{\widehat\theta}.
```

The implementation obtains it by central differences of the analytic household scores. Equation rescaling is used for numerical stability.

### Block groups and design strata

Block groups are treated as primary sampling units for the sandwich calculation. Within design stratum `g`, let `m_g` be the number of observed block groups and `Ψ_gk` the sum of household score vectors in block group `k`. The meat of the sandwich is

```math
B
=
\sum_g\frac{m_g}{m_g-1}
\sum_{k=1}^{m_g}
\left(\Psi_{gk}-\overline\Psi_g\right)
\left(\Psi_{gk}-\overline\Psi_g\right)^{\mathsf T},
\qquad
\widehat{\operatorname{Var}}(\widehat\theta)
=
A^{-1}BA^{-\mathsf T}.
```

This is a stratified, with-replacement, ultimate-cluster approximation. It propagates first-stage parameter uncertainty and covariance with the second-stage scores. It does not include exact multistage finite-population corrections or fully represent the cross-block-group dependence induced by fixed-size subject selection.

Thus **`exact` describes the inclusion probabilities used for point estimation, not exact variance estimation or exact confidence-interval coverage**.

### Interval construction and failure handling

Fixed-effect intervals use the normal approximation. For a positive variance estimate `v̂`, the interval is constructed on the log scale:

```math
\left[
\exp\{\log\widehat v-z_{0.975}\operatorname{SE}(\log\widehat v)\},
\exp\{\log\widehat v+z_{0.975}\operatorname{SE}(\log\widehat v)\}
\right].
```

The exported variance-scale standard error is `v̂ × SE(log v̂)`. This interval construction is also used for the retained residual variance.

The code checks the first-stage score residuals, Jacobian conditioning, number of block groups per stratum, and finiteness of the covariance matrix. At a boundary solution, or when a sandwich check fails, it preserves point estimates but withholds all standard errors and intervals. Missing intervals are reported explicitly rather than replaced by zeros.

## Results from all 1000 Oct2026 samples

The run used `samplemiss_1.csv` through `samplemiss_1000.csv` in `sample_miss00pct`, with exact selection metadata. Samples 201–1,000 were added to the unchanged results for samples 1–200. The population truth was not refitted.

| Component | Census target | Mean estimate | Empirical bias | Empirical SE | Mean model SE | Coverage |
|---|---:|---:|---:|---:|---:|---:|
| Household variance | 22.738064 | 22.549958 | −0.188107 | 5.566258 | 5.189840 | 95.5% |
| Subject variance | 175.174303 | 175.153937 | −0.020366 | 6.261990 | 6.007524 | 93.8% |
| Residual variance | 95.089839 | 95.112469 | +0.022630 | 1.667528 | 1.716098 | 94.6% |

“Empirical SE” is the standard deviation of the estimates across replicates. “Mean model SE” is the arithmetic mean of the exported standard errors. Coverage is the proportion of intervals containing the census target. All three rows have 1,000 point estimates and 1,000 valid intervals.

For subject variance, the Monte Carlo SE of the empirical bias is **0.198022**. The observed bias is small relative to this uncertainty. Its mean model SE is **4.06% below** the empirical SE. Coverage is 1.2 percentage points below nominal, with a Monte Carlo SE of **0.763 percentage points**. A normal Monte Carlo interval around the observed coverage is approximately **92.3%–95.3%**, so these 1,000 replications alone do not establish a departure from 95% coverage at the usual two-sided 5% level.

The variance-score change addresses the prominent subject-variance bias in this run. Standard-error calibration remains a separate issue. The comparison with the original scenario also changes the standard-error calculation: the original scenario used household-cluster standard errors, whereas the new procedure uses a stacked block-group sandwich that propagates uncertainty across both stages. The coverage difference cannot be attributed solely to the point-estimation change.

Across all fixed effects, coverage ranges from 93.1% to 95.5%; their individual results are included in the [complete summary](summaries/oct2026/coverage_scen3_designscore_exact.csv).

### Numerical and data verification

All 1,000 solutions were interior, all 10,000 coefficient rows had usable intervals, and no requested replicate was missing. The largest scaled KKT residual was **4.063 × 10⁻⁷**, the largest first-stage score discrepancy was **7.555 × 10⁻⁶**, and the largest scaled Jacobian condition number was **13,892.60**. Each sample contained 220 observed block groups and 14,200–16,559 sampled subjects.

Independent R calculations agreed with the Python summary's means, empirical and model SEs, coverage rates, and denominators within 10⁻¹⁰. The verification records and per-replicate audit distinguish successful optimization from usable interval estimation. These checks establish numerical consistency of the saved run; they do not make the approximate sandwich an exact design variance estimator.

## Scope and further evaluation

The current implementation applies to Gaussian random-intercept models with exactly three complete visits per subject and the specified household/person sampling structure. Unequal visit counts, missing outcomes, random slopes, or other response distributions require new score derivations or extensions of the implementation.

The method requires full household sizes, valid household inclusion probabilities, and individual and within-household pair probabilities. Eligible subject counts cannot be reconstructed correctly from nonempty sampled households alone. A change in the sampling algorithm requires a corresponding review of the probability calculations.

The first-stage estimation strategy is retained. Its finite-sample behavior can still affect the second-stage variances, and the two-step estimator is not generally identical to fitting all census likelihood parameters jointly. The present empirical evidence is for the saved Oct2026 complete-data samples; it does not establish performance in other populations or missing-data settings.

A useful next evaluation would examine interval calibration using a variance estimator or replicate-weight procedure that more fully represents the sampling stages. Any such investigation should keep point-estimation bias, variance estimation, and interval construction distinguishable. The current evidence does not justify an empirical adjustment to the subject-variance point estimate.

## Code and supporting records

All paths below are relative to `designscore_diagnostics/`.

| File or folder | Role |
|---|---|
| [code/stata/scen3_designscore.do](code/stata/scen3_designscore.do) | Input validation, retained first stage, Mata invocation, and export |
| [code/stata/variance_designscore.mata](code/stata/variance_designscore.mata) | Household summaries, variance scores, objective, optimization, and sandwich |
| [code/preparation/prepare_designscore_frame.py](code/preparation/prepare_designscore_frame.py) | Full household sizes and household inclusion probabilities |
| [code/preparation/recover_designscore_counts.R](code/preparation/recover_designscore_counts.R) | Verified recovery of exact subject-selection counts |
| [code/evaluation/summarize_designscore.py](code/evaluation/summarize_designscore.py) | Bias, SE, coverage, Monte Carlo SE, and missing-output accounting |
| [data/input_locations.json](data/input_locations.json) | Locations of shared population and sample inputs |
| [summaries/oct2026/coverage_scen3_designscore_exact.csv](summaries/oct2026/coverage_scen3_designscore_exact.csv) | Results for all 1,000 samples |
| [summaries/oct2026/audit_scen3_designscore_exact.csv](summaries/oct2026/audit_scen3_designscore_exact.csv) | Numerical and inference status by replicate |
| [summaries/oct2026/provenance_scen3_designscore_exact.json](summaries/oct2026/provenance_scen3_designscore_exact.json) | Source, input-metadata, truth, and output fingerprints |
| [README.md](README.md) | Folder layout, interfaces, and execution instructions |

Within Mata, `ds_build()` constructs the household summaries; `ds_variance_scores()` and `ds_objective()` implement the displayed equations; `ds_solve()` selects an admissible solution; and `ds_stacked_scores()` and `ds_sandwich()` calculate the covariance estimate.

Historical source paths in the saved provenance refer to the source snapshots under `provenance/pre_reorganization/`. The relocation manifest records their current locations. The original scenario code and its saved summaries remain separate from this diagnostic implementation.

