# Design score method for the Oct2026 mixed model

**Method report — 29 October 2026**

The design-score method changes how scenario 3 estimates household and subject variances. It weights residual squares and cross-products to reproduce the estimating equations for the complete population. In 1,000 Oct2026 samples, it reduced subject-variance bias from **+6.000 to −0.020**; coverage improved from **83.0% to 93.8%**.

## Why the original scenario 3 is biased for Oct2026

The original method estimates fixed effects and residual variance first, then holds them fixed while estimating household and subject variances. Its second fit uses `sub_s1`: subject weights rescaled within each sampled household to an effective sample size. For a household with one sampled person, that weight becomes one regardless of the person's selection probability.

Variance estimation depends on both individual residual squares and products between different subjects. Effective-size weights do not generally give these terms the individual and pair inclusion weights needed to recover the population equations. They also change how household size enters the likelihood. This can shift the variance estimates away from the population target.

The populations differ in ways that can expose this problem. [Dec2024](https://github.com/UNC-CSCC/HCHS-Simulations/blob/e9e8128/code/v3/Dec2024/population/generate_population.R) uses an additive outcome generator with Gaussian random effects. [Oct2026](https://github.com/UNC-CSCC/HCHS-Simulations/blob/e9e8128/code/v3/Oct2026/population/README.md#how-renal-calibration-and-clustering-are-preserved-together) calibrates and reassigns whole three-visit trajectories, changing their associations with covariates and clusters. Its calibration model also differs from scenario 3, including block-group effects and visit-specific residual variances. Good performance in December therefore does not guarantee good performance in October; neither generator exactly matches the fitted scenario 3 model.

The saved checks support the weighting explanation: adding stratum indicators leaves the subject-variance bias near 3.4%, while unweighted fits have much smaller bias. Simply removing weights is not a general solution, however: the unweighted age effect has 17.0% relative bias and only 24.2% coverage. These comparisons implicate the weighting procedure without identifying the exact contribution of each population change. See the [stratum-adjusted](../data/summaries/oct2026/coverage_scen3_strata.csv), [household SRS](../data/summaries/oct2026/coverage_scen3_srs.csv), and [unweighted](../data/summaries/oct2026/coverage_scen3_unw.csv) summaries.

## Model and estimation target

### Overview

The method uses a sample to reconstruct the household and subject variance equations for the complete population. It first retains scenario 3's estimates of the fixed effects and visit-level residual variance. It then averages each person's three residuals and weights their squares and cross-products using individual and pair selection probabilities. The Mata solver uses those totals to estimate household and subject variances jointly, and a stacked sandwich calculation carries uncertainty from both stages into the confidence intervals.

The sections below explain the population target, required sampling information, score derivation, numerical solution, and uncertainty calculation.

### Working model and target

Let household `h` contain `M_h` eligible subjects, indexed by `i`. Each subject has exactly three complete visits, indexed by `t`. The working model is

```math
Y_{hit}=x_{hit}^{T}\beta+u_h+v_{hi}+\varepsilon_{hit},\qquad t=1,2,3.
```

with mutually independent Gaussian random terms

```math
u_h\sim N(0,a),\qquad v_{hi}\sim N(0,b),\qquad \varepsilon_{hit}\sim N(0,e).
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
w_{hi}^{G}=c w_{hi},\qquad
c=\frac{n_v}{\sum_{(h,i,t)\in s_v}w_{hi}}
=\frac{n_s}{\sum_{(h,i)\in s}w_{hi}}.
```

Here `n_v` and `n_s` are the numbers of observed visit records and sampled subjects; `s_v` and `s` denote those two sets. The two expressions for `c` agree because every sampled subject contributes three visits. Stata supports weights at the separate levels of a multilevel model; its weighted-likelihood interface is described in the [official meglm manual](https://www.stata.com/manuals/memeglm.pdf).

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
p_{hi}=\frac{n_g}{N_g},\qquad p_{hij}=\frac{n_g(n_g-1)}{N_g(N_g-1)}.
```

For a pair in different age groups, the separate age-group draws give

```math
p_{hij}=\frac{n_y}{N_y}\frac{n_o}{N_o}.
```

Here `y` and `o` denote the younger and older age groups. These probabilities are conditional on the selected household set. The eligible counts cover **all** selected stage-2 households, including those absent from the final person-level file. Computing the denominators only from households that appear in that file would omit eligible people and give incorrect probabilities.

For the local run, [the recovery script](code/preparation/recover_designscore_counts.R) replayed the sampling draws from saved seeds and verified the sampled subject IDs and file checksums. It recovered selection metadata without replacing the sample files or generating new outcomes. The combined records are in `data/design_frames/oct2026_subject_counts.csv`.

All reported results use the explicit `exact` mode. The optional `bernoulli` mode instead uses nominal person probabilities 0.3575 and 0.55 and products for pairs. That approximation was not used for the 1,000-sample results.

## Derivation of the census variance scores

### Reduce the balanced visit records to subject means

For fixed `β`, define the average residual of subject `i`:

```math
r_{hi}=\frac{1}{3}\sum_{t=1}^{3}(Y_{hit}-x_{hit}^{T}\beta).
```

These residuals subtract the fixed-effect prediction only. They do not subtract predicted household or subject random effects.

Define

```math
d=b+\frac{e}{3},\qquad D_h=d+aM_h,\qquad r_h=(r_{h1},\ldots,r_{hM_h})^{T}.
```

Under the working model, the covariance matrix of the subject means is

```math
V_h=dI_{M_h}+aJ_h.
```

Here `I` is the identity matrix and `J_h` is the `M_h` by `M_h` matrix of ones. The within-subject visit contrasts depend on `e` and `β` but not on `a` or `b`. With `β` and `e` fixed, the subject means therefore contain the part of the Gaussian likelihood needed for the two variance scores.

### Express the likelihood through two household totals

Let

```math
B_h=\sum_{i=1}^{M_h}r_{hi}^{2},\qquad
Q_h=\left(\sum_{i=1}^{M_h}r_{hi}\right)^2=B_h+2\sum_{i<j}r_{hi}r_{hj}.
```

The covariance identities are

```math
|V_h|=d^{M_h-1}D_h,\qquad V_h^{-1}=\frac{1}{d}I_{M_h}-\frac{a}{dD_h}J_h.
```

Up to terms independent of `a` and `b`, the census household log likelihood is

```math
\ell_h(a,b\mid\beta,e)=-\frac{1}{2}
\left[M_h\log(d)+\log\left(\frac{D_h}{d}\right)+\frac{B_h}{d}-\frac{aQ_h}{dD_h}\right].
```

Differentiation gives

```math
U_{a,h}=\frac{\partial\ell_h}{\partial a}
=\frac{1}{2}\left[\frac{Q_h}{D_h^2}-\frac{M_h}{D_h}\right].
```

and

```math
U_{b,h}=\frac{\partial\ell_h}{\partial b}
=\frac{1}{2}\left[\frac{B_h}{d^2}-\frac{a(2d+aM_h)Q_h}{d^2D_h^2}-\frac{M_h}{d}+\frac{aM_h}{dD_h}\right].
```

For any fixed candidate parameter values, these scores are linear in the household totals `B_h`, `Q_h`, and the count appearing in the trace terms. That structure permits direct inverse-probability estimation.

## Constructing the design weighted scores

For the sampled subjects `s_h` in household `h`, construct

```math
\hat{N}_h=\sum_{i\in s_h}\frac{1}{p_{hi}},\qquad
\hat{B}_h=\sum_{i\in s_h}\frac{r_{hi}^2}{p_{hi}}.
```

and

```math
\hat{Q}_h=\hat{B}_h+2\sum_{i<j}\frac{r_{hi}r_{hj}}{p_{hij}}.
```

The pair sum includes distinct sampled subjects `i` and `j` within household `h`. The individual squared terms require first-order inclusion probabilities; the distinct-person cross-products require joint inclusion probabilities. This follows the general design-weighted estimating-equation principle discussed by [Rao, Verret, and Hidiroglou, section 3](https://www150.statcan.gc.ca/n1/pub/12-001-x/2013002/article/11887/section3-eng.htm). The particular equations here are derived from the full Gaussian census likelihood above; they are not a direct implementation of that paper's composite likelihood.

The pair-weighted quantity `Q̂_h` is not the square of an estimated total. For example, under an illustrative independent subject draw with probability one-half, a residual square receives weight two, while a distinct-person product receives weight four. Squaring an inverse-probability estimate of the residual total would incorrectly give the squared terms weight four as well. In a particular sample, `Q̂_h` can also be negative even though its census target `Q_h` is nonnegative.

The estimating equations used by the solver are

```math
\hat{U}_a=\frac{1}{2}\sum_{h\in s_{12}}W_h
\left[\frac{\hat{Q}_h}{D_h^2}-\frac{\hat{N}_h}{D_h}\right]=0.
```

and

```math
\hat{U}_b=\frac{1}{2}\sum_{h\in s_{12}}W_h
\left[\frac{\hat{B}_h}{d^2}-\frac{a(2d+aM_h)\hat{Q}_h}{d^2D_h^2}-\frac{\hat{N}_h}{d}+\frac{a\hat{N}_h}{dD_h}\right]=0.
```

Here `s₁₂` denotes the households selected at the first two sampling stages. Households with no selected person have `N̂_h = B̂_h = Q̂_h = 0`, so they need not appear as rows in the model-fitting dataset.

**Full household size remains in the covariance denominators.** The solver uses `D_h = d + aM_h` from the frame, not `d + aN̂_h`, the observed sample size, or an effective sample size. Replacing the denominator by an estimated count would introduce a nonlinear ratio and break the expectation identity below.

### Why the count terms also use an estimated total

Although the full `M_h` is known from the frame, the fitting dataset contains only households with at least one sampled person. Attaching a fixed `M_h` trace contribution only to those nonempty households would condition that contribution on appearing in the final sample.

Using `N̂_h` makes the entire contribution zero for an empty sampled household while retaining the correct expectation. This allows the equations to be evaluated from observed subjects without explicitly inserting all empty stage-2 households.

### The design expectation identity

Conditional on the selected stage-2 household set, the inverse-probability identities give

```math
E_3(\hat{N}_h\mid s_{12})=M_h,\qquad
E_3(\hat{B}_h\mid s_{12})=B_h,\qquad
E_3(\hat{Q}_h\mid s_{12})=Q_h.
```

For fixed `a`, `b`, `β`, and `e`, the covariance denominators are fixed and the scores are linear in these estimated totals. Applying expectation over subject selection and then household selection yields

```math
E_p[\hat{U}_k]=\sum_{h\in H}U_{k,h},\qquad k=a,b.
```

Here `E_p` denotes expectation over the sampling design and `H` is the full population of households. The result requires the stated inclusion probabilities, positive probabilities for the required individual and pair contributions, and correct frame sizes.

This is a **design-unbiased estimate of the census score at fixed parameter values**. It does not imply that the nonlinear solution has zero finite-sample bias. In the implemented procedure, `β̂` and `ê` are estimated from the same sample, so the fixed-parameter identity also does not make the plug-in estimating functions exactly unbiased. The simulation evaluates these remaining effects.

## Numerical solution in Mata

The two equations are gradients of the integrated objective

```math
\hat{L}(a,b\mid\beta,e)=-\frac{1}{2}\sum_h W_h
\left[\hat{N}_h\log(d)+\frac{\hat{N}_h}{M_h}\log\left(\frac{D_h}{d}\right)+\frac{\hat{B}_h}{d}-\frac{a\hat{Q}_h}{dD_h}\right].
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
K=\max_{j=1,2}\frac{(v_j+e/3)k_j}{S}.
```

For a positive variance, `k_j` is the absolute value of its variance-scale score; at a zero boundary, `k_j = max(0, Û_j)`. Accepted solutions require `K ≤ 10⁻⁶`. At a zero boundary, the maximizing score must be nonpositive in the feasible direction. This check prevents a very small variance from creating a misleadingly small log-scale gradient.

The multiple starts and boundary checks improve numerical robustness, but do not prove a global maximum for every possible dataset.

## Standard errors and confidence intervals

### Propagating uncertainty from both stages

Holding `β̂` and `ê` fixed during point estimation does not make them known without error. The standard-error calculation therefore stacks the first- and second-stage estimating equations.

The parameter vector has 13 entries:

```math
\theta=(\beta^{T},\log(a_1),\log(b_1),\log(e),\log(c),\log(a_2),\log(b_2))^{T}.
```

There are seven fixed effects, three first-stage covariance parameters, the global subject-weight normalizer, and two second-stage variances. The household estimating vector contains the corresponding first-stage scores, the normalizer equation

```math
g_{c,h}=n_h-c\sum_{i\in s_h}w_{hi}.
```

and the two design scores. Here `n_h` is the number of sampled subjects in the household. The normalizer equation is not multiplied by the household weight, matching the global normalization used in the Stata fit.

Let `ψ_h(θ)` be this stacked household estimating vector. The Jacobian is

```math
A(\theta)=-\frac{\partial}{\partial\theta^{T}}\sum_h\psi_h(\theta).
```

The Jacobian is evaluated at the estimated parameter vector `θ̂`. The implementation obtains it by central differences of the analytic household scores. Equation rescaling is used for numerical stability.

### Block groups and design strata

Block groups are treated as primary sampling units for the sandwich calculation. Within design stratum `g`, let `m_g` be the number of observed block groups and `Ψ_gk` the sum of household score vectors in block group `k`. The meat of the sandwich is

```math
B=\sum_g\frac{m_g}{m_g-1}\sum_{k=1}^{m_g}
(\Psi_{gk}-\bar{\Psi}_g)(\Psi_{gk}-\bar{\Psi}_g)^{T}.
```

The estimated covariance matrix of `θ̂` is

```math
C=A^{-1}B(A^{-1})^{T}.
```

This is a stratified, with-replacement, ultimate-cluster approximation. It propagates first-stage parameter uncertainty and covariance with the second-stage scores. It does not include exact multistage finite-population corrections or fully represent the cross-block-group dependence induced by fixed-size subject selection.

Thus **`exact` describes the inclusion probabilities used for point estimation, not exact variance estimation or exact confidence-interval coverage**.

### Interval construction and failure handling

Fixed-effect intervals use the normal approximation. For a positive variance estimate `v̂`, the interval is constructed on the log scale:

```math
\left[\exp\left(\log(\hat{v})-z_{0.975}s_v\right),
\exp\left(\log(\hat{v})+z_{0.975}s_v\right)\right].
```

Here `s_v` is the standard error of the estimated log variance and `z₀.₉₇₅` is the standard normal 97.5th percentile (approximately 1.96). The exported variance-scale standard error is `v̂ × s_v`. This interval construction is also used for the retained residual variance.

The code checks the first-stage score residuals, Jacobian conditioning, number of block groups per stratum, and finiteness of the covariance matrix. At a boundary solution, or when a sandwich check fails, it preserves point estimates but withholds all standard errors and intervals. Missing intervals are reported explicitly rather than replaced by zeros.

## Stata implementation example

### Calling the existing wrapper

The following example uses Stata 19 and the already prepared household frame and subject-count file. Start in the parent `HCHS Simulations` directory:

```stata
version 19
cd "hchs_mixedmodels_simulations/designscore_diagnostics"

* Fit sample 1 if its saved coefficient file does not already exist.
do code/stata/run_designscore_local.do 1 1 exact
```

The batch wrapper skips an existing coefficient file. For a new fit it calls `scen3_designscore.do`, which reads `../../data/Oct2026/sample_miss00pct/samplemiss_1.csv`, uses the metadata in `data/design_frames/`, and writes coefficient and audit CSVs under `data/estimates/oct2026/exact/`, with a log under `logs/oct2026/exact/`.

### How Stata passes the fit to Mata

This excerpt illustrates the two fitting stages after input preparation. The wrapper has already merged `m_h` and `p_h`, created `gl_hhid`, `gl_subid`, `gl_psu`, `age_strat_new`, and the intercept column `ds_one`, and sorted records by household, subject, and visit. It has also created the probability row vector `(p_y, p_o, p_yy, p_oo, p_yo)` as the Stata matrix `probabilities`. Use the wrapper above for a complete run; validation and export steps are omitted here.

```stata
* Load the solver once in a fresh Stata session.
quietly do code/stata/variance_designscore.mata

* Step 1: the original globally scaled weighted fit.
gen double w_visit = 1
egen double sum_w_sub_total = total(w_sub)
scalar global_scale = _N / sum_w_sub_total[1]
gen double w_sub_glob = w_sub * global_scale

meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
    gl_hhid: , pweight(w_bghh) || ///
    gl_subid: , pweight(w_sub_glob) vce(cluster gl_hhid)
matrix rt1 = r(table)

* Pass seven fixed effects, a1, b1, e, and the weight normalizer.
matrix initial = J(1,11,.)
local j = 0
foreach v in _cons age_strat_new x12 x17 x18 x6 y_bmi {
    local ++j
    matrix initial[1,`j'] = _b[`v']
}
local hhcol = colnumb(rt1, "/:var(_cons[gl_hhid])")
local subcol = colnumb(rt1, "/:var(_cons[gl_hhid>gl_subid])")
local ecol = colnumb(rt1, "/:var(e.y_gfr)")
matrix initial[1,8] = rt1[1,`hhcol']
matrix initial[1,9] = rt1[1,`subcol']
matrix initial[1,10] = rt1[1,`ecol']
matrix initial[1,11] = global_scale

* Step 2: solve the design scores and compute stacked uncertainty.
mata: ds_run("y_gfr ds_one age_strat_new x12 x17 x18 x6 y_bmi", ///
    "gl_hhid gl_psu strat m_h w_bghh p_h w_sub age_strat_new", ///
    "initial", "probabilities", "estimates", "audit")

matrix colnames estimates = estimate std_error conf_low conf_high
matrix rownames estimates = intercept age_strat_new x12 x17 x18 ///
    x6 y_bmi var_hhid var_subid s2
matrix list estimates
```

Inside `ds_run()`, `ds_build()` prepares household summaries, `ds_solve()` estimates the two variances, and `ds_sandwich()` computes uncertainty for an interior solution. The output retains the first-stage fixed effects and residual variance, substitutes the design-score household and subject variances, and reports standard errors from the stacked calculation when its checks pass.

## Results from 1000 Oct2026 samples

All 1,000 fits produced usable estimates and intervals. The principal comparison is:

| Population and method | Subject-variance bias | Relative bias | Coverage |
|---|---:|---:|---:|
| Dec2024 original scenario 3 | −0.135 | −0.115% | 94.6% |
| Oct2026 original scenario 3 | +6.000 | +3.425% | 83.0% |
| Oct2026 design score | −0.020 | −0.0116% | 93.8% |

Sources: the [Dec2024 original summary](../data/summaries/dec2024/coverage_scen3.csv), [Oct2026 original summary](../data/summaries/oct2026/coverage_scen3.csv), and [design-score summary](summaries/oct2026/coverage_scen3_designscore_exact.csv).

The design-score variance estimates were:

| Component | Population target | Mean estimate | Bias | Coverage |
|---|---:|---:|---:|---:|
| Household variance | 22.738 | 22.550 | −0.188 | 95.5% |
| Subject variance | 175.174 | 175.154 | −0.020 | 93.8% |
| Residual variance | 95.090 | 95.112 | +0.023 | 94.6% |

Subject-variance bias is small relative to its Monte Carlo standard error of 0.198. Coverage remains approximate; 93.8% alone does not establish a departure from 95% with 1,000 replicates. The coverage improvement also reflects a changed standard-error calculation, not just the change in point estimation. These comparisons use saved summaries, rather than a verified paired analysis of replicate inputs.

## Scope and implementation

The current implementation uses exactly three complete visits per person and deliberately ignores the MAR visit flag. Applying that flag requires extending the residual moments, first-stage estimation, and uncertainty calculation to account for visit nonresponse. Baseline-covariate missingness and other covariance models also require further work.

The implementation is in [scen3_designscore.do](code/stata/scen3_designscore.do) and [variance_designscore.mata](code/stata/variance_designscore.mata). The folder [README](README.md) contains execution instructions and file locations; the [audit](summaries/oct2026/audit_scen3_designscore_exact.csv) contains detailed numerical checks. The original scenario 3 code remains unchanged.
