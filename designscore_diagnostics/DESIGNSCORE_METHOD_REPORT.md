# Design score method for the Oct2026 mixed model

**Method report — 29 October 2026**

The design-score method changes how scenario 3 estimates household and subject variances. It weights residual squares and cross-products to reproduce the estimating equations for the complete population. In 1,000 Oct2026 samples, it reduced subject-variance bias from **+6.000 to −0.020**; coverage improved from **83.0% to 93.8%**.

## Why the original scenario 3 is biased for Oct2026

The original method estimates fixed effects and residual variance first, then holds them fixed while estimating household and subject variances. Its second fit uses `sub_s1`: subject weights rescaled within each sampled household to an effective sample size. For a household with one sampled person, that weight becomes one regardless of the person's selection probability.

Variance estimation depends on both individual residual squares and products between different subjects. Effective-size weights do not generally give these terms the individual and pair inclusion weights needed to recover the population equations. They also change how household size enters the likelihood. This can shift the variance estimates away from the population target.

The populations differ in ways that can expose this problem. [Dec2024](https://github.com/UNC-CSCC/HCHS-Simulations/blob/e9e8128/code/v3/Dec2024/population/generate_population.R) uses an additive outcome generator with Gaussian random effects. [Oct2026](https://github.com/UNC-CSCC/HCHS-Simulations/blob/e9e8128/code/v3/Oct2026/population/README.md#how-renal-calibration-and-clustering-are-preserved-together) calibrates and reassigns whole three-visit trajectories, changing their associations with covariates and clusters. Its calibration model also differs from scenario 3, including block-group effects and visit-specific residual variances. Good performance in December therefore does not guarantee good performance in October; neither generator exactly matches the fitted scenario 3 model.

The saved checks support the weighting explanation: adding stratum indicators leaves the subject-variance bias near 3.4%, while unweighted fits have much smaller bias. Simply removing weights is not a general solution, however: the unweighted age effect has 17.0% relative bias and only 24.2% coverage. These comparisons implicate the weighting procedure without identifying the exact contribution of each population change. See the [stratum-adjusted](../data/summaries/oct2026/coverage_scen3_strata.csv), [household SRS](../data/summaries/oct2026/coverage_scen3_srs.csv), and [unweighted](../data/summaries/oct2026/coverage_scen3_unw.csv) summaries.

### Age-specific subject variance in the two populations

A separate census diagnostic fitted both complete populations, each containing **1,266,981 subjects in 633,602 households**, with three visits per subject. Younger means baseline age below 45; older means age 45 or above. The unweighted Gaussian maximum-likelihood model retained scenario 3's fixed-effect specification and a shared household random intercept, while allowing **both subject random-intercept variance and visit-level residual variance to differ by age group**. Fixed effects were jointly reestimated. Mixed-age households retained their shared household effect; no block-group random effect was added. This is a diagnostic extension of the outcome model, not a change to the design-score estimator.

| Population | Subject variance, age <45 | Subject variance, age ≥45 | Older / younger | Relative difference |
|---|---:|---:|---:|---:|
| Dec2024 | 117.425 | 117.138 | 0.998 | −0.24% |
| Oct2026 | 156.312 | 205.267 | 1.313 | +31.32% |

All entries are **variances**, not standard deviations. The corresponding nuisance variance estimates were:

| Population | Shared household variance | Visit residual variance, age <45 | Visit residual variance, age ≥45 |
|---|---:|---:|---:|
| Dec2024 | 38.876 | 57.324 | 56.228 |
| Oct2026 | 22.509 | 86.094 | 109.191 |

The fitted subject variance is essentially equal by age in Dec2024, but about **31% greater among older subjects in Oct2026**, even after allowing older subjects to have greater visit-level noise. A sensitivity fit with a common visit-residual variance gives subject variances of 117.565 versus 116.914 in Dec2024 and 153.313 versus 209.967 in Oct2026; the conclusion is unchanged.

This finding supports an interaction between age-dependent variability, unequal age-group selection probabilities (nominally 0.3575 for younger and 0.55 for older subjects), and scenario 3's within-household weight scaling. The group with higher October variability is sampled more often, and effective-size scaling does not generally restore the population variance equations. These checks do **not** isolate what fraction of the original six-unit subject-variance bias comes from this mechanism. The reported components are fitted parameters under the expanded working model, not identified latent generating variances; these are descriptive census contrasts, with no superpopulation hypothesis test.

The calculation used an exact Gaussian likelihood compressed by household age composition. Direct dense-covariance checks reproduced the likelihood and profiled fixed effects on small observed-data subsets; independent optimization starts agreed after analytic-score refinement. A common-variance fit reproduced the saved Stata census variances within **0.0025%** and fixed effects within **0.000002**. See the [full age-specific estimates and validation records](provenance/age_specific_variance_20261029/age_variance_results.R), [calculation log](provenance/age_specific_variance_20261029/validation.txt), and [reproduction instructions](provenance/age_specific_variance_20261029/README.txt). The diagnostic was run separately on one logical CPU at below-normal priority, without launching or modifying any Stata MI process.

### Why design-score is a more general correction

**For the sampling-and-weighting bias at issue here, design-score is a broader correction than adding age-specific subject variances.** The current estimator still targets one common household variance and one common subject variance in the original census working model. Age-specific variances describe one particular source of outcome heterogeneity and change the model's parameters of interest. They do not, by themselves, ensure that the scaled sample likelihood recovers the corresponding census equations.

Design-score instead corrects the equations directly: residual squares and subject counts receive inverse individual inclusion weights, within-household cross-products receive inverse pair inclusion weights, household contributions receive inverse household inclusion weights, and covariance denominators use full eligible household sizes. At fixed parameter values, their sampling expectation equals the complete-population variance score. This identity can hold even when residual variability is associated with age or other selection-related characteristics; identifying age as the source of heterogeneity is not required. It targets the common-variance census fit even when the working model does not exactly describe the outcome distribution.

The two approaches are complementary. If age-specific variance components became the scientific targets, the design-score equations would need to be extended for that outcome model. The present correction requires valid selection probabilities and household sizes, and its implementation assumes the specified Gaussian model and three complete visits. Design-unbiased scores do not guarantee unbiased finite-sample roots, remove error from the retained first stage, validate an imputation model, or guarantee nominal interval coverage.

## Model and estimation target

### Overview

The method uses a sample to reconstruct the household and subject variance equations for the complete population. It first retains scenario 3's estimates of the fixed effects and visit-level residual variance. It then averages each person's three residuals and weights their squares and cross-products using individual and pair selection probabilities. The Mata solver uses those totals to estimate household and subject variances jointly, and a stacked sandwich calculation carries uncertainty from both stages into the confidence intervals.

The information flow is **subject-average residuals → household estimating equations → block-group score totals → standard errors**. The first two steps produce point estimates; the last two combine information from both estimation stages to quantify uncertainty.

![Subject-average residuals become weighted individual totals and within-household pair products, which enter the explicit household and subject variance equations. Both estimation stages are stacked by household, summed within block groups, and converted to standard errors using a stratified sandwich.](figures/two_stage_designscore.png)

**Figure 1. From subject-average residuals to standard errors.** The upper panels show point estimation: stage 1 supplies fixed effects and visit-residual variance; stage 2 uses weighted counts, residual squares, and products between different subjects in the same household to solve the displayed household-variance and subject-variance equations. The lower panels show uncertainty: weighted household contributions from both stages form a 13-element stack, household stacks are summed within block groups, and their within-stratum variation enters the sandwich. Full eligible household size remains in the covariance denominator. Design strata are the sampling strata, not the younger/older groups. The fixed-parameter expectation identity does not automatically apply at fitted values. [Scalable SVG](figures/two_stage_designscore.svg) · [Figure source](code/evaluation/draw_two_stage_designscore.R).

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

For fixed `β`, define one residual at each of subject `i`'s three visits, then average them:

```math
R_{hit}=Y_{hit}-x_{hit}^{T}\beta,\qquad t=1,2,3,
\qquad
r_{hi}=\frac{R_{hi1}+R_{hi2}+R_{hi3}}{3}.
```

In the fitted estimator, substitute the stage-1 estimate `β̂`. These residuals subtract the fixed-effect prediction only; the household and subject random effects remain in them. In an MI analysis, the calculation uses the current imputation's completed data. The three residuals are three visit-level discrepancies, not three types of random effect.

Under the working model with its mean parameter `β`, the subject average has the decomposition `r_hi = u_h + v_hi + ε̄_hi`. Consequently,

```math
E(r_{hi}^2)=a+b+e/3,\qquad
E(r_{hi}r_{hj})=a\quad(i\ne j\text{ in the same household}).
```

Individual squares contain household variance, subject variance, and averaged visit noise. Cross-products between different subjects contain their shared household component. These model identities explain how the terms distinguish variance components; the finite-population design-expectation result below is a separate statement about the sampling weights.

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

### Individual totals and within-household pairs

**Individual totals** add one contribution per sampled subject within household `h`. For sampled subjects `s_h`, construct

```math
\hat{N}_h=\sum_{i\in s_h}\frac{1}{p_{hi}},\qquad
\hat{B}_h=\sum_{i\in s_h}\frac{r_{hi}^2}{p_{hi}}.
```

Here `p_hi` is the probability that subject `i` is selected, conditional on the selected stage-2 household set. `N̂_h` estimates the eligible subject count; `B̂_h` estimates the total of squared subject-average residuals in the full household. They are household totals formed from individual contributions, not a separate fitted parameter for each person.

**Within-household pairs** are two different sampled subjects `i` and `j` from the same household. They are not two visits from one subject. Their products enter

```math
\hat{Q}_h=\hat{B}_h+2\sum_{\substack{i<j\\i,j\in s_h}}\frac{r_{hi}r_{hj}}{p_{hij}}.
```

The joint probability `p_hij` is the probability that **both subjects** are selected, conditional on the same stage-2 household set. `Q̂_h` estimates the full household's squared residual sum, `Q_h = (Σ_i r_hi)^2`, including both the individual squares and distinct-person products. With one sampled person, the pair sum is empty and `Q̂_h = B̂_h`; households with multiple sampled subjects supply the cross-person information.

The individual squared terms require first-order inclusion probabilities; the distinct-person cross-products require joint inclusion probabilities. This follows the general design-weighted estimating-equation principle discussed by [Rao, Verret, and Hidiroglou, section 3](https://www150.statcan.gc.ca/n1/pub/12-001-x/2013002/article/11887/section3-eng.htm). The particular equations here are derived from the full Gaussian census likelihood above; they are not a direct implementation of that paper's composite likelihood.

The pair-weighted quantity `Q̂_h` is not the square of an estimated total. For example, under an illustrative independent subject draw with probability one-half, a residual square receives weight two, while a distinct-person product receives weight four. Squaring an inverse-probability estimate of the residual total would incorrectly give the squared terms weight four as well. In a particular sample, `Q̂_h` can also be negative even though its census target `Q_h` is nonnegative.

### What solving the two variance equations means

The two unknowns are **`a`, the common household variance, and `b`, the common subject variance**. For stage-2 point fitting, hold `β̂` and `ê` fixed, calculate the residual totals above, and use `d = b + ê/3`, `D_h = d + a M_h`, and `W_h = 1/π_h`, where `π_h` is the household selection probability. For a positive interior solution, the solver finds `a` and `b` satisfying both equations:

```math
\hat{U}_a(a,b)=\frac{1}{2}\sum_{h\in s_{12}}W_h
\left[\frac{\hat{Q}_h}{D_h^2}-\frac{\hat{N}_h}{D_h}\right]=0.
```

and

```math
\hat{U}_b(a,b)=\frac{1}{2}\sum_{h\in s_{12}}W_h
\left[\frac{\hat{B}_h}{d^2}-\frac{a(2d+aM_h)\hat{Q}_h}{d^2D_h^2}-\frac{\hat{N}_h}{d}+\frac{a\hat{N}_h}{dD_h}\right]=0.
```

`Û_a` is the sampling-corrected household-variance score; `Û_b` is the sampling-corrected subject-variance score. They are derivatives of the weighted objective shown below with respect to `a` and `b`. Each compares weighted residual information with what the candidate variances predict. A positive score locally favors increasing its variance while holding the other fixed. Solving both equations balances those discrepancies **across households**; it does not require each household's contribution to equal zero.

Write `Û_{a,h}` and `Û_{b,h}` for the two weighted household summands, including their `W_h/2` factors. Thus `Û_a = Σ_h Û_{a,h}` and `Û_b = Σ_h Û_{b,h}`. These same household contributions become the last two entries of the stacked estimating vector used for standard errors. At a zero-variance boundary, the solver uses the one-sided score conditions described under numerical solution rather than requiring both scores to be zero.

Here `s₁₂` denotes the households selected at the first two sampling stages. Households with no selected person have `N̂_h = B̂_h = Q̂_h = 0`, so they need not appear as rows in the model-fitting dataset. In fixed-parameter derivations, `e` denotes the fixed candidate residual variance; in the actual stage-2 fit it is set to `ê`.

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

Holding `β̂` and `ê` fixed during point estimation does not make them known without error. The standard-error calculation therefore follows the path **subject-average residuals → household estimating equations → block-group score totals → standard errors**, while retaining the full first-stage information.

The parameter vector has 13 entries:

```math
\theta=(\beta^{T},\log(a_1),\log(b_1),\log(e),\log(c),\log(a_2),\log(b_2))^{T}.
```

There are seven fixed effects, three first-stage covariance parameters, the global subject-weight normalizer, and two second-stage variances. The household estimating vector contains the corresponding first-stage scores, the normalizer equation

```math
g_{c,h}=n_h-c\sum_{i\in s_h}w_{hi}.
```

and the two design scores. Here `n_h` is the number of sampled subjects in the household. The normalizer equation is not multiplied by the household weight, matching the global normalization used in the Stata fit.

More explicitly, let `S_h^(1)` contain the ten weighted first-stage equations (seven fixed effects and three variances). Then each household contributes

```math
\psi_h(\theta)=
\begin{pmatrix}
S_h^{(1)}\quad\text{(10 entries)}\\
g_{c,h}\quad\text{(1 entry)}\\
\hat U_{a,h}\quad\text{(1 entry)}\\
\hat U_{b,h}\quad\text{(1 entry)}
\end{pmatrix}.
```

The last two entries use the subject-average residual squares and within-household pair products defined above. The first-stage equations additionally retain within-subject visit deviations and covariate information; averaging residuals does not discard those quantities from the full stack. The household variance, subject variance, and residual variance in the first-stage equations are its own nuisance parameters `a₁`, `b₁`, and `e`, while the last two equations target `a₂` and `b₂`.

Let `ψ_h(θ)` be this stacked household estimating vector. The Jacobian is

```math
A(\theta)=-\frac{\partial}{\partial\theta^{T}}\sum_h\psi_h(\theta).
```

The Jacobian is evaluated at the estimated parameter vector `θ̂`. The implementation obtains it by central differences of the analytic household scores. In this calculation, perturbing `β` recomputes the residual totals, and perturbing `e` changes the stage-2 denominators. These cross-stage derivatives propagate first-stage estimation error into the final variance-component uncertainty. Equation rescaling is used for numerical stability.

### Block groups and design strata

Block groups are treated as primary sampling units for the sandwich calculation. Within design stratum `g`, let `m_g` be the number of observed block groups. First **sum**, rather than average, the already weighted household vectors within each block group `k`, then center the resulting block-group totals within their design stratum:

```math
\Psi_{gk}=\sum_{h\in(g,k)}\psi_h,\qquad
\bar\Psi_g=\frac{1}{m_g}\sum_{k=1}^{m_g}\Psi_{gk},\qquad
\Delta_{gk}=\Psi_{gk}-\bar\Psi_g.
```

The strata here are **sampling-design strata**, not the younger/older groups. All thirteen entries move together as one block-group vector, allowing dependence among subjects and households in that block group and covariance between the two estimation stages. The meat of the sandwich measures the within-stratum variability of these totals:

```math
\mathcal M=\sum_g\frac{m_g}{m_g-1}\sum_{k=1}^{m_g}
\Delta_{gk}\Delta_{gk}^{T}.
```

`𝓜` is a 13-by-13 score covariance estimate, distinct from the scalar residual total `B̂_h` and the full household size `M_h`. The bread `A⁻¹` converts score fluctuations into parameter fluctuations. The estimated covariance matrix of `θ̂` is

```math
\widehat V_\theta=A^{-1}\mathcal M(A^{-1})^{T},\qquad
\operatorname{SE}(\hat\theta_j)=\sqrt{(\widehat V_\theta)_{jj}}.
```

For logged variance parameters, these are log-scale standard errors; the delta method and log-scale intervals below give the reported variance-scale results. The stacked covariance includes cross-stage dependence, so the second-stage standard errors do not treat the first-stage point estimates as known constants.

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

The complete-data wrapper uses exactly three complete visits per person and deliberately ignores the MAR visit flag. The separate paired-MI extension applies the same two-stage estimator to each completed imputation and pools estimates and uncertainty using Rubin's rules; its validity additionally depends on the imputation model. The complete-data results and age-specific census diagnostic in this report are not MI results. Baseline-covariate missingness and other covariance models require further work.

The implementation is in [scen3_designscore.do](code/stata/scen3_designscore.do) and [variance_designscore.mata](code/stata/variance_designscore.mata). The folder [README](README.md) contains execution instructions and file locations; the [audit](summaries/oct2026/audit_scen3_designscore_exact.csv) contains detailed numerical checks. The original scenario 3 code remains unchanged.
