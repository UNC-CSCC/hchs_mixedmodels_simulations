# Scenario 3: design-weighted census variance scores

The original `scen3.do`, original estimates and population truths are unchanged.
Stata 19 is required. The Mata code has been compiled and the local October
samples 1–1,000 have completed with usable estimates and intervals. Local batch
results are stored separately from the original scenarios.

The variant retains the original globally scaled Step 1 and replaces the
constrained `meglm` Step 2 with design-weighted scores for the same Gaussian
census likelihood, conditional on Step 1's fixed effects and residual variance.
It targets `var_hhid` and `var_subid` jointly. It does not impose the observed
October bias as a correction or read population outcome values when fitting.

## Folder organization

All paths below are relative to `hchs_mixedmodels_simulations/designscore_diagnostics/`.
Run the examples from this folder. The Stata argument named `repo_root` and the
Python `--root` option both refer to this diagnostics folder.

```text
code/
  stata/          model, Mata solver, local batch runner
  preparation/    household frame and exact sampling-count recovery
  execution/      resumable local extension from 200 to 1,000 samples
  evaluation/     summary generator
  verification/   checks of saved results and independent R summary verification
  tests/          model-free algebra and summary-denominator tests
  archive/        historical wrappers and checks, retained as run records
data/
  design_frames/  household lookup, exact counts, cached population design columns
  estimates/oct2026/exact/  1,000 coefficient CSVs and 1,000 audit CSVs
  population_coefficients/ unchanged copies of the existing population truths
  run_artifacts/  original smoke-test output
  input_locations.json     references to shared population and sample inputs
logs/
  oct2026/exact/   1,000 per-replicate Stata logs
  local_batches/  batch logs and progress manifests
  local_runs/     sampling-count recovery and supervisor logs
  validation/     checks performed after the reorganization
summaries/oct2026/
  coverage_scen3_designscore_exact.csv
  audit_scen3_designscore_exact.csv
  provenance_scen3_designscore_exact.json
provenance/
  relocation_manifest.json  old/new paths and before/after SHA-256 hashes
  pre_reorganization/       exact source snapshots from before the move
  local_runs/               historical run status and verification records
.statatmp/                  retained scratch files; new runtime scratch
```

The population and 1,000 source sample CSVs remain in the shared HCHS data
folder, at `../../data/Oct2026/`. They were not copied or moved. Their locations
and the production sampling-code location are recorded in `data/input_locations.json`.
Derived data, estimates and logs are kept locally and ignored by Git; code,
summaries and provenance remain available for version control.

The reorganization did not rerun any models or change any estimates. Original
summary/provenance JSON files were retained byte for byte. Their historical
`source_sha256` paths refer to `provenance/pre_reorganization/`; other old paths
can be resolved using `provenance/relocation_manifest.json`. The active Stata
wrapper has path/log-location changes; the Mata solver has only a load-path
comment update. All estimating equations are unchanged. Files in
`code/archive/` are historical records with their original paths, not current
entry points. Future runs use `code/stata/` and `code/execution/`.

## Files and interface

- `code/stata/scen3_designscore.do`: validation, original Step 1, solver invocation and export.
- `code/stata/variance_designscore.mata`: household sufficient statistics, constrained
  numerical solution and stacked sandwich calculation. Loading it only defines
  functions; it does not fit a model.
- `code/preparation/prepare_designscore_frame.py`: builds household frame information
  from design columns in a full saved population; no model fitting.
- `code/tests/test_designscore_algebra.py`: independent, model-free checks of the
  score identities, likelihood derivatives and selection-probability weighting.
  These do not execute or compile the Mata implementation.
- `code/preparation/recover_designscore_counts.R`: replays only the sampling stage
  from saved completion-record seeds; verifies population/sample checksums and
  every sampled subject ID before publishing exact counts. It does not change
  the saved samples or generate new outcome/missingness values.
- `code/stata/run_designscore_local.do`: runs a specified replicate range, saves a progress
  manifest, skips existing coefficient files, records nonconvergence and stops
  the range on other errors. Its arguments are `first last mode [counts_csv]`;
  the optional count file supports separate local batches.
- `code/execution/extend_designscore_1000.py`: resumes the local October study
  from 200 to 1,000 samples. It verifies the previous inputs and estimates,
  recovers exact counts for the additional samples, runs local Stata batches,
  then writes a combined summary. Existing completed fits are reused.
- `code/evaluation/summarize_designscore.py`: summarizes a specified replicate
  range, including missing outputs and interval denominators.

The Stata argument order is:

```text
sim_id population probability_mode [repo_root] [sample_directory] [frame_csv] [counts_csv]
```

`population` is `dec2024` or `oct2026`. `probability_mode` must be supplied as
`exact` or `bernoulli`; the script never silently falls back to an approximation.
Use a fresh Stata session when updating the Mata source: an already loaded
`ds_run()` is reused during repeated calls in the same session.

Defaults when submitted from the `designscore_diagnostics` folder:

```text
repo_root         current working directory (designscore_diagnostics)
sample_directory  ../../data/Oct2026/sample_miss00pct  (or Dec2024)
frame_csv         data/design_frames/oct2026_households.csv
counts_csv        data/design_frames/oct2026_subject_counts.csv
```

Supply the sample directory explicitly on the cluster. Relative custom paths
are resolved against Stata's current directory; paths with spaces need quotes.

## Household frame

Prepare the lookup once, separately from any fit. The October lookup has been
prepared locally from the saved population using this command:

```bash
python code/preparation/prepare_designscore_frame.py \
  --population oct2026 \
  --population-file ../../data/Oct2026/population_3visits_Oct2026.csv \
  --output data/design_frames/oct2026_households.csv
```

The script requires complete subjects sorted by `BGid, hhid, subid, v_num`.
It uses only identifiers, stratum, surname group and age. Required lookup columns:

| Column | Meaning |
|---|---|
| `population` | Release identifier, `dec2024` or `oct2026` |
| `bgid`, `hhid` | Household key |
| `strat` | Original design stratum 1–8 |
| `m_h` | Full eligible household size, before subject sampling |
| `p_h` | Actual BG-selection probability times conditional HH-selection probability |

The supplied builder also saves `n_younger` and `n_older`, plus a JSON source
fingerprint. It implements the existing HCHS stage 1/2 fractions and their
rounding, validates the eight stratum BG counts and refuses to overwrite
existing output. Derived lookups are ignored by Git.

`p_h` is the probability of selecting a household at stages 1/2, before any
subjects are selected. It is not conditional on that household appearing in
the final sample CSV. Step 1 retains the original nominal `w_bghh`; Step 2 uses
`1/p_h` from the lookup, including rounding. This small distinction is intentional.

## Subject-probability modes

### Exact point-estimation probabilities

The `counts_csv` must contain one uniquely identified row per population and
replicate, with these columns:

```text
population,sim_id,eligible_younger,selected_younger,eligible_older,selected_older
```

`eligible_*` counts all people in the selected stage-2 households, including
households with no person in the final sample. `selected_*` counts the people
drawn in stage 3, before missingness. These must come from retained sampling
metadata, not be inferred from only the nonempty sampled households. The
original generator is not modified by this variant, and the current sample
CSVs alone do not provide those exact eligible counts.

The local October samples also have `.complete` sidecars with file-specific
seeds and checksums. Their counts can be recovered and verified without
regenerating samples:

```bash
Rscript code/preparation/recover_designscore_counts.R ../.. data/design_frames/oct2026_subject_counts.csv
```

The helper replays the production sampling draws with each recorded seed,
using cached household candidate lists and the production sampling helper. It
checks the replay against the full production sampling function for the first
replicate in each invocation. It captures the two stage-3 eligible/selected
counts and requires every replayed subject ID to match the saved sample. It
stops on any mismatch. A cached baseline design table is keyed by the checked
population checksum. Optional third and fourth arguments specify the last
and first replicate, respectively.

For an age stratum with `N` eligible and `n` selected people, the script uses
`p=n/N` and `p_pair=n(n-1)/(N(N-1))`. For a mixed-age pair it uses the product
of the two marginal probabilities, because the age-stratum draws are separate.
The implementation requires at least two selected people in each age stratum
and checks selected counts against the sample. This is satisfied by the intended
HCHS samples; samples with zero pair-inclusion probabilities are rejected.

### Explicit Bernoulli approximation

`bernoulli` uses nominal probabilities 0.3575 for younger and 0.55 for older
subjects, and products of these probabilities for pairs. It needs the same full
household lookup but no subject-count metadata. This approximates the actual
fixed-size sampling and its rounding. Every output records `probability_mode`.
Do not combine exact and approximate outputs as one simulation scenario.

Future invocation examples, from the `designscore_diagnostics` folder:

```stata
do code/stata/scen3_designscore.do 1 oct2026 bernoulli
do code/stata/scen3_designscore.do 1 oct2026 exact
```

The second command fails before any model fit if exact metadata is missing.
For cluster paths, pass the diagnostics root and sample directory explicitly:

```stata
do code/stata/scen3_designscore.do 1 oct2026 bernoulli "/path/to/hchs_mixedmodels_simulations/designscore_diagnostics" "/work/users/b/e/beibo/v3data/Oct2026/sample_miss00pct"
```

## Census-score construction

There must be exactly three complete visit records per subject. This remains
the complete-data scenario 3 benchmark: the retained simulated values at visits
flagged missing are used. This is not a missing-data or MI analysis. Positive
baseline missingness, incomplete visits, duplicate keys and missing analysis
values cause failure, rather than silent record deletion.

For subject `i` in household `h`, let `r_hi` be the average of its three residuals
after subtracting Step 1's fixed-effect prediction. Random-effect predictions
are not subtracted. Let `M_h` be the full household size and `W_h=1/p_h`.

The household summaries are:

```text
Nhat_h = sum_i 1/p_i
Bhat_h = sum_i r_hi^2/p_i
Qhat_h = Bhat_h + 2 sum_(i<j) r_hi*r_hj/p_ij
```

With `a=var_hhid`, `b=var_subid`, `e=Step-1 s2`, `d=b+e/3`, and
`D_h=d+a*M_h`, the scores are:

```text
U_a = 0.5 sum_h W_h [Qhat_h/D_h^2 - Nhat_h/D_h]

U_b = 0.5 sum_h W_h [Bhat_h/d^2
                    - a*(2*d+a*M_h)*Qhat_h/(d^2*D_h^2)
                    - Nhat_h/d + a*Nhat_h/(d*D_h)]
```

Using `Nhat_h` in the count terms accounts for households with no selected
subjects, which have zero sample contribution. Actual `M_h` remains in covariance
denominators. A singleton's individual contribution retains its inverse
selection probability; a household's outer selection weight is applied once.

At fixed beta and e, and with correct probabilities, the design expectation of
these scores equals the full-population conditional Gaussian scores. This does
not assert finite-sample unbiased roots or exact equality of the two-step
estimator and a joint census ML fit. The quality of Step 1 and finite-sample
behavior still require evaluation.

The Mata implementation accumulates quadratic sufficient statistics so beta can
vary during sandwich differentiation without re-reading the data. It maximizes
the integrated score objective

```text
L = -0.5 sum_h W_h [Nhat_h*log(d) + (Nhat_h/M_h)*log(D_h/d)
                   + Bhat_h/d - a*Qhat_h/(d*D_h)]
```

using `optimize()` with analytic gradients. Its derivatives are the two scores
above. This is a numerical implementation of the proposed equations, not a
generic pairwise composite-likelihood replacement. It checks three starting
points, both one-variance-zero boundaries and the corner. Raw score/KKT checks
are required, so tiny log variances cannot falsely indicate convergence.

## Standard errors and status

The interior solution uses a 13-parameter stacked sandwich:

1. Seven Step 1 fixed effects; the three Step 1 covariance parameters are also
   included as nuisance parameters.
2. The sample-dependent global weight normalizer, with its own estimating equation.
3. The two new Step 2 variance parameters.

The Jacobian is obtained by central differences of analytic household scores.
The meat aggregates scores by BG within the eight original strata, using a
with-replacement ultimate-cluster approximation. Step 1 uncertainty and its
covariance with Step 2 are propagated. Output SEs for fixed effects and `s2`
also come from this stack, so they need not equal the original household-cluster
SEs even though their point estimates are retained.

**Inference remains approximate in both probability modes.** The sandwich does
not include exact multistage finite-population corrections or every dependence
induced by fixed-size subject sampling across BGs. `exact` describes the point
estimation probabilities, not exact confidence-interval coverage. Every replicate
must pass the numerical score checks; empirical coverage is assessed separately.

Variance intervals use the log scale at interior solutions. Boundary solutions
are exported as point estimates only. Missing SEs are never replaced with zero.
The audit `se_status` codes are:

| Code | Meaning |
|---|---|
| 0 | Stacked sandwich computed |
| 1 | Nonfinite stacked scores |
| 2 | Step 1 scores do not match the fitted point sufficiently closely |
| 3 | Singular or ill-conditioned stacked Jacobian |
| 4 | Fewer than two observed PSUs in a stratum |
| 5 | Nonfinite covariance matrix or negative diagonal |
| 6 | Boundary solution; ordinary interior inference withheld |
| 7 | Nonfinite SE or confidence limit after transformation |

## Outputs and verification

Outputs are isolated under:

```text
data/estimates/<population>/<exact|bernoulli>/
  scen3_designscore_coef_<rep>.csv
  scen3_designscore_audit_<rep>.csv
logs/<population>/<exact|bernoulli>/
  scen3_designscore_<rep>.log
```

Existing coefficient/audit files are protected from replacement. Coefficient
tables retain `term,estimate,std_error,conf_low,conf_high`, adding probability
mode, fit/inference status and SE method. The audit includes the selected
boundary, KKT error, objective, successful starts, Step 1 score check, Jacobian
condition, sample counts, probabilities and Step 1 covariance estimates.

Original coverage scripts deliberately have no preset for this experimental
variant. The dedicated local evaluator distinguishes point-estimate and
valid-interval denominators, reports missing outputs, and keeps probability modes
separate. For a batch of 1,000 samples, run from the `designscore_diagnostics` folder:

```stata
do code/stata/run_designscore_local.do 1 1000 exact
```

```bash
python code/evaluation/summarize_designscore.py --root . --mode exact --first 1 --last 1000
```

The summary, replicate audit and source/output fingerprints are written under
`summaries/oct2026/`, with filenames containing `scen3_designscore_exact`.
`coverage_rate` uses valid intervals as its denominator; `coverage_over_all_requested`
counts missing intervals as failures over all requested replicates. Bias and
coverage Monte Carlo standard errors are also reported.

Model-free algebra checks can be run separately:

```bash
python -m unittest discover -s code/tests -p test_designscore_algebra.py -v
```

These checks alone cannot certify Stata syntax, optimizer behavior, sandwich
accuracy or sampling performance. A Stata smoke run and the per-replicate audit
provide separate runtime checks.

Method references: [design-weighted estimating equations](https://www150.statcan.gc.ca/n1/pub/12-001-x/2013002/article/11887/section3-eng.htm),
[Stata's weighted mixed-model likelihood](https://www.stata.com/manuals/memeglm.pdf),
and [Mata optimize()](https://www.stata.com/manuals/m-5optimize.pdf).

## Full 1,000-sample October run, 29 October 2026

The local extension fitted saved samples 201–1,000 and combined them with the
unchanged estimates for samples 1–200. Stata/BE 19.5 ran sixteen ranges of 50
additional samples, with eight concurrent batches. The estimator source and
population truth were unchanged from the 200-sample run. Checksums confirmed
that the original 200 input samples and coefficient files were unchanged.

All 1,000 fits converged to interior solutions, and all 10,000 coefficient rows
have usable standard errors and confidence intervals. No requested replicate,
point estimate or interval is missing. Exact selection-count records for all
1,000 samples are in `data/design_frames/oct2026_subject_counts.csv`; every
replayed subject ID and sample checksum was verified before fitting.

| Variance component | Population truth | Mean estimate | Empirical bias | Coverage |
|---|---:|---:|---:|---:|
| `var_hhid` | 22.738064 | 22.549958 | −0.188107 | 95.5% |
| `var_subid` | 175.174303 | 175.153937 | −0.020366 | 93.8% |
| `s2` | 95.089839 | 95.112469 | 0.022630 | 94.6% |

The `var_subid` relative bias is −0.0116%, with bias Monte Carlo SE 0.198022.
Its empirical SE is 6.261990, compared with mean model SE 6.007524 (4.06%
lower). Its coverage Monte Carlo SE is 0.763 percentage points. The bias is
near zero in this run; interval coverage is 1.2 percentage points below the
nominal 95%. The sandwich remains the approximate method described above.
The archived original scenario summary reported bias +6.000029 and coverage
83.0%; that comparison does not establish that both runs used identical sample
files.

The maximum scaled KKT error was 4.063e-7, maximum Step 1 score error was
7.555e-6, and maximum Jacobian condition number was 13,892.60. All samples
contained 220 observed BGs, with 14,200–16,559 sampled subjects. Independent R
recalculation of means, empirical/model SEs, coverage rates and denominators
agreed with the Python summary within 1e-10. Final integrity checks verified
all exact count records, numerical diagnostics, interval fields and coefficient
fingerprints, along with the unchanged original `scen3.do`.

The combined outputs are:

- `summaries/oct2026/coverage_scen3_designscore_exact.csv`
- `summaries/oct2026/audit_scen3_designscore_exact.csv`
- `summaries/oct2026/provenance_scen3_designscore_exact.json`

Per-replicate estimates remain under `data/estimates/oct2026/exact/`,
with Stata logs under `logs/oct2026/exact/`.

## Verify saved results without fitting models

```bash
python -B -m unittest discover -s code/tests -p "test_designscore_*.py" -v
python -B code/verification/verify_saved_run.py
Rscript --vanilla code/verification/verify_summary.R
```

To recompute the summary into a separate folder, without changing the saved
summary or fitting any model:

```bash
python -B code/evaluation/summarize_designscore.py --root . --mode exact --first 1 --last 1000 --output-dir .statatmp/summary_check
```
