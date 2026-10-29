Age-specific census variance diagnostic, 29 October 2026

This folder preserves the evidence used in DESIGNSCORE_METHOD_REPORT.md.
The original calculation was run separately from the active Stata MI batch.
The source copies here differ only in repository-root discovery so that they
can run from this archived location. No production estimation code was changed.

Read results_summary.txt for model details and the principal comparison.
age_variance_results.R contains full-precision estimates, fixed effects,
likelihoods, analytic-score convergence diagnostics, and baseline validation.
validation.txt is the original successful calculation log.
input_*.json record the population input paths, file metadata, and counts.

The sufficient_*.R files encode exact cross-products, within-subject contrasts,
and household sums grouped by younger/older household composition. They allow
the Gaussian ML models to be refitted without rereading the large populations.
validation_subset_*.txt contain 40 synthetic subjects per population and are
used only to check the compressed likelihood against direct matrix algebra.
Final model fits use the complete-population sufficient statistics.

To reproduce from designscore_diagnostics, with Python (numpy and pandas) and
base R available:

  python provenance/age_specific_variance_20261029/prepare.py
  Rscript --vanilla provenance/age_specific_variance_20261029/fit.R provenance/age_specific_variance_20261029

The Python preparation step can be skipped when using the saved sufficient
statistics. Preparation reads populations from the parent repository's data
folder and writes only within this evidence folder. On Windows, it limits its
own process to one CPU at below-normal priority. The original R process was
also launched with those limits. Apply the same limits to an R rerun while
another analysis is active. The scripts do not launch or manage Stata.

The result is a descriptive census contrast under an expanded scenario 3
working model. It is not a comparison of estimated subject-specific BLUPs,
nor a claim to identify the original latent generating variances.
