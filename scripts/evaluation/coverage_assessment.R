#!/usr/bin/env Rscript
# Summarize coefficient bias, standard errors, coverage rates, and rejection probabilities.
# Presets: scen9 (dec2024, oct2026) and scen15/scen16 MI variants v1-v3 (oct2026).
# Usage: Rscript coverage_assessment.R <preset>
#   preset: 9_dec2024 | 9_oct2026 | 15_v{1,2,3}_oct2026 | 16_v{1,2,3}_oct2026

repo_root <- "/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"

args <- commandArgs(trailingOnly = TRUE)
scen_num <- if (length(args) > 0) args[[1]] else "9_dec2024"

truth_file <- function(population) {
  file.path(repo_root, "data", "population_coefficients",
            paste0("population_coef_", population, "_hhsub.csv"))
}

presets <- list()

# Scenario 9: two-step pseudo-likelihood on the complete-case samples of each population.
for (population in c("dec2024", "oct2026")) {
  presets[[paste0("9_", population)]] <- list(
    coef_dir = file.path(repo_root, "data", "estimates", population, "scen9"),
    file_pattern = "^scen9_coef_[0-9]+\\.csv$",
    sim_id_pattern = "^scen9_coef_([0-9]+)\\.csv$",
    out_file = file.path(repo_root, "data", "summaries", population, "coverage_scen9.csv"),
    true_file = truth_file(population),
    label = paste0("scen9_", population)
  )
}

# Scenarios 15 (s1) and 16 (s2): MI with scaled weights as imputation covariates, variants v1-v3,
# on the latest (Oct2026) samples.
for (scen_i in c("15", "16")) {
  for (var_i in c("v1", "v2", "v3")) {
    presets[[paste0(scen_i, "_", var_i, "_oct2026")]] <- list(
      coef_dir = file.path(repo_root, "data", "estimates", "oct2026", paste0("scen", scen_i, "_mi"), var_i),
      file_pattern = paste0("^scen", scen_i, "_mi_coef_[0-9]+\\.csv$"),
      sim_id_pattern = paste0("^scen", scen_i, "_mi_coef_([0-9]+)\\.csv$"),
      out_file = file.path(repo_root, "data", "summaries", "oct2026", paste0("coverage_scen", scen_i, "_mi_", var_i, ".csv")),
      true_file = truth_file("oct2026"),
      label = paste0("scen", scen_i, "_mi_", var_i, "_oct2026")
    )
  }
}

preset <- presets[[scen_num]]
if (is.null(preset)) stop(sprintf("Unknown scenario preset '%s'", scen_num))
coef_dir <- preset$coef_dir

read_one <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  sim_id <- sub(preset$sim_id_pattern, "\\1", basename(path))
  df$sim_id <- as.integer(sim_id)

  keep_cols <- c("sim_id", "term", "estimate", "std_error", "conf_low", "conf_high")
  for (col in keep_cols) {
    if (!col %in% names(df)) df[[col]] <- NA_real_
  }
  df[, keep_cols, drop = FALSE]
}

true_coefficients <- utils::read.csv(preset$true_file, stringsAsFactors = FALSE, check.names = FALSE)
names(true_coefficients) <- c("term", "true")
true_coefficients$true <- as.numeric(true_coefficients$true)

coefficient_files <- list.files(coef_dir, pattern = preset$file_pattern, full.names = TRUE)

if (length(coefficient_files) == 0) {
  stop(sprintf("No estimate files found matching pattern '%s' in %s", preset$file_pattern, coef_dir))
}

all_estimates <- do.call(rbind, lapply(coefficient_files, read_one))
all_estimates$estimate <- as.numeric(all_estimates$estimate)
all_estimates$std_error <- as.numeric(all_estimates$std_error)
all_estimates$conf_low <- as.numeric(all_estimates$conf_low)
all_estimates$conf_high <- as.numeric(all_estimates$conf_high)

# Keep terms without a known true value (bias/coverage are NA for those)
analysis_data <- merge(all_estimates, true_coefficients, by = "term", all.x = TRUE)

z975 <- stats::qnorm(0.975)
unique_terms <- sort(unique(analysis_data$term))

term_summaries <- lapply(unique_terms, function(term_i) {
  sub_df <- analysis_data[analysis_data$term == term_i, , drop = FALSE]
  true_val <- sub_df$true[which(!is.na(sub_df$true))[1]]
  mean_est <- mean(sub_df$estimate, na.rm = TRUE)
  emp_bias <- mean(sub_df$estimate - true_val, na.rm = TRUE)
  emp_se <- stats::sd(sub_df$estimate, na.rm = TRUE)
  model_se <- mean(sub_df$std_error, na.rm = TRUE)

  ci_low <- sub_df$conf_low
  ci_high <- sub_df$conf_high
  covered <- ci_low <= true_val & ci_high >= true_val
  reject_h0 <- ci_low > 0 | ci_high < 0

  emp_ci_low <- sub_df$estimate - z975 * emp_se
  emp_ci_high <- sub_df$estimate + z975 * emp_se
  covered_emp <- emp_ci_low <= true_val & emp_ci_high >= true_val

  is_variance_term <- grepl("^(var_|s2$|sd_|lns_|lnsig|sigma|Var\\()", term_i, ignore.case = FALSE)

  data.frame(
    scenario = preset$label,
    term_type = if (is_variance_term) "variance_component" else "fixed_effect",
    term = term_i,
    n_reps = sum(!is.na(sub_df$estimate)),
    true_value = true_val,
    mean_estimate = mean_est,
    empirical_bias = emp_bias,
    relative_bias = if (!is.na(true_val) && abs(true_val) > 0) emp_bias / true_val else NA_real_,
    empirical_se = emp_se,
    model_se = model_se,
    relative_se_diff = if (!is.na(emp_se) && emp_se > 0) {
      (model_se - emp_se) / emp_se
    } else {
      NA_real_
    },
    coverage_rate = mean(covered, na.rm = TRUE),
    coverage_emp_se = mean(covered_emp, na.rm = TRUE),
    p_reject_h0 = mean(reject_h0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

final_summary <- do.call(rbind, term_summaries)

# Sort: fixed effects first, variance components at the bottom
type_order <- c("fixed_effect" = 1L, "variance_component" = 2L)
final_summary <- final_summary[order(type_order[final_summary$term_type], final_summary$term), ]
row.names(final_summary) <- NULL

dir.create(dirname(preset$out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(final_summary, preset$out_file, row.names = FALSE)
message(sprintf("Successfully saved coverage assessment for %s to %s", preset$label, preset$out_file))
