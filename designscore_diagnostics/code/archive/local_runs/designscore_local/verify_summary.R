directory <- "data/estimates/oct2026/scen3_designscore/exact"
files <- file.path(directory, sprintf("scen3_designscore_coef_%d.csv", 1:200))
stopifnot(all(file.exists(files)))
estimates <- do.call(rbind, lapply(files, read.csv, check.names = FALSE))
summary <- read.csv("data/summaries/oct2026/coverage_scen3_designscore_exact.csv")
stopifnot(nrow(estimates) == 2000L, nrow(summary) == 10L,
          all(summary$n_reps == 200), all(summary$n_valid_intervals == 200),
          all(summary$n_missing_points == 0), all(summary$n_missing_intervals == 0))
for (i in seq_len(nrow(summary))) {
  row <- summary[i, ]
  values <- estimates[estimates$term == row$term, ]
  stopifnot(identical(sort(values$sim_id), 1:200), all(values$se_status == 0),
            all(values$inference_usable == 1), all(values$conf_low < values$conf_high))
  expected <- c(mean(values$estimate), sd(values$estimate), mean(values$std_error),
    mean(values$conf_low <= row$true_value & values$conf_high >= row$true_value))
  actual <- unlist(row[c("mean_estimate", "empirical_se", "model_se", "coverage_rate")])
  stopifnot(max(abs(expected-actual)) < 1e-10)
}
cat("Independent R cross-check passed: all 2000 estimates, means, empirical/model SEs, and coverage rates.\n")
