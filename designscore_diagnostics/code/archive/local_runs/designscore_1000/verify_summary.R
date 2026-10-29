directory <- "data/estimates/oct2026/scen3_designscore/exact"
files <- file.path(directory, sprintf("scen3_designscore_coef_%d.csv", 1:1000))
available <- file.exists(files)
estimates <- do.call(rbind, lapply(files[available], read.csv, check.names = FALSE))
summary <- read.csv("data/summaries/oct2026/coverage_scen3_designscore_exact.csv")
stopifnot(nrow(summary) == 10L, all(summary$n_expected == 1000),
          all(summary$first_replicate == 1), all(summary$last_replicate == 1000))
for (i in seq_len(nrow(summary))) {
  row <- summary[i, ]
  values <- estimates[estimates$term == row$term, ]
  stopifnot(!anyDuplicated(values$sim_id), all(values$sim_id %in% 1:1000))
  finite_points <- is.finite(values$estimate)
  good <- values$inference_usable == 1 & is.finite(values$std_error) &
          is.finite(values$conf_low) & is.finite(values$conf_high)
  points <- values$estimate[finite_points]
  covered <- values$conf_low[good] <= row$true_value & values$conf_high[good] >= row$true_value
  expected <- c(mean(points), sd(points), mean(values$std_error[good]), mean(covered),
                sum(covered)/1000)
  actual <- unlist(row[c("mean_estimate", "empirical_se", "model_se", "coverage_rate",
                        "coverage_over_all_requested")])
  stopifnot(isTRUE(all.equal(unname(expected), unname(actual), tolerance=1e-10)),
            row$n_reps == length(points), row$n_valid_intervals == sum(good),
            row$n_missing_points == 1000-length(points),
            row$n_missing_intervals == 1000-sum(good))
}
cat(sprintf("Independent R check passed: %d coefficient rows from %d requested replicate files.\n",
            nrow(estimates), sum(available)))
