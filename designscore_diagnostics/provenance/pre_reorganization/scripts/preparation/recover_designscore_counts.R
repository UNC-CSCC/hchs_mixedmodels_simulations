#!/usr/bin/env Rscript
# Replay ONLY the sampling stage, using saved seeds; never overwrite samples.
# All subject IDs and sample/population checksums must match before publishing.
# Usage: Rscript recover_designscore_counts.R <HCHS-root> <counts-output.csv> [last-replicate] [first-replicate]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 2:4)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
output <- args[2]
if (file.exists(output)) stop("Refusing to overwrite counts: ", output)
generator <- file.path(root, "code/v3/Oct2026/sample/generate_samples.R")
env <- new.env(parent = globalenv())
source(generator, local = env)
config <- make_sample_config(root)
population_file <- config$paths$population_file
population_md5 <- toupper(unname(tools::md5sum(population_file)))
stopifnot(identical(population_md5, config$population$expected_md5))

read_design <- function(file, columns) {
  header <- names(read.csv(file, nrows = 0L, check.names = FALSE))
  if (!all(columns %in% header)) stop("Missing design columns: ", file)
  classes <- ifelse(header %in% columns, NA_character_, "NULL")
  read.csv(file, colClasses = classes, check.names = FALSE)
}
cache <- file.path(dirname(output), paste0("baseline_design_", population_md5, ".rds"))
if (file.exists(cache)) {
  population <- readRDS(cache)
} else {
  message("Reading baseline population design columns")
  population <- read_design(population_file,
    c("strat", "BGid", "hhid", "subid", "v_num", "hisp_strat", "age_strat"))
  # Match the production reader's baseline subject order.
  population <- population[population$v_num == 1L, , drop = FALSE]
  population <- population[order(population$subid), , drop = FALSE]
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(population, cache)
}
stopifnot(!anyDuplicated(population$subid), !anyNA(population))

# Cache the deterministic candidate lists instead of rescanning every household
# for each selected BG. Draw order and sample_values() match the production
# sampler exactly. No weighting, outcome or missingness calculations are needed.
design <- config$sampling
candidates <- lapply(1:8, function(s) sort(unique(population$BGid[population$strat == s])))
households <- population[!duplicated(population$hhid), c("BGid", "hhid", "hisp_strat")]
hh_hisp <- split(households$hhid[households$hisp_strat], households$BGid[households$hisp_strat])
hh_other <- split(households$hhid[!households$hisp_strat], households$BGid[!households$hisp_strat])
sample_values <- env$sample_values
replay_selection <- function(seed) {
  set.seed(seed)
  selected_bgs <- lapply(1:8, function(s)
    sample_values(candidates[[s]], round(design$bg_prob[s] * design$bg_count[s])))
  selected_households <- list()
  for (s in 1:8) for (bg in selected_bgs[[s]]) {
    hisp <- hh_hisp[[as.character(bg)]]
    other <- hh_other[[as.character(bg)]]
    selected_households[[length(selected_households) + 1L]] <- c(
      sample_values(hisp, round(design$hh_prob_hispanic[s] * length(hisp))),
      sample_values(other, round(design$hh_prob_other[s] * length(other))))
  }
  eligible <- population$hhid %in% unique(unlist(selected_households, use.names = FALSE))
  younger <- population$subid[eligible & !population$age_strat]
  older <- population$subid[eligible & population$age_strat]
  ny <- round(design$subject_prob[["younger"]] * length(younger))
  no <- round(design$subject_prob[["older"]] * length(older))
  ids <- c(sample_values(younger, ny), sample_values(older, no))
  list(ids = ids, eligible_younger = length(younger), selected_younger = ny,
       eligible_older = length(older), selected_older = no)
}

directory <- file.path(root, "data/Oct2026/sample_miss00pct")
files <- list.files(directory, "^samplemiss_[0-9]+\\.csv$", full.names = TRUE)
ids <- as.integer(sub("^samplemiss_([0-9]+)\\.csv$", "\\1", basename(files)))
files <- files[order(ids)]; ids <- sort(ids)
if (length(args) >= 3L) {
  first <- if (length(args) == 4L) as.integer(args[4]) else 1L
  last <- as.integer(args[3])
  stopifnot(!is.na(first), !is.na(last), first >= 1L, last >= first)
  selected <- ids >= first & ids <= last
  files <- files[selected]; ids <- ids[selected]
}
if (!length(ids) || anyDuplicated(ids)) stop("Invalid replicate inventory")
results <- vector("list", length(files))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
for (i in seq_along(files)) {
  metadata <- env$sample_read_completion_record(files[i])
  if (is.null(metadata)) stop("Missing completion metadata: ", files[i])
  stopifnot(metadata[["population_md5"]] == population_md5,
            metadata[["population_version"]] == "Oct2026",
            metadata[["replicate"]] == as.character(ids[i]),
            metadata[["missing_percent"]] == "0",
            toupper(unname(tools::md5sum(files[i]))) == metadata[["csv_md5"]])
  seed <- as.integer(metadata[["file_seed"]])
  if (is.na(seed)) stop("No reproducible seed: ", files[i])
  captured <- replay_selection(seed)
  # Independently compare the cached replay with the production function on
  # the first replicate of every invocation, then check all saved sample IDs.
  if (i == 1L) {
    set.seed(seed)
    reference <- env$draw_three_stage_sample(population, config)
    stopifnot(identical(sort(as.integer(reference$subid)), sort(as.integer(captured$ids))))
  }
  sample <- read_design(files[i], c("subid", "v_num", "age_strat"))
  stopifnot(nrow(sample) == as.integer(metadata[["rows"]]))
  observed <- sample[sample$v_num == 1L, , drop = FALSE]
  if (!identical(sort(as.integer(captured$ids)), sort(as.integer(observed$subid)))) {
    stop(sprintf("Replicate %d: replay has %d IDs, saved sample has %d, overlap %d",
      ids[i], length(captured$ids), nrow(observed), length(intersect(captured$ids, observed$subid))))
  }
  stopifnot(!anyDuplicated(observed$subid),
            captured$selected_younger == sum(!observed$age_strat),
            captured$selected_older == sum(observed$age_strat),
            length(captured) == 5L)
  results[[i]] <- data.frame(population = "oct2026", sim_id = ids[i],
    eligible_younger = captured$eligible_younger,
    selected_younger = captured$selected_younger,
    eligible_older = captured$eligible_older,
    selected_older = captured$selected_older,
    file_seed = seed, sample_md5 = metadata[["csv_md5"]],
    population_md5 = population_md5,
    generator_md5 = unname(tools::md5sum(generator)),
    ids_verified = TRUE, r_version = as.character(getRversion()))
  message(sprintf("Verified exact selection counts: %d/%d (replicate %d)",
                  i, length(files), ids[i]))
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, results), output, row.names = FALSE)
message("Saved verified counts to ", output)
