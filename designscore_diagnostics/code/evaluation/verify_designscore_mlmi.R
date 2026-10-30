# Independent base-R recalculation of every pooled estimate and summary.
# Does not fit models or change raw estimates.
args <- commandArgs(TRUE)
root <- normalizePath(args[1], winslash="/", mustWork=TRUE)
label <- if (length(args)>1) args[2] else "mlmi_joint6_m20_run2"
first <- if (length(args)>2) as.integer(args[3]) else 1L
last <- if (length(args)>3) as.integer(args[4]) else 100L
M <- if (length(args)>4) as.integer(args[5]) else 20L
data_dir <- file.path(root,"data/estimates/oct2026",label)
dest <- file.path(root,"summaries/oct2026",paste0(label,"_",first,"_",last))
truth <- read.csv(file.path(root,"data/population_coefficients/population_coef_oct2026_hhsub.csv"))
close <- function(a,b) {
  stopifnot(length(a)==length(b), identical(is.na(a),is.na(b)),
            identical(is.infinite(a),is.infinite(b)))
  keep <- is.finite(a) & is.finite(b)
  stopifnot(all(abs(a[keep]-b[keep])<=1e-9*pmax(1,abs(a[keep]),abs(b[keep]))))
}
all_pooled <- list()
source_map_path <- file.path(dest,"replicate_source_runs.csv")
source_map <- if(file.exists(source_map_path)) read.csv(source_map_path) else NULL
for (sim in first:last) {
  actual_dir <- data_dir
  if(!is.null(source_map)) {
    selected <- source_map$source_run[source_map$sim_id==sim]
    stopifnot(length(selected)==1L)
    actual_dir <- file.path(root,"data/estimates/oct2026",selected)
  }
  raw <- read.csv(file.path(actual_dir,"per_imputation",paste0("mi_estimates_",sim,".csv")))
  pooled <- read.csv(file.path(actual_dir,"pooled",paste0("mi_coef_",sim,".csv")))
  audit <- read.csv(file.path(actual_dir,"audits",paste0("mi_audit_",sim,".csv")))
  stopifnot(nrow(raw)==20L*M, nrow(pooled)==20L, nrow(audit)==M,
            all(raw$inference_usable==1), all(pooled$fit_status=="OK"),
            all(audit$boundary_mode==0), all(audit$se_status==0),
            all(audit$original_converged==1), length(unique(audit$input_signature))==M)
  for (m in seq_len(M)) {
    signatures <- unique(c(raw$input_signature[raw$imputation==m], audit$input_signature[audit$imputation==m]))
    stopifnot(length(signatures)==1L, nzchar(signatures))
  }
  for (i in seq_len(nrow(pooled))) {
    p <- pooled[i,]
    d <- raw[raw$method==p$method & raw$term==p$term,]
    stopifnot(nrow(d)==M, identical(sort(d$imputation),seq_len(M)))
    logarithmic <- p$term %in% c("var_hhid","var_subid","s2")
    q <- if(logarithmic) log(d$estimate) else d$estimate
    u <- if(logarithmic) (d$std_error/d$estimate)^2 else d$std_error^2
    qbar <- mean(q); within <- mean(u); between <- var(q)
    total <- within+(1+1/M)*between
    df <- if(between>0) (M-1)*(1+within/((1+1/M)*between))^2 else Inf
    ci <- qbar+c(-1,1)*qt(.975,df)*sqrt(total)
    est <- if(logarithmic) exp(qbar) else qbar
    if(logarithmic) ci <- exp(ci)
    se <- sqrt(total)*if(logarithmic) est else 1
    close(c(p$estimate,p$std_error,p$conf_low,p$conf_high,p$within_variance,p$between_variance,p$total_variance),
          c(est,se,ci,within,between,total))
    if(is.finite(df)) close(p$mi_df,df)
  }
  all_pooled[[as.character(sim)]] <- pooled
}
all_pooled <- do.call(rbind,all_pooled)
aggregate <- read.csv(file.path(dest,"pooled_estimates.csv"))
key <- function(d) paste(d$sim_id,d$method,d$term,sep="/")
stopifnot(!anyDuplicated(key(aggregate)),setequal(key(aggregate),key(all_pooled)))
aggregate <- aggregate[match(key(all_pooled),key(aggregate)),]
for(k in c("estimate","std_error","conf_low","conf_high","within_variance","between_variance","total_variance")) close(aggregate[[k]],all_pooled[[k]])
for(method in c("original","designscore")) {
  summary <- read.csv(file.path(dest,paste0("coverage_scen3_",method,"_mi.csv")))
  for(i in seq_len(nrow(summary))) {
    s <- summary[i,]; d <- all_pooled[all_pooled$method==method & all_pooled$term==s$term,]
    target <- truth$estimate[truth$term==s$term]
    covered <- d$conf_low<=target & target<=d$conf_high
    stopifnot(s$n_reps==last-first+1L, nrow(d)==s$n_reps)
    close(c(s$mean_estimate,s$empirical_bias,s$empirical_se,s$model_se,s$coverage_rate,s$bias_mcse),
          c(mean(d$estimate),mean(d$estimate)-target,sd(d$estimate),mean(d$std_error),mean(covered),sd(d$estimate)/sqrt(nrow(d))))
  }
}
paired <- read.csv(file.path(dest,"paired_comparison.csv"))
for(i in seq_len(nrow(paired))) {
  p <- paired[i,]; d <- all_pooled[all_pooled$term==p$term,]
  a <- d[d$method=="original",]; b <- d[d$method=="designscore",]
  a <- a[order(a$sim_id),]; b <- b[order(b$sim_id),]
  stopifnot(identical(a$sim_id,b$sim_id))
  delta <- b$estimate-a$estimate
  target <- truth$estimate[truth$term==p$term]
  cd <- as.integer(b$conf_low<=target & target<=b$conf_high)-as.integer(a$conf_low<=target & target<=a$conf_high)
  close(c(p$mean_designscore_minus_original,p$paired_difference_mcse,p$coverage_difference,p$coverage_difference_mcse),
        c(mean(delta),sd(delta)/sqrt(length(delta)),mean(cd),sd(cd)/sqrt(length(cd))))
}
message("PASS: independently verified ",last-first+1L," replicates, ",nrow(all_pooled),
        " pooled rows, both coverage summaries, and paired comparisons.")
