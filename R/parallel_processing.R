#' Run Parallel CpG Classification
#' Classifies multiple CpGs using parallel processing.
#' @param dat_cor A matrix of scaled beta values (CpGs as columns, samples as rows).
#' @param age A numeric vector of ages.
#' @param ages_grid A numeric vector of ages for GAM predictions.
#' @param cores Number of cores to use for parallel processing (default: detectCores() - 2).
#' @return A data frame with classification results for all CpGs.
#' @export
#' @import parallel
run_parallel_classification <- function(dat_cor, age, ages_grid, cores = NULL) {
  if (is.null(cores)) {
    cores <- detectCores() - 2
  }
  cl <- makeCluster(cores)

  clusterEvalQ(cl, {
    library(mgcv)
    library(lmtest)
  })

  results <- parLapply(cl, 1:ncol(dat_cor), function(i) {
    classify_cpg(dat_cor[, i], age, ages_grid)
  })

  stopCluster(cl)

  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(classification = x$classification, lm_pval = x$lm_pval,
               lm_coef = x$lm_coef, dbic_lg = x$dbic_lg, bp_pval = x$bp_pval,
               white_pval = x$white_pval)
  }))

  return(results_df)
}
