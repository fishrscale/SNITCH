#' Run Parallel CpG Classification
#' Classifies multiple CpGs using parallel processing.
#' @param dat_scaled A matrix of scaled beta values (CpGs as columns, samples as rows).
#' @param age A numeric vector of ages.
#' @param ages_grid A numeric vector of ages for GAM predictions.
#' @param covariates Optional data.frame of covariates (samples as rows).
#' @param cores Number of cores to use for parallel processing (default: detectCores() - 2).
#' @return A data frame with classification results for all CpGs.
#' @export
#' @import parallel
#' @importFrom stats p.adjust
run_parallel_classification <- function(dat_scaled, age, ages_grid, covariates = NULL, cores = NULL) {
  if (is.null(cores)) {
    cores <- detectCores() - 2
  }
  cl <- makeCluster(cores)

  clusterEvalQ(cl, {
    library(mgcv)
    library(lmtest)
    library(minerva)
  })
  # Export variables anfd function to workers
  # Export function arguments directly using their names in the function
  clusterExport(cl, varlist = c("dat_scaled", "age", "ages_grid", "covariates", "classify_cpg"), envir = environment())

  results <- parLapply(cl, 1:ncol(dat_scaled), function(i) {
    classify_cpg(dat_scaled[, i], age, ages_grid, colnames(dat_scaled)[i], covariates)
  })

  stopCluster(cl)

  # Convert results to a structured data frame
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(CpG = x$CpG, lm_pval = x$lm_pval,
               lm_coef = x$lm_coef, dbic_lg = x$dbic_lg, bp_pval = x$bp_pval,
               white_pval = x$white_pval,
               Predictions = I(list(x$gam_predictions)),
               gam_edf = x$gam_edf,
               gam_refdf = x$gam_refdf,
               gam_pval = x$gam_pval)
  }))

  results_df$adj_lm <- p.adjust(results_df$lm_pval)
  results_df$adj_white <- p.adjust(results_df$white_pval)
  results_df$classification <- NA
  results_df$adj_gam <- p.adjust(results_df$gam_pval)

  results_df[results_df$adj_lm <= 0.01 & results_df$lm_coef > 0, 'classification'] <- 'LI'
  results_df[results_df$adj_lm <= 0.01 & results_df$lm_coef < 0, 'classification'] <- 'LD'
  results_df[results_df$adj_lm > 0.01 & results_df$adj_white <= 0.01, 'classification'] <- 'VI'
  results_df[results_df$adj_lm > 0.01 & results_df$adj_white > 0.01, 'classification'] <- 'NC'
  results_df[results_df$dbic_lg > 2 & results_df$adj_gam <= 0.01, 'classification'] <- 'NL'

  return(results_df)
}
