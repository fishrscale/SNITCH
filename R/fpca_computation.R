#' Perform FPCA on DNA Methylation Data
#'
#' Runs Functional Principal Component Analysis (FPCA) on smoothed methylation values.
#'
#' @param nl_var_smooth A matrix of smoothed DNA methylation beta values (CpGs as rows, samples as columns).
#' @param ages_grid Numeric vector of age grid values.
#' @param pve Proportion of variance explained threshold for selecting PCs (default = 0.9999).
#' @return A list containing FPCA results including scores, mean function, eigenfunctions, and variance explained.
#' @export
#' @import refund
perform_fpca <- function(nl_var_smooth, ages_grid, pve = 0.9999) {

  n_timepoints <- ncol(nl_var_smooth)
  k <- min(35, floor(0.8 * n_timepoints))

  fpca_result <- refund::fpca.face(Y = nl_var_smooth, argvals = ages_grid, pve = pve, knots = k)

  return(list(
    scores = fpca_result$scores,
    mu = fpca_result$mu,
    efunctions = fpca_result$efunctions,
    evalues = fpca_result$evalues,
    variance_explained = fpca_result$evalues / sum(fpca_result$evalues) * 100
  ))
}
