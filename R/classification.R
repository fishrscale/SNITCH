#' Classify CpG Sites
#' Determines if a CpG site follows a linear, non-linear, or non-correlated trajectory.
#' @param beta_values A numeric vector of methylation beta values for a CpG site.
#' @param age A numeric vector of ages corresponding to the samples.
#' @param ages_grid A numeric vector of ages for GAM predictions.
#' @param cpg_name A value for the name of the variable.
#' @param covariates Optional data.frame of covariates (samples as rows).
#' @return A list with classification results and GAM predictions.
#' @export
#' @import mgcv
#' @import lmtest
#' @importFrom stats BIC as.formula lm median predict rbeta reformulate rnorm
classify_cpg <- function(beta_values, age, ages_grid, cpg_name, covariates = NULL) {
  beta_values <- as.numeric(beta_values)

  # Combine data for modeling
  data <- data.frame(beta = beta_values, Age = age, Age2 = age^2)
  if (!is.null(covariates)) {
    data <- cbind(data, covariates)
  }

  # Build formulas dynamically
  covar_terms <- if (!is.null(covariates)) colnames(covariates) else NULL
  lm_formula <- reformulate(c("Age", covar_terms), response = "beta")
  gam_formula <- reformulate(c("s(Age, k=5)", covar_terms), response = "beta")

  lm_model <- lm(lm_formula, data = data)
  lm_pval <- summary(lm_model)$coefficients[2, 4]
  lm_coef <- summary(lm_model)$coefficients[2, 1]
  bp_pval <- bptest(lm_model, ~ Age, data = data)$p.value
  white_pval <- bptest(lm_model, ~ Age + Age2, data = data)$p.value

  gam_model <- gam(as.formula(gam_formula), data = data, method = "REML", select = FALSE)
  bic_lm <- BIC(lm_model)
  bic_gam <- BIC(gam_model)
  dbic_lg <- bic_lm - bic_gam
  gam_edf <- summary(gam_model)$s.table[1]
  gam_pval <- summary(gam_model)$s.table[4]
  gam_refdf <- summary(gam_model)$s.table[2]


  if (dbic_lg > 2) {
    # Build prediction data frame
    newdata_pred <- data.frame(Age = ages_grid)

    if (!is.null(covariates)) {
      for (covar in colnames(covariates)) {
        if (is.factor(covariates[[covar]])) {
          # Use the first level as the reference
          newdata_pred[[covar]] <- factor(rep(levels(covariates[[covar]])[1], length(ages_grid)),
                                          levels = levels(covariates[[covar]]))
        } else {
          # Use the median (or mean) for numeric covariates
          newdata_pred[[covar]] <- rep(median(covariates[[covar]], na.rm = TRUE), length(ages_grid))
        }
      }
    }

    # Predict using full model formula
    gam_predictions <- predict(gam_model, newdata = newdata_pred)

  } else if (lm_pval > 0.01) {
    gam_predictions <- NA
  } else {
    gam_predictions <- NA
  }


  return(list(CpG = cpg_name, lm_pval = lm_pval, lm_coef = lm_coef,
              dbic_lg = dbic_lg, bp_pval = bp_pval, white_pval = white_pval,
              gam_predictions = gam_predictions, gam_edf = gam_edf,
              gam_refdf = gam_refdf,
              gam_pval = gam_pval))
}
