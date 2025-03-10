#' Classify CpG Sites
#' Determines if a CpG site follows a linear, non-linear, or non-correlated trajectory.
#' @param beta_values A numeric vector of methylation beta values for a CpG site.
#' @param age A numeric vector of ages corresponding to the samples.
#' @param ages_grid A numeric vector of ages for GAM predictions.
#' @return A list with classification results and GAM predictions.
#' @export
#' @import mgcv
#' @import lmtest
classify_cpg <- function(beta_values, age, ages_grid) {
  lm_model <- lm(beta_values ~ age)
  lm_pval <- summary(lm_model)$coefficients[2, 4]
  lm_coef <- summary(lm_model)$coefficients[2, 1]

  if (lm_pval > 0.01) {
    bp_pval <- bptest(lm_model, ~ age)$p.value
    white_pval <- bptest(lm_model, ~ age + I(age^2))$p.value
    return(list(classification = "NC", lm_pval = lm_pval, lm_coef = NA,
                dbic_lg = NA, bp_pval = bp_pval, white_pval = white_pval, gam_predictions = NA))
  }

  gam_model <- gam(beta_values ~ s(age), method = "REML")
  bic_lm <- BIC(lm_model)
  bic_gam <- BIC(gam_model)
  dbic_lg <- bic_lm - bic_gam

  if (dbic_lg > 2) {
    classification <- "NL"
    gam_predictions <- predict(gam_model, newdata = data.frame(age = ages_grid))
  } else {
    classification <- ifelse(lm_coef > 0, "LI", "LD")
    gam_predictions <- NA
  }

  return(list(classification = classification, lm_pval = lm_pval, lm_coef = lm_coef,
              dbic_lg = dbic_lg, bp_pval = NA, white_pval = NA, gam_predictions = gam_predictions))
}
