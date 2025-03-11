#' Classify CpG Sites
#' Determines if a CpG site follows a linear, non-linear, or non-correlated trajectory.
#' @param beta_values A numeric vector of methylation beta values for a CpG site.
#' @param age A numeric vector of ages corresponding to the samples.
#' @param ages_grid A numeric vector of ages for GAM predictions.
#' @param cpg_name A value for the name of the variable.
#' @return A list with classification results and GAM predictions.
#' @export
#' @import mgcv
#' @import lmtest
#' @import minerva
classify_cpg <- function(beta_values, age, ages_grid, cpg_name) {
  beta_values <- as.numeric(beta_values)
  lm_model <- lm(beta_values ~ age)
  lm_pval <- summary(lm_model)$coefficients[2, 4]
  lm_coef <- summary(lm_model)$coefficients[2, 1]
  bp_pval <- bptest(lm_model, ~ age)$p.value
  white_pval <- bptest(lm_model, ~ age + I(age^2))$p.value
  # Step 1: Compute MIC for correlation assessment
  mic_value <- minerva::mine(age, beta_values)$MIC  # Compute MIC

  # Step 2: Check for non-correlation based on MIC
  if (mic_value < 0.3) {  # Non-correlated if MIC is below threshold
    return(list(CpG = cpg_name, classification = "NC", MIC = mic_value, lm_pval = lm_pval, lm_coef = lm_coef,
                dbic_lg = NA, bp_pval = bp_pval, white_pval = white_pval, gam_predictions = NA))
  }

  gam_model <- gam(beta_values ~ s(age), method = "REML")
  bic_lm <- BIC(lm_model)
  bic_gam <- BIC(gam_model)
  dbic_lg <- bic_lm - bic_gam


  if (dbic_lg > 2) {
    classification <- "NL"
    gam_predictions <- predict.gam(gam_model, data.frame(Age = ages_grid))
  } else if (lm_pval > 0.01) {
    classification <- "VI"
    gam_predictions <- NA
  } else {
    classification <- ifelse(lm_coef > 0, "LI", "LD")
    gam_predictions <- NA
  }


  return(list(CpG = cpg_name, classification = classification, MIC = mic_value, lm_pval = lm_pval, lm_coef = lm_coef,
              dbic_lg = dbic_lg, bp_pval = bp_pval, white_pval = white_pval, gam_predictions = gam_predictions))
}
