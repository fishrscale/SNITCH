#' Prepare DNA Methylation Data
#' Transforms and scales DNA methylation data for analysis.
#' @param data A matrix of DNA methylation beta values (CpGs in columns, samples in rows).
#' @param age A numeric vector representing the ages of individuals.
#' @return A list containing the scaled data, the age vector, and an age grid.
#' @export
prepare_data <- function(data, age) {
  if (!is.matrix(data)) stop("Data must be a matrix")
  if (length(age) != nrow(data)) stop("Age vector length must match number of rows in data")

  scaled_data <- scale(data)  # Normalize
  ages_grid <- seq(min(age), max(age))  # Smooth age grid

  return(list(dat_cor = scaled_data, Age = age, ages_grid = ages_grid))
}
