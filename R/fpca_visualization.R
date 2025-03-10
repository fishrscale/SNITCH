#' Generate FPCA Diagnostic Plots
#'
#' Creates and saves plots of FPCA results including the mean function, variance explained, and functional PCs.
#'
#' @param fpca_res A list containing FPCA results as returned by `perform_fpca()`.
#' @param ages_grid Numeric vector of age grid values.
#' @param output_dir Directory where plots should be saved.
#' @export
#' @import ggplot2
#'
plot_fpca_results <- function(fpca_res, ages_grid, output_dir = "./Results") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Mean function plot
  mean_plot <- data.frame(Age = ages_grid, Mean = fpca_res$mu)
  ggplot(mean_plot, aes(x = Age, y = Mean)) +
    geom_line(color = "blue") +
    theme_minimal() +
    labs(title = "Mean Methylation Trajectory", x = "Age", y = "Mean Methylation") +
    ggsave(filename = file.path(output_dir, "SNITCH_fpca_mu.pdf"))

  # Variance explained plot
  variance_df <- data.frame(PC = 1:length(fpca_res$variance_explained),
                            Variance_Explained = fpca_res$variance_explained)
  ggplot(variance_df, aes(x = PC, y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "blue") +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Variance Explained by Principal Components",
         x = "Principal Component", y = "Percentage of Variance Explained") +
    ggsave(filename = file.path(output_dir, "SNITCH_fpca_var.pdf"))
}
