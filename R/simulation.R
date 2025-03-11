#' Generate Simulated Methylation Data
#'
#' Creates a dataset simulating DNA methylation intensity across different ageing-related patterns.
#'
#' @param n_people Number of individuals to simulate trajectories for.
#' @param seed Random seed for reproducibility.
#' @param plot TRUE/FALSE. Visualise the representative functions and the simulated data (default = FALSE).
#' @return A list containing simulated methylation intensities, age data, and group labels.
#' @export
#' @import dplyr
#' @import tibble
#' @import ggplot2
simulate_methylation_data <- function(n_people = 300, seed = 123, plot = FALSE, output_dir = "./") {
  set.seed(seed)

  ages <- sample(1:100, n_people, replace = TRUE)

  # Define methylation intensity functions
  # Define methylation intensity functions
  functions_list <- list(
    "Non-Correlated" = function(age) rep(0.5, length(age)),
    "Linear Increasing" = function(age) age / max(age),
    "Linear Decreasing" = function(age) 1 - (age / max(age)),
    "Quadratic Increasing" = function(age) (age / max(age))^2,
    "Quadratic Decreasing" = function(age) 1 - (age / max(age))^2,
    "Logarithmic Increasing" = function(age) log(age + 1) / log(max(age) + 1),
    "Logarithmic Decreasing" = function(age) 1 - log(age + 1) / log(max(age) + 1),
    "Sigmoid Increasing 25" = function(age) 1 / (1 + exp(-0.1 * (age - 50))),
    "Sigmoid Decreasing 25" = function(age) 1 - 1 / (1 + exp(-0.1 * (age - 50))),
    "Sigmoid Increasing 40" = function(age) 1 / (1 + exp(-0.1 * (age - 80))),
    "Sigmoid Decreasing 40" = function(age) 1 - 1 / (1 + exp(-0.1 * (age - 80))),
    "Sigmoid Increasing 80" = function(age) 1 / (1 + exp(-0.1 * (age - 100))),
    "Sigmoid Decreasing 80" = function(age) 1 - 1 / (1 + exp(-0.1 * (age - 100))),
    "Variance Increasing 25" = function(age) sapply(age, function(a) {
      sd <- ifelse(a < 25, 0.01, 0.01 + (0.3 - 0.01) * ((a - 25) / (max(age) - 25)))
      pmax(pmin(rnorm(1, mean = 0.5, sd = sd), 1), 0)
    }),
    "Non-Monotonic" = function(age) {
      peak_center <- max(age) / 2
      spread <- max(age) / 4
      intensity <- exp(-((age - peak_center)^2) / (2 * spread^2))
      pmax(pmin(intensity, 1), 0)
    }
  )

  groups <- rep(names(functions_list), each = 200)

  simulate_beta <- function(mu, b = 1) {
    epsilon <- 1e-6
    mu <- pmin(pmax(mu, epsilon), 1 - epsilon)
    rbeta(length(mu), mu * b / (1 - mu), b)
  }

  simulated_data <- data.frame(Person = integer(), Site = character(), Age = integer(), Intensity = numeric(), Group = character())
  site_id <- 1

  for (group in unique(groups)) {
    func <- functions_list[[group]]
    for (i in 1:200) {
      mu <- func(ages)
      intensity <- simulate_beta(mu)

      simulated_data <- rbind(simulated_data, data.frame(
        Person = 1:n_people, Site = paste0('Var', site_id), Age = ages, Intensity = intensity, Group = group
      ))
      site_id <- site_id + 1
    }
  }


  if (plot) {
    functions_data <- data.frame(
      Age = rep(ages, length(functions_list)),
      Intensity = unlist(lapply(functions_list, function(f) f(ages))),
      Group = rep(names(functions_list), each = length(ages))
    )

    # Plot representative functions
    ggsave(
      filename = file.path(output_dir, "functions_sim_data.pdf"),
      plot = ggplot(functions_data, aes(x = Age, y = Intensity, color = Group)) +
        geom_line(size = 1) +
        facet_wrap(~ Group, scales = "fixed") +
        labs(
          title = "Methylation Intensity Patterns",
          x = "Age",
          y = "Intensity"
        ) +
        theme_minimal() +
        theme(legend.position = "none"),
      device = "pdf",
      width = 8, height = 6  # Optional: specify dimensions
    )

    # Only plot 2 sites for visibility
    selected_sites <- simulated_data %>%
      group_by(Group) %>%
      filter(Site %in% unique(Site)[1:2]) %>%
      ungroup()

    # Plot simulated data
    ggsave(
      filename = file.path(output_dir, "sample_sim_data.pdf"),
      plot = ggplot(selected_sites, aes(x = Age, y = Intensity, color = Group)) +
        geom_point(alpha = 0.4) +
        geom_smooth(se = FALSE) +
        facet_wrap(~ Group, scales = "fixed") +
        labs(
          title = "Simulated Methylation Intensity by Group",
          x = "Age", y = "Methylation Intensity"
        ) +
        theme_minimal(),
      device = "pdf",
      width = 8, height = 6  # Adjust dimensions as needed
    )
  }

  # Create structured output data frames
  ages_df <- simulated_data %>%
    select(Person, Age) %>%
    distinct()

  groups_df <- simulated_data %>%
    select(Site, Group) %>%
    distinct()

  methylation_df <- simulated_data %>%
    select(Site, Person, Intensity) %>%
    pivot_wider(
      names_from = Person,
      values_from = Intensity
    ) %>%
    column_to_rownames(var = "Site")

  return(list(ages = ages_df, groups = groups_df, meth = methylation_df))
}
