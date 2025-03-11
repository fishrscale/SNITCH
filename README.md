# SNITCH: Semi-supervised Non-linear Identification and Trajectory Clustering for High-dimensional Data

## Overview

SNITCH (**S**emi-supervised **N**on-linear **I**dentification and **T**rajectory **C**lustering for **H**igh-dimensional data) is an R package designed to analyze DNA methylation patterns across age. The goal of SNITCH is to provide a robust framework to:

- Classify CpG sites into linear, non-linear, or non-correlated trajectories using a semi-supervised approach.
- Identify Differentially Methylated Positions (DMPs), Variably Methylated Positions (VMPs), and non-linear DMPs based on age-related methylation changes.
- Perform Functional Principal Component Analysis (FPCA) specifically on non-linear CpGs to capture complex ageing-related methylation patterns.
- Cluster non-linear CpG sites based on their FPCA scores using multiple clustering algorithms (HDBSCAN, k-means, fuzzy clustering) to uncover biologically meaningful subgroups.

SNITCH is designed to provide an efficient, scalable, and interpretable workflow for methylation trajectory analysis, making it ideal for researchers studying ageing, epigenetics, and biomarker discovery.

## Installation
To install the development version of SNITCH, run:

```r
# Install devtools if not installed
install.packages("devtools")

# Install SNITCH from GitHub
devtools::install_github("fishrscale/SNITCH")
```

## Quick Start
### 1️⃣ Simulate DNA Methylation Data
Generate a dataset with simulated methylation values for multiple CpG sites following different ageing-related trajectories.

```r
library(SNITCH)

# Simulate data with 300 individuals
simulated_data <- simulate_methylation_data(n_people = 300, plot = TRUE)

# Extract data frames
ages_df <- simulated_data$ages_df
groups_df <- simulated_data$groups_df
methylation_df <- simulated_data$methylation_df
```
This will generate and save:
- `functions_sim_data.pdf`: Plots of the predefined methylation patterns.
- `sim_data_viz.pdf`: Visualization of the simulated CpG sites.

### 2️⃣ Classify CpG Sites
Determine which CpG sites follow linear, non-linear, or non-correlated trajectories.
```r
classified_cpgs <- classify_cpg_sites(methylation_df, ages_df$Age)
```
This step will categorize CpG sites and prepare for downstream FPCA analysis.

### 3️⃣ Perform Functional Principal Component Analysis (FPCA) on Non-Linear CpGs
```r
# Select only non-linear CpGs for FPCA
non_linear_cpgs <- classified_cpgs[classified_cpgs$classification == "non-linear", ]$CpG
methylation_nl <- methylation_df[rownames(methylation_df) %in% non_linear_cpgs, ]

# Perform FPCA
fpca_results <- perform_fpca(methylation_nl, ages_df$Age)

# Plot FPCA results
plot_fpca_results(fpca_results, ages_df$Age)
```
This step will generate:
- `SNITCH_fpca_mu.pdf`: Mean methylation trajectory.
- `SNITCH_fpca_var.pdf`: Variance explained by principal components.

### 4️⃣ Run Clustering Analysis on Non-Linear CpGs
Cluster non-linear CpGs based on their FPCA scores using different algorithms.

```r
# Plot diagnostic plots to determine optimal clustering parameters
plot_clustering_diagnostics(fpca_results$scores)

# Perform clustering using HDBSCAN
clusters <- perform_clustering(fpca_results$scores, method = "hdbscan", minPts = 15)
```
The diagnostic plots will help determine:
- The optimal number of clusters for k-means (`SNITCH_kmeans_elbow.pdf`).
- The ideal fuzziness parameter for fuzzy clustering (`SNITCH_fuzzy_distance.pdf`).

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on [GitHub](https://github.com/fishrscale/SNITCH).

## License
This package is licensed under the Apache License 2.0. See the LICENSE file for details.

---

🚀 **SNITCH: Bringing insights to DNA methylation analysis!**

