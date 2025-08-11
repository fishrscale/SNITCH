# SNITCH: Semi-supervised Non-linear Identification and Trajectory Clustering for High-dimensional Data

## Cite
If you're using SNITCH in your work, please cite:


## Overview

SNITCH (**S**emi-supervised **N**on-linear **I**dentification and **T**rajectory **C**lustering for **H**igh-dimensional data) is an R package designed to analyze DNA methylation patterns across age. The goal of SNITCH is to provide a robust framework to:

- Classify CpG sites into linear, non-linear, or non-correlated trajectories using a semi-supervised approach.
- Identify Differentially Methylated Positions (DMPs), Variably Methylated Positions (VMPs), and non-linear DMPs based on age-related methylation changes.
- Perform Functional Principal Component Analysis (FPCA) specifically on non-linear CpGs to capture complex ageing-related methylation patterns.
- Cluster non-linear CpG sites based on their FPCA scores using multiple clustering algorithms (HDBSCAN, k-means, fuzzy clustering) to uncover biologically meaningful subgroups.

SNITCH is designed to provide an efficient, scalable, and interpretable workflow for methylation trajectory analysis, making it ideal for researchers studying ageing, epigenetics, and biomarker discovery.

## Installation
To install the development version of SNITCH, run:

```{r}
# Install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install SNITCH from GitHub
devtools::install_github("fishrscale/SNITCH")
```

## Quick Start

### 1️⃣ Simulate DNA Methylation Data

Generate a dataset with simulated methylation values for multiple CpG sites following different ageing-related trajectories.

```{r}
library(SNITCH)

# Simulate data with 300 individuals
simulated_data <- simulate_methylation_data(n_people = 300, plot = TRUE)

# Extract data frames
ages_df        <- simulated_data$ages_df
groups_df      <- simulated_data$groups_df
methylation_df <- simulated_data$methylation_df
```

This will generate and save:
- `functions_sim_data.pdf`: Plots of the predefined methylation patterns.
- `sim_data_viz.pdf`: Visualization of the simulated CpG sites.

### 2️⃣ Classify CpG Sites

Determine which CpG sites follow linear, non-linear, or non-correlated trajectories.

```{r}
classified_cpgs <- classify_cpg_sites(methylation_df, ages_df$Age)
head(classified_cpgs)
```

This step will categorize CpG sites and prepare for downstream FPCA analysis.

### 3️⃣ Perform Functional Principal Component Analysis (FPCA) on Non-Linear CpGs

```{r}
# Select only non-linear CpGs for FPCA
non_linear_cpgs <- classified_cpgs[classified_cpgs$classification == "non-linear", ]$CpG
methylation_nl  <- methylation_df[rownames(methylation_df) %in% non_linear_cpgs, ]

# Perform FPCA
fpca_results <- perform_fpca(methylation_nl, ages_df$Age)

# Plot FPCA results (files are saved by the function)
plot_fpca_results(fpca_results, ages_df$Age)
```

This step will generate:
- `SNITCH_fpca_mu.pdf`: Mean methylation trajectory.
- `SNITCH_fpca_var.pdf`: Variance explained by principal components.

### 4️⃣ Run Clustering Analysis on Non-Linear CpGs

Use your favorite **unsupervised clustering** strategy to group non-linear CpGs using their **FPCA scores**. Below we illustrate three common options: **k-means**, **fuzzy clustering (MFUZZ)**, and **HDBSCAN**, and compare them with **Adjusted Rand Index (ARI)** and **Adjusted Mutual Information (AMI)** against the ground truth labels. Use the diagnostics to pick sensible parameters for your data.

```{r}
library(dbscan)     # HDBSCAN
library(Mfuzz)      # Fuzzy c-means
library(Biobase)    # ExpressionSet
library(factoextra) # k-means diagnostics
library(ggplot2)
library(ggrepel)
library(aricode)    # ARI/AMI
library(dplyr)
library(tibble)

# Ground truth from simulation
truth <- groups_df$Group

# We'll update copies of 'classified_cpgs' so we can compare methods cleanly
class_work <- classified_cpgs
nl_mask    <- class_work$classification == "non-linear"

# Container for method comparison
df_comp <- tibble(Method = character(), ARI = double(), AMI = double())

add_metrics <- function(pred_labels, method_label) {
  tibble(
    Method = method_label,
    ARI    = aricode::ARI(pred_labels, truth),
    AMI    = aricode::AMI(pred_labels, truth)
  )
}

# -------------------------------------------------------------
# A) FUZZY Clustering with MFUZZ
# -------------------------------------------------------------
eset <- new("ExpressionSet", exprs = fpca_results$scores)

# Estimate m and inspect minimum centroid distance
m_opt <- mestimate(eset)
tmp   <- Dmin(eset, m_opt, crange = seq(2, 20, 1), repeats = 3, visu = TRUE)
ggplot2::ggsave("SNITCH_fuzzy_distance.pdf", width = 6, height = 4)

# Choose cluster count (example choice shown here)
c_opt <- 11

# Fit fuzzy c-means
mfuzz_fit <- mfuzz(eset, c = c_opt, m = m_opt)

# Hard-assign each CpG to the cluster with highest membership
fuzzy_assign <- apply(mfuzz_fit$membership, 1, which.max)
pred_fuzzy   <- paste0("NL_", fuzzy_assign)

# Update and score
cw_fuzzy <- class_work
cw_fuzzy$classification[nl_mask] <- pred_fuzzy
df_comp  <- bind_rows(df_comp, add_metrics(cw_fuzzy$classification, "SNITCH + Fuzzy"))

# -------------------------------------------------------------
# B) K-MEANS Clustering
# -------------------------------------------------------------
p_elbow <- fviz_nbclust(fpca_results$scores, kmeans, k.max = 20, method = "wss") +
  labs(title = "Elbow Method for K-Means",
       x = "Number of Clusters (K)",
       y = "Total Within-Cluster Sum of Squares (WCSS)") +
  theme_minimal()
ggplot2::ggsave("SNITCH_kmeans_elbow.pdf", p_elbow, width = 6, height = 4)

k_opt  <- 10
km_fit <- kmeans(as.matrix(fpca_results$scores), centers = k_opt, nstart = 25)
pred_km <- paste0("NL_", km_fit$cluster)

cw_km <- class_work
cw_km$classification[nl_mask] <- pred_km
df_comp <- bind_rows(df_comp, add_metrics(cw_km$classification, "SNITCH + K-Means"))

# -------------------------------------------------------------
# C) HDBSCAN
# -------------------------------------------------------------
hdb_fit <- hdbscan(as.matrix(fpca_results$scores), minPts = 5)
pred_hdb <- paste0("NL_", hdb_fit$cluster)

cw_hdb <- class_work
cw_hdb$classification[nl_mask] <- pred_hdb
df_comp <- bind_rows(df_comp, add_metrics(cw_hdb$classification, "SNITCH + HDBSCAN"))

# -------------------------------------------------------------
# D) Compare methods with ARI/AMI scatter
# -------------------------------------------------------------
p_comp <- ggplot(df_comp, aes(x = ARI, y = AMI, label = Method)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  coord_equal() +
  labs(title = "Clustering Agreement on Simulated Data",
       x = "Adjusted Rand Index (ARI)",
       y = "Adjusted Mutual Information (AMI)") +
  theme_minimal()

print(df_comp)
print(p_comp)
ggplot2::ggsave("SNITCH_clustering_ari_ami.pdf", p_comp, width = 6, height = 4)
```

**What you’ll get**
- `SNITCH_kmeans_elbow.pdf`: k-means elbow diagnostic.
- `SNITCH_fuzzy_distance.pdf`: MFUZZ distance diagnostic.
- `SNITCH_clustering_ari_ami.pdf`: ARI vs AMI scatter comparing methods.
- Console table of ARI/AMI for **SNITCH + Fuzzy**, **SNITCH + K-Means**, and **SNITCH + HDBSCAN**.

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on [GitHub](https://github.com/fishrscale/SNITCH).

## License
This package is licensed under the Apache License 2.0. See the LICENSE file for details.

---

🚀 **SNITCH: Bringing Non-Linear insights to your analysis!**
