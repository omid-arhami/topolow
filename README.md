# Topolow: Force-Directed Euclidean Embedding of Dissimilarity Data

## Overview

`topolow` is an R package that implements a novel, physics-inspired algorithm for **Euclidean embedding** of potentially non-metric, sparse, and noisy dissimilarity data. The algorithm converts dissimilarity matrices into valid Euclidean coordinate systems, making the data compatible with standard statistical and machine learning tools like PCA, clustering, and regression.

**The Problem**: Many datasets contain dissimilarity measurements that violate metric axioms (symmetry, triangle inequality) or are highly sparse with missing values. Standard methods like Multidimensional Scaling (MDS) struggle with such data, leading to poor embeddings or complete failure.

**The Solution**: Topolow uses a physics-inspired approach that models objects as particles connected by springs (for known dissimilarities) and repulsive forces (for missing pairs). This gradient-free optimization is robust to local optima and handles non-metric data naturally.

### Key Advantages

- **Handles non-metric data**: Works with dissimilarities that violate metric axioms
- **Superior performance on sparse data**: Effectively processes datasets with >95% missing values
- **Calculates antigenic velocity vectors (biology)**: Measures the rate and direction of viral evolution, offering early warnings of lineage replacements.
- **Robust optimization**: Gradient-free algorithm avoids local optima
- **Statistical foundation**: Maximum likelihood estimation under Laplace error model
- **Automatic parameter optimization**: Determines optimal dimensions and parameters through cross-validation
- **Handles censored data**: Properly incorporates threshold measurements (e.g., "<1", ">64")
- **Network-based error dampening**: Reduces experimental noise through interconnected force system
- **Complete positioning**: Maps all objects regardless of data sparsity

## Installation

### From CRAN
[![](https://cranlogs.r-pkg.org/badges/topolow)](https://cran.r-project.org/package=topolow)

(Download version 2 from GitHub to use the wizard function Euclidify.)

```r
install.packages("topolow")
```

### From GitHub

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install topolow
devtools::install_github("omid-arhami/topolow")
```

## Quick Start

### Simple Example: Embedding 5 Points

Let's embed 4 points with known coordinates to validate the algorithm:

```r
library(topolow)

# Known coordinates:
# S1 at (0,0), S2 at (3,0), S3 at (4,4), V1 at (2,2), V2 at (0,4)
# We'll provide some distances and let Topolow infer the missing ones
coordinates <- data.frame(
  point = c("S1", "S2", "S3", "V1", "V2"),
  x = c(0, 3, 4, 2, 0),
  y = c(0, 0, 4, 2, 4)
)
# Create a matrix with just the x,y coordinates for distance calculation
coord_matrix <- as.matrix(coordinates[, c("x", "y")])
rownames(coord_matrix) <- coordinates$point
# Calculate the distance matrix using Euclidean distance
dist_mat <- as.matrix(dist(coord_matrix))
# Remove a known distance to see if topolow predicts it accurately
dist_mat["V1", "V2"] <- NA
dist_mat["V2", "V1"] <- NA

# Run Topolow (manual parameters for quick demo)
result <- euclidean_embedding(
  dissimilarity_matrix = dist_mat,
  ndim = 2,
  mapping_max_iter = 1000,
  k0 = 5, 
  cooling_rate = 0.03, 
  c_repulsion = 0.7,
  verbose = TRUE
)

# Check results
print("Original matrix:")
print(dist_mat)
print("Estimated distances:")
print(result$est_distances)

# The missing distance V1-V2 should be approximately 2.83
```

### Automatic Optimization with Euclidify

For real applications, use `Euclidify()` which automatically optimizes all parameters:

```r
# Using automatic parameter optimization
result_auto <- Euclidify(
  dissimilarity_matrix = dist_mat,
  ndim_range = c(2, 4),
  output_dir = tempdir(),  # Required for optimization files
  n_initial_samples = 50,  # Reduced for quick demo
  n_adaptive_samples = 200,
  folds = 4,
  verbose = "standard"
)

# Extract optimized distances
est_dist <- result_auto$est_distances
print(est_dist)

# View optimal parameters found
print(result_auto$optimal_params)
```

## Applications and Examples

### 1. Antigenic Mapping: Viral Evolution with Temporal Visualization

For immunologists studying viral evolution and vaccine effectiveness, Topolow can generate insightful visualizations. A key feature for this application is the calculation of antigenic velocity vectors, which show the rate and direction of antigenic drift for each virus against its recent predecessors. These vectors can reveal evolutionary trends and highlight potential vaccine-escape variants by revealing fast movements.
```r
# Example: H3N2 Influenza A antigenic evolution with temporal mapping and velocity vectors

# Create a more comprehensive antigenic dataset with temporal information
antigen_data <- data.frame(
  virus = c("A/H3N2/HK/1968", "A/H3N2/EN/1972", "A/H3N2/VI/1975", "A/H3N2/TX/1977", 
            "A/H3N2/BK/1979", "A/H3N2/SI/1987", "A/H3N2/BE/1989", "A/H3N2/BE/1992",
            "A/H3N2/WU/1995", "A/H3N2/SY/1997", "A/H3N2/FU/2002", "A/H3N2/WI/2005"),
  serum = rep(c("anti-HK68/1968", "anti-EN72/1972", "anti-VI75/1975", "anti-TX77/1977", "anti-BK79/1979", 
                "anti-SI87/1987"), each = 12),
  titer = c(2560, 1280, 640, 320, 160, 80, 40, "<40", "<40", "<40", "<40", "<40",
            640, 2560, 1280, 640, 320, 160, 80, 40, "<40", "<40", "<40", "<40",
            320, 640, 2560, 1280, 640, 320, 160, 80, 40, "<40", "<40", "<40",
            160, 320, 640, 2560, 1280, 640, 320, 160, 80, 40, "<40", "<40",
            80, 160, 320, 640, 2560, 1280, 640, 320, 160, 80, 40, "<40",
            "<40", 80, 160, 320, 640, 2560, 1280, 640, 320, 160, 80, 40),
  year = rep(c(1968, 1972, 1975, 1977, 1979, 1987, 1989, 1992, 1995, 1997, 2002, 2005), 6)
)

# Convert titers to dissimilarity matrix
results <- process_antigenic_data(
  data = antigen_data,
  antigen_col = "virus",
  serum_col = "serum", 
  value_col = "titer",
  is_similarity = TRUE,  # Titers are similarities
  scale_factor = 10      # Base dilution factor of HI assay
)
antigenic_matrix <- results$matrix

# Create antigenic map with temporal information
antigenic_map <- Euclidify(
  dissimilarity_matrix = antigenic_matrix,
  ndim_range = c(3, 6),
  folds = 10,
  output_dir = tempdir(),
  verbose = "standard"
)

# Prepare data for temporal visualization
positions_df <- data.frame(
  V1 = antigenic_map$positions[, 1],
  V2 = antigenic_map$positions[, 2],
  V3 = antigenic_map$positions[, 3],
  name = rownames(antigenic_map$positions),
  year = as.numeric(sub(".*/([0-9]+).*", "\\1", rownames(antigenic_map$positions))),
  antigen = grepl("^V/", rownames(antigenic_map$positions)),
  antiserum = grepl("^S/", rownames(antigenic_map$positions))
)
  

# Configure visualization aesthetics
aesthetic_config <- new_aesthetic_config(
  point_size = 3.0,
  point_alpha = 0.8,
  point_shapes = c(antigen = 16, antiserum = 5),  # Circle for antigens, diamond for sera
  gradient_colors = list(low = "blue", high = "red"),
  show_labels = TRUE,
  show_title = TRUE,
  title_size = 14,
  axis_text_size = 11,
  show_legend = TRUE,
  legend_position = "right",
  arrow_alpha = 0.7
)

layout_config <- new_layout_config(
  width = 10,
  height = 8,
  save_plot = FALSE,
  arrow_plot_threshold = 0.15,  # Show arrows for significant movements
  show_grid = TRUE,
  grid_type = "major"
)

annotation_config <- new_annotation_config(
  notable_points = c("V/A/H3N2/HK/1968", "V/A/H3N2/FU/2002", "V/A/H3N2/WI/2005")
)

# Create temporal antigenic map with velocity vectors
temporal_plot <- plot_temporal_mapping(
  df_coords = positions_df,
  ndim = 3,
  draw_arrows = TRUE,           # Enable velocity vectors
  annotate_arrows = TRUE,       # Label arrows with strain names
  sigma_t = 2.0,               # Temporal bandwidth (years)
  sigma_x = 1.5,               # Spatial bandwidth (antigenic units)
  aesthetic_config = aesthetic_config,
  layout_config = layout_config,
  annotation_config = annotation_config,
  output_dir = tempdir()
)

print(temporal_plot)

# Create clustered view by antigenic epoch
positions_df$cluster <- cut(positions_df$year, 
                           breaks = c(1965, 1975, 1985, 1995, 2010),
                           labels = c("Pre-1975", "1975-1985", "1985-1995", "Post-1995"))

# Visualize antigenic clusters with evolutionary relationships
cluster_plot <- plot_cluster_mapping(
  df_coords = positions_df,
  ndim = 3,
  draw_arrows = TRUE,
  show_one_arrow_per_cluster = TRUE,  # One representative arrow per cluster
  aesthetic_config = aesthetic_config,
  layout_config = layout_config,
  annotation_config = annotation_config
)

print(cluster_plot)

# Generate 3D visualization for enhanced perspective
if (requireNamespace("rgl", quietly = TRUE)) {
  # Interactive 3D antigenic map
  plot_3d <- plot_3d_mapping(
    positions_df,
    ndim = 3,
    aesthetic_config = aesthetic_config,
  layout_config = layout_config
  )
  
  cat("3D antigenic map created. Use mouse to rotate and zoom.\n")
}
```

### 2. General Data Science: Customer Similarity

```r
# Example: Customer behavior dissimilarity
customer_data <- data.frame(
  customer = rep(paste0("Cust", 1:5), each = 5),
  product = rep(paste0("Prod", 1:5), 5),
  dissimilarity = c(0, 2.1, 3.5, 1.8, 4.2,
                   2.1, 0, 1.9, 3.1, 2.8,
                   3.5, 1.9, 0, 2.4, 3.7,
                   1.8, 3.1, 2.4, 0, 2.9,
                   4.2, 2.8, 3.7, 2.9, 0)
)

# Convert to matrix format
dissim_matrix <- list_to_matrix(
  data = customer_data,
  object_col = "customer",
  reference_col = "product", 
  value_col = "dissimilarity",
  is_similarity = FALSE
)

# Embed in Euclidean space
customer_map <- Euclidify(
  dissimilarity_matrix = dissim_matrix,
  output_dir = tempdir(),
  ndim_range = c(2, 4),
  verbose = "standard"
)

plot(customer_map$positions, main = "Customer Behavior Map")
text(customer_map$positions[,1] + jitter(rep(0, nrow(customer_map$positions)), amount = 0.2), 
     customer_map$positions[,2] + jitter(rep(0, nrow(customer_map$positions)), amount = 0.2), 
     labels = rownames(customer_map$positions), pos = 3, cex = 0.5)
```

### 3. Handling Large and Sparse Data

```r
# Example: Large symmetric sparse matrix with realistic structure
set.seed(12345)  # For reproducibility

# Generate a large, realistic sparse dissimilarity matrix
n_objects <- 50
object_names <- paste0("Object_", sprintf("%02d", 1:n_objects))

# Create base coordinates in 3D space with clustered structure
cluster_centers <- matrix(c(
  c(0, 0, 0),      # Cluster 1
  c(5, 0, 0),      # Cluster 2  
  c(0, 5, 0),      # Cluster 3
  c(5, 5, 0),      # Cluster 4
  c(2.5, 2.5, 3)   # Cluster 5
), ncol = 3, byrow = TRUE)

# Assign objects to clusters
cluster_assignments <- sample(1:5, n_objects, replace = TRUE, 
                             prob = c(0.25, 0.25, 0.20, 0.20, 0.10))

# Generate coordinates with cluster structure + noise
true_coordinates <- matrix(0, n_objects, 3)
for(i in 1:n_objects) {
  cluster_id <- cluster_assignments[i]
  true_coordinates[i, ] <- cluster_centers[cluster_id, ] + rnorm(3, 0, 0.8)
}

rownames(true_coordinates) <- object_names

# Calculate complete Euclidean distance matrix
complete_distances <- as.matrix(dist(true_coordinates))

# Add realistic measurement noise (5% coefficient of variation)
noisy_distances <- complete_distances * (1 + rnorm(n_objects^2, 0, 0.05))
noisy_distances <- pmax(noisy_distances, 0.1)  # Minimum distance threshold

# Make symmetric and zero diagonal
noisy_distances[lower.tri(noisy_distances)] <- t(noisy_distances)[lower.tri(noisy_distances)]
diag(noisy_distances) <- 0

# Introduce structured sparsity (85% missing data)
# Objects in the same cluster are more likely to have measurements
total_pairs <- n_objects * (n_objects - 1) / 2
target_missing_pairs <- round(total_pairs * 0.85)  # 85% sparsity

# Generate upper triangular indices for sampling
upper_tri_indices <- which(upper.tri(noisy_distances), arr.ind = TRUE)

# Create sampling weights: higher probability for within-cluster pairs
sampling_weights <- numeric(nrow(upper_tri_indices))
for(k in 1:nrow(upper_tri_indices)) {
  i <- upper_tri_indices[k, 1]
  j <- upper_tri_indices[k, 2]
  if(cluster_assignments[i] == cluster_assignments[j]) {
    sampling_weights[k] <- 0.3  # Lower chance of being missing (within cluster)
  } else {
    sampling_weights[k] <- 1.0  # Higher chance of being missing (between clusters)
  }
}

# Sample pairs to remove
missing_pair_indices <- sample(
  nrow(upper_tri_indices), 
  target_missing_pairs, 
  prob = sampling_weights
)

# Create sparse matrix
sparse_matrix <- noisy_distances
sparse_matrix[upper_tri_indices[missing_pair_indices, ]] <- NA
sparse_matrix[upper_tri_indices[missing_pair_indices, c(2,1)]] <- NA

rownames(sparse_matrix) <- colnames(sparse_matrix) <- object_names

# Calculate actual sparsity
actual_sparsity <- sum(is.na(sparse_matrix)) / (n_objects * (n_objects-1)) * 100

cat("=== SPARSE DATA EXAMPLE ===\n")
cat("Generated sparse dissimilarity matrix:\n")
cat("- Matrix size:", n_objects, "x", n_objects, "objects\n")
cat("- Actual sparsity:", round(actual_sparsity, 1), "% missing data\n")
cat("- Clustering structure: 5 clusters with", table(cluster_assignments), "objects each\n")
cat("- Noise level: 5% coefficient of variation\n")
cat("- Within-cluster connectivity preserved for realism\n")

cat("=== Check the connectivity of the data graph ===\n")
# Network structure analysis to make sure there are no separate islands in the data
network_analysis <- analyze_network_structure(sparse_matrix)
network_plot <- plot_network_structure(network_analysis)
print(network_plot)

# Demonstrate Topolow's superior sparse data handling
cat("\n=== EMBEDDING SPARSE DATA ===\n")

# Topolow embedding with automatic optimization
sparse_result <- Euclidify(
  dissimilarity_matrix = sparse_matrix,
  ndim_range = c(2, 8),
  output_dir = tempdir(),
  n_initial_samples = 50,
  n_adaptive_samples = 150,
  folds = 20,
  verbose = "standard"
)
print(sparse_result$optimal_params)

# Evaluate embedding quality by visualizing the results 
# and coloring by original cluster for validation
plot_colors <- rainbow(5)[cluster_assignments]

plot(sparse_result$positions[, 1:2], 
     col = plot_colors, 
     pch = 19, 
     cex = 1.2,
     main = paste("Topolow Embedding of Sparse Data\n", 
                 round(actual_sparsity, 1), "% Missing Values"),
     xlab = "Dimension 1", 
     ylab = "Dimension 2")

legend("topleft", 
       legend = paste("Cluster", 1:5), 
       col = rainbow(5), 
       pch = 19, 
       cex = 0.6)

# Add text labels for some points
text(sparse_result$positions[1:10, 1:2], 
     labels = object_names[1:10], 
     pos = 3, 
     cex = 0.6)
```

## How Topolow Works

Topolow employs a novel physical model where:

1. **Objects as particles**: Each object becomes a particle in N-dimensional space
2. **Spring forces**: Pairs with known dissimilarities are connected by springs with rest lengths equal to the measured dissimilarity
3. **Repulsive forces**: Pairs without measurements apply repulsive forces, preventing particle collapse
4. **Mass-weighted motion**: Each particle has effective mass proportional to its number of measurements
5. **Statistical foundation**: Minimizes Mean Absolute Error, equivalent to Maximum Likelihood Estimation under Laplace error model
6. **Gradient-free optimization**: Sequential pairwise updates avoid local optima common in gradient-based methods
7. **Cooling schedule**: Force constants gradually decrease, allowing fine-scale adjustments

**Key Distinction from MDS**: While MDS methods impute missing values and calculate eigenvalues or gradient vectors, Topolow works directly with the structure in the data and uses physics-inspired forces for robust optimization.

## Features

### Core Algorithm
- **Physics-inspired optimization**: Spring-mass system for robust positioning
- **Gradient-free**: Avoids local optima through stochastic pairwise interactions  
- **Non-metric compatibility**: Handles data violating metric axioms
- **Sparse data handling**: No imputation required for missing values
- **Censored data support**: Handles threshold indicators (<, >) as constraints

### Parameter Optimization
- **Automatic dimension selection**: Likelihood-based optimal dimensionality
- **Adaptive Monte Carlo**: Parameter optimization focusing on high-likelihood regions
- **Cross-validation**: Built-in k-fold validation for robust parameter estimates

### Practical Features
- **Parallel processing**: Multi-core support for large datasets
- **Convergence diagnostics**: Automatic convergence detection and monitoring
- **Flexible input formats**: Handles matrices, data frames, similarity/dissimilarity data
- **Comprehensive visualization**: 2D, 3D, temporal, and cluster-based plotting tools, including antigenic velocity vectors to track evolutionary drift.

## Input Data Formats

Topolow accepts data in multiple formats:

### 1. Matrix Format
```r
# Direct dissimilarity matrix
dist_matrix <- matrix(c(0, 1.2, 2.1, 1.2, 0, 1.8, 2.1, 1.8, 0), nrow=3)
```

### 2. Long Format (List)
```r
# Convert long format to matrix
long_data <- data.frame(
  object = c("A", "B", "C"),
  reference = c("X", "Y", "Z"), 
  value = c(1.2, 1.8, 2.1)
)
matrix_data <- list_to_matrix(long_data, "object", "reference", "value")
```

### 3. Threshold Measurements
```r
# Data with detection limits
threshold_matrix <- matrix(c(0, ">64", "<40", ">64", 0, "20", "<40", "20", 0), nrow=3)
```

## Performance Advantages

Based on empirical evaluations in the Bioinformatics paper (Arhami and Rohani, 2025 https://doi.org/10.1093/bioinformatics/btaf372):

- **50-1000% improved accuracy** over MDS on simulated datasets with varying missingness
- **56% improvement** on dengue virus antigenic data
- **41% improvement** on HIV neutralization data  
- **Orders of magnitude better stability** across multiple runs
- **Superior performance maintained** even at 90% data sparsity
- **Robust performance** across wide parameter ranges

## Applications

### Immunology and Virology  
- **Antigenic cartography**: Map viral evolution and vaccine effectiveness
- **Track evolutionary drift**: Use antigenic velocity vectors to identify fast-evolving samples and potential vaccine escape variants
- **Immune repertoire analysis**: Visualize antibody/TCR similarity
- **Vaccine design**: Identify antigenic variants and coverage gaps
- **Pathogen surveillance**: Track emerging variants
- **Cross-reactivity analysis**: Understand immune cross-protection

### General Data Science
- **Customer segmentation**: Embed behavioral dissimilarity data
- **Recommendation systems**: Map user-item preference dissimilarities
- **Network analysis**: Embed graph distances into Euclidean space
- **Dimensionality reduction**: Robust alternative to PCA for non-Euclidean data
- **Anomaly detection**: Identify outliers in dissimilarity space

### Bioinformatics and Computational Biology
- **Protein structure analysis**: Embed structural dissimilarity measures
- **Phylogenetic analysis**: Visualize evolutionary distances
- **Gene expression**: Map correlation-based dissimilarities
- **Drug discovery**: Embed molecular dissimilarity for compound analysis

### Other Domains
- **Psychology**: Embed perceptual or cognitive dissimilarity data
- **Marketing**: Map brand perception and consumer preferences  
- **Geographic analysis**: Handle incomplete distance data
- **Social network analysis**: Embed relationship dissimilarities

## Algorithm Parameters

Key parameters of `euclidean_embedding()` for manual optimization:

- **ndim**: Number of dimensions (typically 2-10)
- **k0**: Initial spring constant (typical range: 0.1-30)
- **cooling_rate**: Parameter decay rate (typical range: 0.0001-0.05)
- **c_repulsion**: Repulsion constant (typical range: 0.0001-0.4)

**Recommendation**: Use `Euclidify()` for automatic parameter optimization if you are not willing to invest time on manual tuning and investigation.

## Comparison with Other Methods

| Method | Non-metric Compatible | Missing Data | Sparse Data | Gradient-free | Stability |
|--------|----------------|-------------|--------------|---------------|-----------|
| Topolow | ✅ | ✅ | ✅ | ✅ | ✅ |
| Classical MDS | ❌ | ❌ (requires imputation) | ❌ | ✅ | ✅ |
| Iterative MDS | ❌ | ❌ (requires imputation) | ❌ | ❌ | ❌ |
| t-SNE | ❌ | ❌ | ❌ | ❌ | ❌ |
| UMAP | ❌ | ❌ | ❌ | ❌ | ⚠️ |

## Using on HPC or SLURM Clusters

topolow can be used on a single HPC system and leverage the larger number of cores by increasing `max_cores` parameters. Distributed processing using SLURM is supported in versions prior to 1.0.0.

### Optional Dependencies

For 3D visualization capabilities, install the `rgl` package:

```r
install.packages("rgl")
```

Note for macOS users: The `rgl` package requires XQuartz. Download from https://www.xquartz.org/, install, and restart your computer.

## Documentation

Full documentation available at:
```r
# View documentation for specific functions
?Euclidify
?euclidean_embedding
?initial_parameter_optimization

# Package overview
help(package = "topolow")
```

## Citation

If you use this package, please cite:

Omid Arhami, Pejman Rohani, Topolow: a mapping algorithm for antigenic cross-reactivity and binding affinity assays, *Bioinformatics*, Volume 41, Issue 7, July 2025, btaf372, https://doi.org/10.1093/bioinformatics/btaf372

```bibtex
@article{10.1093/bioinformatics/btaf372,
    author = {Arhami, Omid and Rohani, Pejman},
    title = {Topolow: a mapping algorithm for antigenic cross-reactivity and binding affinity assays},
    journal = {Bioinformatics},
    volume = {41},
    number = {7},
    pages = {btaf372},
    year = {2025},
    month = {06},
    abstract = {Understanding antigenic evolution through cross-reactivity assays is crucial for tracking rapidly evolving pathogens requiring regular vaccine updates. However, existing cartography methods, commonly based on multidimensional scaling (MDS), face significant challenges with sparse and complex data, producing incomplete and inconsistent maps. There is an urgent need for robust computational methods that can accurately map antigenic relationships from incomplete experimental data while maintaining biological relevance, especially given that more than 95\% of possible measurements could be missing in large-scale studies.We present Topolow, an algorithm that transforms cross-reactivity and binding affinity measurements into accurate positions in a phenotype space. Using a physics-inspired model, Topolow achieved comparable prediction accuracy to MDS for H3N2 influenza and 56\% and 41\% improved accuracy for dengue and HIV, while maintaining complete positioning of all antigens. The method effectively reduces experimental noise and bias, determines optimal dimensionality through likelihood-based estimation, avoiding distortions due to insufficient dimensions, and demonstrates orders of magnitude better stability across multiple runs. We also introduce antigenic velocity vectors, which measure the rate of antigenic advancement of each isolate per unit of time against its temporal and evolutionary related background, revealing the underlying antigenic relationships and cluster transitions.Topolow is implemented in R and freely available at https://doi.org/10.5281/zenodo.15620983 and https://github.com/omid-arhami/topolow.},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btaf372},
    url = {https://doi.org/10.1093/bioinformatics/btaf372},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/41/7/btaf372/63582086/btaf372.pdf},
}
```

Pre-print of the article explaining Euclidean embedding and mathematical properties of the algorithm with evaluations:

Omid Arhami, Pejman Rohani, Topolow: Force-Directed Euclidean Embedding of Dissimilarity Data with Robustness Against Non-Metricity and Sparsity, arXiv:2508.01733, https://doi.org/10.48550/arXiv.2508.01733

```bibtex
@misc{arhami2025topolowforcedirectedeuclideanembedding,
      title={Topolow: Force-Directed Euclidean Embedding of Dissimilarity Data with Robustness Against Non-Metricity and Sparsity}, 
      author={Omid Arhami and Pejman Rohani},
      year={2025},
      eprint={2508.01733},
      archivePrefix={arXiv},
      primaryClass={cs.CG},
      doi={10.48550/arXiv.2508.01733},
      url={https://arxiv.org/abs/2508.01733}, 
}
```

Software DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15620983.svg)](https://doi.org/10.5281/zenodo.15620983)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is distributed under BSD-3-Clause license
YEAR: 2025
COPYRIGHT HOLDER: Omid Arhami

## Contact

- **Maintainer**: Omid Arhami
- **Email**: [omid.arhami@uga.edu](mailto:omid.arhami@uga.edu)
- **GitHub**: @omid-arhami
- **Issues**: https://github.com/omid-arhami/topolow/issues