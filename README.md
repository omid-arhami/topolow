# Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assays

## Overview

`topolow` is an R package that implements a novel, physics-inspired algorithm for antigenic cartography mapping and analysis. The algorithm addresses critical challenges in mapping antigenic relationships from incomplete experimental data, particularly for rapidly evolving pathogens like influenza, SARS-CoV-2, HIV, and dengue viruses.

### Key Advantages

-   **Superior handling of missing data**: Effectively processes datasets with even more than 95% missing values
-   **Complete positioning**: Maps all antigens regardless of dimensionality or data sparsity
-   **Improved accuracy**: Achieves better prediction accuracy than traditional MDS approaches
-   **Stability**: Demonstrates orders of magnitude better consistency across multiple runs
-   **Automatic dimensionality optimization**: Determines optimal mapping dimensions through likelihood-based estimation
-   **Noise reduction**: Effectively reduces experimental noise and bias through network-based error dampening

## Installation

### From GitHub

You can install the development version of topolow directly from GitHub:

``` r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install topolow
devtools::install_github("omid-arhami/topolow")
```

### From Release versions

Alternatively, you can install using the single source file:

1.  Download the latest release
2.  For Windows binary: Install the .zip file
3.  For source package: Install the .tar.gz file

``` r
# For Windows binary
install.packages("path/to/topolow_X.zip", repos = NULL)

# For source package
install.packages("path/to/topolow_X.tar.gz", repos = NULL, type = "source")
```

### Optional Dependencies

For 3D visualization capabilities, install the `rgl` package:

``` r
install.packages("rgl")
```

Note for macOS users: The `rgl` package requires XQuartz to be installed for proper OpenGL support. You can download it from <https://www.xquartz.org/>, then install the downloaded package and restart your computer.

Even without rgl, you can use all core functionality of topolow. The package will automatically fall back to 2D visualizations.

## Quick Start

Here's a simple example to check if Topolow is working and to analytically validate its result.

Let us take 4 points in a 2D space, two reference antigens S/1 and S/2 and two test antigens V/1 and V/2.

S/1 at (0, 0)

S/2 at (3, 0)

V/1 at (2, 2)

V/2 at (0, 4)

The pairwise Euclidean distances between these points are computed as follows:

<!-- For inline equations, use single $ -->

$d(S/1,S/2) = \sqrt{(3-0)^2 + (0-0)^2} = \sqrt{9 + 0} = \sqrt{9} = 3.$

$d(S/1,V/1) = \sqrt{(2-0)^2 + (2-0)^2} = \sqrt{4 + 4} = \sqrt{8} = 2\sqrt{2} \approx 2.828.$

$d(S/1,V/2) = \sqrt{(0-0)^2 + (4-0)^2} = \sqrt{0 + 16} = \sqrt{16} = 4.$

$d(S/2,V/1) = \sqrt{(2-3)^2 + (2-0)^2} = \sqrt{1 + 4} = \sqrt{5} \approx 2.236.$

$d(S/2,V/2) = \sqrt{(0-3)^2 + (4-0)^2} = \sqrt{9 + 16} = \sqrt{25} = 5.$

$d(V/1,V/2) = \sqrt{(0-2)^2 + (4-2)^2} = \sqrt{4 + 4} = \sqrt{8} = 2\sqrt{2} \approx 2.828.$

Imagine we have measured the distances of V/1 against S/1 and S/2, and V/2 against S/1 and S/2. We use Topolow to find the distance between V/1 and V/2 which is missing in the distance matrix (dist_mat in code below). From the analytical calculations we expect d(V/1,V/2) = 2.828.

Remember that this is the simplest example with an analytical solution that lets us verify the result. The true value of using Topolow to find missing distances is when there are many points and many missing distances in the data.

``` r
library(topolow)

# Create a 4×4 simple distance matrix

dist_mat <- matrix(c(
  # S/1  S/2  V/1  V/2
     0,   3,   2.828,   4,    # S/1
     3,   0,  2.236, 5,   # S/2
     2.828,  2.236,   0,   NA,    # V/1
     4, 5 ,  NA,   0     # V/2
), nrow=4)
rownames(dist_mat) <- colnames(dist_mat) <- c("S/1", "S/2", "V/1", "V/2")

# Run TopoLow in 2D
result <- create_topolow_map(dist_mat, ndim=2, mapping_max_iter=1000, 
                             k0=1, cooling_rate=0.0001, c_repulsion=0.001, 
                             write_positions_to_csv = FALSE, verbose = TRUE)

# Investigate the results
print(dist_mat)
print(result$est_distances)
         S/1      S/2      V/1      V/2
S/1 0.000000 3.000027 2.827970 4.000056
S/2 3.000027 0.000000 2.235928 5.000045
V/1 2.827970 2.235928 0.000000 2.828457
V/2 4.000056 5.000045 2.828457 0.000000
```
All of the estimated distances are close to the analytical solution, including model's estimate for the missing distance between V/1 and V/2.


## Reproduction Studies

This package includes computationally intensive examples in the `inst/examples` directory. These examples demonstrate complete use cases in the paper but require computational time and resources.

To run these studies after installing Topolow, you can copy all associated files, subdirectories, and the Rmd files to your machine by calling function below:

``` r
library(topolow)
# Copy to a specific location:
copy_reproduction_examples("path/to/my/examples")
```

Then read through the markdown notebooks and choose which parts you wish to run. There are usually options to use the provided parameters to bypass some parts of the simulations.

Note: Results of time-intensive sections are also provided in csv files and explained at the beginning of each Rmd file.

## How Topolow Works

Topolow employs a novel physical model where:

1.  **Antigens as particles**: Test and reference antigens are represented as particles in an N-dimensional space
2.  **Spring-based connections**: Pairs with known measurements are connected by springs with free lengths equal to their antigenic distance
3.  **Repulsive forces**: Pairs without direct measurements apply repulsive forces to each other, following an inverse square law
4.  **Mass-weighted motion**: Each antigen receives an effective mass proportional to its number of measurements, providing natural regularization
5.  **Cooling schedule**: Spring and repulsion constants gradually decrease during optimization, allowing fine-scale adjustments in final stages

This approach allows Topolow to effectively optimize antigenic positions through a series of one-dimensional calculations, eliminating the need for complex gradient computations required by traditional MDS methods.

## Antigenic Velocity

- **What it is**  
  Computes for each antigen a **velocity vector**  showing the rate and direction of each antigen’s drift. 
  \[
    v_i = \frac{\sum_{j:\,t_j<t_i} K_{ij}\,\frac{x_i - x_j}{t_i - t_j}}
               {\sum_{j:\,t_j<t_i} K_{ij}}
  \]

- **Key parameters**  
  - `sigma_x` (antigenic bandwidth) and `sigma_t` (temporal bandwidth) — default: auto-estimated via Silverman’s rule  
  - `clade_depth` — depth (in tree edges) for phylo-aware clade filtering (Average Leaf-to-Backbone Distance)

## Features

-   **Physics-inspired optimization**: Employs a spring-mass system for robust positioning in high-dimensional spaces
-   **Optimal dimensionality detection**: Automatically determines the best dimensionality through likelihood-based estimation. This is particularly useful for datasets with high levels of missingness and complexity (e.g., due to various serotypes).
-   **Complete antigenic positioning**: Maps all antigens
-   **Noise reduction**: Decreases measurement errors through network-based dampening
-   **Threshold handling**: Properly incorporates low and high reactor thresholds (e.g., \<40) as equality constraints
-   **Cross-validation**: Built-in validation framework for performance assessment
-   **Parallel processing**: Support for multi-core and HPC environments
-   **Visualization tools**: Interactive and publication-ready map generation
-   **Phylogenetically-Aware Clade Detection**: Dynamic depth-based clades (no rooting or branch lengths required) are defined based on Average Leaf-to-Backbone Distance (ALBD) in the tree

## Input Data Format

The algorithm can handle input data in various formats - if the raw input consists of one or multiple long tables with references on columns and challenges on rows, they are converted to the standard matrix form. (See the example scripts in `inst/examples`)

The package accepts distance matrices with the following characteristics:

-   Square symmetric matrices
-   Can contain NA values for missing measurements
-   Can contain threshold indicators (\< or \>) for bounded measurements

## Algorithm Parameters

Key parameters for the TopoLow algorithm:

-   ndim: Number of dimensions (typically 2-20)
-   k0: Initial spring constant (typical range: 0.1-30)
-   cooling_rate: Spring decay rate (typical range: 0.0001-0.1)
-   c_repulsion: Repulsion constant (typical range: 0.00001-0.1)

The optimal values for each data can be determined through adaptive Monte Carlo simulations done by functions `initial_parameter_optimization` and `run_adaptive_sampling`. (See the example scripts in `inst/examples`)

## Performance

Topolow demonstrates significant improvements over traditional MDS approaches:

-   **27 simulated datasets with varying missingness and complexity**: Between 50% to 1000% improved prediction accuracy
-   **H3N2 influenza data (1968 - 2003)**: Similar prediction accuracy to the extensively tested maps in the literature
-   **HIV neutralization data (Subtypes B and C tested)**: 41% improved prediction accuracy
-   **Run-to-run stability**: Orders of magnitude better consistency across multiple runs
-   **Parameter sensitivity**: Performance remains robust across a wide range of parameter values

## Applications

Topolow is particularly valuable for:

-   Understanding antigenic evolution of rapidly evolving viral pathogens
-   Early detection of emerging antigenic variants
-   Predicting antigenic phenotypes for under-characterized strains
-   Amplifying training data for downstream machine learning models
-   Analyzing any continuous and relational phenotype under directional selection pressure

## Using on HPC or SLURM Clusters

When using topolow on HPC systems with SLURM, additional setup might be needed:

1.  Ensure the correct R version is loaded (4.3.2 or newer):

``` bash
module load R/4.4.1
```

2.  Install required dependencies:

``` r
install.packages(c("reshape2", "data.table", "dplyr", "ggplot2"))
```

3.  When submitting SLURM jobs, set the correct R module in the script:

``` r
initial_parameter_optimization(
  # ... other parameters ...
  r_module = "R/4.4.1", # Set this to match your cluster's R module
  use_slurm = TRUE
)
```

## Documentation

See the full documentation of the package and all functionalities in <https://github.com/omid-arhami/topolow/blob/main/build/topolow-manual.pdf>

For detailed documentation of a specific function in Topolow package:

``` r
# View documentation
?function_name
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is protected by a pre-publication license.

-   Researchers can use the software for academic purposes.
-   Redistribution, modification, and commercial use are prohibited before publication.

The license will transition upon publication - see the LICENSE file for details.

## Citation

If you use this package, please cite the article:

Arhami and Rohani, 2025 <doi:10.1101/2025.02.09.637307>

Software doi: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15620983.svg)](https://doi.org/10.5281/zenodo.15620983)

## Contact

-   Maintainer: Omid Arhami
-   Email: [omid.arhami\@uga.edu](mailto:omid.arhami@uga.edu){.email}
-   GitHub: @omid-arhami
