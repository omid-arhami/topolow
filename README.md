# Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assay results

## Overview

`topolow` is an R package that implements a novel algorithm for antigenic cartography mapping and analysis. The package uses a physics-inspired approach combining spring forces and repulsive interactions to find optimal point configurations in high-dimensional spaces.

## Installation

### From GitHub
You can install the development version of topolow directly from GitHub:

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install topolow
devtools::install_github("omid-arhami/topolow")
```

### From Release versions
Alternatively, you can install using the single source file:

1. Download the latest release
2. For Windows binary: Install the .zip file
3. For source package: Install the .tar.gz file

```r
# For Windows binary
install.packages("path/to/topolow_0.1.3.zip", repos = NULL)

# For source package
install.packages("path/to/topolow_0.1.3.tar.gz", repos = NULL, type = "source")
```

### Optional Dependencies

For 3D visualization capabilities, install the `rgl` package:

```r
install.packages("rgl")
```

Note for macOS users: The `rgl` package requires XQuartz to be installed for proper OpenGL support. You can download it from [https://www.xquartz.org/](https://www.xquartz.org/), then install the downloaded package and restart your computer.

Even without rgl, you can use all core functionality of topolow. The package will automatically fall back to 2D visualizations.


## Quick Start

Here's a basic example:

```r
library(topolow)

# Create a simple distance matrix
dist_mat <- matrix(c(0, 2, 3, 2, 0, NA, 3, NA, 0), nrow=3)
rownames(dist_mat) <- colnames(dist_mat) <- c("V/1", "V/2", "S/1")

# Run TopoLow in 2D
result <- create_topolow_map(dist_mat, ndim=2, max_iter=100, 
                      k0=1.0, cooling_rate=0.01, c_repulsion=0.01)

# Investigate the results
print(dist_mat)
print(result$est_distances)
```

## Using on HPC or SLURM Clusters

When using topolow on HPC systems with SLURM, additional setup might be needed:

1. Ensure the correct R version is loaded (4.3.2 or newer):
```bash
module load R/4.4.1
```

2. Install required dependencies:
```r
install.packages(c("reshape2", "data.table", "dplyr", "ggplot2"))
```

3. When submitting SLURM jobs, set the correct R module in the script:
```r
run_parameter_optimization(
  # ... other parameters ...
  r_module = "R/4.4.1", # Set this to match your cluster's R module
  use_slurm = TRUE
)
```

## Documentation

See the full documentation of the package and all functionalities in https://github.com/omid-arhami/topolow/blob/main/build/topolow-manual.pdf

For detailed documentation of a specific function in Topolow package:

```r
# View documentation
?function_name
```

## Reproduction Studies

This package includes computationally intensive examples in the `inst/examples` 
directory. These examples demonstrate complete use cases in the paper but require computational time and resources.

To access these examples after installation:
```r
example_path <- system.file("examples", package="topolow")
# View available examples
list.files(example_path)
```

- `parameter-fitting-h3n2.Rmd`
- `assay_data_analysis_HIV_Subtype.Rmd`
- `methods-comparison-h3n2-hiv.Rmd`
- `synthetic-parameter-fitting.Rmd`
- `synthetic_data_algorithm_comparison.Rmd`
- `procrustes_vignette.Rmd`

To run these studies after installing Topolow, you can copy all associated files, subdirectories, and the Rmd files to your machine by calling function below:

```r
copy_reproduction_examples()
# Or to a specific location:
copy_reproduction_examples("path/to/my/examples")
```

Then read through the markdown notebooks and choose which parts you wish to run. There are usually options to use the provided parameters to bypass some parts of the simulations.

Note: Results of time-intensive sections are provided and explained at the beginning of each Rmd file. 

To generate supplementary figure S-1, the pairwise distances in the original space (10D) versus the pairwise distances after reducing the dimensions to 2, use file `inst/examples/10d_2d_pairwise_distance_comparison_plot.R`

## Features

- Optimized point configuration in high-dimensional spaces
- Handling of missing and thresholded measurements
- Processing of antigenic assay data
- Interactive visualization of antigenic maps
- Cross-validation and error analysis
- Network structure analysis
- Support for parallel processing and high-performance computing environments

## Input Data Format

The algorithm can handle input data in various formats - if the raw input consists of one or multiple long tables with references on columns (rows) and challenges on rows (columns), they are converted to the standard matrix form.

The package accepts distance matrices with the following characteristics:

* Square symmetric matrices
* Can contain NA values for missing measurements
* Can contain threshold indicators (< or >) for bounded measurements
* Row and column names should identify antigens (V/) and antisera (S/)

## Algorithm Parameters

Key parameters for the TopoLow algorithm:

* ndim: Number of dimensions (typically 2-20)
* k0: Initial spring constant (typical range: 0.1-30)
* cooling_rate: Spring decay rate (typical range: 0.0001-0.1)
* c_repulsion: Repulsion constant (typical range: 0.00001-0.1)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is protected by a pre-publication license.

* Researchers can use the software for academic purposes.
* Redistribution, modification, and commercial use are prohibited before publication.

The license will transition upon publication - see the LICENSE file for details.

## Citation

If you use this package, please cite:

Arhami and Rohani, 2025 [doi:to be added]

## Contact

- Maintainer: Omid Arhami
- Email: omid.arhami@uga.edu
- GitHub: @omid-arhami
