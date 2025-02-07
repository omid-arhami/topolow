# topolow: Antigenic Cartography Using Topological Optimization

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
install.packages("path/to/topolow_0.1.0.zip", repos = NULL)

# For source package
install.packages("path/to/topolow_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick Start

Here's a basic example:

```r
library(topolow)

# Create a simple distance matrix
dist_mat <- matrix(c(0, 2, 3, 2, 0, NA, 3, NA, 0), nrow=3)
rownames(dist_mat) <- colnames(dist_mat) <- c("V/1", "V/2", "S/1")

# Run TopoLow in 2D
result <- topolow_full(dist_mat, ndim=2, max_iter=100, 
                      k0=1.0, k_decay=0.01, cqq=0.01)

# Visualize results

plot(result)
```

## Documentation

For detailed documentation:

```r
# View package documentation
?topolow

# List available vignettes
vignette(package = "topolow")

# Read specific vignette
vignette("parameter-fitting", package = "topolow")
```
## Reproduction Studies

Additional vignettes with detailed computational studies are available in the `reproduction_examples/` directory:

- `parameter-fitting-h3n2.Rmd`
- `synthetic-parameter-fitting.Rmd`
- `synthetic_data_algorithm_comparison.Rmd`

To run these studies:

```r
# Clone the repository
git clone https://github.com/omid-arhami/topolow.git

# Run specific reproduction study
rmarkdown::render("reproduction_examples/synthetic_data_algorithm_comparison.Rmd")
```

Note: These studies may take several hours to complete and require significant computational resources.

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
* k_decay: Spring decay rate (typical range: 0.0001-0.1)
* cqq: Repulsion constant (typical range: 0.00001-0.1)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package, please cite:

Arhami and Rohani 2025 [doi:to be added]

## Contact

- Maintainer: Omid Arhami
- Email: omid.arhami@uga.edu
- GitHub: @omid-arhami
