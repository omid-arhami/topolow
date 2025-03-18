# install_github.R
# Installation script for topolow package
# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install required dependencies
required_packages <- c(
  "ggplot2", "dplyr", "data.table", "reshape2", "plotly", "rgl", 
  "Racmacs", "parallel", "coda", "MASS", "vegan", "igraph", 
  "lhs", "umap", "gridExtra", "scales"
)

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Install topolow from GitHub
devtools::install_github("omid-arhami/topolow")

# Test installation
library(topolow)
?topolow  # Should open package documentation
