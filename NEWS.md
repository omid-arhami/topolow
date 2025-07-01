# topolow 0.3.2

*  Initial release to CRAN (revised per CRAN reviewr's instructions).
*  Introduces the Topolow algorithm, a physics-inspired method for antigenic cartography.
*  Provides robust mapping and complete positioning of all antigens, even with highly sparse datasets (>95% missing values).
*  Implements automatic, likelihood-based estimation to determine the optimal dimensionality of the antigenic map.
*  Includes functionality to calculate "antigenic velocity" vectors to quantify the rate and direction of antigenic drift.
*  Features tools for handling and processing cross-reactivity and binding affinity assay data, including those with thresholded values.
*  Demonstrates improved prediction accuracy and run-to-run stability compared to traditional MDS methods.

# topolow 1.0.0 (2025-07-01)

* All exported methods now include `\value` documentation describing the output's class, structure, and meaning.
* Examples for unexported functions have been omitted, and `\dontrun{}` wrappers have been removed5. Slower examples are now wrapped in `\donttest{}` as appropriate.
* Functions no longer write to user directories by default. Functions where writing a file is the main purpose now require the user to specify an output directory.
* The complex distributed processing functionality has been removed, as it was not essential for typical use cases.
* The link to our paper and citation information have been updated.