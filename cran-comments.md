## CRAN Submission: topolow 0.3.2

Dear CRAN maintainers,

This is the initial submission of the `topolow` package to CRAN.

The `topolow` package provides a novel, physics-inspired algorithm for antigenic cartography. [cite_start]It is designed to create complete and accurate antigenic maps from sparse and complex immunological assay data. [cite_start]The package includes tools for automatic dimensionality detection and the calculation of "antigenic velocity" to track pathogen evolution. [cite_start]The theories and methods are described in our paper, Arhami and Rohani (2025) https://www.biorxiv.org/content/10.1101/2025.02.09.637307v1.

## Test Environments
* local: macOS Sonoma 14.5, R 4.4.1
* devtools::check_win_devel()
* R-hub: Windows Server 2022 (R-devel), Ubuntu Linux 20.04 (R-release), Fedora Linux (R-devel)

## R CMD check results

There were no ERRORs or WARNINGs. There was one NOTE regarding file timestamps:

* `checking for future file timestamps ... NOTE unable to verify current time`

This appears to be a system-specific issue related to my local macOS environment and can be safely ignored.

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.

Thank you for your time and consideration.

Sincerely,

Omid Arhami