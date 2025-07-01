## Resubmission: topolow 1.0.0

Dear CRAN maintainers,

This is a revised submission of the `topolow` package. Thank you for your careful evaluation and feedback, which we have addressed in this version.

### Summary of Changes

In response to the previous review, we have made the following changes:

* **Documentation:** All exported methods now include `\value` documentation describing the output's class, structure, and meaning.
* **Examples:** Examples for unexported functions have been omitted, and `\dontrun{}` wrappers have been removed5. Slower examples are now wrapped in `\donttest{}` as appropriate.
* **File Writing:** Functions no longer write to user directories by default. Functions where writing a file is the main purpose now require the user to specify an output directory.
* **Functionality:** The complex distributed processing functionality has been removed, as it was not essential for typical use cases.
* **Citation:** The link to our paper and citation information have been updated.

Due to significant changes and extensive testing, we increased the version number.

* The theories and methods are described in our paper, Arhami and Rohani (2025) https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaf372/8173949.

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