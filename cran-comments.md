# cran-comments.md — topolow 2.2.0

This is an update of the CRAN package `topolow` (currently 2.1.0, published
2026-09-14).

## Test environments

* Local: Windows 11, R 4.5.3 — `R CMD check --as-cran`
* win-builder: R-devel and R-release — PENDING, run before submitting
* macOS builder (R-release) — PENDING, run before submitting
* R-hub: linux, macos, windows — PENDING, run before submitting

(Delete the PENDING markers once those checks have come back clean. Do not
submit with them still in the file.)

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

## What is in this release

The release is additive. It adds one feature — resampling a subset of the points
during parameter optimization and during map building — together with three new
exported functions, two new S3 methods, and new optional arguments on five
existing functions. Every new argument defaults to `NULL` or `FALSE`, no
argument was renamed or removed, and no default value was changed. Code written
against 2.1.0 runs unchanged.

No dependency was added. The new code uses only packages already in `Imports`,
plus `MASS` from `Suggests` behind a `requireNamespace()` guard for one optional
plotting method.

## Notes for the reviewer

### A documented behaviour change in `Euclidify()`

`Euclidify()` accepted an `opt_subsample` argument but did not pass it on to its
adaptive refinement stage, which therefore embedded the full matrix even when
subsampling had been requested. The help page said subsampling applied to "both
initial and adaptive", so the documentation was right and the code was wrong.
The argument is now forwarded. Existing calls that set `opt_subsample` will find
the adaptive stage faster, and its likelihood estimates correspondingly noisier,
than in 2.1.0. This is called out under `## Bug Fixes` in `NEWS.md`.

### Writing to the file system

Unchanged and still the case: no function writes outside `tempdir()` unless the
user supplies `output_dir` explicitly. The new `bootstrap_embedding()` writes
nothing by default; `write_files = TRUE` requires `output_dir` and errors without
it. Every example and test that writes a file writes to `tempdir()`.

### Examples

All examples run; none is wrapped in `\dontrun{}`. `grep -rn "dontrun" man/`
returns nothing. The new functions' examples use a 9-point matrix and a handful
of replicates and complete in well under a second each; the one longer example,
on `bootstrap_blobs()`, is wrapped in `\donttest{}`. `--run-donttest` completes
in 18 seconds for the whole package.

## Downstream dependencies

There are no reverse dependencies on CRAN.
