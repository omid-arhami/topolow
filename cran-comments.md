# cran-comments.md — topolow 2.1.0

This is an update of the CRAN package `topolow` (currently 2.0.1).

## Test environments

* Local: Windows 11, R 4.5.3 — `R CMD check --as-cran`
* win-builder: R-devel and R-release
* macOS builder (R-release)
* R-hub: linux, macos, windows, clang-asan

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

## Notes for the reviewer

### This is the first release of this package that contains compiled code

Versions up to and including 2.0.1 were `NeedsCompilation: no`. In 2.1.0 the
core optimization loop has been moved from R to C++ via Rcpp, so this release is
`NeedsCompilation: yes`.

The backend is Rcpp-only. It calls no BLAS or LAPACK routine and links no
external library, so there is no `src/Makevars` or `src/Makevars.win` on any
platform. Native routines are registered with `R_registerRoutines()` and symbol
search is disabled with `R_useDynamicSymbols(dll, FALSE)`. All sources are ASCII
and compile without warnings under `-Wall -pedantic`.

### The C++ code uses a random number generator that `set.seed()` does not control

Writing R Extensions asks packages to use R's RNG so that `set.seed()` governs
results. This package does that for everything that determines the answer, but
one narrow use of randomness is intentionally kept outside R's stream, and we
would rather flag it than have it found.

The algorithm updates point positions in place (Gauss-Seidel) while sweeping
over all N(N-1)/2 point pairs. On every iteration the C++ code re-permutes the
*order* in which those pairs are visited, using `std::mt19937` seeded from
`std::random_device`.

The shuffle changes only the visiting order. It does not change which pairs are
considered, their observed dissimilarities, their weights, the spring and
repulsion physics, the cooling schedule, or the convergence test — all of which
are deterministic. Because updates are applied in place, a fixed sweep order
would systematically favour points early in that order, which are always moved
against stale neighbour positions; reshuffling each iteration removes that bias
and supplies the perturbation that lets the configuration escape local optima.
This is the layer at which the stochasticity is the method, and it is the same
behaviour the pure-R implementation in 2.0.1 had.

`set.seed()` does control the starting configuration, which is drawn with
`stats::runif()` in R before the C++ call. Run-to-run variation in the final
coordinates is expected and documented: `?euclidean_embedding` has a new
`Reproducibility` section stating exactly what is and is not seeded, and
describing how to obtain a single canonical map (run several times, compare with
the package's own diagnostics, keep the lowest holdout error). `NEWS.md` says
the same. We are happy to change this if you would prefer R's RNG here, but it
would alter published results, so we did not make the change unilaterally.

### `\dontrun{}` has been removed

CRAN asked us to remove `\dontrun{}` wrappers at the 1.0.0 review. Nine had crept
back in by 2.1.0. All nine are gone: the examples that were slow are now
`\donttest{}` and actually run, and the ones that referenced files or objects
that did not exist have been rewritten so they execute. `grep -rn "dontrun" man/`
returns nothing. The `Euclidify()` help page previously had six near-duplicate
example blocks; these are consolidated to three, each with small
`n_initial_samples` and `n_adaptive_samples`, to keep the page inside the
per-example time budget.

### Writing to the file system

No function writes outside `tempdir()` unless the user supplies `output_dir`
explicitly. Every example and test that writes a file writes to `tempdir()`.

## Downstream dependencies

There are no reverse dependencies on CRAN.
