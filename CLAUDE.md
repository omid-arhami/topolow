# CLAUDE.md — topolow

Guidance for Claude Code when working in this repository.

## 1. What this project is

`topolow` is a **CRAN-published scientific R package** (current source version 2.1.0;
the last recorded CRAN submission in `CRAN-SUBMISSION` was 2.0.0). It implements the
Topolow algorithm: force-directed Euclidean embedding of pairwise dissimilarity data
that is sparse, non-metric, and/or censored (thresholded `<` / `>` values). Its flagship
application is antigenic cartography, but the API is deliberately domain-neutral.

The package backs two publications:
- Arhami & Rohani (2025a), *Bioinformatics*, <doi:10.1093/bioinformatics/btaf372>
- Arhami & Rohani (2025b), arXiv:2508.01733 (mathematical properties)

**Three consequences govern every change:**
1. **Published results must stay reproducible.** Changes to the numerical core alter
   published figures and users' maps. Treat `R/core.R` and `src/optimization.cpp` as
   high-risk surfaces (see §7).
2. **CRAN policy is a hard constraint**, not a preference (see §6).
3. **The API is public.** Renaming or removing an exported function or argument is a
   breaking change requiring a deprecation cycle and a `NEWS.md` entry (see §5).

## 2. Repository map

| Path | Contents |
|---|---|
| `R/core.R` | `euclidean_embedding()` (main algorithm), `Euclidify()` (wizard), `create_topolow_map()` (deprecated), `print.topolow`, `summary.topolow` |
| `R/adaptive_sampling.R` | Parameter optimization: `initial_parameter_optimization()` (Latin hypercube), `run_adaptive_sampling()`, `adaptive_MC_sampling()`, `likelihood_function()`, profile likelihood and sensitivity analysis |
| `R/data_preprocessing.R` | Input conversion: `list_to_matrix()`, `table_to_matrix()`, `titers_list_to_matrix()`, `process_antigenic_data()`, `clean_data()`, `prune_sparse_matrix()` |
| `R/visualization.R` | Config constructors (`new_*_config()`), `plot_temporal_mapping()`, `plot_cluster_mapping()`, `plot_3d_mapping()`, `save_plot()`, network and scatter plots |
| `R/diagnostics.R` | MCMC/convergence diagnostics, `analyze_network_structure()`, performance traces |
| `R/euclidify_diagnostics.R` | Diagnostic plots and report for `Euclidify()` output |
| `R/error_metrics.R`, `R/utils.R` | Error comparison, CV folds, connectivity checks, subsampling |
| `R/globals.R` | `utils::globalVariables()` — add NSE column names here to silence R CMD check NOTEs |
| `R/topolow-package.R` | Package doc, `.onLoad()`, package-level `@importFrom` block |
| `src/optimization.cpp` | Rcpp/RcppArmadillo exact O(N^2) optimizer backend |
| `src/RcppExports.cpp`, `R/RcppExports.R` | **Generated** by `Rcpp::compileAttributes()` — never hand-edit |
| `NAMESPACE`, `man/*.Rd` | **Generated** by roxygen2 — never hand-edit |
| `tests/testthat/` | 15 test files, ~2750 lines, testthat edition 3 |
| `data/`, `data-raw/` | Bundled `.rda` datasets and the scripts that build them |
| `inst/examples/`, `vignettes/` | Long-form analyses (`.Rmd`); `inst/examples` is `.Rbuildignore`d |
| `github/workflows/` | Coverage workflow — see §9; this path is **not** `.github/`, so CI is inert |

**Canonical user workflow:** raw assay data → `list_to_matrix()` / `titers_list_to_matrix()`
→ (optional `prune_sparse_matrix()`, `subsample_dissimilarity_matrix()`) →
`Euclidify()` *or* `initial_parameter_optimization()` + `run_adaptive_sampling()` +
`euclidean_embedding()` → `plot_*_mapping()` / diagnostics.

## 3. Commands

R 4.5.3 is at `C:/Program Files/R/R-4.5.3/bin/x64/`. `devtools`, `roxygen2`, `testthat`,
`covr`, `Rcpp`, `RcppArmadillo`, and `rcmdcheck` are installed.

```r
devtools::document()        # regenerate NAMESPACE + man/ after roxygen edits
Rcpp::compileAttributes()   # regenerate RcppExports after editing src/*.cpp
devtools::load_all()        # recompile C++ and load for interactive work
devtools::test()            # full test suite
testthat::test_file("tests/testthat/test-core.R")  # single file — prefer while iterating
devtools::check()           # full R CMD check (SLOW — see below)
covr::package_coverage()    # coverage
```

**Runtime discipline (per the user's global instructions): do not run long jobs
yourself.** `devtools::check()`, `covr::package_coverage()`, vignette builds, and any
`Euclidify()` or `run_adaptive_sampling()` call on real data take minutes to hours.
Write the command out, hand it to the user, and wait for their output. Targeted
`testthat::test_file()` runs on small matrices are fine to run directly.

After editing `src/*.cpp` you **must** recompile (`devtools::load_all()` or
`devtools::install()`) before tests reflect the change — stale `.o` / `.so` artifacts sit
in `src/` and are gitignored.

## 4. Code conventions (match what is already there — do not restyle)

- **Roxygen2 markdown** (`Roxygen: list(markdown = TRUE)`), 2-space indent, `<-` for
  assignment, `snake_case` for functions and arguments.
- **File header**: each `R/*.R` opens with `# Copyright (c) <year> Omid Arhami <email>`
  followed by `# R/<filename>`. Keep this on new files. Existing files are inconsistent
  between `omid.arhami@uga.edu` and `o.arhami@gmail.com` — match the file you are
  editing rather than "fixing" it.
- **Every function gets a roxygen block** with `@description`, `@param` for *all*
  arguments, `@return` describing class and structure, and either `@export` or
  `@keywords internal`. CRAN requires `@return` on every exported function.
- **Argument naming**: `dissimilarity_matrix` (not `distance_matrix` — that spelling is
  the deprecated one), `ndim`, `output_dir`, `max_cores`, `verbose`, `save_plot`.
  Follow the existing vocabulary exactly when adding arguments.
- **Input validation first**: functions open with a block of
  `if (...) stop("<argument> must be <condition>")` in plain lowercase, before any
  computation. Tests assert on these exact message strings (`tests/testthat/test-core.R`),
  so changing a message is a test-visible change.
- **Lifecycle**: use `lifecycle::deprecate_warn()` plus a
  `lifecycle::badge("deprecated")` inline call in the roxygen `@description`. Stable
  functions carry `lifecycle::badge("stable")`.
- **Imports**: add `@importFrom` in the relevant file, then re-document. Never edit
  `NAMESPACE` by hand. Do not add a new package to `Imports` without justifying it to
  the user first — each dependency is a CRAN maintenance liability.
- **NSE columns** used by dplyr/ggplot2/data.table go in `R/globals.R`, or use
  `rlang::.data$`, to keep R CMD check clean.

## 5. Changing the public API

Before renaming or removing anything exported, check `NAMESPACE`, then:
1. Keep the old function/argument working via `lifecycle::deprecate_warn()`.
2. Add a `NEWS.md` entry under a version heading with a **Migration** line showing
   old → new. Follow the existing format (see the 2.0.0 section of `NEWS.md`).
3. Add a test in `tests/testthat/test-deprecated.R`.
4. Bump `Version:` in `DESCRIPTION` per semver — breaking changes require a major bump.

`create_topolow_map()` is the reference example of the full deprecation pattern.

## 6. CRAN constraints (non-negotiable)

- **Never write to the user's filesystem by default.** Functions that produce files take
  an explicit `output_dir` and error if it is missing (`Euclidify()` does this).
  Examples must write only to `tempdir()`.
- **Examples must run fast or be wrapped.** Use `\donttest{}` for slow-but-valid
  examples; `\dontrun{}` only when the example genuinely cannot run.
- **Parallelism**: default to at most `future::availableCores() - 1`; CRAN checks run
  with 2 cores. Never size a cluster from an unconditional `detectCores()`.
- **No `print()` / `cat()` outside explicit `verbose` branches.** Use `message()` or
  `warning()` for diagnostics.
- Add new technical terms to `inst/WORDLIST` when R CMD check flags spelling.
- The target is **0 ERRORs, 0 WARNINGs, 0 NOTEs** from `devtools::check(cran = TRUE)`.

## 7. High-risk surfaces — ask before touching

Stop and ask the user before changing any of the following, and state what you believe
the numerical consequence will be:

- The force model, convergence criterion, or cooling schedule in `src/optimization.cpp`
  or `euclidean_embedding()`.
- The likelihood or weighting used by adaptive sampling (`likelihood_function()`,
  `calculate_weighted_marginals()`, temperature scaling).
- Threshold (`<` / `>`) handling in `vectorized_process_distance_matrix()` — the sign
  conventions are subtle and load-bearing.
- Default parameter values or ranges in `Euclidify()`.
- Anything under `data/` or `data-raw/` (bundled datasets used by examples and tests).
- Regenerating all of `man/` when the actual change touches one function.

For these, "it looks like a bug" is not sufficient — reproduce it with a minimal
matrix, show the before/after numbers, and let the user decide.

## 8. Working agreements

- **Read before editing.** Inspect the function, its roxygen block, its tests, and its
  callers (grep the function name across `R/`, `tests/`, `vignettes/`, `inst/`).
- **Do not guess.** If a parameter's intended semantics, a formula, or a scientific
  convention is not determinable from the repository, ask. Do not invent numerical
  constants, citations, or algorithm details.
- **Minimal, surgical diffs.** No opportunistic refactors, reformatting, or renaming
  outside the requested scope.
- **Tests accompany behavior changes.** Add or update the matching
  `tests/testthat/test-<area>.R`; keep fixtures small (3–10 points) so they run fast.
- **Report honestly.** State exactly which commands were run and their real output. If
  verification could not be run (e.g. a long check), say so and name the residual risk.
- Documentation edits go in the **roxygen block**, then `devtools::document()` — never
  directly in `man/*.Rd`.

## 9. Known repository state (recorded at first setup on this machine, 2026-09-02)

- The working tree shows ~124 modified files that are **file-mode changes only**
  (`100755` → `100644`) from copying the repo to Windows, plus CRLF/LF warnings. Content
  is unchanged. Do not bundle these into a feature commit; ask the user how to handle them.
- **CI is inert**: the coverage workflow lives at `github/workflows/test-coverage.yaml`,
  not `.github/workflows/`. GitHub Actions does not read that path, and there is no
  R CMD check workflow at all. Flag this if the user wants CI; do not silently move it.
- `src/*.o` and `src/topolow.so` are build artifacts, gitignored, and may be stale.
- `.Rapp.history` and `.Rhistory` are tracked-but-empty noise files.

## 10. Git protocol

Commit with the user's explicit identity override — **always this exact form**:

```bash
git -c user.name="Omid Arhami" -c user.email="o.arhami@gmail.com" commit -m "<message>"
```

- Commit messages follow the existing history: a short descriptive summary of the
  user-visible change. No Conventional Commits prefixes, no trailers, no co-author
  lines. Examples from the log: `Return velocity vectors in ndim`,
  `Fix Windows BLAS linking: add src/Makevars.win with LAPACK/BLAS/FLIBS flags`.
- **Only commit when asked.** Do not push, tag, or open PRs without an explicit request.
- Stage precisely with `git add <paths>` — never `git add -A` while the mode-change
  churn described in §9 is outstanding.
- The repository is public on GitHub (`omid-arhami/topolow`). Anything committed is
  published; do not commit data, credentials, or draft manuscript text without asking.
