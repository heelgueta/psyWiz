# Getting psyWiz onto CRAN

A plain-language guide for when the time comes. Nothing here has to happen
before the jamovi module is usable; CRAN is a separate distribution channel
for R users.

## What CRAN is and what it asks for

CRAN is the official R package archive (`install.packages("psyWiz")` would
fetch from it). Volunteers review each submission, mostly through automated
checks. They do **not** judge whether a package is novel or statistically
sound. They require that it:

1. passes `R CMD check --as-cran` with no errors, no warnings and as few
   notes as possible, on Linux, macOS and Windows, with the current and the
   development version of R;
2. documents every exported function and data set;
3. has examples and tests that run quickly (a few seconds per example, the
   whole check within about ten minutes) and never use more than two cores;
4. depends only on packages that are on CRAN (or Bioconductor);
5. has a clear licence, a maintainer with a working e-mail address, and a
   `Description` that says what the package does and cites the methods;
6. does not write to the user's home directory, change global options, or
   install anything.

After acceptance the maintainer must fix problems CRAN reports (for example
when a new lavaan release breaks a test), usually within two to four weeks, or
the package is archived. That ongoing duty is the real cost.

## Where psyWiz stands (checked 2026-11-14)

`R CMD check --as-cran --no-manual` was run with jamovi's R 4.6.0 on Linux:
no errors, tests pass in about a minute. What it reported, and what else is
known:

| Item | Status | What to do |
|---|---|---|
| Help pages for exported functions and `wizdata` | **Missing** (check warning) | Run `roxygen2::roxygenise(roclets = "rd")` on a machine with roxygen2 and commit `man/`. The roxygen comments are already in the source. The jamovi functions (`groupCfaOne()` etc.) get theirs from the generated `R/*.h.R`. |
| Non-ASCII characters in R code | Fixed | Greek letters in table notes are written as `\uXXXX` escapes. Keep it that way; yaml files may use real characters. |
| `LICENSE` file "not mentioned" | Fixed | The file stays in the repository for GitHub but is excluded from the R build (`.Rbuildignore`); `License: GPL-3` in DESCRIPTION is what CRAN reads. |
| Dependencies | OK | jmvcore, R6, lavaan and ggplot2 are all on CRAN. |
| Examples | Present, wrapped in `\donttest{}` | They take a few seconds each. CRAN accepts `\donttest` for slow examples but prefers at least one fast, unwrapped example per function. Add a tiny one (one factor, continuous, ML). |
| Parallel workers | OK | CRAN allows at most two cores. psyWiz caps itself at two workers when `_R_CHECK_LIMIT_CORES_` is set (as it is during CRAN checks), and the tests use `parallel = FALSE` or `ncores = 2`. |
| Test duration | About 60 s | Fine. Keep heavy tests behind `skip_on_cran()`. |
| `Description` field | Needs polish | CRAN wants method references in the form `Author (year) <doi:...>`, for example `Wu and Estabrook (2016) <doi:10.1007/s11336-016-9506-0>`, and software names in single quotes. |
| `data/wizdata.omv` | Excluded from the R build | It is only for jamovi. |
| `jamovi/`, `tools/`, `docs/`, `build/` | Excluded from the R build | Via `.Rbuildignore`. |
| Checks on macOS, Windows, R-devel | **Not done** | `.github/workflows/R-CMD-check.yaml` is written for this but has not run yet. |
| Spelling, URLs | Not checked | `devtools::spell_check()`, `urlchecker::url_check()`. |
| Package name | Not checked | Check with `available::available("psyWiz")` before submitting. Names differing only in case from an existing package are refused. |

## Step by step, when ready

Needs a normal R installation with development tools (not jamovi's R), or
just GitHub Actions for steps 2 to 4.

1. Decide the API is stable. Renaming arguments after release is painful.
2. `roxygen2::roxygenise(roclets = "rd")`, commit `man/`.
3. `devtools::check(cran = TRUE)` locally until clean.
4. Check other platforms: push and read the GitHub Actions results, then
   `devtools::check_win_devel()` (uploads to CRAN's Windows test server and
   e-mails the result).
5. Write `cran-comments.md`: the platforms checked, and one line explaining
   each remaining note ("New submission" is always there the first time).
6. Set `Version: 1.0.0` (or `0.2.0`; avoid development suffixes), update
   `NEWS.md` and the date.
7. `devtools::submit_cran()`, or build the tarball and upload it at
   <https://cran.r-project.org/submit.html>. Confirm via the e-mail link.
8. Expect an answer within a few days, often with small requests (wording of
   the Description, an example to unwrap). Fix, bump the version, resubmit.

## Things that commonly get first submissions rejected

- Title not in Title Case, or starting with "A package for".
- Description starting with the package name or "This package".
- No references for the implemented methods in the Description.
- Examples wrapped in `\dontrun{}` without need (use `\donttest{}`).
- `print()`/`cat()` used for messages that cannot be silenced (psyWiz only
  prints in print methods, which is fine).
- More than two cores used during checks.
- Missing `\value` sections in help pages (each exported function needs one;
  the roxygen comments have `@return`).

## jamovi module and CRAN package from one repository

This combination is established practice: `jmv` itself is on CRAN. Points to
keep in mind:

- The generated files `R/*.h.R` and `R/00jmv.R` must be committed, because
  CRAN builds from the R sources only.
- The jamovi library builds modules from the repository on its own
  infrastructure; CRAN status is irrelevant to it and vice versa.
- Keep one version number for both (DESCRIPTION and `jamovi/0000.yaml`).
