# Developing psyWiz

## What you need

- jamovi (2.4 or later). Its bundled R is enough; a separate R installation
  is optional.
- `jmvtools`, the jamovi module compiler, installed into a library of your own:

  ```r
  install.packages("jmvtools", lib = "~/jmvtools-lib",
                   repos = c("https://repo.jamovi.org", "https://cloud.r-project.org"))
  ```

- Python 3 (only to regenerate the group-analysis yaml files).

lavaan does not need installing: jamovi ships it inside its `jmv` module.

## Workflow

```sh
python3 tools/gen_yaml.py                 # only after editing that script
Rscript tools/build.R <jamovi home> <jmv R library>   # compile, build .jmo, install into jamovi
R CMD INSTALL --library=<test lib> .      # install the R package for testing
Rscript -e 'testthat::test_local()'       # or test_dir("tests/testthat", package = "psyWiz", load_package = "installed")
```

`tools/build.R` regenerates `R/*.h.R` and `R/00jmv.R` from `jamovi/*.yaml`,
writes `psyWiz_<version>.jmo` (git-ignored) and installs the module into the
running jamovi. On most installations plain `jmvtools::install()` does the
same; the script exists for the case below.

### Linux, jamovi as a flatpak (the main development setup)

There is no R outside the sandbox, the sandbox cannot see `/tmp`, and jmvtools
does not find the lavaan that jamovi bundles. Everything therefore runs through
`flatpak run` with an explicit library path:

```sh
JMV=/app/lib/jamovi/modules/jmv/R                 # lavaan lives here
LIBS=$HOME/jmvtools-lib:$HOME/.cache/psywiz/lib:$JMV
RS="flatpak run --env=R_LIBS=$LIBS --command=Rscript org.jamovi.jamovi"
RC="flatpak run --env=R_LIBS=$LIBS --command=R org.jamovi.jamovi"

$RS tools/build.R /app $JMV
mkdir -p ~/.cache/psywiz/lib
$RC CMD INSTALL --no-docs --library=$HOME/.cache/psywiz/lib .
$RS -e 'testthat::test_dir("tests/testthat", package = "psyWiz", load_package = "installed")'
```

Set `NOT_CRAN=true` (`--env=NOT_CRAN=true`) to include the parallel test.

### macOS and Windows

Untested so far. Point `R_HOME`/`PATH` at jamovi's R as described in the
[jamovi developer hub](https://dev.jamovi.org), then
`jmvtools::install()` should work; if it stops with "dependency 'lavaan' is
not available", use `tools/build.R` with the path of jamovi's
`Resources/modules/jmv/R` as second argument.

## Rules

- jamovi marks a module as incompatible (no menu appears) if it was built with
  a different R version than the one in the jamovi it is installed in. Build
  on the jamovi version you target.
- Do not edit `R/*.h.R`, `R/00jmv.R`, or the `groupCfa*.yaml` files by hand.
- An option name in a `.a.yaml` becomes an argument of the generated R
  function; renaming one is a breaking change for R users.
- In headless calls pass variable names as strings (`vars = c("x1", "x2")`).
- Statistics go in `R/engine-*.R` with a test in `tests/testthat/test-engine.R`
  that checks them against plain lavaan. See [DESIGN.md](DESIGN.md).
- User-visible changes go in `NEWS.md`; changes to what is computed also go in
  [METHODS.md](METHODS.md).
- `NAMESPACE` is written by hand. Help pages should be generated from the
  roxygen comments with `roxygen2::roxygenise(roclets = "rd")` (not done yet).

## Regenerating the example data

```sh
Rscript tools/make_examples.R     # writes data/wizdata.rda and data/wizdata.omv
```

## Reporting and fixing bugs

A useful report includes: jamovi and psyWiz versions, the analysis and options
used, and if possible a data set
(simulated is fine) that reproduces the problem. Most fixes start by
reproducing the problem in R with `invarianceWiz()` or `cfaWiz()` and
`keep.fits = TRUE`, then adding a failing test.
