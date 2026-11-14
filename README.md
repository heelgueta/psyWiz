# psyWiz

Psychometric wizardry for confirmatory factor analysis, as an **R package** and a
**jamovi module** in one repository. Built on [lavaan](https://lavaan.ugent.be).

It does two things that are tedious to do by hand:

1. **Measurement invariance across groups, step by step.** Give it a factor
   model and a grouping variable. It fits the whole ladder (configural,
   loadings, intercepts, optionally residuals; for ordinal items the
   threshold-based ladder of Wu & Estabrook, 2016), compares each step with the
   previous one, gives a verdict against cutoffs you can edit, and lists the
   constraints that fit worst when a step fails.
2. **One CFA under several estimators, side by side.** ML, MLR, MLM, ULS,
   ULSMV, DWLS, WLSMV, WLS: fit indices, loadings, factor correlations and
   reliability in adjacent columns, so you can see what changes when Likert
   items are treated as continuous or as ordinal.

Independent models are fitted on several processor cores when that saves time.

> Status: version 0.2, a rewrite. It replaces the unfinished 2024 prototype
> (`cfawiz()`), which is kept on the `psywiz_old` branch. Not yet in the jamovi
> library. Bug reports are welcome: <https://github.com/heelgueta/psywiz/issues>

## In jamovi

The module adds a **psyWiz** menu with three analyses:

| Analysis | What it does |
|---|---|
| One-Factor Group CFA | Invariance of a single factor across groups |
| Correlated-Factors Group CFA | Invariance of several correlated factors across groups |
| Multi-Estimator CFA | One model, several estimators, side by side |

The group analyses can also draw the model: a fixed-layout diagram with one
panel per group, where the width of each arrow shows the size of the loading.
Colour is off by default; one checkbox highlights flagged values and weak
loadings.

Drag items into the factor boxes, set the options, then tick **Run analysis**.
Nothing is fitted before that, so jamovi does not restart a slow analysis every
time you change an option. Changing the decision cutoffs or the output options
afterwards does not refit anything.

Install: build the `.jmo` from this repository (see
[docs/DEVELOPING.md](docs/DEVELOPING.md)) and sideload it in jamovi
(Modules > Sideload). The example data set *Trust and stress in three countries*
is in jamovi's data library after installation.

## In R

```r
# install.packages("remotes")
remotes::install_github("heelgueta/psywiz")
library(psyWiz)
data(wizdata)

model <- list(Trust  = paste0("trust", 1:4),
              Stress = paste0("stress", 1:4))

# measurement invariance across countries
inv <- invarianceWiz(wizdata, model, group = "country")
inv
#> Model fit
#>                          chisq        df    pvalue       cfi       tli     rmsea      srmr
#> Configural              54.368        57     0.574     1.000     1.001     0.000     0.028
#> Thresholds              91.653        89     0.403     0.999     0.999     0.009     0.028
#> Loadings (metric)       94.431       101     0.665     1.000     1.001     0.000     0.029
#> Intercepts (scalar)    201.167       113     <.001     0.982     0.987     0.044     0.030
#>
#> Each step against the previous one
#>                        d.chisq      d.df         p     d.cfi   d.rmsea    d.srmr   verdict
#> Thresholds              38.942        32     0.186    -0.001     0.009     0.000  Supported
#> Loadings (metric)        5.985        12     0.917     0.001    -0.009     0.001  Supported
#> Intercepts (scalar)     84.130        12    !<.001   !-0.018    !0.044     0.001  Not supported
#>
#> Worst-fitting constraints, Intercepts (scalar) (score tests)
#>                       chisq         p
#> stress3 [Brazil]     80.202     <.001
#> ...

# stricter cutoffs, no refitting
invarianceVerdict(inv, criteria = c("cfi", "chisq"), cutoffs = list(cfi = 0.005))

# the same CFA under several estimators ("MULTI" = all eight)
multi <- cfaWiz(wizdata, model, estimators = c("MLR", "ULSMV", "WLSMV"))
print(multi, colored = TRUE)   # colour ramp for fit and loadings, off by default

# diagrams: one panel per group / per estimator
plot(inv)
plot(multi, highlight = TRUE)
```

Both functions return plain data frames (`$fit`, `$comparison`, `$culprits`,
`$loadings`, ...) plus the lavaan syntax that was fitted (`$syntax`). Use
`keep.fits = TRUE` to also get the lavaan objects.

## What it does not do (yet)

- Bifactor and higher-order models, cross-loadings, correlated residuals, or
  free-form lavaan syntax.
- Automatic partial invariance. psyWiz shows which constraints fit worst and
  leaves the decision to you.
- Multi-estimator comparison of the invariance ladder (the two analyses are
  separate for now).
- Latent mean comparisons after invariance has been established.

See [NEWS.md](NEWS.md) for changes and
[docs/METHODS.md](docs/METHODS.md) for exactly which models are fitted and how
the verdicts are reached.

## Documentation

- [docs/METHODS.md](docs/METHODS.md): the statistics, with references.
- [docs/DESIGN.md](docs/DESIGN.md): how the code is organised and why.
- [docs/DEVELOPING.md](docs/DEVELOPING.md): building, testing, contributing.
- [docs/CRAN.md](docs/CRAN.md): what a CRAN submission needs and where psyWiz stands.

## Licence

GPL-3. See [LICENSE](LICENSE).
