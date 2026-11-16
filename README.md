# psyWiz

A **jamovi module** for confirmatory factor analysis and scale psychometrics,
built on [lavaan](https://lavaan.ugent.be). It fills gaps in what jamovi
offers out of the box: CFA that treats ordinal items as ordinal, measurement
invariance step by step with a verdict, and figures in which lengths mean
something.

You drag items into factor boxes. There is no model syntax to write.

> Status: version 0.4, in development. Tested headless on Linux; not yet in
> the jamovi library. Bug reports: <https://github.com/heelgueta/psywiz/issues>

## The menu

| Group | Analysis | What you get by default | On request |
|---|---|---|---|
| Factor Models | **CFA** | Fit, loadings, factor correlations, alpha / omega / AVE, comb plot | Residual correlations; second-order factor; path diagram; factor correlations from composite scores and factor scores; largest modification indices; colour summary |
| | **Group CFA (Invariance)** | The invariance ladder across groups with a verdict per step, and the constraints that fit worst at a failed step | Partial invariance search; robustness of the verdict; loadings by group; figure |
| | **Bifactor CFA** | General plus specific factors: fit, omega hierarchical, explained common variance, antenna plot | Facing combs; path diagram |
| | **Multi-Estimator CFA** | One model under up to eight estimators, side by side | |
| | **CFA vs ESEM** | CFA against a model with cross-loadings; each factor correlation by scoring method | Brush plot |
| Scales | **Scale Check** | Reliability (alpha, lambda-2, omega, AVE) and item statistics | Tau-equivalence and parallelism; Guttman / Mokken scalability and scalogram; comb plot, path diagram |
| | **Item Order (Exploratory)** | For yes/no items: which form a line, share a rung, branch or stand apart; order diagram | Answer patterns |
| | **Composite Check** | For formative indicators: do they act through one composite? | |

The first three rows are the core. Every analysis takes one factor or several
correlated factors, and the simpler ones are deliberately contained in the
more complex ones: the first model of a Group CFA is the plain CFA fitted in
each group.

How it behaves:

- Set everything up, then tick **Run analysis**. Nothing is fitted before
  that, so jamovi does not restart a slow model every time you touch an option.
  Output options and decision cutoffs never refit.
- Items are treated as ordinal (WLSMV, polychoric correlations) or continuous
  (MLR) according to their measure type in jamovi; you can override it.
- Colour is off by default; one **Colour** checkbox marks flagged values.
- Figures are black, white and grey: comb plots (every length is a share of
  variance), path diagrams, a scalogram, antenna plots for bifactor models
  and a brush plot for ESEM.

## Install

Until the module is in the jamovi library: build the `.jmo` from this
repository (see [docs/DEVELOPING.md](docs/DEVELOPING.md)) and sideload it
(Modules > Sideload). The example data set *Trust and stress in three
countries* then appears in jamovi's data library.

## In R

A jamovi module is an R package, so everything is also callable from R. This
side is not the focus yet and its interface may still change; see
[docs/ROADMAP.md](docs/ROADMAP.md).

```r
# install.packages("remotes")
remotes::install_github("heelgueta/psywiz")
library(psyWiz)
data(wizdata)

model <- list(Trust  = paste0("trust", 1:4),
              Stress = paste0("stress", 1:4))

# a plain CFA: fit, loadings, reliability; plot() draws the comb plot
cfa <- factorWiz(wizdata, model)
plot(cfa)

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

# the intercepts step failed: is it a few items? (frees stress3, then it holds)
invariancePartial(inv)

# the same CFA under several estimators ("MULTI" = all eight)
multi <- cfaWiz(wizdata, model, estimators = c("MLR", "ULSMV", "WLSMV"))
print(multi, colored = TRUE)   # colour ramp for fit and loadings, off by default

# figures: one panel per group / per estimator
plot(inv)                       # path diagram
plot(inv, type = "comb")        # comb plot: every length is shared variance

# one scale: reliability, tau-equivalence, Guttman / Mokken scalability
scaleWiz(wizdata, paste0("stress", 1:4), name = "Stress")

# CFA against ESEM, and the factor correlation under five scorings
esemWiz(wizdata, model)

# does the invariance verdict survive other estimators and decision rules?
invarianceMultiverse(wizdata, model, group = "country")

# formative indicators: is one composite enough? (needs two or more outcomes)
# compositeWiz(data, indicators = c("x1", "x2", "x3"), outcomes = c("y1", "y2"))
```

All of them return plain data frames (`$fit`, `$comparison`, `$culprits`,
`$loadings`, ...) plus the lavaan syntax that was fitted (`$syntax`). Use
`keep.fits = TRUE` to also get the lavaan objects.

## What it does not do (yet)

- User-specified cross-loadings, free-form lavaan syntax.
- Residual correlations, second-order and bifactor models across groups
  (they are single-group for now).
- Latent mean comparisons after invariance has been established.
- Release of thresholds in the partial invariance search.

## Documentation

- [NEWS.md](NEWS.md): what changed in each version.
- [docs/ROADMAP.md](docs/ROADMAP.md): where the project came from, where it stands, what is next.
- [docs/DATA.md](docs/DATA.md): the example data, block by block, and what to try on each.
- [docs/METHODS.md](docs/METHODS.md): exactly which models are fitted and how verdicts are reached, with references.
- [docs/DESIGN.md](docs/DESIGN.md): how the code is organised and why.
- [docs/IDEAS.md](docs/IDEAS.md): the less usual ideas, what is new and what is not.
- [docs/DEVELOPING.md](docs/DEVELOPING.md): building, testing, contributing.
- [docs/CRAN.md](docs/CRAN.md): what a CRAN submission would need (later).

## Licence

GPL-3. See [LICENSE](LICENSE).
