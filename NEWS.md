# psyWiz 1.0.0

The module is now the focus; the R interface stays but is not advertised.
The menu was reorganised around a simple core. See `docs/ROADMAP.md`.

## Latent mean comparisons

- New feature in **Group CFA (Invariance)** (`groupCfa`) and R API (`invarianceMeans()`):
  compares factor means across groups relative to a baseline group (defaults to the
  first group) once scalar or partial scalar invariance is supported.
- Reports raw mean differences ($\Delta \kappa$), standard errors derived from the
  parameter covariance matrix via the Delta method, Wald $z$ tests, $p$-values, 95%
  confidence intervals, and standardized effect sizes: Latent Cohen's $d$ using the
  reference group factor standard deviation and pooled factor standard deviation
  (Hancock, 2001).
- In R, `invarianceMeans(x, reference = NULL)` operates on `psw_invariance` and
  `psw_partial` objects, or can be run directly inside `invarianceWiz(..., means = TRUE)`.
- If scalar invariance is not supported across groups, an informative message explains
  why means cannot be compared.
- New **Latent Means Comparison Plot** (`plot(m)` or `plot(x, type = "means")` in R;
  option `latentMeansPlot` in jamovi): an intuitive forest/comparison plot displaying
  mean differences and confidence intervals by factor with reference group baseline.
- Fixed `geom_label()` deprecation warning with `linewidth`.

## Visual Psychometrics & UI Overhaul

- **Clean Publication-Ready Figures:** Removed instructional captions baked inside the plots;
  reading guides now appear cleanly *above* each novel plot in jamovi and can be toggled off
  via **Show plot guide** (on by default).
- **Model Fit on Plots:** All figures (Comb, Path, Antenna, Facing Combs, Brush) support
  optional model fit printing (`below` or `above`), dynamically wrapping to plot width and
  showing the estimator (e.g. `Estimator: MLR`, `WLSMV`), $\chi^2(df)$, $p$, CFI, TLI,
  RMSEA [90% CI], and SRMR.
- **Native Red Table Flagging (No Colour Heatmap):** Removed the separate HTML "Colour Summary"
  table and palette selector; replaced with a single **Flag weak / misfitting values** toggle
  (off by default) that highlights poor fit indices, weak loadings, and ESEM cross-loadings in red
  directly inside the standard tables.
- **Intuitive Menu & Consistent UI Layout:** Reordered the dropdown menu into a clean 6-analysis
  progression (`Scale Check` $\rightarrow$ `CFA` $\rightarrow$ `Measurement Invariance` $\rightarrow$
  `Bifactor CFA` $\rightarrow$ `CFA vs. ESEM` $\rightarrow$ `Multi-Estimator CFA`) with uniform
  within-analysis sections (`Estimation`, `Tables`, `Plots`, analysis extensions, and cutoffs at the bottom).
- **Visual Psychometrics Guide (`docs/PLOTS.md`):** 9 realistic fitted example figures illustrating
  good, poor, and mixed empirical conditions across Comb, Antenna, and Brush plots.
- **Deferred for Future Releases:** Exploratory `Item Order` and `Composite Check` archived out
  of the initial 1.0.0 release.

## Menu

- New **CFA** analysis (`plainCfa` in jamovi, `factorWiz()` in R): one factor
  or several correlated factors in one group. Fit, loadings, factor
  correlations, alpha / omega / AVE and a comb plot by default. On request:
  each factor correlation as obtained from composite scores, corrected
  composites and factor scores; and the largest modification indices
  (residual correlations and cross-loadings).
- **Group CFA (Invariance)** (`groupCfa`) replaces One-Factor Group CFA and
  Correlated-Factors Group CFA. It starts with one factor box; add more for
  correlated factors. *Breaking for R users of the jamovi functions:*
  `groupCfaOne()` and `groupCfaCorr()` are gone; analyses saved in `.omv`
  files with the old names will not reopen.
- **Scale Check** shows reliability and item statistics only. The
  tau-equivalence ladder and the Guttman / Mokken statistics are under
  "More", off by default.
- Two menu groups: *Factor Models* (CFA, Group CFA, Multi-Estimator CFA, CFA
  vs ESEM) and *Scales* (Scale Check, Composite Check).

## Partial invariance

- `invariancePartial()` and the option "Partial invariance" in Group CFA:
  when a step is not supported, the item with the largest score test is left
  free in every group and the step is judged again, keeping at least two
  constrained items per factor. A freed loading frees the intercept, a freed
  intercept frees the residual variance.
- First written by another coding session on 2026-11-15; reviewed and
  rewritten on 2026-11-16. Changed in the review: thresholds are no longer
  released (the model with free thresholds was not identified); the jamovi
  result is cached instead of being refitted on every option change; the
  data are no longer stored in the saved analysis; the summary only claims a
  level of invariance if all steps below it hold, and no longer says that
  observed means become comparable.

## CFA vs ESEM: the ESEM can look worse than the CFA

- Reported by the author on real data (20 ordinal items, 3 factors, WLSMV):
  the ESEM had a higher chi-square and lower CFI than the CFA, which cannot
  be true of the fit itself. Reproduced on simulated data: the uncorrected
  chi-square favours the ESEM (663 vs 981) while the robust (scaled) one says
  the opposite (1467 vs about 1100), because each model gets its own
  correction. The model fit table now also shows the uncorrected chi-square,
  and a note appears when the two disagree. `docs/METHODS.md`, "CFA vs ESEM".

## Guttman scales

- The scalogram statistics are tested against the order in which the items
  were entered (easiest first) by default, with Kendall's tau between that
  order and the observed difficulty. "Sorted by endorsement" gives the old
  behaviour. For items with more than two categories the item steps are
  always sorted. Loevinger's H never depended on the order.
- New figure: scalogram (Scale Check > Guttman Scale; `plot(x, type =
  "scalogram")`).

## Colour

- Still off by default. The checkbox is now called "Colour". In CFA and
  Multi-Estimator CFA it adds a "Colour Summary": fit indices, standardized
  loadings, omega and AVE on a ten-step ramp, as the 2024 `cfawiz()` printed
  them, with a choice of palette: orange to purple (default, colour-blind
  safe), yellow to blue (colour-blind safe), red to blue (the 2024 ramp) or
  greys. Elsewhere it marks flagged values in red, as before.
- Figures are black, white and grey. Colouring comb teeth by loading size
  was tried and dropped the same day.

## Figures (decided by the author on 2026-11-16 from a sheet of prototypes)

- **Comb plot, new default look**: factor name written along the spine (no
  top tooth), factor correlations as arcs joining the spines (thicker =
  stronger), residual correlations in the model as brackets at the tooth
  tips. Option: tooth thickness proportional to the item's observed variance.
- **Second-order comb**: each factor gets a top tooth split into the share
  explained by the second-order factor and the rest.
- **Antenna plot** and **facing combs** for bifactor models; a bifactor path
  diagram.
- **Brush plot** for ESEM: every item on every factor, unique variance in a
  last column.
- The comb plot and the path diagram are separate checkboxes and can be
  shown together. Path diagram options: loadings on the arrows, arrow width
  by loading.
- Removed: the path diagram with variance tiles.
- *Breaking:* options `showDiagram` and `diagramStyle` are replaced by
  `showComb`, `combThick`, `showPath`, `pathNumbers`, `pathWidth`;
  `plot(type = "diagram")` is now `type = "path"`, `type = "tiles"` is gone,
  and `plot()` no longer takes `highlight`.

## Models

- **Bifactor CFA** (new analysis; `bifactorWiz()`): general factor plus
  specific factors, compared with the correlated-factors model; omega,
  omega hierarchical and explained common variance for the total score and
  each subscale (Rodriguez, Reise & Haviland, 2016).
- **CFA**: residual correlations between pairs of items, and a second-order
  factor in place of the factor correlations (three or more factors).
  Omega now counts residual covariances as error and uses the model-implied
  factor variance. Called from R, the jamovi function needs `resCov =
  list()` (a jamovi "Pairs" option has no default, as in `jmv::cfa()`).

## Item Order (exploratory)

- New analysis in *Scales* (`itemOrder`; `orderWiz()` in R) for yes/no items:
  looks for order instead of common factors. It reports cumulative
  (Mokken-type) scales and unscalable items, implications between items
  ("whoever endorses A also endorses B"), items that share a rung, branches
  and items that stand apart, how many respondents fit a single line and how
  many fit the order found, and draws the order as a Hasse diagram.
- Written from a longer brief that asked for the `mokken`, `DAKS` and
  `POSAC` packages. None ships with jamovi, so the methods were implemented
  here in a reduced form; what was left out (significance tests in the item
  selection, POSAC's two-dimensional pattern map, an MDS embedding) is listed
  in `docs/METHODS.md`.
- On the example debt items it recovers the structure they were built with:
  family = friends at the base, bank and moneylender as two branches, card
  on top of bank, mortgage apart; 55% of respondents fit a single line, 91%
  the order.

## Fixed after the first look in the jamovi window (2026-11-16)

- CFA, Group CFA, Scale Check and Bifactor CFA showed only "comparison (==)
  is possible only for atomic and list types". Cause: a result (the path
  diagram image) was named `path`, and Scale Check had a table named
  `items`; both are members of jamovi's own results classes, which the
  generated code then replaced. Renamed to `pathFig` and `itemStats`. The
  yaml generator now refuses reserved names and a test checks every
  analysis. Headless tests could not see this because they never serialise
  the results the way the jamovi window does.
- The "loom plot" is now the **brush plot**: unique variance is a plain last
  column instead of growing leftwards. A clean factor looks like a
  toothbrush.
- The R example that failed the 0.3 check on GitHub was `compositeCheck()`'s:
  it used an object called `data` that does not exist. Examples of the
  jamovi functions are now self-contained.

## Data

- New example data set `wizlab` ("A laboratory of instruments"): ordinal,
  numeric and binary instruments, good and bad, including a bifactor
  questionnaire, a Guttman scale and a partial order. `docs/DATA.md`
  explains every block and what to try on it.

## Internal

- `psw_composite_cors()` and `psw_score_cors()` are shared by `factorWiz()`
  and `esemWiz()`.
- Known: the 0.3.0 check on GitHub failed in an R example ("Argument 'data'
  must be a data frame"); not yet traced. It does not affect the module.

# psyWiz 0.3.0 (2026-11-14)

## New analyses

- **Scale Check** (`scaleWiz()`): for one scale, reliability (alpha,
  lambda-2, omega, ordinal alpha and omega, AVE), nested tests of
  tau-equivalence and parallelism with a plain-language reading for alpha
  and sum scores, item statistics, and Guttman / Mokken scalability
  (Loevinger's H for binary and ordinal items; scalogram statistics on item
  steps).
- **CFA vs ESEM** (`esemWiz()`): the hypothesised model as a CFA and as an
  ESEM with target or geomin rotation; loadings with cross-loadings; and each
  factor correlation under five scoring methods (composites, corrected
  composites, CFA, factor scores, ESEM).
- **Composite Check** (`compositeWiz()`): do formative indicators act
  through one composite? Rank-one test, variance retained per outcome, and a
  cross-validated comparison of separate weights, one estimated composite and
  the unit-weighted sum.
- **Robustness of the verdict** (`invarianceMultiverse()`, option in the
  group analyses): the invariance ladder under four estimator / item-treatment
  combinations and five decision rules, with the share of specifications that
  support each step.

## New figures

- Comb plot (`plot(x, type = "comb")`): every length is shared variance.
- Path diagram with variance tiles (`type = "tiles"`).
- Figure type is an option in the group analyses and in Scale Check.

## Other

- `docs/IDEAS.md`: what was built, what exists elsewhere, what was only
  thought through (response styles, non-linear loadings, PLS, Bayesian SEM).
- jamovi menu regrouped: Single Scale, Measurement Invariance, Model
  Comparison, Formative Indicators.

# psyWiz 0.2.0

Complete rewrite. The 2024 prototype (a single `cfawiz()` function that printed
to the console) is preserved on the `psywiz_old` branch and the
`v0.1.1-legacy` tag.

## New

- `invarianceWiz()`: measurement invariance ladder across groups for
  continuous items (configural, metric, scalar, strict) and for ordinal items
  (configural, thresholds, loadings, intercepts, optionally residuals;
  identification after Wu & Estabrook, 2016).
- `invarianceVerdict()`: flags and a verdict per step from the change in CFI,
  RMSEA, SRMR and/or the chi-square difference test, with editable cutoffs
  (defaults from Chen, 2007). Never refits.
- Score tests that rank the constraints of each step by misfit, computed from
  a single matrix inversion (same values as `lavaan::lavTestScore()`, far
  faster for large models).
- `cfaWiz()`: one single-group CFA under up to eight estimators, side by side.
  Successor of `cfawiz()`; keeps its coloured console tables.
- jamovi module with three analyses: One-Factor Group CFA,
  Correlated-Factors Group CFA and Multi-Estimator CFA. Each has a
  "Run analysis" checkbox, so nothing is fitted while options are being set.
- Independent models are fitted in parallel worker processes when the first
  model takes more than a second.
- Model diagrams with a fixed layout (`plot()` in R; "Model diagram" in the
  jamovi group analyses): one panel per group or estimator, arrow width
  proportional to the standardized loading.
- Colour is off by default everywhere. `print(x, colored = TRUE)`,
  `plot(x, highlight = TRUE)` and the jamovi "Highlight" checkbox turn it on.
- Example data `wizdata` (also in jamovi's data library).

## Changed from the prototype

- Items are assigned to factors explicitly instead of by name prefix.
- Reliability (omega) and AVE are computed per factor, not pooled over all
  loadings.
- One version of each fit index (robust, scaled or standard) is used
  consistently within a table and named in a footnote.
- Licence changed from MIT to GPL-3.

## Known gaps

- Help pages (`man/*.Rd`) are not generated yet: the roxygen comments are in
  the source, but roxygen2 was not available on the development machine.
- Not yet tested inside jamovi on Windows or macOS.
- See "What it does not do (yet)" in the README.
