# psyWiz 0.3.0 (in development)

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
