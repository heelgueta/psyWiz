# psyWiz 0.2.0 (in development)

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
