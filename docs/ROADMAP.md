# Roadmap

Where psyWiz came from, where it stands and what comes next. Update this
file whenever a direction changes; NEWS.md records what shipped, this file
records why and what is intended.

## Where it came from

| When | What |
|---|---|
| 2024 | Prototype: one R function, `cfawiz()`, printing a multi-estimator CFA to the console. Kept on branch `psywiz_old`, tag `v0.1.1-legacy`. |
| 2026-11-14, v0.2 | Rewrite as R package + jamovi module: invariance ladder (continuous and ordinal), verdicts with editable cutoffs, fast score tests, multi-estimator CFA, parallel fitting. |
| 2026-11-14, v0.3 | Breadth: Scale Check, CFA vs ESEM, Composite Check, verdict multiverse, comb plot and variance tiles. Six menu entries. |
| 2026-11-15 | Another coding session added partial invariance (uncommitted). |
| 2026-11-16, v0.4 | Simplification, then figures and models: comb redesign, scalogram, antenna, facing combs, brush plot; residual correlations, second-order factor, Bifactor CFA. The jamovi module is the product; plain CFA added as the entry point; the two group analyses merged; Scale Check reduced to reliability and items by default; partial invariance reviewed, corrected and kept as an option. |

## Where it stands (2026-11-16)

- Eight analyses in two menu groups (see README). The core is CFA, Group CFA
  and Scale Check.
- The author has run CFA vs ESEM in the jamovi window on real data (20 items,
  3 factors; it works, and it exposed the scaled chi-square reversal now
  documented in METHODS). The other analyses have only been checked headless
  (testthat, including the jamovi analyses run through their generated R
  functions).
- GitHub Actions passed for 0.2 on Linux, macOS and Windows and failed for 0.3
  on an R example (not traced).
- Not verified against `mokken` or `semTools`.

## Direction

1. **jamovi module first.** Goal: a module good enough to submit to the
   jamovi library. The R package comes after that, from this same repository
   (a jamovi module is an R package; see docs/DESIGN.md for why there is no
   separate folder).
2. **Simple defaults, depth on request.** Each analysis opens with the
   minimum. Simpler menus are contained in the more complex ones on purpose.
3. **Say what is a search.** Anything guided by the data (modification
   indices, partial invariance) is labelled as such in its own output.

## Before submitting to the jamovi library

1. **Look at every analysis in the jamovi window.** The first look (CFA,
   2026-11-16) found a bug that no headless test could see (a result named
   like one of jamovi's own members; fixed, and now guarded). Still unseen:
   the "Residual Correlations" pairs box, the nested checkboxes in
   "Figures", figure sizes, the Html "Colour Summary", Item Order.
2. **Check the numbers against other software**: Loevinger's H and the item
   selection against `mokken`, implications against `DAKS`, the ordinal
   invariance ladder against `semTools`, the bifactor indices against
   `BifactorIndicesCalculator`. None is installable in the development
   environment; needs a normal R installation.
3. **R side hygiene, because the module is an R package**: keep the GitHub
   check green (0.3 failed on an example that used an undefined `data`;
   fixed); generate the help pages; `label.size` in `geom_label()` is
   deprecated in recent ggplot2.
4. **Other platforms**: build and open the module in jamovi on Windows and
   macOS; confirm lavaan is found (the module relies on the copy that ships
   with jamovi).
5. **Submission procedure** (dev.jamovi.org): not read yet.

Done since the list was first written: committed and pushed (0.4); example
data that exercise everything (`wizlab`, docs/DATA.md).

## After that

- Group CFA with residual correlations, second-order and bifactor models.
- Latent mean comparisons; thresholds in the partial invariance search.
- Dashed bracket for residual correlations suggested by modification indices.
- Response-style "shift and stretch" table, non-linear loading check.
- CRAN (docs/CRAN.md).

## Parked

- lavaan syntax box.
- Warm starts between ladder steps (no gain; broke modification indices).
- Sharing sample statistics across estimators.
