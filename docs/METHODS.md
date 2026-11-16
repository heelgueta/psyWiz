# Methods

What psyWiz fits, in enough detail to reproduce it with plain lavaan. The
lavaan syntax of every model is returned by the R functions (`$syntax`).

## Models

Every model is a confirmatory factor model with simple structure: each item
loads on exactly one factor, factors are freely correlated, residuals are
uncorrelated. Factors are identified by fixing their variance to 1 (and, with
several groups, their mean to 0) in the first group.

## How items are treated

| Setting | Ordinal variable | Continuous variable |
|---|---|---|
| Their measure type (default) | ordinal | continuous |
| Ordinal | ordinal | ordinal (each distinct value is a category) |
| Continuous | continuous (category codes) | continuous |

All items of a model must be treated the same way. Ordinal treatment uses
polychoric correlations and thresholds; it needs every category of every item
to be observed in every group and allows at most 10 categories.

## Invariance ladder

Each step keeps the constraints of the previous steps and is compared with the
previous step. The first group is the reference group.

### Continuous items

| Step | Equal across groups | Freed in groups 2+ |
|---|---|---|
| Configural | nothing | |
| Loadings (metric) | loadings | factor variances |
| Intercepts (scalar) | + intercepts | factor means |
| Residuals (strict, optional) | + residual variances | |

These are the same models as `lavaan::cfa(..., group.equal = ...)`; the test
suite checks that the chi-squares agree.

### Ordinal items

The ladder and identification constraints follow Wu & Estabrook (2016),
because the older practice of constraining loadings before thresholds does not
give an identified, interpretable sequence for ordinal data.

| Step | Equal across groups | Freed in groups 2+ |
|---|---|---|
| Configural | nothing (intercepts 0, scales 1 in all groups) | |
| Thresholds | thresholds | item intercepts and scales |
| Loadings (metric) | + loadings | factor variances |
| Intercepts (scalar) | + intercepts (fixed to 0) | factor means |
| Residuals (strict, optional) | + residual variances (fixed to 1) | |

Notes:

- With three response categories the Thresholds model is statistically
  equivalent to the Configural model (same fit, same degrees of freedom), so
  that step cannot be tested and is reported as "Equivalent to previous model".
- Binary items: scales (residual variances) stay fixed in all groups because
  they are not identified, so threshold, loading and residual differences
  cannot be fully separated. psyWiz says so in a note.
- The default parameterisation is lavaan's "delta". The strict step needs
  "theta", so when it is requested the whole ladder is fitted with "theta".
- Degrees of freedom for p items with T thresholds each, m factors and G
  groups: thresholds step adds p(T-2)(G-1), loadings and intercepts steps add
  (p-m)(G-1) each. The test suite checks this.

## Estimators and missing values

| Items treated as | Default estimator | Others | Missing values (automatic) |
|---|---|---|---|
| continuous | MLR | ML, MLM | FIML for ML and MLR; listwise for MLM |
| ordinal | WLSMV | ULSMV, DWLS, ULS | pairwise |

With no missing values everything uses the complete data.

## Fit indices

lavaan can report standard, scaled and robust versions of the chi-square based
indices. psyWiz reports one version per index, **the same for every model in a
table**, so that differences between models never mix versions, and names the
version in a footnote:

- chi-square: scaled when the estimator provides it, else standard.
- CFI, TLI, RMSEA: robust > scaled > standard for continuous items; scaled >
  standard for ordinal items (lavaan's robust versions for ordinal data are
  experimental and frequently unavailable).

## Comparing steps

- **Chi-square difference test**: `lavaan::lavTestLRT()`, which applies the
  appropriate scaled difference test for robust estimators.
- **Changes in CFI, RMSEA, SRMR**: new step minus previous step.

A value is flagged when CFI drops by more than .010, RMSEA rises by more than
.015, SRMR rises by more than .030 (loadings step) or .010 (other steps), or
the difference test has p < .05. These defaults follow Chen (2007) and are all
editable. By default the verdict uses CFI and RMSEA:

- **Supported**: no selected criterion is flagged.
- **Not supported**: every selected criterion is flagged.
- **Ambiguous**: the selected criteria disagree.

Chen's cutoffs come from simulations with continuous indicators and ML. They
are conventions, not tests; with ordinal estimators treat them as rough guides.

## Worst-fitting constraints

For each step, univariate score (Lagrange multiplier) tests of the constraints
that step added: the expected drop in the model chi-square if that single
constraint were released for one group, relative to the first group. Equality
constraints use `lavaan::lavTestScore()`; parameters fixed to a constant in
groups 2+ (ordinal intercepts and residual variances) use
`lavaan::modindices()`.

lavaan computes these without the robust correction, and with many constraints
some will be large by chance. Use them to rank constraints and to decide where
to look, not as confirmatory tests.

## Partial invariance

Off by default. `invariancePartial()`, or "Partial invariance" in Group CFA,
refits the ladder and lets single items out of the equality constraints:

1. Each step is fitted and judged against the previous model with the
   current decision criteria.
2. If the verdict is "Not supported" or "Ambiguous", the constraints that
   step adds are ranked by their score test (previous section). The item at
   the top is left free in every group, the step is refitted and judged
   again against the same previous model. This repeats until the step is
   supported or no release is allowed.
3. A release is allowed only if at least two items of the factor stay
   constrained (`min.invariant`) and no more than `max.free` (default 2)
   items of that factor are free on that parameter. Two constrained items
   per factor is the usual minimum for the factor to be on a common scale
   (Byrne, Shavelson & Muthen, 1989; Steenkamp & Baumgartner, 1998).
4. Whatever is freed stays free in later steps, and a freed loading also
   frees that item's intercept at the intercepts step; a freed intercept
   frees its residual variance at the strict step.
5. If a step is still not supported when the search stops, the releases made
   are kept, the ladder continues, and the summary reports the highest level
   reached with every step below it supported.

Limits, stated in the output as well:

- Thresholds of ordinal items are never released. An item with free
  thresholds needs its own identification constraints (Wu & Estabrook, 2016),
  which the syntax builder does not write yet. A failed thresholds step ends
  the claim.
- An item is released in all groups at once, although the score test that
  picked it belongs to one group.
- It is a search guided by the data. With many items, some release will
  improve fit by chance, and the tests of the later steps are conditional on
  the search. Partial scalar invariance supports comparing latent means in a
  model that keeps the freed parameters free; observed sum scores still
  contain the freed items.

## Plain CFA

`factorWiz()` / the CFA analysis: the single-group model of section "Models",
factor variances fixed to 1. Reported: one version of each fit index (as in
"Fit indices"), unstandardised and standardised loadings, R-squared
(standardised loading squared), factor correlations.

Reliability per factor: alpha from the items as scored (complete cases of
that factor's items); omega = (sum of loadings)^2 / ((sum of loadings)^2 +
sum of residual variances) from the unstandardised solution. For continuous
items that is the reliability of the sum score under the model; for ordinal
items the solution is on the latent-response scale, so omega describes the
latent responses and is usually higher than what the sum score achieves
(Scale Check reports both). AVE = mean squared standardised loading.

Factor correlations by scoring method (two or more factors), as in "CFA vs
ESEM" without the ESEM column: composite (item mean) scores; composites
corrected for attenuation with the two alphas; the CFA; regression factor
scores and, for continuous items, Bartlett scores.

Largest modification indices: `lavaan::modindices()`, restricted to residual
covariances between two items and loadings of an item on another factor,
sorted by size. The expected value is lavaan's fully standardised expected
parameter change. They are not robust-corrected.

## Multi-estimator comparison

One single-group model per estimator. ML, MLR and MLM treat items as
continuous; ULS, ULSMV, DWLS, WLSMV and WLS treat ordinal items as ordinal and
continuous items as continuous. Each estimator's own best version of each fit
index is reported and named in a footnote. ULS and DWLS without correction
give no valid p value; WLS needs large samples.

Composite reliability (omega) and average variance extracted (AVE) are
computed per factor from the standardised loadings λ:
ω = (Σλ)² / ((Σλ)² + Σ(1 − λ²)), AVE = mean(λ²). Under ordinal treatment they
describe the latent response variables, not the observed sum score.

## Scale Check

For one set of items. Classical statistics use complete cases.

**Reliability.** Cronbach's alpha and Guttman's lambda-2 from the item
covariance matrix; McDonald's omega for the sum score from the congeneric
model on the items as scored, ω = (Σλ)² / ((Σλ)² + Σθ). With ordinal items,
ordinal alpha (alpha on the polychoric matrix) and ordinal omega (from
standardized loadings) are added and labelled as referring to the latent
responses.

**Measurement models.** Congeneric (free loadings), tau-equivalent (equal
unstandardized loadings) and parallel (equal loadings and residual
variances), each compared with the previous one by the same machinery and
decision rules as the invariance ladder. Means are not constrained, so
"tau-equivalent" here is essential tau-equivalence. With ordinal items
(delta parameterisation) equal loadings imply equal residual variances, so
only congeneric and tau-equivalent are fitted. Reading: alpha equals
reliability only under tau-equivalence; unit weights are optimal only for
parallel items.

**Scalability.** Loevinger's H from covariances: H_ij = cov(X_i, X_j) /
cov_max(X_i, X_j), where cov_max is the covariance when both items are sorted
(the largest their marginals allow); item H_i and scale H are ratios of the
corresponding sums. For binary items this is 1 minus observed over expected
Guttman errors; for polytomous items it is Molenaar's (1991) coefficient.
Conventional labels: below .30 not scalable, .30 to .40 weak, .40 to .50
medium, .50 or more strong. Guttman scalogram statistics are computed on item
steps (response at or above each category), ordered by popularity, with
Goodenough-Edwards error counting: CR = 1 - errors / (N x steps), MMR = mean
of max(p, 1 - p), CS = (CR - MMR) / (1 - MMR).

### Order of the items in the Guttman statistics

A Guttman scale is a claim about order: whoever endorses a hard item also
endorses every easier one. CR, MMR and CS count deviations from that
pattern, so they need an order of the items.

- *In the order entered* (default): the order in which the items were given
  is taken as the hypothesised order, easiest first. This tests your
  ordering. Kendall's tau between that order and the observed endorsement
  rates is added (1 = entered exactly from easiest to hardest).
- *Sorted by endorsement*: items sorted from the most to the least endorsed
  in this sample. This is the order that fits the sample best, so CR is at
  its most favourable; it is the usual exploratory scalogram.

For items with more than two categories the statistics are computed on item
steps, and the steps are always sorted by endorsement (an order of items
does not define an order of their steps); only tau uses the entered order.
Loevinger's H does not depend on the order at all.

The scalogram figure shows respondents sorted by total score (averaged in
at most 150 groups) against the items in the chosen order; for binary items
a line marks where a perfect Guttman scale would end.

## CFA vs ESEM

The CFA is the correlated-factors model above. The ESEM lets every item load
on every factor (`efa()` blocks in lavaan) and rotates obliquely: towards a
target with zeros for all cross-loadings (default) or by geomin. ESEM factors
are matched to the hypothesised ones by their mean absolute loading on each
factor's own items and reflected to be positive. The CFA is nested in the
ESEM, so the usual difference test applies.

Factor correlations by method: composite = correlation of item means;
corrected = composite / sqrt(alpha_1 x alpha_2); CFA and ESEM = latent
correlations; factor scores = correlations of `lavaan::lavPredict()` scores
from the CFA (regression method; Bartlett for continuous items).

### When the ESEM looks worse than the CFA

The CFA is the ESEM with the cross-loadings fixed to zero, so the ESEM's
discrepancy from the data can only be equal or smaller. With robust test
statistics (WLSMV, MLR) the table can still show a *higher* chi-square and
*lower* CFI for the ESEM. The reason: the robust chi-square is the
uncorrected one divided by a scaling factor (and shifted) that is estimated
separately for each model, and the two factors can differ a lot. In the
reproduction of 2026-11-16 (20 five-point items, 3 factors, a nuisance
factor on six items, n = 900, WLSMV) the uncorrected chi-squares were 981
(CFA) and 663 (ESEM) but the scaling factor of the ESEM was 0.46, giving a
scaled value of 1467. CFI, TLI and RMSEA are computed from the scaled values
and inherit the reversal.

What psyWiz does: the fit table shows the uncorrected chi-square next to the
corrected one, and when they disagree a note says that the corrected indices
of the two models are not comparable. What to read instead: the difference
test (which is built from the uncorrected statistics and corrected as one
test) and SRMR (which involves no correction). Whether lavaan's correction
for rotated solutions is itself right in such cases was not investigated.

## Robustness of the invariance verdict

`invarianceMultiverse()` refits the ladder under ordinal treatment (WLSMV,
ULSMV) and continuous treatment (MLR, ML) and applies five single-criterion
rules to every step: CFI drop > .010; CFI drop > .002 (Meade et al., 2008);
RMSEA rise > .015; SRMR rise > .030 (loadings) or .010; difference test
p < .05. The summary counts, per step, the specifications in which it is
supported. Shares of .90 or more are labelled robustly supported, .10 or less
robustly not supported, otherwise "depends on choices". The specifications
are correlated, so the share is descriptive. See docs/IDEAS.md, section 2.

## Composite Check

Indicators are standardized; cases with missing values are dropped.
"Separate weights" is ordinary regression of each outcome on all indicators.
"One composite" is rank-one reduced-rank regression: the weights are the
first canonical variate of the indicators with respect to the outcomes, and
each outcome is regressed on that composite. "Unit-weighted" uses the sum of
the standardized indicators. The single-composite hypothesis is tested with
Bartlett's statistic -(N - 1 - (p + q + 1)/2) Σ ln(1 - ρ_i²) over the
second and later canonical correlations, df = (p - 1)(q - 1). Cross-validated
R² is 1 - SSE/SST over held-out cases, k folds, averaged over repetitions
with a fixed seed. See docs/IDEAS.md, section 3.

## Residual correlations and second-order factor (CFA)

Residual correlations: `i ~~ j` for each chosen pair of items; reported as
correlations between the two unique parts. A second-order factor: the
first-order factors load on one more factor and their correlations are no
longer free. With exactly three factors this is the same model as three
correlated factors (same fit), only described differently; from four factors
on it is a restriction. Reliability: omega = T / (T + E), with T = (sum of
unstandardised loadings)^2 x model-implied factor variance and E = sum of
all residual variances and covariances among the factor's items.

## Bifactor CFA

Every item loads on a general factor and on one specific factor; all factors
are uncorrelated (`g ~~ 0*f`, `f ~~ 0*f'`), variances fixed to 1. Each
specific factor needs three or more items. The correlated-factors model is
fitted alongside for comparison; a bifactor model tends to fit better even
when it is not the better description, so the indices matter more than the
fit. With standardised loadings g (general) and s (specific) and unique
variances u = 1 - g^2 - s^2 (Rodriguez, Reise & Haviland, 2016):

- omega, total score: [(sum g)^2 + sum over factors of (sum s)^2] / V, where
  V is that numerator plus the sum of u.
- omega hierarchical: (sum g)^2 / V, the share of the total score's variance
  due to the general factor.
- for a subscale: omega = [(sum g)^2 + (sum s)^2] / V_k and omega
  hierarchical = (sum s)^2 / V_k over its own items: what the subscale score
  measures beyond the general factor.
- ECV = sum g^2 / (sum g^2 + sum s^2); item ECV the same for one item; PUC =
  share of item pairs from different specific factors.

For ordinal items all of this refers to the latent responses.

## Figures

Black, white and grey; fixed layouts; sizes chosen by the code.

- **Comb plot.** One comb per factor, the factor's name along the spine. A
  tooth is an item: dark part = squared standardised loading, pale remainder
  = unique variance, so every length is a share of variance between 0 and 1.
  Arcs join the spines of correlated factors, line width 0.4 + 3.6|r| (not
  drawn for more than four factors). Brackets join the tips of items whose
  residuals correlate in the model. Optional: tooth thickness proportional
  to the item's observed variance (items as scored, per group), relative to
  the largest in the figure. Second-order model: each factor has a top tooth
  whose dark part is its squared loading on the second-order factor.
- **Path diagram.** Factors left, items right, an arrow per loading (width
  0.25 + 1.5|loading| unless switched off).
- **Scalogram.** See "Order of the items in the Guttman statistics".
- **Antenna plot** (bifactor). One row per item; from a mast in the middle,
  the bar to the left is g^2 and the bar to the right is s^2.
- **Facing combs** (bifactor). Each row is the item's whole variance: g^2
  from the left spine, s^2 from the right spine, u in between, where the
  item's name is written.
- **Loom plot** (ESEM). Items as rows, one vertical line per factor; to the
  right of each line the squared standardised loading on that factor (dark
  for the item's own factor, grey for cross-loadings); left of the first
  line the item's unique variance from the model. With correlated factors
  the squared loadings of an item do not add up to its explained variance.

## Colour

Off by default. When on, a value is placed on a ten-step ramp by fixed bands
per statistic (`psw_color_index()` in `R/print.R`): CFI and TLI .70 to .99,
RMSEA .30 to .01, SRMR .35 to .03, loadings and AVE .1 to .9, omega .50 to
.95. The bands are those of the 2024 prototype and are a reading aid, not
cutoffs. The ramp is used in the "Colour Summary" tables only, not in figures. Palettes (`R/colour.R`): orange to purple (ColorBrewer PuOr,
default) and yellow to blue (cividis) are chosen to survive the common forms
of colour blindness; red to blue is the original ramp.

## References

- Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
  measurement invariance. *Structural Equation Modeling, 14*(3), 464-504.
- Franke, G. R., Preacher, K. J., & Rigdon, E. E. (2008). Proportional
  structural effects of formative indicators. *Journal of Business Research,
  61*(12), 1229-1237.
- Graham, J. M. (2006). Congeneric and (essentially) tau-equivalent estimates
  of score reliability. *Educational and Psychological Measurement, 66*(6),
  930-944.
- Marsh, H. W., Morin, A. J. S., Parker, P. D., & Kaur, G. (2014).
  Exploratory structural equation modeling. *Annual Review of Clinical
  Psychology, 10*, 85-110.
- Meade, A. W., Johnson, E. C., & Braddy, P. W. (2008). Power and sensitivity
  of alternative fit indices in tests of measurement invariance. *Journal of
  Applied Psychology, 93*(3), 568-592.
- Molenaar, I. W. (1991). A weighted Loevinger H-coefficient extending Mokken
  scaling to multicategory items. *Kwantitatieve Methoden, 37*, 97-117.
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
  *Journal of Statistical Software, 48*(2), 1-36.
- Steegen, S., Tuerlinckx, F., Gelman, A., & Vanpaemel, W. (2016). Increasing
  transparency through a multiverse analysis. *Perspectives on Psychological
  Science, 11*(5), 702-712.
- Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
  analysis models of different levels of invariance for ordered categorical
  outcomes. *Psychometrika, 81*(4), 1014-1045.
