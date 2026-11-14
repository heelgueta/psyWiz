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
to look, not as confirmatory tests. psyWiz deliberately does not free
parameters automatically.

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

## Figures

All figures use fixed layouts (docs/DESIGN.md). Comb plot: every length is
shared variance on a 0 to 1 scale; top tooth = the factor (1), item tooth =
squared standardized loading on a pale tooth of length 1, thin teeth = squared
correlation with each other factor. Tiles: each item's square has unit area
and contains a filled square of area R².

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
