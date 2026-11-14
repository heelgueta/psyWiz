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

## References

- Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
  measurement invariance. *Structural Equation Modeling, 14*(3), 464-504.
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
  *Journal of Statistical Software, 48*(2), 1-36.
- Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
  analysis models of different levels of invariance for ordered categorical
  outcomes. *Psychometrika, 81*(4), 1014-1045.
