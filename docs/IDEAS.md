# Ideas, experiments and where they stand

A lab notebook for the psychometric ideas explored in psyWiz: what was
built, what already exists elsewhere, what was only thought through, and what
failed. Dated entries; newest thinking at the top of each section. Literature
claims marked *(searched)* come from a quick web search on the date given and
are not a systematic review; treat "we found nothing" as "worth a proper
search before claiming novelty".

## 1. Built in version 0.3 (2026-11-14)

| Idea | Where | Status of the idea in the literature |
|---|---|---|
| Tau-equivalence and parallelism tested as nested models, with a plain-language reading for alpha and sum scores | Scale Check, `scaleWiz()` | Classical (Joreskog 1971; Graham 2006). Rarely run in practice because it needs hand-written constraints. |
| Guttman scalogram statistics for ordinal items via item steps, and Loevinger's H | Scale Check | Classical: Mokken (1971), Molenaar (1991). H *is* the "pseudo-Guttman for ordinal scales". Implemented from the covariance definition; not cross-checked against the `mokken` package (not installable here). |
| ESEM with target rotation next to the CFA | CFA vs ESEM, `esemWiz()` | Established (Asparouhov & Muthen 2009; Marsh et al. 2014). Target rotation is still the standard confirmatory use of ESEM. |
| One factor correlation under five scorings (composites, corrected composites, CFA, regression and Bartlett factor scores, ESEM) | CFA vs ESEM | Each piece is textbook; the side-by-side table is a convenience, and a good teaching device: on the example data the same correlation ranges from -.35 to -.51. |
| Comb plot (redesigned in 0.4, see section 5) | `plot(x, type = "comb")`, figure option in jamovi | The author's idea. Closest precedent found *(searched)*: bar charts of loadings. A figure in which every length is shared variance on one 0-1 scale (factor = 1, item = loading squared, other factor = r squared) was not found. |
| Variance tiles (removed in 0.4) | was `type = "tiles"` | Path diagram where each item has a square (its variance) containing a filled square of area R-squared. Area encodings of variance exist in many forms; not claimed as new. |
| **Verdict multiverse for invariance** | "Robustness of the verdict", `invarianceMultiverse()` | See section 2. |
| **Composite sufficiency check for formative indicators** | Composite Check, `compositeWiz()` | See section 3. |

Added in 0.4 (2026-11-16): partial invariance by score-test search (`invariancePartial()`), which is classical (Byrne et al. 1989; Steenkamp & Baumgartner 1998), and the plain CFA with modification indices and correlations by scoring method. Nothing new as an idea.

## 2. The idea that is ours: a multiverse of invariance verdicts

**Problem.** "Scalar invariance was (not) supported" is reported as a fact
about the data. It is really a fact about the data *plus* choices that could
defensibly have been made otherwise: ordinal or continuous treatment of
Likert items, the estimator, which fit index, which cutoff. Simulation work
keeps showing that cutoffs only hold for the conditions they were derived
under and that the estimator moves the indices *(searched: Groskurth
dissertation; Sokolov working paper; Meade et al. 2008)*. Yet the applied
researcher still picks one path and reports one verdict.

**Idea.** Run all defensible paths and report the verdict as a distribution:
the ladder is refitted under ordinal (WLSMV, ULSMV) and continuous (MLR, ML)
treatment, and five decision rules are applied to every step. A step
supported in 20 of 20 specifications is a robust finding; one supported in 8
of 10 "depends on choices", and the grid shows which choices.

**Is it new?** Multiverse analysis is not (Steegen et al. 2016). Two searches
on 2026-11-14 found no application of it to invariance verdicts and no tool
that does it. It is cheap for us because the engine already fits ladders and
separates fitting from verdicts.

**Limits.** The specifications are not independent, so the share is a
description, not a probability. A multiverse can launder bad options (Del
Giudice & Gangestad 2021): every specification included must be one a
reasonable analyst would defend, which is why WLS and uncorrected ULS/DWLS
are not in it. It does not replace thinking about *why* a step fails.

**Next.** Add the reference group and partial-invariance decisions as
dimensions; a specification-curve figure; a small simulation showing how
often the single-path verdict disagrees with the multiverse majority.

## 3. Formative models: what can actually fail?

The standard complaint is right: internal consistency is meaningless for
formative indicators, and "the composite predicts Y" cannot fail. What can
fail is the claim that the indicators act **through one composite**. Then
each indicator's effects on different outcomes must be proportional (Franke,
Preacher & Rigdon 2008), i.e. the indicator-by-outcome coefficient matrix has
rank one. That is the reduced-rank regression model of econometrics
(Anderson 1951) and has a classical test (Bartlett). The connection between
the proportionality constraint and rank-one regression is our reading; the
2008 paper frames it as a MIMIC constraint.

Composite Check adds what the test alone does not give:

- how much of each outcome's explainable variance one composite retains;
- the same comparison **out of sample**, because separate weights always win
  in sample;
- the unit-weighted sum as a third contender, since that is what people
  actually compute, and it often predicts as well as estimated weights.

Borrowed from elsewhere, not yet built: ecology and economics evaluate
indices axiomatically (monotonicity, no indicator can lower the index when
it improves) and by weight stability under resampling; chemometrics uses
PLS components precisely as cross-validated composites. A bootstrap of the
weights with sign-stability rates would be the natural next table.

## 4. PLS and Bayesian estimation: what do they add for a psychometrician?

**PLS-SEM.** Builds each construct as a weighted sum of its indicators
(a composite) and estimates paths between composites by repeated
regressions. What it gives you: it always converges, works with small
samples and many constructs, handles formative blocks naturally, and is
oriented to prediction. What it costs: composites contain measurement error,
so with truly reflective scales loadings are overestimated and structural
paths underestimated (consistent PLS, "PLSc", corrects this); there is no
model chi-square, so the measurement model is never really put at risk; and
the classical toolkit you already use (CFA fit, invariance ladders) has only
rough equivalents (MICOM). Verdict for scale psychometrics: it does not
cover anything CFA/ESEM does not, *except* genuinely formative or
forecasting-oriented models, which is exactly the territory of section 3.
If you ever need it: the `seminr` or `cSEM` packages.

**Bayesian SEM.** Same models, different estimation: parameters get prior
distributions and the result is a posterior distribution. What it adds that
ML/WLS do not:

1. *Approximate zeros.* Cross-loadings or residual covariances can be given
   small-variance priors around zero instead of being fixed at zero. That is
   a principled middle ground between CFA (too strict) and ESEM (everything
   free), and the same trick gives *approximate measurement invariance*
   (small group differences allowed), which is useful with many groups.
2. *Small samples and awkward models.* No asymptotics; models that will not
   converge or give Heywood cases under ML often behave with mild priors.
3. *Direct uncertainty statements* for any derived quantity: an interval
   for omega, for a difference between two reliabilities, for a factor score.
4. *Ordinal and missing data* handled in one framework.

Costs: slow (minutes to hours), results depend on priors when N is small
(which is when you want it), and fit assessment (posterior predictive p,
BRMSEA) is less familiar. In R: `blavaan`, which needs Stan or JAGS and
therefore cannot ship inside a jamovi module easily. Verdict: worth learning
for approximate invariance and small-N reliability intervals; not needed for
the everyday large-sample scale validation psyWiz targets.

## 5. Thought through, not built

Figures, 2026-11-16. The author proposed a set of variants; prototypes
were drawn with a throwaway script (`tools/plot_proposals.R`, deleted once the
choices were built; it was never committed) and the author chose:

| Idea | Outcome |
|---|---|
| Factor name along the spine | Adopted as the comb default |
| Correlations as arcs between spines | Adopted (up to four factors) |
| Residual correlations as brackets | Adopted for correlations in the model. Not built: a dashed bracket for a pair only *suggested* by a modification index |
| Tooth thickness by item variance | Adopted as an option, off by default |
| Second-order comb | Adopted, with the split top tooth |
| Antenna plot for bifactor models | Adopted, default figure of Bifactor CFA |
| Unexplained variance in the middle, combs from both sides | Built as "facing combs" |
| "Grill" for EFA / ESEM | Built and named **loom plot** (items are the weft, factors the warp). The name is a proposal; rename freely |
| Colouring the figures | Tried (teeth coloured by loading size), dropped |
| Path diagram with variance tiles | Removed |
| Scalogram | Adopted in Scale Check |

Open: arcs for more than four factors (a small matrix instead?); the loom
plot for a plain EFA (psyWiz has no EFA analysis; jamovi has one).

Also wanted: one example data set with several short scales that differ in
item type (binary, ordinal, numeric) and in quality (clean, cross-loading,
residual pair, non-invariant item, dud item, bifactor, second-order), with a
paragraph explaining each.

- **Extreme responding vs. content.** Split each response into direction
  (agree / disagree) and extremity (distance from the midpoint). This is the
  logic of IRTree and multiprocess models (Bockenholt 2012) and of
  threshold-scaling models of response style (Jin & Wang 2014; Lubbe &
  Schuster 2020) *(searched)*. A light version fits our engine: in the
  ordinal invariance ladder, the thresholds step frees an intercept (shift)
  and a scale (stretch) per item in each non-reference group. If the stretch
  is about the same for all items of *different* factors, the group uses the
  response scale more or less extremely, a style, not item bias; item-specific
  deviations are the real non-invariance. A "shift and stretch" table per
  group would make that visible. Caveat: within one factor a common stretch
  cannot be told apart from a real difference in variance.
- **Non-linear loadings.** lavaan cannot estimate them. A cheap, honest
  diagnostic: regress each item on the rest score with a quadratic term (or
  plot the nonparametric item-rest regression, as Mokken analysis does for
  monotonicity) and flag curvature. For ordinal items, part of what looks
  like non-linearity is floor and ceiling, which the threshold model already
  absorbs; the check should be run on the latent-response scale.
- **Categorisation loss.** How much reliability does a 5-point format throw
  away compared with the underlying continuous response? Computable from the
  estimated thresholds; the gap between ordinal omega and sum-score omega in
  Scale Check is a first glimpse of it.
- **Cross-validated choice between sum scores and weighted scores.** The
  chi-square test asks whether tau-equivalence is exactly true; the practical
  question is whether weights estimated in a sample of this size score new
  people better than unit weights. Much is known (Wainer 1976; Grice 2001)
  *(searched)*; a per-dataset cross-validation would be a small addition to
  Scale Check.

## 6. Losses and open issues (2026-11-14)

- Nothing new in this version has been seen in the jamovi window.
- `mokken` and `semTools` are not available in jamovi's R here, so Loevinger's
  H and the Wu-Estabrook ladder are verified by their defining properties and
  degrees of freedom, not against those packages.
- The comb plot shows factor correlations only for up to four factors, and
  relates combs by listing r-squared teeth; a layout in which two combs share
  a physical overlap was sketched and dropped (it only works for two factors).
- ESEM factors are matched to the hypothesised factors by their loadings on
  their own items; with a badly wrong hypothesis that matching is arbitrary.
- Factor scores for ordinal items come from lavaan's empirical Bayes
  predictor; Bartlett scores are continuous-only.
- The robustness check refits four ladders and is slow for large ordinal
  models.
- Bifactor models and latent mean comparisons are still not implemented.
- 0.4: the partial invariance search does not release thresholds.
