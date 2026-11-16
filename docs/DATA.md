# Example data

Two simulated data sets ship with psyWiz. In jamovi they are in the data
library (Open > Data Library > psyWiz); in R, `data(wizdata)` and
`data(wizlab)`. Both are made by `tools/make_examples.R`, where the true
values are written down, so every analysis can be checked against what was
put in.

## wizdata: Trust and stress in three countries

Two correlated factors (r = -.40), four 5-point items each, Chile / Brazil /
Spain, n = 1200. Loadings are equal everywhere; item `stress3` is easier to
endorse in Brazil.

- **Group CFA** with `country`: thresholds and loadings hold, intercepts do
  not, and `stress3` in Brazil tops the worst-fitting constraints. Tick
  "Partial invariance": the intercept of `stress3` is released and the step
  holds.

## wizlab: a laboratory of instruments

n = 1500, three regions (`country`: North, Centre, South). Each block is a
small instrument with a known property.

| Items | Type | Built to be | Try | You should see |
|---|---|---|---|---|
| `calm1`-`calm5` | ordinal 1-5 | a clean scale (loadings .70-.80) | Scale Check; CFA | high omega, even comb teeth |
| `focus1`-`focus6` | ordinal 1-5 | a scale with four different problems | Scale Check; CFA with Calm and Focus | see below |
| `speed1`-`speed4` | numeric | tau-equivalent: equal loadings, unequal error | Scale Check > Tau-equivalence | equal loadings tenable, parallel not |
| `power1`-`power4` | numeric | congeneric, loadings .90 to .35, very unequal variances | Scale Check > Tau-equivalence; comb plot with "thickness by item variance" | tau-equivalence rejected; alpha below omega; teeth of very different thickness |
| `body*`, `mind*`, `social*` (4 each) | ordinal 1-5 | wellbeing: a general factor plus three specific ones | Bifactor CFA; CFA with a second-order factor | antenna plot leaning left (general); `body4` and `mind4` lean right |
| `know1`-`know6` | 0/1 | a near-perfect Guttman scale, entered easiest first | Scale Check > Guttman Scale; Item Order | H near 1, a clean staircase, one line of six rungs |
| `habit1`-`habit5` | 0/1 | not a scale | Scale Check; Item Order | low H, speckled scalogram, items left unconnected |
| `debt_*` (6) | 0/1 | a partial order with a branch and an unrelated item | Item Order | see below |

### The Focus items

- `focus1` also loads on Calm (.35). **CFA vs ESEM** with Calm and Focus: the
  brush plot shows a stray bristle for `focus1` on Calm, and the CFA factor
  correlation is higher than the ESEM one.
- `focus2` is easier in the South. **Group CFA**: intercepts step flagged,
  `focus2` in South on top.
- `focus4` is a dud (loading .15): a stub of a tooth in the comb plot, low
  item-rest correlation, alpha rises without it.
- `focus5` and `focus6` share wording. **CFA > Largest modification
  indices** points at their residual correlation; add the pair under
  "Residual Correlations" and a bracket joins the two teeth.

### The debt items

Built as: family and friends debt go together (the same rung); bank debt
occurs almost only among those with informal debt, and card debt among those
with bank debt; a moneylender is used by people with informal debt who are
shut out of formal credit, whether or not they have bank debt; a mortgage is
unrelated to all of it. So the truth is neither a line nor a set of factors:

```
card -> bank -> family = friends <- lender        mortgage (apart)
```

**Item Order** recovers this diagram. Scale Check on the same six items
gives a mediocre alpha and a speckled scalogram, and a factor analysis
would suggest "one weak factor": neither shows the structure.
