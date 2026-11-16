# Simulates the example data shipped with psyWiz and writes
#   data/wizdata.rda                 (R users)
#   data/wizdata.omv  (jamovi users; needs jmvReadWrite, shipped with jamovi)
# Run from the repository root with jamovi's R, see docs/DEVELOPING.md.
#
# Design: two correlated factors (Trust, Stress), four 5-point Likert items
# each, three countries. Loadings are equal everywhere, but in country
# "Brazil" item stress3 is 'easier' to endorse (shifted intercept), so the
# loadings step should hold and the intercepts step should fail on stress3.
set.seed(20261114)
sim_group <- function(n, country, means = c(0, 0), shift = 0) {
  lam <- c(.80, .75, .70, .65, .80, .75, .70, .65)
  eta <- MASS::mvrnorm(n, means, matrix(c(1, -.4, -.4, 1), 2))
  ystar <- cbind(eta[, 1] %o% lam[1:4], eta[, 2] %o% lam[5:8]) +
    sapply(sqrt(1 - lam^2), function(s) rnorm(n, 0, s))
  ystar[, 7] <- ystar[, 7] + shift
  thr <- c(-1.4, -0.5, 0.4, 1.3)
  y <- apply(ystar, 2, function(v) findInterval(v, thr) + 1L)
  d <- as.data.frame(y)
  names(d) <- c(paste0("trust", 1:4), paste0("stress", 1:4))
  d[] <- lapply(d, function(v) factor(v, levels = 1:5, ordered = TRUE))
  cbind(country = country, d)
}
wizdata <- rbind(sim_group(420, "Chile"),
                 sim_group(380, "Brazil", means = c(-0.3, 0.4), shift = 0.6),
                 sim_group(400, "Spain", means = c(0.2, -0.2)))
wizdata$country <- factor(wizdata$country, levels = c("Chile", "Brazil", "Spain"))
# a little missingness, as in real surveys
for (v in names(wizdata)[-1]) wizdata[[v]][sample(nrow(wizdata), 12)] <- NA
rownames(wizdata) <- NULL
dir.create("data", showWarnings = FALSE)
save(wizdata, file = "data/wizdata.rda", compress = "xz")
if (requireNamespace("jmvReadWrite", quietly = TRUE))
  jmvReadWrite::write_omv(wizdata, "data/wizdata.omv", frcWrt = TRUE)

# ------------------------------------------------------------------------------
# wizlab: one file with many short instruments, good and bad, of every item
# type, so that each analysis and each figure of psyWiz has something to show.
# What every block is built to show is written up in docs/DATA.md; keep the two
# in step. Writes data/wizlab.rda and data/wizlab.omv.
set.seed(20261116)
n <- 1500
country <- factor(rep(c("North", "Centre", "South"), each = n / 3), levels = c("North", "Centre", "South"))
likert <- function(ystar, thr = c(-1.3, -0.4, 0.5, 1.4)) factor(findInterval(ystar, thr) + 1L, levels = 1:5, ordered = TRUE)
item <- function(f, lam, extra = 0) lam * f + extra + rnorm(n, 0, sqrt(max(1 - lam^2, 0.05)))

# -- Calm and Focus: two correlated ordinal scales (r = .5)
cf <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .5, .5, 1), 2))
calm <- sapply(c(.80, .75, .75, .70, .70), function(l) item(cf[, 1], l))          # a clean scale
pair <- rnorm(n)                                                                 # shared wording of focus5 and focus6
focus <- cbind(
  item(cf[, 2], .70, extra = .35 * cf[, 1]),                                      # focus1: cross-loads on Calm
  item(cf[, 2], .75, extra = ifelse(country == "South", .7, 0)),                  # focus2: easier in the South (non-invariant intercept)
  item(cf[, 2], .75),
  item(cf[, 2], .15),                                                             # focus4: a dud
  item(cf[, 2], .65, extra = .45 * pair),                                         # focus5, focus6: residual pair
  item(cf[, 2], .65, extra = .45 * pair))

# -- Speed and Power: numeric scales (r = .4)
sp <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .4, .4, 1), 2))
speed <- sapply(c(.6, .8, 1.0, 1.2), function(e) 50 + 8 * (sp[, 1] + rnorm(n, 0, e)))     # tau-equivalent: equal loadings, unequal error
power <- mapply(function(l, s) 100 + s * (l * sp[, 2] + rnorm(n, 0, sqrt(1 - l^2))),      # congeneric, very unequal loadings and variances
                c(.90, .75, .55, .35), c(5, 20, 10, 40))

# -- Wellbeing: 12 ordinal items with a general factor and three specific ones (bifactor / second-order)
g <- rnorm(n)
wb <- do.call(cbind, lapply(list(c(.70, .65, .60, .30), c(.55, .60, .50, .65), c(.75, .40, .50, .60)), function(gl) {
  s <- rnorm(n)
  sl <- c(.40, .45, .50, .65)
  sapply(1:4, function(j) gl[j] * g + sl[j] * s + rnorm(n, 0, sqrt(1 - gl[j]^2 - sl[j]^2)))
}))

# -- Knowledge: six binary items forming a near-perfect Guttman scale, entered from easiest to hardest
ability <- rnorm(n)
know <- sapply(c(-1.5, -0.8, -0.2, 0.4, 1.0, 1.6), function(b) as.integer(ability + rnorm(n, 0, 0.25) > b))

# -- Habits: five binary items that do not form a scale (weak, one unrelated, one reversed)
h <- rnorm(n)
habit <- cbind(as.integer(.5 * h + rnorm(n) > 0), as.integer(.4 * h + rnorm(n) > .3), as.integer(.3 * h + rnorm(n) > -.4),
               as.integer(rnorm(n) > .2), as.integer(-.4 * h + rnorm(n) > 0))

# -- Debt sources: a partial order, not a line.
#      family = friends (same rung)  <-  bank  <-  card
#                                    <-  lender            (a branch: needs informal debt, unrelated to bank debt)
#      mortgage: unrelated to all of it
strain <- rnorm(n); excluded <- rnorm(n)
informal <- strain > -0.3
flip <- function(x, p = 0.015) as.integer(xor(x, runif(n) < p))
debt <- cbind(family = flip(informal & runif(n) < .97), friends = flip(informal & runif(n) < .97),
              bank = flip(informal & strain > 0.5), card = NA, lender = flip(informal & excluded > 0.7),
              mortgage = as.integer(runif(n) < .25))
debt[, "card"] <- flip(debt[, "bank"] == 1 & runif(n) < .55)

ord <- function(M, stem) { d <- as.data.frame(lapply(seq_len(ncol(M)), function(j) likert(M[, j]))); names(d) <- paste0(stem, seq_len(ncol(M))); d }
num <- function(M, stem, digits = 1) { d <- as.data.frame(round(M, digits)); names(d) <- paste0(stem, seq_len(ncol(M))); d }
wizlab <- cbind(country = country, ord(calm, "calm"), ord(focus, "focus"), num(speed, "speed"), num(power, "power"),
                ord(wb[, 1:4], "body"), ord(wb[, 5:8], "mind"), ord(wb[, 9:12], "social"),
                num(know, "know", 0), num(habit, "habit", 0), setNames(as.data.frame(debt), paste0("debt_", colnames(debt))))
for (v in grep("^(know|habit|debt)", names(wizlab), value = TRUE)) wizlab[[v]] <- as.integer(wizlab[[v]])
# a little missingness in the questionnaire items, as in real surveys
for (v in grep("^(calm|focus|body|mind|social)", names(wizlab), value = TRUE)) wizlab[[v]][sample(n, 10)] <- NA
rownames(wizlab) <- NULL
save(wizlab, file = "data/wizlab.rda", compress = "xz")
if (requireNamespace("jmvReadWrite", quietly = TRUE))
  jmvReadWrite::write_omv(wizlab, "data/wizlab.omv", frcWrt = TRUE)
