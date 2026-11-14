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
  jmvReadWrite::write_omv(wizdata, "data/wizdata.omv")
