# Extracted from test-audit-2026-08-03.R:295

# test -------------------------------------------------------------------------
skip_on_cran()
set.seed(307)
sizes <- c(rep(5L, 20), rep(800L, 20))
ab <- expand.grid(a = factor(1:8), b = factor(1:5))
u <- rnorm(40, 0, 0.15)
d <- do.call(rbind, lapply(1:40, function(j)
    data.frame(a = ab$a[j], b = ab$b[j], x = rnorm(sizes[j]), j = j)))
d$y <- 1 + 0.5 * d$x + u[d$j] + rnorm(nrow(d), 0, 1)
m <- suppressWarnings(fit_maihda(y ~ x + (1 | a:b), data = d[, c("a", "b", "x", "y")]))
raw <- maihda_re_normality_stat(lme4::ranef(m$model)$stratum[["(Intercept)"]])
expect_gt(raw$excess_kurtosis,
            maihda_adequacy_thresholds()$re_excess_kurtosis)
