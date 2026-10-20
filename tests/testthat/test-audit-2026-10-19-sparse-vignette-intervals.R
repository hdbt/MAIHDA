# Audit 2026-10-19 (sparse-vignette-intervals): three claims of
# vignette("bayesian_sparse_maihda").
# - F1, CONFIRMED: it said the lme4 shares carry no interval ("neither estimate carries
#   an interval", "provides a single value with no uncertainty", "at the boundary,
#   reports it without a standard error"). maihda(bootstrap = TRUE) gives them a
#   parametric-bootstrap interval (bootstrap_cc()). On the vignette's fits, 1,000 draws
#   after set.seed(1): Gaussian 23.4% [0.0%, 89.2%], binary 3.1% [0.0%, 100.0%] (brms
#   under the vignette's priors: 13% [0.6%, 53%] and 5% [0.0%, 41%]). The draws are
#   simulated from the fit, so none has a gender effect; every Gaussian refit puts
#   gender's variance below 1e-4 (1,000 of 1,000), the binary refits in 74%. In 200 data
#   sets from data-raw/maihda_sparse_data.R (200 draws each) the Gaussian interval
#   covered the realised share (the mean-square split of the drawn strata) in all 200,
#   median interval [0.000, 0.918], the 145 fits with a dimension's variance at zero
#   included; the binary one (60 sets, 100 draws, fitted on per-stratum counts) in 59,
#   median [0, 1], the miss a fit with the interaction variance at zero ([0, 0.20]).
#   The vignette now says the default output has no interval, what bootstrap = TRUE
#   gives here, and that its draws come from the fit.
# - F2, CONFIRMED: "the prior regularises the between-stratum standard deviations away
#   from the boundary (so the fit is not singular)". The normal(0, 0.5) prior has its
#   mode at zero: gender's SD has posterior median 0.19 under it, 0.47 under brms's
#   default and 0.92 under a flat prior on (0, 10), and its posterior mode is 0 under all
#   three (exact Gaussian marginal likelihood). The vignette now says the SDs are off
#   zero because each is a posterior median.
# - F3, CONFIRMED: the lme4 chunk printed lme4's singular flags. The Gaussian fit's REML
#   optimum has gender's variance exactly at zero (criterion about 20.4 theta^2), but
#   nloptwrap stops anywhere up to theta about 1.5e-3 and lme4 flags only below 1e-4:
#   flagged in 29 of 40 row orders and 21 of 40 one-ulp changes of y, share 0.234
#   throughout (0.2336 to 0.2344). The chunk now prints the variances, and the sentence after
#   it opens with them (its pin is in test-audit-2026-10-19-sparse-vignette-lme4-boundary.R).
#
# The bootstrap block is gated by skip_slow_test(), which restores the random-number
# state; no other block draws a random number.

sparse_quiet <- function(expr) suppressMessages(suppressWarnings(expr))

sparse_fit <- function(outcome, family = "gaussian", ...) {
  data(maihda_sparse_data, envir = environment())
  f <- stats::as.formula(paste(outcome, "~ 1 + (1 | gender:ethnicity:education:age_group)"))
  sparse_quiet(maihda(f, data = maihda_sparse_data, decomposition = "crossed-dimensions",
                      engine = "lme4", family = family, ...))
}

test_that("bootstrap = TRUE gives the sparse vignette's lme4 shares a wide interval", {
  skip_on_cran()
  skip_slow_test()
  # The vignette: from 0% to about 90% on the Gaussian outcome, over the whole range on
  # the binary one. Both are 1,000-draw intervals (the default n_boot): Gaussian upper
  # ends 0.89 to 0.95 over seeds 1 to 5. The block uses fewer draws, so its bounds are
  # looser. Gaussian, 200 draws: upper ends 0.88 to 0.92 over seeds 1 to 5, 0.905 to
  # 0.925 for seed 1's draws refitted on permuted rows. Binary: about 4% of the draws
  # have every dimension's variance at zero and the interaction's not, a share of 1,
  # which puts the upper end of a 1,000-draw interval at 1 or near it in 99% of runs;
  # with 100 draws the upper end falls below 0.8 for 5 to 9% of seeds (seed 1: 0.87,
  # and 0.87 to 0.95 on permuted rows).
  set.seed(1)
  m_g <- sparse_fit("y", bootstrap = TRUE, n_boot = 200)
  dec <- m_g$decomposition
  ci <- dec$interaction_share_ci
  expect_true(isTRUE(dec$bootstrap))
  expect_true(all(is.finite(ci)))
  expect_true(ci[1] <= dec$interaction_share && dec$interaction_share <= ci[2])
  expect_lt(ci[1], 0.01)
  expect_gt(ci[2], 0.8)
  expect_lt(ci[2], 0.98)

  # The draws are simulated from the fit, whose gender variance is at zero: no draw has
  # a gender effect, and no refit finds one.
  expect_lt(dec$per_dim[["gender"]], 1e-4)
  model <- m_g$model$model
  sim <- maihda_simulate_lme4(model, nsim = 20)
  refit_gender <- vapply(sim, function(y) {
    maihda_random_variances_lme4(sparse_quiet(maihda_refit_draw(model, y)))[["gender"]]
  }, numeric(1))
  expect_true(all(refit_gender < 1e-4))

  set.seed(1)
  m_b <- sparse_fit("event", family = "binomial", bootstrap = TRUE, n_boot = 100)
  ci <- m_b$decomposition$interaction_share_ci
  expect_true(all(is.finite(ci)))
  expect_lt(ci[1], 0.01)
  expect_gt(ci[2], 0.7)
})

test_that("the normal(0, 0.5) prior pulls gender's standard deviation toward zero", {
  skip_on_cran()
  # The exact marginal likelihood of the vignette's Gaussian crossed-dimensions model,
  # mu integrated out under a flat prior, from the cell means, the cell sizes and the
  # within-cell sum of squares; along gender's standard deviation, the other variances
  # at their REML values.
  data(maihda_sparse_data, envir = environment())
  d <- maihda_sparse_data
  dims <- c("gender", "ethnicity", "education", "age_group")
  key <- interaction(d[dims], drop = TRUE, sep = ":")
  n <- as.numeric(table(key))
  ybar <- as.numeric(tapply(d$y, key, mean))
  cells <- d[match(levels(key), as.character(key)), dims]
  within_ss <- sum((d$y - stats::ave(d$y, key))^2)
  J <- length(n)
  ZZ <- lapply(dims, function(v) tcrossprod(stats::model.matrix(~ cells[[v]] - 1)))
  log_lik <- function(sd_dim, sd_int, sigma) {
    V <- Reduce(`+`, Map(function(z, s) s^2 * z, ZZ, sd_dim)) +
      diag(sd_int^2 + sigma^2 / n, J)
    R <- chol(V)
    a <- backsolve(R, rep(1, J), transpose = TRUE)
    b <- backsolve(R, ybar, transpose = TRUE)
    -0.5 * ((nrow(d) - J) * log(sigma^2) + within_ss / sigma^2) - sum(log(diag(R))) -
      0.5 * log(sum(a^2)) - 0.5 * (sum(b^2) - sum(a * b)^2 / sum(a^2))
  }
  dec <- sparse_fit("y")$decomposition
  sd_dim <- sqrt(dec$per_dim[dims])
  grid <- seq(0, 4, by = 0.001)
  ll <- vapply(grid, function(s) {
    log_lik(replace(sd_dim, 1L, s), sqrt(dec$interaction_var), sqrt(dec$within_var))
  }, numeric(1))
  posterior <- function(log_prior) {
    w <- exp(ll + log_prior - max(ll + log_prior))
    w / sum(w)
  }
  median_of <- function(w) grid[which(cumsum(w) >= 0.5)[1]]
  half_normal <- posterior(stats::dnorm(grid, 0, 0.5, log = TRUE))
  flat <- posterior(0)
  # The likelihood and the posterior under the prior peak at zero ...
  expect_identical(grid[which.max(ll)], 0)
  expect_identical(grid[which.max(half_normal)], 0)
  # ... the posterior medians are above it, and the prior's median is the lower.
  expect_gt(median_of(half_normal), 0.1)
  expect_lt(median_of(half_normal), median_of(flat) - 0.2)
})

test_that("the sparse vignette's text matches its intervals, prior and variances", {
  skip_on_cran()
  root <- testthat::test_path("..", "..")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  p <- file.path(root, "vignettes", "bayesian_sparse_maihda.Rmd")
  skip_if_not(file.exists(p), "the package tree is not available")
  v <- gsub("[[:space:]]+", " ", paste(readLines(p, warn = FALSE), collapse = " "))
  for (old in c("neither estimate carries an interval",
                "unwarranted precision",
                "has no way to express",
                "a single value with no uncertainty",
                "reports it without a standard error",
                "as if it were certain",
                "away from the boundary",
                "so the fit is not singular",
                "m_g$model$diagnostics$singular")) {
    expect_false(grepl(old, v, fixed = TRUE), info = old)
  }
  for (new in c(
    "Neither estimate comes with an interval in this output.",
    paste0("Asked for one (`bootstrap = TRUE`), `maihda()` simulates outcomes from the ",
           "fitted model and refits it; the interval for the interaction share then runs ",
           "from 0% to about 90% on the Gaussian outcome and over the whole range on the ",
           "binary one."),
    "The draws are simulated from the fit itself, so none of them has a gender effect.",
    paste0("No standard deviation is reported at zero, since each is a posterior median, ",
           "which is above zero whatever the prior; the `normal(0, 0.5)` prior has its ",
           "mode at zero and pulls the standard deviations toward it."),
    "Maximum likelihood's default output is a single point;",
    "whereas `lme4` reports a single value unless it is bootstrapped.",
    paste0("commits to a single variance decomposition and gives no interval for it ",
           "unless one is bootstrapped from that decomposition."),
    "c(m_g$decomposition$per_dim, interaction = m_g$decomposition$interaction_var)",
    "c(m_b$decomposition$per_dim, interaction = m_b$decomposition$interaction_var)")) {
    expect_true(grepl(new, v, fixed = TRUE), info = new)
  }
})
