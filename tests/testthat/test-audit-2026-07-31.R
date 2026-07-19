# Regression tests for the 2026-07-31 audit findings.
#
#   1 [High] Any numeric response with exactly two distinct values was classified
#            Bernoulli (maihda_is_binary_vector) and recoded to 0/1, even a valid
#            PROPORTION response (successes / trials) supplied with weights = trial
#            counts -- the aggregated-binomial form documented at
#            ?maihda_discriminatory_accuracy. Proportions 0.25 and 0.5 were recoded
#            to 0 and 1, silently fitting a different (wrong) model (intercept and
#            AUC diverging from the equivalent cbind(success, failure) fit). Numeric
#            binary detection now requires whole-number values, so a proportion in
#            (0, 1) stays an aggregated binomial() fit.
#   2 [High] maihda_bootstrap_ci() required only an absolute floor of 10 finite
#            draws regardless of n_boot, so a 95% interval could be returned from 10
#            successes out of 1000 (99% failure) with only a warning. It now also
#            requires a MAJORITY of the eligible draws to succeed; legitimately
#            excluded draws (e.g. PCV boundary draws) are passed as n_excluded and do
#            not count against the fraction.
#   3 [Med]  The crossed-dimensions MOR solved a pair-mixture CDF for its median via
#            uniroot on [eps, hi]. When the zero-variance pairs carry >= 0.5 of the
#            weight the median |difference| is exactly 0 (MOR = 1), but both
#            endpoints give cdf >= 0.5 so uniroot could not bracket the root and the
#            MOR was reported NA. It now returns 1 whenever the zero-variance pair
#            mass reaches 0.5.
#   4 [Med]  The parametric-bootstrap loops caught thrown errors only, so a refit
#            returned with an optimiser that did not converge still contributed to
#            the interval and to n_boot_ok, which therefore implied a convergence
#            that was never checked. Contributing draws whose optimiser failed
#            (optinfo$conv$opt != 0) are now counted, reported (warning +
#            n_boot_nonconverged), and retained. lme4's post-hoc relative-gradient
#            flag -- a frequent false positive on simulated refits -- is not used.

# ---- Finding 1: proportion response is not recoded to Bernoulli --------------

test_that("a numeric two-level PROPORTION response is not treated as Bernoulli", {
  d <- data.frame(stratum = factor(rep(seq_len(12), each = 12)))
  # exactly two distinct proportion values, both strictly inside (0, 1)
  d$p <- rep(c(0.25, 0.5), length.out = nrow(d))
  d$n <- 20L

  expect_false(MAIHDA:::maihda_response_is_binary(p ~ (1 | stratum), d, weights = d$n))
  # left un-recoded by maihda_prepare_binomial_response()
  prep <- MAIHDA:::maihda_prepare_binomial_response(d, p ~ (1 | stratum), weights = d$n)
  expect_identical(sort(unique(prep$p)), c(0.25, 0.5))
  expect_null(attr(prep, "response_recoding"))

  # genuine integer codings are still Bernoulli (0/1 and a 1/2 coding)
  expect_true(MAIHDA:::maihda_is_binary_vector(c(0, 1, 1, 0)))
  expect_true(MAIHDA:::maihda_is_binary_vector(c(1, 2, 2, 1)))
  # a response of exactly 0/1 (all-success / all-failure rows) stays Bernoulli
  expect_true(MAIHDA:::maihda_is_binary_vector(c(0, 1, 0, 1)))
  # the proportion is not
  expect_false(MAIHDA:::maihda_is_binary_vector(c(0.25, 0.5)))
})

test_that("proportion+weights and cbind(success, failure) fits now agree", {
  skip_on_cran()
  set.seed(731)
  strat <- factor(rep(seq_len(12), each = 12))
  succ <- sample(c(5L, 10L), length(strat), replace = TRUE)
  d <- data.frame(stratum = strat, p = succ / 20, n = 20L,
                  s = succ, f = 20L - succ)

  m_prop <- suppressMessages(
    fit_maihda(p ~ (1 | stratum), data = d, family = "binomial", weights = n))
  m_cbind <- suppressMessages(
    fit_maihda(cbind(s, f) ~ (1 | stratum), data = d, family = "binomial"))

  # Equivalent data -> identical fixed-effect intercept (previously 0.06 vs -0.50)
  expect_equal(unname(lme4::fixef(m_prop$model))[1],
               unname(lme4::fixef(m_cbind$model))[1], tolerance = 1e-4)
  # ... and identical discriminatory accuracy
  expect_equal(maihda_discriminatory_accuracy(m_prop)$auc,
               maihda_discriminatory_accuracy(m_cbind)$auc, tolerance = 1e-4)
})

# ---- Finding 2: bootstrap needs a majority of eligible draws ----------------

test_that("maihda_bootstrap_ci refuses an interval from a small success fraction", {
  # 10 finite draws of 1000 requested: was returned with only a warning.
  vals <- c(seq(1, 10), rep(NA_real_, 990))
  expect_error(
    MAIHDA:::maihda_bootstrap_ci(vals, n_boot = 1000, conf_level = 0.95, what = "VPC"),
    "must succeed"
  )
  # Legitimately excluded draws do not count against the fraction: 30 successes of
  # 50 eligible (100 requested, 50 excluded) is a majority and is allowed.
  ci <- suppressWarnings(
    MAIHDA:::maihda_bootstrap_ci(rep(0.5, 30), n_boot = 100, conf_level = 0.95,
                                 n_excluded = 50))
  expect_length(ci, 2)
  # Excluding so many that fewer than half of the remainder succeed still errors.
  expect_error(
    MAIHDA:::maihda_bootstrap_ci(rep(0.5, 20), n_boot = 100, conf_level = 0.95,
                                 n_excluded = 10),
    "must succeed"
  )
})

# ---- Finding 3: crossed MOR is 1, not NA, when heterogeneity is a minority ---

test_that("crossed MOR returns 1 when zero-variance pairs carry at least half the mass", {
  # Ragged six-cell design A x {1..5} + B x 1 with dim-2 and interaction variance 0:
  # 10 of 15 pairs (differing only in the zero-variance dimension) have zero contrast
  # variance, so the median |difference| is exactly 0 and MOR = exp(0) = 1.
  parts_bug <- list(
    interaction = 0,
    pattern = matrix(c(0, 1,  1, 0,  1, 1), nrow = 3, byrow = TRUE),
    dims = c(0.8, 0), weight = c(10, 1, 4))
  expect_identical(MAIHDA:::maihda_mor_crossed_from_parts(parts_bug), 1)

  # Exactly half the mass at zero variance still gives MOR 1.
  parts_half <- list(
    interaction = 0,
    pattern = matrix(c(0, 1,  1, 0), nrow = 2, byrow = TRUE),
    dims = c(0.8, 0), weight = c(1, 1))
  expect_identical(MAIHDA:::maihda_mor_crossed_from_parts(parts_half), 1)

  # A minority (0.4) of zero-variance mass leaves a finite MOR > 1.
  parts_ok <- list(
    interaction = 0,
    pattern = matrix(c(0, 1,  1, 0,  1, 1), nrow = 3, byrow = TRUE),
    dims = c(0.8, 0), weight = c(6, 3, 6))
  mor_ok <- MAIHDA:::maihda_mor_crossed_from_parts(parts_ok)
  expect_true(is.finite(mor_ok) && mor_ok > 1)

  # All pairs at zero variance is still the trivial MOR 1.
  parts_zero <- list(
    interaction = 0,
    pattern = matrix(c(0, 1), nrow = 1), dims = c(0.8, 0), weight = 3)
  expect_identical(MAIHDA:::maihda_mor_crossed_from_parts(parts_zero), 1)
})

# ---- Finding 4: bootstrap non-convergence is checked, not assumed ------------

test_that("maihda_lme4_optimizer_failed flags a genuine optimiser failure only", {
  set.seed(741)
  d <- data.frame(stratum = factor(rep(seq_len(10), each = 15)))
  d$y <- rnorm(150) + rnorm(10, sd = 0.6)[d$stratum]
  m <- suppressMessages(fit_maihda(y ~ (1 | stratum), data = d))

  expect_false(MAIHDA:::maihda_lme4_optimizer_failed(m$model))
  bad <- m$model
  bad@optinfo$conv$opt <- 1L        # optimiser return code != 0
  expect_true(MAIHDA:::maihda_lme4_optimizer_failed(bad))
  expect_false(MAIHDA:::maihda_lme4_optimizer_failed("not a model"))
})

test_that("PCV/VPC bootstrap results report a non-convergence count", {
  skip_on_cran()
  set.seed(742)
  d <- data.frame(stratum = factor(rep(seq_len(10), each = 20)), x = rnorm(200))
  d$y <- 1 + 0.3 * d$x + rnorm(10, sd = 0.7)[d$stratum] + rnorm(200, sd = 0.4)
  m0 <- suppressMessages(fit_maihda(y ~ (1 | stratum), data = d))
  m1 <- suppressMessages(fit_maihda(y ~ x + (1 | stratum), data = d))

  pcv <- suppressWarnings(suppressMessages(
    calculate_pcv(m0, m1, bootstrap = TRUE, n_boot = 40)))
  expect_true(is.numeric(pcv$n_boot_nonconverged))
  expect_true(pcv$n_boot_nonconverged >= 0 &&
                pcv$n_boot_nonconverged <= pcv$n_boot_ok)

  s <- suppressWarnings(suppressMessages(summary(m1, bootstrap = TRUE, n_boot = 40)))
  expect_true(is.numeric(s$vpc$n_boot_nonconverged))
  expect_true(s$vpc$n_boot_nonconverged >= 0)
})
