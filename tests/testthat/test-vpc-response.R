# A binomial MAIHDA fit on synthetic data with a controllable between-stratum signal.
maihda_vpcr_fit <- function(seed = 7, n = 1500, sd_u = 1.2) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  lp <- stats::rnorm(nlevels(sk), sd = sd_u)[sk]
  d$y <- stats::rbinom(n, 1, stats::plogis(lp))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "binomial")
  ))
}

test_that("maihda_vpc_response returns a VPC in [0, 1] and is seed-reproducible", {
  m <- maihda_vpcr_fit()
  r1 <- maihda_vpc_response(m, n_sim = 20000, seed = 42)
  r2 <- maihda_vpc_response(m, n_sim = 20000, seed = 42)

  expect_s3_class(r1, "maihda_vpc_response")
  expect_identical(r1$scale, "response")
  expect_true(r1$estimate >= 0 && r1$estimate <= 1)
  expect_equal(r1$estimate, r2$estimate)  # same seed -> identical draws -> identical VPC
})

test_that("maihda_vpc_response matches an independent Goldstein/Browne/Rasbash simulation", {
  m <- maihda_vpcr_fit()
  r <- maihda_vpc_response(m, n_sim = 50000, seed = 1)

  # Reimplement the computation independently with the same seed/draws.
  s2 <- MAIHDA:::extract_between_variance(m)
  lp <- mean(stats::predict(m$model, re.form = NA, type = "link"))
  set.seed(1)
  u <- stats::rnorm(50000, 0, sqrt(s2))
  p <- stats::plogis(lp + u)
  ref <- stats::var(p) / (stats::var(p) + mean(p * (1 - p)))

  expect_equal(r$estimate, ref, tolerance = 1e-6)
})

test_that("response-scale VPC increases with between-stratum variance", {
  low  <- maihda_vpcr_fit(seed = 3, sd_u = 0.3)
  high <- maihda_vpcr_fit(seed = 3, sd_u = 1.6)
  v_low  <- maihda_vpc_response(low,  n_sim = 50000, seed = 5)$estimate
  v_high <- maihda_vpc_response(high, n_sim = 50000, seed = 5)$estimate
  expect_gt(v_high, v_low)
})

test_that("maihda_vpc_response rejects non-binomial models and invalid n_sim", {
  g <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(120), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    family = "gaussian")))
  expect_error(maihda_vpc_response(g$model), "binomial")

  m <- maihda_vpcr_fit()
  expect_error(maihda_vpc_response(m, n_sim = 10), "n_sim")
  # A fractional count is rejected (rnorm() would silently truncate it while the
  # recorded n_sim kept the fraction).
  expect_error(maihda_vpc_response(m, n_sim = 100.5), "whole number")

  # A valid whole number is stored as an integer, matching the realised draw count.
  r <- maihda_vpc_response(m, n_sim = 100, seed = 1)
  expect_identical(r$n_sim, 100L)
})

test_that("response-scale VPC integrates over non-stratum random effects", {
  # Regression: for a site + stratum model the simulation drew ONLY the stratum
  # effect, so the site variance appeared in neither the numerator's scope nor
  # the denominator -- overstating the stratum share (audit repro: 0.0490
  # reported vs 0.0126 from a coherent nested Monte Carlo).
  skip_on_cran()
  set.seed(707)
  d <- expand.grid(stratum = factor(1:8), site = factor(1:10), rep = 1:25)
  su <- stats::rnorm(8, sd = 0.4)[d$stratum]
  so <- stats::rnorm(10, sd = 1.2)[d$site]
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.5 + su + so))
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | site) + (1 | stratum), data = d, family = "binomial")))

  v <- maihda_vpc_response(m, n_sim = 20000, seed = 42)
  expect_gt(v$var_other, 0)
  expect_true(v$estimate >= 0 && v$estimate <= 1)

  # Independent brute-force nested Monte Carlo of the stratum share.
  vc <- lme4::VarCorr(m$model)
  vs <- as.numeric(as.matrix(vc$stratum)[1, 1])
  vo <- as.numeric(as.matrix(vc$site)[1, 1])
  eta0 <- mean(stats::predict(m$model, re.form = NA, type = "link"))
  set.seed(99)
  us <- stats::rnorm(4000, 0, sqrt(vs))
  uo <- stats::rnorm(2000, 0, sqrt(vo))
  P <- stats::plogis(outer(us, uo, `+`) + eta0)
  truth <- stats::var(rowMeans(P)) /
    (stats::var(as.vector(P)) + mean(P * (1 - P)))
  expect_equal(v$estimate, truth, tolerance = 0.1)

  # The old stratum-only computation is strictly larger (the defect).
  p_only <- stats::plogis(eta0 + us)
  old <- stats::var(p_only) / (stats::var(p_only) + mean(p_only * (1 - p_only)))
  expect_gt(old, v$estimate)

  # print() explains the integration.
  expect_output(print(v), "non-stratum random effects")

  # A single-stratum fit is unchanged: var_other 0, same seeded value as the
  # stratum-only formula.
  m0 <- maihda_vpcr_fit()
  v0 <- maihda_vpc_response(m0, n_sim = 20000, seed = 42)
  expect_identical(v0$var_other, 0)
  s2 <- MAIHDA:::extract_between_variance(m0)
  lp <- mean(stats::predict(m0$model, re.form = NA, type = "link"))
  set.seed(42)
  u <- stats::rnorm(20000, 0, sqrt(s2))
  p <- stats::plogis(lp + u)
  expect_equal(v0$estimate,
               stats::var(p) / (stats::var(p) + mean(p * (1 - p))),
               tolerance = 1e-12)
})
