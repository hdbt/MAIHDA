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
