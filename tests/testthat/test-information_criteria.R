# Information criteria (maihda_ic) and their integration into compare_maihda().

make_ic_data <- function(seed = 4242, n = 300) {
  set.seed(seed)
  d <- data.frame(
    g1 = sample(c("F", "M"), n, replace = TRUE),
    g2 = sample(c("A", "B", "C"), n, replace = TRUE),
    x = rnorm(n),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$g1, d$g2, drop = TRUE)
  u <- rnorm(nlevels(stratum), sd = 0.7)[stratum]
  d$y <- 1 + 0.4 * d$x + u + rnorm(n, sd = 0.5)
  # A binary outcome for the glmer path.
  d$bin <- rbinom(n, 1, plogis(-0.2 + u))
  d
}

test_that("maihda_ic on two lme4 Gaussian models reports AIC/BIC with ML refit and delta", {
  d <- make_ic_data()
  null_model <- fit_maihda(y ~ 1 + (1 | g1:g2), data = d)
  adj_model  <- fit_maihda(y ~ x + (1 | g1:g2), data = d)

  ic <- maihda_ic(null_model, adj_model, model_names = c("Null", "Adjusted"))
  expect_s3_class(ic, "maihda_ic")
  expect_equal(ic$model, c("Null", "Adjusted"))
  expect_true(all(c("AIC", "BIC", "df", "logLik", "n") %in% names(ic)))
  # All-NA Bayesian columns are dropped for an all-lme4 comparison.
  expect_false(any(c("WAIC", "LOOIC") %in% names(ic)))
  expect_true(all(is.finite(ic$AIC)) && all(is.finite(ic$BIC)))

  # REML lmer fits are refitted with ML so AIC/BIC compare across fixed effects.
  expect_true(all(grepl("refit from REML", ic$estimator)))

  # delta is the gap from the best model on AIC; the best row is 0.
  expect_true("delta" %in% names(ic))
  expect_equal(min(ic$delta), 0)
  expect_identical(attr(ic, "ic_primary"), "AIC")
})

test_that("a single lme4 model reports the REML estimator and no delta column", {
  d <- make_ic_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d)

  ic <- maihda_ic(m)
  expect_equal(nrow(ic), 1L)
  expect_identical(ic$estimator, "REML")   # lmer default, reported as fitted
  expect_true(is.finite(ic$AIC) && is.finite(ic$BIC))
  expect_false("delta" %in% names(ic))      # delta only with >1 model
})

test_that("a binomial glmer model reports the ML estimator", {
  d <- make_ic_data()
  m <- fit_maihda(bin ~ x + (1 | g1:g2), data = d, family = "binomial")

  ic <- maihda_ic(m)
  expect_identical(ic$estimator, "ML")       # glmer is ML, no REML refit needed
  expect_true(is.finite(ic$AIC))
})

test_that("maihda_ic expands a maihda_analysis into Null and Adjusted rows", {
  d <- make_ic_data()
  a <- suppressMessages(maihda(y ~ x + g1 + g2 + (1 | g1:g2), data = d))

  ic <- maihda_ic(a)
  expect_equal(nrow(ic), 2L)
  expect_match(ic$model[1], "Null")
  expect_match(ic$model[2], "Adjusted")
  expect_true(all(is.finite(ic$AIC)))
})

test_that("maihda_ic validates its arguments", {
  d <- make_ic_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d)

  expect_error(maihda_ic(), "at least one")
  expect_error(maihda_ic(42), "maihda_model or maihda_analysis")
  expect_error(maihda_ic(m, m, model_names = "only-one"), "must match")
})

test_that("compare_maihda appends information-criteria columns when ic = TRUE", {
  d <- make_ic_data()
  m1 <- fit_maihda(y ~ 1 + (1 | g1:g2), data = d)
  m2 <- fit_maihda(y ~ x + (1 | g1:g2), data = d)

  cmp <- compare_maihda(m1, m2, model_names = c("Null", "Adjusted"))
  expect_s3_class(cmp, "maihda_comparison")
  expect_true(all(c("AIC", "BIC") %in% names(cmp)))
  expect_true(all(is.finite(cmp$AIC)))
  # The VPC columns and plotting are unaffected by the extra IC columns.
  expect_true("vpc" %in% names(cmp))
  expect_s3_class(plot(cmp), "ggplot")

  # ic = FALSE restores the lean VPC-only table.
  cmp_lean <- compare_maihda(m1, m2, ic = FALSE)
  expect_false(any(c("AIC", "BIC", "WAIC", "LOOIC") %in% names(cmp_lean)))

  expect_error(compare_maihda(m1, m2, ic = NA), "must be TRUE or FALSE")
})

test_that("maihda_ic reports AIC/BIC for the ordinal (clmm) engine", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  set.seed(909)
  n <- 600
  d <- data.frame(
    g1 = sample(c("F", "M"), n, replace = TRUE),
    g2 = sample(c("A", "B", "C"), n, replace = TRUE),
    x = rnorm(n),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$g1, d$g2, drop = TRUE)
  u <- rnorm(nlevels(stratum), sd = 0.6)[stratum]
  lat <- u + 0.3 * d$x + rlogis(n)
  d$y <- factor(cut(lat, c(-Inf, -1, 0.5, 2, Inf), labels = 1:4), ordered = TRUE)

  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, family = "ordinal")))
  ic <- maihda_ic(m)
  expect_identical(ic$estimator, "ML")
  expect_true(is.finite(ic$AIC) && is.finite(ic$BIC))
})

test_that("design-weighted (wemix) fits report NA criteria without warning", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  set.seed(515)
  n <- 400
  d <- data.frame(
    g1 = sample(c("F", "M"), n, replace = TRUE),
    g2 = sample(c("A", "B", "C"), n, replace = TRUE),
    x = rnorm(n),
    w = runif(n, 0.5, 2),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$g1, d$g2, drop = TRUE)
  u <- rnorm(nlevels(stratum), sd = 0.5)[stratum]
  d$y <- 1 + 0.3 * d$x + u + rnorm(n, sd = 0.5)

  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | g1:g2), data = d,
               engine = "wemix", sampling_weights = "w")))
  ic <- maihda_ic(m)
  expect_identical(ic$estimator, "pseudo-ML (weighted)")
  # AIC/BIC are NA for pseudo-likelihood, so the all-NA columns are dropped.
  expect_false(any(c("AIC", "BIC") %in% names(ic)))
  # n must still be the analytic sample size: WeMixResults has no nobs() method,
  # so maihda_ic falls back to nrow(model$data) rather than reporting NA.
  expect_true(is.finite(ic$n))
  expect_identical(as.integer(ic$n), nrow(m$data))

  # compare_maihda(ic = TRUE) on a wemix pair must not add a warning.
  expect_no_warning(suppressMessages(compare_maihda(m, m)))
})

test_that("maihda_ic_spans_scales flags mixed likelihood/Bayesian IC columns", {
  expect_true(maihda_ic_spans_scales(c("AIC", "BIC", "WAIC", "LOOIC")))
  expect_true(maihda_ic_spans_scales(c("AIC", "LOOIC")))   # one of each is enough
  expect_true(maihda_ic_spans_scales(c("BIC", "WAIC")))
  expect_false(maihda_ic_spans_scales(c("AIC", "BIC")))    # likelihood only
  expect_false(maihda_ic_spans_scales(c("WAIC", "LOOIC"))) # Bayesian only
  expect_false(maihda_ic_spans_scales(character(0)))
})

test_that("compare_maihda warns when appended IC columns mix likelihood and Bayesian scales", {
  d <- make_ic_data()
  # Same outcome, family, sample and strata, so the differences block stays
  # silent and the only thing left to flag is the mixed IC scale.
  m_lik   <- fit_maihda(y ~ x + (1 | g1:g2), data = d)
  m_bayes <- fit_maihda(y ~ 1 + (1 | g1:g2), data = d)

  # Stand in for a brms fit on the second model: its IC row carries WAIC/LOOIC
  # instead of AIC/BIC, so the appended columns span both scales -- without any
  # actual Stan toolchain. compare_maihda() reads the IC rows in argument order.
  calls <- 0L
  local_mocked_bindings(
    maihda_ic_one = function(model, ml = FALSE) {
      calls <<- calls + 1L
      if (calls == 1L) {
        data.frame(n = 300L, estimator = "ML (refit from REML)", df = 4,
                   logLik = -400, AIC = 808, BIC = 823,
                   WAIC = NA_real_, LOOIC = NA_real_, stringsAsFactors = FALSE)
      } else {
        data.frame(n = 300L, estimator = "Bayesian", df = NA_real_,
                   logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                   WAIC = 805, LOOIC = 807, stringsAsFactors = FALSE)
      }
    }
  )

  expect_warning(
    cmp <- compare_maihda(m_lik, m_bayes, model_names = c("lme4", "brms")),
    "mix scales|different scales"
  )
  # The table still carries all four criteria side by side -- the warning flags
  # the mixing rather than dropping columns.
  expect_true(all(c("AIC", "BIC", "WAIC", "LOOIC") %in% names(cmp)))
})

test_that("an all-lme4 comparison appends AIC/BIC without a scale-mixing warning", {
  d <- make_ic_data()
  m1 <- fit_maihda(y ~ x + (1 | g1:g2), data = d)
  m2 <- fit_maihda(y ~ 1 + (1 | g1:g2), data = d)

  expect_no_warning(cmp <- compare_maihda(m1, m2))
  expect_true(all(c("AIC", "BIC") %in% names(cmp)))
  expect_false(any(c("WAIC", "LOOIC") %in% names(cmp)))
})
