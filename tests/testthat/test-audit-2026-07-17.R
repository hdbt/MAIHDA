# Regression tests for the 2026-07-13 (High/Medium/Medium-low) audit findings.
#
#   #1 High    the external offset= was forwarded to lme4 but omitted from the
#              analytic-sample masks, so a row lme4 later drops for offset = NA was
#              still seen by automatic family detection (flipping e.g. gaussian for
#              binomial) and by longitudinal time-centering (mis-anchoring it).
#   #2 Medium  a per-group PCV decomposition that failed (adjusted fit / calculate_pcv
#              error) was silently turned into pcv = NA with status still "ok" and no
#              warning; the failure is now surfaced (pcv_status + aggregated warning).
#   #3 Med-low direct maihda_ic() calls did not guard against a table mixing likelihood
#              AIC/BIC with Bayesian WAIC/LOOIC, so an lme4-vs-brms comparison ranked
#              on AIC with the Bayesian row left an NA delta and no warning.

# ---- #1 High: external offset= folded into the analytic masks ---------------

test_that("maihda_row_mask drops rows whose offset is NA", {
  d <- data.frame(y = 1:5, x = 1:5)
  expect_equal(MAIHDA:::maihda_row_mask(d, offset = c(1, NA, 3, 4, NA)),
               c(TRUE, FALSE, TRUE, TRUE, FALSE))
  # A length-mismatched offset (e.g. a scalar) carries no per-row NA semantics and
  # is ignored, mirroring the weights guard.
  expect_equal(MAIHDA:::maihda_row_mask(d, offset = c(1, NA)), rep(TRUE, 5))
  # Offset-NA and weight-NA drops combine.
  expect_equal(
    MAIHDA:::maihda_row_mask(d, weights = c(1, 1, NA, 1, 1), offset = c(1, NA, 1, 1, 1)),
    c(TRUE, FALSE, FALSE, TRUE, TRUE))
})

test_that("an offset=NA row does not flip automatic family detection", {
  set.seed(101)
  n <- 400
  d <- data.frame(
    g1 = sample(c("F", "M"), n, replace = TRUE),
    g2 = sample(c("A", "B", "C"), n, replace = TRUE),
    x  = rnorm(n),
    stringsAsFactors = FALSE
  )
  d$y   <- rbinom(n, 1, plogis(-0.2 + 0.6 * d$x))   # genuine 0/1 outcome
  d$off <- rnorm(n, 0, 0.3)
  d$y[1]   <- 2L                                     # out-of-sample value ...
  d$off[1] <- NA_real_                               # ... on the row lme4 drops

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, offset = off)))
  # Detection now sees the analytic sample (the offset-NA row removed), so the
  # outcome is recognised as binary -- not gaussian because of the stray 2.
  expect_identical(m$family$family, "binomial")
  expect_true(inherits(m$model, "glmerMod"))
  expect_equal(MAIHDA:::maihda_nobs(m$model), 399)

  # Prefiltering the offset-NA row by hand must give the identical family.
  mb <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d[!is.na(d$off), ], offset = off)))
  expect_identical(mb$family$family, "binomial")
})

test_that("an offset=NA baseline occasion does not mis-anchor longitudinal centering", {
  set.seed(202)
  nid <- 120
  g1 <- sample(c("F", "M"), nid, replace = TRUE)
  g2 <- sample(c("A", "B", "C"), nid, replace = TRUE)
  long <- do.call(rbind, lapply(seq_len(nid), function(i)
    data.frame(id = i, g1 = g1[i], g2 = g2[i], time = c(10, 11, 12),
               stringsAsFactors = FALSE)))
  long$y   <- rnorm(nrow(long))
  long$off <- rnorm(nrow(long), 0, 0.3)
  # A stray time = 0 occasion with a missing offset (dropped by lme4). The fitted
  # times start at 10, so the centre must be 10, not 0.
  extra <- data.frame(id = 1, g1 = g1[1], g2 = g2[1], time = 0,
                      y = 0.4, off = NA_real_, stringsAsFactors = FALSE)
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ 1 + (1 | g1:g2), data = rbind(long, extra),
               id = "id", time = "time", offset = off)))
  expect_equal(m$longitudinal_info$time_center, 10)
})

# ---- #2 Medium: per-group PCV failure is surfaced, not silent ---------------

test_that("a failed per-group PCV is flagged (pcv_status + warning), not silent", {
  set.seed(303)
  mk <- function(country, genders, nn) data.frame(
    country = country,
    gender  = sample(genders, nn, replace = TRUE),
    ses     = sample(c("low", "mid", "high"), nn, replace = TRUE),
    stringsAsFactors = FALSE)
  # g1 has a single gender level, so its adjusted model's gender main effect is a
  # one-level factor lmer rejects -> the adjusted fit (hence the PCV) fails.
  dg <- rbind(mk("g1", "F", 160), mk("g2", c("F", "M"), 240))
  strat <- interaction(dg$gender, dg$ses, drop = TRUE)
  dg$y <- 1 + rnorm(nlevels(strat), 0, 0.5)[strat] + rnorm(nrow(dg), 0, 1)

  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(compare_maihda_groups(y ~ (1 | gender:ses), data = dg,
                                           group = "country")),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })

  expect_true("pcv_status" %in% names(res))
  g1 <- res[res$group == "g1", ]
  g2 <- res[res$group == "g2", ]
  expect_identical(g1$pcv_status, "failed")
  expect_true(is.na(g1$pcv))
  # The group's own status stays "ok": its NULL VPC model succeeded; only the
  # decomposition failed, and that is what pcv_status records.
  expect_identical(g1$status, "ok")
  expect_identical(g2$pcv_status, "ok")
  expect_false(is.na(g2$pcv))
  # The failure is named in an aggregated warning rather than swallowed.
  expect_true(any(grepl("per-group PCV decomposition failed", warns)))
  expect_true(any(grepl("\\bg1\\b", warns)))
})

test_that("pcv_status is dropped when no two-model decomposition runs", {
  set.seed(404)
  d <- data.frame(
    country = sample(c("g1", "g2"), 400, replace = TRUE),
    gender  = sample(c("F", "M"), 400, replace = TRUE),
    ses     = sample(c("low", "mid", "high"), 400, replace = TRUE),
    stringsAsFactors = FALSE)
  d$y <- rnorm(400)
  # crossed-dimensions mode does not compute the two-model PCV, so neither pcv nor
  # pcv_status should appear.
  res <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ (1 | gender:ses), data = d, group = "country",
                          decomposition = "crossed-dimensions")))
  expect_false("pcv" %in% names(res))
  expect_false("pcv_status" %in% names(res))
})

# ---- #3 Med-low: maihda_ic() guards against mixed IC scales ------------------

make_audit_ic_data <- function(seed = 4242, n = 300) {
  set.seed(seed)
  d <- data.frame(
    g1 = sample(c("F", "M"), n, replace = TRUE),
    g2 = sample(c("A", "B", "C"), n, replace = TRUE),
    x  = rnorm(n),
    stringsAsFactors = FALSE)
  stratum <- interaction(d$g1, d$g2, drop = TRUE)
  d$y <- 1 + 0.3 * d$x + rnorm(nlevels(stratum), sd = 0.5)[stratum] + rnorm(n, sd = 0.5)
  d
}

test_that("maihda_ic() warns and omits the delta for a mixed-scale (lme4 vs brms) table", {
  d <- make_audit_ic_data()
  # Same outcome, family and sample -> the outcome/family/sample block stays silent,
  # so the ONLY thing left to flag is the mixed likelihood/Bayesian IC scale.
  m1 <- fit_maihda(y ~ x + (1 | g1:g2), data = d)
  m2 <- fit_maihda(y ~ 1 + (1 | g1:g2), data = d)

  # Stand in for a brms fit on the second model without a Stan toolchain: its IC row
  # carries WAIC/LOOIC instead of AIC/BIC. maihda_ic() reads the rows in order.
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
    ic <- maihda_ic(m1, m2, model_names = c("lme4", "brms")),
    "likelihood/Bayesian divide|different scales|Bayesian WAIC"
  )
  # The delta is withheld (it would rank AIC vs a bare NA), but the per-model
  # criteria are still reported on both scales.
  expect_false("delta" %in% names(ic))
  expect_true(all(c("AIC", "BIC", "WAIC", "LOOIC") %in% names(ic)))
})

test_that("an all-lme4 maihda_ic() table keeps its delta and warns about no scale mix", {
  d <- make_audit_ic_data()
  m1 <- fit_maihda(y ~ x + (1 | g1:g2), data = d)
  m2 <- fit_maihda(y ~ 1 + (1 | g1:g2), data = d)
  expect_no_warning(ic <- suppressMessages(maihda_ic(m1, m2)))
  expect_true("delta" %in% names(ic))
  expect_false(any(c("WAIC", "LOOIC") %in% names(ic)))
})
