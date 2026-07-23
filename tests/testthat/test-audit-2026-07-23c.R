# Regression tests for the 2026-07-23 statistical audit, third file of the day
# (test-audit-2026-07-23.R and -07-23b.R belong to earlier passes).
#
# Finding (CONFIRMED): predict_maihda(type = "strata") silently dropped random
# slopes. maihda_stratum_ranef_lme4() / maihda_stratum_ranef_brms() selected
# only the (Intercept) column of the stratum random-effect block, so a
# NON-longitudinal fit with a hand-written slope -- fit_maihda(y ~ x +
# (1 + x | stratum)) -- returned a complete-looking table (estimate, SE, 95%
# interval) that was in fact the stratum effect at zero of the slope variable:
# an extrapolation for uncentered covariates (repro: returned sd 2.77 vs 0.76
# for the effects at mean(x); Spearman rank corr 0.33, signs flipped), while
# summary() on the SAME fit already errored via the intercept-only guard.
# The brms extractor additionally fell back to the lone effect column when no
# intercept was present, presenting a slope-only (0 + x | stratum) model's
# SLOPE as the stratum effect. Both extractors now reject non-intercept effect
# columns in the requested group with a directed error; the LONGITUDINAL
# summary -- whose growth block carries time slopes by design and reports them
# via its trajectory table -- opts out explicitly (allow_slope_columns = TRUE),
# and longitudinal type = "strata" predictions keep routing to the trajectory
# path before the extractor. The guard is scoped to the requested group: a
# slope on some OTHER random effect leaves the stratum's intercept BLUP a
# complete description of the stratum effect. wemix and ordinal engines reject
# slope formulas at fit time and were never exposed.

maihda_audit_0723c_data <- function() {
  set.seed(2318)
  K <- 6; n_per <- 40
  stratum <- factor(rep(sprintf("s%d", 1:K), each = n_per))
  x <- rnorm(K * n_per, mean = 10, sd = 3)
  u1 <- rnorm(K, 0, 0.3)
  u0 <- rnorm(K, 0, 0.5) - 10 * u1
  k <- as.integer(stratum)
  data.frame(
    y = 2 + 0.1 * x + u0[k] + u1[k] * x + rnorm(K * n_per, 0, 1),
    x = x,
    stratum = stratum,
    site = factor(rep_len(sprintf("c%d", 1:8), K * n_per))
  )
}

test_that("audit 2026-07-23c: stratum ranef slope stopper (pure)", {
  # No offending columns: silent pass-through.
  expect_invisible(maihda_stop_stratum_ranef_slopes("stratum", character(0)))
  expect_true(maihda_stop_stratum_ranef_slopes("stratum", character(0)))

  # Offending columns: directed error naming the group, the columns, and the
  # longitudinal route.
  expect_error(maihda_stop_stratum_ranef_slopes("stratum", c("x", "z")),
               "random slopes \\(x, z\\)")
  expect_error(maihda_stop_stratum_ranef_slopes("stratum", "x"),
               "fit_maihda\\(id = , time = \\)")
})

test_that("audit 2026-07-23c: lme4 stratum extractor rejects slope columns", {
  d <- maihda_audit_0723c_data()

  # Compound intercept + slope: the old extractor silently returned the
  # intercept column only.
  m_slope <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (1 + x | stratum), data = d)))
  expect_error(maihda_stratum_ranef_lme4(m_slope), "random slopes")

  # Slope-only block: no intercept to fall back on; same directed error.
  m_slope_only <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (0 + x | stratum), data = d)))
  expect_error(maihda_stratum_ranef_lme4(m_slope_only), "random slopes")

  # The longitudinal summary's explicit opt-out still reads the intercept
  # column (the trajectory table carries the slopes separately).
  re <- lme4::ranef(m_slope, condVar = TRUE)$stratum
  out <- maihda_stratum_ranef_lme4(m_slope, allow_slope_columns = TRUE)
  expect_identical(out$stratum, rownames(re))
  expect_equal(out$random_effect, re[["(Intercept)"]])
  expect_true(all(is.finite(out$se)) && all(out$se > 0))

  # Scoped to the requested group: a slope on ANOTHER group leaves the
  # intercept-only stratum block extractable, while that other group errors.
  m_other <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (1 | stratum) + (1 + x | site), data = d)))
  out_s <- maihda_stratum_ranef_lme4(m_other, group = "stratum")
  expect_equal(nrow(out_s), nlevels(d$stratum))
  expect_error(maihda_stratum_ranef_lme4(m_other, group = "site"),
               "random slopes")

  # Intercept-only fits are untouched.
  m_ok <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (1 | stratum), data = d)))
  re_ok <- lme4::ranef(m_ok, condVar = TRUE)$stratum
  out_ok <- maihda_stratum_ranef_lme4(m_ok)
  expect_equal(out_ok$random_effect, re_ok[["(Intercept)"]])
})

test_that("audit 2026-07-23c: predict_maihda(type = 'strata') errors on a non-longitudinal slope fit", {
  d <- maihda_audit_0723c_data()

  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 + x | stratum), data = d)))
  expect_null(fit$longitudinal_info)

  # The reported defect: this returned an intercept-only table with no signal.
  expect_error(predict_maihda(fit, type = "strata"), "random slopes")

  # Individual predictions on the same fit are unaffected (lme4's predict
  # handles the slopes itself).
  p_ind <- predict_maihda(fit, type = "individual")
  expect_length(p_ind, nrow(d))
  expect_true(all(is.finite(p_ind)))

  # Intercept-only fits keep the documented table shape.
  fit0 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d)))
  p0 <- predict_maihda(fit0, type = "strata")
  expect_named(p0, c("stratum", "predicted", "se", "lower_95", "upper_95"))
  expect_equal(nrow(p0), nlevels(d$stratum))
})

test_that("audit 2026-07-23c: longitudinal trajectory predictions and summaries are untouched", {
  set.seed(915)
  P <- 40; W <- 3; G <- 4
  pid <- rep(sprintf("p%02d", 1:P), each = W)
  tt <- rep(0:(W - 1), P)
  strat_p <- rep(sprintf("g%d", 1 + (0:(P - 1)) %% G), each = W)
  v0 <- rnorm(G, 0, 0.4); v1 <- rnorm(G, 0, 0.2)
  pi0 <- rnorm(P, 0, 0.3); pi1 <- rnorm(P, 0, 0.1)
  gi <- as.integer(factor(strat_p)); pii <- as.integer(factor(pid))
  dl <- data.frame(
    y = 1 + 0.2 * tt + v0[gi] + v1[gi] * tt + pi0[pii] + pi1[pii] * tt +
      rnorm(P * W, 0, 0.5),
    t = tt, pid = pid, stratum = factor(strat_p)
  )

  fitl <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = dl, id = "pid", time = "t")))

  # The growth block carries a time slope by design; type = "strata" routes to
  # the trajectory table, not the guarded extractor.
  pl <- predict_maihda(fitl, type = "strata")
  expect_true(all(c("stratum", "baseline", "intercept", "slope") %in% names(pl)))
  expect_equal(nrow(pl), G)

  # The longitudinal summary reads the intercept column through the explicit
  # opt-out; the guard must not break it.
  sl <- suppressWarnings(summary(fitl))
  expect_s3_class(sl, "maihda_summary")
  expect_true(is.data.frame(sl$stratum_estimates) &&
                nrow(sl$stratum_estimates) == G)
})

test_that("audit 2026-07-23c: brms stratum extractor rejects slope columns", {
  skip_if_not_installed("brms")

  make_arr <- function(effects) {
    array(
      rep(c(0.2, -0.2, 0.05, 0.06, 0.1, -0.3, 0.3, -0.1),
          times = length(effects)),
      dim = c(2, 4, length(effects)),
      dimnames = list(c("1", "2"),
                      c("Estimate", "Est.Error", "Q2.5", "Q97.5"),
                      effects)
    )
  }
  fake <- structure(list(), class = "brmsfit")

  with_ranef <- function(arr, expr) {
    testthat::local_mocked_bindings(
      ranef = function(object, ...) list(stratum = arr),
      .package = "brms")
    force(expr)
  }

  # Intercept + slope: directed error instead of the silent intercept column.
  arr2 <- make_arr(c("Intercept", "x"))
  expect_error(with_ranef(arr2, maihda_stratum_ranef_brms(fake)),
               "random slopes")

  # Slope-only single column: the old lone-effect fallback presented the SLOPE
  # as the stratum effect; now the same directed error.
  arr_slope <- make_arr("x")
  expect_error(with_ranef(arr_slope, maihda_stratum_ranef_brms(fake)),
               "random slopes")

  # The longitudinal opt-out still reads the intercept column.
  out <- with_ranef(arr2,
                    maihda_stratum_ranef_brms(fake, allow_slope_columns = TRUE))
  expect_equal(out$random_effect, unname(arr2[, "Estimate", "Intercept"]))

  # Intercept-only blocks are untouched.
  arr1 <- make_arr("Intercept")
  out1 <- with_ranef(arr1, maihda_stratum_ranef_brms(fake))
  expect_equal(out1$random_effect, unname(arr1[, "Estimate", "Intercept"]))
  expect_equal(out1$se, unname(arr1[, "Est.Error", "Intercept"]))
})
