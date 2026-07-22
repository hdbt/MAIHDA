# Regression tests for the second 2026-07-22 statistical audit (thirty-ninth
# pass; test-audit-2026-07-22.R belongs to the earlier pass on f993586).
#
# Finding #1 (CONFIRMED): longitudinal validation permitted unidentified growth
# models -- one duplicated id was all it required. A quadratic on all-equal
# times fit "successfully" and reported the optimizer's STARTING VALUES as
# slope variances (the slope columns of Z are all zero, so the likelihood is
# flat in them); a quadratic on two waves gave a ~1300-fold optimizer-dependent
# swing in the between-stratum variance at unobserved times (VarS(0.5) 0.0074
# default vs 9.70 Nelder-Mead on identical data); replicate-only data (ids
# repeating at a single time) reported a fabricated person slope variance. All
# with at most a "boundary (singular) fit" note. maihda_check_longitudinal_times()
# now stops on (i) fewer than time_degree + 1 distinct observed times and
# (ii) no person measured at two distinct times, on the full input
# (maihda_validate_longitudinal) and again on the analytic sample (fit_maihda).
#
# Finding #2 (CONFIRMED, code half): sampling-weighted brms fits are
# pseudo-posteriors (weighted likelihood -- fit_maihda says so at fit time),
# but maihda_ic() reported their WAIC/LOOIC as plain Bayesian criteria while
# the wemix pseudo-likelihood correctly got NA AIC/BIC. maihda_ic_one() now
# reports NA with estimator "Bayesian (weighted pseudo-posterior)"; the
# unweighted brms path is unchanged. The conditional predictive target of
# WAIC/LOOIC (new observations within represented strata, not new strata) is
# documented in the maihda_ic Details.

test_that("audit 2026-07-22b #1: growth polynomial needs time_degree + 1 distinct times", {
  # The finding's exact scenario: quadratic growth, every time value equal.
  d0 <- data.frame(pid = rep(1:8, each = 3), wave = 0,
                   y = rnorm(24))
  expect_error(
    maihda_check_longitudinal_times(d0, "pid", "wave", 2L),
    "1 distinct observed value.*at least 3")
  # Even linear growth is unidentified on a single observed time.
  expect_error(
    maihda_check_longitudinal_times(d0, "pid", "wave", 1L),
    "no time variation")
  # Two waves cannot identify a quadratic: I(t^2) is collinear with t on {0, 1}.
  d2 <- data.frame(pid = rep(1:8, each = 2), wave = rep(0:1, 8), y = rnorm(16))
  expect_error(
    maihda_check_longitudinal_times(d2, "pid", "wave", 2L),
    "2 distinct observed value")
  # ...but exactly time_degree + 1 distinct times passes (boundary case).
  expect_silent(maihda_check_longitudinal_times(d2, "pid", "wave", 1L))
  d3 <- data.frame(pid = rep(1:8, each = 3), wave = rep(0:2, 8), y = rnorm(24))
  expect_silent(maihda_check_longitudinal_times(d3, "pid", "wave", 2L))
})

test_that("audit 2026-07-22b #1: replicate-only ids (no within-person time variation) are rejected", {
  # Ids repeat -- the old duplicated-id gate passes -- but every person sits at
  # ONE time. Between-person spread keeps 3 distinct times, so the global
  # distinct-times check passes and the within-person check must catch it.
  d <- data.frame(pid = rep(1:9, each = 3),
                  wave = rep(c(0, 1, 2), each = 3, length.out = 27),
                  y = rnorm(27))
  stopifnot(length(unique(d$wave)) == 3)
  expect_error(
    maihda_check_longitudinal_times(d, "pid", "wave", 1L),
    "within-person time variation.*replicate measurements")
  # One person with two distinct times clears the gate.
  d_ok <- d
  d_ok$wave[3] <- 1   # pid 1: waves 0, 0, 1
  expect_silent(maihda_check_longitudinal_times(d_ok, "pid", "wave", 1L))
  # NA ids and non-finite times are ignored, never counted as variation.
  d_na <- d_ok
  d_na$wave[3] <- NA   # pid 1's second distinct time was an NA all along
  d_na$pid[4] <- NA    # and an NA id cannot clear the within-person gate
  expect_error(maihda_check_longitudinal_times(d_na, "pid", "wave", 1L),
               "within-person")
})

test_that("audit 2026-07-22b #1: fit_maihda rejects unidentified growth models end-to-end", {
  set.seed(1)
  # Validation path: the finding's exact all-equal-times quadratic.
  d <- data.frame(pid = rep(1:12, each = 3), wave = 0,
                  g1 = rep(c("a", "b"), each = 18),
                  g2 = rep(rep(c("x", "y"), 6), each = 3),  # constant within person
                  y = rnorm(36))
  expect_error(
    fit_maihda(y ~ (1 | g1:g2), data = d, id = "pid", time = "wave",
               time_degree = 2),
    "unidentified")

  # Analytic-sample re-check: the INPUT passes validation (waves 0 and 1, two
  # distinct times per person), but every wave-1 outcome is missing, so the
  # fitted rows collapse to a single time. Persons keep two wave-0 rows each,
  # so the older duplicated-id re-check still passes; only the distinct-times
  # re-check can catch it.
  d2 <- data.frame(pid = rep(1:12, each = 3),
                   wave = rep(c(0, 0, 1), 12),
                   g1 = rep(c("a", "b"), each = 18),
                   g2 = rep(rep(c("x", "y"), 6), each = 3),  # constant within person
                   y = rnorm(36))
  d2$y[d2$wave == 1] <- NA
  expect_silent(maihda_validate_longitudinal("pid", "wave", 1, d2))
  expect_error(
    fit_maihda(y ~ (1 | g1:g2), data = d2, id = "pid", time = "wave"),
    "unidentified")

  # Control: a healthy 3-occasion linear growth fit still goes through.
  d3 <- d2
  d3$y <- rnorm(36)   # no missingness; waves 0, 0, 1 per person
  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | g1:g2), data = d3, id = "pid", time = "wave")))
  expect_s3_class(fit, "maihda_model")
})

test_that("audit 2026-07-22b #2: weighted-brms pseudo-posterior gets NA WAIC/LOOIC", {
  # A sampling-weighted brms fit is a pseudo-posterior; its weighted pointwise
  # log-likelihoods are not log predictive densities, so maihda_ic() must treat
  # it like the wemix pseudo-likelihood (NA criteria), not as plain "Bayesian".
  # The branch keys on the wrapper's $sampling_weights, so a mock fitted object
  # exercises it without a Stan fit. The mock carries $data because a loaded
  # brms namespace dispatches nobs.brmsfit on it (real fits always have $data).
  fake_fit <- structure(list(data = data.frame(y = 1:4, x = 1:4)),
                        class = "brmsfit")
  weighted <- structure(
    list(model = fake_fit, engine = "brms",
         formula = y ~ x + (1 | stratum),
         data = data.frame(y = 1:4, x = 1:4, w = 1,
                           stratum = c("a", "a", "b", "b")),
         sampling_weights = "w"),
    class = "maihda_model")
  row <- maihda_ic_one(weighted)
  expect_identical(row$estimator, "Bayesian (weighted pseudo-posterior)")
  expect_true(is.na(row$WAIC))
  expect_true(is.na(row$LOOIC))

  # And the user-facing table drops the all-NA criterion columns cleanly.
  tab <- maihda_ic(weighted)
  expect_s3_class(tab, "maihda_ic")
  expect_false(any(c("AIC", "BIC", "WAIC", "LOOIC") %in% names(tab)))
  expect_identical(tab$estimator, "Bayesian (weighted pseudo-posterior)")
})

test_that("audit 2026-07-22b #2: unweighted brms path is unchanged (confinement)", {
  skip_if_not_installed("brms")
  # Mock the IC evaluator so no real draws are needed: if the unweighted branch
  # runs, the mocked estimates flow into the row; if the WEIGHTED branch were
  # (wrongly) still computing criteria, the same mock would fill them too.
  est <- matrix(c(111.1, 1.0, 222.2, 2.0), nrow = 2, byrow = TRUE,
                dimnames = list(c("waic", "looic"), c("Estimate", "SE")))
  testthat::local_mocked_bindings(
    maihda_ic_quiet_but_warn = function(expr, what) list(estimates = est)
  )
  fake_fit <- structure(list(data = data.frame(y = 1:4, x = 1:4)),
                        class = "brmsfit")
  unweighted <- structure(
    list(model = fake_fit, engine = "brms",
         formula = y ~ x + (1 | stratum),
         data = data.frame(y = 1:4, x = 1:4)),
    class = "maihda_model")
  row <- maihda_ic_one(unweighted)
  expect_identical(row$estimator, "Bayesian")
  expect_equal(row$WAIC, 111.1)
  expect_equal(row$LOOIC, 222.2)

  # The weighted branch must short-circuit BEFORE any criterion evaluation:
  # with the same mock in place, the criteria stay NA.
  weighted <- unweighted
  weighted$sampling_weights <- "w"
  row_w <- maihda_ic_one(weighted)
  expect_identical(row_w$estimator, "Bayesian (weighted pseudo-posterior)")
  expect_true(is.na(row_w$WAIC) && is.na(row_w$LOOIC))
})
