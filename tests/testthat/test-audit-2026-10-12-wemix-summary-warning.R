# Audit 2026-10-12 (packet 2026-10-12-wemix-summary-warning): summary() on a
# design-weighted logistic fit (engine = "wemix", family = "binomial") passed on
# "diag(V) had non-positive or NA entries; the non-finite result may be dubious"
# whenever the between-stratum variance was small -- CONFIRMED. Whatever calls summary()
# passed it on too, once per summarised fit: tidy(), glance(), maihda_table(), the
# plots, maihda() (two from a two-model analysis) and compare_maihda_groups() (two from
# two groups) among them.
#
# The warning is WeMix's own (WeMix 4.0.3). summary.WeMixResults() rebuilds the sandwich
# covariance of the estimates, where makeSandwich() sets the standard error of a
# variance to NA when log(variance) < -3.6, and then calls stats::cov2cor() on the
# block belonging to the variances to fill its correlation columns. With one random
# intercept that block is 1 x 1, there is no correlation, and the result is thrown
# away: the warning reports a missing diagonal entry of a matrix nothing reads.
# Measured on both sides of exp(-3.6) = 0.02732, over 23 binomial fits with variances
# from 0.003 to 0.84 (unit and non-unit weights): the 12 below it warned, the 11 above
# did not. None of 14 Gaussian fits warned, zero-variance ones included: on that path
# cov2cor() is given the random-effect variance itself, which comes back as about 1e-22
# and not as 0. Under options(warn = 2) the warning became an error inside the summary's
# own tryCatch and the model_summary slot came back NULL.
#
# The package reads nothing from WeMix's summary object; it hands it on as the
# model_summary slot. maihda_wemix_summary() drops that one warning, recognised by the
# call that raised it (R translates the text), and passes every other. The object
# returned and every reported number are as they were.
#
# Warnings are COUNTED here, never suppressed around the call under test. The warning
# had been met before: test-audit-2026-09-13d.R changed its fixture to get away from it.

# Every warning an expression raises: how many, and the function each came from.
ws12_warnings <- function(expr) {
  from <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    cl <- conditionCall(w)
    from <<- c(from, if (is.call(cl)) paste(deparse(cl[[1L]]), collapse = "") else "")
    invokeRestart("muffleWarning")
  })
  list(value = value, n = length(from), from = from)
}

# Eight strata of `n` rows with the given numbers of successes, unit weights. No random
# numbers.
ws12_data <- function(successes, n) {
  g <- expand.grid(d1 = 0:1, d2 = 0:1, d3 = 0:1)
  g$stratum <- seq_len(nrow(g))
  d <- g[rep(1:8, each = n), ]
  d$y <- unlist(lapply(1:8, function(j) rep(c(1L, 0L), c(successes[j], n - successes[j]))))
  d$w <- 1
  d
}

# The reported data. WeMix puts the between-stratum variance at 0.0222, below exp(-3.6).
ws12_reported <- function() ws12_data(c(12, 16, 16, 20, 16, 20, 20, 24), 40)

ws12_null <- y ~ 1 + (1 | stratum)

ws12_fit <- function(formula, data, family) {
  suppressMessages(fit_maihda(formula, data = data, engine = "wemix", family = family,
                              sampling_weights = "w"))
}

# A stand-in for a fit whose summary() raises R's own cov2cor() warning twice, through a
# bare and a namespace-qualified call, and then an unrelated warning.
ws12_mock_fit <- function(fail = FALSE) {
  registerS3method("summary", "ws12_mock_fit", function(object, ...) {
    if (isTRUE(object$fail)) stop("the summary failed")
    cov2cor(matrix(NA_real_, 1, 1))
    stats::cov2cor(matrix(0, 1, 1))
    warning("an unrelated warning")
    "the summary"
  }, envir = asNamespace("base"))
  structure(list(fail = fail), class = "ws12_mock_fit")
}

test_that("the cov2cor() warning is told by the call that raised it, not by its text", {
  words <- "diag(V) had non-positive or NA entries; the non-finite result may be dubious"
  expect_true(maihda_wemix_cov2cor_warning(simpleWarning(words, quote(cov2cor(vc)))))
  # The text is not read: R translates it.
  expect_true(maihda_wemix_cov2cor_warning(
    simpleWarning("the same warning in another language", quote(cov2cor(vc)))))
  # A namespace-qualified call names the same function.
  expect_true(maihda_wemix_cov2cor_warning(
    simpleWarning(words, quote(stats::cov2cor(vc)))))
  expect_true(maihda_wemix_cov2cor_warning(
    simpleWarning(words, quote(stats:::cov2cor(vc)))))

  # The same words from anywhere else are not it, nor is a warning without a call.
  expect_false(maihda_wemix_cov2cor_warning(simpleWarning(words, quote(sqrt(x)))))
  expect_false(maihda_wemix_cov2cor_warning(simpleWarning(words)))
  expect_false(maihda_wemix_cov2cor_warning(simpleWarning(words, quote(cov2cor2(vc)))))
  expect_false(maihda_wemix_cov2cor_warning(
    simpleWarning(words, quote((function(V) V)(vc)))))
  # Nor is a function of that name reached through an object.
  expect_false(maihda_wemix_cov2cor_warning(simpleWarning(words, quote(obj$cov2cor(vc)))))

  # R's own conditions, as this R raises them.
  missing_diag <- tryCatch(stats::cov2cor(matrix(NA_real_, 1, 1)), warning = function(w) w)
  expect_s3_class(missing_diag, "warning")
  expect_true(maihda_wemix_cov2cor_warning(missing_diag))
  zero_diag <- tryCatch(stats::cov2cor(matrix(0, 1, 1)), warning = function(w) w)
  expect_s3_class(zero_diag, "warning")
  expect_true(maihda_wemix_cov2cor_warning(zero_diag))
  expect_false(maihda_wemix_cov2cor_warning(tryCatch(sqrt(-1), warning = function(w) w)))
  expect_false(maihda_wemix_cov2cor_warning(
    tryCatch(as.integer("a"), warning = function(w) w)))
})

test_that("maihda_wemix_summary() drops the cov2cor() warning and passes every other", {
  fit <- ws12_mock_fit()
  # With nothing in the way the stand-in warns three times.
  bare <- ws12_warnings(summary(fit))
  expect_identical(bare$from, c("cov2cor", "stats::cov2cor", "summary.ws12_mock_fit"))

  got <- ws12_warnings(maihda_wemix_summary(fit))
  expect_identical(got$from, "summary.ws12_mock_fit")
  expect_identical(got$value, "the summary")

  # An error leaves as it was raised; summary() makes a NULL slot of it.
  expect_error(maihda_wemix_summary(ws12_mock_fit(fail = TRUE)), "the summary failed")
})

test_that("summary() of a binomial wemix fit with a small variance raises no warning", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  fit <- ws12_warnings(ws12_fit(ws12_null, ws12_reported(), "binomial"))
  m <- fit$value
  expect_identical(fit$n, 0L)
  v <- as.numeric(m$model$vars)[1]
  expect_equal(v, 0.02221119, tolerance = 1e-3)
  expect_lt(log(v), -3.6)

  # WeMix's own summary of this fit, with no package code in between: one warning, from
  # cov2cor(), and no standard error for the variance. If these two expectations start
  # to fail, WeMix no longer warns here and maihda_wemix_summary() has nothing left to
  # drop.
  own <- ws12_warnings(summary(m$model))
  expect_identical(own$from, "cov2cor")
  expect_true(is.na(own$value$varsmat[["Std. Error"]]))

  s <- ws12_warnings(summary(m))
  expect_identical(s$n, 0L)
  s <- s$value

  # The slot still holds WeMix's summary, as WeMix built it.
  expect_s3_class(s$model_summary, "summaryWeMixResults")
  expect_identical(s$model_summary, own$value)

  # The reported numbers.
  vc <- s$variance_components
  expect_identical(vc$variance[vc$component == "Between-stratum (random)"], v)
  expect_equal(vc$variance[vc$component == "Within-stratum (residual)"], pi^2 / 3)
  expect_equal(s$vpc$estimate, 0.006706116, tolerance = 1e-3)
  expect_equal(s$fixed_effects$estimate, -0.2017761, tolerance = 1e-3)
  expect_equal(s$fixed_effects$se, 0.1333183, tolerance = 1e-3)
  expect_equal(s$discriminatory_accuracy$auc, 0.5946970, tolerance = 1e-6)
  expect_equal(s$discriminatory_accuracy$mor, 1.152761, tolerance = 1e-3)
  expect_equal(nrow(s$stratum_estimates), 8L)

  # Everything that calls summary() on the fit.
  expect_identical(ws12_warnings(utils::capture.output(print(s)))$n, 0L)
  expect_identical(ws12_warnings(generics::tidy(m))$n, 0L)
  expect_identical(ws12_warnings(generics::glance(m))$n, 0L)
  expect_identical(ws12_warnings(maihda_table(m))$n, 0L)
  expect_identical(ws12_warnings(plot(m, type = "vpc"))$n, 0L)

  # With warnings turned into errors the slot is not lost.
  strict <- local({
    op <- options(warn = 2)
    on.exit(options(op))
    summary(m)
  })
  expect_identical(strict$model_summary, own$value)
})

test_that("a wemix summary that WeMix built without a warning is as it was", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  # A binomial fit whose variance is above exp(-3.6): WeMix reports a standard error
  # for it and its summary is silent.
  m <- ws12_fit(ws12_null, ws12_data(c(41, 56, 35, 49, 47, 38, 52, 43), 100), "binomial")
  v <- as.numeric(m$model$vars)[1]
  expect_equal(v, 0.03353942, tolerance = 1e-3)
  expect_gt(log(v), -3.6)
  own <- ws12_warnings(summary(m$model))
  expect_identical(own$n, 0L)
  expect_true(is.finite(own$value$varsmat[["Std. Error"]]))
  s <- ws12_warnings(summary(m))
  expect_identical(s$n, 0L)
  expect_identical(s$value$model_summary, own$value)
  expect_identical(
    s$value$variance_components$variance[
      s$value$variance_components$component == "Between-stratum (random)"], v)

  # A Gaussian fit with every stratum mean equal: the variance is zero to rounding and
  # the summary is silent.
  g <- ws12_data(rep(0, 8), 40)
  g$y <- rep(stats::qnorm(stats::ppoints(40)), times = 8)
  mg <- ws12_fit(ws12_null, g, "gaussian")
  expect_lt(as.numeric(mg$model$vars)[1], 1e-10)
  sg <- ws12_warnings(summary(mg))
  expect_identical(sg$n, 0L)
  expect_s3_class(sg$value$model_summary, "summaryWeMixResults")
  expect_lt(sg$value$vpc$estimate, 1e-10)
})

test_that("maihda() and compare_maihda_groups() pass no warning on from wemix summaries", {
  skip_on_cran()
  skip_slow_test()
  skip_if_not_installed("WeMix")

  # Two groups holding the reported data: a null model per group, both variances below
  # exp(-3.6).
  d <- ws12_reported()
  two <- rbind(transform(d, grp = "a"), transform(d, grp = "b"))
  cmp <- ws12_warnings(suppressMessages(compare_maihda_groups(
    ws12_null, data = two, group = "grp", engine = "wemix", family = "binomial",
    sampling_weights = "w")))
  expect_identical(cmp$n, 0L)
  expect_identical(cmp$value$status, c("ok", "ok"))
  expect_equal(cmp$value$var_between, rep(0.02221119, 2), tolerance = 1e-3)
  expect_equal(cmp$value$vpc, rep(0.006706116, 2), tolerance = 1e-3)

  # A two-model analysis whose null and adjusted variances are both below it.
  d3 <- ws12_data(c(31, 25, 25, 35, 25, 35, 35, 29), 60)
  d3$stratum <- NULL
  an <- ws12_warnings(suppressMessages(maihda(
    y ~ 1 + (1 | d1:d2:d3), data = d3, engine = "wemix", family = "binomial",
    sampling_weights = "w")))
  expect_identical(an$n, 0L)
  an <- an$value
  v0 <- as.numeric(an$model$model$vars)[1]
  v1 <- as.numeric(an$model_adjusted$model$vars)[1]
  expect_equal(v0, 0.0182236, tolerance = 1e-3)
  expect_equal(v1, 0.004707988, tolerance = 1e-3)
  expect_lt(log(v0), -3.6)
  expect_lt(log(v1), -3.6)
  expect_equal(an$pcv$pcv, 0.741654, tolerance = 1e-3)
  expect_s3_class(an$summary$model_summary, "summaryWeMixResults")
  expect_s3_class(an$summary_adjusted$model_summary, "summaryWeMixResults")
})
