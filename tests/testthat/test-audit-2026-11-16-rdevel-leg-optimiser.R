# Audit 2026-11-16, one finding from CI: the maximum-likelihood refit behind the AIC and BIC of
# REML lmer fits passed on a false alarm of lme4's optimiser.
#
# maihda_ic() with more than one model, and compare_maihda(), refit a REML lmer fit by
# lme4::refitML(). Its optimiser, minqa's bobyqa, can end at the maximum with an exit code ("a
# trust region step failed to reduce q") where a variance is zero or small: in 2 to 6% of fits
# with a zero or small between-stratum variance, and each of 56 such refits was within 5e-13 of
# the maximum-likelihood deviance. lme4's warning for the code reached the user. Which fits it
# strikes depends on how lme4 and minqa were compiled, so a block of test-compare_maihda.R that
# counts the warnings of compare_maihda() failed on GitHub's ubuntu-26.04 runners, whose gcc
# fuses a multiplication and an addition into one operation, with the seed it passes on
# elsewhere; on Windows another seed of the same design gives the same two warnings.
#
# FIX: maihda_ic_one() makes the refit through maihda_ic_refit_ml(). A refit whose optimiser
# reports a code is made once more by lmer()'s default optimiser, which replaces it when its
# deviance is lower by more than 1e-4, and lme4's warning is passed on only when the second
# refit reports a code as well or fails. A refit without a code is used as lme4 makes it.
#
# Random numbers are drawn only inside ic16_keep_seed(), which puts the suite's stream back.

ic16_keep_seed <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  force(expr)
}

ic16_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression with the text of every warning and message it raised.
ic16_heard <- function(expr) {
  warnings <- character(0)
  messages <- character(0)
  value <- withCallingHandlers(expr,
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    })
  list(value = value, warnings = warnings, messages = messages)
}

# Eight strata of 20 with a between-stratum SD of 0.6: the REML and the maximum-likelihood
# estimates differ (theta 0.564 against 0.523), so a refit stopped after 5 evaluations is short
# of the maximum, by 0.004.
ic16_strata_data <- function() {
  ic16_keep_seed({
    set.seed(1)
    st <- rep(seq_len(8), each = 20)
    d <- data.frame(stratum = factor(st), x = stats::rnorm(160))
    d$y <- 1 + 0.3 * d$x + stats::rnorm(8, 0, 0.6)[st] + stats::rnorm(160)
    d
  })
}

# The design of the block of test-compare_maihda.R that failed: 200 rows, an outcome of pure
# noise and strata of two dimensions, so that the between-stratum variance is at or near zero.
ic16_noise_data <- function(seed) {
  ic16_keep_seed({
    set.seed(seed)
    n <- 200
    d <- data.frame(a = sample(c("p", "q"), n, replace = TRUE),
                    b = sample(c("x", "y"), n, replace = TRUE),
                    cc = sample(c("m", "n"), n, replace = TRUE),
                    age = stats::rnorm(n))
    d$y <- stats::rnorm(n)
    d
  })
}

# `refit` as lme4 returns a refit whose optimiser reported `code`: the code in its record, and
# lme4's warning for it raised.
ic16_flagged <- function(refit, optimiser, code = 3L) {
  refit@optinfo$conv$opt <- code
  warning("convergence code ", code, " from ", optimiser, ": ", optimiser,
          " -- a trust region step failed to reduce q")
  refit
}

code_warning <- "convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q"

# ---- the check ------------------------------------------------------------------------------

test_that("a refit flagged at its maximum is checked by a second refit and returned without lme4's warning", {
  skip_on_cran()
  d <- ic16_strata_data()
  fm <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))$model
  real <- lme4::refitML
  own <- real(fm)
  calls <- character(0)
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    calls <<- c(calls, optimizer)
    refit <- real(x, optimizer = optimizer, ...)
    if (identical(optimizer, "bobyqa")) ic16_flagged(refit, optimizer) else refit
  }, .package = "lme4")
  h <- ic16_heard(maihda_ic_refit_ml(fm))
  # lme4's refit, then lmer's default optimiser from the same estimates.
  expect_identical(calls, c("bobyqa", "nloptwrap"))
  expect_identical(h$warnings, character(0))
  # The refit returned is lme4's own, to the last digit: no number changes.
  expect_identical(h$value@optinfo$optimizer, "bobyqa")
  expect_identical(stats::deviance(h$value), stats::deviance(own))
  expect_identical(lme4::getME(h$value, "theta"), lme4::getME(own, "theta"))
})

test_that("the first refit stays unless the second is lower by more than the tolerance of 1e-4", {
  skip_on_cran()
  d <- ic16_strata_data()
  fm <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))$model
  real <- lme4::refitML
  by_nloptwrap <- stats::deviance(ic16_quiet(real(fm, optimizer = "nloptwrap")))
  # lme4's refit stopped after `evaluations` evaluations, which lme4 flags itself, with a
  # message of its own.
  stopped_after <- function(evaluations) {
    testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
      if (identical(optimizer, "bobyqa")) {
        message("a message of the first refit")
        real(x, control = list(maxfun = evaluations), ...)
      } else {
        real(x, optimizer = optimizer, ...)
      }
    }, .package = "lme4")
    ic16_heard(maihda_ic_refit_ml(fm))
  }
  # The premise: after 7 evaluations it is 3.6e-5 above the second refit's deviance, within
  # the tolerance; after 6, 2.1e-4 above, beyond it.
  seven <- ic16_quiet(real(fm, control = list(maxfun = 7)))
  six <- ic16_quiet(real(fm, control = list(maxfun = 6)))
  expect_gt(stats::deviance(seven) - by_nloptwrap, 1e-5)
  expect_lt(stats::deviance(seven) - by_nloptwrap, 9e-5)
  expect_gt(stats::deviance(six) - by_nloptwrap, 1.5e-4)
  expect_lt(stats::deviance(six) - by_nloptwrap, 3e-4)
  # Within the tolerance the first refit is returned.
  h <- stopped_after(7)
  expect_identical(h$value@optinfo$optimizer, "bobyqa")
  expect_identical(stats::deviance(h$value), stats::deviance(seven))
  # Its exit code is not passed on; what minqa said about the limit is, and its message.
  expect_length(h$warnings, 1L)
  expect_match(h$warnings, "maxfun", fixed = TRUE)
  expect_false(maihda_is_optimiser_code_warning(simpleWarning(h$warnings)))
  expect_identical(h$messages, "a message of the first refit\n")
  # Beyond it the second refit is returned, and nothing of the first is signalled.
  h <- stopped_after(6)
  expect_identical(h$value@optinfo$optimizer, "nloptwrap")
  expect_identical(stats::deviance(h$value), by_nloptwrap)
  expect_identical(h$warnings, character(0))
  expect_identical(h$messages, character(0))
  # The tolerance is the one ?maihda_ic gives.
  expect_identical(formals(maihda_ic_refit_ml)$tol, 1e-4)
})

test_that("a flagged refit short of its maximum gives way to the second refit", {
  skip_on_cran()
  d <- ic16_strata_data()
  fm <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))$model
  direct <- ic16_quiet(lme4::lmer(y ~ x + (1 | stratum), data = d, REML = FALSE))
  real <- lme4::refitML
  # The premise: stopped after 5 evaluations, lme4's refit is flagged by lme4 itself, 0.004
  # above the maximum-likelihood deviance.
  short <- ic16_heard(real(fm, control = list(maxfun = 5)))
  expect_true(maihda_lme4_optimizer_failed(short$value))
  expect_true(any(vapply(short$warnings, function(w) {
    maihda_is_optimiser_code_warning(simpleWarning(w))
  }, logical(1))))
  expect_gt(stats::deviance(short$value) - stats::deviance(direct), 1e-3)
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    if (identical(optimizer, "bobyqa")) {
      real(x, control = list(maxfun = 5), ...)
    } else {
      real(x, optimizer = optimizer, ...)
    }
  }, .package = "lme4")
  h <- ic16_heard(maihda_ic_refit_ml(fm))
  # The second refit is the lower by more than the tolerance and ended without a code: it is
  # returned, and nothing of the first refit is signalled.
  expect_identical(h$value@optinfo$optimizer, "nloptwrap")
  expect_false(maihda_lme4_optimizer_failed(h$value))
  expect_identical(h$warnings, character(0))
  expect_lt(abs(stats::deviance(h$value) - stats::deviance(direct)), 1e-6)
  # Through maihda_ic_one(): the criteria are the maximum-likelihood fit's.
  model <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))
  row <- ic16_heard(maihda_ic_one(model, ml = TRUE))
  expect_identical(row$warnings, character(0))
  expect_identical(row$value$estimator, "ML (refit from REML)")
  expect_lt(abs(row$value$AIC - stats::AIC(direct)), 1e-6)
})

test_that("lme4's warning is passed on when the second refit is flagged too or fails", {
  skip_on_cran()
  d <- ic16_strata_data()
  fm <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))$model
  real <- lme4::refitML
  # Both refits flagged at the same deviance: the first is returned with its warning.
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    ic16_flagged(real(x, optimizer = optimizer, ...), optimizer)
  }, .package = "lme4")
  h <- ic16_heard(maihda_ic_refit_ml(fm))
  expect_identical(h$warnings, code_warning)
  expect_identical(h$value@optinfo$optimizer, "bobyqa")
  # The second refit fails: the first as lme4 returned it.
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    if (!identical(optimizer, "bobyqa")) stop("the second refit fails")
    ic16_flagged(real(x, optimizer = optimizer, ...), optimizer)
  }, .package = "lme4")
  h <- ic16_heard(maihda_ic_refit_ml(fm))
  expect_identical(h$warnings, code_warning)
  expect_identical(h$value@optinfo$optimizer, "bobyqa")
  # Both flagged and the second the lower: the second is returned, with its own warning.
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    if (identical(optimizer, "bobyqa")) {
      real(x, control = list(maxfun = 5), ...)
    } else {
      ic16_flagged(real(x, optimizer = optimizer, ...), optimizer, code = 5L)
    }
  }, .package = "lme4")
  h <- ic16_heard(maihda_ic_refit_ml(fm))
  expect_identical(h$value@optinfo$optimizer, "nloptwrap")
  expect_identical(
    h$warnings,
    "convergence code 5 from nloptwrap: nloptwrap -- a trust region step failed to reduce q")
  # An error of the first refit passes through, and maihda_ic_one() keeps the REML fit.
  testthat::local_mocked_bindings(refitML = function(x, ...) stop("the refit fails"),
                                  .package = "lme4")
  expect_error(maihda_ic_refit_ml(fm), "the refit fails", fixed = TRUE)
  model <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))
  expect_identical(maihda_ic_one(model, ml = TRUE)$estimator, "REML")
  # What a refit raised before it stopped with an error is passed on with the error.
  testthat::local_mocked_bindings(refitML = function(x, ...) {
    warning("a warning before the refit fails")
    message("a message before the refit fails")
    stop("the refit fails")
  }, .package = "lme4")
  h <- ic16_heard(tryCatch(maihda_ic_refit_ml(fm), error = function(e) conditionMessage(e)))
  expect_identical(h$value, "the refit fails")
  expect_identical(h$warnings, "a warning before the refit fails")
  expect_identical(h$messages, "a message before the refit fails\n")
})

test_that("a refit without a code is returned as lme4 returns it, and only the exit code is dropped from a checked one", {
  skip_on_cran()
  d <- ic16_strata_data()
  fm <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))$model
  real <- lme4::refitML
  calls <- character(0)
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    calls <<- c(calls, optimizer)
    warning("a warning of the refit")
    message("a message of the refit")
    real(x, optimizer = optimizer, ...)
  }, .package = "lme4")
  h <- ic16_heard(maihda_ic_refit_ml(fm))
  # No code: no second refit, and what the refit said is passed on.
  expect_identical(calls, "bobyqa")
  expect_identical(h$warnings, "a warning of the refit")
  expect_identical(h$messages, "a message of the refit\n")
  # A checked refit keeps its other warnings and its messages.
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    refit <- real(x, optimizer = optimizer, ...)
    if (!identical(optimizer, "bobyqa")) {
      return(refit)
    }
    warning("a warning of the refit")
    message("a message of the refit")
    ic16_flagged(refit, optimizer)
  }, .package = "lme4")
  h <- ic16_heard(maihda_ic_refit_ml(fm))
  expect_identical(h$warnings, "a warning of the refit")
  expect_identical(h$messages, "a message of the refit\n")
  # The warning that is dropped is lme4's for an optimiser's exit code, and nothing else.
  is_code <- function(text, class = simpleWarning) maihda_is_optimiser_code_warning(class(text))
  expect_true(is_code(code_warning))
  expect_true(is_code("convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded"))
  expect_true(is_code("convergence code -4 from nloptwrap: NLOPT_ROUNDOFF_LIMITED: Roundoff errors led to a breakdown of the optimization algorithm."))
  expect_true(is_code("convergence code 4 from Nelder_Mead"))
  expect_false(is_code("Model failed to converge with max|grad| = 0.0025 (tol = 0.002, component 1)"))
  expect_false(is_code("the convergence code 3 from bobyqa was reported"))
  expect_false(is_code(code_warning, class = simpleMessage))
})

# ---- the criteria ---------------------------------------------------------------------------

test_that("maihda_ic() and compare_maihda() read the criteria from the checked refit", {
  skip_on_cran()
  d <- ic16_strata_data()
  m0 <- ic16_quiet(fit_maihda(y ~ 1 + (1 | stratum), data = d))
  m1 <- ic16_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))
  reference <- ic16_quiet(maihda_ic(m0, m1))
  real <- lme4::refitML
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    refit <- real(x, optimizer = optimizer, ...)
    if (identical(optimizer, "bobyqa")) ic16_flagged(refit, optimizer) else refit
  }, .package = "lme4")
  ic <- ic16_heard(maihda_ic(m0, m1))
  expect_identical(ic$warnings, character(0))
  expect_identical(ic$value$estimator, rep("ML (refit from REML)", 2L))
  expect_identical(ic$value$AIC, reference$AIC)
  expect_identical(ic$value$BIC, reference$BIC)
  expect_true("delta" %in% names(ic$value))
  cmp <- ic16_heard(compare_maihda(m0, m1))
  expect_identical(cmp$warnings, character(0))
  expect_identical(cmp$value$AIC, reference$AIC)
  # Where a warning is an error, the refit is still used: no warning is raised for it.
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  expect_identical(maihda_ic_one(m1, ml = TRUE)$estimator, "ML (refit from REML)")
})

test_that("no exit code of lme4's optimiser reaches the user from fits with a variance at zero", {
  skip_on_cran()
  # Five data sets of the design that failed. On every build of lme4 and minqa measured, the
  # refit of at least one of the ten fits is flagged (seed 1 where the compiler does not fuse
  # multiply-adds, seed 2 where it does, seed 4 where it does in one of the two packages), so
  # the check runs here whatever the build.
  for (seed in 1:5) {
    d <- ic16_noise_data(seed)
    m_ab <- ic16_quiet(fit_maihda(y ~ age + (1 | a:b), data = d))
    m_ac <- ic16_quiet(fit_maihda(y ~ age + (1 | a:cc), data = d))
    cmp <- ic16_heard(compare_maihda(m_ab, m_ac))
    # The one warning is compare_maihda()'s own, about the strata.
    expect_length(cmp$warnings, 1L)
    expect_match(cmp$warnings, "stratum definitions", fixed = TRUE)
    ic <- ic16_heard(maihda_ic(m_ab, m_ac))
    expect_false(any(grepl("convergence code", ic$warnings, fixed = TRUE)))
    # The criteria are those of the models fitted by maximum likelihood.
    direct <- vapply(list(m_ab, m_ac), function(m) {
      stats::AIC(ic16_quiet(lme4::lmer(stats::formula(m$model), data = m$model@frame, REML = FALSE)))
    }, numeric(1))
    expect_lt(max(abs(cmp$value$AIC - direct)), 1e-6)
    expect_lt(max(abs(ic$value$AIC - direct)), 1e-6)
  }
})
