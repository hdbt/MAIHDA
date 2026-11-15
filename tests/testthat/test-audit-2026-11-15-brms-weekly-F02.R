# Audit 2026-11-15, batch brms-weekly, finding F02: maihda_ic() and
# compare_maihda(ic = TRUE) returned a brms fit's row without WAIC and LOOIC, and
# said nothing, when brms could not compute them. CONFIRMED.
#
# brms::waic() and brms::loo() both end in brms's hash_response(), which starts with
# require_package("digest"); digest is a package brms suggests, so a library can hold
# brms without it. maihda_ic_one() caught the error of either call and left NA: the
# table came back with model, n and estimator only, and compare_maihda() with its
# VPC columns only. Seen on two real Stan fits with digest made unloadable (brms
# 2.23.0), and on any brmsfit brms stops on for another reason.
#
# FIX: maihda_ic_one() keeps the error of each criterion brms stopped on (attribute
# "ic_unavailable" of its row), and maihda_ic() and compare_maihda() each raise one
# warning a call from it (maihda_ic_warn_unavailable()): the criteria, the models and
# the error as brms gave it. The table is as before. digest is in Suggests.
#
# The fits are brm(empty = TRUE): genuine brmsfit objects with no Stan compile. brms
# stops on them by itself (they hold no draws), and brms::waic() / brms::loo() are
# replaced where a block needs them to stop as they do without digest, or to work.

# Run `expr` and put R's random-number state back: setup-seed.R fixes one stream for
# the whole suite, and this file must leave it where it found it.
f02_keep_seed <- function(expr) {
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

f02_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression and the text of every warning it raised.
f02_warnings <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

# What a call stops with where a warning is an error (options(warn = 2)); NA when it
# returns.
f02_error_at_warn2 <- function(expr) {
  old <- options(warn = 2)
  on.exit(options(old))
  tryCatch({
    expr
    NA_character_
  }, error = function(e) conditionMessage(e))
}

# The warnings about criteria brms could not compute, among those raised.
f02_unavailable <- function(warnings) {
  warnings[grepl("brms could not compute", warnings, fixed = TRUE)]
}

# 6 strata x 10 rows, no random numbers.
f02_data <- function() {
  d <- data.frame(stratum = rep(seq_len(6), each = 10))
  d$x <- rep(c(-1.5, -0.9, -0.4, -0.1, 0.2, 0.3, 0.6, 0.8, 1.1, 1.7), 6)
  d$y <- 1 + 0.5 * d$x + c(-0.9, -0.5, -0.1, 0.2, 0.5, 0.8)[d$stratum] +
    0.6 * sin(1.7 * seq_len(nrow(d)))
  d
}

# Two brms fits without draws (null and adjusted) and the same two models in lme4,
# built once.
f02_models <- local({
  made <- NULL
  function() {
    if (is.null(made)) {
      d <- f02_data()
      brms_fit <- function(f) {
        f02_keep_seed(f02_quiet(fit_maihda(f, data = d, engine = "brms", empty = TRUE)))
      }
      made <<- list(
        b0 = brms_fit(y ~ 1 + (1 | stratum)),
        b1 = brms_fit(y ~ x + (1 | stratum)),
        l0 = f02_keep_seed(f02_quiet(fit_maihda(y ~ 1 + (1 | stratum), data = d))),
        l1 = f02_keep_seed(f02_quiet(fit_maihda(y ~ x + (1 | stratum), data = d))))
    }
    made
  }
})

# brms's words when 'digest' is not installed (brms:::require_package()).
f02_digest <- "Please install the 'digest' package."

# A stand-in for brms::waic() or brms::loo() that stops as brms does there.
f02_stops <- function(message = f02_digest) {
  function(x, ...) stop(message, call. = FALSE)
}

# A stand-in that returns what maihda_ic_one() reads off the criterion, one value a
# call in the order given.
f02_returns <- function(criterion, values) {
  i <- 0L
  function(x, ...) {
    i <<- i + 1L
    list(estimates = matrix(c(values[[i]], 1), nrow = 1,
                            dimnames = list(criterion, c("Estimate", "SE"))))
  }
}

test_that("maihda_ic() warns once, with the models and brms's error, when neither criterion computes", {
  skip_on_cran()
  skip_if_not_installed("brms")
  m <- f02_models()
  local_mocked_bindings(waic = f02_stops(), loo = f02_stops(), .package = "brms")

  one <- f02_warnings(maihda_ic(m$b0))
  expect_length(one$warnings, 1L)
  expect_identical(
    one$warnings,
    paste0("maihda_ic(): brms could not compute WAIC and LOOIC for 'Model1' ",
           "(missing from the table). The error was: ", f02_digest))
  # The table is the one returned before, criteria left out.
  expect_s3_class(one$value, "maihda_ic")
  expect_identical(names(one$value), c("model", "n", "estimator"))
  expect_identical(one$value$estimator, "Bayesian")
  expect_identical(one$value$n, 60L)
  expect_true(is.na(attr(one$value, "ic_primary")))

  # Two models, two criteria each: still one warning.
  two <- f02_warnings(maihda_ic(m$b0, m$b1, model_names = c("Null", "Adjusted")))
  expect_length(two$warnings, 1L)
  expect_identical(
    two$warnings,
    paste0("maihda_ic(): brms could not compute WAIC and LOOIC for 'Null' and ",
           "'Adjusted' (missing from the table). The error was: ", f02_digest))
  expect_identical(names(two$value), c("model", "n", "estimator"))
  expect_identical(two$value$model, c("Null", "Adjusted"))

  # A maihda() analysis is named by its rows.
  analysis <- structure(list(model = m$b0, model_adjusted = m$b1),
                        class = "maihda_analysis")
  rows <- f02_warnings(maihda_ic(analysis))
  expect_length(rows$warnings, 1L)
  expect_match(rows$warnings, "for 'Model1 (Null)' and 'Model1 (Adjusted)' (missing",
               fixed = TRUE)

  # Where a warning is an error, the call stops with this one.
  expect_match(
    f02_error_at_warn2(maihda_ic(m$b0)),
    paste0("maihda_ic(): brms could not compute WAIC and LOOIC for 'Model1' ",
           "(missing from the table). The error was: ", f02_digest),
    fixed = TRUE)
})

test_that("only the criterion brms stopped on is named, and the other is reported", {
  skip_on_cran()
  skip_if_not_installed("brms")
  m <- f02_models()
  local_mocked_bindings(waic = f02_returns("waic", c(210, 200)),
                        loo = f02_stops("no LOO today"), .package = "brms")
  got <- f02_warnings(maihda_ic(m$b0, m$b1))
  expect_identical(
    got$warnings,
    paste0("maihda_ic(): brms could not compute LOOIC for 'Model1' and 'Model2' ",
           "(missing from the table). The error was: no LOO today"))
  expect_identical(names(got$value), c("model", "n", "estimator", "WAIC", "delta"))
  expect_equal(got$value$WAIC, c(210, 200))
  expect_identical(attr(got$value, "ic_primary"), "WAIC")
  expect_equal(got$value$delta, c(10, 0))

  # One model of two: the column stays, with NA in that model's row.
  calls <- 0L
  local_mocked_bindings(
    waic = f02_returns("waic", c(210, 200)),
    loo = function(x, ...) {
      calls <<- calls + 1L
      if (calls == 2L) stop(f02_digest, call. = FALSE)
      list(estimates = matrix(c(212, 1), nrow = 1,
                              dimnames = list("looic", c("Estimate", "SE"))))
    },
    .package = "brms")
  got <- f02_warnings(maihda_ic(m$b0, m$b1))
  expect_identical(
    got$warnings,
    paste0("maihda_ic(): brms could not compute LOOIC for 'Model2' ",
           "(missing from the table). The error was: ", f02_digest))
  expect_equal(got$value$LOOIC, c(212, NA))
  expect_equal(got$value$WAIC, c(210, 200))
})

test_that("models stopped by different errors are listed apart, in one warning", {
  row <- function(...) {
    out <- data.frame(n = 1L)
    stopped <- c(...)
    if (length(stopped) > 0L) attr(out, "ic_unavailable") <- stopped
    out
  }
  told <- function(rows, names) {
    f02_warnings(maihda_ic_warn_unavailable(rows, names, "maihda_ic"))
  }

  got <- told(list(row(WAIC = "first", LOOIC = "first"), row(),
                   row(LOOIC = "second"), row(WAIC = "first", LOOIC = "first")),
              c("a", "b", "c", "d"))
  expect_true(got$value)
  expect_identical(
    got$warnings,
    paste0("maihda_ic(): brms could not compute WAIC and LOOIC for 'a' and 'd'; ",
           "LOOIC for 'c' (missing from the table). The errors were, in that ",
           "order: first | second"))

  # One model, a different error for each criterion.
  got <- told(list(row(WAIC = "first", LOOIC = "second")), "a")
  expect_identical(
    got$warnings,
    paste0("maihda_ic(): brms could not compute WAIC for 'a'; LOOIC for 'a' ",
           "(missing from the table). The errors were, in that order: first | second"))

  # The same criterion stopped by different errors: a clause for each model.
  got <- told(list(row(LOOIC = "first"), row(LOOIC = "second")), c("a", "b"))
  expect_identical(
    got$warnings,
    paste0("maihda_ic(): brms could not compute LOOIC for 'a'; LOOIC for 'b' ",
           "(missing from the table). The errors were, in that order: first | second"))

  # One error behind different criteria is quoted once.
  got <- told(list(row(WAIC = "first", LOOIC = "first"), row(LOOIC = "first")),
              c("a", "b"))
  expect_identical(
    got$warnings,
    paste0("maihda_ic(): brms could not compute WAIC and LOOIC for 'a'; LOOIC for 'b' ",
           "(missing from the table). The error was: first"))

  # An error over several lines is quoted on one, and one without words is named.
  got <- told(list(row(WAIC = "  line one\n   line two \n"), row(LOOIC = "")),
              c("a", "b"))
  expect_match(got$warnings, "in that order: line one line two | an error without a message",
               fixed = TRUE)

  # Rows that carry nothing: no warning, and the answer says so.
  got <- told(list(row(), row()), c("a", "b"))
  expect_false(got$value)
  expect_length(got$warnings, 0L)
  expect_false(f02_warnings(maihda_ic_warn_unavailable(NULL, character(0), "x"))$value)
})

test_that("a brmsfit without draws: the warning quotes what brms itself says", {
  skip_on_cran()
  skip_if_not_installed("brms")
  m <- f02_models()
  # Nothing is replaced here: brms stops on a fit that holds no draws.
  says <- function(expr) {
    gsub("\\s+", " ", trimws(tryCatch({
      f02_quiet(expr)
      NA_character_
    }, error = function(e) conditionMessage(e))))
  }
  waic_says <- says(brms::waic(m$b0$model))
  loo_says <- says(brms::loo(m$b0$model))
  expect_false(is.na(waic_says))
  expect_false(is.na(loo_says))

  row <- maihda_ic_one(m$b0)
  expect_true(is.na(row$WAIC) && is.na(row$LOOIC))
  expect_identical(row$estimator, "Bayesian")
  expect_identical(gsub("\\s+", " ", trimws(attr(row, "ic_unavailable", exact = TRUE))),
                   c(WAIC = waic_says, LOOIC = loo_says))

  got <- f02_warnings(maihda_ic(m$b0))
  expect_length(got$warnings, 1L)
  expect_true(grepl(waic_says, got$warnings, fixed = TRUE))
  expect_true(grepl("for 'Model1' (missing from the table)", got$warnings, fixed = TRUE))
  expect_identical(names(got$value), c("model", "n", "estimator"))
})

test_that("compare_maihda() warns once and returns its VPC columns", {
  skip_on_cran()
  skip_if_not_installed("brms")
  m <- f02_models()
  # summary() needs draws, so the VPCs are those of the two lme4 fits, and the
  # criteria rows are maihda_ic_one()'s own for the two brms fits, in that order.
  one_row <- maihda_ic_one
  brms_fits <- list(m$b0, m$b1)
  compare <- function(..., collect = f02_warnings) {
    at <- 0L
    local_mocked_bindings(maihda_ic_one = function(model, ml = FALSE) {
      at <<- at + 1L
      one_row(brms_fits[[at]], ml = ml)
    })
    collect(compare_maihda(m$l0, m$l1, ...))
  }
  plain <- f02_warnings(compare_maihda(m$l0, m$l1, model_names = c("Null", "Adjusted"),
                                       ic = FALSE))
  expect_length(plain$warnings, 0L)

  local_mocked_bindings(waic = f02_stops(), loo = f02_stops(), .package = "brms")
  got <- compare(model_names = c("Null", "Adjusted"))
  expect_identical(
    got$warnings,
    paste0("compare_maihda(): brms could not compute WAIC and LOOIC for 'Null' and ",
           "'Adjusted' (missing from the table). The error was: ", f02_digest))
  expect_s3_class(got$value, "maihda_comparison")
  expect_identical(got$value, plain$value)

  # One criterion computes: its column is appended, and the VPCs are the same.
  local_mocked_bindings(waic = f02_returns("waic", c(210, 200)), loo = f02_stops(),
                        .package = "brms")
  got <- compare()
  expect_identical(
    got$warnings,
    paste0("compare_maihda(): brms could not compute LOOIC for 'Model1' and 'Model2' ",
           "(missing from the table). The error was: ", f02_digest))
  expect_identical(names(got$value), c("model", "vpc", "WAIC"))
  expect_identical(got$value$vpc, plain$value$vpc)
  expect_equal(got$value$WAIC, c(210, 200))

  # ic = FALSE asks brms for nothing.
  local_mocked_bindings(waic = f02_stops(), loo = f02_stops(), .package = "brms")
  expect_length(compare(ic = FALSE)$warnings, 0L)

  # Where a warning is an error, the call stops with this one: it is raised outside
  # the guard around the criteria, which would catch it there and return the VPC
  # columns as if nothing were missing.
  expect_match(
    compare(collect = f02_error_at_warn2),
    paste0("compare_maihda(): brms could not compute WAIC and LOOIC for 'Model1' and ",
           "'Model2' (missing from the table). The error was: ", f02_digest),
    fixed = TRUE)
})

test_that("no warning where the criteria compute, or where NA is the answer", {
  skip_on_cran()
  skip_if_not_installed("brms")
  m <- f02_models()

  # brms criteria that compute (three calls each: two models, then one row).
  local_mocked_bindings(waic = f02_returns("waic", c(210, 200, 210)),
                        loo = f02_returns("looic", c(212, 203, 212)), .package = "brms")
  got <- f02_warnings(maihda_ic(m$b0, m$b1))
  expect_length(got$warnings, 0L)
  expect_identical(names(got$value),
                   c("model", "n", "estimator", "WAIC", "LOOIC", "delta"))
  expect_equal(got$value$LOOIC, c(212, 203))
  expect_equal(got$value$delta, c(9, 0))
  expect_null(attr(maihda_ic_one(m$b0), "ic_unavailable", exact = TRUE))

  # From here on brms stops whenever it is asked.
  local_mocked_bindings(waic = f02_stops(), loo = f02_stops(), .package = "brms")

  # A sampling-weighted brms fit is not asked: its NA is the estimator's.
  weighted <- m$b0
  weighted$sampling_weights <- "w"
  got <- f02_warnings(maihda_ic(weighted))
  expect_length(got$warnings, 0L)
  expect_identical(got$value$estimator, "Bayesian (weighted pseudo-posterior)")
  expect_null(attr(maihda_ic_one(weighted), "ic_unavailable", exact = TRUE))

  # lme4 fits: the same tables, no warning.
  got <- f02_warnings(maihda_ic(m$l0, m$l1))
  expect_length(got$warnings, 0L)
  expect_true(all(is.finite(got$value$AIC)))
  expect_null(attr(maihda_ic_one(m$l0, ml = TRUE), "ic_unavailable", exact = TRUE))
  got <- f02_warnings(compare_maihda(m$l0, m$l1))
  expect_length(got$warnings, 0L)
  expect_identical(names(got$value), c("model", "vpc", "AIC", "BIC"))

  # Beside an lme4 fit, only the brms model is named.
  got <- f02_warnings(maihda_ic(m$l0, m$b0))
  expect_length(got$warnings, 1L)
  expect_match(got$warnings, "WAIC and LOOIC for 'Model2' (missing", fixed = TRUE)
  expect_true(is.finite(got$value$AIC[1]))
})
