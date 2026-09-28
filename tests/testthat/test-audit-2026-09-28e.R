# Audit 2026-09-28e: the follow-up to F02 (test-audit-2026-09-28d.R). A brms trials()
# term that calls a function outside base -- a global function, or one from an
# attached package such as coalesce() under library(dplyr) -- lost its trial counts
# everywhere the package reads them itself.
#
# maihda_trials_from_formula() evaluated the term against base functions only
# (enclos = baseenv(), whose parent is the EMPTY environment), while brms finds the
# functions of an addition term through the global environment and the attached
# packages -- and refuses any VARIABLE it cannot find in its data. So brms fitted
# y | trials(plus_one(n)) and the package read no trials at all. Measured on two real
# Stan fits of one model, trials(plus_one(n)) and the same trials precomputed as a
# column m, with IDENTICAL draws: predict_maihda(scale = "response") returned expected
# success counts (0.39 - 12) where the column spelling gives probabilities
# (0.049 - 0.80), on newdata too; the prediction weights were all 1 instead of the
# trials (7 - 15); plot_obs_vs_shrunken() drew mean success counts on its observed
# axis (0.88 - 6.8 against 0.084 - 0.61) and its shrunken axis, plot_predicted_strata(),
# maihda_table(), maihda_interactions(scale = "response") and
# plot_effect_decomposition() moved with the unit weights; maihda_describe() -- on the
# fit and on the formula -- counted one trial a row, stratum proportions up to 6.23
# too high; and the AUC and the deviation panel refused outright.
#
# The same base-only evaluation read the RESPONSE of such a fit for the deviation
# panel, so a response calling a global function -- a28e_floor(raw) | trials(n), which
# brms fits -- left the panel no row to score (all 12 strata without a residual, and a
# warning). The AUC, predict_maihda(), maihda_describe() and plot_obs_vs_shrunken()
# take the evaluated response from brms's own frame and were already right.
#
# FIX: maihda_eval_brms_term() evaluates a term as brms does -- variables from the data
# alone, so a missing column is still no trials rather than a same-named object from
# the session; functions base first, then the global environment and the attached
# packages -- and serves the trials() readers and the panel's successes alike. Every
# output above then equals the column spelling's exactly.

a28e_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
a28e_or_null <- function(expr) tryCatch(expr, error = function(e) NULL)

# Assign `value` as `name` in the global environment -- where brms must find a function
# its formula calls -- skipping the test when the name is already taken there. Returns
# the undo, for the test to register with on.exit().
a28e_global <- function(name, value) {
  skip_if(exists(name, envir = globalenv(), inherits = FALSE),
          paste("the global environment already has an object named", name))
  assign(name, value, envir = globalenv())
  function() {
    if (exists(name, envir = globalenv(), inherits = FALSE)) {
      rm(list = name, envir = globalenv())
    }
  }
}

a28e_data <- function() {
  set.seed(20260928)
  d <- data.frame(stratum = rep(seq_len(12), each = 25))
  d$x <- stats::rnorm(nrow(d))
  u <- stats::rnorm(12, 0, 0.7)
  d$n <- sample(6:14, nrow(d), replace = TRUE)
  d$y <- stats::rbinom(nrow(d), d$n, stats::plogis(-0.3 + 0.4 * d$x + u[d$stratum]))
  d$m <- d$n + 1L
  # The same successes as a column to transform: floor(raw) == y.
  d$raw <- d$y + 0.25
  d
}

test_that("a trials() term's functions are found as brms finds them", {
  # Assigns into the global environment, as brms requires, so not on CRAN.
  skip_on_cran()
  skip_if_not("package:stats" %in% search(), "stats is not attached")
  undo <- a28e_global("a28e_plus_one", function(z) z + 1L)
  on.exit(undo(), add = TRUE)
  d <- data.frame(y = 1:4, n = c(2, 4, 6, 8))
  # A global function, as brms fits it.
  expect_equal(a28e_or_null(maihda_eval_brms_term(quote(a28e_plus_one(n)), d)), c(3, 5, 7, 9))
  expect_equal(maihda_trials_from_formula(y | trials(a28e_plus_one(n)) ~ x, d), c(3, 5, 7, 9))
  # A function from an attached package: median() is stats, not base.
  expect_equal(maihda_trials_from_formula(y | trials(pmax(n, median(n))) ~ x, d),
               c(5, 5, 6, 8))
  # Controls: a column, a constant, base arithmetic.
  expect_equal(maihda_trials_from_formula(y | trials(n) ~ x, d), c(2, 4, 6, 8))
  expect_equal(maihda_trials_from_formula(y | trials(12) ~ x, d), rep(12, 4))
  expect_equal(maihda_trials_from_formula(y | trials(2 * n) ~ x, d), c(4, 8, 12, 16))
})

test_that("a term's variables still come from the data alone", {
  skip_on_cran()
  undo1 <- a28e_global("a28e_plus_one", function(z) z + 1L)
  on.exit(undo1(), add = TRUE)
  undo2 <- a28e_global("a28e_stray", c(10, 20, 30, 40))
  on.exit(undo2(), add = TRUE)
  d <- data.frame(y = 1:4, n = c(2, 4, 6, 8))
  # A stray global object is never read as a missing column ...
  expect_null(maihda_trials_from_formula(y | trials(a28e_stray) ~ x, d))
  expect_null(maihda_trials_from_formula(y | trials(a28e_stray + n) ~ x, d))
  # ... nor is a global FUNCTION, when the term uses its name as a variable.
  expect_null(maihda_trials_from_formula(y | trials(a28e_plus_one) ~ x, d))
  # A column of the data shadows a global object of the same name.
  d$a28e_stray <- c(1, 1, 1, 1)
  expect_equal(maihda_trials_from_formula(y | trials(a28e_stray) ~ x, d), rep(1, 4))
})

test_that("base functions come first, as they do for brms", {
  # brms fits trials(round(n)) with base::round under the same mask (measured).
  skip_on_cran()
  undo <- a28e_global("round", function(x, digits = 0) 999)
  on.exit(undo(), add = TRUE)
  d <- data.frame(y = 1:4, n = c(2.4, 4.6, 6, 8))
  expect_equal(maihda_trials_from_formula(y | trials(round(n)) ~ x, d), c(2, 5, 6, 8))
})

test_that("a brms fit with such a term keeps its trials everywhere it is read (Stan-free)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  undo <- a28e_global("a28e_plus_one", function(z) z + 1L)
  on.exit(undo(), add = TRUE)
  d <- a28e_data()
  fit <- function(f) a28e_or_null(a28e_quiet(fit_maihda(f, data = d, engine = "brms",
                                                        family = "binomial", empty = TRUE)))
  mF <- fit(y | trials(a28e_plus_one(n)) ~ x + (1 | stratum))
  mM <- fit(y | trials(m) ~ x + (1 | stratum))
  expect_s3_class(mF, "maihda_model")
  expect_equal(a28e_or_null(maihda_brms_trial_counts(mF)), as.numeric(d$m))
  expect_equal(a28e_or_null(maihda_prediction_weights(mF)), as.numeric(d$m))
  outcome <- function(x, ...) {
    desc <- a28e_or_null(a28e_quiet(maihda_describe(x, ...)))
    if (is.null(desc)) return(NULL)
    desc$strata[, c("outcome_events", "outcome_trials", "outcome_proportion")]
  }
  expect_equal(outcome(mF), outcome(mM))
  expect_equal(outcome(y | trials(a28e_plus_one(n)) ~ x + (1 | stratum), data = d,
                       family = "binomial"),
               outcome(y | trials(m) ~ x + (1 | stratum), data = d, family = "binomial"))

  # The response side: the deviation panel read no successes for a response calling
  # a global function, so it scored no row at all.
  undo2 <- a28e_global("a28e_floor", function(z) floor(z))
  on.exit(undo2(), add = TRUE)
  mR <- fit(a28e_floor(raw) | trials(n) ~ x + (1 | stratum))
  expect_s3_class(mR, "maihda_model")
  expect_equal(a28e_or_null(maihda_prediction_panel_brms_successes(mR$model, mR$data)),
               as.numeric(d$y))
  expect_equal(a28e_or_null(maihda_prediction_panel_brms_successes(mR$model, d)),
               as.numeric(d$y))
})

test_that("real brms fits: the function spelling reads as the column spelling (opt-in Stan)", {
  # Compiles two Stan models, so OPT-IN (set MAIHDA_TEST_BRMS=true). The two fits are one
  # model with identical draws -- the response and the trials each calling a global
  # function, against both precomputed -- so every output must agree exactly.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  undo <- a28e_global("a28e_plus_one", function(z) z + 1L)
  on.exit(undo(), add = TRUE)
  undo2 <- a28e_global("a28e_floor", function(z) floor(z))
  on.exit(undo2(), add = TRUE)
  d <- a28e_data()
  fit <- function(f) a28e_quiet(fit_maihda(f, data = d, engine = "brms", family = "binomial",
                                           chains = 1, iter = 600, warmup = 300,
                                           refresh = 0, seed = 7))
  mF <- fit(a28e_floor(raw) | trials(a28e_plus_one(n)) ~ x + (1 | stratum))
  mM <- fit(y | trials(m) ~ x + (1 | stratum))
  expect_identical(as.matrix(mF$model), as.matrix(mM$model))
  sF <- a28e_quiet(summary(mF))
  sM <- a28e_quiet(summary(mM))

  pF <- as.numeric(a28e_quiet(predict_maihda(mF, scale = "response")))
  expect_equal(pF, as.numeric(a28e_quiet(predict_maihda(mM, scale = "response"))),
               tolerance = 1e-12)
  expect_true(all(pF > 0 & pF < 1))
  nd <- d[1:40, ]
  expect_equal(as.numeric(a28e_quiet(predict_maihda(mF, newdata = nd, scale = "response"))),
               as.numeric(a28e_quiet(predict_maihda(mM, newdata = nd, scale = "response"))),
               tolerance = 1e-12)
  expect_equal(a28e_or_null(a28e_quiet(maihda_discriminatory_accuracy(mF))$auc),
               a28e_quiet(maihda_discriminatory_accuracy(mM))$auc, tolerance = 1e-12)
  oF <- a28e_quiet(plot_obs_vs_shrunken(mF, sF))$data
  oM <- a28e_quiet(plot_obs_vs_shrunken(mM, sM))$data
  expect_equal(oF$observed, oM$observed, tolerance = 1e-12)
  expect_equal(oF$shrunken, oM$shrunken, tolerance = 1e-12)
  panel <- function(m) a28e_or_null(a28e_quiet(
    plot_prediction_deviation_panels(m, type = "binomial"))[[2]]$data)
  gF <- panel(mF)
  gM <- panel(mM)
  expect_equal(gF$fitted, gM$fitted, tolerance = 1e-12)
  expect_equal(gF$abs_res_dev, gM$abs_res_dev, tolerance = 1e-12)
  expect_true(all(is.finite(gF$abs_res_dev)))
  expect_equal(a28e_quiet(maihda_table(mF))$strata$predicted,
               a28e_quiet(maihda_table(mM))$strata$predicted, tolerance = 1e-12)
})
