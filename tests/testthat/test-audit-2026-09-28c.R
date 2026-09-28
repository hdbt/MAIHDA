# Audit 2026-09-28 (c): the binomial prediction-deviation panel scored no row of a fit
# whose response is a factor with three or more levels, on the model's own rows --
# a REGRESSION in the then-uncommitted work of test-audit-2026-09-27b.R/-27c.R/-28b.R,
# found by merging the case-level shape fix (test-audit-2026-09-28.R) into it.
#
# R's binomial family takes a factor response's first level as the failure and EVERY
# other level as the event, so glm(), glmer() and fit_maihda(family = "binomial") all
# fit a three-level factor. The panel codes supplied rows as the model coded its
# outcome (maihda_prediction_panel_fit_01(), which knows that rule), but coded the
# fitted rows by their own values with maihda_binomial_observed_01(), on the premise
# that they "carry the model's own coding already" -- which codes only two values. So
# an unweighted fit of a three-level factor had every row unscored on its own rows:
# no residual, no label, and a warning that no row could be scored -- 0 of 150 rows
# for a glm or a glmer, 0 of 10 and 0 of 6 strata for a glmer with a stratum and a
# fit_maihda() model -- while the same rows passed as `data` were scored in full, and
# HEAD had scored them with the fit's own deviance residuals. A weighted fit escaped
# (its rows are scored from successes out of trials), and so did every two-valued
# outcome.
#
# FIX: the fitted rows are coded as the model coded its outcome too. The residuals are
# the fit's own |deviance residuals| again (to 3e-16), the same five rows are labelled
# as at HEAD, and the plotted rows now agree whether or not they are passed as `data`
# -- the Wrong / Correct rings included, which the fitted rows gain. Every panel of a
# two-valued, aggregated or proportion outcome is unchanged: 19 of 25 probe panels
# identical() before and after, the other 6 being the multi-level factor fits.

c28_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression and the warnings it raised, muffled.
c28_capture <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr,
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  list(value = value, warnings = w)
}

c28_data <- function() {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(150), grp = factor(rep(1:10, each = 15)))
  d$y3 <- factor(sample(c("lo", "mid", "hi"), 150, replace = TRUE),
                 levels = c("lo", "mid", "hi"))
  d$y7 <- factor(sample(letters[1:7], 150, replace = TRUE))
  d$y <- stats::rbinom(150, 1, stats::plogis(0.8 * d$x))
  d$yr <- factor(ifelse(d$y == 1, "yes", "no"), levels = c("yes", "no"))
  d$w <- sample(1:4, 150, replace = TRUE)
  d
}

# |deviance residual| of each panel row (case level, by id), from the fit itself.
c28_row_oracle <- function(fit, pd) {
  abs(as.numeric(stats::residuals(fit, type = "deviance")))[pd$id]
}

# The prior-weighted stratum mean of the fit's own |deviance residuals|.
c28_stratum_oracle <- function(fit, strat, pd) {
  r <- abs(as.numeric(stats::residuals(fit, type = "deviance")))
  w <- as.numeric(stats::weights(fit, type = "prior"))
  m <- vapply(split(seq_along(r), as.character(strat)),
              function(i) sum(w[i] * r[i]) / sum(w[i]), numeric(1))
  as.numeric(m[as.character(pd$stratum)])
}

c28_labels <- function(p) {
  lab <- Filter(function(l) inherits(l$geom, "GeomLabelRepel"), p[[2]]$layers)[[1]]$data
  if ("stratum" %in% names(lab)) as.character(lab$stratum) else lab$id
}

c28_no_unscored_warning <- function(w) {
  expect_identical(grep("can be scored", w, value = TRUE), character(0))
}

test_that("a glm fit of a three-level factor is scored on its own rows", {
  d <- c28_data()
  g <- stats::glm(y3 ~ x, data = d, family = stats::binomial)
  r <- c28_capture(plot_prediction_deviation_panels(g))
  c28_no_unscored_warning(r$warnings)
  pd <- r$value[[2]]$data
  expect_true(all(is.finite(pd$abs_res_dev)))
  expect_equal(as.numeric(pd$abs_res_dev), c28_row_oracle(g, pd), tolerance = 1e-12)
  # Coded as glm coded it: the first level against the other two.
  expect_identical(pd$obs_outcome_01, as.integer(d$y3 != "lo")[pd$id])
  # The labelled rows are the fit's five worst.
  worst <- order(-abs(as.numeric(stats::residuals(g, type = "deviance"))))[1:5]
  expect_setequal(c28_labels(r$value), worst)

  # The same rows given as `data` -- the fitted frame, or the original columns --
  # are the same panel, rings included.
  for (dd in list(stats::model.frame(g), d)) {
    ps <- c28_quiet(plot_prediction_deviation_panels(g, data = dd))[[2]]$data
    ps <- ps[order(ps$id), ]
    po <- pd[order(pd$id), ]
    expect_equal(as.numeric(ps$abs_res_dev), as.numeric(po$abs_res_dev), tolerance = 1e-12)
    expect_identical(ps$obs_outcome_01, po$obs_outcome_01)
    expect_identical(ps$wrong, po$wrong)
  }

  # Seven levels are coded the same way.
  g7 <- stats::glm(y7 ~ x, data = d, family = stats::binomial)
  r7 <- c28_capture(plot_prediction_deviation_panels(g7))
  c28_no_unscored_warning(r7$warnings)
  p7 <- r7$value[[2]]$data
  expect_equal(as.numeric(p7$abs_res_dev), c28_row_oracle(g7, p7), tolerance = 1e-12)
  expect_identical(p7$obs_outcome_01, as.integer(d$y7 != "a")[p7$id])
})

test_that("glmer fits of one, at case and stratum level, and a fit_maihda() model", {
  d <- c28_data()
  m <- c28_quiet(lme4::glmer(y3 ~ x + (1 | grp), data = d, family = stats::binomial))
  r <- c28_capture(plot_prediction_deviation_panels(m))
  c28_no_unscored_warning(r$warnings)
  pd <- r$value[[2]]$data
  expect_equal(as.numeric(pd$abs_res_dev), c28_row_oracle(m, pd), tolerance = 1e-12)

  ds <- d
  names(ds)[names(ds) == "grp"] <- "stratum"
  ms <- c28_quiet(lme4::glmer(y3 ~ x + (1 | stratum), data = ds, family = stats::binomial))
  r <- c28_capture(plot_prediction_deviation_panels(ms))
  c28_no_unscored_warning(r$warnings)
  pd <- r$value[[2]]$data
  expect_identical(nrow(pd), 10L)
  expect_equal(as.numeric(pd$abs_res_dev), c28_stratum_oracle(ms, ds$stratum, pd),
               tolerance = 1e-12)

  dm <- d
  dm$g1 <- factor(sample(c("A", "B"), 150, replace = TRUE))
  dm$g2 <- factor(sample(c("C", "D", "E"), 150, replace = TRUE))
  fm <- c28_quiet(fit_maihda(y3 ~ x + (1 | g1:g2), data = dm, family = "binomial"))
  r <- c28_capture(plot_prediction_deviation_panels(fm))
  c28_no_unscored_warning(r$warnings)
  pd <- r$value[[2]]$data
  expect_true(all(is.finite(pd$abs_res_dev)))
  expect_equal(as.numeric(pd$abs_res_dev), c28_stratum_oracle(fm$model, fm$data$stratum, pd),
               tolerance = 1e-12)
})

test_that("a weighted fit of one keeps its residuals and gains its coding", {
  d <- c28_data()
  g <- stats::glm(y3 ~ x, data = d, family = stats::binomial, weights = w)
  pd <- c28_quiet(plot_prediction_deviation_panels(g))[[2]]$data
  expect_equal(as.numeric(pd$abs_res_dev), c28_row_oracle(g, pd), tolerance = 1e-12)
  expect_identical(pd$obs_outcome_01, as.integer(d$y3 != "lo")[pd$id])
  ps <- c28_quiet(plot_prediction_deviation_panels(g, data = d))[[2]]$data
  expect_identical(ps$wrong[order(ps$id)], pd$wrong[order(pd$id)])
})

test_that("two-valued outcomes are coded as before (unchanged)", {
  d <- c28_data()
  # A factor declared "yes" before "no": glm takes "yes" as the failure.
  g <- stats::glm(yr ~ x, data = d, family = stats::binomial)
  pd <- c28_quiet(plot_prediction_deviation_panels(g))[[2]]$data
  expect_identical(pd$obs_outcome_01, as.integer(d$yr != "yes")[pd$id])
  expect_equal(as.numeric(pd$abs_res_dev), c28_row_oracle(g, pd), tolerance = 1e-12)
  g01 <- stats::glm(y ~ x, data = d, family = stats::binomial)
  p01 <- c28_quiet(plot_prediction_deviation_panels(g01))[[2]]$data
  expect_identical(p01$obs_outcome_01, as.integer(d$y)[p01$id])
  expect_equal(as.numeric(p01$abs_res_dev), c28_row_oracle(g01, p01), tolerance = 1e-12)
})
