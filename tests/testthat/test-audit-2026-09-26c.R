# Audit 2026-09-26 -- the prediction-deviation panel must not ask WeMix's predict()
# to name a transformed response.
#
# FINDING (external audit). plot_prediction_deviation_panels(fit, type = "auto") on a
# wemix fit errored with
#
#   'language' object cannot be coerced to type 'symbol'
#
# and under plot(type = "all") maihda_try_optional() turned that into "Plot panel
# 'prediction_deviation' could not be computed and was omitted: ..." and dropped the
# panel without a value ever being drawn.
#
# ROOT CAUSE, upstream in WeMix (WeMix:::predict.WeMixResults):
#
#   form     <- as.formula(formula(getCall(object)))
#   response <- as.name(form[[2]])        # <- form[[2]] is log(y), a call
#   if (!deparse(response) %in% colnames(newdata)) newdata[[response]] <- 0
#
# as.name() takes a string or a symbol. Handed a CALL it raises exactly the error
# above, so every wemix fit whose response is an expression -- log(y), I(y > 0),
# sqrt(y) -- died there, while the SAME fit spelled with a pre-computed column
# (ly ~ ...) sailed through, because a bare response is already a symbol. One fit, two
# spellings, two answers; the asymmetry measured on the pre-fix tree was
#
#   engine  response     prediction_deviation panel
#   wemix   ly (bare)    built
#   wemix   log(y)       ERROR (and silently dropped by plot(type = "all"))
#   lme4    ly / log(y)  built
#
# lme4 was never exposed: maihda_prediction_panel_fitted() has a merMod branch that
# predicts the fit's OWN rows WITHOUT newdata, so predict.merMod is the only method
# reached. wemix had no such branch and fell through to the generic
# predict(model, newdata = data) below it.
#
# FIX. A WeMixResults branch that rebuilds the linear predictor from the wrapper via
# maihda_wemix_linpred() -- the route predict_maihda() already takes for this engine,
# which also re-uses the FITTED transformation basis and factor coding -- with
# predict() (no newdata, hence no as.name()) left for a bare fit's own rows, and a
# directed error for the one case that has neither: a bare WeMix fit, prediction data,
# and a transformed response.
#
# CONFINED: every case that worked before is bit-identical after (wemix bare Gaussian
# per-stratum sum 4.8942437404, wemix binomial on the response scale 2.974711829268,
# both bare-fit routes, and every lme4 route), measured on a clean 1f58723 worktree
# and on the fixed tree.
#
# NOT this finding (measured, left alone): maihda_prediction_panel_auto_type() reads
# maihda_family(), which is NULL for a WeMixResults, so type = "auto" routes a wemix
# BINOMIAL fit through the Gaussian branch and plots stratum log-odds under "Fitted
# Value" labels (range [-0.415, 0.788]) where type = "binomial" gives probabilities
# ([0.405, 0.669], against lme4's [0.407, 0.620]). Different root, unchanged here --
# the last test below pins that routing so a fix to it is a deliberate change.

test_that("the panel builds for a wemix fit with a transformed response (2026-09-26)", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  set.seed(1)
  n <- 600
  d <- data.frame(
    g1 = factor(sample(c("a", "b"), n, TRUE)),
    g2 = factor(sample(c("p", "q", "r"), n, TRUE)),
    x = rnorm(n)
  )
  d <- make_strata(d, c("g1", "g2"))$data
  st <- as.integer(factor(d$stratum))
  d$y <- exp(0.4 * d$x + rnorm(length(unique(st)), 0, 0.5)[st] + rnorm(n, 0, 0.4) + 1)
  d$w <- runif(n, 0.5, 2)
  d$ly <- log(d$y)

  f_tx <- suppressMessages(
    fit_maihda(log(y) ~ x + (1 | stratum), data = d, sampling_weights = "w"))
  f_bare <- suppressMessages(
    fit_maihda(ly ~ x + (1 | stratum), data = d, sampling_weights = "w"))

  panel_fitted <- function(p) as.numeric(p[[2]]$data$fitted)

  # The crash itself: an error, not a plot, before the fix.
  p_tx <- plot_prediction_deviation_panels(f_tx, type = "auto")
  expect_s3_class(p_tx, "patchwork")

  # The two spellings of one fit now give ONE answer, not a plot and an error.
  expect_equal(panel_fitted(p_tx),
               panel_fitted(plot_prediction_deviation_panels(f_bare, type = "auto")))

  # Supplied prediction data takes the same route and agrees with the fitted rows.
  expect_equal(panel_fitted(plot_prediction_deviation_panels(f_tx, data = f_tx$data,
                                                            type = "auto")),
               panel_fitted(p_tx))

  # The values are the wrapper's own linear predictor, aggregated per stratum with the
  # prediction weights -- i.e. what predict_maihda() reports for this engine, not
  # something the panel invented.
  eta <- as.numeric(maihda_wemix_linpred(f_tx, include_re = TRUE))
  w <- maihda_prediction_weights(f_tx)
  by_stratum <- vapply(split(seq_len(nrow(f_tx$data)), f_tx$data$stratum),
                       function(i) sum(eta[i] * w[i]) / sum(w[i]), numeric(1))
  expect_equal(sort(panel_fitted(p_tx)), sort(unname(by_stratum)))

  # plot(type = "all") no longer reports the panel as uncomputable and drops it.
  msgs <- character()
  withCallingHandlers(
    invisible(plot(f_tx, type = "all")),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(grep("prediction_deviation", msgs), 0L)

  # A bare WeMix fit has no wrapper to rebuild the design from, and model.frame() of
  # one does not evaluate either, so this one case stays an error -- but a directed
  # one, naming the response and the way out, not WeMix's as.name() message.
  expect_error(
    plot_prediction_deviation_panels(f_tx$model, data = f_tx$data, type = "auto"),
    "response is an expression \\(log\\(y\\)\\)")
  expect_error(
    plot_prediction_deviation_panels(f_tx$model, data = f_tx$data, type = "auto"),
    "fit_maihda\\(\\) model object")

  # A bare fit with a bare response keeps its old predict.WeMixResults() route and its
  # old values, measured on the pre-fix tree. They differ from the wrapper route above
  # because a bare fit exposes no prior weights (maihda_fit_rows_weights() comes back
  # empty for a WeMixResults), so its per-stratum means are unweighted -- pre-existing,
  # and identical on both trees.
  expect_equal(
    sort(panel_fitted(plot_prediction_deviation_panels(f_bare$model, data = f_bare$data,
                                                      type = "auto"))),
    c(0.494053053453, 0.560482248280, 0.724870746282,
      0.751920815704, 1.089107271270, 1.209119769493),
    tolerance = 1e-8)

  # Routing the panel through maihda_wemix_linpred() exposed a latent zero-length path
  # in it: `newdata` without the grouping column made re[as.character(NULL)] a
  # zero-length vector, and `eta + u` then SILENTLY shrank n predictions to none, which
  # reached the caller as dplyr's "`fitted` must be size 600 or 1, not 0". The helper
  # now names the missing column. predict_maihda() never hit this -- it rebuilds
  # `stratum` from the dimension columns first, which is why the path stayed latent.
  nd <- f_bare$data
  nd$stratum <- NULL
  expect_error(maihda_wemix_linpred(f_bare, newdata = nd, include_re = TRUE),
               "must carry the 'stratum' column")
  expect_error(plot_prediction_deviation_panels(f_bare, data = nd, type = "auto"),
               "must carry the 'stratum' column")
  # include_re = FALSE needs no grouping column and is unaffected.
  expect_length(maihda_wemix_linpred(f_bare, newdata = nd, include_re = FALSE), nrow(nd))
  # predict_maihda() still predicts every row of such a frame, as it did before.
  expect_length(predict_maihda(f_bare, newdata = nd), nrow(nd))
})

test_that("a wemix binomial panel keeps its response-scale values (2026-09-26)", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  set.seed(1)
  n <- 600
  d <- data.frame(
    g1 = factor(sample(c("a", "b"), n, TRUE)),
    g2 = factor(sample(c("p", "q", "r"), n, TRUE)),
    x = rnorm(n)
  )
  d <- make_strata(d, c("g1", "g2"))$data
  st <- as.integer(factor(d$stratum))
  # Same RNG stream as the Gaussian test above, so the fit below is the one measured.
  d$y <- exp(0.4 * d$x + rnorm(length(unique(st)), 0, 0.5)[st] + rnorm(n, 0, 0.4) + 1)
  d$w <- runif(n, 0.5, 2)
  d$ly <- log(d$y)
  d$b <- rbinom(n, 1, plogis(0.5 * d$x + rnorm(length(unique(st)), 0, 0.4)[st]))

  fb <- suppressMessages(
    fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
               family = binomial()))
  panel_fitted <- function(p) as.numeric(p[[2]]$data$fitted)

  # type = "binomial" summarises on the PROBABILITY scale, and the new wemix branch
  # reproduces the old predict(type = "response") values to the last bit.
  pb <- panel_fitted(plot_prediction_deviation_panels(fb, type = "binomial"))
  expect_true(all(pb >= 0 & pb <= 1))
  expect_equal(sum(pb), 2.974711829268, tolerance = 1e-9)
  expect_equal(panel_fitted(plot_prediction_deviation_panels(fb, data = fb$data,
                                                            type = "binomial")),
               pb)

  # Those probabilities are plogis() of the wrapper's linear predictor, per stratum.
  eta <- as.numeric(maihda_wemix_linpred(fb, include_re = TRUE))
  w <- maihda_prediction_weights(fb)
  by_stratum <- vapply(split(seq_len(nrow(fb$data)), fb$data$stratum),
                       function(i) sum(stats::plogis(eta[i]) * w[i]) / sum(w[i]),
                       numeric(1))
  expect_equal(sort(pb), sort(unname(by_stratum)))

  # PRE-EXISTING and deliberately unchanged: maihda_family() is NULL for a
  # WeMixResults, so "auto" still calls this binomial fit Gaussian and plots log-odds.
  # If a later fix routes it to "binomial", THESE are the expectations to update.
  expect_identical(maihda_prediction_panel_auto_type(fb$model), "gaussian")
  expect_true(is.null(maihda_family(fb$model)))
  expect_identical(maihda_model_family(fb)$family, "binomial")
  expect_true(any(panel_fitted(plot_prediction_deviation_panels(fb, type = "auto")) < 0))
})
