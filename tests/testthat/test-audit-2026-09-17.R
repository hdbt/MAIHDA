# Audit 2026-09-17: the binomial prediction-deviation panel on an aggregated binomial
# fit -- CONFIRMED (a follow-up the 2026-09-15 surprise-panel pass spawned), and
# wider.
#
# brms's fitted() and posterior_epred() return the expected success COUNT
# (trials * p) for a `y | trials(n)` fit. predict_maihda() divides by the trials
# (maihda_brms_response_to_prob()); plot_prediction_deviation_panels() did not, so on
# a 12-stratum fit with 3 to 30 trials per row it drew stratum "probabilities" of 2.9
# to 16.9 (true 0.14 to 0.80) with a dashed mean of 8.17, and clamped every stratum
# interval to [1, 1]. The same fit's strata all had a |deviance residual| of 0, the
# aggregated response being treated as unobserved, so the five labelled "worst-fit"
# strata were arbitrary. Wider: the lme4 twin, a fit_maihda() cbind(successes,
# failures) fit, could not be drawn at all -- its stored data carries the response as
# a two-column matrix under the name "cbind(s, f)", which the panel took for the 0/1
# outcome ("`obs_outcome` must be size 240 or 1, not 480") -- and plot(type = "all")
# warned and dropped the panel. It worked only when the original data was passed.
#
# FIX (the residual and lme4 choices per the user): the panel reads the fit's trial
# counts off its trials() term on the plotted rows and divides the row estimates,
# their interval half-widths and the per-draw expected values by them; a bare brmsfit,
# and prediction data the stored weights do not match, are weighted by those counts
# as the maihda_model's own rows are; the strata are ranked by the
# binomial deviance residual of successes out of trials at the posterior-mean
# probability, which is lme4's residuals(type = "deviance") for a cbind() fit (sign
# dropped); and a matrix response column counts as unobserved, so a cbind() fit
# falls back to its model residuals as it did with the original data. Bernoulli
# fits, brms and lme4, are untouched.

d17_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression and the warnings it raised, muffled.
d17_capture <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

# 12 strata x 20 rows of successes out of 3 to 30 trials.
d17_data <- function() {
  set.seed(917)
  d <- data.frame(stratum = factor(rep(seq_len(12), each = 20)))
  d$x <- stats::rnorm(nrow(d))
  d$n <- sample(3:30, nrow(d), replace = TRUE)
  p <- stats::plogis(-0.4 + 0.5 * d$x + stats::rnorm(12, sd = 0.8)[d$stratum])
  d$s <- stats::rbinom(nrow(d), d$n, p)
  d$f <- d$n - d$s
  d
}

# |binomial deviance residual| of s successes out of n trials at probability p,
# written out independently of the package.
d17_dev_resid <- function(s, n, p) {
  f <- n - s
  a <- ifelse(s > 0, s * log(s / (n * p)), 0)
  b <- ifelse(f > 0, f * log(f / (n * (1 - p))), 0)
  sqrt(pmax(2 * (a + b), 0))
}

# Weighted mean of x within each stratum, named by stratum.
d17_by_stratum <- function(x, stratum, w) {
  key <- as.character(stratum)
  vapply(split(seq_along(key), key), function(i) sum(x[i] * w[i]) / sum(w[i]), numeric(1))
}

test_that("the aggregated binomial deviance residual is lme4's cbind() residual", {
  skip_on_cran()
  d <- d17_data()
  m <- d17_quiet(fit_maihda(cbind(s, f) ~ x + (1 | stratum), data = d, family = "binomial"))
  p <- as.numeric(stats::predict(m$model, type = "response"))
  expect_equal(maihda_binomial_abs_deviance_residual_agg(d$s, d$n, p),
               abs(as.numeric(stats::residuals(m$model, type = "deviance"))),
               tolerance = 1e-10)
  expect_equal(maihda_binomial_abs_deviance_residual_agg(d$s, d$n, p),
               d17_dev_resid(d$s, d$n, p), tolerance = 1e-12)
})

test_that("the aggregated residual is the 0/1 residual at one trial and NA where unscorable", {
  p <- c(0.2, 0.7, 0.4, 0.5, 0.5, 0.5, 0.3)
  s <- c(0, 1, 1, 0, 3, NA, 4)
  n <- c(1, 1, 1, 1, 0, 5, 3)
  out <- maihda_binomial_abs_deviance_residual_agg(s, n, p)
  expect_equal(out[1:4], maihda_binomial_abs_deviance_residual(c(0L, 1L, 1L, 0L), p[1:4]))
  # No trials, a missing count, and more successes than trials cannot be scored.
  expect_identical(is.na(out), c(rep(FALSE, 4), TRUE, TRUE, TRUE))
  # All successes or none: finite, the other term counting as 0.
  expect_equal(maihda_binomial_abs_deviance_residual_agg(c(0, 5), c(5, 5), c(0.2, 0.9)),
               c(sqrt(-2 * 5 * log(0.8)), sqrt(-2 * 5 * log(0.9))))
  expect_identical(maihda_prediction_panel_per_trial(c(6, 6, 6, 6), c(3, 0, NA, -1)),
                   c(2, NA, NA, NA))
})

test_that("an lme4 cbind() fit draws its panel on its own stored data", {
  skip_on_cran()
  d <- d17_data()
  m <- d17_quiet(fit_maihda(cbind(s, f) ~ x + (1 | stratum), data = d, family = "binomial"))
  # The stored data carries the response as the two-column matrix the panel tripped on.
  expect_true(is.matrix(m$data[["cbind(s, f)"]]))

  r <- d17_capture(plot_prediction_deviation_panels(m))
  expect_s3_class(r$value, "patchwork")
  expect_identical(grep("plot_prediction_deviation_panels", r$warnings, value = TRUE),
                   character(0))
  pd <- r$value[[2]]$data
  w <- as.numeric(maihda_prediction_weights(m))
  expect_equal(w, d$n)
  fit_ref <- d17_by_stratum(as.numeric(stats::predict(m$model, type = "response")), d$stratum, w)
  res_ref <- d17_by_stratum(abs(as.numeric(stats::residuals(m$model, type = "deviance"))),
                            d$stratum, w)
  key <- as.character(pd$stratum)
  expect_equal(pd$fitted, as.numeric(fit_ref[key]), tolerance = 1e-10)
  expect_equal(pd$abs_res_dev, as.numeric(res_ref[key]), tolerance = 1e-10)
  expect_true(all(pd$abs_res_dev > 0))

  # The same numbers as the panel drawn from the original data, which always worked.
  po <- d17_quiet(plot_prediction_deviation_panels(m, data = d))[[2]]$data
  for (col in c("stratum", "fitted", "abs_res_dev", "rank", "direction")) {
    expect_equal(pd[[col]], po[[col]], tolerance = 1e-12, label = col)
  }

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ra <- d17_capture(plots <- plot(m, type = "all"))
  expect_s3_class(plots$prediction_deviation, "patchwork")
  expect_identical(grep("prediction_deviation", ra$warnings, value = TRUE), character(0))
})

test_that("the trial counts are read off a brms trials() term, alone or beside weights()", {
  skip_if_not_installed("brms")
  d <- data.frame(s = c(1, 2, 3), n = c(4, 5, 6), w = c(1, 2, 1))
  mock <- function(f) structure(list(formula = brms::bf(f)), class = "brmsfit")
  expect_identical(maihda_prediction_panel_brms_trials(mock(s | trials(n) ~ 1), d), c(4, 5, 6))
  expect_identical(maihda_prediction_panel_brms_trials(mock(s | trials(20) ~ 1), d), rep(20, 3))
  expect_identical(
    maihda_prediction_panel_brms_trials(mock(s | trials(n) + weights(w) ~ 1), d), c(4, 5, 6))
  expect_identical(
    maihda_prediction_panel_brms_trials(mock(s | weights(w) + trials(n) ~ 1), d), c(4, 5, 6))
  expect_identical(maihda_prediction_panel_brms_successes(mock(s | trials(n) ~ 1), d), c(1, 2, 3))
  # No trials() term, or not a brms fit: no trial counts, so nothing is divided.
  expect_null(maihda_prediction_panel_brms_trials(mock(s | weights(w) ~ 1), d))
  expect_null(maihda_prediction_panel_brms_trials(mock(s ~ 1), d))
  expect_null(maihda_prediction_panel_brms_trials(stats::lm(s ~ n, data = d), d))
  # A trials() term the plotted rows cannot supply is an error, never a count axis.
  expect_error(maihda_prediction_panel_brms_trials(mock(s | trials(n) ~ 1), d[c("s", "w")]),
               "Could not evaluate the trials\\(\\) term")
  expect_null(maihda_prediction_panel_brms_successes(mock(s | trials(n) ~ 1), d[c("n", "w")]))
})

# A brms `y | trials(n)` stand-in: stats::fitted() and brms::posterior_epred() dispatch
# to methods registered on a private subclass and return what brms returns for such a
# fit -- success COUNTS, the per-row probability `p` (and its `draws`) times the trials
# `mult(newdata)` -- indexed by the .row column. Its predict() fails, as it must not be
# what the panel reads.
d17_mock_brms <- function(formula, data, p, draws, mult) {
  registerS3method("fitted", "d17_mock_brmsfit", function(object, newdata = NULL, ...) {
    rows <- newdata$.row
    k <- object$mult(newdata)
    q <- apply(object$draws[, rows, drop = FALSE], 2, stats::quantile, c(0.025, 0.975))
    cbind(Estimate = object$p[rows] * k, Est.Error = 0.1 * k,
          Q2.5 = q[1, ] * k, Q97.5 = q[2, ] * k)
  }, envir = asNamespace("stats"))
  registerS3method("posterior_epred", "d17_mock_brmsfit", function(object, newdata = NULL, ...) {
    rows <- newdata$.row
    sweep(object$draws[, rows, drop = FALSE], 2, object$mult(newdata), "*")
  }, envir = asNamespace("brms"))
  data$.row <- seq_len(nrow(data))
  model <- structure(
    list(formula = brms::bf(formula), data = data, p = p, draws = draws, mult = mult,
         family = structure(list(family = "binomial", link = "logit"), class = "family")),
    class = c("d17_mock_brmsfit", "brmsfit"))
  structure(
    list(model = model, engine = "brms", formula = formula, data = data,
         original_data = data, family = list(family = "binomial", link = "logit"),
         strata_vars = "stratum", strata_info = NULL, context_vars = NULL,
         sampling_weights = NULL, longitudinal_info = NULL, response_recoding = NULL),
    class = "maihda_model")
}

# Draws of each row's probability around `p`, and the stratum references the panel
# must reproduce: trial-weighted means of p, of the per-draw probabilities (then 2.5%
# and 97.5% quantiles), and of the deviance residual at p.
d17_mock_case <- function(d, trials) {
  set.seed(4)
  p <- stats::plogis(-0.4 + 0.5 * d$x + stats::rnorm(12, sd = 0.8)[d$stratum])
  draws <- t(vapply(seq_len(200), function(i) {
    stats::plogis(stats::qlogis(p) + stats::rnorm(1, sd = 0.2) + stats::rnorm(length(p), sd = 0.1))
  }, numeric(length(p))))
  key <- as.character(d$stratum)
  per_draw <- vapply(split(seq_along(key), key), function(i)
    as.vector(draws[, i, drop = FALSE] %*% (trials[i] / sum(trials[i]))), numeric(nrow(draws)))
  list(p = p, draws = draws,
       fitted = d17_by_stratum(p, d$stratum, trials),
       lower = apply(per_draw, 2, stats::quantile, 0.025, names = FALSE),
       upper = apply(per_draw, 2, stats::quantile, 0.975, names = FALSE),
       resid = d17_by_stratum(d17_dev_resid(d$s, trials, p), d$stratum, trials))
}

d17_check_panel <- function(pd, ref) {
  key <- as.character(pd$stratum)
  expect_equal(pd$fitted, as.numeric(ref$fitted[key]), tolerance = 1e-12)
  expect_equal(pd$ci_lower, as.numeric(ref$lower[key]), tolerance = 1e-12)
  expect_equal(pd$ci_upper, as.numeric(ref$upper[key]), tolerance = 1e-12)
  expect_equal(pd$abs_res_dev, as.numeric(ref$resid[key]), tolerance = 1e-12)
  expect_true(all(pd$ci_upper > pd$ci_lower))
}

test_that("a brms trials() panel is drawn on the probability scale (Stan-free)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- d17_data()
  ref <- d17_mock_case(d, d$n)
  m <- d17_mock_brms(s | trials(n) ~ x + (1 | stratum), d, ref$p, ref$draws,
                     function(nd) nd$n)
  p <- d17_quiet(plot_prediction_deviation_panels(m, type = "binomial"))
  pd <- p[[2]]$data
  d17_check_panel(pd, ref)
  expect_true(all(p[[1]]$data$fitted > 0 & p[[1]]$data$fitted < 1))
  # The five labelled strata are the five worst-fitting ones.
  lab <- Filter(function(l) inherits(l$geom, "GeomLabelRepel"), p[[2]]$layers)[[1]]$data
  expect_identical(sort(as.character(lab$stratum)),
                   sort(names(sort(ref$resid, decreasing = TRUE))[1:5]))

  # A bare brmsfit is weighted by its trials too, and gives the same panel.
  pb <- d17_quiet(plot_prediction_deviation_panels(m$model, type = "binomial"))[[2]]$data
  d17_check_panel(pb, ref)

  # So is other prediction data, whose rows the fit's stored weights do not match.
  keep <- which(d$x > 0)
  ds <- d[keep, ]
  per_draw <- vapply(split(seq_along(keep), as.character(ds$stratum)), function(i)
    as.vector(ref$draws[, keep[i], drop = FALSE] %*% (ds$n[i] / sum(ds$n[i]))),
    numeric(nrow(ref$draws)))
  sub_ref <- list(fitted = d17_by_stratum(ref$p[keep], ds$stratum, ds$n),
                  lower = apply(per_draw, 2, stats::quantile, 0.025, names = FALSE),
                  upper = apply(per_draw, 2, stats::quantile, 0.975, names = FALSE),
                  resid = d17_by_stratum(d17_dev_resid(ds$s, ds$n, ref$p[keep]),
                                         ds$stratum, ds$n))
  expect_length(unique(ds$stratum), 12L)
  ps <- d17_quiet(plot_prediction_deviation_panels(m, data = m$data[keep, ],
                                                   type = "binomial"))[[2]]$data
  d17_check_panel(ps, sub_ref)
})

test_that("a constant trials() fit is converted the same way (Stan-free)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- d17_data()
  d$s <- pmin(d$s, 20)
  ref <- d17_mock_case(d, rep(20, nrow(d)))
  m <- d17_mock_brms(s | trials(20) ~ x + (1 | stratum), d, ref$p, ref$draws,
                     function(nd) rep(20, nrow(nd)))
  pd <- d17_quiet(plot_prediction_deviation_panels(m, type = "binomial"))[[2]]$data
  d17_check_panel(pd, ref)
})

test_that("brms trials() panels score the successes out of trials on a real fit", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). The references are
  # computed from the same fit's own draws, so sampler quality is irrelevant.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  d <- d17_data()
  m <- d17_quiet(fit_maihda(s | trials(n) ~ x + (1 | stratum), data = d, engine = "brms",
                            family = "binomial", chains = 1, iter = 600, warmup = 300,
                            refresh = 0, seed = 1))
  n <- m$data$n
  w <- as.numeric(maihda_prediction_weights(m))
  expect_equal(w, n)
  p_hat <- as.numeric(stats::fitted(m$model, newdata = m$data, summary = TRUE)[, "Estimate"]) / n
  ep <- brms::posterior_epred(m$model, newdata = m$data)
  key <- as.character(m$data$stratum)
  per_draw <- vapply(split(seq_along(key), key), function(i)
    as.vector((ep[, i, drop = FALSE] / rep(n[i], each = nrow(ep))) %*% (w[i] / sum(w[i]))),
    numeric(nrow(ep)))
  ref <- list(fitted = d17_by_stratum(p_hat, key, w),
              lower = apply(per_draw, 2, stats::quantile, 0.025, names = FALSE),
              upper = apply(per_draw, 2, stats::quantile, 0.975, names = FALSE),
              resid = d17_by_stratum(d17_dev_resid(m$data$s, n, p_hat), key, w))
  names(ref$lower) <- names(ref$upper) <- colnames(per_draw)

  r <- d17_capture(plot_prediction_deviation_panels(m))
  expect_identical(grep("plot_prediction_deviation_panels", r$warnings, value = TRUE),
                   character(0))
  pd <- r$value[[2]]$data
  d17_check_panel(pd, ref)
  # Row by row, the panel's probabilities are predict_maihda()'s.
  rows <- maihda_prediction_panel_fitted(m$model, m$data, "binomial", trials = n)$fit
  expect_equal(rows, as.numeric(predict_maihda(m, type = "individual", scale = "response")),
               tolerance = 1e-12)

  # The lme4 cbind() twin of the same data lands on the same scale.
  ml <- d17_quiet(fit_maihda(cbind(s, f) ~ x + (1 | stratum), data = d, family = "binomial"))
  pl <- d17_quiet(plot_prediction_deviation_panels(ml))[[2]]$data
  expect_lt(max(abs(pl$fitted - pd$fitted[match(pl$stratum, pd$stratum)])), 0.02)
  expect_lt(max(abs(pl$abs_res_dev - pd$abs_res_dev[match(pl$stratum, pd$stratum)])), 0.05)
})
