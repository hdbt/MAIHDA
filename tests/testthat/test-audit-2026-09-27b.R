# Audit 2026-09-27b: the prediction-deviation panel on rows other than the fitted
# ones -- external findings F05 and F04, both CONFIRMED and wider.
#
# F05. The binomial panel's deviance residuals were wrong on two branches of one
# helper. (1) When the outcome of `data` could not be read as one 0/1 value per row --
# a cbind(successes, failures) response, a proportion -- it took the FITTED rows'
# residuals whenever `data` had as many rows, without aligning them: on this file's
# fixture, reordering the fitted rows put a cbind() fit's residuals up to 2.83 off on a
# row and 0.53 on a stratum mean, and a subset of them got a residual of 0 everywhere.
# (2) A 0/1 outcome with prior weights was scored as a single trial whatever its
# weight, although R's binomial family reads a weight of 25 as 25 trials: every
# residual was exactly five times too small, on the DEFAULT panel too, and the same
# model spelled cbind() gave another panel. Wider: a missing outcome in `data` was
# scored 0, a perfect fit, so plotted on the original data a stratum's mean residual
# shrank by exactly its share of missing outcomes; and glm() fits were affected alike.
#
# F04. The stratum summaries' weights were the fitted rows' weights, paired with the
# plotted rows by position whenever the counts matched: on this file's fixture,
# reordering the fitted rows moved a precision-weighted Gaussian stratum mean by 2.57
# and a cbind() stratum probability by 0.082, a subset was silently unweighted, and new
# rows of the same count took the training rows' weights instead of their own. Found
# while checking the fix: a weighted MASS::polr or bare ordinal::clmm, whose weights
# stats::weights() does not return, was unweighted on every route, and reading the
# supplied rows' weights alone would have made its two routes disagree.
#
# FIX (the choices are the user's): a row's residual is R's binomial deviance residual,
# every prior weight counting as that many trials -- the engine's own
# residuals(type = "deviance") on the fitted rows, and the same quantity rebuilt from
# the plotted rows' own successes and trials elsewhere; a row that cannot be scored is
# NA, left out of the stratum means and the labels (drawn hollow), with a warning when
# trial counts cannot be recovered or no row can be scored. Each plotted row's weight is
# its own -- a column of `data`, a bare fit's `weights =` expression, or the fitted row
# it is shown to be by its row name AND its prediction -- and when any row's weight
# cannot be found every row is weighted equally, with a warning; a polr / clm / clmm fit
# whose weights() has no method is weighted by its frame's "(weights)" on its own rows
# too (the user's choice). The fitted-rows panel of every other fit is unchanged except
# the 0/1-with-weights residual: 34 panels were identical() to HEAD's.

a27b_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression and the panel's own warnings (lme4's se.fit notice and
# ggplot2's are not this panel's business).
a27b_capture <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value,
       warnings = grep("plot_prediction_deviation_panels", w, value = TRUE))
}

# The stratum-level panel data, one row per stratum in stratum order.
a27b_panel <- function(p) {
  x <- p[[2]]$data
  x[order(x$stratum), , drop = FALSE]
}

# |binomial deviance residual| of s successes out of n trials at probability p,
# written out independently of the package.
a27b_dev <- function(s, n, p) {
  f <- n - s
  a <- ifelse(s > 0, s * log(s / (n * p)), 0)
  b <- ifelse(f > 0, f * log(f / (n * (1 - p))), 0)
  sqrt(pmax(2 * (a + b), 0))
}

# Weighted mean of x by stratum over the rows where x is finite, in stratum order.
a27b_wmean <- function(x, g, w) {
  ok <- is.finite(x)
  out <- tapply(seq_along(x)[ok], g[ok], function(i) sum(x[i] * w[i]) / sum(w[i]))
  as.numeric(out[order(as.numeric(names(out)))])
}

# 12 strata x 15 rows: successes out of 3 to 40 trials, 25-trial proportions, a 0/1
# outcome, and a Gaussian outcome with strongly varying precision weights.
a27b_data <- function() {
  set.seed(20260927)
  d <- expand.grid(g1 = paste0("a", 1:3), g2 = paste0("b", 1:4), rep = 1:15,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  st <- as.integer(interaction(d$g1, d$g2, drop = TRUE))
  d$stratum <- st
  d$x <- stats::rnorm(nrow(d), mean = st / 4, sd = 1.5)
  u <- stats::rnorm(12, 0, 0.8)
  p <- stats::plogis(-0.3 + 0.6 * d$x - st / 8 + u[st])
  d$ntr <- sample(3:40, nrow(d), replace = TRUE)
  d$s <- stats::rbinom(nrow(d), d$ntr, p)
  d$f <- d$ntr - d$s
  d$n25 <- 25
  d$p25 <- stats::rbinom(nrow(d), 25, p) / 25
  d$y01 <- stats::rbinom(nrow(d), 1, p)
  d$wg <- round(stats::runif(nrow(d), 0.2, 20), 2)
  d$yg <- 10 + 3 * d$x + u[st] * 2 + stats::rnorm(nrow(d), 0, 1 / sqrt(d$wg))
  d
}

a27b_perm <- function(d) {
  set.seed(4)
  d[sample(nrow(d)), , drop = FALSE]
}

test_that("a cbind() fit is scored on the supplied rows' own successes and trials", {
  skip_on_cran()
  d <- a27b_data()
  dr <- a27b_perm(d)
  m <- a27b_quiet(fit_maihda(cbind(s, f) ~ x + (1 | stratum), data = d,
                             family = "binomial"))

  # Control (unchanged): the fitted-rows panel is the engine's residuals, trial-weighted.
  pd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(m)))
  expect_equal(pd$abs_res_dev,
               a27b_wmean(abs(as.numeric(stats::residuals(m$model, type = "deviance"))),
                          d$stratum, d$ntr), tolerance = 1e-10)

  # Reordering the fitted rows changes nothing -- the residuals (F05) and the trial
  # weights of the stratum probabilities (F04) -- and there is nothing to warn about.
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = dr))
  expect_identical(r$warnings, character(0))
  pr <- a27b_panel(r$value)
  expect_equal(pr$abs_res_dev, pd$abs_res_dev, tolerance = 1e-9)
  expect_equal(pr$fitted, pd$fitted, tolerance = 1e-12)

  # The bare glmer, handed only the columns it needs of the same reordered rows, agrees
  # (row by row: see the glm test below).
  pb <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(
    m$model, data = dr[, c("x", "s", "f", "stratum")], type = "binomial")))
  expect_equal(pb$abs_res_dev, pd$abs_res_dev, tolerance = 1e-9)
  expect_equal(pb$fitted, pd$fitted, tolerance = 1e-12)

  # A subset has fewer rows: scored and weighted on those rows, not set to 0.
  keep <- which(d$x > stats::median(d$x))
  ds <- d[keep, ]
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = ds))
  expect_identical(r$warnings, character(0))
  psub <- a27b_panel(r$value)
  p_s <- as.numeric(stats::predict(m$model, newdata = ds, type = "response"))
  expect_equal(psub$abs_res_dev,
               a27b_wmean(a27b_dev(ds$s, ds$ntr, p_s), ds$stratum, ds$ntr)[
                 match(psub$stratum, sort(unique(ds$stratum)))],
               tolerance = 1e-10)
  expect_equal(psub$fitted,
               a27b_wmean(p_s, ds$stratum, ds$ntr)[
                 match(psub$stratum, sort(unique(ds$stratum)))],
               tolerance = 1e-12)
  expect_true(all(psub$abs_res_dev > 0))
})

test_that("a glm() row's residual is its own successes out of its own trials", {
  skip_on_cran()
  d <- a27b_data()
  dr <- a27b_perm(d)[, c("x", "s", "f", "ntr")]
  g <- stats::glm(cbind(s, f) ~ x, data = d, family = stats::binomial())
  pc <- a27b_quiet(plot_prediction_deviation_panels(g, data = dr))[[2]]$data
  pc <- pc[order(pc$id), , drop = FALSE]
  p <- as.numeric(stats::predict(g, newdata = dr, type = "response"))
  expect_equal(as.numeric(pc$abs_res_dev), a27b_dev(dr$s, dr$ntr, p), tolerance = 1e-10)
  # ... which on the fitted rows is the engine's own residual.
  p0 <- a27b_quiet(plot_prediction_deviation_panels(g))[[2]]$data
  p0 <- p0[order(p0$id), , drop = FALSE]
  expect_equal(as.numeric(p0$abs_res_dev),
               abs(as.numeric(stats::residuals(g, type = "deviance"))), tolerance = 1e-12)
})

test_that("a weighted proportion is scored with its trials on reordered rows", {
  skip_on_cran()
  d <- a27b_data()
  dr <- a27b_perm(d)
  # The wrapper stores its weights as values: the reordered rows are recognised as
  # fitted rows by their row names and predictions.
  m <- a27b_quiet(fit_maihda(p25 ~ x + (1 | stratum), data = d, family = "binomial",
                             weights = n25))
  pd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(m)))
  expect_equal(pd$abs_res_dev,
               a27b_wmean(abs(as.numeric(stats::residuals(m$model, type = "deviance"))),
                          d$stratum, d$n25), tolerance = 1e-10)
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = dr))
  expect_identical(r$warnings, character(0))
  expect_equal(a27b_panel(r$value)$abs_res_dev, pd$abs_res_dev, tolerance = 1e-9)

  # A bare glmer names its weights column, which is read off the supplied rows.
  g <- a27b_quiet(lme4::glmer(p25 ~ x + (1 | stratum), data = d,
                              family = stats::binomial(), weights = n25))
  gd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(g)))
  r <- a27b_capture(plot_prediction_deviation_panels(g, data = dr))
  expect_identical(r$warnings, character(0))
  expect_equal(a27b_panel(r$value)$abs_res_dev, gd$abs_res_dev, tolerance = 1e-9)
})

test_that("a 0/1 outcome with prior weights counts each weight as trials", {
  skip_on_cran()
  d <- a27b_data()
  m <- a27b_quiet(fit_maihda(y01 ~ x + (1 | stratum), data = d, family = "binomial",
                             weights = n25))
  pd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(m)))
  back <- abs(as.numeric(stats::residuals(m$model, type = "deviance")))
  expect_equal(pd$abs_res_dev, a27b_wmean(back, d$stratum, d$n25), tolerance = 1e-10)

  # Five times the one-trial residual at 25 trials a row: the sqrt(trials) factor.
  p <- as.numeric(stats::fitted(m$model))
  one <- maihda_binomial_abs_deviance_residual(d$y01, p)
  expect_equal(back, 5 * one, tolerance = 1e-10)

  # One model, two spellings, one panel.
  d$ys <- 25 * d$y01
  d$yf <- 25 - d$ys
  mc <- a27b_quiet(fit_maihda(cbind(ys, yf) ~ x + (1 | stratum), data = d,
                              family = "binomial"))
  expect_equal(unname(lme4::fixef(mc$model)), unname(lme4::fixef(m$model)),
               tolerance = 1e-8)
  pc <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(mc)))
  expect_equal(pd$abs_res_dev, pc$abs_res_dev, tolerance = 1e-6)
  expect_equal(pd$fitted, pc$fitted, tolerance = 1e-6)

  # And on reordered rows.
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = a27b_perm(d)))
  expect_identical(r$warnings, character(0))
  expect_equal(a27b_panel(r$value)$abs_res_dev, pd$abs_res_dev, tolerance = 1e-9)
})

test_that("a non-integer weight scales the residual as R's binomial family does", {
  skip_on_cran()
  # The user's choice: R's deviance residual for every prior weight, not only whole
  # numbers -- the residual residuals(type = "deviance") reports.
  d <- a27b_data()
  d$wq <- d$wg / 4
  g <- suppressWarnings(stats::glm(y01 ~ x, data = d, family = stats::binomial(),
                                   weights = wq))
  p0 <- a27b_quiet(plot_prediction_deviation_panels(g, data = NULL))[[2]]$data
  p0 <- p0[order(p0$id), , drop = FALSE]
  back <- abs(as.numeric(stats::residuals(g, type = "deviance")))
  expect_equal(as.numeric(p0$abs_res_dev), back, tolerance = 1e-12)
  expect_equal(back, sqrt(d$wq) *
                 maihda_binomial_abs_deviance_residual(d$y01, as.numeric(stats::fitted(g))),
               tolerance = 1e-10)
  # The weights = wq column is read off new rows too.
  dn <- d[rev(seq_len(nrow(d))), c("x", "y01", "wq")]
  pn <- a27b_quiet(plot_prediction_deviation_panels(g, data = dn))[[2]]$data
  pn <- pn[order(pn$id), , drop = FALSE]
  pp <- as.numeric(stats::predict(g, newdata = dn, type = "response"))
  expect_equal(as.numeric(pn$abs_res_dev), a27b_dev(dn$wq * dn$y01, dn$wq, pp),
               tolerance = 1e-10)
})

test_that("a row that cannot be scored is left out, not counted as a perfect fit", {
  expect_identical(maihda_binomial_abs_deviance_residual(c(1L, NA, 0L), c(0.5, 0.5, NA)),
                   c(sqrt(-2 * log(0.5)), NA_real_, NA_real_))

  skip_on_cran()
  d <- a27b_data()
  set.seed(3)
  d$y01[sample(nrow(d), 27)] <- NA
  m <- a27b_quiet(fit_maihda(y01 ~ x + (1 | stratum), data = d, family = "binomial"))
  pd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(m)))
  # The original data, missing outcomes and all: the same stratum residuals as the
  # fitted rows, silently (a missing outcome is not a problem to report).
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = d))
  expect_identical(r$warnings, character(0))
  expect_equal(a27b_panel(r$value)$abs_res_dev, pd$abs_res_dev, tolerance = 1e-12)
})

test_that("rows whose trial counts cannot be found are unscored and warned about", {
  skip_on_cran()
  d <- a27b_data()
  m <- a27b_quiet(fit_maihda(p25 ~ x + (1 | stratum), data = d, family = "binomial",
                             weights = n25))
  # New rows numbered like the fitted ones: the names match, the predictions do not,
  # so the names are refuted and the wrapper's stored weights -- the trial counts --
  # cannot be carried to them.
  set.seed(9)
  dn <- d
  dn$x <- stats::rnorm(nrow(d))
  rownames(dn) <- NULL
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = dn))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "prior weights are not known for 180 of the 180 rows",
               fixed = TRUE)
  expect_match(r$warnings, "deviance residuals are undefined", fixed = TRUE)
  pd <- a27b_panel(r$value)
  expect_true(all(is.nan(pd$abs_res_dev)))
  # No label on an unscored stratum; every stratum still drawn, hollow.
  layers <- r$value[[2]]$layers
  labels <- Filter(function(l) inherits(l$geom, "GeomLabelRepel") ||
                     inherits(l$geom, "GeomLabel"), layers)
  expect_length(labels, 1L)
  expect_identical(nrow(labels[[1]]$data), 0L)
  hollow <- Filter(function(l) inherits(l$geom, "GeomPoint") &&
                     identical(l$aes_params$shape, 1), layers)
  expect_length(hollow, 1L)
  expect_identical(nrow(hollow[[1]]$data), 12L)
})

test_that("data without the outcome keeps its points and says why nothing is scored", {
  skip_on_cran()
  d <- a27b_data()
  m <- a27b_quiet(fit_maihda(cbind(s, f) ~ x + (1 | stratum), data = d,
                             family = "binomial"))
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = d[, c("x", "stratum")]))
  expect_true(any(grepl("does not contain the model's observed outcome", r$warnings,
                        fixed = TRUE)))
  # The trial counts come from the outcome, so the probabilities are unweighted and
  # the warning says that too.
  expect_true(any(grepl("binomial trial counts are not known", r$warnings, fixed = TRUE)))
  pd <- a27b_panel(r$value)
  p <- as.numeric(stats::predict(m$model, newdata = d, type = "response"))
  expect_equal(pd$fitted, a27b_wmean(p, d$stratum, rep(1, nrow(d))), tolerance = 1e-12)
  hollow <- Filter(function(l) inherits(l$geom, "GeomPoint") &&
                     identical(l$aes_params$shape, 1), r$value[[2]]$layers)
  expect_identical(nrow(hollow[[1]]$data), 12L)
})

test_that("precision weights follow their rows: permutations, subsets, new rows", {
  skip_on_cran()
  d <- a27b_data()
  dr <- a27b_perm(d)
  m <- a27b_quiet(fit_maihda(yg ~ x + (1 | stratum), data = d, weights = wg))
  g <- a27b_quiet(lme4::lmer(yg ~ x + (1 | stratum), data = d, weights = wg))
  fv <- as.numeric(stats::fitted(g))
  keep <- which(d$x > stats::median(d$x))
  ds <- d[keep, ]
  for (fit in list(m, g)) {
    pd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(fit, type = "gaussian")))
    expect_equal(pd$fitted, a27b_wmean(fv, d$stratum, d$wg), tolerance = 1e-8)

    r <- a27b_capture(plot_prediction_deviation_panels(fit, data = dr, type = "gaussian"))
    expect_identical(r$warnings, character(0))
    expect_equal(a27b_panel(r$value)$fitted, pd$fitted, tolerance = 1e-10)

    r <- a27b_capture(plot_prediction_deviation_panels(fit, data = ds, type = "gaussian"))
    expect_identical(r$warnings, character(0))
    ps <- a27b_panel(r$value)
    expect_equal(ps$fitted, a27b_wmean(fv[keep], ds$stratum, ds$wg)[
      match(ps$stratum, sort(unique(ds$stratum)))], tolerance = 1e-8)
  }

  # New rows of the same count. A bare lmer names its weights column, so the new rows
  # carry their OWN weights; the wrapper cannot know them and says so.
  set.seed(99)
  dn <- d
  dn$x <- stats::rnorm(nrow(d), mean = d$stratum / 4, sd = 1.5)
  dn$wg <- round(stats::runif(nrow(d), 0.2, 20), 2)
  rownames(dn) <- NULL
  pn <- as.numeric(stats::predict(g, newdata = dn))
  r <- a27b_capture(plot_prediction_deviation_panels(g, data = dn, type = "gaussian"))
  expect_identical(r$warnings, character(0))
  expect_equal(a27b_panel(r$value)$fitted, a27b_wmean(pn, dn$stratum, dn$wg),
               tolerance = 1e-8)
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = dn, type = "gaussian"))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "prior weights are not known for 180 of the 180 rows",
               fixed = TRUE)
  expect_equal(a27b_panel(r$value)$fitted, a27b_wmean(pn, dn$stratum, rep(1, nrow(dn))),
               tolerance = 1e-8)
})

test_that("a weights column is never taken from outside `data`", {
  skip_on_cran()
  d <- a27b_data()
  g <- a27b_quiet(lme4::lmer(yg ~ x + (1 | stratum), data = d, weights = wg))
  # New rows without the weights column, beside a same-named object where the
  # formula can see it: the object has nothing to do with these rows, so the weights
  # are unknown and the rows are weighted equally.
  set.seed(7)
  dn <- d[, c("x", "yg", "stratum")]
  dn$x <- stats::rnorm(nrow(d))
  rownames(dn) <- NULL
  env <- environment(stats::formula(g))
  assign("wg", seq_len(nrow(d)), envir = env)
  on.exit(rm("wg", envir = env), add = TRUE)
  pn <- as.numeric(stats::predict(g, newdata = dn))
  r <- a27b_capture(plot_prediction_deviation_panels(g, data = dn, type = "gaussian"))
  expect_length(r$warnings, 1L)
  expect_equal(a27b_panel(r$value)$fitted, a27b_wmean(pn, dn$stratum, rep(1, nrow(dn))),
               tolerance = 1e-8)
})

test_that("a design-weighted fit's rows keep their sampling weights when reordered", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- a27b_data()
  set.seed(5)
  d$sw <- stats::runif(nrow(d), 0.2, 5)
  m <- a27b_quiet(fit_maihda(y01 ~ x + (1 | stratum), data = d, sampling_weights = "sw",
                             family = stats::binomial()))
  pd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(m)))
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = a27b_perm(m$data)))
  expect_identical(r$warnings, character(0))
  pr <- a27b_panel(r$value)
  expect_equal(pr$fitted, pd$fitted, tolerance = 1e-10)
  expect_equal(pr$abs_res_dev, pd$abs_res_dev, tolerance = 1e-10)
})

test_that("a real brms trials() fit keeps each row's trials when its rows are reordered", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). Every quantity is a
  # deterministic function of the fit's stored draws, so sampler quality is irrelevant.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  d <- a27b_data()
  m <- a27b_quiet(fit_maihda(s | trials(ntr) ~ x + (1 | stratum), data = d,
                             engine = "brms", family = "binomial", chains = 1,
                             iter = 600, warmup = 300, refresh = 0, seed = 1))
  pd <- a27b_panel(a27b_quiet(plot_prediction_deviation_panels(m)))
  r <- a27b_capture(plot_prediction_deviation_panels(m, data = a27b_perm(d)))
  expect_identical(r$warnings, character(0))
  pr <- a27b_panel(r$value)
  for (col in c("fitted", "abs_res_dev", "ci_lower", "ci_upper")) {
    expect_equal(pr[[col]], pd[[col]], tolerance = 1e-10, label = col)
  }
})

test_that("a weighted polr is weighted by its case weights on every route", {
  # stats::weights() has no polr method, so HEAD left a weighted polr unweighted on its
  # own rows; reading its supplied rows' weights = column without also reading its
  # frame's "(weights)" made one fit draw two panels (the user's choice: weight both).
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(11)
  d <- data.frame(stratum = factor(rep(1:6, each = 40)))
  d$x <- stats::rnorm(nrow(d), as.integer(d$stratum) / 3)
  d$w <- round(stats::runif(nrow(d), 0.3, 8), 2)
  lat <- 0.8 * d$x + stats::rnorm(6)[as.integer(d$stratum)] + stats::rlogis(nrow(d))
  d$y <- ordered(cut(lat, c(-Inf, 0, 1.5, Inf), labels = c("lo", "mid", "hi")))
  po <- a27b_quiet(MASS::polr(y ~ x + stratum, data = d, weights = w, Hess = TRUE))
  expect_null(stats::weights(po))
  score <- as.numeric(as.matrix(stats::predict(po, newdata = d, type = "probs")) %*% (1:3))
  oracle <- a27b_wmean(score, as.integer(d$stratum), d$w)
  # Not vacuous: the weights move these means.
  expect_gt(max(abs(oracle - a27b_wmean(score, as.integer(d$stratum), rep(1, nrow(d))))),
            1e-3)
  for (dd in list(NULL, d, a27b_perm(d))) {
    r <- a27b_capture(plot_prediction_deviation_panels(po, data = dd, type = "ordinal",
                                                       ordinal_mode = "expected_score"))
    expect_identical(r$warnings, character(0))
    expect_equal(a27b_panel(r$value)$fitted, oracle, tolerance = 1e-10)
  }
})

test_that("a fitted row is recognised by its name AND its prediction", {
  frame <- data.frame(stratum = c(1, 1, 2, 2), row.names = c("a", "b", "c", "d"))
  fit_v <- c(0.1, 0.2, 0.3, 0.4)
  dat <- frame[c("d", "b"), , drop = FALSE]
  expect_identical(maihda_prediction_panel_row_identity(frame, dat, fit_v, c(0.4, 0.2)),
                   c(4L, 2L))
  # A row whose name is no fitted row is NA ...
  dat2 <- rbind(dat, data.frame(stratum = 2, row.names = "z"))
  expect_identical(maihda_prediction_panel_row_identity(frame, dat2, fit_v,
                                                        c(0.4, 0.2, 0.9)),
                   c(4L, 2L, NA))
  # ... and one named row that is not that fitted row refutes every name.
  expect_null(maihda_prediction_panel_row_identity(frame, dat, fit_v, c(0.4, 0.25)))
  dat3 <- dat
  dat3$stratum[1] <- 1
  expect_null(maihda_prediction_panel_row_identity(frame, dat3, fit_v, c(0.4, 0.2)))
  expect_null(maihda_prediction_panel_row_identity(frame, dat, fit_v, c(NA, 0.2)))
})

test_that("controls: an unweighted Bernoulli panel is unchanged on its own rows", {
  skip_on_cran()
  d <- a27b_data()
  m <- a27b_quiet(fit_maihda(y01 ~ x + (1 | stratum), data = d, family = "binomial"))
  p <- as.numeric(stats::fitted(m$model))
  ref <- a27b_wmean(maihda_binomial_abs_deviance_residual(d$y01, p), d$stratum,
                    rep(1, nrow(d)))
  r0 <- a27b_capture(plot_prediction_deviation_panels(m))
  r1 <- a27b_capture(plot_prediction_deviation_panels(m, data = d))
  expect_identical(c(r0$warnings, r1$warnings), character(0))
  expect_equal(a27b_panel(r0$value)$abs_res_dev, ref, tolerance = 1e-12)
  expect_equal(a27b_panel(r1$value)$abs_res_dev, ref, tolerance = 1e-12)
})
