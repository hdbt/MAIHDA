# Regression tests for the 2026-07-13 (P1/P1/P1/P2/P2/P3) audit findings.
#
#   #1 P1  External/formula offsets were dropped (or errored) in downstream prediction
#          paths (stratum expected counts, longitudinal count VPC, fixed trajectory,
#          prediction-deviation panels) -- fitting honoured the offset but predict() on
#          the stored model frame did not.
#   #2 P1  compare_maihda_groups() built shared numeric strata on ALL input rows, not
#          the pooled analytic sample, so subset= disagreed with prefiltered input.
#   #3 P1  compare_maihda_groups()/maihda() omitted the external offset from family
#          detection and the per-group analytic n.
#   #4 P2  A failed lme4::refitML() was silently kept as REML but still labelled ML.
#   #5 P2  The adjusted-model boundary flag was gated on estimation = "ML".
#   #6 P3  Subsetting a maihda_ic result dropped ic_primary but kept the class, so
#          print.maihda_ic() errored on a zero-length condition.

# ---- #1 P1: offsets retained in prediction paths ----------------------------

test_that("maihda_fitted_offset recovers external and formula offsets", {
  set.seed(1); n <- 250
  d <- data.frame(g1 = factor(sample(c("F", "M"), n, TRUE)),
                  g2 = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = rnorm(n), logE = log(runif(n, 1, 8)))
  d$y <- rpois(n, exp(-0.3 + 0.2 * d$x + d$logE))
  m_ext  <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, family = "poisson", offset = logE)))
  m_form <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + offset(logE) + (1 | g1:g2), data = d, family = "poisson")))
  m_none <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, family = "poisson")))

  off_ext  <- MAIHDA:::maihda_fitted_offset(m_ext$model)
  off_form <- MAIHDA:::maihda_fitted_offset(m_form$model)
  expect_equal(off_ext, d$logE, ignore_attr = TRUE)
  expect_equal(off_form, d$logE, ignore_attr = TRUE)
  expect_null(MAIHDA:::maihda_fitted_offset(m_none$model))
})

test_that("maihda_lme4_fixed_link matches predict(re.form=NA) and adds the offset", {
  set.seed(2); n <- 300
  d <- data.frame(g1 = factor(sample(c("F", "M"), n, TRUE)),
                  g2 = factor(sample(c("A", "B"), n, TRUE)),
                  x = rnorm(n), z = factor(sample(c("a", "b", "c"), n, TRUE)),
                  logE = log(runif(n, 1, 8)))
  d$y <- rpois(n, exp(-0.2 + 0.25 * d$x + 0.1 * (d$z == "b") + d$logE))
  m_none <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + z + (1 | g1:g2), data = d, family = "poisson")))
  nd <- m_none$data; nd$x <- nd$x + 0.3
  # No offset: the helper reproduces predict(newdata, re.form = NA) exactly.
  expect_equal(
    MAIHDA:::maihda_lme4_fixed_link(m_none$model, nd),
    as.numeric(stats::predict(m_none$model, newdata = nd, re.form = NA, type = "link")),
    tolerance = 1e-10)
  # Supplying an offset just shifts the linear predictor by that vector.
  off <- runif(nrow(nd))
  expect_equal(
    MAIHDA:::maihda_lme4_fixed_link(m_none$model, nd, offset = off),
    MAIHDA:::maihda_lme4_fixed_link(m_none$model, nd) + off, tolerance = 1e-12)
})

test_that("lme4 stratum expected counts retain the offset (external and formula agree)", {
  set.seed(3); n <- 600
  d <- data.frame(g1 = factor(sample(c("F", "M"), n, TRUE)),
                  g2 = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = rnorm(n), logE = log(runif(n, 1, 20)))
  strat <- interaction(d$g1, d$g2, drop = TRUE)
  d$y <- rpois(n, exp(-0.3 + 0.25 * d$x + rnorm(nlevels(strat), 0, 0.4)[strat] + d$logE))
  m_ext  <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, family = "poisson", offset = logE)))
  m_form <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + offset(logE) + (1 | g1:g2), data = d, family = "poisson")))
  s_ext  <- suppressWarnings(suppressMessages(summary(m_ext)))
  s_form <- suppressWarnings(suppressMessages(summary(m_form)))

  # The formula-offset path previously ERRORED here (the stored model frame names the
  # offset column "offset(logE)", not "logE"); it must now succeed.
  pe <- MAIHDA:::maihda_stratum_predictions_lme4(m_ext, s_ext, "response")
  pf <- MAIHDA:::maihda_stratum_predictions_lme4(m_form, s_form, "response")

  # External and formula offsets fit identically, so the per-stratum expected counts
  # must match, and match the offset-inclusive fitted() means (the external path used
  # to drop the offset, biasing the counts).
  i <- match(pe$stratum, pf$stratum)
  expect_equal(pe$predicted_row, pf$predicted_row[i], tolerance = 1e-6)
  gt <- tapply(stats::fitted(m_ext$model), as.character(m_ext$data$stratum), mean)
  j <- match(names(gt), pe$stratum)
  expect_equal(pe$predicted_row[j], as.numeric(gt), tolerance = 1e-6)
})

test_that("prediction-deviation panels retain the offset and reject external-offset newdata", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("patchwork")
  skip_if_not_installed("dplyr"); skip_if_not_installed("ggrepel")
  skip_if_not_installed("tidyr")
  set.seed(4); n <- 400
  d <- data.frame(g1 = factor(sample(c("F", "M"), n, TRUE)),
                  g2 = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = rnorm(n), logE = log(runif(n, 1, 10)))
  d$y <- rpois(n, exp(-0.2 + 0.2 * d$x + d$logE))
  m_ext  <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, family = "poisson", offset = logE)))
  m_form <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + offset(logE) + (1 | g1:g2), data = d, family = "poisson")))
  # Fitted-row panels (data = NULL) must build for BOTH offset styles (formula offset
  # used to error, external offset dropped the offset).
  expect_s3_class(
    suppressWarnings(plot_prediction_deviation_panels(m_ext, type = "poisson")),
    "patchwork")
  expect_s3_class(
    suppressWarnings(plot_prediction_deviation_panels(m_form, type = "poisson")),
    "patchwork")
  # A genuine external newdata cannot reconstruct an external offset -> rejected.
  expect_error(
    plot_prediction_deviation_panels(m_ext, data = d, type = "poisson"),
    "external offset")
})

test_that("longitudinal count VPC retains the offset (external and formula agree)", {
  set.seed(5); nid <- 130
  g1 <- sample(c("F", "M"), nid, TRUE); g2 <- sample(c("A", "B", "C"), nid, TRUE)
  base <- do.call(rbind, lapply(seq_len(nid), function(i)
    data.frame(id = i, g1 = g1[i], g2 = g2[i], time = 0:2,
               logE = log(runif(3, 1, 6)), stringsAsFactors = FALSE)))
  base$y <- rpois(nrow(base), exp(-0.2 + 0.15 * base$time + base$logE))
  m_ext <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ 1 + (1 | g1:g2), data = base, family = "poisson",
               id = "id", time = "time", offset = logE)))
  m_form <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ 1 + offset(logE) + (1 | g1:g2), data = base, family = "poisson",
               id = "id", time = "time")))
  s_ext  <- suppressWarnings(suppressMessages(summary(m_ext)))
  s_form <- suppressWarnings(suppressMessages(summary(m_form)))
  # The count resid-grid path (maihda_longitudinal_resid_grid_lme4) runs for both, and
  # identical fits give an identical time-varying VPC (the external path used to drop
  # the offset from the grid, biasing VPC(t); formula would have errored).
  expect_s3_class(s_ext$longitudinal$vpc_t, "data.frame")
  expect_equal(s_ext$longitudinal$vpc_t, s_form$longitudinal$vpc_t, tolerance = 1e-6)

  # 1c: the fixed-part trajectory (mean covariate profile + mean offset) is finite for
  # both offset styles and agrees (previously dropped/errored on the offset).
  grid <- sort(unique(base$time))
  te <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m_ext, grid)
  tf <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m_form, grid)
  expect_true(all(is.finite(te)))
  expect_equal(te, tf, tolerance = 1e-6)
})

# ---- #2 P1: shared binning uses the analytic sample -------------------------

test_that("compare_maihda_groups shared binning: subset= matches prefiltered input", {
  set.seed(11); N <- 900
  dd <- data.frame(country = sample(c("A", "B"), N, TRUE),
                   gender = sample(c("F", "M"), N, TRUE),
                   income = rnorm(N, 50, 15))
  istr <- interaction(dd$gender, cut(dd$income, 3), drop = TRUE)
  dd$y <- rnorm(nlevels(istr), 0, 0.5)[istr] + rnorm(N, 0, 1)
  keep <- dd$income < 60
  r_sub <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ (1 | gender:income), data = dd, group = "country",
                          subset = keep)))
  r_pre <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ (1 | gender:income), data = dd[keep, ], group = "country")))
  expect_equal(r_sub$vpc, r_pre$vpc, tolerance = 1e-8)
  expect_equal(r_sub$var_between, r_pre$var_between, tolerance = 1e-8)
  expect_equal(r_sub$n, r_pre$n)
})

# ---- #3 P1: offset folded into compare/maihda family + analytic n -----------

test_that("an offset=NA row does not flip family detection in maihda()/compare", {
  set.seed(303); n <- 400
  d <- data.frame(country = sample(c("A", "B"), n, TRUE),
                  g1 = sample(c("F", "M"), n, TRUE),
                  g2 = sample(c("X", "Y"), n, TRUE), x = rnorm(n))
  d$yf <- factor(sample(c("low", "mid"), n, TRUE),
                 levels = c("low", "mid", "high"), ordered = TRUE)
  d$off <- rnorm(n, 0, 0.3)
  d$yf[1] <- "high"; d$off[1] <- NA  # 3rd category only on the offset-NA row

  # maihda(): the analytic sample (offset-NA row dropped) is a 2-level ordered factor
  # -> binomial. Without folding the offset in, the workflow saw 3 levels and selected
  # the ordinal engine, then errored against fit_maihda()'s binary analytic sample.
  m <- suppressWarnings(suppressMessages(
    maihda(yf ~ x + (1 | g1:g2), data = d, offset = off)))
  expect_identical(m$model$family$family, "binomial")

  # compare_maihda_groups(): same detection, and the per-group analytic n excludes the
  # offset-NA row.
  cmp <- suppressWarnings(suppressMessages(
    compare_maihda_groups(yf ~ x + (1 | g1:g2), data = d, group = "country",
                          offset = off)))
  expect_identical(attr(cmp, "family"), "binomial")
  expect_equal(sum(cmp$n), n - 1L)
})

# ---- #4 P2: a failed ML refit is surfaced, not silently kept as REML ---------

test_that("a failed lme4::refitML() warns and is flagged, not labelled ML", {
  set.seed(41); n <- 400
  d <- data.frame(g1 = sample(c("F", "M"), n, TRUE),
                  g2 = sample(c("A", "B", "C"), n, TRUE), x = rnorm(n))
  strat <- interaction(d$g1, d$g2, drop = TRUE)
  d$y <- 1 + 0.3 * d$x + rnorm(nlevels(strat), 0, 0.8)[strat] + rnorm(n, 0, 0.6)
  m_null <- suppressWarnings(suppressMessages(fit_maihda(y ~ (1 | g1:g2), data = d)))
  m_adj  <- suppressWarnings(suppressMessages(fit_maihda(y ~ x + (1 | g1:g2), data = d)))

  testthat::local_mocked_bindings(
    refitML = function(object, ...) stop("mocked refitML failure"),
    .package = "lme4")

  expect_warning(
    r <- calculate_pcv(m_null, m_adj, estimation = "ML"),
    "refitML\\(\\) failed|not a pure ML")
  expect_true(isTRUE(r$ml_refit_failed))

  # maihda_ic() must not rank AIC across a mixed REML/ML basis when a refit fell back.
  expect_warning(
    ic <- suppressMessages(maihda_ic(m_null, m_adj)),
    "estimation basis|ML refit failed")
  expect_false("delta" %in% names(ic))
})

# ---- #5 P2: adjusted boundary flagged under any estimation basis -------------

test_that("calculate_pcv flags an adjusted boundary fit under the default (fitted) basis", {
  # Small-n additive strata (no interaction): the null (1 | stratum) has clear between-
  # stratum variance, but the additive full model leaves ~0 residual interaction
  # variance -> a singular (boundary) adjusted fit with PCV ~ 100%.
  set.seed(7101); n <- 240
  d <- data.frame(gender = sample(c("m", "f"), n, TRUE),
                  race = sample(c("a", "b", "c"), n, TRUE), age = rnorm(n))
  d <- make_strata(d, c("gender", "race"))$data
  d$y <- 1 + ifelse(d$gender == "m", 0.8, 0) + c(a = 0, b = 0.9, c = -0.4)[d$race] +
    0.3 * d$age + rnorm(n)
  m_null <- suppressWarnings(suppressMessages(fit_maihda(y ~ 1 + (1 | stratum), data = d)))
  m_full <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ gender + race + age + (1 | stratum), data = d)))

  # A singular adjusted fit is indistinguishable from genuinely additive strata (the
  # common MAIHDA case), so calculate_pcv() does NOT warn under any estimation basis
  # (the previous warning was gated on "ML", leaving the default basis with no flag at
  # all). It records the status on the result and surfaces it honestly only in print().
  expect_no_warning(r <- calculate_pcv(m_null, m_full))
  expect_true(isTRUE(r$adjusted_at_boundary))
  expect_output(print(r), "singularity boundary")

  # The status also propagates to stepwise_pcv() and pcv_importance() as SILENT
  # attributes surfaced in their print methods (never as a per-call warning).
  sp <- suppressWarnings(suppressMessages(stepwise_pcv(d, "y", c("gender", "race", "age"))))
  expect_true(length(attr(sp, "boundary_steps")) >= 1)
  expect_no_warning(pim <- suppressMessages(pcv_importance(d, "y", c("gender", "race", "age"))))
  expect_true(isTRUE(pim$full_at_boundary))
})

# ---- #6 P3: subsetting a maihda_ic result preserves metadata -----------------

test_that("[.maihda_ic preserves ic_primary so print() still works after subset", {
  set.seed(61); n <- 300
  d <- data.frame(g1 = sample(c("F", "M"), n, TRUE),
                  g2 = sample(c("A", "B", "C"), n, TRUE), x = rnorm(n))
  strat <- interaction(d$g1, d$g2, drop = TRUE)
  d$y <- 1 + 0.3 * d$x + rnorm(nlevels(strat), 0, 0.5)[strat] + rnorm(n, 0, 0.5)
  m1 <- suppressWarnings(suppressMessages(fit_maihda(y ~ x + (1 | g1:g2), data = d)))
  m2 <- suppressWarnings(suppressMessages(fit_maihda(y ~ 1 + (1 | g1:g2), data = d)))
  ic <- suppressMessages(maihda_ic(m1, m2))

  sub <- ic[, c("model", "estimator", "delta")]
  expect_s3_class(sub, "maihda_ic")
  expect_false(is.null(attr(sub, "ic_primary")))
  # Previously printed the table and THEN errored on a zero-length `!is.na(NULL)`.
  expect_output(print(sub), "MAIHDA Information Criteria")
  expect_error(print(sub), NA)
})
