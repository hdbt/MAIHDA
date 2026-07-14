# Regression tests for the 2026-07-19 (P1/P1/P2) audit findings on commit 32bd266.
#
#   #1 P1  fit_maihda() accepted duplicate intercept-only stratum terms, e.g.
#          y ~ (1 | stratum) + (1 | stratum). lme4 fits these as two components
#          ("stratum", "stratum.1") and splits the between-stratum variance
#          arbitrarily; summary() then counted only the first as between-stratum
#          variance and misclassified the rest, giving a non-identifiable VPC.
#   #2 P1  PCV boundary detection (maihda_pcv_null_at_boundary) was gated to
#          engine == "lme4", so a wemix/ordinal fit with a tiny-but-positive
#          between-stratum variance (e.g. 3e-9 from clmm, 5e-24 from WeMix) was
#          not flagged: it passed the > 0 denominator guard and could yield an
#          exploding / nonsensical PCV, and adjusted_at_boundary was never TRUE.
#   #3 P2  With estimation = "ML", maihda_pcv_refit_ml() skips refitML() when a
#          model's stratum variance is at the boundary (leaving it on REML) but
#          recorded nothing, so a partly-REML comparison was still labelled a
#          pure, correction-free "ML" comparison.

# ---- #1 P1: duplicate stratum random-effect terms are rejected --------------

test_that("fit_maihda() rejects duplicate intercept-only (1 | stratum) terms", {
  set.seed(1)
  ns <- 30; npg <- 20
  d <- data.frame(y = rnorm(ns * npg),
                  stratum = factor(rep(seq_len(ns), each = npg)))
  expect_error(
    fit_maihda(y ~ (1 | stratum) + (1 | stratum), data = d),
    "non-identifiable")
})

test_that("fit_maihda() rejects a duplicated intersectional shorthand", {
  set.seed(2)
  d <- data.frame(y = rnorm(600),
                  gender = factor(rep(c("M", "F"), 300)),
                  race   = factor(rep(1:6, 100)))
  # Two identical shorthand terms would build one stratum but leave a duplicated
  # (1 | stratum) grouping; rejected (here by the pre-existing single-term guard).
  expect_error(
    fit_maihda(y ~ (1 | gender:race) + (1 | gender:race), data = d))
})

test_that("fit_maihda() still accepts a single stratum / shorthand term", {
  set.seed(3)
  d <- data.frame(y = rnorm(600),
                  gender = factor(rep(c("M", "F"), 300)),
                  race   = factor(rep(1:6, 100)))
  expect_no_error(suppressMessages(fit_maihda(y ~ (1 | gender:race), data = d)))
  # A single stratum intercept plus a genuine random slope is a different (and
  # already-rejected) structure -- the duplicate guard must not pre-empt that
  # more specific intercept-only error path here; the single-term fit is enough.
})

# ---- #2 P1: engine-agnostic boundary detection ------------------------------

test_that("maihda_variance_at_boundary applies lme4's relative SD tolerance", {
  # sd_stratum / sd_residual < 1e-4 -> boundary
  expect_true(MAIHDA:::maihda_variance_at_boundary(1e-9, 1))     # ratio ~3.2e-5
  expect_true(MAIHDA:::maihda_variance_at_boundary(0, 1))        # exact zero
  expect_true(MAIHDA:::maihda_variance_at_boundary(-1, 1))       # negative
  expect_false(MAIHDA:::maihda_variance_at_boundary(0.01, 1))    # ratio 0.1
  expect_false(MAIHDA:::maihda_variance_at_boundary(1, (pi^2)/3))# healthy latent
  # No usable residual scale -> absolute near-zero fallback.
  expect_true(MAIHDA:::maihda_variance_at_boundary(1e-30, NA))
  expect_false(MAIHDA:::maihda_variance_at_boundary(0.5, NA))
})

test_that("ordinal PCV rejects a boundary (near-zero) between-stratum denominator", {
  skip_if_not_installed("ordinal")
  set.seed(101)
  G <- 40; npg <- 40
  si <- rep(seq_len(G), each = npg)
  do <- data.frame(gender = factor(rep(c("M", "F"), length.out = G)[si]),
                   race   = factor(rep(1:20, length.out = G)[si]),
                   age    = rnorm(G * npg))
  eta <- 0.6 * do$age + rlogis(G * npg)                # NO stratum effect
  do$y <- ordered(cut(eta, quantile(eta, c(0, .33, .66, 1)),
                      include.lowest = TRUE, labels = c("lo", "mid", "hi")))
  mo0 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | gender:race), data = do, engine = "ordinal")))
  mo1 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race), data = do, engine = "ordinal")))
  expect_lt(MAIHDA:::maihda_clmm_variances(mo0)$stratum, 1e-4)
  expect_true(MAIHDA:::maihda_pcv_null_at_boundary(mo0))
  # Previously returned a nonsensical PCV (~ -110%) from ~0 / ~0; now a clean stop.
  expect_error(calculate_pcv(mo0, mo1), "zero boundary")
})

test_that("a healthy ordinal fit is NOT flagged at the boundary", {
  skip_if_not_installed("ordinal")
  set.seed(202)
  G <- 40; npg <- 40
  si <- rep(seq_len(G), each = npg)
  u  <- rnorm(G, 0, 1.3)                                # real between-stratum effect
  age <- rnorm(G * npg)
  eta <- 0.5 * age + u[si] + rlogis(G * npg)
  do <- data.frame(gender = factor(rep(c("M", "F"), length.out = G)[si]),
                   race   = factor(rep(1:20, length.out = G)[si]),
                   y = ordered(cut(eta, quantile(eta, c(0, .33, .66, 1)),
                                   include.lowest = TRUE, labels = c("lo", "mid", "hi"))))
  mo <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | gender:race), data = do, engine = "ordinal")))
  expect_false(MAIHDA:::maihda_pcv_null_at_boundary(mo))
})

test_that("wemix PCV flags a boundary between-stratum variance", {
  skip_if_not_installed("WeMix")
  set.seed(3)
  G <- 30; npg <- 80
  si <- rep(seq_len(G), each = npg)
  dw <- data.frame(gender = factor(rep(c("M", "F"), length.out = G)[si]),
                   race   = factor(rep(seq_len(G / 2), length.out = G)[si]),
                   age    = rnorm(G * npg),
                   w      = runif(G * npg, 0.5, 2))
  dw$y <- 1 + 0.5 * dw$age + rnorm(G * npg, 0, 1)       # NO stratum effect
  mw <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race), data = dw, sampling_weights = "w")))
  # WeMix returns a tiny positive (~1e-20 here), not exact zero; still boundary.
  expect_lt(MAIHDA:::maihda_wemix_variances(mw)$stratum, 1e-4)
  expect_true(MAIHDA:::maihda_pcv_null_at_boundary(mw))
})

# ---- #3 P2: honest provenance when the ML refit is skipped at the boundary ---

test_that("maihda_pcv_estimation_used resolves the basis actually used", {
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("fitted", FALSE), "fitted")
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("fitted", TRUE),  "fitted")
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("ML", FALSE),     "ML")
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("ML", TRUE),      "mixed")
})

# A stratum-level covariate x drives ALL between-stratum variation, so the null
# has real between-stratum variance but the adjusted model's is ~0 (boundary).
maihda_boundary_adjusted_df <- function(seed = 42) {
  set.seed(seed)
  G <- 60; npg <- 25
  xstr <- rnorm(G)
  si <- rep(seq_len(G), each = npg)
  data.frame(gender = factor(rep(c("M", "F"), length.out = G)[si]),
             race   = factor(rep(1:10, length.out = G)[si]),
             edu    = factor(rep(1:3, length.out = G)[si]),
             x = xstr[si],                               # constant within stratum
             y = 2 * xstr[si] + rnorm(G * npg, 0, 1))
}

test_that("estimation='ML' with a boundary adjusted model reports a mixed basis", {
  df <- maihda_boundary_adjusted_df()
  m0 <- fit_maihda(y ~ (1 | gender:race:edu), data = df)
  m1 <- fit_maihda(y ~ x + (1 | gender:race:edu), data = df)
  p <- calculate_pcv(m0, m1, estimation = "ML")
  expect_identical(p$estimation, "ML")          # requested basis preserved
  expect_identical(p$estimation_used, "mixed")  # basis actually used
  expect_true(p$adjusted_at_boundary)
  out <- paste(capture.output(print(p)), collapse = "\n")
  expect_match(out, "Variance basis: mixed")
  expect_false(grepl("correction-free", out))   # no false pure-ML claim
})

test_that("a non-boundary estimation='ML' comparison stays a pure ML basis", {
  set.seed(9)
  G <- 60; npg <- 25
  u <- rnorm(G, 0, 1.5)
  si <- rep(seq_len(G), each = npg)
  df <- data.frame(gender = factor(rep(c("M", "F"), length.out = G)[si]),
                   race   = factor(rep(1:10, length.out = G)[si]),
                   edu    = factor(rep(1:3, length.out = G)[si]),
                   x = rnorm(G * npg))
  df$y <- 1 + u[si] + 0.3 * df$x + rnorm(G * npg, 0, 1)
  m0 <- fit_maihda(y ~ (1 | gender:race:edu), data = df)
  m1 <- fit_maihda(y ~ x + (1 | gender:race:edu), data = df)
  p <- calculate_pcv(m0, m1, estimation = "ML")
  expect_identical(p$estimation_used, "ML")
  out <- paste(capture.output(print(p)), collapse = "\n")
  expect_match(out, "Variance basis: ML-refit")
  expect_match(out, "correction-free")
})

# A stratum-level covariate x drives ALL between-stratum variation, so a model that
# includes x has an EXACT REML-zero between-stratum variance (a singular fit at the
# boundary as fitted, before any ML refit). Under estimation = "ML" that model's ML
# refit is therefore SKIPPED (it keeps its REML fit), making the series' basis "mixed"
# -- as opposed to a model the ML refit merely pushes toward the boundary, which is a
# genuine (pure-ML) refit. x2 is an individual-level nuisance that leaves between-
# stratum variance untouched (so a subset with x2 but not x is not at the boundary).
maihda_boundary_series_df <- function(seed = 5) {
  set.seed(seed)
  G <- 60; npg <- 20
  base <- data.frame(gender = factor(rep(c("M", "F"), length.out = G)),
                     race   = factor(rep(1:10, length.out = G)),
                     edu    = factor(rep(1:3, length.out = G)))
  d <- base[rep(seq_len(G), each = npg), ]
  d <- make_strata(d, c("gender", "race", "edu"))$data
  slev <- sort(unique(as.character(d$stratum)))
  xmap <- stats::setNames(rnorm(length(slev)), slev)
  d$x  <- xmap[as.character(d$stratum)]   # constant within stratum
  d$x2 <- rnorm(nrow(d))                  # individual-level nuisance
  d$y  <- 2 * d$x + rnorm(nrow(d), 0, 1)
  rownames(d) <- NULL
  d
}

test_that("stepwise_pcv reports a mixed basis when a step hits the boundary under ML", {
  d <- maihda_boundary_series_df()
  sw <- suppressWarnings(stepwise_pcv(d, "y", c("x2", "x"), estimation = "ML"))
  # The step adding x drives the between-stratum variance to the exact boundary.
  expect_true(length(attr(sw, "boundary_steps")) > 0)
  expect_identical(attr(sw, "estimation"), "ML")
  expect_identical(attr(sw, "estimation_used"), "mixed")
  out <- paste(capture.output(print(sw)), collapse = "\n")
  expect_match(out, "Variance basis: mixed")
})

test_that("pcv_importance reports a mixed basis when a subset is at the boundary under ML", {
  d <- maihda_boundary_series_df()
  imp <- suppressWarnings(pcv_importance(d, "y", c("x", "x2"), estimation = "ML"))
  expect_true(isTRUE(imp$full_at_boundary))
  expect_identical(imp$estimation, "ML")
  expect_identical(imp$estimation_used, "mixed")
  out <- paste(capture.output(print(imp)), collapse = "\n")
  expect_match(out, "Variance basis: mixed")
})
