# Audit 2026-09-13 (third pass): "A7. Symmetric ordinal thresholds do not
# generally restore marginal proportional odds" -- CONFIRMED, on both counts.
#
# ?maihda_proportional_odds_test (and the parallel comment above
# maihda_ordinal_po_stat()) said that an exactly symmetric threshold
# configuration is an exception, "where the marginal slopes coincide and the
# fixed-only statistic is valid", because logit E[plogis(eta - u)] is odd for
# symmetric u. Oddness makes the derivative even, which equates the marginal
# slopes at cut points -c and +c only where their arguments are opposite, i.e.
# where the location x'beta is 0. It says nothing about the rest of the covariate
# range, and nothing about observations that share a stratum.
#
# Measured for this pass (quadrature, the expected-data LRT, and simulation):
#   * cut points -1/+1, unit coefficient, tau = 1: marginal slopes at x = 1 are
#     -0.8816705 and -0.8264839 (the audit's figures, reproduced to 4e-8 by
#     Gauss-Hermite and by integrate(), the derivative checked by finite
#     differences); equal at x = 0, mirror images elsewhere. The old
#     comment's own "slope spread 0.027 for (-1.5, 0, 1.5)" was a symmetric
#     configuration with unequal slopes.
#   * Threshold symmetry is neither sufficient nor necessary for a zero
#     population gap between the best nominal and proportional fixed-only fits.
#     Cut points -1/+1 with a 0/1 covariate (unit coefficient, tau = 1): gap
#     3.4e-05 per observation, which a non-central chi-squared approximation turns
#     into about 44% rejection at n = 96,000; cut points -0.5/1.5 with the same
#     covariate: gap 0. What matters is the design as a whole, and whether cut
#     points look symmetric is a matter of coding (block 2).
#   * Even a design symmetric as a whole -- three categories cut at -1 and +1, a
#     standard normal covariate with coefficient 0.8, zero population gap -- is
#     not rescued once the random effect is realised in finitely many strata:
#     with 12 strata the chi-squared reference rejected data simulated from the
#     correctly specified model at 4.9% / 6.2% / 10.0% for n = 1,440 / 24,000 /
#     96,000 at tau = 0.5, and 6.6% / 25.0% / 53.2% at tau = 1 (2,000, 1,000 and
#     600 datasets per cell). A fresh tau = 0 control at n = 96,000 rejected 4.8%.
#
# FIX: the exception is gone from both places, the neighbouring sentence no
# longer claims the slopes differ "for any non-zero stratum variance", and the
# help page now says why symmetric thresholds do not rescue the chi-squared
# reference. Code is unchanged. A rejection rate is a long simulation, not a unit
# test, so the blocks below pin the wording and the two identities the new text
# relies on; the audit log keeps the runs.

test_that("the proportional-odds help page claims no symmetric-threshold exception", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree (test_local / load_all);
  # R CMD check copies only tests/, so there is nothing to read there.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  # Whitespace-normalised, because Rd re-wraps and a phrase can straddle a line.
  po <- gsub("[[:space:]]+", " ", paste(
    readLines(file.path(man, "maihda_proportional_odds_test.Rd"), warn = FALSE),
    collapse = " "))

  # The replacement says what symmetry does and does not give ...
  expect_true(grepl("Symmetric thresholds do not rescue the chi-squared reference.",
                    po, fixed = TRUE))
  expect_true(grepl(paste0("which equates the marginal slopes at thresholds \\eqn{-c} ",
                           "and \\eqn{+c} where the location \\eqn{x'\\beta}{x'beta} ",
                           "is zero, but not over the range of the covariates"),
                    po, fixed = TRUE))
  expect_true(grepl(paste0("Whether fitted thresholds look symmetric depends on how ",
                           "the covariates are coded"), po, fixed = TRUE))
  expect_true(grepl("without changing the fit or the statistic", po, fixed = TRUE))
  expect_true(grepl(paste0("no threshold configuration makes observations that share ",
                           "a stratum independent"), po, fixed = TRUE))
  expect_true(grepl("The chi-squared reference also treats the observations as independent",
                    po, fixed = TRUE))
  expect_true(grepl("Even when both the thresholds and the covariate are symmetric about zero",
                    po, fixed = TRUE))
  expect_true(grepl("generally differ across thresholds once the stratum variance is non-zero",
                    po, fixed = TRUE))

  # ... and no spelling of the old claim survives. Ban the words rather than the
  # sentences, so a reworded relapse is caught too; the one legitimate "valid" is
  # removed first.
  expect_false(grepl("exception", po, ignore.case = TRUE))
  expect_false(grepl("coincide", po, ignore.case = TRUE))
  expect_false(grepl("share a slope", po, fixed = TRUE))
  expect_false(grepl("valid", gsub("is not valid here", "", po, fixed = TRUE),
                     ignore.case = TRUE))
  expect_false(grepl("for any non-zero stratum variance", po, fixed = TRUE))
  expect_false(grepl("without bound", po, fixed = TRUE))
  expect_false(grepl("arises in practice", po, fixed = TRUE))
})

test_that("threshold symmetry is a property of the coding, not of the fit or the statistic", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  set.seed(7)
  strata <- expand.grid(a = factor(1:4), b = factor(1:3))
  d <- strata[rep(seq_len(nrow(strata)), each = 100), ]
  idx <- interaction(d$a, d$b, drop = TRUE)
  u <- stats::rnorm(nlevels(idx), 0, 1)[as.integer(idx)]
  d$g01 <- stats::rbinom(nrow(d), 1, 0.5)
  lat <- d$g01 + u + stats::rlogis(nrow(d))
  d$y <- factor(cut(lat, c(-Inf, -1, 1, Inf), labels = FALSE), levels = 1:3,
                ordered = TRUE)
  fit <- function(term) suppressMessages(suppressWarnings(fit_maihda(
    stats::as.formula(paste("y ~", term, "+ (1 | a:b)")), data = d,
    family = "ordinal")))

  m01 <- fit("g01")
  th01 <- maihda_clmm_cutpoints(m01$model)
  beta <- unname(m01$model$beta)
  # Shift the covariate so that the fitted cut points come out symmetric: adding
  # s to it adds beta * s to every cut point.
  s <- -mean(th01) / beta
  d$gs <- d$g01 + s
  d$gc <- d$g01 - 0.5
  ms <- fit("gs")
  mc <- fit("gc")
  ths <- maihda_clmm_cutpoints(ms$model)
  thc <- maihda_clmm_cutpoints(mc$model)

  # fixture premise: the 0/1 and centred codings are clearly asymmetric here
  expect_gt(abs(mean(th01)), 0.1)
  expect_gt(abs(mean(thc)), 0.1)
  # the shifted coding is symmetric, and each shift moves every cut point by beta * s
  expect_lt(abs(mean(ths)), 1e-3)
  expect_equal(unname(ths - th01), rep(beta * s, 2), tolerance = 1e-3)
  expect_equal(unname(thc - th01), rep(-0.5 * beta, 2), tolerance = 1e-3)

  # ... while the model and the fixed-only statistic do not move at all.
  for (m in list(ms, mc)) {
    expect_equal(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(m01$model)),
                 tolerance = 1e-7)
    expect_equal(unname(m$model$beta), beta, tolerance = 1e-4)
  }
  po01 <- maihda_ordinal_po_stat(m01$model)
  expect_false(is.null(po01))
  expect_equal(maihda_ordinal_po_stat(ms$model)$lrt, po01$lrt, tolerance = 1e-8)
  expect_equal(maihda_ordinal_po_stat(mc$model)$lrt, po01$lrt, tolerance = 1e-8)
})

test_that("symmetric cut points share a marginal slope where the location is zero, not beyond", {
  gh <- maihda_gauss_hermite_normal(80)
  # d/dx logit E[plogis(cut - beta * x - u)], u ~ N(0, tau^2), by quadrature
  slope <- function(cut, x, beta = 1, tau = 1) {
    eta <- cut - beta * x - tau * gh$nodes
    p <- sum(gh$weights * stats::plogis(eta))
    -beta * sum(gh$weights * stats::dlogis(eta)) / (p * (1 - p))
  }
  # the local identity the old help page generalised ...
  expect_equal(slope(-1, 0), slope(1, 0), tolerance = 1e-12)
  expect_equal(slope(-1, 1), slope(1, -1), tolerance = 1e-12)
  # ... and the audit's counterexample one unit away from it
  expect_equal(c(slope(-1, 1), slope(1, 1)), c(-0.8816705, -0.8264839), tolerance = 1e-6)
  expect_gt(abs(slope(-1, 1) - slope(1, 1)), 0.05)
  # no stratum variance, no attenuation: every slope is -beta
  expect_equal(c(slope(-1, 1, tau = 0), slope(1, 1, tau = 0)), c(-1, -1), tolerance = 1e-12)
})
