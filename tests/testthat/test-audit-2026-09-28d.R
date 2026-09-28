# Audit 2026-09-28d: external finding F02 -- the analytic-sample checks behind
# calculate_pcv(), compare_maihda() and maihda_ic() never saw a brms fit's binomial
# trial counts. CONFIRMED, and wider: a brms weights() term was invisible too.
#
# brms keeps a response's addition terms out of the terms of its model frame, so
# model.response() on a `y | trials(n) ~ ...` fit returns the successes alone (the
# frame's terms are `y ~ y + n + x + stratum`), and a brmsfit has no weights() method,
# so the prior-weight fingerprint read "unit". Two fits of the same successes and row
# ids out of 12 and out of 24 trials therefore passed as one analytic sample. Measured
# on real Stan fits (brms 2.23, 2 chains of 1500 iterations, standata() confirming 12
# vs 24): calculate_pcv() returned PCV 0.389 (0.9517 -> 0.5813) where the same
# adjustment on the 12-trial data gives -0.025, maihda_ic() reported a LOOIC delta of
# 110.6, and compare_maihda() warned about nothing. Varying trials, sampling-weighted fits and reordered rows all
# passed the same way. The same blind spot passed two brms `y | weights(w)` fits with
# different weights, which lme4's weights= check refuses.
#
# FIX: maihda_brms_addition_values() reads a trials() or weights() term's values as
# brms passed them to Stan (standata(), falling back to the stored frame). The
# response fingerprint encodes a brms trials() response as successes and failures --
# the pair an lme4 cbind(successes, failures) response already is -- and the
# prior-weight fingerprint takes a brms weights() term (not the package's own
# weights(.maihda_sw), which the sampling-weight check keys). Every fit is built with
# brm(empty = TRUE): a genuine brmsfit and frame with no Stan compile, and the
# fingerprints never read draws.

a28d_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression, or NULL where it errors, so an assertion it feeds fails
# rather than aborting the block.
a28d_or_null <- function(expr) tryCatch(expr, error = function(e) NULL)

# The PCV sample checks' refusal for two models, or "" when they pass.
a28d_pcv_refusal <- function(a, b) {
  tryCatch({
    validate_pcv_models(a, b)
    ""
  }, error = function(e) conditionMessage(e))
}

# 12 strata x 25 rows: successes out of 12 trials, a Gaussian outcome, weights.
a28d_data <- function() {
  set.seed(20260928)
  d <- data.frame(stratum = rep(seq_len(12), each = 25))
  d$x <- stats::rnorm(nrow(d))
  u <- stats::rnorm(12, 0, 0.7)
  d$y <- stats::rbinom(nrow(d), 12L, stats::plogis(-0.3 + 0.4 * d$x + u[d$stratum]))
  d$n <- 12L
  d$g <- 1 + 0.5 * d$x + u[d$stratum] + stats::rnorm(nrow(d))
  d$w1 <- stats::runif(nrow(d), 0.5, 2)
  d$w2 <- stats::runif(nrow(d), 0.5, 2)
  d$sw <- stats::runif(nrow(d), 0.5, 3)
  d$yb <- stats::rbinom(nrow(d), 1L, stats::plogis(0.2 * d$x + u[d$stratum]))
  d
}

a28d_brms <- function(formula, data, family, ...) {
  a28d_or_null(a28d_quiet(fit_maihda(formula, data = data, engine = "brms",
                                     family = family, empty = TRUE, ...)))
}

a28d_fp <- function(m) a28d_or_null(maihda_wrapper_response_fingerprint(m))
a28d_wfp <- function(m) a28d_or_null(maihda_weight_fingerprint(m$model))

test_that("premise: brms reports the successes alone as the response of a trials() fit", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- a28d_data()
  dB <- d
  dB$n <- 24L
  mA <- a28d_brms(y | trials(n) ~ x + (1 | stratum), d, "binomial")
  mB <- a28d_brms(y | trials(n) ~ x + (1 | stratum), dB, "binomial")
  frame <- maihda_model_frame(mA$model)
  expect_identical(unname(stats::model.response(frame)), d$y)
  # brms's own data: the same successes, out of different trials.
  expect_identical(brms::standata(mA$model)$Y, brms::standata(mB$model)$Y)
  expect_equal(as.numeric(brms::standata(mA$model)$trials), rep(12, 300))
  expect_equal(as.numeric(brms::standata(mB$model)$trials), rep(24, 300))
})

test_that("the same successes out of different trials are a different analytic sample", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- a28d_data()
  f <- y | trials(n) ~ x + (1 | stratum)
  mA <- a28d_brms(f, d, "binomial")
  dB <- d
  dB$n <- 24L
  mB <- a28d_brms(f, dB, "binomial")
  set.seed(5)
  dV <- d
  dV$n <- sample(12:20, nrow(d), replace = TRUE)
  mV <- a28d_brms(f, dV, "binomial")

  expect_false(identical(a28d_fp(mA), a28d_fp(mB)))
  expect_false(identical(a28d_fp(mA), a28d_fp(mV)))
  expect_match(maihda_ic_delta_issues(list(mA, mB)), "analytic sample", all = FALSE)
  expect_match(maihda_ic_delta_issues(list(mA, mV)), "analytic sample", all = FALSE)
  expect_error(calculate_pcv(mA, mB), "outcome values differ")
  expect_error(calculate_pcv(mA, mV), "outcome values differ")
  # The weights are untouched: it is the outcome that differs.
  expect_identical(a28d_wfp(mA), "unit")
  expect_identical(a28d_wfp(mB), "unit")
})

test_that("design weights and row order do not hide, or invent, a trials difference", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- a28d_data()
  f <- y | trials(n) ~ x + (1 | stratum)
  dB <- d
  dB$n <- 24L
  sA <- a28d_brms(f, d, "binomial", sampling_weights = "sw")
  sB <- a28d_brms(f, dB, "binomial", sampling_weights = "sw")
  expect_false(identical(a28d_fp(sA), a28d_fp(sB)))
  expect_error(calculate_pcv(sA, sB), "outcome values differ")

  mA <- a28d_brms(f, d, "binomial")
  set.seed(3)
  perm <- sample(nrow(d))
  mP <- a28d_brms(f, d[perm, ], "binomial")
  mPB <- a28d_brms(f, dB[perm, ], "binomial")
  expect_false(identical(row.names(mA$data), row.names(mP$data)))
  expect_identical(a28d_fp(mA), a28d_fp(mP))
  expect_length(maihda_ic_delta_issues(list(mA, mP)), 0L)
  expect_identical(a28d_pcv_refusal(mA, mP), "")
  expect_false(identical(a28d_fp(mA), a28d_fp(mPB)))
  expect_match(a28d_pcv_refusal(mA, mPB), "outcome values differ")
  # A null and an adjusted model on the same data: the comparison MAIHDA exists for.
  m0 <- a28d_brms(y | trials(n) ~ 1 + (1 | stratum), d, "binomial")
  expect_length(maihda_ic_delta_issues(list(m0, mA)), 0L)
  expect_identical(a28d_pcv_refusal(m0, mA), "")
})

test_that("equivalent encodings of one binomial dataset share a fingerprint", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- a28d_data()
  mA <- a28d_brms(y | trials(n) ~ x + (1 | stratum), d, "binomial")
  # The lme4 spelling of the same observations: successes and failures.
  mL <- a28d_quiet(fit_maihda(cbind(y, n - y) ~ x + (1 | stratum), data = d,
                              family = "binomial"))
  expect_identical(a28d_fp(mA), a28d_fp(mL))
  dd <- d
  dd$n <- as.numeric(dd$n)
  expect_identical(a28d_fp(mA), a28d_fp(a28d_brms(y | trials(n) ~ x + (1 | stratum), dd,
                                                  "binomial")))
  expect_identical(a28d_fp(mA), a28d_fp(a28d_brms(y | trials(12) ~ x + (1 | stratum), d,
                                                  "binomial")))
})

test_that("trials() are read as brms evaluated them, a global function included", {
  skip_on_cran()
  skip_if_not_installed("brms")
  # brms evaluates an addition term with a lookup that reaches the global environment;
  # the package's own evaluation on the stored frame sees base functions only (so it
  # cannot pick up a stray global object), and on its own would return no trials at
  # all here. The function must live in the global environment for brms to fit it.
  skip_if(exists("a28d_plus_one", envir = globalenv(), inherits = FALSE),
          "the global environment already has an a28d_plus_one")
  assign("a28d_plus_one", function(z) z + 1L, envir = globalenv())
  on.exit(if (exists("a28d_plus_one", envir = globalenv(), inherits = FALSE)) {
    rm("a28d_plus_one", envir = globalenv())
  }, add = TRUE)

  d <- a28d_data()
  mA <- a28d_brms(y | trials(n) ~ x + (1 | stratum), d, "binomial")
  mF <- a28d_brms(y | trials(a28d_plus_one(n)) ~ x + (1 | stratum), d, "binomial")
  d13 <- d
  d13$m <- d13$n + 1L
  m13 <- a28d_brms(y | trials(m) ~ x + (1 | stratum), d13, "binomial")
  expect_s3_class(mF, "maihda_model")
  expect_null(maihda_trials_from_formula(mF$model$formula, mF$data))
  expect_equal(a28d_or_null(maihda_brms_addition_values(mF$model, "trials", mF$data)),
               rep(13, 300))
  expect_identical(a28d_fp(mF), a28d_fp(m13))
  expect_false(identical(a28d_fp(mF), a28d_fp(mA)))
})

test_that("a brms weights() term enters the prior-weight check, as lme4 weights= do", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- a28d_data()
  f <- g | weights(w1) ~ x + (1 | stratum)
  wA <- a28d_brms(f, d, "gaussian")
  wB <- a28d_brms(f, transform(d, w1 = w2), "gaussian")
  expect_false(identical(a28d_wfp(wA), a28d_wfp(wB)))
  expect_match(maihda_ic_delta_issues(list(wA, wB)), "prior weights", all = FALSE)
  expect_error(calculate_pcv(wA, wB), "same prior weights")
  # lme4 has always refused the same pair.
  lA <- a28d_quiet(fit_maihda(g ~ x + (1 | stratum), data = d, weights = w1))
  lB <- a28d_quiet(fit_maihda(g ~ x + (1 | stratum), data = transform(d, w1 = w2),
                              weights = w1))
  expect_error(calculate_pcv(lA, lB), "same prior weights")

  # The same weights, a null and an adjusted model: comparable.
  w0 <- a28d_brms(g | weights(w1) ~ 1 + (1 | stratum), d, "gaussian")
  expect_identical(a28d_wfp(w0), a28d_wfp(wA))
  expect_identical(a28d_pcv_refusal(w0, wA), "")
  # weights(scale = TRUE) are the weights brms fitted: w and 2w are one likelihood.
  fS <- g | weights(w1, scale = TRUE) ~ x + (1 | stratum)
  sA <- a28d_brms(fS, d, "gaussian")
  expect_identical(a28d_wfp(sA), a28d_wfp(a28d_brms(fS, transform(d, w1 = 2 * w1),
                                                     "gaussian")))
  expect_false(identical(a28d_wfp(sA), a28d_wfp(a28d_brms(fS, transform(d, w1 = w2),
                                                           "gaussian"))))
  d$one <- 1
  expect_identical(a28d_wfp(a28d_brms(g | weights(one) ~ x + (1 | stratum), d,
                                      "gaussian")), "unit")
  # Beside trials(), each difference is caught by its own check.
  fT <- y | trials(n) + weights(w1) ~ x + (1 | stratum)
  tA <- a28d_brms(fT, d, "binomial")
  tW <- a28d_brms(fT, transform(d, w1 = w2), "binomial")
  tN <- a28d_brms(fT, transform(d, n = 24L), "binomial")
  expect_identical(a28d_fp(tA), a28d_fp(tW))
  expect_match(a28d_pcv_refusal(tA, tW), "same prior weights")
  expect_identical(a28d_wfp(tA), a28d_wfp(tN))
  expect_match(a28d_pcv_refusal(tA, tN), "outcome values differ")
})

test_that("controls: sampling weights keep their own check; no trials() term, no change", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- a28d_data()
  f <- y | trials(n) ~ x + (1 | stratum)
  sA <- a28d_brms(f, d, "binomial", sampling_weights = "sw")
  sB <- a28d_brms(f, transform(d, sw = w2), "binomial", sampling_weights = "sw")
  # weights(.maihda_sw) are sampling weights, keyed once, by their own check.
  expect_identical(a28d_wfp(sA), "unit")
  expect_match(a28d_pcv_refusal(sA, sB), "same sampling weights")
  expect_identical(maihda_ic_delta_issues(list(sA, sB)), "sampling weights")
  # A Bernoulli brms fit keeps the successes-only fingerprint.
  mb <- a28d_brms(yb ~ x + (1 | stratum), d, "binomial")
  frame <- maihda_model_frame(mb$model)
  y01 <- maihda_order_by_ids(unname(stats::model.response(frame)), row.names(frame))
  expect_identical(a28d_fp(mb), paste(formatC(y01, format = "g", digits = 12),
                                      collapse = "\r"))
  expect_null(maihda_brms_addition_values(mb$model, "trials", frame))
  expect_null(maihda_brms_addition_values(mb$model, "weights", frame))
})

test_that("the addition-term parser reads either term off a formula or a brmsformula", {
  skip_if_not_installed("brms")
  expect_identical(maihda_brms_addition_arg(y | trials(n) ~ x, "trials"), as.name("n"))
  expect_identical(maihda_brms_addition_arg(y | trials(n) + weights(w) ~ x, "weights"),
                   as.name("w"))
  expect_identical(maihda_brms_addition_arg(y | weights(w) + trials(n) ~ x, "trials"),
                   as.name("n"))
  expect_identical(maihda_brms_addition_arg(y | weights(w, scale = TRUE) ~ x, "weights"),
                   as.name("w"))
  expect_identical(maihda_brms_addition_arg(brms::bf(y | trials(12) ~ x), "trials"), 12)
  expect_null(maihda_brms_addition_arg(y | weights(w) ~ x, "trials"))
  expect_null(maihda_brms_addition_arg(y ~ x, "trials"))
  expect_null(maihda_brms_addition_arg(~ x, "weights"))
  expect_identical(maihda_find_trials_expr(quote(trials(n) + weights(w))), as.name("n"))
  # Only a brms fit carries its terms this way.
  expect_null(maihda_brms_addition_values(stats::lm(dist ~ speed, cars), "weights", cars))

  # Where brms::standata() cannot run (here a stand-in brmsfit with no fit behind it),
  # the term is evaluated on the stored frame instead.
  frame <- data.frame(y = c(1, 2, 3), n = c(4, 5, 6), w = c(0.5, 1, 2))
  stub <- function(f) structure(list(formula = brms::bf(f), data = frame), class = "brmsfit")
  expect_error(suppressWarnings(brms::standata(stub(y | trials(n) ~ 1))))
  expect_equal(maihda_brms_addition_values(stub(y | trials(n) ~ 1), "trials", frame),
               c(4, 5, 6))
  expect_equal(maihda_brms_addition_values(stub(y | trials(5) ~ 1), "trials", frame),
               c(5, 5, 5))
  expect_equal(maihda_brms_addition_values(stub(y | trials(n) + weights(w) ~ 1), "weights",
                                           frame), c(0.5, 1, 2))
  expect_null(maihda_brms_addition_values(stub(y | weights(.maihda_sw) ~ 1), "weights",
                                          frame))
  expect_null(maihda_brms_addition_values(stub(y | trials(no_such) ~ 1), "trials", frame))
})

test_that("real brms fits: different trials are refused end to end (opt-in Stan)", {
  # Compiles two Stan models, so OPT-IN (set MAIHDA_TEST_BRMS=true). The refusals do
  # not depend on the draws, so sampler quality is irrelevant.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  d <- a28d_data()
  dB <- d
  dB$n <- 24L
  fit <- function(f, dd) a28d_quiet(fit_maihda(f, data = dd, engine = "brms",
                                               family = "binomial", chains = 1,
                                               iter = 600, warmup = 300, refresh = 0,
                                               seed = 1))
  m0 <- fit(y | trials(n) ~ 1 + (1 | stratum), d)
  mB <- fit(y | trials(n) ~ x + (1 | stratum), dB)
  # The trials the check reads are the data's own, row for row.
  expect_equal(maihda_brms_addition_values(m0$model, "trials", m0$data), as.numeric(d$n))
  expect_equal(maihda_brms_addition_values(mB$model, "trials", mB$data), as.numeric(dB$n))
  expect_error(calculate_pcv(m0, mB), "outcome values differ")
  # Collect every warning: short chains add brms's own Rhat warning to both calls.
  warnings_of <- function(expr) {
    w <- character(0)
    value <- withCallingHandlers(expr, warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = w)
  }
  ic <- warnings_of(maihda_ic(m0, mB))
  expect_false("delta" %in% names(ic$value))
  expect_true(any(grepl("analytic sample", ic$warnings, fixed = TRUE)))
  cmp <- warnings_of(compare_maihda(m0, mB, ic = FALSE))
  expect_true(any(grepl("analytic sample", cmp$warnings, fixed = TRUE)))
  expect_length(maihda_ic_delta_issues(list(m0, m0)), 0L)
})
