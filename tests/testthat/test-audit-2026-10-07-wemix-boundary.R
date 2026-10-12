# Audit 2026-10-07 (packet 2026-10-07-wemix-boundary): fit_maihda(engine = "wemix",
# family = "binomial") stopped with "Matrix seems negative semi-definite" whenever the
# between-stratum variance was estimated at zero -- CONFIRMED, and WeMix's own:
# WeMix::mix() called directly on the same rows raises it.
#
# WeMix 4.0.3 optimises a binomial fit in the Cholesky factor, where a zero variance is
# an ordinary stationary point, and its Newton loop converges there (to 1.8e-15 on the
# reported data). It then takes a numerical Hessian in the VARIANCE, by central
# differences of half-width (1 + v) * eps^(1/4) = 1.22e-4: the lower point is a negative
# variance, and Matrix::nearPD() stops on it. So every estimate below 0.000122 fails
# (measured: fails at 0.000084, returns at 0.00045), and the additive case -- an
# adjusted model that explains all of the between-stratum variance, where lme4 returns
# 0 and a PCV of 1 -- cannot be returned by this engine. The Gaussian path
# differentiates in the Cholesky factor and returns such fits (10 of 10 measured,
# flagged singular, PCV 1).
#
# Per the user the fit stays an error, one that now says what happened.
# maihda_fit_wemix() recognises WeMix's error by its erroring call (nearPD, on the
# binomial path) and not by its text, which Matrix translates ("Matrix scheint negativ
# semidefinit zu sein"), and stops with the diagnosis, what a zero variance means and
# the engine that returns a fit there. Every other WeMix error and every returned fit
# is as WeMix produced it. That the error implies a zero variance was checked on 40
# fits (36 with non-unit weights) against glmer with prior weights and 13-point
# quadrature, which maximises the same criterion: WeMix errored in the 16 where that
# variance is at most 0.000084 and returned in the 24 where it is at least 0.00045.

wb07_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The error a call stops with, or NULL when it returns -- so a call that stops (or does
# not) fails each assertion it feeds rather than aborting the block.
wb07_error <- function(expr) tryCatch({ expr; NULL }, error = function(e) e)
wb07_text <- function(e) if (inherits(e, "error")) conditionMessage(e) else ""
wb07_head <- function(e) {
  cl <- if (inherits(e, "error")) conditionCall(e)
  if (is.call(cl)) cl[[1L]]
}

# Eight strata of `n` rows whose success counts are additive on the probability scale
# (12, 16, 16, 20, 16, 20, 20, 24 of 40). The additive logit model leaves far less
# between-stratum variation than binomial noise, so the adjusted model's variance is
# zero while the null model keeps a positive one. No random numbers.
wb07_data <- function(n = 40) {
  g <- expand.grid(d1 = 0:1, d2 = 0:1, d3 = 0:1)
  g$s <- round(n * (0.3 + 0.1 * (g$d1 + g$d2 + g$d3)))
  g$stratum <- seq_len(nrow(g))
  d <- g[rep(1:8, each = n), c("d1", "d2", "d3", "stratum")]
  d$y <- unlist(lapply(1:8, function(j) rep(c(1L, 0L), c(g$s[j], n - g$s[j]))))
  for (v in c("d1", "d2", "d3")) d[[v]] <- factor(d[[v]])
  d$w <- 1
  d
}

wb07_adj <- y ~ d1 + d2 + d3 + (1 | stratum)
wb07_diagnosis <- "estimated the between-stratum variance of this binomial model at zero"

test_that("the boundary error is told by its erroring call and the family, not its text", {
  near_pd <- quote(nearPD(X, posd.tol = tol * sqrt(.Machine$double.eps)))
  binom <- stats::binomial(link = "logit")
  words <- "Matrix seems negative semi-definite"

  # WeMix's condition as it arrives, in English and as Matrix translates it.
  expect_true(maihda_wemix_boundary_error(simpleError(words, near_pd), binom))
  expect_true(maihda_wemix_boundary_error(
    simpleError("Matrix scheint negativ semidefinit zu sein", near_pd), binom))
  # A namespace-qualified call names the same function.
  qualified <- as.call(list(call("::", as.name("Matrix"), as.name("nearPD")),
                            as.name("X")))
  expect_true(maihda_wemix_boundary_error(simpleError("any text", qualified), binom))

  # The same words from anywhere else are not it, nor is an error without a call.
  expect_false(maihda_wemix_boundary_error(
    simpleError(words, quote(solve.default(a))), binom))
  expect_false(maihda_wemix_boundary_error(simpleError(words), binom))
  expect_false(maihda_wemix_boundary_error(simpleError(words, quote(nearPD2(X))), binom))
  # The linear path returns a zero-variance fit, so nearPD() stopping there is
  # another matter.
  expect_false(maihda_wemix_boundary_error(simpleError(words, near_pd),
                                           stats::gaussian()))

  # The message quotes whatever WeMix said.
  expect_match(maihda_wemix_boundary_message(simpleError("WeMix's own words")),
               "(\"WeMix's own words\")", fixed = TRUE)
})

test_that("a binomial wemix fit at a zero variance stops with the diagnosis", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- wb07_data()

  # The state under test, shown without WeMix: the adjusted model's variance is zero.
  ref <- wb07_quiet(lme4::glmer(wb07_adj, data = d, family = stats::binomial()))
  expect_lt(as.numeric(lme4::VarCorr(ref)$stratum), 1e-8)

  # WeMix's own error on these rows, with no package code in between. If this stops
  # failing, WeMix returns zero-variance binomial fits and the catch in
  # maihda_fit_wemix() is no longer reached.
  d2 <- d
  d2$l2 <- 1
  own <- wb07_error(wb07_quiet(WeMix::mix(
    wb07_adj, data = d2, weights = c("w", "l2"),
    family = stats::binomial(link = "logit"))))
  expect_s3_class(own, "error")
  expect_identical(wb07_head(own), as.name("nearPD"))

  err <- wb07_error(fit_maihda(wb07_adj, data = d, engine = "wemix",
                               family = "binomial", sampling_weights = "w"))
  expect_s3_class(err, "error")
  expect_match(wb07_text(err), wb07_diagnosis, fixed = TRUE)
  expect_match(wb07_text(err), "that is a PCV of 1", fixed = TRUE)
  expect_match(wb07_text(err), "engine = \"brms\" with 'sampling_weights'",
               fixed = TRUE)
  # WeMix's words are quoted, in whatever language Matrix put them.
  expect_match(wb07_text(err), paste0("(\"", wb07_text(own), "\")"), fixed = TRUE)

  # Non-unit weights that leave every stratum's weighted proportion where it was: the
  # same zero variance, the same diagnosis.
  d$w <- rep(c(2, 1), length.out = nrow(d))
  err_w <- wb07_error(fit_maihda(wb07_adj, data = d, engine = "wemix",
                                 family = "binomial", sampling_weights = "w"))
  expect_match(wb07_text(err_w), wb07_diagnosis, fixed = TRUE)
})

test_that("maihda() passes the diagnosis on when its adjusted model is the one at zero", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- wb07_data()
  d$stratum <- NULL

  err <- wb07_error(wb07_quiet(maihda(
    y ~ 1 + (1 | d1:d2:d3), data = d, engine = "wemix", family = "binomial",
    sampling_weights = "w")))
  expect_match(wb07_text(err), wb07_diagnosis, fixed = TRUE)
})

test_that("a binomial fit WeMix returns is passed on as WeMix made it", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- wb07_data()
  f_null <- y ~ 1 + (1 | stratum)

  fit <- wb07_quiet(fit_maihda(f_null, data = d, engine = "wemix",
                               family = "binomial", sampling_weights = "w"))
  d$l2 <- 1
  own <- wb07_quiet(WeMix::mix(f_null, data = d, weights = c("w", "l2"),
                               family = stats::binomial(link = "logit")))
  for (field in c("coef", "vars", "SE", "lnl", "theta", "invHessian", "cov_mat")) {
    expect_equal(fit$model[[field]], own[[field]], tolerance = 1e-8, info = field)
  }
  # The fixture's null model keeps a positive variance: only its adjusted model is
  # at zero.
  expect_gt(as.numeric(fit$model$vars[1]), 0.01)
})

test_that("other WeMix errors and the Gaussian path are left alone", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- wb07_data()

  # An error of WeMix's that is not the boundary one leaves in WeMix's words.
  other <- wb07_error(fit_maihda(wb07_adj, data = d, engine = "wemix",
                                 family = "binomial", sampling_weights = "w",
                                 nQuad = 0))
  expect_match(wb07_text(other), "nQuad", fixed = TRUE)
  expect_false(grepl("between-stratum variance", wb07_text(other), fixed = TRUE))

  # A Gaussian outcome whose stratum means follow the additive model exactly: the
  # adjusted variance is zero, and the fit returns, flagged as a boundary fit.
  d$z <- 0.5 * (d$d1 == "1") + 0.45 * (d$d2 == "1") + 0.55 * (d$d3 == "1") +
    rep(seq(-1, 1, length.out = 40), times = 8)
  gauss <- tryCatch(
    wb07_quiet(fit_maihda(z ~ d1 + d2 + d3 + (1 | stratum), data = d,
                          engine = "wemix", sampling_weights = "w")),
    error = function(e) e)
  expect_s3_class(gauss, "maihda_model")
  expect_true(isTRUE(as.numeric(gauss$model$vars[1]) < 1e-8))
  expect_true(isTRUE(gauss$diagnostics$singular))
})
