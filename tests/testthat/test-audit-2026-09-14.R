# Audit 2026-09-14: fit_maihda(engine = "wemix") forwarded WeMix::mix()'s
# center_grand / center_group, whose centring the package's predictions never applied
# -- CONFIRMED, wider in the routes that reach it, narrower in two consumers.
#
# mix() centres the covariates those arguments name INSIDE the fit -- subtracting their
# level-1 (sampling) weighted grand mean, or the weighted mean within each stratum --
# before the formula is evaluated, and keeps no centring constants.
# maihda_wemix_linpred() rebuilt the design from the RAW covariates and multiplied it
# by the centred-scale coefficients. On this file's fixture, against the unmodified
# code:
#   * predictions and stratum tables were off by exactly slope * centring mean: 1.706
#     for center_grand = ~ x and 2.262 for center_group = list(stratum = ~ x), against
#     WeMix's own fitted values (y - resid, aligned by position);
#   * binomial stratum standard errors were off by up to 32% (grand) and 38% (group);
#     the binomial AUC (0.679 against 0.697) and plot(type = "effect_decomp") were
#     wrong under GROUP centring only -- a grand constant shifts every row alike,
#     which the rank-based AUC and the within-stratum decomposition ignore;
#   * partial names (center_gra =, center_grou =) reached the same centring, and so
#     did maihda(), which forwards ... to both of its fits.
# Grand centring is a pure reparametrisation (fitted values unchanged); group
# centring is a different model (stratum variance 0.149 against 0.181). Either way,
# the covariate centred in `data` reproduces the centred fit exactly.
#
# FIX (chosen over recording the centring constants): maihda_fit_wemix() rejects both
# arguments and points to centring in `data`. Forwarded names are first resolved to
# the formals mix() binds them to (maihda_resolve_dot_names(), via match.call(); the
# audit 2026-09-14b pass generalised it from WeMix to every engine), which also closes
# the same hole in the max_iteration check:
# max_iter = 0 reached mix() as max_iteration = 0, and a binomial fit returned its
# unoptimised starting values (intercept -1.813 against the converged -1.659). A
# saved fit whose WeMix call records centring is refused by maihda_wemix_linpred()
# and by the effect_decomp plot, whose tryCatch would otherwise turn the refusal into
# NA points. A saved maihda() analysis keeps summaries computed when it was fitted:
# for a binomial fit its stored stratum SEs and intervals and the interaction table
# were wrong (grand centring: 0.060 adjusted, 0.114 null, 0.091 interactions, against
# a pre-centred twin; group centring also moved the stored AUC, 0.671 against 0.691)
# while the VPC, PCV and every stored Gaussian number were right, so print(),
# summary(), tidy(), glance() and plot() of such a binomial analysis warn (the user's
# choice over refusing or only documenting it). Fits without centring are
# bit-identical to the previous code.

a14_data <- function() {
  set.seed(20260913)
  n <- 600
  d <- data.frame(gender = factor(sample(c("f", "m"), n, TRUE)),
                  race = factor(sample(c("a", "b", "c"), n, TRUE)))
  d <- suppressWarnings(suppressMessages(make_strata(d, c("gender", "race"))))$data
  sid <- as.integer(factor(d$stratum))
  # x has mean about 3 and a different mean in each stratum, so grand and group
  # centring differ and both move the intercept.
  d$x <- stats::rnorm(n, mean = 3) + 0.4 * (sid - 3.5)
  u <- stats::rnorm(6, sd = 0.6)[sid]
  d$y <- 1 + 0.5 * d$x + u + stats::rnorm(n)
  d$yb <- as.integer(-1.5 + 0.5 * d$x + 1.5 * u + stats::rlogis(n) > 0)
  d$w <- stats::runif(n, 0.5, 2)
  d
}

a14_fit <- function(...) suppressWarnings(suppressMessages(fit_maihda(...)))

a14_wemix <- function(formula, data, ...) {
  a14_fit(formula, data = data, engine = "wemix", sampling_weights = "w", ...)
}

# WeMix's own fitted values. $resid follows the data rows by POSITION; its names are
# not the data's row names.
a14_oracle <- function(fit) fit$data$y - unname(fit$model$resid)

# A model saved before the fix: the package's object around a mix() fit that used
# centring, called the way fit_maihda() calls mix() -- do.call() with evaluated values.
a14_saved <- function(base, family = NULL, ...) {
  args <- list(formula = base$formula, data = base$data,
               weights = c("w", .maihda_wemix_l2_col), ...)
  if (!is.null(family)) {
    args$family <- family
  }
  base$model <- suppressWarnings(do.call(WeMix::mix, args))
  base
}

test_that("wemix rejects WeMix's centring arguments, under partial names too", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- a14_data()
  cases <- list(
    list(dots = list(center_grand = ~ x), arg = "'center_grand'"),
    list(dots = list(center_group = list(stratum = ~ x)), arg = "'center_group'"),
    list(dots = list(center_gra = ~ x), arg = "'center_grand'"),
    list(dots = list(center_grou = list(stratum = ~ x)), arg = "'center_group'"))
  for (case in cases) {
    expect_error(do.call(a14_wemix, c(list(y ~ x + (1 | stratum), d), case$dots)),
                 paste("does not accept", case$arg), fixed = TRUE,
                 info = names(case$dots))
  }
  # The advice names the sampling-weight column.
  expect_error(a14_wemix(y ~ x + (1 | stratum), d, center_grand = ~ x),
               "weighted.mean(data$x, data$w)", fixed = TRUE)
  # maihda() forwards ... to both of its fits.
  d2 <- d
  d2$stratum <- NULL
  expect_error(
    suppressWarnings(suppressMessages(maihda(y ~ x + (1 | gender:race), data = d2,
                                             engine = "wemix", sampling_weights = "w",
                                             center_grand = ~ x))),
    "does not accept 'center_grand'", fixed = TRUE)
  # An ambiguous prefix binds to no single formal, so R's own matching refuses it.
  expect_error(a14_wemix(y ~ x + (1 | stratum), d, center_gr = ~ x))
  # NULL is mix()'s own default -- no centring -- and fits as if absent.
  plain <- a14_wemix(y ~ x + (1 | stratum), d)
  with_null <- a14_wemix(y ~ x + (1 | stratum), d, center_grand = NULL)
  expect_identical(with_null$model$coef, plain$model$coef)
  expect_lt(max(abs(predict_maihda(with_null, scale = "link") - a14_oracle(with_null))),
            1e-10)
})

test_that("centring in data reproduces WeMix's centred fit, as the error advises", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  # Documents the advice rather than the defect: it holds with or without the fix.
  d <- a14_data()
  plain <- a14_wemix(y ~ x + (1 | stratum), d)
  expect_identical(nrow(plain$data), nrow(d))
  direct <- function(...) {
    suppressWarnings(do.call(WeMix::mix, list(
      formula = y ~ x + (1 | stratum), data = plain$data,
      weights = c("w", .maihda_wemix_l2_col), ...)))
  }
  grand <- direct(center_grand = ~ x)
  group <- direct(center_group = list(stratum = ~ x))
  d$x_grand <- d$x - stats::weighted.mean(d$x, d$w)
  d$x_group <- d$x - stats::ave(d$x * d$w, d$stratum, FUN = sum) /
    stats::ave(d$w, d$stratum, FUN = sum)
  pre_grand <- a14_wemix(y ~ x_grand + (1 | stratum), d)
  pre_group <- a14_wemix(y ~ x_group + (1 | stratum), d)
  # Two SEPARATE WeMix optimisations of the same model, so how many digits they share
  # is the platform's business, not the package's: CI's macOS build lands 4.2e-7 apart
  # in a variance component (7.6e-08 on 0.181) where Windows and Linux agree to 1e-10.
  # 1e-5 covers that and still convicts the alternative this documents -- internal
  # centring differing from centring in data would move a coefficient by slope * mean,
  # 1.706 here against an estimate of 2.409, some 70%.
  expect_equal(unname(pre_grand$model$coef), unname(grand$coef), tolerance = 1e-5)
  expect_equal(unname(pre_grand$model$vars), unname(grand$vars), tolerance = 1e-5)
  expect_equal(unname(pre_group$model$coef), unname(group$coef), tolerance = 1e-5)
  expect_equal(unname(pre_group$model$vars), unname(group$vars), tolerance = 1e-5)
  expect_lt(max(abs(predict_maihda(pre_grand, scale = "link") - a14_oracle(pre_grand))),
            1e-10)
  expect_lt(max(abs(predict_maihda(pre_group, scale = "link") - a14_oracle(pre_group))),
            1e-10)
  # Grand centring reparametrises the uncentred model; group centring changes it.
  expect_equal(a14_oracle(pre_grand), a14_oracle(plain), tolerance = 1e-6)
  expect_gt(max(abs(a14_oracle(pre_group) - a14_oracle(plain))), 0.01)
})

test_that("a saved fit that used WeMix's centring refuses to predict", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- a14_data()
  plain <- a14_wemix(y ~ x + (1 | stratum), d)
  saved_grand <- a14_saved(plain, center_grand = ~ x)
  saved_group <- a14_saved(plain, center_group = list(stratum = ~ x))
  # mix() records a partial name under the full one.
  saved_partial <- a14_saved(plain, center_gra = ~ x)
  msg <- "internal centring"
  expect_error(predict_maihda(saved_grand), msg, fixed = TRUE)
  expect_error(predict_maihda(saved_partial, scale = "link"), msg, fixed = TRUE)
  group_summary <- summary(saved_group)
  expect_error(maihda_stratum_predictions_wemix(saved_group, group_summary, scale = "link"),
               msg, fixed = TRUE)
  expect_error(plot(saved_group, type = "effect_decomp"), msg, fixed = TRUE)
  # A binomial fit's stratum standard errors are built from the same predictions.
  plain_b <- a14_fit(yb ~ x + (1 | stratum), data = d, engine = "wemix",
                     sampling_weights = "w", family = "binomial")
  saved_b <- a14_saved(plain_b, family = stats::binomial(), center_grand = ~ x)
  expect_error(summary(saved_b), msg, fixed = TRUE)
  # A recorded NULL is no centring: that fit still predicts.
  saved_null <- a14_saved(plain, center_grand = NULL)
  expect_lt(max(abs(predict_maihda(saved_null, scale = "link") - a14_oracle(saved_null))),
            1e-10)
})

test_that("a saved binomial analysis whose WeMix fits used centring warns when shown", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- a14_data()
  d$stratum <- NULL
  quiet_maihda <- function(...) suppressWarnings(suppressMessages(maihda(...)))
  saved_analysis <- function(analysis, family, ...) {
    analysis$model <- a14_saved(analysis$model, family = family, ...)
    analysis$model_adjusted <- a14_saved(analysis$model_adjusted, family = family, ...)
    analysis
  }
  binomial <- quiet_maihda(yb ~ x + (1 | gender:race), data = d, engine = "wemix",
                           sampling_weights = "w", family = "binomial")
  expect_no_warning(utils::capture.output(print(binomial)))
  expect_no_warning(glance(binomial))

  # Its stored stratum SEs, intervals and interaction tests came from the uncentred
  # predictions; the methods that show them say so.
  msg <- "internal centring"
  grand <- saved_analysis(binomial, stats::binomial(), center_grand = ~ x)
  expect_warning(utils::capture.output(print(grand)), msg, fixed = TRUE)
  expect_warning(summary(grand, which = "adjusted"), msg, fixed = TRUE)
  expect_warning(tidy(grand), msg, fixed = TRUE)
  expect_warning(glance(grand), "intervals and interaction tests were", fixed = TRUE)
  expect_warning(plot(grand, type = "vpc"), msg, fixed = TRUE)
  # Group centring moved the AUC too; a grand constant leaves its ranks alone.
  group <- saved_analysis(binomial, stats::binomial(), center_group = list(stratum = ~ x))
  expect_warning(glance(group), "interaction tests and AUC were", fixed = TRUE)

  # A Gaussian analysis's stored summaries never used those predictions.
  gaussian <- quiet_maihda(y ~ x + (1 | gender:race), data = d, engine = "wemix",
                           sampling_weights = "w")
  gaussian <- saved_analysis(gaussian, NULL, center_group = list(stratum = ~ x))
  expect_no_warning(utils::capture.output(print(gaussian)))
  expect_no_warning(glance(gaussian))
})

test_that("a partially named max_iteration is checked like the full name", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- a14_data()
  expect_error(a14_wemix(y ~ x + (1 | stratum), d, max_iter = 0), "whole number",
               fixed = TRUE)
  fit3 <- a14_wemix(y ~ x + (1 | stratum), d, max_iter = 3)
  # Checked, then forwarded under the full name as the integer the check coerces to.
  expect_identical(fit3$model$call$max_iteration, 3L)
})

test_that("forwarded names resolve to the formals mix() binds them to", {
  skip_if_not_installed("WeMix")
  supplied <- c("formula", "data", "weights")
  resolved <- maihda_resolve_dot_names(
    list(center_gra = 1, max_iter = 2, nQuad = 9, TRUE), WeMix::mix, supplied)
  expect_identical(names(resolved), c("center_grand", "max_iteration", "nQuad", "cWeights"))
  expect_identical(unname(resolved), list(1, 2, 9, TRUE))
  # A name binding to no formal (family is supplied by a binomial fit) or to several
  # is left as given, for mix() to refuse.
  expect_identical(names(maihda_resolve_dot_names(list(fam = 1), WeMix::mix,
                                                  c(supplied, "family"))), "fam")
  expect_identical(names(maihda_resolve_dot_names(list(center_gr = 1), WeMix::mix,
                                                  supplied)), "center_gr")
  expect_identical(maihda_resolve_dot_names(list(), WeMix::mix, supplied), list())
})
