# Audit 2026-09-17 (b): engine = "wemix" accepted a formula offset() that
# WeMix::mix() never fits -- CONFIRMED, and wider than reported (binomial too).
#
# WeMix 4.0.3's mix() takes the response from its model frame and the fixed design
# from lme4's X. Neither its linear nor its adaptive (binomial) likelihood has an
# offset term; only the unweighted lmer()/glmer() fit it draws starting values from
# sees one. An offset fit therefore returned the coefficients and variance components
# of the model WITHOUT the offset (to 1.4e-12 and 2.4e-8 on the auditor's data), while
# maihda_wemix_linpred() added the offset to every prediction and the convergence
# check still reported TRUE. On the auditor's Gaussian data, whose offset varies
# between strata: stratum variance 3.8246 against 0.2835 for the equivalent
# response-minus-offset fit, VPC 0.7993 against 0.2279, predictions off by up to
# 14.54. At unit weights, where the pseudo-likelihood is the likelihood,
# lmer(REML = FALSE) with the offset agrees with the response-minus-offset fit
# (predictions to 3.9e-8) and the package's offset fit agrees with lmer WITHOUT the
# offset. Binomial, which the audit did not test: the offset fit equals the
# offset-free fit (to 1.9e-8), and at unit weights glmer(nAGQ = 13) with the offset
# gives slope 0.547 and stratum variance 0.336 where the offset fit gave 1.330 and
# 0.979, the offset-free values, with probabilities off by up to 0.22.
#
# test-design-weights.R rebuilt the offset fit's prediction as X %*% beta + offset + u
# from the fit's own coefficients, so it compared the predictions with the same wrong
# fit; it now checks the refusal.
#
# FIX (the user's choice over fitting a Gaussian response minus the offset
# internally): the wemix engine refuses a formula offset() in fit_maihda(), and up
# front in maihda() and compare_maihda_groups() -- the latter would otherwise report a
# failed fit, with a warning, for every group. The message carries the exact Gaussian
# recipe and points a binary outcome to brms, whose Stan code with sampling weights
# includes the offset. The test is the terms() offset attribute that model.offset()
# reads. Fits saved before the fix are left as they were (the user's choice); fits
# without an offset are identical to the previous code.

b17_msg <- "does not support offset() terms"

# The auditor's design: 24 strata of 30, an offset linear in x plus a stratum shift.
b17_data <- function() {
  set.seed(682)
  d <- data.frame(stratum = factor(rep(1:24, each = 30)))
  d$x <- stats::rnorm(nrow(d))
  d$w <- stats::runif(nrow(d), 0.5, 2)
  d$off <- 3 + 2 * d$x + stats::rnorm(24, sd = 1.8)[d$stratum]
  d$y_adj <- 1 + 0.5 * d$x + stats::rnorm(24, sd = 0.5)[d$stratum] +
    stats::rnorm(nrow(d))
  d$y <- d$y_adj + d$off
  d
}

# Shorthand strata, a binary outcome and a comparison group, for the entry points.
b17_dim_data <- function() {
  set.seed(170917)
  n <- 360
  d <- data.frame(a = factor(sample(c("p", "q"), n, TRUE)),
                  b = factor(sample(c("r", "s", "t"), n, TRUE)),
                  grp = sample(c("g1", "g2"), n, TRUE))
  d$x <- stats::rnorm(n)
  d$w <- stats::runif(n, 0.5, 2)
  d$off <- stats::rnorm(n)
  d$o1 <- stats::rnorm(n)
  d$o2 <- stats::rnorm(n)
  d$y <- 1 + 0.3 * d$x + d$off + stats::rnorm(n)
  d$yb <- stats::rbinom(n, 1, 0.4)
  d
}

b17_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("WeMix::mix() fits the same model with and without a formula offset", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- b17_data()
  d$l2 <- 1
  mix <- function(f) b17_quiet(WeMix::mix(f, data = d, weights = c("w", "l2")))
  stratum_var <- function(m) m$varDF$vcov[m$varDF$grp == "stratum"][1]
  with_off <- mix(y ~ x + offset(off) + (1 | stratum))
  without <- mix(y ~ x + (1 | stratum))
  # If these start to differ, WeMix fits offsets now and the refusal in
  # maihda_wemix_check_offset() can be revisited.
  expect_equal(with_off$coef, without$coef, tolerance = 1e-6)
  expect_equal(stratum_var(with_off), stratum_var(without), tolerance = 1e-6)

  # The model the offset formula names -- its response minus the offset -- is
  # materially different: the offset varies between strata.
  adjusted <- mix(y_adj ~ x + (1 | stratum))
  resid_var <- function(m) m$varDF$vcov[m$varDF$grp == "Residual"][1]
  vpc <- function(m) stratum_var(m) / (stratum_var(m) + resid_var(m))
  expect_equal(vpc(without), 0.7993171, tolerance = 1e-5)
  expect_equal(vpc(adjusted), 0.2279454, tolerance = 1e-5)

  # A linear offset shifts the coefficients by exactly its own: 3 and 2.
  d$y_lin <- d$y_adj + 3 + 2 * d$x
  lin_off <- mix(y_lin ~ x + offset(3 + 2 * x) + (1 | stratum))
  expect_equal(unname(lin_off$coef - adjusted$coef), c(3, 2), tolerance = 1e-6)
})

test_that("the recipe in the refusal fits the offset model", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- b17_data()
  d$w1 <- 1
  msg <- tryCatch({
    b17_quiet(fit_maihda(y ~ x + offset(off) + (1 | stratum), d, engine = "wemix",
                         sampling_weights = "w1"))
    "fitted without an error"
  }, error = conditionMessage)
  expect_match(msg, b17_msg, fixed = TRUE)
  # Run the recipe line exactly as the message prints it.
  recipe <- regmatches(msg, regexpr("data\\$[^\n]+", msg))
  expect_identical(recipe, "data$y_adj <- with(data, y - off)")
  data <- d
  data$y_adj <- NULL
  eval(parse(text = recipe))
  expect_equal(data$y_adj, d$y - d$off)

  # At unit weights the pseudo-likelihood is the likelihood, so lmer(REML = FALSE)
  # with the offset is an independent reference for the recipe's fit (d$y_adj is
  # the column the recipe builds, as checked above).
  fit <- b17_quiet(fit_maihda(y_adj ~ x + (1 | stratum), d, engine = "wemix",
                              sampling_weights = "w1"))
  ref <- lme4::lmer(y ~ x + offset(off) + (1 | stratum), d, REML = FALSE)
  ref_var <- as.data.frame(lme4::VarCorr(ref))$vcov
  expect_equal(unname(fit$model$coef), unname(lme4::fixef(ref)), tolerance = 1e-5)
  expect_equal(maihda_wemix_variances(fit)$stratum, ref_var[1], tolerance = 1e-4)
  expect_equal(maihda_wemix_variances(fit)$residual, ref_var[2], tolerance = 1e-4)
  # Its predictions plus the offset are the offset model's.
  expect_equal(unname(predict_maihda(fit, scale = "link") + d$off),
               unname(stats::predict(ref)), tolerance = 1e-5)
})

test_that("fit_maihda(engine = 'wemix') refuses a formula offset before fitting", {
  d <- b17_dim_data()
  # The refusal must come before the engine: reaching it is a failure.
  local_mocked_bindings(maihda_fit_wemix = function(...) stop("reached the engine"))
  refuse <- function(f, ...) {
    expect_error(
      b17_quiet(fit_maihda(f, d, engine = "wemix", sampling_weights = "w", ...)),
      b17_msg, fixed = TRUE)
  }
  refuse(y ~ x + offset(off) + (1 | a:b))
  refuse(y ~ x + (1 | a:b) + offset(off))
  refuse(y ~ offset(off) + x + (1 | a:b))
  refuse(y ~ 0 + x + offset(off) + (1 | a:b))
  refuse(y ~ offset(off) + (1 | a:b))
  refuse(y ~ x + (offset(o1) + offset(o2)) + (1 | a:b))
  refuse(y ~ x + offset(0.5 * x) + (1 | a:b))
  refuse(yb ~ x + offset(off) + (1 | a:b), family = "binomial")
  # The default engine chosen by sampling_weights refuses too.
  expect_error(b17_quiet(fit_maihda(y ~ x + offset(off) + (1 | a:b), d,
                                    sampling_weights = "w")),
               b17_msg, fixed = TRUE)
  # A pre-built stratum column.
  d2 <- b17_quiet(make_strata(d, c("a", "b")))$data
  expect_error(b17_quiet(fit_maihda(y ~ x + offset(off) + (1 | stratum), d2,
                                    engine = "wemix", sampling_weights = "w")),
               b17_msg, fixed = TRUE)
})

test_that("maihda() and compare_maihda_groups() refuse it once, before any fit", {
  d <- b17_dim_data()
  local_mocked_bindings(fit_maihda = function(...) stop("a fit was attempted"))
  f <- y ~ x + offset(off) + (1 | a:b)
  expect_error(b17_quiet(maihda(f, d, engine = "wemix", sampling_weights = "w")),
               b17_msg, fixed = TRUE)
  expect_error(b17_quiet(maihda(f, d, sampling_weights = "w")), b17_msg, fixed = TRUE)
  expect_error(b17_quiet(maihda(f, d, group = "grp", engine = "wemix",
                                sampling_weights = "w")),
               b17_msg, fixed = TRUE)

  # compare_maihda_groups() turns a failing group fit into a warning and a row, so
  # count the warnings rather than suppress them: the refusal must be one error.
  warns <- character(0)
  err <- tryCatch({
    withCallingHandlers(
      suppressMessages(compare_maihda_groups(f, d, group = "grp", engine = "wemix",
                                             sampling_weights = "w")),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    "returned without an error"
  }, error = conditionMessage)
  expect_match(err, b17_msg, fixed = TRUE)
  expect_length(warns, 0)
})

test_that("the refusal names every offset and a runnable recipe", {
  check <- function(f) tryCatch(maihda_wemix_check_offset(f), error = conditionMessage)
  two <- check(y ~ x + (offset(o1) + offset(o2)) + (1 | stratum))
  expect_match(two, "(offset(o1), offset(o2))", fixed = TRUE)
  expect_match(two, "data$y_adj <- with(data, y - (o1 + o2))", fixed = TRUE)
  expect_match(check(y ~ x + offset(3 + 2 * x) + (1 | stratum)),
               "with(data, y - (3 + 2 * x))", fixed = TRUE)
  expect_match(check(y ~ x + offset(-o1) + (1 | stratum)),
               "with(data, y - (-o1))", fixed = TRUE)
  expect_match(check(log(y) ~ x + offset(log(e)) + (1 | stratum)),
               "data$y_adj <- with(data, log(y) - log(e))", fixed = TRUE)
  expect_match(check(yb ~ offset(off) + (1 | stratum)), "data$yb_adj <-", fixed = TRUE)
  expect_match(check(`my y` ~ offset(`log exposure`) + (1 | stratum)),
               "data$my.y_adj <- with(data, `my y` - `log exposure`)", fixed = TRUE)
  # terms() also reports an offset written after a minus sign, and model.offset()
  # would add it, so it is refused as well.
  expect_match(check(y ~ x + (1 | stratum) - offset(off)), b17_msg, fixed = TRUE)
  expect_match(check(y ~ . + offset(off) + (1 | stratum)), b17_msg, fixed = TRUE)
  expect_match(check(y ~ x + offset() + (1 | stratum)), b17_msg, fixed = TRUE)

  # Each recipe computes what its formula's offset adds.
  d <- b17_dim_data()
  for (f in list(y ~ x + (offset(o1) + offset(o2)) + (1 | stratum),
                 y ~ x + offset(-o1) + (1 | stratum),
                 y ~ x + offset(3 + 2 * x) + (1 | stratum))) {
    recipe <- regmatches(check(f), regexpr("data\\$[^\n]+", check(f)))
    data <- d
    eval(parse(text = recipe))
    off <- stats::model.offset(stats::model.frame(maihda_nobars(f), d))
    expect_equal(data$y_adj, d$y - off)
  }
  # Non-syntactic names keep their backticks, so that recipe runs too.
  d[["my y"]] <- d$y
  d[["log exposure"]] <- d$o1
  f <- `my y` ~ offset(`log exposure`) + (1 | stratum)
  recipe <- regmatches(check(f), regexpr("data\\$[^\n]+", check(f)))
  data <- d
  eval(parse(text = recipe))
  expect_equal(data$my.y_adj, d$y - d$o1)

  # No offset, or only a name that looks like one: nothing to refuse.
  for (f in list(y ~ x + (1 | stratum),
                 y ~ (1 | stratum),
                 y ~ x + I(offset(off)) + (1 | stratum),
                 y ~ x + offset_var + (1 | stratum),
                 y ~ x + log(offset) + (1 | stratum))) {
    expect_silent(maihda_wemix_check_offset(f))
  }
  expect_silent(maihda_wemix_check_offset("y ~ x + offset(off)"))
})

test_that("the other engines still fit a formula offset", {
  d <- b17_dim_data()
  m <- b17_quiet(fit_maihda(y ~ x + offset(off) + (1 | a:b), d))
  ref <- b17_quiet(lme4::lmer(y ~ x + offset(off) + (1 | stratum),
                              transform(d, stratum = m$data$stratum)))
  expect_equal(lme4::fixef(m$model), lme4::fixef(ref))
  skip_if_not_installed("ordinal")
  d$yo <- factor(sample(1:3, nrow(d), TRUE), ordered = TRUE)
  mo <- b17_quiet(fit_maihda(yo ~ x + offset(off) + (1 | a:b), d, family = "ordinal"))
  expect_identical(mo$engine, "ordinal")
})
