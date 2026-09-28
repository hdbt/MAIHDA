# Audit 2026-09-28b: a brms fit_maihda() model whose two-level outcome was recoded to
# 0/1 refused rows that carried the outcome as the user fitted it -- the finding the
# 2026-09-27c pass measured and spawned. PARTIAL: confirmed for a factor, character or
# logical outcome; refuted for a 1/2 coding.
#
# fit_maihda() recodes a two-level outcome to 0/1 before brms sees it
# ($response_recoding), and brms (2.23) checks that every model variable numeric in its
# fitted data is numeric in `newdata`, the outcome included although no prediction
# reads it. So predict_maihda(newdata = ) on either scale and
# plot_prediction_deviation_panels(data = ) stopped with "Variable 'yf' was originally
# numeric but is not in 'newdata'" -- measured on real fits, on HEAD and the tree
# before this fix alike, for a factor, a logical and a sampling-weighted factor
# outcome, and for character labels handed to the factor fit. A 1/2 coding was never
# refused: the check is on type, not values, and its predictions equalled the
# outcome-dropped ones exactly. Nor was an outcome column holding an NA -- brms skips
# the check for such a column (brms:::validate_newdata()) -- so the refusal hit fully
# observed outcomes only. Dropping the outcome changes no prediction (0 difference, and
# the stored frame handed back gives predictions identical() to the fitted rows'):
# brms's prediction methods fill a missing response themselves. Also measured: brms
# 2.23 does not need a weighted fit's .maihda_sw column for any of these predictions
# (it fills a missing weights column with 1), so the panel lacked nothing there.
#
# FIX: maihda_brms_newdata() hands brms a copy of supplied rows without the recoded
# outcome column, for predict_maihda() and the panel alike; the panel keeps the
# original column for its own outcome coding. Fitted-row routes are unchanged.

a28b_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression and the panel's own warnings.
a28b_capture <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value,
       warnings = grep("plot_prediction_deviation_panels", w, value = TRUE))
}

# The value of an expression, or NULL where it errors -- so a refusal fails the
# assertion it feeds rather than aborting the block, and every assertion reports.
a28b_or_null <- function(expr) tryCatch(expr, error = function(e) NULL)

a28b_bern <- function(y, p) ifelse(y == 1, sqrt(-2 * log(p)), sqrt(-2 * log1p(-p)))
a28b_mean_by <- function(x, g) as.numeric(tapply(x, g, mean))
a28b_strata <- function(p) {
  x <- p[[2]]$data
  x[order(x$stratum), , drop = FALSE]
}

# 12 strata x 25 rows; a control/case factor outcome (control declared first) and its
# logical twin. `.row` lets the Stan-free fit below find each row's probability.
a28b_data <- function() {
  set.seed(22)
  e <- data.frame(stratum = rep(1:12, each = 25))
  e$x <- stats::rnorm(nrow(e))
  e$yf <- factor(ifelse(stats::rbinom(nrow(e), 1, stats::plogis(0.6 * e$x +
                   stats::rnorm(12, 0, 0.7)[e$stratum])) == 1, "case", "control"),
                 levels = c("control", "case"))
  e$ylg <- e$yf == "case"
  e$.row <- seq_len(nrow(e))
  e
}

# brms 2.23's newdata check (brms:::validate_newdata(), and as measured on real fits): a
# model variable numeric in the fitted data must be numeric in newdata, the outcome
# included -- unless the newdata column holds an NA, which switches the check off for
# it. No weights() column is required.
a28b_brms_check <- function(object, newdata) {
  for (nm in intersect(names(object$data), names(newdata))) {
    if (is.numeric(object$data[[nm]]) && !anyNA(newdata[[nm]]) &&
        !is.numeric(newdata[[nm]])) {
      stop("Variable '", nm, "' was originally numeric but is not in 'newdata'.",
           call. = FALSE)
    }
  }
}

# A fit_maihda() brms Bernoulli model, Stan-free: its stored data holds the outcome
# recoded to 0/1, as fit_maihda() hands it to brms, and fitted() / posterior_epred() /
# posterior_linpred() -- registered on a private subclass -- apply brms's check and
# return fixed probabilities (and draws) by the .row column.
a28b_mock <- function(e, outcome = "yf", levels = c("control", "case")) {
  cls <- "a28b_mock_brmsfit"
  rows <- function(object, newdata) {
    if (is.null(newdata)) newdata <- object$data
    a28b_brms_check(object, newdata)
    newdata$.row
  }
  registerS3method("fitted", cls, function(object, newdata = NULL, ...) {
    p <- object$p[rows(object, newdata)]
    cbind(Estimate = p, Est.Error = 0.05, Q2.5 = pmax(p - 0.1, 0),
          Q97.5 = pmin(p + 0.1, 1))
  }, envir = asNamespace("stats"))
  registerS3method("posterior_epred", cls, function(object, newdata = NULL, ...) {
    object$draws[, rows(object, newdata), drop = FALSE]
  }, envir = asNamespace("brms"))
  registerS3method("posterior_linpred", cls, function(object, newdata = NULL, ...) {
    stats::qlogis(object$draws[, rows(object, newdata), drop = FALSE])
  }, envir = asNamespace("brms"))
  set.seed(4)
  p <- stats::plogis(-0.2 + 0.6 * e$x + stats::rnorm(12, 0, 0.7)[e$stratum])
  draws <- t(vapply(seq_len(100), function(i) {
    stats::plogis(stats::qlogis(p) + stats::rnorm(1, sd = 0.2))
  }, numeric(length(p))))
  fit_data <- e
  fit_data[[outcome]] <- as.integer(as.character(e[[outcome]]) == levels[2])
  f <- stats::as.formula(paste(outcome, "~ x + (1 | stratum)"))
  model <- structure(
    list(formula = brms::bf(f), data = fit_data, p = p, draws = draws,
         family = structure(list(family = "bernoulli", link = "logit"), class = "family")),
    class = c(cls, "brmsfit"))
  structure(
    list(model = model, engine = "brms", formula = f, data = fit_data, original_data = e,
         family = list(family = "bernoulli", link = "logit"), strata_vars = NULL,
         strata_info = NULL, context_vars = NULL, sampling_weights = NULL,
         longitudinal_info = NULL,
         response_recoding = data.frame(level = levels, value = 0:1,
                                        role = c("reference", "event"),
                                        stringsAsFactors = FALSE)),
    class = "maihda_model")
}

test_that("predict_maihda() hands brms supplied rows without the recoded outcome", {
  skip_if_not_installed("brms")
  e <- a28b_data()
  m <- a28b_mock(e)
  link <- as.numeric(colMeans(stats::qlogis(m$model$draws)))
  pred <- function(obj, nd, ...) {
    a28b_or_null(as.numeric(predict_maihda(obj, newdata = nd, type = "individual", ...)))
  }
  ec <- e
  ec$yf <- as.character(ec$yf)
  for (nd in list(e, ec)) {
    expect_equal(pred(m, nd), m$model$p, tolerance = 1e-12)
    expect_equal(pred(m, nd, scale = "link"), link, tolerance = 1e-12)
  }
  ml <- a28b_mock(e, outcome = "ylg", levels = c("FALSE", "TRUE"))
  expect_equal(pred(ml, e), ml$model$p, tolerance = 1e-12)
})

test_that("the deviation panel predicts brms rows and codes their outcome as fitted", {
  skip_if_not_installed("brms")
  e <- a28b_data()
  m <- a28b_mock(e)
  y <- as.integer(e$yf == "case")
  r <- a28b_or_null(a28b_capture(plot_prediction_deviation_panels(m, data = e,
                                                                  type = "binomial")))
  expect_false(is.null(r))
  expect_identical(r$warnings, character(0))
  x <- if (!is.null(r)) a28b_strata(r$value)
  expect_equal(x$fitted, a28b_mean_by(m$model$p, e$stratum), tolerance = 1e-12)
  expect_equal(x$abs_res_dev, a28b_mean_by(a28b_bern(y, m$model$p), e$stratum),
               tolerance = 1e-12)
})

test_that("controls: the model's own rows, its stored frame, rows missing an outcome", {
  # Right before the fix too: these rows carry the outcome already recoded to 0/1, or
  # (the last) an outcome column brms does not check.
  skip_if_not_installed("brms")
  e <- a28b_data()
  m <- a28b_mock(e)
  y <- as.integer(e$yf == "case")
  expect_equal(as.numeric(predict_maihda(m, type = "individual")), m$model$p,
               tolerance = 1e-12)
  expect_equal(as.numeric(predict_maihda(m, newdata = m$data, type = "individual")),
               m$model$p, tolerance = 1e-12)
  d <- a28b_strata(a28b_quiet(plot_prediction_deviation_panels(m, type = "binomial")))
  expect_equal(d$fitted, a28b_mean_by(m$model$p, e$stratum), tolerance = 1e-12)
  expect_equal(d$abs_res_dev, a28b_mean_by(a28b_bern(y, m$model$p), e$stratum),
               tolerance = 1e-12)
  # brms does not check a column holding an NA, so rows whose outcome has a missing
  # value were never refused; they are still predicted in full.
  ena <- e
  ena$yf[c(3, 40)] <- NA
  expect_equal(as.numeric(predict_maihda(m, newdata = ena, type = "individual")),
               m$model$p, tolerance = 1e-12)
})

test_that("only the recoded outcome, and only on supplied rows, is dropped", {
  skip_if_not_installed("brms")
  e <- a28b_data()
  m <- a28b_mock(e)
  expect_false("yf" %in% names(maihda_brms_newdata(m, e)))
  expect_identical(maihda_brms_newdata(m, e, drop_outcome = FALSE), e)
  # An outcome the fit did not recode, and a model that is not a fit_maihda() one, are
  # left alone.
  m0 <- m
  m0$response_recoding <- NULL
  expect_identical(maihda_brms_newdata(m0, e), e)
  expect_identical(maihda_brms_newdata(m$model, e), e)
  expect_identical(maihda_brms_newdata(m, e[, c("x", "stratum", ".row")]),
                   e[, c("x", "stratum", ".row")])
})

test_that("real brms fits predict rows carrying the outcome as fitted", {
  # Compiles one Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). Every reference is
  # the same fit's own predictions, so sampler quality is irrelevant.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  e <- a28b_data()
  e$.row <- NULL
  m <- a28b_quiet(fit_maihda(yf ~ x + (1 | stratum), data = e, engine = "brms",
                             family = "binomial", chains = 1, iter = 600, warmup = 300,
                             refresh = 0, seed = 1))
  expect_true(is.numeric(m$data$yf))
  pred <- function(...) as.numeric(a28b_quiet(predict_maihda(m, type = "individual", ...)))
  ref <- pred(newdata = e[, c("x", "stratum")])
  ec <- e
  ec$yf <- as.character(ec$yf)
  # A missing outcome switches brms's check off for the column: accepted before too.
  ena <- e
  ena$yf[c(3, 40)] <- NA
  for (nd in list(e, ec, ena)) {
    expect_equal(pred(newdata = nd), ref, tolerance = 1e-12)
  }
  expect_equal(pred(), ref, tolerance = 1e-12)
  d <- a28b_strata(a28b_quiet(plot_prediction_deviation_panels(m)))
  r <- a28b_capture(plot_prediction_deviation_panels(m, data = e))
  expect_identical(r$warnings, character(0))
  x <- a28b_strata(r$value)
  expect_equal(x$fitted, d$fitted, tolerance = 1e-12)
  expect_equal(x$abs_res_dev, d$abs_res_dev, tolerance = 1e-12)
})
