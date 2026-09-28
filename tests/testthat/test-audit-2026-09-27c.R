# Audit 2026-09-27c: the binomial prediction-deviation panel codes a supplied outcome
# as the MODEL coded it -- a finding the 2026-09-27b pass measured and spawned.
#
# The panel took each supplied row's 0/1 outcome from maihda_binomial_observed_01(),
# which re-derives the coding from the plotted rows themselves: a factor's own level
# order, sorted values otherwise. So (1) the fitted rows handed back with their factor
# levels reversed had every residual scored against the wrong outcome and every
# "Wrong" / "Correct" swapped -- 1.54 off on a row of a glm and all 200 of its
# classifications flipped, 0.65 on a stratum mean of glmer and of a fit_maihda() model,
# 0.72 on a wemix one, identically on HEAD; and (2) rows holding only one of the two
# values could not be coded at all: scored 0 on HEAD, a perfect fit, and left unscored
# with a warning since the 2026-09-27b pass.
#
# FIX: maihda_prediction_panel_fit_01() keys the coding on the fit -- a fit_maihda()
# model's $response_recoding (original label -> 0/1), or the fitted response itself (a
# factor's first level is the failure, as in R's binomial family and in brms; TRUE is
# the event; 0/1 numbers are themselves; for a wemix model, its outcome expression
# evaluated on its stored rows) -- for supplied data only; the fitted rows already
# carry the model's coding and are read as before. A numeric 0/1 column that can only
# be the recoded one handed back (object$data) is taken as it stands: its values are
# not all original labels, or the original cannot have been numeric (its recorded
# levels are not ascending numbers -- my first rule missed this and inverted
# object$data for factor(y, levels = c(1, 0)), pinned below). A label the fit never
# saw is left out with a warning.

a27c_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression and the panel's own warnings.
a27c_capture <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value,
       warnings = grep("plot_prediction_deviation_panels", w, value = TRUE))
}

# |Bernoulli deviance residual| of a 0/1 outcome at probability p, by hand.
a27c_bern <- function(y, p) ifelse(y == 1, sqrt(-2 * log(p)), sqrt(-2 * log1p(-p)))

# Case-level panel rows in the order of `data`; stratum-level rows in stratum order.
a27c_cases <- function(p) {
  x <- p[[2]]$data
  x[order(x$id), , drop = FALSE]
}
a27c_strata <- function(p) {
  x <- p[[2]]$data
  x[order(x$stratum), , drop = FALSE]
}
a27c_mean_by <- function(x, g) as.numeric(tapply(x, g, mean))

# 200 rows with a No/Yes factor outcome (No declared first, so No is the failure).
a27c_glm_data <- function() {
  set.seed(21)
  d <- data.frame(x = stats::rnorm(200))
  d$yf <- factor(ifelse(stats::rbinom(200, 1, stats::plogis(0.8 * d$x)) == 1, "Yes", "No"),
                 levels = c("No", "Yes"))
  d
}

# 12 strata x 25 rows with a control/case factor outcome (control declared first).
a27c_strata_data <- function() {
  set.seed(22)
  e <- data.frame(stratum = rep(1:12, each = 25))
  e$x <- stats::rnorm(nrow(e))
  e$yf <- factor(ifelse(stats::rbinom(nrow(e), 1, stats::plogis(0.6 * e$x +
                   stats::rnorm(12, 0, 0.7)[e$stratum])) == 1, "case", "control"),
                 levels = c("control", "case"))
  e
}

a27c_reverse <- function(d) {
  d$yf <- factor(as.character(d$yf), levels = rev(levels(d$yf)))
  d
}

test_that("a glm's rows with their factor levels reversed are scored as fitted", {
  skip_on_cran()
  d <- a27c_glm_data()
  g <- stats::glm(yf ~ x, data = d, family = stats::binomial())
  p <- as.numeric(stats::fitted(g))
  y <- as.integer(d$yf == "Yes")
  r <- a27c_capture(plot_prediction_deviation_panels(g, data = a27c_reverse(d)[, c("x", "yf")]))
  expect_identical(r$warnings, character(0))
  x <- a27c_cases(r$value)
  expect_equal(x$abs_res_dev, a27c_bern(y, p), tolerance = 1e-12)
  expect_identical(as.character(x$wrong),
                   ifelse((p > 0.5 & y == 0) | (p < 0.5 & y == 1), "Wrong", "Correct"))
  # Control (right on HEAD too): the labels supplied as character.
  dc <- d
  dc$yf <- as.character(dc$yf)
  xc <- a27c_cases(a27c_quiet(plot_prediction_deviation_panels(g, data = dc[, c("x", "yf")])))
  expect_equal(xc$abs_res_dev, a27c_bern(y, p), tolerance = 1e-12)
})

test_that("rows holding a single value are scored, not left out", {
  skip_on_cran()
  d <- a27c_glm_data()
  g <- stats::glm(yf ~ x, data = d, family = stats::binomial())
  p <- as.numeric(stats::fitted(g))
  yes <- which(d$yf == "Yes")
  r <- a27c_capture(plot_prediction_deviation_panels(g, data = d[yes, c("x", "yf")]))
  expect_identical(r$warnings, character(0))
  expect_equal(a27c_cases(r$value)$abs_res_dev, a27c_bern(rep(1, length(yes)), p[yes]),
               tolerance = 1e-12)
})

test_that("a raw fit's logical or 0/1 response codes single-valued rows too", {
  skip_on_cran()
  d <- a27c_glm_data()
  d$yl <- d$yf == "Yes"
  d$y01 <- as.integer(d$yl)
  yes <- which(d$yl)
  for (f in list(yl ~ x, y01 ~ x, I(y01 > 0) ~ x)) {
    g <- stats::glm(f, data = d, family = stats::binomial())
    p <- as.numeric(stats::fitted(g))
    r <- a27c_capture(plot_prediction_deviation_panels(g, data = d[yes, ]))
    lab <- paste(deparse(f), collapse = " ")
    expect_identical(r$warnings, character(0), label = lab)
    expect_equal(a27c_cases(r$value)$abs_res_dev, a27c_bern(rep(1, length(yes)), p[yes]),
                 tolerance = 1e-12, label = lab)
  }
})

test_that("glmer and fit_maihda() models code supplied rows as they were fitted", {
  skip_on_cran()
  e <- a27c_strata_data()
  y <- as.integer(e$yf == "case")
  g <- a27c_quiet(lme4::glmer(yf ~ x + (1 | stratum), data = e, family = stats::binomial()))
  pg <- as.numeric(stats::fitted(g))
  expect_equal(a27c_strata(a27c_quiet(plot_prediction_deviation_panels(g, data = a27c_reverse(e))))$abs_res_dev,
               a27c_mean_by(a27c_bern(y, pg), e$stratum), tolerance = 1e-12)

  m <- a27c_quiet(fit_maihda(yf ~ x + (1 | stratum), data = e, family = "binomial"))
  expect_identical(m$response_recoding$level[m$response_recoding$value == 1], "case")
  pm <- as.numeric(stats::fitted(m$model))
  expect_equal(a27c_strata(a27c_quiet(plot_prediction_deviation_panels(m, data = a27c_reverse(e))))$abs_res_dev,
               a27c_mean_by(a27c_bern(y, pm), e$stratum), tolerance = 1e-12)
  ic <- which(e$yf == "case")
  r <- a27c_capture(plot_prediction_deviation_panels(m, data = e[ic, ]))
  expect_identical(r$warnings, character(0))
  expect_equal(a27c_strata(r$value)$abs_res_dev,
               a27c_mean_by(a27c_bern(rep(1, length(ic)), pm[ic]), e$stratum[ic]),
               tolerance = 1e-12)
  # Control (right on HEAD too): the stored, recoded frame handed back.
  expect_equal(a27c_strata(a27c_quiet(plot_prediction_deviation_panels(m, data = m$data)))$abs_res_dev,
               a27c_strata(a27c_quiet(plot_prediction_deviation_panels(m)))$abs_res_dev,
               tolerance = 1e-12)
})

test_that("a 1/2-coded outcome is read in either coding, one value or two", {
  skip_on_cran()
  e <- a27c_strata_data()
  e$y12 <- ifelse(e$yf == "case", 2, 1)
  m <- a27c_quiet(fit_maihda(y12 ~ x + (1 | stratum), data = e, family = "binomial"))
  p <- as.numeric(stats::fitted(m$model))
  y <- as.integer(e$y12 == 2)
  panel <- function(dd) a27c_strata(a27c_quiet(plot_prediction_deviation_panels(m, data = dd)))$abs_res_dev
  # Control (right on HEAD too): the original data, both values present.
  expect_equal(panel(e), a27c_mean_by(a27c_bern(y, p), e$stratum), tolerance = 1e-12)
  i2 <- which(e$y12 == 2)
  i1 <- which(e$y12 == 1)
  expect_equal(panel(e[i2, ]), a27c_mean_by(a27c_bern(rep(1, length(i2)), p[i2]), e$stratum[i2]),
               tolerance = 1e-12)
  expect_equal(panel(e[i1, ]), a27c_mean_by(a27c_bern(rep(0, length(i1)), p[i1]), e$stratum[i1]),
               tolerance = 1e-12)
  # The recoded 0/1 frame handed back, whole and as its reference rows.
  expect_equal(panel(m$data), a27c_mean_by(a27c_bern(y, p), e$stratum), tolerance = 1e-12)
  expect_equal(panel(m$data[m$data$y12 == 0, ]),
               a27c_mean_by(a27c_bern(rep(0, length(i1)), p[i1]), e$stratum[i1]),
               tolerance = 1e-12)
})

test_that("word labels and logical outcomes holding one value are scored", {
  skip_on_cran()
  e <- a27c_strata_data()
  ic <- which(e$yf == "case")
  e$ych <- ifelse(e$yf == "case", "yes", "no")
  e$ylg <- e$yf == "case"
  for (resp in c("ych", "ylg")) {
    m <- a27c_quiet(fit_maihda(stats::reformulate(c("x", "(1 | stratum)"), response = resp),
                               data = e, family = "binomial"))
    p <- as.numeric(stats::fitted(m$model))
    r <- a27c_capture(plot_prediction_deviation_panels(m, data = e[ic, ]))
    expect_identical(r$warnings, character(0), label = resp)
    expect_equal(a27c_strata(r$value)$abs_res_dev,
                 a27c_mean_by(a27c_bern(rep(1, length(ic)), p[ic]), e$stratum[ic]),
                 tolerance = 1e-12, label = resp)
  }
})

test_that("a wemix binomial fit codes supplied rows as it was fitted", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  e <- a27c_strata_data()
  set.seed(23)
  e$sw <- stats::runif(nrow(e), 0.5, 2)
  m <- a27c_quiet(fit_maihda(yf ~ x + (1 | stratum), data = e, sampling_weights = "sw",
                             family = stats::binomial()))
  pd <- a27c_strata(a27c_quiet(plot_prediction_deviation_panels(m)))
  pr <- a27c_strata(a27c_quiet(plot_prediction_deviation_panels(m, data = a27c_reverse(e))))
  expect_equal(pr$abs_res_dev, pd$abs_res_dev, tolerance = 1e-12)
  expect_equal(pr$fitted, pd$fitted, tolerance = 1e-12)

  # An outcome written as an expression is not recoded, so it carries no
  # $response_recoding, and the wemix stored rows are plain data: the fitted response is
  # the expression evaluated on them. Rows holding only events are scored.
  set.seed(31)
  e$ly <- stats::plogis(0.8 * e$x + stats::rnorm(nrow(e)))
  mx <- a27c_quiet(fit_maihda(I(ly > 0.5) ~ x + (1 | stratum), data = e,
                              sampling_weights = "sw", family = stats::binomial()))
  expect_null(mx$response_recoding)
  ev <- which(e$ly > 0.5)
  p <- as.numeric(maihda_prediction_panel_fitted(mx$model, e, "binomial", maihda_obj = mx)$fit)
  r <- a27c_capture(plot_prediction_deviation_panels(mx, data = e[ev, ]))
  expect_identical(r$warnings, character(0))
  ref <- vapply(split(seq_along(ev), e$stratum[ev]), function(i)
    sum(a27c_bern(rep(1, length(i)), p[ev][i]) * e$sw[ev][i]) / sum(e$sw[ev][i]), numeric(1))
  expect_equal(a27c_strata(r$value)$abs_res_dev, as.numeric(ref), tolerance = 1e-12)
})

test_that("the coding rule, keyed on the fit", {
  rec <- function(levels) {
    list(response_recoding = data.frame(level = levels, value = 0:1,
                                        role = c("reference", "event"),
                                        stringsAsFactors = FALSE))
  }
  code <- function(obj, x) as.vector(maihda_prediction_panel_fit_01(obj, NULL, x))
  fac <- rec(c("control", "case"))
  x <- factor(c("case", "control", NA, "case"), levels = c("case", "control"))
  expect_identical(code(fac, x), c(1L, 0L, NA, 1L))
  # A label the fit never saw is uncoded, and reported; a missing value is neither.
  other <- maihda_prediction_panel_fit_01(fac, NULL, c("case", "other", NA))
  expect_identical(as.vector(other), c(1L, NA, NA))
  expect_identical(attr(other, "unmatched"), 2L)
  # The recoded column handed back: numeric, all 0/1, not original labels.
  expect_identical(code(fac, c(0L, 1L, 1L)), c(0L, 1L, 1L))
  # A 1/2 coding: the original values, the recoded ones, and one value of each.
  r12 <- rec(c("1", "2"))
  expect_identical(code(r12, c(1, 2, NA)), c(0L, 1L, NA))
  expect_identical(code(r12, c(0, 1)), c(0L, 1L))
  expect_identical(code(r12, c(2, 2)), c(1L, 1L))
  expect_identical(code(r12, c(0, 0)), c(0L, 0L))
  # Only 1s is both codings at once; it is read as the original data `data` is
  # documented to be.
  expect_identical(code(r12, c(1, 1)), c(0L, 0L))
  lgl <- rec(c("FALSE", "TRUE"))
  expect_identical(code(lgl, c(TRUE, TRUE)), c(1L, 1L))
  # Labels "1"/"0" in that order can only be a factor's (a numeric original is
  # recorded ascending): its labels are read as labels, and numbers are the recoded
  # column.
  r10 <- rec(c("1", "0"))
  expect_identical(code(r10, factor(c("0", "1"), levels = c("1", "0"))), c(1L, 0L))
  expect_identical(code(r10, c(0L, 1L, 1L)), c(0L, 1L, 1L))
  # TRUE / FALSE for an outcome coded 0/1 are 1 / 0 (the old derivation's reading too);
  # for word labels they are no label at all.
  r01 <- rec(c("0", "1"))
  expect_identical(code(r01, c(TRUE, FALSE, NA)), c(1L, 0L, NA))
  expect_identical(code(fac, c(TRUE, FALSE)), c(NA_integer_, NA_integer_))
})

test_that("a factor labelled 0/1 in reverse order reads both codings", {
  skip_on_cran()
  e <- a27c_strata_data()
  e$f10 <- factor(ifelse(e$yf == "case", "0", "1"), levels = c("1", "0"))
  m <- a27c_quiet(fit_maihda(f10 ~ x + (1 | stratum), data = e, family = "binomial"))
  expect_identical(m$response_recoding$level[m$response_recoding$value == 1], "0")
  p <- as.numeric(stats::fitted(m$model))
  ref <- a27c_mean_by(a27c_bern(as.integer(e$f10 == "0"), p), e$stratum)
  panel <- function(dd) a27c_strata(a27c_quiet(plot_prediction_deviation_panels(m, data = dd)))$abs_res_dev
  expect_equal(panel(e), ref, tolerance = 1e-12)
  # The stored, recoded frame: its 0/1 are the recoded values, not the labels "0"/"1".
  expect_equal(panel(m$data), ref, tolerance = 1e-12)
})

test_that("a 0/1 outcome supplied as logical is read as 1 / 0", {
  # Right on HEAD too; my first keyed rule matched TRUE / FALSE against the recorded
  # labels "0" / "1" and left every row unscored.
  skip_on_cran()
  e <- a27c_strata_data()
  e$y01 <- as.integer(e$yf == "case")
  m <- a27c_quiet(fit_maihda(y01 ~ x + (1 | stratum), data = e, family = "binomial"))
  p <- as.numeric(stats::fitted(m$model))
  el <- e
  el$y01 <- el$y01 == 1
  r <- a27c_capture(plot_prediction_deviation_panels(m, data = el))
  expect_identical(r$warnings, character(0))
  expect_equal(a27c_strata(r$value)$abs_res_dev, a27c_mean_by(a27c_bern(e$y01, p), e$stratum),
               tolerance = 1e-12)
})

test_that("an outcome the fit never saw is left out, and said to be", {
  skip_on_cran()
  d <- a27c_glm_data()
  g <- stats::glm(yf ~ x, data = d, family = stats::binomial())
  p <- as.numeric(stats::fitted(g))
  y <- as.integer(d$yf == "Yes")
  dm <- d[, c("x", "yf")]
  dm$yf <- as.character(dm$yf)
  dm$yf[1:3] <- "Maybe"
  r <- a27c_capture(plot_prediction_deviation_panels(g, data = dm))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "3 row(s) of `data` is not one the model was fitted on (e.g. 'Maybe')",
               fixed = TRUE)
  x <- a27c_cases(r$value)
  expect_true(all(is.na(x$abs_res_dev[1:3])))
  expect_equal(x$abs_res_dev[-(1:3)], a27c_bern(y, p)[-(1:3)], tolerance = 1e-12)
  # Every label unknown (the right words in the wrong case): one warning, naming why.
  dl <- d[, c("x", "yf")]
  dl$yf <- tolower(as.character(dl$yf))
  r <- a27c_capture(plot_prediction_deviation_panels(g, data = dl))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "no row of `data` can be scored", fixed = TRUE)
  expect_match(r$warnings, "is not one the model was fitted on", fixed = TRUE)

  # A weighted fit still scores a proportion as successes out of its trials, so a
  # value its 0/1 outcome never took is not reported as unscorable.
  d$y01 <- y
  d$w <- 10
  gw <- stats::glm(y01 ~ x, data = d, family = stats::binomial(), weights = w)
  pw <- as.numeric(stats::fitted(gw))
  dp <- d[, c("x", "y01", "w")]
  dp$y01[1] <- 0.4
  r <- a27c_capture(plot_prediction_deviation_panels(gw, data = dp))
  expect_identical(r$warnings, character(0))
  x <- a27c_cases(r$value)
  expect_equal(x$abs_res_dev[1],
               sqrt(2 * (4 * log(4 / (10 * pw[1])) + 6 * log(6 / (10 * (1 - pw[1]))))),
               tolerance = 1e-12)
  expect_equal(x$abs_res_dev[-1], sqrt(10) * a27c_bern(y, pw)[-1], tolerance = 1e-12)
})

test_that("real brms fits code supplied rows as they were fitted", {
  # Compiles two Stan models, so OPT-IN (set MAIHDA_TEST_BRMS=true). The references are
  # the same fits' own panels and draws, so sampler quality is irrelevant.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  e <- a27c_strata_data()
  # The wrapper hands brms the RECODED 0/1 column, and brms's own validate_newdata()
  # then refuses rows carrying the original factor (a separate, pre-existing defect),
  # so the wrapper is exercised on its stored rows: those holding only events, which
  # the plotted rows' own values could not code.
  m <- a27c_quiet(fit_maihda(yf ~ x + (1 | stratum), data = e, engine = "brms",
                             family = "binomial", chains = 1, iter = 600, warmup = 300,
                             refresh = 0, seed = 1))
  pm <- as.numeric(stats::fitted(m$model, newdata = m$data, summary = TRUE)[, "Estimate"])
  ev <- which(m$data$yf == 1)
  r <- a27c_capture(plot_prediction_deviation_panels(m, data = m$data[ev, ]))
  expect_identical(r$warnings, character(0))
  expect_equal(a27c_strata(r$value)$abs_res_dev,
               a27c_mean_by(a27c_bern(rep(1, length(ev)), pm[ev]), m$data$stratum[ev]),
               tolerance = 1e-10)
  # A bare brmsfit keeps the factor, coded by position in the fit's levels.
  rb <- a27c_quiet(brms::brm(yf ~ x + (1 | stratum), data = e, family = brms::bernoulli(),
                             chains = 1, iter = 600, warmup = 300, refresh = 0, seed = 1))
  p <- as.numeric(stats::fitted(rb, newdata = e, summary = TRUE)[, "Estimate"])
  y <- as.integer(e$yf == "case")
  pb <- a27c_strata(a27c_quiet(plot_prediction_deviation_panels(rb, data = a27c_reverse(e),
                                                                type = "binomial")))
  expect_equal(pb$abs_res_dev, a27c_mean_by(a27c_bern(y, p), e$stratum), tolerance = 1e-10)
  ic <- which(e$yf == "case")
  ps <- a27c_strata(a27c_quiet(plot_prediction_deviation_panels(rb, data = e[ic, ],
                                                                type = "binomial")))
  expect_equal(ps$abs_res_dev, a27c_mean_by(a27c_bern(rep(1, length(ic)), p[ic]), e$stratum[ic]),
               tolerance = 1e-10)
})
