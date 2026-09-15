# Audit 2026-09-15 (second pass), finding 5: the ordinal surprise panel matched each
# row's observed category LABEL against the probability columns' NAMES -- CONFIRMED,
# and wider.
#
# plot_prediction_deviation_panels(ordinal_mode = "surprise") read P(observed
# category) with match(label, colnames(probs)). The rebuilt clmm matrix names its
# columns "1".."K", so only an outcome coded 1..K was scored right: low < mid < high
# matched nothing (every stratum surprise NaN, an empty panel), 0 < 1 < 2 lost
# category "0" and scored "1" and "2" against the category before (stratum surprise
# 0.396 off on the auditor's data), and 3 < 2 < 1 was complete and silently wrong
# (1.785 off). Wider than reported: formula() of a brmsfit is a brmsformula, whose
# as.character()[2] is "list()", so a brms fit never found its response -- surprise
# undefined for every coding, and in the binomial branch, which reads the response the
# same way, every brms stratum's deviance residual was 0 and the labelled strata
# arbitrary; brms category probabilities were tabulated from simulated responses and
# changed between calls; and a label used as a column key was overwritten by the
# panel's own columns (a polr category "n" drawn as the stratum row count, "weight"
# as the prior weight).
#
# FIX: the observed category is located by its POSITION among the model's fitted
# categories (clmm y.levels, the brms stored response levels, polr's own column
# names); the probability columns are keyed prob_1..prob_K and the labels only name
# the legend; the response name unwraps a brmsformula and its weights() / thres()
# additions (a trials() count stays unobserved); brms probabilities are fitted()'s
# posterior means; and rows that cannot be scored are warned about. lme4 / glm
# binomial, Gaussian, Poisson and the clmm and polr expected-score panels were
# identical() to the old ones; the surprise values of 1..K clmm and polr fits were
# identical, only their probability column keys renamed. The brms expected-score
# panel changed by design (posterior means instead of draw tabulations).
#
# Why the suite missed it: the one test that drew a clmm surprise panel
# (test-ordinal-engine.R) used labels 1:4 and only checked that it did not error.

d15b_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression and the warnings it raised, muffled.
d15b_capture <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

# Only this function's own warnings: the assertions do not demand silence from code
# the fix does not touch.
d15b_panel_warnings <- function(w) {
  grep("plot_prediction_deviation_panels()", w, fixed = TRUE, value = TRUE)
}

# 12 strata x 25 rows with a 3-category outcome. The categories depend on the seed
# alone, so every coding below is the same outcome under different labels.
d15b_data <- function(labels) {
  set.seed(915)
  d <- data.frame(stratum = factor(rep(seq_len(12), each = 25)))
  d$x <- stats::rnorm(nrow(d))
  lat <- 0.5 * d$x + stats::rnorm(12, sd = 0.7)[d$stratum] + stats::rlogis(nrow(d))
  k <- as.integer(cut(lat, c(-Inf, -0.8, 0.8, Inf)))
  d$yo <- factor(labels[k], levels = labels, ordered = TRUE)
  d
}

d15b_fit <- function(labels) {
  d15b_quiet(fit_maihda(yo ~ x + (1 | stratum), data = d15b_data(labels),
                        engine = "ordinal", family = "ordinal"))
}

# Reference stratum surprise: the mean over each stratum's rows of -log P(observed
# category), P being fitted() of the clmm -- ordinal's own probability of each row's
# observed category, with no package code involved.
d15b_oracle <- function(fit, rows = seq_len(nrow(fit$data))) {
  tapply(-log(stats::fitted(fit$model))[rows],
         as.character(fit$data$stratum[rows]), mean)
}

# The surprise panel's per-stratum values, named by stratum.
d15b_surprise <- function(p) {
  pd <- p[[2]]$data
  stats::setNames(pd$surprise, as.character(pd$stratum))
}

# The strata the surprise panel labels, i.e. its five most surprising.
d15b_labelled <- function(p) {
  lay <- Filter(function(l) inherits(l$geom, "GeomLabelRepel"), p[[2]]$layers)
  sort(as.character(lay[[1]]$data$stratum))
}

test_that("the clmm surprise panel scores the observed category under any response labels", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # low < mid < high matched no column; 0 < 1 < 2 was scored one category off;
  # 3 < 2 < 1 was complete and silently wrong. 1 < 2 < 3, the one coding the label
  # match got right, is the control.
  for (labels in list(c("low", "mid", "high"), c("0", "1", "2"), c("3", "2", "1"),
                      c("1", "2", "3"))) {
    tag <- paste(labels, collapse = " < ")
    fit <- d15b_fit(labels)
    p <- d15b_quiet(plot_prediction_deviation_panels(fit))
    got <- d15b_surprise(p)
    ref <- d15b_oracle(fit)
    expect_length(got, 12L)
    expect_true(all(is.finite(got)), label = tag)
    expect_equal(as.numeric(got[names(ref)]), as.numeric(ref), tolerance = 1e-10, label = tag)
    expect_identical(d15b_labelled(p), sort(names(sort(ref, decreasing = TRUE))[1:5]),
                     label = tag)
    # The legend names the categories, not the column positions.
    expect_identical(levels(p[[1]]$data$Category), labels, label = tag)
  }
})

test_that("prediction data that lacks a category or re-declares its levels is scored by the fitted categories", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  fit <- d15b_fit(c("low", "mid", "high"))
  keep <- which(fit$data$yo != "mid")
  ref <- d15b_oracle(fit, keep)
  sub <- fit$data[keep, ]
  dropped <- sub
  dropped$yo <- droplevels(dropped$yo)          # "high" is now the 2nd level
  relevelled <- sub
  relevelled$yo <- factor(as.character(sub$yo), levels = c("high", "low", "mid"))
  chr <- sub
  chr$yo <- as.character(sub$yo)
  for (nd in list(sub, dropped, relevelled, chr)) {
    r <- d15b_capture(plot_prediction_deviation_panels(fit, data = nd))
    got <- d15b_surprise(r$value)
    expect_equal(as.numeric(got[names(ref)]), as.numeric(ref), tolerance = 1e-10)
    expect_identical(d15b_panel_warnings(r$warnings), character(0))
  }
})

test_that("a category named like a panel column keeps its own probabilities", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  # As column keys, the polr category "n" was overwritten by each stratum's row
  # count and "weight" by the prior weight.
  for (labels in list(c("y", "n", "m"), c("low", "weight", "high"))) {
    d <- d15b_data(labels)
    m <- MASS::polr(yo ~ x, data = d, Hess = TRUE)
    pr <- stats::predict(m, newdata = d, type = "probs")
    p <- plot_prediction_deviation_panels(m, d, type = "ordinal")
    long <- p[[1]]$data
    for (k in seq_along(labels)) {
      rows <- long$Category == labels[k]
      ref <- tapply(pr[, k], as.character(d$stratum), mean)
      expect_equal(long$Probability[rows], as.numeric(ref[as.character(long$stratum[rows])]),
                   tolerance = 1e-10, label = paste0("P(", labels[k], ")"))
    }
    p_obs <- pr[cbind(seq_len(nrow(d)), as.integer(d$yo))]
    ref_s <- tapply(-log(p_obs), as.character(d$stratum), mean)
    got <- d15b_surprise(p)
    expect_equal(as.numeric(got[names(ref_s)]), as.numeric(ref_s), tolerance = 1e-10)
  }
})

test_that("rows whose observed category cannot be scored are warned about", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  fit <- d15b_fit(c("low", "mid", "high"))
  n <- nrow(fit$data)

  # Part of the response recoded to numbers: those rows are counted in a warning, and
  # every other row is still scored exactly.
  recoded <- fit$data
  recoded$yo <- as.character(recoded$yo)
  bad <- seq(2, n, by = 5)
  recoded$yo[bad] <- "2"
  r <- d15b_capture(plot_prediction_deviation_panels(fit, data = recoded))
  w <- d15b_panel_warnings(r$warnings)
  expect_length(w, 1L)
  expect_match(w, sprintf(paste0("%d row(s) of `data` have an observed category that is ",
                                 "not one of the model's fitted categories (low, mid, high)"),
                          length(bad)), fixed = TRUE)
  ref <- d15b_oracle(fit, setdiff(seq_len(n), bad))
  got <- d15b_surprise(r$value)
  expect_equal(as.numeric(got[names(ref)]), as.numeric(ref), tolerance = 1e-10)

  # Nothing to score: the whole response recoded, or no response column at all.
  all_recoded <- fit$data
  all_recoded$yo <- as.integer(all_recoded$yo)
  w <- d15b_panel_warnings(d15b_capture(
    plot_prediction_deviation_panels(fit, data = all_recoded))$warnings)
  expect_length(w, 1L)
  expect_match(w, "no row can be scored", fixed = TRUE)
  expect_match(w, "fitted categories (low, mid, high)", fixed = TRUE)
  no_resp <- fit$data
  no_resp$yo <- NULL
  w <- d15b_panel_warnings(d15b_capture(
    plot_prediction_deviation_panels(fit, data = no_resp))$warnings)
  expect_length(w, 1L)
  expect_match(w, "no row can be scored", fixed = TRUE)
  expect_match(w, "no column 'yo'", fixed = TRUE)

  # A missing response is left out without a warning, and the expected-score mode
  # needs no response at all.
  missing <- fit$data
  gone <- seq(1, n, by = 7)
  missing$yo[gone] <- NA
  r <- d15b_capture(plot_prediction_deviation_panels(fit, data = missing))
  expect_identical(d15b_panel_warnings(r$warnings), character(0))
  ref <- d15b_oracle(fit, setdiff(seq_len(n), gone))
  got <- d15b_surprise(r$value)
  expect_equal(as.numeric(got[names(ref)]), as.numeric(ref), tolerance = 1e-10)
  r <- d15b_capture(plot_prediction_deviation_panels(fit, data = no_resp,
                                                     ordinal_mode = "expected_score"))
  expect_identical(d15b_panel_warnings(r$warnings), character(0))
})

test_that("the observed probability is read by category position, not column name", {
  # Columns named by position, as the rebuilt clmm matrix names them.
  pm <- matrix(c(0.2, 0.5, 0.3,
                 0.6, 0.3, 0.1,
                 0.1, 0.1, 0.8,
                 0.3, 0.3, 0.4), ncol = 3, byrow = TRUE,
               dimnames = list(NULL, c("1", "2", "3")))
  cats <- c("0", "1", "2")
  expect_identical(maihda_prediction_panel_observed_prob(pm, c("0", "2", "1", NA), cats),
                   c(0.2, 0.1, 0.1, NA))
  expect_warning(
    out <- maihda_prediction_panel_observed_prob(pm, c("0", "3", "1", "2"), cats),
    "1 row(s) of `data` have an observed category that is not one of the model's fitted categories (0, 1, 2), e.g. '3'",
    fixed = TRUE)
  expect_identical(out, c(0.2, NA, 0.1, 0.4))
  expect_warning(maihda_prediction_panel_observed_prob(pm, rep(NA, 4), cats),
                 "every observed response in `data` is missing", fixed = TRUE)
  expect_warning(maihda_prediction_panel_observed_prob(pm, rep(NA, 4), cats,
                                                       resp_name = "y", resp_found = FALSE),
                 "`data` has no column 'y'", fixed = TRUE)
  expect_warning(maihda_prediction_panel_observed_prob(pm, rep(NA, 4), cats,
                                                       resp_found = FALSE),
                 "the model's response could not be identified", fixed = TRUE)
})

test_that("the fitted categories come from the fit, one per probability column", {
  probs <- data.frame(`1` = 0.2, `2` = 0.3, `3` = 0.5, check.names = FALSE)
  expect_identical(
    maihda_prediction_panel_ordinal_categories(
      structure(list(y.levels = c("low", "mid", "high")), class = "clmm"), probs),
    c("low", "mid", "high"))
  # polr and other predict(type = "probs") methods: the columns name the levels.
  expect_identical(
    maihda_prediction_panel_ordinal_categories(
      structure(list(), class = "polr"),
      data.frame(a = 0.2, b = 0.3, c = 0.5)),
    c("a", "b", "c"))
  expect_error(
    maihda_prediction_panel_ordinal_categories(
      structure(list(y.levels = c("a", "b")), class = "clmm"), probs),
    "Could not match the ordinal model's 3 category probabilities")
  # No y.levels: an error, not the column names 1..K the label match used to read.
  expect_error(
    maihda_prediction_panel_ordinal_categories(structure(list(), class = "clmm"), probs),
    "Could not match")

  skip_if_not_installed("brms")
  mock <- function(y) {
    structure(list(formula = brms::bf(yo ~ x + (1 | stratum)),
                   data = data.frame(yo = y, x = 0)), class = "brmsfit")
  }
  # brms numbers an ordered factor by its levels, an integer response by its value.
  expect_identical(
    maihda_prediction_panel_ordinal_categories(
      mock(factor(c("3", "1"), levels = c("3", "2", "1"), ordered = TRUE)), probs),
    c("3", "2", "1"))
  expect_identical(maihda_prediction_panel_ordinal_categories(mock(c(1L, 3L)), probs),
                   c("1", "2", "3"))
  expect_error(
    maihda_prediction_panel_ordinal_categories(
      mock(factor("a", levels = c("a", "b", "c", "d"), ordered = TRUE)), probs),
    "Could not match")
  no_resp <- structure(list(formula = brms::bf(yo ~ x), data = data.frame(x = 0)),
                       class = "brmsfit")
  expect_error(maihda_prediction_panel_ordinal_categories(no_resp, probs), "Could not match")
})

test_that("the response name is read through a brmsformula and its addition terms", {
  d <- d15b_data(c("a", "b", "c"))
  # Every other model keeps as.character(formula)[2].
  expect_identical(maihda_prediction_panel_response_name(stats::lm(x ~ yo, data = d)), "x")
  g <- suppressWarnings(stats::glm(cbind(c(1, 2, 3), c(3, 2, 1)) ~ c(1, 2, 3),
                                   family = stats::binomial))
  expect_identical(maihda_prediction_panel_response_name(g), "cbind(c(1, 2, 3), c(3, 2, 1))")

  skip_if_not_installed("brms")
  mock <- function(f) structure(list(formula = f), class = "brmsfit")
  expect_identical(maihda_prediction_panel_response_name(
    mock(brms::bf(yo ~ x + (1 | stratum)))), "yo")
  expect_identical(maihda_prediction_panel_response_name(
    mock(brms::bf(yo | weights(.maihda_sw) ~ x + (1 | stratum)))), "yo")
  expect_identical(maihda_prediction_panel_response_name(
    mock(brms::bf(yo | thres(3) ~ x))), "yo")
  # A trials() count is not an observation-level outcome and stays unobserved.
  expect_null(maihda_prediction_panel_response_name(mock(brms::bf(s | trials(n) ~ x))))
  expect_null(maihda_prediction_panel_response_name(
    mock(brms::bf(s | weights(w) + trials(n) ~ x))))
})

# A brmsfit stand-in: stats::fitted() dispatches to a method registered on a private
# subclass and returns the summary shape brms returns, built from `estimate` (an
# ordinal nobs x K probability matrix, or a binomial probability vector) and indexed
# by the .row column of newdata. brms's own predict() and posterior_epred() fail on it.
d15b_mock_brms <- function(formula, data, family, estimate) {
  registerS3method("fitted", "d15b_mock_brmsfit", function(object, newdata = NULL, ...) {
    rows <- if (is.null(newdata)) object$data$.row else newdata$.row
    stats_names <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    est <- object$estimate
    if (is.matrix(est)) {
      lev <- levels(object$data$yo)
      out <- array(0.01, dim = c(length(rows), 4L, ncol(est)),
                   dimnames = list(NULL, stats_names, paste0("P(Y = ", lev, ")")))
      out[, "Estimate", ] <- est[rows, , drop = FALSE]
    } else {
      out <- matrix(0.01, nrow = length(rows), ncol = 4L,
                    dimnames = list(NULL, stats_names))
      out[, "Estimate"] <- est[rows]
    }
    out
  }, envir = asNamespace("stats"))
  data$.row <- seq_len(nrow(data))
  structure(list(formula = brms::bf(formula), data = data,
                 family = structure(list(family = family, link = "logit"), class = "family"),
                 estimate = estimate),
            class = c("d15b_mock_brmsfit", "brmsfit"))
}

test_that("the brms surprise panel finds the response and reads fitted() probabilities (Stan-free)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- d15b_data(c("0", "1", "2"))
  set.seed(3)
  est <- matrix(stats::runif(nrow(d) * 3, 0.05, 1), ncol = 3)
  est <- est / rowSums(est)
  # The weights() addition term is the shape of every weighted brms MAIHDA formula.
  m <- d15b_mock_brms(yo | weights(w) ~ x + (1 | stratum), d, "cumulative", est)
  p <- d15b_quiet(plot_prediction_deviation_panels(m, data = m$data, type = "ordinal"))
  p_obs <- est[cbind(seq_len(nrow(d)), as.integer(d$yo))]
  ref <- tapply(-log(p_obs), as.character(d$stratum), mean)
  got <- d15b_surprise(p)
  expect_true(all(is.finite(got)))
  expect_equal(as.numeric(got[names(ref)]), as.numeric(ref), tolerance = 1e-12)
  expect_identical(levels(p[[1]]$data$Category), c("0", "1", "2"))
})

test_that("the brms binomial panel reads the observed outcome for its deviance residuals (Stan-free)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- d15b_data(c("0", "1", "2"))
  d$yb <- as.integer(d$yo != "0")
  set.seed(4)
  pb <- stats::plogis(stats::rnorm(nrow(d)))
  mb <- d15b_mock_brms(yb ~ x + (1 | stratum), d, "bernoulli", pb)
  pd <- d15b_quiet(plot_prediction_deviation_panels(mb, data = mb$data, type = "binomial"))[[2]]$data
  dev <- sqrt(ifelse(d$yb == 1, -2 * log(pb), -2 * log1p(-pb)))
  ref_b <- tapply(dev, as.character(d$stratum), mean)
  expect_equal(pd$abs_res_dev, as.numeric(ref_b[as.character(pd$stratum)]), tolerance = 1e-12)
})

test_that("brms surprise and binomial panels score the observed response on real fits", {
  # Compiles two Stan models, so OPT-IN (set MAIHDA_TEST_BRMS=true). Sampler quality is
  # irrelevant: the reference is computed from the same fit's own draws.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  fit <- d15b_quiet(fit_maihda(yo ~ x + (1 | stratum), data = d15b_data(c("3", "2", "1")),
                               engine = "brms", family = "ordinal", chains = 1,
                               iter = 600, warmup = 300, refresh = 0, seed = 11))
  ep <- brms::posterior_epred(fit$model, newdata = fit$data)
  p_obs <- apply(ep, c(2, 3), mean)[cbind(seq_len(nrow(fit$data)), as.integer(fit$data$yo))]
  ref <- tapply(-log(p_obs), as.character(fit$data$stratum), mean)
  r <- d15b_capture(plot_prediction_deviation_panels(fit))
  expect_identical(d15b_panel_warnings(r$warnings), character(0))
  got <- d15b_surprise(r$value)
  expect_equal(as.numeric(got[names(ref)]), as.numeric(ref), tolerance = 1e-10)
  expect_identical(levels(r$value[[1]]$data$Category), c("3", "2", "1"))
  # Posterior means, not a tabulation of simulated responses: every call agrees.
  expect_identical(d15b_quiet(plot_prediction_deviation_panels(fit))[[2]]$data,
                   r$value[[2]]$data)

  db <- d15b_data(c("a", "b", "c"))
  db$yb <- as.integer(db$yo != "a")
  fb <- d15b_quiet(fit_maihda(yb ~ x + (1 | stratum), data = db, engine = "brms",
                              family = "binomial", chains = 1, iter = 600, warmup = 300,
                              refresh = 0, seed = 12))
  p_hat <- as.numeric(stats::fitted(fb$model, newdata = fb$data, summary = TRUE)[, "Estimate"])
  dev <- sqrt(ifelse(fb$data$yb == 1, -2 * log(p_hat), -2 * log1p(-p_hat)))
  ref_b <- tapply(dev, as.character(fb$data$stratum), mean)
  pd <- d15b_quiet(plot_prediction_deviation_panels(fb))[[2]]$data
  expect_equal(pd$abs_res_dev, as.numeric(ref_b[as.character(pd$stratum)]), tolerance = 1e-10)
})
