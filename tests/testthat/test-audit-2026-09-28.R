# Audit 2026-09-28: the case-level binomial prediction-deviation panel dropped the
# point of every case without an observed 0/1 outcome -- CONFIRMED (self-found during
# the 2026-09-27 F05/F04 pass; identical on HEAD 7f11dcf, so pre-existing), and wider.
#
# With no `stratum` column the panel draws one point per case and maps its shape to
# the observed outcome, as.factor() of the response column. An aggregated binomial has
# no such outcome: a cbind(successes, failures) response is a two-column matrix, so
# obs_outcome was factor(rep(NA, n)), and a brms `y | trials(n)` fit names no response
# column at all. ggplot's shape scale gives NA a missing shape and the point is removed
# at draw time -- all 200 points of a glm cbind() fit, of its glmer twin with a
# grouping factor not called `stratum`, and all 240 of a brms trials() stand-in, with
# nothing on the panel but the segments, bars and labels. PARTIAL on one mechanism: a
# proportion with trial weights is NOT all-NA -- it is a factor of the proportions,
# 64 levels on the fixture below, and the shape palette holds six, so 187 of 200 points
# were dropped rather than all. WIDER: a Bernoulli fit plotted on `data` that lacks the
# outcome column lost every point the same way, rows whose outcome is missing lost
# theirs, and a glm factor outcome with seven categories lost 21 of 150 past the sixth
# shape. The deviance residuals and labels were right throughout; only points vanished.
#
# FIX: the shape shows the observed outcome only where it is categorical (whole
# numbers when numeric) with no more values than the palette holds; every other case
# takes a fixed shape -- the default when no case shows an outcome, a hollow ring beside
# cases that do -- and the "Observed" legend title is dropped when nothing maps to it.
# A proportion therefore draws on one shape even when it has few values. Bernoulli
# panels are untouched: on 19 Bernoulli and control panels the ggplot_build() data are
# identical() to a git-archive export of HEAD.

p28_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The points ggplotGrob() draws in the second panel's plotting area -- the "points"
# grobs of the panel, so legend keys are not counted -- with each grob's plotting
# symbols, and every warning or message drawing it raised.
p28_drawn <- function(p) {
  conds <- character(0)
  g <- withCallingHandlers(ggplot2::ggplotGrob(p[[2]]),
    warning = function(e) {
      conds <<- c(conds, conditionMessage(e))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      conds <<- c(conds, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  panel <- g$grobs[[which(grepl("^panel", g$layout$name))[1]]]
  pts <- list()
  walk <- function(x) {
    if (inherits(x, "points")) pts[[length(pts) + 1L]] <<- x
    if (inherits(x, "gTree")) for (ch in x$children) walk(ch)
  }
  walk(panel)
  list(n = vapply(pts, function(x) length(x$x), 1L),
       pch = lapply(pts, function(x) x$pch),
       conditions = conds)
}

# Every case is drawn exactly once, plus one red ring per misclassified case, and
# drawing removed nothing.
p28_expect_all_drawn <- function(p, n_cases) {
  pd <- p[[2]]$data
  r <- p28_drawn(p)
  n_wrong <- if ("wrong" %in% names(pd)) sum(pd$wrong == "Wrong", na.rm = TRUE) else 0L
  expect_identical(nrow(pd), as.integer(n_cases))
  expect_true(all(is.finite(pd$fitted)))
  expect_identical(sum(r$n), as.integer(n_cases + n_wrong))
  expect_identical(grep("Removed|shape palette", r$conditions, value = TRUE), character(0))
  invisible(r)
}

# The finding's fixture: 200 rows of successes out of 2 to 20 trials.
p28_counts <- function() {
  set.seed(21)
  d <- data.frame(x = stats::rnorm(200))
  d$n <- sample(2:20, 200, replace = TRUE)
  d$s <- stats::rbinom(200, d$n, stats::plogis(0.3 * d$x))
  d$f <- d$n - d$s
  d$p <- d$s / d$n
  d
}

test_that("a cbind() fit draws a point for every case", {
  d <- p28_counts()
  g <- stats::glm(cbind(s, f) ~ x, data = d, family = stats::binomial)
  p <- p28_quiet(plot_prediction_deviation_panels(g))
  r <- p28_expect_all_drawn(p, 200)
  # One default shape, and no "Observed" legend with nothing mapped to it.
  expect_true(all(unlist(r$pch) == 19))
  expect_null(p[[2]]$labels$shape)
  # The premise: the residuals were right, only the points were missing.
  expect_equal(as.numeric(sort(p[[2]]$data$abs_res_dev)),
               as.numeric(sort(abs(stats::residuals(g, type = "deviance")))), tolerance = 1e-12)

  # The same fit given its data, and the glmer twin at case level.
  p28_expect_all_drawn(p28_quiet(plot_prediction_deviation_panels(g, data = d)), 200)
  d$grp <- factor(rep(1:10, each = 20))
  gm <- p28_quiet(lme4::glmer(cbind(s, f) ~ x + (1 | grp), data = d,
                              family = stats::binomial))
  p28_expect_all_drawn(p28_quiet(plot_prediction_deviation_panels(gm)), 200)
})

test_that("a proportion with trial weights draws every case on one shape", {
  d <- p28_counts()
  expect_length(unique(d$p), 64L)
  g <- stats::glm(p ~ x, data = d, family = stats::binomial, weights = n)
  p <- p28_quiet(plot_prediction_deviation_panels(g))
  r <- p28_expect_all_drawn(p, 200)
  expect_true(all(unlist(r$pch) == 19))
  expect_null(p[[2]]$labels$shape)

  # With two trials a row's proportion is 0, 0.5 or 1: few enough values for the
  # palette, but still an aggregated outcome with no shape to show.
  set.seed(5)
  d2 <- data.frame(x = stats::rnorm(200), n = 2)
  d2$p <- stats::rbinom(200, 2, stats::plogis(0.3 * d2$x)) / 2
  g2 <- stats::glm(p ~ x, data = d2, family = stats::binomial, weights = n)
  p2 <- p28_quiet(plot_prediction_deviation_panels(g2))
  r2 <- p28_expect_all_drawn(p2, 200)
  expect_true(all(unlist(r2$pch) == 19))
  expect_null(p2[[2]]$labels$shape)
})

test_that("a brms trials() fit draws a point for every case (Stan-free)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  # brms's fitted() returns expected success counts for a trials() fit; a private
  # subclass answers it with the per-row probability times the row's trials.
  registerS3method("fitted", "p28_mock_brmsfit", function(object, newdata = NULL, ...) {
    rows <- newdata$.row
    k <- newdata$n
    cbind(Estimate = object$p[rows] * k, Est.Error = 0.1 * k,
          Q2.5 = 0.9 * object$p[rows] * k, Q97.5 = pmin(1.1 * object$p[rows], 1) * k)
  }, envir = asNamespace("stats"))
  set.seed(917)
  d <- data.frame(x = stats::rnorm(240))
  d$n <- sample(3:30, 240, replace = TRUE)
  p_true <- stats::plogis(-0.4 + 0.5 * d$x)
  d$s <- stats::rbinom(240, d$n, p_true)
  d$.row <- seq_len(240)
  m <- structure(
    list(formula = brms::bf(s | trials(n) ~ x), data = d, p = p_true,
         family = structure(list(family = "binomial", link = "logit"), class = "family")),
    class = c("p28_mock_brmsfit", "brmsfit"))
  p <- p28_quiet(plot_prediction_deviation_panels(m, data = d, type = "binomial"))
  r <- p28_expect_all_drawn(p, 240)
  expect_true(all(unlist(r$pch) == 19))
  expect_null(p[[2]]$labels$shape)
  # Drawn where they belong: per-trial probabilities, and successes out of trials.
  pd <- p[[2]]$data
  expect_equal(pd$fitted, p_true[pd$id], tolerance = 1e-12)
  expect_equal(pd$abs_res_dev,
               maihda_binomial_abs_deviance_residual_agg(d$s[pd$id], d$n[pd$id], p_true[pd$id]),
               tolerance = 1e-12)
})

test_that("cases without an observed outcome are drawn: absent, or missing beside observed ones", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(150))
  d$y <- stats::rbinom(150, 1, stats::plogis(0.8 * d$x))
  g <- stats::glm(y ~ x, data = d, family = stats::binomial)

  # Prediction rows that do not carry the outcome: none wears an outcome's shape.
  p <- p28_quiet(plot_prediction_deviation_panels(g, data = d["x"]))
  r <- p28_expect_all_drawn(p, 150)
  expect_false(any(unlist(r$pch) %in% c(16, 17)))
  expect_null(p[[2]]$labels$shape)

  # Ten missing outcomes: those cases are drawn as hollow rings, which no observed
  # case uses, and the observed cases keep their outcome shapes and legend.
  dn <- d
  dn$y[1:10] <- NA
  p <- p28_quiet(plot_prediction_deviation_panels(g, data = dn))
  p28_expect_all_drawn(p, 150)
  expect_identical(p[[2]]$labels$shape, "Observed")
  b <- ggplot2::ggplot_build(p[[2]])
  point_layers <- which(vapply(p[[2]]$layers, function(l) inherits(l$geom, "GeomPoint"), TRUE))
  hollow <- Filter(function(i) identical(p[[2]]$layers[[i]]$aes_params$shape, 1) &&
                     is.null(p[[2]]$layers[[i]]$aes_params$colour), point_layers)
  expect_length(hollow, 1L)
  expect_setequal(p[[2]]$layers[[hollow]]$data$id, 1:10)
  shaped <- point_layers[1]
  expect_identical(sort(p[[2]]$layers[[shaped]]$data$id), 11:150)
  obs <- dn$y[p[[2]]$layers[[shaped]]$data$id]
  expect_identical(as.numeric(b$data[[shaped]]$shape), ifelse(obs == 1, 17, 16))
})

test_that("Bernoulli panels keep their outcome shapes (unchanged)", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(150))
  d$y <- stats::rbinom(150, 1, stats::plogis(0.8 * d$x))
  g <- stats::glm(y ~ x, data = d, family = stats::binomial)

  # The first point layer maps every case's outcome to its shape (the other layer is
  # the red ring of the misclassified).
  first_points <- function(p) {
    i <- which(vapply(p[[2]]$layers, function(l) inherits(l$geom, "GeomPoint"), TRUE))[1]
    list(layer = p[[2]]$layers[[i]], built = ggplot2::ggplot_build(p[[2]])$data[[i]])
  }
  p <- p28_quiet(plot_prediction_deviation_panels(g))
  p28_expect_all_drawn(p, 150)
  expect_identical(p[[2]]$labels$shape, "Observed")
  main <- first_points(p)
  expect_true(inherits(main$layer$mapping$shape, "quosure"))
  expect_identical(nrow(main$built), 150L)
  obs <- d$y[p[[2]]$data$id]
  expect_identical(as.numeric(main$built$shape), ifelse(obs == 1, 17, 16))

  # Rows showing one outcome value keep their shape and legend too.
  ps <- p28_quiet(plot_prediction_deviation_panels(g, data = d[d$y == 1, ]))
  p28_expect_all_drawn(ps, sum(d$y == 1))
  expect_identical(ps[[2]]$labels$shape, "Observed")
  main <- first_points(ps)
  expect_identical(nrow(main$built), sum(d$y == 1))
  expect_true(all(main$built$shape == 16))

  # A glm factor outcome with three categories keeps one shape per category; one with
  # seven, more than the palette holds, is drawn on one shape rather than losing cases.
  # (Read off the outcome layer: where these fits' cases are coded as glm codes them --
  # the first level against the rest -- the misclassified also get their red rings.)
  d$y3 <- factor(sample(c("lo", "mid", "hi"), 150, replace = TRUE), levels = c("lo", "mid", "hi"))
  p3 <- p28_quiet(plot_prediction_deviation_panels(stats::glm(y3 ~ x, data = d, family = stats::binomial)))
  p28_expect_all_drawn(p3, 150)
  main <- first_points(p3)
  expect_identical(nrow(main$built), 150L)
  expect_setequal(main$built$shape, c(15, 16, 17))
  d$y7 <- factor(sample(letters[1:7], 150, replace = TRUE))
  p7 <- p28_quiet(plot_prediction_deviation_panels(stats::glm(y7 ~ x, data = d, family = stats::binomial)))
  p28_expect_all_drawn(p7, 150)
  main <- first_points(p7)
  expect_identical(nrow(main$built), 150L)
  expect_true(all(main$built$shape == 19))
})

test_that("the shape rule: categorical outcomes the palette can hold, nothing else", {
  expect_true(maihda_prediction_panel_shape_outcome(c(0, 1, 1, NA)))
  expect_true(maihda_prediction_panel_shape_outcome(c(1, 1)))
  expect_true(maihda_prediction_panel_shape_outcome(c(1, 2)))
  expect_true(maihda_prediction_panel_shape_outcome(c(TRUE, FALSE)))
  expect_true(maihda_prediction_panel_shape_outcome(c("no", "yes")))
  expect_true(maihda_prediction_panel_shape_outcome(factor(c("a", "b", "c"))))
  expect_true(maihda_prediction_panel_shape_outcome(I(c(TRUE, FALSE))))
  expect_false(maihda_prediction_panel_shape_outcome(NULL))
  expect_false(maihda_prediction_panel_shape_outcome(cbind(s = 1:3, f = 3:1)))
  expect_false(maihda_prediction_panel_shape_outcome(c(0, 0.5, 1)))
  expect_false(maihda_prediction_panel_shape_outcome(c(0.25, 0.5)))
  expect_false(maihda_prediction_panel_shape_outcome(letters[1:7]))
  expect_false(maihda_prediction_panel_shape_outcome(c(0, 1, Inf)))
})
