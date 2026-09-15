# Audit 2026-09-15 (third pass), finding 2: an lme4 fit carrying BOTH a formula
# offset() and an external offset = ... accepted newdata and silently left the external
# offset out of the prediction -- CONFIRMED.
#
# predict_maihda() refuses newdata for a fit with an external offset, because
# predict.merMod re-evaluates only formula offset() terms on newdata and the external
# offset survives solely as its fitted values. The guard, maihda_lme4_has_external_offset(),
# found the model frame's "(offset)" column but then returned FALSE whenever the formula
# also had an offset() term, on the premise that lme4 names the column "(offset)" for
# either kind. It does not: lme4 builds its frame with stats::model.frame(), which gives
# "(offset)" to the offset= argument alone and keeps a formula offset() as "offset(o1)"
# (measured on lme4 2.0.1 fits and in base R), so the extra check separated nothing and
# only switched the guard off on a fit with both. On the training rows
# predict(m, newdata = d) then differed from predict(m) by exactly o2 (up to 0.9997,
# within 4.4e-16 of o2; a factor exp(o2), up to 2.72, on the Poisson response scale, and
# the same under glmer.nb), and a five-row request returned the offset-less values. The
# raw-merMod twin the prediction-deviation panels use, maihda_mermod_has_external_offset(),
# already tested the column's presence alone.
#
# FIX: the presence of "(offset)" is the whole test. Unchanged: training predictions (no
# newdata), newdata predictions of a formula-offset-only fit, type = "strata", and the
# longitudinal grids, which add the external and formula parts themselves.

d15c_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

d15c_data <- function() {
  set.seed(20260915)
  n <- 400
  d <- data.frame(g1 = factor(sample(c("F", "M"), n, TRUE)),
                  g2 = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = stats::rnorm(n), o1 = stats::runif(n), o2 = stats::runif(n))
  d <- make_strata(d, vars = c("g1", "g2"))$data
  st <- as.character(d$stratum)
  u <- stats::setNames(stats::rnorm(length(unique(st)), sd = 0.5), sort(unique(st)))
  d$y <- 1 + 0.5 * d$x + d$o1 + d$o2 + u[st] + stats::rnorm(n, sd = 0.3)
  d$cnt <- stats::rpois(n, exp(-0.5 + 0.3 * d$x + d$o1 + d$o2 + u[st]))
  d
}

test_that("a fit with a formula offset() and an external offset refuses newdata", {
  d <- d15c_data()
  m <- d15c_quiet(fit_maihda(y ~ x + offset(o1) + (1 | stratum), data = d, offset = o2))
  # The frame keeps each offset under its own column.
  expect_true(all(c("offset(o1)", "(offset)") %in% names(m$data)))
  expect_true(maihda_lme4_has_external_offset(m))

  expect_error(predict(m, newdata = d, scale = "link"), "external offset")
  expect_error(predict(m, newdata = d, scale = "response"), "external offset")
  expect_error(predict_maihda(m, newdata = d[1:5, ]), "external offset")
  expect_error(predict(m, newdata = d[1:5, ], allow_new_levels = TRUE), "external offset")
  # The prediction-deviation panels' own guard already refused it.
  expect_error(maihda_prediction_panel_fitted(m$model, d, "gaussian"), "external offset")
})

test_that("the newdata guard follows the external offset alone, whatever the formula carries", {
  d <- d15c_data()
  fits <- list(
    none = d15c_quiet(fit_maihda(y ~ x + (1 | stratum), data = d)),
    form = d15c_quiet(fit_maihda(y ~ x + offset(o1) + (1 | stratum), data = d)),
    ext  = d15c_quiet(fit_maihda(y ~ x + (1 | stratum), data = d, offset = o2)),
    both = d15c_quiet(fit_maihda(y ~ x + offset(o1) + (1 | stratum), data = d, offset = o2))
  )
  wrapper <- vapply(fits, maihda_lme4_has_external_offset, logical(1))
  expect_identical(wrapper, c(none = FALSE, form = FALSE, ext = TRUE, both = TRUE))
  # The same answer as the raw-merMod helper, and exactly where lme4 writes "(offset)".
  expect_identical(
    vapply(fits, function(m) maihda_mermod_has_external_offset(m$model), logical(1)),
    wrapper)
  expect_identical(vapply(fits, function(m) "(offset)" %in% names(m$data), logical(1)),
                   wrapper)

  # Negative controls: without an external offset newdata still predicts, and on the
  # training rows it reproduces the training predictions; an external offset alone is
  # still refused.
  for (k in c("none", "form")) {
    expect_equal(as.numeric(predict(fits[[k]], newdata = d)),
                 as.numeric(predict(fits[[k]])), tolerance = 1e-10)
  }
  expect_error(predict(fits$ext, newdata = d), "external offset")
})

test_that("training and stratum predictions of a both-offsets fit are unchanged", {
  # Negative control: passes before and after the fix.
  d <- d15c_data()
  m <- d15c_quiet(fit_maihda(y ~ x + offset(o1) + (1 | stratum), data = d, offset = o2))
  beta <- lme4::fixef(m$model)
  X <- stats::model.matrix(~ x, d)
  b <- lme4::ranef(m$model)$stratum[as.character(d$stratum), "(Intercept)"]
  expect_equal(as.numeric(predict(m)),
               as.numeric(X %*% beta[colnames(X)]) + b + d$o1 + d$o2, tolerance = 1e-10)
  st <- predict(m, newdata = d[1:5, ], type = "strata")
  expect_s3_class(st, "data.frame")
  expect_setequal(as.character(st$stratum), unique(as.character(d$stratum[1:5])))
})

test_that("count fits with both offsets refuse newdata (glmer and glmer.nb)", {
  skip_on_cran()
  d <- d15c_data()
  pb <- d15c_quiet(fit_maihda(cnt ~ x + offset(o1) + (1 | stratum), data = d,
                              family = "poisson", offset = o2))
  expect_true(maihda_lme4_has_external_offset(pb))
  expect_error(predict(pb, newdata = d, scale = "link"), "external offset")
  expect_error(predict(pb, newdata = d, scale = "response"), "external offset")
  expect_equal(as.numeric(predict(pb, scale = "response")),
               as.numeric(stats::fitted(pb$model)), tolerance = 1e-10)

  nb <- d15c_quiet(fit_maihda(cnt ~ x + offset(o1) + (1 | stratum), data = d,
                              family = "negbinomial", offset = o2))
  expect_true(maihda_lme4_has_external_offset(nb))
  expect_error(predict(nb, newdata = d, scale = "response"), "external offset")
})

test_that("a longitudinal count fit with both offsets refuses newdata and matches one total offset", {
  skip_on_cran()
  set.seed(5)
  nid <- 130
  g1 <- sample(c("F", "M"), nid, TRUE)
  g2 <- sample(c("A", "B", "C"), nid, TRUE)
  base <- do.call(rbind, lapply(seq_len(nid), function(i)
    data.frame(id = i, g1 = g1[i], g2 = g2[i], time = 0:2,
               o1 = log(stats::runif(3, 1, 6)), o2 = log(stats::runif(3, 1, 3)),
               stringsAsFactors = FALSE)))
  base$o12 <- base$o1 + base$o2
  # Real stratum and person effects: without them the fit's variance components are 0,
  # VPC(t) is about 1e-12 whatever the offsets are, and the VPC(t) comparisons below
  # would pass however the grid treated either offset.
  s <- interaction(base$g1, base$g2, drop = TRUE)
  u_s <- stats::rnorm(nlevels(s), sd = 0.5)[s]
  u_i <- stats::rnorm(nid, sd = 0.3)[base$id]
  base$y <- stats::rpois(nrow(base), exp(-0.2 + 0.15 * base$time + base$o12 + u_s + u_i))
  fl <- function(f, ...) d15c_quiet(fit_maihda(f, data = base, family = "poisson",
                                               id = "id", time = "time", ...))
  l_both <- fl(y ~ 1 + offset(o1) + (1 | g1:g2), offset = o2)
  l_ext  <- fl(y ~ 1 + (1 | g1:g2), offset = o12)
  l_form <- fl(y ~ 1 + offset(o1) + offset(o2) + (1 | g1:g2))

  expect_error(predict(l_both, newdata = base[1:6, ]), "external offset")
  expect_equal(as.numeric(predict(l_form, newdata = l_form$original_data)),
               as.numeric(predict(l_form)), tolerance = 1e-8)

  # Negative controls: the VPC(t) grid and the fixed trajectory add both parts, so the
  # same total offset gives the same answer however it was supplied.
  s_both <- d15c_quiet(summary(l_both))
  expect_true(all(s_both$longitudinal$vpc_t$estimate > 0.2))
  expect_equal(s_both$longitudinal$vpc_t, d15c_quiet(summary(l_ext))$longitudinal$vpc_t,
               tolerance = 1e-6)
  expect_equal(s_both$longitudinal$vpc_t, d15c_quiet(summary(l_form))$longitudinal$vpc_t,
               tolerance = 1e-6)
  grid <- sort(unique(base$time))
  expect_equal(maihda_longitudinal_fixed_trajectory(l_both, grid),
               maihda_longitudinal_fixed_trajectory(l_form, grid), tolerance = 1e-6)
})
