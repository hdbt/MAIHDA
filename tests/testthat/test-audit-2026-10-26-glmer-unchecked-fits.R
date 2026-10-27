# Audit pass 2026-10-26: glmer fits that lme4 runs no convergence check on.
#
# glmer() computes the derivatives its convergence checks are taken from only for a fit of fewer
# rows and parameters than glmerControl()'s check.conv.nobsmax and check.conv.nparmax (10,000 and
# 20 by default, the fixed effects counted with the variance parameters); without them it runs no
# check at all, not even its singular test,
# so the record had nothing to go on: on simulated binomial and Poisson designs of 10,000 to
# 20,016 rows, glmer's own fits with a random slope (what fit_maihda() returns from a start, which
# skips its refits) were recorded as converged 0.017 to 14.3 log-likelihood units below the
# maximum.
# Per the user, a glmer fit outside growth models with a random slope that lme4 did not check, and
# that lme4 does not call singular, is measured as a fit whose checks fired is (the distance of
# test-audit-2026-10-20-glmer-slope-scale.R), and counts as not converged 0.1 standard errors or
# more from a stationary point; where that measure cannot be taken, far from the maximum, it is
# searched as a fit lme4 calls singular is (test-audit-2026-10-20-singular-flag-units.R).
# Intercept-only and growth fits are not measured.
#
# Random numbers are drawn only inside with_seed(), which puts the suite's stream back.

with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  force(expr)
}

quiet <- function(e) suppressWarnings(suppressMessages(e))

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 0.3) of `n` rows, each row with a
# birth year between 1940 and 2000 and a stratum slope deviation of SD 0.01 per year (deviations
# centred at 1970), and a binary outcome yb.
year_data <- function(seed, n = 80) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, 0.01)
    st <- rep(seq_len(18), each = n)
    x <- sample(1940:2000, length(st), replace = TRUE)
    eta <- -0.2 + add[st] + (0.02 + v[st]) * (x - 1970)
    data.frame(yb = stats::rbinom(length(st), 1, stats::plogis(eta)), x = x, stratum = factor(st))
  })
}

slope <- yb ~ x + (1 + x | stratum)

# A glmer fit placed at (theta, beta): optimisers that return their start. `...` goes to
# glmerControl().
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}
placed <- function(point, data, formula = slope, ...) {
  quiet(lme4::glmer(formula, data = data, family = stats::binomial,
                    start = list(theta = point$theta, fixef = point$beta),
                    control = lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE,
                                                 ...)))
}

# The maximum of glmer's deviance function searched with no bound on the Cholesky diagonal (PIRLS
# to 1e-12) and points along one direction of the centred and scaled parameters below it in
# deviance:
# - year_data(1), 1,440 rows: the maximum, and points 0.0081, 0.0121 and 0.03 below it (those of
#   test-audit-2026-10-20-glmer-slope-scale.R);
# - year_data(1, 556), 10,008 rows: a deviance of 13251.5916131 at the maximum, and a point 0.05
#   below it; glmer's own fit stops 40.4 below it, singular.
max_1 <- list(theta = c(31.821741, -0.016132383, 0.00015508396), beta = c(-41.57003, 0.021048822))
below_0.0081 <- list(theta = c(31.519655, -0.015981713, 0.00015503863),
                     beta = c(-40.957266, 0.020736257))
below_0.0121 <- list(theta = c(31.453792, -0.015948879, 0.00015502309),
                     beta = c(-40.821168, 0.020666834))
below_0.03 <- list(theta = c(31.248715, -0.015846687, 0.00015496039),
                   beta = c(-40.391206, 0.020447514))
max_big <- list(theta = c(21.525687, -0.010989632, 0.00016047682), beta = c(-44.279504, 0.022454854))
below_big_0.05 <- list(theta = c(21.994717, -0.011223628, 0.00016476934),
                       beta = c(-43.859924, 0.022243428))

# ---- which fits lme4 does not check -------------------------------------------------------

test_that("glmer runs no convergence check on a fit of 10,000 rows or with 20 parameters, singular or not", {
  skip_on_cran()
  big <- year_data(1, 556)
  m <- placed(below_big_0.05, big)
  expect_false(lme4::isSingular(m))
  expect_null(m@optinfo$derivs)
  expect_true(maihda_lme4_unchecked(m))
  # A singular glmer fit of 10,008 rows has no entry either: glmer skips its singular test too.
  singular <- placed(list(theta = c(below_big_0.05$theta[1:2], 0), beta = below_big_0.05$beta), big)
  expect_true(lme4::isSingular(singular))
  expect_true(maihda_lme4_unchecked(singular))
  # 9,990 rows: lme4 checks the same kind of point, and its checks fire.
  small <- placed(below_big_0.05, year_data(1, 555))
  expect_gt(length(small@optinfo$conv$lme4$messages), 0L)
  expect_false(maihda_lme4_unchecked(small))
  # Twenty parameters on 1,440 rows, the fixed effects counted with the variance parameters, and
  # lme4's own switches.
  d <- year_data(1)
  expect_false(maihda_lme4_unchecked(placed(below_0.03, d)))
  expect_true(maihda_lme4_unchecked(placed(below_0.03, d, calc.derivs = FALSE)))
  # A fit lme4 checked and found clean carries an empty entry, which is not counted.
  clean <- quiet(lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = lme4::cbpp, family = stats::binomial))
  expect_false(is.null(clean@optinfo$conv$lme4))
  expect_length(clean@optinfo$conv$lme4, 0L)
  expect_false(maihda_lme4_unchecked(clean))
  expect_true(maihda_lme4_unchecked(placed(below_0.03, d, check.conv.nobsmax = 100)))
  for (i in 1:4) d[[c("z", "w", "v", "u")[i]]] <- with_seed(4 + i, stats::rnorm(1440))
  identity_5 <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1)
  twenty <- placed(list(theta = identity_5, beta = numeric(5)), d,
                   yb ~ z + w + v + u + (1 + z + w + v + u | stratum))
  expect_identical(length(lme4::getME(twenty, "theta")) + length(lme4::fixef(twenty)), 20L)
  expect_true(maihda_lme4_unchecked(twenty))
  nineteen <- placed(list(theta = identity_5, beta = numeric(4)), d,
                     yb ~ z + w + v + (1 + z + w + v + u | stratum))
  expect_false(maihda_lme4_unchecked(nineteen))
})

# ---- the record ---------------------------------------------------------------------------

test_that("an unchecked glmer slope fit 0.1 standard errors or more from a stationary point is recorded as not converged", {
  skip_on_cran()
  d <- year_data(1)
  m <- placed(below_0.03, d, calc.derivs = FALSE)
  expect_false(lme4::isSingular(m))
  dg <- maihda_fit_diagnostics(m)
  expect_false(dg$converged)
  # 0.03 below in deviance is sqrt(0.03) standard errors from the maximum.
  expect_lt(abs(dg$stationary_distance / sqrt(0.03) - 1), 0.01)
  expect_true(dg$lme4_unchecked)
  expect_identical(dg$messages, paste0(
    "lme4 ran no convergence check on this fit (by default it runs none from 10,000 rows or 20 ",
    "parameters on, fixed effects included)."))
  expect_identical(maihda_format_fit_diagnostics(dg), c(
    "Convergence problems of the lme4 fit:",
    paste0("  - lme4 ran no convergence check on this fit (by default it runs none from 10,000 rows ",
           "or 20 parameters on, fixed effects included)."),
    "In centred and scaled units the fit's gradient and Hessian put it 0.17 standard",
    "  errors from a stationary point."))
  # From both sides of the line, and at the maximum: converged, nothing recorded, nothing printed.
  expect_false(maihda_fit_diagnostics(placed(below_0.0121, d, calc.derivs = FALSE))$converged)
  for (point in list(below_0.0081, max_1)) {
    inside <- maihda_fit_diagnostics(placed(point, d, calc.derivs = FALSE))
    expect_true(inside$converged)
    expect_null(inside$stationary_distance)
    expect_null(inside$lme4_unchecked)
    expect_identical(maihda_format_fit_diagnostics(inside), character(0))
  }
})

test_that("the record's line for an unchecked glmer fit is a tenth of a standard error", {
  skip_on_cran()
  m <- placed(below_0.03, year_data(1), calc.derivs = FALSE)
  at <- function(distance) {
    local_mocked_bindings(maihda_lme4_stationary_distance = function(model) distance)
    maihda_fit_diagnostics(m)
  }
  on_line <- at(0.1)
  expect_false(on_line$converged)
  expect_identical(on_line$stationary_distance, 0.1)
  below <- at(0.0999)
  expect_true(below$converged)
  expect_null(below$stationary_distance)
  # Where the measure reads NA the search decides (the next block).
})

test_that("only an unchecked glmer fit with a random slope outside growth models that lme4 does not call singular is measured", {
  skip_on_cran()
  measured <- 0L
  local_mocked_bindings(maihda_lme4_stationary_distance = function(model) {
    measured <<- measured + 1L
    0.5
  })
  d <- year_data(1)
  m <- placed(below_0.03, d, calc.derivs = FALSE)
  expect_false(maihda_fit_diagnostics(m)$converged)
  expect_identical(measured, 1L)
  # A growth model, a singular fit, a fit without a random slope, an optimiser failure: not.
  expect_true(maihda_fit_diagnostics(m, longitudinal_info = list(id = "id"))$converged)
  singular <- placed(list(theta = c(below_0.03$theta[1:2], 0), beta = below_0.03$beta), d,
                     calc.derivs = FALSE)
  expect_true(lme4::isSingular(singular))
  maihda_fit_diagnostics(singular)
  intercept_only <- quiet(lme4::glmer(yb ~ x + (1 | stratum), data = d, family = stats::binomial,
                                      control = lme4::glmerControl(calc.derivs = FALSE)))
  expect_true(maihda_lme4_unchecked(intercept_only))
  expect_true(maihda_fit_diagnostics(intercept_only)$converged)
  failed <- m
  failed@optinfo$conv$opt <- 1L
  expect_false(maihda_fit_diagnostics(failed)$converged)
  expect_identical(measured, 1L)
  # A fit lme4 checked is measured only when its checks fire, as before, and not recorded here.
  fired <- maihda_fit_diagnostics(placed(below_0.03, d))
  expect_identical(measured, 2L)
  expect_null(fired$lme4_unchecked)
  # An lmer fit lme4 did not check is searched, not measured this way.
  d$y <- d$yb + with_seed(3, stats::rnorm(1440, 0, 0.5))
  lm_fit <- quiet(lme4::lmer(y ~ x + (1 | stratum), data = d,
                             control = lme4::lmerControl(calc.derivs = FALSE)))
  expect_false(lme4::isSingular(lm_fit))
  expect_true(maihda_lme4_unchecked(lm_fit))
  expect_null(maihda_fit_diagnostics(lm_fit)$lme4_unchecked)
  expect_identical(measured, 2L)
})

test_that("where the measure cannot be taken, an unchecked glmer fit is searched as a singular one is", {
  skip_on_cran()
  d <- year_data(1)
  # 19 below the maximum in deviance, the slope's entry at 0.0005: the measure's Hessian is not
  # positive definite there, and a search from the fit reaches 0.01 or more below it.
  far <- placed(list(theta = c(31.821741, -0.016132383, 0.0005), beta = max_1$beta), d,
                calc.derivs = FALSE)
  gap <- -2 * as.numeric(stats::logLik(far)) + 2 * as.numeric(stats::logLik(placed(max_1, d)))
  expect_gt(gap, 15)
  expect_false(lme4::isSingular(far))
  expect_true(is.na(maihda_lme4_stationary_distance(far)))
  dg <- maihda_fit_diagnostics(far)
  expect_false(dg$converged)
  expect_true(dg$lme4_unchecked)
  expect_gte(dg$search_gain, 0.01)
  expect_lte(dg$search_gain, gap + 1e-6)
  expect_null(dg$stationary_distance)
  expect_true(startsWith(dg$messages, paste0(
    "lme4 ran no convergence check on this fit (by default it runs none from 10,000 rows or 20 ",
    "parameters on, fixed effects included), and a search of the same model with its random-effect ",
    "terms and fixed effects centred and scaled reaches a deviance ")))
  expect_identical(maihda_format_fit_diagnostics(dg)[1], "Convergence problems of the lme4 fit:")
  # The gate, with the measure and the search stood in for.
  searches <- 0L
  measure <- NA_real_
  local_mocked_bindings(
    maihda_lme4_stationary_distance = function(model) measure,
    maihda_glmer_search_gain = function(model, ...) {
      searches <<- searches + 1L
      0.5
    })
  m <- placed(below_0.03, d, calc.derivs = FALSE)
  flagged <- maihda_fit_diagnostics(m)
  expect_false(flagged$converged)
  expect_identical(flagged$search_gain, 0.5)
  expect_identical(searches, 1L)
  # Refitted by fit_maihda() from such searches: not searched again.
  expect_true(maihda_fit_diagnostics(m, searched = TRUE)$converged)
  # A measure that can be taken decides alone.
  measure <- 0.05
  expect_true(maihda_fit_diagnostics(m)$converged)
  measure <- 0.5
  expect_null(maihda_fit_diagnostics(m)$search_gain)
  expect_identical(searches, 1L)
  # A search that reaches nothing 0.01 lower leaves the fit converged.
  measure <- NA_real_
  local_mocked_bindings(maihda_glmer_search_gain = function(model, ...) 0.005)
  expect_true(maihda_fit_diagnostics(m)$converged)
})

# ---- live fits ----------------------------------------------------------------------------

test_that("a glmer fit of 10,000 rows made from a start is measured, and print() says why", {
  skip_slow_test()
  skip_on_cran()
  big <- year_data(1, 556)
  control <- lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE)
  m <- quiet(fit_maihda(slope, data = big, family = "binomial",
                        start = list(theta = below_big_0.05$theta, fixef = below_big_0.05$beta),
                        control = control))
  expect_null(m$model@optinfo$conv$lme4)
  expect_false(m$diagnostics$converged)
  expect_lt(abs(m$diagnostics$stationary_distance / sqrt(0.05) - 1), 0.03)
  shown <- utils::capture.output(print(m))
  expect_true(any(grepl("lme4 ran no convergence check on this fit", shown, fixed = TRUE)))
  expect_true(any(grepl("Convergence problems of the lme4 fit:", shown, fixed = TRUE)))
  # At the maximum: converged, nothing recorded.
  at <- quiet(fit_maihda(slope, data = big, family = "binomial",
                         start = list(theta = max_big$theta, fixef = max_big$beta), control = control))
  expect_true(at$diagnostics$converged)
  expect_null(at$diagnostics$lme4_unchecked)
})
