# Audit pass 2026-10-20: glmer fits with a random slope on a covariate far from zero or in large
# units recorded as not converged at the maximum.
#
# lme4 takes its convergence checks of a glmer fit by central differences with a step of 1e-4 in
# the variance parameters and fixed effects as the formula has them. The re-check of
# test-audit-2026-10-14-lme4-covariate-scale.R took its own differences with the fixed effects
# centred and scaled but theta as fitted, where a step along a random slope's entry moves each
# observation's random-effect variance by an amount that grows with the square of its covariate.
# On simulated binomial fits with a slope on a birth year (18 strata of 80 rows) 6 of 12 fits
# returned by fit_maihda() were recorded as not converged, every one within 0.033 standard errors
# of the maximum found by an unbounded search of glmer's own deviance function; the re-check gave
# NA on four and 0.41 / 0.35 on two fits within 0.0001 of it, and with a slope and no intercept
# in the block it overruled lme4 at points 0.11 to 0.15 standard errors short. Per the user, the
# re-check now takes theta with each random-effect block's terms centred and scaled (in a block
# without an intercept, divided by their root mean square), as an lmer fit's is
# (test-audit-2026-10-20-lmer-optimum-messages.R), with the differences free to cross a variance
# bound; a fit with a block lme4 holds as a diagonal, diag(1 + w | g), is measured as the same
# model written with || (test-audit-2026-10-25-diag-block-scale.R). On 1,245 points placed at
# known distances from the maximum (19 designs:
# binomial and Poisson, a covariate on a scale of 1, an age, a birth year, a calendar year, an
# income, a year in thousands of years, blocks without an intercept, growth models) the distance
# it gives was 0.92 to 1.05 times the true one (1% to 99%).
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

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 0.3) of 80 rows, each row with a
# birth year between 1940 and 2000 and a stratum slope deviation of SD `slope_sd` per year
# (deviations centred at 1970): yb a binary outcome, or with `family = "poisson"` yc a count.
year_data <- function(seed, family = "binomial", slope_sd = 0.01) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, slope_sd)
    st <- rep(seq_len(18), each = 80)
    x <- sample(1940:2000, length(st), replace = TRUE)
    if (family == "binomial") {
      eta <- -0.2 + add[st] + (0.02 + v[st]) * (x - 1970)
      data.frame(yb = stats::rbinom(length(st), 1, stats::plogis(eta)), x = x, stratum = factor(st))
    } else {
      eta <- 0.5 + add[st] + (0.01 + v[st] / 2) * (x - 1970)
      data.frame(yc = stats::rpois(length(st), exp(eta)), x = x, stratum = factor(st))
    }
  })
}

slope <- yb ~ x + (1 + x | stratum)
noint <- yb ~ x + (0 + x | stratum)

# A glmer fit placed at (theta, beta): optimisers that return their start, so lme4's checks are
# taken at that point exactly.
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}
placed <- function(point, data, formula = slope, family = stats::binomial) {
  quiet(lme4::glmer(formula, data = data, family = family,
                    start = list(theta = point$theta, fixef = point$beta),
                    control = lme4::glmerControl(optimizer = list(stay, stay),
                                                 nAGQ0initStep = FALSE)))
}
deviance_of <- function(m) -2 * as.numeric(stats::logLik(m))

# The measure with the fixed effects centred and scaled and theta as fitted.
fitted_theta_distance <- function(m) {
  dev <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
  theta <- lme4::getME(m, "theta")
  beta <- lme4::fixef(m)
  k <- length(theta)
  A <- maihda_centre_scale_map(lme4::getME(m, "X"))
  f <- function(p) dev(c(p[seq_len(k)], drop(A %*% p[-seq_len(k)])))
  at <- c(theta, solve(A, beta))
  lower <- c(lme4::getME(m, "lower"), rep(-Inf, length(beta)))
  maihda_newton_distance(maihda_central_derivatives(f, at, lower = lower, fx = f(at)))
}

# Seed 1: the maximum of glmer's deviance function searched with no bound on the Cholesky
# diagonal (PIRLS to 1e-12), where lme4 reports max|grad| = 0.59, and points along one direction
# of the centred and scaled parameters 0.0081, 0.0121 and 0.03 below it in deviance.
max_1 <- list(theta = c(31.821741, -0.016132383, 0.00015508396), beta = c(-41.57003, 0.021048822))
below_0.0081 <- list(theta = c(31.519655, -0.015981713, 0.00015503863),
                     beta = c(-40.957266, 0.020736257))
below_0.0121 <- list(theta = c(31.453792, -0.015948879, 0.00015502309),
                     beta = c(-40.821168, 0.020666834))
below_0.03 <- list(theta = c(31.248715, -0.015846687, 0.00015496039),
                   beta = c(-40.391206, 0.020447514))

# The point carried to the same model on the covariate (x - shift) * k: the block (Intercept), x
# and the fixed effects become T b with T = [1, shift; 0, 1 / k].
carry <- function(point, shift = 0, k = 1) {
  Tm <- matrix(c(1, 0, shift, 1 / k), 2)
  L <- matrix(c(point$theta[1], point$theta[2], 0, point$theta[3]), 2)
  L2 <- maihda_chol_psd(Tm %*% tcrossprod(L) %*% t(Tm))
  list(theta = L2[lower.tri(L2, diag = TRUE)], beta = drop(Tm %*% point$beta))
}

# ---- the distance --------------------------------------------------------------------------

test_that("a glmer slope fit is measured in its random-effect terms centred and scaled, whatever the covariate's spelling", {
  skip_on_cran()
  skip_slow_test()
  d <- year_data(1)
  at <- placed(max_1, d)
  expect_lt(maihda_lme4_stationary_distance(at), 1e-4)
  # Theta as fitted, at the same point.
  expect_false(isTRUE(fitted_theta_distance(at) < 0.1))
  # The year as x - 1900, and in decades: the same model, the same distance.
  d$x1900 <- d$x - 1900
  d$dec <- d$x / 10
  spellings <- function(point) {
    list(placed(point, d),
         placed(carry(point, shift = 1900), d, yb ~ x1900 + (1 + x1900 | stratum)),
         placed(carry(point, k = 1 / 10), d, yb ~ dec + (1 + dec | stratum)))
  }
  near <- vapply(spellings(max_1), maihda_lme4_stationary_distance, numeric(1))
  expect_true(all(near < 1e-4))
  far <- vapply(spellings(below_0.03), maihda_lme4_stationary_distance, numeric(1))
  expect_lt(max(abs(far / far[1] - 1)), 1e-3)
})

test_that("the distance is the deviance's shortfall below the maximum", {
  skip_on_cran()
  d <- year_data(1)
  top <- deviance_of(placed(max_1, d))
  for (case in list(list(point = below_0.0081, gap = 0.0081),
                    list(point = below_0.0121, gap = 0.0121),
                    list(point = below_0.03, gap = 0.03))) {
    m <- placed(case$point, d)
    # The deviance glmer reports, from its iterations to lme4's default tolerance, is within
    # about 3e-5 of the converged one at these points.
    gap <- deviance_of(m) - top
    expect_lt(abs(gap - case$gap), 1e-4)
    # A ratio, not expect_equal(tolerance =): under this suite's testthat edition a scalar
    # passes that whenever the absolute difference is below the tolerance.
    expect_lt(abs(maihda_lme4_stationary_distance(m) / sqrt(gap) - 1), 0.02)
  }
})

# Seeds 1 and 3 fitted with the slope alone in the block, (0 + x | stratum): the maximum of
# each, and a point of seed 3 0.0225 below its maximum.
noint_max_1 <- list(theta = 0.00014645191, beta = c(-41.534808, 0.021030495))
noint_max_3 <- list(theta = 0.00011099502, beta = c(-25.435037, 0.012794443))
noint_below_3 <- list(theta = 0.00011662035, beta = c(-25.669047, 0.012914299))

test_that("a slope without an intercept is measured with its block divided by its root mean square", {
  skip_on_cran()
  at <- placed(noint_max_1, year_data(1), noint)
  expect_gt(length(at@optinfo$conv$lme4$messages), 0L)
  expect_lt(maihda_lme4_stationary_distance(at), 1e-4)
  expect_true(maihda_fit_diagnostics(at)$converged)
  # Theta as fitted is 1.5e-4 per year, so a step of 1e-4 along it is most of its value.
  expect_gt(fitted_theta_distance(at), 0.1)
  # 0.0225 below the maximum, 0.15 standard errors: not converged, and how far. Theta as fitted
  # puts it within 0.1.
  d3 <- year_data(3)
  top <- deviance_of(placed(noint_max_3, d3, noint))
  short <- placed(noint_below_3, d3, noint)
  gap <- deviance_of(short) - top
  expect_lt(abs(gap - 0.0225), 1e-4)
  expect_gt(length(short@optinfo$conv$lme4$messages), 0L)
  expect_lt(fitted_theta_distance(short), 0.1)
  ds <- maihda_fit_diagnostics(short)
  expect_false(ds$converged)
  expect_lt(abs(ds$stationary_distance / sqrt(gap) - 1), 0.02)
  expect_null(ds$lme4_overruled)
})

test_that("the differences may cross a variance bound", {
  skip_on_cran()
  # The year in thousands of years from 1970, the strata sharing one slope: the maximum is a
  # block of rank 1. At a point with its slope diagonal at 5e-5 in the centred and scaled
  # terms, within a step of its bound, lme4's own entry is 0.0028, so lme4 checks such a fit.
  d <- year_data(1, slope_sd = 0)
  d$x <- (d$x - 1970) / 1000
  near <- placed(list(theta = c(0.35255329, -8.6059568, 0.0028308795),
                      beta = c(-0.064069481, 20.946614)), d)
  expect_false(lme4::isSingular(near))
  expect_lt(maihda_lme4_stationary_distance(near), 0.01)
})

# Strata as in year_data() with a fixed effect of the birth year x and a stratum slope of SD 0.4 on
# w, a standard normal covariate.
w_data <- function(seed) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, 0.4)
    st <- rep(seq_len(18), each = 80)
    x <- sample(1940:2000, length(st), replace = TRUE)
    w <- stats::rnorm(length(st))
    eta <- -0.2 + add[st] + 0.02 * (x - 1970) + (0.5 + v[st]) * w
    data.frame(yb = stats::rbinom(length(st), 1, stats::plogis(eta)), x = x, w = w,
               stratum = factor(st))
  })
}

test_that("a fit with a block lme4 holds as a diagonal is measured as the same model written with ||", {
  skip_on_cran()
  skip_slow_test()
  skip_if(utils::packageVersion("lme4") < "2.0.0")
  d <- w_data(1)
  m <- quiet(lme4::glmer(yb ~ x + w + diag(1 + w | stratum), data = d, family = stats::binomial))
  # Theta holds the block's two standard deviations, not a Cholesky factor, and lme4's checks
  # fire on the birth year's coefficient.
  expect_length(lme4::getME(m, "theta"), 2L)
  expect_gt(length(m@optinfo$conv$lme4$messages), 0L)
  expect_true(maihda_fit_diagnostics(m)$converged)
  # The same model with the block split, (1 + w || stratum): the same distance.
  split <- quiet(lme4::glmer(yb ~ x + w + (1 + w || stratum), data = d, family = stats::binomial))
  expect_lt(abs(maihda_lme4_stationary_distance(split) - maihda_lme4_stationary_distance(m)), 1e-4)
})

test_that("an intercept-only glmer fit is measured as before", {
  skip_on_cran()
  # A block of the intercept alone has nothing to centre or scale.
  for (counts in c(FALSE, TRUE)) {
    m <- if (counts) {
      quiet(lme4::glmer(yc ~ x + (1 | stratum), data = year_data(1, "poisson"),
                        family = stats::poisson))
    } else {
      quiet(lme4::glmer(yb ~ x + (1 | stratum), data = year_data(1), family = stats::binomial))
    }
    expect_identical(maihda_lme4_stationary_distance(m), fitted_theta_distance(m))
  }
})

# ---- the record ---------------------------------------------------------------------------

test_that("a glmer slope fit at its maximum is converged although lme4's checks fire", {
  skip_on_cran()
  m <- placed(max_1, year_data(1))
  said <- as.character(m@optinfo$conv$lme4$messages)
  # The premise: lme4's checks fire at the maximum, and its optimiser reports success.
  expect_gt(length(said), 0L)
  expect_true(any(grepl("failed to converge", said, fixed = TRUE)))
  expect_identical(as.integer(m@optinfo$conv$opt), 0L)
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$converged)
  expect_identical(dg$messages, character(0))
  expect_identical(dg$lme4_overruled, said)
  expect_lt(dg$stationary_distance, 1e-4)
  expect_identical(maihda_format_fit_diagnostics(dg),
                   c("lme4 warned that the model failed to converge; in centred and scaled units the",
                     "  fit is less than 0.0001 standard errors from the maximum."))
})

test_that("the line is a tenth of a standard error, from both sides", {
  skip_on_cran()
  d <- year_data(1)
  inside <- maihda_fit_diagnostics(placed(below_0.0081, d))
  expect_lt(inside$stationary_distance, 0.1)
  expect_true(inside$converged)
  m <- placed(below_0.0121, d)
  outside <- maihda_fit_diagnostics(m)
  expect_gt(outside$stationary_distance, 0.1)
  expect_false(outside$converged)
  expect_identical(outside$messages, as.character(m@optinfo$conv$lme4$messages))
  expect_null(outside$lme4_overruled)
})

test_that("a Poisson slope fit at its maximum is converged", {
  skip_on_cran()
  # Seed 3 of the count design at its maximum, where lme4 reports max|grad| = 3.9.
  m <- placed(list(theta = c(7.3168189, -0.003708262, 0.00014962102),
                   beta = c(-16.074386, 0.0084192517)),
              year_data(3, "poisson"), yc ~ x + (1 + x | stratum), stats::poisson)
  expect_gt(length(m@optinfo$conv$lme4$messages), 0L)
  expect_gt(fitted_theta_distance(m), 1)
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$converged)
  expect_lt(dg$stationary_distance, 0.001)
})

test_that("fit_maihda() records a birth-year slope fit at its maximum as converged", {
  skip_on_cran()
  skip_slow_test()
  m <- quiet(fit_maihda(slope, data = year_data(1), family = "binomial"))
  dg <- m$diagnostics
  # lme4's checks fire on the fit fit_maihda() keeps, which is at the maximum.
  expect_gt(length(dg$lme4_overruled), 0L)
  expect_lt(deviance_of(m$model) - deviance_of(placed(max_1, year_data(1))), 1e-4)
  expect_true(dg$converged)
  expect_lt(dg$stationary_distance, 0.01)
  shown <- utils::capture.output(print(m))
  expect_true(any(grepl("lme4 warned that the model failed to converge", shown, fixed = TRUE)))
  expect_false(any(grepl("Convergence warnings reported by lme4", shown, fixed = TRUE)))
})

test_that("a fit at nAGQ = 0 is not measured, and measuring leaves the fit as lme4 made it", {
  skip_on_cran()
  d <- year_data(1)
  m0 <- quiet(lme4::glmer(slope, data = d, family = stats::binomial, nAGQ = 0))
  expect_true(is.na(maihda_lme4_stationary_distance(m0)))
  m <- placed(max_1, d)
  # Copies, not references: lme4 writes into some of the fit's vectors in place.
  state <- function(m) unserialize(serialize(list(
    theta = lme4::getME(m, "theta"), lambda = lme4::getME(m, "Lambdat")@x + 0,
    b = lme4::getME(m, "b"), u = lme4::getME(m, "u"), fixef = lme4::fixef(m),
    ranef = lme4::ranef(m), loglik = stats::logLik(m), predict = stats::predict(m),
    mu = m@resp$mu + 0), NULL))
  before <- state(m)
  first <- maihda_lme4_stationary_distance(m)
  expect_identical(state(m), before)
  dg <- maihda_fit_diagnostics(m)
  expect_identical(state(m), before)
  expect_identical(dg$stationary_distance, first)
})
