# Audit pass 2026-10-25: lmer fits that lme4 runs no convergence check on.
#
# lmer() computes the derivatives its convergence checks are taken from only for a fit of fewer
# rows and variance parameters than lmerControl()'s check.conv.nobsmax and check.conv.nparmax
# (10,000 and 10 by default) that it does not call singular. For any other fit it records no
# convergence check at all, so the record had nothing to go on: on simulated designs of 10,000 to
# 30,000 rows, lmer's own fits with a random slope (what fit_maihda() returns when given lme4's
# starting values as its start, which skips its refits) stopped up to 1,236 above the REML optimum
# and were recorded as converged. Per the user, a cross-sectional lmer fit with a random slope that
# lme4 did not check, and that fit_maihda() did not refit, is searched as a fit lme4 calls singular
# is, and counts as not converged when a search reaches a criterion 0.01 or more below it.
# Intercept-only and growth fits are not searched.
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
# calendar year between 2010 and 2020, a stratum slope deviation of SD 0.05 per year around 2015
# and a residual SD of 1.
calendar_data <- function(seed, n) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, 0.05)
    st <- rep(seq_len(18), each = n)
    x <- sample(2010:2020, length(st), replace = TRUE)
    y <- 1 + add[st] + (0.1 + v[st]) * (x - 2015) + stats::rnorm(length(st))
    data.frame(y = y, x = x, stratum = factor(st))
  })
}

# Eighteen strata of 40 rows with a birth year between 1940 and 2000, a stratum slope deviation of
# SD 0.005 per year and a residual SD of 1.
year_data <- function(seed) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    u <- stats::rnorm(18, 0, 0)
    v <- stats::rnorm(18, 0, 0.005)
    st <- rep(seq_len(18), each = 40)
    x <- sample(1940:2000, 720, replace = TRUE)
    y <- 1 + add[st] + u[st] + (0.3 + v[st]) * (x - 1970) + stats::rnorm(720)
    data.frame(y = y, x = x, stratum = factor(st))
  })
}

slope <- y ~ x + (1 + x | stratum)

# REML optima from lmer's deviance function searched with no bound on the Cholesky diagonal:
# - calendar_data(12, 556), 10,008 rows: 28469.4621737 at theta (79.034429, -0.0391618528,
#   0.000228845); lmer's own fit stops 106.45 above at (1.779906, -0.000868207, 0.0002464972),
#   the strata's slope variance nearly all missed, and is not singular.
# - year_data(3), 720 rows: 1994.85420 at (18.347146, -0.0093008294, 0.00017077300); a point
#   along one direction of the standardised theta lies 0.03 above it.
optimum_12 <- 28469.4621737
at_optimum_12 <- c(79.034429, -0.0391618528, 0.000228845)
stuck_12 <- c(1.779906, -0.000868207, 0.0002464972)
optimum_3 <- 1994.85420
at_optimum_3 <- c(18.347146, -0.0093008294, 0.00017077300)
above_0.03 <- c(18.987570, -0.0096323374, 0.00017433472)

# An lmer fit placed at `theta`: an optimiser that returns its start, and no move of an entry
# near zero onto zero (lme4's boundary.tol). `...` goes to lmerControl().
stay <- function(fn, par, lower, upper, control) {
  list(par = par, fval = fn(par), conv = 0, feval = 1, message = "stay")
}
placed <- function(theta, data, formula = slope, REML = TRUE, ...) {
  quiet(lme4::lmer(formula, data = data, REML = REML, start = theta,
                   control = lme4::lmerControl(optimizer = stay, boundary.tol = 0, ...)))
}

# ---- which fits lme4 does not check -------------------------------------------------------

test_that("lme4 runs no convergence check on an lmer fit of 10,000 rows or with 10 variance parameters", {
  skip_on_cran()
  # 10,008 rows: no derivatives and no convergence entry for a fit lme4 does not call singular.
  big <- placed(stuck_12, calendar_data(12, 556))
  expect_false(lme4::isSingular(big))
  expect_null(big@optinfo$derivs)
  expect_null(big@optinfo$conv$lme4)
  expect_true(maihda_lme4_unchecked(big))
  # With the derivatives asked for, lme4 computes them and still checks nothing above 10,000 rows.
  asked <- placed(stuck_12, calendar_data(12, 556), calc.derivs = TRUE)
  expect_false(is.null(asked@optinfo$derivs))
  expect_true(maihda_lme4_unchecked(asked))
  # 9,990 rows of the same design: lme4 checks the same point and its checks fire.
  small <- placed(stuck_12, calendar_data(12, 555))
  expect_false(lme4::isSingular(small))
  expect_gt(length(small@optinfo$conv$lme4$messages), 0L)
  expect_false(maihda_lme4_unchecked(small))
  # Ten variance parameters on 720 rows, and lme4's own switches.
  d <- year_data(3)
  d$z <- with_seed(5, stats::rnorm(720))
  d$w <- with_seed(6, stats::rnorm(720))
  ten <- placed(c(1, 0, 0, 1, 0, 1), d, y ~ x + z + w + (1 + z + w | stratum))
  expect_false(maihda_lme4_unchecked(ten))
  d$v <- with_seed(7, stats::rnorm(720))
  ten <- placed(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 1), d, y ~ x + z + w + (1 + z + w + v | stratum))
  expect_length(lme4::getME(ten, "theta"), 10L)
  expect_true(maihda_lme4_unchecked(ten))
  expect_true(maihda_lme4_unchecked(placed(above_0.03, d, calc.derivs = FALSE)))
  expect_true(maihda_lme4_unchecked(placed(above_0.03, d, check.conv.nobsmax = 100)))
  expect_false(maihda_lme4_unchecked(placed(above_0.03, d)))
  # A fit lme4 calls singular carries its message, and so does a glmer fit lme4 checked.
  singular <- placed(c(above_0.03[1], above_0.03[2], 0), d)
  expect_true(lme4::isSingular(singular))
  expect_false(maihda_lme4_unchecked(singular))
  d$yb <- as.integer(d$y > stats::median(d$y))
  expect_false(maihda_lme4_unchecked(quiet(lme4::glmer(yb ~ x + (1 | stratum), data = d,
                                                       family = stats::binomial))))
  # A random slope is any term other than an intercept.
  expect_true(maihda_lme4_has_slope(big))
  expect_true(maihda_lme4_has_slope(placed(c(1, 1), d, y ~ x + (1 | stratum) + (0 + x | stratum))))
  expect_false(maihda_lme4_has_slope(placed(1, d, y ~ x + (1 | stratum))))
})

# ---- the record ---------------------------------------------------------------------------

test_that("an unchecked lmer fit with a random slope that a search reaches 0.01 below is recorded as not converged, by how much", {
  skip_on_cran()
  d <- year_data(3)
  m <- placed(above_0.03, d, calc.derivs = FALSE)
  gap <- lme4::REMLcrit(m) - optimum_3
  expect_lt(abs(gap - 0.03), 1e-4)
  dg <- maihda_fit_diagnostics(m)
  expect_false(dg$converged)
  expect_lt(abs(dg$search_gain / gap - 1), 1e-3)
  expect_false(isTRUE(dg$singular))
  expect_identical(dg$messages, paste0(
    "lme4 ran no convergence check on this fit (by default it runs none from 10,000 rows or 10 ",
    "variance parameters on), and a search of the same model with each random-effect block's ",
    "terms centred and scaled reaches a REML criterion 0.03 lower: the fit stopped short of the ",
    "optimum."))
  lines <- maihda_format_fit_diagnostics(dg)
  expect_identical(lines[1], "Convergence problems of the lme4 fit:")
  expect_false(any(grepl("Singular fit", lines, fixed = TRUE)))
  # At the optimum nothing is lower: converged, no field, nothing printed.
  at <- maihda_fit_diagnostics(placed(at_optimum_3, d, calc.derivs = FALSE))
  expect_true(at$converged)
  expect_null(at$search_gain)
  expect_identical(maihda_format_fit_diagnostics(at), character(0))
  # By maximum likelihood the fit's criterion is its deviance.
  ml <- maihda_fit_diagnostics(placed(above_0.03, d, REML = FALSE, calc.derivs = FALSE))
  expect_false(ml$converged)
  expect_true(grepl("reaches a deviance 0.15 lower", ml$messages, fixed = TRUE))
  # A fit lme4 calls singular keeps its own sentence.
  expect_true(grepl("^lme4 runs no convergence check on a singular fit, ",
                    maihda_lmer_search_message(placed(c(above_0.03[1:2], 0), d), 0.5)))
})

test_that("the record's line for an unchecked fit is 0.01", {
  skip_on_cran()
  m <- placed(above_0.03, year_data(3), calc.derivs = FALSE)
  at_gain <- function(gain) {
    local_mocked_bindings(maihda_lmer_search_gain = function(model) gain)
    maihda_fit_diagnostics(m)
  }
  expect_false(at_gain(0.01)$converged)
  expect_identical(at_gain(0.01)$search_gain, 0.01)
  below <- at_gain(0.01 - 1e-7)
  expect_true(below$converged)
  expect_null(below$search_gain)
  expect_true(at_gain(NA_real_)$converged)
})

test_that("only an unchecked lmer fit with a random slope outside growth models, not refitted by fit_maihda(), is searched", {
  skip_on_cran()
  searched <- 0L
  local_mocked_bindings(maihda_lmer_search_gain = function(model) {
    searched <<- searched + 1L
    0.5
  })
  d <- year_data(3)
  m <- placed(above_0.03, d, calc.derivs = FALSE)
  expect_false(maihda_fit_diagnostics(m)$converged)
  expect_identical(searched, 1L)
  # Refitted by fit_maihda() from the same searches, or a growth model: not searched.
  expect_true(maihda_fit_diagnostics(m, searched = TRUE)$converged)
  expect_true(maihda_fit_diagnostics(m, longitudinal_info = list(id = "id"))$converged)
  # Checked by lme4, whether its checks find nothing or fire: not searched.
  sleep <- quiet(lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy))
  expect_identical(sleep@optinfo$conv$lme4, list())
  expect_true(maihda_fit_diagnostics(sleep)$converged)
  fired <- placed(at_optimum_3, d)
  expect_gt(length(fired@optinfo$conv$lme4$messages), 0L)
  expect_true(maihda_fit_diagnostics(fired)$converged)
  expect_identical(searched, 1L)
  # Unchecked without a random slope: not searched.
  expect_true(maihda_fit_diagnostics(placed(1, d, y ~ x + (1 | stratum), calc.derivs = FALSE))$converged)
  # An optimiser failure or a boundary trap has decided already.
  failed <- m
  failed@optinfo$conv$opt <- 1L
  expect_false(maihda_fit_diagnostics(failed)$converged)
  expect_identical(searched, 1L)
  # Ten variance parameters below 10,000 rows: searched.
  d$z <- with_seed(5, stats::rnorm(720))
  d$w <- with_seed(6, stats::rnorm(720))
  d$v <- with_seed(7, stats::rnorm(720))
  ten <- placed(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 1), d, y ~ x + z + w + (1 + z + w + v | stratum))
  expect_false(maihda_fit_diagnostics(ten)$converged)
  expect_identical(searched, 2L)
  local_mocked_bindings(maihda_lme4_boundary_traps = function(model, ...) {
    data.frame(block = "stratum", term = "(Intercept)", index = 1L, drop = 1, n_terms = 2L,
               stringsAsFactors = FALSE)
  })
  expect_false(maihda_fit_diagnostics(m)$converged)
  expect_identical(searched, 2L)
})

# ---- live fits ----------------------------------------------------------------------------

test_that("a fit of 10,000 rows made from a start is searched, and fit_maihda()'s refitted fit is not", {
  skip_slow_test()
  skip_on_cran()
  d <- calendar_data(12, 556)
  gap <- lme4::REMLcrit(placed(stuck_12, d)) - optimum_12
  expect_gt(gap, 100)
  # fit_maihda() given a start makes the fit lmer makes from it, here placed where lmer's own fit
  # stops; nothing in lme4 looks at it.
  m <- quiet(fit_maihda(slope, data = d, start = stuck_12,
                        control = lme4::lmerControl(optimizer = stay, boundary.tol = 0)))
  expect_null(m$model@optinfo$conv$lme4)
  expect_false(m$diagnostics$converged)
  expect_lt(abs(m$diagnostics$search_gain / gap - 1), 1e-3)
  shown <- utils::capture.output(print(m))
  expect_true(any(grepl("lme4 ran no convergence check on this fit", shown, fixed = TRUE)))
  # At the optimum: converged, no field.
  at <- maihda_fit_diagnostics(placed(at_optimum_12, d))
  expect_true(at$converged)
  expect_null(at$search_gain)
  # Refitted by fit_maihda(), the fit reaches the optimum and is not searched again.
  calls <- 0L
  real <- maihda_lmer_search_gain
  local_mocked_bindings(maihda_lmer_search_gain = function(model) {
    calls <<- calls + 1L
    real(model)
  })
  own <- quiet(fit_maihda(slope, data = d))
  expect_lt(abs(lme4::REMLcrit(own$model) - optimum_12), 1e-4)
  expect_true(own$diagnostics$converged)
  expect_identical(calls, 0L)
  # On 9,990 rows lme4 checks the same point itself, and the record keeps its verdict.
  small <- maihda_fit_diagnostics(placed(stuck_12, calendar_data(12, 555)))
  expect_false(small$converged)
  expect_null(small$search_gain)
  expect_identical(calls, 0L)
})
