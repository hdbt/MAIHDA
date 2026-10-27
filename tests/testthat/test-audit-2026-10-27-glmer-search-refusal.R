# Audit pass 2026-10-27: glmer fits with a random slope outside growth models whose search was
# refused, left short of the maximum.
#
# fit_maihda() searches such a fit in re-expressed units on a deviance function rebuilt from a
# copy of the fit. That function starts its penalised least squares from the same point at every
# evaluation, and lme4 takes the Laplace deviance with the log-determinant of the weights before
# the last step, so at lme4's relative tolerance of 1e-7 its value at the fit's estimates can lie
# well above the converged one; lme4's own -2 log-likelihood, from the last step of its fit, can
# lie more than 1e-5 from it too. The re-expression was refused when the function at the fit
# missed lme4's -2 log-likelihood by more than 1e-5 relative: on simulated binomial null models
# with strata spread on the logit scale with an SD of 2 or more, fits were refused and left up to
# 13.6 below the maximum, some recorded as converged. Per the user, the searches now run that
# function to a relative tolerance of 1e-10, a fit is re-expressed when the deviance at its own
# re-expressed estimates matches that function at its estimates and that lies within 1e-3
# (relative) of lme4's -2 log-likelihood, and the singular-fit record's search does the same.
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
deviance_of <- function(m) -2 * as.numeric(stats::logLik(m))

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 2 on the logit scale) of 40 rows, a
# stratum slope on an N(0, 1) covariate (mean 0.3, SD 0.2), and a 0 / 1 outcome.
strong_strata <- function() {
  with_seed(5, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 2)[lev$d1] + stats::rnorm(3, 0, 2)[lev$d2] +
      stats::rnorm(2, 0, 2)[lev$d3]
    v <- stats::rnorm(18, 0, 0.2)
    st <- rep(seq_len(18), each = 40)
    x <- stats::rnorm(length(st))
    eta <- -0.2 + add[st] + (0.3 + v[st]) * x
    data.frame(y = stats::rbinom(length(st), 1, stats::plogis(eta)), x = x, stratum = factor(st))
  })
}

slope <- y ~ x + (1 + x | stratum)

# glmer's own fit, rounded to seven digits: at the maximum of the Laplace likelihood, which an
# unbounded search of glmer's deviance function (PIRLS to 1e-12) puts at 491.965239, with the
# slope's variance after the intercept at zero.
at_max <- list(theta = c(4.11175, 0.1080736, 0), beta = c(-1.756927, 0.235049))
deviance_max <- 491.965239
# The same with the intercept-slope entry's sign reversed: 1.0027 below the maximum.
flipped <- list(theta = c(4.11175, -0.1080736, 0), beta = c(-1.756927, 0.235049))
# The same with the intercept's entry at 0.8 of its value: 1.0308 below the maximum, where the
# fit's own -2 log-likelihood lies 2.7e-5 relative above the deviance function's converged value.
low <- list(theta = c(3.2894, 0.1080736, 0), beta = c(-1.756927, 0.235049))

# A glmer fit placed at a point: optimisers that return their start, and no move of an entry near
# zero onto zero (lme4's boundary.tol).
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}
stay_control <- function() {
  lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE, boundary.tol = 0)
}
placed <- function(point, data) {
  quiet(lme4::glmer(slope, data = data, family = stats::binomial,
                    start = list(theta = point$theta, fixef = point$beta), control = stay_control()))
}

# The fit's deviance function (from a detached copy) at its estimates, its penalised least squares
# run to `tolerance` (lme4's 1e-7 when NULL).
deviance_at <- function(m, tolerance = NULL) {
  f <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
  if (!is.null(tolerance)) environment(f)$tolPwrss <- tolerance
  f(c(unname(lme4::getME(m, "theta")), unname(lme4::fixef(m))))
}

# ---- the premise ---------------------------------------------------------------------------

test_that("the deviance function at lme4's tolerance misses a strong-strata fit's -2 log-likelihood", {
  skip_on_cran()
  d <- strong_strata()
  m <- placed(at_max, d)
  converged <- deviance_at(m, 1e-12)
  expect_lt(abs(converged / deviance_max - 1), 1e-8)
  # At lme4's tolerance, 2.5e-5 relative above lme4's -2 log-likelihood; at 1e-10, the converged
  # deviance, 6.7e-7 relative from lme4's.
  expect_gt(deviance_at(m) / deviance_of(m) - 1, 1e-5)
  expect_lt(abs(deviance_at(m, 1e-10) / converged - 1), 1e-9)
  expect_lt(abs(deviance_of(m) / converged - 1), 2e-6)
  # lme4's -2 log-likelihood of the fit placed lower lies 2.7e-5 relative above the converged one.
  m_low <- placed(low, d)
  expect_gt(deviance_of(m_low) / deviance_at(m_low, 1e-12) - 1, 1e-5)
})

# ---- the re-expression ---------------------------------------------------------------------

test_that("a fit is re-expressed when its deviance matches the function at its estimates", {
  skip_on_cran()
  d <- strong_strata()
  # Made at lme4's tolerance, though the function there misses lme4's -2 log-likelihood by more
  # than 1e-5; its deviance at the fit is the function's at the fit's estimates.
  m <- placed(flipped, d)
  re <- maihda_glmer_reexpressed(m)
  expect_false(is.null(re))
  expect_lt(abs(re$at / deviance_at(m) - 1), 1e-10)
  # Made at 1e-10 where lme4's -2 log-likelihood lies 2.7e-5 relative from the converged
  # deviance, which is the deviance at the fit.
  m_low <- placed(low, d)
  re <- maihda_glmer_reexpressed(m_low, tolerance = 1e-10)
  expect_false(is.null(re))
  expect_lt(abs(re$at / deviance_at(m_low, 1e-12) - 1), 1e-10)
  expect_lt(abs(re$deviance(re$p) / re$at - 1), 1e-12)
})

test_that("a fit is not re-expressed where the function is not the fit's or does not carry back", {
  skip_on_cran()
  d <- strong_strata()
  m <- placed(at_max, d)
  converged <- deviance_at(m, 1e-12)
  # lme4's -2 log-likelihood 5e-4 relative from the function at the fit's estimates: re-expressed;
  # 2e-3: not, as a deviance function that is not the fit's.
  local_mocked_bindings(maihda_glmer_criterion = function(model) converged * (1 + 5e-4))
  expect_false(is.null(maihda_glmer_reexpressed(m, tolerance = 1e-10)))
  local_mocked_bindings(maihda_glmer_criterion = function(model) converged * (1 + 2e-3))
  expect_null(maihda_glmer_reexpressed(m, tolerance = 1e-10))
  local_mocked_bindings(maihda_glmer_criterion = function(model) converged * (1 - 2e-3))
  expect_null(maihda_glmer_reexpressed(m, tolerance = 1e-10))
})

test_that("a fit is not re-expressed where its estimates do not carry back to themselves", {
  skip_on_cran()
  d <- strong_strata()
  m <- placed(flipped, d)
  local_mocked_bindings(maihda_glmer_carry_theta = function(theta, cnms, Ts) theta * 1.1)
  expect_null(maihda_glmer_reexpressed(m, tolerance = 1e-10))
})

# ---- the search refit ----------------------------------------------------------------------

test_that("maihda_fit_glmer_bound() runs the searches' deviance function to 1e-10", {
  seen <- NULL
  env <- new.env()
  env$fake_fit <- function() list(crit = 10)
  local_mocked_bindings(
    maihda_glmer_mirror_refit = function(...) NULL,
    maihda_glmer_search_refit = function(kept, fit_call, fit_env, tol = 1e-4, tolerance = NULL,
                                         gain = NULL) {
      seen <<- list(crit = kept$value$crit, tol = tol, tolerance = tolerance, gain = gain)
      NULL
    }
  )
  expect_identical(maihda_fit_glmer_bound(quote(fake_fit()), env)$crit, 10)
  # The tolerance given, lme4's PIRLS to 1e-10, and a refit for any gain above `tol`.
  expect_identical(seen, list(crit = 10, tol = 1e-4, tolerance = 1e-10, gain = NULL))
})

test_that("a fit stopped short whose function at lme4's tolerance missed its -2 log-likelihood is refitted", {
  skip_on_cran()
  skip_slow_test()
  env <- new.env()
  env$data <- strong_strata()
  env$start <- list(theta = flipped$theta, fixef = flipped$beta)
  env$control <- stay_control()
  call <- quote(lme4::glmer(y ~ x + (1 + x | stratum), data = data, family = stats::binomial,
                            start = start, control = control))
  # The first fit is placed 1.0027 below the maximum; the search from it reaches the maximum,
  # which the refit from its lower end keeps.
  m <- quiet(maihda_fit_glmer_bound(call, env))
  expect_lt(deviance_at(m, 1e-12) - deviance_max, 1e-4)
  expect_identical(m@optinfo$optimizer, "bobyqa")
})

# ---- the singular-fit record ---------------------------------------------------------------

test_that("the singular-fit record's search runs its deviance function to 1e-10", {
  skip_on_cran()
  d <- strong_strata()
  # A fit given its own start is not refitted from the searches, so the record searches it when
  # lme4 calls it singular: 1.0027 below the maximum, the search reaches 0.01 lower or more.
  fit <- quiet(fit_maihda(slope, data = d, family = "binomial",
                          start = list(theta = flipped$theta, fixef = flipped$beta),
                          control = stay_control()))
  expect_true(fit$diagnostics$singular)
  expect_false(fit$diagnostics$converged)
  expect_gte(fit$diagnostics$search_gain, 0.01)
  expect_lte(fit$diagnostics$search_gain, 1.0027 + 1e-4)
  # The record's deviance function at 1e-10, and its check of the function against lme4's
  # -2 log-likelihood at 1e-3: a fit whose -2 log-likelihood lies 2.7e-5 from the converged
  # deviance is searched, one at the maximum gains nothing.
  expect_gte(maihda_glmer_search_gain(placed(low, d)), 0.01)
  m <- placed(at_max, d)
  expect_identical(maihda_glmer_search_gain(m), 0)
  used <- NULL
  real <- maihda_glmer_search_with
  local_mocked_bindings(maihda_glmer_search_with = function(model, deviance, within) {
    used <<- environment(deviance)$tolPwrss
    real(model, deviance, within)
  })
  maihda_glmer_search_gain(m)
  expect_identical(used, 1e-10)
  # The fit's -2 log-likelihood 2e-3 relative from the function at its estimates: no search.
  off <- m
  off@devcomp$cmp[["dev"]] <- deviance_at(m, 1e-12) * (1 + 2e-3)
  expect_lt(abs(deviance_of(off) / deviance_at(m, 1e-12) - 1 - 2e-3), 1e-9)
  expect_true(is.na(maihda_glmer_search_gain(off)))
})
