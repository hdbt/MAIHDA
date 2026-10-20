# Audit pass 2026-10-19: the refit of a glmer growth fit across a variance bound, which could
# stop short of the maximum.
#
# A glmer growth fit kept at a bound of lme4's parametrisation -- a block's intercept
# Cholesky entry at 0 beside a non-zero entry below it -- is refitted once from the other
# side of that bound by glmer's second stage alone. Made by Nelder-Mead, that refit could
# stop short of a maximum where the stratum block has rank one (its intercept and slope
# perfectly correlated), along which the likelihood is flat: 0.023 below on the seed-30
# design in weeks with the rows in another order, recorded as not converged, and 0.005 below,
# recorded as converged, from that fit's estimates. The refit is now made by bobyqa, as the
# refit of other glmer fits is, and a refit that is kept is refitted once more: across the
# bound again when it stopped at one, as bobyqa can, and otherwise from where it stopped,
# as bobyqa can stop on the bound of a rank-one block beside a maximum just inside it.
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

heard <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(expr,
    message = function(m) { seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { seen <<- c(seen, paste("warning:", conditionMessage(w))); invokeRestart("muffleWarning") })
  list(value = value, seen = seen)
}

# test-audit-2026-10-18-glmer-growth-optimum.R's design: 16 factorial strata of 6 persons,
# waves 0..5 times k, stratum and person intercepts and slopes, a Poisson count n.
growth_design <- function(seed, k = 1) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
    add0 <- 0.3 * (g$A - 0.5) + 0.2 * (g$B - 0.5) - 0.25 * (g$C - 0.5) + 0.15 * (g$D - 0.5)
    add1 <- 0.04 * (g$A - 0.5) - 0.03 * (g$B - 0.5) + 0.05 * (g$C - 0.5) - 0.02 * (g$D - 0.5)
    u <- cbind(stats::rnorm(16, 0, 0.2), stats::rnorm(16, 0, 0.06), stats::rnorm(16, 0, 0))
    p <- cbind(stats::rnorm(96, 0, 0.5), stats::rnorm(96, 0, 0.15), stats::rnorm(96, 0, 0))
    id <- rep(seq_len(96), each = 6)
    st <- rep(seq_len(16), each = 6)[id]
    t <- rep(0:5, times = 96)
    eta <- add0[st] + u[st, 1] + (add1[st] + u[st, 2]) * t + p[id, 1] + p[id, 2] * t
    y <- stats::rnorm(576, 0, 0.5)
    data.frame(id = id, wave = t * k,
               A = factor(g$A[st]), B = factor(g$B[st]), C = factor(g$C[st]), D = factor(g$D[st]),
               n = stats::rpois(576, exp(0.5 + 0.1 * t + 0.6 * eta)))
  })
}

# A stand-in fitting function: a first fit, then a restart per start, each with a criterion.
fake_fits <- function(first_crit, restart_crit) {
  env <- new.env()
  env$fake_fit <- function(start = NULL, control = NULL) {
    if (is.null(start)) {
      list(crit = first_crit)
    } else {
      list(crit = restart_crit[[as.character(start$theta[["a"]])]])
    }
  }
  env
}

# fit_maihda() on a Poisson growth design with its first fit placed at `point` (theta and
# fixef) by an optimiser that returns its start, and no restarts, so that the fit kept for the
# refit is that point on every platform. Returns the fit, the fit the first refit across a
# bound started from, and the refit the second refit started from.
fit_from_point <- function(d, point) {
  stay <- function(par, fn, lower, upper, control = list(), ...) {
    list(par = par, fval = fn(par), conv = 0, message = "the start", feval = 1)
  }
  kept <- NULL
  refit <- NULL
  real_growth <- maihda_fit_glmer_growth
  real_mirror <- maihda_glmer_mirror_refit
  real_pass <- maihda_glmer_refit_pass
  local_mocked_bindings(
    maihda_glmer_growth_starts = function(...) list(),
    maihda_fit_glmer_growth = function(fit_call, fit_env, id, ...) {
      first_call <- fit_call
      first_call$start <- point
      first_call$control <- lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE)
      real_growth(fit_call, fit_env, id, first = maihda_eval_holding_conditions(first_call, fit_env), ...)
    },
    maihda_glmer_mirror_refit = function(kept_fit, ...) {
      if (is.null(kept)) kept <<- kept_fit$value
      real_mirror(kept_fit, ...)
    },
    maihda_glmer_refit_pass = function(refit_fit, ...) {
      refit <<- refit_fit$value
      real_pass(refit_fit, ...)
    }
  )
  m <- quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                        family = "poisson"))
  list(m = m, kept = kept, refit = refit)
}

# The stratum block of a growth fit's theta as its Cholesky entries, and its correlation.
stratum_chol <- function(model) {
  th <- lme4::getME(model, "theta")
  unname(th[c("stratum.(Intercept)", "stratum..maihda_ctime.(Intercept)", "stratum..maihda_ctime")])
}
chol_cor <- function(L) L[1] * L[2] / sqrt(L[1]^2 * (L[2]^2 + L[3]^2))

# ---- the refits, with stand-in fits -------------------------------------------------------

test_that("a growth fit kept at a bound is refitted across it by bobyqa, then once more", {
  across <- NULL
  passed <- NULL
  local_mocked_bindings(
    maihda_glmer_growth_starts = function(theta, cnms, id, scale = 0.5) list(c(a = 0.5), c(a = 2)),
    maihda_glmer_criterion = function(model) model$crit,
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4, optimizer = NULL) {
      across <<- list(kept$value$crit, tol, optimizer)
      if (kept$value$crit > 7) list(value = list(crit = 6), conditions = list(simpleMessage("across\n"))) else NULL
    },
    maihda_glmer_refit_pass = function(refit, fit_call, fit_env, tol = 1e-4) {
      passed <<- list(refit$value$crit, tol)
      if (refit$value$crit > 5.5) list(value = list(crit = 5), conditions = list(simpleMessage("pass\n"))) else NULL
    }
  )
  # From the fit the restarts kept, by bobyqa, with the tolerance given; then the pass from
  # the refit, and only the pass speaks.
  env <- fake_fits(10, list("0.5" = 9.5, "2" = 8))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid", tol = 0.5))
  expect_identical(across, list(8, 0.5, "bobyqa"))
  expect_identical(passed, list(6, 0.5))
  expect_identical(h$value$crit, 5)
  expect_identical(h$seen, "pass\n")
  # A bootstrap draw's refit given as the first fit: the same refits.
  across <- NULL
  passed <- NULL
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid",
                                     first = list(value = list(crit = 7.5), conditions = list())))
  expect_identical(across, list(7.5, 1e-4, "bobyqa"))
  expect_identical(h$value$crit, 5)
  # No refit across the bound: no pass, and the kept fit stands.
  across <- NULL
  passed <- NULL
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(7, list("0.5" = 9, "2" = 8)), "pid"))
  expect_identical(across[[1]], 7)
  expect_null(passed)
  expect_identical(h$value$crit, 7)
  # A pass that is not returned leaves the refit, and a pass that fails is passed over.
  local_mocked_bindings(maihda_glmer_refit_pass = function(refit, fit_call, fit_env, tol = 1e-4) NULL)
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))
  expect_identical(h$value$crit, 6)
  expect_identical(h$seen, "across\n")
  local_mocked_bindings(maihda_glmer_refit_pass = function(...) stop("the pass fails"))
  expect_identical(heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))$value$crit, 6)
  expect_identical(formals(maihda_fit_glmer_growth)$tol, 1e-4)
})

test_that("the pass is glmer's second stage from the refit's estimates, with the refit's control", {
  skip_on_cran()
  env <- new.env()
  env$data <- growth_design(30, 7)
  env$data$t <- env$data$wave / 10
  call <- quote(lme4::glmer(n ~ t + (t | id) + (t | A:B:C:D), data = data,
                            family = stats::poisson))
  ctrl <- lme4::glmerControl(optimizer = "bobyqa", nAGQ0initStep = FALSE,
                             optCtrl = list(maxfun = 1e5))
  assign(".maihda_arg_refit_control", ctrl, envir = env)
  # A refit stopped beside the maximum: its stratum block at (0.0538, -0.0904, 0.0135).
  at <- call
  at$start <- list(theta = c(0.412, -0.0278, 0.1366, 0.0538, -0.0904, 0.0135),
                   fixef = c(0.4187, 0.0909))
  at$control <- lme4::glmerControl(optimizer = list(
    function(par, fn, lower, upper, control = list(), ...) list(par = par, fval = fn(par), conv = 0, feval = 1),
    function(par, fn, lower, upper, control = list(), ...) list(par = par, fval = fn(par), conv = 0, feval = 1)),
    nAGQ0initStep = FALSE)
  refit <- quiet(maihda_eval_holding_conditions(at, env))
  seen <- NULL
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen <<- call
    real(call, env)
  })
  out <- quiet(maihda_glmer_refit_pass(refit, call, env))
  expect_equal(seen$start$theta, unname(lme4::getME(refit$value, "theta")))
  expect_identical(seen$start$fixef, lme4::fixef(refit$value))
  expect_identical(seen$control, as.name(".maihda_arg_refit_control"))
  # It moves on to the maximum, -1021.05075, and is returned.
  expect_false(is.null(out))
  expect_lt(abs(as.numeric(stats::logLik(out$value)) - (-1021.05075)), 0.005)
  # A pass that is not lower by more than the tolerance is not returned.
  expect_null(quiet(maihda_glmer_refit_pass(out, call, env)))
  expect_null(quiet(maihda_glmer_refit_pass(refit, call, env, tol = 1)))
  # A refit stopped at a bound again is refitted across it once more, by bobyqa, with the
  # tolerance given, and not from its own estimates.
  at$start <- list(theta = c(0.412581, -0.030526, 0.139232, 0, 0.055931, 0.051393),
                   fixef = c(0.420019, 0.090793))
  at_bound <- quiet(real(at, env))
  across <- NULL
  local_mocked_bindings(
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4, optimizer = NULL) {
      across <<- list(kept$value, fit_call, tol, optimizer)
      "across"
    },
    maihda_eval_holding_conditions = function(call, env) stop("a pass from its own estimates")
  )
  expect_identical(maihda_glmer_refit_pass(at_bound, call, env, tol = 0.5), "across")
  expect_identical(across, list(at_bound$value, call, 0.5, "bobyqa"))
})

# ---- live fits ----------------------------------------------------------------------------

test_that("the refit from a stratum bound reaches a maximum where the block has rank one", {
  skip_on_cran()
  # Seed 30 in weeks, from the stratum-bound fit the restarts kept with the rows in another
  # order (with_seed(20, sample.int(576))), its estimates rounded to six digits. At the maximum
  # of -1021.05075 the stratum intercept and slope are perfectly correlated. From this point
  # glmer's Nelder-Mead stopped beside it, at a correlation of -0.989 and 0.005 below,
  # recorded as converged; bobyqa reaches it.
  point <- list(theta = c(0.412583, -0.030527, 0.139233, 0, -0.002480, 0.075917),
                fixef = c(0.420017, 0.090792))
  f <- fit_from_point(growth_design(30, 7), point)
  # The premise: the fit kept for the refit is that point, at the stratum bound.
  expect_identical(unname(lme4::getME(f$kept, "theta")), point$theta)
  expect_lt(chol_cor(stratum_chol(f$m$model)), -0.999)
  expect_lt(abs(as.numeric(stats::logLik(f$m$model)) - (-1021.05075)), 0.01)
  expect_true(f$m$diagnostics$converged)
})

test_that("the refit from a person-level bound reaches the maximum where Nelder-Mead stops short", {
  skip_on_cran()
  # Seed 3 in waves, from its maximum of -1082.00427 with the person intercept's Cholesky
  # entry set to 0 (the other estimates rounded to six digits). From the other side of that
  # bound glmer's Nelder-Mead stopped 0.032 below the maximum, and a second Nelder-Mead pass
  # stayed there, with the rows in any of four orders; bobyqa reaches it.
  point <- list(theta = c(0, -0.037739, 0.127576, 0.0319, -0.029645, 3.1e-05),
                fixef = c(0.588719, 0.073881))
  f <- fit_from_point(growth_design(3), point)
  th <- lme4::getME(f$kept, "theta")
  expect_identical(unname(th), point$theta)
  expect_identical(vapply(maihda_lme4_trap_columns(th, lme4::getME(f$kept, "cnms")), `[[`,
                          character(1), "block"), "id")
  expect_lt(abs(as.numeric(stats::logLik(f$m$model)) - (-1082.00427)), 0.005)
  expect_true(f$m$diagnostics$converged)
})

test_that("a refit that stops at a bound again is refitted across it once more", {
  skip_on_cran()
  # Seed 45 in weeks, from glmer's own first fit (estimates rounded to six digits), the
  # person intercept's Cholesky entry at 0. With the rows as they are, and in two of three
  # other orders, the bobyqa refit from across that bound stops with the stratum intercept's
  # entry at 0 instead, 0.136 below the maximum of -1081.44047, recorded as not converged;
  # the second refit, across that bound, reaches it.
  point <- list(theta = c(0, -0.167315, 0.083878, 0.006261, 0.02037, 3.7e-05),
                fixef = c(0.600034, 0.111866))
  f <- fit_from_point(growth_design(45, 7), point)
  expect_identical(unname(lme4::getME(f$kept, "theta")), point$theta)
  expect_lt(abs(as.numeric(stats::logLik(f$m$model)) - (-1081.44047)), 0.005)
  expect_length(maihda_lme4_trap_columns(lme4::getME(f$m$model, "theta"),
                                         lme4::getME(f$m$model, "cnms")), 0L)
  expect_true(f$m$diagnostics$converged)
})
