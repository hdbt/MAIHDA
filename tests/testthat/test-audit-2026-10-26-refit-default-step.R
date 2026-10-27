# Audit pass 2026-10-26: glmer fits that fit_maihda() refits by bobyqa outside growth models.
#
# lme4 records the optimiser that made a fit, and lme4::refit() of a fit without a control runs that
# optimiser at its package's defaults: for bobyqa a first step of min(0.95, 0.2 max |par|), which
# moves every entry within it of its bound up to it before it starts. From the refits of
# maihda_fit_glmer_bound() (across a bound, from the searches), recorded as bobyqa, refit() ended up
# to 188 above the fit in -2 log-likelihood. Per the user, a fit that either refit replaced is
# refitted once more from its estimates by the second stage of the call's own control (Nelder-Mead
# by default), with first steps of 0.1% of each estimate where Nelder-Mead's own fail or it stops
# without converging, so the fit records that optimiser as glmer's own fits do.
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

# Twelve groups of 25 rows, a covariate on a scale of 1 and a Poisson count: a fit glmer makes
# in a fraction of a second.
small_counts <- function() {
  with_seed(11, {
    g <- rep(seq_len(12), each = 25)
    x <- stats::rnorm(300)
    u <- stats::rnorm(12, 0, 0.4)
    v <- stats::rnorm(12, 0, 0.2)
    data.frame(y = stats::rpois(300, exp(0.3 + u[g] + (0.2 + v[g]) * x)), x = x, g = factor(g))
  })
}
small_call <- quote(lme4::glmer(y ~ x + (1 + x | g), data = data, family = stats::poisson))

# A value / conditions pair of a fit made by bobyqa in both stages, as a refit of
# maihda_fit_glmer_bound() records it.
bobyqa_fit <- function(env) {
  call <- small_call
  call$control <- lme4::glmerControl(optimizer = c("bobyqa", "bobyqa"))
  quiet(maihda_eval_holding_conditions(call, env))
}
small_env <- function() {
  env <- new.env()
  env$data <- small_counts()
  env
}

# ---- the pass ------------------------------------------------------------------------------

test_that("a refitted fit is refitted once more by the second stage of the call's own control", {
  skip_on_cran()
  env <- small_env()
  kept <- bobyqa_fit(env)
  expect_identical(kept$value@optinfo$optimizer, "bobyqa")
  seen <- list()
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen[[length(seen) + 1L]] <<- call
    real(call, env)
  })
  out <- quiet(maihda_glmer_finish(kept, small_call, env))
  # From the fit's estimates, glmer's second stage alone, the call's own optimiser pair.
  expect_length(seen, 1L)
  expect_identical(seen[[1]]$start, list(theta = unname(lme4::getME(kept$value, "theta")),
                                         fixef = lme4::fixef(kept$value)))
  expect_identical(seen[[1]]$control, as.name(".maihda_arg_finish_control"))
  ctrl <- get(".maihda_arg_finish_control", envir = env)
  expect_identical(ctrl$optimizer, c("bobyqa", "Nelder_Mead"))
  expect_false(ctrl$nAGQ0initStep)
  expect_identical(ctrl$optCtrl, list())
  # The fit records Nelder-Mead, as glmer's own fit does, and is no worse on one deviance
  # function; the two fits' criteria differ by the precision of lme4's PIRLS.
  expect_identical(out$value@optinfo$optimizer, "Nelder_Mead")
  expect_lt(maihda_glmer_rise(out$value, kept$value), 1e-4)
  expect_lt(abs(deviance_of(out$value) - deviance_of(kept$value)), 1e-2)
  own <- quiet(eval(small_call, env))
  expect_identical(out$value@optinfo$optimizer, own@optinfo$optimizer)
  # A second-stage optimiser of the user's own is the one the pass takes.
  call <- small_call
  call$control <- quote(user_control)
  env$user_control <- lme4::glmerControl(optimizer = c("bobyqa", "nloptwrap"), tolPwrss = 1e-8)
  out <- quiet(maihda_glmer_finish(kept, call, env))
  expect_identical(out$value@optinfo$optimizer, "nloptwrap")
  expect_identical(get(".maihda_arg_finish_control", envir = env)$tolPwrss, 1e-8)
})

test_that("a fit that records the call's own optimiser, or a control naming a function, has no pass", {
  skip_on_cran()
  env <- small_env()
  kept <- bobyqa_fit(env)
  calls <- 0L
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    calls <<- calls + 1L
    real(call, env)
  })
  call <- small_call
  call$control <- quote(user_control)
  env$user_control <- lme4::glmerControl(optimizer = "bobyqa")
  expect_null(maihda_glmer_finish(kept, call, env))
  stay <- function(par, fn, lower, upper, control = list(), ...) list(par = par, fval = fn(par), conv = 0)
  env$user_control <- lme4::glmerControl(optimizer = list(stay, stay))
  expect_null(maihda_glmer_finish(kept, call, env))
  expect_identical(calls, 0L)
})

test_that("where Nelder-Mead's own first steps fail, the pass takes first steps of 0.1% of each estimate", {
  skip_on_cran()
  env <- small_env()
  kept <- bobyqa_fit(env)
  # The same fit with its intercept-slope entry placed at 0.
  stay <- function(par, fn, lower, upper, control = list(), ...) {
    list(par = par, fval = fn(par), conv = 0, feval = 1)
  }
  zero <- small_call
  zero$start <- list(theta = unname(lme4::getME(kept$value, "theta")) * c(1, 0, 1),
                     fixef = lme4::fixef(kept$value))
  zero$control <- lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE,
                                     boundary.tol = 0)
  at_zero <- quiet(maihda_eval_holding_conditions(zero, env))
  controls <- list()
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    controls[[length(controls) + 1L]] <<- get(".maihda_arg_finish_control", envir = env)$optCtrl
    if (length(controls) == 1L) stop("Downdated VtV is not positive definite")
    real(call, env)
  })
  out <- quiet(maihda_glmer_finish(kept, small_call, env))
  expect_length(controls, 2L)
  expect_identical(controls[[1]], list())
  par <- c(unname(lme4::getME(kept$value, "theta")), unname(lme4::fixef(kept$value)))
  expect_identical(controls[[2]]$xst, 1e-3 * pmax(abs(par), 1e-3))
  expect_identical(controls[[2]]$xt, controls[[2]]$xst * 5e-4)
  expect_identical(out$value@optinfo$optimizer, "Nelder_Mead")
  expect_lt(maihda_glmer_rise(out$value, kept$value), 1e-4)
  # An estimate at 0 takes a first step of 1e-6.
  controls <- list()
  quiet(maihda_glmer_finish(at_zero, small_call, env))
  expect_length(controls, 2L)
  expect_identical(controls[[2]]$xst[2], 1e-6)
  # Another optimiser of the user's own that fails is not retried, and the refit stands.
  controls <- list()
  call <- small_call
  call$control <- quote(user_control)
  env$user_control <- lme4::glmerControl(optimizer = c("bobyqa", "nloptwrap"))
  expect_null(quiet(maihda_glmer_finish(kept, call, env)))
  expect_length(controls, 1L)
})

test_that("a pass whose optimiser stops without converging counts as failed", {
  skip_on_cran()
  env <- small_env()
  kept <- bobyqa_fit(env)
  calls <- 0L
  stopped <- 1L
  real <- maihda_eval_holding_conditions
  # As an evaluation limit of the user's own stops Nelder-Mead: convergence code 4.
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    calls <<- calls + 1L
    out <- real(call, env)
    if (calls <= stopped) out$value@optinfo$conv$opt <- 4L
    out
  })
  # Nelder-Mead stopped: the pass is made again with small first steps, and kept.
  out <- quiet(maihda_glmer_finish(kept, small_call, env))
  expect_identical(calls, 2L)
  expect_false(is.null(get(".maihda_arg_finish_control", envir = env)$optCtrl$xst))
  expect_identical(as.integer(out$value@optinfo$conv$opt), 0L)
  # Stopped both times: the refit stands.
  calls <- 0L
  stopped <- 2L
  expect_null(quiet(maihda_glmer_finish(kept, small_call, env)))
  expect_identical(calls, 2L)
  # Another optimiser of the user's own that stops is not retried.
  calls <- 0L
  stopped <- 1L
  call <- small_call
  call$control <- quote(user_control)
  env$user_control <- lme4::glmerControl(optimizer = c("bobyqa", "nloptwrap"))
  expect_null(quiet(maihda_glmer_finish(kept, call, env)))
  expect_identical(calls, 1L)
})

test_that("the pass is kept unless it ends higher than the fit by more than the tolerance", {
  skip_on_cran()
  env <- small_env()
  kept <- bobyqa_fit(env)
  rise <- 0
  asked <- list()
  local_mocked_bindings(maihda_glmer_rise = function(to, from) {
    asked[[length(asked) + 1L]] <<- c(to = to@optinfo$optimizer, from = from@optinfo$optimizer)
    rise
  })
  rise <- 5e-5
  expect_false(is.null(quiet(maihda_glmer_finish(kept, small_call, env))))
  # The pass is measured against the fit, the pass first.
  expect_identical(asked[[1]], c(to = "Nelder_Mead", from = "bobyqa"))
  rise <- 2e-4
  expect_null(quiet(maihda_glmer_finish(kept, small_call, env)))
  expect_false(is.null(quiet(maihda_glmer_finish(kept, small_call, env, tol = 1e-3))))
  rise <- NA_real_
  expect_null(quiet(maihda_glmer_finish(kept, small_call, env)))
})

test_that("the rise is measured on one deviance function, the criteria where it cannot be", {
  skip_on_cran()
  env <- small_env()
  kept <- bobyqa_fit(env)
  out <- quiet(maihda_glmer_finish(kept, small_call, env))
  # Both ends of the pass at one maximum: no rise on the pass's own deviance function, whatever
  # the two fits' criteria, which carry the precision of lme4's PIRLS.
  expect_lt(abs(maihda_glmer_rise(out$value, kept$value)), 1e-4)
  expect_lt(abs(maihda_glmer_rise(kept$value, out$value)), 1e-4)
  # A fit placed below the maximum, its slope's own SD halved, rises from it.
  stay <- function(par, fn, lower, upper, control = list(), ...) {
    list(par = par, fval = fn(par), conv = 0, feval = 1)
  }
  low <- small_call
  low$start <- list(theta = unname(lme4::getME(kept$value, "theta")) * c(1, 1, 0.5),
                    fixef = lme4::fixef(kept$value))
  low$control <- lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE,
                                    boundary.tol = 0)
  placed <- quiet(eval(low, env))
  expect_gt(maihda_glmer_rise(placed, kept$value), 0.01)
  expect_lt(maihda_glmer_rise(kept$value, placed), -0.01)
  # The evaluation leaves the fits as they were.
  before <- lme4::getME(kept$value, "Lambdat")
  maihda_glmer_rise(kept$value, placed)
  expect_identical(lme4::getME(kept$value, "Lambdat"), before)
  # Without a deviance function, the criteria.
  local_mocked_bindings(maihda_lme4_detached_copy = function(model) stop("no copy"))
  expect_identical(maihda_glmer_rise(out$value, kept$value),
                   maihda_glmer_criterion(out$value) - maihda_glmer_criterion(kept$value))
})

# ---- in maihda_fit_glmer_bound() -----------------------------------------------------------

test_that("maihda_fit_glmer_bound() passes over a fit only when a refit replaced it", {
  env <- new.env()
  env$fake_fit <- function(crit) {
    message("first fit says")
    list(crit = crit)
  }
  pair <- function(crit, said) list(value = list(crit = crit), conditions = list(simpleMessage(said)))
  mirrored <- NULL
  searched <- NULL
  finished <- NULL
  given <- list()
  local_mocked_bindings(
    maihda_glmer_mirror_refit = function(kept, ...) mirrored,
    maihda_glmer_search_refit = function(kept, ...) searched,
    maihda_glmer_finish = function(kept, fit_call, fit_env, tol = 1e-4) {
      given[[length(given) + 1L]] <<- list(kept = kept$value, tol = tol)
      if (is.function(finished)) finished() else finished
    }
  )
  fit <- function(...) {
    seen <- character(0)
    value <- withCallingHandlers(maihda_fit_glmer_bound(quote(fake_fit(10)), env, ...),
      message = function(m) { seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage") })
    list(value = value, seen = seen)
  }
  # Neither refit kept: lme4's own fit, no pass.
  h <- fit()
  expect_identical(h$value$crit, 10)
  expect_length(given, 0L)
  # A refit across a bound kept: the pass is made from it, with the tolerance given.
  mirrored <- pair(8, "across\n")
  h <- fit(tol = 0.5)
  expect_identical(given[[1]], list(kept = list(crit = 8), tol = 0.5))
  # No pass returned: the refit stands, its own messages heard.
  expect_identical(h$value$crit, 8)
  expect_identical(h$seen, "across\n")
  # A pass returned replaces it, and speaks alone.
  finished <- pair(7.9, "finished\n")
  h <- fit()
  expect_identical(h$value$crit, 7.9)
  expect_identical(h$seen, "finished\n")
  # From the searches' refit when it is kept; a pass that fails is passed over.
  given <- list()
  mirrored <- NULL
  searched <- pair(6, "searched\n")
  finished <- function() stop("the pass fails")
  h <- fit()
  expect_identical(given[[1]]$kept, list(crit = 6))
  expect_identical(h$value$crit, 6)
  expect_identical(h$seen, "searched\n")
})

# ---- a live fit ----------------------------------------------------------------------------

# Eighteen strata of 80 rows, each row with a birth year between 1940 and 2000 and a Poisson count
# whose stratum slope deviates around 1970: the data of test-audit-2026-10-25-glmer-zero-variance-
# stop.R, where fit_maihda() refits glmer's fit from its searches to the maximum, 4935.556312 in
# -2 log-likelihood.
year_counts <- function(seed) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, 0.01)
    st <- rep(seq_len(18), each = 80)
    x <- sample(1940:2000, length(st), replace = TRUE)
    eta <- 0.5 + add[st] + (0.01 + v[st] / 2) * (x - 1970)
    data.frame(y = stats::rpois(length(st), exp(eta)), x = x, stratum = factor(st))
  })
}

test_that("lme4::refit() of a fit fit_maihda() refitted from its searches stays at the fit", {
  skip_on_cran()
  skip_slow_test()
  m <- quiet(fit_maihda(y ~ x + (1 + x | stratum), data = year_counts(1), family = "poisson"))
  expect_lt(abs(deviance_of(m$model) - 4935.556312), 1e-4)
  # Recorded as glmer's own fit records its second stage; refit() without a control runs it.
  expect_identical(m$model@optinfo$optimizer, "Nelder_Mead")
  expect_identical(m$model@call$control, quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  expect_lt(abs(deviance_of(quiet(lme4::refit(m$model))) - deviance_of(m$model)), 1e-3)
  expect_lt(abs(deviance_of(quiet(lme4::refit(m$model, control = lme4::glmerControl()))) -
                  deviance_of(m$model)), 1e-3)
  expect_true(m$diagnostics$converged)
  expect_output(print(m))
  expect_equal(deviance_of(quiet(stats::update(m$model))), deviance_of(m$model), tolerance = 1e-8)
})
