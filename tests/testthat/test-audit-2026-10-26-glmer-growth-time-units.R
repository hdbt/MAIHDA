# Audit pass 2026-10-26: glmer growth fits that stop short of the maximum, quadratic ones above
# all, some recorded as converged.
#
# Growth time is fitted in units that put its span between 1 and 10, so a unit such as days does
# not by itself make glmer stop short. Quadratic growth does: a quadratic glmer growth fit gets no
# restarts, its maximum often has the quadratic term's variance at or near zero, and glmer's
# second stage stopped 47 of 240 simulated fits more than 0.01 below the maximum after the
# package's refits, by up to 21 in -2 log-likelihood; the 11 recorded as converged were singular,
# which lme4 does not check. Per the user, every glmer growth fit, and each bootstrap draw and
# null-restricted refit of one, is now searched as other glmer fits with a random slope are, the
# searches' deviance function run to a relative tolerance of 1e-10 -- at lme4's 1e-7 it starts
# its penalised least squares from the same point at every evaluation and stopped 0.022 above a
# Bernoulli growth fit's own -2 log-likelihood, where the search refuses the fit -- and
# refitted from the lower end when that lies at least 0.001 below: refits gaining less had
# taken half the time the searches added.
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

# Sixteen factorial strata (A..D) of 6 persons, `waves` occasions at wave * k, stratum and person
# intercepts, slopes and quadratic terms with SDs `s` and `p` (times 0.6 on the linear
# predictor), and a Poisson count n or a 0 / 1 outcome y.
growth_data <- function(seed, waves, k, outcome, s = c(0.2, 0.06, 0.012), p = c(0.5, 0.15, 0.02)) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
    add0 <- drop(as.matrix(g - 0.5) %*% c(0.3, 0.2, -0.25, 0.15))
    add1 <- drop(as.matrix(g - 0.5) %*% c(0.04, -0.03, 0.05, -0.02))
    u <- sapply(1:3, function(j) stats::rnorm(16, 0, s[j]))
    pe <- sapply(1:3, function(j) stats::rnorm(96, 0, p[j]))
    id <- rep(seq_len(96), each = waves)
    st <- rep(seq_len(16), each = 6)[id]
    t <- rep(seq_len(waves) - 1, times = 96)
    eta <- add0[st] + u[st, 1] + (add1[st] + u[st, 2]) * t + u[st, 3] * t^2 +
      pe[id, 1] + pe[id, 2] * t + pe[id, 3] * t^2
    d <- data.frame(id = id, wave = t * k, A = factor(g$A[st]), B = factor(g$B[st]),
                    C = factor(g$C[st]), D = factor(g$D[st]))
    if (outcome == "count") {
      d$n <- stats::rpois(nrow(d), exp(0.5 + 0.1 * t - 0.01 * t^2 + 0.6 * eta))
    } else {
      d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(-0.3 + 0.1 * t - 0.01 * t^2 + 0.6 * eta))
    }
    d
  })
}

deviance_of <- function(m) -2 * as.numeric(stats::logLik(m))

# Optimisers that return their start, so a fit is made at a given point on every platform.
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0, message = "the start", feval = 1)
}
stay_control <- function() {
  lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE, boundary.tol = 0)
}

# fit_maihda() on a quadratic growth design with its first fit placed at `point` (theta and fixef
# on the internal time). Returns the fit and the fit maihda_glmer_search_refit() was given.
fit_from_point <- function(d, formula, family, point) {
  given <- NULL
  real_growth <- maihda_fit_glmer_growth
  real_search <- maihda_glmer_search_refit
  local_mocked_bindings(
    maihda_fit_glmer_growth = function(fit_call, fit_env, id, ...) {
      first_call <- fit_call
      first_call$start <- point
      first_call$control <- stay_control()
      real_growth(fit_call, fit_env, id, first = maihda_eval_holding_conditions(first_call, fit_env),
                  ...)
    },
    maihda_glmer_search_refit = function(kept, ...) {
      given <<- kept$value
      real_search(kept, ...)
    }
  )
  fit <- quiet(fit_maihda(formula, data = d, id = "id", time = "wave", family = family,
                          time_degree = 2))
  list(fit = fit, given = given)
}

# Seed 6 of the Poisson design in days (6 waves at 0, 365, ..., 1825): glmer's fit stops with the
# person block's intercept SD near zero and the stratum block's last entry below 1e-4, singular,
# 9.31 below the maximum of the Laplace likelihood, which an unbounded search of glmer's deviance
# function (PIRLS to 1e-12) from many starts, and glmmTMB, put at 2042.071867. Its estimates,
# rounded to seven digits; the fixed effects named, as those of a refit from a search are, so
# that a reduced model's refit takes the ones it has.
stuck_days <- list(
  theta = c(4.112030e-03, 4.259621e-01, -3.482405e-01, 3.308587e-01, -5.933422e-03, 1.422582e-03,
            5.746344e-02, 1.900423e-02, 8.364204e-02, 1.681879e-01, -9.445602e-02, 7.092443e-05),
  fixef = c("(Intercept)" = 0.4491436, ".maihda_ctime" = 0.6360170,
            "I(.maihda_ctime^2)" = -0.3096026))
deviance_max_days <- 2042.071867

# Seed 1 of the 0 / 1 design with 10 waves in days: glmer's fit stops 0.72 below the maximum
# (1274.712102), recorded as not converged; its own -2 log-likelihood there, at the end of its
# optimisation, is 1275.432582. Its estimates, rounded to seven digits.
stuck_binary <- list(
  theta = c(0.5826776, -0.1972564, -0.1441493, 0.7696946, -0.2160046, 0.04805291, 0.3684006,
            -0.827052, 0.1727089, 0.2053519, -0.03951521, 0.008997533),
  fixef = c(-0.1853592, 0.2882030, -0.0790065))

# ---- the refits, with stand-in fits -------------------------------------------------------

# A stand-in fitting function: every fit has the criterion `crit`.
fake_fits <- function(crit) {
  env <- new.env()
  env$fake_fit <- function(start = NULL, control = NULL) list(crit = crit)
  env
}

test_that("a glmer growth fit is refitted from searches after the refits across a bound", {
  seen <- NULL
  local_mocked_bindings(
    maihda_glmer_growth_starts = function(...) list(),
    maihda_glmer_criterion = function(model) model$crit,
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4, optimizer = NULL) {
      if (kept$value$crit > 7) list(value = list(crit = 6), conditions = list(simpleMessage("across\n"))) else NULL
    },
    maihda_glmer_refit_pass = function(refit, fit_call, fit_env, tol = 1e-4) NULL,
    maihda_glmer_search_refit = function(kept, fit_call, fit_env, tol = 1e-4, tolerance = NULL,
                                         gain = NULL) {
      seen <<- list(crit = kept$value$crit, call = fit_call, tol = tol, tolerance = tolerance,
                    gain = gain)
      if (kept$value$crit > 5) list(value = list(crit = 5), conditions = list(simpleMessage("searched\n"))) else NULL
    }
  )
  # From the fit the refit across the bound kept, with the fit's own call, the tolerance given,
  # the deviance function to 1e-10 and a refit only for a gain of 0.001 or more; only the refit
  # from the search speaks.
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(10), "pid", tol = 0.5))
  expect_identical(seen, list(crit = 6, call = quote(fake_fit()), tol = 0.5, tolerance = 1e-10,
                              gain = 1e-3))
  expect_identical(h$value$crit, 5)
  expect_identical(h$seen, "searched\n")
  # A bootstrap draw's refit given as the first fit: the same refit.
  seen <- NULL
  env <- fake_fits(10)
  first <- maihda_eval_holding_conditions(quote(fake_fit()), env)
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid", first = first))
  expect_identical(seen$crit, 6)
  expect_identical(h$value$crit, 5)
  # Nothing lower: the kept fit.
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(4), "pid"))
  expect_identical(seen$crit, 4)
  expect_identical(h$value$crit, 4)
  # A search refit that fails is passed over.
  local_mocked_bindings(maihda_glmer_search_refit = function(...) stop("no search"))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(10), "pid"))
  expect_identical(h$value$crit, 6)
  expect_identical(h$seen, "across\n")
})

test_that("the search refit hands its tolerance to both searches", {
  got <- list()
  local_mocked_bindings(
    maihda_lme4_has_slope = function(model) TRUE,
    maihda_glmer_search = function(model, from, optimiser = NULL, tolerance = NULL) {
      got[[from]] <<- if (is.null(tolerance)) "lme4's" else tolerance
      NULL
    }
  )
  kept <- list(value = "a fit", conditions = list())
  expect_null(maihda_glmer_search_refit(kept, quote(f()), new.env(), tolerance = 1e-10))
  expect_identical(got, list(standardised = 1e-10, estimates = 1e-10))
  got <- list()
  expect_null(maihda_glmer_search_refit(kept, quote(f()), new.env()))
  expect_identical(got, list(standardised = "lme4's", estimates = "lme4's"))
})

test_that("the search refit is made only for the gain asked for", {
  local_mocked_bindings(
    maihda_lme4_has_slope = function(model) TRUE,
    maihda_glmer_search = function(model, from, optimiser = NULL, tolerance = NULL) {
      list(start = list(theta = 1, fixef = c(a = 1)), deviance = 9.9995, gain = 5e-4, lower = 0)
    }
  )
  kept <- list(value = "a fit", conditions = list())
  env <- new.env()
  env$f <- function(...) stop("the refit is made")
  # Lower by 5e-4: more than tol, so refitted where no gain is set, as for other glmer fits; not
  # for a gain of 0.001.
  expect_error(maihda_glmer_search_refit(kept, quote(f()), env), "the refit is made")
  expect_null(maihda_glmer_search_refit(kept, quote(f()), env, gain = 1e-3))
  expect_error(maihda_glmer_search_refit(kept, quote(f()), env, gain = 5e-4), "the refit is made")
})

# ---- the calls a draw and a null-restricted refit replay ------------------------------------

test_that("a growth fit's call is replayed without the control of a refit from a search", {
  call <- quote(lme4::glmer(n ~ t + (t | id), data = d, family = stats::poisson,
                            start = list(theta = 1), control = .maihda_arg_search_control))
  env <- new.env()
  fresh <- maihda_growth_fresh_call(call, env)
  expect_null(fresh$start)
  expect_null(fresh$control)
  env$.maihda_arg_control <- lme4::glmerControl()
  expect_identical(maihda_growth_fresh_call(call, env)$control, as.name(".maihda_arg_control"))
  # The control of a refit across a bound as before; a control of the user's own is kept.
  call$control <- as.name(".maihda_arg_refit_control")
  expect_identical(maihda_growth_fresh_call(call, env)$control, as.name(".maihda_arg_control"))
  call$control <- quote(my_control)
  expect_identical(maihda_growth_fresh_call(call, env)$control, quote(my_control))
})

test_that("a glmer growth fit's null-restricted refit leaves the fit's own controls alone", {
  skip_on_cran()
  d <- growth_data(6, 6, 365, "count")
  m <- quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                        family = "poisson", time_degree = 2, start = stuck_days,
                        control = stay_control()))$model
  env <- environment(stats::formula(m))
  assign(".maihda_arg_search_control", "the fit's", envir = env)
  used <- NULL
  local_mocked_bindings(
    maihda_refit_growth_call = function(model, fit_call, fit_env, spec, first = NULL) {
      used <<- fit_env
      assign(".maihda_arg_search_control", "the null fit's", envir = fit_env)
      stop("a stand-in refit")
    }
  )
  red_form <- stats::update(stats::formula(m), . ~ . - I(.maihda_ctime^2))
  maihda_refit_reduced_growth(m, red_form, list(id = "id"))
  # Refitted in a child of the fit's environment, which keeps the control the fit's call names.
  expect_false(identical(used, env))
  expect_identical(parent.env(used), env)
  expect_identical(get(".maihda_arg_search_control", envir = env), "the fit's")
})

test_that("a null-restricted refit from a search's start takes the user's control, not the search's", {
  skip_on_cran()
  d <- growth_data(6, 6, 365, "count")
  m <- quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                        family = "poisson", time_degree = 2, start = stuck_days,
                        control = stay_control()))$model
  # As if kept from a search: its call names the search's control beside its start.
  env <- environment(stats::formula(m))
  assign(".maihda_arg_search_control", stay_control(), envir = env)
  m@call$control <- as.name(".maihda_arg_search_control")
  local_mocked_bindings(maihda_refit_growth_call = function(...) stop("a stand-in refit"))
  red <- maihda_refit_reduced_growth(m, stats::update(stats::formula(m), . ~ . - I(.maihda_ctime^2)),
                                     list(id = "id"))
  expect_identical(stats::getCall(red)$control, as.name(".maihda_arg_control"))
  expect_false(is.null(stats::getCall(red)$start))
})

# ---- the searches' deviance function ------------------------------------------------------

test_that("the searches' deviance function is run to the tolerance asked for", {
  skip_on_cran()
  d <- growth_data(1, 10, 365, "binary")
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                        family = "binomial", time_degree = 2, start = stuck_binary,
                        control = stay_control()))$model
  at <- function(tolerance = NULL) {
    devfun <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
    if (!is.null(tolerance)) environment(devfun)$tolPwrss <- tolerance
    devfun(c(unname(lme4::getME(m, "theta")), unname(lme4::fixef(m))))
  }
  # The premise: from the point it starts at, the deviance function at lme4's tolerance stops
  # 0.022 above the -2 log-likelihood at the fit's estimates, which 1e-10 reaches (1e-12 agrees).
  # lme4 2.1-0 and later take the Laplace deviance at the converged step, and stop within 1e-5.
  precise <- at(1e-10)
  if (utils::packageVersion("lme4") < "2.1-0") {
    expect_gt(at() - precise, 0.015)
  } else {
    expect_lt(abs(at() - precise), 1e-4)
  }
  expect_lt(abs(at(1e-12) - precise), 1e-6)
  # A fit whose -2 log-likelihood is that value, as glmer's own fit's is at the end of its
  # optimisation, is re-expressed at either tolerance, its deviance at the fit the function's own
  # at the fit's estimates (test-audit-2026-10-27-glmer-search-refusal.R): at lme4's tolerance
  # 0.022 above, at 1e-10 that value.
  local_mocked_bindings(maihda_glmer_criterion = function(model) precise)
  expect_lt(abs(maihda_glmer_reexpressed(m)$at - at()), 1e-9)
  re <- maihda_glmer_reexpressed(m, tolerance = 1e-10)
  expect_lt(abs(re$at - precise), 1e-6)
  stand_in <- function(par, fn, lower, upper) list(par = par, fval = fn(par))
  expect_lt(abs(maihda_glmer_search(m, "estimates", optimiser = stand_in)$gain), 1e-6)
  expect_lt(abs(maihda_glmer_search(m, "estimates", optimiser = stand_in, tolerance = 1e-10)$gain),
            1e-6)
})

# ---- fits ---------------------------------------------------------------------------------

test_that("a quadratic glmer growth fit stopped short and called converged reaches the maximum", {
  skip_on_cran()
  skip_slow_test()
  d <- growth_data(6, 6, 365, "count")
  form <- n ~ wave + (1 | A:B:C:D)
  # Without the refit from the searches: glmer's stop, singular, lme4's only message, and
  # recorded as converged, 9.31 below the maximum. Days are fitted on time divided by 1000.
  before <- with_mocked_bindings(fit_from_point(d, form, "poisson", stuck_days)$fit,
                                 maihda_glmer_search_refit = function(...) NULL)
  expect_equal(range(before$model@frame$.maihda_ctime), c(0, 1.825))
  expect_true(lme4::isSingular(before$model))
  expect_identical(as.character(before$model@optinfo$conv$lme4$messages),
                   "boundary (singular) fit: see help('isSingular')")
  expect_true(before$diagnostics$converged)
  expect_gt(deviance_of(before$model) - deviance_max_days, 9.3)
  # With it: from that fit, the maximum (lme4's -2 log-likelihood about 3e-4 above the value its
  # deviance function reaches at 1e-12), made by glmer's second stage from the search.
  out <- fit_from_point(d, form, "poisson", stuck_days)
  expect_gt(deviance_of(out$given) - deviance_max_days, 9.3)
  expect_lt(deviance_of(out$fit$model) - deviance_max_days, 1e-3)
  expect_identical(stats::getCall(out$fit$model)$control, as.name(".maihda_arg_search_control"))
  expect_true(out$fit$diagnostics$converged)
})

test_that("a linear glmer growth fit at its maximum is left as it is", {
  skip_on_cran()
  skip_slow_test()
  # The 0 / 1 design in waves without quadratic terms: the searches gain about 2e-8.
  d <- growth_data(1, 6, 1, "binary", s = c(0.2, 0.06, 0), p = c(0.5, 0.15, 0))
  fit <- function() {
    quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                     family = "binomial"))$model
  }
  searched <- 0L
  real_search <- maihda_glmer_search_refit
  with_search <- with_mocked_bindings(fit(), maihda_glmer_search_refit = function(...) {
    searched <<- searched + 1L
    real_search(...)
  })
  without <- with_mocked_bindings(fit(), maihda_glmer_search_refit = function(...) NULL)
  expect_identical(searched, 1L)
  expect_identical(lme4::getME(with_search, "theta"), lme4::getME(without, "theta"))
  expect_identical(lme4::fixef(with_search), lme4::fixef(without))
  expect_identical(deparse(stats::getCall(with_search)), deparse(stats::getCall(without)))
  expect_null(stats::getCall(with_search)$control)
})
