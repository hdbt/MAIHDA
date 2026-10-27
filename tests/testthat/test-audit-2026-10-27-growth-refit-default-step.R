# Audit pass 2026-10-27: glmer growth fits that a refit replaced recorded bobyqa, so lme4::refit()
# of them, with no control, ran bobyqa at minqa's defaults: on 29 simulated quadratic and linear
# growth fits so recorded it ended up to 1.8 above 2 of them, and on 87 responses simulated from
# them it stopped at minqa's limit of evaluations on 19, 13 of them within 0.001 of the maximum,
# where the refit of a bootstrap draw keeps the stop and counts it as not converged. Per the user,
# such a fit is now refitted once more from its estimates by the second stage of the call's own
# control (maihda_glmer_finish()), as other glmer fits a refit replaced are, and records that
# optimiser; the calls a bootstrap draw and a null-restricted refit replay take the user's own
# control whatever control the fit's call names.
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

# ---- the growth fitter, with stand-in fits ------------------------------------------------

# A stand-in fitting function: every fit has the criterion `crit`.
fake_fits <- function(crit) {
  env <- new.env()
  env$fake_fit <- function(start = NULL, control = NULL) list(crit = crit)
  env
}

# Stand-ins for the refits: the refit across a bound replaces a fit above `mirror_above` with one
# of criterion 6, the search refit a fit above `search_above` with one of criterion 5, and the
# finish records the fit it was given; each refit says what it is.
stand_in_refits <- function(mirror_above = Inf, search_above = Inf, finish = "finished") {
  finished <- NULL
  local_mocked_bindings(
    maihda_glmer_growth_starts = function(...) list(),
    maihda_glmer_criterion = function(model) model$crit,
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4, optimizer = NULL) {
      if (kept$value$crit > mirror_above) {
        list(value = list(crit = 6), conditions = list(simpleMessage("across\n")))
      }
    },
    maihda_glmer_refit_pass = function(refit, fit_call, fit_env, tol = 1e-4) NULL,
    maihda_glmer_search_refit = function(kept, fit_call, fit_env, tol = 1e-4, tolerance = NULL,
                                         gain = NULL) {
      if (kept$value$crit > search_above) {
        list(value = list(crit = 5), conditions = list(simpleMessage("searched\n")))
      }
    },
    maihda_glmer_finish = function(kept, fit_call, fit_env, tol = 1e-4) {
      finished <<- list(crit = kept$value$crit, call = fit_call, env = fit_env, tol = tol)
      switch(finish,
        finished = list(value = list(crit = kept$value$crit, finished = TRUE),
                        conditions = list(simpleMessage("finished\n"))),
        none = NULL,
        fails = stop("no pass"))
    },
    .env = parent.frame()
  )
  function() finished
}

test_that("a glmer growth fit that a refit replaced is refitted once more by the call's own second stage", {
  # Replaced by the refit across a bound and then by the refit from a search: the last pass is
  # made from the search's fit, with the fit's own call, environment and tolerance, and only the
  # pass speaks.
  got <- stand_in_refits(mirror_above = 7, search_above = 5.5)
  env <- fake_fits(10)
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid", tol = 0.5))
  expect_identical(got(), list(crit = 5, call = quote(fake_fit()), env = env, tol = 0.5))
  expect_identical(h$value, list(crit = 5, finished = TRUE))
  expect_identical(h$seen, "finished\n")
})

test_that("each refit of a glmer growth fit leads to the last pass, and a fit kept as made does not", {
  # Replaced by the refit across a bound alone.
  got <- stand_in_refits(mirror_above = 7)
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(10), "pid"))
  expect_identical(got()$crit, 6)
  expect_identical(h$value, list(crit = 6, finished = TRUE))
  # Replaced by the refit from a search alone.
  got <- stand_in_refits(search_above = 5.5)
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(10), "pid"))
  expect_identical(got()$crit, 5)
  expect_identical(h$value, list(crit = 5, finished = TRUE))
  # A bootstrap draw's refit given as the first fit and replaced: the same.
  env <- fake_fits(10)
  first <- maihda_eval_holding_conditions(quote(fake_fit()), env)
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid", first = first))
  expect_identical(h$value, list(crit = 5, finished = TRUE))
  # Nothing replaced it: no pass, and the fit as it was made.
  got <- stand_in_refits(mirror_above = 7, search_above = 5.5)
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(4), "pid"))
  expect_null(got())
  expect_identical(h$value, list(crit = 4))
  expect_identical(h$seen, character(0))
})

test_that("a last pass that is not kept or fails leaves the refit", {
  got <- stand_in_refits(search_above = 5.5, finish = "none")
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(10), "pid"))
  expect_identical(got()$crit, 5)
  expect_identical(h$value, list(crit = 5))
  expect_identical(h$seen, "searched\n")
  got <- stand_in_refits(mirror_above = 7, finish = "fails")
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), fake_fits(10), "pid"))
  expect_identical(h$value, list(crit = 6))
  expect_identical(h$seen, "across\n")
})

# ---- the last pass's rise ------------------------------------------------------------------

test_that("the rise of a last pass is measured at each fit's own negative-binomial dispersion", {
  skip_on_cran()
  # glmer.nb() searches the dispersion again in the pass: two fits of one model at dispersions
  # 2 (near the data's) and 20.
  d <- with_seed(3, {
    g <- factor(rep(1:12, each = 25))
    u <- stats::rnorm(12, 0, 0.5)
    data.frame(g = g, y = stats::rnbinom(300, mu = exp(1 + u[as.integer(g)]), size = 2))
  })
  fit_at <- function(theta, optimizer = "Nelder_Mead") {
    quiet(lme4::glmer(y ~ 1 + (1 | g), data = d, family = MASS::negative.binomial(theta),
                      control = lme4::glmerControl(optimizer = c("bobyqa", optimizer))))
  }
  near <- fit_at(2)
  far <- fit_at(20)
  at <- function(m) c(unname(lme4::getME(m, "theta")), unname(lme4::fixef(m)))
  own <- function(m) {
    devfun <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
    environment(devfun)$tolPwrss <- 1e-12
    devfun(at(m))
  }
  # Each at its own dispersion, far ends well above near; on far's deviance function alone near's
  # estimates, taken at far's dispersion, lie above far's own, which would call the pass lower.
  rise <- maihda_glmer_rise(far, near)
  expect_gt(rise, 1)
  expect_lt(abs(rise - (own(far) - own(near))), 1e-6)
  devfun <- lme4::getME(maihda_lme4_detached_copy(far), "devfun")
  expect_lt(devfun(at(far)) - devfun(at(near)), 0)
  # Two fits of one dispersion: on the first fit's deviance function, as for glmer() fits.
  other <- fit_at(2, optimizer = "bobyqa")
  devfun <- lme4::getME(maihda_lme4_detached_copy(near), "devfun")
  expect_identical(maihda_glmer_rise(near, other), devfun(at(near)) - devfun(at(other)))
})

# ---- the calls a draw and a null-restricted refit replay ------------------------------------

test_that("a growth fit's call is replayed with the user's own control, whatever control it names", {
  env <- new.env()
  call <- quote(lme4::glmer(n ~ t + (t | id), data = d, family = stats::poisson,
                            start = list(theta = 1), control = .maihda_arg_finish_control))
  # The last pass's control, and any other a refit names, by name or written out: none without a
  # control of the user's own, the user's with one.
  for (control in list(as.name(".maihda_arg_finish_control"), as.name(".maihda_arg_refit_control"),
                       as.name(".maihda_arg_search_control"),
                       quote(lme4::glmerControl(nAGQ0initStep = FALSE)))) {
    call$control <- control
    fresh <- maihda_growth_fresh_call(call, env)
    expect_null(fresh$start)
    expect_null(fresh$control)
    env$.maihda_arg_control <- lme4::glmerControl()
    expect_identical(maihda_growth_fresh_call(call, env)$control, as.name(".maihda_arg_control"))
    rm(".maihda_arg_control", envir = env)
  }
  # A call with the user's own control or none keeps it.
  call$control <- as.name(".maihda_arg_control")
  env$.maihda_arg_control <- lme4::glmerControl()
  expect_identical(maihda_growth_fresh_call(call, env)$control, as.name(".maihda_arg_control"))
  call$control <- NULL
  rm(".maihda_arg_control", envir = env)
  expect_identical(maihda_growth_fresh_call(call, env), local({ call$start <- NULL; call }))
})

# ---- a fit ----------------------------------------------------------------------------------

# The 0 / 1 design with 10 waves in days of the growth search's tests: under lme4 before 2.1-0
# glmer's quadratic growth fit stops 0.72 below the maximum (1274.712102), where the refit from a
# search takes it; lme4 2.1-0 fits it to the maximum.
growth_binary <- function(seed = 1, waves = 10, k = 365, s = c(0.2, 0.06, 0.012),
                          p = c(0.5, 0.15, 0.02)) {
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
    d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(-0.3 + 0.1 * t - 0.01 * t^2 + 0.6 * eta))
    d
  })
}
deviance_max_binary <- 1274.712102
# Where glmer's fit of it stops under lme4 before 2.1-0, rounded to seven digits (the stuck_binary
# of test-audit-2026-10-26-glmer-growth-time-units.R).
stuck_binary <- list(
  theta = c(0.5826776, -0.1972564, -0.1441493, 0.7696946, -0.2160046, 0.04805291, 0.3684006,
            -0.827052, 0.1727089, 0.2053519, -0.03951521, 0.008997533),
  fixef = c(-0.1853592, 0.2882030, -0.0790065))

# fit_maihda() on a quadratic growth design with glmer's first fit placed at `point` (theta and
# fixef on the internal time) by optimisers that return their start; the refits run as they do.
fit_from_point <- function(d, formula, family, point) {
  stay <- function(par, fn, lower, upper, control = list(), ...) {
    list(par = par, fval = fn(par), conv = 0, message = "the start", feval = 1)
  }
  real_growth <- maihda_fit_glmer_growth
  local_mocked_bindings(
    maihda_fit_glmer_growth = function(fit_call, fit_env, id, ...) {
      first_call <- fit_call
      first_call$start <- point
      first_call$control <- lme4::glmerControl(optimizer = list(stay, stay),
                                               nAGQ0initStep = FALSE, boundary.tol = 0)
      real_growth(fit_call, fit_env, id, first = maihda_eval_holding_conditions(first_call, fit_env),
                  ...)
    }
  )
  quiet(fit_maihda(formula, data = d, id = "id", time = "wave", family = family, time_degree = 2))
}

test_that("a null-restricted refit carries the user's control in place of the last pass's", {
  skip_on_cran()
  # A fit made at a given point by optimisers that return their start.
  stay <- function(par, fn, lower, upper, control = list(), ...) {
    list(par = par, fval = fn(par), conv = 0, message = "the start", feval = 1)
  }
  stay_control <- function() {
    lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE, boundary.tol = 0)
  }
  d <- growth_binary(waves = 4, k = 1)
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                        family = "binomial",
                        start = list(theta = c(0.5, 0, 0.1, 0.2, 0, 0.05),
                                     fixef = c("(Intercept)" = -0.3, wave = 0.1)),
                        control = stay_control()))$model
  env <- environment(stats::formula(m))
  local_mocked_bindings(maihda_refit_growth_call = function(...) stop("a stand-in refit"))
  carried_control <- function(control) {
    m@call$control <- control
    red <- maihda_refit_reduced_growth(m, stats::update(stats::formula(m), . ~ . - wave),
                                       list(id = "id"))
    expect_false(is.null(stats::getCall(red)$start))
    stats::getCall(red)$control
  }
  # The last pass's control, by name or written out as fit_maihda() writes it: the user's own.
  assign(".maihda_arg_finish_control", stay_control(), envir = env)
  expect_identical(carried_control(as.name(".maihda_arg_finish_control")),
                   as.name(".maihda_arg_control"))
  expect_identical(carried_control(quote(lme4::glmerControl(nAGQ0initStep = FALSE))),
                   as.name(".maihda_arg_control"))
  # The control of the refit across a bound is carried as it is.
  assign(".maihda_arg_refit_control", stay_control(), envir = env)
  expect_identical(carried_control(as.name(".maihda_arg_refit_control")),
                   as.name(".maihda_arg_refit_control"))
})

test_that("a glmer growth fit refitted from a search records Nelder-Mead and lme4::refit() stays at it", {
  skip_on_cran()
  skip_slow_test()
  d <- growth_binary()
  searched <- NULL
  real_finish <- maihda_glmer_finish
  # glmer's own fit before lme4 2.1-0; from 2.1-0, which fits the design to its maximum, the first
  # fit placed where the earlier glmer stops.
  own_stop <- utils::packageVersion("lme4") < "2.1-0"
  m <- with_mocked_bindings(
    if (own_stop) {
      quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                       family = "binomial", time_degree = 2))$model
    } else {
      fit_from_point(d, y ~ wave + (1 | A:B:C:D), "binomial", stuck_binary)$model
    },
    maihda_glmer_finish = function(kept, ...) {
      searched <<- kept$value
      real_finish(kept, ...)
    })
  # Fits compared on one deviance function, its penalised least squares run to 1e-10: on this
  # design lme4's own -2 log-likelihoods of two fits at one point differ by up to about 0.02.
  devfun <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
  environment(devfun)$tolPwrss <- 1e-10
  at <- function(fit) devfun(c(unname(lme4::getME(fit, "theta")), unname(lme4::fixef(fit))))
  # The fit the last pass is given, the refit from the search, is at the maximum and records
  # bobyqa; from the search that glmer's own stop leads to, lme4::refit() runs bobyqa at minqa's
  # defaults and ends 0.27 above it (from the placed stop's search, it stays).
  expect_identical(searched@optinfo$optimizer, "bobyqa")
  expect_identical(stats::getCall(searched)$control, as.name(".maihda_arg_search_control"))
  expect_lt(abs(at(searched) - deviance_max_binary), 1e-3)
  if (own_stop) {
    expect_gt(at(quiet(lme4::refit(searched))) - at(searched), 0.1)
  }
  # After it: the same maximum, recorded as glmer's own fits are, and lme4::refit() stays there
  # (within 4e-4).
  expect_identical(m@optinfo$optimizer, "Nelder_Mead")
  expect_identical(stats::getCall(m)$control, quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  expect_lt(abs(at(m) - at(searched)), 1e-4)
  expect_lt(at(quiet(lme4::refit(m))) - at(m), 0.01)
  # A bootstrap draw's call replays the fit's call without the last pass's control.
  replay <- maihda_growth_draw_call(m, stats::model.response(m@frame))
  expect_null(replay$call$control)
  expect_null(replay$call$start)
})
