# Audit 2026-11-16, one finding: the maximum-likelihood refit of a REML lmer fit could end far
# above the maximum without an exit code or a warning.
#
# maihda_ic() with more than one model, compare_maihda() and the PCVs of estimation = "ML" read
# a REML lmer fit's maximum-likelihood fit from lme4::refitML(), which runs minqa's bobyqa from
# the REML estimates. bobyqa takes one first step for all parameters, min(0.95, 0.2 * max|theta|),
# and before it starts it moves every parameter closer to its bound than that step up to the
# step. With a random slope on a covariate far from zero or in large units the intercept's entry
# of theta is large and the slopes' are small, so the refit starts hundreds of deviance units
# above the REML estimates it was given and ends with "Normal exit": 220 of 280 fits with a slope
# on an age, a birth or calendar year or an income ended more than 1e-4 above the maximum, by up
# to 165, and AIC or BIC preferred the other model in 15 of 30 pairs with slopes on an age and a
# birth year. Of 224 growth fits 8 ended short, by up to 0.40, some with the adjusted model's
# stratum intercept variance at zero, which is a PCV of exactly 1 where the maximum-likelihood
# fits gave 0.9997 to -0.18. Of 786 fits without such a slope none did.
#
# FIX: maihda_lmer_refit_ml() makes the refit. A fit with a random-effect term other than an
# intercept is also fitted by maximum likelihood from its own call, through the fitter
# fit_maihda() used, and that fit replaces lme4's refit when it is a fit of the same model (the
# same response, offset, weights and design matrices) and its deviance is lower by more than
# 1e-4; of what lme4 raises for it only an optimiser's exit code is passed on.
# maihda_pcv_refit_ml() and maihda_longitudinal_refit_ml() make their refits through it too.
#
# The optima below are those of the maximum-likelihood deviance written out as dense algebra,
# with nothing of lme4, minimised from 25 starts with the random-effect terms centred and scaled
# (every start within 1e-6 of the value given); for the growth pair, of lme4's deviance function
# searched the same way from 13 starts and evaluated by that algebra at its minimum.
#
# Random numbers are drawn only inside mr16_keep_seed(), which puts the suite's stream back.

mr16_keep_seed <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  force(expr)
}

mr16_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression with the text of every warning and message it raised.
mr16_heard <- function(expr) {
  warnings <- character(0)
  messages <- character(0)
  value <- withCallingHandlers(expr,
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    })
  list(value = value, warnings = warnings, messages = messages)
}

# Eighteen strata (3 x 3 x 2 dimensions a, b, c with additive effects, and a stratum deviation
# of SD 0.3) of `n` rows, residual SD 1, with a covariate x whose slope differs by stratum:
# "z" standard normal, "income" 18,000 to 80,000 in currency units, "age" 18 to 80 beside a
# birth year `yr` 1940 to 2000 with a slope of its own.
mr16_data <- function(seed, x = c("z", "income", "age"), n = 40) {
  x <- match.arg(x)
  mr16_keep_seed({
    set.seed(seed)
    lev <- expand.grid(a = 1:3, b = 1:3, c = 1:2)
    st <- rep(seq_len(18), each = n)
    rows <- length(st)
    main <- stats::rnorm(3, 0, 0.3)[lev$a] + stats::rnorm(3, 0, 0.3)[lev$b] +
      stats::rnorm(2, 0, 0.3)[lev$c]
    d <- data.frame(a = factor(lev$a[st]), b = factor(lev$b[st]), c = factor(lev$c[st]),
                    stratum = factor(st))
    y <- 1 + main[st] + stats::rnorm(18, 0, 0.3)[st]
    if (x == "z") {
      d$x <- stats::rnorm(rows)
      y <- y + (0.5 + stats::rnorm(18, 0, 0.15)[st]) * d$x
    } else if (x == "income") {
      d$x <- round(stats::runif(rows, 18000, 80000))
      y <- y + (3e-05 + stats::rnorm(18, 0, 1e-05)[st]) * (d$x - 49000)
    } else {
      d$x <- sample(18:80, rows, replace = TRUE)
      d$yr <- sample(1940:2000, rows, replace = TRUE)
      y <- y + (0.03 + stats::rnorm(18, 0, 0.01)[st]) * (d$x - 49) +
        (0.02 + stats::rnorm(18, 0, 0.008)[st]) * (d$yr - 1970)
    }
    d$y <- y + stats::rnorm(rows)
    d
  })
}

# Long data: 10 persons in each of 8 strata (dimensions a and b), seen at times 0 to 3, with
# person-level intercepts and slopes (SDs 1 and 0.2) and stratum-level ones (0.3 and 0.08).
mr16_long <- function(seed) {
  mr16_keep_seed({
    set.seed(seed)
    lev <- expand.grid(a = 1:4, b = 1:2)
    main <- stats::rnorm(4, 0, 0.3)[lev$a] + stats::rnorm(2, 0, 0.3)[lev$b]
    pst <- rep(seq_len(8), each = 10)
    rows <- expand.grid(k = 1:4, pid = seq_len(80))
    st <- pst[rows$pid]
    t <- rows$k - 1
    u0 <- stats::rnorm(8, 0, 0.3)
    u1 <- stats::rnorm(8, 0, 0.08)
    b0 <- stats::rnorm(80, 0, 1)
    b1 <- stats::rnorm(80, 0, 0.2)
    y <- 1 + main[st] + u0[st] + b0[rows$pid] + (0.2 + u1[st] + b1[rows$pid]) * t +
      stats::rnorm(nrow(rows), 0, 0.7)
    data.frame(y = y, time = t, pid = factor(rows$pid), a = factor(lev$a[st]),
               b = factor(lev$b[st]), stratum = factor(st))
  })
}

# The fits of a fixture, made once for the file by the first block that asks for them.
mr16_made <- new.env()
mr16_fits <- function(kind) {
  if (is.null(mr16_made[[kind]])) {
    mr16_made[[kind]] <- mr16_quiet(switch(kind,
      z = {
        d <- mr16_data(3, "z")
        list(data = d, null = fit_maihda(y ~ x + (1 + x | stratum), data = d),
             adj = fit_maihda(y ~ x + a + b + c + (1 + x | stratum), data = d))
      },
      plain = {
        d <- mr16_data(3, "z")
        list(data = d, null = fit_maihda(y ~ 1 + (1 | stratum), data = d),
             adj = fit_maihda(y ~ a + b + c + (1 | stratum), data = d))
      },
      income = {
        d <- mr16_data(1, "income")
        list(data = d, null = fit_maihda(y ~ x + (1 + x | stratum), data = d),
             adj = fit_maihda(y ~ x + a + b + c + (1 + x | stratum), data = d))
      },
      long = {
        d <- mr16_long(4)
        list(data = d,
             null = fit_maihda(y ~ time + (1 | stratum), data = d, id = "pid", time = "time"),
             adj = fit_maihda(y ~ time + a + b + (1 | stratum), data = d, id = "pid",
                              time = "time"))
      }))
  }
  mr16_made[[kind]]
}

# Maximum-likelihood deviances at the optimum, null model then adjusted model.
optimum_z <- c(2058.849737, 2036.774338)
optimum_income <- c(2192.369172, 2165.243399)
optimum_age <- c(3196.889305, 3174.801661)
optimum_long <- c(848.560036, 835.608676)

# The model of the lmer fit `fit` as a maximum-likelihood fit placed at `theta`: the fit's own
# call with an optimiser that returns its start, which lme4 records with exit code 0.
mr16_stay <- function(fn, par, lower, upper, control) {
  list(par = par, fval = fn(par), conv = 0, feval = 1, message = "stay")
}
mr16_placed <- function(fit, theta) {
  call <- stats::getCall(fit)
  call$REML <- FALSE
  call$start <- theta
  call$control <- lme4::lmerControl(optimizer = mr16_stay, boundary.tol = 0)
  mr16_quiet(eval(call, environment(stats::formula(fit))))
}
# The REML estimates of `fit` with every variance parameter at 0.95, where bobyqa's first step
# puts those that lie closer to zero when an entry of theta is 4.75 or more.
mr16_moved <- function(fit) {
  theta <- unname(lme4::getME(fit, "theta"))
  theta[lme4::getME(fit, "lower") == 0] <- 0.95
  theta
}
# A stand-in for lme4::refitML() that returns `refit` for every model.
mr16_refit_is <- function(refit) function(x, optimizer = "bobyqa", ...) refit

# ---- the refit --------------------------------------------------------------------------------

test_that("a refit that ends short of the maximum without a code gives way to the model's own maximum-likelihood fit", {
  skip_on_cran()
  fits <- mr16_fits("z")
  fm <- fits$null$model
  short <- mr16_placed(fm, mr16_moved(fm))
  # The premise: lme4 records the placed fit as a refit that ended normally, far above the
  # maximum.
  expect_false(lme4::isREML(short))
  expect_false(maihda_lme4_optimizer_failed(short))
  expect_gt(stats::deviance(short) - optimum_z[1], 10)
  calls <- character(0)
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    calls <<- c(calls, optimizer)
    short
  }, .package = "lme4")
  # The criteria are those of the maximum, read without a warning or a message.
  row <- mr16_heard(maihda_ic_one(fits$null, ml = TRUE))
  expect_identical(row$warnings, character(0))
  expect_identical(row$messages, character(0))
  expect_identical(row$value$estimator, "ML (refit from REML)")
  expect_lt(abs(-2 * row$value$logLik - optimum_z[1]), 1e-3)
  expect_lt(abs(row$value$AIC - (optimum_z[1] + 2 * 6)), 1e-3)
  expect_lt(abs(row$value$BIC - (optimum_z[1] + log(720) * 6)), 1e-3)
  # No code, so no second refit by lmer's default optimiser: one call of lme4::refitML().
  expect_identical(calls, "bobyqa")
  h <- mr16_heard(maihda_lmer_refit_ml(fm))
  expect_identical(h$warnings, character(0))
  expect_identical(h$messages, character(0))
  expect_s4_class(h$value, "lmerMod")
  expect_false(lme4::isREML(h$value))
  expect_lt(abs(stats::deviance(h$value) - optimum_z[1]), 1e-3)
  # It is the fit fit_maihda() makes of the model by maximum likelihood.
  direct <- mr16_quiet(fit_maihda(y ~ x + (1 + x | stratum), data = fits$data, REML = FALSE))
  expect_lt(abs(stats::deviance(h$value) - stats::deviance(direct$model)), 1e-6)
})

test_that("lme4's refit stays unless the model's own fit is lower by more than the tolerance of 1e-4", {
  skip_on_cran()
  fits <- mr16_fits("z")
  fm <- fits$null$model
  direct <- mr16_quiet(fit_maihda(y ~ x + (1 + x | stratum), data = fits$data,
                                  REML = FALSE))$model
  at <- unname(lme4::getME(direct, "theta"))
  deviance_at <- mr16_quiet(stats::update(fm, REML = FALSE, devFunOnly = TRUE))
  # A refit `by` above the model's own fit: its first parameter raised until the deviance is.
  above <- function(by) {
    rise <- function(e) deviance_at(at * c(1 + e, 1, 1)) - stats::deviance(direct) - by
    mr16_placed(fm, at * c(1 + stats::uniroot(rise, c(0, 1), tol = 1e-12)$root, 1, 1))
  }
  within <- above(3e-5)
  beyond <- above(3e-4)
  expect_lt(abs(stats::deviance(within) - stats::deviance(direct) - 3e-5), 1e-6)
  expect_lt(abs(stats::deviance(beyond) - stats::deviance(direct) - 3e-4), 1e-6)
  # Within the tolerance lme4's refit is returned as it is.
  testthat::local_mocked_bindings(refitML = mr16_refit_is(within), .package = "lme4")
  expect_identical(maihda_ic_one(fits$null, ml = TRUE)$AIC, as.numeric(stats::AIC(within)))
  kept <- maihda_lmer_refit_ml(fm)
  expect_identical(stats::deviance(kept), stats::deviance(within))
  expect_identical(lme4::getME(kept, "theta"), lme4::getME(within, "theta"))
  # Beyond it the model's own fit is.
  testthat::local_mocked_bindings(refitML = mr16_refit_is(beyond), .package = "lme4")
  expect_lt(abs(maihda_ic_one(fits$null, ml = TRUE)$AIC - stats::AIC(direct)), 1e-6)
  expect_lt(abs(stats::deviance(maihda_lmer_refit_ml(fm)) - stats::deviance(direct)), 1e-6)
  # The tolerance is the one ?maihda_ic gives.
  expect_identical(formals(maihda_lmer_refit_ml)$tol, 1e-4)
  expect_identical(formals(maihda_lmer_checked_refit)$tol, 1e-4)
})

test_that("a refit at the maximum comes back as lme4 made it", {
  skip_on_cran()
  fits <- mr16_fits("z")
  for (model in list(fits$null, fits$adj)) {
    fm <- model$model
    own <- lme4::refitML(fm)
    from_call <- mr16_quiet(maihda_lmer_ml_from_call(fm))
    # Where lme4's refit is within the tolerance of the model's own fit, as it is for a slope on
    # a standard-normal covariate on every build measured, nothing moves.
    skip_if(stats::deviance(own) - stats::deviance(from_call) > 5e-5,
            "lme4's refit of this fit is not at the maximum on this build")
    h <- mr16_heard(maihda_lmer_refit_ml(fm))
    expect_identical(h$warnings, character(0))
    expect_identical(h$messages, character(0))
    expect_identical(stats::deviance(h$value), stats::deviance(own))
    expect_identical(lme4::getME(h$value, "theta"), lme4::getME(own, "theta"))
    expect_identical(lme4::fixef(h$value), lme4::fixef(own))
    row <- maihda_ic_one(model, ml = TRUE)
    expect_identical(row$AIC, as.numeric(stats::AIC(own)))
    expect_identical(row$BIC, as.numeric(stats::BIC(own)))
  }
  expect_lt(abs(stats::deviance(lme4::refitML(fits$null$model)) - optimum_z[1]), 1e-3)
  expect_lt(abs(stats::deviance(lme4::refitML(fits$adj$model)) - optimum_z[2]), 1e-3)
})

test_that("only a fit with a random slope is fitted from its own call, and once", {
  skip_on_cran()
  real <- maihda_lmer_ml_from_call
  made <- 0L
  testthat::local_mocked_bindings(maihda_lmer_ml_from_call = function(fit, longitudinal = NULL) {
    made <<- made + 1L
    real(fit, longitudinal)
  }, .package = "MAIHDA")
  # Random intercepts only: lme4's refit alone, to the last digit.
  plain <- mr16_fits("plain")
  expect_false(maihda_lme4_has_slope(plain$adj$model))
  refit <- maihda_lmer_refit_ml(plain$adj$model)
  expect_identical(made, 0L)
  own <- lme4::refitML(plain$adj$model)
  expect_identical(stats::deviance(refit), stats::deviance(own))
  expect_identical(lme4::getME(refit, "theta"), lme4::getME(own, "theta"))
  expect_identical(maihda_ic_one(plain$adj, ml = TRUE)$AIC, as.numeric(stats::AIC(own)))
  expect_identical(made, 0L)
  # A random slope: one fit from the call for each refit.
  fits <- mr16_fits("z")
  expect_true(maihda_lme4_has_slope(fits$adj$model))
  mr16_quiet(maihda_lmer_refit_ml(fits$adj$model))
  expect_identical(made, 1L)
  mr16_quiet(maihda_ic(fits$null, fits$adj))
  expect_identical(made, 3L)
})

test_that("the model's own fit is made as fit_maihda() made the fit", {
  skip_on_cran()
  real_slopes <- maihda_fit_lmer_slopes
  real_growth <- maihda_refit_growth_call
  seen <- list()
  testthat::local_mocked_bindings(
    maihda_fit_lmer_slopes = function(fit_call, fit_env, ...) {
      seen[[length(seen) + 1L]] <<- list(by = "slopes", call = fit_call)
      real_slopes(fit_call, fit_env, ...)
    },
    maihda_refit_growth_call = function(model, fit_call, fit_env, spec, first = NULL) {
      seen[[length(seen) + 1L]] <<- list(by = "growth", call = fit_call, spec = spec,
                                         first = first)
      real_growth(model, fit_call, fit_env, spec, first = first)
    }, .package = "MAIHDA")
  fits <- mr16_fits("z")
  fm <- fits$null$model
  # A cross-sectional fit: the refits of maihda_fit_lmer_slopes(), by maximum likelihood, in
  # the environment of the fit's formula.
  own <- mr16_quiet(maihda_lmer_ml_from_call(fm))
  expect_length(seen, 1L)
  expect_identical(seen[[1]]$by, "slopes")
  expect_identical(seen[[1]]$call$REML, FALSE)
  expect_false(lme4::isREML(own))
  expect_identical(stats::nobs(own), stats::nobs(fm))
  # The numeric start of a refit fit_maihda() kept is not carried into it.
  restarted <- fm
  restarted@call$start <- unname(lme4::getME(fm, "theta"))
  mr16_quiet(maihda_lmer_ml_from_call(restarted))
  expect_length(seen, 2L)
  expect_identical(seen[[2]]$by, "slopes")
  expect_null(seen[[2]]$call$start)
  expect_identical(seen[[2]]$call$REML, FALSE)
  # A fit made from starting values of the user's own, which fit_maihda() does not refit, is
  # fitted from them alone.
  started <- mr16_quiet(fit_maihda(y ~ x + (1 + x | stratum), data = fits$data,
                                   start = c(1, 0, 1)))
  expect_true(is.name(stats::getCall(started$model)$start))
  own <- mr16_quiet(maihda_lmer_ml_from_call(started$model))
  expect_length(seen, 2L)
  expect_false(lme4::isREML(own))
  expect_true(is.name(stats::getCall(own)$start))
  # A growth model: through the growth fitter, from lme4's starting values (the start of a
  # restart fit_maihda() kept is not carried into it), with the tag of the model; the criteria
  # and the PCVs of estimation = "ML" pass the tag on.
  long <- mr16_fits("long")
  restarted <- long$null$model
  restarted@call$start <- unname(lme4::getME(restarted, "theta"))
  own <- mr16_quiet(maihda_lmer_ml_from_call(restarted, long$null$longitudinal_info))
  expect_length(seen, 3L)
  expect_identical(seen[[3]]$by, "growth")
  expect_identical(seen[[3]]$call$REML, FALSE)
  expect_null(seen[[3]]$call$start)
  expect_null(seen[[3]]$first)
  expect_identical(seen[[3]]$spec$id, "pid")
  expect_identical(seen[[3]]$spec$time_degree, 1L)
  expect_false(lme4::isREML(own))
  expect_identical(stats::nobs(own), stats::nobs(long$null$model))
  mr16_quiet(maihda_ic_one(long$adj, ml = TRUE))
  expect_length(seen, 4L)
  expect_identical(seen[[4]]$by, "growth")
  mr16_quiet(maihda_longitudinal_pcv(long$null, long$adj, estimation = "ML"))
  expect_identical(vapply(seen[5:6], `[[`, character(1), "by"), c("growth", "growth"))
  expect_length(seen, 6L)
  # A growth fit made from starting values of the user's own is fitted from them alone.
  started <- mr16_quiet(fit_maihda(y ~ time + (1 | stratum), data = long$data, id = "pid",
                                   time = "time", start = c(1, 0, 1, 1, 0, 1)))
  expect_true(is.name(stats::getCall(started$model)$start))
  own <- mr16_quiet(maihda_lmer_ml_from_call(started$model, started$longitudinal_info))
  expect_length(seen, 6L)
  expect_false(lme4::isREML(own))
  expect_true(is.name(stats::getCall(own)$start))
  # So is one whose tag gives no growth structure, which fit_maihda() fits by lmer() alone.
  tag <- long$null$longitudinal_info
  tag$time_degree <- NA
  expect_null(maihda_lme4_growth_spec(long$null$model, tag))
  own <- mr16_quiet(maihda_lmer_ml_from_call(long$null$model, tag))
  expect_length(seen, 6L)
  expect_false(lme4::isREML(own))
  expect_null(stats::getCall(own)$start)
})

test_that("of what lme4 says about the model's own fit only an optimiser's exit code is passed on", {
  skip_on_cran()
  fits <- mr16_fits("z")
  fm <- fits$null$model
  direct <- mr16_quiet(fit_maihda(y ~ x + (1 + x | stratum), data = fits$data,
                                  REML = FALSE))$model
  short <- mr16_placed(fm, mr16_moved(fm))
  code <- "convergence code 5 from nloptwrap: NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached."
  speaks <- function(fit) {
    function(fit_arg, longitudinal = NULL) {
      message("boundary (singular) fit: see help('isSingular')")
      warning("Model failed to converge with max|grad| = 2.63106 (tol = 0.002, component 1)")
      warning(code)
      fit
    }
  }
  # The model's own fit is returned: its exit code is passed on, lme4's checks are not.
  testthat::local_mocked_bindings(refitML = mr16_refit_is(short), .package = "lme4")
  testthat::local_mocked_bindings(maihda_lmer_ml_from_call = speaks(direct), .package = "MAIHDA")
  h <- mr16_heard(maihda_lmer_refit_ml(fm))
  expect_identical(stats::deviance(h$value), stats::deviance(direct))
  expect_identical(h$warnings, code)
  expect_identical(h$messages, character(0))
  row <- mr16_heard(maihda_ic_one(fits$null, ml = TRUE))
  expect_identical(row$warnings, code)
  expect_identical(row$messages, character(0))
  # lme4's refit is returned: nothing the model's own fit raised is passed on.
  testthat::local_mocked_bindings(refitML = mr16_refit_is(direct), .package = "lme4")
  testthat::local_mocked_bindings(maihda_lmer_ml_from_call = speaks(short), .package = "MAIHDA")
  h <- mr16_heard(maihda_lmer_refit_ml(fm))
  expect_identical(stats::deviance(h$value), stats::deviance(direct))
  expect_identical(h$warnings, character(0))
  expect_identical(h$messages, character(0))
  # What lme4's refit raised is passed on with it, as for a model without a slope.
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    warning("a warning of the refit")
    message("a message of the refit")
    direct
  }, .package = "lme4")
  h <- mr16_heard(maihda_lmer_refit_ml(fm))
  expect_identical(h$warnings, "a warning of the refit")
  expect_identical(h$messages, "a message of the refit\n")
  # And dropped when the model's own fit is returned in its place.
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    warning("a warning of the refit")
    short
  }, .package = "lme4")
  testthat::local_mocked_bindings(maihda_lmer_ml_from_call = function(fit, longitudinal = NULL) {
    direct
  }, .package = "MAIHDA")
  h <- mr16_heard(maihda_lmer_refit_ml(fm))
  expect_identical(stats::deviance(h$value), stats::deviance(direct))
  expect_identical(h$warnings, character(0))
})

test_that("a fit from the call that fails, or is not this model's maximum-likelihood fit, leaves lme4's refit", {
  skip_on_cran()
  fits <- mr16_fits("z")
  fm <- fits$null$model
  direct <- mr16_quiet(fit_maihda(y ~ x + (1 + x | stratum), data = fits$data,
                                  REML = FALSE))$model
  short <- mr16_placed(fm, mr16_moved(fm))
  testthat::local_mocked_bindings(refitML = mr16_refit_is(short), .package = "lme4")
  returns <- function(fit) {
    testthat::local_mocked_bindings(maihda_lmer_ml_from_call = function(fit_arg,
                                                                        longitudinal = NULL) {
      if (is.null(fit)) stop("the call cannot be evaluated") else fit
    }, .package = "MAIHDA")
    mr16_heard(maihda_lmer_refit_ml(fm))
  }
  # The premise: this model's own maximum-likelihood fit replaces the short refit.
  expect_identical(stats::deviance(returns(direct)$value), stats::deviance(direct))
  # An error: the refit, and nothing is raised.
  h <- returns(NULL)
  expect_identical(stats::deviance(h$value), stats::deviance(short))
  expect_identical(h$warnings, character(0))
  # A REML fit, whose criterion is not a deviance.
  expect_lt(lme4::REMLcrit(fm), stats::deviance(short))
  expect_identical(stats::deviance(returns(fm)$value), stats::deviance(short))
  # Another response, other fixed effects, other random-effect terms: each a lower deviance.
  d <- fits$data
  d$y <- d$y / 10
  other_response <- mr16_quiet(lme4::lmer(y ~ x + (1 + x | stratum), data = d, REML = FALSE))
  other_fixed <- mr16_quiet(lme4::lmer(y ~ x + a + (1 + x | stratum), data = fits$data,
                                       REML = FALSE))
  other_random <- mr16_quiet(lme4::lmer(y ~ x + (1 + x | stratum) + (1 | a), data = fits$data,
                                        REML = FALSE))
  for (other in list(other_response, other_fixed)) {
    expect_lt(stats::deviance(other), stats::deviance(short))
    expect_identical(stats::deviance(returns(other)$value), stats::deviance(short))
  }
  expect_identical(stats::deviance(returns(other_random)$value), stats::deviance(short))
  # The same formula on another design, with weights, with an offset, with the strata in
  # another order: each has the deviance of this model's own fit and the same response, names
  # and terms, and is another model's fit all the same.
  d <- fits$data
  d$x <- 2 * d$x
  other_design <- mr16_quiet(lme4::lmer(y ~ x + (1 + x | stratum), data = d, REML = FALSE))
  d <- fits$data
  d$w <- 4
  other_weights <- mr16_quiet(lme4::lmer(y ~ x + (1 + x | stratum), data = d, weights = w,
                                         REML = FALSE))
  d$off <- 1
  other_offset <- mr16_quiet(lme4::lmer(y ~ x + (1 + x | stratum), data = d, offset = off,
                                        REML = FALSE))
  d$stratum <- factor(d$stratum, levels = rev(levels(d$stratum)))
  other_order <- mr16_quiet(lme4::lmer(y ~ x + (1 + x | stratum), data = d, REML = FALSE))
  for (other in list(other_design, other_weights, other_offset, other_order)) {
    expect_lt(abs(stats::deviance(other) - stats::deviance(direct)), 1e-4)
    expect_identical(names(lme4::fixef(other)), names(lme4::fixef(direct)))
    expect_identical(lme4::getME(other, "cnms"), lme4::getME(direct, "cnms"))
    expect_false(maihda_lmer_same_model(other, short))
    expect_identical(stats::deviance(returns(other)$value), stats::deviance(short))
  }
  # A slope covariate that is in no fixed effect, doubled: the same fixed-effect design, the
  # same random-effect pattern, other random-effect values.
  slope_only <- mr16_quiet(lme4::lmer(y ~ 1 + (1 + x | stratum), data = fits$data,
                                      REML = FALSE))
  d <- fits$data
  d$x <- 2 * d$x
  other_slope <- mr16_quiet(lme4::lmer(y ~ 1 + (1 + x | stratum), data = d, REML = FALSE))
  expect_true(maihda_lmer_same_model(slope_only, slope_only))
  expect_false(maihda_lmer_same_model(other_slope, slope_only))
  expect_true(maihda_lmer_same_model(direct, short))
  expect_true(maihda_lmer_same_model(direct, lme4::refitML(fm)))
  expect_false(maihda_lmer_same_model(fm, short))
  expect_false(maihda_lmer_same_model(other_response, short))
  expect_false(maihda_lmer_same_model(other_fixed, short))
  expect_false(maihda_lmer_same_model(other_random, short))
  expect_false(maihda_lmer_same_model(NULL, short))
})

test_that("the model's own fit is made on the fit's rows, with its weights and its offset", {
  skip_on_cran()
  d <- mr16_data(3, "z")
  d$w <- rep(c(1, 2, 0.5, 1.5), length.out = nrow(d))
  d$off <- 0.1 * as.integer(d$a)
  m <- mr16_quiet(fit_maihda(y ~ x + (1 + x | stratum), data = d, weights = w,
                             subset = x > -1.2, offset = off))
  fm <- m$model
  # The premise: fewer rows than the data, unequal weights, an offset that is not zero.
  expect_lt(stats::nobs(fm), nrow(d))
  expect_gt(stats::sd(stats::weights(fm)), 0)
  expect_gt(stats::sd(lme4::getME(fm, "offset")), 0)
  own <- mr16_quiet(maihda_lmer_ml_from_call(fm))
  expect_false(lme4::isREML(own))
  expect_identical(stats::nobs(own), stats::nobs(fm))
  expect_identical(unname(stats::weights(own)), unname(stats::weights(fm)))
  expect_identical(unname(lme4::getME(own, "offset")), unname(lme4::getME(fm, "offset")))
  expect_identical(unname(lme4::getME(own, "y")), unname(lme4::getME(fm, "y")))
  expect_true(maihda_lmer_same_model(own, lme4::refitML(fm)))
  # A refit that ends short gives way to it, and the criteria are those of the weighted model
  # with its offset, fitted by maximum likelihood on those rows.
  direct <- mr16_quiet(fit_maihda(y ~ x + (1 + x | stratum), data = d, weights = w,
                                  subset = x > -1.2, offset = off, REML = FALSE))$model
  short <- mr16_placed(fm, mr16_moved(fm))
  expect_identical(stats::nobs(short), stats::nobs(fm))
  expect_gt(stats::deviance(short) - stats::deviance(direct), 10)
  testthat::local_mocked_bindings(refitML = mr16_refit_is(short), .package = "lme4")
  row <- mr16_heard(maihda_ic_one(m, ml = TRUE))
  expect_identical(row$warnings, character(0))
  expect_lt(abs(row$value$AIC - stats::AIC(direct)), 1e-6)
  expect_identical(row$value$n, stats::nobs(fm))
  expect_lt(abs(stats::deviance(maihda_lmer_refit_ml(fm)) - stats::deviance(direct)), 1e-6)
})

test_that("a vector outside the data that the call reads, changed since the model was fitted, does not change the criteria", {
  skip_on_cran()
  d <- mr16_data(3, "z")
  # An offset the outcome does not follow, held in this environment and not in the data.
  v <- rep(c(-2, 2, 0, 1), length.out = nrow(d))
  m <- mr16_quiet(fit_maihda(y ~ x + offset(v) + (1 + x | stratum), data = d))
  fm <- m$model
  refit <- lme4::refitML(fm)
  expect_identical(unname(lme4::getME(refit, "offset")), v)
  expect_true(maihda_lmer_same_model(mr16_quiet(maihda_lmer_ml_from_call(fm)), refit))
  # The vector changes: the call now gives the model without the offset, with a far lower
  # deviance. lme4's refit carries the offset the model was fitted with.
  v <- rep(0, nrow(d))
  changed <- mr16_quiet(maihda_lmer_ml_from_call(fm))
  expect_false(lme4::isREML(changed))
  expect_lt(stats::deviance(changed), stats::deviance(refit) - 100)
  expect_false(maihda_lmer_same_model(changed, refit))
  row <- mr16_heard(maihda_ic_one(m, ml = TRUE))
  expect_identical(row$warnings, character(0))
  expect_identical(row$value$AIC, as.numeric(stats::AIC(refit)))
  expect_identical(row$value$logLik, as.numeric(stats::logLik(refit)))
  kept <- maihda_lmer_refit_ml(fm)
  expect_identical(stats::deviance(kept), stats::deviance(refit))
  expect_identical(unname(lme4::getME(kept, "offset")), rep(c(-2, 2, 0, 1), length.out = nrow(d)))
})

test_that("a refit with an exit code is checked by lmer's default optimiser and by the model's own fit", {
  skip_on_cran()
  fits <- mr16_fits("z")
  fm <- fits$null$model
  real <- lme4::refitML
  own <- real(fm)
  calls <- character(0)
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    calls <<- c(calls, optimizer)
    refit <- real(x, optimizer = optimizer, ...)
    if (identical(optimizer, "bobyqa")) {
      refit@optinfo$conv$opt <- 3L
      warning("convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q")
    }
    refit
  }, .package = "lme4")
  made <- 0L
  real_from_call <- maihda_lmer_ml_from_call
  testthat::local_mocked_bindings(maihda_lmer_ml_from_call = function(fit, longitudinal = NULL) {
    made <<- made + 1L
    real_from_call(fit, longitudinal)
  }, .package = "MAIHDA")
  h <- mr16_heard(maihda_lmer_refit_ml(fm))
  expect_identical(calls, c("bobyqa", "nloptwrap"))
  expect_identical(made, 1L)
  # The flagged refit stands at the maximum: it is returned, without lme4's warning.
  expect_identical(h$warnings, character(0))
  expect_identical(stats::deviance(h$value), stats::deviance(own))
  expect_identical(lme4::getME(h$value, "theta"), lme4::getME(own, "theta"))
})

# ---- the criteria ---------------------------------------------------------------------------

test_that("maihda_ic() reads the criteria of models with a slope on an income from their maximum-likelihood fits", {
  skip_on_cran()
  # On the builds measured lme4's refit of these two fits ends 0.027 and 1.95 above the maximum,
  # the adjusted model's with a stratum intercept variance of 0.34 where the maximum has 0.68.
  fits <- mr16_fits("income")
  ic <- mr16_heard(maihda_ic(fits$null, fits$adj))
  expect_identical(ic$warnings, character(0))
  expect_identical(ic$messages, character(0))
  expect_identical(ic$value$estimator, rep("ML (refit from REML)", 2L))
  expect_lt(max(abs(-2 * ic$value$logLik - optimum_income)), 0.01)
  expect_lt(max(abs(ic$value$AIC - (optimum_income + 2 * c(6, 11)))), 0.01)
  expect_lt(max(abs(ic$value$BIC - (optimum_income + log(720) * c(6, 11)))), 0.01)
  # The adjusted model is the better one, by the difference of the criteria at the maxima.
  expect_lt(max(abs(ic$value$delta - c(optimum_income[1] - optimum_income[2] - 10, 0))), 0.01)
  # The variances are those of the maximum.
  adjusted <- mr16_quiet(maihda_lmer_refit_ml(fits$adj$model))
  expect_lt(abs(lme4::VarCorr(adjusted)$stratum[1, 1] - 0.6765), 0.01)
})

test_that("maihda_ic() reads the criteria of models with slopes on an age and a birth year from their maximum-likelihood fits", {
  skip_on_cran()
  # On the builds measured lme4's refit of these two fits ends 156.4 and 159.4 above the maximum.
  d <- mr16_data(1, "age", n = 60)
  m0 <- mr16_quiet(fit_maihda(y ~ x + yr + (1 + x + yr | stratum), data = d))
  m1 <- mr16_quiet(fit_maihda(y ~ x + yr + a + b + c + (1 + x + yr | stratum), data = d))
  ic <- mr16_heard(maihda_ic(m0, m1))
  expect_identical(ic$warnings, character(0))
  expect_identical(ic$messages, character(0))
  expect_lt(max(abs(-2 * ic$value$logLik - optimum_age)), 0.01)
  expect_lt(max(abs(ic$value$AIC - (optimum_age + 2 * c(10, 15)))), 0.01)
  expect_lt(max(abs(ic$value$delta - c(optimum_age[1] - optimum_age[2] - 10, 0))), 0.01)
})

test_that("maihda_ic() and compare_maihda() read a growth pair's criteria from its maximum-likelihood fits", {
  skip_on_cran()
  # On the builds measured lme4's refit of the adjusted model ends 0.112 above the maximum,
  # with the stratum intercept variance at zero.
  long <- mr16_fits("long")
  ic <- mr16_heard(maihda_ic(long$null, long$adj))
  expect_identical(ic$warnings, character(0))
  expect_lt(max(abs(-2 * ic$value$logLik - optimum_long)), 0.01)
  expect_lt(max(abs(ic$value$AIC - (optimum_long + 2 * c(9, 13)))), 0.01)
  cmp <- mr16_quiet(compare_maihda(long$null, long$adj))
  expect_identical(cmp$AIC, ic$value$AIC)
  expect_identical(cmp$BIC, ic$value$BIC)
})

# ---- the PCVs of estimation = "ML" ----------------------------------------------------------

test_that("the growth PCV of estimation = \"ML\" is read from the maximum-likelihood fits", {
  skip_on_cran()
  long <- mr16_fits("long")
  real <- lme4::refitML
  # lme4's refit of the adjusted model with the stratum intercept variance at zero, where it
  # ends on the builds measured; its refit of the null model as it is.
  at_zero <- function(fit) {
    theta <- unname(lme4::getME(real(fit), "theta"))
    stratum <- which(names(lme4::getME(fit, "cnms")) == "stratum")
    expect_identical(stratum, 2L)
    theta[4] <- 0
    mr16_placed(fit, theta)
  }
  zero <- at_zero(long$adj$model)
  expect_identical(lme4::VarCorr(zero)$stratum[1, 1], 0)
  expect_gt(stats::deviance(zero) - optimum_long[2], 0.01)
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    if (length(lme4::fixef(x)) == length(lme4::fixef(long$adj$model))) zero else real(x)
  }, .package = "lme4")
  pcv <- mr16_heard(maihda_longitudinal_pcv(long$null, long$adj, estimation = "ML"))
  expect_identical(pcv$warnings, character(0))
  expect_identical(pcv$value$estimation_used, "ML")
  expect_true(pcv$value$ml_refit)
  # Not the PCV of exactly 1 that a variance of zero gives, but that of the two models fitted
  # by maximum likelihood.
  expect_gt(pcv$value$var_baseline_adjusted, 1e-3)
  expect_lt(pcv$value$pcv_intercept, 0.999)
  d <- long$data
  direct <- mr16_quiet(maihda_longitudinal_pcv(
    fit_maihda(y ~ time + (1 | stratum), data = d, id = "pid", time = "time", REML = FALSE),
    fit_maihda(y ~ time + a + b + (1 | stratum), data = d, id = "pid", time = "time",
               REML = FALSE)))
  expect_lt(abs(pcv$value$pcv_intercept - direct$pcv_intercept), 1e-3)
  expect_lt(abs(pcv$value$var_baseline_adjusted - direct$var_baseline_adjusted), 1e-4)
  expect_lt(max(abs(pcv$value$pcv_t$pcv - direct$pcv_t$pcv)), 1e-3)
  # The stored models keep their REML fits.
  expect_true(lme4::isREML(long$adj$model))
})

test_that("the growth PCV of estimation = \"ML\" does not report a variance of zero that lme4's refit leaves", {
  skip_on_cran()
  # Nothing replaced: on the builds measured lme4's own refit of the adjusted model leaves the
  # stratum intercept variance at zero, 0.112 above the maximum, where it is 0.0026.
  long <- mr16_fits("long")
  pcv <- mr16_quiet(maihda_longitudinal_pcv(long$null, long$adj, estimation = "ML"))
  expect_gt(pcv$var_baseline_adjusted, 1e-3)
  expect_lt(pcv$pcv_intercept, 0.999)
  fitted <- mr16_quiet(maihda_longitudinal_pcv(long$null, long$adj))
  expect_identical(fitted$estimation_used, "fitted")
  expect_false(fitted$ml_refit)
})

test_that("calculate_pcv(estimation = \"ML\") does not pass on lme4's exit code for a refit at its maximum", {
  skip_on_cran()
  plain <- mr16_fits("plain")
  reference <- mr16_heard(calculate_pcv(plain$null, plain$adj, estimation = "ML"))
  expect_identical(reference$warnings, character(0))
  expect_identical(reference$value$estimation_used, "ML")
  real <- lme4::refitML
  calls <- character(0)
  testthat::local_mocked_bindings(refitML = function(x, optimizer = "bobyqa", ...) {
    calls <<- c(calls, optimizer)
    refit <- real(x, optimizer = optimizer, ...)
    if (identical(optimizer, "bobyqa")) {
      refit@optinfo$conv$opt <- 3L
      warning("convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q")
    }
    refit
  }, .package = "lme4")
  pcv <- mr16_heard(calculate_pcv(plain$null, plain$adj, estimation = "ML"))
  expect_identical(calls, rep(c("bobyqa", "nloptwrap"), 2L))
  expect_identical(pcv$warnings, character(0))
  expect_identical(pcv$value$pcv, reference$value$pcv)
  expect_identical(pcv$value$var_model1, reference$value$var_model1)
  expect_identical(pcv$value$var_model2, reference$value$var_model2)
  expect_identical(pcv$value$estimation_used, "ML")
  # The variances are those of lme4's refits.
  expect_identical(pcv$value$var_model2,
                   as.numeric(lme4::VarCorr(real(plain$adj$model))$stratum[1, 1]))
  # A refit that fails is still said to have failed.
  testthat::local_mocked_bindings(refitML = function(x, ...) stop("the refit fails"),
                                  .package = "lme4")
  failed <- mr16_heard(calculate_pcv(plain$null, plain$adj, estimation = "ML"))
  expect_true(any(grepl("lme4::refitML() failed", failed$warnings, fixed = TRUE)))
  expect_identical(failed$value$estimation_used, "mixed")
})
