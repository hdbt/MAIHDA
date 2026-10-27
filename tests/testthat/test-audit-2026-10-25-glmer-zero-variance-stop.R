# Audit pass 2026-10-25: glmer fits with a random slope that stop short of the maximum, often with
# what is left of the slope's variance at zero, recorded as converged.
#
# glmer's second stage, Nelder-Mead, takes its first steps 0.02 from where its first stage stopped
# in every entry of theta. A slope on a covariate far from zero or in large units, such as a birth
# or calendar year, has entries of about 1e-4 in its column, which such steps cannot resolve, and
# the fit stops where the first stage left them: often with what is left of the slope's variance
# after the intercept at zero, where lme4 calls the fit singular and runs no convergence check, and
# at times short of the maximum without a variance at zero. Per the user, fit_maihda() now searches
# a glmer fit with a random slope with each random-effect block's terms re-expressed on their own
# and the fixed effects centred and scaled, from lme4's starting values and from the fit's
# estimates, and refits from the lower end by glmer's second stage with bobyqa, its first step
# below the end's smallest positive variance entry.
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

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 0.3) of 80 rows, each row with a
# birth year between 1940 and 2000 and a count whose stratum slope deviates with SD 0.005 per year
# around 1970: the Poisson data of test-audit-2026-10-20-glmer-slope-scale.R.
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

slope <- y ~ x + (1 + x | stratum)

# Seed 1: glmer's own fit stops with the slope's variance after the intercept at zero; the maximum
# of the Laplace likelihood, from an unbounded search of glmer's deviance function (PIRLS to 1e-12)
# from many starts, lies 2.4634 lower in -2 log-likelihood. Both rounded to seven digits; the gap
# is the same on five orders of the rows.
stuck_1 <- list(theta = c(6.415436, -0.003434893, 0), beta = c(-21.6468, 0.01131645))
max_1 <- list(theta = c(9.80502, -0.005099931, 0.0001291164), beta = c(-21.60469, 0.01129371))
deviance_max_1 <- 4935.556312

# A glmer fit placed at a point: optimisers that return their start, and no move of an entry near
# zero onto zero (lme4's boundary.tol).
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}
placed <- function(point, data, formula = slope) {
  quiet(lme4::glmer(formula, data = data, family = stats::poisson,
                    start = list(theta = point$theta, fixef = point$beta),
                    control = lme4::glmerControl(optimizer = list(stay, stay),
                                                 nAGQ0initStep = FALSE, boundary.tol = 0)))
}
deviance_of <- function(m) -2 * as.numeric(stats::logLik(m))

# ---- the units of the search ---------------------------------------------------------------

test_that("a block's factor is carried to its re-expressed terms without forming the covariance", {
  cnms <- list(stratum = c("(Intercept)", "x"))
  Tm <- matrix(c(1, 0, -1970 / 17.6, 1 / 17.6), 2)
  factor_of <- function(theta) matrix(c(theta[1], theta[2], 0, theta[3]), 2)
  theta <- c(6.4, -0.0034, 0.00013)
  carried <- maihda_glmer_carry_theta(theta, cnms, list(Tm))
  # The covariance of T b, a non-negative diagonal, and back again.
  expect_equal(tcrossprod(factor_of(carried)), Tm %*% tcrossprod(factor_of(theta)) %*% t(Tm))
  expect_true(all(carried[c(1, 3)] >= 0))
  expect_equal(maihda_glmer_carry_theta(carried, cnms, list(solve(Tm))), theta, tolerance = 1e-10)
  # A block of rank one whose intercept is near zero at the mean year, carried to the year centred
  # and scaled, keeps its zero.
  rank_one <- c(1970 * 0.0034 + 2e-5, -0.0034, 0)
  expect_identical(maihda_glmer_carry_theta(rank_one, cnms, list(solve(Tm)))[3], 0)
  # A theta that is not a full factor per block, as lme4 holds a diag() block.
  expect_null(maihda_glmer_carry_theta(c(1, 2), cnms, list(Tm)))
})

test_that("the re-expressed glmer fit is the same model", {
  skip_on_cran()
  d <- year_counts(1)
  m <- placed(stuck_1, d)
  re <- maihda_glmer_reexpressed(m)
  # At the fit, the fit's own -2 log-likelihood (to the precision of lme4's PIRLS, which stops at
  # a relative change of 1e-7); elsewhere, its deviance function at the point carried back.
  expect_lt(abs(re$deviance(re$p) - deviance_of(m)), 1e-3)
  devfun <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
  p <- re$p + c(0.01, -0.02, 0.03, 0.05, -0.01)
  back <- re$carry(p)
  expect_equal(re$deviance(p), devfun(c(back$theta, unname(back$fixef))))
  # The block by the year's mean and SD, the coefficients of the year centred and scaled.
  expect_equal(re$scales[[1]], matrix(c(1, 0, -mean(d$x) / stats::sd(d$x), 1 / stats::sd(d$x)), 2))
  expect_equal(re$carry(re$p)$fixef, lme4::fixef(m))
  expect_equal(re$carry(re$p)$theta, unname(lme4::getME(m, "theta")), tolerance = 1e-8)
  expect_identical(re$lower, c(0, -Inf, 0, -Inf, -Inf))
  # Not made at nAGQ 0, whose fixed effects were not optimised with the variances, nor where the
  # deviance function does not give the fit's own -2 log-likelihood.
  m0 <- quiet(lme4::glmer(slope, data = d, family = stats::poisson, nAGQ = 0))
  expect_null(maihda_glmer_reexpressed(m0))
  local_mocked_bindings(maihda_glmer_criterion = function(model) {
    -2 * as.numeric(stats::logLik(model)) + 10
  })
  expect_null(maihda_glmer_reexpressed(m))
})

# ---- the searches --------------------------------------------------------------------------

test_that("the searches start from lme4's starting values and from the fit's estimates", {
  skip_on_cran()
  d <- year_counts(1)
  m <- placed(stuck_1, d)
  re <- maihda_glmer_reexpressed(m)
  seen <- list()
  stand_in <- function(par, fn, lower, upper) {
    seen[[length(seen) + 1L]] <<- list(par = par, lower = lower)
    list(par = par, fval = fn(par))
  }
  maihda_glmer_search(m, "standardised", optimiser = stand_in)
  # theta at lme4's start (1 on the diagonal), the fixed effects the fit's, within lme4's bounds;
  # then once more from where it stopped.
  expect_length(seen, 2L)
  expect_equal(seen[[1]]$par, c(1, 0, 1, re$p[4:5]))
  expect_identical(seen[[1]]$lower, re$lower)
  expect_equal(seen[[2]]$par, seen[[1]]$par)
  seen <- list()
  found <- maihda_glmer_search(m, "estimates", optimiser = stand_in)
  expect_equal(seen[[1]]$par, re$p)
  expect_lt(abs(found$gain), 1e-8)
  expect_identical(found$lower, c(0, -Inf, 0))
  # A stop across from a bound: once more from the other side.
  seen <- list()
  trapped <- function(par, fn, lower, upper) {
    seen[[length(seen) + 1L]] <<- list(par = par)
    if (length(seen) == 1L) par[1] <- 0
    list(par = par, fval = fn(par))
  }
  maihda_glmer_search(m, "estimates", optimiser = trapped)
  expect_length(seen, 3L)
  expect_identical(seen[[2]]$par[1], 0.05)
  expect_equal(seen[[2]]$par[2], -re$p[2])
  # Intercepts only: nothing to re-express, so no search from lme4's starting values.
  io <- quiet(lme4::glmer(y ~ x + (1 | stratum), data = d, family = stats::poisson))
  expect_null(maihda_glmer_search(io, "standardised"))
})

test_that("searched from glmer's stuck fit, the deviance reaches the maximum", {
  skip_on_cran()
  d <- year_counts(1)
  m <- placed(stuck_1, d)
  # The premise: lme4 calls the fit singular, says nothing else, and it is 2.46 below the maximum.
  expect_true(lme4::isSingular(m))
  expect_identical(as.character(m@optinfo$conv$lme4$messages),
                   "boundary (singular) fit: see help('isSingular')")
  expect_gt(deviance_of(m) - deviance_max_1, 2.46)
  for (from in c("standardised", "estimates")) {
    found <- maihda_glmer_search(m, from)
    expect_gt(found$gain, 2.46)
    expect_lt(found$deviance - deviance_max_1, 1e-3)
  }
  # At the maximum a search gains nothing.
  expect_lt(maihda_glmer_search(placed(max_1, d), "estimates")$gain, 1e-4)
})

# The counts of test-audit-2026-10-26-tied-block-order.R's glmer check: 600 rows, a slope on a
# birth year in the block of g1, beside g2 with as many levels and gbig with 30. lme4 lists the
# term columns of such tied blocks, and of the two blocks (1 + yr || g1) makes behind gbig's, in
# another order than their variance parameters.
tied_counts <- function() {
  with_seed(11, {
    n <- 600
    d <- data.frame(g1 = factor(sample(8, n, TRUE)), g2 = factor(sample(8, n, TRUE)),
                    gbig = factor(sample(30, n, TRUE)), x = round(stats::rnorm(n), 2))
    zlev <- stats::rnorm(8)
    d$z <- zlev[d$g1]
    d$yr <- round(stats::runif(n, 1940, 2000))
    u1 <- stats::rnorm(8, 0, 0.5)
    v1 <- stats::rnorm(8, 0, 0.3)
    u2 <- stats::rnorm(8, 0, 0.4)
    ub <- stats::rnorm(30, 0, 0.3)
    d$y <- 1 + 0.5 * d$x + 0.3 * d$z + u1[d$g1] + v1[d$g1] * d$x + u2[d$g2] + ub[d$gbig] +
      stats::rnorm(n)
    d$cnt <- stats::rpois(n, exp(0.2 + 0.2 * d$z + 0.4 * u1[d$g1] + 0.004 * (d$yr - 1970) +
                                   stats::rnorm(8, 0, 0.003)[d$g1] * (d$yr - 1970) +
                                   0.3 * u2[d$g2] + 0.3 * ub[d$gbig]))
    d
  })
}

test_that("a glmer fit whose blocks lme4 lists in another order is searched as the sorted spelling is", {
  skip_on_cran()
  d <- tied_counts()
  blocks <- function(model) {
    cnms <- lme4::getME(model, "cnms")
    paste(names(cnms), vapply(cnms, paste, character(1), collapse = " + "))
  }
  # Each pair placed at lme4's own fit of its first spelling (to seven digits), theta in the
  # order of each spelling's blocks.
  pairs <- list(
    list(spelt = cnt ~ yr + z + (1 + yr | g1) + (1 | g2) + (1 | gbig),
         sorted = cnt ~ yr + z + (1 | gbig) + (1 | g2) + (1 + yr | g1),
         theta = c(0, 0.1587222, 9.296477, -0.004574353, 0),
         theta_sorted = c(0, 0.1587222, 9.296477, -0.004574353, 0),
         beta = c(-3.587002, 0.001972961, 0.3291751)),
    list(spelt = cnt ~ yr + z + (1 + yr || g1) + (1 | gbig),
         sorted = cnt ~ yr + z + (1 | gbig) + (1 | g1) + (0 + yr | g1),
         theta = c(0.05508204, 0, 0.262959),
         theta_sorted = c(0.05508204, 0.262959, 0),
         beta = c(-1.025796, 0.0006771691, 0.3432832)))
  for (pair in pairs) {
    m <- placed(list(theta = pair$theta, beta = pair$beta), d, formula = pair$spelt)
    s <- placed(list(theta = pair$theta_sorted, beta = pair$beta), d, formula = pair$sorted)
    expect_lt(abs(deviance_of(m) - deviance_of(s)), 1e-6)
    # Paired by position, the term columns do not fit the blocks.
    expect_null(maihda_lmer_block_scales(lme4::getME(m, "mmList"), lme4::getME(m, "cnms")))
    # Each block re-expressed as in the sorted spelling, and both searches made.
    re <- maihda_glmer_reexpressed(m)
    expect_false(is.null(re))
    expect_equal(unname(re$scales),
                 unname(maihda_glmer_reexpressed(s)$scales[match(blocks(m), blocks(s))]))
    stand_in <- function(par, fn, lower, upper) list(par = par, fval = fn(par))
    expect_false(is.null(maihda_glmer_search(m, "standardised", optimiser = stand_in)))
    expect_lt(abs(maihda_glmer_search(m, "estimates")$deviance -
                    maihda_glmer_search(s, "estimates")$deviance), 1e-6)
  }
})

# ---- the refit -----------------------------------------------------------------------------

test_that("the refit's first step lies below the smallest positive variance entry where the search ended", {
  lower <- c(0, -Inf, 0)
  # Half the smallest positive diagonal entry: bobyqa would move a slope's entry of 1.3e-4 on a
  # birth year to glmer's own first step of 2e-4 before it starts.
  expect_equal(maihda_glmer_first_step(c(9.8, -0.0051, 0.000129), lower), 0.0000645)
  # Entries below the diagonal and zeros do not count; at most 2e-4.
  expect_identical(maihda_glmer_first_step(c(0.5, -1e-6, 0), lower), 2e-4)
  expect_identical(maihda_glmer_first_step(c(0.8, 0.1, 0.3), lower), 2e-4)
  expect_identical(maihda_glmer_first_step(c(0, 0.1, 0), lower), 2e-4)
})

test_that("a fit with a random slope is refitted from the lower end of its searches, kept only when lower by more than the tolerance", {
  env <- new.env()
  env$crit <- c(first = 10, refit = 9)
  env$fake_fit <- function(start = NULL, control = NULL) {
    message(if (is.null(start)) "first fit says" else "refit says")
    if (identical(start$fixef, c(b = -1))) stop("this refit fails")
    list(crit = env$crit[[if (is.null(start)) "first" else "refit"]], start = start, control = control)
  }
  end_at <- function(theta, deviance, gain, fixef = c(b = 1)) {
    list(start = list(theta = theta, fixef = fixef), deviance = deviance, gain = gain,
         lower = c(0, -Inf, 0))
  }
  ends <- list(standardised = end_at(c(0.5, 0.1, 3e-4), 9.2, 0.8),
               estimates = end_at(c(0.6, 0.1, 1e-4), 9.1, 0.9))
  asked <- character(0)
  slope_held <- TRUE
  local_mocked_bindings(
    maihda_glmer_mirror_refit = function(...) NULL,
    maihda_lme4_has_slope = function(model) slope_held,
    maihda_glmer_criterion = function(model) model$crit,
    maihda_glmer_search = function(model, from, ...) {
      asked <<- c(asked, from)
      ends[[from]]
    }
  )
  fit <- function() heard(maihda_fit_glmer_bound(quote(fake_fit()), env))
  # Both searches; the refit from the lower end, by glmer's second stage alone with bobyqa, its
  # first step half the end's smallest positive variance entry.
  h <- fit()
  expect_identical(asked, c("standardised", "estimates"))
  expect_identical(h$value$start, list(theta = c(0.6, 0.1, 1e-4), fixef = c(b = 1)))
  expect_identical(h$value$control$optimizer, c("bobyqa", "bobyqa"))
  expect_false(h$value$control$nAGQ0initStep)
  expect_equal(h$value$control$optCtrl, list(rhobeg = 5e-5, rhoend = 5e-8, maxfun = 1e5))
  expect_identical(h$seen, "refit says\n")
  # With a control of the user's own: its first stage and other settings kept, its optimiser's
  # settings dropped.
  own <- quiet(maihda_fit_glmer_bound(quote(fake_fit(control = lme4::glmerControl(
    optimizer = c("Nelder_Mead", "Nelder_Mead"), optCtrl = list(maxfun = 3e4), tolPwrss = 1e-9))), env))
  expect_identical(own$control$optimizer, c("Nelder_Mead", "bobyqa"))
  expect_equal(own$control$optCtrl, list(rhobeg = 5e-5, rhoend = 5e-8, maxfun = 1e5))
  expect_identical(own$control$tolPwrss, 1e-9)
  # A refit lower by less than 1e-4, a search that gains no more than that, a refit that fails:
  # the first fit stays.
  env$crit[["refit"]] <- 10 - 5e-5
  expect_null(fit()$value$start)
  env$crit[["refit"]] <- 9
  ends <- list(standardised = end_at(c(0.5, 0.1, 3e-4), 9.99995, 5e-5), estimates = NULL)
  expect_null(fit()$value$start)
  ends <- list(standardised = end_at(c(0.5, 0.1, 3e-4), 9, 1, fixef = c(b = -1)), estimates = NULL)
  h <- fit()
  expect_null(h$value$start)
  expect_identical(h$seen, "first fit says\n")
  # Intercepts only: no search.
  asked <- character(0)
  slope_held <- FALSE
  expect_null(fit()$value$start)
  expect_identical(asked, character(0))
})

test_that("the record searches a singular glmer slope fit unless fit_maihda() refitted it from its searches", {
  skip_on_cran()
  m <- placed(stuck_1, year_counts(1))
  calls <- 0L
  real <- maihda_glmer_search_gain
  local_mocked_bindings(maihda_glmer_search_gain = function(model, ...) {
    calls <<- calls + 1L
    real(model, ...)
  })
  # Made otherwise, as from a start of the user's own: searched, and found short.
  open <- maihda_fit_diagnostics(m)
  expect_identical(calls, 1L)
  expect_false(open$converged)
  expect_gt(open$search_gain, 0.01)
  # Refitted by fit_maihda() from its searches already: not searched again.
  done <- maihda_fit_diagnostics(m, searched = TRUE)
  expect_identical(calls, 1L)
  expect_null(done$search_gain)
})

# ---- live fits -----------------------------------------------------------------------------

test_that("fit_maihda() reaches the maximum of a birth-year slope fit glmer leaves at a variance of zero", {
  skip_on_cran()
  skip_slow_test()
  d <- year_counts(1)
  own <- quiet(lme4::glmer(slope, data = d, family = stats::poisson))
  expect_gt(deviance_of(own) - deviance_max_1, 2.46)
  told <- NULL
  real <- maihda_fit_diagnostics
  local_mocked_bindings(maihda_fit_diagnostics = function(model, longitudinal_info = NULL, searched = FALSE) {
    told <<- searched
    real(model, longitudinal_info, searched)
  })
  m <- quiet(fit_maihda(slope, data = d, family = "poisson"))
  # The record is told the fit comes from the searches.
  expect_true(told)
  expect_lt(abs(deviance_of(m$model) - deviance_max_1), 1e-4)
  expect_false(m$diagnostics$singular)
  expect_true(m$diagnostics$converged)
  # The fit kept is a glmer fit of the formula, refitted by bobyqa from the search's end and then
  # by glmer's own second stage from there: it prints, update() makes it again, and refit() with
  # glmer's default control stays there.
  expect_identical(m$model@call$control, as.name(".maihda_arg_finish_control"))
  expect_identical(m$model@optinfo$optimizer, "Nelder_Mead")
  expect_output(print(m))
  expect_output(print(m$model))
  expect_equal(deviance_of(quiet(stats::update(m$model))), deviance_of(m$model), tolerance = 1e-8)
  refitted <- quiet(lme4::refit(m$model, control = lme4::glmerControl()))
  expect_lt(abs(deviance_of(refitted) - deviance_of(m$model)), 1e-3)
})

test_that("a negative-binomial fit is refitted the same way", {
  skip_on_cran()
  skip_slow_test()
  # Negative-binomial counts (size 20) of the same design with a calendar year between 2010 and
  # 2020: glmer.nb's own fit stops singular at 5223.46 in -2 log-likelihood, 8.4 above the fit
  # fit_maihda() refits to.
  d <- with_seed(6, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, 0.01)
    st <- rep(seq_len(18), each = 80)
    x <- sample(2010:2020, length(st), replace = TRUE)
    xc <- (x - 2015) / 3.16 * 17.6
    eta <- 0.5 + add[st] + (0.01 + v[st] / 2) * xc
    data.frame(y = stats::rnbinom(length(st), mu = exp(eta), size = 20), x = x, stratum = factor(st))
  })
  m <- quiet(fit_maihda(slope, data = d, family = "negbinomial"))
  expect_lt(deviance_of(m$model), 5216)
  expect_identical(m$model@optinfo$optimizer, "Nelder_Mead")
  expect_output(print(m))
})

test_that("an intercept-only glmer fit is not searched", {
  skip_on_cran()
  searched <- 0L
  told <- NULL
  real <- maihda_fit_diagnostics
  local_mocked_bindings(
    maihda_glmer_search = function(model, from, ...) {
      searched <<- searched + 1L
      NULL
    },
    maihda_fit_diagnostics = function(model, longitudinal_info = NULL, searched = FALSE) {
      told <<- searched
      real(model, longitudinal_info, searched)
    }
  )
  d <- year_counts(1)
  m <- quiet(fit_maihda(y ~ x + (1 | stratum), data = d, family = "poisson"))
  direct <- quiet(lme4::glmer(y ~ x + (1 | stratum), data = d, family = stats::poisson))
  expect_identical(searched, 0L)
  expect_false(told)
  expect_identical(lme4::getME(m$model, "theta"), lme4::getME(direct, "theta"))
  expect_identical(lme4::fixef(m$model), lme4::fixef(direct))
})
