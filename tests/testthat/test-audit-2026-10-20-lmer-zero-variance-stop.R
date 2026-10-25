# Audit pass 2026-10-20: lmer fits stopped with a variance at zero, or beside one, short of the
# REML optimum, which lme4 calls singular and does not check.
#
# In a random-effect column with nothing below it -- a block's last term, or a block of one term
# such as (0 + x | ctx) or the slope of (1 + x || stratum) -- the criterion depends on the Cholesky
# diagonal entry only through its square, so at zero it is flat to first order, and lmer's
# optimiser can stop there while the criterion falls as that variance moves off zero. With a slope
# covariate far from zero (a calendar or birth year) in a block without an intercept, lmer also
# stops on the ridge between that block and the intercept's. lme4 calls such fits singular and runs
# no convergence check on them. Now fit_maihda() refits from a start with each random-effect
# block's terms re-expressed on their own (a term in a block without an intercept divided by its
# root mean square, not by its standard deviation, which leaves a calendar year at about 630), and
# from a variance at zero moved to the lowest criterion along it, from the kept fit and from lmer's
# own; and the record searches a singular lmer fit it did not refit that way, not converged when a
# search reaches a criterion 0.01 or more below it.
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

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 0.3) of `n` rows, each row with a
# calendar year (2010..2020), a birth year (1940..2000) or an age (18..80) `x`, a stratum slope
# deviation of SD `slope_sd` per unit of x around the middle of its range, and a residual SD of 1;
# with `z`, a second covariate of that kind with slope deviations of SD 0.005; with
# `ctx_slope_sd`, one of six contexts with slope deviations of that SD.
strata_data <- function(seed, x = "cal", slope_sd = 0.05, beta = 0.1, n = 40, z = NULL,
                        ctx_slope_sd = NULL) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, slope_sd)
    st <- rep(seq_len(18), each = n)
    draw <- function(type) switch(type, cal = sample(2010:2020, length(st), replace = TRUE),
                                  year = sample(1940:2000, length(st), replace = TRUE),
                                  age = sample(18:80, length(st), replace = TRUE))
    centre <- c(cal = 2015, year = 1970, age = 49)
    xs <- draw(x)
    y <- 1 + add[st] + (beta + v[st]) * (xs - centre[[x]])
    out <- data.frame(x = xs, d1 = factor(lev$d1[st]), d2 = factor(lev$d2[st]),
                      d3 = factor(lev$d3[st]), stratum = factor(st))
    if (!is.null(z)) {
      out$z <- draw(z)
      y <- y + (0.2 + stats::rnorm(18, 0, 0.005)[st]) * (out$z - centre[[z]])
    }
    if (!is.null(ctx_slope_sd)) {
      ctx <- sample(1:6, length(st), replace = TRUE)
      y <- y + stats::rnorm(6, 0, ctx_slope_sd)[ctx] * (xs - centre[[x]])
      out$ctx <- factor(ctx)
    }
    out$y <- y + stats::rnorm(length(st))
    out
  })
}

uncorrelated <- y ~ x + (1 + x || stratum)
by_context <- y ~ x + (1 + x | stratum) + (0 + x | ctx)
slope_only <- y ~ x + (1 | stratum) + (0 + x | stratum)
three <- y ~ x + z + d1 + d2 + d3 + (1 + x + z | stratum)

# REML optima from lmer's deviance function searched with no bound on the Cholesky diagonal
# (nlminb from thirty and more starts), and where lmer's own fit stops:
# - uncorrelated, seed 11: 2107.81576875 at theta (3.3286e-05, 0.000211307), the slope carrying
#   the strata's variance; lmer stops at (0.425751, 0), the intercept carrying it, 0.0294 above,
#   where the criterion falls as the slope's variance moves off zero.
# - uncorrelated, seed 1: 2179.85291366 at (0.300031, 0), where it rises off zero.
# - by_context, seed 9: 2112.83311693; lmer stops 12.8 above with the stratum slope's variance
#   (what is left of it after the intercept) and the context's at zero.
# - slope_only, seed 1 (a birth year): 2182.20733012; lmer stops 0.051 above at (0.338996,
#   1.03135e-05), and a refit from a start with the year divided by its standard deviation at
#   (0.337278, 1.77235e-05), the slope 0.035 in its block's root-mean-square units: on the ridge,
#   not at zero.
# - three, seed 2 (age and birth year, 60 rows a stratum): 3118.90149579, the last column at zero
#   there too; lmer stops 6.0 above.
# - by_context with the dimensions as fixed effects, seed 7: 2150.74922730; lmer stops 5.8 above.
optimum_unc_11 <- 2107.81576875
stuck_unc_11 <- c(0.425751, 0)
at_unc_11 <- c(3.3286e-05, 0.000211307)
boundary_unc_1 <- c(0.300031, 0)
optimum_ctx_9 <- 2112.83311693
optimum_slope_only_1 <- 2182.20733012
stuck_slope_only_1 <- c(0.337278, 1.77235e-05)
optimum_three_2 <- 3118.90149579
optimum_ctx_adjusted_7 <- 2150.74922730

# An lmer fit placed at `theta`: an optimiser that returns its start, and no move of an entry
# near zero onto zero (lme4's boundary.tol).
stay <- function(fn, par, lower, upper, control) {
  list(par = par, fval = fn(par), conv = 0, feval = 1, message = "stay")
}
placed <- function(theta, data, formula, REML = TRUE) {
  quiet(lme4::lmer(formula, data = data, REML = REML, start = theta,
                   control = lme4::lmerControl(optimizer = stay, boundary.tol = 0)))
}

# ---- the start, block by block ------------------------------------------------------------

test_that("the standardised start re-expresses each block on its own, a block without an intercept by its root mean square", {
  skip_on_cran()
  d <- strata_data(11)
  own <- quiet(lme4::lmer(uncorrelated, data = d))
  calls <- list()
  stand_in <- function(par, fn, lower, upper) {
    calls[[length(calls) + 1L]] <<- list(par = par, lower = lower, value = fn(c(0.3, 0.2)))
    list(par = c(0.3, 0.2), fval = fn(c(0.3, 0.2)))
  }
  start <- maihda_lmer_standardised_start(own, optimiser = stand_in)
  dev <- lme4::lmer(uncorrelated, data = d, devFunOnly = TRUE)
  rms <- sqrt(mean(d$x^2))
  # From lme4's starting values within its bounds; the intercept block as it is, the slope block
  # divided by the root mean square of the year.
  expect_length(calls, 1L)
  expect_equal(calls[[1]]$par, c(1, 1))
  expect_equal(calls[[1]]$lower, c(0, 0))
  expect_equal(calls[[1]]$value, dev(c(0.3, 0.2 / rms)))
  expect_equal(start, list(c(0.3, 0.2 / rms)))
  # Divided by its standard deviation the year stays far from zero, and the start is another.
  expect_gt(mean(d$x) / stats::sd(d$x), 600)
  expect_gt(abs(calls[[1]]$value - dev(c(0.3, 0.2 / stats::sd(d$x)))), 1)
  # Missing values enter neither the root mean square nor a block's mean and SD.
  x <- c(2011, NA, 2014, 2016, 2019, 2020)
  Ts <- maihda_lmer_block_scales(list(cbind(`(Intercept)` = 1, x = x), cbind(x = x)),
                                 list(stratum = c("(Intercept)", "x"), ctx = "x"))
  xs <- x[!is.na(x)]
  expect_equal(Ts[[1]], matrix(c(1, 0, -mean(xs) / stats::sd(xs), 1 / stats::sd(xs)), 2))
  expect_equal(Ts[[2]], matrix(1 / sqrt(mean(xs^2))))
})

# ---- the start from a variance at zero ----------------------------------------------------

test_that("a variance at zero where the criterion falls off zero is moved to the lowest criterion along it", {
  skip_on_cran()
  d <- strata_data(11)
  stuck <- placed(stuck_unc_11, d, uncorrelated)
  dev <- lme4::lmer(uncorrelated, data = d, devFunOnly = TRUE)
  rms <- sqrt(mean(d$x^2))
  along <- function(t) dev(c(stuck_unc_11[1], t / rms))
  # The premise: the slope's variance at zero, 0.029 above the optimum, where the criterion falls
  # as it moves off zero.
  expect_lt(abs(dev(stuck_unc_11) - optimum_unc_11 - 0.0294), 1e-3)
  expect_lt(along(1e-3), along(0))
  handed <- NULL
  stand_in <- function(par, fn, lower, upper) {
    handed <<- par
    list(par = par, fval = fn(par))
  }
  maihda_lmer_zero_column_start(stuck, optimiser = stand_in)
  # The intercept is held; the slope, re-expressed, is moved to the lowest criterion along it.
  expect_equal(handed[1], stuck_unc_11[1])
  expect_gt(handed[2], 1e-3)
  expect_lt(along(handed[2]), along(1e-3))
  expect_lte(along(handed[2]), min(along(0.98 * handed[2]), along(1.02 * handed[2])))
  # A variance near zero, not at it, is moved as well.
  handed <- NULL
  maihda_lmer_zero_column_start(placed(c(stuck_unc_11[1], 5e-4 / rms), d, uncorrelated),
                                optimiser = stand_in)
  expect_true(isTRUE(handed[2] > 1e-3))
  # lmer refitted from that start, as maihda_fit_lmer_slopes() refits it, reaches the optimum.
  start <- maihda_lmer_zero_column_start(stuck)
  expect_length(start, 1L)
  refit <- quiet(lme4::lmer(uncorrelated, data = d, start = start[[1]]))
  expect_lt(abs(lme4::REMLcrit(refit) - optimum_unc_11), 1e-5)
  # With a context slope the start itself reaches it, from lmer's own fit 12.8 above with the
  # variances of the stratum slope (after the intercept) and of the context's at zero.
  dc <- strata_data(9, ctx_slope_sd = 0.05)
  stuck_ctx <- c(0.641034, -0.000493873, 0, 0)
  dev_ctx <- lme4::lmer(by_context, data = dc, devFunOnly = TRUE)
  expect_gt(dev_ctx(stuck_ctx) - optimum_ctx_9, 12)
  from_ctx <- maihda_lmer_zero_column_start(placed(stuck_ctx, dc, by_context))
  expect_lt(dev_ctx(from_ctx[[1]]) - optimum_ctx_9, 1e-5)
  # Nothing to move at the optimum, where the intercept's variance is near zero but the criterion
  # rises off it, nor where the slope's variance is zero at the optimum itself.
  expect_identical(maihda_lmer_zero_column_start(placed(at_unc_11, d, uncorrelated)), list())
  expect_identical(maihda_lmer_zero_column_start(placed(boundary_unc_1, strata_data(1), uncorrelated)),
                   list())
  # A deviance function that does not give the fit's own criterion.
  local_mocked_bindings(maihda_lmer_criterion = function(model) lme4::REMLcrit(model) + 1)
  expect_identical(maihda_lmer_zero_column_start(stuck), list())
})

test_that("the kept fit and lmer's own are refitted from a variance at zero, a refit kept only when lower by more than the tolerance", {
  env <- new.env()
  env$crit <- c(first = 10, standardised = 9, zero_standardised = 9.5, zero_first = 8)
  env$fake_fit <- function(start = "first") {
    message("fit from ", start)
    if (identical(start, "boom")) stop("this refit fails")
    list(crit = env$crit[[start]], start = start)
  }
  zero_starts <- list(standardised = "zero_standardised", first = "zero_first")
  asked <- character(0)
  local_mocked_bindings(
    maihda_lmer_refit_starts = function(model) list("standardised"),
    maihda_lmer_criterion = function(model) model$crit,
    maihda_lmer_zero_column_start = function(model) {
      asked <<- c(asked, model$start)
      if (model$start %in% names(zero_starts)) list(zero_starts[[model$start]]) else list()
    }
  )
  fit <- function() heard(maihda_fit_lmer_slopes(quote(fake_fit()), env))
  # The standardised refit is kept; from its zero the refit is higher, from lmer's own lower.
  h <- fit()
  expect_identical(asked, c("standardised", "first"))
  expect_identical(h$value$start, "zero_first")
  expect_identical(h$seen, "fit from zero_first\n")
  # Lower by less than 1e-4: the kept fit stays.
  asked <- character(0)
  env$crit[["zero_first"]] <- 9 - 5e-5
  expect_identical(fit()$value$start, "standardised")
  # lmer's own fit kept: it is refitted from its zero once.
  asked <- character(0)
  env$crit <- c(first = 10, standardised = 11, zero_first = 9)
  expect_identical(fit()$value$start, "zero_first")
  expect_identical(asked, "first")
  # No variance at zero to move, or a zero refit that fails: no change.
  zero_starts <- list()
  expect_identical(fit()$value$start, "first")
  zero_starts <- list(first = "boom")
  expect_identical(fit()$value$start, "first")
})

# ---- live fits ----------------------------------------------------------------------------

test_that("fit_maihda() reaches the optimum of fits that stop at a variance at zero or on the ridge beside it", {
  skip_slow_test()
  skip_on_cran()
  cases <- list(
    list(formula = uncorrelated, data = strata_data(11), optimum = optimum_unc_11),
    list(formula = by_context, data = strata_data(9, ctx_slope_sd = 0.05), optimum = optimum_ctx_9),
    list(formula = slope_only, data = strata_data(1, x = "year", slope_sd = 0.005, beta = 0.3),
         optimum = optimum_slope_only_1),
    list(formula = three, data = strata_data(2, x = "age", slope_sd = 0.005, beta = 0.3, n = 60, z = "year"),
         optimum = optimum_three_2))
  for (case in cases) {
    # The premise: lmer's own fit stops short (by 0.029, 12.8, 0.051 and 6.0).
    own <- quiet(lme4::lmer(case$formula, data = case$data))
    expect_gt(lme4::REMLcrit(own) - case$optimum, 0.02)
    m <- quiet(fit_maihda(case$formula, data = case$data))
    expect_lt(abs(lme4::REMLcrit(m$model) - case$optimum), 1e-4)
    expect_true(m$diagnostics$converged)
    expect_null(m$diagnostics$search_gain)
  }
})

test_that("fit_maihda() refits from lmer's own variance at zero where the kept fit has none to move", {
  skip_on_cran()
  d <- strata_data(7, ctx_slope_sd = 0.05)
  adjusted <- y ~ x + d1 + d2 + d3 + (1 + x | stratum) + (0 + x | ctx)
  own <- quiet(lme4::lmer(adjusted, data = d))
  expect_gt(lme4::REMLcrit(own) - optimum_ctx_adjusted_7, 5)
  # Here the start re-expressed block by block stops 0.015 above the optimum, where its last
  # column's variance hardly falls off zero; from lmer's own fit the zero-column refit ends within
  # 0.002 (0.00004 to 0.0013 on four row orders).
  m <- quiet(fit_maihda(adjusted, data = d))
  expect_lt(lme4::REMLcrit(m$model) - optimum_ctx_adjusted_7, 0.01)
})

# ---- the record ---------------------------------------------------------------------------

test_that("a singular lmer fit that a search reaches 0.01 below is recorded as not converged, by how much", {
  skip_on_cran()
  d <- strata_data(11)
  stuck <- placed(stuck_unc_11, d, uncorrelated)
  # lme4 calls it singular and says nothing else.
  expect_true(lme4::isSingular(stuck))
  expect_identical(as.character(stuck@optinfo$conv$lme4$messages),
                   "boundary (singular) fit: see help('isSingular')")
  gap <- lme4::REMLcrit(stuck) - optimum_unc_11
  expect_lt(abs(maihda_lmer_search_gain(stuck) / gap - 1), 1e-3)
  dg <- maihda_fit_diagnostics(stuck)
  expect_false(dg$converged)
  expect_lt(abs(dg$search_gain / gap - 1), 1e-3)
  expect_identical(dg$messages, paste0(
    "lme4 runs no convergence check on a singular fit, and a search of the same model with each ",
    "random-effect block's terms centred and scaled reaches a REML criterion 0.029 lower: the fit ",
    "stopped short of the optimum."))
  lines <- maihda_format_fit_diagnostics(dg)
  expect_true("Convergence problems of the lme4 fit:" %in% lines)
  # At the optimum, singular too, nothing is lower: converged, no field.
  at <- maihda_fit_diagnostics(placed(at_unc_11, d, uncorrelated))
  expect_true(lme4::isSingular(placed(at_unc_11, d, uncorrelated)))
  expect_true(at$converged)
  expect_null(at$search_gain)
  expect_lte(maihda_lmer_search_gain(placed(at_unc_11, d, uncorrelated)), 1e-6)
  # On the ridge: singular to lme4 with the slope off zero once re-expressed.
  ridge <- maihda_fit_diagnostics(placed(stuck_slope_only_1,
                                         strata_data(1, x = "year", slope_sd = 0.005, beta = 0.3),
                                         slope_only))
  expect_false(ridge$converged)
  expect_gt(ridge$search_gain, 0.05)
  # By maximum likelihood the fit's criterion is its deviance.
  ml <- maihda_fit_diagnostics(placed(stuck_unc_11, d, uncorrelated, REML = FALSE))
  expect_false(ml$converged)
  expect_true(grepl("reaches a deviance 0.056 lower", ml$messages, fixed = TRUE))
})

test_that("without a random slope the record's search starts from the variance at zero", {
  skip_on_cran()
  two <- y ~ 1 + (1 | stratum) + (1 | ctx)
  d <- with_seed(4, {
    st <- rep(1:18, each = 40)
    ctx <- sample(1:6, 720, replace = TRUE)
    y <- stats::rnorm(18, 0, 0.4)[st] + stats::rnorm(6, 0, 0.5)[ctx] + stats::rnorm(720)
    data.frame(y = y, stratum = factor(st), ctx = factor(ctx))
  })
  dev <- lme4::lmer(two, data = d, devFunOnly = TRUE)
  # The context's variance held at zero and the strata's at its best given that: a saddle, far
  # above the optimum.
  face <- stats::optimize(function(t) dev(c(t, 0)), c(0, 5))$minimum
  gap <- dev(c(face, 0)) - stats::nlminb(c(face, 0.3), dev)$objective
  expect_gt(gap, 100)
  m <- placed(c(face, 0), d, two)
  # Intercepts only: nothing to re-express, so no search from lme4's starting values.
  expect_identical(maihda_lmer_standardised_start(m), list())
  dg <- maihda_fit_diagnostics(m)
  expect_false(dg$converged)
  expect_lt(abs(dg$search_gain / gap - 1), 1e-6)
})

test_that("the record's line is 0.01, and a fit not at the maximum is not said to leave the VPC unaffected", {
  skip_on_cran()
  stuck <- placed(stuck_unc_11, strata_data(11), uncorrelated)
  at_gain <- function(gain) {
    local_mocked_bindings(maihda_lmer_search_gain = function(model) gain)
    maihda_fit_diagnostics(stuck)
  }
  expect_false(at_gain(0.01)$converged)
  expect_identical(at_gain(0.01)$search_gain, 0.01)
  below <- at_gain(0.01 - 1e-7)
  expect_true(below$converged)
  expect_null(below$search_gain)
  expect_true(at_gain(NA_real_)$converged)
  # A singular block that is not the stratum's: the stratum variance is said to be unaffected only
  # where nothing found the fit short of the maximum.
  base <- list(engine = "lme4", singular = TRUE, singular_terms = "ctx", stratum_singular = FALSE,
               converged = TRUE, messages = character(0))
  unaffected <- "  The 'stratum' block is NOT at the boundary, so the between-stratum variance and"
  expect_true(unaffected %in% maihda_format_fit_diagnostics(base))
  short <- utils::modifyList(base, list(converged = FALSE, search_gain = 0.2, messages = "short"))
  lines <- maihda_format_fit_diagnostics(short)
  expect_false(unaffected %in% lines)
  expect_identical(lines[3:4], c("Convergence problems of the lme4 fit:", "  - short"))
})

test_that("only a singular lmer fit outside growth models, not found short already and not refitted by fit_maihda(), is searched", {
  skip_on_cran()
  searched <- 0L
  local_mocked_bindings(maihda_lmer_search_gain = function(model) {
    searched <<- searched + 1L
    0.5
  })
  stuck <- placed(stuck_unc_11, strata_data(11), uncorrelated)
  expect_false(maihda_fit_diagnostics(stuck)$converged)
  expect_identical(searched, 1L)
  # Refitted by fit_maihda() from the same searches, or a growth model: not searched.
  expect_true(maihda_fit_diagnostics(stuck, searched = TRUE)$converged)
  expect_true(maihda_fit_diagnostics(stuck, longitudinal_info = list(id = "id"))$converged)
  # Not singular, which lme4 checks itself.
  sleep <- quiet(lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy))
  expect_false(lme4::isSingular(sleep))
  expect_true(maihda_fit_diagnostics(sleep)$converged)
  # A glmer fit.
  d <- strata_data(11)
  d$yb <- as.integer(d$y > stats::median(d$y))
  gl <- quiet(lme4::glmer(yb ~ x + (1 | stratum), data = d, family = stats::binomial))
  maihda_fit_diagnostics(gl)
  # Already found short at a boundary trap.
  local_mocked_bindings(maihda_lme4_boundary_traps = function(model, ...) {
    data.frame(block = "stratum", term = "(Intercept)", index = 1L, drop = 1, n_terms = 1L,
               stringsAsFactors = FALSE)
  })
  expect_false(maihda_fit_diagnostics(stuck)$converged)
  expect_identical(searched, 1L)
})

test_that("fit_maihda() searches a fit made from the user's own start and none it refitted", {
  skip_on_cran()
  d <- strata_data(11)
  # lmer started at the saddle stays there.
  m <- quiet(fit_maihda(uncorrelated, data = d, start = stuck_unc_11))
  expect_gt(lme4::REMLcrit(m$model) - optimum_unc_11, 0.02)
  expect_false(m$diagnostics$converged)
  expect_gt(m$diagnostics$search_gain, 0.02)
  calls <- 0L
  real <- maihda_lmer_search_gain
  local_mocked_bindings(maihda_lmer_search_gain = function(model) {
    calls <<- calls + 1L
    real(model)
  })
  m2 <- quiet(fit_maihda(uncorrelated, data = d))
  expect_identical(calls, 0L)
  expect_true(lme4::isSingular(m2$model))
  expect_true(m2$diagnostics$converged)
})
