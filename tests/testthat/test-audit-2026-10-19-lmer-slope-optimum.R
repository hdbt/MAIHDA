# Audit pass 2026-10-19: lmer fits with a random slope, outside growth models, that stop short
# of the REML optimum.
#
# A random-effect block with a slope, such as (1 + x | stratum), fitted by lmer can stop with
# its intercept Cholesky entry at its bound of 0 beside a non-zero intercept-slope entry, as
# the glmer fits of the same packet's predecessor did; and with a slope covariate far from zero
# or in large units (age in years, a birth year) lmer stops short of the optimum without the
# intercept at its bound. Neither was repaired nor recorded. Now such a fit is refitted from
# across the bound and from a start carried back from its random-effect terms centred and
# scaled, the lowest criterion kept, and the convergence record flags a cross-sectional lmer
# fit still at such a bound.
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

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 0.3) of 40 rows, a stratum slope
# deviation of SD 0.15 per unit of x ~ N(0, 1), or of SD 0.005 per year of age (18 to 80,
# deviations centred at 49), and a residual SD of 1: the adjusted model's intercept variance
# is small next to its slope variance.
slope_data <- function(seed, age = FALSE) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    u <- stats::rnorm(18, 0, 0)
    v <- stats::rnorm(18, 0, if (age) 0.005 else 0.15)
    st <- rep(seq_len(18), each = 40)
    x <- if (age) sample(18:80, 720, replace = TRUE) else stats::rnorm(720)
    y <- 1 + add[st] + u[st] + (0.3 + v[st]) * (x - if (age) 49 else 0) + stats::rnorm(720)
    data.frame(y = y, x = x, d1 = factor(lev$d1[st]), d2 = factor(lev$d2[st]),
               d3 = factor(lev$d3[st]), stratum = factor(st))
  })
}

adjusted <- y ~ x + d1 + d2 + d3 + (1 + x | stratum)

# Seed 3: lmer's own fit stops with the intercept entry at 0 (theta 0, 0.09122, 0.1524) at a
# REML criterion of 2041.02859; started there it stays. The optimum, from lmer's deviance
# function searched with no bound on the Cholesky diagonal, is 2040.83268 (a rank-1 block,
# correlation -1). Fitted by maximum likelihood: 2016.31505 at (0, 0.02171, 0.1644); the
# optimum is 2016.18939.
trap_3 <- c(0, 0.09122, 0.1524)
optimum_3 <- 2040.83268
trap_3_ml <- c(0, 0.02171, 0.1644)
optimum_3_ml <- 2016.18939
# Seed 7 with age in years: lmer's own fit stops 0.371 short of the optimum, 2145.58274, with
# the intercept entry off its bound (theta 0.2900, -0.005057, 0: the block at a correlation of
# -1), in each of ten row orders. The optimum is the same with age written as I(x - 18).
optimum_7_age <- 2145.58274

# ---- which fits ---------------------------------------------------------------------------

test_that("lmer fits outside growth models with a random slope are refitted, unless the user gives a start", {
  lmer <- quote(lme4::lmer)
  slope <- y ~ x + (1 + x | stratum)
  expect_true(maihda_lmer_slopes_apply(lmer, NULL, list(), slope))
  expect_true(maihda_lmer_slopes_apply(lmer, NULL, list(REML = FALSE, control = 1), slope))
  expect_true(maihda_lmer_slopes_apply(lmer, NULL, list(), y ~ x + (0 + x | stratum)))
  expect_true(maihda_lmer_slopes_apply(lmer, NULL, list(), y ~ x + (x || stratum)))
  expect_true(maihda_lmer_slopes_apply(lmer, NULL, list(), y ~ x + (1 | stratum) + (x | school)))
  # Intercept-only random effects take lmer's own fit alone.
  expect_false(maihda_lmer_slopes_apply(lmer, NULL, list(), y ~ x + (1 | stratum)))
  expect_false(maihda_lmer_slopes_apply(lmer, NULL, list(), y ~ x + (1 | stratum) + (1 | school)))
  # A start of the user's own, under any spelling lmer() binds to `start`.
  expect_false(maihda_lmer_slopes_apply(lmer, NULL, list(start = 1), slope))
  expect_false(maihda_lmer_slopes_apply(lmer, NULL, list(st = 1), slope))
  # Not glmer, which has its own refit, and not a growth model, which has its own restarts.
  expect_false(maihda_lmer_slopes_apply(quote(lme4::glmer), NULL, list(), slope))
  expect_false(maihda_lmer_slopes_apply(quote(lme4::glmer.nb), NULL, list(), slope))
  lng <- list(id = "id", time = "wave", time_degree = 1L, stratum_slope = FALSE)
  expect_false(maihda_lmer_slopes_apply(lmer, lng, list(), slope))
})

test_that("fit_maihda() sends lmer fits with a random slope through the refit, and no other fit", {
  skip_on_cran()
  calls <- 0L
  real <- maihda_fit_lmer_slopes
  local_mocked_bindings(maihda_fit_lmer_slopes = function(...) { calls <<- calls + 1L; real(...) })
  d <- slope_data(3)
  quiet(fit_maihda(adjusted, data = d))
  quiet(fit_maihda(adjusted, data = d, REML = FALSE))
  expect_identical(calls, 2L)
  quiet(fit_maihda(y ~ x + d1 + d2 + d3 + (1 | stratum), data = d))
  quiet(fit_maihda(adjusted, data = d, start = trap_3))
  quiet(fit_maihda(round(exp(y / 4)) ~ x + (1 + x | stratum), data = d, family = "poisson"))
  data(maihda_long_data, package = "MAIHDA", envir = environment())
  ld <- maihda_long_data[maihda_long_data$id %in% unique(maihda_long_data$id)[1:60], ]
  quiet(fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = ld, id = "id",
                   time = "wave", stratum_slope = FALSE))
  expect_identical(calls, 2L)
})

# ---- the standardised spelling ------------------------------------------------------------

test_that("random-effect terms are scaled, and centred where every block holding them has an intercept", {
  x <- c(1, 2, 4, 8, 16, 32)
  w <- c(10, 20, 30, 40, 50, 70)
  fb <- c(0, 1, 0, 1, 0, 1)
  values <- list(stratum = cbind(`(Intercept)` = 1, x = x), school = cbind(w = w),
                 other = cbind(`(Intercept)` = 1, fb = fb, k = 5))
  cnms <- lapply(values, colnames)
  s <- maihda_lmer_term_scales(values, cnms)
  # w sits in a block without an intercept, so it is only scaled; a constant term is left out.
  expect_equal(s$centre, c(x = mean(x), w = 0, fb = 0.5))
  expect_equal(s$scale, c(x = stats::sd(x), w = stats::sd(w), fb = stats::sd(fb)))
  # A term in a block without an intercept as well is only scaled.
  both <- list(a = cbind(`(Intercept)` = 1, x = x), g = cbind(x = x))
  expect_identical(maihda_lmer_term_scales(both, lapply(both, colnames))$centre, c(x = 0))
  # Nothing to standardise, or values that do not match the blocks.
  expect_null(maihda_lmer_term_scales(list(cbind(`(Intercept)` = rep(1, 4))), list(stratum = "(Intercept)")))
  expect_null(maihda_lmer_term_scales(list(cbind(`(Intercept)` = 1, z = 1:4)),
                                      list(stratum = c("(Intercept)", "x"))))
  expect_null(maihda_lmer_term_scales(values[1:2], cnms))
  # Missing values do not enter the mean or the SD.
  values$stratum[2, "x"] <- NA
  s3 <- maihda_lmer_term_scales(values, cnms)
  expect_equal(s3$centre[["x"]], mean(x[-2]))
  expect_equal(s3$scale[["x"]], stats::sd(x[-2]))
})

test_that("the theta carried back gives every observation the random-effect variance it had", {
  # A block (Intercept), x, w fitted on z = (x - 49) / 18 and v = (w - 2) / 3, and an
  # intercept-only block, which is unchanged.
  L <- matrix(0, 3, 3)
  L[lower.tri(L, diag = TRUE)] <- c(0.5, 0.2, -0.1, 0.3, 0.05, 0.2)
  cnms <- list(stratum = c("(Intercept)", "x", "w"), school = "(Intercept)")
  theta_z <- c(L[lower.tri(L, diag = TRUE)], 0.7)
  centre <- c(x = 49, w = 2)
  scale <- c(x = 18, w = 3)
  back <- maihda_lmer_unstandardise_theta(theta_z, cnms, centre, scale)
  Lx <- matrix(0, 3, 3)
  Lx[lower.tri(Lx, diag = TRUE)] <- back[1:6]
  expect_true(all(diag(Lx) >= 0))
  for (p in list(c(18, -1), c(49, 2), c(80, 5), c(0, 0))) {
    xrow <- c(1, p)
    zrow <- c(1, (p - centre) / scale)
    expect_equal(drop(xrow %*% tcrossprod(Lx) %*% xrow), drop(zrow %*% tcrossprod(L) %*% zrow))
  }
  expect_identical(back[7], 0.7)
  # A rank-1 block (correlation -1) stays rank 1.
  r1 <- maihda_lmer_unstandardise_theta(c(0.4, -0.2, 0), list(stratum = c("(Intercept)", "x")),
                                        c(x = 49), c(x = 18))
  expect_lt(abs(r1[3]), 1e-8)
  expect_equal(drop(c(1, 80) %*% tcrossprod(matrix(c(r1[1], r1[2], 0, r1[3]), 2)) %*% c(1, 80)),
               (0.4 - 0.2 * (80 - 49) / 18)^2)
  # A block without an intercept is only rescaled; a theta of another layout gives NULL.
  expect_equal(maihda_lmer_unstandardise_theta(0.6, list(g = "x"), c(x = 0), c(x = 2)), 0.3)
  expect_null(maihda_lmer_unstandardise_theta(c(0.1, 0.2), cnms, centre, scale))
})

# ---- the refit, with stand-in fits --------------------------------------------------------

test_that("a refit replaces lmer's fit only when lower by more than the tolerance, and the kept fit speaks alone", {
  env <- new.env()
  env$crit <- c(first = 10, across = 9, standardised = 9.5)
  env$fake_fit <- function(start = "first") {
    message("fit from ", start)
    if (identical(start, "boom")) stop("this refit fails")
    list(crit = env$crit[[start]], start = start)
  }
  starts <- list("across", "standardised")
  local_mocked_bindings(
    maihda_lmer_refit_starts = function(model) starts,
    maihda_lmer_criterion = function(model) model$crit
  )
  fit <- function(...) heard(maihda_fit_lmer_slopes(quote(fake_fit()), env, ...))
  h <- fit()
  expect_identical(h$value$start, "across")
  expect_identical(h$seen, "fit from across\n")
  env$crit <- c(first = 10, across = 9.6, standardised = 9.5)
  expect_identical(fit()$value$start, "standardised")
  # Lower by less than 1e-4: lmer's own fit stays; by more: the refit replaces it.
  env$crit <- c(first = 10, across = 10 - 5e-5, standardised = 10 - 2e-5)
  h <- fit()
  expect_identical(h$value$start, "first")
  expect_identical(h$seen, "fit from first\n")
  env$crit <- c(first = 10, across = 10 - 2e-4, standardised = 11)
  expect_identical(fit()$value$start, "across")
  expect_identical(fit(tol = 0.5)$value$start, "first")
  expect_identical(formals(maihda_fit_lmer_slopes)$tol, 1e-4)
  # A refit that fails is passed over, and so are starts that cannot be made.
  env$crit <- c(first = 10, across = 9, standardised = 9.5)
  starts <- list("boom", "standardised")
  expect_identical(fit()$value$start, "standardised")
  local_mocked_bindings(maihda_lmer_refit_starts = function(model) stop("no starts"))
  expect_identical(fit()$value$start, "first")
  # A first fit that fails fails.
  env$fake_fit <- function(start = "first") stop("the first fit fails")
  expect_error(maihda_fit_lmer_slopes(quote(fake_fit()), env), "the first fit fails")
})

test_that("the starts are the fit across its bound and the fit with its terms standardised", {
  skip_on_cran()
  d <- slope_data(3)
  first <- quiet(lme4::lmer(adjusted, data = d, start = trap_3))
  th <- unname(lme4::getME(first, "theta"))
  expect_lt(th[1], 1e-4)
  starts <- maihda_lmer_refit_starts(first)
  expect_length(starts, 2L)
  expect_equal(starts[[1]], c(0.05, -th[2], th[3]))
  # The standardised fit's start, carried back, is the optimum already.
  dev <- lme4::lmer(adjusted, data = d, devFunOnly = TRUE)
  expect_lt(abs(dev(starts[[2]]) - optimum_3), 0.005)
  # A fit at an interior optimum is not taken across anything, and its standardised start
  # is the same point: the spellings are one model.
  inner <- slope_data(1)
  at <- quiet(lme4::lmer(y ~ x + (1 + x | stratum), data = inner))
  s <- maihda_lmer_refit_starts(at)
  expect_length(s, 1L)
  dev_null <- lme4::lmer(y ~ x + (1 + x | stratum), data = inner, devFunOnly = TRUE)
  expect_lt(abs(dev_null(s[[1]]) - lme4::REMLcrit(at)), 1e-5)
  # A slope term written as an expression is standardised as it enters the model, on the
  # fitted rows, with the fixed part as the formula has it.
  age <- slope_data(7, age = TRUE)
  spelled <- y ~ I(x - 18) + d1 + d2 + d3 + (1 + I(x - 18) | stratum)
  own <- quiet(lme4::lmer(spelled, data = age))
  st <- maihda_lmer_standardised_start(own)
  dev_age <- lme4::lmer(spelled, data = age, devFunOnly = TRUE)
  expect_lt(abs(dev_age(st[[1]]) - optimum_7_age), 0.005)
  # With a second grouping beside the strata, every block is standardised.
  d$school <- factor(rep(1:12, length.out = nrow(d)))
  two <- quiet(lme4::lmer(y ~ x + d1 + d2 + d3 + (1 + x | stratum) + (1 | school), data = d))
  st2 <- maihda_lmer_standardised_start(two)
  expect_length(st2, 1L)
  dev2 <- lme4::lmer(y ~ x + d1 + d2 + d3 + (1 + x | stratum) + (1 | school), data = d,
                     devFunOnly = TRUE)
  expect_lte(dev2(st2[[1]]), lme4::REMLcrit(two) + 1e-6)
  # No term to standardise: no start. A deviance function that does not give the fit's own
  # criterion: no start either.
  expect_identical(maihda_lmer_standardised_start(quiet(lme4::lmer(y ~ x + (1 | stratum), data = d))),
                   list())
  local_mocked_bindings(maihda_lmer_criterion = function(model) lme4::REMLcrit(model) + 1)
  expect_identical(maihda_lmer_standardised_start(at), list())
})

test_that("the standardised minimum is taken across its bound only when that is lower", {
  skip_on_cran()
  first <- quiet(lme4::lmer(adjusted, data = slope_data(3)))
  cnms <- lme4::getME(first, "cnms")
  std <- maihda_lmer_term_scales(lme4::getME(first, "mmList"), cnms)
  carry <- function(th) maihda_lmer_unstandardise_theta(th, cnms, std$centre, std$scale)
  # A stand-in optimiser: its first minimum stops at the bound, its second, from across the
  # bound, is `across`.
  calls <- list()
  across_value <- 9
  stand_in <- function(par, fn, lower, upper) {
    calls[[length(calls) + 1L]] <<- list(par = par, lower = lower, upper = upper)
    if (length(calls) == 1L) list(par = c(0, 0.2, 0.1), fval = 10) else
      list(par = c(0.3, -0.2, 0.1), fval = across_value)
  }
  expect_equal(maihda_lmer_standardised_start(first, optimiser = stand_in), list(carry(c(0.3, -0.2, 0.1))))
  # lme4's starting values and bounds, then the first minimum taken across its bound.
  expect_equal(calls[[1]]$par, c(1, 0, 1))
  expect_equal(calls[[1]]$lower, unname(lme4::getME(first, "lower")))
  expect_equal(calls[[2]]$par, c(0.05, -0.2, 0.1))
  # Across the bound is higher: the first minimum is carried back.
  calls <- list()
  across_value <- 11
  expect_equal(maihda_lmer_standardised_start(first, optimiser = stand_in), list(carry(c(0, 0.2, 0.1))))
  # A first minimum off the bound is not taken across.
  calls <- list()
  off <- function(par, fn, lower, upper) {
    calls[[length(calls) + 1L]] <<- par
    list(par = c(0.3, 0.2, 0.1), fval = 10)
  }
  expect_equal(maihda_lmer_standardised_start(first, optimiser = off), list(carry(c(0.3, 0.2, 0.1))))
  expect_length(calls, 1L)
  # The criterion it is handed is the fit's own deviance at the carried-back theta.
  probe <- NULL
  maihda_lmer_standardised_start(first, optimiser = function(par, fn, lower, upper) {
    probe <<- fn(c(0.3, -0.2, 0.1))
    list(par = par, fval = fn(par))
  })
  dev <- lme4::lmer(adjusted, data = slope_data(3), devFunOnly = TRUE)
  expect_equal(probe, dev(carry(c(0.3, -0.2, 0.1))), tolerance = 1e-8)
})

# ---- live fits ----------------------------------------------------------------------------

test_that("an lmer fit stopped at the bound is refitted to the REML optimum", {
  skip_on_cran()
  env <- new.env()
  env$data <- slope_data(3)
  trapped_call <- quote(lme4::lmer(y ~ x + d1 + d2 + d3 + (1 + x | stratum), data = data))
  trapped_call$start <- trap_3
  m <- quiet(maihda_fit_lmer_slopes(trapped_call, env))
  expect_lt(abs(lme4::REMLcrit(m) - optimum_3), 0.005)
  expect_gt(lme4::getME(m, "theta")[[1]], 0.01)
  # Through fit_maihda() from lmer's own start, by REML and by maximum likelihood.
  fm <- quiet(fit_maihda(adjusted, data = env$data))
  expect_lt(abs(lme4::REMLcrit(fm$model) - optimum_3), 0.005)
  expect_true(fm$diagnostics$converged)
  expect_null(fm$diagnostics$boundary_traps)
  ml <- quiet(fit_maihda(adjusted, data = env$data, REML = FALSE))
  expect_false(lme4::isREML(ml$model))
  expect_lt(abs(stats::deviance(ml$model) - optimum_3_ml), 0.005)
})

test_that("an lmer fit with a slope on age in years reaches the REML optimum", {
  skip_on_cran()
  d <- slope_data(7, age = TRUE)
  m <- quiet(fit_maihda(adjusted, data = d))
  expect_lt(abs(lme4::REMLcrit(m$model) - optimum_7_age), 0.005)
  # It is stored as a fit of the formula on age in years, with the data as given.
  expect_identical(unname(lme4::getME(m$model, "X")[, "x"]), as.numeric(d$x))
  expect_identical(m$data$x, d$x)
  # Written as an expression, and beside a transformed fixed term, on all the rows.
  spelled <- quiet(fit_maihda(y ~ I(x - 18) + d1 + d2 + d3 + (1 + I(x - 18) | stratum), data = d))
  expect_lt(abs(lme4::REMLcrit(spelled$model) - optimum_7_age), 0.005)
  logged <- quiet(fit_maihda(y ~ log(x) + d1 + d2 + d3 + (1 + x | stratum), data = d))
  own <- quiet(lme4::lmer(y ~ log(x) + d1 + d2 + d3 + (1 + x | stratum), data = d))
  expect_identical(stats::nobs(logged$model), nrow(d))
  expect_lte(lme4::REMLcrit(logged$model), lme4::REMLcrit(own) + 1e-8)
})

test_that("the refits keep weights, an offset and a subset given as arguments", {
  skip_on_cran()
  d <- slope_data(3)
  d$w <- rep(c(1, 2), length.out = nrow(d))
  d$off <- round(sin(seq_len(nrow(d))), 3)
  m <- quiet(fit_maihda(adjusted, data = d, weights = w, offset = off, subset = x > -2))
  keep <- d$x > -2
  expect_identical(stats::nobs(m$model), sum(keep))
  expect_identical(unname(stats::weights(m$model)), d$w[keep])
  expect_identical(unname(lme4::getME(m$model, "offset")), d$off[keep])
  # The kept fit is the one lmer makes with the same arguments from its estimates, and it is
  # no worse than lmer's own.
  direct <- quiet(lme4::lmer(adjusted, data = d, weights = w, offset = off, subset = x > -2,
                             start = unname(lme4::getME(m$model, "theta"))))
  expect_lt(abs(lme4::REMLcrit(m$model) - lme4::REMLcrit(direct)), 1e-6)
  own <- quiet(lme4::lmer(adjusted, data = d, weights = w, offset = off, subset = x > -2))
  expect_lte(lme4::REMLcrit(m$model), lme4::REMLcrit(own) + 1e-8)
})

test_that("an intercept-only lmer fit is made as before, and the refits draw no random numbers", {
  skip_on_cran()
  d <- slope_data(3)
  held <- 0L
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    held <<- held + 1L
    real(call, env)
  })
  m <- quiet(fit_maihda(y ~ x + d1 + d2 + d3 + (1 | stratum), data = d))
  expect_identical(held, 0L)
  direct <- quiet(lme4::lmer(y ~ x + d1 + d2 + d3 + (1 | stratum), data = d))
  expect_identical(lme4::getME(m$model, "theta"), lme4::getME(direct, "theta"))
  expect_identical(lme4::fixef(m$model), lme4::fixef(direct))
  # lme4's own message still reaches the caller: strata that hold the same values have no
  # variance between them.
  same <- data.frame(y = rep(slope_data(1)$y[1:30], times = 6), stratum = factor(rep(1:6, each = 30)))
  h <- heard(fit_maihda(y ~ 1 + (1 | stratum), data = same))
  expect_true(any(grepl("boundary (singular) fit", h$seen, fixed = TRUE)))
  with_seed(1, {
    before <- get(".Random.seed", envir = globalenv())
    quiet(fit_maihda(adjusted, data = d))
    expect_identical(get(".Random.seed", envir = globalenv()), before)
  })
})

# ---- the record ---------------------------------------------------------------------------

test_that("the record flags a cross-sectional lmer fit left at such a bound, not an lmer growth fit", {
  skip_on_cran()
  d <- slope_data(3)
  # The user's own start is fitted as given: lmer's trap.
  m <- quiet(fit_maihda(adjusted, data = d, start = trap_3))
  expect_gt(lme4::REMLcrit(m$model), optimum_3 + 0.15)
  dg <- m$diagnostics
  expect_false(dg$converged)
  expect_identical(dg$boundary_traps$block, "stratum")
  expect_identical(dg$boundary_traps$term, "(Intercept)")
  sentence <- paste0("The intercept variance of the 'stratum' block is at zero, where the ",
                     "likelihood still rises with the sign of its covariance with x reversed: ",
                     "the fit stopped at a bound of lme4's parametrisation, not at the maximum.")
  expect_true(sentence %in% dg$messages)
  printed <- utils::capture.output(print(m))
  expect_true(any(grepl("Convergence problems of the lme4 fit:", printed, fixed = TRUE)))
  # The drop is the REML criterion's, a step of 1e-3 across the bound.
  th <- unname(lme4::getME(m$model, "theta"))
  dev <- lme4::lmer(adjusted, data = d, devFunOnly = TRUE)
  expect_equal(dg$boundary_traps$drop, dev(th) - dev(c(1e-3, -th[2], th[3])), tolerance = 1e-8)
  expect_gt(dg$boundary_traps$drop, .maihda_lme4_trap_drop)
  # Measuring leaves the fit where it is.
  before <- lme4::getME(m$model, "theta")
  maihda_lme4_boundary_traps(m$model)
  expect_identical(lme4::getME(m$model, "theta"), before)
  # Fitted by maximum likelihood, on the deviance.
  ml <- quiet(fit_maihda(adjusted, data = d, REML = FALSE, start = trap_3_ml))
  expect_false(ml$diagnostics$converged)
  expect_identical(ml$diagnostics$boundary_traps$term, "(Intercept)")
  # An lmer growth fit is not read by the record.
  expect_null(maihda_fit_diagnostics(m$model, list(id = "id", time = "wave"))$boundary_traps)
  # At the optimum there is nothing to flag.
  ok <- quiet(fit_maihda(adjusted, data = d))
  expect_identical(nrow(maihda_lme4_boundary_traps(ok$model)), 0L)
  expect_true(ok$diagnostics$converged)
})
