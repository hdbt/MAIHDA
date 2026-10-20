# Audit pass 2026-10-20: lmer fits with a random slope on a covariate far from zero or in large
# units recorded as not converged at the REML optimum.
#
# lme4 takes its convergence checks of an lmer fit in theta as the formula has it, by central
# differences with a step of 1e-4 and absolute lines (0.002 on a scaled gradient, 1e6 on an
# eigenvalue of the Hessian). For a slope on a birth year or an age in years a step along a
# slope entry moves each observation's random-effect variance by an amount that grows with the
# square of its covariate, so the checks fire at the optimum itself: on simulated designs
# (18 strata of 40 rows) 17 of 60 birth-year, 47 of 180 age-in-years and 24 of 60
# income-in-thousands fits returned by fit_maihda() were recorded as not converged, every one
# within 0.001 of the REML optimum found by an unbounded search of lmer's own deviance
# function. Per the user, when lme4's checks fire on an lmer fit whose optimiser reports
# success, the fit is measured as a glmer fit is, with each random-effect block's terms centred
# and scaled (in a block without an intercept, divided by their root mean square): within 0.1
# standard errors of a stationary point (sqrt(g' H^-1 g / 2) of the REML criterion) it counts
# as converged, lme4's messages are kept in `lme4_overruled`, and print() says what lme4 warned
# of and how far the fit is. This holds for every lmer fit, growth fits included; a fit lme4
# does not check is not measured.
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

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 0.3) of 40 rows, each row with a
# birth year between 1940 and 2000, a stratum slope deviation of SD 0.005 per year (deviations
# centred at 1970) and a residual SD of 1. With `slope_sd = 0` and the year written as
# (x - 1970) / 1000 the strata share one slope.
year_data <- function(seed, slope_sd = 0.005, small = FALSE) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    if (slope_sd > 0) {
      u <- stats::rnorm(18, 0, 0)
      v <- stats::rnorm(18, 0, slope_sd)
    } else {
      u <- v <- numeric(18)
    }
    st <- rep(seq_len(18), each = 40)
    x <- sample(1940:2000, 720, replace = TRUE)
    y <- 1 + add[st] + u[st] + (0.3 + v[st]) * (x - 1970) + stats::rnorm(720)
    data.frame(y = y, x = if (small) (x - 1970) / 1000 else x, stratum = factor(st))
  })
}

# Eighteen strata of 40 rows with a calendar year between 2010 and 2020, slope deviations of
# SD 0.05 per year by stratum and by one of six contexts, and a residual SD of 1.
calendar_data <- function(seed) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    st <- rep(seq_len(18), each = 40)
    x <- sample(2010:2020, 720, replace = TRUE)
    v <- stats::rnorm(18, 0, 0.05)
    ctx <- sample(1:6, 720, replace = TRUE)
    w <- stats::rnorm(6, 0, 0.05)
    y <- 1 + add[st] + (0.1 + v[st] + w[ctx]) * (x - 2015) + stats::rnorm(720)
    data.frame(y = y, x = x, stratum = factor(st), ctx = factor(ctx))
  })
}

slope <- y ~ x + (1 + x | stratum)

# An lmer fit placed at `theta`: an optimiser that returns its start, so lme4's checks are
# taken at that point exactly.
stay <- function(fn, par, lower, upper, control) {
  list(par = par, fval = fn(par), conv = 0, feval = 1, message = "stay")
}
placed <- function(theta, data = year_data(3), formula = slope, REML = TRUE) {
  quiet(lme4::lmer(formula, data = data, REML = REML, start = theta,
                   control = lme4::lmerControl(optimizer = stay)))
}

# Seed 3: the REML optimum, from lmer's deviance function searched with no bound on the
# Cholesky diagonal, is 1994.85420 at theta (18.347146, -0.0093008294, 0.00017077300); there
# lme4 reports max|grad| = 1.28 and a very large eigenvalue (2e9). Points along one direction
# of the standardised theta lie 0.0081, 0.0121 and 0.03 above it.
optimum_3 <- c(18.347146, -0.0093008294, 0.00017077300)
above_0.0081 <- c(18.674103, -0.0094701332, 0.00017264300)
above_0.0121 <- c(18.748465, -0.0095086216, 0.00017305280)
above_0.03 <- c(18.987570, -0.0096323374, 0.00017433472)
# Fitted by maximum likelihood, the optimum is a deviance of 1981.71431 at theta (17.377846,
# -0.0088087318, 0.00016502863), where the REML criterion is 1994.907; lme4 fires there too.
optimum_3_ml <- c(17.377846, -0.0088087318, 0.00016502863)

# A block (Intercept), x carried to the same block on the covariate T applied to it.
carry_block <- function(theta, Tm) {
  L <- matrix(c(theta[1], theta[2], 0, theta[3]), 2)
  L2 <- maihda_chol_psd(Tm %*% tcrossprod(L) %*% t(Tm))
  L2[lower.tri(L2, diag = TRUE)]
}

# ---- the distance --------------------------------------------------------------------------

test_that("the distance is sqrt(g' H^-1 g / 2), and NA without a positive-definite Hessian", {
  d <- list(gradient = c(0.3, -0.4), hessian = diag(c(2, 0.5)))
  expect_equal(maihda_newton_distance(d), sqrt((0.3^2 / 2 + 0.4^2 / 0.5) / 2))
  expect_identical(maihda_newton_distance(list(gradient = c(0, 0), hessian = diag(2))), 0)
  # Curved down along the second parameter, with a gradient along the first only.
  expect_true(is.na(maihda_newton_distance(list(gradient = c(1, 0), hessian = diag(c(2, -1))))))
  expect_true(is.na(maihda_newton_distance(list(gradient = c(1, 0), hessian = diag(c(2, 0))))))
  expect_true(is.na(maihda_newton_distance(list(gradient = c(NA, 0), hessian = diag(2)))))
  expect_true(is.na(maihda_newton_distance(list(gradient = c(1, 0), hessian = matrix(Inf, 2, 2)))))
  expect_true(is.na(maihda_newton_distance(NULL)))
})

test_that("an lmer fit is measured in its standardised terms, the same whatever the covariate's spelling", {
  skip_on_cran()
  d <- year_data(3)
  at <- placed(optimum_3, d)
  expect_lt(maihda_lmer_stationary_distance(at), 1e-4)
  # The year as I(x - 1900) moves the intercept to 1900; in decades the slope is per ten years.
  d$dec <- d$x / 10
  in_1900 <- placed(carry_block(optimum_3, matrix(c(1, 0, 1900, 1), 2)), d,
                    y ~ I(x - 1900) + (1 + I(x - 1900) | stratum))
  in_decades <- placed(carry_block(optimum_3, diag(c(1, 10))), d, y ~ dec + (1 + dec | stratum))
  expect_lt(abs(lme4::REMLcrit(in_1900) - lme4::REMLcrit(at)), 1e-6)
  distance <- vapply(list(at, in_1900, in_decades), maihda_lmer_stationary_distance, numeric(1))
  expect_true(all(distance < 1e-4))
  # Where the fit is 0.03 above the optimum, the three spellings agree on how far it is.
  far <- vapply(list(placed(above_0.03, d),
                     placed(carry_block(above_0.03, matrix(c(1, 0, 1900, 1), 2)), d,
                            y ~ I(x - 1900) + (1 + I(x - 1900) | stratum)),
                     placed(carry_block(above_0.03, diag(c(1, 10))), d, y ~ dec + (1 + dec | stratum))),
                maihda_lmer_stationary_distance, numeric(1))
  expect_lt(max(abs(far / far[1] - 1)), 1e-4)
  # lme4's own derivatives, in theta as each formula has it, say something else on each.
  raw <- function(m) {
    dev <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
    th <- unname(lme4::getME(m, "theta"))
    maihda_newton_distance(maihda_central_derivatives(dev, th, fx = dev(th)))
  }
  expect_false(isTRUE(raw(at) < 1e-4))
})

test_that("the distance is the criterion's shortfall where the fit is short of the optimum", {
  skip_on_cran()
  d <- year_data(3)
  opt <- lme4::REMLcrit(placed(optimum_3, d))
  for (case in list(list(theta = above_0.0081, gap = 0.0081), list(theta = above_0.0121, gap = 0.0121),
                    list(theta = above_0.03, gap = 0.03))) {
    m <- placed(case$theta, d)
    gap <- lme4::REMLcrit(m) - opt
    expect_lt(abs(gap / case$gap - 1), 1e-4)
    # A ratio, not expect_equal(tolerance =): under this suite's testthat edition a scalar
    # passes that whenever the absolute difference is below the tolerance.
    expect_lt(abs(maihda_lmer_stationary_distance(m) / sqrt(gap) - 1), 0.05)
  }
})

test_that("the differences may cross a variance bound", {
  skip_on_cran()
  # The strata share one slope, so the optimum is a rank-1 block: what is left of the slope's
  # variance after the intercept is zero. At a point with that standardised diagonal at 5e-5,
  # within a step of its bound, lme4's own entry is 0.0029 (the year in thousands of years), so
  # lme4 checks such a fit.
  d <- year_data(1, slope_sd = 0, small = TRUE)
  near <- placed(c(0.35650, 1.2490, 0.0029270), d)
  expect_false(lme4::isSingular(near))
  expect_lt(maihda_lmer_stationary_distance(near), 0.01)
})

test_that("each block's terms are centred and scaled in that block, by its root mean square without an intercept", {
  x <- c(2011, 2013, 2014, 2016, 2019, 2020)
  values <- list(stratum = cbind(`(Intercept)` = 1, x = x), stratum.1 = cbind(x = x),
                 other = cbind(`(Intercept)` = 1, k = rep(5, 6)))
  cnms <- lapply(values, colnames)
  Ts <- maihda_lmer_block_scales(values, cnms)
  expect_equal(Ts[[1]], matrix(c(1, 0, -mean(x) / stats::sd(x), 1 / stats::sd(x)), 2))
  expect_equal(Ts[[2]], matrix(1 / sqrt(mean(x^2)), 1))
  # A constant term is left as it is.
  expect_identical(Ts[[3]], diag(1, 2))
  # Values that do not match the blocks.
  expect_null(maihda_lmer_block_scales(values[1:2], cnms))
  expect_null(maihda_lmer_block_scales(list(cbind(`(Intercept)` = 1, z = x)), list(stratum = c("(Intercept)", "x"))))
  # The theta carried there and back is the theta, and gives each row the variance it had.
  theta <- c(0.4, -0.02, 0.003, 0.0004, 0.7, 0.2, 0.1)
  there <- maihda_lmer_transform_theta(theta, cnms, lapply(Ts, solve))
  expect_equal(maihda_lmer_transform_theta(there, cnms, Ts), theta)
  L <- matrix(c(theta[1], theta[2], 0, theta[3]), 2)
  La <- matrix(c(there[1], there[2], 0, there[3]), 2)
  for (xi in c(2010, 2015, 2021)) {
    row <- c(1, xi)
    z <- c(1, (xi - mean(x)) / stats::sd(x))
    expect_equal(drop(z %*% tcrossprod(La) %*% z), drop(row %*% tcrossprod(L) %*% row))
  }
  expect_equal(there[4]^2 * (2015 / sqrt(mean(x^2)))^2, theta[4]^2 * 2015^2)
  expect_null(maihda_lmer_transform_theta(theta[-1], cnms, Ts))
})

test_that("a slope in a block without an intercept is measured with that block scaled alone", {
  skip_on_cran()
  d <- calendar_data(3)
  uncorrelated <- y ~ x + (1 + x || stratum)
  # fit_maihda()'s fit of seed 3, 4.4e-8 above the optimum, where lme4 fires.
  at <- placed(c(0.00035376665, 0.00014952541), d, uncorrelated)
  expect_false(lme4::isSingular(at))
  expect_gt(length(at@optinfo$conv$lme4$messages), 0L)
  dg <- maihda_fit_diagnostics(at)
  expect_true(dg$converged)
  expect_lt(dg$stationary_distance, 0.01)
  # 0.015 above the optimum, where lme4 fires too: not converged, and how far.
  short <- placed(c(0.0050111125, 0.00015361898), d, uncorrelated)
  gap <- lme4::REMLcrit(short) - lme4::REMLcrit(at)
  expect_lt(abs(gap / 0.015 - 1), 1e-3)
  expect_gt(length(short@optinfo$conv$lme4$messages), 0L)
  ds <- maihda_fit_diagnostics(short)
  expect_false(ds$converged)
  expect_lt(abs(ds$stationary_distance / sqrt(gap) - 1), 0.05)
  # The year in a block with an intercept and in one without: seed 8 at its optimum.
  mixed <- placed(c(147.58778, -0.073194586, 0.00012321111, 2.6887156e-05), calendar_data(8),
                  y ~ x + (1 + x | stratum) + (0 + x | ctx))
  expect_lt(maihda_lmer_stationary_distance(mixed), 0.001)
})

# ---- the record ---------------------------------------------------------------------------

test_that("an lmer fit at the REML optimum is converged although lme4's checks fire", {
  skip_on_cran()
  m <- placed(optimum_3)
  said <- m@optinfo$conv$lme4$messages
  # The premise: lme4's checks fire at the optimum, and its optimiser reports success.
  expect_gt(length(said), 0L)
  expect_true(any(grepl("failed to converge", said, fixed = TRUE)))
  expect_identical(as.integer(m@optinfo$conv$opt), 0L)
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$converged)
  expect_identical(dg$messages, character(0))
  expect_identical(dg$lme4_overruled, as.character(said))
  expect_lt(dg$stationary_distance, 1e-4)
  expect_identical(maihda_format_fit_diagnostics(dg),
                   c("lme4 warned that the model failed to converge; in centred and scaled units the",
                     "  fit is less than 0.0001 standard errors from the maximum."))
})

test_that("a fit by maximum likelihood is measured on its deviance", {
  skip_on_cran()
  m <- placed(optimum_3_ml, REML = FALSE)
  expect_false(lme4::isREML(m))
  expect_gt(length(m@optinfo$conv$lme4$messages), 0L)
  expect_gt(lme4::REMLcrit(m) - stats::deviance(m), 10)
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$converged)
  expect_lt(dg$stationary_distance, 1e-4)
  expect_identical(dg$lme4_overruled, as.character(m@optinfo$conv$lme4$messages))
})

test_that("a fit short of the optimum stays not converged, and says how far", {
  skip_on_cran()
  m <- placed(above_0.03)
  said <- as.character(m@optinfo$conv$lme4$messages)
  expect_gt(length(said), 0L)
  dg <- maihda_fit_diagnostics(m)
  expect_false(dg$converged)
  expect_identical(dg$messages, said)
  expect_null(dg$lme4_overruled)
  expect_gt(dg$stationary_distance, 0.15)
  lines <- maihda_format_fit_diagnostics(dg)
  expect_identical(lines[1], "Convergence warnings reported by lme4:")
  expect_identical(utils::tail(lines, 2L),
                   c("In centred and scaled units the fit's gradient and Hessian put it 0.18 standard",
                     "  errors from a stationary point."))
})

test_that("the line is a tenth of a standard error, from both sides", {
  skip_on_cran()
  d <- year_data(3)
  inside <- maihda_fit_diagnostics(placed(above_0.0081, d))
  expect_lt(inside$stationary_distance, 0.1)
  expect_true(inside$converged)
  outside <- maihda_fit_diagnostics(placed(above_0.0121, d))
  expect_gt(outside$stationary_distance, 0.1)
  expect_false(outside$converged)
  # At the line itself.
  m <- placed(optimum_3, d)
  at <- function(distance) {
    local_mocked_bindings(maihda_lmer_stationary_distance = function(model) distance)
    maihda_fit_diagnostics(m)
  }
  expect_true(at(0.0999)$converged)
  on_line <- at(0.1)
  expect_false(on_line$converged)
  expect_identical(on_line$stationary_distance, 0.1)
  expect_null(on_line$lme4_overruled)
})

test_that("fit_maihda() records a birth-year slope fit at its optimum as converged", {
  skip_on_cran()
  m <- quiet(fit_maihda(slope, data = year_data(3)))
  dg <- m$diagnostics
  # lme4's checks fire on the fit fit_maihda() keeps, which is at the optimum.
  expect_gt(length(dg$lme4_overruled), 0L)
  expect_lt(lme4::REMLcrit(m$model) - 1994.85420, 1e-4)
  expect_true(dg$converged)
  expect_lt(dg$stationary_distance, 0.01)
  shown <- utils::capture.output(print(m))
  expect_true(any(grepl("lme4 warned that the model failed to converge", shown, fixed = TRUE)))
  expect_false(any(grepl("Convergence warnings reported by lme4", shown, fixed = TRUE)))
})

test_that("growth and intercept-only lmer fits are measured too, and a fit lme4 does not flag is not", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA", envir = environment())
  ld <- maihda_long_data[maihda_long_data$id %in% unique(maihda_long_data$id)[1:80], ]
  g <- quiet(fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = ld,
                        id = "id", time = "wave"))
  # lme4 does not flag this growth fit: nothing is measured.
  expect_length(g$model@optinfo$conv$lme4$messages, 0L)
  expect_true(g$diagnostics$converged)
  expect_null(g$diagnostics$stationary_distance)
  expect_null(g$diagnostics$lme4_overruled)
  # Had lme4 flagged it, the fit would be measured, in its growth terms centred and scaled.
  said <- list("Model is nearly unidentifiable: very large eigenvalue\n - Rescale variables?")
  flagged <- g$model
  flagged@optinfo$conv$lme4$messages <- said
  dg <- maihda_fit_diagnostics(flagged, g$longitudinal_info)
  expect_true(dg$converged)
  expect_lt(dg$stationary_distance, 0.01)
  expect_identical(dg$lme4_overruled, unlist(said))
  # An intercept-only block has nothing to standardise and is measured in theta as fitted.
  one <- quiet(lme4::lmer(y ~ x + (1 | stratum), data = year_data(3)))
  one@optinfo$conv$lme4$messages <- said
  d1 <- maihda_fit_diagnostics(one)
  expect_true(d1$converged)
  expect_lt(d1$stationary_distance, 0.01)
})

test_that("an optimiser failure, a deviance that does not reproduce the fit and a saddle keep lme4's verdict", {
  skip_on_cran()
  m <- placed(optimum_3)
  said <- as.character(m@optinfo$conv$lme4$messages)
  # The optimiser's own failure is not measured again.
  failed <- m
  failed@optinfo$conv$opt <- 1L
  failed@optinfo$message <- "stopped at the evaluation limit"
  df <- maihda_fit_diagnostics(failed)
  expect_false(df$converged)
  expect_true("stopped at the evaluation limit" %in% df$messages)
  expect_null(df$stationary_distance)
  expect_null(df$lme4_overruled)
  # A deviance that does not reproduce the fit's criterion, and one with the fit's criterion
  # at the fit but curved down along one direction.
  deviance <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
  expect_lt(maihda_lmer_distance_with(m, deviance), 1e-4)
  expect_true(is.na(maihda_lmer_distance_with(m, function(th) deviance(th) + 1)))
  target <- lme4::REMLcrit(m)
  own <- unname(lme4::getME(m, "theta"))
  saddle <- function(th) target + sum(c(1, -1, 1) * ((th - own) / c(1, 0.001, 0.001))^2)
  expect_true(is.na(maihda_lmer_distance_with(m, saddle)))
  # A measure that cannot be taken keeps lme4's verdict.
  local_mocked_bindings(maihda_lmer_stationary_distance = function(model) NA_real_)
  dn <- maihda_fit_diagnostics(m)
  expect_false(dn$converged)
  expect_identical(dn$messages, said)
  expect_null(dn$stationary_distance)
  expect_null(dn$lme4_overruled)
})

test_that("only lmer fits are measured this way, and measuring leaves the fit as lmer made it", {
  skip_on_cran()
  expect_true(is.na(maihda_lmer_stationary_distance(list())))
  d <- year_data(3)
  d$yb <- as.integer(d$y > stats::median(d$y))
  gl <- quiet(lme4::glmer(yb ~ x + (1 | stratum), data = d, family = stats::binomial))
  expect_true(is.na(maihda_lmer_stationary_distance(gl)))
  m <- placed(optimum_3, d)
  expect_true(is.na(maihda_lme4_stationary_distance(m)))
  # Copies, not references: lme4 writes into some of the fit's vectors in place.
  state <- function(m) unserialize(serialize(list(
    theta = lme4::getME(m, "theta"), lambda = lme4::getME(m, "Lambdat")@x + 0,
    b = lme4::getME(m, "b"), u = lme4::getME(m, "u"), fixef = lme4::fixef(m),
    ranef = lme4::ranef(m), sigma = stats::sigma(m), crit = lme4::REMLcrit(m),
    predict = stats::predict(m)), NULL))
  before <- state(m)
  first <- maihda_lmer_stationary_distance(m)
  expect_identical(state(m), before)
  dg <- maihda_fit_diagnostics(m)
  expect_identical(state(m), before)
  expect_identical(dg$stationary_distance, first)
})
