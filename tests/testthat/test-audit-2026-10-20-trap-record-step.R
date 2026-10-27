# Audit pass 2026-10-20: the convergence record's fall across a variance bound.
#
# The record measured a fit stopped at a Cholesky bound by the fall of its deviance with the
# diagonal entry at a fixed step of 0.001. Next to the entries of a slope on age in years or a
# covariate in thousands that step is far too long: the deviance falls only over a much
# shorter one and has risen by 0.001, and the fall along one column's entries can be small
# where the likelihood rises by much more as the block's other entries move with them. Fits
# short of the maximum were recorded as converged. The record now also takes steps of 0.01
# and 0.1, and a Newton step in the block's entries with its terms centred and scaled, and
# its sentence says what showed the likelihood rising.
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

# Optimisers that keep their start: a fit held at a given point.
stay_lmer <- function(fn, par, lower, upper, control) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}
stay_glmer <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0L, message = "start kept")
}
hold_lmer <- lme4::lmerControl(optimizer = stay_lmer, calc.derivs = FALSE)
hold_glmer <- lme4::glmerControl(nAGQ0initStep = FALSE, optimizer = stay_glmer, calc.derivs = FALSE)

# Twelve strata (A 3 x B 4) of a Gaussian outcome with stratum slopes on age in years (20 to 70,
# uncentred) and on z, fanning out from age 0; with z = FALSE the block is (1 + age | stratum).
age_data <- function(seed, z = TRUE) {
  with_seed(seed, {
    n <- 1500
    A <- sample(0:2, n, TRUE)
    B <- sample(0:3, n, TRUE)
    st <- 4L * A + B + 1L
    age <- round(stats::runif(n, 20, 70), 1)
    if (z) {
      zz <- stats::rbinom(n, 1, 0.5)
      ua <- stats::rnorm(12, 0, 0.02)
      uz <- stats::rnorm(12, 0, 0.3)
      y <- 1 + 0.02 * age + 0.3 * zz + c(0, 0.3, 0.5)[A + 1] + c(0, -0.2, 0.1, 0.3)[B + 1] +
        ua[st] * age + uz[st] * zz + stats::rnorm(n)
      data.frame(y = y, age = age, z = zz, A = factor(A), B = factor(B), stratum = factor(st))
    } else {
      ua <- stats::rnorm(12, 0, 0.02)
      y <- 1 + 0.02 * age + c(0, 0.3, 0.5)[A + 1] + c(0, -0.2, 0.1, 0.3)[B + 1] +
        ua[st] * age + stats::rnorm(n)
      data.frame(y = y, age = age, A = factor(A), B = factor(B), stratum = factor(st))
    }
  })
}
age3 <- y ~ age + z + A + B + (1 + age + z | stratum)
age2 <- y ~ age + A + B + (1 + age | stratum)

# The same strata with a slope on x from 1,000 to 9,000 and on z ~ N(0, 1): a Poisson outcome
# (count = TRUE) or a Gaussian one.
thousands_data <- function(seed, count = TRUE) {
  with_seed(seed, {
    n <- if (count) 2000 else 1500
    A <- sample(0:2, n, TRUE)
    B <- sample(0:3, n, TRUE)
    st <- 4L * A + B + 1L
    x <- round(stats::runif(n, 1000, 9000))
    z <- round(stats::rnorm(n), 2)
    if (count) {
      ux <- stats::rnorm(12, 0, 4e-5)
      uz <- stats::rnorm(12, 0, 0.2)
      eta <- 0.2 + 5e-5 * x + 0.1 * z + c(0, 0.3, 0.5)[A + 1] + c(0, -0.2, 0.1, 0.3)[B + 1] +
        ux[st] * x + uz[st] * z
      y <- stats::rpois(n, exp(eta))
    } else {
      ux <- stats::rnorm(12, 0, 1e-4)
      uz <- stats::rnorm(12, 0, 0.3)
      y <- 1 + 1e-4 * x + 0.3 * z + c(0, 0.3, 0.5)[A + 1] + c(0, -0.2, 0.1, 0.3)[B + 1] +
        ux[st] * x + uz[st] * z + stats::rnorm(n)
    }
    data.frame(y = y, x = x, z = z, A = factor(A), B = factor(B), stratum = factor(st))
  })
}
thousands3 <- y ~ x + z + A + B + (1 + x + z | stratum)

# Points found by holding one column of the block at its bound and optimising the rest, and
# the REML optimum of each data set from a search of lmer's deviance function with no bound on
# the Cholesky diagonal. age_data(4): the age column at its bound, 0.1834 below the optimum
# 4304.96131 in REML log-likelihood; age_data(5, z = FALSE): the intercept column, 0.1556 below
# 4353.135147. lme4's own Poisson fit of thousands_data(1) stops with the x column at its bound,
# 0.270 log-likelihood units below the maximum, deviance 6905.36164.
trap_age <- c(0.03802615, 0.01748130, -0.11735360, 0, -0.15898450, 0.23796110)
optimum_age <- 4304.96131
# age_data(18): the age column at its bound, 0.1698 below the optimum 4289.875741, where the
# criterion falls along the first-order direction and rises by the reflection.
trap_first <- c(0.2301560, 0.01160693, -0.0205690, 0, -0.1926840, 0.0001716914)
optimum_first <- 4289.875741
# age_data(28): the intercept column at its bound, 0.4384 below the optimum 4328.702263, where
# the first-order direction falls further than the Newton step and the reflection rises.
trap_largest_first <- c(0, -0.01413018, 0.03969768, 0.0178475, 0.2368465, 0.0006430777)
optimum_largest_first <- 4328.702263
# age_data(39): the age column at its bound, 0.607 below the optimum 4347.818748, where the
# reflection falls by less than counts (0.00074) and the Newton step by 0.11.
trap_small_reflection <- c(0.06549441, -0.01889301, -0.09306663, 0, 0.1311113, 0.1888992)
optimum_small_reflection <- 4347.818748
# lmer's own fit of thousands_data(27, count = FALSE): 6.85 below the optimum 4421.387935, the
# x column's diagonal entry 5.6e-5 at its natural size (most of the x slope's variance).
stop_thousands <- c(0.3640868, -3.068656e-05, 0.2053429, 5.588063e-05, 0.1370348, 0.4953119)
optimum_thousands <- 4421.387935
trap_saddle <- c(0, -0.002481465, 0.022990900)
optimum_saddle <- 4353.135147
trap_count <- c(4.979712e-02, -4.226155e-05, 2.781052e-02, 0, -3.848508e-03, 1.832879e-01)
fixef_count <- c(2.055697e-01, 4.107582e-05, 1.515375e-02, 2.505560e-01, 4.139357e-01,
                 -2.446245e-01, 1.853708e-01, 2.935087e-01)
maximum_count <- 6905.36164

moved <- paste0("and the likelihood still rises as the block's variances and covariances move: ",
                "the fit stopped at a bound of lme4's parametrisation, not at the maximum.")

# ---- the Newton step --------------------------------------------------------------------

test_that("a fixed step can overshoot where a Newton step in the block's own units falls", {
  # A deviance of the block's covariance in units where x is divided by 1000: the x column of
  # (1 + x + z | g) at its bound beside z's entry, far from the minimum (0).
  cnms <- list(g = c("(Intercept)", "x", "z"))
  s <- 1000
  scale_x <- diag(c(1, 1 / s, 1))
  target <- matrix(c(0.30, 0.05, -0.04, 0.05, 0.20, 0.06, -0.04, 0.06, 0.25), 3)
  cov_of <- function(th) {
    L <- matrix(0, 3, 3)
    L[lower.tri(L, diag = TRUE)] <- th
    tcrossprod(L)
  }
  dev <- function(th) {
    S <- solve(scale_x) %*% cov_of(th) %*% t(solve(scale_x))
    1000 * sum((S - target)^2)
  }
  trap <- c(0.5, 0.3 / s, -0.1, 0, 0.4, 0.3)
  expect_equal(dev(trap), 51.1, tolerance = 1e-8)
  # Every step of the record's rays rises; a much shorter one falls.
  for (step in c(1e-3, 0.01, 0.1)) {
    expect_true(all(maihda_lme4_column_falls(trap, cnms, 4L, dev, step) < -1))
  }
  expect_gt(maihda_lme4_column_falls(trap, cnms, 4L, dev, 1e-5)[["first_order"]], 1)
  fall <- maihda_lme4_block_newton_fall(trap, cnms, 1L, dev, list(scale_x))
  expect_gt(fall, 5)
  expect_lte(fall, dev(trap))
  # At the minimum the step finds nothing; NA where the deviance cannot be evaluated or is flat.
  L <- t(chol(scale_x %*% target %*% t(scale_x)))
  minimum <- L[lower.tri(L, diag = TRUE)]
  expect_lt(abs(maihda_lme4_block_newton_fall(minimum, cnms, 1L, dev, list(scale_x))), 1e-12)
  expect_identical(maihda_lme4_block_newton_fall(trap, cnms, 1L, function(th) NA_real_,
                                                 list(scale_x)), NA_real_)
  expect_identical(maihda_lme4_block_newton_fall(trap, cnms, 1L, function(th) 3, list(scale_x)),
                   NA_real_)
  # Only the named block's entries move: a deviance of the other block alone gives the first
  # nothing to step along.
  two <- list(g = cnms$g, h = c("(Intercept)", "w"))
  th2 <- c(trap, 0.3, 0.1, 0.2)
  of_h <- function(th) sum((cov_of(c(th[7], th[8], 0, th[9], 0, 0))[1:2, 1:2] - diag(c(0.2, 0.1)))^2)
  expect_identical(maihda_lme4_block_newton_fall(th2, two, 1L, of_h, list(scale_x, diag(2))),
                   NA_real_)
  expect_gt(maihda_lme4_block_newton_fall(th2, two, 2L, of_h, list(scale_x, diag(2))), 0.01)
})

test_that("the record's two falls across the bound are the reflection's and the first-order one", {
  cnms <- list(g = c("(Intercept)", "x"))
  theta <- c(0, 0.3, 0.4)
  dev <- function(th) 12 * th[1] * th[2] + 2.4e4 * (th[1] * th[2])^2
  falls <- maihda_lme4_column_falls(theta, cnms, 1L, dev, 1e-3)
  expect_identical(names(falls), c("reflection", "first_order"))
  expect_identical(falls[["reflection"]], dev(theta) - dev(maihda_lme4_mirror_theta(theta, cnms, 1L, 1e-3)))
  expect_identical(maihda_lme4_column_drop(theta, cnms, 1L, dev, 1e-3), max(falls))
  expect_identical(unname(maihda_lme4_column_falls(theta, cnms, 1L, function(th) NA_real_)),
                   c(NA_real_, NA_real_))
})

# ---- the record on fits held at a bound -------------------------------------------------

test_that("a fit held where only the Newton step shows the likelihood rising is not converged", {
  skip_on_cran()
  d <- age_data(4)
  m <- quiet(fit_maihda(age3, data = d, start = trap_age, control = hold_lmer))
  expect_identical(unname(lme4::getME(m$model, "theta")), trap_age)
  gap <- lme4::REMLcrit(m$model) - optimum_age
  expect_gt(gap, 0.36)
  dg <- m$diagnostics
  expect_false(dg$converged)
  tr <- dg$boundary_traps
  expect_identical(tr$block, "stratum")
  expect_identical(tr$term, "age")
  expect_identical(tr$n_terms, 3L)
  expect_true(tr$at_bound)
  # Across the bound the criterion rises at every step; the Newton step reaches most of the gap.
  expect_lt(tr$across, 0)
  expect_lt(tr$reflected, 0)
  expect_gt(tr$drop, 0.8 * gap)
  expect_lte(tr$drop, gap + 1e-6)
  expect_true(paste("In the 'stratum' block the variance of age left after the intercept is zero,",
                    moved) %in% dg$messages)
  # In the formula's own units the same step finds nothing.
  cnms <- lme4::getME(m$model, "cnms")
  raw <- maihda_lme4_block_newton_fall(trap_age, cnms, 1L, maihda_lme4_theta_deviance(m$model),
                                       list(diag(3)))
  expect_lt(raw, .maihda_lme4_trap_drop)
})

test_that("where the first-order step shows the rise and the reflection does not, the sentence says so", {
  skip_on_cran()
  d <- age_data(18)
  m <- quiet(fit_maihda(age3, data = d, start = trap_first, control = hold_lmer))
  expect_identical(unname(lme4::getME(m$model, "theta")), trap_first)
  dg <- m$diagnostics
  expect_false(dg$converged)
  tr <- dg$boundary_traps
  expect_identical(tr$term, "age")
  expect_lt(tr$reflected, 0)
  expect_gt(tr$across, 0.05)
  expect_lte(tr$drop, lme4::REMLcrit(m$model) - optimum_first + 1e-6)
  expect_true(paste0("In the 'stratum' block the variance of age left after the intercept is zero, ",
                     "where the likelihood still rises when that variance and its covariance with ",
                     "z move off zero: the fit stopped at a bound of lme4's parametrisation, not at ",
                     "the maximum.") %in% dg$messages)
})

test_that("the record's fall is the first-order one where that falls furthest", {
  skip_on_cran()
  m <- quiet(fit_maihda(age3, data = age_data(28), start = trap_largest_first, control = hold_lmer))
  expect_identical(unname(lme4::getME(m$model, "theta")), trap_largest_first)
  tr <- m$diagnostics$boundary_traps
  expect_identical(tr$term, "(Intercept)")
  cnms <- lme4::getME(m$model, "cnms")
  newton <- maihda_lme4_block_newton_fall(trap_largest_first, cnms, 1L,
                                          maihda_lme4_theta_deviance(m$model),
                                          maihda_lmer_block_scales(lme4::getME(m$model, "mmList"), cnms))
  # The step across falls 0.796 on every platform; the Newton step's fall depends on it (0.738
  # on Windows, 0.791 on macOS).
  expect_gt(tr$across, newton)
  expect_lt(tr$reflected, 0)
  expect_identical(tr$drop, tr$across)
  expect_lte(tr$drop, lme4::REMLcrit(m$model) - optimum_largest_first + 1e-6)
})

test_that("a reflection that falls by less than counts does not choose the sentence", {
  skip_on_cran()
  m <- quiet(fit_maihda(age3, data = age_data(39), start = trap_small_reflection, control = hold_lmer))
  expect_identical(unname(lme4::getME(m$model, "theta")), trap_small_reflection)
  dg <- m$diagnostics
  expect_false(dg$converged)
  tr <- dg$boundary_traps
  expect_true(tr$at_bound)
  expect_gt(tr$reflected, 0)
  expect_lte(tr$reflected, .maihda_lme4_trap_drop)
  expect_lte(tr$across, .maihda_lme4_trap_drop)
  expect_gt(tr$drop, 0.05)
  expect_lte(tr$drop, lme4::REMLcrit(m$model) - optimum_small_reflection + 1e-6)
  expect_true(paste("In the 'stratum' block the variance of age left after the intercept is zero,",
                    moved) %in% dg$messages)
})

test_that("a short fit whose slope entry is below 1e-4 at its natural size names no bound", {
  skip_on_cran()
  m <- quiet(fit_maihda(thousands3, data = thousands_data(27, count = FALSE), start = stop_thousands,
                        control = hold_lmer))
  expect_identical(unname(lme4::getME(m$model, "theta")), stop_thousands)
  gap <- lme4::REMLcrit(m$model) - optimum_thousands
  dg <- m$diagnostics
  expect_false(dg$converged)
  tr <- dg$boundary_traps
  expect_identical(tr$term, "x")
  expect_false(tr$at_bound)
  # No step across the bound falls; the Newton step does, once halved seven times.
  expect_lt(tr$across, -1)
  expect_gt(tr$drop, 1)
  expect_lte(tr$drop, gap + 1e-6)
  expect_true(paste0("The 'stratum' block is not at the maximum: the likelihood still rises as ",
                     "its variances and covariances move.") %in% dg$messages)
  expect_false(any(grepl("is zero", dg$messages, fixed = TRUE)))
})

test_that("a column is at its bound when its own entry carries almost none of its term's variance", {
  cnms <- list(g = c("(Intercept)", "x", "z"))
  # theta: (L11, L21, L31, L22, L32, L33); the x column's entry L22 against its row (L21, L22).
  at <- function(theta, index, scales = NULL) maihda_lme4_column_at_bound(theta, cnms, index, scales)
  expect_true(at(c(0.3, 0.02, 0.1, 0, 0.2, 0.3), 4L))
  expect_true(at(c(0.3, 0.02, 0.1, 0.0002, 0.2, 0.3), 4L))
  expect_false(at(c(0.3, 0.02, 0.1, 0.0003, 0.2, 0.3), 4L))
  expect_false(at(c(0.3, -3e-5, 0.1, 5.6e-5, 0.2, 0.3), 4L))
  # A row of zeros is a variance of zero.
  expect_true(at(c(0.3, 0, 0.1, 0, 0.2, 0.3), 4L))
  # The first column: lme4's 1e-4 for an intercept; for a term in a block without an
  # intercept, that entry times the root mean square of the term's column.
  expect_true(at(c(5e-5, 0.02, 0.1, 0.01, 0.2, 0.3), 1L))
  slopes <- list(g = c("x", "z"))
  rms <- diag(c(1 / 5000, 1))
  expect_false(maihda_lme4_column_at_bound(c(4e-5, 0.2, 0.3), slopes, 1L, list(rms)))
  expect_true(maihda_lme4_column_at_bound(c(1e-9, 0.2, 0.3), slopes, 1L, list(rms)))
  expect_true(maihda_lme4_column_at_bound(c(4e-5, 0.2, 0.3), slopes, 1L))
  # The sentence for a column not at its bound names the block only.
  expect_identical(maihda_lme4_trap_message("stratum", "x", 3L, c("(Intercept)", "x", "z"), "interior"),
                   paste0("The 'stratum' block is not at the maximum: the likelihood still rises as ",
                          "its variances and covariances move."))
})

test_that("a saddle across the bound is found by the longer steps", {
  skip_on_cran()
  d <- age_data(5, z = FALSE)
  m <- quiet(fit_maihda(age2, data = d, start = trap_saddle, control = hold_lmer))
  expect_identical(unname(lme4::getME(m$model, "theta")), trap_saddle)
  gap <- lme4::REMLcrit(m$model) - optimum_saddle
  dg <- m$diagnostics
  expect_false(dg$converged)
  tr <- dg$boundary_traps
  expect_identical(tr$term, "(Intercept)")
  # At 0.001 the criterion falls by less than counts; by the reflection at 0.1 it falls by more.
  cnms <- lme4::getME(m$model, "cnms")
  dev <- maihda_lme4_theta_deviance(m$model)
  expect_lt(maihda_lme4_column_drop(trap_saddle, cnms, 1L, dev, 1e-3), .maihda_lme4_trap_drop)
  falls <- vapply(c(1e-3, 0.01, 0.1), function(s) {
    maihda_lme4_column_falls(trap_saddle, cnms, 1L, dev, s)
  }, numeric(2))
  expect_identical(tr$reflected, max(falls["reflection", ]))
  expect_identical(tr$across, max(falls))
  expect_gt(tr$reflected, 0.1)
  expect_gte(tr$drop, tr$across)
  expect_lte(tr$drop, gap + 1e-6)
  reflected_sentence <- paste0("The intercept variance of the 'stratum' block is at zero, where the ",
                               "likelihood still rises with the sign of its covariance with age ",
                               "reversed: the fit stopped at a bound of lme4's parametrisation, not ",
                               "at the maximum.")
  expect_true(reflected_sentence %in% dg$messages)
  # A record without the falls' columns reads as one the reflection showed at a bound.
  local_mocked_bindings(maihda_lme4_boundary_traps = function(model, ...) {
    data.frame(block = "stratum", term = "(Intercept)", index = 1L, drop = 1, n_terms = 2L,
               stringsAsFactors = FALSE)
  })
  old <- maihda_fit_diagnostics(m$model)
  expect_false(old$converged)
  expect_true(reflected_sentence %in% old$messages)
})

test_that("a Poisson fit held at lme4's stop with x in thousands is not converged", {
  skip_on_cran()
  d <- thousands_data(1)
  m <- quiet(fit_maihda(thousands3, data = d, family = "poisson",
                        start = list(theta = trap_count, fixef = fixef_count), control = hold_glmer))
  expect_identical(unname(lme4::getME(m$model, "theta")), trap_count)
  deviance_gap <- -2 * as.numeric(stats::logLik(m$model)) - maximum_count
  expect_gt(deviance_gap, 0.5)
  dg <- m$diagnostics
  expect_false(dg$converged)
  tr <- dg$boundary_traps
  expect_identical(tr$term, "x")
  expect_true(tr$at_bound)
  expect_lt(tr$across, -1)
  expect_gt(tr$drop, 0.5 * deviance_gap)
  expect_lt(tr$drop, deviance_gap + 0.01)
  expect_true(paste("In the 'stratum' block the variance of x left after the intercept is zero,",
                    moved) %in% dg$messages)
})

test_that("a fit at its optimum with an entry of a slope in thousands below 1e-4 stays converged", {
  skip_on_cran()
  m <- quiet(fit_maihda(thousands3, data = thousands_data(1, count = FALSE)))
  # The x column's diagonal sits below lme4's 1e-4 at its natural size, so the record measures it.
  expect_lt(lme4::getME(m$model, "theta")[[4]], 1e-4)
  tr <- maihda_lme4_boundary_traps(m$model)
  expect_identical(tr$term, "x")
  expect_lt(tr$drop, 1e-5)
  expect_true(m$diagnostics$converged)
  expect_null(m$diagnostics$boundary_traps)
  # A fit without such a column has the record's columns and no rows.
  plain <- quiet(fit_maihda(y ~ A + B + (1 | stratum), data = thousands_data(1, count = FALSE)))
  none <- maihda_lme4_boundary_traps(plain$model)
  expect_identical(names(none), c("block", "term", "index", "drop", "across", "reflected", "at_bound",
                                 "n_terms"))
  expect_identical(nrow(none), 0L)
})

# ---- the sentence -----------------------------------------------------------------------

test_that("the sentence says what showed the likelihood rising", {
  tail <- ": the fit stopped at a bound of lme4's parametrisation, not at the maximum."
  two <- c("(Intercept)", "x")
  three <- c("(Intercept)", "age", "z")
  # The reflection: the sentences as before.
  expect_identical(maihda_lme4_trap_message("stratum", "(Intercept)", 2L, two, "reflection"),
                   maihda_lme4_trap_message("stratum", "(Intercept)", 2L, two))
  expect_identical(maihda_lme4_trap_message("stratum", "age", 3L, three),
                   paste0("In the 'stratum' block the variance of age left after the intercept is ",
                          "zero, where the likelihood still rises with the sign of its covariance ",
                          "with z reversed", tail))
  # A first-order step across the bound.
  expect_identical(maihda_lme4_trap_message("stratum", "(Intercept)", 2L, two, "across"),
                   paste0("The intercept variance of the 'stratum' block is at zero, where the ",
                          "likelihood still rises when that variance and its covariance with x ",
                          "move off zero", tail))
  expect_identical(maihda_lme4_trap_message("stratum", "age", 3L, three, "across"),
                   paste0("In the 'stratum' block the variance of age left after the intercept is ",
                          "zero, where the likelihood still rises when that variance and its ",
                          "covariance with z move off zero", tail))
  expect_identical(maihda_lme4_trap_message("stratum", "(Intercept)", 3L, three, "across"),
                   maihda_lme4_trap_message("stratum", "(Intercept)", 3L, three, "reflection"))
  # Only the Newton step.
  expect_identical(maihda_lme4_trap_message("stratum", "(Intercept)", 2L, two, "newton"),
                   paste("The intercept variance of the 'stratum' block is at zero,", moved))
  expect_identical(maihda_lme4_trap_message("g", "x", 3L, c("x", "z", "w"), "newton"),
                   paste("The variance of x in the 'g' block is at zero,", moved))
  # Growth blocks.
  expect_identical(maihda_lme4_trap_message("id", "(Intercept)", 2L, NULL, "across"),
                   paste0("The intercept variance of the 'id' block is at zero, where the ",
                          "likelihood still rises when that variance and its covariance with the ",
                          "slope move off zero", tail))
  expect_identical(maihda_lme4_trap_message("id", "(Intercept)", 3L, NULL, "across"),
                   maihda_lme4_trap_message("id", "(Intercept)", 3L))
  expect_identical(maihda_lme4_trap_message("id", "(Intercept)", 2L, NULL, "newton"),
                   paste("The intercept variance of the 'id' block is at zero,", moved))
  expect_identical(maihda_lme4_trap_message("id", ".maihda_ctime", 3L, NULL, "across"),
                   paste0("In the 'id' block the slope is perfectly correlated with the intercept, ",
                          "where the likelihood still rises when what is left of its variance and ",
                          "its covariance with the curvature move off zero", tail))
  expect_identical(maihda_lme4_trap_message("id", ".maihda_ctime", 3L, NULL, "newton"),
                   paste("In the 'id' block the slope is perfectly correlated with the intercept,",
                         moved))
})
