# Audit pass 2026-10-25: fits with a diag() random-effect block (lme4 2.0 on) measured in the
# units of the formula, or not at all.
#
# The convergence re-checks re-express each random-effect block with its terms centred and
# scaled (test-audit-2026-10-20-lmer-optimum-messages.R, test-audit-2026-10-20-glmer-slope-scale.R).
# A diag(1 + x | g) block holds its terms' standard deviations in theta, not a Cholesky factor,
# so the re-expression could not read it: a glmer fit was measured in theta as fitted, where a
# step along a birth-year slope's entry is most of its value (4 of 8 simulated Poisson fits at
# the maximum were recorded as not converged, read at 0.5 to 1.5 standard errors or NA), and an lmer
# fit was not measured at all, keeping lme4's verdict wherever its checks fired. Per the user,
# the terms of a diagonal block are divided by their root mean square and not centred (the block
# cannot hold the covariances centring adds), on both engines: the same model and the same
# measure as the block written with ||. On 270 points placed at known distances (7 designs) the
# measure equalled that of the || spelling to 3e-7 and was 0.89 to 1.05 times the true distance
# (1% to 99%). Other structures (cs(), ar1()) are left as they were.
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

# Eighteen additive strata (3 x 3 x 2 dimensions, effects SD 0.3) of 80 rows, each row with a
# birth year between 1940 and 2000 and a stratum slope deviation of SD 0.01 per year (centred at
# 1970): yc a count.
year_counts <- function(seed) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, 0.01)
    st <- rep(seq_len(18), each = 80)
    x <- sample(1940:2000, length(st), replace = TRUE)
    eta <- 0.5 + add[st] + (0.01 + v[st] / 2) * (x - 1970)
    data.frame(yc = stats::rpois(length(st), exp(eta)), x = x, stratum = factor(st))
  })
}

# Eighteen strata of 60 rows with an age x between 18 and 80 and a standard normal w, the
# stratum intercepts (SD 0.4) and age slopes (SD 0.012 per year) independent: y Gaussian or,
# with `family = "poisson"`, a count.
age_data <- function(seed, family = "gaussian") {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    st <- rep(seq_len(18), each = 60)
    x <- sample(18:80, length(st), replace = TRUE)
    w <- stats::rnorm(length(st))
    u <- stats::rnorm(18, 0, 0.4)
    v <- stats::rnorm(18, 0, 0.012)
    eta <- add[st] + u[st] + (0.01 + v[st]) * x + 0.3 * w
    d <- data.frame(x = x, w = w, stratum = factor(st))
    if (family == "gaussian") {
      d$y <- 1 + eta + stats::rnorm(length(st))
    } else {
      d$y <- stats::rpois(length(st), exp(0.3 + 0.4 * eta))
    }
    d
  })
}

# Fits placed at a point: optimisers that return their start, so lme4's checks are taken at
# that point exactly (boundary.tol = 0 keeps lme4 from moving an entry near its bound).
placed_lmer <- function(formula, data, theta) {
  quiet(lme4::lmer(formula, data = data, start = theta, control = lme4::lmerControl(
    optimizer = function(fn, par, lower, upper, control) {
      list(par = par, fval = fn(par), conv = 0, feval = 1)
    }, boundary.tol = 0)))
}
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}
placed_glmer <- function(formula, data, theta, beta) {
  quiet(lme4::glmer(formula, data = data, family = stats::poisson,
                    start = list(theta = theta, fixef = beta),
                    control = lme4::glmerControl(optimizer = list(stay, stay),
                                                 nAGQ0initStep = FALSE, boundary.tol = 0)))
}
deviance_of <- function(m) -2 * as.numeric(stats::logLik(m))

# The glmer measure with the fixed effects centred and scaled and theta as fitted.
fitted_theta_distance <- function(m) {
  dev <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
  theta <- lme4::getME(m, "theta")
  beta <- lme4::fixef(m)
  k <- length(theta)
  A <- maihda_centre_scale_map(lme4::getME(m, "X"))
  f <- function(p) dev(c(p[seq_len(k)], drop(A %*% p[-seq_len(k)])))
  at <- c(theta, solve(A, beta))
  lower <- c(lme4::getME(m, "lower"), rep(-Inf, length(beta)))
  maihda_newton_distance(maihda_central_derivatives(f, at, lower = lower, fx = f(at)))
}

# ---- the layout and the re-expression ---------------------------------------------------------

test_that("a block held as a diagonal is read from theta's bounds", {
  cnms <- list(a = c("(Intercept)", "x"), b = "(Intercept)", c = c("(Intercept)", "x", "z"))
  # Full Cholesky factors: 0 on the diagonal entries, no bound below.
  expect_identical(maihda_lme4_diagonal_blocks(cnms, c(0, -Inf, 0, 0, 0, -Inf, -Inf, 0, -Inf, 0)),
                   c(FALSE, FALSE, FALSE))
  # Diagonal blocks hold one standard deviation per term, each bounded at 0.
  expect_identical(maihda_lme4_diagonal_blocks(cnms, c(0, 0, 0, 0, 0, 0)), c(TRUE, FALSE, TRUE))
  expect_identical(maihda_lme4_diagonal_blocks(cnms[c(3, 1)], c(0, 0, 0, 0, -Inf, 0)),
                   c(TRUE, FALSE))
  expect_identical(maihda_lme4_diagonal_blocks(cnms["b"], c(stratum = 0)), FALSE)
  # cs() and ar1() hold a correlation with its own bound; bounds that do not fit the blocks.
  expect_null(maihda_lme4_diagonal_blocks(cnms["c"], c(0, 0, 0, -0.5)))
  expect_null(maihda_lme4_diagonal_blocks(cnms["c"], c(0, -1)))
  # A two-term ar1() block: a standard deviation and a correlation, as many bounds as terms.
  expect_null(maihda_lme4_diagonal_blocks(list(g = c("f1", "f2")), c(0, -1)))
  expect_null(maihda_lme4_diagonal_blocks(cnms["a"], c(0, -Inf, 0, 0)))
  expect_null(maihda_lme4_diagonal_blocks(cnms["a"], 0))
})

test_that("a diagonal block's terms are only scaled, its standard deviations carried there and back", {
  x <- c(1955, 1962, 1970, 1981, 1990, 1998)
  values <- list(stratum = cbind(`(Intercept)` = 1, x = x), ctx = cbind(`(Intercept)` = 1, x = x))
  cnms <- lapply(values, colnames)
  Ts <- maihda_lmer_block_scales(values, cnms, diagonal = c(TRUE, FALSE))
  # Divided by its root mean square, not centred, in the diagonal block; centred and scaled in
  # the other.
  expect_equal(Ts[[1]], diag(c(1, 1 / sqrt(mean(x^2)))))
  expect_equal(Ts[[2]], matrix(c(1, 0, -mean(x) / stats::sd(x), 1 / stats::sd(x)), 2))
  expect_identical(maihda_lmer_block_scales(values, cnms), maihda_lmer_block_scales(values, cnms, c(FALSE, FALSE)))
  # theta: two standard deviations for the diagonal block, a factor for the other.
  theta <- c(0.6, 0.004, 0.5, -0.0002, 0.0003)
  there <- maihda_lmer_transform_theta(theta, cnms, lapply(Ts, solve), c(TRUE, FALSE))
  expect_equal(there[1:2], c(0.6, 0.004 * sqrt(mean(x^2))))
  expect_equal(maihda_lmer_transform_theta(there, cnms, Ts, c(TRUE, FALSE)), theta)
  # Each row's random-effect variance in the diagonal block is unchanged.
  for (xi in c(1950, 1975, 2001)) {
    expect_equal(there[1]^2 + there[2]^2 * (xi / sqrt(mean(x^2)))^2, theta[1]^2 + theta[2]^2 * xi^2)
  }
  # A standard deviation's sign does not matter; a T that adds covariances cannot be held.
  expect_equal(maihda_lmer_transform_theta(c(-0.6, 0.004, theta[3:5]), cnms, Ts, c(TRUE, FALSE))[1],
               0.6)
  expect_null(maihda_lmer_transform_theta(theta, cnms, list(Ts[[2]], Ts[[2]]), c(TRUE, FALSE)))
  # Theta of the wrong length for the layout.
  expect_null(maihda_lmer_transform_theta(theta, cnms, Ts, c(FALSE, FALSE)))
})

# ---- glmer ---------------------------------------------------------------------------------

test_that("a glmer fit with a diag() slope on a birth year is measured as the same model written with ||", {
  skip_on_cran()
  skip_if(utils::packageVersion("lme4") < "2.0.0")
  d <- year_counts(1)
  # Seed 1's fit by glmer(), at the maximum: lme4 reports max|grad| = 2.1 and very large
  # eigenvalues there.
  point <- list(theta = c(0.00010182583, 0.00018368866), beta = c(-23.759184, 0.012387786))
  m <- placed_glmer(yc ~ x + diag(1 + x | stratum), d, point$theta, point$beta)
  said <- as.character(m@optinfo$conv$lme4$messages)
  expect_false(lme4::isSingular(m))
  expect_true(any(grepl("failed to converge", said, fixed = TRUE)))
  split <- placed_glmer(yc ~ x + (1 + x || stratum), d, point$theta, point$beta)
  expect_lt(abs(deviance_of(m) - deviance_of(split)), 1e-6)
  distance <- maihda_lme4_stationary_distance(m)
  expect_lt(distance, 0.01)
  expect_lt(abs(distance - maihda_lme4_stationary_distance(split)), 1e-6)
  # Theta as fitted reads 1.5 standard errors.
  expect_gt(fitted_theta_distance(m), 0.1)
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$converged)
  expect_identical(dg$messages, character(0))
  expect_identical(dg$lme4_overruled, said)
})

test_that("a three-term diag() block is measured as its terms written with ||", {
  skip_on_cran()
  skip_if(utils::packageVersion("lme4") < "2.0.0")
  d <- age_data(1, "poisson")
  # Seed 1's maximum, where lme4's checks fire.
  theta <- c(0.17494762, 0.003897966, 0.056043894)
  beta <- c(0.32120264, 0.0064098783, 0.11790069)
  m <- placed_glmer(y ~ x + w + diag(1 + x + w | stratum), d, theta, beta)
  expect_gt(length(m@optinfo$conv$lme4$messages), 0L)
  split <- placed_glmer(y ~ x + w + (1 + x + w || stratum), d, theta, beta)
  distance <- maihda_lme4_stationary_distance(m)
  expect_lt(distance, 0.01)
  expect_lt(abs(distance - maihda_lme4_stationary_distance(split)), 1e-6)
})

# ---- lmer ----------------------------------------------------------------------------------

test_that("an lmer fit with a diag() block is measured, from both sides of the line", {
  skip_on_cran()
  skip_if(utils::packageVersion("lme4") < "2.0.0")
  d <- age_data(1)
  # The REML optimum of seed 1, and points along one direction of the scaled standard
  # deviations 0.0001, 0.0081, 0.0121 and 0.0225 above it in the criterion; lme4's checks fire
  # at each of those points (max|grad| from 0.06 to 0.93).
  optimum <- c(0.35497227, 0.011526904)
  points <- list(list(theta = c(0.35586141, 0.011509911), gap = 0.0001, converged = TRUE),
                 list(theta = c(0.36296939, 0.011374062), gap = 0.0081, converged = TRUE),
                 list(theta = c(0.36474499, 0.011340126), gap = 0.0121, converged = FALSE),
                 list(theta = c(0.36829456, 0.011272286), gap = 0.0225, converged = FALSE))
  top <- lme4::REMLcrit(placed_lmer(y ~ x + diag(1 + x | stratum), d, optimum))
  for (p in points) {
    m <- placed_lmer(y ~ x + diag(1 + x | stratum), d, p$theta)
    gap <- lme4::REMLcrit(m) - top
    expect_lt(abs(gap - p$gap), 1e-6)
    said <- as.character(m@optinfo$conv$lme4$messages)
    expect_gt(length(said), 0L)
    distance <- maihda_lmer_stationary_distance(m)
    # A ratio, not expect_equal(tolerance =): under this suite's testthat edition a scalar
    # passes that whenever the absolute difference is below the tolerance.
    expect_lt(abs(distance / sqrt(gap) - 1), 0.02)
    split <- placed_lmer(y ~ x + (1 + x || stratum), d, p$theta)
    expect_lt(abs(distance - maihda_lmer_stationary_distance(split)), 1e-6)
    dg <- maihda_fit_diagnostics(m)
    expect_identical(dg$converged, p$converged)
    if (p$converged) {
      expect_identical(dg$lme4_overruled, said)
    } else {
      expect_identical(dg$messages, said)
    }
  }
})

test_that("a cs() block is left as it was", {
  skip_on_cran()
  skip_if(utils::packageVersion("lme4") < "2.0.0")
  m <- quiet(lme4::lmer(y ~ x + cs(1 + x + w | stratum), data = age_data(1)))
  # lme4 reports theta as a factor of 6 entries for 4 parameters.
  expect_false(length(lme4::getME(m, "theta")) == length(lme4::getME(m, "lower")))
  expect_null(maihda_lme4_diagonal_blocks(lme4::getME(m, "cnms"), lme4::getME(m, "lower")))
  expect_true(is.na(maihda_lmer_stationary_distance(m)))
})
