# Audit pass 2026-10-20: lme4's singular-fit test read in the formula's units on a random slope
# far from zero or in large units.
#
# lme4 calls a fit singular when a diagonal entry of a random-effect block's relative Cholesky
# factor, theta as the formula has it, is below 1e-4. For a slope that entry is the slope's
# standard deviation given the earlier terms at zero, per unit of its covariate: on a birth year
# the intercept at year 0 makes it tiny though the block is not at the boundary, and lme4 then
# runs no convergence check. On simulated designs (18 strata of 40 rows, a stratum slope SD of
# 0.005 per year) 14 of the 88 birth-year fits lme4 calls singular at the REML optimum are
# interior with each block's terms centred and scaled, and print() said "Singular fit ... may be
# unreliable" for them. Per the user, an lmer fit outside a growth model takes lme4's test on
# each block with a slope re-expressed as maihda_lmer_stationary_distance() re-expresses it
# (centred and divided by the standard deviation with an intercept, divided by the root mean
# square without), with the line at 0.01 there: lme4's 1e-4 is below where optimisers stop at a
# boundary optimum in those terms (up to 9.3e-4 in lmer fits, 6.1e-3 in glmer fits), and no
# interior optimum simulated came below 0.0265. Blocks with an intercept only and growth fits
# keep lme4's test; lme4's own verdict is `lme4_singular` where the two differ, and print() says
# so in one sentence. The scaled block's factor comes from an unpivoted QR decomposition rather
# than from its covariance factored again, which left the slope's entry of a block of rank one
# at up to 5.8e-4.
#
# glmer fits: lme4 runs no convergence check on one it calls singular, and of 119 such
# simulated fits with a birth-year, age, N(0, 1) or currency slope, 37 were 0.01 to 23 deviance
# short of the maximum, 36 of them recorded converged. Per the user, a glmer fit with nAGQ of 1
# or more outside a growth model whose random effects hold a slope and that lme4 calls singular
# is searched (bobyqa from the fit, with its random-effect terms and fixed effects centred and
# scaled, stopped 0.01 below it) and counts as not converged when the search falls that far, as
# lmer fits do (`search_gain`); its singular verdict is then taken as an lmer fit's. A glmer fit
# with nAGQ 0, which cannot be searched, and one with intercepts only keep lme4's test and are
# not searched, as are fits with a diag(), cs() or ar1() block, whose theta lme4 bounds as
# standard deviations and correlations rather than as a Cholesky factor.
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
# centred at 1970) and a residual SD of 1.
birth_year_data <- function(seed, int_sd = 0) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    dims <- list(stats::rnorm(3, 0, 0.3), stats::rnorm(3, 0, 0.3), stats::rnorm(2, 0, 0.3))
    add <- dims[[1]][lev$d1] + dims[[2]][lev$d2] + dims[[3]][lev$d3]
    u <- stats::rnorm(18, 0, int_sd)
    v <- stats::rnorm(18, 0, 0.005)
    st <- rep(seq_len(18), each = 40)
    x <- sample(1940:2000, 720, replace = TRUE)
    y <- 1 + add[st] + u[st] + (0.3 + v[st]) * (x - 1970) + stats::rnorm(720)
    d <- data.frame(y = y, x = x, lev[st, ], stratum = factor(st))
    for (nm in c("d1", "d2", "d3")) d[[nm]] <- factor(d[[nm]])
    rownames(d) <- NULL
    d
  })
}

adjusted <- y ~ x + d1 + d2 + d3 + (1 + x | stratum)

# A binary outcome over 18 additive strata (dimension effects SD 0.4) with a stratum intercept
# deviation of SD 0.15 and a slope deviation of SD `slope_sd` per unit of x: x a birth year
# between 1940 and 2000 (deviations centred at 1970) or an income in currency units, mean
# 50,000 and SD 15,000 (centred at 50,000).
binary_data <- function(seed, n = 80, x_type = "year", slope_sd = 0.01, beta_x = 0.02) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.4)[lev$d1] + stats::rnorm(3, 0, 0.4)[lev$d2] +
      stats::rnorm(2, 0, 0.4)[lev$d3]
    u <- stats::rnorm(18, 0, 0.15)
    v <- stats::rnorm(18, 0, slope_sd)
    st <- rep(seq_len(18), each = n)
    x <- if (x_type == "year") sample(1940:2000, 18 * n, replace = TRUE) else
      stats::rnorm(18 * n, 50000, 15000)
    centre <- if (x_type == "year") 1970 else 50000
    y <- stats::rbinom(18 * n, 1, stats::plogis(-0.3 + add[st] + u[st] + (beta_x + v[st]) * (x - centre)))
    d <- data.frame(y = y, x = x, lev[st, ], stratum = factor(st))
    for (nm in c("d1", "d2", "d3")) d[[nm]] <- factor(d[[nm]])
    rownames(d) <- NULL
    d
  })
}

# A glmer fit placed at `theta` and `fixef` (an optimiser that returns its start), nAGQ 1 unless
# given.
keep <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0L, message = "start kept")
}
placed_glmer <- function(theta, fixef, data, formula = y ~ x + (1 + x | stratum), nAGQ = 1L) {
  control <- if (nAGQ == 0L) {
    lme4::glmerControl(optimizer = keep, boundary.tol = 0)
  } else {
    lme4::glmerControl(nAGQ0initStep = FALSE, optimizer = keep, boundary.tol = 0)
  }
  start <- if (nAGQ == 0L) list(theta = theta) else list(theta = theta, fixef = fixef)
  quiet(lme4::glmer(formula, data = data, family = stats::binomial, nAGQ = nAGQ, start = start,
                    control = control))
}

# Birth-year seed 7, y ~ x + (1 + x | stratum): fit_maihda()'s binomial fit, at a deviance of
# 1810.3243, 3.15 above the maximum (1807.1707 at theta (24.877306, -0.012721077, 0.00032340709)
# found by unbounded searches of the deviance function). lme4's entry for the slope there is
# 1.6e-6, the block's scaled entries are 0.66 and 3.6e-4.
year_7_short <- list(theta = c(8.39611037, -0.00459586462, 1.60042606e-06),
                     fixef = c(-35.859405873, 0.017979860362))
# Currency seed 3 (200 rows per stratum): fit_maihda()'s fit, at the maximum (an unbounded
# search falls 2.4e-8 below it), with lme4's entry for the slope at 6.0e-6 and the block's
# scaled entries 0.39 and 0.18.
currency_3_max <- list(theta = c(0.788190538, -1.086026624e-05, 5.961999791e-06),
                       fixef = c(-1.5176879118, 2.22108310965e-05))
# Birth-year seed 2, y ~ x + d1 + d2 + d3 + (1 + x | stratum): fit_maihda()'s binomial fit, at a
# deviance of 1896.39238, 0.842 above the maximum. lme4's entry for the slope is 1.4e-5, the
# block's scaled entries are 0.073, -0.153 and 0.063.
year_2_adjusted <- list(theta = c(18.87465413, -0.009548702190, 1.401297664e-05),
                        fixef = c(-34.4610912772801, 0.0173054942321, -0.0607937976128,
                                  0.5292727647418, 0.3367018978494, 0.2327724356633,
                                  -0.3746221751159))

# Twenty-four additive strata (2 x 3 x 4 dimensions, effects SD 0.3) of 50 rows with a count over
# an exposure E between 0.5 and 3, a stratum intercept deviation of SD 0.25 and a slope deviation
# of SD 0.006 per year of birth (1935 to 2004, deviations centred at 1970).
poisson_year_data <- function(seed) {
  with_seed(seed, {
    lev <- expand.grid(a = 1:2, b = 1:3, c = 1:4)
    st <- rep(seq_len(24), each = 50)
    d <- lev[st, ]
    d$stratum <- factor(st)
    for (nm in c("a", "b", "c")) d[[nm]] <- factor(d[[nm]])
    rownames(d) <- NULL
    add <- stats::rnorm(2, 0, 0.3)[lev$a] + stats::rnorm(3, 0, 0.3)[lev$b] +
      stats::rnorm(4, 0, 0.3)[lev$c]
    u <- stats::rnorm(24, 0, 0.25)
    v <- stats::rnorm(24, 0, 0.006)
    d$x <- sample(1935:2004, 1200, replace = TRUE)
    d$E <- stats::runif(1200, 0.5, 3)
    d$y <- stats::rpois(1200, exp(-0.6 + add[st] + u[st] + (0.012 + v[st]) * (d$x - 1970) + log(d$E)))
    d
  })
}
# Seed 6, y ~ x + a + b + c + (1 + x | stratum) with offset = log(E): fit_maihda()'s fit, the
# slope's entry at 0, at a deviance of 3321.2326, 0.234 above the maximum.
poisson_6 <- list(theta = c(5.418014162233, -0.002854870055, 0),
                  fixef = c(-21.6881484085672, 0.0109862555141, -0.0860131066162, -0.2327694405517,
                            -0.4185356242199, -0.2909808214500, 0.1498404366988, 0.0265816254644))

# An lmer fit placed at `theta`: an optimiser that returns its start, and no move of an entry
# within lme4's boundary.tol of zero onto zero.
stay <- function(fn, par, lower, upper, control) {
  list(par = par, fval = fn(par), conv = 0, feval = 1, message = "stay")
}
placed <- function(theta, data, formula = adjusted) {
  quiet(lme4::lmer(formula, data = data, REML = TRUE, start = theta,
                   control = lme4::lmerControl(optimizer = stay, boundary.tol = 0)))
}

# Seed 2: the REML optimum, from lmer's deviance function searched with no bound on the Cholesky
# diagonal (and in each block's centred and scaled terms), is 2058.639975 at this theta. lme4's
# entry for the slope there is 3.1e-5; with the year centred at its mean and divided by its SD
# the block's factor is (0.080468, 0.066857, 0.077714).
optimum_2 <- c(11.666161, -0.0059493842, 3.1108755e-05)
scaled_2 <- c(0.080467934, 0.066857424, 0.077713621)

# The theta of a single (Intercept), x block for a factor given in its centred and scaled terms.
from_scaled <- function(phi, model) {
  Tm <- maihda_lmer_block_scales(lme4::getME(model, "mmList"), lme4::getME(model, "cnms"))[[1]]
  L <- matrix(0, 2, 2)
  L[lower.tri(L, diag = TRUE)] <- phi
  B <- maihda_reexpress_block_factor(L, Tm)
  B[lower.tri(B, diag = TRUE)]
}

test_that("a block's factor is carried to other terms without forming its covariance", {
  L <- matrix(c(2, -0.5, 0.3, 0, 1.5, -0.2, 0, 0, 0.7), 3)
  M <- matrix(c(1, 0, 0, -3, 0.5, 0, 2, 0, 0.25), 3)
  B <- maihda_reexpress_block_factor(L, M)
  expect_equal(tcrossprod(B), M %*% tcrossprod(L) %*% t(M), tolerance = 1e-12)
  expect_identical(B[upper.tri(B)], c(0, 0, 0))
  expect_true(all(diag(B) >= 0))
  # Of rank one in a birth year's units, the intercept at the mean 1.5e-4: the slope's entry
  # stays at zero.
  T2 <- matrix(c(1, 0, -1969.67 / 17.6, 1 / 17.6), 2)
  R1 <- matrix(c(1.5e-4 + 1969.67 * 0.012, -0.012, 0, 0), 2)
  B1 <- maihda_reexpress_block_factor(R1, solve(T2))
  expect_identical(B1[2, 2], 0)
  expect_equal(B1[1, 1], 1.5e-4, tolerance = 1e-8)
})

test_that("lme4's test fires on a birth-year block at an interior REML optimum", {
  at <- placed(optimum_2, birth_year_data(2))
  expect_lt(abs(maihda_lmer_criterion(at) - 2058.639975), 1e-5)
  expect_lt(maihda_lmer_stationary_distance(at), 0.001)
  expect_true(lme4::isSingular(at))
  expect_identical(maihda_singular_terms_lme4(at), "stratum")
  cnms <- lme4::getME(at, "cnms")
  sc <- maihda_lmer_block_scales(lme4::getME(at, "mmList"), cnms)
  expect_equal(maihda_lmer_transform_theta(optimum_2, cnms, lapply(sc, solve)), scaled_2,
               tolerance = 1e-6)
  expect_identical(maihda_re_boundary_components_scaled_lme4(at),
                   list(stratum = c("(Intercept)" = FALSE, x = FALSE)))
})

test_that("the record calls that fit not singular and print() says what lme4 reported", {
  at <- placed(optimum_2, birth_year_data(2))
  expect_true(any(grepl("boundary (singular)", at@optinfo$conv$lme4$messages, fixed = TRUE)))
  dg <- maihda_fit_diagnostics(at)
  expect_false(dg$singular)
  expect_identical(dg$singular_terms, character(0))
  expect_true(is.na(dg$stratum_singular))
  expect_true(dg$lme4_singular)
  # lme4's message is not a convergence warning.
  expect_true(dg$converged)
  expect_identical(dg$messages, character(0))
  lines <- maihda_format_fit_diagnostics(dg)
  expect_identical(lines, c(
    "lme4 reports a singular fit, taking its test in the units of the formula; with",
    "  each random-effects block's terms scaled, and centred where the block has an",
    "  intercept, no block is at the boundary."))
  expect_identical(maihda_app_diagnostics_items(lines), list(list(
    text = paste("lme4 reports a singular fit, taking its test in the units of the formula;",
                 "with each random-effects block's terms scaled, and centred where the block",
                 "has an intercept, no block is at the boundary."),
    sub = character(0))))
})

test_that("the line in the centred and scaled terms is 0.01, and a block of rank one is singular", {
  skip_on_cran()
  d <- birth_year_data(2)
  at <- placed(optimum_2, d)
  # The slope's scaled entry at 0.012, then at 0.008, 2.5e-4 and 5e-5 (lme4's entry from 7.3e-6
  # down), the rest of the block as at the optimum.
  above <- placed(from_scaled(replace(scaled_2, 3, 0.012), at), d)
  below <- lapply(c(0.008, 2.5e-4, 5e-5), function(s) placed(from_scaled(replace(scaled_2, 3, s), at), d))
  rank_one <- placed(replace(optimum_2, 3, 0), d)
  for (m in c(list(above, rank_one), below)) expect_true(lme4::isSingular(m))
  dg <- maihda_fit_diagnostics(above)
  expect_false(dg$singular)
  expect_true(dg$lme4_singular)
  for (m in c(below, list(rank_one))) {
    dg <- maihda_fit_diagnostics(m)
    expect_true(dg$singular)
    expect_identical(dg$singular_terms, "stratum")
    expect_true(dg$stratum_singular)
    expect_null(dg$lme4_singular)
    lines <- maihda_format_fit_diagnostics(dg)
    expect_true(any(grepl("Singular fit: the random-effects block for 'stratum'", lines, fixed = TRUE)))
    expect_false(any(grepl("lme4 reports a singular fit", lines, fixed = TRUE)))
  }
})

test_that("a block of rank one is singular in its scaled terms, at the intercept's entry too", {
  skip_on_cran()
  d <- birth_year_data(2)
  at <- placed(optimum_2, d)
  xbar <- mean(d$x)
  # Of rank one in the formula's units, with the intercept at the mean at 0.05 and a slope SD of
  # 0.01 or 0.012 per year: the slope's entry alone.
  for (l21 in c(-0.01, -0.012)) {
    m <- placed(c(0.05 - xbar * l21, l21, 0), d)
    expect_true(lme4::isSingular(m))
    expect_identical(maihda_re_boundary_components_scaled_lme4(m),
                     list(stratum = c("(Intercept)" = FALSE, x = TRUE)))
    dg <- maihda_fit_diagnostics(m)
    expect_true(dg$singular)
    expect_identical(dg$singular_terms, "stratum")
    expect_null(dg$lme4_singular)
  }
  # The intercept at the mean at 0.008, the slope as at the optimum: the intercept's entry alone.
  m <- placed(from_scaled(replace(scaled_2, 1, 0.008), at), d)
  expect_identical(maihda_re_boundary_components_scaled_lme4(m),
                   list(stratum = c("(Intercept)" = TRUE, x = FALSE)))
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$singular)
  expect_true(dg$stratum_singular)
})

test_that("a block lme4's test misses in large units is singular in its scaled terms", {
  # The year in thousands of years from 1970: lme4's entry for the slope is 0.0029, its entry
  # per standard deviation of the year 5.0e-5.
  d <- birth_year_data(2)
  d$x <- (d$x - 1970) / 1000
  at <- placed(c(0.3, 0.5, 0.0029), d)
  expect_false(lme4::isSingular(at))
  dg <- maihda_fit_diagnostics(at)
  expect_true(dg$singular)
  expect_identical(dg$singular_terms, "stratum")
  expect_true(dg$stratum_singular)
  expect_false(dg$lme4_singular)
  lines <- maihda_format_fit_diagnostics(dg)
  expect_true(any(grepl("Singular fit: the random-effects block for 'stratum'", lines, fixed = TRUE)))
  expect_false(any(grepl("lme4 reports a singular fit", lines, fixed = TRUE)))
})

test_that("a block without an intercept is scaled by its root mean square, each block on its own", {
  skip_on_cran()
  d <- birth_year_data(2)
  # (1 | stratum) + (0 + x | stratum): lme4's entry for the slope is 3.8e-5 per year, 0.075
  # times the root mean square of the year.
  slope_only <- y ~ x + d1 + d2 + d3 + (1 | stratum) + (0 + x | stratum)
  at <- placed(c(0.08, 3.8e-5), d, slope_only)
  expect_true(lme4::isSingular(at))
  cnms <- lme4::getME(at, "cnms")
  sc <- maihda_lmer_block_scales(lme4::getME(at, "mmList"), cnms)
  expect_equal(maihda_lmer_transform_theta(c(0.08, 3.8e-5), cnms, lapply(sc, solve)),
               c(0.08, 3.8e-5 * sqrt(mean(d$x^2))), tolerance = 1e-10)
  expect_identical(lapply(maihda_re_boundary_components_scaled_lme4(at), unname),
                   list(stratum = FALSE, stratum = FALSE))
  dg <- maihda_fit_diagnostics(at)
  expect_false(dg$singular)
  expect_true(dg$lme4_singular)
  # The intercept's block at 0.005 keeps lme4's line.
  at <- placed(c(0.005, 3.8e-5), d, slope_only)
  expect_identical(lapply(maihda_re_boundary_components_scaled_lme4(at), unname),
                   list(stratum = FALSE, stratum = FALSE))
  expect_false(maihda_fit_diagnostics(at)$singular)
  # (1 + x || stratum) with the intercept's block at zero: that block is singular either way.
  uncorrelated <- y ~ x + d1 + d2 + d3 + (1 + x || stratum)
  at <- placed(c(0, 3.8e-5), d, uncorrelated)
  expect_identical(lapply(maihda_re_boundary_components_scaled_lme4(at), unname),
                   list(stratum = TRUE, stratum = FALSE))
  dg <- maihda_fit_diagnostics(at)
  expect_true(dg$singular)
  expect_identical(dg$singular_terms, "stratum")
  expect_null(dg$lme4_singular)
})

test_that("a context block keeps lme4's entries beside a stratum block in its scaled terms", {
  # Six contexts with no effect of their own: the REML optimum has the context SD at zero and
  # the stratum block of seed 2's optimum.
  d <- birth_year_data(2)
  d$ctx <- with_seed(2, factor(sample(1:6, nrow(d), replace = TRUE)))
  context_formula <- y ~ x + d1 + d2 + d3 + (1 + x | stratum) + (1 | ctx)
  at <- placed(c(optimum_2, 0), d, context_formula)
  expect_lt(abs(maihda_lmer_criterion(at) - 2058.639975), 1e-5)
  expect_identical(maihda_singular_terms_lme4(at), c("stratum", "ctx"))
  dg <- maihda_fit_diagnostics(at)
  expect_true(dg$singular)
  expect_identical(dg$singular_terms, "ctx")
  expect_false(dg$stratum_singular)
  expect_null(dg$lme4_singular)
  lines <- maihda_format_fit_diagnostics(dg)
  expect_true(any(grepl("the random-effects block for 'ctx' is at", lines, fixed = TRUE)))
  expect_true(any(grepl("The 'stratum' block is NOT at the boundary", lines, fixed = TRUE)))
})

test_that("intercept-only fits, glmer fits with nAGQ 0 and growth fits keep lme4's test", {
  skip_on_cran()
  same <- data.frame(g = factor(rep(1:10, each = 10)), y = rep(seq(-1, 1, length.out = 10), times = 10))
  m <- quiet(lme4::lmer(y ~ 1 + (1 | g), data = same))
  expect_true(lme4::isSingular(m))
  expect_null(maihda_re_boundary_components_scaled_lme4(m))
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$singular)
  expect_identical(dg$singular_terms, "g")
  expect_null(dg$lme4_singular)
  # Seed 2's fit read as a growth fit: lme4's verdict.
  at <- placed(optimum_2, birth_year_data(2))
  dg <- maihda_fit_diagnostics(at, longitudinal_info = list(id = "id", time = "x"))
  expect_true(dg$singular)
  expect_identical(dg$singular_terms, "stratum")
  expect_null(dg$lme4_singular)
  # A glmer fit with nAGQ 0 at seed 2's adjusted theta, lme4's entry for the slope at 1.4e-5 and
  # the block's scaled entries 0.06 or more: lme4's verdict, and no search.
  g0 <- placed_glmer(year_2_adjusted$theta, NULL, binary_data(2), formula = adjusted, nAGQ = 0L)
  expect_true(lme4::isSingular(g0))
  expect_false(maihda_glmer_quadrature(g0))
  expect_true(is.na(maihda_glmer_search_gain(g0)))
  dg <- maihda_fit_diagnostics(g0)
  expect_true(dg$singular)
  expect_identical(dg$singular_terms, "stratum")
  expect_null(dg$lme4_singular)
  expect_null(dg$search_gain)
})

test_that("a singular glmer fit short of the maximum is found by the search", {
  skip_on_cran()
  g <- placed_glmer(year_7_short$theta, year_7_short$fixef, binary_data(7))
  expect_lt(abs(-2 * as.numeric(stats::logLik(g)) - 1810.3243), 1e-3)
  expect_true(lme4::isSingular(g))
  expect_true(maihda_glmer_quadrature(g))
  gain <- maihda_glmer_search_gain(g)
  expect_gte(gain, 0.01)
  # The search stops once it is `within` below the fit, short of the maximum.
  expect_lt(gain, 3)
  expect_gte(maihda_glmer_search_gain(g, within = 0.5), 0.5)
  dg <- maihda_fit_diagnostics(g)
  # At the boundary in its scaled terms too.
  expect_true(dg$singular)
  expect_identical(dg$singular_terms, "stratum")
  expect_null(dg$lme4_singular)
  expect_false(dg$converged)
  expect_equal(dg$search_gain, gain, tolerance = 1e-6)
  expect_identical(dg$messages, maihda_glmer_search_message(gain))
  lines <- maihda_format_fit_diagnostics(dg)
  expect_true(any(grepl("Singular fit: the random-effects block for 'stratum'", lines, fixed = TRUE)))
  expect_true(any(lines == "Convergence problems of the lme4 fit:"))
  expect_true(any(grepl("the fit stopped short of the maximum.", lines, fixed = TRUE)))
  # Read as a growth fit, the same point keeps lme4's verdict and is not searched.
  dg <- maihda_fit_diagnostics(g, longitudinal_info = list(id = "id", time = "x"))
  expect_true(dg$singular)
  expect_null(dg$lme4_singular)
  expect_null(dg$search_gain)
  # A fit its optimiser did not vouch for is not searched either.
  failed <- function(par, fn, lower, upper, control = list(), ...) {
    list(par = par, fval = fn(par), conv = 1L, message = "stopped")
  }
  g1 <- quiet(lme4::glmer(y ~ x + (1 + x | stratum), data = binary_data(7), family = stats::binomial,
                          start = list(theta = year_7_short$theta, fixef = year_7_short$fixef),
                          control = lme4::glmerControl(nAGQ0initStep = FALSE, optimizer = failed,
                                                       boundary.tol = 0)))
  dg <- maihda_fit_diagnostics(g1)
  expect_false(dg$converged)
  expect_null(dg$search_gain)
})

test_that("a glmer fit lme4 calls singular is searched when its scaled terms are off the boundary", {
  skip_on_cran()
  g <- placed_glmer(year_2_adjusted$theta, year_2_adjusted$fixef, binary_data(2), formula = adjusted)
  expect_lt(abs(-2 * as.numeric(stats::logLik(g)) - 1896.39238), 1e-4)
  expect_true(lme4::isSingular(g))
  gain <- maihda_glmer_search_gain(g)
  expect_gte(gain, 0.01)
  expect_lte(gain, 0.843)
  dg <- maihda_fit_diagnostics(g)
  # Without the search it would read clean.
  expect_false(dg$singular)
  expect_true(dg$lme4_singular)
  expect_false(dg$converged)
  expect_equal(dg$search_gain, gain, tolerance = 1e-6)
  expect_identical(dg$messages, maihda_glmer_search_message(gain))
  lines <- maihda_format_fit_diagnostics(dg)
  expect_true(any(grepl("lme4 reports a singular fit", lines, fixed = TRUE)))
  expect_false(any(grepl("Singular fit:", lines, fixed = TRUE)))
  expect_true(any(lines == "Convergence problems of the lme4 fit:"))
  expect_true(any(grepl("the fit stopped short of the maximum.", lines, fixed = TRUE)))
})

test_that("a Poisson fit with an offset is searched with its fixed effects centred and scaled", {
  skip_on_cran()
  d <- poisson_year_data(6)
  g <- quiet(lme4::glmer(y ~ x + a + b + c + (1 + x | stratum), data = d, family = stats::poisson,
                         offset = log(E), start = list(theta = poisson_6$theta, fixef = poisson_6$fixef),
                         control = lme4::glmerControl(nAGQ0initStep = FALSE, optimizer = keep,
                                                      boundary.tol = 0)))
  expect_lt(abs(-2 * as.numeric(stats::logLik(g)) - 3321.2326), 1e-4)
  expect_true(lme4::isSingular(g))
  offset <- g@frame[["(offset)"]] + 0
  effects <- lme4::ranef(g)
  # With the fixed effects left in the formula's units, the same search falls 0.0006.
  gain <- maihda_glmer_search_gain(g)
  expect_gte(gain, 0.01)
  expect_lte(gain, 0.235)
  # The search leaves the fit's offset and random effects as they were.
  expect_identical(g@frame[["(offset)"]], offset)
  expect_identical(lme4::ranef(g), effects)
  dg <- maihda_fit_diagnostics(g)
  expect_true(dg$singular)
  expect_false(dg$converged)
  expect_equal(dg$search_gain, gain, tolerance = 1e-6)
})

test_that("a glmer fit 0.02 short at a block of rank one is found, the search stopped 0.01 below it", {
  skip_on_cran()
  # Birth-year seed 9, fit_maihda()'s fit: the slope's entry at 0, 0.0207 above the maximum.
  g <- placed_glmer(c(23.457235047, -0.0121062044, 0), c(-36.473532157, 0.018315511174),
                    binary_data(9))
  gain <- maihda_glmer_search_gain(g)
  expect_gte(gain, 0.01)
  expect_lt(gain, 0.0208)
  dg <- maihda_fit_diagnostics(g)
  expect_true(dg$singular)
  expect_null(dg$lme4_singular)
  expect_false(dg$converged)
  expect_equal(dg$search_gain, gain, tolerance = 1e-6)
})

test_that("a singular glmer fit at the maximum is not singular in its scaled terms nor found short", {
  skip_on_cran()
  g <- placed_glmer(currency_3_max$theta, currency_3_max$fixef,
                    binary_data(3, n = 200, x_type = "currency", slope_sd = 1e-5, beta_x = 2e-5))
  expect_true(lme4::isSingular(g))
  gain <- maihda_glmer_search_gain(g)
  expect_gte(gain, 0)
  expect_lt(gain, 0.01)
  dg <- maihda_fit_diagnostics(g)
  expect_false(dg$singular)
  expect_true(dg$lme4_singular)
  expect_true(dg$converged)
  expect_null(dg$search_gain)
  expect_identical(dg$messages, character(0))
})

test_that("a glmer fit lme4 does not call singular is left to lme4's checks", {
  skip_on_cran()
  # Seed 7's fit with the slope's entry at 2e-4, 0.20 or more short of the maximum, placed with
  # lme4's checks switched off.
  g <- quiet(lme4::glmer(y ~ x + (1 + x | stratum), data = binary_data(7), family = stats::binomial,
                         start = list(theta = replace(year_7_short$theta, 3, 2e-4),
                                      fixef = year_7_short$fixef),
                         control = lme4::glmerControl(nAGQ0initStep = FALSE, optimizer = keep,
                                                      boundary.tol = 0, check.conv.grad = "ignore",
                                                      check.conv.hess = "ignore")))
  expect_false(lme4::isSingular(g))
  expect_gte(maihda_glmer_search_gain(g), 0.01)
  dg <- maihda_fit_diagnostics(g)
  expect_false(dg$singular)
  expect_true(dg$converged)
  expect_null(dg$search_gain)
})

test_that("fits with a diag() or cs() block keep lme4's test and are not searched", {
  skip_on_cran()
  d <- birth_year_data(2)
  expect_true(maihda_lme4_full_factors(placed(optimum_2, d)))
  # A cs() block's theta is its Cholesky factor, but its bounds are those of its standard
  # deviations and correlation, here 0.06, 1e-5 per year and 0.3; a diag() block's are its
  # standard deviations, 0.07 and 0.
  cs_block <- placed(c(0.06, 1e-5, 0.3), d, y ~ x + d1 + d2 + d3 + cs(1 + x | stratum))
  diag_block <- placed(c(0.07, 0), d, y ~ x + d1 + d2 + d3 + diag(1 + x | stratum))
  for (m in list(cs_block, diag_block)) {
    expect_true(lme4::isSingular(m))
    expect_false(maihda_lme4_full_factors(m))
    expect_null(maihda_re_boundary_components_scaled_lme4(m))
    dg <- maihda_fit_diagnostics(m)
    expect_true(dg$singular)
    expect_null(dg$lme4_singular)
  }
})

test_that("a glmer fit with a diag() block is not searched", {
  # Standard deviations of 1e-5 and 3e-4 per year.
  g <- placed_glmer(c(1e-5, 3e-4), year_7_short$fixef, binary_data(7),
                    formula = y ~ x + diag(1 + x | stratum))
  expect_true(lme4::isSingular(g))
  expect_false(maihda_lme4_full_factors(g))
  expect_true(is.na(maihda_glmer_search_gain(g)))
  dg <- maihda_fit_diagnostics(g)
  expect_true(dg$singular)
  expect_null(dg$lme4_singular)
  expect_null(dg$search_gain)
})

test_that("a glmer fit with intercepts only is not searched", {
  skip_on_cran()
  # The stratum variance at zero, 97 deviance short of its maximum (theta 0.658).
  d <- binary_data(7)
  g <- placed_glmer(0, c(-36.674537775, 0.0183899957887), d, formula = y ~ x + (1 | stratum))
  expect_true(lme4::isSingular(g))
  dg <- maihda_fit_diagnostics(g)
  expect_true(dg$singular)
  expect_true(dg$converged)
  expect_null(dg$search_gain)
})

test_that("the glmer search's sentence and its re-expression of theta", {
  expect_identical(maihda_glmer_search_message(0.02367),
                   paste("lme4 runs no convergence check on a singular fit, and a search of the same",
                         "model with its random-effect terms and fixed effects centred and scaled",
                         "reaches a deviance 0.024 lower: the fit stopped short of the maximum."))
  cnms <- list(stratum = c("(Intercept)", "x"), ctx = "(Intercept)")
  theta <- c(1.2, -0.3, 0.4, 0.25)
  M <- list(matrix(c(1, 0, 3, 2), 2), matrix(2, 1, 1))
  there <- maihda_reexpress_theta(theta, cnms, M)
  expect_equal(there[4], 0.5, tolerance = 1e-12)
  L <- matrix(c(1.2, -0.3, 0, 0.4), 2)
  B <- matrix(c(there[1], there[2], 0, there[3]), 2)
  expect_equal(tcrossprod(B), M[[1]] %*% tcrossprod(L) %*% t(M[[1]]), tolerance = 1e-12)
  expect_equal(maihda_reexpress_theta(there, cnms, lapply(M, solve)), theta, tolerance = 1e-12)
  # Not a glmer fit with its fixed effects among the optimised parameters.
  expect_false(maihda_glmer_quadrature(placed(optimum_2, birth_year_data(2))))
  expect_true(is.na(maihda_glmer_search_gain(placed(optimum_2, birth_year_data(2)))))
})

test_that("fit_maihda()'s fit of seed 2 reaches the optimum and prints no singular block", {
  skip_on_cran()
  m <- quiet(fit_maihda(adjusted, data = birth_year_data(2)))
  expect_lt(abs(maihda_lmer_criterion(m$model) - 2058.639975), 1e-4)
  expect_false(m$diagnostics$singular)
  expect_true(m$diagnostics$lme4_singular)
  expect_true(m$diagnostics$converged)
  out <- utils::capture.output(print(m))
  expect_false(any(grepl("Singular fit", out, fixed = TRUE)))
  expect_true(any(grepl("lme4 reports a singular fit", out, fixed = TRUE)))
})
