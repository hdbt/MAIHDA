# Audit pass 2026-10-27: glmer fits whose lme4 convergence checks fire at the maximum, left
# unmeasured and recorded as not converged, and glmer fits stopped at a variance bound that the
# convergence record did not measure.
#
# When lme4's checks fire on a glmer fit whose optimiser reported success, the fit is measured in
# centred and scaled units (maihda_lme4_stationary_distance()) and counts as converged within 0.1
# standard errors of a stationary point. The measure takes its derivatives on the fit's deviance
# function rebuilt from a copy of the fit, which starts its penalised least squares from the same
# point at every evaluation; at lme4's relative tolerance of 1e-7 its value at the fit's estimates
# lies up to about 1e-4 relative above the converged deviance, and lme4's own -2 log-likelihood,
# from the last step of its fit, lies more than 1e-5 from that at times. The fit was measured only
# when the function missed lme4's -2 log-likelihood by at most 1e-5 relative: on simulated 0 / 1
# outcomes with strata spread on the logit scale with an SD of 2.5 and a birth year, 17 of 135
# fits were left unmeasured and recorded as not converged, all 0.0003 to 0.066 standard errors
# from the maximum. The boundary-trap record (maihda_lme4_theta_deviance()) took its falls on such
# a function with the same check: where another random-effect block is large, as persons with an
# intercept SD of 2.5, 7 of 26 simulated stops at a variance bound went unmeasured. Per the user,
# both functions now run to a relative tolerance of 1e-10, and a glmer fit is measured when the
# function at its theta matches the function at its own variance parameters (lme4's `par`; 1e-5
# relative) and that lies within 1e-3 (relative) of lme4's -2 log-likelihood. A block whose
# function takes other parameters than theta, as cs() and ar1() do, fails the first.
# lme4 2.1-0 and later take the Laplace deviance at the converged step, where the function at
# lme4's tolerance and lme4's own -2 log-likelihood agree to about 1e-8 relative: the premises
# below hold before 2.1-0, and the fixes hold with either.
#
# Random numbers are drawn only inside with_seed(), which puts the suite's stream back.

# Whether lme4 takes its Laplace deviance after the last step of its penalised least squares.
lagged_laplace <- utils::packageVersion("lme4") < "2.1-0"

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

# Twelve strata of 60 rows, their effects drawn with an SD of 2.5 on the logit scale around a
# baseline of -2.5, a birth year from 1940 to 2000 with a log-odds slope of 0.02, and a 0 / 1
# outcome. lme4's checks fire on the null model's fit with the year, at its maximum.
year_strata <- function(seed) {
  with_seed(seed, {
    u <- stats::rnorm(12, 0, 2.5)
    st <- rep(seq_len(12), each = 60)
    year <- sample(1940:2000, length(st), replace = TRUE)
    eta <- -2.5 + u[st] + 0.02 * (year - 1970)
    data.frame(y = stats::rbinom(length(st), 1, stats::plogis(eta)), year = year,
               stratum = factor(st))
  })
}

null_year <- y ~ year + (1 | stratum)

# glmer's own fits of those data, rounded to seven digits, each within 0.05 standard errors of its
# maximum. With lme4 before 2.1-0: for seed 6, lme4's -2 log-likelihood lies 2.0e-5 relative above
# the converged deviance; for seeds 31 and 40, the deviance function at lme4's tolerance lies
# 1.4e-5 and 1.5e-5 above lme4's -2 log-likelihood, and taken on it the distance is 0.0073 and
# 0.0039 where the converged function gives 0.0192 and 0.0124.
points <- list(
  "6" = list(theta = 2.852984, beta = c(-41.92978, 0.01993696)),
  "31" = list(theta = 2.442014, beta = c(-28.18353, 0.01320105)),
  "40" = list(theta = 3.37779, beta = c(-62.74365, 0.02941082))
)

# A glmer fit placed at a point: optimisers that return their start.
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}
stay_control <- function(...) {
  lme4::glmerControl(optimizer = list(stay, stay), nAGQ0initStep = FALSE, boundary.tol = 0, ...)
}
placed <- function(seed) {
  point <- points[[as.character(seed)]]
  quiet(lme4::glmer(null_year, data = year_strata(seed), family = stats::binomial,
                    start = list(theta = point$theta, fixef = point$beta), control = stay_control()))
}

# The fit's deviance function (from a detached copy), its penalised least squares run to
# `tolerance` (lme4's 1e-7 when NULL).
devfun_at <- function(m, tolerance = NULL) {
  f <- lme4::getME(maihda_lme4_detached_copy(m), "devfun")
  if (!is.null(tolerance)) environment(f)$tolPwrss <- tolerance
  f
}
estimates <- function(m) c(unname(lme4::getME(m, "theta")), unname(lme4::fixef(m)))
relative_gap <- function(a, b) abs(a - b) / abs(b)
stored_deviance <- function(m) -2 * as.numeric(stats::logLik(m))

test_that("a fit whose checks fire at the maximum is measured and recorded as converged", {
  skip_on_cran()
  for (seed in c(6, 40)) {
    m <- placed(seed)
    said <- as.character(m@optinfo$conv$lme4$messages)
    # The premise: lme4's checks fire, its optimiser reports success, and the deviance function
    # at lme4's tolerance misses lme4's -2 log-likelihood by more than 1e-5 relative.
    expect_true(any(grepl("failed to converge", said, fixed = TRUE)))
    expect_identical(as.integer(m@optinfo$conv$opt), 0L)
    if (lagged_laplace) {
      expect_gt(relative_gap(devfun_at(m)(estimates(m)), stored_deviance(m)), 1e-5)
    }
    distance <- maihda_lme4_stationary_distance(m)
    expect_true(is.finite(distance))
    expect_lt(distance, 0.05)
    dg <- maihda_fit_diagnostics(m)
    expect_true(dg$converged)
    expect_identical(dg$messages, character(0))
    expect_identical(dg$lme4_overruled, said)
    expect_identical(dg$stationary_distance, distance)
  }
})

test_that("the distance is taken on the deviance function run to 1e-10", {
  skip_on_cran()
  for (seed in c(31, 40)) {
    m <- placed(seed)
    # Taken on the function run to 1e-12 with the same checks.
    exact <- maihda_lme4_distance_with(m, devfun_at(m, 1e-12))
    expect_true(is.finite(exact))
    expect_lt(abs(maihda_lme4_stationary_distance(m) - exact), 1e-3)
    # At lme4's tolerance the function's error moves the distance by more.
    if (lagged_laplace) {
      expect_gt(abs(maihda_lme4_distance_with(m, devfun_at(m)) - exact), 0.005)
    }
  }
})

test_that("lme4's -2 log-likelihood is a guard, not the reference", {
  skip_on_cran()
  m <- placed(6)
  f <- devfun_at(m, 1e-10)
  # The premise: lme4's -2 log-likelihood lies more than 1e-5 relative from the converged
  # deviance at its own estimates.
  if (lagged_laplace) {
    expect_gt(relative_gap(stored_deviance(m), devfun_at(m, 1e-12)(estimates(m))), 1e-5)
  }
  expect_true(is.finite(maihda_lme4_distance_with(m, f)))
  # A function within 1e-3 relative of lme4's -2 log-likelihood is taken as the fit's; one further
  # away is not, as a function rebuilt without the fit's offset or weights is not.
  expect_true(is.finite(maihda_lme4_distance_with(m, function(p) f(p) * (1 + 5e-4))))
  expect_true(is.na(maihda_lme4_distance_with(m, function(p) f(p) * (1 + 2e-3))))
})

test_that("the deviance at the re-expressed estimates must match the function at the fit's", {
  skip_on_cran()
  m <- placed(40)
  f <- devfun_at(m, 1e-10)
  expect_true(is.finite(maihda_lme4_distance_with(m, f)))
  # A re-expression that does not carry the fit's estimates back to themselves.
  carried <- maihda_lmer_transform_theta
  local_mocked_bindings(maihda_lmer_transform_theta = function(theta, cnms, Ts, diagonal = NULL) {
    carried(theta, cnms, Ts, diagonal) * 1.01
  })
  expect_true(is.na(maihda_lme4_distance_with(m, f)))
})

test_that("theta is checked against the parameters the deviance function takes", {
  skip_on_cran()
  skip_if(utils::packageVersion("lme4") < "2.0.0")
  # Sixteen strata of 60 rows, effects SD 1.5, a slope on an N(0, 1) covariate varying with an SD
  # of 0.5, and a birth year.
  d <- with_seed(1101, {
    st <- rep(seq_len(16), each = 60)
    x <- stats::rnorm(length(st))
    u <- stats::rnorm(16, 0, 1.5)
    v <- stats::rnorm(16, 0, 0.5)
    year <- sample(1930:2005, length(st), replace = TRUE)
    eta <- -1 + u[st] + (0.3 + v[st]) * x + 0.02 * (year - 1967)
    data.frame(y = stats::rbinom(length(st), 1, stats::plogis(eta)), x = x, year = year,
               stratum = factor(st))
  })
  fit <- function(formula) quiet(lme4::glmer(formula, data = d, family = stats::binomial))
  # diag(hom = TRUE) holds one standard deviation, which is theta: measured, as theta as fitted.
  hom <- fit(y ~ year + x + diag(1 + x | stratum, hom = TRUE))
  expect_null(maihda_lme4_diagonal_blocks(lme4::getME(hom, "cnms"), lme4::getME(hom, "lower")))
  expect_equal(maihda_lme4_parameters(hom), unname(lme4::getME(hom, "theta")))
  expect_lt(maihda_lme4_stationary_distance(hom), 0.1)
  # cs() holds standard deviations and a correlation, which are not theta: at theta the function
  # misses the fit's -2 log-likelihood, and the fit is not measured.
  cs_fit <- fit(y ~ year + x + cs(1 + x | stratum))
  expect_false(isTRUE(all.equal(maihda_lme4_parameters(cs_fit),
                                unname(lme4::getME(cs_fit, "theta")))))
  expect_gt(relative_gap(devfun_at(cs_fit, 1e-10)(estimates(cs_fit)), stored_deviance(cs_fit)),
            1e-3)
  expect_true(is.na(maihda_lme4_stationary_distance(cs_fit)))
  # Told that its function takes other parameters than theta, the diag(hom = TRUE) fit is not
  # measured either.
  expect_false(is.null(maihda_lme4_theta_deviance(hom)))
  own <- maihda_lme4_parameters
  local_mocked_bindings(maihda_lme4_parameters = function(model) own(model) * 1.5)
  expect_true(is.na(maihda_lme4_stationary_distance(hom)))
  expect_null(maihda_lme4_theta_deviance(hom))
})

test_that("fit_maihda() records such a fit as converged and says how far it is", {
  skip_on_cran()
  # Given glmer's estimates as its start and optimisers that keep them, fit_maihda() makes no refit.
  point <- points[["6"]]
  m <- quiet(fit_maihda(null_year, data = year_strata(6), family = "binomial",
                        start = list(theta = point$theta, fixef = point$beta),
                        control = stay_control()))
  dg <- m$diagnostics
  expect_gt(length(dg$lme4_overruled), 0L)
  expect_true(dg$converged)
  expect_lt(dg$stationary_distance, 0.1)
  shown <- maihda_format_fit_diagnostics(dg)
  expect_identical(shown[1],
                   "lme4 warned that the model failed to converge; in centred and scaled units the")
  expect_false(any(grepl("Convergence warnings reported by lme4", shown, fixed = TRUE)))
})

# 500 persons with an intercept SD of 2.5 seen at five times, additive strata (3 x 4 dimensions)
# with a stratum slope on time (mean 0.2, SD 0.3), and a 0 / 1 outcome.
person_waves <- function() {
  with_seed(61001, {
    nid <- 500
    a <- sample(1:3, nid, TRUE)
    b <- sample(1:4, nid, TRUE)
    st <- (a - 1) * 4 + b
    ui <- stats::rnorm(nid, 0, 2.5)
    us <- stats::rnorm(12, 0, 0.1)
    vs <- stats::rnorm(12, 0, 0.3)
    id <- rep(seq_len(nid), each = 5)
    time <- rep(0:4, nid)
    eta <- -1 + 0.4 * a[id] - 0.3 * b[id] + ui[id] + us[st[id]] + (0.2 + vs[st[id]]) * time
    data.frame(y = stats::rbinom(length(id), 1, stats::plogis(eta)), time = time,
               a = factor(a[id]), b = factor(b[id]), stratum = factor(st[id]), id = factor(id))
  })
}
person_slope <- y ~ time + a + b + (1 + time | stratum) + (1 | id)

# glmer's own fit of the adjusted model, rounded to seven digits: the stratum block's intercept
# variance at zero beside a slope entry of 0.25, 0.56 below what the likelihood reaches across
# that bound.
trapped <- list(theta = c(2.199874, 0, 0.2524547, 0.2010752),
                beta = c(-1.0088, 0.1956842, 0.1755558, 1.013829, -0.2314996, -0.3379305,
                         -0.9104165))

test_that("a stop at a variance bound beside a large random-effect block is measured", {
  skip_on_cran()
  m <- quiet(lme4::glmer(person_slope, data = person_waves(), family = stats::binomial,
                         start = list(theta = trapped$theta, fixef = trapped$beta),
                         control = stay_control(calc.derivs = FALSE)))
  # The premise: a column at its bound, and the deviance function at lme4's tolerance missing
  # lme4's -2 log-likelihood by more than 1e-5 relative.
  expect_length(maihda_lme4_trap_columns(lme4::getME(m, "theta"), lme4::getME(m, "cnms")), 1L)
  if (lagged_laplace) {
    expect_gt(relative_gap(devfun_at(m)(estimates(m)), stored_deviance(m)), 1e-5)
  }
  expect_false(is.null(maihda_lme4_theta_deviance(m)))
  # The same falls taken on the function run to 1e-12 give 0.5578; at lme4's tolerance before
  # 2.1-0 they read 0.5651.
  traps <- maihda_lme4_boundary_traps(m)
  expect_lt(abs(traps$drop - 0.5578), 1e-3)
  # As a growth fit, which no search covers, it is recorded as not converged.
  dg <- maihda_fit_diagnostics(m, longitudinal_info = list(id = "id"))
  expect_false(dg$converged)
  expect_identical(nrow(dg$boundary_traps), 1L)
  # A stored -2 log-likelihood 1% from the deviance function's: a function that is not the fit's.
  shifted <- m
  shifted@devcomp$cmp[["dev"]] <- shifted@devcomp$cmp[["dev"]] * 1.01
  expect_null(maihda_lme4_theta_deviance(shifted))
})
