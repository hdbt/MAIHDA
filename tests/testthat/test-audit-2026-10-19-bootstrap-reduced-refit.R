# Audit pass 2026-10-19: the null-restricted refits of summary(df_method = "bootstrap") on
# lme4 growth fits.
#
# For each fixed-effect coefficient the bootstrap refits the model with that coefficient
# constrained to zero (maihda_refit_reduced(), else maihda_restrict_fixef()) and simulates
# its null distribution from that fit. For a growth fit the refit started from lme4's own
# start -- or from the start (and, for glmer, the control) update() carried over from the
# fit kept -- without the refits fit_maihda() gives the model: on a balanced 16-stratum
# design 14 of 110 null refits of quadratic growth fits stopped more than 0.01 above the
# lowest criterion found (up to 1.95 REML units; linear growth 2 of 60). Given the model's
# longitudinal tag, both fitters now refit a growth model as fit_maihda() fits it, from
# lme4's start and with the user's own control, and the formula route keeps the lower of
# that fit and the former refit.
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

# A balanced growth design: 16 factorial strata of 6 persons, waves 0..5, with stratum and
# person intercepts, slopes and (deg = 2) curvatures and a Gaussian outcome y.
growth_design <- function(seed, deg = 2) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
    add0 <- 0.3 * (g$A - 0.5) + 0.2 * (g$B - 0.5) - 0.25 * (g$C - 0.5) + 0.15 * (g$D - 0.5)
    add1 <- 0.04 * (g$A - 0.5) - 0.03 * (g$B - 0.5) + 0.05 * (g$C - 0.5) - 0.02 * (g$D - 0.5)
    u <- cbind(stats::rnorm(16, 0, 0.2), stats::rnorm(16, 0, 0.06),
               stats::rnorm(16, 0, if (deg >= 2) 0.02 else 0))
    p <- cbind(stats::rnorm(96, 0, 0.5), stats::rnorm(96, 0, 0.15),
               stats::rnorm(96, 0, if (deg >= 2) 0.06 else 0))
    id <- rep(seq_len(96), each = 6)
    st <- rep(seq_len(16), each = 6)[id]
    t <- rep(0:5, times = 96)
    eta <- add0[st] + u[st, 1] + (add1[st] + u[st, 2]) * t + u[st, 3] * t^2 +
      p[id, 1] + p[id, 2] * t + p[id, 3] * t^2
    data.frame(id = id, wave = t,
               A = factor(g$A[st]), B = factor(g$B[st]), C = factor(g$C[st]), D = factor(g$D[st]),
               y = 2 + 0.25 * t - 0.02 * t^2 + eta + stats::rnorm(576, 0, 0.5))
  })
}

adjusted_quadratic <- y ~ wave + A + B + C + D + wave:A + wave:B + wave:C + wave:D +
  I(wave^2):A + I(wave^2):B + I(wave^2):C + I(wave^2):D + (1 | A:B:C:D)

fit_quadratic <- function(seed, formula = adjusted_quadratic, ...) {
  quiet(fit_maihda(formula, data = growth_design(seed), id = "id", time = "wave",
                   time_degree = 2, ...))
}

# The null-restricted fit maihda_bootstrap_fixef() builds for the one-column term `term` of
# the fit_maihda() object `fit`, given its longitudinal tag.
reduced_fit <- function(fit, term) {
  model <- fit$model
  X <- lme4::getME(model, "X")
  labels <- attr(stats::terms(stats::formula(model, fixed.only = TRUE)), "term.labels")
  j <- which(attr(X, "assign") == match(term, labels))
  red <- maihda_refit_reduced(model, stats::update(stats::formula(model), paste(". ~ . -", term)),
                              fit$longitudinal_info)
  if (!is.null(red) &&
      !maihda_same_column_space(lme4::getME(red, "X"), X[, -j, drop = FALSE])) {
    red <- NULL
  }
  if (is.null(red)) red <- maihda_restrict_fixef(model, j, fit$longitudinal_info)
  red
}

spec_of <- function(fit) maihda_lme4_growth_spec(fit$model, fit$longitudinal_info)
spec_core <- function(spec) spec[c("id", "time", "time_degree", "stratum_slope")]

# A cross-sectional fit with random slopes on wave for the strata and for the persons:
# random-effect terms of a growth model's form, in a fit made without id and time.
cross_sectional_slopes <- function() {
  d <- growth_design(2, deg = 1)
  d$stratum <- interaction(d$A, d$B, d$C, d$D)
  quiet(fit_maihda(y ~ wave + A + (1 + wave | stratum) + (1 + wave | id), data = d))
}

# ---- the pieces ----------------------------------------------------------------------------

test_that("a growth fit's structure is read off the tag fit_maihda() put on it", {
  skip_on_cran()
  q <- fit_quadratic(2, y ~ wave + (1 | A:B:C:D))
  expect_identical(spec_core(spec_of(q)),
                   list(id = "id", time = "wave", time_degree = 2L, stratum_slope = TRUE))
  expect_identical(spec_of(q)$longitudinal, q$longitudinal_info)
  ns <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(2, deg = 1), id = "id",
                         time = "wave", stratum_slope = FALSE))
  expect_identical(spec_core(spec_of(ns)),
                   list(id = "id", time = "wave", time_degree = 1L, stratum_slope = FALSE))
  # On an axis that does not start at 0 the model's time is the internal column.
  d10 <- growth_design(2)
  d10$wave <- d10$wave + 10
  cq <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d10, id = "id", time = "wave",
                         time_degree = 2))
  expect_identical(spec_of(cq)$time, ".maihda_ctime")
  # A restart's numeric start in the call does not hide the structure.
  r1 <- fit_quadratic(1)
  expect_true(is.numeric(stats::getCall(r1$model)$start))
  expect_identical(spec_of(r1)$time_degree, 2L)
  # Not a growth fit, or fitted from the user's own start: NULL.
  w0 <- growth_design(2)
  xs <- quiet(fit_maihda(y ~ A + (1 | A:B:C:D), data = w0[w0$wave == 0, ]))
  expect_null(spec_of(xs))
  # A cross-sectional fit is not a growth fit, whatever the form of its random-effect terms.
  cs <- cross_sectional_slopes()
  expect_identical(sort(names(lme4::getME(cs$model, "cnms"))), c("id", "stratum"))
  expect_true(all(vapply(lme4::getME(cs$model, "cnms"), identical, logical(1),
                         c("(Intercept)", "wave"))))
  expect_null(cs$longitudinal_info)
  expect_null(spec_of(cs))
  us <- fit_quadratic(2, y ~ wave + (1 | A:B:C:D), start = rep(c(1, 0, 0, 1, 0, 1), 2))
  expect_true(is.name(stats::getCall(us$model)$start))
  expect_null(spec_of(us))
  # A time column whose name is not syntactic: the spec gives the column's own name.
  dt <- growth_design(2, deg = 1)
  names(dt)[names(dt) == "wave"] <- "time point"
  nt <- quiet(fit_maihda(y ~ A + (1 | A:B:C:D), data = dt, id = "id", time = "time point"))
  expect_identical(spec_of(nt)$time, "time point")
})

test_that("a growth fit's null refits go through fit_maihda()'s fitter, from lme4's start", {
  skip_on_cran()
  # Seed 1: the adjusted fit kept is a restart, so its call carries a numeric start.
  m <- fit_quadratic(1)
  expect_true(is.numeric(stats::getCall(m$model)$start))
  seen <- list()
  real <- maihda_fit_lmer_growth
  local_mocked_bindings(maihda_fit_lmer_growth = function(fit_call, fit_env, id, ...) {
    seen[[length(seen) + 1L]] <<- list(start = fit_call$start, data = fit_call$data, id = id,
                                       args = list(...))
    real(fit_call, fit_env, id, ...)
  })
  # The formula route (update()), without the start of the kept restart.
  red <- reduced_fit(m, "wave:A")
  expect_length(seen, 1L)
  expect_null(seen[[1]]$start)
  expect_identical(seen[[1]]$id, "id")
  expect_identical(seen[[1]]$args, list(time = "wave", time_degree = 2L, first = NULL))
  expect_s4_class(red, "lmerMod")
  # Dropping `wave` from the formula restricts nothing (R recodes wave:A and its siblings to
  # span the same space), so that fit is set aside and the design route
  # (maihda_restrict_fixef()) refits, with the data where the restarts find them.
  seen <- list()
  red <- reduced_fit(m, "wave")
  expect_length(seen, 2L)
  expect_true(is.name(seen[[2]]$data))
  expect_null(seen[[2]]$start)
  expect_identical(seen[[2]]$args, list(time = "wave", time_degree = 2L, first = NULL))
  expect_identical(ncol(lme4::getME(red, "X")), ncol(lme4::getME(m$model, "X")) - 1L)
  # A linear model with a stratum slope takes the restarts too; one without takes none,
  # as at fit time.
  d1 <- growth_design(2, deg = 1)
  lin <- quiet(fit_maihda(y ~ wave + A + (1 | A:B:C:D), data = d1, id = "id", time = "wave"))
  seen <- list()
  red <- maihda_refit_reduced(lin$model, stats::update(stats::formula(lin$model), ". ~ . - A"),
                              lin$longitudinal_info)
  expect_length(seen, 1L)
  expect_identical(seen[[1]]$args, list(time = "wave", time_degree = 1L, first = NULL))
  # Without the tag, the same fit's null is update()'s alone.
  seen <- list()
  red <- maihda_refit_reduced(lin$model, stats::update(stats::formula(lin$model), ". ~ . - A"))
  expect_length(seen, 0L)
  expect_identical(lme4::REMLcrit(red), lme4::REMLcrit(quiet(stats::update(
    lin$model, formula. = stats::update(stats::formula(lin$model), ". ~ . - A")))))
  flat <- quiet(fit_maihda(y ~ wave + A + (1 | A:B:C:D), data = d1, id = "id", time = "wave",
                           stratum_slope = FALSE))
  seen <- list()
  red <- maihda_refit_reduced(flat$model, stats::update(stats::formula(flat$model), ". ~ . - A"),
                              flat$longitudinal_info)
  expect_length(seen, 0L)
  expect_s4_class(red, "lmerMod")
  # A dropped covariate with missing values: update() refits on more rows, so the null is
  # refitted on the model's own frame -- through the restarts as well.
  d1$x <- with_seed(3, stats::rnorm(nrow(d1)))
  d1$x[c(5, 77, 300)] <- NA
  mx <- quiet(fit_maihda(y ~ wave + x + (1 | A:B:C:D), data = d1, id = "id", time = "wave"))
  seen <- list()
  red <- maihda_refit_reduced(mx$model, stats::update(stats::formula(mx$model), ". ~ . - x"),
                              mx$longitudinal_info)
  expect_length(seen, 2L)
  expect_true(is.name(seen[[2]]$data))
  expect_identical(stats::nobs(red), stats::nobs(mx$model))
})

test_that("the former refit is kept when it is the lower", {
  skip_on_cran()
  # Seed 1's kept fit is a restart, so the formula route also refits from its start. A
  # stand-in growth fitter that stops at lme4's starting values loses to it.
  m <- fit_quadratic(1)
  local_mocked_bindings(maihda_fit_lmer_growth = function(fit_call, fit_env, id, ...) {
    fit_call$control <- lme4::lmerControl(optimizer = NULL)
    eval(fit_call, fit_env)
  })
  red <- reduced_fit(m, "wave:A")
  expect_true(is.numeric(stats::getCall(red)$start))
  stuck <- quiet(stats::update(m$model, formula. = stats::formula(red), start = NULL,
                               control = lme4::lmerControl(optimizer = NULL)))
  expect_lt(lme4::REMLcrit(red), lme4::REMLcrit(stuck) - 1)
})

test_that("a glmer null refit does not keep the control of glmer's refit across a bound", {
  skip_on_cran()
  skip_slow_test()
  # A Poisson quadratic growth fit whose kept fit is glmer's refit across a variance bound
  # followed by its last pass, both made with glmer's second stage only (the call gives the
  # pass's control, which skips glmer's first stage). Refitted with the refit's
  # control and from lme4's start, the null for B1 stopped 80 log-likelihood units on the
  # -2 scale above the lowest found, 999.9206 (bobyqa, Nelder_Mead and nlminb from lme4's
  # start and a polish from the refit's estimates).
  d <- with_seed(791, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1)
    add0 <- as.vector(as.matrix(g) %*% c(0.3, 0.05, -0.2))
    add1 <- as.vector(as.matrix(g) %*% c(0.05, 0.005, -0.04))
    u <- cbind(stats::rnorm(8, 0, 0.2), stats::rnorm(8, 0, 0.06), stats::rnorm(8, 0, 0.02))
    p <- cbind(stats::rnorm(64, 0, 0.4), stats::rnorm(64, 0, 0.12), stats::rnorm(64, 0, 0.03))
    id <- rep(seq_len(64), each = 5)
    st <- rep(rep(seq_len(8), each = 8), each = 5)
    t <- rep(0:4, times = 64)
    eta <- add0[st] + u[st, 1] + (add1[st] + u[st, 2]) * t + u[st, 3] * t^2 +
      p[id, 1] + p[id, 2] * t + p[id, 3] * t^2
    y <- 1.5 + 0.2 * t - 0.015 * t^2 + eta + stats::rnorm(length(id), 0, 0.6)
    data.frame(id = id, wave = t, A = factor(g$A[st]), B = factor(g$B[st]), C = factor(g$C[st]),
               cnt = stats::rpois(length(id), exp(0.3 + 0.35 * (eta - mean(eta)) / stats::sd(eta) + 0.05 * t)),
               y = y)
  })
  m <- quiet(fit_maihda(cnt ~ wave + A + B + C + B:wave + (1 | A:B:C), data = d, id = "id",
                        time = "wave", time_degree = 2, family = "poisson"))
  expect_identical(stats::getCall(m$model)$control, quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  controls <- list()
  real <- maihda_fit_glmer_growth
  local_mocked_bindings(maihda_fit_glmer_growth = function(fit_call, fit_env, id, ...) {
    controls[[length(controls) + 1L]] <<- list(control = fit_call$control, start = fit_call$start)
    real(fit_call, fit_env, id, ...)
  })
  red <- maihda_refit_reduced(m$model, stats::update(stats::formula(m$model), ". ~ . - B"),
                              m$longitudinal_info)
  expect_length(controls, 1L)
  expect_null(controls[[1]]$control)
  expect_null(controls[[1]]$start)
  expect_lt(-2 * as.numeric(stats::logLik(red)) - 999.9206, 0.05)
  # With a control of the user's own, the null is refitted with that control.
  mu <- quiet(fit_maihda(cnt ~ wave + A + B + C + B:wave + (1 | A:B:C), data = d, id = "id",
                         time = "wave", time_degree = 2, family = "poisson",
                         control = lme4::glmerControl()))
  expect_identical(stats::getCall(mu$model)$control, quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  controls <- list()
  red <- maihda_refit_reduced(mu$model, stats::update(stats::formula(mu$model), ". ~ . - B"),
                              mu$longitudinal_info)
  expect_identical(controls[[1]]$control, as.name(".maihda_arg_control"))
})

test_that("other fits keep their refits as they were", {
  skip_on_cran()
  local_mocked_bindings(
    maihda_fit_lmer_growth = function(...) stop("restarts used"),
    maihda_fit_glmer_growth = function(...) stop("restarts used"))
  # A cross-sectional fit: update() from lme4's start, unchanged.
  w0 <- growth_design(2)
  xs <- quiet(fit_maihda(y ~ A + B + (1 | A:B:C:D), data = w0[w0$wave == 0, ]))
  red_form <- stats::update(stats::formula(xs$model), ". ~ . - A")
  red <- maihda_refit_reduced(xs$model, red_form, xs$longitudinal_info)
  ref <- quiet(stats::update(xs$model, formula. = red_form))
  expect_identical(lme4::getME(red, "theta"), lme4::getME(ref, "theta"))
  expect_identical(lme4::REMLcrit(red), lme4::REMLcrit(ref))
  # A cross-sectional fit with random slopes of a growth model's form: refitted as before,
  # by update() and on its design columns.
  cs <- cross_sectional_slopes()
  red_form <- stats::update(stats::formula(cs$model), ". ~ . - A")
  red <- maihda_refit_reduced(cs$model, red_form, cs$longitudinal_info)
  expect_identical(lme4::REMLcrit(red),
                   lme4::REMLcrit(quiet(stats::update(cs$model, formula. = red_form))))
  expect_s4_class(maihda_restrict_fixef(cs$model, 2L, cs$longitudinal_info), "lmerMod")
  # A growth fit made from the user's own start: refitted from it alone.
  us <- quiet(fit_maihda(y ~ wave + A + (1 | A:B:C:D), data = growth_design(2), id = "id",
                         time = "wave", time_degree = 2, start = rep(c(1, 0, 0, 1, 0, 1), 2)))
  red <- maihda_refit_reduced(us$model, stats::update(stats::formula(us$model), ". ~ . - A"),
                              us$longitudinal_info)
  expect_s4_class(red, "lmerMod")
  expect_true(is.name(stats::getCall(red)$start))
})

test_that("a glmer growth fit's null refits go through its own fitter", {
  skip_on_cran()
  skip_slow_test()
  data(maihda_long_data, envir = environment())
  d <- maihda_long_data[maihda_long_data$id %in% unique(maihda_long_data$id)[1:150], ]
  m <- quiet(fit_maihda(low_wellbeing ~ wave + gender + (1 | gender:ethnicity:education), data = d,
                        id = "id", time = "wave", family = "binomial"))
  calls <- 0L
  real <- maihda_fit_glmer_growth
  local_mocked_bindings(maihda_fit_glmer_growth = function(fit_call, fit_env, id, ...) {
    calls <<- calls + 1L
    real(fit_call, fit_env, id, ...)
  })
  labels <- attr(stats::terms(stats::formula(m$model, fixed.only = TRUE)), "term.labels")
  term <- labels[startsWith(labels, "gender")][1]
  red <- maihda_refit_reduced(m$model, stats::update(stats::formula(m$model), paste(". ~ . -", term)),
                              m$longitudinal_info)
  expect_identical(calls, 1L)
  expect_s4_class(red, "glmerMod")
})

# ---- the fits it changes -------------------------------------------------------------------

test_that("null refits stuck above their optimum now reach it", {
  skip_on_cran()
  skip_slow_test()
  # The lowest criteria found for these reduced models by bobyqa and nlminb from lme4's
  # start and by eight random bobyqa starts; the null refits stopped 1.54, 1.77 and 0.57
  # above them (seed 14 from the start of the restart its full fit kept).
  cases <- list(list(4, "I(wave^2)", 1305.405999), list(14, "A", 1403.331923),
                list(7, "wave:A", 1279.321657))
  for (cs in cases) {
    m <- fit_quadratic(cs[[1]])
    red <- reduced_fit(m, cs[[2]])
    expect_lt(abs(lme4::REMLcrit(red) - cs[[3]]), 5e-3)
  }
})
