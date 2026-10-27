# Audit pass 2026-10-19: the refits of the bootstrap draws of lme4 growth fits.
#
# summary(bootstrap = TRUE) and summary(df_method = "bootstrap") refitted each simulated
# response once, from the fitted model's estimates (lme4::refit(); for glmer, Nelder-Mead on
# all the parameters), without the refits fit_maihda() gives a growth model. On a Poisson
# growth design of 16 strata of 6 persons, 12 to 19 of 100 draws per time unit stopped more
# than 0.01 short of their maximum (up to 1.9 in log-likelihood), none recorded as failed,
# and the VPC(t) band's upper limit ran up to 25 per cent low at early times; an lmer draw
# stopped on the wrong side of the stratum block's boundary. Each draw of a growth fit is now
# refitted from the fitted model's estimates and then, with that refit as the first fit, as
# fit_maihda() fits the model. (The null-restricted fits of df_method = "bootstrap" are
# tested in test-audit-2026-10-19-bootstrap-reduced-refit.R.)
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

# The balanced growth design of the earlier growth passes: 16 factorial strata of 6 persons,
# waves 0..5 times k, stratum and person intercepts, slopes and (deg = 2) curvatures, a
# Gaussian outcome y and a Poisson count n.
growth_design <- function(seed, k = 1, deg = 1) {
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
    data.frame(id = id, wave = t * k,
               A = factor(g$A[st]), B = factor(g$B[st]), C = factor(g$C[st]), D = factor(g$D[st]),
               y = 2 + 0.25 * t - 0.02 * t^2 + eta + stats::rnorm(576, 0, 0.5),
               n = stats::rpois(576, exp(0.5 + 0.1 * t + 0.6 * eta)))
  })
}

# A small count design: 4 strata of 6 persons, waves 0..3.
small_counts <- function(seed) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1)
    u <- stats::rnorm(4, 0, 0.3)
    p <- cbind(stats::rnorm(24, 0, 0.5), stats::rnorm(24, 0, 0.15))
    id <- rep(1:24, each = 4)
    st <- rep(1:4, each = 6)[id]
    t <- rep(0:3, 24)
    mu <- exp(1 + 0.1 * t + u[st] + p[id, 1] + p[id, 2] * t)
    data.frame(id = id, wave = t, A = factor(g$A[st]), B = factor(g$B[st]),
               expo = rep(c(1, 2), length.out = 96),
               n = stats::rnbinom(96, mu = mu, size = 3),
               s = pmin(stats::rpois(96, mu), 5L), tr = 5L)
  })
}

# ---- the replay of a draw ------------------------------------------------------------------

test_that("a draw is replayed through the fit's own call, on its own rows", {
  skip_on_cran()
  # Time in weeks (internally scaled), missing values in the response and in time, and
  # precision weights.
  d <- growth_design(1, 7)
  d$y[c(3, 50, 200)] <- NA
  d$wave[7] <- NA
  d$w <- rep(c(1, 2, 0.5), length.out = nrow(d))
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave", weights = w))
  y <- with_seed(5, maihda_simulate_lme4(m$model, nsim = 1))[[1]]
  rp <- maihda_growth_draw_call(m$model, y)
  expect_identical(rp$call$formula[[2]], as.name(".maihda_boot_y"))
  expect_null(rp$call$start)
  expect_identical(rp$call$weights, as.name(".maihda_arg_weights"))
  expect_true(exists(".maihda_arg_weights", envir = rp$env))
  # The data are the fit's own, with the draw on the model frame's rows.
  data <- get("data", envir = rp$env, inherits = FALSE)
  rows <- match(rownames(m$model@frame), rownames(data))
  expect_equal(data$.maihda_boot_y[rows], as.numeric(y))
  expect_true(all(is.na(data$.maihda_boot_y[-rows])))
  # Fitted from the refit's own estimates, the replayed call ends where the refit does.
  ref <- quiet(maihda_refit_draw(m$model, y))
  rp$call$start <- lme4::getME(ref, "theta")
  fit <- quiet(eval(rp$call, rp$env))
  expect_identical(stats::nobs(fit), stats::nobs(m$model))
  expect_equal(unname(stats::weights(fit)), unname(stats::weights(m$model)))
  expect_lt(abs(maihda_lmer_criterion(fit) - maihda_lmer_criterion(ref)), 1e-6)
  # Under na.exclude a draw padded back to the input rows is unpadded first.
  mx <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                         na.action = stats::na.exclude))
  yx <- with_seed(5, maihda_simulate_lme4(mx$model, nsim = 1))[[1]]
  padded <- stats::napredict(attr(mx$model@frame, "na.action"), as.numeric(yx))
  expect_gt(length(padded), length(yx))
  rx <- maihda_growth_draw_call(mx$model, padded)
  expect_false(is.null(rx))
  dx <- get("data", envir = rx$env, inherits = FALSE)
  expect_equal(dx$.maihda_boot_y[match(rownames(mx$model@frame), rownames(dx))], as.numeric(yx))
  fx <- quiet(eval(rx$call, rx$env))
  expect_identical(stats::nobs(fx), stats::nobs(mx$model))
})

test_that("a count draw is replayed with its offset, a count-matrix response and theta fixed", {
  skip_on_cran()
  skip_slow_test()
  d <- small_counts(11)
  replay_frame <- function(m, y) {
    rp <- maihda_growth_draw_call(m$model, y)
    ref <- quiet(maihda_refit_draw(m$model, y))
    rp$call$start <- list(theta = lme4::getME(ref, "theta"), fixef = lme4::fixef(ref))
    fit <- quiet(eval(rp$call, rp$env))
    list(rp = rp, fit = fit, ref = ref)
  }
  po <- quiet(fit_maihda(n ~ wave + offset(log(expo)) + (1 | A:B), data = d, id = "id",
                         time = "wave", family = "poisson"))
  r <- replay_frame(po, with_seed(3, maihda_simulate_lme4(po$model, nsim = 1))[[1]])
  expect_equal(lme4::getME(r$fit, "offset"), lme4::getME(r$ref, "offset"))
  expect_equal(unname(r$fit@resp$y), unname(r$ref@resp$y))
  bi <- quiet(fit_maihda(cbind(s, tr - s) ~ wave + (1 | A:B), data = d, id = "id", time = "wave",
                         family = "binomial"))
  yb <- with_seed(3, maihda_simulate_lme4(bi$model, nsim = 1))[[1]]
  expect_true(is.matrix(yb))
  r <- replay_frame(bi, yb)
  expect_equal(unname(r$fit@resp$y), unname(r$ref@resp$y))
  expect_equal(unname(r$fit@resp$weights), unname(r$ref@resp$weights))
  # A factor response written as an expression: its draw is a factor, replayed as 0 / 1.
  d$f <- factor(ifelse(d$s > 2, "hi", "lo"))
  bf <- quiet(fit_maihda(relevel(f, ref = "lo") ~ wave + (1 | A:B), data = d, id = "id",
                         time = "wave", family = "binomial"))
  yf <- with_seed(3, maihda_simulate_lme4(bf$model, nsim = 1))[[1]]
  r <- replay_frame(bf, yf)
  expect_equal(unname(r$fit@resp$y), unname(r$ref@resp$y))
  # glmer.nb() records its fit as glmer() with the fitted theta: the replay holds it fixed,
  # as lme4::refit() does.
  nb <- quiet(fit_maihda(n ~ wave + (1 | A:B), data = d, id = "id", time = "wave",
                         family = "negbinomial"))
  rp <- maihda_growth_draw_call(nb$model, with_seed(3, maihda_simulate_lme4(nb$model, nsim = 1))[[1]])
  expect_identical(rp$call[[1]], quote(lme4::glmer))
  expect_identical(eval(rp$call$family)$family, stats::family(nb$model)$family)
})

test_that("the call loses the start of a kept restart and the control of a refit across a bound", {
  env <- new.env()
  cl <- quote(lme4::glmer(n ~ t + (t | id), data = data, start = list(theta = 1),
                          control = .maihda_arg_refit_control))
  fresh <- maihda_growth_fresh_call(cl, env)
  expect_null(fresh$start)
  expect_null(fresh$control)
  assign(".maihda_arg_control", lme4::glmerControl(), envir = env)
  expect_identical(maihda_growth_fresh_call(cl, env)$control, as.name(".maihda_arg_control"))
  # The user's own control is left as it is.
  cl$control <- as.name(".maihda_arg_control")
  expect_identical(maihda_growth_fresh_call(cl, env)$control, as.name(".maihda_arg_control"))
})

# ---- the refit as the first fit ------------------------------------------------------------

# A stand-in fitting function for the restarts: each start returns a criterion and speaks.
fake_fits <- function(restart_crit) {
  env <- new.env()
  env$calls <- 0L
  env$fake_fit <- function(start = NULL, control = NULL) {
    env$calls <- env$calls + 1L
    if (is.null(start)) stop("fitted from lme4's start: the first fit was not taken")
    a <- if (is.list(start)) start$theta[["a"]] else start[["a"]]
    message("restart ", a)
    list(crit = restart_crit[[as.character(a)]])
  }
  env
}
heard <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(expr,
    message = function(m) { seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { seen <<- c(seen, paste("warning:", conditionMessage(w))); invokeRestart("muffleWarning") })
  list(value = value, seen = seen)
}
first_fit <- function(crit) list(value = list(crit = crit), conditions = list(simpleMessage("refit\n")))

test_that("the growth fitters take a draw's refit as their first fit", {
  local_mocked_bindings(
    maihda_glmer_growth_starts = function(theta, cnms, id, scale = 0.5) list(c(a = 0.5), c(a = 2)),
    maihda_growth_restart_starts = function(theta, cnms, id, scale = 0.5) list(c(a = 0.5), c(a = 2)),
    maihda_glmer_criterion = function(model) model$crit,
    maihda_lmer_criterion = function(model) model$crit,
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4, optimizer = NULL) NULL)
  # glmer: only the two restarts are fitted; a lower restart replaces the refit, and only the
  # kept fit speaks.
  env <- fake_fits(list("0.5" = 9.5, "2" = 8))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid", first = first_fit(10)))
  expect_identical(env$calls, 2L)
  expect_identical(h$value$crit, 8)
  expect_identical(h$seen, "restart 2\n")
  # The refit stays when no restart is lower by more than 1e-4, and it speaks.
  env <- fake_fits(list("0.5" = 10 - 5e-5, "2" = 11))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid", first = first_fit(10)))
  expect_identical(h$value$crit, 10)
  expect_identical(h$seen, "refit\n")
  # lmer, linear growth: the same.
  env <- fake_fits(list("0.5" = 7, "2" = 8))
  h <- heard(maihda_fit_lmer_growth(quote(fake_fit()), env, "pid", first = first_fit(10)))
  expect_identical(env$calls, 2L)
  expect_identical(h$value$crit, 7)
})

test_that("where fit_maihda() makes no refits, the draw's refit is returned as it is", {
  skip_on_cran()
  # An lmer growth fit without a stratum slope gets no restarts at fit time.
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(2), id = "id",
                        time = "wave", stratum_slope = FALSE))
  spec <- maihda_lme4_growth_spec(m$model, m$longitudinal_info)
  expect_false(is.null(spec))
  y <- with_seed(7, maihda_simulate_lme4(m$model, nsim = 1))[[1]]
  ref <- quiet(maihda_refit_draw(m$model, y))
  h <- heard(maihda_refit_growth_draw(m$model, y, spec))
  expect_identical(lme4::getME(h$value, "theta"), lme4::getME(ref, "theta"))
  expect_identical(lme4::fixef(h$value), lme4::fixef(ref))
  # No growth spec -- a cross-sectional fit, or a growth fit from the user's own start: the
  # refit alone, and nothing replayed.
  local_mocked_bindings(maihda_growth_draw_call = function(model, newresp) stop("replayed"))
  w0 <- growth_design(2)
  xs <- quiet(fit_maihda(y ~ A + (1 | A:B:C:D), data = w0[w0$wave == 0, ]))
  expect_null(maihda_lme4_growth_spec(xs$model, xs$longitudinal_info))
  yx <- with_seed(7, maihda_simulate_lme4(xs$model, nsim = 1))[[1]]
  expect_identical(lme4::fixef(maihda_refit_growth_draw(xs$model, yx, NULL)),
                   lme4::fixef(quiet(maihda_refit_draw(xs$model, yx))))
  us <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(2), id = "id",
                         time = "wave", start = c(0.5, 0, 0.2, 0.3, 0, 0.1)))
  expect_null(maihda_lme4_growth_spec(us$model, us$longitudinal_info))
})

# ---- live draws ----------------------------------------------------------------------------

test_that("a Poisson growth draw stopped short reaches its maximum; one at it keeps its refit", {
  skip_on_cran()
  skip_slow_test()
  # Time in days, 100 draws simulated with seed 11 from the fit's estimates as lme4 2.1-0 makes
  # them on Windows, written out: the fit's near-zero entry (5e-4) moves with the platform and the
  # lme4 version, and 2e-3 there changes 97 of draw 99's 576 counts. Refitted once from those
  # estimates, draw 99 stopped 1.88 below its maximum in log-likelihood, -2 logLik 2208.2049 (the
  # best of the refit, a second refit, the refit across a bound, the two restarts and the same
  # fits in waves and weeks; 61 searches from random starts found none lower); draw 4 was at its
  # maximum, 2156.0861. Where the fit's estimates are others, the refits start elsewhere.
  m <- quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = growth_design(2, 365), id = "id",
                        time = "wave", family = "poisson"))
  expect_identical(m$longitudinal_info$time_scale, 1000)
  estimates <- list(theta = c(0.323128458, 0.0345750579, 0.219021977, 0.0501899232,
                              -0.200851056, 0.000544115668),
                    beta = c(0.4183902, 0.383956752))
  sim <- with_seed(11, stats::simulate(m$model, nsim = 100, newparams = estimates))
  same_start <- isTRUE(all.equal(unname(lme4::getME(m$model, "theta")), estimates$theta,
                                 tolerance = 1e-8))
  spec <- maihda_lme4_growth_spec(m$model, m$longitudinal_info)
  old <- quiet(maihda_refit_draw(m$model, sim[[99]]))
  new <- quiet(maihda_refit_growth_draw(m$model, sim[[99]], spec))
  if (same_start) {
    expect_gt(maihda_glmer_criterion(old) - 2208.2049, 3.5)
  }
  expect_lt(abs(maihda_glmer_criterion(new) - 2208.2049), 0.01)
  # Draw 4: the refit is kept, bit for bit, at its maximum; a refit replaces it only from lower.
  old <- quiet(maihda_refit_draw(m$model, sim[[4]]))
  new <- quiet(maihda_refit_growth_draw(m$model, sim[[4]], spec))
  if (same_start || identical(lme4::getME(new, "theta"), lme4::getME(old, "theta"))) {
    expect_identical(lme4::getME(new, "theta"), lme4::getME(old, "theta"))
    expect_identical(lme4::fixef(new), lme4::fixef(old))
  } else {
    expect_lt(maihda_glmer_criterion(new), maihda_glmer_criterion(old))
    expect_lt(abs(maihda_glmer_criterion(new) - 2156.0861), 0.01)
  }
})

test_that("an lmer growth draw on the wrong side of the stratum block's bound reaches the optimum", {
  skip_on_cran()
  # Where lme4's own refit stops depends on the optimiser's path, so this runs on one leg.
  skip_slow_test()
  # Waves, 100 draws simulated with seed 11. Refitted once, draw 98 stopped with the stratum
  # intercept's Cholesky entry at 0 and the intercept-slope correlation's sign reversed,
  # 0.062 above the REML optimum 1195.1589 (multi-start bobyqa on the draw's REML criterion);
  # a second refit from there stays put.
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(1), id = "id",
                        time = "wave"))
  sim <- with_seed(11, maihda_simulate_lme4(m$model, nsim = 100))
  spec <- maihda_lme4_growth_spec(m$model, m$longitudinal_info)
  old <- quiet(maihda_refit_draw(m$model, sim[[98]]))
  expect_gt(lme4::REMLcrit(old) - 1195.1589, 0.05)
  expect_lt(lme4::getME(old, "theta")[["stratum.(Intercept)"]], 1e-4)
  new <- quiet(maihda_refit_growth_draw(m$model, sim[[98]], spec))
  expect_lt(abs(lme4::REMLcrit(new) - 1195.1589), 1e-3)
  expect_gt(lme4::getME(new, "theta")[["stratum.(Intercept)"]], 0.01)
  # A draw at its optimum keeps its refit, bit for bit.
  old <- quiet(maihda_refit_draw(m$model, sim[[1]]))
  new <- quiet(maihda_refit_growth_draw(m$model, sim[[1]], spec))
  expect_identical(lme4::getME(new, "theta"), lme4::getME(old, "theta"))
  expect_identical(lme4::fixef(new), lme4::fixef(old))
})

test_that("a quadratic lmer draw gets the centred-time start on its own response", {
  skip_on_cran()
  # The start of maihda_growth_centred_start() is fitted on the call's data, found by name in
  # the environment the call is evaluated in: here the replayed draw's.
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(4, deg = 2), id = "id",
                        time = "wave", time_degree = 2))
  y <- with_seed(2, maihda_simulate_lme4(m$model, nsim = 1))[[1]]
  spec <- maihda_lme4_growth_spec(m$model, m$longitudinal_info)
  seen <- NULL
  real <- maihda_growth_centred_start
  local_mocked_bindings(maihda_growth_centred_start = function(fit_call, fit_env, time, time_degree) {
    seen <<- list(response = get(as.character(fit_call$data), envir = fit_env, inherits = FALSE)$.maihda_boot_y,
                  time = time, time_degree = time_degree)
    real(fit_call, fit_env, time, time_degree)
  })
  ref <- quiet(maihda_refit_draw(m$model, y))
  new <- quiet(maihda_refit_growth_draw(m$model, y, spec))
  expect_equal(seen$response[!is.na(seen$response)], as.numeric(y))
  expect_identical(seen$time, "wave")
  expect_identical(seen$time_degree, 2L)
  expect_lte(lme4::REMLcrit(new), lme4::REMLcrit(ref))
})

# ---- summary() -----------------------------------------------------------------------------

test_that("summary()'s two bootstraps refit a growth fit's draws this way, other fits as before", {
  skip_on_cran()
  skip_slow_test()
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(3), id = "id",
                        time = "wave"))
  specs <- list()
  tags <- list()
  real_draw <- maihda_refit_growth_draw
  real_reduced <- maihda_refit_reduced
  # Recorded as one-element lists: a NULL spec or tag is a call too.
  local_mocked_bindings(
    maihda_refit_growth_draw = function(model, newresp, spec) {
      specs <<- c(specs, list(spec))
      real_draw(model, newresp, spec)
    },
    maihda_refit_reduced = function(model, red_form, longitudinal = NULL) {
      tags <<- c(tags, list(longitudinal))
      real_reduced(model, red_form, longitudinal)
    })
  with_seed(1, quiet(summary(m, bootstrap = TRUE, n_boot = 10)))
  expect_length(specs, 10L)
  expect_true(all(vapply(specs, function(s) identical(s$longitudinal, m$longitudinal_info), logical(1))))
  specs <- list()
  with_seed(1, quiet(summary(m, df_method = "bootstrap", n_boot = 10)))
  expect_length(specs, 10L)
  expect_true(all(vapply(specs, function(s) identical(s$id, "id"), logical(1))))
  expect_length(tags, 1L)
  expect_identical(tags[[1]], m$longitudinal_info)
  # A cross-sectional fit: no spec, and the null fit gets no tag.
  w0 <- growth_design(3)
  xs <- quiet(fit_maihda(y ~ A + (1 | A:B:C:D), data = w0[w0$wave == 0, ]))
  specs <- list()
  tags <- list()
  with_seed(1, quiet(summary(xs, df_method = "bootstrap", n_boot = 10)))
  expect_length(specs, 10L)
  expect_true(all(vapply(specs, is.null, logical(1))))
  expect_null(tags[[1]])
})

test_that("summary(df_method = \"bootstrap\") simulates a growth fit's null from the refitted null", {
  skip_on_cran()
  skip_slow_test()
  # Poisson, time in days: refitted by update() alone, the time coefficient's null stopped at
  # a variance bound, -2 logLik 2221.169; refitted as fit_maihda() fits the model, 2217.038,
  # the lowest found. (test-audit-2026-10-19-bootstrap-reduced-refit.R tests the null fits.)
  m <- quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = growth_design(1, 365), id = "id",
                        time = "wave", family = "poisson"))
  nulls <- list()
  local_mocked_bindings(maihda_simulate_lme4 = function(model, nsim) {
    nulls <<- c(nulls, list(model))
    stop("enough")
  })
  expect_error(quiet(summary(m, df_method = "bootstrap", n_boot = 10)), "could not be simulated")
  expect_length(nulls, 1L)
  expect_lt(abs(maihda_glmer_criterion(nulls[[1]]) - 2217.038), 0.005)
})
