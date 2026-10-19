# Audit pass 2026-10-18: lme4 growth fits with a quadratic or higher time polynomial that
# stop above the REML optimum.
#
# On balanced 16-stratum designs lmer's default optimiser stopped above the exact REML
# optimum in 10 to 20 per cent of quadratic growth fits, at any time scale and with or
# without a stratum slope, and the fit was recorded as converged: mostly with the
# person-level block collapsed to rank 2 where the optimum's is full rank, which the
# stratum-level restarts of a linear model (made on the first fit's person-level block)
# often cannot leave. fit_maihda() now refits such a model from lmer's fit of the same
# model on time centred at the middle of its range, carried back to the model's own time,
# and, when that refit and lmer's own fit differ by more than 0.01, from a start with both
# blocks diagonal; the fit with the lowest criterion is kept.
#
# The exact optimum of the balanced design is known in closed form (Anderson, Anderson &
# Olkin 1986, in canonical coordinates of the within- and between-stratum scatter of the
# persons' least-squares coefficients); the criteria and PCVs pinned below are its values,
# checked to 5e-3 REML units (lmer's own fits missed them by 0.19 to 5.8) so that an
# optimiser stopping a little short on another platform does not fail them.
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

# A balanced growth design: 16 factorial strata of 6 persons, waves 0..5 times k, with
# stratum and person intercepts, slopes and curvatures and a Gaussian outcome y.
growth_design <- function(seed, k = 1) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
    add0 <- 0.3 * (g$A - 0.5) + 0.2 * (g$B - 0.5) - 0.25 * (g$C - 0.5) + 0.15 * (g$D - 0.5)
    add1 <- 0.04 * (g$A - 0.5) - 0.03 * (g$B - 0.5) + 0.05 * (g$C - 0.5) - 0.02 * (g$D - 0.5)
    u <- cbind(stats::rnorm(16, 0, 0.2), stats::rnorm(16, 0, 0.06), stats::rnorm(16, 0, 0.02))
    p <- cbind(stats::rnorm(96, 0, 0.5), stats::rnorm(96, 0, 0.15), stats::rnorm(96, 0, 0.06))
    id <- rep(seq_len(96), each = 6)
    st <- rep(seq_len(16), each = 6)[id]
    t <- rep(0:5, times = 96)
    eta <- add0[st] + u[st, 1] + (add1[st] + u[st, 2]) * t + u[st, 3] * t^2 +
      p[id, 1] + p[id, 2] * t + p[id, 3] * t^2
    data.frame(id = id, wave = t * k,
               A = factor(g$A[st]), B = factor(g$B[st]), C = factor(g$C[st]), D = factor(g$D[st]),
               y = 2 + 0.25 * t - 0.02 * t^2 + eta + stats::rnorm(576, 0, 0.5))
  })
}

quad_cnms <- list(id = c("(Intercept)", "wave", "I(wave^2)"),
                  stratum = c("(Intercept)", "wave", "I(wave^2)"))
quad_theta <- stats::setNames(c(0.9, -0.05, 0.02, 0.3, -0.06, 0.13, 0.15, -0.08, 0.009, 0.18, -0.002, 0.01),
                              paste0("theta", 1:12))

block_of <- function(theta, cnms, g) {
  sizes <- vapply(cnms, function(cn) length(cn) * (length(cn) + 1) / 2, numeric(1))
  end <- cumsum(sizes)[[g]]
  q <- length(cnms[[g]])
  L <- matrix(0, q, q)
  L[lower.tri(L, diag = TRUE)] <- theta[seq.int(end - sizes[[g]] + 1, end)]
  L
}

# ---- the pieces ----------------------------------------------------------------------------

test_that("a positive semi-definite block is factored as lme4 holds it", {
  S <- matrix(c(4, 2, 1, 2, 3, 0.5, 1, 0.5, 2), 3)
  expect_equal(maihda_chol_psd(S), t(chol(S)), tolerance = 1e-12)
  # Rank 2: the third column is left at zero and the factor reproduces the block.
  v1 <- c(1, 0.5, -0.2)
  v2 <- c(0, 0.3, 0.6)
  S2 <- tcrossprod(v1) + tcrossprod(v2)
  L2 <- maihda_chol_psd(S2)
  expect_equal(tcrossprod(L2), S2, tolerance = 1e-12)
  expect_identical(L2[3, 3], 0)
  expect_true(all(diag(L2) >= 0))
  expect_identical(L2[upper.tri(L2)], c(0, 0, 0))
  # A zero first variance leaves the first column at zero.
  S3 <- matrix(c(0, 0, 0, 0, 1, 0.5, 0, 0.5, 0.25), 3)
  L3 <- maihda_chol_psd(S3)
  expect_identical(L3[, 1], c(0, 0, 0))
  expect_equal(tcrossprod(L3), S3, tolerance = 1e-12)
  expect_identical(maihda_chol_psd(matrix(0, 2, 2)), matrix(0, 2, 2))
  expect_equal(maihda_chol_psd(matrix(0.49, 1, 1)), matrix(0.7, 1, 1))
})

test_that("a theta fitted on shifted time is carried back by U Sigma U'", {
  sh <- 2.5
  back <- maihda_growth_unshift_theta(quad_theta, quad_cnms, "wave", 2, sh)
  # Random effects over (1, t - c, (t - c)^2) are U b over (1, t, t^2).
  U <- matrix(c(1, 0, 0, -sh, 1, 0, sh^2, -2 * sh, 1), 3)
  for (g in c("id", "stratum")) {
    Lc <- block_of(quad_theta, quad_cnms, g)
    Lb <- block_of(back, quad_cnms, g)
    expect_equal(tcrossprod(Lb), U %*% tcrossprod(Lc) %*% t(U), tolerance = 1e-12)
    expect_true(all(diag(Lb) >= 0))
    # The variance a block implies at each time is the same function of time.
    tt <- c(0, 1.7, 5)
    expect_equal(maihda_var_at_time(tcrossprod(Lb), tt),
                 maihda_var_at_time(tcrossprod(Lc), tt - sh), tolerance = 1e-10)
  }
  expect_identical(names(back), names(quad_theta))
  # A rank-deficient block stays rank-deficient: its zero pivot is not turned into a
  # rounding residue.
  r2 <- quad_theta
  r2[12] <- 0
  expect_identical(unname(maihda_growth_unshift_theta(r2, quad_cnms, "wave", 2, 1.7)[12]), 0)
  # Back and forth.
  expect_equal(maihda_growth_unshift_theta(back, quad_cnms, "wave", 2, -sh), quad_theta,
               tolerance = 1e-10)
  # A cubic block takes the binomial coefficients: the variance function is kept.
  c_cnms <- list(id = c("(Intercept)", "t", "I(t^2)", "I(t^3)"))
  Lc <- matrix(0, 4, 4)
  Lc[lower.tri(Lc, diag = TRUE)] <- c(0.8, -0.1, 0.05, 0.01, 0.3, -0.04, 0.02, 0.1, -0.01, 0.03)
  cb <- maihda_growth_unshift_theta(Lc[lower.tri(Lc, diag = TRUE)], c_cnms, "t", 3, 1.5)
  Lb <- block_of(cb, c_cnms, "id")
  expect_equal(maihda_var_at_time(tcrossprod(Lb), c(0, 0.4, 2, 3)),
               maihda_var_at_time(tcrossprod(Lc), c(0, 0.4, 2, 3) - 1.5), tolerance = 1e-10)
  # An intercept-only stratum block (stratum_slope = FALSE) is unchanged.
  ns_cnms <- list(id = quad_cnms$id, stratum = "(Intercept)")
  ns_theta <- c(quad_theta[1:6], 0.4)
  ns <- maihda_growth_unshift_theta(ns_theta, ns_cnms, "wave", 2, sh)
  expect_identical(ns[[7]], 0.4)
  expect_equal(ns[1:6], back[1:6])
  # Blocks that are not the polynomial's leading terms in order, and bad input.
  none <- function(...) expect_null(maihda_growth_unshift_theta(...))
  none(quad_theta, quad_cnms, "t", 2, sh)
  none(quad_theta, list(id = quad_cnms$id[c(1, 3, 2)], stratum = quad_cnms$stratum), "wave", 2, sh)
  none(c(quad_theta[1:6], 1, 0, 1), list(id = quad_cnms$id, stratum = c("(Intercept)", "I(wave^2)")),
       "wave", 2, sh)
  none(quad_theta, quad_cnms, "wave", 1, sh)
  none(quad_theta[-1], quad_cnms, "wave", 2, sh)
  none(quad_theta, quad_cnms, "wave", 2, NA_real_)
  none(NULL, quad_cnms, "wave", 2, sh)
})

test_that("the diagonal start puts every correlation at zero", {
  st <- maihda_growth_diagonal_start(quad_theta, quad_cnms, id = "id")
  expect_length(st, 1L)
  sds <- sqrt(rowSums(block_of(quad_theta, quad_cnms, "id")^2))
  expect_equal(unname(st[[1]][1:6]), c(sds[1], 0, 0, sds[2], 0, sds[3]))
  expect_equal(unname(st[[1]][7:12]), 0.5 * c(sds[1], 0, 0, sds[2], 0, sds[3]))
  expect_identical(names(st[[1]]), names(quad_theta))
  # Whatever the order of the terms; an intercept-only stratum block takes the intercept's.
  rv <- maihda_growth_diagonal_start(quad_theta[c(7:12, 1:6)], rev(quad_cnms), id = "id")
  expect_equal(rv[[1]], st[[1]][c(7:12, 1:6)])
  ns <- maihda_growth_diagonal_start(c(quad_theta[1:6], 0.4),
                                     list(id = quad_cnms$id, stratum = "(Intercept)"), id = "id")
  expect_equal(unname(ns[[1]]), c(sds[1], 0, 0, sds[2], 0, sds[3], 0.5 * sds[1]))
  none <- function(...) expect_identical(maihda_growth_diagonal_start(...), list())
  z <- quad_theta; z[1] <- 0; z[2:3] <- 0; none(z, quad_cnms, id = "id")
  z <- quad_theta; z[6] <- NA; none(z, quad_cnms, id = "id")
  none(quad_theta, quad_cnms, id = "pid")
  none(quad_theta[1:6], quad_cnms["id"], id = "id")
  none(c(0.9, 0.1, 0.3, quad_theta[7:12]), list(id = c("(Intercept)", "wave"), stratum = quad_cnms$stratum),
       id = "id")
  none(quad_theta[-1], quad_cnms, id = "id")
  none(quad_theta, unname(quad_cnms), id = "id")
})

test_that("the centred start is the same model fitted on time centred at mid-range, carried back", {
  skip_on_cran()
  d <- growth_design(1)
  d$stratum <- interaction(d$A, d$B, d$C, d$D, drop = TRUE)
  # Half the persons miss the last wave: the times run from 0 to 5, mean 2.27.
  d <- d[!(d$id <= 48 & d$wave == 5), ]
  f <- y ~ wave + I(wave^2) + (wave + I(wave^2) | id) + (wave + I(wave^2) | stratum)
  env <- new.env()
  env$data <- d
  fit_call <- quote(lme4::lmer(formula = NULL, data = data))
  fit_call$formula <- f
  st <- maihda_growth_centred_start(fit_call, env, "wave", 2)
  expect_length(st, 1L)
  dc <- d
  dc$wave <- dc$wave - 2.5
  fc <- lme4::lmer(f, data = dc)
  cn <- lme4::getME(fc, "cnms")
  expect_identical(st[[1]], maihda_growth_unshift_theta(lme4::getME(fc, "theta"), cn, "wave", 2, 2.5))
  expect_identical(env$data, d)
  # Shifting time leaves the REML criterion where it was: the start is that fit's optimum.
  dev0 <- lme4::lmer(f, data = d, devFunOnly = TRUE)
  dev1 <- lme4::lmer(f, data = dc, devFunOnly = TRUE)
  expect_lt(abs(dev0(st[[1]]) - dev1(lme4::getME(fc, "theta"))), 1e-6)
  # No start when the data or the column cannot be read, or time is centred already.
  call2 <- quote(lme4::lmer(formula = NULL, data = data[1:300, ]))
  expect_identical(maihda_growth_centred_start(call2, env, "wave", 2), list())
  expect_identical(maihda_growth_centred_start(fit_call, env, "no_such_column", 2), list())
  expect_identical(maihda_growth_centred_start(fit_call, env, c("wave", "y"), 2), list())
  env_na <- new.env()
  env_na$data <- transform(d, wave = NA_real_)
  expect_identical(maihda_growth_centred_start(fit_call, env_na, "wave", 2), list())
  env_c <- new.env()
  env_c$data <- dc
  expect_identical(maihda_growth_centred_start(fit_call, env_c, "wave", 2), list())
})

# ---- which fit is kept ---------------------------------------------------------------------

test_that("a quadratic fit is refitted from the centred start, and from the diagonal one when the two differ by > 0.01", {
  env <- new.env()
  crits <- list()
  env$fake_fit <- function(start = NULL) {
    nm <- if (is.null(start)) "first" else names(start)
    message("fit ", nm)
    if (is.function(crits[[nm]])) crits[[nm]]() else list(crit = crits[[nm]])
  }
  seen <- character(0)
  local_mocked_bindings(
    maihda_growth_centred_start = function(fit_call, fit_env, time, time_degree) {
      seen <<- c(seen, paste("centred", time, time_degree))
      list(c(centred = 1))
    },
    maihda_growth_diagonal_start = function(theta, cnms, id, scale = 0.5) {
      seen <<- c(seen, paste("diagonal", id))
      list(c(diagonal = 1))
    },
    maihda_growth_restart_starts = function(theta, cnms, id, scale = 0.5) {
      seen <<- c(seen, "signs")
      list()
    },
    maihda_lmer_criterion = function(model) model$crit
  )
  run <- function(first, centred, diagonal) {
    crits <<- list(first = first, centred = centred, diagonal = diagonal)
    seen <<- character(0)
    msgs <- character(0)
    value <- withCallingHandlers(
      maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid", time = "wave", time_degree = 2),
      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
    list(crit = value$crit, seen = seen, msgs = msgs)
  }
  # The centred refit agrees with lmer's fit: lmer's fit is kept and no third fit is made.
  r <- run(10, 10 - 5e-5, 1)
  expect_identical(r$crit, 10)
  expect_identical(r$seen, "centred wave 2")
  expect_identical(r$msgs, "fit first\n")
  # It is better by 0.005: it replaces lmer's fit, and a difference that small does not
  # call for the diagonal start.
  r <- run(10, 9.995, 1)
  expect_identical(r$crit, 9.995)
  expect_identical(r$seen, "centred wave 2")
  expect_identical(r$msgs, "fit centred\n")
  # Worse by 0.005: lmer's fit is kept, again without the diagonal start.
  r <- run(10, 10.005, 1)
  expect_identical(r$crit, 10)
  expect_identical(r$seen, "centred wave 2")
  # Better by 0.02: the diagonal start is tried as well, and the centred refit is kept
  # when it is the lowest.
  r <- run(10, 9.98, 11)
  expect_identical(r$crit, 9.98)
  expect_identical(r$seen, c("centred wave 2", "diagonal pid"))
  expect_identical(r$msgs, "fit centred\n")
  # Worse by 0.02, or by more: the diagonal start is tried and wins when it is the lowest.
  r <- run(10, 10.02, 9.9)
  expect_identical(r$crit, 9.9)
  expect_identical(r$msgs, "fit diagonal\n")
  r <- run(10, 12, 9)
  expect_identical(r$crit, 9)
  r <- run(10, 9, 8)
  expect_identical(r$crit, 8)
  # A gain below 1e-4 does not replace lmer's fit.
  r <- run(10, 11, 10 - 5e-5)
  expect_identical(r$crit, 10)
  expect_identical(r$seen, c("centred wave 2", "diagonal pid"))
  # A centred refit that fails, or a centred start that fails, sends the fit on to the
  # diagonal start.
  r <- run(10, function() stop("no fit"), 9.5)
  expect_identical(r$crit, 9.5)
  local_mocked_bindings(maihda_growth_centred_start = function(...) stop("no start"))
  r <- run(10, 1, 9.5)
  expect_identical(r$crit, 9.5)
  expect_identical(r$seen, "diagonal pid")
  local_mocked_bindings(maihda_growth_centred_start = function(...) list())
  r <- run(10, 1, 9.5)
  expect_identical(r$crit, 9.5)
  # A linear model takes the sign starts and nothing else.
  seen <- character(0)
  r <- withCallingHandlers(
    maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid", time = "wave", time_degree = 1),
    message = function(m) invokeRestart("muffleMessage"))
  expect_identical(seen, "signs")
})

test_that("lmer growth fits of degree 2 or more are restarted, with or without a stratum slope", {
  lng <- list(id = "id", stratum_slope = TRUE, time_degree = 2)
  expect_true(maihda_growth_restarts_apply(quote(lme4::lmer), lng, list()))
  expect_true(maihda_growth_restarts_apply(quote(lme4::lmer), utils::modifyList(lng, list(time_degree = 3L)), list()))
  expect_true(maihda_growth_restarts_apply(quote(lme4::lmer), utils::modifyList(lng, list(stratum_slope = FALSE)),
                                           list()))
  expect_true(maihda_growth_restarts_apply(quote(lme4::lmer), lng, list(REML = FALSE)))
  expect_false(maihda_growth_restarts_apply(quote(lme4::lmer), lng, list(start = 1)))
  expect_false(maihda_growth_restarts_apply(quote(lme4::glmer), lng, list()))
  # A linear model without a stratum slope still is not.
  expect_false(maihda_growth_restarts_apply(quote(lme4::lmer),
                                            list(id = "id", stratum_slope = FALSE, time_degree = 1), list()))
  expect_false(maihda_growth_restarts_apply(quote(lme4::lmer), utils::modifyList(lng, list(time_degree = NA)),
                                            list()))
  expect_false(maihda_growth_restarts_apply(quote(lme4::lmer), list(id = "id"), list()))
})

test_that("fit_maihda() hands the restarts the model's own time and degree", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- maihda_long_data[maihda_long_data$id %in% unique(maihda_long_data$id)[1:120], ]
  args <- list()
  real <- maihda_fit_lmer_growth
  local_mocked_bindings(maihda_fit_lmer_growth = function(fit_call, fit_env, id, ...) {
    args[[length(args) + 1L]] <<- list(...)
    real(fit_call, fit_env, id, ...)
  })
  f <- wellbeing ~ wave + (1 | gender:ethnicity:education)
  quiet(fit_maihda(f, data = d, id = "id", time = "wave", time_degree = 2))
  quiet(fit_maihda(f, data = d, id = "id", time = "wave", time_degree = 2, stratum_slope = FALSE))
  d10 <- d
  d10$wave <- d10$wave + 10
  quiet(fit_maihda(f, data = d10, id = "id", time = "wave", time_degree = 2))
  expect_length(args, 3L)
  expect_identical(args[[1]], list(time = "wave", time_degree = 2L))
  expect_identical(args[[2]], list(time = "wave", time_degree = 2L))
  # A centred axis: the restarts work on the internal time column.
  expect_identical(args[[3]], list(time = ".maihda_ctime", time_degree = 2L))
  quiet(fit_maihda(f, data = d, id = "id", time = "wave", time_degree = 2,
                   start = rep(c(1, 0, 0, 1, 0, 1), 2)))
  quiet(fit_maihda(f, data = d, id = "id", time = "wave", stratum_slope = FALSE))
  expect_length(args, 3L)
})

# ---- the fits it changes -------------------------------------------------------------------

test_that("quadratic fits stuck in either block reach the exact optimum", {
  skip_on_cran()
  fit <- function(seed) {
    quiet(maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(seed), id = "id", time = "wave",
                 time_degree = 2, decomposition = "longitudinal"))
  }
  crit <- function(a) c(lme4::REMLcrit(a$model$model), lme4::REMLcrit(a$model_adjusted$model))
  # Seed 1: lmer's adjusted fit stopped 0.537 above the optimum with the person-level
  # block collapsed to rank 2 (PCV 0.950 / -0.182).
  a1 <- fit(1)
  expect_lt(max(abs(crit(a1) - c(1364.094411, 1385.235936))), 5e-3)
  expect_false(is.null(a1$model_adjusted$model@call$start))
  expect_lt(abs(a1$pcv$pcv_intercept - 0.968129), 5e-3)
  expect_lt(abs(a1$pcv$pcv_slope - 0.004328), 0.05)
  # Seed 6: lmer's null fit stopped 5.81 above it with the stratum intercept SD at exactly
  # 0, so PCV_intercept was NA.
  a6 <- fit(6)
  expect_lt(max(abs(crit(a6) - c(1313.788951, 1336.808919))), 5e-3)
  expect_gt(a6$pcv$Sigma_stratum_null[1, 1], 0.05)
  expect_lt(abs(a6$pcv$pcv_intercept - 0.706629), 5e-3)
  # Seed 7: both fits stopped short (0.22 and 0.19), with the stratum block on a boundary
  # the optimum does not have; PCV_slope was NA.
  a7 <- fit(7)
  expect_lt(max(abs(crit(a7) - c(1257.644313, 1281.390388))), 5e-3)
  expect_true(is.finite(a7$pcv$pcv_slope))
})

test_that("a quadratic fit without a stratum slope reaches the lowest criterion found", {
  skip_on_cran()
  a <- quiet(maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(5), id = "id", time = "wave",
                    time_degree = 2, stratum_slope = FALSE, decomposition = "longitudinal"))
  # lmer stopped 0.214 and 0.419 above the lowest criteria found by seven optimisers and
  # random starts, 1310.261059 and 1334.339547.
  expect_lt(abs(lme4::REMLcrit(a$model$model) - 1310.261059), 5e-3)
  expect_lt(abs(lme4::REMLcrit(a$model_adjusted$model) - 1334.339547), 5e-3)
})

test_that("lmer's own quadratic fit is returned unless a restart is lower by more than 1e-4", {
  skip_on_cran()
  fit <- function() {
    quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(2), id = "id", time = "wave",
                     time_degree = 2))
  }
  m <- fit()
  # lmer's own fit, through the same call without the restarts. Here lmer reaches the
  # optimum itself and its fit comes back unchanged; an optimiser that stops a little short
  # of it is replaced by a fit lower by more than 1e-4.
  local_mocked_bindings(maihda_growth_restarts_apply = function(...) FALSE)
  own <- fit()
  if (is.null(m$model@call$start)) {
    expect_identical(lme4::getME(m$model, "theta"), lme4::getME(own$model, "theta"))
    expect_identical(lme4::REMLcrit(m$model), lme4::REMLcrit(own$model))
  } else {
    expect_lt(lme4::REMLcrit(m$model), lme4::REMLcrit(own$model) - 1e-4)
  }
  expect_lt(lme4::REMLcrit(m$model) - 1319.370830, 5e-3)
})
