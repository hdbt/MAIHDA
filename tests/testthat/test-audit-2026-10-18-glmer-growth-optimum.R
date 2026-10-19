# Audit pass 2026-10-18: glmer growth fits that stop below the maximum of their likelihood.
#
# glmer() fits a growth model in two stages, bobyqa on the variance parameters (nAGQ = 0)
# and then Nelder-Mead on all parameters from there. The first stage often stops with a
# block's intercept Cholesky entry at its bound of 0 beside a non-zero intercept-slope
# entry; the covariance matrix there is the same with that entry's sign reversed, and the
# likelihood is higher on that other side, which the second stage cannot reach across the
# bound. Counted in fine units such as days, lme4's fixed optimiser steps do not suit the
# slope entries at all. On a Poisson growth design (16 strata of 6 persons, waves 0..5)
# about a third of 120 fits ended more than 0.01 below the maximum, by up to 32, nearly all
# recorded as converged. A glmer growth fit is now made on the internally scaled time, a
# linear one with a stratum slope is refitted from two more starts, any one with a block's
# intercept variance at zero is refitted once from the other side of that bound, and the
# convergence record flags a growth fit stopped at such a bound where the likelihood still
# rises.
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

# The finding's design: 16 factorial strata of 6 persons, waves 0..5 times k, stratum and
# person intercepts and slopes, a Poisson count n.
growth_design <- function(seed, k = 1) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
    add0 <- 0.3 * (g$A - 0.5) + 0.2 * (g$B - 0.5) - 0.25 * (g$C - 0.5) + 0.15 * (g$D - 0.5)
    add1 <- 0.04 * (g$A - 0.5) - 0.03 * (g$B - 0.5) + 0.05 * (g$C - 0.5) - 0.02 * (g$D - 0.5)
    u <- cbind(stats::rnorm(16, 0, 0.2), stats::rnorm(16, 0, 0.06), stats::rnorm(16, 0, 0))
    p <- cbind(stats::rnorm(96, 0, 0.5), stats::rnorm(96, 0, 0.15), stats::rnorm(96, 0, 0))
    id <- rep(seq_len(96), each = 6)
    st <- rep(seq_len(16), each = 6)[id]
    t <- rep(0:5, times = 96)
    eta <- add0[st] + u[st, 1] + (add1[st] + u[st, 2]) * t + p[id, 1] + p[id, 2] * t
    y <- stats::rnorm(576, 0, 0.5)
    data.frame(id = id, wave = t * k,
               A = factor(g$A[st]), B = factor(g$B[st]), C = factor(g$C[st]), D = factor(g$D[st]),
               n = stats::rpois(576, exp(0.5 + 0.1 * t + 0.6 * eta)))
  })
}

poisson_growth <- function(d, ...) {
  quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave",
                   family = "poisson", ...))
}

# Two 2 x 2 growth blocks as lme4 orders them: person ("pid") first, then stratum.
growth_cnms <- list(pid = c("(Intercept)", "wave"), stratum = c("(Intercept)", "wave"))
growth_theta <- c(pid.1 = 0.4, pid.2 = 0.1, pid.3 = 0.2, st.1 = 0.3, st.2 = -0.05, st.3 = 0.02)

# ---- which fits ---------------------------------------------------------------------------

test_that("every lme4 growth fit is made on the internally scaled time, a brms one is not", {
  for (fam in list("gaussian", "poisson", "binomial", "negbinomial", stats::poisson,
                   stats::gaussian(link = "log"), MASS::negative.binomial(2))) {
    expect_true(maihda_longitudinal_scales_time("lme4", fam))
  }
  expect_false(maihda_longitudinal_scales_time("brms", "poisson"))
  expect_false(maihda_longitudinal_scales_time("brms", "gaussian"))
  expect_false(maihda_longitudinal_scales_time("lme4", "no such family"))
})

test_that("glmer and glmer.nb growth fits are refitted, unless the user gives a start", {
  lng <- list(id = "id", time = "wave", time_degree = 1L, stratum_slope = TRUE)
  expect_true(maihda_glmer_growth_applies(quote(lme4::glmer), lng, list()))
  expect_true(maihda_glmer_growth_applies(quote(lme4::glmer.nb), lng, list()))
  expect_true(maihda_glmer_growth_applies(quote(lme4::glmer), lng, list(control = 1)))
  # Every growth structure: the restarts need two 2 x 2 blocks, the step across a bound
  # does not.
  expect_true(maihda_glmer_growth_applies(quote(lme4::glmer),
                                          utils::modifyList(lng, list(stratum_slope = FALSE)), list()))
  expect_true(maihda_glmer_growth_applies(quote(lme4::glmer),
                                          utils::modifyList(lng, list(time_degree = 2L)), list()))
  expect_false(maihda_glmer_growth_applies(quote(lme4::glmer), lng, list(start = 1)))
  expect_false(maihda_glmer_growth_applies(quote(lme4::glmer), NULL, list()))
  expect_false(maihda_glmer_growth_applies(quote(lme4::lmer), lng, list()))
  expect_false(maihda_glmer_growth_applies(quote(brms::brm), lng, list()))
  # The lmer restarts stay lmer's own.
  expect_false(maihda_growth_restarts_apply(quote(lme4::glmer), lng, list()))
})

# ---- the starts and the criterion ---------------------------------------------------------

test_that("the restarts reset the person block and start the stratum block at half of it", {
  starts <- maihda_glmer_growth_starts(growth_theta, growth_cnms, id = "pid")
  expect_length(starts, 2L)
  sd_int <- 0.4
  sd_slope <- sqrt(0.1^2 + 0.2^2)
  expect_equal(unname(starts[[1]]), c(sd_int, 0, sd_slope, 0.5 * c(sd_int, sd_slope, sd_slope)))
  expect_equal(unname(starts[[2]]), c(sd_int, 0, sd_slope, 0.5 * c(sd_int, -sd_slope, sd_slope)))
  expect_identical(names(starts[[1]]), names(growth_theta))
  # The slices follow the order of the terms.
  rs <- maihda_glmer_growth_starts(growth_theta[c(4:6, 1:3)], rev(growth_cnms), id = "pid")
  expect_equal(unname(rs[[2]]), c(0.5 * c(sd_int, -sd_slope, sd_slope), sd_int, 0, sd_slope))
  # A person SD at its bound (1e-4 or less) starts at 1, glmer's own start.
  at_bound <- growth_theta
  at_bound[1] <- 0
  expect_equal(unname(maihda_glmer_growth_starts(at_bound, growth_cnms, "pid")[[1]][1:4]),
               c(1, 0, sd_slope, 0.5))
  at_bound[1] <- 1e-4
  expect_equal(unname(maihda_glmer_growth_starts(at_bound, growth_cnms, "pid")[[1]][1]), 1)
  at_bound[1] <- 2e-4
  expect_equal(unname(maihda_glmer_growth_starts(at_bound, growth_cnms, "pid")[[1]][1]), 2e-4)
  flat <- growth_theta
  flat[2:3] <- 0
  expect_equal(unname(maihda_glmer_growth_starts(flat, growth_cnms, "pid")[[2]]),
               c(0.4, 0, 1, 0.5 * c(0.4, -1, 1)))
  # No restarts unless both blocks are 2 x 2.
  none <- function(theta, cnms, id = "pid") {
    expect_identical(maihda_glmer_growth_starts(theta, cnms, id = id), list())
  }
  none(c(growth_theta[1:3], st = 0.3), list(pid = growth_cnms$pid, stratum = "(Intercept)"))
  none(c(growth_theta, 0.1, 0.1, 0.1), list(pid = c("(Intercept)", "wave", "w2"),
                                            stratum = c("(Intercept)", "wave")))
  none(growth_theta, growth_cnms, id = "person")
  none(growth_theta[-1], growth_cnms)
  none(growth_theta, unname(growth_cnms))
  none("a", growth_cnms)
})

test_that("the criterion of a glmer fit is -2 log-likelihood, not its deviance()", {
  skip_on_cran()
  d <- growth_design(31)
  m <- quiet(lme4::glmer(n ~ wave + (1 | id), data = d, family = stats::poisson))
  expect_identical(maihda_glmer_criterion(m), -2 * as.numeric(stats::logLik(m)))
  # deviance() of a glmer fit is the conditional deviance, a different number.
  expect_gt(abs(stats::deviance(m) - maihda_glmer_criterion(m)), 1)
})

# ---- the refits, with stand-in fits -------------------------------------------------------

# A stand-in fitting function: lme4's first fit, then a restart per start, each returning a
# criterion and speaking a condition.
fake_fits <- function(first_crit, restart_crit) {
  env <- new.env()
  env$fake_fit <- function(start = NULL, control = NULL) {
    if (is.null(start)) {
      message("first fit says")
      list(crit = first_crit, start = NULL)
    } else {
      message("restart ", start$theta[["a"]])
      list(crit = restart_crit[[as.character(start$theta[["a"]])]], start = start)
    }
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

test_that("the fit with the lowest criterion is kept, and only it speaks", {
  local_mocked_bindings(
    maihda_glmer_growth_starts = function(theta, cnms, id, scale = 0.5) list(c(a = 0.5), c(a = 2)),
    maihda_glmer_criterion = function(model) model$crit,
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4) NULL
  )
  env <- fake_fits(10, list("0.5" = 9.5, "2" = 8))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, id = "pid"))
  expect_identical(h$value$crit, 8)
  expect_identical(h$seen, "restart 2\n")
  # The restart is made with the start as glmer's theta.
  expect_identical(h$value$start, list(theta = c(a = 2)))
  # The first restart wins when it is the lowest.
  env <- fake_fits(10, list("0.5" = 7, "2" = 8))
  expect_identical(heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))$value$crit, 7)
  # A restart lower by 5e-5 does not replace the first fit; one lower by 0.01 does.
  env <- fake_fits(10, list("0.5" = 10 - 5e-5, "2" = 10 + 1))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))
  expect_identical(h$value$crit, 10)
  expect_identical(h$seen, "first fit says\n")
  env <- fake_fits(10, list("0.5" = 10 - 0.01, "2" = 10 + 1))
  expect_identical(heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))$value$crit, 10 - 0.01)
  # A failing restart is passed over; a failing first fit fails.
  env <- fake_fits(10, list("0.5" = 9, "2" = 8))
  failing <- env$fake_fit
  env$fake_fit <- function(start = NULL, control = NULL) {
    if (!is.null(start) && start$theta[["a"]] == 2) stop("no")
    failing(start)
  }
  expect_identical(heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))$value$crit, 9)
  env$fake_fit <- function(start = NULL, control = NULL) stop("the first fit fails")
  expect_error(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"), "the first fit fails")
})

test_that("the kept fit is refitted across a bound, and the refit kept only when lower", {
  across <- NULL
  local_mocked_bindings(
    maihda_glmer_growth_starts = function(theta, cnms, id, scale = 0.5) list(c(a = 0.5), c(a = 2)),
    maihda_glmer_criterion = function(model) model$crit,
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4) {
      across <<- kept$value$crit
      if (kept$value$crit < 9) list(value = list(crit = 6), conditions = list(simpleMessage("across\n"))) else NULL
    }
  )
  # The step is taken from the fit the restarts kept, and its refit speaks alone.
  env <- fake_fits(10, list("0.5" = 9.5, "2" = 8))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))
  expect_identical(across, 8)
  expect_identical(h$value$crit, 6)
  expect_identical(h$seen, "across\n")
  # No refit returned: the kept fit stands.
  env <- fake_fits(10, list("0.5" = 9.5, "2" = 9.2))
  h <- heard(maihda_fit_glmer_growth(quote(fake_fit()), env, "pid"))
  expect_identical(h$value$crit, 9.2)
  expect_identical(h$seen, "restart 2\n")
})

# ---- the columns at a bound ---------------------------------------------------------------

test_that("the columns of a Cholesky factor that can sit at a boundary trap", {
  cols <- maihda_lme4_theta_columns(growth_cnms)
  expect_identical(vapply(cols, `[[`, integer(1), "index"), c(1L, 4L))
  expect_identical(lapply(cols, `[[`, "below"), list(2L, 5L))
  expect_identical(vapply(cols, `[[`, character(1), "block"), c("pid", "stratum"))
  quad <- list(pid = c("(Intercept)", "t", "t2"), stratum = "(Intercept)")
  qc <- maihda_lme4_theta_columns(quad)
  expect_identical(vapply(qc, `[[`, integer(1), "index"), c(1L, 4L))
  expect_identical(lapply(qc, `[[`, "below"), list(2:3, 5L))
  expect_identical(vapply(qc, `[[`, character(1), "term"), c("(Intercept)", "t"))
  expect_length(maihda_lme4_theta_columns(list(stratum = "(Intercept)")), 0L)
  # At a bound: the diagonal entry below 1e-4, an entry below it not 0.
  th <- growth_theta
  expect_length(maihda_lme4_trap_columns(th, growth_cnms), 0L)
  th[1] <- 0
  expect_identical(vapply(maihda_lme4_trap_columns(th, growth_cnms), `[[`, integer(1), "index"), 1L)
  th[1] <- 9e-5
  expect_length(maihda_lme4_trap_columns(th, growth_cnms), 1L)
  th[1] <- 1e-4
  expect_length(maihda_lme4_trap_columns(th, growth_cnms), 0L)
  th[c(1, 2)] <- 0
  expect_length(maihda_lme4_trap_columns(th, growth_cnms), 0L)
  th[c(4, 5)] <- c(0, -1e-9)
  expect_identical(vapply(maihda_lme4_trap_columns(th, growth_cnms), `[[`, integer(1), "index"), 4L)
  # Taken across the bound: the diagonal entry set, the entries below negated.
  m <- maihda_lme4_mirror_theta(growth_theta, growth_cnms, 4L, 0.05)
  expect_equal(unname(m), c(0.4, 0.1, 0.2, 0.05, 0.05, 0.02))
  # A quadratic person block: its first two columns have entries below the diagonal.
  m <- maihda_lme4_mirror_theta(c(0, 0.3, -0.2, 0.5, 0.6, 0.7, 0.1), quad, c(1L, 4L), 0.01)
  expect_equal(m, c(0.01, -0.3, 0.2, 0.01, -0.6, 0.7, 0.1))
})

# ---- live fits ----------------------------------------------------------------------------

test_that("a glmer growth fit reaches the maximum in every unit of time", {
  skip_on_cran()
  skip_slow_test()
  # Seed 3: lme4's own fit stops 6.37 below the maximum in waves (the person intercept SD at
  # exactly 0) and 24.8 below it in days. The maximum, from the best of about 40 fits by
  # several optimisers and starts, is -1082.00427.
  fits <- lapply(c(waves = 1, weeks = 7, days = 365), function(k) poisson_growth(growth_design(3, k)))
  for (u in names(fits)) {
    ll <- as.numeric(stats::logLik(fits[[u]]$model))
    expect_lt(abs(ll - (-1082.00427)), 0.005, label = u)
    expect_true(fits[[u]]$diagnostics$converged, label = u)
    expect_null(fits[[u]]$diagnostics$boundary_traps)
  }
  expect_identical(vapply(fits, function(f) f$longitudinal_info$time_scale, numeric(1)),
                   c(waves = 1, weeks = 10, days = 1000))
  expect_gt(lme4::getME(fits$days$model, "theta")[["id.(Intercept)"]], 0.3)
  # The time-varying VPC is the same in every unit.
  vpc <- lapply(fits, function(f) quiet(summary(f))$longitudinal$vpc_t$estimate)
  expect_equal(vpc$weeks, vpc$waves, tolerance = 1e-3)
  expect_equal(vpc$days, vpc$waves, tolerance = 1e-3)
})

test_that("a fit the restarts leave at a stratum bound is refitted from across it", {
  skip_on_cran()
  # Seed 30 in weeks: both restarts leave the stratum intercept SD at 0, 0.174 below the
  # maximum of -1021.05075. They tie to 4e-7, and for the one kept here a step of 0.001
  # across the bound raises the deviance; the refit from across the bound reaches the
  # maximum whichever is kept.
  m <- poisson_growth(growth_design(30, 7))
  expect_lt(abs(as.numeric(stats::logLik(m$model)) - (-1021.05075)), 0.005)
  expect_gt(lme4::getME(m$model, "theta")[["stratum.(Intercept)"]], 0.01)
  expect_true(m$diagnostics$converged)
})

test_that("a growth fit without a stratum slope is refitted across a bound too", {
  skip_on_cran()
  # stratum_slope = FALSE in days: on the internal time lme4's fit stops with the person
  # intercept SD at 0, 10.9 below the maximum of -1075.5359; there are no restarts for this
  # structure, and the refit from across the bound reaches it.
  m <- poisson_growth(growth_design(5, 365), stratum_slope = FALSE)
  expect_identical(m$longitudinal_info$time_scale, 1000)
  expect_lt(abs(as.numeric(stats::logLik(m$model)) - (-1075.5359)), 0.005)
  expect_gt(lme4::getME(m$model, "theta")[["id.(Intercept)"]], 0.3)
})

test_that("the record flags a fit stopped at a bound where the likelihood still rises", {
  skip_on_cran()
  skip_slow_test()
  # The user's own start is fitted as given: here glmer's own trap in waves.
  d <- growth_design(3)
  trap <- c(0, -0.1156, 0.0812, 0.1206, -0.0185, 0)
  m <- poisson_growth(d, start = list(theta = trap))
  expect_lt(as.numeric(stats::logLik(m$model)), -1088)
  dg <- m$diagnostics
  expect_false(dg$converged)
  expect_identical(dg$boundary_traps$block, "id")
  expect_identical(dg$boundary_traps$term, "(Intercept)")
  expect_gt(dg$boundary_traps$drop, 0.01)
  expect_true(any(grepl("The intercept variance of the 'id' block is at zero", dg$messages, fixed = TRUE)))
  out <- maihda_format_fit_diagnostics(dg)
  expect_true("Convergence problems of the lme4 fit:" %in% out)
  expect_false("Convergence warnings reported by lme4:" %in% out)
  # The same fit read as one outside growth models is flagged too, its sentence naming the
  # block's terms.
  cs <- maihda_fit_diagnostics(m$model)
  expect_false(cs$converged)
  expect_identical(cs$boundary_traps$block, "id")
  expect_true(any(grepl("The intercept variance of the 'id' block is at zero", cs$messages,
                        fixed = TRUE)))
  expect_false(any(grepl("with the slope reversed", cs$messages, fixed = TRUE)))
  expect_false(maihda_fit_diagnostics(m$model, m$longitudinal_info)$converged)
  # The traps themselves: the drop across the bound, and nothing at the maximum.
  tr <- maihda_lme4_boundary_traps(m$model)
  expect_identical(nrow(tr), 1L)
  expect_identical(tr$index, 1L)
  expect_identical(tr$n_terms, 2L)
  ok <- poisson_growth(d)
  expect_identical(nrow(maihda_lme4_boundary_traps(ok$model)), 0L)
  expect_true(ok$diagnostics$converged)
  # Not a glmer fit with its fixed effects optimised: nothing to look at.
  lm_fit <- quiet(lme4::lmer(n ~ wave + (wave | id), data = d))
  expect_identical(nrow(maihda_lme4_boundary_traps(lm_fit)), 0L)
  # An nAGQ = 0 fit at the same bound: its fixed effects are not among the optimised
  # parameters, so it is not measured.
  agq0 <- quiet(lme4::glmer(n ~ wave + (wave | id) + (wave | A:B:C:D), data = d,
                            family = stats::poisson, nAGQ = 0, start = list(theta = trap)))
  expect_gt(length(maihda_lme4_trap_columns(lme4::getME(agq0, "theta"),
                                            lme4::getME(agq0, "cnms"))), 0L)
  expect_identical(nrow(maihda_lme4_boundary_traps(agq0)), 0L)
})

test_that("the record's sentence says what is at the bound", {
  # A block's intercept column: its variance is at zero.
  two <- maihda_lme4_trap_message("id", "(Intercept)", 2L)
  expect_match(two, "The intercept variance of the 'id' block is at zero", fixed = TRUE)
  expect_match(two, "covariance with the slope reversed", fixed = TRUE)
  expect_match(maihda_lme4_trap_message("stratum", "(Intercept)", 3L),
               "covariance with the time terms reversed", fixed = TRUE)
  # The linear column of a quadratic block: the slope is a multiple of the intercept, its
  # variance is not zero.
  slope <- maihda_lme4_trap_message("id", ".maihda_ctime", 3L)
  expect_match(slope, "In the 'id' block the slope is perfectly correlated with the intercept",
               fixed = TRUE)
  expect_false(grepl("is at zero", slope, fixed = TRUE))
  expect_false(grepl("slope variance", slope, fixed = TRUE))
  for (s in c(two, slope)) {
    expect_match(s, "not at the maximum.", fixed = TRUE)
  }
  # The note on a boundary outside the stratum block calls the between-stratum variance
  # unaffected only when the fit is not stopped at a trap: a fit short of the maximum can
  # have every variance off.
  dg <- list(engine = "lme4", singular = TRUE, singular_terms = "id", stratum_singular = FALSE,
             converged = TRUE, messages = character(0))
  expect_true(any(grepl("unaffected", maihda_format_fit_diagnostics(dg), fixed = TRUE)))
  dg$converged <- FALSE
  dg$messages <- two
  dg$boundary_traps <- data.frame(block = "id", term = "(Intercept)", index = 1L, drop = 0.06,
                                  n_terms = 2L)
  out <- maihda_format_fit_diagnostics(dg)
  expect_false(any(grepl("unaffected", out, fixed = TRUE)))
  expect_true("Convergence problems of the lme4 fit:" %in% out)
  expect_true(paste0("  - ", two) %in% out)
})

test_that("the refit across a bound starts there, skips the first stage, and must be lower", {
  skip_on_cran()
  skip_slow_test()
  env <- new.env()
  env$data <- growth_design(3)
  call <- quote(lme4::glmer(n ~ wave + (wave | id) + (wave | A:B:C:D), data = data,
                            family = stats::poisson))
  trap <- c(0, -0.1156, 0.0812, 0.1206, -0.0185, 0)
  trapped_call <- call
  trapped_call$start <- list(theta = trap)
  kept <- maihda_eval_holding_conditions(trapped_call, env)
  th <- lme4::getME(kept$value, "theta")
  expect_lt(abs(unname(th[1])), 1e-4)
  seen <- NULL
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen <<- call
    real(call, env)
  })
  out <- maihda_glmer_mirror_refit(kept, call, env)
  # From the kept fit's estimates with the person intercept entry across its bound...
  expect_equal(seen$start$theta, c(0.05, -th[[2]], unname(th[3:6])))
  expect_identical(seen$start$fixef, lme4::fixef(kept$value))
  # ... by glmer's second stage only ...
  expect_false(eval(seen$control, env)$nAGQ0initStep)
  # ... and kept: it gains most of the 6.37 the trap was short of the maximum.
  expect_gt(as.numeric(stats::logLik(out$value)) - as.numeric(stats::logLik(kept$value)), 6)
  # A user's control is kept, with the first stage switched off.
  env$ctrl <- lme4::glmerControl(optimizer = "Nelder_Mead")
  with_ctrl <- call
  with_ctrl$control <- quote(ctrl)
  maihda_glmer_mirror_refit(kept, with_ctrl, env)
  expect_identical(eval(seen$control, env)$optimizer, c("Nelder_Mead", "Nelder_Mead"))
  expect_false(eval(seen$control, env)$nAGQ0initStep)
  # A refit that is not lower by more than 1e-4 is not returned.
  local_mocked_bindings(maihda_glmer_criterion = function(model) {
    if (identical(model, kept$value)) 2000 else 2000 - 5e-5
  })
  expect_null(maihda_glmer_mirror_refit(kept, call, env))
  # Nothing at a bound: no refit.
  inside_call <- call
  inside_call$start <- list(theta = c(0.358, -0.0377, 0.1276, 0.0319, -0.0296, 0.001))
  inside <- real(inside_call, env)
  expect_length(maihda_lme4_trap_columns(lme4::getME(inside$value, "theta"),
                                         lme4::getME(inside$value, "cnms")), 0L)
  seen <- NULL
  expect_null(maihda_glmer_mirror_refit(inside, call, env))
  expect_null(seen)
})

test_that("a growth block at zero is not flagged when the likelihood does not rise across it", {
  skip_on_cran()
  # A stratum block estimated at 0: its intercept entry sits at the bound beside an entry
  # below it of about 1e-6, and across the bound the deviance does not fall.
  d <- growth_design(8)
  m <- poisson_growth(d)
  tr <- maihda_lme4_boundary_traps(m$model)
  expect_true(all(is.na(tr$drop) | tr$drop <= .maihda_lme4_trap_drop))
  expect_null(m$diagnostics$boundary_traps)
})

test_that("fit_maihda() refits glmer growth models and nothing else", {
  skip_on_cran()
  skip_slow_test()
  calls <- 0L
  real <- maihda_fit_glmer_growth
  local_mocked_bindings(maihda_fit_glmer_growth = function(...) { calls <<- calls + 1L; real(...) })
  d <- growth_design(31)
  poisson_growth(d)
  expect_identical(calls, 1L)
  poisson_growth(d, stratum_slope = FALSE)
  expect_identical(calls, 2L)
  poisson_growth(d, start = list(theta = c(0.3, 0, 0.1, 0.2, 0, 0.05)))
  quiet(fit_maihda(n ~ (1 | A:B:C:D), data = d, family = "poisson"))
  quiet(fit_maihda(log1p(n) ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave"))
  expect_identical(calls, 2L)
})

test_that("a negative-binomial growth fit is refitted through glmer.nb", {
  skip_on_cran()
  skip_slow_test()
  d <- with_seed(11, {
    g <- expand.grid(A = 0:1, B = 0:1)
    u <- stats::rnorm(4, 0, 0.3)
    p <- cbind(stats::rnorm(24, 0, 0.5), stats::rnorm(24, 0, 0.15))
    id <- rep(1:24, each = 4)
    st <- rep(1:4, each = 6)[id]
    t <- rep(0:3, 24)
    data.frame(id = id, wave = t * 365, A = factor(g$A[st]), B = factor(g$B[st]),
               n = stats::rnbinom(96, mu = exp(1 + 0.1 * t + u[st] + p[id, 1] + p[id, 2] * t), size = 3))
  })
  seen <- list()
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen[[length(seen) + 1L]] <<- call
    real(call, env)
  })
  m <- quiet(fit_maihda(n ~ wave + (1 | A:B), data = d, id = "id", time = "wave",
                        family = "negbinomial"))
  expect_identical(m$longitudinal_info$time_scale, 1000)
  expect_match(stats::family(m$model)$family, "Negative Binomial", fixed = TRUE)
  heads <- vapply(seen, function(cl) deparse(cl[[1]]), character(1))
  expect_true(all(heads == "lme4::glmer.nb"))
  expect_gte(length(seen), 3L)
  expect_null(seen[[1]]$start)
  expect_identical(names(seen[[2]]$start), "theta")
  expect_identical(names(seen[[3]]$start), "theta")
})
