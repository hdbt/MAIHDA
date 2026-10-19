# Audit pass 2026-10-14: lme4 growth fits that stop above the REML optimum.
#
# lmer's default optimiser stopped above the REML optimum on growth models with a linear
# stratum slope: on the wrong side of the stratum block's boundary (the stratum
# intercept-slope correlation with the opposite sign to the optimum's, often the intercept
# SD at exactly 0, so PCV_intercept printed as 100%), and, with time counted in days, far
# from any optimum in every fit measured (maihda_long_data in days: PCV_slope -443962%).
# fit_maihda() now refits such a model from two starts scaled by the person-level block,
# with both signs of the correlation, and keeps the fit with the lowest criterion.
#
# The rank-deficiency verdict of a linear stratum block compared its slope Cholesky entry,
# a slope SD per unit of time, with 1e-4 of the residual SD, so it changed with the unit of
# time; it now compares the slope's SD conditional on the intercept with its own SD.
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

# A growth design with 12 factorial strata, 12 persons each, waves 0..4, fixed additive
# parts in level and slope and a stratum interaction whose slope variance is small next to
# the sampling variance of a stratum's slope -- where lmer's own fit gets stuck.
growth_data <- function(seed) {
  with_seed(seed, {
    g <- expand.grid(A = 1:2, B = 1:2, C = 1:3)
    cc <- sqrt(0.08)
    add0 <- c(-0.3, 0.3)[g$A] + c(-0.2, 0.2)[g$B] + c(-cc, 0, cc)[g$C]
    add1 <- sqrt(0.00625 / 0.2) * (c(-0.2, 0.2)[g$A] + c(-0.3, 0.3)[g$B] + c(0, cc, -cc)[g$C])
    u0 <- stats::rnorm(12, 0, 0.2)
    u1 <- stats::rnorm(12, 0, sqrt(0.00125))
    p0 <- stats::rnorm(144, 0, 0.6)
    p1 <- stats::rnorm(144, 0, 0.2)
    id <- rep(seq_len(144), each = 5)
    wave <- rep(0:4, times = 144)
    st <- rep(seq_len(12), each = 12)[id]
    e <- stats::rnorm(720, 0, 0.5)
    data.frame(id = id, wave = wave,
               A = factor(c("a1", "a2")[g$A[st]]), B = factor(c("b1", "b2")[g$B[st]]),
               C = factor(c("c1", "c2", "c3")[g$C[st]]),
               y = 5 + 0.3 * wave + add0[st] + u0[st] + (add1[st] + u1[st]) * wave +
                 p0[id] + p1[id] * wave + e)
  })
}

growth_cnms <- list(pid = c("(Intercept)", "wave"), stratum = c("(Intercept)", "wave"))
growth_theta <- c(pid.1 = 0.8, pid.2 = 0.06, pid.3 = 0.2, st.1 = 0, st.2 = -0.12, st.3 = 9e-5)

# ---- the starts --------------------------------------------------------------------------

test_that("the restarts start the stratum block at half the person-level SDs, with both signs", {
  starts <- maihda_growth_restart_starts(growth_theta, growth_cnms, id = "pid")
  expect_length(starts, 2L)
  sd_int <- 0.8
  sd_slope <- sqrt(0.06^2 + 0.2^2)
  expect_equal(unname(starts[[1]][4:6]), 0.5 * c(sd_int, sd_slope, sd_slope))
  expect_equal(unname(starts[[2]][4:6]), 0.5 * c(sd_int, -sd_slope, sd_slope))
  # The person-level block starts where lmer's own fit ended, and the names stay.
  for (s in starts) {
    expect_identical(s[1:3], growth_theta[1:3])
    expect_identical(names(s), names(growth_theta))
  }
  # The slices follow the order of the terms, whatever it is.
  rev_cnms <- rev(growth_cnms)
  rev_theta <- growth_theta[c(4:6, 1:3)]
  rs <- maihda_growth_restart_starts(rev_theta, rev_cnms, id = "pid")
  expect_equal(unname(rs[[2]][1:3]), 0.5 * c(sd_int, -sd_slope, sd_slope))
  expect_identical(rs[[2]][4:6], rev_theta[4:6])
  # Another grouping factor between them is left alone.
  ctx_cnms <- list(pid = growth_cnms$pid, ctx = "(Intercept)", stratum = growth_cnms$stratum)
  ctx_theta <- c(growth_theta[1:3], ctx = 0.3, growth_theta[4:6])
  cs <- maihda_growth_restart_starts(ctx_theta, ctx_cnms, id = "pid")
  expect_identical(unname(cs[[1]][4]), 0.3)
  expect_equal(unname(cs[[1]][5:7]), 0.5 * c(sd_int, sd_slope, sd_slope))
  # The starts are in the units of time of the fit: a slope per day is a slope per year
  # divided by 365, and so is the start.
  days <- growth_theta
  days[c(2, 3, 5, 6)] <- days[c(2, 3, 5, 6)] / 365
  ds <- maihda_growth_restart_starts(days, growth_cnms, id = "pid")
  expect_equal(unname(ds[[1]][5:6]), unname(starts[[1]][5:6]) / 365)
  expect_equal(unname(ds[[1]][4]), unname(starts[[1]][4]))
})

test_that("no restart is offered where there is nothing to scale it by", {
  none <- function(theta, cnms, id = "pid") {
    expect_identical(maihda_growth_restart_starts(theta, cnms, id = id), list())
  }
  # An intercept-only stratum level, a quadratic stratum block.
  none(c(growth_theta[1:3], st = 0.3), list(pid = growth_cnms$pid, stratum = "(Intercept)"))
  q_cnms <- list(pid = c("(Intercept)", "t", "I(t^2)"), stratum = c("(Intercept)", "t", "I(t^2)"))
  none(rep(0.3, 12), q_cnms)
  # No person-level block of that name, or one without a slope.
  none(growth_theta, growth_cnms, id = "person")
  none(c(pid = 0.8, growth_theta[4:6]), list(pid = "(Intercept)", stratum = growth_cnms$stratum))
  # A person-level SD of zero or a non-finite one.
  z <- growth_theta; z[1] <- 0; none(z, growth_cnms)
  z <- growth_theta; z[2:3] <- 0; none(z, growth_cnms)
  z <- growth_theta; z[3] <- NA; none(z, growth_cnms)
  # theta and the terms do not match, or are not there.
  none(growth_theta[-1], growth_cnms)
  none(NULL, growth_cnms)
  none(growth_theta, unname(growth_cnms))
  none(growth_theta, NULL)
})

# ---- which fit is kept --------------------------------------------------------------------

test_that("lmer's own fit is kept unless another is lower by more than 1e-4", {
  expect_identical(maihda_pick_lowest(c(10, 10 - 5e-5, 10 + 1)), 1L)
  expect_identical(maihda_pick_lowest(c(10, 10 - 2e-4)), 2L)
  expect_identical(maihda_pick_lowest(c(10, 9.5, 8)), 3L)
  expect_identical(maihda_pick_lowest(c(10, NA, 9)), 3L)
  expect_identical(maihda_pick_lowest(c(10, NA, NA)), 1L)
  expect_identical(maihda_pick_lowest(10), 1L)
  expect_identical(maihda_pick_lowest(c(NA, 12, 11)), 3L)
  expect_identical(maihda_pick_lowest(c(NA, NA)), 1L)
  expect_identical(maihda_pick_lowest(numeric(0)), NA_integer_)
  expect_identical(maihda_pick_lowest(c(10, 9.9, 9.95), tol = 0.2), 1L)
})

test_that("only the kept fit's warnings and messages are signalled", {
  env <- new.env()
  env$fake_fit <- function(start = NULL) {
    if (is.null(start)) {
      warning("first fit warns")
      message("first fit says")
      list(crit = 10)
    } else {
      message("restart ", start[["a"]])
      list(crit = 10 - start[["a"]])
    }
  }
  local_mocked_bindings(
    maihda_growth_restart_starts = function(theta, cnms, id, scale = 0.5) list(c(a = 0.5), c(a = 2)),
    maihda_lmer_criterion = function(model) model$crit
  )
  # Every condition that reaches the caller, in order.
  heard <- function(expr) {
    seen <- character(0)
    value <- withCallingHandlers(expr,
      warning = function(w) { seen <<- c(seen, paste("warning:", conditionMessage(w))); invokeRestart("muffleWarning") },
      message = function(m) { seen <<- c(seen, paste("message:", conditionMessage(m))); invokeRestart("muffleMessage") })
    list(value = value, seen = seen)
  }
  h <- heard(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid"))
  expect_identical(h$value$crit, 8)
  expect_identical(h$seen, "message: restart 2\n")

  # When no restart is better, lmer's own fit and its own conditions come back.
  env$fake_fit <- function(start = NULL) {
    if (is.null(start)) {
      warning("first fit warns")
      list(crit = 1)
    } else {
      warning("restart warns")
      list(crit = 5)
    }
  }
  h <- heard(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid"))
  expect_identical(h$value$crit, 1)
  expect_identical(h$seen, "warning: first fit warns")

  # A restart better by less than 1e-4 is the same optimum reached by another path:
  # lmer's own fit stays. One better by 0.01 -- the size of most real gains -- is kept.
  env$fake_fit <- function(start = NULL) {
    if (is.null(start)) list(crit = 10) else list(crit = 10 - 5e-5 * start[["a"]])
  }
  expect_identical(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid")$crit, 10)
  env$fake_fit <- function(start = NULL) {
    if (is.null(start)) list(crit = 10) else list(crit = 10 - 0.005 * start[["a"]])
  }
  expect_equal(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid")$crit, 9.99)

  # Both starts are tried, whichever of them wins.
  local_mocked_bindings(
    maihda_growth_restart_starts = function(theta, cnms, id, scale = 0.5) list(c(a = 2), c(a = 0.5)))
  env$fake_fit <- function(start = NULL) if (is.null(start)) list(crit = 10) else list(crit = 10 - start[["a"]])
  expect_identical(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid")$crit, 8)

  # A restart that fails is passed over; an error of lmer's own fit is not caught.
  env$fake_fit <- function(start = NULL) {
    if (is.null(start)) list(crit = 3) else if (start[["a"]] < 1) stop("no fit") else list(crit = 1)
  }
  expect_identical(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid")$crit, 1)
  env$fake_fit <- function(start = NULL) stop("the first fit fails")
  expect_error(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid"), "the first fit fails")
})

test_that("fits are compared on the criterion lmer minimised", {
  skip_on_cran()
  reml <- suppressMessages(lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy))
  ml <- suppressMessages(lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy,
                                    REML = FALSE))
  expect_identical(maihda_lmer_criterion(reml), lme4::REMLcrit(reml))
  expect_identical(maihda_lmer_criterion(ml), stats::deviance(ml))
  # The two criteria of one REML fit are not the same number.
  expect_gt(abs(lme4::REMLcrit(reml) - stats::deviance(reml, REML = FALSE)), 1)
})

test_that("warnings and messages held back are signalled unchanged and in order", {
  held <- maihda_eval_holding_conditions(
    quote({ message("one"); warning("two"); message("three"); 42 }), new.env())
  expect_identical(held$value, 42)
  expect_identical(vapply(held$conditions, function(x) class(x)[2], ""),
                   c("message", "warning", "message"))
  seen <- character(0)
  withCallingHandlers(
    maihda_signal_conditions(held$conditions),
    warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_identical(seen, c("one\n", "two", "three\n"))
})

# ---- which fits are restarted ---------------------------------------------------------------

test_that("lmer growth fits with a linear stratum slope and no user start are restarted", {
  lng <- list(id = "id", stratum_slope = TRUE, time_degree = 1)
  expect_true(maihda_growth_restarts_apply(quote(lme4::lmer), lng, list()))
  expect_true(maihda_growth_restarts_apply(quote(lme4::lmer), lng, list(REML = FALSE)))
  expect_false(maihda_growth_restarts_apply(quote(lme4::lmer), lng, list(start = 1)))
  expect_false(maihda_growth_restarts_apply(quote(lme4::glmer), lng, list()))
  expect_false(maihda_growth_restarts_apply(quote(lme4::lmer), NULL, list()))
  expect_false(maihda_growth_restarts_apply(quote(lme4::lmer),
                                            utils::modifyList(lng, list(stratum_slope = FALSE)), list()))
  # A quadratic one is restarted too, from other starts (test-audit-2026-10-18-lme4-growth-
  # quadratic-restarts.R).
  expect_true(maihda_growth_restarts_apply(quote(lme4::lmer),
                                           utils::modifyList(lng, list(time_degree = 2)), list()))
})

test_that("fit_maihda() takes the restart path for exactly those fits", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- maihda_long_data[maihda_long_data$id %in% unique(maihda_long_data$id)[1:120], ]
  calls <- 0L
  real <- maihda_fit_lmer_growth
  local_mocked_bindings(maihda_fit_lmer_growth = function(...) { calls <<- calls + 1L; real(...) })
  f <- wellbeing ~ wave + (1 | gender:ethnicity:education)
  q <- function(e) suppressWarnings(suppressMessages(e))
  q(fit_maihda(f, data = d, id = "id", time = "wave"))
  expect_identical(calls, 1L)
  q(fit_maihda(f, data = d, id = "id", time = "wave", REML = FALSE))
  expect_identical(calls, 2L)
  q(fit_maihda(f, data = d, id = "id", time = "wave", start = c(0.8, 0, 0.2, 1, 0.2, 0.25)))
  q(fit_maihda(f, data = d, id = "id", time = "wave", stratum_slope = FALSE))
  q(fit_maihda(low_wellbeing ~ wave + (1 | gender:ethnicity:education), data = d, id = "id",
               time = "wave", family = "binomial"))
  q(fit_maihda(wellbeing ~ age + (1 | gender:ethnicity:education), data = d))
  expect_identical(calls, 2L)
  q(fit_maihda(f, data = d, id = "id", time = "wave", time_degree = 2))
  expect_identical(calls, 3L)
})

# ---- the fits it changes ------------------------------------------------------------------

test_that("a fit stopped on the wrong side of the boundary is refitted to the optimum", {
  skip_on_cran()
  d <- growth_data(10075)
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ wave + (1 | A:B:C), data = d, id = "id", time = "wave",
           decomposition = "longitudinal")))
  # The exact REML optimum of this balanced design gives 1531.844 for the adjusted model
  # and a PCV of 0.73845 / 0.67861; lmer's own fit stopped at 1535.574 with the stratum
  # intercept SD at exactly 0, which printed PCV_intercept as 100%.
  expect_equal(lme4::REMLcrit(a$model_adjusted$model), 1531.8440, tolerance = 1e-7)
  expect_equal(a$pcv$pcv_intercept, 0.73845, tolerance = 1e-3)
  expect_equal(a$pcv$pcv_slope, 0.67861, tolerance = 1e-3)
  expect_gt(lme4::getME(a$model_adjusted$model, "theta")[["stratum.(Intercept)"]], 0.3)
})

test_that("the PCV of maihda_long_data does not depend on the unit of time", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  f <- wellbeing ~ wave + (1 | gender:ethnicity:education)
  fit_in <- function(k) {
    d <- maihda_long_data
    d$wave <- d$wave * k
    suppressWarnings(suppressMessages(
      maihda(f, data = d, id = "id", time = "wave", decomposition = "longitudinal")))
  }
  waves <- fit_in(1)
  days <- fit_in(365)
  # In days lmer's own fits stopped far from the optimum (PCV_intercept -0.15, PCV_slope
  # -4440); the restarts reach it.
  expect_equal(days$pcv$pcv_intercept, waves$pcv$pcv_intercept, tolerance = 1e-3)
  expect_equal(days$pcv$pcv_slope, waves$pcv$pcv_slope, tolerance = 1e-3)
  expect_equal(waves$pcv$pcv_intercept, 0.9186, tolerance = 1e-3)
  # A variance per day squared is the variance per wave divided by 365^2.
  expect_equal(days$pcv$var_slope_null * 365^2, waves$pcv$var_slope_null, tolerance = 1e-3)
  for (fld in c("null_rank_deficient", "null_slope_at_boundary", "adjusted_rank_deficient")) {
    expect_identical(days$pcv[[fld]], waves$pcv[[fld]])
  }
})

# ---- the rank-deficiency verdict ----------------------------------------------------------

test_that("a linear block's boundary verdict is a ratio of SDs, not a slope SD per unit of time", {
  v <- function(l) unname(maihda_linear_block_boundary(l))
  # A perfect correlation (lme4's slope entry at exactly 0), and one lmer stopped just short
  # of: the slope's conditional SD is 0 and 0.3% of its own SD.
  expect_identical(v(c(0.4, 0.1, 0)), c(FALSE, TRUE))
  expect_identical(v(c(0.4, 0.1, 3e-4)), c(FALSE, TRUE))
  # An interior block: conditional SD 85% of the slope SD.
  expect_identical(v(c(0.4, 0.05, 0.08)), c(FALSE, FALSE))
  # The same blocks with time in finer units: the slope entries shrink, the verdict does
  # not -- where the old rule (an entry below 1e-4) called the interior block singular.
  for (k in c(12, 365, 3650)) {
    expect_identical(v(c(0.4, 0.1 / k, 3e-4 / k)), c(FALSE, TRUE))
    expect_identical(v(c(0.4, 0.05 / k, 0.08 / k)), c(FALSE, FALSE))
  }
  expect_lt(0.08 / 3650, 1e-4)
  # The cut-off: a conditional SD of 1% of the slope SD.
  expect_identical(v(c(0.4, 0.1, 0.1 * 0.0099 / sqrt(1 - 0.0099^2))), c(FALSE, TRUE))
  expect_identical(v(c(0.4, 0.1, 0.1 * 0.0101 / sqrt(1 - 0.0101^2))), c(FALSE, FALSE))
  # With the intercept SD at zero the slope variance is identified on its own: the slope
  # is at the boundary only if its SD is exactly zero.
  expect_identical(v(c(0, -0.12, 9e-5)), c(TRUE, FALSE))
  expect_identical(v(c(5e-5, 0.12, 0)), c(TRUE, FALSE))
  expect_identical(v(c(0, 0, 0)), c(TRUE, TRUE))
  expect_identical(v(c(0.4, 0, 0)), c(FALSE, TRUE))
  expect_named(maihda_linear_block_boundary(c(0.4, 0.1, 0)), c("intercept", "slope"))
})

test_that("the verdict on a fitted block reads the stratum term, and keeps lme4's rule for other sizes", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  f <- wellbeing ~ wave + (1 | gender:ethnicity:education)
  q <- function(e) suppressWarnings(suppressMessages(e))
  m <- q(fit_maihda(f, data = maihda_long_data, id = "id", time = "wave"))
  l <- maihda_re_theta_block_lme4(m$model, "stratum")
  expect_length(l, 3L)
  expect_identical(l, unname(lme4::getME(m$model, "theta")[c("stratum.(Intercept)",
                                                             "stratum.wave.(Intercept)",
                                                             "stratum.wave")]))
  expect_identical(maihda_stratum_block_boundary_lme4(m$model),
                   stats::setNames(maihda_linear_block_boundary(l), c("(Intercept)", "wave")))
  # Pin a perfect correlation on the fitted block: the slope component is at the boundary.
  pin <- function(l) {
    x <- m$model
    x@theta[match(c("stratum.(Intercept)", "stratum.wave.(Intercept)", "stratum.wave"),
                  names(lme4::getME(m$model, "theta")))] <- l
    x
  }
  pinned <- pin(c(0.4, 0.1, 0))
  expect_true(maihda_stratum_slope_at_boundary_lme4(pinned))
  expect_true(maihda_stratum_growth_rank_deficient_lme4(pinned))
  expect_false(maihda_stratum_growth_rank_deficient_lme4(m$model))
  # Where the two rules part: a slope entry lme4 left at 3e-4 next to a slope SD of 0.1
  # (conditional SD 0.3% of it) is at the boundary; an interior block with time in tenths
  # of a day (slope entries far below 1e-4) is not.
  expect_true(maihda_stratum_slope_at_boundary_lme4(pin(c(0.4, 0.1, 3e-4))))
  expect_true(maihda_stratum_growth_rank_deficient_lme4(pin(c(0.4, 0.1, 3e-4))))
  expect_false(maihda_stratum_slope_at_boundary_lme4(pin(c(0.4, 0.05, 0.08) / c(1, 3650, 3650))))
  expect_false(maihda_stratum_growth_rank_deficient_lme4(pin(c(0.4, 0.05, 0.08) / c(1, 3650, 3650))))
  # With the intercept SD at zero the block is rank-deficient but the slope is not at the
  # boundary.
  expect_true(maihda_stratum_growth_rank_deficient_lme4(pin(c(0, -0.12, 9e-5))))
  expect_false(maihda_stratum_slope_at_boundary_lme4(pin(c(0, -0.12, 9e-5))))
  # A quadratic block keeps lme4's rule, entry by entry.
  m2 <- q(fit_maihda(f, data = maihda_long_data, id = "id", time = "wave", time_degree = 2))
  expect_identical(maihda_stratum_block_boundary_lme4(m2$model),
                   maihda_re_boundary_components_lme4(m2$model)[["stratum"]])
  expect_null(maihda_re_theta_block_lme4(m$model, "no such term"))
  expect_null(maihda_stratum_block_boundary_lme4(list()))
})
