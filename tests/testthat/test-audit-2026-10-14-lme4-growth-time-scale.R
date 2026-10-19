# Audit pass 2026-10-14: lme4 growth fits on time counted in fine units.
#
# lme4 optimises a growth model's random-effect covariance on its relative Cholesky scale,
# where a slope entry is per unit of time: counted in days those entries are tiny next to
# the intercept's, and lmer stopped far from the optimum, mostly recording the fit as
# converged. The restarts of maihda_fit_lmer_growth() brought the REML fits of
# maihda_long_data in days to the waves answer, but not the paths that refit the stored
# fit: estimation = "ML" (refitML() stopped 155 and 174 deviance units above the ML
# optimum; PCV 0.92 / 0.35 against 0.96 / 0.92 in waves) and maihda_ic() (78 and 87
# log-likelihood units off). Quadratic growth failed from a span of about 12, in months.
# A Gaussian growth fit lme4 makes with lmer() is now made on time divided by the power of
# ten that puts its span in [1, 10); the stored fit stays on that time and every quantity
# in units of time is converted back. brms growth fits are not scaled; glmer growth fits
# are, with refits of their own (test-audit-2026-10-18-glmer-growth-optimum.R).
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
# stratum and person intercepts, slopes and (deg = 2) curvatures, a Gaussian outcome y and
# a Poisson count n.
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

# ---- the rule ----------------------------------------------------------------------------

test_that("the internal time scale is the power of ten that puts the span in [1, 10)", {
  skip_slow_test()
  sc <- function(...) maihda_longitudinal_scale(c(...))
  # A span already between 1 and 10 is left alone, wherever the axis starts.
  for (s in c(1, 4, 5, 9.99)) expect_identical(sc(0, s), 1)
  expect_identical(sc(10, 14), 1)
  expect_identical(sc(0, 10), 10)
  expect_identical(sc(0, 48), 10)
  expect_identical(sc(0, 1460), 1000)
  expect_identical(sc(0, 1825), 1000)
  expect_identical(sc(2000, 2030), 10)
  expect_identical(sc(0, 0.4), 0.1)
  expect_identical(sc(0, 0.999), 0.1)
  # The span itself settles a span next to a power of ten.
  expect_identical(sc(0, 1000), 1000)
  expect_identical(sc(0, 1000 - 1e-9), 100)
  for (s in 10^seq(-4, 6, by = 0.01)) {
    k <- sc(0, s)
    expect_true(s / k >= 1 && s / k < 10, info = s)
  }
  # Missing values are ignored; no finite span is no scaling.
  expect_identical(sc(0, NA, 1460), 1000)
  expect_identical(sc(3, 3), 1)
  expect_identical(sc(NA, 5), 1)
  expect_identical(maihda_longitudinal_scale(numeric(0)), 1)
  # Stored objects without a scale read as unscaled.
  expect_identical(maihda_lng_time_scale(list(time = "wave")), 1)
  expect_identical(maihda_lng_time_scale(list(time_scale = 1000)), 1000)
  expect_identical(maihda_lng_internal_time(list(time_center = 10, time_scale = 1000), 1010), 1)
  expect_identical(maihda_lng_internal_time(list(time_center = 10), 12), 2)
})

# ---- the conversions ----------------------------------------------------------------------

test_that("a covariance block on internal time is divided by the scale to the power of time", {
  S <- matrix(c(4, 1, 0.5, 1, 2, 0.3, 0.5, 0.3, 1), 3,
              dimnames = list(c("(Intercept)", "t", "I(t^2)"), c("(Intercept)", "t", "I(t^2)")))
  k <- 1000
  expect_equal(maihda_unscale_block(S, k), S / k^outer(0:2, 0:2, `+`))
  expect_identical(dimnames(maihda_unscale_block(S, k)), dimnames(S))
  expect_identical(maihda_unscale_block(S, 1), S)
  # a(t)' Sigma a(t) in user time equals the internal block at the internal time.
  a <- function(t) t^(0:2)
  tu <- 730
  expect_equal(maihda_var_at_time(maihda_unscale_block(S, k), tu),
               maihda_var_at_time(S, tu / k))
  expect_equal(maihda_slope_var_at_time(maihda_unscale_block(S, k), tu),
               maihda_slope_var_at_time(S, tu / k) / k^2)
})

test_that("fixed-effect rows built on the internal time are converted, and named on an axis from 0", {
  lng <- list(time = "wave", time_term = ".maihda_ctime", time_center = 0, time_scale = 1000,
              time_degree = 2)
  terms <- c("(Intercept)", ".maihda_ctime", "I(.maihda_ctime^2)", "A1", "A1:.maihda_ctime",
             "A1:I(.maihda_ctime^2)", "x")
  expect_identical(maihda_lng_coef_time_power(terms, lng), c(0, 1, 2, 0, 1, 2, 0))
  fe <- maihda_fixed_effects_table(terms, estimate = c(1, 2, 3, 4, 5, 6, 7),
                                   se = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.6))
  obj <- list(longitudinal_info = lng)
  out <- maihda_lng_unscale_fixed(fe, obj)
  f <- 1000^-c(0, 1, 2, 0, 1, 2, 0)
  expect_equal(out$estimate, fe$estimate * f)
  expect_equal(out$se, fe$se * f)
  expect_equal(out$lower, fe$lower * f)
  expect_equal(out$upper, fe$upper * f)
  expect_identical(out$statistic, fe$statistic)
  expect_identical(out$p_value, fe$p_value)
  expect_identical(out$term, c("(Intercept)", "wave", "I(wave^2)", "A1", "A1:wave",
                               "A1:I(wave^2)", "x"))
  # A non-syntactic time name is quoted as maihda_time_terms() quotes it.
  lng_q <- lng; lng_q$time <- "time point"
  expect_identical(maihda_lng_unscale_fixed(fe, list(longitudinal_info = lng_q))$term[2:3],
                   c("`time point`", "I(`time point`^2)"))
  # A centered fit keeps the internal name, its lower-order rows being at the baseline.
  lng_c <- lng; lng_c$time_center <- 600
  out_c <- maihda_lng_unscale_fixed(fe, list(longitudinal_info = lng_c))
  expect_identical(out_c$term, terms)
  expect_equal(out_c$estimate, fe$estimate * f)
  # Unscaled: untouched.
  lng_1 <- lng; lng_1$time_scale <- 1
  expect_identical(maihda_lng_unscale_fixed(fe, list(longitudinal_info = lng_1)), fe)
  expect_identical(maihda_lng_unscale_fixed(fe, list(longitudinal_info = NULL)), fe)
})

test_that("a time put into the model frame is divided by the scale", {
  nd <- data.frame(.maihda_ctime = c(0, 0.365), wave = c(600, 965))
  moved <- maihda_longitudinal_set_time(nd, ".maihda_ctime", 730, orig_time = "wave",
                                        center = 600, scale = 1000)
  expect_identical(moved$.maihda_ctime, c(0.73, 0.73))
  expect_identical(moved$wave, c(1330, 1330))
  expect_identical(maihda_longitudinal_set_time(nd, ".maihda_ctime", 2, orig_time = "wave",
                                                center = 600),
                   data.frame(.maihda_ctime = c(2, 2), wave = c(602, 602)))
})

test_that("the reserved internal column is guarded when the time will be scaled", {
  d <- data.frame(pid = rep(1:4, each = 2), day = rep(c(0, 730), 4), y = 1:8,
                  .maihda_ctime = 1)
  expect_error(maihda_validate_longitudinal("pid", "day", 1, d, formula = y ~ (1 | g)),
               "reserved")
  # A glmer fit is scaled as well; a brms fit is not, and this axis starts at 0: nothing
  # would be written.
  expect_error(maihda_validate_longitudinal("pid", "day", 1, d, formula = y ~ (1 | g),
                                            family = "poisson"), "reserved")
  expect_silent(maihda_validate_longitudinal("pid", "day", 1, d, engine = "brms",
                                             formula = y ~ (1 | g)))
  # A span between 1 and 10 from 0 writes no column either.
  d$day <- rep(c(0, 5), 4)
  expect_silent(maihda_validate_longitudinal("pid", "day", 1, d, formula = y ~ (1 | g)))
})

# ---- fits ---------------------------------------------------------------------------------

test_that("maihda_long_data gives the same results in days as in waves", {
  skip_on_cran()
  skip_slow_test()
  data(maihda_long_data, envir = environment())
  f <- wellbeing ~ wave + (1 | gender:ethnicity:education)
  fit_in <- function(k) {
    d <- maihda_long_data
    d$wave <- d$wave * k
    quiet(maihda(f, data = d, id = "id", time = "wave", decomposition = "longitudinal"))
  }
  waves <- fit_in(1)
  suppressWarnings(expect_message(
    fit_maihda(f, data = transform(maihda_long_data, wave = wave * 365), id = "id",
               time = "wave"),
    "the time axis 'wave' spans 1460 time units, so the growth terms are fit on internally scaled time (wave / 1000)",
    fixed = TRUE))
  days <- fit_in(365)
  # Waves span 4: no internal time. Days span 1460: fitted on wave / 1000.
  expect_identical(waves$model$longitudinal_info$time_scale, 1)
  expect_identical(waves$model$longitudinal_info$time_term, "wave")
  expect_false(".maihda_ctime" %in% names(waves$model$data))
  expect_identical(days$model$longitudinal_info$time_scale, 1000)
  expect_identical(days$model$longitudinal_info$time_term, ".maihda_ctime")
  expect_identical(days$model$data$wave, maihda_long_data$wave * 365)

  # The paths that refit the stored fit. An ML deviance does not depend on the unit of
  # time: in days refitML() stopped 155.5 and 173.7 units above it.
  dev_ml <- function(a) c(deviance(lme4::refitML(a$model$model)),
                          deviance(lme4::refitML(a$model_adjusted$model)))
  expect_equal(dev_ml(days), dev_ml(waves), tolerance = 1e-7)
  ml_w <- maihda_longitudinal_pcv(waves$model, waves$model_adjusted, estimation = "ML")
  ml_d <- maihda_longitudinal_pcv(days$model, days$model_adjusted, estimation = "ML")
  expect_equal(ml_w$pcv_intercept, 0.95550, tolerance = 1e-4)
  expect_equal(ml_d$pcv_intercept, ml_w$pcv_intercept, tolerance = 1e-4)
  expect_equal(ml_d$pcv_slope, ml_w$pcv_slope, tolerance = 1e-4)
  ic_w <- quiet(maihda_ic(waves))
  ic_d <- quiet(maihda_ic(days))
  expect_equal(ic_d$logLik, ic_w$logLik, tolerance = 1e-8)
  expect_equal(ic_d$delta, ic_w$delta, tolerance = 1e-5)
  # One model alone is reported on its REML likelihood, which depends on the scale of the
  # time column: a fit of days has log(365) less than the same fit of waves, as an
  # unscaled fit of days has, not the log(1000) more of the internal time.
  ic1_w <- quiet(maihda_ic(waves$model))
  ic1_d <- quiet(maihda_ic(days$model))
  expect_identical(ic1_d$estimator, "REML")
  expect_equal(ic1_d$logLik, ic1_w$logLik - log(365), tolerance = 1e-8)
  expect_equal(ic1_d$AIC, ic1_w$AIC + 2 * log(365), tolerance = 1e-8)
  # The residual-autocorrelation diagnostic gives its time gap in the user's units.
  expect_identical(waves$model$diagnostics$adequacy$autocorrelation$gap, 1)
  expect_identical(days$model$diagnostics$adequacy$autocorrelation$gap, 365)
  expect_equal(days$model$diagnostics$adequacy$autocorrelation$acf1,
               waves$model$diagnostics$adequacy$autocorrelation$acf1, tolerance = 1e-3)

  # Everything in units of time comes back per day.
  expect_equal(days$pcv$pcv_intercept, waves$pcv$pcv_intercept, tolerance = 1e-3)
  expect_equal(days$pcv$pcv_slope, waves$pcv$pcv_slope, tolerance = 1e-3)
  expect_equal(days$pcv$var_slope_null * 365^2, waves$pcv$var_slope_null, tolerance = 1e-3)
  expect_equal(unname(days$pcv$Sigma_stratum_null * 365^outer(0:1, 0:1, `+`)),
               unname(waves$pcv$Sigma_stratum_null), tolerance = 1e-3)
  sw <- summary(waves$model)
  sd <- summary(days$model)
  expect_identical(sd$fixed_effects$term, sw$fixed_effects$term)
  per <- ifelse(sd$fixed_effects$term == "wave", 365, 1)
  expect_equal(sd$fixed_effects$estimate * per, sw$fixed_effects$estimate, tolerance = 1e-3)
  expect_equal(sd$fixed_effects$se * per, sw$fixed_effects$se, tolerance = 1e-3)
  expect_equal(sd$fixed_effects$p_value, sw$fixed_effects$p_value, tolerance = 1e-3)
  vc_per <- ifelse(grepl("slope \\(wave\\)", sd$variance_components$component), 365^2,
                   ifelse(grepl("covariance", sd$variance_components$component), 365, 1))
  expect_equal(sd$variance_components$variance * vc_per, sw$variance_components$variance,
               tolerance = 1e-3)
  expect_identical(sd$longitudinal$vpc_t$time, sw$longitudinal$vpc_t$time * 365)
  expect_equal(sd$longitudinal$vpc_t$estimate, sw$longitudinal$vpc_t$estimate, tolerance = 1e-3)
  expect_equal(sd$longitudinal$vpc_slope, sw$longitudinal$vpc_slope, tolerance = 1e-3)
  st_w <- predict_maihda(waves$model, type = "strata")
  st_d <- predict_maihda(days$model, type = "strata")
  expect_equal(st_d$slope * 365, st_w$slope, tolerance = 1e-3)
  expect_equal(st_d$baseline, st_w$baseline, tolerance = 1e-3)
  nd <- maihda_long_data[c(1:5, 2001:2005), ]
  nd_days <- transform(nd, wave = wave * 365)
  expect_equal(predict_maihda(days$model, newdata = nd_days),
               predict_maihda(waves$model, newdata = nd), tolerance = 1e-4)
  tr_w <- plot(waves$model, type = "trajectories")$data
  tr_d <- plot(days$model, type = "trajectories")$data
  expect_identical(tr_d$time, tr_w$time * 365)
  expect_equal(tr_d$value, tr_w$value, tolerance = 1e-4)
  expect_match(paste(utils::capture.output(print(days$model)), collapse = "\n"),
               "growth terms on internal time wave / 1000", fixed = TRUE)
})

test_that("quadratic growth in months gives the waves answer", {
  skip_on_cran()
  fit_in <- function(k) {
    quiet(maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(52, k = k, deg = 2),
                 id = "id", time = "wave", time_degree = 2, decomposition = "longitudinal"))
  }
  waves <- fit_in(1)
  months <- fit_in(12)
  expect_identical(months$model$longitudinal_info$time_scale, 10)
  # The ML deviance does not depend on the unit of time; in months lmer stopped hundreds
  # of units above it.
  dev_ml <- function(a) c(deviance(lme4::refitML(a$model$model)),
                          deviance(lme4::refitML(a$model_adjusted$model)))
  expect_equal(dev_ml(months), dev_ml(waves), tolerance = 1e-6)
  # The REML fits agree to the optimiser's precision on a quadratic block, which is
  # coarser than for a linear one (PCV_slope 0.489 and 0.486 here).
  expect_equal(months$pcv$pcv_intercept, waves$pcv$pcv_intercept, tolerance = 1e-2)
  expect_equal(months$pcv$pcv_slope, waves$pcv$pcv_slope, tolerance = 1e-2)
  sw <- summary(waves$model)
  sm <- summary(months$model)
  expect_equal(sm$longitudinal$vpc_t$estimate, sw$longitudinal$vpc_t$estimate, tolerance = 1e-2)
  expect_identical(sm$fixed_effects$term, c("(Intercept)", "wave", "I(wave^2)"))
  expect_equal(sm$fixed_effects$estimate * c(1, 12, 144), sw$fixed_effects$estimate,
               tolerance = 1e-3)
  expect_equal(unname(sm$longitudinal$Sigma_stratum * 12^outer(0:2, 0:2, `+`)),
               unname(sw$longitudinal$Sigma_stratum), tolerance = 1e-2)
  # Each stratum's coefficient of time^d is per month^d.
  st_w <- predict_maihda(waves$model, type = "strata")
  st_m <- predict_maihda(months$model, type = "strata")
  expect_equal(st_m$slope * 12, st_w$slope, tolerance = 1e-2)
  expect_equal(st_m$slope2 * 144, st_w$slope2, tolerance = 1e-2)
  tr_w <- plot(waves$model, type = "trajectories")$data
  tr_m <- plot(months$model, type = "trajectories")$data
  expect_equal(tr_m$value, tr_w$value, tolerance = 1e-3)
})

test_that("a glmer growth fit is scaled, and a scaled count fit converts back", {
  skip_on_cran()
  skip_slow_test()
  fit_in <- function(k) {
    quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = growth_design(31, k = k),
                     id = "id", time = "wave", family = "poisson"))
  }
  # glmer growth fits are made on the internal time too, and refitted as
  # maihda_fit_glmer_growth() says.
  days <- fit_in(365)
  expect_identical(days$longitudinal_info$time_scale, 1000)
  expect_identical(days$longitudinal_info$time_term, ".maihda_ctime")
  expect_true(".maihda_ctime" %in% names(days$data))
  expect_true(maihda_longitudinal_scales_time("lme4", "gaussian"))
  expect_true(maihda_longitudinal_scales_time("lme4", stats::gaussian))
  expect_true(maihda_longitudinal_scales_time("lme4", "poisson"))
  expect_true(maihda_longitudinal_scales_time("lme4", "binomial"))
  expect_true(maihda_longitudinal_scales_time("lme4", stats::gaussian(link = "log")))
  expect_false(maihda_longitudinal_scales_time("brms", "gaussian"))
  expect_false(maihda_longitudinal_scales_time("lme4", "no such family"))
  # A count model's level-1 variance is evaluated at each reporting time on the internal
  # time, and its likelihood is the same.
  waves <- fit_in(1)
  sw <- summary(waves)
  sd <- summary(days)
  expect_equal(sd$longitudinal$var_resid_t, sw$longitudinal$var_resid_t, tolerance = 1e-3)
  expect_equal(sd$longitudinal$vpc_t$estimate, sw$longitudinal$vpc_t$estimate, tolerance = 1e-3)
  expect_lt(abs(as.numeric(stats::logLik(days$model)) - as.numeric(stats::logLik(waves$model))),
            0.005)
})

test_that("the bootstrap VPC band of a scaled fit is in the user's time", {
  skip_on_cran()
  skip_slow_test()
  fit_in <- function(k) {
    quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = growth_design(52, k = k), id = "id",
                     time = "wave"))
  }
  waves <- fit_in(1)
  days <- fit_in(365)
  # The refits of the bootstrap are of the stored fit, on the internal time; their blocks
  # are converted back like the fit's own. Same simulated responses on both.
  sw <- with_seed(8, quiet(summary(waves, bootstrap = TRUE, n_boot = 20)))
  sd <- with_seed(8, quiet(summary(days, bootstrap = TRUE, n_boot = 20)))
  expect_equal(sd$longitudinal$vpc_t[, c("estimate", "lower", "upper")],
               sw$longitudinal$vpc_t[, c("estimate", "lower", "upper")], tolerance = 1e-3)
  expect_equal(sd$longitudinal$vpc_slope_ci, sw$longitudinal$vpc_slope_ci, tolerance = 1e-3)
  expect_equal(c(sd$vpc$ci_lower, sd$vpc$ci_upper), c(sw$vpc$ci_lower, sw$vpc$ci_upper),
               tolerance = 1e-3)
  # A count model's bootstrap draws evaluate their level-1 variance on the internal time,
  # as the estimate does: with every draw the fitted model itself, the band is the
  # estimate. (Refits of a glmer growth model differ between time scales by more than the
  # conversions could show.)
  local_mocked_bindings(maihda_refit_draw = function(model, newresp) model)
  pd <- quiet(fit_maihda(n ~ wave + (1 | A:B:C:D), data = growth_design(31, k = 365),
                         id = "id", time = "wave", family = "poisson"))
  expect_identical(pd$longitudinal_info$time_scale, 1000)
  bd <- with_seed(9, quiet(summary(pd, bootstrap = TRUE, n_boot = 10)))
  expect_equal(bd$longitudinal$vpc_t$lower, bd$longitudinal$vpc_t$estimate, tolerance = 1e-12)
  expect_equal(bd$longitudinal$vpc_t$upper, bd$longitudinal$vpc_t$estimate, tolerance = 1e-12)
})

test_that("the scale is read off the analytic rows", {
  skip_on_cran()
  # A late occasion with no outcome widens the span of the input to 50, but not of the
  # rows fitted (0 to 5): no scaling, as no centring.
  d <- growth_design(52)
  extra <- d[d$id == 1, ][1, ]
  extra$wave <- 50
  extra$y <- NA
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = rbind(d, extra), id = "id",
                        time = "wave"))
  expect_identical(m$longitudinal_info$time_scale, 1)
  expect_identical(m$longitudinal_info$time_term, "wave")
  # Fitted rows spanning 50 are scaled.
  d50 <- growth_design(52, k = 10)
  m50 <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d50, id = "id", time = "wave"))
  expect_identical(m50$longitudinal_info$time_scale, 10)
})

test_that("a centered and scaled fit converts back and keeps the internal name", {
  skip_on_cran()
  fit_in <- function(k, at) {
    d <- growth_design(52, k = k)
    d$wave <- d$wave + at
    quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave"))
  }
  waves <- fit_in(1, 0)
  months <- fit_in(12, 600)
  lng <- months$longitudinal_info
  expect_identical(c(lng$time_center, lng$time_scale), c(600, 10))
  expect_identical(lng$ref_time, 600)
  sw <- summary(waves)
  sm <- summary(months)
  expect_identical(sm$fixed_effects$term, c("(Intercept)", ".maihda_ctime"))
  expect_equal(sm$fixed_effects$estimate * c(1, 12), sw$fixed_effects$estimate, tolerance = 1e-3)
  expect_identical(sm$longitudinal$vpc_t$time, sw$longitudinal$vpc_t$time * 12 + 600)
  expect_equal(sm$longitudinal$vpc_t$estimate, sw$longitudinal$vpc_t$estimate, tolerance = 1e-3)
  expect_match(paste(utils::capture.output(print(months)), collapse = "\n"),
               "growth terms on internal time (wave - 600) / 10", fixed = TRUE)
})

test_that("each fitted row keeps its own time, also where dividing and multiplying back does not", {
  skip_on_cran()
  # Days in steps of 1001: 1001 / 1000 * 1000 is not 1001 in floating point.
  d <- growth_design(52, k = 1001)
  expect_false(all(d$wave / 1000 * 1000 == d$wave))
  m <- quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id", time = "wave"))
  expect_identical(m$longitudinal_info$time_scale, 1000)
  expect_identical(m$data$wave, d$wave)
  expect_identical(summary(m)$longitudinal$time_grid, sort(unique(d$wave)))
})

test_that("a brms growth fit is not scaled", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  d <- growth_design(52, k = 365)
  # brms draws from R's random numbers while fitting; with_seed() puts the stream back.
  m <- with_seed(3, quiet(fit_maihda(y ~ wave + (1 | A:B:C:D), data = d, id = "id",
                                     time = "wave", engine = "brms", chains = 1,
                                     iter = 400, warmup = 200, refresh = 0, seed = 3)))
  expect_identical(m$longitudinal_info$time_scale, 1)
  expect_identical(m$longitudinal_info$time_term, "wave")
  expect_true("wave" %in% rownames(brms::fixef(m$model)))
})
