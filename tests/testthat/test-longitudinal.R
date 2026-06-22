# Longitudinal / growth-curve MAIHDA: time-varying VPC and additive-vs-
# multiplicative PCV (Bell, Evans, Holman & Leckie 2024).

# ---- pure helpers (no model fit) -------------------------------------------

test_that("maihda_var_at_time evaluates a(t)' Sigma a(t)", {
  Sigma <- matrix(c(2, 0.1, 0.1, 0.5), nrow = 2)  # v0=2, cov=0.1, v1=0.5
  # at t = 3: 2 + 2*3*0.1 + 9*0.5 = 7.1
  expect_equal(maihda_var_at_time(Sigma, 3), 7.1)
  expect_equal(maihda_var_at_time(Sigma, 0), 2)            # intercept variance
  expect_equal(maihda_var_at_time(Sigma, c(0, 3)), c(2, 7.1))
})

test_that("maihda_longitudinal_formula builds the 3-level growth structure", {
  f <- maihda_longitudinal_formula(y ~ x + (1 | stratum), id = "pid",
                                   time = "wave", time_degree = 1)
  bars <- vapply(reformulas::findbars(f),
                 function(b) paste(deparse(b), collapse = ""), character(1))
  expect_true(any(grepl("wave \\| pid", bars)))
  expect_true(any(grepl("wave \\| stratum", bars)))
  # the time term enters the fixed part
  expect_true("wave" %in% attr(stats::terms(reformulas::nobars(f)), "term.labels"))

  f2 <- maihda_longitudinal_formula(y ~ (1 | stratum), id = "pid",
                                    time = "t", time_degree = 2)
  bars2 <- vapply(reformulas::findbars(f2),
                  function(b) paste(deparse(b), collapse = ""), character(1))
  expect_true(any(grepl("I\\(t\\^2\\) \\| stratum", bars2)))
})

test_that("maihda_re_cov_draws_brms builds the 2x2 block from draws (Stan-free)", {
  # Hand-built posterior draws: SD and correlation columns in brms' naming.
  draws <- data.frame(
    sd_stratum__Intercept = c(2, 4),
    sd_stratum__wave = c(0.5, 1),
    cor_stratum__Intercept__wave = c(0.5, -0.5)
  )
  blk <- maihda_re_cov_draws_brms(draws, "stratum", "wave")
  expect_equal(blk$v0, c(4, 16))                 # sd0^2
  expect_equal(blk$v1, c(0.25, 1))               # sd1^2
  expect_equal(blk$cov, c(0.5 * 2 * 0.5, -0.5 * 4 * 1))  # cor * sd0 * sd1
  expect_error(maihda_re_cov_draws_brms(draws, "id", "wave"), "Could not find")
})

test_that("longitudinal components table is honest about intercept-vs-baseline and lists every covariance", {
  # Quadratic (3x3) stratum block; linear (2x2) individual block. Pure helper.
  Ss <- matrix(c(2,    0.3,  0.05,
                 0.3,  0.5,  0.02,
                 0.05, 0.02, 0.1), nrow = 3, byrow = TRUE)
  Si <- matrix(c(1, 0.1, 0.1, 0.4), nrow = 2)
  tab <- maihda_longitudinal_components_table(Ss, Si, var_resid = 0.7,
                                              time = "wave", id = "pid")

  # The intercept variance is the time-0 quantity, NOT the baseline (ref_time).
  expect_true("Between-stratum: intercept (time = 0)" %in% tab$component)
  expect_false(any(grepl("baseline", tab$component)))

  # A quadratic block contributes ALL THREE off-diagonal covariances (not just
  # intercept-slope), each carrying the corresponding Sigma cell.
  expect_equal(tab$variance[tab$component == "Between-stratum: intercept-slope covariance"],
               Ss[1, 2])
  expect_equal(tab$variance[tab$component == "Between-stratum: intercept-slope^2 covariance"],
               Ss[1, 3])
  expect_equal(tab$variance[tab$component == "Between-stratum: slope-slope^2 covariance"],
               Ss[2, 3])

  # The linear (2x2) individual block still yields exactly one covariance row.
  expect_equal(
    sum(grepl("^Between-individual \\(pid\\): .*covariance$", tab$component)), 1L)
})

test_that("maihda_validate_longitudinal enforces its contract", {
  d <- data.frame(pid = rep(1:3, each = 2), t = rep(0:1, 3), y = rnorm(6))
  expect_error(maihda_validate_longitudinal(NULL, "t", 1, d), "needs 'id'")
  expect_error(maihda_validate_longitudinal("pid", "missing", 1, d), "not found")
  # not longitudinal: every id unique
  d2 <- data.frame(pid = 1:6, t = 0:5, y = rnorm(6))
  expect_error(maihda_validate_longitudinal("pid", "t", 1, d2), "not look longitudinal")
  # unsupported engine / weights / context
  expect_error(maihda_validate_longitudinal("pid", "t", 1, d, engine = "wemix"),
               "lme4")
  expect_error(maihda_validate_longitudinal("pid", "t", 1, d,
                                            sampling_weights = "w"), "design-weighted")
  expect_error(maihda_validate_longitudinal("pid", "t", 2, d, engine = "brms"),
               "linear growth only")
})

# ---- fitted-model tests (lme4) ---------------------------------------------

skip_on_cran()

data(maihda_long_data, package = "MAIHDA")

m_g <- fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                  data = maihda_long_data, id = "id", time = "wave")

# Shared null/adjusted longitudinal decomposition, fit once and reused.
a_g <- maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
              data = maihda_long_data, id = "id", time = "wave",
              decomposition = "longitudinal")

test_that("fit tags the model and builds the growth formula", {
  expect_s3_class(m_g, "maihda_model")
  expect_false(is.null(m_g$longitudinal_info))
  expect_identical(m_g$longitudinal_info$time, "wave")
  bars <- vapply(reformulas::findbars(m_g$formula),
                 function(b) paste(deparse(b), collapse = ""), character(1))
  expect_true(any(grepl("wave \\| id", bars)))
  expect_true(any(grepl("wave \\| stratum", bars)))
})

test_that("summary reports a time-varying VPC", {
  s <- summary(m_g)
  expect_false(is.null(s$longitudinal))
  vt <- s$longitudinal$vpc_t
  expect_true(all(c("time", "estimate") %in% names(vt)))
  expect_true(all(vt$estimate >= 0 & vt$estimate <= 1))
  # the headline VPC equals VPC at the reference (baseline) time
  ref_row <- vt$estimate[vt$time == s$longitudinal$ref_time]
  expect_equal(s$vpc$estimate, ref_row, tolerance = 1e-8)
  # stratum slope variance is identified (> 0): the injected trajectory differences
  expect_gt(s$longitudinal$Sigma_stratum[2, 2], 0)
})

test_that("longitudinal PCV recovers a mostly-additive trajectory split", {
  a <- a_g
  expect_identical(a$mode, "longitudinal")
  expect_s3_class(a$pcv, "maihda_long_pcv")
  # both PCVs are genuine proportions strictly inside (0, 1) by construction
  expect_gt(a$pcv$pcv_intercept, 0.5)
  expect_lt(a$pcv$pcv_intercept, 1)
  expect_gt(a$pcv$pcv_slope, 0.5)       # trajectories mostly additive
  expect_lt(a$pcv$pcv_slope, 1)         # but a multiplicative residual survives
  # the adjusted model retains some stratum slope variance (the interaction)
  expect_gt(a$pcv$Sigma_stratum_adjusted[2, 2], 0)
})

test_that("longitudinal PCV baseline is the variance at ref_time, not raw time 0", {
  # Shift time off zero (waves 10..14): the baseline PCV must be the PCV of the
  # between-stratum variance AT the observed baseline (ref_time = 10), evaluated
  # via a(t)'Sigma a(t), not the raw time-0 intercept-variance cell Sn[1, 1].
  d <- maihda_long_data
  d$wave <- d$wave + 10
  # Fitting on raw time far from zero stresses lme4's optimizer (the time-0
  # intercept variance is a far extrapolation); the convergence notice is
  # immaterial here -- the assertions below are algebraic identities on whatever
  # covariance block is returned.
  a <- suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = d, id = "id", time = "wave", decomposition = "longitudinal"))
  pcv <- a$pcv
  expect_identical(pcv$ref_time, 10)

  Sn <- pcv$Sigma_stratum_null
  Sa <- pcv$Sigma_stratum_adjusted
  v_base_n <- maihda_var_at_time(Sn, 10)
  v_base_a <- maihda_var_at_time(Sa, 10)
  expect_equal(pcv$var_baseline_null, v_base_n)
  expect_equal(pcv$var_baseline_adjusted, v_base_a)
  expect_equal(pcv$pcv_intercept, (v_base_n - v_base_a) / v_base_n)

  # It must NOT equal the (meaningless) raw time-0 cell PCV when time is off zero.
  raw_cell_pcv <- (Sn[1, 1] - Sa[1, 1]) / Sn[1, 1]
  expect_false(isTRUE(all.equal(pcv$pcv_intercept, raw_cell_pcv)))

  # The print method reports the baseline at ref_time (= 10), not time 0.
  expect_output(print(pcv), "Baseline \\(wave = 10\\)")
})

test_that("longitudinal models refuse cross-sectional scalar rankings/plots", {
  # A growth model's stratum estimand is a trajectory, so the scalar BLUP views
  # are not defined and must redirect to the trajectory tools.
  for (ty in c("predicted", "obs_vs_shrunken", "risk_vs_effect",
               "effect_decomp", "prediction_deviation", "ternary")) {
    expect_error(plot(m_g, type = ty), "longitudinal MAIHDA")
  }
  expect_error(maihda_strata_ranking(m_g, summary(m_g)), "longitudinal MAIHDA")

  # maihda_table omits the ranking and explains why (rather than silently showing
  # a misleading cross-sectional rank).
  tab <- maihda_table(a_g)
  expect_null(tab$strata)
  expect_match(tab$strata_note, "trajector")
  expect_output(print(tab), "Strata are trajectories")
})

test_that("predict(type = 'strata') returns trajectory parameters", {
  ps <- predict_maihda(m_g, type = "strata")
  expect_true(all(c("stratum", "baseline", "intercept", "slope") %in% names(ps)))
  expect_equal(nrow(ps), nrow(m_g$strata_info))
  # baseline = a(ref_time)' coef = intercept + slope*ref_time for a linear model.
  ref <- m_g$longitudinal_info$ref_time
  expect_equal(ps$baseline, ps$intercept + ps$slope * ref, tolerance = 1e-8)
})

test_that("predict(type = 'strata') baseline differs from the raw intercept off zero", {
  # With waves shifted to 10.., the baseline deviation (at ref_time = 10) is NOT
  # the raw time-0 intercept -- the column must reflect ref_time, not time 0.
  d <- maihda_long_data
  d$wave <- d$wave + 10
  m <- suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = d, id = "id", time = "wave"))
  ps <- predict_maihda(m, type = "strata")
  expect_identical(m$longitudinal_info$ref_time, 10)
  expect_equal(ps$baseline, ps$intercept + ps$slope * 10, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(ps$baseline, ps$intercept)))
})

test_that("maihda_table reports the baseline between-stratum variance for a longitudinal fit", {
  tab <- maihda_table(a_g)
  bv <- tab$models$estimate[tab$models$statistic == "Between-stratum variance"]
  expect_length(bv, 1)
  expect_true(is.finite(bv))   # was NA before: the cross-sectional row label never matched
  # It is the between-stratum variance at ref_time, matching the VPC anchor.
  s <- a_g$summary
  expect_equal(bv,
               as.numeric(maihda_var_at_time(s$longitudinal$Sigma_stratum,
                                             s$longitudinal$ref_time)),
               tolerance = 1e-8)
})

test_that("plots return ggplot objects", {
  expect_s3_class(plot(m_g, type = "vpc_trajectory"), "ggplot")
  expect_s3_class(plot(m_g, type = "trajectories"), "ggplot")
  expect_s3_class(plot(a_g, type = "pcv_trajectory"), "ggplot")
  expect_s3_class(plot(a_g, type = "vpc_trajectory"), "ggplot")
  expect_s3_class(plot(a_g, type = "trajectories"), "ggplot")
})

test_that("print methods cover the longitudinal branches", {
  expect_output(print(m_g), "Longitudinal")
  expect_output(print(summary(m_g)), "baseline")
  expect_output(print(a_g), "longitudinal")
  expect_output(print(a_g$pcv), "PCV")
})

test_that("plot(type = 'all') dispatches the longitudinal trajectory views", {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  pm <- plot(m_g)            # longitudinal model -> vpc_trajectory + trajectories
  pa <- plot(a_g)            # longitudinal analysis -> the three trajectory views
  expect_type(pm, "list")
  expect_true("vpc_trajectory" %in% names(pm))
  expect_type(pa, "list")
  expect_true("pcv_trajectory" %in% names(pa))
})

test_that("bootstrap gives a VPC-trajectory ribbon", {
  s <- summary(m_g, bootstrap = TRUE, n_boot = 10)
  expect_true(is.finite(s$vpc$ci_lower) && is.finite(s$vpc$ci_upper))
  expect_true(any(is.finite(s$longitudinal$vpc_t$lower)))
})

test_that("maihda strips user-written dimension main effects from the null", {
  a2 <- maihda(
    wellbeing ~ wave + gender + ethnicity + education +
      (1 | gender:ethnicity:education),
    data = maihda_long_data, id = "id", time = "wave",
    decomposition = "longitudinal")
  expect_identical(a2$mode, "longitudinal")
  expect_s3_class(a2$pcv, "maihda_long_pcv")
  # the null model carries no dimension main effects (they belong to the adjusted)
  expect_false(any(c("gender", "ethnicity", "education") %in%
                     attr(stats::terms(reformulas::nobars(a2$formula)), "term.labels")))
})

test_that("compare_maihda warns and vpc_trajectory errors off a longitudinal model", {
  expect_warning(compare_maihda(m_g, m_g), "time-varying")
  strata <- make_strata(maihda_long_data,
                        vars = c("gender", "ethnicity", "education"))
  m_cs <- fit_maihda(wellbeing ~ 1 + (1 | stratum), data = strata$data)
  expect_error(plot(m_cs, type = "vpc_trajectory"), "longitudinal")
})

test_that("brms longitudinal path gives a time-varying VPC with credible bands", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:120])
  m <- suppressWarnings(suppressMessages(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = d,
               id = "id", time = "wave", engine = "brms",
               chains = 2, iter = 600, refresh = 0, seed = 1)))
  s <- summary(m)
  expect_false(is.null(s$longitudinal))
  expect_identical(s$vpc$method, "posterior")
  expect_true(all(is.finite(s$longitudinal$vpc_t$estimate)))
  expect_true(any(is.finite(s$longitudinal$vpc_t$lower)))   # credible band

  a <- suppressWarnings(suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = d,
           id = "id", time = "wave", engine = "brms",
           decomposition = "longitudinal", chains = 2, iter = 600,
           refresh = 0, seed = 1)))
  expect_s3_class(a$pcv, "maihda_long_pcv")
  expect_true(is.finite(a$pcv$pcv_slope))
})

test_that("maihda_ic works on a longitudinal lme4 fit", {
  ic <- maihda_ic(m_g)
  expect_true(is.data.frame(ic) || is.list(ic))
})

test_that("binomial longitudinal fit gives a latent-scale time-varying VPC", {
  mb <- fit_maihda(low_wellbeing ~ wave + (1 | gender:ethnicity:education),
                   data = maihda_long_data, id = "id", time = "wave",
                   family = "binomial")
  s <- summary(mb)
  expect_false(is.null(s$longitudinal))
  expect_true(all(is.finite(s$longitudinal$vpc_t$estimate)))
})

# ---- edge cases: higher time_degree, singular fits, explicit-times PCV ------
# Shared fits for the edge-case block; each is reused across the tests below.

# Quadratic (time_degree = 2) growth: a 3x3 (intercept, slope, slope^2)
# covariance block at each level instead of the linear 2x2.
m_q <- suppressWarnings(suppressMessages(
  fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
             data = maihda_long_data, id = "id", time = "wave", time_degree = 2)))
a_q <- suppressWarnings(suppressMessages(
  maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
         data = maihda_long_data, id = "id", time = "wave",
         time_degree = 2, decomposition = "longitudinal")))

# A deliberately data-starved fit (24 individuals spread across the 12 strata)
# so lme4 drives one or more random-effect variances to the boundary: the
# time-varying VPC must still be computable on a singular fit.
m_sing <- suppressWarnings(suppressMessages(
  fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
             data = subset(maihda_long_data,
                           id %in% unique(maihda_long_data$id)[1:24]),
             id = "id", time = "wave")))

test_that("quadratic growth fits an I(time^2) structure and a 3x3 covariance block", {
  expect_identical(m_q$longitudinal_info$time_degree, 2L)
  bars <- vapply(reformulas::findbars(m_q$formula),
                 function(b) paste(deparse(b), collapse = ""), character(1))
  expect_true(any(grepl("I\\(wave\\^2\\) \\| id", bars)))
  expect_true(any(grepl("I\\(wave\\^2\\) \\| stratum", bars)))

  s <- summary(m_q)
  # The stratum block is now (intercept, slope, slope^2): 3x3, ordered to a(t).
  expect_equal(dim(s$longitudinal$Sigma_stratum), c(3L, 3L))
  expect_true(all(s$longitudinal$vpc_t$estimate >= 0 &
                    s$longitudinal$vpc_t$estimate <= 1))
  # The components table carries the quadratic (slope^2) variance row plus every
  # off-diagonal covariance of the 3x3 block (not just intercept-slope).
  vc <- s$variance_components
  expect_identical(attr(vc, "kind"), "longitudinal")
  expect_true("Between-stratum: slope^2 (wave)" %in% vc$component)
  expect_true("Between-stratum: intercept-slope^2 covariance" %in% vc$component)
  expect_true("Between-stratum: slope-slope^2 covariance" %in% vc$component)
})

test_that("quadratic predict(type = 'strata') exposes the slope^2 term", {
  ps <- predict_maihda(m_q, type = "strata")
  expect_true(all(c("baseline", "intercept", "slope", "slope2") %in% names(ps)))
  expect_equal(nrow(ps), nrow(m_q$strata_info))
  # baseline = a(ref_time)' coef = intercept + slope*ref + slope2*ref^2.
  ref <- m_q$longitudinal_info$ref_time
  expect_equal(ps$baseline,
               ps$intercept + ps$slope * ref + ps$slope2 * ref^2,
               tolerance = 1e-8)
})

test_that("quadratic longitudinal PCV splits a 3x3 trajectory block", {
  expect_identical(a_q$mode, "longitudinal")
  expect_s3_class(a_q$pcv, "maihda_long_pcv")
  expect_identical(a_q$pcv$time_degree, 2L)
  expect_equal(dim(a_q$pcv$Sigma_stratum_null), c(3L, 3L))
  # Both the baseline and slope PCVs are genuine proportions inside (0, 1).
  expect_gt(a_q$pcv$pcv_intercept, 0); expect_lt(a_q$pcv$pcv_intercept, 1)
  expect_gt(a_q$pcv$pcv_slope, 0);     expect_lt(a_q$pcv$pcv_slope, 1)
  expect_s3_class(plot(a_q, type = "pcv_trajectory"), "ggplot")
})

test_that("the time-varying VPC stays finite on a singular fit", {
  expect_true(lme4::isSingular(m_sing$model))   # a boundary fit by construction
  s <- summary(m_sing)
  expect_false(is.null(s$longitudinal))
  expect_true(is.finite(s$vpc$estimate))
  expect_true(all(is.finite(s$longitudinal$vpc_t$estimate)))
  expect_true(all(s$longitudinal$vpc_t$estimate >= 0 &
                    s$longitudinal$vpc_t$estimate <= 1))
  # Per-stratum trajectory parameters are still recoverable on the boundary fit.
  ps <- predict_maihda(m_sing, type = "strata")
  expect_true(all(is.finite(ps$baseline)))
})

test_that("a VPC-trajectory bootstrap survives singular bootstrap refits", {
  # On a boundary fit many parametric-bootstrap refits are themselves singular;
  # the per-draw tryCatch must keep the reference-time CI finite regardless.
  set.seed(101)
  s <- summary(m_sing, bootstrap = TRUE, n_boot = 20)
  expect_true(is.finite(s$vpc$ci_lower) && is.finite(s$vpc$ci_upper))
  expect_lte(s$vpc$ci_lower, s$vpc$ci_upper)
})

test_that("longitudinal summary anchors the headline VPC at a nonzero ref_time", {
  # Waves shifted to 10..14: the headline VPC must be VPC(ref_time = 10), i.e.
  # the partition of a(10)' Sigma a(10) -- NOT the raw time-0 value.
  d <- maihda_long_data
  d$wave <- d$wave + 10
  m <- suppressWarnings(suppressMessages(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = d, id = "id", time = "wave")))
  s <- summary(m)
  expect_identical(s$longitudinal$ref_time, 10)
  vt <- s$longitudinal$vpc_t
  expect_equal(s$vpc$estimate, vt$estimate[vt$time == 10], tolerance = 1e-8)
  # equals the partition recomputed by hand at ref_time = 10
  Ss <- s$longitudinal$Sigma_stratum; Si <- s$longitudinal$Sigma_id
  vr <- s$longitudinal$var_resid
  vs <- maihda_var_at_time(Ss, 10); vi <- maihda_var_at_time(Si, 10)
  expect_equal(s$vpc$estimate, vs / (vs + vi + vr), tolerance = 1e-8)
})

test_that("maihda_longitudinal_pcv honours an explicit times grid", {
  # The time-specific PCV is reported on whatever times the caller supplies (not
  # only the default reporting grid); each row is (Var_null - Var_adj) / Var_null.
  times <- c(0, 2, 4)
  pcv <- maihda_longitudinal_pcv(a_g$model, a_g$model_adjusted, times = times)
  expect_equal(pcv$pcv_t$time, times)
  expect_true(all(c("var_null", "var_adjusted", "pcv") %in% names(pcv$pcv_t)))
  # Recompute the time-specific PCV straight from the covariance blocks.
  vn <- maihda_var_at_time(pcv$Sigma_stratum_null, times)
  va <- maihda_var_at_time(pcv$Sigma_stratum_adjusted, times)
  expect_equal(pcv$pcv_t$var_null, vn, tolerance = 1e-10)
  expect_equal(pcv$pcv_t$var_adjusted, va, tolerance = 1e-10)
  expect_equal(pcv$pcv_t$pcv, (vn - va) / vn, tolerance = 1e-10)
})

test_that("longitudinal plot builders carry the summary/PCV data through to the ggplot", {
  s <- summary(m_g)
  pv <- plot_vpc_trajectory(s)
  # The VPC-curve layer plots exactly the summary's vpc_t rows.
  expect_equal(pv$data$time, s$longitudinal$vpc_t$time)
  expect_equal(pv$data$estimate, s$longitudinal$vpc_t$estimate)

  pp <- plot_pcv_trajectory(a_g$pcv)
  expect_equal(pp$data$pcv, a_g$pcv$pcv_t$pcv)

  # One trajectory line per stratum over the reporting grid: grid x strata rows.
  pt <- plot_stratum_trajectories(m_g, s)
  expect_equal(nrow(pt$data),
               length(s$longitudinal$time_grid) * nrow(m_g$strata_info))
  expect_equal(length(unique(pt$data$stratum)), nrow(m_g$strata_info))
})

# ---- guards ----------------------------------------------------------------

test_that("scalar between-variance helpers reject a longitudinal model", {
  expect_error(extract_between_variance(m_g), "time-varying")
  m2 <- fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                   data = maihda_long_data, id = "id", time = "wave")
  expect_error(calculate_pvc(m_g, m2), "time-varying")
})

test_that("a non-longitudinal random slope is still rejected by summary", {
  strata <- make_strata(maihda_long_data,
                        vars = c("gender", "ethnicity", "education"))
  d <- strata$data
  # No id/time: this is NOT tagged longitudinal, so the intercept-only guard fires.
  m_bad <- fit_maihda(wellbeing ~ wave + (wave | stratum), data = d)
  expect_null(m_bad$longitudinal_info)
  expect_error(summary(m_bad), "intercept-only")
})

test_that("maihda() rejects incompatible longitudinal combinations", {
  expect_error(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "crossed-dimensions"),
    "longitudinal")
  expect_error(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave", group = "gender"),
    "does not support")
  # time without id
  expect_error(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, time = "wave",
           decomposition = "longitudinal"),
    "requires both")
  # per-group longitudinal comparison is out of scope
  expect_error(
    compare_maihda_groups(wellbeing ~ wave + (1 | gender:ethnicity:education),
                          data = maihda_long_data, group = "education",
                          decomposition = "longitudinal"),
    "out of scope")
})

# ---- dataset ---------------------------------------------------------------

test_that("maihda_long_data is long-format with repeated measures", {
  data(maihda_long_data, package = "MAIHDA")
  expect_true(all(c("id", "wave", "gender", "ethnicity", "education",
                    "wellbeing", "low_wellbeing") %in% names(maihda_long_data)))
  expect_gt(anyDuplicated(maihda_long_data$id), 0)
  expect_true(is.numeric(maihda_long_data$wave))
})
