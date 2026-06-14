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

test_that("predict(type = 'strata') returns trajectory parameters", {
  ps <- predict_maihda(m_g, type = "strata")
  expect_true(all(c("stratum", "intercept", "slope") %in% names(ps)))
  expect_equal(nrow(ps), nrow(m_g$strata_info))
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
