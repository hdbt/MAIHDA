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
  a <- maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
              data = maihda_long_data, id = "id", time = "wave",
              decomposition = "longitudinal")
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
  a <- maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
              data = maihda_long_data, id = "id", time = "wave",
              decomposition = "longitudinal")
  expect_s3_class(plot(a, type = "pcv_trajectory"), "ggplot")
  expect_s3_class(plot(a, type = "vpc_trajectory"), "ggplot")
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
})

# ---- dataset ---------------------------------------------------------------

test_that("maihda_long_data is long-format with repeated measures", {
  data(maihda_long_data, package = "MAIHDA")
  expect_true(all(c("id", "wave", "gender", "ethnicity", "education",
                    "wellbeing", "low_wellbeing") %in% names(maihda_long_data)))
  expect_gt(anyDuplicated(maihda_long_data$id), 0)
  expect_true(is.numeric(maihda_long_data$wave))
})
