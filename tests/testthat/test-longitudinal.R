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

test_that("maihda_slope_var_at_time evaluates b(t)' Sigma b(t), b = da/dt", {
  # Linear block: the instantaneous-slope variance is Sigma[2,2] at every t.
  S2 <- matrix(c(2, 0.1, 0.1, 0.5), nrow = 2)
  expect_equal(maihda_slope_var_at_time(S2, c(0, 3, 10)), rep(0.5, 3))
  # Quadratic block: b(t) = (0, 1, 2t), so the slope variance is time-varying.
  S3 <- matrix(c(2,    0.3,  0.05,
                 0.3,  0.5,  0.02,
                 0.05, 0.02, 0.1), nrow = 3, byrow = TRUE)
  expect_equal(maihda_slope_var_at_time(S3, 0), 0.5)
  # at t = 1: S22 + 2*(2*S23) + 4*S33 = 0.5 + 4*0.02 + 4*0.1
  expect_equal(maihda_slope_var_at_time(S3, 1), 0.5 + 4 * 0.02 + 4 * 0.1)
})

test_that("internal time centering: helpers and reserved-name guards", {
  # Centering offset: the minimum finite time; 0 when the axis starts at 0.
  expect_identical(maihda_longitudinal_center(c(0, 1, 2)), 0)
  expect_identical(maihda_longitudinal_center(c(10, 11, NA, 14)), 10)
  expect_identical(maihda_longitudinal_center(c(-2, 0, 2)), -2)
  expect_identical(maihda_longitudinal_center(numeric(0)), 0)

  # NULL-safe accessors (objects stored by pre-centering package versions).
  expect_identical(maihda_lng_time_term(list(time = "wave")), "wave")
  expect_identical(maihda_lng_time_term(list(time = "wave",
                                             time_term = ".maihda_ctime")),
                   ".maihda_ctime")
  expect_identical(maihda_lng_time_center(list(time = "wave")), 0)
  expect_identical(maihda_lng_time_center(list(time_center = 10)), 10)

  d <- data.frame(pid = rep(1:4, each = 2), wave = rep(10:11, 4), y = rnorm(8),
                  g = rep(c("a", "b"), 4), h = rep(c("x", "y"), each = 4),
                  .maihda_ctime = 1)
  # A fresh user call whose data carries the reserved column is rejected (the
  # centered fit would overwrite it) ...
  expect_error(maihda_validate_longitudinal("pid", "wave", 1, d,
                                            formula = y ~ (1 | g:h)),
               "reserved")
  # ... but a package-derived refit -- whose formula already references the
  # internal column -- passes (the column is the package's own).
  expect_silent(maihda_validate_longitudinal("pid", "wave", 1, d,
                                             formula = y ~ .maihda_ctime + (1 | g:h)))
  # id/time may not use the reserved name.
  expect_error(maihda_validate_longitudinal("pid", ".maihda_ctime", 1, d),
               "reserved")
})

test_that("maihda_longitudinal_formula replaces raw time terms under centering", {
  f <- maihda_longitudinal_formula(y ~ wave + x + (1 | stratum), id = "pid",
                                   time = ".maihda_ctime", time_degree = 1,
                                   orig_time = "wave")
  labs <- attr(stats::terms(reformulas::nobars(f)), "term.labels")
  expect_true(".maihda_ctime" %in% labs)
  expect_false("wave" %in% labs)   # replaced, not kept alongside (collinear)
  expect_true("x" %in% labs)
  bars <- vapply(reformulas::findbars(f),
                 function(b) paste(deparse(b), collapse = ""), character(1))
  expect_true(any(grepl("maihda_ctime \\| pid", bars)))
  expect_true(any(grepl("maihda_ctime \\| stratum", bars)))
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

test_that("maihda_time_terms quotes non-syntactic names in polynomial terms", {
  # A syntactic time name is unchanged (maihda_quote_name() is a no-op), matching
  # the (Intercept), time, I(time^2), ... rownames lme4 records in VarCorr.
  expect_identical(maihda_time_terms("wave", 1L), "wave")
  expect_identical(maihda_time_terms("wave", 2L), c("wave", "I(wave^2)"))

  # A non-syntactic time name must be back-quoted INSIDE I() too, otherwise the
  # higher-degree term is unparseable formula text like I(time point^2).
  tt <- maihda_time_terms("time point", 2L)
  expect_identical(tt, c("`time point`", "I(`time point`^2)"))
  # Each term parses on its own and the whole growth formula is well-formed.
  expect_silent(lapply(tt, function(s) stats::as.formula(paste("~", s))))
  f <- maihda_longitudinal_formula(y ~ (1 | stratum), id = "pid",
                                   time = "time point", time_degree = 2L)
  expect_s3_class(f, "formula")
  expect_true("I(`time point`^2)" %in%
                attr(stats::terms(reformulas::nobars(f)), "term.labels"))
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

  # The intercept variance is labelled with the coefficient origin (raw time 0
  # here), NOT called the baseline (ref_time).
  expect_true("Between-stratum: intercept (time = 0)" %in% tab$component)
  expect_false(any(grepl("baseline", tab$component)))

  # Under internal centering the label carries the centering offset instead.
  tab10 <- maihda_longitudinal_components_table(Ss, Si, var_resid = 0.7,
                                                time = "wave", id = "pid",
                                                center = 10)
  expect_true("Between-stratum: intercept (time = 10)" %in% tab10$component)

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

test_that("longitudinal ids reused across strata are rejected", {
  # Ids numbered within a site/group (person "1" in every site) used to be
  # silently merged into ONE level-2 unit by (time | id), pooling different
  # people's trajectories. An id spanning more than one stratum is the
  # observable symptom and is now an error.
  set.seed(77)
  person <- rbind(
    data.frame(pid = 1:10, gender = "F", edu = rep(c("lo", "hi"), 5)),
    data.frame(pid = 1:10, gender = "M", edu = rep(c("lo", "hi"), 5))
  )
  d <- person[rep(seq_len(nrow(person)), each = 3), ]
  d$wave <- rep(0:2, nrow(person))
  d$y <- rnorm(nrow(d))

  expect_error(
    fit_maihda(y ~ wave + (1 | gender:edu), data = d, id = "pid", time = "wave"),
    "more than one stratum"
  )

  # The direct guard: unique-per-person ids pass; NA rows are ignored.
  d_ok <- d
  d_ok$pid <- paste(d_ok$gender, d_ok$pid)   # globally unique
  d_ok$stratum <- interaction(d_ok$gender, d_ok$edu)
  expect_silent(maihda_check_longitudinal_ids(d_ok, "pid"))
  d_na <- d_ok
  d_na$pid[1] <- NA
  expect_silent(maihda_check_longitudinal_ids(d_na, "pid"))

  # ...and the reused ids fail with the offending values named.
  d_bad <- d
  d_bad$stratum <- interaction(d_bad$gender, d_bad$edu)
  expect_error(maihda_check_longitudinal_ids(d_bad, "pid"),
               "10 id value\\(s\\)")
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

# Shifted-axis twins of m_g / a_g (waves moved to 10..): exercise the internal
# time centering, which must reproduce the zero-anchored results exactly -- the
# centered design is numerically identical to the 0-anchored coding. Before
# centering, offset-axis fits could silently converge to a false optimum.
d_shift <- maihda_long_data
d_shift$wave <- d_shift$wave + 10
m_s <- suppressMessages(
  fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
             data = d_shift, id = "id", time = "wave"))
a_s <- suppressMessages(
  maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
         data = d_shift, id = "id", time = "wave",
         decomposition = "longitudinal"))

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

test_that("longitudinal PCV: default 'fitted' uses REML, estimation = 'ML' refits", {
  # The null and adjusted growth models differ in fixed effects (the dimensions'
  # main effects + dim:time), across which REML applies a model-specific correction.
  # estimation selects the basis (see calculate_pcv()): the default "fitted" keeps
  # each fit's own REML covariance block; "ML" refits both first. The stored fits
  # (and summary()'s time-varying VPC) always keep REML.
  expect_true(lme4::isREML(a_g$model$model))
  expect_true(lme4::isREML(a_g$model_adjusted$model))

  ref_c <- a_g$model$longitudinal_info$ref_time -
    maihda_lng_time_center(a_g$model$longitudinal_info)

  # Default (fitted): NO ML refit; the baseline variances equal the REML blocks.
  expect_false(isTRUE(a_g$pcv$ml_refit))
  vn_reml <- maihda_var_at_time(maihda_re_block(a_g$model, "stratum"), ref_c)
  va_reml <- maihda_var_at_time(maihda_re_block(a_g$model_adjusted, "stratum"), ref_c)
  expect_equal(a_g$pcv$var_baseline_null, vn_reml, tolerance = 1e-6)
  expect_equal(a_g$pcv$var_baseline_adjusted, va_reml, tolerance = 1e-6)
  expect_equal(a_g$pcv$pcv_intercept, (vn_reml - va_reml) / vn_reml, tolerance = 1e-6)

  # estimation = "ML": refit both growth fits with ML before differencing the blocks.
  a_mlbasis <- suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal", estimation = "ML"))
  expect_true(isTRUE(a_mlbasis$pcv$ml_refit))
  null_ml <- a_g$model;          null_ml$model <- lme4::refitML(a_g$model$model)
  adj_ml  <- a_g$model_adjusted; adj_ml$model  <- lme4::refitML(a_g$model_adjusted$model)
  vn <- maihda_var_at_time(maihda_re_block(null_ml, "stratum"), ref_c)
  va <- maihda_var_at_time(maihda_re_block(adj_ml, "stratum"), ref_c)
  expect_equal(a_mlbasis$pcv$var_baseline_null, vn, tolerance = 1e-6)
  expect_equal(a_mlbasis$pcv$var_baseline_adjusted, va, tolerance = 1e-6)
  expect_equal(a_mlbasis$pcv$pcv_intercept, (vn - va) / vn, tolerance = 1e-6)

  # The ML basis genuinely differs from the REML (fitted) baseline it replaces.
  expect_false(isTRUE(all.equal(a_mlbasis$pcv$var_baseline_null, vn_reml,
                                tolerance = 1e-6)))

  # An explicitly ML-FITTED pair (REML = FALSE) needs no refit under either basis and
  # reproduces the ML-basis decomposition.
  a_mlfit <- suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal", REML = FALSE))
  expect_false(isTRUE(a_mlfit$pcv$ml_refit))   # nothing needed refitting
  expect_equal(a_mlfit$pcv$pcv_intercept, a_mlbasis$pcv$pcv_intercept, tolerance = 1e-4)
  expect_equal(a_mlfit$pcv$pcv_slope, a_mlbasis$pcv$pcv_slope, tolerance = 1e-4)

  # The print method discloses the ML basis when (and only when) a refit ran.
  expect_output(print(a_mlbasis$pcv), "refitted with maximum likelihood")
  expect_false(any(grepl("refitted with maximum likelihood",
                         capture.output(print(a_g$pcv)))))
})

test_that("maihda_interactions refuses a longitudinal analysis or model", {
  # A direct call must NOT fall through to the scalar crossed-dimensions diagnostic:
  # each stratum's interaction is a TRAJECTORY (intercept + slope), so the single
  # per-stratum BLUP would drop the slope and return a cross-sectional value. Both
  # the longitudinal analysis (mode = "longitudinal") and a bare longitudinal model
  # (carrying longitudinal_info) error, mirroring the automatic-attachment skip.
  expect_error(maihda_interactions(a_g), "longitudinal MAIHDA")
  expect_error(maihda_interactions(m_g), "longitudinal MAIHDA")
})

test_that("a shifted time axis reproduces the zero-anchored fit (internal centering)", {
  # Waves moved to 10..: the growth terms are fit on internally centered time,
  # so the fit is the SAME optimization problem as the 0-anchored coding. Before
  # centering, this offset-axis fit silently converged far below the true
  # optimum (no lme4 warning), corrupting the VPC and the baseline variance.
  expect_identical(m_s$longitudinal_info$time_term, ".maihda_ctime")
  expect_equal(m_s$longitudinal_info$time_center, 10)
  expect_equal(m_s$longitudinal_info$ref_time, 10)
  expect_equal(as.numeric(logLik(m_s$model)), as.numeric(logLik(m_g$model)),
               tolerance = 1e-6)

  s <- summary(m_s); s0 <- summary(m_g)
  expect_equal(s$longitudinal$time_center, 10)
  expect_equal(s$vpc$estimate, s0$vpc$estimate, tolerance = 1e-6)
  # The reporting grid stays on the ORIGINAL axis; the VPC curve matches the
  # zero-anchored fit's point for point.
  expect_equal(s$longitudinal$vpc_t$time, s0$longitudinal$vpc_t$time + 10)
  expect_equal(s$longitudinal$vpc_t$estimate, s0$longitudinal$vpc_t$estimate,
               tolerance = 1e-6)
  # The headline VPC anchors at ref_time = 10 on the original axis; recomputed by
  # hand the covariance blocks are in CENTERED coordinates, so the baseline is
  # a(0)' Sigma a(0).
  vt <- s$longitudinal$vpc_t
  expect_equal(s$vpc$estimate, vt$estimate[vt$time == 10], tolerance = 1e-8)
  Ss <- s$longitudinal$Sigma_stratum; Si <- s$longitudinal$Sigma_id
  vr <- s$longitudinal$var_resid
  vs <- maihda_var_at_time(Ss, 0); vi <- maihda_var_at_time(Si, 0)
  expect_equal(s$vpc$estimate, vs / (vs + vi + vr), tolerance = 1e-8)
  # The components table names the coefficient origin (the centering offset).
  expect_true(any(grepl("intercept \\(time = 10\\)",
                        s$variance_components$component)))
  # The stored analytic data carries BOTH the original and the centered column.
  expect_true(all(c("wave", ".maihda_ctime") %in% names(m_s$data)))
  expect_equal(m_s$data$wave, m_s$data$.maihda_ctime + 10)
})

test_that("the longitudinal PCV is invariant to the time anchoring", {
  pcv <- a_s$pcv
  expect_equal(pcv$ref_time, 10)
  expect_equal(pcv$time_center, 10)

  # Matches the zero-anchored decomposition exactly: pcv_slope is the
  # instantaneous-slope variance PCV at the baseline (for linear growth the
  # slope variance is the same at every time), and pcv_intercept the baseline
  # variance PCV.
  expect_equal(pcv$pcv_intercept, a_g$pcv$pcv_intercept, tolerance = 1e-6)
  expect_equal(pcv$pcv_slope, a_g$pcv$pcv_slope, tolerance = 1e-6)
  expect_equal(pcv$var_baseline_null, a_g$pcv$var_baseline_null, tolerance = 1e-6)
  expect_equal(pcv$var_slope_null, a_g$pcv$var_slope_null, tolerance = 1e-6)
  expect_equal(unname(pcv$Sigma_stratum_null), unname(a_g$pcv$Sigma_stratum_null),
               tolerance = 1e-6)

  # Blocks are in centered coordinates: the baseline is the intercept cell.
  expect_equal(pcv$var_baseline_null,
               maihda_var_at_time(pcv$Sigma_stratum_null, 0))
  # The time-specific PCV grid stays on the original axis.
  expect_equal(pcv$pcv_t$time, a_g$pcv$pcv_t$time + 10)
  expect_equal(pcv$pcv_t$pcv, a_g$pcv$pcv_t$pcv, tolerance = 1e-6)

  # The print method reports the baseline on the original axis.
  expect_output(print(pcv), "Baseline \\(wave = 10\\)")
})

test_that("longitudinal models refuse cross-sectional scalar rankings/plots", {
  # A growth model's stratum estimand is a trajectory, so the scalar BLUP views
  # are not defined and must redirect to the trajectory tools.
  for (ty in c("predicted", "obs_vs_shrunken",
               "effect_decomp", "prediction_deviation")) {
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

test_that("predict(type = 'strata') trajectory parameters are anchoring-invariant", {
  # With waves shifted to 10.., the coefficients are in centered coordinates
  # anchored at the baseline: the baseline deviation IS the intercept (deviation
  # at the centering origin = ref_time), and both match the zero-anchored fit.
  ps <- predict_maihda(m_s, type = "strata")
  ps0 <- predict_maihda(m_g, type = "strata")
  expect_equal(m_s$longitudinal_info$ref_time, 10)
  expect_equal(ps$baseline, ps$intercept, tolerance = 1e-8)
  expect_equal(ps$baseline, ps0$baseline, tolerance = 1e-6)
  expect_equal(ps$slope, ps0$slope, tolerance = 1e-6)
})

test_that("individual predictions accept newdata on the original time axis", {
  # Caller newdata carries only the original time column; the centered column is
  # rebuilt internally, so predictions match those on the fitted rows -- and the
  # shifted fit predicts the same values as the zero-anchored fit.
  nd <- m_s$data[1:8, setdiff(names(m_s$data), ".maihda_ctime"), drop = FALSE]
  p_nd <- predict_maihda(m_s, newdata = nd, type = "individual")
  p_all <- predict_maihda(m_s, type = "individual")
  expect_equal(unname(p_nd), unname(p_all[1:8]), tolerance = 1e-8)
  p0 <- predict_maihda(m_g, type = "individual")
  expect_equal(unname(p_all), unname(p0), tolerance = 1e-6)
})

test_that("longitudinal ref_time/time_range come from the fitted frame, not dropped rows", {
  # lme4 drops rows with a missing outcome. If an entire baseline wave's outcomes
  # are missing, the reported VPC/PCV baseline must move to the first wave that
  # actually survives the fit -- staying at the pre-fit min(time) would anchor the
  # headline VPC/PCV at a time NOT represented in the fitted sample (an
  # extrapolation), because the time column itself is non-missing even where the
  # outcome is NA.
  d <- maihda_long_data
  waves <- sort(unique(d$wave))
  d$wellbeing[d$wave == waves[1]] <- NA   # wipe the entire baseline wave's outcome

  m <- suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = d, id = "id", time = "wave"))

  fitted_waves <- sort(unique(m$data$wave))
  expect_false(waves[1] %in% fitted_waves)              # baseline really dropped
  expect_gt(m$longitudinal_info$ref_time, waves[1])     # moved off the stale min
  expect_identical(m$longitudinal_info$ref_time, fitted_waves[1])
  expect_identical(m$longitudinal_info$time_range, range(fitted_waves))

  # The headline VPC summary anchors at the fitted baseline, not raw min(time).
  s <- summary(m)
  expect_identical(s$longitudinal$ref_time, fitted_waves[1])
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

test_that("trajectories select = 'deviation' keeps the most divergent strata", {
  s <- summary(m_g)
  all_strata <- unique(as.character(
    plot(m_g, type = "trajectories", n_strata = NULL)$data$stratum))
  cap <- 3L
  skip_if(length(all_strata) <= cap, "not enough strata to exercise the cap")

  p_dev <- plot(m_g, type = "trajectories", n_strata = cap, select = "deviation")
  p_ord <- plot(m_g, type = "trajectories", n_strata = cap, select = "order")
  dev <- unique(as.character(p_dev$data$stratum))
  expect_length(dev, cap)
  expect_length(unique(as.character(p_ord$data$stratum)), cap)
  expect_match(p_dev$labels$caption, "most extreme by trajectory deviation")

  # the survivors are exactly the most divergent by peak |random deviation| over
  # the grid: every kept stratum's peak is at least every dropped stratum's peak.
  re <- MAIHDA:::maihda_longitudinal_stratum_re(m_g)
  grid <- s$longitudinal$time_grid
  peak <- vapply(seq_len(nrow(re)), function(i) {
    a <- vapply(grid, function(t) sum(re$coef[[i]] * t^(0:(length(re$coef[[i]]) - 1))),
                numeric(1))
    max(abs(a))
  }, numeric(1))
  names(peak) <- as.character(re$stratum)
  expect_gte(min(peak[dev]), max(peak[setdiff(names(peak), dev)]))
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
  # A refit on a simulated draw occasionally trips lme4's convergence-warning
  # tolerance by a hair (max|grad| marginally over 0.002) -- expected Monte-Carlo
  # chatter in a 10-draw bootstrap, not the behaviour under test (the ribbon).
  s <- suppressWarnings(summary(m_g, bootstrap = TRUE, n_boot = 10))
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

test_that("count longitudinal fit evaluates the level-1 variance at time-varying marginal means", {
  # A Poisson growth model's marginal expected count carries a TIME-VARYING
  # lognormal correction: lambda(t) = exp(x'beta(t) + v(t)/2) with v(t) the stratum
  # plus individual growth-block variance a(t)' Sigma a(t) at time t. The level-1
  # variance log1p(1/lambda(t)) therefore changes over the reporting grid, so the
  # VPC(t) residual must be evaluated AT each grid time -- not one sample-wide
  # average reused everywhere. Recompute the grid residual independently, setting
  # the model's time term to each grid value (marginalizing over the rows) with the
  # correction v(t)/2 from the summary's blocks.
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:150])
  set.seed(909)
  d$events <- rpois(nrow(d), lambda = exp(0.6 + 0.15 * d$wave))
  mp <- suppressWarnings(suppressMessages(
    fit_maihda(events ~ wave + (1 | gender:ethnicity:education), data = d,
               id = "id", time = "wave", family = "poisson")))
  s <- suppressWarnings(summary(mp))
  expect_false(is.null(s$longitudinal))
  expect_true(all(is.finite(s$longitudinal$vpc_t$estimate)))

  lng <- s$longitudinal
  tt <- lng$time_term
  center <- lng$time_center
  frame <- MAIHDA:::maihda_model_frame(mp$model)
  resid_at <- function(t_orig) {
    tc <- t_orig - center
    nd <- frame
    nd[[tt]] <- tc
    eta <- as.numeric(stats::predict(mp$model, newdata = nd, re.form = NA,
                                     type = "link"))
    v_t <- maihda_var_at_time(lng$Sigma_stratum, tc) +
      maihda_var_at_time(lng$Sigma_id, tc)
    lambda <- pmax(exp(eta + v_t / 2), .Machine$double.eps)
    mean(log1p(1 / lambda))
  }

  # The per-grid residual matches the independent recomputation, and the headline
  # scalar is the reference-time value (NOT the old global average).
  expected_grid <- vapply(lng$time_grid, resid_at, numeric(1))
  expect_equal(lng$var_resid_t, expected_grid, tolerance = 1e-8)
  expect_equal(lng$var_resid, resid_at(lng$ref_time), tolerance = 1e-8)

  # The whole point of the fix: the residual is genuinely time-varying here (the
  # marginal count rises with wave, so log1p(1/lambda) falls), not a constant.
  expect_gt(diff(range(lng$var_resid_t)), 1e-6)
  # The old behavior -- one sample-wide average over rows at their OWN times --
  # differs from the corrected reference-time residual, so this is a real change.
  tv <- mp$data[[lng$time]] - center
  v_rows <- maihda_var_at_time(lng$Sigma_stratum, tv) +
    maihda_var_at_time(lng$Sigma_id, tv)
  eta_own <- as.numeric(stats::predict(mp$model, re.form = NA, type = "link"))
  old_global <- mean(log1p(1 / pmax(exp(eta_own + v_rows / 2), .Machine$double.eps)))
  expect_false(isTRUE(all.equal(lng$var_resid, old_global)))
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
  s <- suppressWarnings(summary(m_sing, bootstrap = TRUE, n_boot = 20))
  expect_true(is.finite(s$vpc$ci_lower) && is.finite(s$vpc$ci_upper))
  expect_lte(s$vpc$ci_lower, s$vpc$ci_upper)
})

test_that("quadratic growth is anchoring-invariant (regression: silent false convergence)", {
  # THE original defect: before internal centering, a quadratic growth fit on an
  # offset time axis (waves 10..) silently converged ~128 log-likelihood units
  # below the true optimum WITHOUT any lme4 convergence warning, reporting a
  # baseline between-stratum variance orders of magnitude off; and pcv_slope
  # read the raw Sigma[2,2] cell -- the slope variance at raw time 0, an
  # extrapolation that is NOT invariant to the time coding for degree >= 2.
  a_qs <- suppressWarnings(suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = d_shift, id = "id", time = "wave",
           time_degree = 2, decomposition = "longitudinal")))
  expect_equal(as.numeric(logLik(a_qs$model$model)),
               as.numeric(logLik(a_q$model$model)), tolerance = 1e-6)
  expect_equal(a_qs$pcv$pcv_intercept, a_q$pcv$pcv_intercept, tolerance = 1e-5)
  expect_equal(a_qs$pcv$pcv_slope, a_q$pcv$pcv_slope, tolerance = 1e-5)
  expect_equal(a_qs$pcv$var_baseline_null, a_q$pcv$var_baseline_null,
               tolerance = 1e-5)
  expect_equal(a_qs$pcv$var_slope_null, a_q$pcv$var_slope_null, tolerance = 1e-5)
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
  expect_error(calculate_pcv(m_g, m2), "time-varying")
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
