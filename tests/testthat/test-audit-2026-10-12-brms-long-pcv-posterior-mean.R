# Audit 2026-10-12 (brms-long-pcv-posterior-mean): the longitudinal PCV of
# engine = "brms" fits compares the posterior MEDIANS of the per-draw between-stratum
# variances.
#
# PCV_intercept, PCV_slope and the time-specific PCV used to be read off each growth
# model's posterior-MEAN covariance block, so every cell was one minus a ratio of
# posterior means from two separately fitted models. The posterior of a between-stratum
# variance is right-skewed when strata are few, the adjusted model's more than the
# null's, so that ratio understates the PCV. Measured against a known truth of 0.83
# with 12 strata: 0.73 from posterior means, 0.80 from per-draw medians, 0.82 from
# lme4 (REML).
#
# The blocks below pin the summary on hand-built draws (no Stan): the three reported
# pieces, the variance at each time as the median of the per-draw variance (not a
# quadratic form of summarised cells), the stored blocks staying posterior means,
# centred time, an intercept-only stratum level and the print note; then that an lme4
# pair is read off its fitted blocks as before, and the same identities on a real brms
# pair when MAIHDA_TEST_BRMS=true.

# Evaluate `expr` and put the caller's RNG state back, so this file leaves the random
# stream the later test files read as it found it (see setup-seed.R). brms draws its
# sampler seed from R's RNG.
long_pcv_keep_rng <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  force(expr)
}

# Stratum SD and correlation draws of a growth block, built from fixed normal scores
# (no random numbers). An odd number, so each median is one of the draws. The SDs are
# log-normal in shape, the adjusted model's more spread out than the null's, so every
# variance has its mean well above its median. The slope SD is the intercept SD's
# shape in another order and the correlation swings in sign, so the three cells do not
# rise and fall together: the median of the variance at a later time is then not the
# quadratic form of the three cell medians.
long_pcv_draws <- function(which = c("null", "adjusted"), time_term = "wave") {
  which <- match.arg(which)
  n <- 201L
  i <- seq_len(n)
  z <- stats::qnorm((i - 0.5) / n)
  if (identical(which, "null")) {
    sd0 <- 0.60 * exp(0.45 * z)
    sd1 <- 0.18 * exp(0.45 * z[(i * 64L) %% n + 1L])
    rho <- 0.9 * cos(2.1 * i)
  } else {
    sd0 <- 0.25 * exp(0.75 * z)
    sd1 <- 0.08 * exp(0.75 * z[(i * 50L) %% n + 1L])
    rho <- 0.9 * sin(1.3 * i)
  }
  out <- data.frame(b_Intercept = 0, sd0 = sd0, sd1 = sd1, rho = rho)
  names(out)[2:4] <- c("sd_stratum__Intercept", paste0("sd_stratum__", time_term),
                       paste0("cor_stratum__Intercept__", time_term))
  out
}

# The per-draw pieces of a draws frame, written out as the package builds them.
long_pcv_pieces <- function(draws) {
  sd0 <- draws[[2]]; sd1 <- draws[[3]]; rho <- draws[[4]]
  list(v0 = sd0^2, v1 = sd1^2, cov = rho * sd0 * sd1)
}
long_pcv_var_at <- function(p, t) p$v0 + 2 * t * p$cov + t^2 * p$v1

# A maihda_model shell around a brmsfit stub that carries its own draws.
long_pcv_stub <- function(draws, times = 0:4, time_term = "wave", center = 0,
                          ref_time = min(times), stratum_slope = TRUE) {
  structure(
    list(model = structure(list(draws = draws), class = "brmsfit"),
         engine = "brms",
         longitudinal_info = list(id = "id", time = "wave", time_degree = 1L,
                                  time_term = time_term, time_center = center,
                                  ref_time = ref_time, stratum_slope = stratum_slope),
         data = data.frame(wave = rep(times, 2))),
    class = "maihda_model"
  )
}

# The draws are read off the stub; nothing here needs brms or Stan.
long_pcv_mock_draws <- function(env = parent.frame()) {
  local_mocked_bindings(
    maihda_posterior_draws_brms = function(model) model$draws,
    .env = env
  )
}

test_that("a brms growth pair's PCV compares per-draw posterior medians", {
  long_pcv_mock_draws()
  m0 <- long_pcv_stub(long_pcv_draws("null"))
  m1 <- long_pcv_stub(long_pcv_draws("adjusted"))
  pn <- long_pcv_pieces(m0$model$draws)
  pa <- long_pcv_pieces(m1$model$draws)

  res <- maihda_longitudinal_pcv(m0, m1)
  expect_s3_class(res, "maihda_long_pcv")
  expect_identical(res$var_baseline_null, stats::median(pn$v0))
  expect_identical(res$var_baseline_adjusted, stats::median(pa$v0))
  expect_identical(res$var_slope_null, stats::median(pn$v1))
  expect_identical(res$var_slope_adjusted, stats::median(pa$v1))
  expect_equal(res$pcv_intercept, 1 - stats::median(pa$v0) / stats::median(pn$v0),
               tolerance = 1e-12)
  expect_equal(res$pcv_slope, 1 - stats::median(pa$v1) / stats::median(pn$v1),
               tolerance = 1e-12)
  expect_identical(res$estimation_used, "posterior")
  expect_identical(res$variance_summary, "posterior median")
  expect_identical(res$Sigma_summary, "posterior mean")

  # Not one minus the ratio of posterior means, which is far lower on these draws.
  pcv_int_means <- 1 - mean(pa$v0) / mean(pn$v0)
  pcv_slope_means <- 1 - mean(pa$v1) / mean(pn$v1)
  expect_lt(pcv_int_means, res$pcv_intercept - 0.1)
  expect_lt(pcv_slope_means, res$pcv_slope - 0.1)

  # estimation = "ML" refits nothing on a brms pair.
  res_ml <- maihda_longitudinal_pcv(m0, m1, estimation = "ML")
  keep <- setdiff(names(res), "estimation")
  expect_identical(res_ml[keep], res[keep])
  expect_false(res_ml$ml_refit)
})

test_that("the variance at each time is the median of the per-draw variance", {
  long_pcv_mock_draws()
  m0 <- long_pcv_stub(long_pcv_draws("null"))
  m1 <- long_pcv_stub(long_pcv_draws("adjusted"))
  pn <- long_pcv_pieces(m0$model$draws)
  pa <- long_pcv_pieces(m1$model$draws)
  med_at <- function(p, t) vapply(t, function(ti) stats::median(long_pcv_var_at(p, ti)),
                                  numeric(1))

  res <- maihda_longitudinal_pcv(m0, m1)
  expect_identical(res$pcv_t$time, 0:4)
  expect_equal(res$pcv_t$var_null, med_at(pn, 0:4), tolerance = 1e-12)
  expect_equal(res$pcv_t$var_adjusted, med_at(pa, 0:4), tolerance = 1e-12)
  expect_equal(res$pcv_t$pcv, 1 - med_at(pa, 0:4) / med_at(pn, 0:4), tolerance = 1e-12)
  # The baseline row is PCV_intercept.
  expect_equal(res$pcv_t$pcv[1], res$pcv_intercept, tolerance = 1e-12)

  # Times off the observed grid.
  res2 <- maihda_longitudinal_pcv(m0, m1, times = c(0.5, 2, 7))
  expect_equal(res2$pcv_t$var_null, med_at(pn, c(0.5, 2, 7)), tolerance = 1e-12)
  expect_equal(res2$pcv_t$pcv, 1 - med_at(pa, c(0.5, 2, 7)) / med_at(pn, c(0.5, 2, 7)),
               tolerance = 1e-12)

  # Fixture premise: neither a quadratic form of the three cell medians nor one of the
  # posterior-mean block gives these variances away from the baseline.
  cell_med <- function(p, t) {
    stats::median(p$v0) + 2 * t * stats::median(p$cov) + t^2 * stats::median(p$v1)
  }
  cell_mean <- function(p, t) mean(p$v0) + 2 * t * mean(p$cov) + t^2 * mean(p$v1)
  later <- 2:5                                   # rows for t = 1..4
  expect_true(all(abs(res$pcv_t$var_adjusted[later] - cell_med(pa, 1:4)) >
                    0.05 * res$pcv_t$var_adjusted[later]))
  expect_true(all(abs(res$pcv_t$var_adjusted[later] - cell_mean(pa, 1:4)) >
                    0.05 * res$pcv_t$var_adjusted[later]))
  pcv_cell_med <- 1 - cell_med(pa, 1:4) / cell_med(pn, 1:4)
  expect_true(all(abs(res$pcv_t$pcv[later] - pcv_cell_med) > 0.01))
})

test_that("the stored stratum blocks stay posterior means and are labelled", {
  long_pcv_mock_draws()
  m0 <- long_pcv_stub(long_pcv_draws("null"))
  m1 <- long_pcv_stub(long_pcv_draws("adjusted"))
  pn <- long_pcv_pieces(m0$model$draws)
  pa <- long_pcv_pieces(m1$model$draws)
  block <- function(p) {
    matrix(c(mean(p$v0), mean(p$cov), mean(p$cov), mean(p$v1)), nrow = 2,
           dimnames = list(c("(Intercept)", "wave"), c("(Intercept)", "wave")))
  }

  res <- maihda_longitudinal_pcv(m0, m1)
  expect_identical(res$Sigma_stratum_null, block(pn))
  expect_identical(res$Sigma_stratum_adjusted, block(pa))
  # The blocks summary() tabulates for such a fit.
  expect_identical(res$Sigma_stratum_null, maihda_re_block(m0, "stratum"))
  expect_identical(res$Sigma_stratum_adjusted, maihda_re_block(m1, "stratum"))
  # So the reported variances are not read off them.
  expect_gt(res$Sigma_stratum_adjusted[1, 1], 2 * res$var_baseline_adjusted)
  expect_gt(res$Sigma_stratum_adjusted[2, 2], 2 * res$var_slope_adjusted)
})

test_that("centred time and a baseline away from the coefficient origin are honoured", {
  long_pcv_mock_draws()
  ct <- ".maihda_ctime"
  m0 <- long_pcv_stub(long_pcv_draws("null", ct), times = 10:14, time_term = ct,
                      center = 10)
  m1 <- long_pcv_stub(long_pcv_draws("adjusted", ct), times = 10:14, time_term = ct,
                      center = 10)
  pn <- long_pcv_pieces(m0$model$draws)
  pa <- long_pcv_pieces(m1$model$draws)
  med_at <- function(p, t) vapply(t, function(ti) stats::median(long_pcv_var_at(p, ti)),
                                  numeric(1))

  res <- maihda_longitudinal_pcv(m0, m1)
  expect_identical(res$pcv_t$time, 10:14)
  expect_identical(res$time_center, 10)
  # Original times 10..14 are 0..4 on the coefficient scale.
  expect_equal(res$pcv_t$var_null, med_at(pn, 0:4), tolerance = 1e-12)
  expect_equal(res$pcv_t$var_adjusted, med_at(pa, 0:4), tolerance = 1e-12)
  expect_identical(res$var_baseline_null, stats::median(pn$v0))
  expect_identical(rownames(res$Sigma_stratum_null), c("(Intercept)", ct))

  # A fitted sample whose first occasion is one unit after the centring origin: the
  # baseline variance is the per-draw variance one unit along, the slope variance is
  # the same at every time.
  m0b <- m0; m0b$longitudinal_info$ref_time <- 11
  m1b <- m1; m1b$longitudinal_info$ref_time <- 11
  resb <- maihda_longitudinal_pcv(m0b, m1b)
  expect_equal(resb$var_baseline_null, med_at(pn, 1), tolerance = 1e-12)
  expect_equal(resb$var_baseline_adjusted, med_at(pa, 1), tolerance = 1e-12)
  expect_equal(resb$pcv_intercept, 1 - med_at(pa, 1) / med_at(pn, 1), tolerance = 1e-12)
  expect_identical(resb$var_slope_null, res$var_slope_null)
  expect_identical(resb$pcv_slope, res$pcv_slope)
  expect_false(isTRUE(all.equal(resb$var_baseline_null, res$var_baseline_null)))
})

test_that("an intercept-only stratum level gives one time-constant median variance", {
  long_pcv_mock_draws()
  d0 <- long_pcv_draws("null")[, c("b_Intercept", "sd_stratum__Intercept")]
  d1 <- long_pcv_draws("adjusted")[, c("b_Intercept", "sd_stratum__Intercept")]
  m0 <- long_pcv_stub(d0, stratum_slope = FALSE)
  m1 <- long_pcv_stub(d1, stratum_slope = FALSE)
  v0n <- d0$sd_stratum__Intercept^2
  v0a <- d1$sd_stratum__Intercept^2

  res <- maihda_longitudinal_pcv(m0, m1)
  expect_identical(res$var_baseline_null, stats::median(v0n))
  expect_identical(res$var_baseline_adjusted, stats::median(v0a))
  expect_equal(res$pcv_intercept, 1 - stats::median(v0a) / stats::median(v0n),
               tolerance = 1e-12)
  expect_true(is.na(res$pcv_slope))
  expect_true(is.na(res$var_slope_null))
  expect_false(res$stratum_slope)
  expect_equal(res$pcv_t$var_null, rep(stats::median(v0n), 5), tolerance = 1e-12)
  expect_equal(res$pcv_t$pcv, rep(res$pcv_intercept, 5), tolerance = 1e-12)
  expect_identical(dim(res$Sigma_stratum_null), c(1L, 1L))
  expect_identical(res$Sigma_stratum_null[1, 1], mean(v0n))
  expect_identical(res$variance_summary, "posterior median")

  out <- capture.output(print(res))
  expect_true(any(grepl("random intercept only", out, fixed = TRUE)))
  expect_false(any(grepl("variance at baseline:", out, fixed = TRUE)))
  expect_true(any(grepl("posterior median of its per-draw value", out, fixed = TRUE)))

  # A pair whose ADJUSTED model alone lacks the stratum slope (maihda() never builds
  # one): whether a slope cell is reported follows the null block, as for an lme4
  # pair, and the adjusted model then has no slope variance left.
  mixed <- maihda_longitudinal_pcv(long_pcv_stub(long_pcv_draws("null")), m1)
  expect_identical(mixed$var_slope_adjusted, 0)
  expect_identical(mixed$pcv_slope, 1)
  expect_identical(mixed$var_slope_null,
                   stats::median(long_pcv_pieces(long_pcv_draws("null"))$v1))
})

test_that("non-finite draws are left out and a null variance at zero gives NA", {
  long_pcv_mock_draws()
  d0 <- long_pcv_draws("null")
  d1 <- long_pcv_draws("adjusted")
  pa <- long_pcv_pieces(d1)

  # One unusable draw in the null model: the median is taken over the other 200.
  d0_bad <- d0
  d0_bad$sd_stratum__Intercept[7] <- Inf
  res <- maihda_longitudinal_pcv(long_pcv_stub(d0_bad), long_pcv_stub(d1))
  pn_ok <- long_pcv_pieces(d0[-7, ])
  expect_identical(res$var_baseline_null, stats::median(pn_ok$v0))
  expect_true(is.finite(res$pcv_intercept))
  expect_true(all(is.finite(res$pcv_t$pcv)))

  # No between-stratum variance in the null model at all: nothing to take a share of.
  d0_zero <- d0
  d0_zero$sd_stratum__Intercept <- 0
  d0_zero$sd_stratum__wave <- 0
  res0 <- maihda_longitudinal_pcv(long_pcv_stub(d0_zero), long_pcv_stub(d1))
  expect_identical(res0$var_baseline_null, 0)
  expect_true(is.na(res0$pcv_intercept))
  expect_true(is.na(res0$pcv_slope))
  expect_true(all(is.na(res0$pcv_t$pcv)))
  expect_identical(res0$var_baseline_adjusted, stats::median(pa$v0))
})

test_that("print says which summary the variances are, only for an object that says so", {
  long_pcv_mock_draws()
  m0 <- long_pcv_stub(long_pcv_draws("null"))
  m1 <- long_pcv_stub(long_pcv_draws("adjusted"))
  pn <- long_pcv_pieces(m0$model$draws)
  pa <- long_pcv_pieces(m1$model$draws)
  res <- maihda_longitudinal_pcv(m0, m1)

  out <- capture.output(print(res))
  expect_true(any(grepl(sprintf("variance: %.4f (null) -> %.4f (adjusted)",
                                stats::median(pn$v0), stats::median(pa$v0)),
                        out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("variance at baseline: %.4f (null) -> %.4f (adjusted)",
                                stats::median(pn$v1), stats::median(pa$v1)),
                        out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("PCV_intercept: %.1f%%", 100 * res$pcv_intercept),
                        out, fixed = TRUE)))
  note <- grepl("posterior median of its per-draw value", out, fixed = TRUE)
  expect_identical(sum(note), 1L)
  expect_true(any(grepl("posterior means of the covariance blocks", out, fixed = TRUE)))
  expect_true(any(grepl("PCV runs low", out, fixed = TRUE)))

  # An object stored without the field holds posterior-mean variances: no such claim.
  old <- res
  old$variance_summary <- NULL
  old$Sigma_summary <- NULL
  out_old <- capture.output(print(old))
  expect_false(any(grepl("posterior median", out_old, fixed = TRUE)))
  expect_identical(out_old, out[!grepl("brms fits:|posterior means of the covariance|PCV runs low",
                                       out)])
})

test_that("an lme4 growth pair is read off its fitted blocks, as before", {
  skip_on_cran()
  a <- long_pcv_keep_rng(suppressWarnings(suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = maihda_long_data,
           id = "id", time = "wave", decomposition = "longitudinal"))))
  p <- a$pcv
  Sn <- maihda_re_block(a$model, "stratum")
  Sa <- maihda_re_block(a$model_adjusted, "stratum")
  ref_c <- p$ref_time - p$time_center

  expect_identical(p$Sigma_stratum_null, Sn)
  expect_identical(p$Sigma_stratum_adjusted, Sa)
  expect_identical(p$var_baseline_null, maihda_var_at_time(Sn, ref_c))
  expect_identical(p$var_baseline_adjusted, maihda_var_at_time(Sa, ref_c))
  expect_identical(p$var_slope_null, maihda_slope_var_at_time(Sn, ref_c))
  expect_identical(p$var_slope_adjusted, maihda_slope_var_at_time(Sa, ref_c))
  expect_identical(p$pcv_t$var_null, maihda_var_at_time(Sn, p$pcv_t$time - p$time_center))
  expect_identical(p$pcv_t$var_adjusted,
                   maihda_var_at_time(Sa, p$pcv_t$time - p$time_center))
  expect_identical(p$pcv_intercept,
                   (p$var_baseline_null - p$var_baseline_adjusted) / p$var_baseline_null)
  expect_identical(p$variance_summary, "estimate")
  expect_identical(p$Sigma_summary, "estimate")
  expect_identical(p$estimation_used, "fitted")
  expect_false(any(grepl("posterior", capture.output(print(p)), fixed = TRUE)))
})

test_that("a real brms growth pair reports per-draw medians (opt-in Stan)", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # Short chains on a subset: every check below is an identity between the package's
  # numbers and the fits' own draws, which holds whether or not the chains have
  # converged, so convergence is not asserted.
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:60])
  a <- long_pcv_keep_rng(suppressWarnings(suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = d,
           id = "id", time = "wave", engine = "brms", decomposition = "longitudinal",
           chains = 2, iter = 600, warmup = 300, refresh = 0))))
  p <- a$pcv
  expect_s3_class(p, "maihda_long_pcv")
  tt <- maihda_lng_time_term(a$model$longitudinal_info)
  dn <- maihda_re_cov_draws_brms(maihda_posterior_draws_brms(a$model$model), "stratum", tt)
  da <- maihda_re_cov_draws_brms(maihda_posterior_draws_brms(a$model_adjusted$model),
                                 "stratum", tt)
  t_c <- p$pcv_t$time - p$time_center
  med_at <- function(b, t) vapply(t, function(ti) {
    stats::median(b$v0 + 2 * ti * b$cov + ti^2 * b$v1)
  }, numeric(1))

  expect_equal(p$var_baseline_null, med_at(dn, p$ref_time - p$time_center),
               tolerance = 1e-12)
  expect_equal(p$var_baseline_adjusted, med_at(da, p$ref_time - p$time_center),
               tolerance = 1e-12)
  expect_equal(p$var_slope_null, stats::median(dn$v1), tolerance = 1e-12)
  expect_equal(p$var_slope_adjusted, stats::median(da$v1), tolerance = 1e-12)
  expect_equal(p$pcv_slope, 1 - stats::median(da$v1) / stats::median(dn$v1),
               tolerance = 1e-12)
  expect_equal(p$pcv_t$pcv, 1 - med_at(da, t_c) / med_at(dn, t_c), tolerance = 1e-12)
  expect_identical(p$variance_summary, "posterior median")

  # The stored blocks are the posterior means summary() tabulates, and summary()'s
  # own block is that same matrix.
  expect_equal(p$Sigma_stratum_null[1, 1], mean(dn$v0), tolerance = 1e-12)
  expect_equal(p$Sigma_stratum_adjusted[2, 2], mean(da$v1), tolerance = 1e-12)
  expect_identical(p$Sigma_stratum_null, a$summary$longitudinal$Sigma_stratum)
  expect_identical(p$Sigma_stratum_adjusted, a$summary_adjusted$longitudinal$Sigma_stratum)
  expect_true(any(grepl("posterior median of its per-draw value",
                        capture.output(print(p)), fixed = TRUE)))
})
