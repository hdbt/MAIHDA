# Audit 2026-09-04 (second pass of that date; the first, on maihda_describe()
# trial counts, is test-audit-2026-09-04.R).
#
# FINDING: the Bayesian trajectory VPCs (Bell et al. 2024, eq. 5) were
# built from the POSTERIOR-MEAN covariance blocks and returned as bare point
# estimates, while every neighbouring quantity on the same brms summary -- vpc_t,
# the headline VPC, the contextual partition -- is summarised per draw with a
# credible interval. A ratio of posterior means is a different estimator from the
# posterior median of the ratio: the between-stratum variance posterior is
# right-skewed whenever the strata are few (a MAIHDA has few by construction), so
# the plug-in runs high. Measured on maihda_long_data (150 ids, 12 strata, 2000
# draws, max Rhat 1.007): 0.6145 against a posterior median of 0.5765, a 6.6%
# overstatement, and the credible interval that was never reported at all spans
# [0.336, 0.823]. The gap widens as the strata get fewer.
#
# Fixed by summarising both trajectory VPCs per draw (posterior median + credible
# interval) on the brms path, and by adding bootstrap intervals -- point estimates
# unchanged -- on the lme4 path, which already pays for the refits.

# A deterministic stand-in for a posterior: lognormal quantiles on a regular grid,
# so there is no RNG and the skew is exactly the one the finding names. The
# between-stratum variance is diffuse and right-skewed; the between-individual
# variance is well identified and near-symmetric.
skewed_draws <- function(n = 400, meanlog = log(0.3), sdlog = 0.9) {
  stats::qlnorm(seq(0.5 / n, 1 - 0.5 / n, length.out = n),
                meanlog = meanlog, sdlog = sdlog)
}

test_that("trajectory VPC draws give the posterior MEDIAN of the ratio, not the ratio of means", {
  vs <- skewed_draws()
  vi <- stats::qlnorm(seq(0.5 / 400, 1 - 0.5 / 400, length.out = 400),
                      meanlog = log(0.25), sdlog = 0.12)
  sig_s <- list(v0 = vs, v1 = rep(0.04, 400), cov = rep(0, 400))
  sig_i <- list(v0 = vi, v1 = rep(0.02, 400), cov = rep(0, 400))

  out <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, ref_c = 0,
                                                  has_slope = TRUE)

  ratios <- vs / (vs + vi)
  expect_equal(out$vpc_intercept, stats::median(ratios))
  expect_equal(out$vpc_intercept_ci,
               stats::quantile(ratios, c(0.025, 0.975), names = FALSE))
  # The point estimate is the 50th percentile of the same posterior the interval
  # is read off, so it lies inside it by construction -- an invariant the plug-in
  # would not have, since a ratio of means need not sit within the quantiles of
  # the ratio.
  expect_gte(out$vpc_intercept, out$vpc_intercept_ci[1])
  expect_lte(out$vpc_intercept, out$vpc_intercept_ci[2])

  # THE TEETH. The old estimator -- what maihda_longitudinal_trajectory_vpc()
  # returns from the posterior-MEAN blocks -- must be materially different, or this
  # test would pass against the bug it exists to catch.
  plug_in <- mean(vs) / (mean(vs) + mean(vi))
  expect_gt(plug_in - out$vpc_intercept, 0.05)
  # And it is the same number the point-estimate sibling produces from the
  # posterior-mean blocks, i.e. exactly the pre-fix output.
  Sigma_s <- matrix(c(mean(vs), 0, 0, 0.04), 2, 2)
  Sigma_i <- matrix(c(mean(vi), 0, 0, 0.02), 2, 2)
  expect_equal(maihda_longitudinal_trajectory_vpc(Sigma_s, Sigma_i, 0)$vpc_intercept,
               plug_in)
})

test_that("the per-draw intercept VPC uses the full growth block at a non-zero reference time", {
  n <- 200
  v0s <- skewed_draws(n)
  sig_s <- list(v0 = v0s, v1 = seq(0.02, 0.06, length.out = n),
                cov = seq(-0.01, 0.03, length.out = n))
  sig_i <- list(v0 = rep(0.25, n), v1 = rep(0.02, n), cov = rep(0.005, n))
  t_c <- 2

  out <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, ref_c = t_c,
                                                  has_slope = TRUE)
  at <- function(b) b$v0 + 2 * t_c * b$cov + t_c^2 * b$v1
  expect_equal(out$vpc_intercept, stats::median(at(sig_s) / (at(sig_s) + at(sig_i))))
  # The covariance term is load-bearing away from the origin: dropping it changes
  # the answer, so this is not silently a v0-only calculation.
  no_cov <- sig_s$v0 + t_c^2 * sig_s$v1
  expect_false(isTRUE(all.equal(
    out$vpc_intercept,
    stats::median(no_cov / (no_cov + sig_i$v0 + t_c^2 * sig_i$v1)))))
})

test_that("a degenerate posterior reproduces the plug-in exactly", {
  # With no posterior uncertainty the two estimators must coincide -- the new code
  # computes the same estimand, not a different one.
  n <- 50
  sig_s <- list(v0 = rep(0.36, n), v1 = rep(0.04, n), cov = rep(0.02, n))
  sig_i <- list(v0 = rep(0.24, n), v1 = rep(0.01, n), cov = rep(0.00, n))
  Sigma_s <- matrix(c(0.36, 0.02, 0.02, 0.04), 2, 2)
  Sigma_i <- matrix(c(0.24, 0.00, 0.00, 0.01), 2, 2)
  for (t_c in c(0, 1.5)) {
    d <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, t_c, has_slope = TRUE)
    p <- maihda_longitudinal_trajectory_vpc(Sigma_s, Sigma_i, t_c)
    expect_equal(d$vpc_intercept, p$vpc_intercept)
    expect_equal(d$vpc_slope, p$vpc_slope)
    expect_equal(d$vpc_intercept_ci, rep(p$vpc_intercept, 2))
  }
})

test_that("no stratum slope gives NA, not the 0 the degenerate block would return", {
  n <- 100
  # maihda_re_cov_draws_brms() reports an intercept-only level as v1 = cov = 0
  # rather than as a shorter block, so the slope ratio would silently compute
  # 0 / (0 + v1_id) = 0 unless has_slope is honoured.
  sig_s <- list(v0 = skewed_draws(n), v1 = rep(0, n), cov = rep(0, n))
  sig_i <- list(v0 = rep(0.25, n), v1 = rep(0.02, n), cov = rep(0, n))

  out <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, 0, has_slope = FALSE)
  expect_true(is.na(out$vpc_slope))
  expect_true(all(is.na(out$vpc_slope_ci)))
  expect_length(out$vpc_slope_ci, 2L)
  # The intercept share is still estimated.
  expect_true(is.finite(out$vpc_intercept))

  # Guard the trap explicitly: with has_slope = TRUE the same degenerate block
  # WOULD return 0, which is why the flag cannot be inferred from the draws.
  wrong <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, 0, has_slope = TRUE)
  expect_equal(wrong$vpc_slope, 0)
})

test_that("draws with no total variance drop out instead of counting as zero", {
  n <- 100
  zero <- rep(0, 20)
  sig_s <- list(v0 = c(zero, rep(0.3, n - 20)), v1 = rep(0, n), cov = rep(0, n))
  sig_i <- list(v0 = c(zero, rep(0.3, n - 20)), v1 = rep(0, n), cov = rep(0, n))
  out <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, 0, has_slope = FALSE)
  # 0/0 is undefined, not 0: the 20 empty draws are dropped, so the median is the
  # 0.5 the remaining draws imply rather than being dragged toward 0.
  expect_equal(out$vpc_intercept, 0.5)

  # An entirely degenerate posterior has no defined share at all.
  allzero <- list(v0 = rep(0, n), v1 = rep(0, n), cov = rep(0, n))
  none <- maihda_longitudinal_trajectory_vpc_draws(allzero, allzero, 0,
                                                   has_slope = TRUE)
  expect_true(is.na(none$vpc_intercept))
  expect_true(all(is.na(none$vpc_intercept_ci)))
  expect_true(is.na(none$vpc_slope))
})

test_that("conf_level is honoured by the trajectory VPC credible interval", {
  vs <- skewed_draws(500)
  sig_s <- list(v0 = vs, v1 = rep(0.04, 500), cov = rep(0, 500))
  sig_i <- list(v0 = rep(0.25, 500), v1 = rep(0.02, 500), cov = rep(0, 500))
  wide <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, 0, TRUE,
                                                   conf_level = 0.95)
  narrow <- maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, 0, TRUE,
                                                     conf_level = 0.50)
  expect_gt(diff(wide$vpc_intercept_ci), diff(narrow$vpc_intercept_ci))
  expect_equal(narrow$vpc_intercept_ci,
               stats::quantile(vs / (vs + sig_i$v0), c(0.25, 0.75), names = FALSE))
  # Summarising thousands of draws must not warn (a quantile() on a vector
  # carrying dropped NA draws is the obvious way this would start to).
  expect_silent(maihda_longitudinal_trajectory_vpc_draws(sig_s, sig_i, 0, TRUE))
})

# ---------------------------------------------------------------------------
# lme4 wiring: intervals added, point estimates untouched.
# ---------------------------------------------------------------------------

test_that("lme4 trajectory VPCs gain a bootstrap interval without moving the estimate", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:200])
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = d, id = "id", time = "wave")))

  s0 <- suppressWarnings(summary(m))
  expect_true(all(is.na(s0$longitudinal$vpc_intercept_ci)))
  expect_true(all(is.na(s0$longitudinal$vpc_slope_ci)))
  expect_true(is.na(s0$longitudinal$trajectory_vpc_method))
  # Still the plug-in from the fitted blocks (the 2026-08-25 contract).
  Ss <- s0$longitudinal$Sigma_stratum; Si <- s0$longitudinal$Sigma_id
  expect_equal(s0$longitudinal$vpc_intercept, Ss[1, 1] / (Ss[1, 1] + Si[1, 1]))

  set.seed(4)
  s1 <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 60))
  # The interval is new; the point estimate must be BIT-IDENTICAL, since the
  # bootstrap only reuses refits the VPC(t) band was already paying for.
  expect_identical(s1$longitudinal$vpc_intercept, s0$longitudinal$vpc_intercept)
  expect_identical(s1$longitudinal$vpc_slope, s0$longitudinal$vpc_slope)
  expect_true(all(is.finite(s1$longitudinal$vpc_intercept_ci)))
  expect_true(all(is.finite(s1$longitudinal$vpc_slope_ci)))
  expect_lt(s1$longitudinal$vpc_intercept_ci[1], s1$longitudinal$vpc_intercept_ci[2])
  expect_identical(s1$longitudinal$trajectory_vpc_method, "bootstrap")

  out <- paste(utils::capture.output(print(s1)), collapse = "\n")
  expect_match(out, "Trajectory VPCs (Bell et al. 2024, eq. 5", fixed = TRUE)
  expect_match(out, sprintf("[%.4f, %.4f]", s1$longitudinal$vpc_intercept_ci[1],
                            s1$longitudinal$vpc_intercept_ci[2]), fixed = TRUE)
  expect_match(out, "(Bootstrap 95% CI)", fixed = TRUE)

  # Without a bootstrap the line stays bare -- no empty brackets, no label.
  out0 <- paste(utils::capture.output(print(s0)), collapse = "\n")
  expect_match(out0, "Trajectory VPCs (Bell et al. 2024, eq. 5", fixed = TRUE)
  expect_false(grepl("(Bootstrap 95% CI)", out0, fixed = TRUE))
  expect_false(grepl("credible interval", out0, fixed = TRUE))
})

test_that("a bootstrapped fit with no stratum slope keeps an NA slope interval", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:200])
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = d, id = "id", time = "wave", stratum_slope = FALSE)))
  set.seed(5)
  s <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 40))
  expect_true(is.na(s$longitudinal$vpc_slope))
  expect_true(all(is.na(s$longitudinal$vpc_slope_ci)))
  expect_true(all(is.finite(s$longitudinal$vpc_intercept_ci)))
  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  # The 2026-08-25 contract: "not estimated", never a bracketed zero.
  expect_match(out, "Slope: NA", fixed = TRUE)
  expect_false(grepl("Slope: NA [", out, fixed = TRUE))
})

test_that("conf_level reaches the trajectory interval and its printed label", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:200])
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = d, id = "id", time = "wave")))
  set.seed(6)
  s95 <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 60))
  set.seed(6)
  s80 <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 60,
                                  conf_level = 0.80))
  # Same seed, so the replicates are the same and only the quantiles move.
  expect_identical(s80$longitudinal$vpc_intercept, s95$longitudinal$vpc_intercept)
  expect_lt(diff(s80$longitudinal$vpc_intercept_ci),
            diff(s95$longitudinal$vpc_intercept_ci))
  expect_lt(diff(s80$longitudinal$vpc_slope_ci),
            diff(s95$longitudinal$vpc_slope_ci))
  out <- paste(utils::capture.output(print(s80)), collapse = "\n")
  expect_match(out, "(Bootstrap 80% CI)", fixed = TRUE)
})

test_that("the trajectory intervals survive into a maihda() longitudinal analysis", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:150])
  set.seed(7)
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = d,
           id = "id", time = "wave", decomposition = "longitudinal",
           bootstrap = TRUE, n_boot = 40, interactions = FALSE)))
  # summary.maihda_analysis() hands back the summary maihda() already computed,
  # so the new fields reach the headline entry point only if they were stored.
  s <- summary(a)
  expect_true(all(is.finite(s$longitudinal$vpc_intercept_ci)))
  expect_identical(s$longitudinal$trajectory_vpc_method, "bootstrap")
  expect_equal(s$longitudinal$vpc_intercept,
               summary(a$model)$longitudinal$vpc_intercept)
  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "(Bootstrap 95% CI)", fixed = TRUE)
})

test_that("print tolerates a summary that predates the interval fields", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:150])
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = d, id = "id", time = "wave")))
  s <- suppressWarnings(summary(m))

  traj_line <- function(obj) {
    o <- utils::capture.output(print(obj))
    o[grep("Trajectory VPCs", o, fixed = TRUE) + 1L]
  }

  # A summary object saved before this change has none of the three new fields.
  # It must still print, with the bare pre-change line and no interval label.
  old <- s
  old$longitudinal$vpc_intercept_ci <- NULL
  old$longitudinal$vpc_slope_ci <- NULL
  old$longitudinal$trajectory_vpc_method <- NULL
  expect_silent(out_old <- utils::capture.output(print(old)))
  expect_false(grepl("[", traj_line(old), fixed = TRUE))
  expect_false(any(grepl("credible interval", out_old, fixed = TRUE)))
  expect_false(any(grepl("Bootstrap 95% CI", out_old, fixed = TRUE)))

  # Malformed intervals must be ignored rather than half-printed: the guard tests
  # length AND finiteness, not merely non-NULL.
  bad_len <- s; bad_len$longitudinal$vpc_intercept_ci <- 0.5
  expect_false(grepl("[", traj_line(bad_len), fixed = TRUE))
  half_na <- s; half_na$longitudinal$vpc_intercept_ci <- c(0.3, NA_real_)
  expect_false(grepl("[", traj_line(half_na), fixed = TRUE))
})

# ---------------------------------------------------------------------------
# brms end-to-end (opt-in: needs a Stan toolchain).
# ---------------------------------------------------------------------------

test_that("brms trajectory VPCs are posterior medians of per-draw ratios with a CrI", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  data(maihda_long_data, envir = environment())

  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:150])
  m <- suppressWarnings(suppressMessages(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = d,
               id = "id", time = "wave", engine = "brms",
               chains = 2, iter = 2000, warmup = 1000, refresh = 0, seed = 7,
               control = list(adapt_delta = 0.95))))
  s <- summary(m)
  lg <- s$longitudinal
  expect_identical(lg$trajectory_vpc_method, "posterior")
  expect_true(all(is.finite(lg$vpc_intercept_ci)))
  expect_true(all(is.finite(lg$vpc_slope_ci)))

  # Recomputed from the raw posterior, independently of the summary path.
  lng <- m$longitudinal_info
  tt <- maihda_lng_time_term(lng)
  ref_c <- lng$ref_time - maihda_lng_time_center(lng)
  pd <- as.data.frame(brms::as_draws_df(m$model))
  blk <- function(g) {
    s0 <- pd[[paste0("sd_", g, "__Intercept")]]
    s1 <- pd[[paste0("sd_", g, "__", tt)]]
    r <- pd[[paste0("cor_", g, "__Intercept__", tt)]]
    list(v0 = s0^2, v1 = s1^2, cov = r * s0 * s1)
  }
  bs <- blk("stratum"); bi <- blk(lng$id)
  at <- function(b) b$v0 + 2 * ref_c * b$cov + ref_c^2 * b$v1
  r_int <- at(bs) / (at(bs) + at(bi))
  r_slope <- bs$v1 / (bs$v1 + bi$v1)
  expect_equal(lg$vpc_intercept, stats::median(r_int))
  expect_equal(lg$vpc_intercept_ci,
               stats::quantile(r_int, c(0.025, 0.975), names = FALSE))
  expect_equal(lg$vpc_slope, stats::median(r_slope))
  expect_equal(lg$vpc_slope_ci,
               stats::quantile(r_slope, c(0.025, 0.975), names = FALSE))
  # Both estimates sit inside their own credible intervals (see the unit test).
  expect_gte(lg$vpc_intercept, lg$vpc_intercept_ci[1])
  expect_lte(lg$vpc_intercept, lg$vpc_intercept_ci[2])
  expect_gte(lg$vpc_slope, lg$vpc_slope_ci[1])
  expect_lte(lg$vpc_slope, lg$vpc_slope_ci[2])

  # It is NOT the old ratio-of-posterior-means, which on this fit runs ~6% high.
  Ss <- maihda_re_block_brms(m$model, "stratum", tt, 1L)
  Si <- maihda_re_block_brms(m$model, lng$id, tt, 1L)
  plug_in <- maihda_longitudinal_trajectory_vpc(Ss, Si, ref_c)
  expect_gt(plug_in$vpc_intercept, lg$vpc_intercept)

  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "(Posterior 95% credible interval)", fixed = TRUE)
})
