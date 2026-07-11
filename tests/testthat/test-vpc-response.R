# A binomial MAIHDA fit on synthetic data with a controllable between-stratum signal.
maihda_vpcr_fit <- function(seed = 7, n = 1500, sd_u = 1.2) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  lp <- stats::rnorm(nlevels(sk), sd = sd_u)[sk]
  d$y <- stats::rbinom(n, 1, stats::plogis(lp))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "binomial")
  ))
}

test_that("maihda_gauss_hermite_normal integrates standard-normal moments exactly", {
  gh <- MAIHDA:::maihda_gauss_hermite_normal(40)
  x <- gh$nodes
  w <- gh$weights

  expect_length(x, 40L)
  expect_equal(sum(w), 1, tolerance = 1e-12)        # weights are a probability measure
  expect_equal(sum(w * x), 0, tolerance = 1e-10)    # E[Z]   = 0
  expect_equal(sum(w * x^2), 1, tolerance = 1e-10)  # E[Z^2] = 1
  expect_equal(sum(w * x^4), 3, tolerance = 1e-9)   # E[Z^4] = 3
  expect_equal(sum(w * x^6), 15, tolerance = 1e-8)  # E[Z^6] = 15

  # A smooth bounded integrand (an inverse-logit, exactly the response-scale VPC
  # kernel) matches a fine Monte-Carlo reference.
  set.seed(1)
  zz <- stats::rnorm(2e6)
  expect_equal(sum(w * stats::plogis(0.3 + 1.2 * x)),
               mean(stats::plogis(0.3 + 1.2 * zz)), tolerance = 1e-3)

  expect_error(MAIHDA:::maihda_gauss_hermite_normal(1), "nodes")
})

test_that("maihda_vpc_response returns a VPC in [0, 1] and is seed-reproducible", {
  m <- maihda_vpcr_fit()
  r1 <- maihda_vpc_response(m, n_sim = 20000, seed = 42)
  r2 <- maihda_vpc_response(m, n_sim = 20000, seed = 42)

  expect_s3_class(r1, "maihda_vpc_response")
  expect_identical(r1$scale, "response")
  expect_true(r1$estimate >= 0 && r1$estimate <= 1)
  expect_equal(r1$estimate, r2$estimate)  # same seed -> identical draws -> identical VPC
})

test_that("maihda_vpc_response matches an independent Goldstein/Browne/Rasbash simulation", {
  m <- maihda_vpcr_fit()
  r <- maihda_vpc_response(m, n_sim = 50000, seed = 1)

  # Reimplement the computation independently with the same seed/draws.
  s2 <- MAIHDA:::extract_between_variance(m)
  lp <- mean(stats::predict(m$model, re.form = NA, type = "link"))
  set.seed(1)
  u <- stats::rnorm(50000, 0, sqrt(s2))
  p <- stats::plogis(lp + u)
  ref <- stats::var(p) / (stats::var(p) + mean(p * (1 - p)))

  expect_equal(r$estimate, ref, tolerance = 1e-6)
})

test_that("response-scale VPC increases with between-stratum variance", {
  low  <- maihda_vpcr_fit(seed = 3, sd_u = 0.3)
  high <- maihda_vpcr_fit(seed = 3, sd_u = 1.6)
  v_low  <- maihda_vpc_response(low,  n_sim = 50000, seed = 5)$estimate
  v_high <- maihda_vpc_response(high, n_sim = 50000, seed = 5)$estimate
  expect_gt(v_high, v_low)
})

test_that("maihda_vpc_response rejects non-binomial models and invalid n_sim", {
  g <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(120), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    family = "gaussian")))
  expect_error(maihda_vpc_response(g$model), "binomial")

  m <- maihda_vpcr_fit()
  expect_error(maihda_vpc_response(m, n_sim = 10), "n_sim")
  # A fractional count is rejected (rnorm() would silently truncate it while the
  # recorded n_sim kept the fraction).
  expect_error(maihda_vpc_response(m, n_sim = 100.5), "whole number")

  # A valid whole number is stored as an integer, matching the realised draw count.
  r <- maihda_vpc_response(m, n_sim = 100, seed = 1)
  expect_identical(r$n_sim, 100L)
})

test_that("maihda_vpc_response rejects random-slope models with a clear error", {
  # A manually-specified (x | stratum) binomial model carries a slope variance and
  # an intercept-slope covariance that the scalar simulation does not integrate over,
  # so its between-stratum variance is not a single number. It must error explicitly
  # (matching the scalar summary/VPC path) rather than reading only the intercept
  # variance and silently returning a wrong -- or NA -- response-scale VPC.
  skip_on_cran()
  set.seed(313)
  n <- 900
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE),
    x      = stats::rnorm(n)
  )
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  d$y <- stats::rbinom(n, 1, stats::plogis(0.3 * d$x +
           stats::rnorm(nlevels(factor(d$stratum)), sd = 0.6)[factor(d$stratum)]))
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (x | stratum), data = d, family = "binomial")))

  expect_error(maihda_vpc_response(m), "intercept-only random effects")

  # A slope on a NON-stratum grouping is rejected the same way (the simulation
  # would silently drop the site slope variance from var_other).
  d2 <- expand.grid(stratum = factor(1:6), site = factor(1:8), rep = 1:20)
  d2$x <- stats::rnorm(nrow(d2))
  d2$y <- stats::rbinom(nrow(d2), 1, stats::plogis(-0.2 + 0.4 * d2$x))
  m2 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum) + (x | site), data = d2, family = "binomial")))
  expect_error(maihda_vpc_response(m2), "intercept-only random effects")
})

test_that("response-scale VPC integrates over non-stratum random effects", {
  # Regression: for a site + stratum model the simulation drew ONLY the stratum
  # effect, so the site variance appeared in neither the numerator's scope nor
  # the denominator -- overstating the stratum share (audit repro: 0.0490
  # reported vs 0.0126 from a coherent nested Monte Carlo).
  skip_on_cran()
  set.seed(707)
  d <- expand.grid(stratum = factor(1:8), site = factor(1:10), rep = 1:25)
  su <- stats::rnorm(8, sd = 0.4)[d$stratum]
  so <- stats::rnorm(10, sd = 1.2)[d$site]
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.5 + su + so))
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | site) + (1 | stratum), data = d, family = "binomial")))

  v <- maihda_vpc_response(m, n_sim = 20000, seed = 42)
  expect_gt(v$var_other, 0)
  expect_true(v$estimate >= 0 && v$estimate <= 1)

  # Independent brute-force nested Monte Carlo of the stratum share.
  vc <- lme4::VarCorr(m$model)
  vs <- as.numeric(as.matrix(vc$stratum)[1, 1])
  vo <- as.numeric(as.matrix(vc$site)[1, 1])
  eta0 <- mean(stats::predict(m$model, re.form = NA, type = "link"))
  set.seed(99)
  us <- stats::rnorm(4000, 0, sqrt(vs))
  uo <- stats::rnorm(2000, 0, sqrt(vo))
  P <- stats::plogis(outer(us, uo, `+`) + eta0)
  truth <- stats::var(rowMeans(P)) /
    (stats::var(as.vector(P)) + mean(P * (1 - P)))
  expect_equal(v$estimate, truth, tolerance = 0.1)

  # The old stratum-only computation is strictly larger (the defect).
  p_only <- stats::plogis(eta0 + us)
  old <- stats::var(p_only) / (stats::var(p_only) + mean(p_only * (1 - p_only)))
  expect_gt(old, v$estimate)

  # print() explains the integration.
  expect_output(print(v), "non-stratum random effects")

  # A single-stratum fit is unchanged: var_other 0, same seeded value as the
  # stratum-only formula.
  m0 <- maihda_vpcr_fit()
  v0 <- maihda_vpc_response(m0, n_sim = 20000, seed = 42)
  expect_identical(v0$var_other, 0)
  expect_identical(v0$inner_method, "none")
  expect_true(is.na(v0$inner_nodes))
  s2 <- MAIHDA:::extract_between_variance(m0)
  lp <- mean(stats::predict(m0$model, re.form = NA, type = "link"))
  set.seed(42)
  u <- stats::rnorm(20000, 0, sqrt(s2))
  p <- stats::plogis(lp + u)
  expect_equal(v0$estimate,
               stats::var(p) / (stats::var(p) + mean(p * (1 - p))),
               tolerance = 1e-12)
})

test_that("response-scale VPC integrates the non-stratum effect by exact quadrature", {
  # The inner integral over the non-stratum effect is now DETERMINISTIC
  # Gauss-Hermite quadrature (no hidden fixed 500-draw inner Monte-Carlo sample),
  # so the estimate reproduces a reference that reuses the SAME outer draws over
  # the SAME quadrature nodes to floating-point tolerance -- the property the old
  # inner sample could not offer (its error did not shrink as n_sim grew).
  skip_on_cran()
  set.seed(707)
  d <- expand.grid(stratum = factor(1:8), site = factor(1:10), rep = 1:25)
  su <- stats::rnorm(8, sd = 0.4)[d$stratum]
  so <- stats::rnorm(10, sd = 1.2)[d$site]
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.5 + su + so))
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | site) + (1 | stratum), data = d, family = "binomial")))

  v <- maihda_vpc_response(m, n_sim = 20000, seed = 42)
  expect_gt(v$var_other, 0)
  expect_identical(v$inner_method, "gauss-hermite")
  expect_gt(v$inner_nodes, 1L)
  expect_output(print(v), "quadrature")

  # Reconstruct the exact computation: same seed/outer draws, same GH nodes.
  s2b <- MAIHDA:::extract_between_variance(m)
  vc <- lme4::VarCorr(m$model)
  vo <- as.numeric(as.matrix(vc$site)[1, 1])
  eta0 <- mean(stats::predict(m$model, re.form = NA, type = "link"))
  set.seed(42)
  us <- stats::rnorm(20000, 0, sqrt(s2b))
  gh <- MAIHDA:::maihda_gauss_hermite_normal(80)
  z <- gh$nodes * sqrt(vo)
  wq <- gh$weights
  P <- stats::plogis(outer(us, z, `+`) + eta0)
  ms <- as.vector(P %*% wq)
  ref <- stats::var(ms) /
    (max(mean(as.vector((P * P) %*% wq)) - mean(ms)^2, 0) +
       mean(as.vector((P * (1 - P)) %*% wq)))
  expect_equal(v$estimate, ref, tolerance = 1e-9)
})

test_that("maihda_vpc_response sums dimension + interaction variance for a crossed-dimensions fit", {
  # Audit finding: the response-scale VPC previously read only the interaction
  # ("stratum") component as between-stratum and dumped the additive dimension REs
  # into var_other, understating the intersectional VPC by orders of magnitude
  # (here ~0.005 vs ~0.22). For a crossed-dimensions fit the between-stratum
  # variance is the SUM of the dimension REs + the interaction -- independent
  # crossed effects that sum at the intersection level, exactly as the MOR
  # (maihda_mor_between_variance) and the crossed-dimensions VPC (maihda_cc_partition).
  skip_on_cran()
  set.seed(4242)
  n <- 6000
  d <- data.frame(
    a = sample(paste0("a", 1:4), n, replace = TRUE),
    b = sample(paste0("b", 1:4), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  ua <- stats::setNames(stats::rnorm(4, sd = 0.9), paste0("a", 1:4))
  ub <- stats::setNames(stats::rnorm(4, sd = 1.1), paste0("b", 1:4))
  stratum <- interaction(d$a, d$b, drop = TRUE)
  uint <- stats::rnorm(nlevels(stratum), sd = 0.25)[stratum]
  d$y <- stats::rbinom(n, 1, stats::plogis(-0.3 + ua[d$a] + ub[d$b] + uint))

  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ (1 | a:b), data = d, decomposition = "crossed-dimensions",
           family = "binomial")
  ))
  m <- cc$model
  expect_false(is.null(m$cc_info))

  av <- MAIHDA:::maihda_random_variances_lme4(m$model)
  groups <- unique(c(m$cc_info$interaction_group, unname(m$cc_info$dim_groups)))
  expected_between <- sum(av[groups])

  r <- maihda_vpc_response(m, n_sim = 20000, seed = 1)
  # var_between is the full intersectional variance (dims + interaction), NOT just
  # the interaction component; and nothing (no context here) leaks into var_other.
  expect_equal(unname(r$var_between), unname(expected_between), tolerance = 1e-8)
  expect_gt(r$var_between, 2 * av[["stratum"]])   # additive dims dominate
  expect_equal(r$var_other, 0, tolerance = 1e-10)

  # The VPC is a real intersectional share, far above the interaction-only value
  # the old partition produced (< 0.02 on this fit).
  expect_gt(r$estimate, 0.1)
  expect_lt(r$estimate, 1)
})
