# Regression tests for the 2026-07-13 statistical audit (six substantive defects
# plus three lower-severity issues). Each test pins the corrected behaviour so the
# specific defect cannot silently return.

test_that("predict_maihda keeps an external offset for training predictions (finding 1)", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  d <- data.frame(
    a = factor(sample(c("p", "q"), n, TRUE)),
    b = factor(sample(c("x", "y"), n, TRUE)),
    age = rnorm(n),
    exposure = runif(n, 1, 40)
  )
  d$y <- rpois(n, exp(-1 + 0.3 * d$age + log(d$exposure)))

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | a:b), data = d, family = poisson,
               offset = log(exposure))))
  # Training predictions (no newdata) must include the external offset, matching
  # the engine's own fitted values -- the pre-fix wrapper dropped it.
  expect_equal(as.numeric(predict_maihda(m, type = "response")),
               as.numeric(stats::fitted(m$model)),
               tolerance = 1e-8)
  # An external offset cannot be reconstructed for genuine newdata: reject, don't
  # silently return an offset-less prediction.
  expect_error(predict_maihda(m, newdata = d, type = "response"),
               "external offset")
  # A FORMULA offset is re-evaluated by predict.merMod and keeps working on newdata.
  mf <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + offset(log(exposure)) + (1 | a:b), data = d,
               family = poisson)))
  expect_equal(as.numeric(predict_maihda(mf, newdata = d, type = "response")),
               as.numeric(stats::fitted(mf$model)),
               tolerance = 1e-6)
})

test_that("crossed-dimensions rejects an undeclared extra random effect (finding 2)", {
  skip_on_cran()
  set.seed(3)
  n <- 900
  d <- data.frame(
    a = factor(sample(c("p", "q"), n, TRUE)),
    b = factor(sample(c("x", "y", "z"), n, TRUE)),
    site = factor(sample(paste0("s", 1:6), n, TRUE))
  )
  d$y <- rnorm(n)
  ds <- make_strata(d, vars = c("a", "b"))$data

  # A (1 | site) written in the formula would be silently dropped by the crossed
  # builder; it must error instead.
  expect_error(
    suppressWarnings(suppressMessages(
      maihda(y ~ (1 | stratum) + (1 | site), data = ds,
             decomposition = "crossed-dimensions"))),
    "crossed-dimensions")
  # Supplying the same grouping through context = composes correctly.
  mcc <- suppressWarnings(suppressMessages(
    maihda(y ~ (1 | stratum), data = ds, context = "site",
           decomposition = "crossed-dimensions")))
  expect_true("site" %in% names(lme4::VarCorr(mcc$model$model)))
})

test_that("longitudinal centering uses the analytic (subset) sample (finding 3)", {
  skip_on_cran()
  set.seed(11)
  np <- 250
  waves <- c(0, 50, 100, 150, 200, 250)
  d <- expand.grid(pid = 1:np, t = waves)
  d$a <- factor(d$pid %% 2)
  d$b <- factor(d$pid %% 3)
  u0 <- rnorm(np, 0, 0.6)[d$pid]
  u1 <- rnorm(np, 0, 0.003)[d$pid]
  d$y <- 2 + 0.01 * d$t + u0 + u1 * d$t + rnorm(nrow(d), 0, 0.8)
  d$pid <- factor(d$pid)

  # subset keeps waves 100..250 (4 analytic waves), so the analytic sample begins
  # at 100 while the full data begins at 0.
  m_sub <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | a:b), data = d, id = "pid", time = "t", subset = t >= 100)))
  m_pre <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | a:b), data = d[d$t >= 100, ], id = "pid", time = "t")))

  # The centre is the analytic minimum (100), not the full-data minimum (0), so a
  # subset behaves like fitting the pre-filtered data.
  expect_equal(m_sub$longitudinal_info$time_center, 100)
  expect_equal(m_sub$longitudinal_info$time_center,
               m_pre$longitudinal_info$time_center)
})

test_that("weighted PCV setup drops invalid-weight rows before family detection (finding 4)", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  set.seed(21)
  n <- 400
  d <- data.frame(
    a = factor(sample(c("p", "q"), n, TRUE)),
    b = factor(sample(c("x", "y"), n, TRUE)),
    w = runif(n, 0.5, 2)
  )
  d$y <- rbinom(n, 1, 0.4)
  # One contaminating row: outcome 2 with a zero sampling weight (wemix drops it).
  d <- rbind(d, data.frame(a = "p", b = "x", w = 0, y = 2))
  ds <- make_strata(d, vars = c("a", "b"))$data

  setup <- suppressWarnings(suppressMessages(
    maihda_pcv_attribution_setup(ds, "y", c("a", "b"), engine = "lme4",
      family = "gaussian", context = NULL, sampling_weights = "w",
      engine_missing = TRUE, family_missing = TRUE)))

  expect_identical(setup$family, "binomial")   # the stray 2 no longer forces Gaussian
  expect_equal(nrow(setup$data), n)            # n_obs not inflated by the invalid row
  expect_false(2 %in% setup$data$y)
})

test_that("autobin cut-points use the subset sample (finding 5)", {
  skip_on_cran()
  set.seed(5)
  n <- 600
  d <- data.frame(
    income = c(rnorm(n * 2 / 3, 20, 3), rnorm(n / 3, 45, 8)),  # bimodal
    sex = factor(sample(c("F", "M"), n, TRUE))
  )
  d$y <- rnorm(n)
  keep <- d$income < 30   # removes the high-income cluster

  b_sub <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | income:sex), data = d,
               subset = keep)))$strata_autobin_info$income$breaks
  b_pre <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | income:sex),
               data = d[keep, ])))$strata_autobin_info$income$breaks

  expect_equal(b_sub, b_pre, tolerance = 1e-9)
})

test_that("calculate_pcv/compare_maihda flag a changed random-effects structure (finding 6)", {
  skip_on_cran()
  set.seed(9)
  n <- 900
  d <- data.frame(
    a = factor(sample(c("p", "q"), n, TRUE)),
    b = factor(sample(c("x", "y", "z"), n, TRUE)),
    site = factor(sample(paste0("s", 1:6), n, TRUE)),
    x = rnorm(n)
  )
  strat <- interaction(d$a, d$b, drop = TRUE)
  strat_eff <- rnorm(nlevels(strat), 0, 0.7)[strat]   # real between-stratum signal
  site_eff <- rnorm(6, 0, 0.5)[d$site]
  d$y <- 0.5 * d$x + strat_eff + site_eff + rnorm(n, 0, 0.8)
  ds <- make_strata(d, vars = c("a", "b"))$data

  m_null <- suppressWarnings(suppressMessages(fit_maihda(y ~ (1 | stratum), data = ds)))
  m_adj <- suppressWarnings(suppressMessages(fit_maihda(y ~ x + (1 | stratum), data = ds)))
  m_site <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum) + (1 | site), data = ds)))

  # {stratum} vs {site, stratum}: a changed decomposition, not covariate adjustment.
  expect_error(calculate_pcv(m_null, m_site), "random-effects grouping structure")
  # The ordinary null-vs-adjusted PCV (same {stratum}) still works.
  expect_error(suppressMessages(calculate_pcv(m_null, m_adj)), NA)
  # compare_maihda() is warning-based; it must surface the structure mismatch.
  expect_warning(suppressMessages(compare_maihda(m_null, m_site)),
                 "random-effects structure")
})

test_that("longitudinal brms count VPC propagates residual uncertainty (finding 7)", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  set.seed(2)
  np <- 100
  waves <- 0:4
  d <- expand.grid(pid = 1:np, t = waves)
  d$a <- factor(d$pid %% 2)
  d$b <- factor(d$pid %% 3)
  u <- rnorm(np, 0, 0.5)[d$pid]
  d$y <- rpois(nrow(d), exp(0.2 + 0.05 * d$t + u))  # low counts
  d$pid <- factor(d$pid)

  m <- suppressWarnings(suppressMessages(fit_maihda(
    y ~ (1 | a:b), data = d, id = "pid", time = "t", family = poisson,
    engine = "brms", chains = 2, iter = 600, refresh = 0, seed = 1)))

  s <- suppressWarnings(suppressMessages(summary(m)))
  vt <- s$longitudinal$vpc_t
  expect_true(all(is.finite(as.matrix(vt))))

  # The count residual now varies across posterior draws (was a posterior-mean
  # scalar), so its uncertainty enters the VPC band.
  draws <- maihda_posterior_draws_brms(m$model)
  tt <- maihda_lng_time_term(m$longitudinal_info)
  sig_s <- maihda_re_cov_draws_brms(draws, "stratum", tt)
  sig_i <- maihda_re_cov_draws_brms(draws, "pid", tt)
  nd <- m$model$data
  nd[[tt]] <- 2
  eta <- as.matrix(brms::posterior_linpred(m$model, newdata = nd, re_formula = NA))
  v_t <- (sig_s$v0 + 4 * sig_s$cov + 4 * sig_s$v1) +
    (sig_i$v0 + 4 * sig_i$cov + 4 * sig_i$v1)
  res_new <- maihda_count_resid_var_from_linpred(
    eta, v_t, w = maihda_fit_prior_weights(m$model))
  expect_length(res_new, nrow(eta))
  expect_gt(stats::sd(res_new), 0)
})

test_that("contextual binary stepwise_pcv reports AUC/MOR (finding 8)", {
  skip_on_cran()
  set.seed(31)
  n <- 1500
  d <- data.frame(
    a = factor(sample(c("p", "q"), n, TRUE)),
    b = factor(sample(c("x", "y", "z"), n, TRUE)),
    site = factor(sample(paste0("s", 1:8), n, TRUE))
  )
  strat <- interaction(d$a, d$b, drop = TRUE)
  strat_eff <- rnorm(nlevels(strat), 0, 0.8)[strat]   # real between-stratum signal
  site_eff <- rnorm(8, 0, 0.4)[d$site]
  d$y <- rbinom(n, 1, plogis(-0.2 + strat_eff + site_eff))
  ds <- make_strata(d, vars = c("a", "b"))$data

  res <- suppressWarnings(suppressMessages(
    stepwise_pcv(ds, "y", c("a", "b"), context = "site")))
  # The intersectional-scope AUC excludes the context effect, matching the
  # net-of-context PCV columns, so it is no longer suppressed.
  expect_true("AUC" %in% names(res))
  expect_true("MOR" %in% names(res))
})

test_that("plot(type='all') warns on a failed panel instead of dropping it (finding 9)", {
  skip_on_cran()
  set.seed(4)
  n <- 600
  d <- data.frame(
    a = factor(sample(c("p", "q"), n, TRUE)),
    b = factor(sample(c("x", "y", "z"), n, TRUE)),
    z = rnorm(n)
  )
  d$y <- 1 + 0.4 * d$z + rnorm(n)
  an <- suppressWarnings(suppressMessages(maihda(y ~ z + (1 | a:b), data = d)))

  # Force one optional panel to fail; the montage must warn (not silently omit).
  local_mocked_bindings(plot_obs_vs_shrunken = function(...) stop("boom"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(suppressMessages(plot(an, type = "all")), "obs_vs_shrunken")
})
