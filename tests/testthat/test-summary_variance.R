test_that("summary handles binomial and gaussian residual variance correctly", {
  # Setup data
  strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))

  # Create a gaussian model
  mod_gauss <- fit_maihda(health_outcome ~ age + (1 | stratum),
                         data = strata_result$data,
                         engine = "lme4", family = "gaussian")

  # Create a binomial model
  strata_result$data$binary_outcome <- ifelse(strata_result$data$health_outcome > mean(strata_result$data$health_outcome), 1, 0)
  mod_binom <- fit_maihda(binary_outcome ~ age + (1 | stratum),
                         data = strata_result$data,
                         engine = "lme4", family = "binomial")

  # Check summaries
  summ_gauss <- summary(mod_gauss)
  summ_binom <- summary(mod_binom)

  # Gaussian residual variance should not be pi^2/3
  gauss_resid_var <- summ_gauss$variance_components$variance[summ_gauss$variance_components$component == "Within-stratum (residual)"]
  expect_true(abs(gauss_resid_var - (pi^2 / 3)) > 0.1) # Shouldn't match pi^2/3

  # Binomial residual variance should be exactly pi^2/3 for logit link
  binom_resid_var <- summ_binom$variance_components$variance[summ_binom$variance_components$component == "Within-stratum (residual)"]
  expect_equal(binom_resid_var, (pi^2) / 3, tolerance = 1e-6)
})

test_that("summary print reports requested bootstrap confidence level", {
  set.seed(1101)
  d <- data.frame(
    stratum = rep(seq_len(8), each = 8),
    x = rnorm(64)
  )
  d$y <- 1 + d$x + rnorm(8, sd = 1)[d$stratum] + rnorm(64, sd = 0.5)

  model <- fit_maihda(y ~ x + (1 | stratum), data = d)
  summ <- summary(model, bootstrap = TRUE, n_boot = 5, conf_level = 0.80)

  expect_equal(summ$vpc$conf_level, 0.80)
  expect_output(print(summ), "Bootstrap 80% CI", fixed = TRUE)
})

test_that("summary errors clearly when brms bootstrap is requested", {
  fake_brms_model <- structure(
    list(engine = "brms", model = NULL),
    class = "maihda_model"
  )

  expect_error(
    summary(fake_brms_model, bootstrap = TRUE),
    "only supported for lme4",
    fixed = TRUE
  )
})

test_that("brms latent residual variance recognises the bernoulli family", {
  # Unit test of the fix without compiling a Stan model: a bernoulli/logit model
  # should use the latent residual variance pi^2/3, not error as unsupported.
  stub_logit <- structure(
    list(family = structure(list(family = "bernoulli", link = "logit"),
                            class = "family")),
    class = "brmsfit"
  )
  expect_equal(MAIHDA:::maihda_residual_variance_brms(stub_logit), (pi^2) / 3,
               tolerance = 1e-8)

  stub_probit <- structure(
    list(family = structure(list(family = "bernoulli", link = "probit"),
                            class = "family")),
    class = "brmsfit"
  )
  expect_equal(MAIHDA:::maihda_residual_variance_brms(stub_probit), 1)
})

# ---- brms VPC/ICC from posterior draws (Stan-free helper tests) ------------

test_that("maihda_vpc_posterior_summary summarises the per-draw VPC", {
  set.seed(42)
  var_stratum <- rgamma(2000, shape = 2, rate = 1)
  var_residual <- rep((pi^2) / 3, 2000)

  res <- MAIHDA:::maihda_vpc_posterior_summary(
    var_stratum, 0, var_residual, conf_level = 0.90
  )

  vpc_draws <- var_stratum / (var_stratum + (pi^2) / 3)
  expect_equal(res$estimate, stats::median(vpc_draws))
  expect_equal(res$ci_lower, stats::quantile(vpc_draws, 0.05, names = FALSE))
  expect_equal(res$ci_upper, stats::quantile(vpc_draws, 0.95, names = FALSE))
  expect_equal(res$conf_level, 0.90)
  expect_equal(res$n_draws, 2000L)
  expect_true(res$ci_lower < res$estimate && res$estimate < res$ci_upper)

  # point = "mean" returns the posterior mean of the per-draw VPC.
  res_mean <- MAIHDA:::maihda_vpc_posterior_summary(
    var_stratum, 0, var_residual, point = "mean"
  )
  expect_equal(res_mean$estimate, mean(vpc_draws))
})

test_that("maihda_vpc_posterior_summary handles degenerate and invalid input", {
  # A VPC that is constant across draws collapses to a point interval.
  const <- MAIHDA:::maihda_vpc_posterior_summary(rep(1, 100), 0, rep(1, 100))
  expect_equal(const$estimate, 0.5)
  expect_equal(const$ci_lower, 0.5)
  expect_equal(const$ci_upper, 0.5)

  # Vectors of incompatible (non-scalar) length are rejected.
  expect_error(
    MAIHDA:::maihda_vpc_posterior_summary(c(1, 2, 3), c(1, 2), 1),
    "common length", fixed = TRUE
  )

  # All-zero variances give 0/0 = NaN draws, which are dropped, leaving nothing.
  expect_error(
    MAIHDA:::maihda_vpc_posterior_summary(0, 0, 0),
    "No finite VPC draws", fixed = TRUE
  )
})

test_that("maihda_random_variance_draws_brms extracts per-draw variances from a draws frame", {
  draws <- data.frame(
    b_Intercept = rnorm(5),
    sd_stratum__Intercept = c(0.5, 1, 1.5, 2, 2.5),
    sd_region__Intercept = c(0.1, 0.2, 0.3, 0.4, 0.5),
    sigma = rep(1, 5)
  )

  rv <- MAIHDA:::maihda_random_variance_draws_brms(draws)
  expect_equal(rv$stratum, c(0.5, 1, 1.5, 2, 2.5)^2)
  expect_equal(rv$total, c(0.5, 1, 1.5, 2, 2.5)^2 + c(0.1, 0.2, 0.3, 0.4, 0.5)^2)
  expect_equal(rv$other, c(0.1, 0.2, 0.3, 0.4, 0.5)^2)
})

test_that("maihda_random_variance_draws_brms reports zero other variance for a single group", {
  rv <- MAIHDA:::maihda_random_variance_draws_brms(
    data.frame(sd_stratum__Intercept = c(1, 2, 3))
  )
  expect_equal(rv$stratum, c(1, 4, 9))
  expect_equal(rv$total, c(1, 4, 9))
  expect_true(all(rv$other == 0))
})

test_that("maihda_random_variance_draws_brms errors on missing or random-slope SD columns", {
  expect_error(
    MAIHDA:::maihda_random_variance_draws_brms(data.frame(b_Intercept = 1:3)),
    "No random-effect standard-deviation draws", fixed = TRUE
  )
  expect_error(
    MAIHDA:::maihda_random_variance_draws_brms(data.frame(sd_region__Intercept = 1:3)),
    "No 'stratum' random-effect SD draws", fixed = TRUE
  )
  expect_error(
    MAIHDA:::maihda_random_variance_draws_brms(
      data.frame(sd_stratum__Intercept = c(1, 2, 3), sd_stratum__x = c(0.1, 0.2, 0.3))
    ),
    "random slopes", fixed = TRUE
  )
})

test_that("maihda_residual_variance_draws_brms returns per-draw residual variance by family", {
  draws <- data.frame(sd_stratum__Intercept = c(1, 2), sigma = c(2, 4))

  gauss <- structure(
    list(family = structure(list(family = "gaussian", link = "identity"),
                            class = "family")),
    class = "brmsfit"
  )
  expect_equal(MAIHDA:::maihda_residual_variance_draws_brms(gauss, draws), c(4, 16))

  logit <- structure(
    list(family = structure(list(family = "bernoulli", link = "logit"),
                            class = "family")),
    class = "brmsfit"
  )
  expect_equal(MAIHDA:::maihda_residual_variance_draws_brms(logit, draws),
               rep((pi^2) / 3, 2))

  probit <- structure(
    list(family = structure(list(family = "binomial", link = "probit"),
                            class = "family")),
    class = "brmsfit"
  )
  expect_equal(MAIHDA:::maihda_residual_variance_draws_brms(probit, draws), c(1, 1))
})

test_that("print.maihda_summary labels the brms posterior credible interval", {
  fake <- structure(
    list(
      vpc = list(estimate = 0.2, ci_lower = 0.1, ci_upper = 0.3,
                 conf_level = 0.95, bootstrap = FALSE, method = "posterior"),
      variance_components = data.frame(
        component = c("Between-stratum (random)", "Within-stratum (residual)", "Total"),
        variance = c(1, 3, 4),
        sd = sqrt(c(1, 3, 4)),
        proportion = c(0.25, 0.75, 1)
      ),
      fixed_effects = data.frame(term = "(Intercept)", estimate = 0),
      stratum_estimates = NULL,
      model_summary = NULL,
      engine = "brms"
    ),
    class = "maihda_summary"
  )

  out <- capture.output(print(fake))
  expect_true(any(grepl("0.2000 [0.1000, 0.3000]", out, fixed = TRUE)))
  expect_true(any(grepl("Posterior 95% credible interval", out, fixed = TRUE)))
})

test_that("brms summary returns a posterior credible interval for the VPC/ICC", {
  # Compiles a Stan model, so skip on CI/CRAN (toolchain is unreliable there);
  # the draws-based VPC logic is covered Stan-free by the helper tests above.
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  set.seed(2024)
  d <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    x = rnorm(200)
  )
  d$y <- 1 + 0.5 * d$x + rnorm(8, sd = 0.8)[d$stratum] + rnorm(200, sd = 1)

  model <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, engine = "brms",
               family = "gaussian", chains = 2, iter = 500, refresh = 0, seed = 1)
  ))

  summ <- suppressWarnings(summary(model))

  expect_false(summ$vpc$bootstrap)
  expect_identical(summ$vpc$method, "posterior")
  expect_equal(summ$vpc$conf_level, 0.95)
  expect_true(is.finite(summ$vpc$estimate))
  expect_true(summ$vpc$ci_lower < summ$vpc$ci_upper)
  expect_true(summ$vpc$ci_lower >= 0 && summ$vpc$ci_upper <= 1)
  expect_true(summ$vpc$estimate >= summ$vpc$ci_lower &&
              summ$vpc$estimate <= summ$vpc$ci_upper)
  expect_output(print(summ), "Posterior 95% credible interval", fixed = TRUE)

  # Lowering the confidence level narrows (nests inside) the 95% interval since
  # both are quantiles of the same posterior VPC draws.
  summ80 <- suppressWarnings(summary(model, conf_level = 0.80))
  expect_equal(summ80$vpc$conf_level, 0.80)
  expect_true(summ80$vpc$ci_lower >= summ$vpc$ci_lower - 1e-8)
  expect_true(summ80$vpc$ci_upper <= summ$vpc$ci_upper + 1e-8)

  # bootstrap = TRUE remains rejected for brms with a clear message.
  expect_error(summary(model, bootstrap = TRUE), "only supported for lme4", fixed = TRUE)
})

test_that("summary validates bootstrap arguments before simulation", {
  fake_model <- structure(
    list(engine = "lme4", model = NULL),
    class = "maihda_model"
  )

  expect_error(
    summary(fake_model, bootstrap = c(TRUE, FALSE)),
    "'bootstrap' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    summary(fake_model, bootstrap = TRUE, n_boot = 0),
    "'n_boot' must be a single positive whole number",
    fixed = TRUE
  )
  expect_error(
    summary(fake_model, bootstrap = TRUE, conf_level = 1),
    "'conf_level' must be a single number between 0 and 1",
    fixed = TRUE
  )
})
