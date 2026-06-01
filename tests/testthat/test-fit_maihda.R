test_that("fit_maihda works with lme4", {
  # Create test data
  set.seed(123)
  data <- data.frame(
    stratum = rep(1:10, each = 10),
    age = rnorm(100),
    outcome = rnorm(100)
  )

  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = data,
                     engine = "lme4")

  # Check structure
  expect_true(inherits(model, "maihda_model"))
  expect_equal(model$engine, "lme4")
  expect_true(inherits(model$model, "lmerMod"))
})

test_that("fit_maihda handles different families", {
  # Create test data for binomial
  set.seed(123)
  data <- data.frame(
    stratum = rep(1:10, each = 10),
    age = rnorm(100),
    outcome = rbinom(100, 1, 0.5)
  )

  # Fit binomial model
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = data,
                     engine = "lme4",
                     family = "binomial")

  expect_true(inherits(model, "maihda_model"))
  expect_equal(model$family$family, "binomial")
})

test_that("fit_maihda accepts family constructor functions", {
  set.seed(124)
  data <- data.frame(
    stratum = rep(1:10, each = 12),
    age = rnorm(120),
    outcome = rbinom(120, 1, 0.5)
  )

  model <- fit_maihda(outcome ~ age + (1 | stratum),
                      data = data,
                      engine = "lme4",
                      family = binomial)

  expect_true(inherits(model, "maihda_model"))
  expect_equal(model$family$family, "binomial")
  expect_error(
    fit_maihda(outcome ~ age + (1 | stratum), data = data, family = list()),
    "family name, family object, or family function",
    fixed = TRUE
  )
})

test_that("fit_maihda creates strata automatically when interaction is passed", {
  set.seed(123)
  data <- data.frame(
    gender = sample(c("M", "F"), 100, replace = TRUE),
    race = sample(c("W", "B"), 100, replace = TRUE),
    age = rnorm(100),
    outcome = rnorm(100)
  )

  # Using older manual strata way to ensure they match
  strata_result <- make_strata(data, c("gender", "race"))
  model1 <- fit_maihda(outcome ~ age + (1 | stratum), data = strata_result$data)

  # Using auto strata way
  model2 <- fit_maihda(outcome ~ age + (1 | gender:race), data = data)

  expect_equal(summary(model1), summary(model2))
  expect_true(!is.null(model2$strata_info))
})

test_that("fit_maihda fits a binary outcome with brms via bernoulli()", {
  # Compiles a Stan model, so skip on CI/CRAN (toolchain is unreliable there);
  # the bernoulli residual-variance fix is covered without Stan in
  # test-summary_variance.R. This runs locally when brms + a compiler are present.
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  set.seed(321)
  d <- data.frame(
    stratum = factor(rep(seq_len(6), each = 30)),
    x = rnorm(180)
  )
  d$y <- rbinom(180, 1, plogis(-0.2 + 0.4 * d$x + rnorm(6, sd = 0.5)[d$stratum]))

  # Tiny chains -> Stan convergence warnings, irrelevant to what we assert
  # (family routing + the deterministic latent residual variance).
  model <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, engine = "brms",
               family = "binomial", chains = 1, iter = 200, refresh = 0)
  ))
  # binomial 0/1 must be routed to bernoulli, and summary must not error
  expect_equal(model$family$family, "bernoulli")
  summ <- suppressWarnings(summary(model))
  resid_var <- summ$variance_components$variance[
    summ$variance_components$component == "Within-stratum (residual)"
  ]
  expect_equal(resid_var, (pi^2) / 3, tolerance = 1e-6)
})

test_that("fit_maihda validates inputs", {
  data <- data.frame(x = 1:10, y = 1:10)

  # Invalid formula
  expect_error(fit_maihda("not a formula", data = data),
               "must be a formula")

  # Invalid data
  expect_error(fit_maihda(y ~ x, data = "not a data frame"),
               "must be a data frame")

  # Invalid engine
  expect_error(fit_maihda(y ~ x, data = data, engine = "invalid"),
               "should be one of")
})
