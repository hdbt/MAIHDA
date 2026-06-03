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

test_that("maihda_response_is_binary distinguishes Bernoulli from aggregated binomial", {
  d <- data.frame(
    y01  = rep(0:1, 10),
    yf   = factor(rep(c("no", "yes"), 10)),
    cnt  = rep(0:4, 4),
    succ = rep(3L, 20),
    fail = rep(7L, 20),
    n    = rep(10L, 20),
    x    = rnorm(20)
  )
  # Two-level vector responses -> Bernoulli
  expect_true(MAIHDA:::maihda_response_is_binary(y01 ~ x, d))
  expect_true(MAIHDA:::maihda_response_is_binary(yf ~ x, d))
  # Counts and aggregated binomial responses must NOT be treated as Bernoulli
  expect_false(MAIHDA:::maihda_response_is_binary(cnt ~ x, d))
  expect_false(MAIHDA:::maihda_response_is_binary(cbind(succ, fail) ~ x, d))
  expect_false(MAIHDA:::maihda_response_is_binary(succ | trials(n) ~ x, d))
})

test_that("maihda_response_is_binary uses the analytic (complete-case) sample", {
  d <- data.frame(
    y = c(rep(0:1, 9), 2L, 2L),   # full column has a spurious third level
    x = c(rnorm(18), NA, NA),     # ...only in rows that are incomplete on x
    stratum = factor(rep(seq_len(2), 10))
  )
  # Full column is not binary, but the analytic sample (complete x) is 0/1.
  expect_false(MAIHDA:::maihda_is_binary_vector(d$y))
  expect_true(MAIHDA:::maihda_response_is_binary(y ~ x + (1 | stratum), d))
})

test_that("maihda_is_colon_interaction accepts only symbols and colon interactions", {
  expect_true(MAIHDA:::maihda_is_colon_interaction(quote(a)))
  expect_true(MAIHDA:::maihda_is_colon_interaction(quote(a:b)))
  expect_true(MAIHDA:::maihda_is_colon_interaction(quote(a:b:c)))
  expect_false(MAIHDA:::maihda_is_colon_interaction(quote(interaction(a, b))))
  expect_false(MAIHDA:::maihda_is_colon_interaction(quote(paste(a, b))))
  expect_false(MAIHDA:::maihda_is_colon_interaction(quote(cut(age, 3))))
  expect_false(MAIHDA:::maihda_is_colon_interaction(quote(f(x):b)))
})

test_that("fit_maihda rejects function-call grouping terms in automatic strata", {
  set.seed(5)
  d <- data.frame(
    g = sample(c("a", "b"), 80, replace = TRUE),
    r = sample(c("x", "y"), 80, replace = TRUE),
    age = rnorm(80)
  )
  d$y <- rnorm(80)

  expect_error(
    fit_maihda(y ~ age + (1 | interaction(g, r)), data = d),
    "make_strata", fixed = TRUE
  )
  # The colon-interaction shorthand still works.
  expect_s3_class(fit_maihda(y ~ age + (1 | g:r), data = d), "maihda_model")
})

test_that("fit_maihda fits a binary outcome with brms via bernoulli()", {
  # Compiles a Stan model, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it
  # (otherwise a plain local test() would hang on Stan compilation). The bernoulli
  # residual-variance fix is covered Stan-free in test-summary_variance.R.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
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
