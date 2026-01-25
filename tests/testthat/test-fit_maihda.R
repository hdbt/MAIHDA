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
