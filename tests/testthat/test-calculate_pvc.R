test_that("calculate_pvc works with basic models", {
  # Create test data with actual stratum effects
  set.seed(123)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 2)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit two models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum), 
                       data = data, 
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PVC
  pvc_result <- calculate_pvc(model1, model2)
  
  # Check structure
  expect_true(inherits(pvc_result, "pvc_result"))
  expect_true(is.numeric(pvc_result$pvc))
  expect_true(is.numeric(pvc_result$var_model1))
  expect_true(is.numeric(pvc_result$var_model2))
  expect_false(pvc_result$bootstrap)
  
  # Check that variances are positive
  expect_true(pvc_result$var_model1 > 0)
  expect_true(pvc_result$var_model2 > 0)
})

test_that("calculate_pvc works with bootstrap", {
  # Create test data with actual stratum effects
  set.seed(456)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1.5)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit two models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PVC with bootstrap (small number for testing)
  pvc_result <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 50)
  
  # Check structure
  expect_true(inherits(pvc_result, "pvc_result"))
  expect_true(pvc_result$bootstrap)
  expect_true(!is.null(pvc_result$ci_lower))
  expect_true(!is.null(pvc_result$ci_upper))
  expect_true(is.numeric(pvc_result$ci_lower))
  expect_true(is.numeric(pvc_result$ci_upper))
  
  # CI should be in reasonable range
  expect_true(pvc_result$ci_lower <= pvc_result$pvc)
  expect_true(pvc_result$ci_upper >= pvc_result$pvc)
})

test_that("calculate_pvc validates inputs", {
  # Create test data with actual stratum effects
  set.seed(789)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Invalid first argument
  expect_error(calculate_pvc("not a model", model1),
               "must be a maihda_model")
  
  # Invalid second argument
  expect_error(calculate_pvc(model1, "not a model"),
               "must be a maihda_model")
  
  # Both arguments invalid
  expect_error(calculate_pvc(data, data),
               "must be a maihda_model")
})

test_that("calculate_pvc handles same model comparison", {
  # Create test data with actual stratum effects
  set.seed(111)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit same model twice
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PVC
  pvc_result <- calculate_pvc(model1, model2)
  
  # PVC should be very close to 0 (same model structure)
  expect_true(abs(pvc_result$pvc) < 0.1)
})

test_that("calculate_pvc calculates correct direction", {
  # Create test data with known structure
  set.seed(222)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 2)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Model without gender (should have more between-stratum variance)
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Model with gender (might reduce between-stratum variance if gender explains some)
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PVC
  pvc_result <- calculate_pvc(model1, model2)
  
  # Variances should be positive
  expect_true(pvc_result$var_model1 > 0)
  expect_true(pvc_result$var_model2 > 0)
  
  # PVC formula: (var1 - var2) / var1
  expected_pvc <- (pvc_result$var_model1 - pvc_result$var_model2) / pvc_result$var_model1
  expect_equal(pvc_result$pvc, expected_pvc)
})

test_that("calculate_pvc print method works", {
  # Create test data with actual stratum effects
  set.seed(333)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1.5)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PVC
  pvc_result <- calculate_pvc(model1, model2)
  
  # Print should work without error
  expect_output(print(pvc_result), "Proportional Change in Variance")
  expect_output(print(pvc_result), "PVC:")
  expect_output(print(pvc_result), "Between-stratum variance:")
})

test_that("calculate_pvc handles binomial models", {
  # Create test data for binomial with stratum effects on logit scale
  set.seed(444)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects on logit scale
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Generate outcome on logit scale with stratum effects
  logit_p <- -0.5 + 0.3 * data$age + stratum_effects[data$stratum]
  prob <- plogis(logit_p)
  data$outcome <- rbinom(nrow(data), 1, prob)
  
  # Fit binomial models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4",
                       family = "binomial")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4",
                       family = "binomial")
  
  # Calculate PVC
  pvc_result <- calculate_pvc(model1, model2)
  
  # Check structure
  expect_true(inherits(pvc_result, "pvc_result"))
  expect_true(is.numeric(pvc_result$pvc))
  expect_true(pvc_result$var_model1 > 0)
  expect_true(pvc_result$var_model2 > 0)
})

test_that("calculate_pvc handles zero variance error", {
  # Create test data with NO stratum effects (will result in singular fit)
  set.seed(555)
  data <- data.frame(
    stratum = rep(1:10, each = 10),
    age = rnorm(100),
    outcome = 5 + 0.5 * rnorm(100)  # No stratum effects
  )
  
  # Fit models (may have singular fit warnings)
  suppressWarnings({
    model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                         data = data,
                         engine = "lme4")
    model2 <- fit_maihda(outcome ~ 1 + (1 | stratum),
                         data = data,
                         engine = "lme4")
  })
  
  # Should error due to zero/negative variance
  expect_error(calculate_pvc(model1, model2),
               "Between-stratum variance")
})
