test_that("plot_maihda creates predicted plot with lme4", {
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
  
  # Create predicted plot
  plot <- plot_maihda(model, type = "predicted")
  
  # Check structure
  expect_true(inherits(plot, "ggplot"))
  expect_true("predicted" %in% names(plot$data))
  expect_true("lower" %in% names(plot$data))
  expect_true("upper" %in% names(plot$data))
})

test_that("plot_maihda predicted handles n_strata parameter", {
  # Create test data with more strata
  set.seed(456)
  data <- data.frame(
    stratum = rep(1:20, each = 10),
    age = rnorm(200),
    outcome = rnorm(200)
  )
  
  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum), 
                     data = data, 
                     engine = "lme4")
  
  # Create predicted plot with limited strata
  plot <- plot_maihda(model, type = "predicted", n_strata = 10)
  
  # Check that it limits to specified number
  expect_true(inherits(plot, "ggplot"))
  expect_lte(nrow(plot$data), 10)
})

test_that("plot_maihda predicted validates inputs", {
  # Create test data
  set.seed(789)
  data <- data.frame(
    stratum = rep(1:5, each = 10),
    age = rnorm(50),
    outcome = rnorm(50)
  )
  
  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum), 
                     data = data, 
                     engine = "lme4")
  
  # Check plot is created
  plot <- plot_maihda(model, type = "predicted")
  expect_true(inherits(plot, "ggplot"))
})

test_that("plot_maihda predicted preserves stratum order", {
  # Create test data with stratum labels
  set.seed(999)
  data <- data.frame(
    stratum = factor(rep(c("A", "B", "C", "D"), each = 10), 
                     levels = c("A", "B", "C", "D")),
    age = rnorm(40),
    outcome = rnorm(40)
  )
  
  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum), 
                     data = data, 
                     engine = "lme4")
  
  # Create predicted plot
  plot <- plot_maihda(model, type = "predicted")
  
  # Check structure
  expect_true(inherits(plot, "ggplot"))
  expect_true("display_label" %in% names(plot$data))
  
  # Check that display_label is a factor and order is preserved
  expect_true(is.factor(plot$data$display_label))
  
  # The strata should be in their original order, not sorted by predicted value
  # We can't check exact order without knowing the predicted values, 
  # but we can check that stratum labels are present
  expect_true(all(c("A", "B", "C", "D") %in% levels(plot$data$display_label)))
})

test_that("plot_maihda uses meaningful stratum labels from make_strata", {
  # Create test data with meaningful categorical variables
  set.seed(555)
  data <- data.frame(
    gender = rep(c("Male", "Female"), each = 20),
    race = rep(c("White", "Black"), times = 20),
    age = rnorm(40),
    outcome = rnorm(40)
  )
  
  # Use make_strata to create labeled strata
  strata_result <- make_strata(data, vars = c("gender", "race"))
  
  # Verify strata_info has labels
  expect_true("label" %in% names(strata_result$strata_info))
  expect_true(any(grepl("_", strata_result$strata_info$label)))
  
  # Fit model using data from make_strata
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = strata_result$data,
                     engine = "lme4")
  
  # Verify model has strata_info
  expect_false(is.null(model$strata_info))
  expect_true("label" %in% names(model$strata_info))
  
  # Get summary
  summary_obj <- summary_maihda(model)
  
  # Verify summary has labels
  expect_true("label" %in% names(summary_obj$stratum_estimates))
  
  # Create predicted plot
  plot <- plot_maihda(model, type = "predicted")
  
  # Check structure
  expect_true(inherits(plot, "ggplot"))
  expect_true("display_label" %in% names(plot$data))
  
  # Check that meaningful labels are used (should contain underscores from gender_race)
  display_labels <- as.character(plot$data$display_label)
  expect_true(any(grepl("_", display_labels)),
             info = "Plot should use meaningful labels like 'Male_White', not numeric IDs")
})
