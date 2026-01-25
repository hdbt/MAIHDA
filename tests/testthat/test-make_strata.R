test_that("make_strata creates strata correctly", {
  # Create test data
  data <- data.frame(
    gender = rep(c("M", "F"), each = 50),
    race = rep(c("White", "Black"), 50),
    outcome = rnorm(100)
  )
  
  # Create strata
  result <- make_strata(data, vars = c("gender", "race"))
  
  # Check structure
  expect_true(inherits(result, "maihda_strata"))
  expect_true("stratum" %in% names(result$data))
  expect_equal(nrow(result$data), 100)
  expect_equal(nrow(result$strata_info), 4)
  
  # Check that all observations have valid strata
  expect_equal(sum(is.na(result$data$stratum)), 0)
})

test_that("make_strata handles minimum count threshold", {
  data <- data.frame(
    gender = c(rep("M", 50), rep("F", 50)),
    race = c(rep("White", 48), rep("Black", 2), rep("White", 48), rep("Black", 2)),
    outcome = rnorm(100)
  )
  
  # Create strata with minimum count
  result <- make_strata(data, vars = c("gender", "race"), min_n = 10)
  
  # Check that small strata are excluded
  expect_equal(nrow(result$strata_info), 2)
  expect_equal(sum(is.na(result$data$stratum)), 4)
})

test_that("make_strata handles errors correctly", {
  data <- data.frame(x = 1:10, y = 1:10)
  
  # Missing variables
  expect_error(make_strata(data, vars = c("a", "b")), 
               "Variables not found")
  
  # Invalid data type
  expect_error(make_strata("not a data frame", vars = "x"),
               "'data' must be a data frame")
  
  # Empty vars
  expect_error(make_strata(data, vars = character(0)),
               "at least one variable name")
})

test_that("make_strata handles missing values correctly", {
  # Create test data with missing values
  data <- data.frame(
    gender = c("M", "F", "M", "F", "M", NA, "F", "M"),
    race = c("White", "Black", "White", NA, "Black", "White", "Black", "White"),
    outcome = rnorm(8)
  )
  
  # Create strata
  result <- make_strata(data, vars = c("gender", "race"))
  
  # Check structure
  expect_true(inherits(result, "maihda_strata"))
  expect_true("stratum" %in% names(result$data))
  expect_equal(nrow(result$data), 8)
  
  # Check that rows with missing values have NA stratum
  # Row 4: race is NA
  expect_true(is.na(result$data$stratum[4]))
  # Row 6: gender is NA
  expect_true(is.na(result$data$stratum[6]))
  
  # Check that rows without missing values have valid strata
  expect_false(is.na(result$data$stratum[1]))
  expect_false(is.na(result$data$stratum[2]))
  expect_false(is.na(result$data$stratum[3]))
  
  # Check that we have 3 strata (M_White, M_Black, F_Black)
  # Note: F_White doesn't exist in the complete cases
  expect_equal(nrow(result$strata_info), 3)
  
  # Check that total observations with missing values are marked as NA
  # Rows 4 and 6 should have NA stratum
  expect_equal(sum(is.na(result$data$stratum)), 2)
})

test_that("make_strata handles all missing values in one variable", {
  # Create test data where one variable is entirely missing
  data <- data.frame(
    gender = c("M", "F", "M", "F"),
    race = c(NA, NA, NA, NA),
    outcome = rnorm(4)
  )
  
  # Create strata
  result <- make_strata(data, vars = c("gender", "race"))
  
  # All observations should have NA stratum
  expect_equal(sum(is.na(result$data$stratum)), 4)
  
  # No valid strata should be created
  expect_equal(nrow(result$strata_info), 0)
})
