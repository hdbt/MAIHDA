test_that("maihda_sim_data is available and has correct structure", {
  data(maihda_sim_data, envir = environment())
  
  # Check structure
  expect_true(is.data.frame(maihda_sim_data))
  expect_equal(nrow(maihda_sim_data), 500)
  expect_equal(ncol(maihda_sim_data), 6)
  
  # Check column names
  expected_cols <- c("id", "gender", "race", "education", "age", "health_outcome")
  expect_equal(names(maihda_sim_data), expected_cols)
  
  # Check data types
  expect_true(is.numeric(maihda_sim_data$id))
  expect_true(is.character(maihda_sim_data$gender))
  expect_true(is.character(maihda_sim_data$race))
  expect_true(is.character(maihda_sim_data$education))
  expect_true(is.numeric(maihda_sim_data$age))
  expect_true(is.numeric(maihda_sim_data$health_outcome))
})
