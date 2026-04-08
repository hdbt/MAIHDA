test_that("maihda_health_data loads and is formatted correctly", {
  expect_true("maihda_health_data" %in% data(package = "MAIHDA")$results[, "Item"])
  expect_s3_class(maihda_health_data, "data.frame")
  expect_equal(ncol(maihda_health_data), 7)
  expect_true(all(c("BMI", "Obese", "Age", "Gender", "Race", "Education", "Poverty") %in% names(maihda_health_data)))
})

test_that("make_strata works on maihda_health_data", {
  strata_result <- make_strata(maihda_health_data, vars = c("Gender", "Race", "Education"))
  expect_s3_class(strata_result, "maihda_strata")
  expect_true("stratum" %in% names(strata_result$data))
  expect_true(nrow(strata_result$strata_info) > 0)
})

test_that("fit_maihda works on maihda_health_data summary", {
  # Subset further just to make the test run fast
  health_subset <- maihda_health_data[1:200, ]
  strata_result <- make_strata(health_subset, vars = c("Gender", "Race"))

  model <- fit_maihda(BMI ~ Age + (1 | stratum), data = strata_result$data, engine = "lme4")
  expect_s3_class(model, "maihda_model")

  summ <- summary_maihda(model)
  expect_s3_class(summ, "maihda_summary")
  expect_true(!is.na(summ$vpc$estimate))
})
