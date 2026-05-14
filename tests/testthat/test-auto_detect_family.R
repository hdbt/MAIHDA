test_that("fit_maihda auto-detects binary outcomes and warns when family missing", {
  data <- maihda_sim_data
  # create binary outcome
  data$bin_out <- ifelse(data$health_outcome > mean(data$health_outcome), 1, 0)

  expect_warning(
    m <- fit_maihda(bin_out ~ age + (1 | gender:race), data = data, engine = "lme4"),
    "The outcome variable appears to be binary. Automatically switching to family = 'binomial'",
    fixed = TRUE
  )

  expect_true(inherits(m, "maihda_model"))
  expect_equal(m$family$family, "binomial")
})
