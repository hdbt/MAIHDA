test_that("Ternary plot functions work", {
  skip_if_not_installed("lme4")

  # Minimal dummy data
  set.seed(123)
  df <- data.frame(
    stratum = rep(letters[1:5], each = 10),
    y = rnorm(50),
    x = rnorm(50)
  )

  library(lme4)
  fit <- lmer(y ~ x + (1 | stratum), data = df)

  model <- list(
    model = fit,
    engine = "lme4",
    data = df
  )
  class(model) <- "maihda_model"

  # Test compute function
  td <- compute_maihda_ternary_data(model)
  expect_s3_class(td, "tbl_df")
  expect_true(all(c("stratum", "additive_prop", "interaction_prop", "uncertainty_prop") %in% names(td)))
  expect_equal(rowSums(td[, c("additive_prop", "interaction_prop", "uncertainty_prop")]), rep(1, 5), tolerance = 1e-4)

  td_ci <- compute_maihda_ternary_data(model, uncertainty_method = "ci_width")
  expect_equal(td_ci$uncertainty, td$uncertainty * 3.92, tolerance = 1e-8)
  expect_error(
    compute_maihda_ternary_data(model, uncertainty_method = "posterior_sd"),
    "only available for brms"
  )

  # Test wrapper function
  skip_if_not_installed("ggtern")
  out <- maihda_ternary_plot(model)
  expect_type(out, "list")
  expect_true(!is.null(out$plot))
  expect_s3_class(out$data, "tbl_df")
})
