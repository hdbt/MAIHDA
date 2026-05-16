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

test_that("response-scale ternary interaction signal uses response-scale differences", {
  skip_if_not_installed("lme4")

  set.seed(124)
  df <- data.frame(
    stratum = factor(rep(seq_len(8), each = 30)),
    x = rnorm(240)
  )
  stratum_u <- rnorm(8, sd = 0.9)[df$stratum]
  df$y <- rbinom(nrow(df), 1, stats::plogis(-0.2 + 0.6 * df$x + stratum_u))

  model <- fit_maihda(y ~ x + (1 | stratum), data = df, family = "binomial")
  td <- suppressWarnings(compute_maihda_ternary_data(model, scale = "response", verbose = FALSE))

  fe_link <- stats::predict(model$model, newdata = model$data, re.form = NA, type = "link")
  u_by_row <- td$u_j[match(as.character(model$data$stratum), as.character(td$stratum))]
  expected <- stats::aggregate(
    list(
      additive_only = stats::plogis(fe_link),
      full_prediction = stats::plogis(fe_link + u_by_row)
    ),
    by = list(stratum = as.character(model$data$stratum)),
    FUN = mean
  )
  expected$interaction_signal <- abs(expected$full_prediction - expected$additive_only)

  idx <- match(as.character(td$stratum), expected$stratum)
  expect_equal(td$additive_only, expected$additive_only[idx], tolerance = 1e-8)
  expect_equal(td$interaction_signal, expected$interaction_signal[idx], tolerance = 1e-8)
  expect_false(isTRUE(all.equal(td$interaction_signal, abs(td$u_j), tolerance = 1e-4)))
})
