test_that("plot_prediction_deviation_panels works for gaussian", {
  set.seed(123)
  df <- data.frame(
    x = rnorm(100),
    y = rnorm(100)
  )
  m <- lm(y ~ x, data = df)

  p <- plot_prediction_deviation_panels(m, df, type = "gaussian")
  expect_s3_class(p, "patchwork")
})

test_that("plot_prediction_deviation_panels works for binomial", {
  set.seed(123)
  df <- data.frame(
    x = rnorm(100),
    y = sample(0:1, 100, replace = TRUE)
  )
  m <- glm(y ~ x, data = df, family = binomial)

  p <- plot_prediction_deviation_panels(m, df, type = "binomial")
  expect_s3_class(p, "patchwork")
})
