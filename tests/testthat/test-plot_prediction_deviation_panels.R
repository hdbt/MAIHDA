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

test_that("plot_prediction_deviation_panels handles factor binomial outcomes without coercion warnings", {
  set.seed(456)
  df <- data.frame(
    x = rnorm(100),
    y = factor(sample(c("No", "Yes"), 100, replace = TRUE), levels = c("No", "Yes"))
  )
  m <- glm(y ~ x, data = df, family = binomial)

  expect_warning(
    p <- plot_prediction_deviation_panels(m, df, type = "binomial"),
    NA
  )
  expect_s3_class(p, "patchwork")
  expect_equal(
    MAIHDA:::maihda_binomial_observed_01(df$y, nrow(df)),
    as.integer(df$y == "Yes")
  )
})

test_that("prediction deviation auto-detects brmsfit families", {
  fake_bernoulli <- structure(
    list(family = structure(list(family = "bernoulli", link = "logit"), class = "family")),
    class = "brmsfit"
  )
  fake_gaussian <- structure(
    list(family = structure(list(family = "gaussian", link = "identity"), class = "family")),
    class = "brmsfit"
  )
  fake_ordinal <- structure(
    list(family = structure(list(family = "cumulative", link = "logit"), class = "family")),
    class = "brmsfit"
  )

  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_bernoulli), "binomial")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_gaussian), "gaussian")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_ordinal), "ordinal")
})
