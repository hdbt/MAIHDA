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

test_that("plot_prediction_deviation_panels aligns gaussian predictions to supplied data", {
  set.seed(124)
  df <- data.frame(
    x = rnorm(100),
    y = rnorm(100)
  )
  df$x[1:3] <- NA_real_
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

test_that("plot_prediction_deviation_panels aligns binomial predictions to supplied data", {
  set.seed(125)
  df <- data.frame(
    x = rnorm(100),
    y = sample(0:1, 100, replace = TRUE)
  )
  df$x[1:3] <- NA_real_
  m <- glm(y ~ x, data = df, family = binomial)

  p <- plot_prediction_deviation_panels(m, df, type = "binomial")

  expect_s3_class(p, "patchwork")
})

test_that("binomial fallback residuals stay on the deviance scale", {
  fitted <- c(0.2, 0.8, NA_real_, 0.5)
  obs <- c(0L, 1L, 1L, 0L)

  resids <- MAIHDA:::maihda_prediction_panel_binomial_residuals(
    structure(list(), class = "no_residual_method"),
    data.frame(id = seq_along(fitted)),
    fitted,
    obs
  )

  expected <- c(
    sqrt(-2 * log1p(-0.2)),
    sqrt(-2 * log(0.8)),
    0,
    sqrt(-2 * log1p(-0.5))
  )
  expect_equal(resids, expected)
})

test_that("binomial residuals align to supplied row order", {
  set.seed(127)
  df <- data.frame(x = rnorm(80))
  df$y <- rbinom(80, 1, plogis(-0.2 + 0.8 * df$x))
  m <- glm(y ~ x, data = df, family = binomial)
  shuffled <- df[sample(seq_len(nrow(df))), , drop = FALSE]

  fitted <- MAIHDA:::maihda_prediction_panel_fitted(m, shuffled, "binomial")$fit
  obs <- MAIHDA:::maihda_binomial_observed_01(shuffled$y, nrow(shuffled))
  resids <- MAIHDA:::maihda_prediction_panel_binomial_residuals(
    m,
    shuffled,
    fitted,
    obs
  )
  expected <- MAIHDA:::maihda_binomial_abs_deviance_residual(obs, fitted)

  expect_equal(resids, expected)
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

test_that("plot_prediction_deviation_panels aligns ordinal predictions to supplied data", {
  skip_if_not_installed("MASS")

  set.seed(126)
  df <- data.frame(
    x = rnorm(100),
    y = ordered(
      sample(c("low", "mid", "high"), 100, replace = TRUE),
      levels = c("low", "mid", "high")
    )
  )
  df$x[1:3] <- NA_real_
  m <- MASS::polr(y ~ x, data = df, Hess = TRUE)

  p <- plot_prediction_deviation_panels(m, df, type = "ordinal")
  p_expected <- plot_prediction_deviation_panels(
    m,
    df,
    type = "ordinal",
    ordinal_mode = "expected_score"
  )

  expect_s3_class(p, "patchwork")
  expect_s3_class(p_expected, "patchwork")
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
  fake_poisson <- structure(
    list(family = structure(list(family = "poisson", link = "log"), class = "family")),
    class = "brmsfit"
  )

  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_bernoulli), "binomial")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_gaussian), "gaussian")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_ordinal), "ordinal")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_poisson), "poisson")
})

test_that("Poisson prediction panels use response-scale (count) fitted values", {
  set.seed(2005)
  df <- data.frame(x = rnorm(150))
  df$y <- rpois(150, lambda = exp(1 + 0.3 * df$x))
  m <- glm(y ~ x, data = df, family = poisson)

  # Poisson is no longer routed through the Gaussian (link-scale) branch.
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(m), "poisson")

  fitted_used <- MAIHDA:::maihda_prediction_panel_fitted(m, df, "poisson")$fit
  expect_equal(fitted_used, unname(predict(m, type = "response")), tolerance = 1e-8)
  # ...and NOT the link (log) scale the old Gaussian routing produced.
  expect_false(isTRUE(all.equal(fitted_used, unname(predict(m, type = "link")))))

  expect_s3_class(plot_prediction_deviation_panels(m, df), "patchwork")
})
