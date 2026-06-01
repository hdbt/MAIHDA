test_that("plot() on a maihda_comparison validates required columns before plotting", {
  bad <- structure(data.frame(vpc = c(0.1, 0.2)),
                   class = c("maihda_comparison", "data.frame"))
  expect_error(plot(bad), "model", fixed = TRUE)

  good <- structure(data.frame(model = c("A", "B"), vpc = c(0.1, 0.2)),
                    class = c("maihda_comparison", "data.frame"))
  expect_s3_class(plot(good), "ggplot")
})

test_that("compare_maihda output is a maihda_comparison and plots via plot()", {
  set.seed(2201)
  d <- data.frame(stratum = rep(seq_len(8), each = 12), x = rnorm(96))
  d$y <- 1 + d$x + rnorm(8, sd = 0.8)[d$stratum] + rnorm(96, sd = 0.4)
  m1 <- fit_maihda(y ~ x + (1 | stratum), data = d)
  m2 <- fit_maihda(y ~ 1 + (1 | stratum), data = d)

  cmp <- compare_maihda(m1, m2)
  expect_s3_class(cmp, "maihda_comparison")
  expect_s3_class(cmp, "data.frame")          # still a data.frame
  expect_s3_class(plot(cmp), "ggplot")
})

test_that("plot_comparison() is deprecated but still works", {
  df <- data.frame(model = c("A", "B"), vpc = c(0.1, 0.2))
  expect_warning(p <- plot_comparison(df), "deprecated")
  expect_s3_class(p, "ggplot")
})

test_that("compare_maihda warns when models use different families", {
  set.seed(2202)
  d <- data.frame(stratum = rep(seq_len(8), each = 25), x = rnorm(200))
  u <- rnorm(8, sd = 0.7)[d$stratum]
  d$y_gauss <- 1 + 0.3 * d$x + u + rnorm(200, sd = 0.4)
  d$y_count <- rpois(200, lambda = exp(0.3 + 0.2 * d$x + u))

  m_gauss <- fit_maihda(y_gauss ~ x + (1 | stratum), data = d)
  m_pois  <- fit_maihda(y_count ~ x + (1 | stratum), data = d, family = "poisson")

  expect_warning(
    compare_maihda(m_gauss, m_pois),
    "different (outcomes|families)"
  )
})
