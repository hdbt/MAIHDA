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

test_that("fit_maihda auto-detects binary on the analytic (complete-case) sample", {
  set.seed(77)
  n <- 200
  d <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    x = rnorm(n)
  )
  # Outcome is 0/1, except two rows valued 2 whose covariate x is missing, so the
  # analytic sample (complete cases) is binary even though the raw column is not.
  y <- rbinom(n, 1, plogis(-0.2 + 0.4 * d$x))
  y[1:2] <- 2L
  d$x[1:2] <- NA_real_
  d$y <- y

  expect_warning(
    m <- fit_maihda(y ~ x + (1 | stratum), data = d),
    "binary", ignore.case = TRUE
  )
  expect_equal(m$family$family, "binomial")
})
