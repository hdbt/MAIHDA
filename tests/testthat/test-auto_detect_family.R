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

test_that("fit_maihda detects binary on the post-transformation analytic frame", {
  set.seed(2002)
  n <- 300
  d <- data.frame(
    stratum = factor(rep(seq_len(10), each = 30)),
    # 15 non-positive x make log(x) NaN, so lme4 drops those rows.
    x = c(rep(-1, 15), runif(n - 15, 0.2, 4))
  )
  y <- rbinom(n, 1, 0.45)
  y[1:15] <- 2L          # a spurious third level, only on the log(x)-dropped rows
  d$y <- y

  # The raw column has three values, but the analytic frame (after log(x) drops
  # the non-positive rows) is 0/1, so the family must switch to binomial.
  expect_warning(
    m <- fit_maihda(y ~ log(x) + (1 | stratum), data = d),
    "binary", ignore.case = TRUE
  )
  expect_equal(m$family$family, "binomial")
})

test_that("the maihda() wrappers pick the engine from the analytic sample, not raw", {
  # An ordered 3-level outcome subset to two observed levels is a BINARY analytic
  # sample. The wrappers must detect that on the post-subset frame (engine = lme4,
  # binomial) instead of pinning engine = "ordinal" from the raw 3-level column,
  # which would error with an engine/family contradiction.
  set.seed(303)
  n <- 600
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("A", "B", "C"), n, replace = TRUE),
    x = rnorm(n)
  )
  lat <- 0.5 * d$x + rnorm(n)
  d$y <- factor(cut(lat, c(-Inf, -0.3, 0.6, Inf), labels = c("lo", "mid", "hi")),
                levels = c("lo", "mid", "hi"), ordered = TRUE)
  d$grp <- sample(c("g1", "g2"), n, replace = TRUE)

  # Baseline: a direct fit resolves binomial/lme4 on the analytic sample.
  fm <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | gender:race), data = d, subset = y %in% c("lo", "mid"))))
  expect_identical(fm$engine, "lme4")
  expect_identical(fm$family$family, "binomial")

  # maihda() must agree -- and the two-model decomposition must complete (the
  # response-referencing subset is re-used by value, immune to the 0/1 recoding).
  mh <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | gender:race), data = d, subset = y %in% c("lo", "mid"))))
  expect_identical(mh$model$engine, "lme4")
  expect_identical(mh$model$family$family, "binomial")
  expect_false(is.null(mh$model_adjusted))
  expect_identical(mh$model_adjusted$engine, "lme4")

  # compare_maihda_groups() must agree and fit every group (no engine/family clash).
  cm <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ x + (1 | gender:race), data = d, group = "grp",
                          subset = y %in% c("lo", "mid"), min_group_n = 5)))
  expect_identical(attr(cm, "engine"), "lme4")
  expect_true(all(cm$status == "ok"))
  expect_true(all(is.finite(cm$vpc)))

  # Regression guard: with the full (unsubset) 3-level ordered outcome the ordinal
  # engine is still selected.
  skip_if_not_installed("ordinal")
  mh_ord <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | gender:race), data = d)))
  expect_identical(mh_ord$model$engine, "ordinal")
})
