make_binary_upset_model <- function(seed = 1, n = 1500) {
  set.seed(seed)
  d <- data.frame(
    parent  = rbinom(n, 1, 0.40),
    chronic = rbinom(n, 1, 0.25),
    low_edu = rbinom(n, 1, 0.35),
    rural   = rbinom(n, 1, 0.30)
  )
  d$outcome <- -1 + 0.6 * d$low_edu + 0.4 * d$chronic +
    0.7 * d$low_edu * d$chronic + rnorm(n)
  s <- make_strata(d, vars = c("parent", "chronic", "low_edu", "rural"))
  fit_maihda(outcome ~ (1 | stratum), data = s$data, engine = "lme4")
}

test_that("type = 'upset' returns a patchwork of binary-dimension strata", {
  model <- make_binary_upset_model()
  p <- plot(model, type = "upset")
  expect_s3_class(p, "patchwork")
})

test_that("upset is offered by the dispatcher and honours n_strata", {
  model <- make_binary_upset_model()
  # No error for a capped view, and it still composes.
  expect_s3_class(plot(model, type = "upset", n_strata = 6), "patchwork")
})

test_that("upset renders multi-level factor dimensions (category matrix)", {
  set.seed(7)
  n <- 1800
  d <- data.frame(
    parent = rbinom(n, 1, 0.4),
    edu3   = factor(sample(c("low", "mid", "high"), n, replace = TRUE)),
    rural  = rbinom(n, 1, 0.3)
  )
  d$outcome <- rnorm(n)
  s <- make_strata(d, vars = c("parent", "edu3", "rural"))
  model <- fit_maihda(outcome ~ (1 | stratum), data = s$data, engine = "lme4")

  # No longer an error: the factor expands to one matrix row per level.
  expect_s3_class(plot(model, type = "upset"), "patchwork")
})

test_that("upset quantity toggles between predicted values and random effects", {
  model <- make_binary_upset_model()
  expect_s3_class(plot(model, type = "upset", quantity = "predicted"), "patchwork")
  expect_s3_class(plot(model, type = "upset", quantity = "interaction"), "patchwork")
  # An unknown quantity is rejected by match.arg.
  expect_error(plot(model, type = "upset", quantity = "nonsense"))
})

test_that("dimension helpers classify indicators vs multi-level factors", {
  expect_true(maihda_dim_is_indicator(c(0L, 1L, 1L, 0L)))
  expect_true(maihda_dim_is_indicator(c(TRUE, FALSE, TRUE)))
  expect_false(maihda_dim_is_indicator(factor(c("female", "male"))))
  expect_false(maihda_dim_is_indicator(c(1L, 2L, 3L)))

  expect_equal(maihda_dim_levels(c(0L, 1L, 1L)), c(0, 1))
  expect_equal(maihda_dim_levels(factor(c("b", "a", "c"),
                                        levels = c("a", "b", "c"))),
               c("a", "b", "c"))
})

test_that("upset errors clearly when the model carries no per-dimension table", {
  # A model built straight from a numeric `stratum` column has no strata_info.
  set.seed(11)
  d <- data.frame(
    stratum = rep(1:8, each = 20),
    outcome = rnorm(160)
  )
  model <- fit_maihda(outcome ~ (1 | stratum), data = d, engine = "lme4")
  expect_error(plot(model, type = "upset"), "make_strata")
})

test_that("the predicted-strata refactor still carries n through prep", {
  model <- make_binary_upset_model()
  prep <- maihda_prepare_predicted_strata(
    model, summary(model), n_strata = 50,
    highlight = NULL, only_flagged = FALSE, select = "order")
  expect_false(isTRUE(prep$no_flagged))
  expect_true(all(c("predicted", "lower", "upper", "n") %in% names(prep$stratum_est)))
  expect_true(is.numeric(prep$fixed_reference))
})

test_that("plot.maihda_analysis dispatches type = 'upset' and forwards rope", {
  set.seed(3)
  n <- 1500
  d <- data.frame(
    a = rbinom(n, 1, 0.45),
    b = rbinom(n, 1, 0.30),
    c = factor(sample(c("lo", "mid", "hi"), n, replace = TRUE))
  )
  d$y <- rbinom(n, 1, plogis(-0.8 + 0.5 * d$a + 0.4 * d$b))
  a <- maihda(y ~ a + b + c + (1 | a:b:c), data = d,
              family = "binomial", interactions = "BH")
  expect_s3_class(plot(a, type = "upset"), "patchwork")
  # ROPE highlighting forwards through the analysis-level method (the screen gets a
  # `decision` column, resolved from the adjusted model).
  expect_s3_class(
    plot(a, type = "upset", highlight_interactions = TRUE,
         highlight_by = "rope", rope = 0.3),
    "patchwork")
})

test_that("maihda_upset_size scales rows with dimensions and cols with n_strata", {
  model <- make_binary_upset_model()           # 4 binary dims -> 4 rows
  n_total <- nrow(model$strata_info)

  sz <- maihda_upset_size(model, n_strata = 50)
  expect_true(is.numeric(sz$width) && is.numeric(sz$height))
  expect_identical(sz$rows, 4L)
  expect_identical(sz$cols, as.integer(min(n_total, 50)))

  # Capping the columns narrows the recommended width (n_total > 6 here).
  sz_small <- maihda_upset_size(model, n_strata = 6)
  expect_identical(sz_small$cols, 6L)
  expect_lt(sz_small$width, sz$width)
})

test_that("maihda_upset_size counts one row per factor level", {
  set.seed(8)
  n <- 1500
  d <- data.frame(
    sex = factor(sample(c("f", "m"), n, TRUE)),        # 2-level factor -> 2 rows
    edu = factor(sample(c("lo", "mid", "hi"), n, TRUE)), # 3-level factor -> 3 rows
    chronic = rbinom(n, 1, 0.3)                          # indicator       -> 1 row
  )
  d$y <- rnorm(n)
  s <- make_strata(d, vars = c("sex", "edu", "chronic"))
  model <- fit_maihda(y ~ (1 | stratum), data = s$data, engine = "lme4")
  expect_identical(maihda_upset_size(model)$rows, 6L)
})

test_that("maihda_upset_size errors without a make_strata table", {
  d <- data.frame(stratum = rep(1:8, each = 20), outcome = rnorm(160))
  model <- fit_maihda(outcome ~ (1 | stratum), data = d, engine = "lme4")
  expect_error(maihda_upset_size(model), "make_strata")
})
