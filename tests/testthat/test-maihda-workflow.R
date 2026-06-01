make_workflow_data <- function(seed = 4001, n = 360) {
  set.seed(seed)
  d <- data.frame(
    country = rep(c("A", "B", "C"), length.out = n),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)
  d$bin <- rbinom(n, 1, plogis(-0.2 + 0.4 * d$age + rnorm(nlevels(sk), sd = 0.6)[sk]))
  d
}

test_that("maihda() returns a consistent bundle; groups NULL without group", {
  d <- make_workflow_data()
  a <- maihda(y ~ age + (1 | gender:race), data = d)

  expect_s3_class(a, "maihda_analysis")
  expect_s3_class(a$model, "maihda_model")
  expect_s3_class(a$summary, "maihda_summary")
  expect_null(a$groups)

  # Overall VPC equals fitting + summarising directly
  m_direct <- fit_maihda(y ~ age + (1 | gender:race), data = d)
  expect_equal(a$summary$vpc$estimate, summary(m_direct)$vpc$estimate, tolerance = 1e-8)
})

test_that("maihda() with group attaches a comparison equal to compare_maihda_groups()", {
  d <- make_workflow_data(4002)
  a <- maihda(y ~ age + (1 | gender:race), data = d, group = "country")

  expect_s3_class(a$groups, "maihda_group_comparison")
  expect_setequal(a$groups$group, c("A", "B", "C"))

  cmp_direct <- compare_maihda_groups(y ~ age + (1 | gender:race), d, group = "country")
  expect_equal(a$groups$vpc, cmp_direct$vpc, tolerance = 1e-8)
  expect_equal(a$groups$var_between, cmp_direct$var_between, tolerance = 1e-8)
})

test_that("maihda() forwards comparison-only args without leaking them into lmer", {
  d <- make_workflow_data(4007)
  # min_group_n and shared_strata are compare_maihda_groups args, not lmer args;
  # they must not be passed through to the model fitter.
  expect_no_error(
    a <- maihda(y ~ age + (1 | gender:race), data = d, group = "country",
                min_group_n = 10, shared_strata = TRUE)
  )
  expect_s3_class(a$groups, "maihda_group_comparison")
})

test_that("maihda() auto-detects a binary outcome consistently for model and groups", {
  d <- make_workflow_data(4003)
  expect_warning(
    a <- maihda(bin ~ age + (1 | gender:race), data = d, group = "country"),
    "binary", ignore.case = TRUE
  )
  expect_equal(a$model$family$family, "binomial")
  # The group comparison used the same resolved family (no silent gaussian fallback)
  expect_equal(attr(a$groups, "family"), "binomial")
})

test_that("print and summary methods work for maihda_analysis", {
  d <- make_workflow_data(4004)
  a <- maihda(y ~ age + (1 | gender:race), data = d, group = "country")

  expect_output(print(a), "MAIHDA Analysis")
  expect_output(print(a), "VPC/ICC")
  expect_output(print(a), "Group comparison by 'country'")

  s <- summary(a)
  expect_s3_class(s, "maihda_summary")
  expect_s3_class(attr(s, "groups"), "maihda_group_comparison")
})

test_that("plot.maihda_analysis dispatches to model and group plots", {
  d <- make_workflow_data(4005)
  a <- maihda(y ~ age + (1 | gender:race), data = d, group = "country")

  expect_s3_class(plot(a, type = "vpc"), "ggplot")
  expect_s3_class(plot(a, type = "predicted"), "ggplot")
  expect_s3_class(plot(a, type = "group_vpc"), "ggplot")
  expect_s3_class(plot(a, type = "group_components"), "ggplot")

  # group plots require a group argument
  a_nogroup <- maihda(y ~ age + (1 | gender:race), data = d)
  expect_error(plot(a_nogroup, type = "group_vpc"), "group", ignore.case = TRUE)
})

test_that("compute_maihda_ternary_data is classed and plots via plot()", {
  skip_if_not_installed("ggtern")
  d <- make_workflow_data(4006)
  model <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  td <- compute_maihda_ternary_data(model, verbose = FALSE)
  expect_s3_class(td, "maihda_ternary")
  expect_s3_class(plot(td), "ggplot")

  expect_warning(plot_maihda_ternary(td), "deprecated")
})
