test_that("compare_maihda_groups VPC matches a manual per-subset fit (shared strata)", {
  set.seed(3001)
  n <- 480
  d <- data.frame(
    country = rep(c("A", "B", "C"), length.out = n),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  strata_key <- interaction(d$gender, d$race, drop = TRUE)
  # Country-specific between-stratum signal so VPCs differ across groups
  u_a <- rnorm(nlevels(strata_key), sd = 1.2)[strata_key]
  u_b <- rnorm(nlevels(strata_key), sd = 0.3)[strata_key]
  base_u <- ifelse(d$country == "A", u_a, u_b)
  d$y <- 1 + 0.4 * d$age + base_u + rnorm(n, sd = 0.5)

  cmp <- compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "country")

  expect_s3_class(cmp, "maihda_group_comparison")
  expect_setequal(cmp$group, c("A", "B", "C"))
  expect_true(all(cmp$status == "ok"))

  # Manual anchor: build shared strata once, subset, re-attach attrs, fit.
  strata <- make_strata(d, vars = c("gender", "race"))
  dd <- strata$data
  attr_names <- c("strata_info", "strata_vars", "strata_sep", "strata_autobin_info")

  for (g in c("A", "B", "C")) {
    sub <- dd[dd$country == g, , drop = FALSE]
    for (a in attr_names) attr(sub, a) <- attr(dd, a)
    manual <- fit_maihda(y ~ age + (1 | stratum), data = sub)
    manual_vpc <- summary(manual)$vpc$estimate
    expect_equal(cmp$vpc[cmp$group == g], manual_vpc, tolerance = 1e-8)
  }
})

test_that("compare_maihda_groups reports variance components and strata counts", {
  set.seed(3002)
  n <- 300
  d <- data.frame(
    region = rep(c("North", "South"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  strata_key <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 2 + 0.3 * d$age + rnorm(nlevels(strata_key), sd = 0.7)[strata_key] +
    rnorm(n, sd = 0.4)

  cmp <- compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "region")

  expect_true(all(c("group", "n", "n_strata", "vpc", "var_between",
                    "var_residual", "status") %in% names(cmp)))
  # Full crossing of 2x2 is present in both regions
  expect_true(all(cmp$n_strata == 4))
  # VPC equals between / (between + residual) for a single random intercept model
  manual_vpc <- cmp$var_between / (cmp$var_between + cmp$var_residual)
  expect_equal(cmp$vpc, manual_vpc, tolerance = 1e-8)
})

test_that("compare_maihda_groups skips groups below min_group_n with a warning", {
  set.seed(3003)
  big <- data.frame(
    grp = "big",
    gender = sample(c("F", "M"), 200, replace = TRUE),
    race = sample(c("X", "Y"), 200, replace = TRUE),
    age = rnorm(200)
  )
  small <- data.frame(
    grp = "small",
    gender = sample(c("F", "M"), 8, replace = TRUE),
    race = sample(c("X", "Y"), 8, replace = TRUE),
    age = rnorm(8)
  )
  d <- rbind(big, small)
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.6)[sk] + rnorm(nrow(d), sd = 0.4)

  expect_warning(
    cmp <- compare_maihda_groups(y ~ age + (1 | gender:race), data = d,
                                 group = "grp", min_group_n = 30),
    "min_group_n"
  )
  expect_true(is.na(cmp$vpc[cmp$group == "small"]))
  expect_match(cmp$status[cmp$group == "small"], "skipped")
  expect_true(cmp$status[cmp$group == "big"] == "ok")
})

test_that("compare_maihda_groups reports VPC 0 (not an error) for a singular group", {
  set.seed(3004)
  n <- 240
  d <- data.frame(
    grp = rep(c("signal", "flat"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  signal_u <- rnorm(nlevels(sk), sd = 1.0)[sk]
  # "flat" group has no between-stratum variance -> singular fit -> VPC 0
  d$y <- ifelse(
    d$grp == "signal",
    1 + 0.3 * d$age + signal_u + rnorm(n, sd = 0.4),
    1 + 0.3 * d$age + rnorm(n, sd = 0.4)
  )

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "grp")
  )
  expect_true(all(is.finite(cmp$vpc)))
  expect_equal(cmp$vpc[cmp$group == "flat"], 0, tolerance = 1e-6)
  expect_true(cmp$vpc[cmp$group == "signal"] > cmp$vpc[cmp$group == "flat"])
})

test_that("compare_maihda_groups reports the analytic n after NA handling", {
  set.seed(3011)
  n <- 160
  d <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    x = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$x + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)
  # Drop 30 covariate values in group A: 80 raw rows -> 50 analytic rows.
  d$x[which(d$grp == "A")[1:30]] <- NA_real_

  cmp <- compare_maihda_groups(y ~ x + (1 | gender:race), data = d, group = "grp")
  expect_equal(cmp$n[cmp$group == "A"], 50L)   # analytic, not the raw 80
  expect_equal(cmp$n[cmp$group == "B"], 80L)
})

test_that("compare_maihda_groups bootstrap returns ordered per-group CIs", {
  set.seed(3005)
  n <- 300
  d <- data.frame(
    country = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.8)[sk] + rnorm(n, sd = 0.4)

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "country",
                          bootstrap = TRUE, n_boot = 25)
  )
  expect_true(all(c("ci_lower", "ci_upper") %in% names(cmp)))
  ok <- cmp$status == "ok"
  expect_true(all(cmp$ci_lower[ok] <= cmp$vpc[ok] + 1e-8))
  expect_true(all(cmp$ci_upper[ok] >= cmp$vpc[ok] - 1e-8))
})

test_that("compare_maihda_groups works with a precomputed stratum column", {
  set.seed(3006)
  n <- 260
  d <- data.frame(
    country = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)

  strata <- make_strata(d, vars = c("gender", "race"))
  cmp <- compare_maihda_groups(y ~ age + (1 | stratum), data = strata$data,
                               group = "country")
  expect_s3_class(cmp, "maihda_group_comparison")
  expect_true(all(cmp$status == "ok"))
})

test_that("compare_maihda_groups validates inputs", {
  d <- data.frame(country = rep(c("A", "B"), 5), gender = "F", race = "X",
                  age = rnorm(10), y = rnorm(10))

  expect_error(compare_maihda_groups("y ~ age + (1 | gender:race)", d, "country"),
               "must be a formula")
  expect_error(compare_maihda_groups(y ~ age + (1 | gender:race), "nope", "country"),
               "must be a data frame")
  expect_error(compare_maihda_groups(y ~ age + (1 | gender:race), d, "missing_col"),
               "Group variable not found")
  expect_error(
    compare_maihda_groups(y ~ age + (1 | gender:race), d, "country",
                          bootstrap = TRUE, engine = "brms"),
    "only supported for the lme4 engine", fixed = TRUE
  )
})

test_that("compare_maihda_groups ignores missing group values", {
  set.seed(3008)
  n <- 300
  d <- data.frame(
    country = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  d$country[c(1, 2, 3)] <- NA  # missing group labels must not become a group
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)

  cmp <- compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "country")
  expect_setequal(cmp$group, c("A", "B"))
  expect_false(any(is.na(cmp$group)))
})

test_that("compare_maihda_groups rejects shorthand formula when a stratum column exists", {
  set.seed(3009)
  n <- 200
  d <- data.frame(
    country = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)
  strata <- make_strata(d, vars = c("gender", "race"))

  expect_error(
    compare_maihda_groups(y ~ age + (1 | gender:race), data = strata$data,
                          group = "country"),
    "already has a 'stratum' column", fixed = TRUE
  )
})

test_that("compare_maihda_groups captures other random-effect variance consistently", {
  set.seed(3010)
  d <- expand.grid(
    country = c("A", "B"), gender = c("F", "M"), race = c("X", "Y"),
    site = factor(1:4), rep = 1:6
  )
  d$age <- rnorm(nrow(d))
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age +
    rnorm(nlevels(sk), sd = 0.7)[sk] +
    rnorm(4, sd = 0.9)[d$site] +
    rnorm(nrow(d), sd = 0.3)
  strata <- make_strata(d, vars = c("gender", "race"))

  cmp <- compare_maihda_groups(y ~ age + (1 | site) + (1 | stratum),
                               data = strata$data, group = "country")
  expect_true("var_other" %in% names(cmp))
  ok <- cmp$status == "ok"
  expect_true(all(cmp$var_other[ok] > 0))
  # VPC must use the full denominator: between / (between + other + residual)
  manual <- cmp$var_between / (cmp$var_between + cmp$var_other + cmp$var_residual)
  expect_equal(cmp$vpc[ok], manual[ok], tolerance = 1e-8)
  expect_s3_class(plot(cmp, type = "components"), "ggplot")
})

test_that("plot() on a maihda_group_comparison returns ggplot objects for both types", {
  set.seed(3007)
  n <- 300
  d <- data.frame(
    country = rep(c("A", "B", "C"), length.out = n),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)

  cmp <- compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "country")

  p_vpc <- plot(cmp, type = "vpc")
  p_comp <- plot(cmp, type = "components")
  expect_s3_class(p_vpc, "ggplot")
  expect_s3_class(p_comp, "ggplot")

  # the method's class guard (reachable via a direct call)
  expect_error(plot.maihda_group_comparison(mtcars), "maihda_group_comparison")

  # deprecated alias still works but warns
  expect_warning(plot_group_comparison(cmp, type = "vpc"), "deprecated")
})
