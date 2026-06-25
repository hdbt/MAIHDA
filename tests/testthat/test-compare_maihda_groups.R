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

test_that("compare_maihda_groups warns and names groups with a singular fit", {
  # "flat" carries the identical response pattern in every stratum -> zero
  # between-stratum variance -> a deterministically singular fit; "signal" has
  # real between-stratum variance and must NOT be named in the warning.
  flat <- data.frame(
    grp = "flat",
    gender = rep(c("F", "M"), each = 20),
    race   = rep(rep(c("X", "Y"), each = 10), 2),
    y = rep(c(-2, -1, 0, 1, 2), times = 8)
  )
  set.seed(3105)
  signal <- data.frame(
    grp = "signal",
    gender = sample(c("F", "M"), 200, replace = TRUE),
    race = sample(c("X", "Y"), 200, replace = TRUE)
  )
  sk <- interaction(signal$gender, signal$race, drop = TRUE)
  signal$y <- rnorm(nlevels(sk), sd = 1.0)[sk] + rnorm(200, sd = 0.4)
  d <- rbind(flat, signal)

  w <- testthat::capture_warnings(
    cmp <- compare_maihda_groups(y ~ 1 + (1 | gender:race), data = d, group = "grp")
  )
  fit_warn <- w[grepl("fit problems", w, fixed = TRUE)]
  expect_length(fit_warn, 1)
  expect_match(fit_warn, "singular fit: flat")
  # The non-singular group must not be named in the diagnostics warning.
  expect_false(grepl("signal", fit_warn))
  # The fit still completed, so status stays "ok"; the warning is the signal.
  expect_equal(cmp$status[cmp$group == "flat"], "ok")
  expect_equal(cmp$status[cmp$group == "signal"], "ok")
})

test_that("compare_maihda_groups skips a group with one analytic stratum (two raw)", {
  set.seed(2302)
  N <- 200
  d <- data.frame(
    grp = rep(c("A", "B"), each = N / 2),
    gender = sample(c("F", "M"), N, replace = TRUE),
    ses = sample(c("lo", "hi"), N, replace = TRUE)
  )
  sk <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- rnorm(nlevels(sk), sd = 0.6)[sk] + rnorm(N, sd = 0.4)
  # Group B has two raw strata (F.lo, F.hi), but every F.hi outcome is missing, so
  # only one stratum survives into the analytic frame. The pre-fit guard must skip
  # it as VPC-undefined rather than letting lme4 fail with "grouping factors must
  # have > 1 sampled level".
  bF <- d$grp == "B"
  d$gender[bF] <- "F"
  d$y[bF & d$ses == "hi"] <- NA_real_

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          min_group_n = 5)
  )
  expect_equal(cmp$n_strata[cmp$group == "B"], 1L)
  expect_match(cmp$status[cmp$group == "B"], "skipped")
  expect_match(cmp$status[cmp$group == "B"], "VPC undefined")
  expect_true(is.na(cmp$vpc[cmp$group == "B"]))
})

test_that("compare_maihda_groups accepts a data-column weights argument", {
  set.seed(4201)
  n <- 320
  d <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    ses = sample(c("lo", "hi"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)
  d$wcol <- runif(n, 0.5, 1.5)

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          min_group_n = 5, weights = wcol)
  )
  expect_true(all(cmp$status == "ok"))
})

test_that("compare_maihda_groups accepts a family function as well as a string/object", {
  set.seed(4202)
  n <- 320
  d <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    ses = sample(c("lo", "hi"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)

  run <- function(fam) {
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          min_group_n = 5, family = fam)
  }
  # As documented ("As in fit_maihda"), a name, a family object, and a bare family
  # function are all accepted. The bare function previously crashed when recording
  # the family attribute ("object of type 'closure' is not subsettable").
  cmp_str  <- run("gaussian")
  cmp_obj  <- run(stats::gaussian())
  cmp_fun  <- run(stats::gaussian)            # regression: must not error

  expect_equal(attr(cmp_str, "family"), "gaussian")
  expect_equal(attr(cmp_obj, "family"), "gaussian")
  expect_equal(attr(cmp_fun, "family"), "gaussian")
  # The three forms describe the same model, so the VPCs must coincide.
  expect_equal(cmp_fun$vpc, cmp_str$vpc, tolerance = 1e-10)
  expect_equal(cmp_obj$vpc, cmp_str$vpc, tolerance = 1e-10)
})

test_that("compare_maihda_groups applies an external weights vector per group", {
  set.seed(4202)
  n <- 320
  d <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    ses = sample(c("lo", "hi"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)
  w_ext <- runif(n, 0.5, 1.5)   # external vector, NOT a column of d

  # External weights used to fail every group with a length mismatch.
  cmp_ext <- suppressWarnings(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          min_group_n = 5, weights = w_ext)
  )
  expect_true(all(cmp_ext$status == "ok"))

  # The external vector must be applied exactly like the same weights as a column.
  d$wcol <- w_ext
  cmp_col <- suppressWarnings(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          min_group_n = 5, weights = wcol)
  )
  expect_equal(cmp_ext$vpc, cmp_col$vpc, tolerance = 1e-8)
  expect_equal(cmp_ext$var_residual, cmp_col$var_residual, tolerance = 1e-8)
})

test_that("compare_maihda_groups slices an external subset to each group's rows", {
  set.seed(4203)
  n <- 240
  d <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    ses = sample(c("lo", "hi"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)

  # Non-uniform external subset: keep 80 of group A (rows 1:80) and 30 of group B
  # (rows 121:150). A naive recycle of the full vector onto each group would fit
  # the wrong rows for group B.
  keep_ext <- rep(FALSE, n)
  keep_ext[1:80] <- TRUE
  keep_ext[121:150] <- TRUE

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          min_group_n = 5, subset = keep_ext)
  )
  expect_equal(cmp$n[cmp$group == "A"], 80L)
  expect_equal(cmp$n[cmp$group == "B"], 30L)
})

test_that("compare_maihda_groups warns when groups have different populated strata", {
  set.seed(2116)
  N <- 400
  g <- data.frame(
    grp = rep(c("A", "B"), each = N / 2),
    gender = sample(c("F", "M"), N, replace = TRUE),
    ses = sample(c("lo", "hi"), N, replace = TRUE)
  )
  # Group B contains only gender == "F", so it populates 2 of the 4 shared strata.
  g$gender[g$grp == "B"] <- "F"
  sk <- interaction(g$gender, g$ses, drop = TRUE)
  g$y <- rnorm(nlevels(sk), sd = 0.6)[sk] + rnorm(N, sd = 0.4)

  w <- testthat::capture_warnings(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = g, group = "grp",
                          min_group_n = 5)
  )
  expect_true(any(grepl("different populated strata", w)))
})

test_that("compare_maihda_groups skips a group whose analytic sample is below min_group_n", {
  set.seed(3111)
  N <- 200
  d <- data.frame(
    grp    = rep(c("ok", "tiny"), each = N / 2),
    gender = sample(c("F", "M"), N, replace = TRUE),
    ses    = sample(c("lo", "hi"), N, replace = TRUE),
    y      = rnorm(N)
  )
  # 'tiny' has 100 raw rows but only 8 non-missing outcomes -> 8 analytic rows.
  tiny_rows <- which(d$grp == "tiny")
  d$y[tiny_rows[seq_len(length(tiny_rows) - 8)]] <- NA_real_

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          min_group_n = 30)
  )
  # Guard is on the analytic sample, so 'tiny' is skipped despite 100 raw rows.
  expect_match(cmp$status[cmp$group == "tiny"], "skipped")
  expect_true(is.na(cmp$vpc[cmp$group == "tiny"]))
  expect_equal(cmp$n[cmp$group == "tiny"], 8L)
  expect_equal(cmp$status[cmp$group == "ok"], "ok")
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
  # Percentile bootstrap intervals are well-ordered and finite, but need not contain
  # the per-group point VPC (skewed/small-sample bootstrap distributions), so assert
  # ordering and finiteness rather than containment.
  expect_true(all(is.finite(cmp$ci_lower[ok])) && all(is.finite(cmp$ci_upper[ok])))
  expect_true(all(cmp$ci_lower[ok] <= cmp$ci_upper[ok]))
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

test_that("compare_maihda_groups rejects function-call grouping terms", {
  set.seed(6)
  d <- data.frame(
    country = rep(c("A", "B"), each = 60),
    g = sample(c("a", "b"), 120, replace = TRUE),
    r = sample(c("x", "y"), 120, replace = TRUE),
    age = rnorm(120)
  )
  d$y <- rnorm(120)
  expect_error(
    compare_maihda_groups(y ~ age + (1 | interaction(g, r)), data = d, group = "country"),
    "make_strata", fixed = TRUE
  )
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
  p_bv <- plot(cmp, type = "between_variance")
  expect_s3_class(p_vpc, "ggplot")
  expect_s3_class(p_comp, "ggplot")
  expect_s3_class(p_bv, "ggplot")
  # between_variance plots the absolute var_between, not the ratio
  expect_identical(p_bv$labels$y, "Between-stratum variance")

  # the method's class guard (reachable via a direct call)
  expect_error(plot.maihda_group_comparison(mtcars), "maihda_group_comparison")

  # deprecated alias still works but warns, including for the new type
  expect_warning(plot_group_comparison(cmp, type = "vpc"), "deprecated")
  expect_warning(plot_group_comparison(cmp, type = "between_variance"), "deprecated")
})

test_that("plot() notes omitted groups in the caption instead of dropping silently", {
  set.seed(3201)
  # Two large groups (estimable) and one tiny group skipped below min_group_n, so the
  # plot must still flag that 'tiny' was dropped rather than omit it without trace.
  big <- data.frame(
    grp = rep(c("A", "B"), each = 200),
    gender = sample(c("F", "M"), 400, replace = TRUE),
    race = sample(c("X", "Y"), 400, replace = TRUE),
    age = rnorm(400)
  )
  tiny <- data.frame(
    grp = "tiny", gender = sample(c("F", "M"), 10, replace = TRUE),
    race = sample(c("X", "Y"), 10, replace = TRUE), age = rnorm(10)
  )
  d <- rbind(big, tiny)
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(nrow(d), sd = 0.4)

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "grp",
                          min_group_n = 30)
  )

  for (ty in c("vpc", "components", "between_variance")) {
    cap <- plot(cmp, type = ty)$labels$caption
    expect_match(cap, "omitted")
    expect_match(cap, "tiny")
  }
})

test_that("compare_maihda_groups reports a per-group PCV with >= 2 dimensions", {
  set.seed(3301)
  n <- 480
  d <- data.frame(
    country = rep(c("A", "B", "C"), length.out = n),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.4 * d$age + rnorm(nlevels(sk), sd = 1.0)[sk] + rnorm(n, sd = 0.5)

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "country")
  )
  expect_true(all(c("pcv", "var_between_adjusted", "var_between_adjusted_ml") %in%
                    names(cmp)))
  ok <- cmp$status == "ok"
  expect_true(all(is.finite(cmp$pcv[ok])))
  expect_true(all(is.finite(cmp$var_between_adjusted[ok])))
  expect_true(all(is.finite(cmp$var_between_adjusted_ml[ok])))
  # The per-group PCV equals the proportional change in between-stratum variance from
  # the null to the adjusted model (a negative PCV is possible -- suppression -- so
  # we check the identity, not a sign/monotonicity). var_between_adjusted is the
  # derived coherence value, so this identity holds against it exactly.
  expect_equal(cmp$pcv[ok],
               (cmp$var_between[ok] - cmp$var_between_adjusted[ok]) / cmp$var_between[ok],
               tolerance = 1e-8)

  # var_between_adjusted_ml is the adjusted fit's ACTUAL between-stratum variance, so
  # it must reproduce the per-group adjusted model fitted on the same sample/strata
  # (the ML-refit var_model2 from calculate_pvc()). It is nonzero here and -- because
  # of the REML-vs-ML gap in the null variance -- not identical to the derived
  # var_between_adjusted coherence value.
  d_ok <- droplevels(d[d$country == "A", ])
  null_a <- fit_maihda(y ~ age + (1 | gender:race), data = d_ok)
  adj_a  <- fit_maihda(y ~ age + gender + race + (1 | gender:race), data = d_ok)
  ml_a   <- calculate_pvc(null_a, adj_a)$var_model2
  expect_equal(cmp$var_between_adjusted_ml[cmp$group == "A"], ml_a, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(
    cmp$var_between_adjusted[cmp$group == "A"],
    cmp$var_between_adjusted_ml[cmp$group == "A"]
  )))

  expect_s3_class(plot(cmp, type = "pcv"), "ggplot")
})

test_that("compare_maihda_groups strips supplied dimension main effects before the per-group PCV", {
  # Regression: passing a fully-specified ADJUSTED formula (the dimensions' additive
  # main effects written out) directly to compare_maihda_groups() must give the SAME
  # per-group decomposition as the bare null shorthand. Previously the supplied
  # formula was fitted as the per-group BASE model, then maihda_adjusted_formula()
  # re-added the same main effects, so the "adjusted" model equalled the base: PCV
  # collapsed to 0/NA and the stratum random effect went singular.
  set.seed(3303)
  n <- 600
  d <- data.frame(
    country = rep(c("A", "B", "C"), length.out = n),
    gender = sample(c("F", "M"), n, replace = TRUE),
    ses = sample(c("lo", "mid", "hi"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- 1 + 0.4 * d$age + rnorm(nlevels(sk), sd = 1.0)[sk] + rnorm(n, sd = 0.5)

  # Bare shorthand: the supplied fit IS the null, the documented base form.
  cmp_null <- suppressWarnings(
    compare_maihda_groups(y ~ age + (1 | gender:ses), data = d, group = "country")
  )
  # Fully-specified adjusted formula: the dimensions' main effects are present and
  # must be stripped so the per-group base stays the null.
  cmp_adj <- suppressWarnings(
    compare_maihda_groups(y ~ age + gender + ses + (1 | gender:ses), data = d,
                          group = "country")
  )

  ok <- cmp_null$status == "ok" & cmp_adj$status == "ok"
  expect_true(any(ok))
  # Identical null base => identical VPC, between-stratum variance, and PCV.
  expect_equal(cmp_adj$vpc[ok], cmp_null$vpc[ok], tolerance = 1e-6)
  expect_equal(cmp_adj$var_between[ok], cmp_null$var_between[ok], tolerance = 1e-6)
  expect_equal(cmp_adj$pcv[ok], cmp_null$pcv[ok], tolerance = 1e-6)
  # Finite and not collapsed to ~0 (the symptom of the unstripped re-added terms).
  expect_true(all(is.finite(cmp_adj$pcv[ok])))
  expect_false(all(abs(cmp_adj$pcv[ok]) < 1e-6))
})

test_that("compare_maihda_groups omits the PCV columns with a single dimension", {
  set.seed(3302)
  n <- 300
  d <- data.frame(
    country = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE)
  )
  d$y <- rnorm(2, sd = 0.6)[as.integer(factor(d$gender))] + rnorm(n, sd = 0.4)

  cmp <- suppressWarnings(
    compare_maihda_groups(y ~ 1 + (1 | gender), data = d, group = "country")
  )
  expect_false("pcv" %in% names(cmp))
  expect_false("var_between_adjusted" %in% names(cmp))
  expect_error(plot(cmp, type = "pcv"), "fewer than two")
})

test_that("plot.maihda_analysis dispatches group_between_variance", {
  set.seed(3202)
  n <- 240
  d <- data.frame(
    country = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    ses = sample(c("lo", "hi"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$ses, drop = TRUE)
  d$math <- 1 + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)

  a <- suppressWarnings(suppressMessages(
    maihda(math ~ 1 + (1 | gender:ses), data = d, group = "country")
  ))
  expect_s3_class(plot(a, type = "group_between_variance"), "ggplot")

  # Without a group, the group types error rather than silently no-op.
  a_nogroup <- suppressWarnings(suppressMessages(maihda(math ~ 1 + (1 | gender:ses), data = d)))
  expect_error(plot(a_nogroup, type = "group_between_variance"),
               "No group comparison")
})

test_that("supplied dimension main effects strip even for backtick-containing names", {
  # When the formula already lists the stratum dimensions as fixed effects, the
  # two-model decomposition strips them before deriving the null/adjusted models.
  # Detection used maihda_quote_name() (proper backtick escaping) but removal used a
  # bare sprintf("`%s`", .), which mis-parses a dimension name that itself contains a
  # backtick -- so as.formula() would error on the update. Both paths must quote the
  # same way. A backtick in a column name is rare but legal (check.names = FALSE).
  set.seed(3101)
  n <- 600
  d <- data.frame(
    grp  = sample(c("P", "Q"), n, replace = TRUE),
    dimA = sample(c("F", "M"), n, replace = TRUE),
    weird = sample(c("lo", "hi"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  names(d)[names(d) == "weird"] <- "a`b"        # backtick-containing legal name
  key <- interaction(d$dimA, d[["a`b"]], drop = TRUE)
  d$y <- 1 + 0.3 * (d$dimA == "M") + rnorm(nlevels(key), sd = 0.6)[key] + rnorm(n, sd = 0.5)

  q <- maihda_quote_name("a`b")
  f <- stats::as.formula(paste("y ~ dimA +", q, "+ (1 | dimA:", q, ")"))

  # Pre-fix this errored in the term-stripping update(); now it completes.
  cmp <- suppressWarnings(suppressMessages(
    compare_maihda_groups(f, data = d, group = "grp", engine = "lme4")))
  expect_s3_class(cmp, "maihda_group_comparison")
  expect_setequal(cmp$group, c("P", "Q"))
})

test_that("compare_maihda_groups rejects a mixed-sign numeric subset like base R", {
  set.seed(9091)
  d <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    gender = rep(c("F", "M"), 20),
    ses = rep(c("lo", "hi"), each = 2, length.out = 40),
    y = rnorm(40)
  )
  # c(-1, 2) mixes positive and negative subscripts: base R and a bare fit error on
  # it, so the comparison must too -- not silently normalise it to an empty (n = 0)
  # mask that reports every group as skipped.
  expect_error(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                          subset = c(-1, 2)),
    "subset")
})

test_that("subsetting a maihda_group_comparison keeps its print metadata", {
  set.seed(7711)
  n <- 300
  d <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    gender = sample(c("F", "M"), n, replace = TRUE),
    ses = sample(c("lo", "hi"), n, replace = TRUE)
  )
  key <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- rnorm(nlevels(key), sd = 0.7)[key] + rnorm(n, sd = 0.4)

  cmp <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp")))

  sub <- cmp[, c("group", "n", "status")]
  expect_s3_class(sub, "maihda_group_comparison")
  # Metadata attributes survive the column subset (regression: `[.data.frame` kept
  # the class but dropped them, blanking the printed header).
  expect_identical(attr(sub, "group_var"), attr(cmp, "group_var"))
  expect_identical(attr(sub, "engine"), attr(cmp, "engine"))
  expect_identical(attr(sub, "family"), attr(cmp, "family"))

  # Printing the column subset still shows the populated header line.
  out <- paste(utils::capture.output(print(sub)), collapse = "\n")
  expect_match(out, "Group variable: grp")

  # A single-column drop returns a bare vector, not a classed object.
  expect_false(is.data.frame(cmp[, "group"]))
})
