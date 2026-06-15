test_that("fit_maihda preserves explicit additional random effects", {
  set.seed(1001)
  d <- expand.grid(
    stratum = factor(seq_len(5)),
    site = factor(seq_len(4)),
    rep = seq_len(8)
  )
  d$x <- rnorm(nrow(d))
  stratum_u <- rnorm(5, sd = 0.8)[d$stratum]
  site_u <- rnorm(4, sd = 1.2)[d$site]
  d$y <- 2 + 0.5 * d$x + stratum_u + site_u + rnorm(nrow(d), sd = 0.3)

  model <- fit_maihda(y ~ x + (1 | site) + (1 | stratum), data = d)

  expect_true(all(c("site", "stratum") %in% names(lme4::VarCorr(model$model))))
  expect_null(model$strata_info)

  summ <- summary(model)
  vc <- lme4::VarCorr(model$model)
  expected <- as.numeric(as.matrix(vc[["stratum"]])["(Intercept)", "(Intercept)"])
  observed <- summ$variance_components$variance[1]
  expect_equal(observed, expected, tolerance = 1e-8)

  site_var <- as.numeric(as.matrix(vc[["site"]])["(Intercept)", "(Intercept)"])
  resid_var <- attr(vc, "sc")^2
  expected_vpc <- expected / (expected + site_var + resid_var)
  expect_equal(summ$vpc$estimate, expected_vpc, tolerance = 1e-8)
  expect_equal(
    summ$variance_components$variance[
      summ$variance_components$component == "Other random effects"
    ],
    site_var,
    tolerance = 1e-8
  )
})

test_that("fit_maihda rejects ambiguous automatic strata formulas", {
  set.seed(1002)
  d <- data.frame(
    gender = sample(c("F", "M"), 100, replace = TRUE),
    race = sample(c("A", "B"), 100, replace = TRUE),
    site = sample(letters[1:4], 100, replace = TRUE),
    x = rnorm(100),
    y = rnorm(100)
  )

  expect_error(
    fit_maihda(y ~ x + (1 | gender:race) + (1 | site), data = d),
    "Automatic strata creation"
  )
})

test_that("automatic strata creation does not bin numeric fixed effects", {
  set.seed(1003)
  n <- 180
  d <- data.frame(
    age = runif(n, 18, 80),
    gender = sample(c("F", "M"), n, replace = TRUE)
  )
  strata <- interaction(cut(d$age, breaks = 3), d$gender, drop = TRUE)
  stratum_u <- rnorm(nlevels(strata), sd = 0.7)[as.integer(strata)]
  d$y <- 1 + 0.2 * d$age + ifelse(d$gender == "M", 0.5, 0) + stratum_u + rnorm(n, sd = 0.4)

  model <- fit_maihda(y ~ age + gender + (1 | age:gender), data = d)

  expect_true(is.numeric(model$data$age))
  expect_true("age" %in% names(lme4::fixef(model$model)))
  expect_false(any(grepl("^ageage_", names(lme4::fixef(model$model)))))
})

test_that("make_strata preserves original numeric grouping variables", {
  d <- data.frame(
    age = seq_len(60),
    gender = rep(c("F", "M"), 30),
    y = rnorm(60)
  )

  strata <- make_strata(d, vars = c("age", "gender"))

  expect_true(is.numeric(strata$data$age))
  expect_true(any(grepl("age_", strata$strata_info$label)))
})

test_that("make_strata does not collapse values that contain the display separator", {
  d <- data.frame(
    a = rep(c("A \u00d7 B", "A"), each = 20),
    b = rep(c("C", "B \u00d7 C"), each = 20),
    x = rnorm(40)
  )

  strata <- make_strata(d, vars = c("a", "b"))

  expect_equal(nrow(strata$strata_info), 2)
  expect_equal(length(unique(stats::na.omit(strata$data$stratum))), 2)
  expect_equal(as.character(strata$strata_info$a), c("A \u00d7 B", "A"))
  expect_equal(as.character(strata$strata_info$b), c("C", "B \u00d7 C"))
})

test_that("maihda_match_strata_rows is collision-safe and keeps row-by-row semantics", {
  lookup <- data.frame(
    a = c("A \u00d7 B", "A", "X"),
    b = c("C", "B \u00d7 C", "Y"),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    a = c("A", "A \u00d7 B", "X", "A", NA, "Z"),
    b = c("B \u00d7 C", "C", "Y", "C", "C", "Y"),
    stringsAsFactors = FALSE
  )
  # First matching lookup row per data row; NA when nothing matches or a value is
  # missing/absent from the lookup. Row 4 (A, C) must NOT collide with the display
  # labels of rows 1/2, and row 6 (Z, ...) is an unknown value.
  expect_equal(
    MAIHDA:::maihda_match_strata_rows(data, lookup, c("a", "b")),
    c(2L, 1L, 3L, NA, NA, NA)
  )
})

test_that("calculate_pvc errors for different analytic samples", {
  set.seed(1004)
  d <- data.frame(
    stratum = factor(rep(seq_len(10), each = 10)),
    x = rnorm(100),
    z = rnorm(100)
  )
  d$z[seq_len(20)] <- NA_real_
  d$y <- 1 + d$x + rnorm(10, sd = 0.8)[d$stratum] + rnorm(100, sd = 0.3)

  model1 <- fit_maihda(y ~ x + (1 | stratum), data = d)
  model2 <- fit_maihda(y ~ x + z + (1 | stratum), data = d)

  expect_error(calculate_pvc(model1, model2), "same analytic sample")
})

test_that("calculate_pvc rejects different outcomes and families", {
  set.seed(1018)
  d <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    x = rnorm(200)
  )
  stratum_u <- rnorm(8, sd = 0.7)[d$stratum]
  d$y_gaussian <- 1 + 0.4 * d$x + stratum_u + rnorm(200, sd = 0.3)
  d$y_count <- rpois(200, lambda = exp(0.2 + 0.3 * d$x + stratum_u))

  gaussian_model <- fit_maihda(y_gaussian ~ x + (1 | stratum), data = d)
  poisson_model <- fit_maihda(y_count ~ x + (1 | stratum), data = d, family = "poisson")
  count_gaussian_model <- fit_maihda(y_count ~ x + (1 | stratum), data = d)

  expect_error(
    calculate_pvc(gaussian_model, poisson_model),
    "same outcome",
    fixed = TRUE
  )
  expect_error(
    calculate_pvc(count_gaussian_model, poisson_model),
    "same model family and link",
    fixed = TRUE
  )
})

test_that("poisson summaries use fitted mean based residual variance", {
  set.seed(1005)
  d <- data.frame(
    stratum = factor(rep(seq_len(12), each = 25)),
    x = rnorm(300)
  )
  stratum_u <- rnorm(12, sd = 0.5)[d$stratum]
  d$y <- rpois(nrow(d), lambda = exp(1 + 0.25 * d$x + stratum_u))

  model <- fit_maihda(y ~ x + (1 | stratum), data = d, family = "poisson")
  summ <- summary(model)

  observed <- summ$variance_components$variance[
    summ$variance_components$component == "Within-stratum (residual)"
  ]
  expected <- mean(log1p(1 / pmax(stats::fitted(model$model), .Machine$double.eps)))
  expect_equal(observed, expected, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(observed, 1)))
})

test_that("binomial predict_maihda defaults to response scale and supports link scale", {
  set.seed(1006)
  d <- data.frame(
    stratum = factor(rep(seq_len(10), each = 25)),
    x = rnorm(250)
  )
  stratum_u <- rnorm(10, sd = 0.8)[d$stratum]
  d$y <- rbinom(nrow(d), 1, plogis(-0.3 + 0.7 * d$x + stratum_u))

  model <- fit_maihda(y ~ x + (1 | stratum), data = d, family = "binomial")
  pred_response <- predict_maihda(model)
  pred_link <- predict_maihda(model, type = "link")

  expect_true(all(pred_response >= 0 & pred_response <= 1))
  expect_equal(pred_response, plogis(pred_link), tolerance = 1e-8)

  p <- plot(model, type = "predicted")
  expect_true(all(p$data$predicted >= 0 & p$data$predicted <= 1, na.rm = TRUE))
})

test_that("fit_maihda recodes two-level responses for glmer", {
  set.seed(1010)
  n <- 240
  d <- data.frame(
    stratum = factor(rep(seq_len(12), each = n / 12)),
    x = rnorm(n)
  )
  stratum_u <- rnorm(12, sd = 0.5)[d$stratum]
  y01 <- rbinom(n, 1, plogis(-0.2 + 0.5 * d$x + stratum_u))

  d$y <- ifelse(y01 == 1, 2, 1)
  expect_warning(
    numeric_model <- fit_maihda(y ~ x + (1 | stratum), data = d),
    "Automatically switching to family = 'binomial'",
    fixed = TRUE
  )
  expect_s3_class(numeric_model, "maihda_model")
  expect_equal(sort(unique(numeric_model$data$y)), c(0L, 1L))

  d$y <- ifelse(y01 == 1, "case", "control")
  expect_warning(
    character_model <- fit_maihda(y ~ x + (1 | stratum), data = d),
    "Automatically switching to family = 'binomial'",
    fixed = TRUE
  )
  expect_s3_class(character_model, "maihda_model")
  expect_equal(sort(unique(character_model$data$y)), c(0L, 1L))
})

test_that("stepwise_pcv quotes non-syntactic variable names", {
  set.seed(1011)
  d <- data.frame(
    check.names = FALSE,
    "health outcome" = rnorm(80),
    "age years" = rnorm(80),
    "race group" = rep(letters[1:8], each = 10)
  )
  strata <- make_strata(d, "race group")

  out <- stepwise_pcv(strata$data, "health outcome", "age years")

  expect_s3_class(out, "maihda_stepwise")
  expect_equal(out$Added_Variable[2], "age years")
})

test_that("stepwise_pcv adds no discriminatory-accuracy columns for a non-binary outcome", {
  # Backward-compat guarantee: a gaussian stepwise table is exactly the historical
  # six columns -- the AUC/MOR trajectory appears only for a binary outcome.
  set.seed(2024)
  d <- data.frame(
    stratum = factor(rep(seq_len(6), each = 20)),
    x = rnorm(120)
  )
  d$y <- 1 + d$x + rnorm(6, sd = 0.8)[d$stratum] + rnorm(120, sd = 0.3)

  out <- stepwise_pcv(d, "y", "x")
  expect_identical(
    names(out),
    c("Step", "Model", "Added_Variable", "Variance", "Step_PCV", "Total_PCV")
  )
})

test_that("stepwise_pcv uses one complete analytic sample across steps", {
  set.seed(1016)
  d <- data.frame(
    stratum = factor(rep(seq_len(6), each = 20)),
    x = rnorm(120),
    z = rnorm(120)
  )
  d$z[seq_len(30)] <- NA_real_
  d$y <- 1 + d$x + rnorm(6, sd = 0.8)[d$stratum] + rnorm(120, sd = 0.3)

  out <- stepwise_pcv(d, "y", c("x", "z"))
  complete_d <- d[stats::complete.cases(d[, c("y", "stratum", "x", "z")]), ]
  # stepwise_pcv refits REML lmer fits with ML before the cross-step variance
  # comparison (see calculate_pvc()), so the reference models must be refit likewise.
  ml <- function(m) MAIHDA:::extract_between_variance(MAIHDA:::maihda_pcv_refit_ml(m))
  null_model <- fit_maihda(y ~ 1 + (1 | stratum), complete_d)
  x_model <- fit_maihda(y ~ x + (1 | stratum), complete_d)
  z_model <- fit_maihda(y ~ x + z + (1 | stratum), complete_d)

  expect_equal(out$Variance[1], ml(null_model), tolerance = 1e-8)
  expect_equal(out$Variance[2], ml(x_model), tolerance = 1e-8)
  expect_equal(out$Variance[3], ml(z_model), tolerance = 1e-8)
})

test_that("stepwise_pcv errors when no complete analytic sample remains", {
  d <- data.frame(
    stratum = factor(c("a", "b")),
    x = c(NA_real_, 1),
    y = c(1, NA_real_)
  )

  expect_error(
    stepwise_pcv(d, "y", "x"),
    "No complete cases remain",
    fixed = TRUE
  )
})

test_that("stepwise_pcv auto-detects a binary outcome when family is the default", {
  set.seed(910)
  base <- data.frame(
    g = sample(c("a", "b"), 240, replace = TRUE),
    r = sample(c("x", "y"), 240, replace = TRUE),
    age = rnorm(240)
  )
  strata <- make_strata(base, vars = c("g", "r"))
  d <- strata$data
  d$y <- rbinom(240, 1, plogis(-0.2 + 0.3 * d$age))

  # Default family must not silently fit a Gaussian PCV on a binary outcome.
  expect_warning(
    out <- stepwise_pcv(d, "y", "age"),
    "binary", ignore.case = TRUE
  )
  expect_s3_class(out, "maihda_stepwise")

  # Explicit gaussian is respected: no binary auto-switch (other warnings, e.g.
  # singular fits, are irrelevant here).
  w <- testthat::capture_warnings(stepwise_pcv(d, "y", "age", family = "gaussian"))
  expect_false(any(grepl("binary", w, ignore.case = TRUE)))
})

test_that("predict_maihda(type = 'strata') respects newdata", {
  set.seed(2301)
  d <- data.frame(
    gender = rep(c("F", "M"), each = 60),
    race = rep(c("A", "B"), times = 60),
    age = rnorm(120)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.4 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(120, sd = 0.3)
  model <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  # With no newdata, every training stratum is returned.
  all_strata <- predict_maihda(model, type = "strata")
  expect_gt(nrow(all_strata), 1L)

  # newdata from a single stratum returns only that stratum, not all of them.
  nd <- data.frame(gender = "F", race = "A", age = rnorm(3))
  got <- predict_maihda(model, newdata = nd, type = "strata")
  expect_equal(nrow(got), 1L)

  # A stratum the model never saw is an error, as for type = "individual".
  nd_unknown <- nd
  nd_unknown$stratum <- "9999"
  expect_error(
    predict_maihda(model, newdata = nd_unknown, type = "strata"),
    "not present in the fitted model", fixed = TRUE
  )

  # A stratum column that is present but entirely missing names no training
  # stratum, so the result is empty -- not silently every training stratum.
  nd_all_na <- nd
  nd_all_na$stratum <- NA_character_
  empty <- predict_maihda(model, newdata = nd_all_na, type = "strata")
  expect_equal(nrow(empty), 0L)
  expect_equal(names(empty), names(all_strata))
})

test_that("predict_maihda rebuilds automatic strata for raw newdata", {
  set.seed(1008)
  d <- data.frame(
    gender = rep(c("F", "M"), each = 40),
    race = rep(c("A", "B"), times = 40),
    age = rnorm(80)
  )
  stratum_key <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.5 * d$age + rnorm(nlevels(stratum_key), sd = 0.8)[stratum_key] + rnorm(80, sd = 0.3)

  model <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  raw_newdata <- d[1:20, c("gender", "race", "age")]
  explicit_newdata <- raw_newdata
  explicit_newdata$stratum <- model$original_data$stratum[1:20]

  expect_equal(
    as.numeric(predict_maihda(model, newdata = raw_newdata)),
    as.numeric(predict_maihda(model, newdata = explicit_newdata)),
    tolerance = 1e-8
  )
})

test_that("predict_maihda rebuilds manual make_strata strata for raw newdata", {
  set.seed(1012)
  d <- data.frame(
    gender = rep(c("F", "M"), each = 50),
    race = rep(c("A", "B"), times = 50),
    age = rnorm(100)
  )
  stratum_key <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.5 * d$age + rnorm(nlevels(stratum_key), sd = 0.8)[stratum_key] + rnorm(100, sd = 0.3)

  strata <- make_strata(d, c("gender", "race"))
  model <- fit_maihda(y ~ age + (1 | stratum), data = strata$data)

  raw_newdata <- d[1:20, c("gender", "race", "age")]
  explicit_newdata <- raw_newdata
  explicit_newdata$stratum <- strata$data$stratum[1:20]

  expect_equal(model$strata_vars, c("gender", "race"))
  expect_equal(
    as.numeric(predict_maihda(model, newdata = raw_newdata)),
    as.numeric(predict_maihda(model, newdata = explicit_newdata)),
    tolerance = 1e-8
  )
})

test_that("predict_maihda rebuilds strata when grouping values contain the separator", {
  set.seed(1019)
  d <- data.frame(
    a = rep(c("A \u00d7 B", "A"), each = 40),
    b = rep(c("C", "B \u00d7 C"), each = 40),
    x = rnorm(80)
  )
  stratum_key <- interaction(d$a, d$b, drop = TRUE)
  d$y <- 1 + 0.3 * d$x + rnorm(nlevels(stratum_key), sd = 0.8)[stratum_key] +
    rnorm(80, sd = 0.2)

  strata <- make_strata(d, c("a", "b"))
  model <- fit_maihda(y ~ x + (1 | stratum), data = strata$data)
  raw_newdata <- d[1:20, c("a", "b", "x")]
  explicit_newdata <- raw_newdata
  explicit_newdata$stratum <- strata$data$stratum[1:20]

  expect_equal(
    as.numeric(predict_maihda(model, newdata = raw_newdata)),
    as.numeric(predict_maihda(model, newdata = explicit_newdata)),
    tolerance = 1e-8
  )
})

test_that("fit_maihda refreshes strata counts after model NA handling", {
  set.seed(1013)
  d <- data.frame(
    group = rep(LETTERS[1:4], each = 15),
    x = rnorm(60)
  )
  d$x[which(d$group == "A")[1:5]] <- NA_real_
  group_effect <- rnorm(4, sd = 0.8)[match(d$group, LETTERS[1:4])]
  d$y <- 1 + 0.5 * d$x + group_effect + rnorm(60, sd = 0.2)

  strata <- make_strata(d, "group")
  model <- fit_maihda(y ~ x + (1 | stratum), data = strata$data)
  analytic_counts <- table(as.character(model$data$stratum), useNA = "no")
  expected_counts <- as.integer(analytic_counts[
    match(as.character(model$strata_info$stratum), names(analytic_counts))
  ])
  expected_counts[is.na(expected_counts)] <- 0L

  expect_equal(model$strata_info$n, expected_counts)
  expect_equal(model$strata_info$n[model$strata_info$group == "A"], 10L)
})

test_that("summary rejects random-slope models for VPC accounting", {
  set.seed(1014)
  d <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    age = rnorm(200)
  )
  intercept_u <- rnorm(8, sd = 0.7)[d$stratum]
  slope_u <- rnorm(8, sd = 0.4)[d$stratum]
  d$y <- 1 + 0.3 * d$age + intercept_u + slope_u * d$age + rnorm(200, sd = 0.3)

  model <- fit_maihda(y ~ age + (age | stratum), data = d)

  expect_error(
    summary(model),
    "intercept-only random effects",
    fixed = TRUE
  )
})

test_that("calculate_pvc rejects random-slope models for PVC accounting", {
  set.seed(1017)
  d <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    age = rnorm(200),
    poverty = rnorm(200)
  )
  intercept_u <- rnorm(8, sd = 0.7)[d$stratum]
  slope_u <- rnorm(8, sd = 0.4)[d$stratum]
  d$y <- 1 + 0.3 * d$age + 0.2 * d$poverty +
    intercept_u + slope_u * d$age + rnorm(200, sd = 0.3)

  model1 <- fit_maihda(y ~ age + (age | stratum), data = d)
  model2 <- fit_maihda(y ~ age + poverty + (age | stratum), data = d)

  expect_error(
    calculate_pvc(model1, model2),
    "intercept-only random effects",
    fixed = TRUE
  )
})

test_that("predict_maihda reuses training bins for numeric automatic strata", {
  set.seed(1009)
  n <- 120
  d <- data.frame(
    age = runif(n, 18, 80),
    gender = sample(c("F", "M"), n, replace = TRUE)
  )
  strata <- interaction(cut(d$age, breaks = 3), d$gender, drop = TRUE)
  d$y <- 2 + 0.1 * d$age + rnorm(nlevels(strata), sd = 0.6)[strata] + rnorm(n, sd = 0.4)

  model <- fit_maihda(y ~ age + gender + (1 | age:gender), data = d)

  raw_newdata <- d[1:30, c("age", "gender")]
  explicit_newdata <- raw_newdata
  explicit_newdata$stratum <- model$original_data$stratum[1:30]

  expect_equal(
    as.numeric(predict_maihda(model, newdata = raw_newdata)),
    as.numeric(predict_maihda(model, newdata = explicit_newdata)),
    tolerance = 1e-8
  )
})

test_that("predict_maihda errors clearly for out-of-range numeric auto-bins", {
  set.seed(1015)
  n <- 120
  d <- data.frame(
    age = runif(n, 18, 80),
    gender = sample(c("F", "M"), n, replace = TRUE)
  )
  strata <- interaction(cut(d$age, breaks = 3), d$gender, drop = TRUE)
  d$y <- 2 + 0.1 * d$age + rnorm(nlevels(strata), sd = 0.6)[strata] + rnorm(n, sd = 0.4)

  model <- fit_maihda(y ~ age + gender + (1 | age:gender), data = d)
  out_of_range <- data.frame(age = c(10, 100), gender = c("F", "M"))

  expect_error(
    predict_maihda(model, newdata = out_of_range),
    "outside the training auto-bin ranges",
    fixed = TRUE
  )
})

test_that("maihda() resolves caller-local weights/subset through the data mask", {
  # Reproduces the wrapper bug: from inside a function, weights/subset that name a
  # caller-local variable were lost because ... forwarding became ..1 promises.
  run <- function() {
    set.seed(2110)
    n <- 240
    d <- data.frame(
      gender = sample(c("F", "M"), n, replace = TRUE),
      race   = sample(c("X", "Y"), n, replace = TRUE),
      x = rnorm(n)
    )
    sk <- interaction(d$gender, d$race, drop = TRUE)
    d$y <- 1 + 0.4 * d$x + rnorm(4, sd = 0.6)[sk] + rnorm(n, sd = 0.3)
    w_local <- runif(n, 0.5, 1.5)       # caller-local, NOT a column of d
    keep_local <- d$x > 0               # caller-local logical

    # The forwarding must reach BOTH the null and the adjusted fit maihda() builds.
    a <- suppressMessages(maihda(y ~ x + (1 | gender:race), data = d, weights = w_local))
    expect_equal(unname(stats::weights(a$model$model, type = "prior")), w_local)
    expect_equal(unname(stats::weights(a$model_adjusted$model, type = "prior")), w_local)

    b <- suppressMessages(maihda(y ~ x + (1 | gender:race), data = d, subset = keep_local))
    expect_equal(as.integer(stats::nobs(b$model$model)), sum(keep_local))
    expect_equal(as.integer(stats::nobs(b$model_adjusted$model)), sum(keep_local))
  }
  run()
})

test_that("compare_maihda_groups resolves a weights column for every group", {
  run <- function() {
    set.seed(2111)
    n <- 240
    d <- data.frame(
      grp = rep(c("A", "B"), each = n / 2),
      gender = sample(c("F", "M"), n, replace = TRUE),
      ses = sample(c("lo", "hi"), n, replace = TRUE)
    )
    sk <- interaction(d$gender, d$ses, drop = TRUE)
    d$y <- rnorm(nlevels(sk), sd = 0.6)[sk] + rnorm(n, sd = 0.4)
    d$wcol <- runif(n, 0.5, 1.5)
    cmp <- suppressWarnings(
      compare_maihda_groups(y ~ 1 + (1 | gender:ses), data = d, group = "grp",
                            min_group_n = 5, weights = wcol)
    )
    expect_true(all(cmp$status == "ok"))
  }
  run()
})

test_that("fit_maihda recodes a binary outcome without breaking a label subset", {
  set.seed(2120)
  n <- 200
  d <- data.frame(stratum = factor(rep(seq_len(8), each = 25)), x = rnorm(n))
  d$y <- ifelse(rbinom(n, 1, plogis(0.3 * d$x)) == 1, "yes", "no")

  # subset references the ORIGINAL labels; recoding to 0/1 must not break it.
  expect_warning(
    m <- fit_maihda(y ~ x + (1 | stratum), data = d, subset = y %in% c("no", "yes")),
    "binary", ignore.case = TRUE
  )
  expect_equal(m$family$family, "binomial")
  expect_equal(as.integer(stats::nobs(m$model)), n)
})

test_that("weighted Gaussian VPC uses the mean conditional residual variance", {
  set.seed(2112)
  n <- 300
  d <- data.frame(stratum = factor(rep(seq_len(10), each = 30)), x = rnorm(n))
  d$y <- 1 + 0.4 * d$x + rnorm(10, sd = 0.8)[d$stratum] + rnorm(n, sd = 0.5)
  d$wt <- runif(n, 0.2, 3)

  m <- fit_maihda(y ~ x + (1 | stratum), data = d, weights = wt)
  s <- summary(m)
  vc <- lme4::VarCorr(m$model)
  sigma2 <- attr(vc, "sc")^2
  var_between <- as.numeric(as.matrix(vc[["stratum"]])["(Intercept)", "(Intercept)"])
  prior_w <- stats::weights(m$model, type = "prior")
  resid_expected <- sigma2 * mean(1 / prior_w)   # mean conditional residual variance

  observed_resid <- s$variance_components$variance[
    s$variance_components$component == "Within-stratum (residual)"
  ]
  expect_equal(observed_resid, resid_expected, tolerance = 1e-8)
  expect_equal(s$vpc$estimate, var_between / (var_between + resid_expected), tolerance = 1e-8)
  # ...and it differs from the naive sigma^2-only VPC the package used to report.
  expect_false(isTRUE(all.equal(s$vpc$estimate,
                                var_between / (var_between + sigma2))))
})

test_that("PVC/comparison flag models that differ only in prior weights", {
  set.seed(2401)
  d <- make_strata(maihda_sim_data, vars = c("gender", "race"))$data
  n <- nrow(d)
  w <- runif(n, 0.2, 5)
  m_unw  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = d)
  m_unw2 <- fit_maihda(health_outcome ~ age + (1 | stratum), data = d)
  m_w    <- fit_maihda(health_outcome ~ age + (1 | stratum), data = d, weights = w)
  m_unit <- fit_maihda(health_outcome ~ age + (1 | stratum), data = d, weights = rep(1, n))

  # Different weights, same rows/outcome/strata: hard error in PVC, warning in compare.
  expect_error(calculate_pvc(m_unw, m_w), "same prior weights", fixed = TRUE)
  expect_true(any(grepl("prior weights",
                        testthat::capture_warnings(compare_maihda(m_unw, m_w)))))

  # Unweighted and explicit unit weights are equivalent (no weight complaint).
  expect_s3_class(calculate_pvc(m_unw, m_unit), "pvc_result")
  expect_false(any(grepl("prior weights",
                         testthat::capture_warnings(compare_maihda(m_unw, m_unit)))))
  expect_s3_class(calculate_pvc(m_unw, m_unw2), "pvc_result")
})

test_that("binary recoding records and reports the level -> 0/1 mapping", {
  set.seed(2402)
  n <- 200
  d <- data.frame(stratum = factor(rep(seq_len(8), each = 25)), x = rnorm(n))
  d$y <- ifelse(rbinom(n, 1, plogis(0.3 * d$x)) == 1, "yes", "no")

  # Explicit family = "binomial" still surfaces the mapping (message + stored).
  expect_message(
    m_char <- fit_maihda(y ~ x + (1 | stratum), data = d, family = "binomial"),
    "recoded to 0/1"
  )
  expect_equal(m_char$response_recoding$level[m_char$response_recoding$value == 1L], "yes")

  # A factor whose levels are reversed maps the OPPOSITE way, and that is recorded.
  d$yf <- factor(d$y, levels = c("yes", "no"))
  m_fac <- suppressMessages(
    fit_maihda(yf ~ x + (1 | stratum), data = d, family = "binomial")
  )
  expect_equal(m_fac$response_recoding$level[m_fac$response_recoding$value == 1L], "no")

  # An already-0/1 numeric outcome is a no-op: no recoding message.
  d$y01 <- rbinom(n, 1, 0.5)
  msgs01 <- character(0)
  withCallingHandlers(
    suppressWarnings(fit_maihda(y01 ~ x + (1 | stratum), data = d, family = "binomial")),
    message = function(mm) {
      msgs01 <<- c(msgs01, conditionMessage(mm)); invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("recoded to 0/1", msgs01)))
})

test_that("stratum predictions honour prior weights; unweighted aggregation unchanged", {
  set.seed(2403)
  N <- 160
  dat <- data.frame(A = rep(c("a", "b"), each = N / 2),
                    B = rep(c("p", "q"), times = N / 2))
  dat$x <- rnorm(N)
  sk <- interaction(dat$A, dat$B, drop = TRUE)
  dat$y <- 2 + 1.5 * dat$x + rnorm(4, sd = 0.7)[sk] + rnorm(N, sd = 0.3)
  dat$wt <- exp(0.9 * dat$x)   # weights correlate with x within stratum

  mw <- fit_maihda(y ~ x + (1 | A:B), data = dat, weights = wt)
  ph <- MAIHDA:::maihda_stratum_predictions_lme4(mw, summary(mw), scale = "response")
  fit_incl <- predict(mw$model)
  prw <- stats::weights(mw$model, type = "prior")
  strat <- as.character(mw$data$stratum)
  wmean <- tapply(seq_along(fit_incl), strat,
                  function(i) stats::weighted.mean(fit_incl[i], prw[i]))
  helper <- stats::setNames(ph$predicted_row, as.character(ph$stratum))
  expect_equal(as.numeric(helper[names(wmean)]), as.numeric(wmean), tolerance = 1e-6)
  expect_true("w_sum" %in% names(ph))

  # Unweighted model: the per-stratum aggregation is the plain mean (unchanged).
  mu <- fit_maihda(y ~ x + (1 | A:B), data = dat)
  pu <- MAIHDA:::maihda_stratum_predictions_lme4(mu, summary(mu), scale = "response")
  umean <- tapply(predict(mu$model), as.character(mu$data$stratum), mean)
  helperu <- stats::setNames(pu$predicted_row, as.character(pu$stratum))
  expect_equal(as.numeric(helperu[names(umean)]), as.numeric(umean), tolerance = 1e-10)
  expect_equal(pu$w_sum, as.numeric(pu$n))
})

test_that("binary detection respects negative subset indices and NA weights", {
  set.seed(2113)
  n <- 150
  d <- data.frame(stratum = factor(rep(seq_len(6), each = 25)), x = rnorm(n))
  y <- rbinom(n, 1, 0.5)
  y[1:6] <- 2L                 # a spurious third value on the first six rows
  d$y <- y

  # Negative subset excludes exactly those rows -> analytic response is 0/1.
  m_sub <- suppressWarnings(fit_maihda(y ~ x + (1 | stratum), data = d, subset = -(1:6)))
  expect_equal(m_sub$family$family, "binomial")

  # NA weights drop the same rows (lme4 removes them), so detection must too.
  w <- rep(1, n)
  w[1:6] <- NA_real_
  m_w <- suppressWarnings(fit_maihda(y ~ x + (1 | stratum), data = d, weights = w))
  expect_equal(m_w$family$family, "binomial")
})

test_that("effect decomposition isolates the stratum random effect from other REs", {
  set.seed(2105)
  d <- expand.grid(stratum = factor(seq_len(5)), site = factor(seq_len(4)), rep = seq_len(12))
  d$x <- rnorm(nrow(d))
  d$y <- 2 + 0.5 * d$x + rnorm(5, sd = 0.8)[d$stratum] +
    rnorm(4, sd = 1.5)[d$site] + rnorm(nrow(d), sd = 0.3)

  model <- fit_maihda(y ~ x + (1 | site) + (1 | stratum), data = d)
  summ <- summary(model)
  p <- MAIHDA:::plot_effect_decomposition(model, summ)

  seg <- NULL
  for (ly in p$layers) {
    if (!is.null(ly$data) && "Component" %in% names(ly$data)) seg <- ly$data
  }
  expect_false(is.null(seg))
  inter <- seg[seg$Component == "Stratum random-effect component", ]
  inter_dev <- inter$y_end - inter$y_start
  # The intersectional component must be exactly the stratum random effects, NOT
  # total-minus-fixed (which would also absorb the (1 | site) variance).
  expect_equal(sort(inter_dev), sort(summ$stratum_estimates$random_effect),
               tolerance = 1e-6)
})

test_that("bootstrap VPC reports Monte Carlo error", {
  set.seed(2114)
  n <- 200
  d <- data.frame(stratum = factor(rep(seq_len(8), each = 25)), x = rnorm(n))
  d$y <- 1 + 0.3 * d$x + rnorm(8, sd = 0.7)[d$stratum] + rnorm(n, sd = 0.3)
  m <- fit_maihda(y ~ x + (1 | stratum), data = d)

  s <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 50))
  expect_true(is.finite(s$vpc$mc_se))
  expect_true(s$vpc$n_boot_ok >= 10 && s$vpc$n_boot_ok <= 50)
  expect_output(print(s), "Monte Carlo SE", fixed = TRUE)
})

test_that("VPC/ICC errors for a Gaussian model with a non-identity link", {
  set.seed(2001)
  d <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    x = rnorm(200)
  )
  stratum_u <- rnorm(8, sd = 0.3)[d$stratum]
  d$y <- exp(0.5 + 0.2 * d$x + stratum_u) * rlnorm(200, 0, 0.15)

  model <- suppressWarnings(
    fit_maihda(y ~ x + (1 | stratum), data = d, family = gaussian(link = "log"))
  )
  # The log link is honoured (routed through glmer), but the VPC mixes a
  # response-scale residual with a link-scale random effect, so it must error
  # rather than report an invalid partition.
  expect_equal(model$family$link, "log")
  expect_error(summary(model), "non-identity link", fixed = TRUE)
  expect_error(
    MAIHDA:::maihda_residual_variance_lme4(model$model),
    "non-identity link", fixed = TRUE
  )
})

test_that("fit_maihda forwards data-masked weights= and subset= to lme4", {
  set.seed(2004)
  d <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    x = rnorm(200)
  )
  d$y <- 1 + 0.4 * d$x + rnorm(8, sd = 0.6)[d$stratum] + rnorm(200, sd = 0.3)
  d$w <- runif(200, 0.5, 1.5)

  # Both reproduced a "..1 used in an incorrect context" error before the call was
  # built explicitly. weights= is a column of `data`; subset= is an expression over it.
  expect_s3_class(
    fit_maihda(y ~ x + (1 | stratum), data = d, weights = w),
    "maihda_model"
  )
  m_sub <- fit_maihda(y ~ x + (1 | stratum), data = d, subset = x > 0)
  expect_s3_class(m_sub, "maihda_model")
  expect_equal(as.integer(stats::nobs(m_sub$model)), sum(d$x > 0))
})

test_that("bootstrap n_boot must clear the 10-refit minimum up front", {
  set.seed(2006)
  d <- data.frame(
    stratum = factor(rep(seq_len(6), each = 20)),
    x = rnorm(120)
  )
  d$y <- 1 + 0.3 * d$x + rnorm(6, sd = 0.6)[d$stratum] + rnorm(120, sd = 0.3)
  model <- fit_maihda(y ~ x + (1 | stratum), data = d)

  expect_error(
    summary(model, bootstrap = TRUE, n_boot = 5),
    "n_boot' must be a single whole number >= 10", fixed = TRUE
  )
})

test_that("ternary additive component is invariant to row order", {
  set.seed(1007)
  d <- data.frame(
    stratum = factor(rep(seq_len(6), each = 20)),
    age = rep(seq(-2, 2, length.out = 20), 6)
  )
  stratum_u <- rnorm(6, sd = 0.5)[d$stratum]
  d$y <- 10 + 3 * d$age + stratum_u + rnorm(nrow(d), sd = 0.2)

  model <- fit_maihda(y ~ age + (1 | stratum), data = d)
  reordered <- model
  reordered$data <- model$data[order(model$data$stratum, -model$data$age), ]

  td1 <- compute_maihda_ternary_data(model, verbose = FALSE)
  td2 <- compute_maihda_ternary_data(reordered, verbose = FALSE)
  comp <- merge(
    as.data.frame(td1[, c("stratum", "additive_only")]),
    as.data.frame(td2[, c("stratum", "additive_only")]),
    by = "stratum",
    suffixes = c("_original", "_reordered")
  )

  expect_equal(comp$additive_only_original, comp$additive_only_reordered, tolerance = 1e-8)
})
