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
  expected <- mean(1 / pmax(stats::fitted(model$model), .Machine$double.eps))
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
  null_model <- fit_maihda(y ~ 1 + (1 | stratum), complete_d)
  x_model <- fit_maihda(y ~ x + (1 | stratum), complete_d)
  z_model <- fit_maihda(y ~ x + z + (1 | stratum), complete_d)

  expect_equal(out$Variance[1], MAIHDA:::extract_between_variance(null_model), tolerance = 1e-8)
  expect_equal(out$Variance[2], MAIHDA:::extract_between_variance(x_model), tolerance = 1e-8)
  expect_equal(out$Variance[3], MAIHDA:::extract_between_variance(z_model), tolerance = 1e-8)
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
