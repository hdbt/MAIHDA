make_table_data <- function(seed = 7001, n = 400) {
  set.seed(seed)
  d <- data.frame(
    country = rep(c("A", "B", "C"), length.out = n),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y", "Z"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.8)[sk] + rnorm(n, sd = 0.5)
  d$bin <- rbinom(n, 1, plogis(-0.2 + 0.4 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk]))
  d
}

test_that("maihda_table() on a two-model analysis assembles the canonical results table", {
  d <- make_table_data()
  a <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race), data = d))

  tab <- maihda_table(a)
  expect_s3_class(tab, "maihda_table")
  expect_s3_class(tab$models, "data.frame")
  expect_identical(tab$model_keys, c("null", "adjusted"))

  # The canonical statistics are present as rows.
  expect_true(all(c("Intercept", "Between-stratum variance", "Between-stratum SD",
                    "VPC/ICC", "PCV (null -> adjusted)") %in% tab$models$statistic))

  # The wide schema: an estimate + interval column per model.
  expect_true(all(c("null", "null_lower", "null_upper",
                    "adjusted", "adjusted_lower", "adjusted_upper") %in% names(tab$models)))

  # Values match the analysis the table is built from (no recomputation).
  vpc_row <- tab$models[tab$models$statistic == "VPC/ICC", ]
  expect_equal(vpc_row$null, a$summary$vpc$estimate, tolerance = 1e-10)

  pcv_row <- tab$models[tab$models$statistic == "PCV (null -> adjusted)", ]
  expect_equal(pcv_row$adjusted, a$pcv$pvc, tolerance = 1e-10)
  expect_true(is.na(pcv_row$null))   # PCV is a null -> adjusted quantity

  # SD is sqrt of the between-stratum variance.
  bv <- tab$models[tab$models$statistic == "Between-stratum variance", "null"]
  sdv <- tab$models[tab$models$statistic == "Between-stratum SD", "null"]
  expect_equal(sdv, sqrt(bv), tolerance = 1e-10)
})

test_that("maihda_table() ranks every stratum by predicted value, descending", {
  d <- make_table_data(7002)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))

  tab <- maihda_table(a)
  expect_s3_class(tab$strata, "data.frame")
  # One row per stratum, ranked.
  expect_equal(nrow(tab$strata), nrow(a$summary$stratum_estimates))
  expect_identical(tab$strata$rank, seq_len(nrow(tab$strata)))
  expect_true(all(diff(tab$strata$predicted) <= 1e-12))   # non-increasing
  expect_true(all(c("rank", "stratum", "label", "n", "predicted",
                    "random_effect") %in% names(tab$strata)))
  expect_equal(sum(tab$strata$n), nrow(a$model$data))
})

test_that("maihda_table() adds AUC and MOR rows for a binary outcome", {
  d <- make_table_data(7003)
  a <- suppressWarnings(suppressMessages(
    maihda(bin ~ age + gender + race + (1 | gender:race), data = d, family = "binomial")))

  tab <- maihda_table(a)
  expect_true(all(c("AUC", "MOR") %in% tab$models$statistic))

  auc_row <- tab$models[tab$models$statistic == "AUC", ]
  expect_equal(auc_row$null, a$summary$discriminatory_accuracy$auc, tolerance = 1e-10)
  # Adjusted-model AUC comes from the adjusted summary.
  expect_equal(auc_row$adjusted, a$summary_adjusted$discriminatory_accuracy$auc,
               tolerance = 1e-10)
})

test_that("maihda_table() accepts a single fitted model (no PCV)", {
  d <- make_table_data(7004)
  m <- suppressMessages(fit_maihda(y ~ age + (1 | gender:race), data = d))

  tab <- maihda_table(m)
  expect_s3_class(tab, "maihda_table")
  expect_identical(tab$model_keys, "estimate")
  expect_false("PCV (null -> adjusted)" %in% tab$models$statistic)
  expect_true("VPC/ICC" %in% tab$models$statistic)
  expect_equal(tab$models[tab$models$statistic == "VPC/ICC", "estimate"],
               summary(m)$vpc$estimate, tolerance = 1e-10)
  expect_s3_class(tab$strata, "data.frame")
})

test_that("maihda_table() reports additive/interaction shares for crossed-dimensions", {
  d <- make_table_data(7005)
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ age + (1 | gender:race), data = d, decomposition = "crossed-dimensions")))

  tab <- maihda_table(a)
  expect_identical(tab$model_keys, "estimate")
  expect_true(all(c("Additive share", "Interaction share") %in% tab$models$statistic))
  expect_false("PCV (null -> adjusted)" %in% tab$models$statistic)
  expect_equal(tab$models[tab$models$statistic == "Additive share", "estimate"],
               a$decomposition$additive_share, tolerance = 1e-10)
})

test_that("maihda_table() carries the VPC/PCV intervals when bootstrapped", {
  d <- make_table_data(7006)
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ age + (1 | gender:race), data = d, bootstrap = TRUE, n_boot = 20)))

  tab <- maihda_table(a)
  vpc_row <- tab$models[tab$models$statistic == "VPC/ICC", ]
  expect_false(is.na(vpc_row$null_lower))
  expect_false(is.na(vpc_row$null_upper))

  pcv_row <- tab$models[tab$models$statistic == "PCV (null -> adjusted)", ]
  expect_false(is.na(pcv_row$adjusted_lower))
})

test_that("maihda_table(which = 'adjusted') ranks strata by the adjusted model", {
  d <- make_table_data(7007)
  a <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race), data = d))

  tab_null <- maihda_table(a, which = "null")
  tab_adj <- maihda_table(a, which = "adjusted")
  expect_identical(tab_null$ranked_by, "null")
  expect_identical(tab_adj$ranked_by, "adjusted")

  # The adjusted-model strata predictions match plot(type = "predicted") on the
  # adjusted model -- i.e. they come from a different fit than the null ranking.
  pred_adj <- maihda_stratum_predictions_lme4(a$model_adjusted, a$summary_adjusted,
                                              scale = "response")
  idx <- match(tab_adj$strata$stratum, pred_adj$stratum)
  expect_equal(sort(tab_adj$strata$predicted), sort(pred_adj$predicted_row[idx]),
               tolerance = 1e-8)
})

test_that("print.maihda_table() runs and shows both tables", {
  d <- make_table_data(7008)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))
  tab <- maihda_table(a, n_strata = 2)

  out <- capture.output(print(tab))
  expect_true(any(grepl("MAIHDA Results Table", out)))
  expect_true(any(grepl("Model results:", out)))
  expect_true(any(grepl("Strata ranked by", out)))
})

test_that("maihda_table() errors on an unsupported input", {
  expect_error(maihda_table(list(a = 1)), "maihda_analysis|maihda_model")
})
