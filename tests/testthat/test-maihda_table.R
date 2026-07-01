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

test_that("maihda_table() reports the context share for a contextual fit", {
  d <- make_table_data(7009)
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ age + (1 | gender:race), data = d, context = "country")))

  tab <- maihda_table(a, n_strata = 50)   # 50 > #strata: exercise the full-list print
  expect_true("Context share (VPC)" %in% tab$models$statistic)
  ctx_row <- tab$models[tab$models$statistic == "Context share (VPC)", ]
  expect_equal(ctx_row$null, a$summary$context$vpc_context_total, tolerance = 1e-10)

  # The between-stratum variance comes from the contextual partition's stratum term.
  bv <- tab$models[tab$models$statistic == "Between-stratum variance", "null"]
  expect_equal(bv, a$summary$context$var_stratum, tolerance = 1e-10)

  out <- capture.output(print(tab))
  expect_true(any(grepl("Context: country", out, fixed = TRUE)))
})

test_that("maihda_table() notes the REML-variance / ML-PCV basis mismatch", {
  d <- make_table_data(7011)
  a <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race), data = d))

  tab <- maihda_table(a)
  # A Gaussian lme4 two-model fit reports REML variance rows but an ML-refitted PCV,
  # so a clarifying note is attached...
  expect_type(tab$models_note, "character")
  expect_match(tab$models_note, "maximum-likelihood")

  # ...and the note is warranted: the ML between-stratum variance the PCV is built
  # from genuinely differs from the displayed REML variance row.
  v_null_disp <- tab$models[tab$models$statistic == "Between-stratum variance", "null"]
  expect_false(isTRUE(all.equal(a$pcv$var_model1, v_null_disp)))

  # The note surfaces in print().
  out <- capture.output(print(tab))
  expect_true(any(grepl("Note:", out, fixed = TRUE)))
})

test_that("maihda_table() attaches no PCV note for single or crossed-dimensions fits", {
  d <- make_table_data(7012)

  m <- suppressMessages(fit_maihda(y ~ age + (1 | gender:race), data = d))
  expect_null(maihda_table(m)$models_note)

  a <- suppressWarnings(suppressMessages(
    maihda(y ~ age + (1 | gender:race), data = d, decomposition = "crossed-dimensions")))
  expect_null(maihda_table(a)$models_note)
})

test_that("maihda_extract_intercept handles every fixed-effects shape", {
  # brms-style matrix with an Intercept row + Estimate column.
  m <- matrix(c(1.5, 0.2, 1.1, 1.9), nrow = 1,
              dimnames = list("Intercept", c("Estimate", "Est.Error", "Q2.5", "Q97.5")))
  expect_equal(maihda_extract_intercept(list(fixed_effects = m)),
               c(1.5, NA_real_, NA_real_))

  # Matrix without an intercept row -> NA triple.
  m2 <- matrix(c(0.3), nrow = 1, dimnames = list("age", "Estimate"))
  expect_true(all(is.na(maihda_extract_intercept(list(fixed_effects = m2)))))

  # NULL fixed effects -> NA triple.
  expect_true(all(is.na(maihda_extract_intercept(list(fixed_effects = NULL)))))

  # Data frame without term/estimate columns (e.g. a degenerate summary) -> NA.
  expect_true(all(is.na(
    maihda_extract_intercept(list(fixed_effects = data.frame(x = 1))))))

  # Data frame with terms but no intercept (the ordinal/threshold case) -> NA.
  fe <- data.frame(term = c("age", "genderM"), estimate = c(0.3, -0.2))
  expect_true(all(is.na(maihda_extract_intercept(list(fixed_effects = fe)))))
})

test_that("maihda_build_results_table fills statistics absent from one model", {
  # The adjusted model lacks AUC the null has -> the missing cell becomes NA.
  ms <- list(
    null = list("VPC/ICC" = c(0.10, NA, NA), "AUC" = c(0.62, NA, NA)),
    adjusted = list("VPC/ICC" = c(0.05, NA, NA))
  )
  df <- maihda_build_results_table(ms, pcv = NULL)
  auc <- df[df$statistic == "AUC", ]
  expect_equal(auc$null, 0.62)
  expect_true(is.na(auc$adjusted))
})

# Mocking internal bindings needs the dev package loaded as a namespace (it is
# under R CMD check / load_all; under covr::file_coverage's bare sourcing the only
# MAIHDA namespace is the installed copy, which may predate these internals). Gate
# on a current internal binding actually being present in the namespace.
maihda_ns_loaded <- function() {
  isTRUE(tryCatch(
    exists("maihda_stratum_predictions_wemix", envir = asNamespace("MAIHDA"),
           inherits = FALSE),
    error = function(e) FALSE))
}

test_that("maihda_strata_ranking dispatches per engine and validates input", {
  skip_if_not(maihda_ns_loaded(), "MAIHDA namespace not loaded")
  fake_pred <- function(object, summary_obj, scale) {
    data.frame(stratum = c("a", "b"), predicted_row = c(2, 1),
               lower_row = c(1.5, 0.5), upper_row = c(2.5, 1.5),
               n = c(10L, 12L), stringsAsFactors = FALSE)
  }
  se <- data.frame(stratum = c("a", "b"), label = c("A", "B"),
                   random_effect = c(0.5, -0.5), lower_95 = c(0.1, -0.9),
                   upper_95 = c(0.9, -0.1), stringsAsFactors = FALSE)
  so <- list(stratum_estimates = se)

  testthat::local_mocked_bindings(
    maihda_stratum_predictions_brms = fake_pred,
    maihda_stratum_predictions_wemix = fake_pred,
    maihda_stratum_predictions_ordinal = fake_pred,
    .package = "MAIHDA"
  )
  for (eng in c("brms", "wemix", "ordinal")) {
    r <- maihda_strata_ranking(list(engine = eng), so, scale = "response")
    expect_equal(r$rank, c(1L, 2L))
    expect_equal(r$stratum[1], "a")   # predicted 2 ranks above 1
  }

  # Unknown engine -> the switch's stop().
  expect_error(maihda_strata_ranking(list(engine = "nope"), so),
               "Unsupported engine")
})

test_that("maihda_strata_ranking errors when there are no stratum estimates", {
  skip_if_not(maihda_ns_loaded(), "MAIHDA namespace not loaded")
  fake_pred <- function(object, summary_obj, scale) {
    data.frame(stratum = "a", predicted_row = 1, lower_row = 0, upper_row = 2, n = 1L)
  }
  testthat::local_mocked_bindings(maihda_stratum_predictions_lme4 = fake_pred,
                                  .package = "MAIHDA")
  expect_error(
    maihda_strata_ranking(list(engine = "lme4"), list(stratum_estimates = NULL)),
    "No stratum estimates")
})

test_that("maihda_table() degrades gracefully when strata ranking fails", {
  skip_if_not(maihda_ns_loaded(), "MAIHDA namespace not loaded")
  d <- make_table_data(7010)
  m <- suppressMessages(fit_maihda(y ~ age + (1 | gender:race), data = d))

  testthat::local_mocked_bindings(maihda_strata_ranking = function(...) stop("boom"),
                                  .package = "MAIHDA")
  tab <- maihda_table(m)
  expect_null(tab$strata)
  # The strata count still falls back to the summary's stratum estimates.
  expect_equal(tab$n_strata_total, nrow(summary(m)$stratum_estimates))
})
