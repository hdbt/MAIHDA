test_that("Shiny app dependency gate includes ternary plotting dependency", {
  expect_true("ggtern" %in% MAIHDA:::maihda_app_required_packages())
})

test_that("maihda_app_fit_models switches a binary factor outcome to binomial", {
  set.seed(11)
  d <- data.frame(
    Obese = factor(sample(c("No", "Yes"), 240, replace = TRUE)),
    Gender = sample(c("F", "M"), 240, replace = TRUE),
    Race = sample(c("A", "B"), 240, replace = TRUE),
    Age = rnorm(240)
  )
  # App default family = "gaussian"; a Yes/No factor must be auto-switched to
  # binomial (gaussian would error: lmer requires a numeric response).
  expect_message(
    res <- MAIHDA:::maihda_app_fit_models(
      d, outcome_var = "Obese", grouping_vars = c("Gender", "Race"),
      additional_covars = "Age", family = "gaussian"
    ),
    "binomial"
  )
  expect_equal(res$model$family$family, "binomial")
  expect_equal(res$null_model$family$family, "binomial")
})

test_that("maihda_app_fit_models switches a numeric 0/1 outcome to binomial", {
  set.seed(2115)
  n <- 600
  d <- data.frame(
    Gender = sample(c("F", "M"), n, replace = TRUE),
    Race = sample(c("A", "B"), n, replace = TRUE)
  )
  sk <- interaction(d$Gender, d$Race, drop = TRUE)
  d$Obese <- rbinom(n, 1, plogis(rnorm(nlevels(sk), sd = 1)[sk]))  # numeric 0/1

  # A no-code user leaving the default gaussian on a numeric 0/1 outcome must not
  # silently get a linear probability model; the app mirrors the core API.
  expect_message(
    res <- MAIHDA:::maihda_app_fit_models(
      d, outcome_var = "Obese", grouping_vars = c("Gender", "Race"),
      family = "gaussian"
    ),
    "binomial"
  )
  expect_equal(res$null_model$family$family, "binomial")
})

test_that("Shiny app dependency gate leaves upload-only readers optional", {
  expect_false("haven" %in% MAIHDA:::maihda_app_required_packages())
})

test_that("Shiny PVC HUD display separates residual variance from unmasking", {
  positive <- MAIHDA:::maihda_app_pvc_display(35)
  expect_equal(positive$label, "Residual Strata Variance")
  expect_equal(positive$value, "65%")
  expect_equal(positive$status, "nonnegative")

  negative <- MAIHDA:::maihda_app_pvc_display(-12.5)
  expect_equal(negative$label, "Unmasked Variance")
  expect_equal(negative$value, "+12.5%")
  expect_equal(negative$remaining_value, "112.5%")
  expect_equal(negative$status, "negative")
})

test_that("Shiny ternary plotly helper includes boundary strata", {
  skip_if_not_installed("plotly")

  td <- data.frame(
    additive_prop = c(0, 1),
    interaction_prop = c(1, 0),
    uncertainty_prop = c(0, 0),
    label = c("A", "B"),
    n = c(10, 25)
  )

  p <- plotly::plotly_build(MAIHDA:::maihda_app_ternary_plotly(td))

  expect_equal(p$x$layout$ternary$aaxis$min, 0)
  expect_equal(p$x$layout$ternary$baxis$min, 0)
  expect_equal(p$x$layout$ternary$caxis$min, 0)
})

test_that("Shiny app fit helper builds the model objects used by the dashboard", {
  dat <- MAIHDA::maihda_sim_data[seq_len(150), ]

  res <- MAIHDA:::maihda_app_fit_models(
    dat = dat,
    outcome_var = "health_outcome",
    grouping_vars = c("gender", "race"),
    additional_covars = "age",
    family = "gaussian",
    use_boot = FALSE,
    n_boot = 10,
    autobin = TRUE,
    engine = "lme4"
  )

  expect_s3_class(res$null_model, "maihda_model")
  expect_s3_class(res$model, "maihda_model")
  expect_s3_class(res$pvc, "pvc_result")
  expect_s3_class(res$stepwise, "maihda_stepwise")
  expect_identical(row.names(res$null_model$data), row.names(res$model$data))
  expect_equal(res$model$strata_vars, c("gender", "race"))
  expect_equal(nrow(res$stepwise), 4)
})

test_that("maihda_app_fit_models degrades gracefully when baseline between-stratum variance is zero", {
  # Four strata that each hold the IDENTICAL multiset of (>2) outcome values: the
  # null model has exactly zero between-stratum variance, so calculate_pvc()
  # errors by design. The fit itself is valid, though, so the helper must still
  # return the models (for the dashboard to show VPC/summaries/plots) and flag the
  # PCV as unavailable rather than aborting the whole analysis.
  base <- data.frame(
    gender = rep(c("F", "M"), each = 6),
    race   = rep(rep(c("A", "B"), each = 3), 2),
    y      = rep(c(2, 4, 6), times = 4)
  )
  d <- base[rep(seq_len(nrow(base)), 10), ]

  res <- expect_no_error(
    suppressWarnings(suppressMessages(
      MAIHDA:::maihda_app_fit_models(
        d, outcome_var = "y", grouping_vars = c("gender", "race"),
        family = "gaussian"
      )
    ))
  )

  expect_s3_class(res$null_model, "maihda_model")
  expect_s3_class(res$model, "maihda_model")
  expect_s3_class(res$pvc, "pvc_result")
  # PCV flagged unavailable, with the underlying reason carried for the UI.
  expect_true(is.na(res$pvc$pvc))
  expect_false(isTRUE(res$pvc$available))
  expect_true(is.character(res$pvc$message) && nzchar(res$pvc$message))
  # Stepwise PCV already tolerates zero variance (returns NA), so it still builds.
  expect_s3_class(res$stepwise, "maihda_stepwise")
})

test_that("maihda_app_fit_models computes VPC/ICC bootstrap intervals when requested", {
  skip_on_cran()
  dat <- MAIHDA::maihda_sim_data[seq_len(150), ]

  res <- suppressWarnings(suppressMessages(
    MAIHDA:::maihda_app_fit_models(
      dat = dat, outcome_var = "health_outcome",
      grouping_vars = c("gender", "race"), additional_covars = "age",
      family = "gaussian", use_boot = TRUE, n_boot = 25, engine = "lme4"
    )
  ))

  # The worker returns finite VPC/ICC intervals for both the null and adjusted
  # models -- this is what the "Compute Bootstrap CIs" control advertises.
  expect_length(res$vpc_ci_null, 2)
  expect_length(res$vpc_ci_adjusted, 2)
  expect_true(all(is.finite(res$vpc_ci_null)))
  expect_lte(res$vpc_ci_null[1], res$vpc_ci_null[2])

  # Attaching an interval yields a summary the VPC interval helpers recognise.
  summ <- MAIHDA:::maihda_app_attach_vpc_ci(summary(res$null_model), res$vpc_ci_null)
  expect_true(MAIHDA:::maihda_vpc_has_interval(summ$vpc))
  expect_identical(summ$vpc$method, "bootstrap")
  expect_true(isTRUE(summ$vpc$bootstrap))
})

test_that("maihda_app helpers leave summaries unchanged without a bootstrap interval", {
  res <- suppressWarnings(suppressMessages(
    MAIHDA:::maihda_app_fit_models(
      dat = MAIHDA::maihda_sim_data[seq_len(120), ],
      outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
      family = "gaussian", use_boot = FALSE, engine = "lme4"
    )
  ))
  # No bootstrap requested => no VPC intervals carried back.
  expect_null(res$vpc_ci_null)
  expect_null(res$vpc_ci_adjusted)

  summ <- summary(res$null_model)
  expect_identical(MAIHDA:::maihda_app_attach_vpc_ci(summ, NULL), summ)
  expect_false(MAIHDA:::maihda_vpc_has_interval(MAIHDA:::maihda_app_attach_vpc_ci(summ, NULL)$vpc))
})

maihda_source_app_for_test <- function() {
  for (pkg in MAIHDA:::maihda_app_required_packages()) {
    skip_if_not_installed(pkg)
  }

  app_file <- system.file("shiny", "app.R", package = "MAIHDA")
  if (app_file == "") {
    app_file <- file.path("inst", "shiny", "app.R")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  app_env <- new.env(parent = globalenv())
  suppressPackageStartupMessages(suppressWarnings(sys.source(app_file, envir = app_env)))
  app_env
}

test_that("Shiny app stores previous future plan for restoration", {
  app_env <- maihda_source_app_for_test()
  expect_true(exists("maihda_app_previous_future_plan", envir = app_env, inherits = FALSE))
})

test_that("Shiny server loads data and derives HUD plot data from real results", {
  app_env <- maihda_source_app_for_test()
  dat <- MAIHDA::maihda_sim_data[seq_len(120), ]
  res <- MAIHDA:::maihda_app_fit_models(
    dat = dat,
    outcome_var = "health_outcome",
    grouping_vars = c("gender", "race"),
    additional_covars = "age",
    family = "gaussian",
    use_boot = FALSE,
    n_boot = 10,
    autobin = TRUE,
    engine = "lme4"
  )

  shiny::testServer(app_env$server, {
    session$setInputs(
      dataset = "sim",
      hud_top_n = 5,
      hud_sort_var = "effect",
      hud_color_var = "deviant"
    )

    expect_equal(nrow(reactive_data()), nrow(MAIHDA::maihda_sim_data))

    model_results(res$model)
    null_summary_results(summary(res$null_model))
    summary_results(summary(res$model))
    pvc_results(res$pvc)

    hud <- hud_plot_data()
    expect_lte(nrow(hud), 5)
    expect_true(all(c("display_label", "deviant", "random_effect") %in% names(hud)))
    expect_false("stratum.1" %in% names(hud))
    expect_true(any(hud$display_label != paste0("Stratum ", hud$stratum)))
  })
})
