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

test_that("Shiny server loads the selected built-in dataset", {
  app_env <- maihda_source_app_for_test()

  shiny::testServer(app_env$server, {
    session$setInputs(dataset = "pisa")
    expect_equal(nrow(reactive_data()), nrow(MAIHDA::maihda_country_data))

    session$setInputs(dataset = "health")
    expect_equal(nrow(reactive_data()), nrow(MAIHDA::maihda_health_data))
  })
})

test_that("A single grouping variable triggers the non-intersectional hint", {
  app_env <- maihda_source_app_for_test()

  shiny::testServer(app_env$server, {
    session$setInputs(dataset = "pisa", group_vars = "gender")
    hint_one <- paste(unlist(output$group_var_hint), collapse = " ")
    expect_match(hint_one, "intersectional", ignore.case = TRUE)

    session$setInputs(group_vars = c("gender", "race"))
    hint_two <- paste(unlist(output$group_var_hint), collapse = " ")
    expect_false(grepl("intersectional", hint_two, ignore.case = TRUE))
  })
})

test_that("Explorer module derives HUD plot data from real results", {
  for (pkg in MAIHDA:::maihda_app_required_packages()) skip_if_not_installed(pkg)
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

  shiny::testServer(
    MAIHDA:::mod_explorer_server,
    args = list(
      model_results = shiny::reactiveVal(res$model),
      null_summary_results = shiny::reactiveVal(summary(res$null_model)),
      summary_results = shiny::reactiveVal(summary(res$model)),
      pvc_results = shiny::reactiveVal(res$pvc),
      group_vars = shiny::reactiveVal(c("gender", "race"))
    ),
    expr = {
      session$setInputs(hud_top_n = 5, hud_sort_var = "effect", hud_color_var = "deviant")

      hud <- hud_plot_data()
      expect_lte(nrow(hud), 5)
      expect_true(all(c("display_label", "deviant", "random_effect") %in% names(hud)))
      expect_false("stratum.1" %in% names(hud))
      expect_true(any(hud$display_label != paste0("Stratum ", hud$stratum)))
      # The per-stratum absolute predicted outcome is aggregated and merged in (it
      # feeds the plot tooltip); a broken aggregate() call would silently drop it.
      expect_true("abs_pred" %in% names(hud))
      expect_true(all(is.finite(hud$abs_pred)))

      # The deviance plot (colour + shape -> plotly) builds without error.
      expect_no_error(output$interactive_plot)

      # The highlighted-strata CSV export runs end-to-end: accessing the download
      # output invokes the handler's write.csv() content path and returns the file.
      csv_path <- output$download_hud_data
      expect_true(file.exists(csv_path))
      exported <- utils::read.csv(csv_path)
      expect_false(any(c("tooltip", "display_label") %in% names(exported)))
      expect_equal(nrow(exported), nrow(hud))
    }
  )
})

test_that("Visualizations module produces a plot object from the fitted model", {
  for (pkg in MAIHDA:::maihda_app_required_packages()) skip_if_not_installed(pkg)
  res <- suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(120), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    family = "gaussian"
  ))

  shiny::testServer(
    MAIHDA:::mod_visualizations_server,
    args = list(model_results = shiny::reactiveVal(res$model)),
    expr = {
      session$setInputs(plot_type = "vpc")
      p <- expect_no_error(current_plot())
      expect_false(is.null(p))
    }
  )
})

test_that("maihda_app_fit_models reports the resolved family and auto-switch flag", {
  set.seed(11)
  d <- data.frame(
    Obese = factor(sample(c("No", "Yes"), 240, replace = TRUE)),
    Gender = sample(c("F", "M"), 240, replace = TRUE),
    Race = sample(c("A", "B"), 240, replace = TRUE)
  )
  # A binary outcome left on the default gaussian is fit as binomial; the result
  # records both the requested and the resolved family plus the switch flag, so
  # the dashboard can tell the user what was actually fit.
  res <- suppressMessages(MAIHDA:::maihda_app_fit_models(
    d, outcome_var = "Obese", grouping_vars = c("Gender", "Race"), family = "gaussian"
  ))
  expect_equal(res$family_requested, "gaussian")
  expect_equal(res$family_used, "binomial")
  expect_true(res$family_autoswitched)

  # A continuous outcome is fit as requested -- no switch.
  res2 <- suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(120), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    family = "gaussian"
  ))
  expect_equal(res2$family_used, "gaussian")
  expect_false(res2$family_autoswitched)
})

test_that("maihda_app_fit_models bins an auto-binned numeric stratum dimension in the adjusted model", {
  skip_on_cran()
  set.seed(1)
  d <- MAIHDA::maihda_health_data
  d <- d[sample(nrow(d), min(500L, nrow(d))), , drop = FALSE]

  # Gender x Race x Age, with Age (numeric) auto-binned into tertiles to define the
  # strata. The adjusted model's main effect for Age must be that SAME tertile factor
  # (.maihda_dim_Age), not a raw linear Age term -- matching core maihda().
  res <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    dat = d, outcome_var = "Obese",
    grouping_vars = c("Gender", "Race", "Age"),
    family = "gaussian", autobin = TRUE, engine = "lme4"
  )))

  fixed_labels <- attr(stats::terms(reformulas::nobars(res$model$formula)), "term.labels")
  expect_true(".maihda_dim_Age" %in% fixed_labels)
  expect_false("Age" %in% fixed_labels)
  expect_true(".maihda_dim_Age" %in% names(res$model$data))
})

test_that("Reproducible code export mirrors the fitted model specification", {
  code <- MAIHDA:::maihda_app_generate_code(
    outcome_var = "math",
    grouping_vars = c("gender", "ses"),
    additional_covars = "escs",
    family = "gaussian",
    autobin = TRUE, use_boot = TRUE, n_boot = 200, seed = 7,
    dataset = "pisa"
  )
  expect_true(grepl("data <- MAIHDA::maihda_country_data", code, fixed = TRUE))
  # Complete-case preprocessing mirrors maihda_app_fit_models() and precedes
  # make_strata() so the auto-bin cut-points match the dashboard's fit.
  expect_true(grepl('model_vars <- c("math", "gender", "ses", "escs")', code, fixed = TRUE))
  expect_true(grepl("data <- data[complete.cases(data[, model_vars, drop = FALSE]), , drop = FALSE]",
                    code, fixed = TRUE))
  expect_true(grepl('make_strata(data, vars = c("gender", "ses"), autobin = TRUE)', code, fixed = TRUE))
  # Keeping make_strata()'s returned data carries the auto-bin recipe forward.
  expect_true(grepl("data <- strata$data", code, fixed = TRUE))
  # The make_strata() complete-case step must come before strata construction.
  expect_lt(regexpr("complete.cases", code, fixed = TRUE),
            regexpr("make_strata", code, fixed = TRUE))
  # With >= 2 grouping vars the script delegates the two-model decomposition to
  # maihda(): the formula's fixed part is the covariate (the null model), and the
  # stratum dimensions' additive main effects (incl. any auto-binned tertile factor)
  # are added by maihda() for the adjusted model -- so covariates are held in BOTH.
  expect_true(grepl('analysis <- maihda(math ~ escs + (1 | stratum), data = data, family = "gaussian", bootstrap = TRUE, n_boot = 200)',
                    code, fixed = TRUE))
  # No hand-rolled raw-term adjusted formula (the source of the binning mismatch).
  expect_false(grepl("math ~ gender + ses + escs + (1 | stratum)", code, fixed = TRUE))
  expect_true(grepl('family = "gaussian"', code, fixed = TRUE))
  expect_true(grepl("set.seed(7)", code, fixed = TRUE))
  # The stepwise PCV must rerun under the same family as the dashboard fit (a
  # Poisson/binomial app fit would otherwise reproduce as Gaussian).
  expect_true(grepl('stepwise_pcv(data, outcome = "math", vars = c("gender", "ses", "escs"), family = "gaussian")',
                    code, fixed = TRUE))
})

test_that("Reproducible code export omits seed/bootstrap when unused and handles uploads", {
  code <- MAIHDA:::maihda_app_generate_code(
    outcome_var = "y", grouping_vars = "g1", additional_covars = character(),
    family = "gaussian", autobin = FALSE, use_boot = FALSE, seed = NULL,
    dataset = "upload", upload_name = "mydata.csv"
  )
  expect_true(grepl('read.csv("mydata.csv")', code, fixed = TRUE))
  expect_true(grepl("autobin = FALSE", code, fixed = TRUE))
  expect_true(grepl("data <- strata$data", code, fixed = TRUE))
  expect_false(grepl("set.seed", code, fixed = TRUE))
  expect_false(grepl("bootstrap = TRUE", code, fixed = TRUE))
  # A single grouping variable is not an intersectional MAIHDA, so there is no
  # additive/intersectional decomposition: fall back to a single strata-only fit.
  expect_false(grepl("analysis <- maihda(", code, fixed = TRUE))
  expect_true(grepl('model <- fit_maihda(y ~ 1 + (1 | stratum), data = data, family = "gaussian")',
                    code, fixed = TRUE))
  expect_true(grepl('stepwise_pcv(data, outcome = "y", vars = c("g1"), family = "gaussian")',
                    code, fixed = TRUE))
})

test_that("Reproducible code export escapes string literals in names and uploads", {
  # Column names / upload filenames are user-controlled and may contain a double
  # quote. A naive paste would break the script (or inject code when sourced), so
  # every emitted literal must be escaped with encodeString(quote = '"').
  evil_outcome <- 'y"); system("oops"); read.csv("x'
  code <- MAIHDA:::maihda_app_generate_code(
    outcome_var = evil_outcome,
    grouping_vars = c('g"1', "g2"),
    additional_covars = character(),
    family = "gaussian", autobin = FALSE, use_boot = FALSE, seed = NULL,
    dataset = "upload", upload_name = 'a"b.csv'
  )

  # The generated script is still valid, self-contained R (no broken/injected
  # literals from the embedded quotes).
  expect_silent(parse(text = code))

  # Each literal appears in its properly escaped form...
  expect_true(grepl(encodeString('a"b.csv', quote = '"'), code, fixed = TRUE))
  expect_true(grepl(encodeString('g"1', quote = '"'), code, fixed = TRUE))
  expect_true(grepl(encodeString(evil_outcome, quote = '"'), code, fixed = TRUE))

  # ...and never in the unescaped form that would close the literal early.
  expect_false(grepl('read.csv("a"b.csv")', code, fixed = TRUE))
})

test_that("A fixed seed makes the bootstrap intervals reproducible across fits", {
  skip_on_cran()
  fit <- function() suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    dat = MAIHDA::maihda_sim_data[seq_len(150), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    additional_covars = "age", family = "gaussian",
    use_boot = TRUE, n_boot = 25, engine = "lme4", seed = 99
  )))
  r1 <- fit()
  r2 <- fit()
  expect_equal(r1$vpc_ci_null, r2$vpc_ci_null)
  expect_equal(r1$vpc_ci_adjusted, r2$vpc_ci_adjusted)
})

test_that("Model-comparison module renders the nested null-vs-adjusted comparison", {
  for (pkg in MAIHDA:::maihda_app_required_packages()) skip_if_not_installed(pkg)
  res <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(150), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    additional_covars = "age", family = "gaussian"
  )))
  cmp <- compare_maihda(res$null_model, res$model,
                        model_names = c("Model 1: Null", "Model 2: Adjusted"))
  expect_s3_class(cmp, "maihda_comparison")
  expect_equal(nrow(cmp), 2)
  expect_true(all(c("model", "vpc") %in% names(cmp)))

  shiny::testServer(
    MAIHDA:::mod_compare_server,
    args = list(
      comparison_results = shiny::reactiveVal(cmp),
      reactive_data = shiny::reactiveVal(MAIHDA::maihda_sim_data),
      fit_params = shiny::reactiveVal(list(
        outcome = "health_outcome", grouping_vars = c("gender", "race"),
        covariates = "age", autobin = TRUE)),
      fitted_family = shiny::reactiveVal("gaussian")
    ),
    expr = {
      # The nested comparison table renders from the precomputed result.
      expect_no_error(output$nested_table)
    }
  )
})

test_that("maihda_app_default_vars gives dataset-aware defaults and a heuristic fallback", {
  pisa <- MAIHDA:::maihda_app_default_vars("pisa", MAIHDA::maihda_country_data)
  expect_equal(pisa$outcome, "math")
  expect_equal(pisa$groups, c("gender", "ses"))

  health <- MAIHDA:::maihda_app_default_vars("health", MAIHDA::maihda_health_data)
  expect_equal(health$outcome, "Obese")
  expect_equal(health$groups, c("Gender", "Race", "Age"))

  # Upload / unknown dataset: first column as outcome, first two categorical-ish
  # columns as strata (continuous columns are not strata candidates).
  up <- data.frame(
    y = rnorm(20),
    region = rep(c("N", "S"), 10),
    sex = rep(c("F", "M"), 10),
    income = rnorm(20)
  )
  res <- MAIHDA:::maihda_app_default_vars("upload", up)
  expect_equal(res$outcome, "y")
  expect_equal(res$groups, c("region", "sex"))

  expect_equal(MAIHDA:::maihda_app_default_vars("upload", data.frame()),
               list(outcome = NULL, groups = character(0)))
})

test_that("Selecting the NHANES dataset yields fittable defaults (outcome + 2 strata)", {
  app_env <- maihda_source_app_for_test()
  shiny::testServer(app_env$server, {
    session$setInputs(dataset = "health")
    expect_equal(nrow(reactive_data()), nrow(MAIHDA::maihda_health_data))
    spec <- MAIHDA:::maihda_app_default_vars("health", reactive_data())
    expect_equal(spec$outcome, "Obese")
    expect_gte(length(spec$groups), 2)
  })
})

test_that("maihda_app_fit_models supports the cross-classified decomposition", {
  dat <- MAIHDA::maihda_sim_data[seq_len(300), ]
  res <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    dat = dat, outcome_var = "health_outcome",
    grouping_vars = c("gender", "race", "education"),
    family = "gaussian", decomposition = "cross-classified"
  )))
  expect_identical(res$decomposition_mode, "cross-classified")
  expect_false(is.null(res$model$cc_info))
  expect_s3_class(res$summary_obj, "maihda_summary")
  expect_null(res$pvc)
  expect_null(res$stepwise)
  d <- res$decomposition
  expect_false(is.null(d))
  expect_equal(d$additive_share + d$interaction_share, 1, tolerance = 1e-8)
  expect_equal(d$additive_var + d$interaction_var, d$between_var, tolerance = 1e-8)
})

test_that("the two-model app fit still tags decomposition_mode", {
  res <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(150), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    family = "gaussian"
  )))
  expect_identical(res$decomposition_mode, "two-model")
  expect_s3_class(res$pvc, "pvc_result")
})

test_that("cross-classified app fit needs at least two grouping variables", {
  dat <- MAIHDA::maihda_sim_data[seq_len(120), ]
  expect_error(
    MAIHDA:::maihda_app_fit_models(
      dat, outcome_var = "health_outcome", grouping_vars = "gender",
      family = "gaussian", decomposition = "cross-classified"),
    "two grouping")
})

test_that("Reproducible code mirrors the cross-classified decomposition", {
  code <- MAIHDA:::maihda_app_generate_code(
    outcome_var = "math", grouping_vars = c("gender", "ses"),
    additional_covars = "escs", family = "gaussian",
    autobin = TRUE, use_boot = FALSE, seed = NULL, dataset = "pisa",
    decomposition = "cross-classified"
  )
  expect_true(grepl('decomposition = "cross-classified"', code, fixed = TRUE))
  expect_true(grepl("analysis$decomposition$additive_share", code, fixed = TRUE))
  # Cross-classified mode has no two-model PCV or stepwise step.
  expect_false(grepl("analysis$pcv", code, fixed = TRUE))
  expect_false(grepl("stepwise_pcv(", code, fixed = TRUE))
})

test_that("Explorer renders in cross-classified mode from the decomposition", {
  for (pkg in MAIHDA:::maihda_app_required_packages()) skip_if_not_installed(pkg)
  dat <- MAIHDA::maihda_sim_data[seq_len(300), ]
  res <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    dat = dat, outcome_var = "health_outcome",
    grouping_vars = c("gender", "race", "education"),
    family = "gaussian", decomposition = "cross-classified")))

  shiny::testServer(
    MAIHDA:::mod_explorer_server,
    args = list(
      model_results = shiny::reactiveVal(res$model),
      null_summary_results = shiny::reactiveVal(res$summary_obj),
      summary_results = shiny::reactiveVal(res$summary_obj),
      pvc_results = shiny::reactiveVal(NULL),
      group_vars = shiny::reactiveVal(c("gender", "race", "education")),
      decomposition_results = shiny::reactiveVal(res$decomposition)
    ),
    expr = {
      session$setInputs(hud_top_n = 5, hud_sort_var = "effect", hud_color_var = "deviant")
      ui_txt <- paste(unlist(output$interactive_explorer_ui), collapse = " ")
      expect_false(grepl("Nothing to explore yet", ui_txt))
      expect_match(ui_txt, "Additive share")
      expect_no_error(hud_plot_data())
    }
  )
})

test_that("Full server renders the cross-classified decomposition on the PCV tab", {
  for (pkg in MAIHDA:::maihda_app_required_packages()) skip_if_not_installed(pkg)
  app_env <- maihda_source_app_for_test()
  dat <- MAIHDA::maihda_sim_data[seq_len(300), ]
  res <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    dat = dat, outcome_var = "health_outcome",
    grouping_vars = c("gender", "race", "education"),
    family = "gaussian", decomposition = "cross-classified")))

  shiny::testServer(app_env$server, {
    session$setInputs(dataset = "pisa")
    model_results(res$model)
    summary_results(res$summary_obj)
    null_summary_results(res$summary_obj)
    decomposition_results(res$decomposition)
    decomposition_mode("cross-classified")
    pvc_results(NULL)
    stepwise_results(NULL)
    session$flushReact()

    pcv_txt <- paste(unlist(output$pvc_summary_ui), collapse = " ")
    expect_match(pcv_txt, "Additive vs. Intersectional Decomposition")
    expect_match(pcv_txt, "Additive share")

    step_txt <- paste(unlist(output$stepwise_pcv_ui), collapse = " ")
    expect_match(step_txt, "Stepwise PCV not used here")
  })
})

test_that("The reproducible-code output mirrors the last fit (Reproduce in R)", {
  app_env <- maihda_source_app_for_test()

  shiny::testServer(app_env$server, {
    session$setInputs(dataset = "pisa")
    # repro_code_text() needs these three reactiveVals non-NULL; the script content
    # is built from fit_params + the resolved family (no fitted model object needed).
    model_results(list(placeholder = TRUE))
    fitted_family("gaussian")
    fit_params(list(
      dataset = "pisa", upload_name = NULL, outcome = "math",
      grouping_vars = c("gender", "ses"), covariates = "escs",
      autobin = TRUE, use_boot = FALSE, n_boot = 100, seed = NULL
    ))

    code <- repro_code_text()
    expect_match(code, "MAIHDA::maihda_country_data", fixed = TRUE)
    # >= 2 grouping vars -> maihda() decomposition; covariate in the null formula,
    # stratum main effects added by maihda() for the adjusted model.
    expect_match(code, 'analysis <- maihda(math ~ escs + (1 | stratum), data = data, family = "gaussian")',
                 fixed = TRUE)
  })
})
