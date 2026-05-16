test_that("Shiny app dependency gate includes ternary plotting dependency", {
  expect_true("ggtern" %in% MAIHDA:::maihda_app_required_packages())
})

test_that("Shiny app dependency gate leaves upload-only readers optional", {
  expect_false("haven" %in% MAIHDA:::maihda_app_required_packages())
})

test_that("Shiny HUD keeps null and adjusted summaries separate", {
  app_file <- system.file("shiny", "app.R", package = "MAIHDA")
  if (app_file == "") {
    app_file <- file.path("inst", "shiny", "app.R")
  }
  app <- readLines(app_file, warn = FALSE)

  expect_true(any(grepl("null_summary_results <- reactiveVal(NULL)", app, fixed = TRUE)))
  expect_true(any(grepl("null_summary_results(summary(res$null_model))", app, fixed = TRUE)))
  expect_true(any(grepl("vpc_val <- round(null_res$vpc$estimate * 100, 2)", app, fixed = TRUE)))
})

test_that("Shiny stepwise PCV hover text uses explicit added-variable column", {
  app_file <- system.file("shiny", "app.R", package = "MAIHDA")
  if (app_file == "") {
    app_file <- file.path("inst", "shiny", "app.R")
  }
  app <- readLines(app_file, warn = FALSE)

  expect_true(any(grepl("res$Added_Variable", app, fixed = TRUE)))
  expect_false(any(grepl("res$Added)", app, fixed = TRUE)))
})
