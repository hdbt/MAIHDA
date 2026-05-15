test_that("Shiny app dependency gate includes ternary plotting dependency", {
  expect_true("ggtern" %in% MAIHDA:::maihda_app_required_packages())
})
