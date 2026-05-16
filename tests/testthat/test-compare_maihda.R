test_that("plot_comparison validates required columns before plotting", {
  expect_error(
    plot_comparison(data.frame(vpc = c(0.1, 0.2))),
    "model",
    fixed = TRUE
  )

  p <- plot_comparison(data.frame(model = c("A", "B"), vpc = c(0.1, 0.2)))

  expect_s3_class(p, "ggplot")
})
