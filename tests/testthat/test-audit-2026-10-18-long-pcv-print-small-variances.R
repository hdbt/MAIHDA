# Audit pass 2026-10-18: variances below 1e-4 in print.maihda_long_pcv().
#
# A slope variance is per unit of time squared. A growth fit of maihda_long_data with
# time in days stores a between-stratum slope variance of 3.215e-07 (null) and 4.318e-08
# (adjusted), and the print's four decimals showed both as 0.0000; in weeks both, in
# months the adjusted one. A variance that is not zero and below 1e-4 now prints to four
# significant digits. Every other value prints with four decimals, as before.
#
# Random numbers are drawn only inside keep_rng(), which puts the suite's stream back.

keep_rng <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  force(expr)
}

quiet <- function(e) suppressWarnings(suppressMessages(e))

# A maihda_long_pcv as maihda_longitudinal_pcv() stores it for a linear stratum block.
# The defaults are maihda_long_data's in days.
small_var_pcv <- function(vb = c(0.407694594, 0.03321365647),
                          vs = c(3.214653219e-07, 4.31753632e-08)) {
  block <- function(v0, v1) {
    matrix(c(v0, 0, 0, v1), 2, 2,
           dimnames = list(c("(Intercept)", "wave"), c("(Intercept)", "wave")))
  }
  structure(
    list(pcv_intercept = 1 - vb[2] / vb[1],
         pcv_slope = if (vs[1] > 0) 1 - vs[2] / vs[1] else NA_real_,
         var_baseline_null = vb[1], var_baseline_adjusted = vb[2],
         var_slope_null = vs[1], var_slope_adjusted = vs[2],
         ref_time = 0, time_center = 0,
         pcv_t = data.frame(time = c(0, 1460), var_null = vb[1], var_adjusted = vb[2],
                            pcv = 1 - vb[2] / vb[1]),
         Sigma_stratum_null = block(vb[1], vs[1]),
         Sigma_stratum_adjusted = block(vb[2], vs[2]),
         time = "wave", time_degree = 1L, ml_refit = FALSE, estimation = "fitted",
         estimation_used = "fitted", null_at_boundary = FALSE,
         adjusted_at_boundary = FALSE, null_rank_deficient = FALSE,
         null_slope_at_boundary = FALSE, adjusted_rank_deficient = FALSE,
         stratum_slope = TRUE, variance_summary = "estimate",
         Sigma_summary = "estimate"),
    class = "maihda_long_pcv")
}

printed_lines <- function(x) cli::ansi_strip(utils::capture.output(print(x)))
# The printed line that starts with `which` ("Baseline" or "Slope").
variance_line <- function(x, which) {
  out <- printed_lines(x)
  out[startsWith(out, which)]
}

test_that("a variance below 1e-4 is rendered to four significant digits", {
  at_least <- c(0.4077447403, 0.03320826683, 0.04284050265, 0.005751945193,
                2.975450511e-04, 2e-04, 1e-04, 1, 123.456)
  expect_identical(maihda_format_variance(at_least), sprintf("%.4f", at_least))
  expect_identical(
    maihda_format_variance(c(3.214653219e-07, 4.31753632e-08, 3.994149333e-05,
                             1.584720313e-05, 2.126840965e-06)),
    c("3.215e-07", "4.318e-08", "3.994e-05", "1.585e-05", "2.127e-06"))
  # The switch is at 1e-4; four significant digits keep their trailing zeros.
  expect_identical(maihda_format_variance(c(1e-04, 9.9e-05, 1e-05)),
                   c("0.0001", "9.900e-05", "1.000e-05"))
  # Zero, missing and infinite values print as they did.
  expect_identical(maihda_format_variance(c(0, NA, Inf)), c("0.0000", "NA", "Inf"))
  expect_identical(maihda_format_variance(NA), "NA")
  expect_identical(maihda_format_variance(c(NA, 3.214653219e-07, NaN)),
                   c("NA", "3.215e-07", "NaN"))
  expect_identical(maihda_format_variance(NULL), character(0))
})

test_that("the longitudinal PCV print shows a slope variance in days", {
  x <- small_var_pcv()
  expect_identical(
    variance_line(x, "Slope"),
    "Slope (wave) variance at baseline: 3.215e-07 (null) -> 4.318e-08 (adjusted)")
  expect_identical(variance_line(x, "Baseline"),
                   "Baseline (wave = 0) variance: 0.4077 (null) -> 0.0332 (adjusted)")
  # In months the null value is above 1e-4 and keeps its four decimals.
  m <- small_var_pcv(vs = c(2.975450511e-04, 3.994149333e-05))
  expect_identical(
    variance_line(m, "Slope"),
    "Slope (wave) variance at baseline: 0.0003 (null) -> 3.994e-05 (adjusted)")
  # A baseline variance below 1e-4, from an outcome in small units, is shown the same way.
  b <- small_var_pcv(vb = c(4.2e-05, 6.1e-06))
  expect_identical(
    variance_line(b, "Baseline"),
    "Baseline (wave = 0) variance: 4.200e-05 (null) -> 6.100e-06 (adjusted)")
  # An exact zero, as for an adjusted model without a stratum slope, stays 0.0000.
  z <- small_var_pcv(vs = c(0.0324, 0))
  expect_identical(
    variance_line(z, "Slope"),
    "Slope (wave) variance at baseline: 0.0324 (null) -> 0.0000 (adjusted)")
  # An object without var_slope_* reads the slope cells of its blocks, formatted alike.
  cells <- x
  cells$var_slope_null <- NULL
  cells$var_slope_adjusted <- NULL
  expect_identical(variance_line(cells, "Slope"), variance_line(x, "Slope"))
})

test_that("variances of at least 1e-4 print with four decimals, as before", {
  # maihda_long_data in waves.
  vb <- c(0.4077447403, 0.03320826683)
  vs <- c(0.04284050265, 0.005751945193)
  w <- small_var_pcv(vb = vb, vs = vs)
  expect_identical(
    variance_line(w, "Baseline"),
    sprintf("Baseline (wave = 0) variance: %.4f (null) -> %.4f (adjusted)", vb[1], vb[2]))
  expect_identical(
    variance_line(w, "Slope"),
    sprintf("Slope (wave) variance at baseline: %.4f (null) -> %.4f (adjusted)",
            vs[1], vs[2]))
  expect_identical(
    variance_line(small_var_pcv(vs = c(0.002, 2e-04)), "Slope"),
    "Slope (wave) variance at baseline: 0.0020 (null) -> 0.0002 (adjusted)")
})

test_that("maihda_long_data in days prints the slope variances it stores", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  d <- maihda_long_data
  d$wave <- d$wave * 365
  a <- keep_rng(quiet(maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                             data = d, id = "id", time = "wave",
                             decomposition = "longitudinal")))
  p <- a$pcv
  expect_lt(p$var_slope_null, 1e-06)
  expect_lt(p$var_slope_adjusted, 1e-06)
  expect_gt(p$var_slope_adjusted, 0)
  slope <- sprintf("Slope (wave) variance at baseline: %s (null) -> %s (adjusted)",
                   sprintf("%.3e", p$var_slope_null), sprintf("%.3e", p$var_slope_adjusted))
  expect_true(slope %in% printed_lines(p))
  # print() on the analysis shows the PCV the same way.
  expect_true(slope %in% printed_lines(a))
  expect_false(any(grepl("0.0000 (null)", printed_lines(a), fixed = TRUE)))
})
