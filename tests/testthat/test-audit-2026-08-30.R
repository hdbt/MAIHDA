# Regression tests for the audit pass of 2026-08-30.
#
# Finding [P2] "automatic strata construction changes the estimand by default"
# asked for make_strata(autobin = ) to default to FALSE. REFUTED: on a purely
# additive DGP with a nonlinear tertile effect, autobin = TRUE reproduces the
# hand-binned reference exactly (PCV 1.000) while the SAME call with
# autobin = FALSE returns PCV 0.109 -- ~89% spurious "interaction" -- and a
# genuinely continuous dimension gives an opaque lme4 "number of levels ... must
# be < number of observations" error instead. The default therefore stands.
#
# What the pass DID change is transparency. The cut-points were announced only by
# make_strata()'s fit-time message(), which is gone once a model is saved and
# reloaded; they are now reported by the print methods as well, so the binning
# recipe stays on the page with the numbers it produced.

test_that("maihda_print_autobin_info prints nothing when nothing was binned", {
  expect_output(maihda_print_autobin_info(NULL), NA)
  expect_output(maihda_print_autobin_info(list()), NA)
  expect_output(maihda_print_autobin_info(list(a = 1)), NA)   # unnamed/odd shape
  # A recipe whose breaks are unusable is skipped, not half-printed with a
  # dangling header.
  expect_output(
    maihda_print_autobin_info(list(v = list(breaks = c(NA, 1, 2, 3),
                                            labels = c("a", "b", "c")))), NA)
  expect_output(
    maihda_print_autobin_info(list(v = list(breaks = numeric(0),
                                            labels = character(0)))), NA)
})

test_that("maihda_print_autobin_info reports the cut-points and their labels", {
  out <- cli::ansi_strip(capture.output(maihda_print_autobin_info(list(
    a = list(breaks = c(1, 2, 3, 4), labels = c("a_Low", "a_Mid", "a_High")),
    b = list(breaks = c(3623.41, 17383.6, 28891.5, 216678),
             labels = c("b_Low", "b_Mid", "b_High"))))))
  expect_true(any(grepl("Auto-binned dimensions", out, fixed = TRUE)))
  expect_true(any(grepl("a: 1, 2, 3, 4 -> a_Low, a_Mid, a_High", out, fixed = TRUE)))
  # Each cut-point is formatted on its own: a vector-wide format() would pad the
  # small values to the widest one's decimals ("17383.60" next to "216678.00").
  expect_true(any(grepl("b: 3623.41, 17383.6, 28891.5, 216678", out, fixed = TRUE)))
  expect_false(any(grepl("17383.60", out, fixed = TRUE)))
})

test_that("maihda_print_autobin_info honours the indent argument", {
  out <- cli::ansi_strip(capture.output(maihda_print_autobin_info(
    list(a = list(breaks = 1:4, labels = c("l", "m", "h"))), indent = "   ")))
  out <- out[nzchar(trimws(out))]
  expect_gt(length(out), 0)
  expect_true(all(startsWith(out, "   ")))
})

test_that("summary.maihda_model carries the auto-bin recipe and prints it", {
  skip_on_cran()
  set.seed(2026)
  n <- 400
  d <- data.frame(income = round(stats::rlnorm(n, 10, 0.6), 2),
                  gender = sample(c("m", "f"), n, TRUE))
  d$y <- stats::rnorm(n) + as.numeric(d$income > stats::median(d$income))

  m <- suppressMessages(maihda(y ~ (1 | income:gender), data = d))
  brk <- m$model$strata_autobin_info$income$breaks
  expect_length(brk, 4)

  s <- suppressMessages(summary(m$model))
  expect_equal(s$strata_autobin_info$income$breaks, brk)

  first <- format(brk[1], digits = 6, trim = TRUE)
  out_s <- cli::ansi_strip(capture.output(print(s)))
  expect_true(any(grepl("Auto-binned dimensions", out_s, fixed = TRUE)))
  expect_true(any(grepl(first, out_s, fixed = TRUE)))
  expect_true(any(grepl("income_Low", out_s, fixed = TRUE)))

  out_a <- cli::ansi_strip(capture.output(print(m)))
  expect_true(any(grepl("Auto-binned dimensions", out_a, fixed = TRUE)))
  expect_true(any(grepl(first, out_a, fixed = TRUE)))
})

test_that("the crossed-dimensions print branch reports the cut-points", {
  skip_on_cran()
  set.seed(11)
  n <- 400
  d <- data.frame(income = round(stats::rlnorm(n, 10, 0.6), 2),
                  gender = sample(c("m", "f"), n, TRUE))
  d$y <- stats::rnorm(n) + as.numeric(d$income > stats::median(d$income))

  m <- suppressMessages(maihda(y ~ (1 | income:gender), data = d,
                               decomposition = "crossed-dimensions"))
  brk <- m$model$strata_autobin_info$income$breaks
  out <- cli::ansi_strip(capture.output(print(m)))
  expect_true(any(grepl("Auto-binned dimensions", out, fixed = TRUE)))
  expect_true(any(grepl(format(brk[2], digits = 6, trim = TRUE), out, fixed = TRUE)))
})

test_that("the longitudinal print branch reports the cut-points", {
  skip_on_cran()
  m <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ (1 | age:gender), data = maihda_long_data,
           id = "id", time = "wave")))
  brk <- m$model$strata_autobin_info$age$breaks
  expect_length(brk, 4)
  out <- cli::ansi_strip(capture.output(print(m)))
  expect_true(any(grepl("Auto-binned dimensions", out, fixed = TRUE)))
  expect_true(any(grepl("age_Low", out, fixed = TRUE)))
})

test_that("a categorical-dimension fit prints no auto-bin block", {
  skip_on_cran()
  m <- suppressMessages(maihda(health_outcome ~ (1 | gender:race), data = maihda_sim_data))
  expect_length(m$model$strata_autobin_info, 0L)
  out <- cli::ansi_strip(capture.output(print(m)))
  expect_false(any(grepl("Auto-binned", out, fixed = TRUE)))
  out_s <- cli::ansi_strip(capture.output(print(suppressMessages(summary(m$model)))))
  expect_false(any(grepl("Auto-binned", out_s, fixed = TRUE)))
})

test_that("the reported cut-points are the ones make_strata() actually used", {
  # The printed recipe must round-trip: cutting the raw column at the reported
  # breaks has to reproduce the stratum-defining bins, or the block would be
  # documenting a different binning from the one behind the VPC/PCV.
  set.seed(5)
  d <- data.frame(income = round(stats::rlnorm(300, 10, 0.6), 2),
                  gender = sample(c("m", "f"), 300, TRUE))
  st <- suppressMessages(make_strata(d, vars = c("income", "gender")))
  info <- st$autobin_info$income
  expect_length(info$breaks, 4)

  # Re-cut the RAW column at the reported breaks and check it reproduces, row by
  # row, the bin each observation was actually assigned -- strata_info[[v]] holds
  # the binned value behind every stratum, so this compares the printed recipe
  # against the strata the model was fitted on rather than against itself.
  rebuilt <- cut(d$income, breaks = info$breaks, include.lowest = TRUE,
                 labels = info$labels)
  assigned <- st$strata_info$income[st$data$stratum]
  expect_false(anyNA(rebuilt))
  expect_false(anyNA(assigned))
  expect_identical(as.character(rebuilt), as.character(assigned))

  # Teeth: the same comparison against SHIFTED breaks must fail, or the check
  # above would pass for any cut-points at all.
  shifted <- info$breaks
  shifted[2] <- shifted[2] + diff(range(d$income)) / 10
  wrong <- cut(d$income, breaks = shifted, include.lowest = TRUE,
               labels = info$labels)
  expect_false(identical(as.character(wrong), as.character(assigned)))

  # ... and the labels the block prints are the ones in the stratum labels.
  expect_true(all(vapply(levels(rebuilt), function(l) {
    any(grepl(l, st$strata_info$label, fixed = TRUE))
  }, logical(1))))
})
