# Regression tests for the 2026-07-13 (P1/P1/P2) audit findings.
#
#   P1 #1  an explicit ordinal model whose top category sat only on rows with a
#          missing predictor silently collapsed to a binomial clmm fit: the
#          category count was validated on the FULL data (before the analytic-
#          sample drop), so the vanished category slipped through.
#   P1 #2  the longitudinal analytic-sample mask used complete.cases() on the RAW
#          columns, so a non-finite transformed predictor (log(x) of x <= 0) was
#          retained -- mis-anchoring the time centre and raising a false
#          "id ... appear in more than one stratum" error.
#   P2 #3  estimation = "fitted" is the default, but pcv_importance(),
#          stepwise_pcv() and compare_maihda_groups() did not retain or print the
#          chosen estimator, losing that provenance from serialized results.

# ---- shared fixtures --------------------------------------------------------

maihda_test_long_df <- function(seed = 2) {
  set.seed(seed)
  np <- 30
  waves <- c(10, 11, 12, 13)                 # 4 occasions -> identifiable slope
  d <- data.frame(pid  = rep(seq_len(np), each = length(waves)),
                  wave = rep(waves, np))
  d$gender <- rep(sample(c("F", "M"), np, replace = TRUE), each = length(waves))
  d$edu    <- rep(sample(c("lo", "hi"), np, replace = TRUE), each = length(waves))
  d$x <- runif(nrow(d), 0.5, 3)              # positive -> log(x) finite
  d$y <- rnorm(nrow(d))
  d
}

maihda_test_pcv_df <- function(seed = 99) {
  set.seed(seed)
  N <- 18 * 22
  d <- data.frame(g   = factor(sample(letters[1:3], N, replace = TRUE)),
                  r   = factor(sample(LETTERS[1:6], N, replace = TRUE)),
                  age = rnorm(N))
  d <- make_strata(d, c("g", "r"))$data
  slev <- unique(as.character(d$stratum))
  eff <- stats::setNames(rnorm(length(slev), sd = 2), slev)
  d$y <- 3 + 0.5 * d$age + eff[as.character(d$stratum)] + rnorm(nrow(d))
  d
}

# ---- P1 #1: ordinal category count is checked on the ANALYTIC sample --------

test_that("maihda_ordinal_assert_min_levels errors below 3 analytic categories", {
  f3 <- factor(c("a", "b", "c", "a", "b", "c"), ordered = TRUE)
  expect_true(MAIHDA:::maihda_ordinal_assert_min_levels(f3, "y"))
  f2 <- droplevels(factor(c("a", "b", "a", "b"),
                          levels = c("a", "b", "c"), ordered = TRUE))
  expect_error(MAIHDA:::maihda_ordinal_assert_min_levels(f2, "y"),
               "at least 3 response categories", fixed = TRUE)
})

test_that("an ordinal top category confined to missing-predictor rows is rejected", {
  skip_if_not_installed("ordinal")
  set.seed(1)
  n <- 240
  d <- data.frame(gender = sample(c("F", "M"), n, replace = TRUE),
                  edu    = sample(c("lo", "hi"), n, replace = TRUE),
                  x      = rnorm(n))
  yy <- sample(c("1", "2"), n, replace = TRUE)
  idx3 <- sample(seq_len(n), 6)
  yy[idx3] <- "3"
  d$x[idx3] <- NA                            # the ONLY level-3 rows miss x
  d$y <- factor(yy, levels = c("1", "2", "3"), ordered = TRUE)
  expect_error(
    suppressWarnings(suppressMessages(
      fit_maihda(y ~ x + (1 | gender:edu), data = d,
                 family = "ordinal", engine = "ordinal"))),
    "at least 3 response categories", fixed = TRUE)
})

test_that("the brms ordinal path rejects the same case before any Stan compile", {
  skip_if_not_installed("brms")
  set.seed(1)
  n <- 240
  d <- data.frame(gender = sample(c("F", "M"), n, replace = TRUE),
                  edu    = sample(c("lo", "hi"), n, replace = TRUE),
                  x      = rnorm(n))
  yy <- sample(c("1", "2"), n, replace = TRUE)
  idx3 <- sample(seq_len(n), 6)
  yy[idx3] <- "3"
  d$x[idx3] <- NA
  d$y <- factor(yy, levels = c("1", "2", "3"), ordered = TRUE)
  # The analytic-sample re-check fires in argument handling, before brms::brm().
  expect_error(
    suppressWarnings(suppressMessages(
      fit_maihda(y ~ x + (1 | gender:edu), data = d,
                 family = "ordinal", engine = "brms"))),
    "at least 3 response categories", fixed = TRUE)
})

test_that("a clean 3-category ordinal fit is unaffected", {
  skip_if_not_installed("ordinal")
  set.seed(2)
  n <- 300
  d <- data.frame(gender = sample(c("F", "M"), n, replace = TRUE),
                  edu    = sample(c("lo", "hi"), n, replace = TRUE),
                  x      = rnorm(n))
  d$y <- factor(sample(c("1", "2", "3"), n, replace = TRUE),
                levels = c("1", "2", "3"), ordered = TRUE)
  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | gender:edu), data = d,
               family = "ordinal", engine = "ordinal")))
  expect_equal(length(fit$model$alpha), 2)   # 3 categories -> 2 thresholds
  expect_equal(nlevels(fit$data$y), 3)
})

# ---- P1 #2: the longitudinal analytic mask is transformation-aware ----------

test_that("a non-finite transformed predictor does not mis-anchor the time centre", {
  skip_on_cran()
  d <- maihda_test_long_df()
  bad <- d[d$pid == 1, ][1, ]
  bad$wave <- 0; bad$x <- -1                 # log(-1) = NaN, dropped by lmer
  d_bad <- rbind(d, bad)
  fit_bad <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ log(x) + (1 | gender:edu), data = d_bad,
               id = "pid", time = "wave")))
  # The excluded wave-0 row must not drag the centre down from 10 to 0.
  expect_equal(fit_bad$longitudinal_info$time_center, 10)
})

test_that("an excluded transformed-NA row raises no false cross-stratum id error", {
  skip_on_cran()
  d <- maihda_test_long_df()
  bad <- d[d$pid == 1, ][1, ]
  bad$wave <- 14; bad$x <- -1
  bad$gender <- if (bad$gender == "F") "M" else "F"   # a DIFFERENT stratum
  d_bad <- rbind(d, bad)
  expect_error(
    suppressWarnings(suppressMessages(
      fit_maihda(y ~ log(x) + (1 | gender:edu), data = d_bad,
                 id = "pid", time = "wave"))),
    NA)                                       # must NOT error
})

# ---- P2 #3: the estimation basis is retained and printed --------------------

test_that("pcv_importance retains and prints the estimation basis", {
  skip_on_cran()
  d <- maihda_test_pcv_df()
  imp <- pcv_importance(d, "y", c("g", "r"))
  expect_identical(imp$estimation, "fitted")
  expect_match(paste(capture.output(print(imp)), collapse = "\n"),
               "Variance basis: as fitted")
  imp_ml <- pcv_importance(d, "y", c("g", "r"), estimation = "ML")
  expect_identical(imp_ml$estimation, "ML")
  expect_match(paste(capture.output(print(imp_ml)), collapse = "\n"),
               "Variance basis: ML-refit")
})

test_that("stepwise_pcv records the estimation basis (attribute + print)", {
  skip_on_cran()
  d <- maihda_test_pcv_df()
  sw <- stepwise_pcv(d, "y", c("g", "r"))
  expect_identical(attr(sw, "estimation"), "fitted")
  expect_match(paste(capture.output(print(sw)), collapse = "\n"),
               "Variance basis: as fitted")
  sw_ml <- stepwise_pcv(d, "y", c("g", "r"), estimation = "ML")
  expect_identical(attr(sw_ml, "estimation"), "ML")
})

test_that("compare_maihda_groups records the estimation basis; it survives subsetting", {
  skip_on_cran()
  set.seed(7)
  n <- 60 * 8
  d <- data.frame(gender = sample(c("F", "M"), n, replace = TRUE),
                  race   = sample(c("A", "B"), n, replace = TRUE),
                  grp    = sample(c("g1", "g2"), n, replace = TRUE),
                  age    = rnorm(n))
  d$y <- rnorm(n)
  cmp <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ age + (1 | gender:race), data = d, group = "grp")))
  expect_identical(attr(cmp, "estimation"), "fitted")
  # The [ method must carry the attribute across a row subset.
  expect_identical(attr(cmp[1, ], "estimation"), "fitted")
  expect_match(paste(capture.output(print(cmp)), collapse = "\n"),
               "Variance basis: as fitted")
})
