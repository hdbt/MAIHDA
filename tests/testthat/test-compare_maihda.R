test_that("plot() on a maihda_comparison validates required columns before plotting", {
  bad <- structure(data.frame(vpc = c(0.1, 0.2)),
                   class = c("maihda_comparison", "data.frame"))
  expect_error(plot(bad), "model", fixed = TRUE)

  good <- structure(data.frame(model = c("A", "B"), vpc = c(0.1, 0.2)),
                    class = c("maihda_comparison", "data.frame"))
  expect_s3_class(plot(good), "ggplot")
})

test_that("compare_maihda output is a maihda_comparison and plots via plot()", {
  set.seed(2201)
  d <- data.frame(stratum = rep(seq_len(8), each = 12), x = rnorm(96))
  d$y <- 1 + d$x + rnorm(8, sd = 0.8)[d$stratum] + rnorm(96, sd = 0.4)
  m1 <- fit_maihda(y ~ x + (1 | stratum), data = d)
  m2 <- fit_maihda(y ~ 1 + (1 | stratum), data = d)

  cmp <- compare_maihda(m1, m2)
  expect_s3_class(cmp, "maihda_comparison")
  expect_s3_class(cmp, "data.frame")          # still a data.frame
  expect_s3_class(plot(cmp), "ggplot")
})

test_that("plot_comparison() is deprecated but still works", {
  df <- data.frame(model = c("A", "B"), vpc = c(0.1, 0.2))
  expect_warning(p <- plot_comparison(df), "deprecated")
  expect_s3_class(p, "ggplot")
})

test_that("compare_maihda warns when models use different families", {
  set.seed(2202)
  d <- data.frame(stratum = rep(seq_len(8), each = 25), x = rnorm(200))
  u <- rnorm(8, sd = 0.7)[d$stratum]
  d$y_gauss <- 1 + 0.3 * d$x + u + rnorm(200, sd = 0.4)
  d$y_count <- rpois(200, lambda = exp(0.3 + 0.2 * d$x + u))

  m_gauss <- fit_maihda(y_gauss ~ x + (1 | stratum), data = d)
  m_pois  <- fit_maihda(y_count ~ x + (1 | stratum), data = d, family = "poisson")

  # These models differ in BOTH outcome and family: still a single warning.
  # The Poisson model's marginal count is below 2, so summarising it also emits
  # the count level-1 approximation note; drop that one -- the COMPARABILITY
  # check itself must stay a single warning.
  w <- testthat::capture_warnings(compare_maihda(m_gauss, m_pois))
  w <- grep("count VPC/ICC level-1 variance", w, value = TRUE, invert = TRUE)
  expect_length(w, 1)
  expect_match(w, "differ in")
  expect_match(w, "outcomes")
  expect_match(w, "families/links")
})

test_that("compare_maihda warns when models use different analytic samples", {
  set.seed(2203)
  d <- data.frame(stratum = rep(seq_len(8), each = 25), x = rnorm(200), z = rnorm(200))
  u <- rnorm(8, sd = 0.7)[d$stratum]
  d$y <- 1 + 0.3 * d$x + u + rnorm(200, sd = 0.4)
  d$z[seq_len(40)] <- NA_real_   # z drops 40 rows from the second model's sample

  m_full <- fit_maihda(y ~ x + (1 | stratum), data = d)          # n = 200
  m_sub  <- fit_maihda(y ~ x + z + (1 | stratum), data = d)       # n = 160

  w <- testthat::capture_warnings(compare_maihda(m_full, m_sub))
  expect_length(w, 1)
  expect_match(w, "analytic sample")
})

test_that("compare_maihda warns when strata definitions differ despite shared IDs", {
  set.seed(88)
  n <- 200
  d <- data.frame(
    a = sample(c("p", "q"), n, replace = TRUE),
    b = sample(c("x", "y"), n, replace = TRUE),
    cc = sample(c("m", "n"), n, replace = TRUE),
    age = rnorm(n)
  )
  d$y <- rnorm(n)
  # Both models have strata numbered 1..4, but defined from different variables.
  m_ab <- fit_maihda(y ~ age + (1 | a:b), data = d)
  m_ac <- fit_maihda(y ~ age + (1 | a:cc), data = d)

  w <- testthat::capture_warnings(compare_maihda(m_ab, m_ac))
  expect_length(w, 1)
  expect_match(w, "stratum definitions")
})

test_that("compare_maihda warns for disjoint analytic samples with the same n and strata", {
  set.seed(2024)
  big <- data.frame(
    stratum = factor(rep(seq_len(8), times = 40)),  # all 8 strata throughout
    x = rnorm(320)
  )
  big$y <- 1 + 0.3 * big$x + rnorm(8, sd = 0.7)[big$stratum] + rnorm(320, sd = 0.4)
  d1 <- big[1:160, ]      # rows 1..160
  d2 <- big[161:320, ]    # disjoint rows; same n = 160 and same 8 strata

  m1 <- fit_maihda(y ~ x + (1 | stratum), data = d1)
  m2 <- fit_maihda(y ~ x + (1 | stratum), data = d2)

  w <- testthat::capture_warnings(compare_maihda(m1, m2))
  expect_length(w, 1)
  expect_match(w, "analytic sample")
})

test_that("compare_maihda includes interval columns only when an interval exists", {
  set.seed(99)
  d <- data.frame(stratum = rep(seq_len(8), each = 12), x = rnorm(96))
  d$y <- 1 + d$x + rnorm(8, sd = 0.9)[d$stratum] + rnorm(96, sd = 0.4)
  m1 <- fit_maihda(y ~ x + (1 | stratum), data = d)
  m2 <- fit_maihda(y ~ 1 + (1 | stratum), data = d)

  # lme4 with bootstrap -> interval columns present
  cmp_ci <- suppressWarnings(compare_maihda(m1, m2, bootstrap = TRUE, n_boot = 25))
  expect_true(all(c("ci_lower", "ci_upper") %in% names(cmp_ci)))

  # lme4 without bootstrap -> no interval available -> columns dropped
  cmp_noci <- suppressWarnings(compare_maihda(m1, m2))
  expect_false("ci_lower" %in% names(cmp_noci))
})
