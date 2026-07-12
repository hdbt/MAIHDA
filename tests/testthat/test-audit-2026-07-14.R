# Regression tests for the 2026-07-14 statistical audit (two P1, two P2, one P3
# doc-only). Each test pins the corrected behaviour so the specific defect cannot
# silently return.

test_that("calculate_pcv rejects crossed-dimensions models (finding 1)", {
  skip_on_cran()
  set.seed(7001)
  n <- 1200
  d <- data.frame(
    a  = sample(c("a1", "a2", "a3"), n, TRUE),
    b  = sample(c("b1", "b2", "b3"), n, TRUE),
    cc = sample(c("c1", "c2", "c3"), n, TRUE),
    x  = rnorm(n), stringsAsFactors = FALSE)
  ua <- stats::setNames(rnorm(3, sd = 0.9), c("a1", "a2", "a3"))
  ub <- stats::setNames(rnorm(3, sd = 0.7), c("b1", "b2", "b3"))
  uc <- stats::setNames(rnorm(3, sd = 0.5), c("c1", "c2", "c3"))
  st <- interaction(d$a, d$b, d$cc, drop = TRUE)
  d$y <- 2 + 0.3 * d$x + ua[d$a] + ub[d$b] + uc[d$cc] +
    rnorm(nlevels(st), sd = 0.8)[st] + rnorm(n)

  ccm <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d,
           decomposition = "crossed-dimensions")))$model
  expect_false(is.null(ccm$cc_info))

  # The between-stratum variance of a crossed-dimensions fit is additive + interaction
  # (summary$decomposition$between_var), NOT the single intersection ("stratum")
  # component. The two-model PCV path returns only the intersection variance, so it
  # must reject a cc fit rather than silently divide by the wrong (and possibly
  # sign-reversing) denominator.
  expect_error(calculate_pcv(ccm, ccm), "crossed-dimensions")
  expect_error(extract_between_variance(ccm), "crossed-dimensions")

  # A canonical (single-stratum) two-model PCV is unaffected.
  an <- suppressWarnings(suppressMessages(maihda(y ~ x + (1 | a:b:cc), data = d)))
  expect_true(is.finite(an$pcv$pcv))
})

test_that("auto-bin cut-points ignore rows the fit drops (finding 2)", {
  skip_on_cran()
  set.seed(101)
  n <- 600
  d <- data.frame(income = round(runif(n, 1, 120), 2),
                  sex = sample(c("M", "F"), n, TRUE), x = rnorm(n),
                  stringsAsFactors = FALSE)
  d$y <- 1 + 0.02 * d$income + 0.3 * d$x + rnorm(n)
  # Extreme, non-missing values of the binning variable carried by rows the fit
  # drops. Pre-fix these leaked into the tertile quantiles and moved OTHER rows'
  # strata; the cut-points must instead match the pre-filtered analytic sample.
  bad_income <- c(500, 700, 900, 1060, 480, 620, 350, 800)
  brk <- function(fit) fit$strata_autobin_info$income$breaks

  # (a) rows excluded by a MISSING OUTCOME
  bad <- data.frame(income = bad_income, sex = "M", x = rnorm(8), y = NA_real_)
  fa <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | income:sex), rbind(d, bad))))
  fc <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | income:sex), d)))
  expect_equal(nobs(fa$model), nobs(fc$model))
  expect_equal(brk(fa), brk(fc))

  # (b) rows excluded by an INVALID (zero) SAMPLING WEIGHT
  d2 <- d
  d2$w <- runif(n, 0.5, 2)
  badw <- data.frame(income = bad_income, sex = "M", x = rnorm(8),
                     y = rnorm(8), w = 0)
  fb <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | income:sex), rbind(d2, badw), sampling_weights = "w")))
  fbc <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | income:sex), d2, sampling_weights = "w")))
  expect_equal(brk(fb), brk(fbc))

  # (c) rows excluded by a MISSING PRECISION WEIGHT
  badpw <- data.frame(income = bad_income, sex = "M", x = rnorm(8),
                      y = rnorm(8), w = NA_real_)
  fp <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | income:sex), rbind(d2, badpw), weights = w)))
  fpc <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | income:sex), d2, weights = w)))
  expect_equal(brk(fp), brk(fpc))
})

test_that("longitudinal id/stratum validation uses the analytic sample (finding 3)", {
  skip_on_cran()
  set.seed(5)
  np <- 40
  d <- expand.grid(pid = 1:np, wave = 0:3)
  d$gender <- ifelse(d$pid %% 2 == 0, "F", "M")     # person-level (constant within pid)
  d$ses    <- ifelse(d$pid <= 20, "low", "high")    # person-level
  d$y <- rnorm(nrow(d)) + 0.1 * d$wave
  d$pid <- as.character(d$pid)
  # A contaminating row assigning person "1" (normally M) to a SECOND stratum.
  bad <- data.frame(pid = "1", wave = 4, gender = "F", ses = "low",
                    y = 0, stringsAsFactors = FALSE)
  dall <- rbind(d, bad)
  keep <- c(rep(TRUE, nrow(d)), FALSE)              # subset EXCLUDES the bad row

  # With the offending row excluded by `subset`, the cross-stratum id check must see
  # the analytic sample and pass -- identical to fitting the pre-filtered data.
  f_sub <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | gender:ses), dall, id = "pid", time = "wave", subset = keep)))
  expect_s3_class(f_sub, "maihda_model")
  f_pre <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | gender:ses), dall[keep, ], id = "pid", time = "wave")))
  expect_equal(nobs(f_sub$model), nobs(f_pre$model))

  # The guard is intact: the SAME contamination, NOT excluded, still errors.
  expect_error(suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | gender:ses), dall, id = "pid", time = "wave"))),
    "more than one stratum")
})

test_that("brms group-variance parser handles group names containing __ (finding 4)", {
  set.seed(1)
  nd <- 40
  # brms names a group-level SD column sd_<group>__<coef>; a group literally named
  # site__id yields sd_site__id__Intercept. The parser must recover "site__id", not
  # the "site" left by splitting on the FIRST "__".
  draws <- data.frame(
    `sd_site__id__Intercept` = abs(rnorm(nd, 0.7, 0.05)),
    `sd_stratum__Intercept`  = abs(rnorm(nd, 0.9, 0.05)),
    sigma = abs(rnorm(nd, 1, 0.02)), check.names = FALSE)
  gv <- maihda_group_variance_draws_brms(draws)
  expect_setequal(names(gv), c("site__id", "stratum"))
  expect_equal(gv[["site__id"]], draws[["sd_site__id__Intercept"]]^2)

  # A single-underscore dimension name is unaffected.
  draws1 <- data.frame(
    `sd_.maihda_dim_age__Intercept` = abs(rnorm(nd, 0.5, 0.05)),
    `sd_stratum__Intercept`         = abs(rnorm(nd, 0.8, 0.05)),
    sigma = abs(rnorm(nd, 1, 0.02)), check.names = FALSE)
  expect_true(".maihda_dim_age" %in% names(maihda_group_variance_draws_brms(draws1)))

  # Random slopes on a __-named group are still rejected (two coefs for one group).
  draws2 <- data.frame(
    `sd_site__id__Intercept` = abs(rnorm(nd, 0.7, 0.05)),
    `sd_site__id__wave`      = abs(rnorm(nd, 0.2, 0.05)),
    sigma = abs(rnorm(nd, 1, 0.02)), check.names = FALSE)
  expect_error(maihda_group_variance_draws_brms(draws2), "site__id")
})
