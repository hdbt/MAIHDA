# Tidy / glance (broom interface) methods for the MAIHDA classes. These repackage
# the summary slots, so the tests assert the broom-shaped column contracts and the
# headline glance() across the lme4, wemix, ordinal and (opt-in) brms engines.

make_tidy_data <- function(seed = 5151, n = 600) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE),
    age    = rnorm(n),
    stringsAsFactors = FALSE
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  u  <- rnorm(nlevels(sk), sd = 0.7)[sk]
  d$y    <- 1 + 0.3 * d$age + u + rnorm(n, sd = 0.5)
  d$ybin <- rbinom(n, 1, stats::plogis(-0.2 + 0.4 * d$age + u))
  lat    <- u + 0.3 * d$age + rlogis(n)
  d$yord <- factor(cut(lat, c(-Inf, -0.8, 0.6, Inf), labels = 1:3), ordered = TRUE)
  d$w    <- runif(n, 0.5, 3)
  d
}

# gender (2) x race (3), all combinations observed at n = 600
N_STRATA <- 6L

test_that("tidy.maihda_model returns broom-shaped strata / variance / fixed tibbles", {
  d <- make_tidy_data()
  m <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  st <- tidy(m)
  expect_s3_class(st, "tbl_df")
  expect_identical(names(st),
                   c("stratum", "label", "estimate", "std.error", "conf.low", "conf.high"))
  expect_equal(nrow(st), N_STRATA)
  expect_true(is.numeric(st$estimate) && all(is.finite(st$estimate)))

  vc <- tidy(m, component = "variance")
  expect_identical(names(vc), c("component", "variance", "sd", "proportion"))
  expect_true("Total" %in% vc$component)
  expect_equal(vc$proportion[vc$component == "Total"], 1)

  fe <- tidy(m, component = "fixed")
  expect_identical(names(fe),
                   c("term", "estimate", "std.error", "statistic", "p.value",
                     "conf.low", "conf.high"))
  expect_true(all(c("(Intercept)", "age") %in% fe$term))
})

# Regression (2026-08-09): the lme4 summary used to store term/estimate only, so
# tidy(component = "fixed") returned an all-NA std.error and no interval.
test_that("tidy(component = 'fixed') carries lme4 SEs, statistics and Wald intervals", {
  d <- make_tidy_data()
  m <- fit_maihda(y ~ age + (1 | gender:race), data = d)
  fe <- tidy(m, component = "fixed")

  cf <- stats::coef(summary(m$model))
  expect_equal(fe$estimate, unname(cf[, "Estimate"]), tolerance = 1e-10)
  expect_equal(fe$std.error, unname(cf[, "Std. Error"]), tolerance = 1e-10)
  expect_true(all(is.finite(fe$std.error)) && all(fe$std.error > 0))

  # Wald z, its two-sided normal p, and the matching 95% interval.
  expect_equal(fe$statistic, fe$estimate / fe$std.error, tolerance = 1e-10)
  expect_equal(fe$p.value, 2 * stats::pnorm(-abs(fe$statistic)), tolerance = 1e-10)
  expect_equal(fe$conf.low,
               fe$estimate - stats::qnorm(0.975) * fe$std.error, tolerance = 1e-10)
  expect_equal(fe$conf.high,
               fe$estimate + stats::qnorm(0.975) * fe$std.error, tolerance = 1e-10)
  # age has a real effect (0.3) at n = 600 -> interval excludes zero.
  expect_true(fe$conf.low[fe$term == "age"] > 0)
  expect_lt(fe$p.value[fe$term == "age"], 0.001)
})

test_that("summary(conf_level =) sets the fixed-effect interval level", {
  d <- make_tidy_data()
  m <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  fe99 <- tidy(summary(m, conf_level = 0.99), component = "fixed")
  fe95 <- tidy(m, component = "fixed")
  expect_equal(fe99$conf.low,
               fe99$estimate - stats::qnorm(0.995) * fe99$std.error,
               tolerance = 1e-10)
  expect_true(all(fe99$conf.low < fe95$conf.low))
  expect_true(all(fe99$conf.high > fe95$conf.high))
  expect_equal(fe99$p.value, fe95$p.value, tolerance = 1e-10)  # level-free

  expect_error(summary(m, conf_level = 1.5), "conf_level")
})

test_that("glance.maihda_model is a one-row tibble with the headline columns", {
  d <- make_tidy_data()
  m <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  g <- glance(m)
  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1L)
  expect_identical(
    names(g),
    c("vpc", "vpc.conf.low", "vpc.conf.high", "additive.share", "interaction.share",
      "auc", "mor", "n_strata", "nobs", "engine", "family")
  )
  expect_equal(g$vpc, summary(m)$vpc$estimate, tolerance = 1e-8)
  expect_equal(g$n_strata, N_STRATA)
  expect_equal(g$nobs, nrow(d))
  expect_identical(g$engine, "lme4")
  expect_true(g$vpc > 0 && g$vpc < 1)
  expect_true(is.na(g$auc))   # gaussian -> no discriminatory accuracy
})

test_that("glance surfaces discriminatory accuracy for a binomial fit", {
  d <- make_tidy_data()
  m <- suppressWarnings(
    fit_maihda(ybin ~ age + (1 | gender:race), data = d, family = "binomial"))

  g <- glance(m)
  expect_true(is.finite(g$auc) && g$auc > 0.4 && g$auc <= 1)
  expect_true(is.finite(g$mor))   # logit link -> MOR defined
})

test_that("glance.maihda_analysis adds PCV, adjusted-model and mode columns", {
  d <- make_tidy_data()
  a <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race), data = d))

  g <- glance(a)
  expect_equal(nrow(g), 1L)
  expect_identical(
    names(g),
    c("vpc", "vpc.conf.low", "vpc.conf.high",
      "pcv", "pcv.conf.low", "pcv.conf.high",
      "additive.share", "interaction.share",
      "auc", "auc.adjusted", "mor",
      "n_strata", "nobs", "engine", "family", "mode")
  )
  expect_equal(g$pcv, a$pcv$pcv, tolerance = 1e-8)
  expect_identical(g$mode, "two-model")
  expect_equal(g$nobs, nrow(d))
  expect_equal(g$n_strata, N_STRATA)
})

test_that("tidy.maihda_analysis selects the null vs adjusted model", {
  d <- make_tidy_data()
  a <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race), data = d))

  fe_null <- tidy(a, component = "fixed", which = "null")
  fe_adj  <- tidy(a, component = "fixed", which = "adjusted")
  expect_true("age" %in% fe_null$term)
  expect_gt(nrow(fe_adj), nrow(fe_null))   # adjusted adds gender/race main effects
  expect_equal(nrow(tidy(a)), N_STRATA)    # strata estimates from the null model
})

test_that("tidy.maihda_analysis errors when the requested summary is absent", {
  d <- make_tidy_data()
  a <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race), data = d))
  a$summary_adjusted <- NULL
  expect_error(tidy(a, which = "adjusted"), "No 'adjusted' summary")
})

test_that("tidy/glance work for the wemix engine (design-weighted)", {
  skip_if_not_installed("WeMix")
  d <- make_tidy_data()
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race), data = d,
               engine = "wemix", sampling_weights = "w")))

  fe <- tidy(m, component = "fixed")
  expect_true(all(c("term", "estimate", "std.error") %in% names(fe)))
  expect_true(any(is.finite(fe$std.error)))   # WeMix reports sandwich (robust) SEs
  expect_true(all(is.finite(fe$conf.low)) && all(fe$conf.low < fe$conf.high))

  g <- glance(m)
  expect_identical(g$engine, "wemix")
  expect_equal(g$n_strata, N_STRATA)
  expect_true(g$vpc > 0 && g$vpc < 1)
})

test_that("tidy/glance work for the ordinal engine", {
  skip_if_not_installed("ordinal")
  d <- make_tidy_data()
  m <- suppressWarnings(suppressMessages(
    fit_maihda(yord ~ age + (1 | gender:race), data = d, family = "ordinal")))

  st <- tidy(m)
  expect_equal(nrow(st), N_STRATA)

  # Location coefficients only -- a cumulative model's "intercepts" are the
  # thresholds, reported by summary()$thresholds.
  fe <- tidy(m, component = "fixed")
  expect_identical(fe$term, "age")
  expect_true(all(is.finite(fe$std.error)) && all(is.finite(fe$conf.low)))
  expect_true(all(fe$conf.low < fe$estimate & fe$estimate < fe$conf.high))

  g <- glance(m)
  expect_identical(g$engine, "ordinal")
  expect_true(g$vpc > 0 && g$vpc < 1)
})

test_that("tidy/glance work for the brms engine", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  d <- make_tidy_data()
  # The old chains = 1, iter = 600 fit reached max Rhat 1.075 with a bulk ESS of
  # 36, and passed no seed at all -- so the Stan seed came from the ambient RNG
  # and the block was not reproducible either. brms warned out of the glance()
  # and tidy() calls below while every assertion passed. Only 6 strata here, so
  # the intercept/sd_stratum funnel needs adapt_delta above the default to clear
  # the zero-divergence half of the package's convergence bar.
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race), data = d, engine = "brms",
               chains = 4, iter = 2000, warmup = 1000, refresh = 0, seed = 101,
               control = list(adapt_delta = 0.95))))
  # Fixture guard -- see the note in test-ordinal-engine.R.
  expect_true(isTRUE(m$diagnostics$converged))

  g <- glance(m)
  expect_identical(g$engine, "brms")
  expect_true(is.finite(g$vpc.conf.low) && is.finite(g$vpc.conf.high))  # posterior CI
  # The band must actually cover the truth, not merely be finite. make_tidy_data
  # draws 6 stratum effects at sd 0.7, realizing sd 0.6107 (var 0.3730) against a
  # residual variance of 0.5^2, so the VPC is 0.3730 / (0.3730 + 0.25) = 0.599.
  expect_true(g$vpc.conf.low < 0.599 && 0.599 < g$vpc.conf.high)
  expect_true(g$vpc >= g$vpc.conf.low && g$vpc <= g$vpc.conf.high)

  fe <- tidy(m, component = "fixed")
  expect_true(all(c("conf.low", "conf.high") %in% names(fe)))
  expect_true(all(is.finite(fe$conf.low)) && all(is.finite(fe$conf.high)))
  # any(is.finite(.)) passed on a single finite cell. Assert the estimate the
  # posterior is actually meant to recover: age enters at 0.3.
  age <- fe[fe$term == "age", , drop = FALSE]
  expect_equal(nrow(age), 1L)
  expect_true(age$conf.low < 0.3 && 0.3 < age$conf.high)
  expect_true(age$conf.low < age$estimate && age$estimate < age$conf.high)
})
