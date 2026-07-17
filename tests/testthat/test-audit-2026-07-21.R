# Regression tests for the 2026-07-21 audit findings on commit 684b6fa.
#
#   #1 P1  compare_maihda_groups() and maihda() re-derived the analytic sample
#          WITHOUT fit_maihda()'s rule that zero / negative / non-finite lme4
#          precision weights are mapped to NA (the rows lme4 drops). A zero-weight
#          row carrying an out-of-sample outcome value therefore flipped the
#          detected family/engine, and a group's analytic n was counted including
#          rows the fit drops (passing min_group_n). Now all three functions share
#          maihda_normalize_precision_weights() before detection / binning / guards.
#   #2 P2  maihda_resolve_strata_formula() rejected duplicate stratum intercepts
#          only when EVERY stratum term was exactly a constant intercept, so the
#          compound (1 | stratum) + (1 + x | stratum) -- and the implicit-intercept
#          (1 | stratum) + (x | stratum) -- were accepted although both terms carry
#          an intercept (lme4 splits the between-stratum variance across
#          'stratum'/'stratum.1'). Now each term is tested for an intercept via its
#          formula intercept attribute (maihda_re_lhs_has_intercept), catching the
#          compound/implicit forms while still allowing (0 + x | stratum).
#   #3 P2  Design-weighted brms fits pruned invalid-weight / incomplete rows and
#          stored that filtered frame as original_data, so maihda_describe() could
#          not report the original total / missingness / excluded-weight counts and
#          made total == analytic. Now the pre-prune frame is kept as original_data.

# ---- #1 P1: high-level APIs exclude invalid precision-weight rows -------------

test_that("compare_maihda_groups() family detection ignores invalid precision-weight rows", {
  skip_on_cran()
  set.seed(101)
  n <- 400
  d <- data.frame(
    a = sample(c("a1", "a2", "a3"), n, TRUE),
    b = sample(c("b1", "b2"), n, TRUE),
    grp = sample(c("g1", "g2"), n, TRUE),
    x = rnorm(n))
  d$y <- rbinom(n, 1, 0.5)             # 0/1 outcome
  d$w <- 1
  d$w[1] <- 0; d$y[1] <- 2             # a zero-weight row carries an out-of-sample value

  # fit_maihda() drops the zero-weight row, sees {0, 1}, and detects binomial.
  fam_fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | a:b), data = d, weights = w)$family$family))
  expect_identical(fam_fit, "binomial")

  # compare_maihda_groups() must resolve the SAME family. Previously it kept the raw
  # weights, saw {0, 1, 2}, and stayed gaussian for every group.
  cmp <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ x + (1 | a:b), data = d, group = "grp",
                          weights = w, min_group_n = 5)))
  expect_identical(attr(cmp, "family"), "binomial")
})

test_that("compare_maihda_groups() min_group_n counts the valid-weight analytic sample", {
  skip_on_cran()
  mk <- function(grp, n_valid, n_zero) {
    m <- n_valid + n_zero
    data.frame(a = rep(c("a1", "a2"), length.out = m),
               b = rep(c("b1", "b2"), length.out = m),
               grp = grp, x = rnorm(m), y = rnorm(m),
               w = c(rep(1, n_valid), rep(0, n_zero)))
  }
  set.seed(7)
  d <- rbind(mk("small", 8, 6),        # 8 valid + 6 zero-weight = 14 raw rows
             mk("big", 60, 0))
  res <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ x + (1 | a:b), data = d, group = "grp",
                          family = "gaussian", weights = w, min_group_n = 10)))
  small <- res[res$group == "small", ]
  # 8 analytic rows < min_group_n = 10, so the group is skipped (previously the raw
  # count of 14 passed the guard and the group was fitted with only 8 rows).
  expect_match(small$status, "skipped")
  expect_equal(small$n, 8)
})

test_that("maihda() ordinal-engine detection ignores invalid precision-weight rows", {
  skip_on_cran()
  set.seed(202)
  n <- 300
  d <- data.frame(a = sample(c("a1", "a2", "a3"), n, TRUE),
                  b = sample(c("b1", "b2"), n, TRUE),
                  x = rnorm(n))
  lvl <- ifelse(runif(n) < 0.5, "lo", "mid")
  d$w <- 1
  lvl[1:5] <- "hi"; d$w[1:5] <- 0      # a third ordered level only on zero-weight rows
  d$y <- factor(lvl, levels = c("lo", "mid", "hi"), ordered = TRUE)

  # On the valid-weight analytic sample the outcome has two observed levels, so the
  # binomial (lme4) model is correct. Previously maihda() detected three levels on
  # the raw column, selected engine = "ordinal", and then errored on the binary
  # analytic sample.
  m <- suppressWarnings(suppressMessages(
    maihda(y ~ (1 | a:b), data = d, weights = w)))
  expect_identical(m$model$engine, "lme4")
  expect_identical(m$model$family$family, "binomial")
})

# ---- #2 P2: compound / implicit duplicate stratum intercepts are rejected -----

test_that("maihda_re_lhs_has_intercept detects explicit and implicit intercepts", {
  f <- MAIHDA:::maihda_re_lhs_has_intercept
  expect_true(f(quote(1)))
  expect_true(f(quote(0 + 1)))
  expect_true(f(quote(1 + x)))
  expect_true(f(quote(x)))             # implicit: lme4 reads (x | g) as (1 + x | g)
  expect_false(f(quote(0 + x)))        # slope only
  expect_false(f(quote(x - 1)))        # slope only
})

test_that("fit_maihda() rejects a compound duplicate stratum intercept, allows an uncorrelated slope", {
  skip_on_cran()
  set.seed(303)
  ns <- 30; per <- 15
  strat <- rep(sprintf("s%02d", seq_len(ns)), each = per)
  N <- length(strat)
  sidx <- as.integer(factor(strat))
  u0 <- rnorm(ns, 0, 1)[sidx]
  u1 <- rnorm(ns, 0, 0.5)[sidx]
  d <- data.frame(stratum = strat, x = rnorm(N))
  d$y <- 2 + u0 + u1 * d$x + 0.5 * d$x + rnorm(N)

  # (1 + x | stratum) carries an intercept, so together with (1 | stratum) there are
  # two stratum intercepts -- lme4 splits the between-stratum variance arbitrarily.
  expect_error(
    suppressMessages(fit_maihda(y ~ x + (1 | stratum) + (1 + x | stratum), data = d)),
    "non-identifiable")
  # (x | stratum) has an IMPLICIT intercept (lme4 expands it to 1 + x) -> also rejected.
  expect_error(
    suppressMessages(fit_maihda(y ~ x + (1 | stratum) + (x | stratum), data = d)),
    "non-identifiable")
  # (0 + x | stratum) is slope-only (no intercept) -> identifiable, still accepted.
  expect_no_error(
    suppressWarnings(suppressMessages(
      fit_maihda(y ~ x + (1 | stratum) + (0 + x | stratum), data = d))))
})

# ---- #3 P2: weighted brms fits keep the full pre-fit data as original_data ----

test_that("design-weighted brms fits keep the full pre-fit data as original_data", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(nzchar(Sys.getenv("MAIHDA_TEST_BRMS")),
              "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  set.seed(404)
  n <- 80
  d <- data.frame(a = sample(c("a1", "a2", "a3"), n, TRUE),
                  b = sample(c("b1", "b2"), n, TRUE),
                  x = rnorm(n), sw = runif(n, 0.5, 2))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  d$sw[1:3] <- 0                       # invalid sampling weights (dropped by brms)
  d$x[10:11] <- NA                     # rows incomplete on a covariate (dropped)
  n_analytic <- sum(is.finite(d$sw) & d$sw > 0 & !is.na(d$x))

  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | a:b), data = d, engine = "brms",
               sampling_weights = "sw", chains = 1, iter = 200,
               refresh = 0, seed = 1)))

  # original_data is the FULL pre-fit frame (so maihda_describe() can report the
  # original total, missingness, and excluded-weight counts); $data is the analytic
  # sample the engine actually fitted. Previously both were the pruned frame.
  expect_equal(nrow(fit$original_data), n)
  expect_equal(nrow(fit$data), n_analytic)
  expect_gt(nrow(fit$original_data), nrow(fit$data))
})
