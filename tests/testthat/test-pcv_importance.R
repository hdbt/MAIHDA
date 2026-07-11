# pcv_importance(): order-invariant PCV attribution (Shapley / dominance),
# issue #63. The acceptance criteria pinned here: the efficiency identity
# (contributions sum to the full-model Total PCV), order invariance of the
# Shapley/dominance results, parity with stepwise_pcv() (same shared
# complete-case sample, same auto-binned dimension reconstruction, same
# engine/family handshakes), Monte-Carlo convergence to the exact values, and
# signed (suppression-capable) contributions.

# Gaussian data with real additive dimension effects plus a covariate; 2 x 3
# strata keep every lmer fit fast.
make_imp_data <- function(n = 240, seed = 7101) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("m", "f"), n, replace = TRUE),
    race   = sample(c("a", "b", "c"), n, replace = TRUE),
    age    = rnorm(n)
  )
  s <- make_strata(d, c("gender", "race"))
  d <- s$data
  d$y <- 1 +
    ifelse(d$gender == "m", 0.8, 0) +
    c(a = 0, b = 0.9, c = -0.4)[d$race] +
    0.3 * d$age +
    rnorm(n, sd = 1)
  d
}

# Reference between-stratum variance on pcv_importance()'s default basis, which is
# now estimation = "fitted" (each fit's own REML variance for a Gaussian lmer fit),
# so the reference is read WITHOUT the ML refit.
fitted_var <- function(m) {
  MAIHDA:::extract_between_variance(m)
}

test_that("exact Shapley contributions satisfy the efficiency identity", {
  d <- make_imp_data()
  imp <- suppressMessages(pcv_importance(d, "y", c("gender", "race", "age")))

  expect_s3_class(imp, "maihda_pcv_importance")
  expect_identical(imp$method, "shapley")
  expect_identical(imp$approx, "exact")
  expect_equal(nrow(imp$importance), 3L)
  # Efficiency: the contributions sum EXACTLY to the full-model Total PCV.
  expect_equal(sum(imp$importance$Contribution), imp$total_pcv,
               tolerance = 1e-10)
  # 2^3 = 8 models fit, including the null.
  expect_identical(imp$n_fits, 8L)

  # The total matches the direct null-vs-full calculation on the same fitted
  # (REML) variance basis used by calculate_pcv()/stepwise_pcv() by default.
  null_m <- fit_maihda(y ~ 1 + (1 | stratum), d)
  full_m <- suppressMessages(fit_maihda(y ~ gender + race + age + (1 | stratum), d))
  v0 <- fitted_var(null_m)
  vf <- fitted_var(full_m)
  expect_equal(imp$total_pcv, (v0 - vf) / v0, tolerance = 1e-8)
  expect_equal(imp$null_variance, v0, tolerance = 1e-8)
  expect_equal(imp$full_variance, vf, tolerance = 1e-8)

  # Shares are the signed contributions over the total.
  expect_equal(imp$importance$Share,
               imp$importance$Contribution / imp$total_pcv, tolerance = 1e-12)
})

test_that("Shapley contributions are invariant to the order of vars", {
  d <- make_imp_data()
  imp1 <- suppressMessages(pcv_importance(d, "y", c("gender", "race", "age")))
  imp2 <- suppressMessages(pcv_importance(d, "y", c("age", "race", "gender")))

  c1 <- setNames(imp1$importance$Contribution, imp1$importance$Variable)
  c2 <- setNames(imp2$importance$Contribution, imp2$importance$Variable)
  # Same subsets are refit from a differently ordered formula, so allow
  # optimizer-level noise only.
  expect_equal(c1[names(c2)], c2, tolerance = 1e-6)
  expect_equal(imp1$total_pcv, imp2$total_pcv, tolerance = 1e-8)

  # The sequential method, by design, is NOT order-invariant on this data
  # (gender and race are correlated with each other through the strata).
  seq1 <- suppressMessages(
    pcv_importance(d, "y", c("gender", "race"), method = "sequential"))
  seq2 <- suppressMessages(
    pcv_importance(d, "y", c("race", "gender"), method = "sequential"))
  s1 <- setNames(seq1$importance$Contribution, seq1$importance$Variable)
  s2 <- setNames(seq2$importance$Contribution, seq2$importance$Variable)
  expect_false(isTRUE(all.equal(s1[names(s2)], s2, tolerance = 1e-4)))
  # ... but both orders still sum to the same total (efficiency).
  expect_equal(sum(s1), sum(s2), tolerance = 1e-8)
})

test_that("sequential contributions are the increments of stepwise_pcv()'s Total_PCV", {
  d <- make_imp_data()
  vars <- c("gender", "race", "age")
  imp <- suppressMessages(pcv_importance(d, "y", vars, method = "sequential"))
  sw <- suppressMessages(stepwise_pcv(d, "y", vars))

  # Identical fits through the shared setup helper: exact agreement.
  expect_equal(imp$importance$Contribution, diff(c(0, sw$Total_PCV[-1])),
               tolerance = 1e-10)
  expect_equal(imp$total_pcv, sw$Total_PCV[nrow(sw)], tolerance = 1e-10)
  # Only the k path models plus the null are fit.
  expect_identical(imp$n_fits, 4L)
  # The path variances agree with the stepwise table.
  path <- imp$subsets[match(c("gender", "gender + race", "gender + race + age"),
                            imp$subsets$Variables), ]
  expect_equal(path$Variance, sw$Variance[-1], tolerance = 1e-10)
})

test_that("dominance analysis returns general == Shapley plus dominance tables", {
  d <- make_imp_data()
  vars <- c("gender", "race", "age")
  dm <- suppressMessages(pcv_importance(d, "y", vars, method = "dominance"))
  sh <- suppressMessages(pcv_importance(d, "y", vars))

  # General dominance coincides with the Shapley values (same cached fits, so
  # the agreement is exact up to float regrouping).
  expect_equal(dm$importance$Contribution, sh$importance$Contribution,
               tolerance = 1e-12)
  expect_equal(sum(dm$importance$Contribution), dm$total_pcv, tolerance = 1e-10)

  # Conditional dominance: k x k matrix, one column per adjustment-set size,
  # whose row means are the general dominance weights.
  expect_identical(dim(dm$conditional), c(3L, 3L))
  expect_identical(rownames(dm$conditional), vars)
  expect_equal(unname(rowMeans(dm$conditional)), dm$importance$Contribution,
               tolerance = 1e-12)

  # Complete dominance: logical with NA diagonal; race carries the largest
  # simulated effect and completely dominates the weak covariate age.
  expect_identical(dim(dm$complete_dominance), c(3L, 3L))
  expect_true(all(is.na(diag(dm$complete_dominance))))
  expect_true(dm$complete_dominance["race", "age"])
  expect_false(dm$complete_dominance["age", "race"])
})

test_that("Monte-Carlo approximation converges to the exact Shapley values", {
  d <- make_imp_data()
  vars <- c("gender", "race", "age")
  exact <- suppressMessages(pcv_importance(d, "y", vars))
  set.seed(4242)
  mc <- suppressMessages(
    pcv_importance(d, "y", vars, approx = "montecarlo", n_perm = 400))

  expect_identical(mc$approx, "montecarlo")
  # The permutation count is reported, and per-variable MC standard errors
  # are attached.
  expect_identical(mc$n_perm, 400L)
  expect_true("MC_SE" %in% names(mc$importance))
  expect_true(all(is.finite(mc$importance$MC_SE)))

  # Convergence to the exact values (observed MC SE here is ~0.002).
  expect_equal(mc$importance$Contribution, exact$importance$Contribution,
               tolerance = 0.02)
  # Every permutation's marginals telescope to v(N), so the efficiency
  # identity holds EXACTLY for any n_perm -- not just in the limit.
  expect_equal(sum(mc$importance$Contribution), mc$total_pcv,
               tolerance = 1e-10)
  # All subset fits are cached: no more than 2^3 = 8 distinct models.
  expect_lte(mc$n_fits, 8L)
})

test_that("pcv_importance works for a binary outcome and auto-detects the family", {
  d <- make_imp_data(n = 300)
  set.seed(11)
  d$yb <- rbinom(nrow(d), 1,
                 plogis(ifelse(d$gender == "m", 0.9, 0) +
                          c(a = 0, b = 1, c = -0.5)[d$race] - 0.5))

  # Default family must not silently fit a Gaussian attribution on 0/1.
  expect_warning(
    imp <- suppressMessages(pcv_importance(d, "yb", c("gender", "race"))),
    "binary", ignore.case = TRUE
  )
  expect_s3_class(imp, "maihda_pcv_importance")
  expect_match(imp$family, "binomial")
  expect_equal(sum(imp$importance$Contribution), imp$total_pcv,
               tolerance = 1e-10)

  # Explicit binomial family: no auto-detect warning.
  w <- testthat::capture_warnings(suppressMessages(
    pcv_importance(d, "yb", c("gender", "race"), family = "binomial")))
  expect_false(any(grepl("binary", w, ignore.case = TRUE)))
})

test_that("pcv_importance shares one complete-case sample across every subset fit", {
  d <- make_imp_data()
  d$age[seq_len(40)] <- NA_real_

  imp <- suppressMessages(pcv_importance(d, "y", c("gender", "age")))
  complete_d <- d[stats::complete.cases(d[, c("y", "stratum", "gender", "age")]), ]
  expect_identical(imp$n_obs, nrow(complete_d))

  # Every subset -- including gender-only, which has no missing values of its
  # own -- is fit on the SHARED filtered sample.
  gender_m <- fit_maihda(y ~ gender + (1 | stratum), complete_d)
  expect_equal(imp$subsets$Variance[imp$subsets$Variables == "gender"],
               fitted_var(gender_m), tolerance = 1e-8)
})

test_that("pcv_importance reconstructs auto-binned dimensions like stepwise_pcv", {
  set.seed(2026)
  n <- 240
  d <- data.frame(
    gender = sample(c("m", "f"), n, replace = TRUE),
    income = rnorm(n, 50, 12)  # many-valued numeric -> auto-binned to tertiles
  )
  s <- make_strata(d, c("gender", "income"))
  d <- s$data
  d$y <- 1 + ifelse(d$gender == "m", 0.7, 0) +
    0.04 * d$income + rnorm(n, sd = 1)

  expect_true("income" %in% names(attr(d, "strata_autobin_info")))

  imp <- suppressMessages(
    pcv_importance(d, "y", c("gender", "income"), method = "sequential"))
  sw <- suppressMessages(stepwise_pcv(d, "y", c("gender", "income")))

  # Identical variances step-for-step proves income entered as the SAME
  # reconstructed tertile factor, not a raw linear term.
  expect_equal(imp$importance$Contribution, diff(c(0, sw$Total_PCV[-1])),
               tolerance = 1e-10)
  path <- imp$subsets[match(c("gender", "gender + income"),
                            imp$subsets$Variables), ]
  expect_equal(path$Variance, sw$Variance[-1], tolerance = 1e-10)

  # And the Shapley attribution runs on the same reconstructed terms.
  sh <- suppressMessages(pcv_importance(d, "y", c("gender", "income")))
  expect_equal(sum(sh$importance$Contribution), sh$total_pcv, tolerance = 1e-10)
})

test_that("bootstrap intervals cover each contribution and report the draw count", {
  d <- make_imp_data()
  set.seed(515)
  bs <- suppressMessages(
    pcv_importance(d, "y", c("gender", "race"), bootstrap = TRUE, n_boot = 40))

  expect_true(bs$bootstrap)
  expect_true(all(c("CI_lower", "CI_upper") %in% names(bs$importance)))
  expect_true(all(is.finite(bs$importance$CI_lower)))
  expect_true(all(bs$importance$CI_lower <= bs$importance$CI_upper))
  expect_identical(bs$conf_level, 0.95)
  expect_true(is.numeric(bs$n_boot_ok) && bs$n_boot_ok >= 10)
})

test_that("pcv_importance rejects invalid inputs and infeasible requests", {
  d <- make_imp_data()

  # No stratum column.
  expect_error(pcv_importance(data.frame(y = 1:5, x = 1:5), "y", "x"),
               "stratum", fixed = TRUE)
  # Unknown variable.
  expect_error(suppressMessages(pcv_importance(d, "y", c("gender", "nope"))),
               "Variables not found in data")
  # Duplicated variables.
  expect_error(pcv_importance(d, "y", c("gender", "gender")),
               "duplicated")
  # Empty vars.
  expect_error(pcv_importance(d, "y", character(0)),
               "at least one predictor")

  # Dominance needs every subset model: no Monte-Carlo route.
  expect_error(
    pcv_importance(d, "y", c("gender", "race"), method = "dominance",
                   approx = "montecarlo"),
    "dominance analysis needs every subset model")

  # Exact attribution past 15 variables is rejected up front (before any fit).
  d16 <- d
  vars16 <- sprintf("v%02d", seq_len(16))
  for (v in vars16) d16[[v]] <- rnorm(nrow(d16))
  expect_error(
    pcv_importance(d16, "y", vars16, approx = "exact"),
    "not feasible")

  # Invalid n_perm.
  expect_error(
    pcv_importance(d, "y", c("gender", "race"), approx = "montecarlo",
                   n_perm = 3),
    "'n_perm' must be a single whole number between 10 and 1e6", fixed = TRUE)

  # Bootstrap is lme4-only: sampling weights (-> wemix) and the ordinal
  # engine are rejected before any model is fit.
  d$w <- runif(nrow(d), 0.5, 2)
  expect_error(
    suppressMessages(pcv_importance(d, "y", c("gender", "race"),
                                    sampling_weights = "w", bootstrap = TRUE)),
    "only available for\\s+lme4|only available for lme4")
  expect_error(
    suppressMessages(pcv_importance(d, "y", c("gender", "race"),
                                    family = "ordinal", bootstrap = TRUE)),
    "lme4")

  # Bootstrap does not combine with the Monte-Carlo approximation.
  expect_error(
    pcv_importance(d, "y", c("gender", "race"), approx = "montecarlo",
                   bootstrap = TRUE),
    "exact")
})

test_that("print and plot methods work for every attribution flavour", {
  d <- make_imp_data()
  imp <- suppressMessages(pcv_importance(d, "y", c("gender", "race", "age")))
  dm <- suppressMessages(
    pcv_importance(d, "y", c("gender", "race"), method = "dominance"))
  sq <- suppressMessages(
    pcv_importance(d, "y", c("gender", "race"), method = "sequential"))
  set.seed(3)
  mc <- suppressMessages(
    pcv_importance(d, "y", c("gender", "race"), approx = "montecarlo",
                   n_perm = 60))

  expect_output(print(imp), "Shapley values \\(exact\\)")
  expect_output(print(imp), "Total PCV")
  expect_output(print(imp), "efficiency")
  expect_output(print(dm), "Conditional dominance")
  expect_output(print(dm), "Complete dominance")
  expect_output(print(sq), "Order-dependent")
  expect_output(print(mc), "Monte-Carlo, 60 permutations")
  expect_output(print(mc), "MC_SE")

  expect_s3_class(plot(imp), "ggplot")
  expect_s3_class(plot(dm), "ggplot")
  expect_s3_class(plot(sq), "ggplot")

  # Signed contributions: this data yields a (small) negative age
  # contribution; the suppression note is printed and nothing is normalised
  # to non-negative.
  if (any(imp$importance$Contribution < 0)) {
    expect_output(print(imp), "suppression")
  }
})

test_that("a single-variable attribution reduces to the total PCV", {
  d <- make_imp_data()
  imp <- suppressMessages(pcv_importance(d, "y", "gender"))
  expect_equal(imp$importance$Contribution, imp$total_pcv, tolerance = 1e-10)
  expect_identical(imp$n_fits, 2L)
})

test_that("pcv_importance(context = ) attributes the PCV net of the context", {
  set.seed(881)
  n <- 300
  d <- data.frame(
    gender = sample(c("m", "f"), n, replace = TRUE),
    race   = sample(c("a", "b", "c"), n, replace = TRUE),
    site   = sample(sprintf("s%02d", 1:12), n, replace = TRUE)
  )
  s <- make_strata(d, c("gender", "race"))
  d <- s$data
  d$y <- 1 + ifelse(d$gender == "m", 0.7, 0) +
    c(a = 0, b = 0.8, c = -0.3)[d$race] +
    rnorm(12, sd = 0.6)[as.integer(factor(d$site))] +
    rnorm(n, sd = 1)

  imp <- suppressWarnings(suppressMessages(
    pcv_importance(d, "y", c("gender", "race"), context = "site")))
  expect_identical(imp$context, "site")
  expect_true(all(is.finite(imp$importance$Contribution)))
  expect_equal(sum(imp$importance$Contribution), imp$total_pcv,
               tolerance = 1e-10)

  # The full-model between-stratum variance is read NET OF the (1 | site)
  # intercept held in every subset fit -- it must match a direct contextual
  # fit, not the context-free one.
  m_ctx <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ gender + race + (1 | stratum), data = d, context = "site")))
  expect_equal(imp$full_variance, fitted_var(m_ctx), tolerance = 1e-6)

  # The wemix/ordinal engines fit no crossed random effects: rejected with
  # the same guard as stepwise_pcv().
  d$w <- runif(n, 0.5, 2)
  expect_error(
    suppressMessages(pcv_importance(d, "y", c("gender", "race"),
                                    context = "site", sampling_weights = "w")),
    "does not support 'context'")
})

test_that("pcv_importance rejects an effectively singular null model", {
  # Audit follow-up: the strict null_variance <= 0 guard missed a boundary-level
  # positive null (e.g. 1.56e-17), and v_of() then divided every contribution by
  # that degenerate denominator (a reported Total PCV around -2.77e14). The lme4
  # relative singularity guard -- applied before the <= 0 check -- now rejects it,
  # matching calculate_pcv() and stepwise_pcv().
  skip_on_cran()
  set.seed(4)
  n <- 320
  d <- data.frame(
    g = sample(c("F", "M"), n, replace = TRUE),
    r = sample(c("A", "B"), n, replace = TRUE),
    x = rnorm(n)
  )
  sk <- interaction(d$g, d$r, drop = TRUE)
  d$y <- 0.3 * d$x + rnorm(nlevels(sk), sd = 0.12)[sk] + rnorm(n)
  st <- make_strata(d, c("g", "r"))
  d$stratum <- st$data$stratum

  expect_error(
    suppressWarnings(suppressMessages(pcv_importance(d, "y", "x"))),
    "zero boundary"
  )
})
