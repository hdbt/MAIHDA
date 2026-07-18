# Regression tests for the 2026-07-28 audit findings.
#
#   1 [P2] stepwise_pcv() gated Step_PCV on `prev_var > 0` alone. The null
#          denominator is hard-stopped at the singularity boundary, but a STEP
#          model at the boundary (already recorded in boundary_steps) still became
#          the NEXT step's denominator. lme4 returns an exact zero there, so the
#          `> 0` test caught it by accident; the wemix and ordinal optimizers
#          return a tiny POSITIVE value that sails through it. A reproduced wemix
#          series ran 3.756 -> 1.43e-21 -> 3.97e-16 and reported a Step_PCV of
#          -276564 (-27656406%) -- purely the ratio of two optimizer artefacts.
#          Step_PCV is now NA whenever the PRECEDING model is at the boundary,
#          the affected steps are listed in an "undefined_step_pcv" attribute, and
#          print() explains it. Total_PCV divides by the (guaranteed healthy) null
#          variance and is unchanged.
#   2 [P3] The contextual branch of maihda_vpc_response() built its numerator from
#          stats::var() (denominator n-1) but its total from pooled population
#          moments E[p^2] - E[p]^2 (denominator n), so the two disagreed on the
#          between piece and the VPC was inflated by exactly n/(n-1) -- 1.01% at
#          the allowed minimum n_sim = 100. The total is now assembled by the law
#          of total variance from the SAME sample variance, which also makes the
#          contextual branch exactly continuous with the var_other == 0 branch.

# ---- Finding 1: a boundary step must not become the next step's denominator ----

# A healthy series: every model keeps a strictly positive, non-boundary between-
# stratum variance, and x1/x2 are individual-level covariates that do not consume
# it. Boundary status is then injected below, so the test exercises the GATE rather
# than an optimizer's seed-fragile boundary behaviour.
maihda_a728_series_df <- function(seed = 21) {
  set.seed(seed)
  G <- 40; npg <- 25
  base <- data.frame(gender = factor(rep(c("M", "F"), length.out = G)),
                     race   = factor(rep(1:20, length.out = G)))
  d <- base[rep(seq_len(G), each = npg), ]
  d <- make_strata(d, c("gender", "race"))$data
  slev <- sort(unique(as.character(d$stratum)))
  u <- stats::setNames(stats::rnorm(length(slev), 0, 1.2),
                       slev)[as.character(d$stratum)]
  d$x1 <- stats::rnorm(nrow(d))
  d$x2 <- stats::rnorm(nrow(d))
  d$y  <- 1 + u + 0.4 * d$x1 + 0.3 * d$x2 + stats::rnorm(nrow(d))
  rownames(d) <- NULL
  d
}

test_that("a healthy stepwise series keeps every Step_PCV (no spurious NA)", {
  d <- maihda_a728_series_df()
  sw <- suppressWarnings(suppressMessages(stepwise_pcv(d, "y", c("x1", "x2"))))
  expect_true(all(is.finite(sw$Step_PCV)))
  expect_true(all(is.finite(sw$Total_PCV)))
  expect_identical(attr(sw, "undefined_step_pcv"), character(0))
  expect_identical(attr(sw, "boundary_steps"), character(0))
  out <- paste(capture.output(print(sw)), collapse = "\n")
  expect_false(grepl("Step_PCV is NA", out, fixed = TRUE))
})

test_that("Step_PCV is NA when the PRECEDING model is at the boundary", {
  d <- maihda_a728_series_df()

  # maihda_pcv_null_at_boundary() is called exactly once per model in
  # stepwise_pcv() -- null first, then one per step. Report the null as healthy
  # (a boundary null is a hard stop) and both step models as at the boundary.
  # Their variances stay strictly positive and far from zero, so `prev_var > 0`
  # is TRUE throughout: only the boundary FLAG can produce the NA.
  calls <- 0L
  fake <- function(model) {
    calls <<- calls + 1L
    calls > 1L
  }
  orig <- MAIHDA:::maihda_pcv_null_at_boundary
  assignInNamespace("maihda_pcv_null_at_boundary", fake, ns = "MAIHDA")
  sw <- tryCatch(
    suppressWarnings(suppressMessages(stepwise_pcv(d, "y", c("x1", "x2")))),
    finally = assignInNamespace("maihda_pcv_null_at_boundary", orig, ns = "MAIHDA"))

  expect_identical(calls, 3L)                    # null + 2 steps
  # The old code would have divided by this: strictly positive, not an exact zero.
  expect_gt(sw$Variance[2], 0.5)

  # Step 1 divides by the null variance, which the hard stop guarantees is healthy.
  expect_true(is.finite(sw$Step_PCV[2]))
  # Step 2 divides by step 1's boundary variance -> undefined.
  expect_true(is.na(sw$Step_PCV[3]))
  # Total_PCV always divides by the null variance, so it survives for both.
  expect_true(all(is.finite(sw$Total_PCV)))

  expect_identical(attr(sw, "boundary_steps"), c("x1", "x2"))
  expect_identical(attr(sw, "undefined_step_pcv"), "x2")

  out <- paste(capture.output(print(sw)), collapse = "\n")
  expect_match(out, "Step_PCV is NA after adding: x2", fixed = TRUE)
  expect_match(out, "Total_PCV is unaffected", fixed = TRUE)
})

test_that("the boundary gate survives a real wemix tiny-positive series", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  # `x` is constant within stratum and drives ALL between-stratum variation, so
  # every model that includes it sits on the boundary; `x2` is individual-level
  # noise that leaves it there. The wemix optimizer reports that boundary as a tiny
  # POSITIVE variance rather than lme4's exact zero.
  set.seed(5)
  G <- 60; npg <- 20
  base <- data.frame(gender = factor(rep(c("M", "F"), length.out = G)),
                     race   = factor(rep(1:10, length.out = G)),
                     edu    = factor(rep(1:3, length.out = G)))
  d <- base[rep(seq_len(G), each = npg), ]
  d <- make_strata(d, c("gender", "race", "edu"))$data
  slev <- sort(unique(as.character(d$stratum)))
  xmap <- stats::setNames(stats::rnorm(length(slev)), slev)
  d$x  <- xmap[as.character(d$stratum)]
  d$x2 <- stats::rnorm(nrow(d))
  d$y  <- 2 * d$x + stats::rnorm(nrow(d), 0, 1)
  d$w  <- stats::runif(nrow(d), 0.5, 2)
  rownames(d) <- NULL

  sw <- suppressWarnings(suppressMessages(
    stepwise_pcv(d, "y", c("x", "x2"), sampling_weights = "w")))

  # Boundary fits are optimizer- and seed-fragile; only assert the gate when this
  # platform's WeMix actually produced the tiny-positive series being guarded.
  if (length(attr(sw, "boundary_steps")) < 2 || !isTRUE(sw$Variance[2] > 0)) {
    skip("WeMix did not produce the tiny-positive boundary series on this platform")
  }

  expect_lt(sw$Variance[2], 1e-8)     # boundary...
  expect_gt(sw$Variance[2], 0)        # ...but strictly positive: `> 0` cannot catch it
  expect_true(is.na(sw$Step_PCV[3]))  # was -276564 (a ratio of two artefacts)
  expect_true(is.finite(sw$Total_PCV[3]))
  expect_identical(attr(sw, "undefined_step_pcv"), "x2")
})

# ---- Finding 2: consistent moments in the contextual response-scale VPC --------

maihda_a728_ctx_fit <- function(seed = 11) {
  set.seed(seed)
  S <- 12; K <- 20; npg <- 6
  gr <- expand.grid(school = factor(seq_len(K)), s = seq_len(S))
  d <- gr[rep(seq_len(nrow(gr)), each = npg), ]
  d$gender <- factor(c("F", "M")[(d$s - 1) %% 2 + 1])
  d$race   <- factor(((d$s - 1) %/% 2) + 1)
  d <- make_strata(d, c("gender", "race"))$data
  u_s <- stats::rnorm(S, 0, 0.8)[d$s]
  u_k <- stats::rnorm(K, 0, 0.7)[as.integer(d$school)]
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.4 + u_s + u_k))
  rownames(d) <- NULL
  suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum) + (1 | school), data = d, family = "binomial")))
}

test_that("the contextual response-scale VPC uses one basis for the between variance", {
  skip_on_cran()
  m <- maihda_a728_ctx_fit()
  n_sim <- 100L                                  # the allowed minimum: worst case
  r <- maihda_vpc_response(m, n_sim = n_sim, seed = 1)
  expect_gt(r$var_other, 0)
  expect_identical(r$inner_method, "gauss-hermite")

  # Replay the exact computation (same seed, same quadrature nodes).
  linkinv <- stats::family(m$model)$linkinv
  set.seed(1)
  u <- stats::rnorm(n_sim, 0, sqrt(r$var_between))
  gh <- MAIHDA:::maihda_gauss_hermite_normal(MAIHDA:::maihda_vpc_response_inner_nodes)
  P  <- linkinv(outer(u, gh$nodes * sqrt(r$var_other), `+`) + r$lp_fixed)
  ms <- as.vector(P %*% gh$weights)                      # E[p | u_i]
  si <- as.vector((P * P) %*% gh$weights)                # E[p^2 | u_i]
  vw <- mean(as.vector((P * (1 - P)) %*% gh$weights))    # E[p(1-p)]

  # Law of total variance: the denominator's between piece IS the numerator.
  ltv <- stats::var(ms) / (stats::var(ms) + mean(pmax(si - ms^2, 0)) + vw)
  # The old form: n-1 numerator over an n population-moment total.
  old <- stats::var(ms) / (max(mean(si) - mean(ms)^2, 0) + vw)

  expect_equal(r$estimate, ltv, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(r$estimate, old, tolerance = 1e-9)))
  expect_gt(old, r$estimate)                     # the mismatch inflated the VPC

  # The whole discrepancy is the n vs n-1 divisor on the between term.
  expect_equal(stats::var(ms) / (mean(ms^2) - mean(ms)^2),
               n_sim / (n_sim - 1), tolerance = 1e-9)
  expect_true(r$estimate >= 0 && r$estimate <= 1)
})

test_that("the contextual VPC formula is continuous with the var_other == 0 branch", {
  skip_on_cran()
  # Same data, no non-stratum random effect: the inner conditional variance is
  # identically zero, so the law-of-total-variance form must collapse EXACTLY onto
  # the stratum-only branch. The population-moment total did not (it differed by
  # the same n/(n-1) factor), leaving a discontinuity at var_other -> 0.
  set.seed(11)
  S <- 12; K <- 20; npg <- 6
  gr <- expand.grid(school = factor(seq_len(K)), s = seq_len(S))
  d <- gr[rep(seq_len(nrow(gr)), each = npg), ]
  d$gender <- factor(c("F", "M")[(d$s - 1) %% 2 + 1])
  d$race   <- factor(((d$s - 1) %/% 2) + 1)
  d <- make_strata(d, c("gender", "race"))$data
  u_s <- stats::rnorm(S, 0, 0.8)[d$s]
  u_k <- stats::rnorm(K, 0, 0.7)[as.integer(d$school)]
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.4 + u_s + u_k))
  rownames(d) <- NULL

  m0 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "binomial")))
  r0 <- maihda_vpc_response(m0, n_sim = 100L, seed = 1)
  expect_identical(r0$var_other, 0)

  set.seed(1)
  u <- stats::rnorm(100L, 0, sqrt(r0$var_between))
  p <- stats::plogis(r0$lp_fixed + u)
  expect_equal(r0$estimate,
               stats::var(p) / (stats::var(p) + mean(p * (1 - p))),
               tolerance = 1e-12)
})
