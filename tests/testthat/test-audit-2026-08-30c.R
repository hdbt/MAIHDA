# Audit 2026-08-30 (third pass of the day): every MAIHDA parametric bootstrap
# called stats::simulate() on the lme4 fit directly. lme4's simulate.merMod()
# draws the Gaussian conditional residual as sigma * rnorm(n) for EVERY row
# (lme4:::.simulateFun, the isLMM() branch), ignoring the prior weights -- even
# though an lmer fitted with weights w_i has conditional residual variance
# sigma^2 / w_i, which is what maihda_gaussian_residual_variance_lme4() encodes
# for the point VPC and what lme4::refit() assumes when it refits the draw.
# Simulating equal-variance noise and refitting under 1/w_i weights inflated the
# bootstrap residual variance by about mean(w), pushing every VPC/PCV interval
# off its own point estimate. Simulation is now centralised in
# maihda_simulate_lme4(), which adds rnorm(sd = sigma / sqrt(w)) itself.
#
# Note a scope correction to the report: only ALL-UNIT weights are safe. A
# CONSTANT weight c != 1 was biased by exactly c, not spared.

# Small weighted Gaussian MAIHDA fit; `w` is the prior-weight vector.
audit_0830c_fit <- function(w) {
  set.seed(404)
  n <- 240L
  d <- data.frame(
    a = factor(rep(1:3, length.out = n)),
    b = factor(rep(1:4, each = 3, length.out = n)),
    x = stats::rnorm(n)
  )
  d$stratum <- interaction(d$a, d$b, drop = TRUE)
  # 1/sqrt(0) is Inf and rnorm(sd = Inf) is NA, which would make lmer DROP the
  # zero-weight row instead of keeping it -- exactly the case being tested.
  # Identical to 1/sqrt(w) for every positive weight.
  d$y <- 2 + 0.5 * d$x + stats::rnorm(nlevels(d$stratum))[d$stratum] +
    stats::rnorm(n, sd = 1 / sqrt(ifelse(w > 0, w, 1)))
  d$w <- w
  suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (1 | stratum), data = d, weights = w)
  ))
}

test_that("only a weighted Gaussian LMM is diverted off lme4's simulate()", {
  n <- 240L
  ones <- rep(1, n)

  # Unweighted and all-weights-1 LMMs stay on lme4's own path.
  expect_null(MAIHDA:::maihda_lme4_simulation_weights(audit_0830c_fit(ones)))

  # Varying weights, and a CONSTANT non-unit weight, are both corrected.
  w_vary <- rep(c(1, 25), length.out = n)
  expect_equal(MAIHDA:::maihda_lme4_simulation_weights(audit_0830c_fit(w_vary)),
               w_vary)
  w_const <- rep(9, n)
  expect_equal(MAIHDA:::maihda_lme4_simulation_weights(audit_0830c_fit(w_const)),
               w_const)

  # A GLMM is left alone: lme4:::.simulateFun() already hands `wts` to the
  # family's simulation function, so correcting it here would double-count.
  set.seed(9)
  db <- data.frame(g = factor(rep(1:8, each = 20)), x = stats::rnorm(160))
  db$succ <- stats::rbinom(160, size = 5, prob = 0.4)
  db$fail <- 5L - db$succ
  gm <- suppressWarnings(suppressMessages(
    lme4::glmer(cbind(succ, fail) ~ x + (1 | g), data = db, family = stats::binomial())
  ))
  expect_null(MAIHDA:::maihda_lme4_simulation_weights(gm))
})

test_that("lme4 still honours cond.sim = FALSE (the fix rests on it)", {
  # maihda_simulate_lme4() asks lme4 for a residual-free predictor and adds its own
  # noise. cond.sim reaches lme4:::.simulateFun() through simulate.merMod()'s `...`,
  # so if a future lme4 drops or renames it the argument would be silently swallowed
  # and the residual counted TWICE. Pin the invariant: with no conditional noise the
  # draw is constant within a stratum (it is fixed effects + one random intercept),
  # whereas the ordinary draw is not.
  m <- audit_0830c_fit(rep(c(1, 25), length.out = 240L))
  g <- lme4::getME(m, "flist")[[1]]

  set.seed(41)
  quiet <- stats::simulate(m, nsim = 1, cond.sim = FALSE)[[1]]
  noisy <- stats::simulate(m, nsim = 1)[[1]]
  resid_free <- quiet - stats::predict(m, re.form = NA)
  expect_lt(max(tapply(resid_free, g, stats::sd), na.rm = TRUE), 1e-8)
  expect_gt(max(tapply(noisy - stats::predict(m, re.form = NA), g, stats::sd),
                na.rm = TRUE), 1e-3)
})

test_that("a single unusable weight does not disable the whole correction", {
  # lmer accepts a zero weight, and the first cut of this fix returned NULL for the
  # entire fit when it saw one -- silently restoring the equal-variance bug for
  # every other row. Zero-weight rows draw at lme4's sigma (they contribute nothing
  # to the weighted likelihood); the rest keep sigma / sqrt(w_i).
  w <- rep(c(1, 25), length.out = 240L)
  w[1] <- 0
  m <- audit_0830c_fit(w)

  expect_equal(MAIHDA:::maihda_lme4_simulation_weights(m), w)

  sig <- stats::sigma(m)
  set.seed(51); got <- MAIHDA:::maihda_simulate_lme4(m, nsim = 3)
  set.seed(51)
  want <- stats::simulate(m, nsim = 3, cond.sim = FALSE)
  sd_i <- ifelse(w > 0, sig / sqrt(ifelse(w > 0, w, 1)), sig)
  for (i in seq_len(3)) {
    want[[i]] <- want[[i]] + stats::rnorm(length(w), mean = 0, sd = sd_i)
  }
  expect_equal(got, want)
  expect_true(all(vapply(got, function(z) all(is.finite(z)), logical(1))))

  # A fit whose usable weights are all 1 stays on lme4's path even when a zero is
  # present: the residual sd is sigma on every row either way, so lme4 is already
  # right and the draws stay bit-identical. (The "no usable weight at all" guard is
  # defensive only -- an all-zero-weight lmer cannot be fitted at all, it fails in
  # the deviance function with a non-positive-definite downdate.)
  w1 <- rep(1, 240L); w1[1] <- 0
  m1 <- audit_0830c_fit(w1)
  expect_null(MAIHDA:::maihda_lme4_simulation_weights(m1))
  set.seed(61); via_helper <- MAIHDA:::maihda_simulate_lme4(m1, nsim = 2)
  set.seed(61); via_lme4 <- stats::simulate(m1, nsim = 2)
  expect_identical(via_helper, via_lme4)
})

test_that("unweighted simulation is bit-identical to stats::simulate()", {
  # The fix must not perturb the overwhelmingly common unweighted path, RNG
  # stream included.
  m <- audit_0830c_fit(rep(1, 240L))
  set.seed(7); via_helper <- MAIHDA:::maihda_simulate_lme4(m, nsim = 3)
  set.seed(7); via_lme4 <- stats::simulate(m, nsim = 3)
  expect_identical(via_helper, via_lme4)
})

test_that("weighted draws carry residual variance sigma^2 / w_i", {
  w <- rep(c(1, 25), length.out = 240L)
  m <- audit_0830c_fit(w)
  sig <- stats::sigma(m)

  # Exact, deterministic reconstruction of the intended draw: the residual-free
  # predictor plus rnorm(sd = sigma / sqrt(w)). Under the old code the added
  # noise was sigma * rnorm(n) for every row, so this identity fails.
  set.seed(31); got <- MAIHDA:::maihda_simulate_lme4(m, nsim = 4)
  set.seed(31)
  want <- stats::simulate(m, nsim = 4, cond.sim = FALSE)
  for (i in seq_len(4)) {
    want[[i]] <- want[[i]] + stats::rnorm(length(w), mean = 0, sd = sig / sqrt(w))
  }
  expect_equal(got, want)

  # ... and the realised noise really is heteroscedastic by weight. The old
  # behaviour put this ratio at 1; lme4's semantics require 25.
  set.seed(32); full <- MAIHDA:::maihda_simulate_lme4(m, nsim = 400)
  set.seed(32); base <- stats::simulate(m, nsim = 400, cond.sim = FALSE)
  noise <- as.matrix(full) - as.matrix(base)
  v <- tapply(seq_along(w), w, function(i) mean(apply(noise[i, , drop = FALSE], 1, stats::var)))
  expect_equal(unname(v[["1"]] / v[["25"]]), 25, tolerance = 0.15)
  expect_equal(unname(v[["1"]]), sig^2, tolerance = 0.15)
  expect_equal(unname(v[["25"]]), sig^2 / 25, tolerance = 0.15)
})

test_that("refitting a draw recovers sigma^2 instead of inflating it by mean(w)", {
  skip_on_cran()
  # The statistically material consequence: refit() keeps the 1/w_i weights, so
  # equal-variance draws came back with sigma^2 scaled up by about mean(w)
  # (~13 for 1-vs-25 weights, and exactly c for a constant weight c).
  for (w in list(rep(c(1, 25), length.out = 240L), rep(9, 240L))) {
    m <- audit_0830c_fit(w)
    set.seed(5)
    sim <- MAIHDA:::maihda_simulate_lme4(m, nsim = 12)
    ratio <- mean(vapply(seq_len(12), function(i) {
      stats::sigma(suppressWarnings(suppressMessages(
        lme4::refit(m, newresp = sim[[i]]))))^2
    }, numeric(1))) / stats::sigma(m)^2
    # Old code landed near mean(w); the corrected draw returns to ~1.
    expect_equal(ratio, 1, tolerance = 0.25)
    expect_lt(ratio, 0.5 * mean(w))
  }
})

test_that("the weighted bootstrap VPC interval surrounds its point estimate", {
  skip_on_cran()
  w <- rep(c(1, 25), length.out = 240L)
  m <- audit_0830c_fit(w)
  point <- MAIHDA:::maihda_stratum_variance_lme4(m) /
    (MAIHDA:::maihda_stratum_variance_lme4(m) +
       MAIHDA:::maihda_gaussian_residual_variance_lme4(m))

  set.seed(13)
  sim <- MAIHDA:::maihda_simulate_lme4(m, nsim = 60)
  vpc <- vapply(seq_len(60), function(i) {
    bm <- suppressWarnings(suppressMessages(lme4::refit(m, newresp = sim[[i]])))
    tau <- MAIHDA:::maihda_stratum_variance_lme4(bm)
    tau / (tau + MAIHDA:::maihda_gaussian_residual_variance_lme4(bm))
  }, numeric(1))

  ci <- unname(stats::quantile(vpc, c(0.025, 0.975), na.rm = TRUE))
  # The defect drove the whole interval below the point estimate.
  expect_gte(point, ci[1])
  expect_lte(point, ci[2])
})
