# Weighted count-model VPC/PCV: the Poisson / negative-binomial latent-scale
# level-1 variance must average its per-row terms by the model's prior weights,
# so a frequency-weighted fit reproduces the plain mean over the equivalent
# duplicated-row data. Regression guard for the bug where the average ignored the
# weights (a targeted weighted Poisson fit gave a 154% larger residual variance).

# ---- pure averaging helper (no fit required) ---------------------------------

test_that("maihda_weighted_obs_mean reproduces the duplicated-row mean", {
  set.seed(11)
  v <- rnorm(40)
  w <- sample(1:4, 40, replace = TRUE)

  # The defining equivalence: a frequency-weighted mean equals the plain mean over
  # the data expanded by those integer weights.
  expect_equal(maihda_weighted_obs_mean(v, w), mean(rep(v, w)))

  # Falls back to the unweighted mean when weights are absent, all 1, or the
  # wrong length (so the unweighted VPC path is byte-for-byte unchanged).
  expect_equal(maihda_weighted_obs_mean(v, NULL), mean(v))
  expect_equal(maihda_weighted_obs_mean(v, rep(1, 40)), mean(v))
  expect_equal(maihda_weighted_obs_mean(v, c(1, 2)), mean(v))

  # Non-finite / non-positive weights are dropped rather than poisoning the mean.
  expect_equal(maihda_weighted_obs_mean(c(1, 2, 3, 4), c(2, NA, 0, 1)),
               (2 * 1 + 1 * 4) / (2 + 1))
})

# ---- weighted-vs-expanded equivalence on a real glmer fit --------------------

test_that("weighted Poisson residual variance matches the duplicated-row fit", {
  skip_on_cran()

  set.seed(303)
  n <- 300
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("A", "B", "C"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$gender, d$race, drop = TRUE)
  u <- rnorm(nlevels(stratum), sd = 0.5)[stratum]
  d$y <- rpois(n, lambda = exp(0.5 + 0.3 * (d$gender == "M") + u))
  # Integer frequency weights give the exact duplicated-row equivalence.
  d$w <- sample(1:3, n, replace = TRUE)

  expanded <- d[rep(seq_len(n), d$w), c("gender", "race", "y")]
  rownames(expanded) <- NULL

  fml <- y ~ gender + (1 | gender:race)
  m_wt  <- suppressWarnings(fit_maihda(fml, data = d, family = "poisson", weights = w))
  m_exp <- suppressWarnings(fit_maihda(fml, data = expanded, family = "poisson"))

  # The weighted fit really carries lme4 prior weights for the helper to read.
  expect_equal(as.numeric(stats::weights(m_wt$model, type = "prior")), as.numeric(d$w))

  rv_wt  <- maihda_residual_variance_lme4(m_wt$model)
  rv_exp <- maihda_residual_variance_lme4(m_exp$model)

  # Frequency-weighted average == plain mean over the equivalent expanded rows.
  expect_equal(rv_wt, rv_exp, tolerance = 1e-3)

  # The pre-fix behaviour (ignore the weights, plain mean over the un-expanded
  # rows) is materially different -- the bug this regression guards against. The
  # negative-binomial path averages through the same helper, so it is covered too.
  mu_wt <- pmax(as.numeric(stats::fitted(m_wt$model)), .Machine$double.eps)
  expect_false(isTRUE(all.equal(rv_wt, mean(log1p(1 / mu_wt)), tolerance = 1e-3)))
})
