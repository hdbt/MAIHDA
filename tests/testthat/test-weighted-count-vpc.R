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
  # (The per-row term is computed at the same marginal means the package uses,
  # so the contrast isolates the weighting alone.)
  mu_wt <- MAIHDA:::maihda_count_marginal_mu_lme4(m_wt$model)
  expect_false(isTRUE(all.equal(rv_wt, mean(log1p(1 / mu_wt)), tolerance = 1e-3)))
})

# ---- brms sampling-weighted count VPC reads the .maihda_sw column ------------

test_that("brms Poisson count VPC averages the latent variance by the sampling weights", {
  # Regression guard: brms exposes no weights.brmsfit, so maihda_fit_prior_weights()
  # saw nothing through stats::weights() and the population count-family VPC fell back
  # to an UNWEIGHTED mean of the per-row latent variance, despite the design-weighted
  # claim. The sampling weights live in the reserved .maihda_sw data column, so they
  # must be read from there and folded into the average.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  set.seed(515)
  n <- 300
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  u <- stats::rnorm(nlevels(sk), sd = 0.6)[sk]
  d$y <- stats::rpois(n, lambda = exp(0.3 + 0.5 * (d$gender == "M") + u))
  # Weights strongly tied to the linear predictor so the weighted and unweighted
  # means of the per-row latent variance genuinely differ.
  d$w <- ifelse(d$gender == "M", 3, 0.4) * stats::runif(n, 0.8, 1.2)

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ gender + (1 | gender:race), data = d,
               engine = "brms", family = "poisson", sampling_weights = "w",
               chains = 1, iter = 300, refresh = 0)
  ))

  # The raw brmsfit really exposes the normalized sampling weights for the helper.
  w_read <- maihda_fit_prior_weights(m$model)
  expect_false(is.null(w_read))
  expect_equal(length(w_read), maihda_nobs(m$model))

  # Per-row latent variance at the same marginal expected counts the package
  # uses, so the weighted-vs-unweighted contrast isolates the weighting alone.
  draws <- as.data.frame(m$model)
  mu <- MAIHDA:::maihda_count_marginal_mu_brms(m$model, draws)
  rv_unwt <- mean(log1p(1 / mu))
  rv_wt   <- maihda_weighted_obs_mean(log1p(1 / mu), w_read)

  rv_path <- maihda_residual_variance_draws_brms(m$model, draws)
  expect_equal(rv_path[1], rv_wt)                          # uses the weighted average
  expect_false(isTRUE(all.equal(rv_path[1], rv_unwt)))     # not the old unweighted one
})
