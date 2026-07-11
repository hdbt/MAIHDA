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

test_that("maihda_count_resid_var_from_linpred folds per-draw marginal means", {
  # Pure, Stan-free core of the per-draw count VPC residual variance: given the
  # fixed-part linear-predictor DRAWS (ndraws x nobs) and each draw's total
  # random-intercept variance, it forms lambda_{d,i} = exp(eta + v_d/2) per draw and
  # returns the (weighted) mean of log1p(1/lambda [+ 1/shape]) over observations.
  eta <- matrix(c(0, 0.5, 1.0,
                  0.2, 0.2, 0.2), nrow = 2, byrow = TRUE)   # 2 draws x 3 obs
  vtot <- c(0.4, 1.0)
  mu <- exp(sweep(eta, 1L, vtot / 2, "+"))

  # Poisson (no extra term), unweighted.
  expect_equal(MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot),
               rowMeans(log1p(1 / mu)))

  # Negative-binomial: the per-draw 1/shape term is added to 1/lambda per draw.
  shape <- c(5, 2)
  expect_equal(
    MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot, extra = 1 / shape),
    rowMeans(log1p(sweep(1 / mu, 1L, 1 / shape, "+")))
  )

  # Prior weights average the per-row terms within each draw.
  w <- c(3, 1, 1)
  expect_equal(
    MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot, w = w),
    vapply(seq_len(2), function(d) sum(w * log1p(1 / mu[d, ])) / sum(w), numeric(1))
  )

  # The marginal mean genuinely varies across draws (fixed part + RE variance),
  # so the residual variance is NOT constant -- the defect this fix removes.
  expect_gt(stats::sd(MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot)), 0)

  # Length guards.
  expect_error(MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot[1]),
               "one value per posterior draw")
})

test_that("per-draw count residual variance gates on an intercept-only structure", {
  # A random-slope (longitudinal-style) structure must fall back to the plug-in
  # path (NULL here) -- the fast per-draw shortcut uses the total random-INTERCEPT
  # variance and does not carry the row-varying slope design, so it must not be
  # applied to a growth-model count fit. Gated before any brms call, so Stan-free.
  slope_model <- list(formula = y ~ x + (x | id))
  draws <- data.frame(sd_id__Intercept = c(1, 2), sd_id__x = c(0.3, 0.4))
  expect_null(
    MAIHDA:::maihda_count_residual_variance_draws_brms_perdraw(slope_model, draws))

  # Intercept-only but with no sd_* draws to build the RE variance from -> NULL,
  # so the caller's posterior-mean plug-in still applies.
  ic_model <- list(formula = y ~ (1 | stratum))
  expect_null(
    MAIHDA:::maihda_count_residual_variance_draws_brms_perdraw(
      ic_model, data.frame(b_Intercept = c(0, 1))))
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

  # Per-draw latent variance now uses the per-draw marginal mean lambda_{1,i} =
  # exp(eta_{1,i} + v_1/2) (fixed part + total RE variance for draw 1), so the
  # weighted-vs-unweighted contrast is reconstructed against THAT draw, not the
  # posterior-mean plug-in the old path used.
  draws <- as.data.frame(m$model)
  eta_link <- brms::posterior_linpred(m$model, re_formula = NA)   # ndraws x nobs
  sd_cols <- grep("^sd_", names(draws), value = TRUE)
  vtot <- Reduce(`+`, lapply(sd_cols, function(cn) draws[[cn]]^2))
  mu1 <- exp(eta_link[1, ] + vtot[1] / 2)
  mu1[mu1 < .Machine$double.eps] <- .Machine$double.eps
  rv_unwt_1 <- mean(log1p(1 / mu1))
  rv_wt_1   <- maihda_weighted_obs_mean(log1p(1 / mu1), w_read)

  rv_path <- maihda_residual_variance_draws_brms(m$model, draws)
  expect_length(rv_path, nrow(draws))
  expect_equal(rv_path[1], rv_wt_1)                        # per-draw + weighted
  expect_false(isTRUE(all.equal(rv_path[1], rv_unwt_1)))   # weighting still matters
  # The per-draw marginal means propagate fixed-effect + RE-variance uncertainty,
  # so the residual variance varies across draws (the old plug-in was constant).
  expect_gt(stats::sd(rv_path), 0)

  # Public paths over the same fit -- regression guards for the draws-only
  # posterior_linpred (real brms ignores summary = TRUE): summary() of a brms
  # count fit, link-scale prediction, and the stratum predictions behind
  # plot(type = "predicted") / maihda_table() must all work against real brms.
  s <- suppressWarnings(summary(m))
  expect_true(is.finite(s$vpc$estimate) && s$vpc$estimate > 0 && s$vpc$estimate < 1)
  eta <- predict_maihda(m, type = "individual", scale = "link")
  expect_length(eta, maihda_nobs(m$model))
  expect_true(all(is.finite(eta)))
  sp <- MAIHDA:::maihda_stratum_predictions_brms(m, s, scale = "response")
  expect_true(all(is.finite(sp$predicted_row)))
})
