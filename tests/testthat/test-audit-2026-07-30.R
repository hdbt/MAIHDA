# Regression tests for the 2026-07-30 audit findings.
#
#   1 [High] The count-family level-1 (residual) variance averaged the TRANSFORM
#            over observations, mean_i log1p(1 / lambda_i), while citing Stryhn et
#            al. (2006) and Nakagawa, Johnson & Schielzeth (2017), whose estimator
#            is defined at a single global mean count: log1p(1 / mean_i lambda_i).
#            log1p(1 / x) is convex, so by Jensen's inequality the package returned
#            a systematically LARGER residual variance and a SMALLER VPC whenever
#            the fitted means varied -- on an adjusted Poisson probe 0.4197 against
#            0.1516, carrying VPC 0.531 against 0.758. The two agree exactly when
#            every lambda_i is equal, which is why a null model without an offset
#            (the MAIHDA headline) was and remains correct: its value reproduces
#            insight::get_variance()'s var.distribution to 12 digits. The reduction
#            now happens on the COUNTS, in one shared helper, for both engines,
#            both count families, the per-draw path, and the longitudinal VPC(t).
#   2 [High] maihda_da_aggregated_counts() read ANY integer-valued prior weights
#            with a value above 1 as aggregated binomial trial counts, on the
#            premise that "a genuine Bernoulli fit has unit prior weights". But
#            fit_maihda() documents `weights =` as lme4 PRECISION weights, so an
#            integer precision-weight vector was silently converted into
#            case/control frequency mass: a 600-row Bernoulli fit with weights 1:5
#            reported 727 cases and 1069 controls against the observed 251 / 349,
#            and an AUC of 0.6689 against the observation-level 0.6898 -- while the
#            SAME weights plus 1e-6 took the ordinary path. Aggregation is now
#            inferred structurally: a cbind() matrix response, or R's proportion
#            response + trials-as-weights idiom, neither of which a Bernoulli fit
#            can present.

# ---- Finding 1: the count level-1 variance reduces on the counts -------------

test_that("maihda_count_level1_variance reduces before the transform", {
  # Pure helper, no fit required. lambda is the MEAN count, so the result is the
  # transform of the mean -- not the mean of the transforms.
  mu <- c(0.2, 1, 5, 40)
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu), log1p(1 / mean(mu)))
  expect_false(isTRUE(all.equal(MAIHDA:::maihda_count_level1_variance(mu),
                                mean(log1p(1 / mu)))))

  # Jensen's inequality, in the direction the finding describes: the old
  # observation-averaged form is strictly larger whenever the counts vary...
  expect_gt(mean(log1p(1 / mu)), MAIHDA:::maihda_count_level1_variance(mu))
  # ...and exactly equal when they do not (the null-model case).
  flat <- rep(3.7, 5)
  expect_equal(MAIHDA:::maihda_count_level1_variance(flat), mean(log1p(1 / flat)))

  # The negative-binomial 1/theta term is added to 1/lambda, and vanishes as
  # theta -> Inf (reducing to the Stryhn Poisson form).
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, theta = 2),
               log1p(1 / mean(mu) + 1 / 2))
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, theta = Inf),
               MAIHDA:::maihda_count_level1_variance(mu))

  # Prior weights weight the mean COUNT, and reproduce the plain mean over the
  # equivalent duplicated-row data.
  w <- c(3, 1, 2, 1)
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, w = w),
               log1p(1 / mean(rep(mu, w))))
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, w = rep(1, 4)),
               MAIHDA:::maihda_count_level1_variance(mu))
})

test_that("the per-draw count residual variance uses each draw's mean count", {
  # Stan-free core of the brms per-draw path.
  eta  <- matrix(c(0, 0.5, 1.0,
                   0.2, 0.2, 0.2), nrow = 2, byrow = TRUE)
  vtot <- c(0.4, 1.0)
  mu   <- exp(sweep(eta, 1L, vtot / 2, "+"))

  expect_equal(MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot),
               log1p(1 / rowMeans(mu)))
  # Draw 1 has varying counts, so it must differ from the old row-average form;
  # draw 2 is flat, so it must agree -- the two halves of the same mechanism.
  old <- rowMeans(log1p(1 / mu))
  got <- MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot)
  expect_false(isTRUE(all.equal(got[1], old[1])))
  expect_equal(got[2], old[2])

  shape <- c(5, 2)
  expect_equal(
    MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot, extra = 1 / shape),
    log1p(1 / rowMeans(mu) + 1 / shape))

  w <- c(3, 1, 1)
  expect_equal(
    MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot, w = w),
    log1p(1 / vapply(seq_len(2), function(d) sum(w * mu[d, ]) / sum(w), numeric(1))))
})

test_that("Poisson VPC uses the global mean count, leaving the null model intact", {
  skip_on_cran()
  set.seed(3)
  n <- 1200
  d <- data.frame(gender = sample(c("M", "F"), n, replace = TRUE),
                  race   = sample(c("A", "B", "C"), n, replace = TRUE),
                  stringsAsFactors = FALSE)
  sk <- interaction(d$gender, d$race, drop = TRUE)
  u  <- stats::rnorm(nlevels(sk), sd = 0.7)[sk]
  d$x <- stats::rnorm(n)
  # A strong fixed effect spreads the fitted counts, which is the only condition
  # under which the two forms differ at all.
  d$y <- stats::rpois(n, exp(0.7 + u + 1.2 * d$x))
  d$stratum <- make_strata(d, vars = c("gender", "race"))$data$stratum

  adj <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, family = "poisson")))
  mu  <- MAIHDA:::maihda_count_marginal_mu_lme4(adj$model)
  rv  <- MAIHDA:::maihda_residual_variance_lme4(adj$model)

  expect_equal(rv, log1p(1 / mean(mu)))
  expect_false(isTRUE(all.equal(rv, mean(log1p(1 / mu)), tolerance = 1e-6)))
  # The counts really do vary here, so the contrast is not vacuous.
  expect_gt(stats::sd(mu), 0)
  expect_gt(mean(log1p(1 / mu)), rv)

  # The null model's counts are constant, so its residual variance -- and hence
  # the headline VPC -- is byte-identical to the pre-fix value and to Nakagawa's
  # exp(beta0 + sigma^2 / 2) plug-in.
  null <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "poisson")))
  mu0 <- MAIHDA:::maihda_count_marginal_mu_lme4(null$model)
  rv0 <- MAIHDA:::maihda_residual_variance_lme4(null$model)
  expect_length(unique(round(mu0, 12)), 1L)
  expect_equal(rv0, mean(log1p(1 / mu0)))
  b0 <- unname(lme4::fixef(null$model)[1])
  s2 <- as.numeric(lme4::VarCorr(null$model)$stratum[1, 1])
  expect_equal(rv0, log1p(1 / exp(b0 + s2 / 2)), tolerance = 1e-10)
})

test_that("an offset makes the null model's counts vary, and it follows the fix", {
  skip_on_cran()
  set.seed(4)
  n <- 800
  d <- data.frame(gender = sample(c("M", "F"), n, replace = TRUE),
                  race   = sample(c("A", "B", "C"), n, replace = TRUE),
                  stringsAsFactors = FALSE)
  sk <- interaction(d$gender, d$race, drop = TRUE)
  u  <- stats::rnorm(nlevels(sk), sd = 0.6)[sk]
  d$logE <- log(stats::runif(n, 0.5, 4))
  d$y <- stats::rpois(n, exp(0.8 + u + d$logE))
  d$stratum <- make_strata(d, vars = c("gender", "race"))$data$stratum

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ offset(logE) + (1 | stratum), data = d, family = "poisson")))
  mu <- MAIHDA:::maihda_count_marginal_mu_lme4(m$model)
  rv <- MAIHDA:::maihda_residual_variance_lme4(m$model)

  # Even with no covariates the offset spreads lambda_i, so this null-shaped model
  # is one of the cases the fix moves.
  expect_gt(stats::sd(mu), 0)
  expect_equal(rv, log1p(1 / mean(mu)))
  expect_false(isTRUE(all.equal(rv, mean(log1p(1 / mu)), tolerance = 1e-6)))
})

test_that("the negative-binomial level-1 variance follows the same reduction", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(9)
  n <- 900
  d <- data.frame(gender = sample(c("M", "F"), n, replace = TRUE),
                  race   = sample(c("A", "B", "C"), n, replace = TRUE),
                  stringsAsFactors = FALSE)
  sk <- interaction(d$gender, d$race, drop = TRUE)
  u  <- stats::rnorm(nlevels(sk), sd = 0.6)[sk]
  d$x <- stats::rnorm(n)
  d$y <- MASS::rnegbin(n, mu = exp(0.6 + u + 0.9 * d$x), theta = 2)
  d$stratum <- make_strata(d, vars = c("gender", "race"))$data$stratum

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, family = "negbinomial")))
  mu    <- MAIHDA:::maihda_count_marginal_mu_lme4(m$model)
  theta <- MAIHDA:::maihda_negbin_theta_lme4(m$model)
  rv    <- MAIHDA:::maihda_residual_variance_lme4(m$model)

  expect_equal(rv, log1p(1 / mean(mu) + 1 / theta))
  expect_false(isTRUE(all.equal(rv, mean(log1p(1 / mu + 1 / theta)),
                                tolerance = 1e-6)))
})

test_that("the longitudinal count VPC(t) grid reduces on the counts at each time", {
  skip_on_cran()
  set.seed(31)
  nid <- 160
  d <- data.frame(id   = rep(seq_len(nid), each = 5),
                  wave = rep(0:4, nid),
                  g    = rep(sample(c("M", "F"), nid, replace = TRUE), each = 5),
                  r    = rep(sample(c("A", "B"), nid, replace = TRUE), each = 5),
                  stringsAsFactors = FALSE)
  sk <- interaction(d$g, d$r, drop = TRUE)
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rpois(nrow(d),
                      exp(0.4 + stats::rnorm(nlevels(sk), sd = 0.5)[sk] +
                            0.18 * d$wave + 1.1 * d$x))
  d$stratum <- make_strata(d, vars = c("g", "r"))$data$stratum

  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, family = "poisson",
               id = "id", time = "wave")))
  lng    <- fit$longitudinal_info
  center <- lng$time_center
  grid_c <- MAIHDA:::maihda_longitudinal_time_grid(d$wave) - center
  Ss <- MAIHDA:::maihda_re_block_lme4(fit$model, "stratum", lng$time_term,
                                      lng$time_degree)
  Si <- MAIHDA:::maihda_re_block_lme4(fit$model, lng$id, lng$time_term,
                                      lng$time_degree)
  vt <- MAIHDA:::maihda_var_at_time(Ss, grid_c) + MAIHDA:::maihda_var_at_time(Si, grid_c)

  got <- MAIHDA:::maihda_longitudinal_resid_grid_lme4(
    fit$model, lng$time_term, grid_c, vt,
    orig_time = lng$time, center = center)

  frame <- MAIHDA:::maihda_model_frame(fit$model)
  w     <- MAIHDA:::maihda_fit_prior_weights(fit$model)
  mu_at <- function(j) {
    nd  <- MAIHDA:::maihda_longitudinal_set_time(frame, lng$time_term, grid_c[j],
                                                 orig_time = lng$time,
                                                 center = center)
    eta <- MAIHDA:::maihda_lme4_fixed_link(fit$model, nd, offset = NULL)
    pmax(exp(eta + vt[j] / 2), .Machine$double.eps)
  }
  expected <- vapply(seq_along(grid_c),
                     function(j) log1p(1 / MAIHDA:::maihda_weighted_obs_mean(mu_at(j), w)),
                     numeric(1))
  old <- vapply(seq_along(grid_c),
                function(j) MAIHDA:::maihda_weighted_obs_mean(log1p(1 / mu_at(j)), w),
                numeric(1))

  expect_equal(got, expected)
  # The covariate spreads the counts within every wave, so each grid time moves.
  expect_true(all(abs(got - old) > 1e-6))
  expect_true(all(old > got))
})

# ---- Finding 2: aggregation is a structural property, not a weight value -----

test_that("integer precision weights are not read as binomial trial counts", {
  skip_on_cran()
  set.seed(11)
  n <- 600
  d <- data.frame(gender = sample(c("M", "F"), n, replace = TRUE),
                  race   = sample(c("A", "B", "C"), n, replace = TRUE),
                  stringsAsFactors = FALSE)
  sk <- interaction(d$gender, d$race, drop = TRUE)
  u  <- stats::rnorm(nlevels(sk), sd = 1.1)[sk]
  d$x <- stats::rnorm(n)
  d$y <- stats::rbinom(n, 1, stats::plogis(-0.4 + u + 0.8 * d$x))
  d$stratum <- make_strata(d, vars = c("gender", "race"))$data$stratum

  auc_for <- function(wvec) {
    dd <- d; dd$w <- wvec
    m <- suppressWarnings(suppressMessages(
      fit_maihda(y ~ (1 | stratum), data = dd, family = "binomial", weights = w)))
    da <- maihda_discriminatory_accuracy(m)
    prob <- predict_maihda(m, type = "individual", scale = "response")
    # Observed case/control counts, not weight mass.
    expect_equal(da$n_case, sum(d$y == 1))
    expect_equal(da$n_control, sum(d$y == 0))
    # The ordinary observation-level concordance, flagged as ignoring the weights.
    expect_equal(da$auc, maihda_auc(prob, d$y))
    expect_true(isTRUE(da$precision_weights_ignored))
    expect_false(isTRUE(da$weighted))
    expect_null(MAIHDA:::maihda_da_aggregated_counts(m))
    da$auc
  }

  set.seed(11)
  w_int <- sample(1:5, n, replace = TRUE)
  a_int <- auc_for(w_int)
  # The discontinuity the finding turned on: an arbitrarily small perturbation of
  # the weights used to switch the estimand. It no longer does.
  a_eps <- auc_for(w_int + 1e-6)
  expect_equal(a_int, a_eps, tolerance = 1e-5)
  # A uniform integer weight used to double both reported totals.
  auc_for(rep(2, n))
  auc_for(w_int + 0.5)
})

test_that("genuine aggregated binomial fits are still detected structurally", {
  skip_on_cran()
  set.seed(7)
  d <- data.frame(gender = rep(c("M", "F"), each = 3),
                  race   = rep(c("A", "B", "C"), 2),
                  stringsAsFactors = FALSE)
  d <- d[rep(seq_len(6), each = 8), ]
  d$trials <- sample(5:40, nrow(d), replace = TRUE)
  sk <- interaction(d$gender, d$race, drop = TRUE)
  u  <- stats::rnorm(nlevels(sk), sd = 0.9)[sk]
  d$succ <- stats::rbinom(nrow(d), d$trials, stats::plogis(-0.3 + u))
  d$fail <- d$trials - d$succ
  d$prop <- d$succ / d$trials
  d$stratum <- make_strata(d, vars = c("gender", "race"))$data$stratum
  # The fallback's discriminator: a genuine proportion response.
  expect_true(any(d$prop > 0 & d$prop < 1))

  # (a) cbind(): read exactly off the matrix response.
  m_cb <- suppressWarnings(suppressMessages(
    fit_maihda(cbind(succ, fail) ~ (1 | stratum), data = d, family = "binomial")))
  cb <- MAIHDA:::maihda_da_aggregated_counts(m_cb)
  expect_equal(cb$successes, as.numeric(d$succ))
  expect_equal(cb$trials, as.numeric(d$trials))
  da_cb <- maihda_discriminatory_accuracy(m_cb)
  expect_equal(da_cb$n_case, sum(d$succ))
  expect_equal(da_cb$n_control, sum(d$fail))

  # (b) R's other aggregated idiom -- proportion response, trials as prior weights
  # (?glm). This is what the prior-weights fallback exists for, and it must give
  # the same answer as the cbind() spelling of the same data.
  m_pr <- suppressWarnings(suppressMessages(
    fit_maihda(prop ~ (1 | stratum), data = d, family = "binomial",
               weights = trials)))
  pr <- MAIHDA:::maihda_da_aggregated_counts(m_pr)
  expect_equal(pr$successes, as.numeric(d$succ))
  expect_equal(pr$trials, as.numeric(d$trials))
  da_pr <- maihda_discriminatory_accuracy(m_pr)
  expect_equal(da_pr$n_case, sum(d$succ))
  expect_equal(da_pr$n_control, sum(d$fail))
  expect_equal(da_pr$auc, da_cb$auc)

  # (c) An aggregated fit carrying extra precision weights: the matrix response is
  # checked FIRST, so the trial counts stay exact even though the prior weights
  # are now trials * precision weight.
  d$pw <- rep(c(1, 2), length.out = nrow(d))
  m_cw <- suppressWarnings(suppressMessages(
    fit_maihda(cbind(succ, fail) ~ (1 | stratum), data = d, family = "binomial",
               weights = pw)))
  cw <- MAIHDA:::maihda_da_aggregated_counts(m_cw)
  expect_equal(cw$trials, as.numeric(d$trials))
  expect_equal(cw$successes, as.numeric(d$succ))
})
