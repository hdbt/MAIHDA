test_that("plot_prediction_deviation_panels works for gaussian", {
  set.seed(123)
  df <- data.frame(
    x = rnorm(100),
    y = rnorm(100)
  )
  m <- lm(y ~ x, data = df)

  p <- plot_prediction_deviation_panels(m, df, type = "gaussian")
  expect_s3_class(p, "patchwork")
})

test_that("plot_prediction_deviation_panels aligns gaussian predictions to supplied data", {
  set.seed(124)
  df <- data.frame(
    x = rnorm(100),
    y = rnorm(100)
  )
  df$x[1:3] <- NA_real_
  m <- lm(y ~ x, data = df)

  p <- plot_prediction_deviation_panels(m, df, type = "gaussian")

  expect_s3_class(p, "patchwork")
})

test_that("plot_prediction_deviation_panels works for binomial", {
  set.seed(123)
  df <- data.frame(
    x = rnorm(100),
    y = sample(0:1, 100, replace = TRUE)
  )
  m <- glm(y ~ x, data = df, family = binomial)

  p <- plot_prediction_deviation_panels(m, df, type = "binomial")
  expect_s3_class(p, "patchwork")
})

test_that("plot_prediction_deviation_panels aligns binomial predictions to supplied data", {
  set.seed(125)
  df <- data.frame(
    x = rnorm(100),
    y = sample(0:1, 100, replace = TRUE)
  )
  df$x[1:3] <- NA_real_
  m <- glm(y ~ x, data = df, family = binomial)

  p <- plot_prediction_deviation_panels(m, df, type = "binomial")

  expect_s3_class(p, "patchwork")
})

test_that("binomial fallback residuals stay on the deviance scale", {
  fitted <- c(0.2, 0.8, NA_real_, 0.5)
  obs <- c(0L, 1L, 1L, 0L)

  resids <- MAIHDA:::maihda_prediction_panel_binomial_residuals(
    structure(list(), class = "no_residual_method"),
    data.frame(id = seq_along(fitted)),
    fitted,
    obs
  )

  expected <- c(
    sqrt(-2 * log1p(-0.2)),
    sqrt(-2 * log(0.8)),
    0,
    sqrt(-2 * log1p(-0.5))
  )
  expect_equal(resids, expected)
})

test_that("binomial residuals align to supplied row order", {
  set.seed(127)
  df <- data.frame(x = rnorm(80))
  df$y <- rbinom(80, 1, plogis(-0.2 + 0.8 * df$x))
  m <- glm(y ~ x, data = df, family = binomial)
  shuffled <- df[sample(seq_len(nrow(df))), , drop = FALSE]

  fitted <- MAIHDA:::maihda_prediction_panel_fitted(m, shuffled, "binomial")$fit
  obs <- MAIHDA:::maihda_binomial_observed_01(shuffled$y, nrow(shuffled))
  resids <- MAIHDA:::maihda_prediction_panel_binomial_residuals(
    m,
    shuffled,
    fitted,
    obs
  )
  expected <- MAIHDA:::maihda_binomial_abs_deviance_residual(obs, fitted)

  expect_equal(resids, expected)
})

test_that("plot_prediction_deviation_panels handles factor binomial outcomes without coercion warnings", {
  set.seed(456)
  df <- data.frame(
    x = rnorm(100),
    y = factor(sample(c("No", "Yes"), 100, replace = TRUE), levels = c("No", "Yes"))
  )
  m <- glm(y ~ x, data = df, family = binomial)

  expect_warning(
    p <- plot_prediction_deviation_panels(m, df, type = "binomial"),
    NA
  )
  expect_s3_class(p, "patchwork")
  expect_equal(
    MAIHDA:::maihda_binomial_observed_01(df$y, nrow(df)),
    as.integer(df$y == "Yes")
  )
})

test_that("plot_prediction_deviation_panels aligns ordinal predictions to supplied data", {
  skip_if_not_installed("MASS")

  set.seed(126)
  df <- data.frame(
    x = rnorm(100),
    y = ordered(
      sample(c("low", "mid", "high"), 100, replace = TRUE),
      levels = c("low", "mid", "high")
    )
  )
  df$x[1:3] <- NA_real_
  m <- MASS::polr(y ~ x, data = df, Hess = TRUE)

  p <- plot_prediction_deviation_panels(m, df, type = "ordinal")
  p_expected <- plot_prediction_deviation_panels(
    m,
    df,
    type = "ordinal",
    ordinal_mode = "expected_score"
  )

  expect_s3_class(p, "patchwork")
  expect_s3_class(p_expected, "patchwork")
})

test_that("prediction deviation auto-detects brmsfit families", {
  fake_bernoulli <- structure(
    list(family = structure(list(family = "bernoulli", link = "logit"), class = "family")),
    class = "brmsfit"
  )
  fake_gaussian <- structure(
    list(family = structure(list(family = "gaussian", link = "identity"), class = "family")),
    class = "brmsfit"
  )
  fake_ordinal <- structure(
    list(family = structure(list(family = "cumulative", link = "logit"), class = "family")),
    class = "brmsfit"
  )
  fake_poisson <- structure(
    list(family = structure(list(family = "poisson", link = "log"), class = "family")),
    class = "brmsfit"
  )

  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_bernoulli), "binomial")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_gaussian), "gaussian")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_ordinal), "ordinal")
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(fake_poisson), "poisson")
})

test_that("ordinal stratum surprise averages per-observation surprise (log loss)", {
  skip_if_not_installed("MASS")
  set.seed(2108)
  n <- 300
  d <- data.frame(stratum = factor(rep(seq_len(6), each = 50)), x = rnorm(n))
  eta <- 0.6 * d$x + rnorm(6, sd = 0.8)[d$stratum]
  d$y <- factor(findInterval(eta + rlogis(n), c(-0.8, 0.8)),
                labels = c("low", "mid", "high"), ordered = TRUE)
  m <- MASS::polr(y ~ x, data = d, Hess = TRUE)

  p <- plot_prediction_deviation_panels(m, d, type = "ordinal",
                                        ordinal_mode = "surprise")
  plotted <- p[[2]]$data            # the surprise panel's data, per stratum

  # Reference: per-observation surprise -log(P(observed)), averaged per stratum.
  probs <- predict(m, type = "probs")
  obs <- as.character(d$y)
  p_obs <- vapply(seq_len(nrow(d)),
                  function(i) probs[i, match(obs[i], colnames(probs))], numeric(1))
  ref <- tapply(-log(p_obs), d$stratum, mean)

  m_plot <- stats::setNames(plotted$surprise, as.character(plotted$stratum))
  common <- intersect(names(m_plot), names(ref))
  p_mean <- tapply(p_obs, d$stratum, mean)
  expect_gt(length(common), 0)
  expect_equal(as.numeric(m_plot[common]), as.numeric(ref[common]), tolerance = 1e-6)
  # The mean-of-logs differs from the (incorrect) log-of-mean it replaced.
  expect_false(isTRUE(all.equal(as.numeric(m_plot[common]),
                                as.numeric(-log(p_mean[common])))))
})

test_that("maihda_stratum_interval_from_epred aggregates draws then takes quantiles", {
  # 4 posterior draws x 4 observations; two strata of two rows each. The correct
  # stratum interval aggregates predictions WITHIN each draw, then takes posterior
  # quantiles -- it is not a function of any row-level SE.
  ep <- rbind(
    c(0.10, 0.20, 0.60, 0.80),
    c(0.15, 0.25, 0.55, 0.85),
    c(0.20, 0.10, 0.65, 0.75),
    c(0.05, 0.30, 0.50, 0.90)
  )
  strat <- c("A", "A", "B", "B")

  out <- MAIHDA:::maihda_stratum_interval_from_epred(ep, strat, level = 0.5)
  expect_equal(out$stratum, c("A", "B"))

  mA <- rowMeans(ep[, 1:2])
  mB <- rowMeans(ep[, 3:4])
  expect_equal(out$fitted, c(mean(mA), mean(mB)))
  qA <- stats::quantile(mA, c(0.25, 0.75), names = FALSE)
  qB <- stats::quantile(mB, c(0.25, 0.75), names = FALSE)
  expect_equal(out$ci_lower, c(qA[1], qB[1]))
  expect_equal(out$ci_upper, c(qA[2], qB[2]))

  # Weights bias the per-draw stratum mean (a weighted mean over the stratum rows).
  w <- c(3, 1, 1, 1)
  outw <- MAIHDA:::maihda_stratum_interval_from_epred(ep, strat, w = w, level = 0.5)
  mAw <- as.vector(ep[, 1:2] %*% (c(3, 1) / 4))
  expect_equal(outw$fitted[1], mean(mAw))

  expect_error(MAIHDA:::maihda_stratum_interval_from_epred(ep, strat[1:3]),
               "one label per column")
})

test_that("frequentist (lme4) stratum panels omit error bars rather than averaging SEs", {
  set.seed(404)
  n <- 900
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y  <- stats::rbinom(n, 1, stats::plogis(stats::rnorm(nlevels(sk), sd = 0.8)[sk]))
  d$yc <- stats::rnorm(n, stats::rnorm(nlevels(sk), sd = 0.8)[sk])
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum

  mb <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "binomial")))
  # suppressWarnings: predict.merMod(se.fit = TRUE) emits an "approximation"
  # notice; the point of the fix is precisely that this row-level SE is NOT
  # aggregated into a (statistically invalid) stratum SE -- the stratum interval
  # is omitted for a frequentist fit, so ci is NA rather than fabricated.
  pb <- suppressWarnings(plot_prediction_deviation_panels(mb, type = "binomial"))
  expect_s3_class(pb, "patchwork")
  expect_true(all(is.na(pb[[2]]$data$ci_lower)))
  expect_true(all(is.na(pb[[2]]$data$ci_upper)))

  mg <- suppressWarnings(suppressMessages(
    fit_maihda(yc ~ (1 | stratum), data = d, family = "gaussian")))
  pg <- suppressWarnings(plot_prediction_deviation_panels(mg, type = "gaussian"))
  expect_true(all(is.na(pg[[2]]$data$ci_lower)))
  expect_true(all(is.na(pg[[2]]$data$ci_upper)))
})

test_that("Poisson prediction panels use response-scale (count) fitted values", {
  set.seed(2005)
  df <- data.frame(x = rnorm(150))
  df$y <- rpois(150, lambda = exp(1 + 0.3 * df$x))
  m <- glm(y ~ x, data = df, family = poisson)

  # Poisson is no longer routed through the Gaussian (link-scale) branch.
  expect_equal(MAIHDA:::maihda_prediction_panel_auto_type(m), "poisson")

  fitted_used <- MAIHDA:::maihda_prediction_panel_fitted(m, df, "poisson")$fit
  expect_equal(fitted_used, unname(predict(m, type = "response")), tolerance = 1e-8)
  # ...and NOT the link (log) scale the old Gaussian routing produced.
  expect_false(isTRUE(all.equal(fitted_used, unname(predict(m, type = "link")))))

  expect_s3_class(plot_prediction_deviation_panels(m, df), "patchwork")
})
