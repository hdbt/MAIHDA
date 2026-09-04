# Audit 2026-09-02: a GLMM's fixed-effect p-values were a Wald z, which is badly
# anticonservative with few strata. Measured on a binomial MAIHDA
# (ev ~ d1 + d2 + d3 + (1 | stratum), all dimension betas 0, n = 120 per stratum,
# stratum SD 0.5, 8 strata, nominal 5%): the z rejected at 21.8% over 400
# replicates x 3 dimension coefficients.
#
# Five candidate repairs were measured and REJECTED. The run each number comes
# from is named, because they are NOT all from the same replicates:
#   containment t alone (the Gaussian pass's fix)  10.3%  (500 reps x 3 coefs)
#   bootstrap SD + containment t                   10.0%  (400 reps x 3 coefs)
#   studentized bootstrap-t from the fitted model  13.3%  (400 reps x 3 coefs)
#   percentile bootstrap from the fitted model     18.9%  (400 reps x 3 coefs)
#   REML-style SE rescaling + containment t         3.4%  (500 reps x 3 coefs)
# The reason none works is that the standard error itself is wrong, not just its
# reference distribution: measured Wald SE / true sampling SD = 0.701, and the
# bootstrap SD of the same coefficient = 0.703. glmer is ML-only, so its stratum
# variance carries the downward bias REML removes, and in 10-13% of adjusted fits
# it lands exactly on the boundary -- where SE / truth falls to 0.387 and the z
# rejects at 49%. An UNRESTRICTED bootstrap plugs that same tau_hat in and
# inherits the whole error.
#
# The last candidate is the most instructive failure: 3.4% overall looks like a
# conservative fix, but conditioning on the fit shows 10.5% on the singular tenth
# against 2.6% on the rest -- two errors cancelling, not calibration.
#
# What works is a null-RESTRICTED bootstrap: simulate from the model WITHOUT the
# tested term, whose stratum variance is estimated from more information (mean
# tau_hat 0.341 vs 0.286, singular 3.7% vs 10.0%), refit the full model on each
# draw, and refer the observed Wald statistic to the resulting |t*| distribution.
# summary(df_method = "bootstrap") does that. On the same 8-strata design, scored
# on one coefficient over 400 replicates at B = 99 (MC SE 0.011), it rejected at
# 7.0% where the z rejected at 19.3% -- the z figure differs from the 21.8% above
# because that run scored three coefficients under a different seed.
# The residual error is concentrated in the tenth of fits whose stratum variance
# is singular, where it still rejects at 43%; on a non-singular fit it is 3.1%.
# A singular fit is already reported by print() as unreliable.
#
# Power is not the price: at beta = 0.8 the bootstrap rejects at 0.425 [0.356,
# 0.494] (200 reps) against a ceiling of 0.466 for an EXACTLY 5% test on this
# design (noncentrality 0.8 / 0.427). The z's 0.790 is above that ceiling, which
# is only reachable at its broken level.
#
# The tests below pin the CONSTRUCTION, not the level: a level is a 90-minute
# simulation, not a unit test.

maihda_boot_fe_data <- function(n_per = 40, seed = 3) {
  set.seed(seed)
  grid <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:2))
  d <- grid[rep(seq_len(nrow(grid)), each = n_per), , drop = FALSE]
  d$stratum <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  u <- stats::rnorm(nlevels(d$stratum), 0, 0.5)
  d$y <- u[as.integer(d$stratum)] + stats::rnorm(nrow(d))
  d$ev <- stats::rbinom(nrow(d), 1,
                        stats::plogis(stats::qlogis(0.3) + u[as.integer(d$stratum)]))
  rownames(d) <- NULL
  d
}

test_that("df_method = 'bootstrap' gives a GLMM a null-restricted reference", {
  skip_on_cran()
  d <- maihda_boot_fe_data()
  fit <- fit_maihda(ev ~ d1 + d2 + d3 + (1 | stratum), data = d,
                    family = stats::binomial())
  set.seed(21)
  s <- suppressWarnings(summary(fit, df_method = "bootstrap", n_boot = 25))
  fe <- s$fixed_effects
  tested <- fe$term != "(Intercept)"

  expect_identical(s$df_method, "bootstrap")
  # There is no t here, so no degrees of freedom are invented for the df column.
  expect_true(all(is.na(fe$df)))
  expect_true(all(is.finite(fe$p_value[tested])))
  expect_true(all(is.finite(fe$lower[tested]) & is.finite(fe$upper[tested])))
  # The intercept has no reduced model to simulate from, so it is NA rather than
  # a Wald z that is miscalibrated in exactly this regime.
  expect_true(is.na(fe$p_value[!tested]))
  expect_true(is.na(fe$lower[!tested]) && is.na(fe$upper[!tested]))
  # The estimate and the model's own standard error are untouched: only the
  # reference distribution changes.
  expect_equal(fe$estimate, unname(lme4::fixef(fit$model)))
  expect_equal(fe$se, unname(sqrt(diag(as.matrix(stats::vcov(fit$model))))))
  expect_equal(fe$statistic, fe$estimate / fe$se)
  # A bootstrap p is never exactly zero: it is (1 + count) / (B + 1).
  expect_true(all(fe$p_value[tested] >= 1 / 26))
  expect_identical(names(fe), names(summary(fit)$fixed_effects))
})

test_that("the reference is built from the model WITHOUT the tested term", {
  skip_on_cran()
  # Pinned on a Gaussian fit, whose refits are cheap enough to replay the whole
  # bootstrap by hand. maihda_bootstrap_fixef() does no RNG work before the
  # per-term simulate() calls, and takes the terms in order, so one seed
  # reproduces its draws exactly. This is what separates the shipped method from
  # an unrestricted bootstrap, which was measured not to work.
  d <- maihda_boot_fe_data()
  m <- lme4::lmer(y ~ d1 + d2 + d3 + (1 | stratum), data = d)
  B <- 40L
  set.seed(77)
  fx <- suppressWarnings(maihda_bootstrap_fixef(m, n_boot = B, conf_level = 0.95))

  X <- lme4::getME(m, "X")
  asg <- attr(X, "assign")
  tl <- attr(stats::terms(stats::formula(m, fixed.only = TRUE)), "term.labels")
  est <- unname(lme4::fixef(m))
  se <- unname(sqrt(diag(as.matrix(stats::vcov(m)))))

  set.seed(77)
  for (k in seq_along(tl)) {
    red <- stats::update(m, formula. = stats::update(stats::formula(m),
                                                     paste(". ~ . -", tl[k])))
    sim <- maihda_simulate_lme4(red, nsim = B)
    cols <- which(asg == k)
    t_star <- matrix(NA_real_, B, length(cols))
    for (i in seq_len(B)) {
      bm <- lme4::refit(m, newresp = sim[[i]])
      t_star[i, ] <- lme4::fixef(bm)[cols] /
        sqrt(diag(as.matrix(stats::vcov(bm))))[cols]
    }
    for (ci in seq_along(cols)) {
      j <- cols[ci]
      ts <- abs(t_star[, ci])
      expect_equal(fx$p_value[j], min(1, (1 + sum(ts >= abs(est[j] / se[j]))) / (B + 1)))
      r <- ceiling(0.05 * (B + 1)) - 1
      crit <- sort(ts, decreasing = TRUE)[r]
      expect_equal(fx$lower[j], est[j] - crit * se[j])
      expect_equal(fx$upper[j], est[j] + crit * se[j])
    }
  }
  # The interval is a bootstrapped critical value times the model's own standard
  # error, so it is symmetric about the estimate and agrees with the p-value:
  # zero falls outside exactly when the statistic beats the critical value.
  tested <- which(asg > 0L)
  expect_equal(fx$upper[tested] - fx$estimate[tested],
               fx$estimate[tested] - fx$lower[tested])
  expect_equal((fx$lower[tested] > 0 | fx$upper[tested] < 0),
               fx$p_value[tested] < 0.05)
})

test_that("the reduced fit keeps the weights and the offset", {
  skip_on_cran()
  # The null distribution is simulated FROM the reduced fit, so a reduced fit
  # that silently dropped the weights or the offset would quietly change the
  # reference. Checked with the data removed, which forces the model-frame
  # fallback -- the path where lme4's non-standard evaluation of `weights` bites.
  set.seed(4)
  g <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:2))
  u <- stats::rnorm(8, 0, 0.5)
  wt <- g[rep(seq_len(8), each = 30), ]
  wt$stratum <- interaction(wt$d1, wt$d2, wt$d3, drop = TRUE)
  wt$w <- rep(c(1, 4, 9), length.out = nrow(wt))
  wt$y <- u[as.integer(wt$stratum)] + stats::rnorm(nrow(wt), 0, 1 / sqrt(wt$w))
  rownames(wt) <- NULL
  m <- lme4::lmer(y ~ d1 + d2 + d3 + (1 | stratum), wt, weights = w)

  red_form <- stats::update(stats::formula(m), ". ~ . - d1")
  red <- maihda_refit_reduced(m, red_form)
  expect_s4_class(red, "merMod")
  expect_equal(unname(stats::weights(red, type = "prior")),
               unname(stats::weights(m, type = "prior")))

  rm(wt)                                  # update() can no longer be used at all
  red2 <- maihda_refit_reduced(m, red_form)
  expect_s4_class(red2, "merMod")
  expect_equal(unname(stats::weights(red2, type = "prior")),
               unname(stats::weights(m, type = "prior")))
  expect_equal(lme4::fixef(red2), lme4::fixef(red), tolerance = 1e-6)

  # The offset half of this test's name: a Poisson rate model, in all three
  # spellings, must reach the same reduced fit with its exposure intact.
  cnt <- g[rep(seq_len(8), each = 8), ]
  cnt$stratum <- interaction(cnt$d1, cnt$d2, cnt$d3, drop = TRUE)
  cnt$expo <- 40
  cnt$logexp <- log(cnt$expo)
  cnt$y <- stats::rpois(nrow(cnt),
                        cnt$expo * exp(-2 + u[as.integer(cnt$stratum)]))
  rownames(cnt) <- NULL
  mo <- lme4::glmer(y ~ d1 + d2 + d3 + (1 | stratum) + offset(logexp),
                    cnt, stats::poisson())
  ro <- maihda_refit_reduced(mo, stats::update(stats::formula(mo), ". ~ . - d1"))
  expect_s4_class(ro, "merMod")
  expect_equal(unname(lme4::getME(ro, "offset")), unname(lme4::getME(mo, "offset")))

  ma <- lme4::glmer(y ~ d1 + d2 + d3 + (1 | stratum), cnt, stats::poisson(),
                    offset = logexp)
  ra <- maihda_refit_reduced(ma, stats::update(stats::formula(ma), ". ~ . - d1"))
  expect_s4_class(ra, "merMod")
  expect_equal(unname(lme4::getME(ra, "offset")), unname(lme4::getME(ma, "offset")))
  expect_equal(unname(lme4::fixef(ra)), unname(lme4::fixef(ro)), tolerance = 1e-6)
})

test_that("the reduced fit covers the full model's rows even if the data changed", {
  skip_on_cran()
  # stats::update() re-evaluates the ORIGINAL call, so it refits on whatever the
  # data argument now names. Left unchecked, a data object that has changed since
  # the fit yields a reduced model for a DIFFERENT dataset -- whose simulated
  # responses are then the wrong length for refit(), or worse, the right length
  # and the wrong rows. The same mismatch arises when the dropped term carried
  # missing values. Both must fall through to the stored model frame.
  d <- maihda_boot_fe_data()
  dd <- d
  m <- lme4::lmer(y ~ d1 + d2 + d3 + (1 | stratum), data = dd)
  n_full <- stats::nobs(m)
  dd <- dd[seq_len(200), ]                      # same name, fewer rows

  red <- maihda_refit_reduced(m, stats::update(stats::formula(m), ". ~ . - d1"))
  expect_s4_class(red, "merMod")
  expect_equal(stats::nobs(red), n_full)
  expect_false("d1" %in% names(lme4::fixef(red)))
  # And the whole bootstrap still runs rather than erroring or silently using
  # the truncated data.
  set.seed(4)
  fx <- suppressWarnings(maihda_bootstrap_fixef(m, n_boot = 25, conf_level = 0.95))
  expect_true(all(is.finite(fx$p_value[fx$term != "(Intercept)"])))
})

test_that("a formula that transforms a column is not refitted from the frame", {
  skip_on_cran()
  # The model frame stores EVALUATED columns under their deparsed names, so
  # refitting log(x) against it would take the log a second time. The guard must
  # refuse rather than silently double-transform.
  d <- maihda_boot_fe_data()
  d$x <- exp(stats::rnorm(nrow(d)))
  m <- lme4::lmer(y ~ d1 + log(x) + (1 | stratum), data = d)
  fr <- m@frame
  expect_true("log(x)" %in% names(fr))
  expect_false(maihda_frame_refit_safe(y ~ log(x) + (1 | stratum), fr))
  expect_true(maihda_frame_refit_safe(y ~ d1 + (1 | stratum), fr))
})

test_that("summary() on an analysis refuses arguments it cannot apply", {
  skip_on_cran()
  # maihda() computes its summaries once. Silently dropping df_method would hand
  # back an uncorrected Wald reference while looking as though the requested one
  # had been applied -- the exact failure mode this pass exists to remove.
  d <- maihda_boot_fe_data()
  d$stratum <- NULL
  a <- suppressMessages(suppressWarnings(maihda(y ~ (1 | d1:d2:d3), data = d)))

  expect_error(summary(a, df_method = "bootstrap"), "cannot apply")
  expect_error(summary(a, df_method = "bootstrap"), "df_method")
  expect_error(summary(a, bootstrap = TRUE), "cannot apply")
  # Ordinary use is untouched.
  expect_s3_class(summary(a), "maihda_summary")
  expect_s3_class(summary(a, which = "adjusted"), "maihda_summary")
})

test_that("the fixed-effect bootstrap survives every supported response shape", {
  skip_on_cran()
  # An aggregated binomial reaches refit() as a two-column response, and a
  # weighted Gaussian goes down maihda_simulate_lme4()'s separate weights branch.
  # Neither is exercised by the plain Bernoulli tests above.
  set.seed(4)
  g <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:2))
  agg <- g[rep(seq_len(8), each = 5), ]
  agg$stratum <- interaction(agg$d1, agg$d2, agg$d3, drop = TRUE)
  u <- stats::rnorm(8, 0, 0.5)
  agg$trials <- 40
  agg$succ <- stats::rbinom(nrow(agg), agg$trials,
                            stats::plogis(stats::qlogis(0.3) + u[as.integer(agg$stratum)]))
  agg$fail <- agg$trials - agg$succ
  agg$prop <- agg$succ / agg$trials
  rownames(agg) <- NULL

  fits <- list(
    cbind_binom = lme4::glmer(cbind(succ, fail) ~ d1 + d2 + d3 + (1 | stratum),
                              agg, stats::binomial()),
    prop_weights = lme4::glmer(prop ~ d1 + d2 + d3 + (1 | stratum), agg,
                               stats::binomial(), weights = trials)
  )
  out <- lapply(fits, function(m) {
    set.seed(8)
    suppressWarnings(maihda_bootstrap_fixef(m, n_boot = 25, conf_level = 0.95))
  })
  for (r in out) {
    k <- r$term != "(Intercept)"
    expect_true(all(is.finite(r$p_value[k])))
    expect_true(all(r$lower[k] < r$estimate[k] & r$estimate[k] < r$upper[k]))
  }
  # The two spellings of the same aggregated binomial must give the same answer.
  expect_equal(out$cbind_binom$p_value, out$prop_weights$p_value)
  expect_equal(out$cbind_binom$lower, out$prop_weights$lower)
})

test_that("a fixed-effect bootstrap validates n_boot without bootstrap = TRUE", {
  d <- maihda_boot_fe_data()
  fit <- fit_maihda(y ~ d1 + d2 + d3 + (1 | stratum), data = d)
  # bootstrap = FALSE, so before this pass n_boot was never looked at on this path.
  expect_error(summary(fit, df_method = "bootstrap", n_boot = 3), "n_boot")
  expect_error(summary(fit, df_method = "bootstrap", conf_level = 1.4), "conf_level")
  expect_warning(summary(fit, df_method = "bootstrap", n_boot = 20), "n_boot = 20")
})

test_that("df_method = 'bootstrap' is refused where it cannot be built", {
  skip_on_cran()
  d <- maihda_boot_fe_data()
  # An intercept-only fixed part has no term to drop, so there is no null model
  # to simulate from -- say so rather than returning a table of NAs.
  m <- lme4::lmer(y ~ (1 | stratum), data = d)
  expect_error(maihda_bootstrap_fixef(m, n_boot = 25, conf_level = 0.95),
               "no fixed-effect term")

  skip_if_not_installed("WeMix")
  d$sw <- 1
  fit <- suppressWarnings(fit_maihda(y ~ d1 + d2 + d3 + (1 | stratum), data = d,
                                     engine = "wemix", sampling_weights = "sw"))
  # Refused before n_boot is validated, so the message is about the engine.
  expect_error(summary(fit, df_method = "bootstrap", n_boot = 20), "lme4 engine")
})

test_that("the fixed-effect bootstrap leaves the fitted model untouched", {
  d <- maihda_boot_fe_data()
  m <- lme4::lmer(y ~ d1 + d2 + d3 + (1 | stratum), data = d)
  before <- list(fixef = lme4::fixef(m), theta = lme4::getME(m, "theta"),
                 sigma = stats::sigma(m), vcov = as.matrix(stats::vcov(m)),
                 logLik = as.numeric(stats::logLik(m)))
  set.seed(9)
  invisible(suppressWarnings(maihda_bootstrap_fixef(m, n_boot = 25, conf_level = 0.95)))

  expect_equal(lme4::fixef(m), before$fixef, tolerance = 0)
  expect_equal(lme4::getME(m, "theta"), before$theta, tolerance = 0)
  expect_equal(stats::sigma(m), before$sigma, tolerance = 0)
  expect_equal(as.matrix(stats::vcov(m)), before$vcov, tolerance = 0)
  expect_equal(as.numeric(stats::logLik(m)), before$logLik, tolerance = 0)
})

test_that("the defaults are unchanged: Gaussian keeps the t, a GLMM keeps the z", {
  d <- maihda_boot_fe_data()
  gauss <- summary(fit_maihda(y ~ d1 + d2 + d3 + (1 | stratum), data = d))
  expect_identical(gauss$df_method, "between-within")
  expect_equal(gauss$fixed_effects$df, rep(4, 4))

  glmm_fit <- fit_maihda(ev ~ d1 + d2 + d3 + (1 | stratum), data = d,
                         family = stats::binomial())
  glmm <- summary(glmm_fit)
  expect_identical(glmm$df_method, "normal")
  expect_true(all(is.na(glmm$fixed_effects$df)))
  # Still lme4's own reported p-value when nothing is asked for.
  expect_equal(glmm$fixed_effects$p_value,
               unname(stats::coef(summary(glmm_fit$model))[, "Pr(>|z|)"]),
               tolerance = 1e-8)
  # And the containment rule still declines to answer for a glmerMod: a t cannot
  # repair a standard error, which is what is wrong there.
  expect_null(maihda_containment_df(glmm_fit$model))
})
