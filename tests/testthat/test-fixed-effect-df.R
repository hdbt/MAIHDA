# Audit 2026-09-01: Gaussian fixed-effect p-values were synthesized from a Wald z
# with no denominator degrees of freedom, which is anticonservative for terms
# constant within a stratum -- the intercept and an adjusted model's dimension
# main effects, whose effective sample size is the number of strata, not n.
# Reproduced at 400 replicates: a nominal 5% z test on those terms rejected at
# 7-11% with 12 strata and 15-16% with 8. summary() now refers a Gaussian lme4
# fit's Wald statistic to a t on containment (between-within) degrees of freedom,
# which needs nothing beyond lme4. The three tests below pin how closely that
# tracks Kenward-Roger across the balance regimes: exactly when balanced, to ~1%
# under mild imbalance, and one design-level value inside KR's per-term span when
# severely unbalanced.

maihda_df_test_data <- function(n_per = 60, seed = 1) {
  set.seed(seed)
  grid <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:3))
  d <- grid[rep(seq_len(nrow(grid)), each = n_per), , drop = FALSE]
  d$stratum <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nlevels(d$stratum), 0, 0.4)[as.integer(d$stratum)] +
    0.3 * d$x + stats::rnorm(nrow(d))
  rownames(d) <- NULL
  d
}

test_that("a Gaussian adjusted fit tests dimension main effects against the strata", {
  d <- maihda_df_test_data()
  fit <- fit_maihda(y ~ d1 + d2 + d3 + x + (1 | stratum), data = d)
  fe <- summary(fit)$fixed_effects

  expect_true("df" %in% names(fe))
  # 12 strata minus the 5 columns constant within a stratum (intercept + 1 + 1 + 2).
  strat_terms <- c("(Intercept)", "d12", "d22", "d32", "d33")
  expect_equal(fe$df[match(strat_terms, fe$term)], rep(7, length(strat_terms)))
  # A covariate varying within stratum keeps essentially all of n.
  expect_equal(fe$df[fe$term == "x"], nrow(d) - 12 - 1)
})

test_that("the t reference widens intervals and p-values relative to the z", {
  d <- maihda_df_test_data()
  fit <- fit_maihda(y ~ d1 + d2 + d3 + x + (1 | stratum), data = d)
  t_fe <- summary(fit)$fixed_effects
  z_fe <- summary(fit, df_method = "normal")$fixed_effects

  expect_equal(t_fe$estimate, z_fe$estimate)
  expect_equal(t_fe$se, z_fe$se)
  expect_equal(t_fe$statistic, z_fe$statistic)
  # Never anticonservative relative to the old behaviour: every p-value is at
  # least as large and every interval at least as wide.
  expect_true(all(t_fe$p_value >= z_fe$p_value - 1e-12))
  expect_true(all((t_fe$upper - t_fe$lower) >= (z_fe$upper - z_fe$lower) - 1e-12))
  # And materially so on a stratum-level term.
  i <- which(t_fe$term == "d33")
  expect_gt((t_fe$upper[i] - t_fe$lower[i]) / (z_fe$upper[i] - z_fe$lower[i]), 1.15)

  expect_true(all(is.na(z_fe$df)))
  expect_identical(summary(fit, df_method = "normal")$df_method, "normal")
  expect_identical(summary(fit)$df_method, "between-within")
})

# pbkrtest sits in Suggests, not Imports, so the package takes on no hard
# dependency for it. It has to be declared: "checking for unstated dependencies
# in 'tests'" walks tests/ recursively, tests/testthat/ included, and an
# undeclared '::' here is a WARNING that fails R-CMD-check under
# error-on: warning. These cross-checks are the only external validation that
# the containment rule is statistically right rather than merely
# self-consistent, so they run wherever pbkrtest is installed and skip
# elsewhere.
maihda_kr_df <- function(m) {
  est <- lme4::fixef(m)
  adj <- pbkrtest::vcovAdj(m)
  vapply(seq_along(est), function(j) {
    L <- rep(0, length(est)); L[j] <- 1
    pbkrtest::Lb_ddf(L, stats::vcov(m), adj)
  }, numeric(1))
}

test_that("containment df match Kenward-Roger for a BALANCED adjusted MAIHDA", {
  skip_if_not_installed("pbkrtest")
  d <- maihda_df_test_data()                       # equal n per stratum
  m <- fit_maihda(y ~ d1 + d2 + d3 + x + (1 | stratum), data = d)$model

  ours <- maihda_containment_df(m)
  # Kenward-Roger differentiates numerically, so it lands on 7.00 rather than a
  # bit-exact 7; the agreement is to its own tolerance, not to machine epsilon.
  expect_equal(unname(ours), maihda_kr_df(m), tolerance = 5e-3)
})

test_that("containment df still track Kenward-Roger under MILD imbalance", {
  skip_if_not_installed("pbkrtest")
  set.seed(5)
  grid <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:3))
  n_per <- sample(30:300, 12, replace = TRUE)      # ~6x largest/smallest
  d <- grid[rep(seq_len(12), times = n_per), , drop = FALSE]
  d$stratum <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(12, 0, 0.4)[as.integer(d$stratum)] + 0.3 * d$x +
    stats::rnorm(nrow(d))
  rownames(d) <- NULL

  m <- fit_maihda(y ~ d1 + d2 + d3 + x + (1 | stratum), data = d)$model
  ours <- unname(maihda_containment_df(m))
  kr <- maihda_kr_df(m)
  expect_lt(max(abs(ours - kr) / kr), 0.05)        # observed: 0.009
})

test_that("under SEVERE imbalance the single containment value sits inside the KR span", {
  skip_if_not_installed("pbkrtest")
  # maihda_health_data's strata run from 1 to 349 observations (349x). Containment
  # is a design-level rule, so it reports ONE value for every stratum-level term
  # where Kenward-Roger varies them term by term with the imbalance. It stays
  # inside their span and far below n -- it does not degenerate back towards the
  # z this fix replaced -- but it is not KR, and the docs do not claim it is.
  data("maihda_health_data")
  m <- lme4::lmer(BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
                  data = maihda_health_data)
  ours <- maihda_containment_df(m)
  kr <- maihda_kr_df(m)
  strat <- names(ours) != "Age"

  expect_length(unique(ours[strat]), 1L)
  one <- unique(ours[strat])
  expect_gt(one, min(kr[strat]))
  expect_lt(one, max(kr[strat]))
  expect_lt(one, nrow(maihda_health_data) / 50)    # nowhere near the z's effective n
})

test_that("a NEARLY stratum-constant covariate is a known gap, not an oversight", {
  # Containment changes the reference distribution, never the standard error, so
  # a covariate that is nearly (not exactly) constant within a stratum falls in a
  # gap: it varies within, so nothing absorbs it and it gets residual df, yet its
  # information really is between strata. Measured at 12 strata over 1500-2000
  # replicates, true effect zero: rejection 7-8% at a nominal 5% around ICC 0.9,
  # where Kenward-Roger gives 4.7-5.0%. At ICC 0.5 or below there is no gap
  # (4.9-5.6%), and it closes as strata grow: 7.5% at 8 strata, 5.7% at 30, 5.0% at 50.
  # This test pins the STRUCTURE that produces it, so the boundary stays visible.
  set.seed(12)
  grid <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:3))
  d <- grid[rep(seq_len(12), each = 100), , drop = FALSE]
  d$stratum <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  # ICC ~0.9: almost all of m's variance is between strata, but not all.
  d$m <- stats::rnorm(12, 0, 1)[as.integer(d$stratum)] +
    stats::rnorm(nrow(d), 0, 0.3)
  d$y <- stats::rnorm(12, 0, 0.5)[as.integer(d$stratum)] + stats::rnorm(nrow(d))
  rownames(d) <- NULL

  m <- lme4::lmer(y ~ d1 + d2 + d3 + m + (1 | stratum), data = d)
  df <- maihda_containment_df(m)

  # It is NOT absorbed -- it gets the residual df, the same reference the z used.
  expect_gt(df[["m"]], 1000)
  expect_equal(unname(df[["m"]]), nrow(d) - 12 - 1)
  # while the genuinely stratum-constant terms are still caught.
  expect_equal(unname(df[["d12"]]), 7)

  # Make an exactly-constant version and confirm containment DOES catch that,
  # which is what makes this a discontinuity rather than a blanket failure.
  d$m_exact <- stats::rnorm(12, 0, 1)[as.integer(d$stratum)]
  m2 <- lme4::lmer(y ~ d1 + d2 + d3 + m_exact + (1 | stratum), data = d)
  expect_equal(unname(maihda_containment_df(m2)[["m_exact"]]), 6)
})

test_that("a term with a random slope is tested against its grouping, not n", {
  set.seed(4)
  n_id <- 200; n_t <- 5
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = n_t)),
    time = rep(seq_len(n_t) - 1, n_id)
  )
  d$g1 <- factor(rep(sample(1:2, n_id, TRUE), each = n_t))
  d$g2 <- factor(rep(sample(1:3, n_id, TRUE), each = n_t))
  d$stratum <- interaction(d$g1, d$g2, drop = TRUE)
  d$y <- stats::rnorm(nlevels(d$stratum), 0, 0.3)[as.integer(d$stratum)] +
    stats::rnorm(n_id, 0, 0.5)[as.integer(d$id)] + 0.2 * d$time +
    stats::rnorm(nrow(d))

  m <- suppressWarnings(lme4::lmer(y ~ time + g1 + (time | stratum) + (time | id),
                                   data = d))
  df <- maihda_containment_df(m)
  # `time` varies within every grouping, so a plain between-within split would
  # hand it n - levels df; it carries a random slope at the stratum level, so it
  # is absorbed there and belongs to the 6 strata.
  expect_lt(df[["time"]], 10)
  expect_gte(df[["time"]], 1)
  expect_lt(df[["g12"]], 10)
})

test_that("a level contributing no variance does not cap the degrees of freedom", {
  set.seed(7)
  d <- maihda_df_test_data(n_per = 80, seed = 7)
  # A context grouping with 4 units and no signal: its variance component is
  # driven to (near) zero, so it must not cap the intercept at 3 df.
  d$region <- factor(rep_len(paste0("R", 1:4), nrow(d)))
  m <- suppressWarnings(
    lme4::lmer(y ~ d1 + x + (1 | stratum) + (1 | region), data = d)
  )
  skip_if(lme4::VarCorr(m)$region[1] > 1e-6, "region variance not at the boundary")

  df <- maihda_containment_df(m)
  expect_gt(df[["(Intercept)"]], 3)
})

test_that("a singular stratum variance still tests against the strata, not n", {
  # The materiality filter drops a grouping that contributes nothing to a
  # coefficient's variance. It must not do so when that grouping is the ONLY one
  # absorbing the term: a boundary sigma^2_stratum = 0 would otherwise hand a
  # stratum-level term all of n (df 4 -> 631 between two similar datasets).
  set.seed(2)
  grid <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:2))
  d <- grid[rep(seq_len(8), each = 40), ]
  d$stratum <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  d$y <- stats::rnorm(nrow(d))          # no stratum variance at all
  rownames(d) <- NULL

  m <- suppressWarnings(lme4::lmer(y ~ d1 + d2 + d3 + (1 | stratum), data = d))
  skip_if_not(lme4::isSingular(m), "stratum variance did not hit the boundary")

  df <- maihda_containment_df(m)
  expect_equal(unname(df[c("(Intercept)", "d12", "d22", "d32")]), rep(4, 4))
})

test_that("engines without a finite-sample t keep the Wald z", {
  d <- maihda_df_test_data()
  d$event <- stats::rbinom(nrow(d), 1, stats::plogis(-0.5 + 0.2 * d$x))
  fit <- fit_maihda(event ~ x + (1 | stratum), data = d, family = binomial())
  s <- summary(fit)

  expect_true(all(is.na(s$fixed_effects$df)))
  expect_identical(s$df_method, "normal")
  # The z p-value is the one lme4's own summary reports.
  expect_equal(s$fixed_effects$p_value,
               unname(stats::coef(summary(fit$model))[, "Pr(>|z|)"]),
               tolerance = 1e-8)
})

test_that("tidy() carries the degrees of freedom into the broom shape", {
  d <- maihda_df_test_data()
  fit <- fit_maihda(y ~ d1 + d2 + d3 + x + (1 | stratum), data = d)
  td <- tidy(fit, component = "fixed")

  expect_true(all(c("term", "estimate", "std.error", "statistic", "df",
                    "p.value", "conf.low", "conf.high") %in% names(td)))
  expect_equal(td$df, summary(fit)$fixed_effects$df)
})

test_that("computing the degrees of freedom leaves the fitted model untouched", {
  # The variance shares replay lme4's predictor object at perturbed theta. That
  # object has reference semantics, so a shallow copy would silently corrupt the
  # user's fit -- every quantity below is compared at zero tolerance.
  d <- maihda_df_test_data()
  m <- lme4::lmer(y ~ d1 + d2 + d3 + x + (1 | stratum), data = d)
  before <- list(vcov = as.matrix(stats::vcov(m)), fixef = lme4::fixef(m),
                 theta = lme4::getME(m, "theta"), sigma = stats::sigma(m),
                 ranef = as.matrix(lme4::ranef(m)[[1]]),
                 logLik = as.numeric(stats::logLik(m)),
                 pp_theta = m@pp$theta)

  invisible(maihda_containment_df(m))
  invisible(maihda_variance_shares(m))

  expect_equal(as.matrix(stats::vcov(m)), before$vcov, tolerance = 0)
  expect_equal(lme4::fixef(m), before$fixef, tolerance = 0)
  expect_equal(lme4::getME(m, "theta"), before$theta, tolerance = 0)
  expect_equal(stats::sigma(m), before$sigma, tolerance = 0)
  expect_equal(as.matrix(lme4::ranef(m)[[1]]), before$ranef, tolerance = 0)
  expect_equal(as.numeric(stats::logLik(m)), before$logLik, tolerance = 0)
  expect_equal(m@pp$theta, before$pp_theta, tolerance = 0)

  expect_identical(maihda_containment_df(m), maihda_containment_df(m))
})

test_that("summary() rejects an unknown df_method", {
  d <- maihda_df_test_data()
  fit <- fit_maihda(y ~ d1 + x + (1 | stratum), data = d)
  # match.arg's message is translated, so match on the untranslated option list.
  expect_error(summary(fit, df_method = "kenward-roger"), "between-within")
})
