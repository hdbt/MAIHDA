# Regression tests for the 2026-07-23 statistical audit (fortieth pass;
# test-audit-2026-07-23.R belongs to an earlier pass).
#
# Finding #1 (CONFIRMED): the contextual lme4 summary read per-group variances
# through maihda_random_variances_lme4(), which returns only the INTERCEPT
# diagonal of each random-effect covariance block and carried no intercept-only
# validation. A contextual fit with (1 + x | stratum) therefore reported a
# partition whose total silently omitted the fitted slope variance and the
# intercept-slope covariance (repro: slope var 0.318 absent from a reported
# total of 0.812), while the SAME model without context = errored on the
# ordinary path's guard; a slope-only context (0 + x | site) returned an
# all-NA partition instead of an error. The extractor now enforces the same
# intercept-only contract as maihda_total_random_variance_lme4(), which also
# guards the crossed/contextual bootstraps, the crossed MOR, and the stepwise
# context share built on the same vector. The brms twin
# (maihda_group_variance_draws_brms) rejected compound slopes (two sd_ columns
# per group) but accepted a slope-only single column (sd_site__x), silently
# squaring SLOPE draws as the group's intercept variance; it now requires the
# single coefficient to be the Intercept.
#
# Finding #2 (CONFIRMED): maihda_cross_classified_formula() validated only the
# GROUPING side of each random-effect bar, then stripped all bars and inserted
# canonical intercepts, so y ~ x + (1 + x | stratum) was silently rewritten to
# y ~ x + (1 | dim_a) + (1 | dim_b) + (1 | stratum) -- no warning, wrong model
# -- through both maihda(decomposition = "crossed-dimensions") and the
# per-group crossed workflow. The builder now rejects a non-intercept
# left-hand side on an allowed grouping factor with a directed error, mirroring
# the existing extra-random-effect guard in the same function.

test_that("audit 2026-07-23b #1: lme4 per-group variance extractor rejects random slopes", {
  set.seed(2023)
  n <- 300
  d <- data.frame(
    g1 = sample(c("m", "f"), n, TRUE),
    g2 = sample(c("a", "b", "c"), n, TRUE),
    site = factor(sample(1:8, n, TRUE)),
    x = rnorm(n)
  )
  d$stratum <- interaction(d$g1, d$g2, sep = "_")
  d$y <- 1 + 0.4 * d$x + rnorm(n)

  # Compound intercept + slope on the stratum: the old extractor returned the
  # intercept diagonal only, silently dropping the slope variance/covariance.
  m_slope <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (1 + x | stratum) + (1 | site), data = d)))
  expect_error(maihda_random_variances_lme4(m_slope),
               "intercept-only random effects")

  # Slope-only term: the old extractor returned NA for the group, surfacing
  # downstream as an all-NA contextual partition rather than an error.
  m_slope_only <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (0 + x | site) + (1 | stratum), data = d)))
  expect_error(maihda_random_variances_lme4(m_slope_only),
               "intercept-only random effects")

  # Intercept-only models are untouched: full named vector, values matching the
  # VarCorr intercept diagonals.
  m_ok <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ x + (1 | stratum) + (1 | site), data = d)))
  v <- maihda_random_variances_lme4(m_ok)
  expect_setequal(names(v), c("stratum", "site"))
  vc <- lme4::VarCorr(m_ok)
  expect_equal(unname(v["stratum"]), as.numeric(vc$stratum[1, 1]))
  expect_equal(unname(v["site"]), as.numeric(vc$site[1, 1]))
})

test_that("audit 2026-07-23b #1: contextual summary errors on random slopes instead of silently dropping them", {
  set.seed(2024)
  n <- 400
  d <- data.frame(
    g1 = sample(c("m", "f"), n, TRUE),
    g2 = sample(c("lo", "hi"), n, TRUE),
    site = factor(sample(1:8, n, TRUE)),
    x = rnorm(n)
  )
  d$stratum <- interaction(d$g1, d$g2, sep = "_")
  K <- nlevels(d$stratum)
  u0 <- rnorm(K, 0, 0.5)
  u1 <- rnorm(K, 0, 0.7)
  vs <- rnorm(8, 0, 0.5)
  d$y <- 1 + 0.4 * d$x + u0[as.integer(d$stratum)] +
    u1[as.integer(d$stratum)] * d$x + vs[as.integer(d$site)] +
    rnorm(n, 0, 0.5)

  # The finding's headline case: compound slope + context. Before the fix this
  # returned a partition omitting the slope variance and covariance; the
  # ordinary (no-context) summary of the same model already errored.
  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 + x | stratum), data = d, context = "site")))
  expect_error(summary(fit), "intercept-only random effects")

  # Slope-only context term: before the fix, an all-NA partition.
  fit_so <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (0 + x | site) + (1 | stratum), data = d,
               context = "site")))
  expect_error(summary(fit_so), "intercept-only random effects")
  expect_error(summary(fit_so), "site")

  # Confinement: a valid intercept-only contextual fit still partitions, and
  # the reported total is exactly the sum of the VarCorr intercept variances
  # plus the residual (nothing dropped, nothing added).
  fit_ok <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, context = "site")))
  sm <- summary(fit_ok)
  comp <- sm$variance_components
  expect_true(all(is.finite(comp$variance)))
  vc <- as.data.frame(lme4::VarCorr(fit_ok$model))
  expected_total <- sum(vc$vcov)
  expect_equal(comp$variance[comp$component == "Total"], expected_total,
               tolerance = 1e-8)
  expect_equal(sm$context$vpc_stratum,
               unname(sm$context$var_stratum / expected_total),
               tolerance = 1e-8)
})

test_that("audit 2026-07-23b #1: brms per-group draws reject slope-only single-coefficient groups", {
  # Slope-only (0 + x | site): ONE sd column whose coefficient is not the
  # Intercept. The old parser squared the slope draws as site's intercept
  # variance; the compound two-column case was already rejected.
  draws_slope <- data.frame(
    sd_stratum__Intercept = c(0.5, 0.6),
    sd_site__x = c(0.3, 0.4),
    sigma = c(1, 1), check.names = FALSE)
  expect_error(maihda_group_variance_draws_brms(draws_slope), "slope-only")
  expect_error(maihda_group_variance_draws_brms(draws_slope), "sd_site__x")

  # The compound (multi-coefficient) rejection survives the loop restructure.
  draws_multi <- data.frame(
    sd_stratum__Intercept = c(0.5, 0.6),
    sd_stratum__x = c(0.2, 0.3),
    sigma = c(1, 1), check.names = FALSE)
  expect_error(maihda_group_variance_draws_brms(draws_multi),
               "multiple sd_stratum")

  # Intercept-only draws still parse, including a group name containing "__".
  draws_ok <- data.frame(
    `sd_site__id__Intercept` = c(0.7, 0.7),
    `sd_stratum__Intercept` = c(0.9, 0.8),
    sigma = c(1, 1), check.names = FALSE)
  gv <- maihda_group_variance_draws_brms(draws_ok)
  expect_setequal(names(gv), c("site__id", "stratum"))
  expect_equal(gv$stratum, draws_ok[["sd_stratum__Intercept"]]^2)

  # A column with no "__" separator has no coefficient name to check and keeps
  # the historical permissive fallback (never emitted by real brms fits).
  draws_nosep <- data.frame(
    sd_legacy = c(0.5, 0.5),
    sigma = c(1, 1), check.names = FALSE)
  gv2 <- maihda_group_variance_draws_brms(draws_nosep)
  expect_equal(names(gv2), "legacy")
})

test_that("audit 2026-07-23b #2: crossed builder rejects random slopes on allowed factors", {
  set.seed(2025)
  n <- 120
  d <- data.frame(
    a = factor(sample(c("F", "M"), n, TRUE)),
    b = factor(sample(c("lo", "hi"), n, TRUE)),
    x = rnorm(n)
  )
  d$stratum <- interaction(d$a, d$b, sep = "_")
  d$y <- rnorm(n)

  # Every slope spelling on an ALLOWED grouping factor: before the fix each was
  # silently rewritten to (1 | a) + (1 | b) + (1 | stratum).
  expect_error(
    maihda_cross_classified_formula(y ~ x + (1 + x | stratum), c("a", "b"),
                                    list(), d),
    "random-slope term")
  expect_error(
    maihda_cross_classified_formula(y ~ x + (0 + x | stratum), c("a", "b"),
                                    list(), d),
    "random-slope term")
  # (x | stratum) carries an implicit intercept alongside the slope.
  expect_error(
    maihda_cross_classified_formula(y ~ x + (x | stratum), c("a", "b"),
                                    list(), d),
    "random-slope term")
  # A slope bar next to a clean intercept bar, and a slope on a DIMENSION
  # factor rather than the intersection, are both caught.
  expect_error(
    maihda_cross_classified_formula(
      y ~ x + (1 | stratum) + (0 + x | stratum), c("a", "b"), list(), d),
    "random-slope term")
  expect_error(
    maihda_cross_classified_formula(
      y ~ x + (1 | stratum) + (0 + x | a), c("a", "b"), list(), d),
    "random-slope term")

  # Intercept-only input still builds the canonical crossed formula.
  cc <- maihda_cross_classified_formula(y ~ x + (1 | stratum), c("a", "b"),
                                        list(), d)
  rhs <- paste(deparse(cc$formula), collapse = " ")
  expect_true(grepl("(1 | a)", rhs, fixed = TRUE))
  expect_true(grepl("(1 | b)", rhs, fixed = TRUE))
  expect_true(grepl("(1 | stratum)", rhs, fixed = TRUE))

  # The pre-existing guard for an UNDECLARED grouping factor keeps its own
  # (different) error.
  expect_error(
    maihda_cross_classified_formula(
      y ~ x + (1 | stratum) + (1 | clinic), c("a", "b"), list(), d),
    "extra random effect")
})

test_that("audit 2026-07-23b #2: maihda(crossed-dimensions) surfaces the slope error end-to-end", {
  set.seed(2026)
  n <- 200
  d0 <- data.frame(
    a = factor(sample(c("F", "M"), n, TRUE)),
    b = factor(sample(c("lo", "hi"), n, TRUE)),
    x = rnorm(n)
  )
  s <- make_strata(d0, vars = c("a", "b"))
  d <- s$data
  d$y <- 1 + 0.4 * d$x + rnorm(n)

  expect_error(
    suppressWarnings(suppressMessages(
      maihda(y ~ x + (1 + x | stratum), data = d,
             decomposition = "crossed-dimensions"))),
    "random-slope term")

  # Confinement: the intercept-only crossed decomposition still runs and keeps
  # a coherent partition (shares summing to 1).
  res <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | stratum), data = d,
           decomposition = "crossed-dimensions")))
  dec <- res$summary$decomposition
  expect_equal(dec$additive_share + dec$interaction_share, 1, tolerance = 1e-8)
})
