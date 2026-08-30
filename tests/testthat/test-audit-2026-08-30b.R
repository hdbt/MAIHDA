# Audit 2026-08-30 (second pass of the day): the count VPC hard-coded the
# lognormal level-1 variance approximation, with no way to reach the delta or
# trigamma alternatives of Nakagawa, Johnson & Schielzeth (2017) and no report of
# which one produced the number -- which matters because the three diverge by a
# factor of six below a marginal count of 2, the regime a rare-outcome count
# MAIHDA lives in. fit_maihda(count_approximation = ) now selects; the default is
# unchanged.

test_that("the three approximations match Nakagawa et al. table 1 definitions", {
  # cv2 is Var(y | u) / mu^2: 1/lambda (Poisson), 1/lambda + 1/theta (negbin).
  lambda <- c(0.1, 0.5, 1, 2, 5, 20)
  cv2 <- 1 / lambda

  expect_equal(MAIHDA:::maihda_count_level1_transform(cv2, "lognormal"),
               log1p(1 / lambda))
  expect_equal(MAIHDA:::maihda_count_level1_transform(cv2, "delta"), 1 / lambda)
  # The trigamma form is evaluated at the RECIPROCAL of cv2 -- for the Poisson
  # that is lambda itself, exactly as insight::get_variance() computes it.
  expect_equal(MAIHDA:::maihda_count_level1_transform(cv2, "trigamma"),
               trigamma(lambda))

  # Negative binomial: trigamma((1/lambda + 1/theta)^-1), NOT trigamma(lambda).
  theta <- 2
  cv2_nb <- 1 / lambda + 1 / theta
  expect_equal(MAIHDA:::maihda_count_level1_transform(cv2_nb, "trigamma"),
               trigamma(1 / (1 / lambda + 1 / theta)))
  # ... which reduces onto the Poisson row as theta -> Inf.
  expect_equal(MAIHDA:::maihda_count_level1_transform(1 / lambda + 1 / 1e12,
                                                      "trigamma"),
               trigamma(lambda), tolerance = 1e-6)

  # Strict ordering for every lambda > 0: trigamma > delta > lognormal, so the
  # VPC is largest under the lognormal default and smallest under trigamma.
  expect_true(all(trigamma(lambda) > 1 / lambda))
  expect_true(all(1 / lambda > log1p(1 / lambda)))

  # The divergence the finding reported, pinned at the two ends.
  expect_equal(log1p(1 / 0.5) / trigamma(0.5), 0.2226254, tolerance = 1e-6)
  expect_equal(log1p(1 / 5) / trigamma(5), 0.8237806, tolerance = 1e-6)
})

test_that("maihda_count_level1_variance defaults to the pre-existing lognormal form", {
  mu <- c(0.2, 0.4, 1.5, 3)
  # Unchanged from before the argument existed: reduce to the mean count, THEN
  # transform (audit 2026-07-30).
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu), log1p(1 / mean(mu)))
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, theta = 2),
               log1p(1 / mean(mu) + 1 / 2))
  # Explicit "lognormal" is the same value; the alternatives are not.
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, approximation = "lognormal"),
               log1p(1 / mean(mu)))
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, approximation = "delta"),
               1 / mean(mu))
  expect_equal(MAIHDA:::maihda_count_level1_variance(mu, approximation = "trigamma"),
               trigamma(mean(mu)))
  # Weights still weight the mean COUNT under every approximation.
  w <- c(1, 1, 3, 5)
  expect_equal(
    MAIHDA:::maihda_count_level1_variance(mu, w = w, approximation = "trigamma"),
    trigamma(sum(w * mu) / sum(w)))
})

test_that("count_approximation is validated and defaults are backward compatible", {
  expect_identical(MAIHDA:::maihda_check_count_approximation(NULL), "lognormal")
  # A caller forwarding the whole c("lognormal", "delta", "trigamma") default.
  expect_identical(
    MAIHDA:::maihda_check_count_approximation(
      MAIHDA:::maihda_count_approximation_choices), "lognormal")
  expect_identical(MAIHDA:::maihda_check_count_approximation("trigamma"), "trigamma")
  expect_error(MAIHDA:::maihda_check_count_approximation("trigama"),
               "must be one of")
  expect_error(MAIHDA:::maihda_check_count_approximation(c("delta", "trigamma")),
               "must be one of")
  expect_error(MAIHDA:::maihda_check_count_approximation(3), "must be one of")
  # A model object fitted before the argument existed reads as the default,
  # rather than erroring on the missing slot.
  expect_identical(MAIHDA:::maihda_count_approximation(list(engine = "lme4")),
                   "lognormal")
  expect_identical(
    MAIHDA:::maihda_count_approximation(list(count_approximation = "delta")),
    "delta")
})

test_that("the per-draw brms count path honours the approximation too", {
  # Pure inputs, no Stan: eta_link is a draws x obs matrix.
  eta <- matrix(log(c(0.3, 0.5, 0.8, 1.2)), nrow = 2, ncol = 4, byrow = TRUE)
  vtot <- c(0.2, 0.4)
  lam <- rowMeans(exp(sweep(eta, 1L, vtot / 2, "+")))

  expect_equal(MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot),
               log1p(1 / lam))
  expect_equal(
    MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot,
                                                 approximation = "trigamma"),
    trigamma(lam))
  shape <- c(2, 3)
  expect_equal(
    MAIHDA:::maihda_count_resid_var_from_linpred(eta, vtot, extra = 1 / shape,
                                                 approximation = "trigamma"),
    trigamma(1 / (1 / lam + 1 / shape)))
})

# ---- end to end -----------------------------------------------------------

make_low_count_data <- function() {
  set.seed(20260830)
  d <- expand.grid(sex = c("F", "M"), edu = c("low", "mid", "high"),
                   race = c("A", "B", "C"), rep = seq_len(60))
  st <- interaction(d$sex, d$edu, d$race, drop = TRUE)
  u <- stats::rnorm(nlevels(st), 0, sqrt(0.25))
  d$y <- stats::rpois(nrow(d), exp(log(0.3) + u[as.integer(st)]))
  d$g <- stats::rnorm(nrow(d))
  d
}

test_that("the fit-time choice reaches the VPC, the components table and the note", {
  skip_on_cran()
  d <- make_low_count_data()
  f <- y ~ 1 + (1 | sex:edu:race)

  fits <- lapply(c("lognormal", "delta", "trigamma"), function(a) {
    m <- fit_maihda(f, data = d, family = poisson(), count_approximation = a)
    expect_identical(m$count_approximation, a)
    suppressWarnings(summary(m))
  })
  names(fits) <- c("lognormal", "delta", "trigamma")

  vpc <- vapply(fits, function(s) s$vpc$estimate, numeric(1))
  lvl <- vapply(fits, function(s) s$count_vpc$level1_variance, numeric(1))

  # Strict ordering (see the transform test): trigamma has the largest level-1
  # variance and therefore the smallest VPC.
  expect_lt(lvl[["lognormal"]], lvl[["delta"]])
  expect_lt(lvl[["delta"]], lvl[["trigamma"]])
  expect_gt(vpc[["lognormal"]], vpc[["delta"]])
  expect_gt(vpc[["delta"]], vpc[["trigamma"]])
  # Not a rounding difference: at this marginal count the choice moves the VPC by
  # more than a factor of five.
  expect_gt(vpc[["lognormal"]] / vpc[["trigamma"]], 5)

  # The components table carries the SAME level-1 variance the note reports, so
  # the printed decomposition and the note can never disagree.
  for (nm in names(fits)) {
    s <- fits[[nm]]
    resid_row <- s$variance_components$variance[
      s$variance_components$component == "Within-stratum (residual)"]
    expect_equal(resid_row, s$count_vpc$level1_variance)
    expect_identical(s$count_vpc$approximation, nm)
    expect_true(s$count_vpc$low_count)
    # All three alternatives are reported, whichever was used, and the used one
    # is one of them. (On lme4 the plug-in IS the number the VPC used; on brms
    # the summary works draw by draw and reports E[sigma^2_e], which differs
    # from sigma^2_e(E[lambda]) by Jensen -- which is why print() reports the
    # method and lambda rather than restating the variance beside the table.)
    expect_named(s$count_vpc$alternatives, c("lognormal", "delta", "trigamma"))
    expect_equal(unname(s$count_vpc$alternatives[[nm]]),
                 s$count_vpc$level1_variance)
  }

  # The default is still the lognormal number -- no silent change to any
  # published count VPC.
  m_default <- fit_maihda(f, data = d, family = poisson())
  expect_identical(m_default$count_approximation, "lognormal")
  expect_equal(suppressWarnings(summary(m_default))$vpc$estimate,
               vpc[["lognormal"]])
})

test_that("summary warns below the lambda = 2 threshold and is quiet above it", {
  skip_on_cran()
  d <- make_low_count_data()
  m_low <- fit_maihda(y ~ 1 + (1 | sex:edu:race), data = d, family = poisson())
  expect_warning(summary(m_low), "count VPC/ICC level-1 variance")
  expect_warning(summary(m_low), "agree only above 2")

  # A high-count fit on the same design must NOT warn.
  set.seed(99)
  d2 <- d
  st <- interaction(d2$sex, d2$edu, d2$race, drop = TRUE)
  u <- stats::rnorm(nlevels(st), 0, 0.5)
  d2$y <- stats::rpois(nrow(d2), exp(log(25) + u[as.integer(st)]))
  m_high <- fit_maihda(y ~ 1 + (1 | sex:edu:race), data = d2, family = poisson())
  s_high <- expect_silent(summary(m_high))
  expect_false(s_high$count_vpc$low_count)
  expect_gt(s_high$count_vpc$lambda, 2)
})

test_that("the null-model level-1 variance equals insight's null plug-in exactly", {
  skip_on_cran()
  # ?fit_maihda claims MAIHDA and insight::get_variance() agree EXACTLY on a null
  # model and diverge only on adjusted ones. Pin the exact half here so the claim
  # cannot rot: insight evaluates log1p(1/lambda) at a SINGLE lambda taken from an
  # intercept-only null model, exp(beta0 + sigma^2/2), while this package averages
  # the row-wise marginal counts exp(x_i'beta + v_i/2). With no covariate every
  # row-wise lambda IS that plug-in, so the mean of them is too.
  #
  # insight's estimator is reimplemented here rather than called: insight is not a
  # dependency (not even in Suggests), and the invariant is about OUR lambda, not
  # about their package being installed.
  d <- make_low_count_data()
  d$g <- interaction(d$sex, d$edu, d$race, drop = TRUE)
  m <- suppressWarnings(
    lme4::glmer(y ~ 1 + (1 | g), data = d, family = stats::poisson()))

  pkg <- MAIHDA:::maihda_residual_variance_lme4(m)
  b0 <- unname(lme4::fixef(m)[["(Intercept)"]])
  s2 <- unname(as.numeric(lme4::VarCorr(m)$g[1, 1]))
  insight_style <- log1p(1 / exp(b0 + s2 / 2))
  expect_equal(pkg, insight_style, tolerance = 1e-12)

  # ... and the divergence on an adjusted fit is real and SIGNED: Jensen makes the
  # averaged row-wise lambda the larger, so the package level-1 variance is the
  # smaller (and its VPC the larger). This is the half the docs must not call
  # equivalence.
  set.seed(4)
  d$x <- stats::rnorm(nrow(d))
  d$y2 <- stats::rpois(nrow(d), exp(log(2.5) + 1.5 * d$x +
                                      stats::rnorm(nlevels(d$g), 0, 0.5)[d$g]))
  m2 <- suppressWarnings(
    lme4::glmer(y2 ~ x + (1 | g), data = d, family = stats::poisson()))
  m2_null <- suppressWarnings(
    lme4::glmer(y2 ~ 1 + (1 | g), data = d, family = stats::poisson()))
  pkg2 <- MAIHDA:::maihda_residual_variance_lme4(m2)
  ins2 <- log1p(1 / exp(unname(lme4::fixef(m2_null)[["(Intercept)"]]) +
                          unname(as.numeric(lme4::VarCorr(m2_null)$g[1, 1])) / 2))
  expect_lt(pkg2, ins2)
  expect_false(isTRUE(all.equal(pkg2, ins2, tolerance = 1e-3)))
})

test_that("maihda() forwards count_approximation to both derived fits", {
  skip_on_cran()
  d <- make_low_count_data()
  # maihda() has no count_approximation formal -- it reaches fit_maihda through
  # `...`/dots_eval, and every derived null/adjusted fit must get the same one,
  # or a null and an adjusted VPC in the same analysis would sit on different
  # level-1 scales and the printed comparison would be meaningless.
  a <- suppressWarnings(
    maihda(y ~ sex + edu + race + (1 | sex:edu:race), data = d,
           family = "poisson", count_approximation = "trigamma",
           interactions = FALSE))
  expect_identical(a$model$count_approximation, "trigamma")
  expect_identical(a$model_adjusted$count_approximation, "trigamma")
  expect_identical(a$summary$count_vpc$approximation, "trigamma")
  expect_identical(a$summary_adjusted$count_vpc$approximation, "trigamma")

  b <- suppressWarnings(
    maihda(y ~ sex + edu + race + (1 | sex:edu:race), data = d,
           family = "poisson", interactions = FALSE))
  expect_identical(b$model$count_approximation, "lognormal")
  expect_gt(b$summary$vpc$estimate, a$summary$vpc$estimate)
  # The PCV is a ratio of BETWEEN-stratum variances, so it must not move at all.
  expect_equal(a$pcv$pcv, b$pcv$pcv)
})

test_that("a non-count family is untouched and carries no note", {
  skip_on_cran()
  d <- make_low_count_data()
  f <- g ~ 1 + (1 | sex:edu:race)
  s_default <- summary(fit_maihda(f, data = d))
  s_trigamma <- summary(fit_maihda(f, data = d, count_approximation = "trigamma"))
  # Bit-identical: the argument is inert for a Gaussian fit.
  expect_identical(s_default$vpc$estimate, s_trigamma$vpc$estimate)
  expect_null(s_default$count_vpc)
  expect_null(s_trigamma$count_vpc)
})

test_that("the bootstrap interval uses the model's own approximation", {
  skip_on_cran()
  d <- make_low_count_data()
  m <- fit_maihda(y ~ 1 + (1 | sex:edu:race), data = d, family = poisson())
  # maihda_bootstrap_ci() enforces an absolute floor of 10 successful refits.
  set.seed(7)
  ln <- bootstrap_vpc(m$model, m$data, m$formula, n_boot = 12, conf_level = 0.95,
                      approximation = "lognormal")
  set.seed(7)
  tg <- bootstrap_vpc(m$model, m$data, m$formula, n_boot = 12, conf_level = 0.95,
                      approximation = "trigamma")
  # Same simulated responses, same refits: the ONLY difference is the level-1
  # transform. Every draw's lognormal VPC exceeds its trigamma VPC (the level-1
  # ordering is strict), and element-wise dominance carries to the order
  # statistics, so BOTH interval endpoints must move. maihda_bootstrap_ci()
  # returns a bare c(lower, upper). If the parameter never reached the refit
  # loop these would be identical.
  expect_gt(as.numeric(ln)[1], as.numeric(tg)[1])
  expect_gt(as.numeric(ln)[2], as.numeric(tg)[2])
})

test_that("the low-count flag uses the EFFECTIVE count, so overdispersion counts", {
  # A negative binomial can sit far above lambda = 2 and still be deep inside the
  # divergent regime, because (1/lambda + 1/theta)^-1 is what the approximations
  # are indexed by. A bare lambda > 2 test would wave this through.
  lambda <- 12
  theta <- 0.5
  eff <- 1 / (1 / lambda + 1 / theta)
  expect_lt(eff, 2)
  expect_gt(lambda, 2)
  ln <- MAIHDA:::maihda_count_level1_transform(1 / lambda + 1 / theta, "lognormal")
  tg <- MAIHDA:::maihda_count_level1_transform(1 / lambda + 1 / theta, "trigamma")
  expect_lt(ln / tg, 0.25)
})
