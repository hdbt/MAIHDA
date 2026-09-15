# Audit 2026-09-15 (e): the lme4 likelihoods maihda_ic() reads -- a follow-up to
# external audit-2026-09-14 finding 7 (test-audit-2026-09-15d.R). Three defects, each
# measured against a direct integration of the marginal likelihood:
#
# 1. A glmer fit of a family with a scale parameter (a Gaussian with a non-identity
#    link, Gamma, inverse Gaussian) reports a logLik() that is not the marginal
#    likelihood. Maximising the integrated likelihood over every parameter, lme4's
#    value was ABOVE the achievable maximum by 15.4 / 6.4 / 12.8, and nested AIC
#    differences were off by 45.2 / 4.2 / 4.0. lme4's help pages do not mention it.
#    Per the user such fits report NA criteria.
# 2. Under adaptive quadrature (nAGQ > 1) glmer's logLik() leaves out the saturated
#    log-likelihood sum(log p(y | mu = y)): 675.8 on this 720-row Poisson outcome,
#    640.3 for the aggregated binomial one, 923.2 for a negative binomial with theta
#    1.6, 0 for a Bernoulli outcome. A Laplace and a quadrature fit of one model were
#    1351.6 AIC apart, and the Poisson-vs-negative-binomial delta that finding 7 made
#    reportable would have read 751.4 instead of 256.5 under quadrature (theta fixed
#    at 1.6). Per the user the term is restored, from lme4's own response module.
# 3. glmer.nb() chooses theta on that incomplete likelihood, so family = "negbinomial"
#    with nAGQ > 1 silently returned theta 0.087 (1.62 at nAGQ = 1) and a stratum SD
#    of 0, lme4 itself giving no warning. Per the user fit_maihda() refuses it.

e_data <- function() {
  set.seed(20260914)
  d <- data.frame(stratum = factor(rep(seq_len(24), each = 30)))
  n <- nrow(d)
  d$x <- rnorm(n)
  set.seed(17)
  d$cnt <- rnbinom(n, mu = exp(0.4 + 0.3 * d$x + rnorm(24, sd = 0.4)[d$stratum]), size = 2)
  set.seed(99)
  d$yb <- rbinom(n, 1, plogis(-0.3 + 0.6 * d$x + rnorm(24, sd = 0.6)[d$stratum]))
  set.seed(21)
  d$trials <- sample(3:8, n, TRUE)
  d$succ <- rbinom(n, d$trials, plogis(-0.2 + 0.5 * d$x + rnorm(24, sd = 0.5)[d$stratum]))
  d$fail <- d$trials - d$succ
  set.seed(7)
  d$yg <- exp(1 + 0.3 * d$x + rnorm(24, sd = 0.3)[d$stratum]) + rnorm(n, sd = 0.5)
  d
}

e_fit <- function(...) suppressWarnings(suppressMessages(fit_maihda(...)))

# Marginal log-likelihood of a random-intercept GLMM by adaptive Gauss-Hermite
# quadrature, stratum by stratum, at the fit's own estimates.
e_marginal_loglik <- function(y, eta, g, tau, logf, n_nodes = 30L) {
  i <- seq_len(n_nodes - 1L)
  jacobi <- matrix(0, n_nodes, n_nodes)
  jacobi[cbind(i, i + 1L)] <- sqrt(i / 2)
  jacobi[cbind(i + 1L, i)] <- sqrt(i / 2)
  eig <- eigen(jacobi, symmetric = TRUE)
  nodes0 <- eig$values
  weights0 <- sqrt(pi) * eig$vectors[1, ]^2
  total <- 0
  for (idx in split(seq_along(y), g)) {
    logh <- function(u) {
      sum(logf(y[idx], eta[idx] + u)) + stats::dnorm(u, 0, tau, log = TRUE)
    }
    opt <- stats::optimize(logh, c(-15, 15), maximum = TRUE, tol = 1e-10)
    h <- 1e-4
    curv <- (logh(opt$maximum + h) - 2 * opt$objective + logh(opt$maximum - h)) / h^2
    s <- 1 / sqrt(-curv)
    nodes <- opt$maximum + sqrt(2) * s * nodes0
    vals <- vapply(nodes, logh, numeric(1)) - opt$objective + nodes0^2
    total <- total + opt$objective + log(sqrt(2) * s * sum(weights0 * exp(vals)))
  }
  total
}

test_that("a glmer fit of a family with a scale parameter reports no information criteria", {
  skip_on_cran()
  d <- e_data()
  label <- "ML (glmer scale family: no likelihood)"
  m_log1 <- e_fit(yg ~ x + (1 | stratum), data = d, family = stats::gaussian("log"))
  m_log0 <- e_fit(yg ~ 1 + (1 | stratum), data = d, family = stats::gaussian("log"))
  m_gamma <- e_fit(yg ~ x + (1 | stratum), data = d, family = stats::Gamma("log"))

  for (m in list(m_log1, m_gamma)) {
    expect_true(maihda_glmer_has_scale(m$model))
    expect_no_warning(ic <- maihda_ic(m))
    expect_identical(ic$estimator, label)
    expect_false(any(c("df", "logLik", "AIC", "BIC") %in% names(ic)))
  }
  # Two such fits: nothing to rank, so no delta and nothing to warn about.
  expect_no_warning(ic <- maihda_ic(m_log0, m_log1))
  expect_false(any(c("AIC", "delta") %in% names(ic)))

  # Beside another model its row is NA while the other keeps lme4's criteria.
  m_pois <- e_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson")
  ic <- suppressWarnings(maihda_ic(m_log1, m_pois))
  expect_true(is.na(ic$AIC[1]) && is.na(ic$logLik[1]) && is.na(ic$df[1]))
  expect_identical(ic$AIC[2], stats::AIC(m_pois$model))

  # Families without a scale parameter, and lmer() (a scale, but not a GLMM), are
  # unaffected.
  m_lmer <- e_fit(yg ~ x + (1 | stratum), data = d, family = "gaussian")
  expect_false(maihda_glmer_has_scale(m_pois$model))
  expect_false(maihda_glmer_has_scale(m_lmer$model))
  expect_identical(maihda_ic(m_lmer)$AIC, stats::AIC(m_lmer$model))
})

test_that("maihda_ic() restores the saturated log-likelihood glmer leaves out under nAGQ > 1", {
  skip_on_cran()
  d <- e_data()
  m1 <- e_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson")
  m5 <- e_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson", nAGQ = 5)
  saturated <- sum(stats::dpois(d$cnt, d$cnt, log = TRUE))

  ic5 <- maihda_ic(m5)
  expect_equal(ic5$logLik, as.numeric(stats::logLik(m5$model)) + saturated, tolerance = 1e-12)
  expect_equal(ic5$AIC, stats::AIC(m5$model) - 2 * saturated, tolerance = 1e-12)
  expect_equal(ic5$BIC, stats::BIC(m5$model) - 2 * saturated, tolerance = 1e-12)

  # The restored value is the likelihood: it matches the integrated oracle to
  # quadrature error, where lme4's own value is hundreds of units away.
  tau <- sqrt(as.numeric(lme4::VarCorr(m5$model)$stratum[1, 1]))
  eta <- as.numeric(stats::predict(m5$model, re.form = NA, type = "link"))
  oracle <- e_marginal_loglik(d$cnt, eta, d$stratum, tau,
                              function(y, e) stats::dpois(y, exp(e), log = TRUE))
  expect_lt(abs(ic5$logLik - oracle), 1e-3)
  expect_gt(abs(as.numeric(stats::logLik(m5$model)) - oracle), 600)

  # A Laplace and a quadrature fit of the same model now agree to approximation error.
  expect_no_warning(ic <- maihda_ic(m1, m5))
  expect_lt(max(ic$delta), 0.1)

  # Aggregated binomial: the saturated term of the binomial mass function.
  mb5 <- e_fit(cbind(succ, fail) ~ x + (1 | stratum), data = d, family = "binomial", nAGQ = 5)
  expect_equal(maihda_ic(mb5)$logLik,
               as.numeric(stats::logLik(mb5$model)) +
                 sum(stats::dbinom(d$succ, d$trials, d$succ / d$trials, log = TRUE)),
               tolerance = 1e-12)

  # Bernoulli outcomes have a saturated likelihood of 1: criteria exactly lme4's.
  mbb5 <- e_fit(yb ~ x + (1 | stratum), data = d, family = "binomial", nAGQ = 5)
  expect_identical(maihda_glmer_saturated_loglik(mbb5$model), 0)
  ic_b <- maihda_ic(mbb5)
  expect_identical(ic_b$logLik, as.numeric(stats::logLik(mbb5$model)))
  expect_identical(ic_b$AIC, stats::AIC(mbb5$model))
  expect_identical(ic_b$BIC, stats::BIC(mbb5$model))
})

test_that("a Poisson-vs-negative-binomial delta under quadrature matches the Laplace one", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  d <- e_data()
  laplace <- maihda_ic(
    e_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson"),
    e_fit(cnt ~ x + (1 | stratum), data = d, family = MASS::negative.binomial(1.6)))
  expect_no_warning(quadrature <- maihda_ic(
    e_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson", nAGQ = 5),
    e_fit(cnt ~ x + (1 | stratum), data = d, family = MASS::negative.binomial(1.6), nAGQ = 5)))
  expect_gt(laplace$delta[1], 100)
  expect_lt(abs(quadrature$delta[1] - laplace$delta[1]), 0.05)
  expect_identical(as.numeric(quadrature$df), c(3, 3))
})

test_that("fit_maihda() refuses an estimated negative-binomial theta with nAGQ > 1", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  d <- e_data()
  msg <- "cannot estimate the negative-binomial theta with nAGQ > 1"
  expect_error(fit_maihda(cnt ~ x + (1 | stratum), data = d, family = "negbinomial", nAGQ = 5),
               msg)
  expect_error(fit_maihda(cnt ~ x + (1 | stratum), data = d, family = "negbinomial", nAGQ = 2L),
               msg)
  # A partial spelling binds to glmer()'s nAGQ, so it is refused too.
  expect_error(fit_maihda(cnt ~ x + (1 | stratum), data = d, family = "negbinomial", nAG = 5),
               msg)
  expect_error(fit_maihda(cnt ~ x + (1 | stratum), data = d, family = "negbinomial",
                          nAGQ = 1, nAG = 5),
               "supplied more than once")

  da <- d
  da$ga <- c("a", "b")[(as.integer(da$stratum) - 1L) %/% 12L + 1L]
  da$gb <- LETTERS[(as.integer(da$stratum) - 1L) %% 12L + 1L]
  da$stratum <- NULL
  expect_error(suppressMessages(maihda(cnt ~ x + ga + gb + (1 | ga:gb), data = da,
                                       family = "negbinomial", nAGQ = 5)),
               msg)

  # The Laplace and nAGQ = 0 bases estimate theta correctly, and a fixed theta has
  # nothing to estimate.
  for (q in c(0L, 1L)) {
    m <- e_fit(cnt ~ x + (1 | stratum), data = d, family = "negbinomial", nAGQ = q)
    expect_equal(maihda_negbin_theta_lme4(m$model), 1.6164, tolerance = 1e-2)
  }
  expect_no_error(e_fit(cnt ~ x + (1 | stratum), data = d,
                        family = MASS::negative.binomial(1.6), nAGQ = 5))
})
