# Audit 2026-09-15 (d): external audit-2026-09-14 finding 7 [P2] -- "Information
# criteria are incorrectly declared incomparable across families or links".
#
# maihda_ic() withheld the delta whenever two models' "family(link)" keys differed:
# the strict equality the VPC and PCV need, whose level-1 variance depends on the
# family and link. An information criterion compares log-likelihoods of the same
# observations, so it needs a common response measure and complete normalising
# constants, not a common family. The auditor's F5 reproduced exactly (Poisson AIC
# 2932.0088 against negative binomial 2677.5431, lme4's own delta 254.4656516
# withheld), binomial logit / probit / cloglog fits were withheld the same way, and
# so was the ordinal twin the report did not name, clmm logit vs probit. Both count
# likelihoods matched a direct integration of the marginal likelihood to Laplace
# error (0.018, 0.014) against a log(y!) constant of 1236, so they are comparable.
# A delta is now reported within three response classes -- counts (Poisson, negative
# binomial), binomial / Bernoulli outcomes, ordered categories (cumulative) -- and
# every other family or link difference still withholds it: a Gaussian density is
# not comparable with a probability, and glmer's non-identity Gaussian likelihood
# is not on lmer()'s scale (it exceeded the maximum marginal log-likelihood by 15.4
# on the fit this pass measured; test-audit-2026-09-15e.R now reports no criteria
# for such fits, and restores the term lme4 leaves out under nAGQ > 1, which would
# otherwise have made a Poisson-vs-negative-binomial delta under quadrature wrong).
#
# Reporting those deltas exposed a second, pre-existing defect: lme4's parameter
# count adds one for theta in EVERY negative-binomial family, including a fixed
# MASS::negative.binomial(theta) passed to glmer(), which estimates no theta. Such a
# fit's AIC was 2 too high and its BIC log(n) too high, so a delta against a Poisson
# or estimated-theta fit would have carried a parameter that does not exist. Per the
# user, maihda_ic() counts only the estimated parameters, as glm() does for that
# family.

f7_data <- function(seed = 20260915) {
  set.seed(seed)
  d <- data.frame(stratum = factor(rep(seq_len(24), each = 30)))
  n <- nrow(d)
  d$x <- rnorm(n)
  u <- rnorm(24, sd = 0.5)[d$stratum]
  d$cnt <- rnbinom(n, mu = exp(0.4 + 0.3 * d$x + u), size = 2)
  d$yb <- rbinom(n, 1, plogis(-0.3 + 0.6 * d$x + u))
  lat <- 0.5 * d$x + u + rlogis(n)
  d$yo <- factor(cut(lat, c(-Inf, -0.7, 0.9, Inf), labels = c("lo", "mid", "hi")),
                 levels = c("lo", "mid", "hi"), ordered = TRUE)
  d
}

f7_fit <- function(...) suppressWarnings(suppressMessages(fit_maihda(...)))

# Marginal log-likelihood of a random-intercept GLMM, integrated stratum by stratum
# with adaptive Gauss-Hermite quadrature at the fit's own estimates: an oracle that
# owes nothing to lme4's likelihood code.
f7_marginal_loglik <- function(y, eta, g, tau, logf, n_nodes = 30L) {
  i <- seq_len(n_nodes - 1L)
  jacobi <- matrix(0, n_nodes, n_nodes)
  jacobi[cbind(i, i + 1L)] <- sqrt(i / 2)
  jacobi[cbind(i + 1L, i)] <- sqrt(i / 2)
  eig <- eigen(jacobi, symmetric = TRUE)
  nodes0 <- eig$values
  weights0 <- sqrt(pi) * eig$vectors[1, ]^2
  total <- 0
  for (lev in unique(g)) {
    idx <- which(g == lev)
    logh <- function(u) {
      sum(logf(y[idx], eta[idx] + u)) + stats::dnorm(u, 0, tau, log = TRUE)
    }
    opt <- stats::optimize(logh, c(-20, 20), maximum = TRUE, tol = 1e-10)
    h <- 1e-4
    curv <- (logh(opt$maximum + h) - 2 * opt$objective + logh(opt$maximum - h)) / h^2
    s <- 1 / sqrt(-curv)
    nodes <- opt$maximum + sqrt(2) * s * nodes0
    vals <- vapply(nodes, logh, numeric(1)) - opt$objective + nodes0^2
    total <- total + opt$objective + log(sqrt(2) * s * sum(weights0 * exp(vals)))
  }
  total
}

test_that("maihda_ic() reports the delta between Poisson and negative-binomial fits of the same counts", {
  skip_on_cran()
  d <- f7_data()
  mp <- f7_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson")
  mn <- f7_fit(cnt ~ x + (1 | stratum), data = d, family = "negbinomial")

  expect_no_warning(ic <- maihda_ic(mp, mn))
  expect_true("delta" %in% names(ic))
  aic <- c(stats::AIC(mp$model), stats::AIC(mn$model))
  expect_identical(ic$delta, aic - min(aic))
  expect_gt(max(ic$delta), 10)

  # The printed rule no longer says the family must match.
  expect_output(print(ic), "unless all are count \\(Poisson / negative binomial\\)")

  # The premise: both are complete mass functions of the same counts. Each lme4
  # log-likelihood matches the integrated oracle to Laplace error, far inside the
  # log(y!) constant a likelihood on another normalisation would differ by.
  tau <- function(fm) sqrt(as.numeric(lme4::VarCorr(fm)$stratum[1, 1]))
  eta <- function(fm) as.numeric(stats::predict(fm, re.form = NA, type = "link"))
  theta <- maihda_negbin_theta_lme4(mn$model)
  ll_p <- f7_marginal_loglik(d$cnt, eta(mp$model), d$stratum, tau(mp$model),
                             function(y, e) stats::dpois(y, exp(e), log = TRUE))
  ll_n <- f7_marginal_loglik(d$cnt, eta(mn$model), d$stratum, tau(mn$model),
                             function(y, e) stats::dnbinom(y, size = theta, mu = exp(e), log = TRUE))
  expect_gt(sum(lgamma(d$cnt + 1)), 100)
  expect_lt(abs(ic$logLik[1] - ll_p), 0.1)
  expect_lt(abs(ic$logLik[2] - ll_n), 0.1)
})

test_that("maihda_ic() reports the delta between binomial fits that differ only in their link", {
  skip_on_cran()
  d <- f7_data()
  fits <- lapply(c("logit", "probit", "cloglog"), function(link) {
    f7_fit(yb ~ x + (1 | stratum), data = d, family = stats::binomial(link))
  })

  expect_no_warning(ic <- do.call(maihda_ic, fits))
  expect_true("delta" %in% names(ic))
  aic <- vapply(fits, function(m) stats::AIC(m$model), numeric(1))
  expect_identical(ic$delta, aic - min(aic))
})

test_that("maihda_ic() reports the delta between cumulative fits that differ only in their link", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- f7_data()
  m_logit <- f7_fit(yo ~ x + (1 | stratum), data = d, engine = "ordinal",
                    family = maihda_cumulative("logit"))
  m_probit <- f7_fit(yo ~ x + (1 | stratum), data = d, engine = "ordinal",
                     family = maihda_cumulative("probit"))

  expect_no_warning(ic <- maihda_ic(m_logit, m_probit))
  expect_true("delta" %in% names(ic))
  aic <- c(stats::AIC(m_logit$model), stats::AIC(m_probit$model))
  expect_identical(ic$delta, aic - min(aic))
})

test_that("a fixed negative-binomial theta is not counted as an estimated parameter", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  d <- f7_data()
  n <- nrow(d)
  models <- list(
    f7_fit(cnt ~ x + (1 | stratum), data = d, family = MASS::negative.binomial(2)),
    f7_fit(cnt ~ x + (1 | stratum), data = d, family = MASS::negative.binomial(5)),
    f7_fit(cnt ~ x + (1 | stratum), data = d, family = "negbinomial"),
    f7_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson")
  )

  expect_no_warning(ic <- do.call(maihda_ic, models))
  # Estimated parameters counted by hand: two fixed effects and the stratum SD, plus
  # theta only where glmer.nb() estimates it.
  k <- c(3, 3, 4, 3)
  ll <- vapply(models, function(m) as.numeric(stats::logLik(m$model)), numeric(1))
  hand_aic <- -2 * ll + 2 * k
  expect_identical(as.numeric(ic$df), k)
  expect_identical(ic$logLik, ll)
  expect_equal(ic$AIC, hand_aic, tolerance = 1e-12)
  expect_equal(ic$BIC, -2 * ll + log(n) * k, tolerance = 1e-12)
  expect_equal(ic$delta, hand_aic - min(hand_aic), tolerance = 1e-12)
  # The estimated-theta and Poisson rows are lme4's own criteria, untouched.
  expect_identical(ic$AIC[3:4], c(stats::AIC(models[[3]]$model), stats::AIC(models[[4]]$model)))
  expect_identical(ic$BIC[3:4], c(stats::BIC(models[[3]]$model), stats::BIC(models[[4]]$model)))

  # compare_maihda() appends the same criteria (its warnings concern the VPCs).
  cmp <- suppressWarnings(compare_maihda(models[[1]], models[[3]]))
  expect_identical(cmp$AIC, ic$AIC[c(1, 3)])
  expect_identical(cmp$BIC, ic$BIC[c(1, 3)])

  # Defensive: the correction applies only while lme4's count is one above the
  # estimated parameters, so an lme4 that stops counting a fixed theta is not
  # corrected twice. A Poisson fit recorded as fixed-theta has no extra count.
  stub <- models[[4]]
  stub$family <- MASS::negative.binomial(2)
  expect_true(maihda_negbin_theta_is_fixed(stub))
  expect_identical(maihda_ic_one(stub)$df, attr(stats::logLik(models[[4]]$model), "df"))
})

test_that("the delta is still withheld across response classes and for continuous families", {
  skip_on_cran()
  d <- f7_data()
  pairs <- list(
    # a Gaussian density and a Poisson probability of the same counts
    list(f7_fit(cnt ~ x + (1 | stratum), data = d, family = "gaussian"),
         f7_fit(cnt ~ x + (1 | stratum), data = d, family = "poisson")),
    # a binomial and a Poisson fit of the same 0/1 outcome
    list(f7_fit(yb ~ x + (1 | stratum), data = d, family = "binomial"),
         f7_fit(yb ~ x + (1 | stratum), data = d, family = "poisson"))
  )
  # (A continuous family whose link differs is pinned with brms-style stubs below:
  # a glmer fit of a Gaussian with a non-identity link reports no criteria at all,
  # see test-audit-2026-09-15e.R.)
  for (pair in pairs) {
    expect_warning(ic <- maihda_ic(pair[[1]], pair[[2]]), "families/links")
    expect_false("delta" %in% names(ic))
    expect_true(all(is.finite(ic$AIC)))
  }
})

test_that("brms family names join the same classes, and no class crosses the likelihood/Bayesian divide", {
  d <- data.frame(stratum = factor(rep(1:5, each = 10)), x = seq_len(50) / 10,
                  cnt = rep(0:4, 10), yb = rep(0:1, 25))
  stub <- function(family, link, response) {
    structure(list(model = structure(list(), class = "f7_stub_fit"),
                   family = list(family = family, link = link),
                   formula = stats::as.formula(paste(response, "~ x + (1 | stratum)")),
                   data = d),
              class = "maihda_model")
  }
  ic_row <- function(aic = NA_real_, looic = NA_real_) {
    data.frame(n = 50L, estimator = if (is.na(looic)) "ML" else "Bayesian",
               df = NA_real_, logLik = NA_real_, AIC = aic, BIC = aic + 5,
               WAIC = looic - 1, LOOIC = looic, stringsAsFactors = FALSE)
  }
  rows <- list()
  local_mocked_bindings(maihda_ic_one = function(model, ml = FALSE) {
    out <- rows[[1]]
    rows <<- rows[-1]
    out
  })

  rows <- list(ic_row(looic = 210), ic_row(looic = 200))
  expect_no_warning(ic <- maihda_ic(stub("poisson", "log", "cnt"), stub("negbinomial", "log", "cnt")))
  expect_identical(attr(ic, "ic_primary"), "LOOIC")
  expect_equal(ic$delta, c(10, 0))

  rows <- list(ic_row(looic = 150), ic_row(looic = 153))
  expect_no_warning(ic <- maihda_ic(stub("bernoulli", "logit", "yb"), stub("bernoulli", "probit", "yb")))
  expect_equal(ic$delta, c(0, 3))

  rows <- list(ic_row(looic = 210), ic_row(looic = 200))
  expect_warning(ic <- maihda_ic(stub("gaussian", "identity", "cnt"), stub("poisson", "log", "cnt")),
                 "families/links")
  expect_false("delta" %in% names(ic))

  # A continuous family keeps strict equality: the same Gaussian outcome under
  # another link withholds the delta.
  rows <- list(ic_row(looic = 210), ic_row(looic = 200))
  expect_warning(ic <- maihda_ic(stub("gaussian", "identity", "cnt"), stub("gaussian", "log", "cnt")),
                 "families/links")
  expect_false("delta" %in% names(ic))

  # Compatible families on different scales: the scale check still withholds it.
  rows <- list(ic_row(aic = 205), ic_row(looic = 200))
  expect_warning(ic <- maihda_ic(stub("poisson", "log", "cnt"), stub("negbinomial", "log", "cnt")),
                 "likelihood/Bayesian divide")
  expect_false("delta" %in% names(ic))
})

test_that("maihda_ic_response_class() groups only the families with a common response measure", {
  stub <- function(family, link = "log") {
    structure(list(model = structure(list(), class = "f7_stub_fit"),
                   family = list(family = family, link = link)),
              class = "maihda_model")
  }
  expect_identical(maihda_ic_response_class(stub("poisson")), "count")
  expect_identical(maihda_ic_response_class(stub("negbinomial")), "count")
  expect_identical(maihda_ic_response_class(stub("Negative Binomial(2)")), "count")
  expect_identical(maihda_ic_response_class(stub("binomial", "probit")), "binomial")
  expect_identical(maihda_ic_response_class(stub("bernoulli", "logit")), "binomial")
  expect_identical(maihda_ic_response_class(stub("cumulative", "probit")), "ordered categories")
  expect_identical(maihda_ic_response_class(stub("ordinal", "logit")), "ordered categories")
  for (fam in c("gaussian", "Gamma", "inverse.gaussian", "zero_inflated_poisson",
                "quasipoisson", "sratio")) {
    expect_identical(maihda_ic_response_class(stub(fam)), NA_character_)
  }
  no_family <- structure(list(model = structure(list(), class = "f7_stub_fit")),
                         class = "maihda_model")
  expect_identical(maihda_ic_response_class(no_family), NA_character_)
})
