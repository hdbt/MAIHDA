# Count-scale VPC of Leckie, Browne, Goldstein, Merlo & Austin (2020), "Partitioning
# variation in multilevel models for count data", Psychological Methods 25(6),
# 787-801: the closed form maihda_vpc_response() returns for a Poisson or
# negative-binomial fit.
#
# For u ~ N(0, s2), a conditional mean mu = exp(eta + u) and
# Var(y | u) = mu + alpha * mu^2, with lambda = exp(eta + s2 / 2) the marginal
# expected count:
#   level 2 = lambda^2 * (exp(s2) - 1)
#   level 1 = lambda + alpha * lambda^2 * exp(s2)
#   VPC     = level 2 / (level 2 + level 1)          (their eq. 8, 13 and 14)
# The oracle values are the paper's own Model 1 and hand evaluations of those
# expressions; a brute-force block checks the expressions themselves against a
# simulation of the model.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this
# file leaves the random stream the later test files read as it found it (see
# setup-seed.R).
count_vpc_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

# The expressions written out directly: the independent reference.
count_vpc_by_hand <- function(lambda, s2, alpha = 0) {
  level2 <- lambda^2 * (exp(s2) - 1)
  level2 / (level2 + lambda + alpha * lambda^2 * exp(s2))
}

count_vpc_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# Twelve strata (gender x race x edu) with a stratum effect of SD 0.5 on the log
# scale. `y` has no covariate and no offset, `yx` has both, `ynb` is overdispersed.
count_vpc_data <- function() {
  count_vpc_with_seed(20261011, {
    n <- 600
    d <- data.frame(
      gender = sample(c("F", "M"), n, replace = TRUE),
      race   = sample(c("A", "B", "C"), n, replace = TRUE),
      edu    = sample(c("lo", "hi"), n, replace = TRUE),
      age    = stats::rnorm(n),
      expo   = stats::runif(n, 0.5, 4),
      school = factor(sample(1:10, n, replace = TRUE)),
      stringsAsFactors = FALSE
    )
    sk <- interaction(d$gender, d$race, d$edu, drop = TRUE)
    u <- stats::rnorm(nlevels(sk), sd = 0.5)[sk]
    d$y   <- stats::rpois(n, exp(1.2 + u))
    d$yx  <- stats::rpois(n, d$expo * exp(0.3 + 0.25 * d$age + u))
    d$ynb <- stats::rnbinom(n, mu = exp(1.2 + 0.2 * d$age + u), size = 2)
    d
  })
}

# Each fixture is fitted once, on first use, and shared by the blocks below.
count_vpc_fits <- new.env()
count_vpc_fit <- function(name) {
  if (is.null(count_vpc_fits[[name]])) {
    d <- count_vpc_data()
    count_vpc_fits[[name]] <- count_vpc_quiet(switch(
      name,
      null = fit_maihda(y ~ (1 | gender:race:edu), data = d, family = "poisson"),
      adjusted = fit_maihda(yx ~ age + offset(log(expo)) + (1 | gender:race:edu),
                            data = d, family = "poisson"),
      negbin = fit_maihda(ynb ~ age + (1 | gender:race:edu), data = d,
                          family = "negbinomial"),
      context = fit_maihda(y ~ (1 | gender:race:edu), data = d, family = "poisson",
                           context = "school"),
      stop("unknown fixture: ", name)
    ))
  }
  count_vpc_fits[[name]]
}

# ---- the closed form (no fit required) ---------------------------------------

test_that("the closed form reproduces the published and hand-computed values", {
  # Model 1 of the paper: Poisson, intercept 2.085, cluster variance 0.100. It
  # reports a marginal expectation of 8.46, a level-2 variance of 7.52, a level-1
  # variance of 8.46 and a VPC of 0.47.
  m1 <- maihda_count_vpc_exact(exp(2.085 + 0.100 / 2), 0.100)
  expect_equal(round(m1$lambda, 2), 8.46)
  expect_equal(round(m1$level2, 2), 7.52)
  expect_equal(round(m1$level1, 2), 8.46)
  expect_equal(round(m1$vpc, 3), 0.471)
  expect_equal(round(m1$vpc, 2), 0.47)

  # Poisson, variance 0.2, at a low, a moderate and a high expected count.
  pois <- maihda_count_vpc_exact(c(0.34, 2, 10), 0.2)
  expect_equal(round(pois$vpc, 4), c(0.0700, 0.3069, 0.6889))

  # Negative binomial in its mean-dispersion form, Var(y | u) = mu + alpha * mu^2.
  expect_equal(round(maihda_count_vpc_exact(2, 0.2, alpha = 0.5)$vpc, 4), 0.1662)
  expect_equal(round(maihda_count_vpc_exact(10, 1, alpha = 0.5)$vpc, 4), 0.5408)

  # The two components are the VPC's numerator and the rest of its denominator,
  # and they add up to the marginal variance of eq. 13-14,
  # lambda + lambda^2 * (exp(s2) * (1 + alpha) - 1).
  nb <- maihda_count_vpc_exact(c(0.5, 2, 10), 0.7, alpha = 0.4)
  expect_equal(nb$vpc, nb$level2 / (nb$level2 + nb$level1))
  expect_equal(nb$level2 + nb$level1,
               nb$lambda + nb$lambda^2 * (exp(0.7) * (1 + 0.4) - 1))
  expect_equal(nb$vpc, count_vpc_by_hand(c(0.5, 2, 10), 0.7, alpha = 0.4))
})

test_that("the closed form is a partition of squared coefficients of variation", {
  # With CVb = exp(s2) - 1 between strata and CVw = 1 / lambda within, the exact
  # VPC is CVb / (CVb + CVw). The lognormal latent-scale VPC is the same partition
  # after log(1 + .) of each part, and the delta one puts s2 in place of CVb. The
  # latent-scale values are those of fit_maihda(count_approximation = ), which
  # this feature leaves as they are.
  lambda <- c(0.34, 2, 10)
  s2 <- 0.2
  cvb <- exp(s2) - 1
  cvw <- 1 / lambda
  expect_equal(maihda_count_vpc_exact(lambda, s2)$vpc, cvb / (cvb + cvw))

  latent <- function(approximation) {
    s2 / (s2 + maihda_count_level1_transform(cvw, approximation))
  }
  expect_equal(latent("lognormal"), log1p(cvb) / (log1p(cvb) + log1p(cvw)))
  expect_equal(latent("delta"), s2 / (s2 + cvw))
  expect_equal(round(latent("lognormal"), 4), c(0.1273, 0.3303, 0.6773))
  expect_equal(round(latent("delta"), 4), c(0.0637, 0.2857, 0.6667))
  expect_equal(round(latent("trigamma"), 4), c(0.0201, 0.2367, 0.6554))
})

test_that("the closed form has the right limits and stays finite", {
  # No between-stratum variance, no share; no expected count, no share.
  expect_identical(maihda_count_vpc_exact(c(0.5, 3), 0)$vpc, c(0, 0))
  expect_identical(maihda_count_vpc_exact(0, 0.4)$vpc, 0)
  # alpha = 0 is the Poisson model.
  expect_identical(maihda_count_vpc_exact(3, 0.4, alpha = 0),
                   maihda_count_vpc_exact(3, 0.4))
  # Overdispersion is level-1 variance, so it lowers the share; a larger expected
  # count or between-stratum variance raises it.
  expect_lt(maihda_count_vpc_exact(3, 0.4, alpha = 0.5)$vpc,
            maihda_count_vpc_exact(3, 0.4)$vpc)
  expect_true(all(diff(maihda_count_vpc_exact(c(0.1, 1, 5, 50), 0.4)$vpc) > 0))
  expect_true(all(diff(vapply(c(0.05, 0.2, 1, 3), function(v) {
    maihda_count_vpc_exact(2, v, alpha = 0.3)$vpc
  }, numeric(1))) > 0))

  # The ratio is formed without lambda^2 or exp(s2), so a huge count or variance
  # reaches the limit -- 1 for a Poisson, (1 - exp(-s2)) / (1 - exp(-s2) + alpha)
  # for NB2 -- rather than Inf / Inf.
  expect_equal(maihda_count_vpc_exact(1e200, 0.3)$vpc, 1)
  expect_equal(maihda_count_vpc_exact(5, 800)$vpc, 1)
  expect_equal(maihda_count_vpc_exact(5, 800, alpha = 0.5)$vpc, 1 / 1.5)
  expect_equal(maihda_count_vpc_exact(1e200, 0.3, alpha = 0.5)$vpc,
               (1 - exp(-0.3)) / (1 - exp(-0.3) + 0.5))

  # Vectorised over lambda; one variance and one overdispersion.
  expect_length(maihda_count_vpc_exact(c(1, 2, 3), 0.2)$vpc, 3L)
  expect_error(maihda_count_vpc_exact(c(1, NA), 0.2), "lambda")
  expect_error(maihda_count_vpc_exact(-1, 0.2), "lambda")
  expect_error(maihda_count_vpc_exact(Inf, 0.2), "lambda")
  expect_error(maihda_count_vpc_exact(1, -0.1), "var_between")
  expect_error(maihda_count_vpc_exact(1, c(0.1, 0.2)), "var_between")
  expect_error(maihda_count_vpc_exact(1, 0.2, alpha = -1), "alpha")
  expect_error(maihda_count_vpc_exact(1, 0.2, alpha = Inf), "alpha")
})

test_that("the closed form agrees with a brute-force simulation of the model", {
  skip_on_cran()
  # Draw the stratum effect, form the conditional mean, and take the two terms of
  # the law of total variance as they stand: Var(E[y | u]) and E[Var(y | u)].
  brute <- function(lambda, s2, alpha, n = 2e6) {
    u <- stats::rnorm(n, 0, sqrt(s2))
    mu <- exp(log(lambda) - s2 / 2 + u)
    stats::var(mu) / (stats::var(mu) + mean(mu + alpha * mu^2))
  }
  cases <- list(c(0.34, 0.2, 0), c(2, 0.2, 0), c(10, 0.2, 0),
                c(2, 0.2, 0.5), c(10, 1, 0.5))
  simulated <- count_vpc_with_seed(20261011, vapply(cases, function(cs) {
    brute(cs[1], cs[2], cs[3])
  }, numeric(1)))
  exact <- vapply(cases, function(cs) {
    maihda_count_vpc_exact(cs[1], cs[2], alpha = cs[3])$vpc
  }, numeric(1))
  # Each case on its own. The nearest wrong formulas -- s2 for exp(s2) - 1, or the
  # NB2 term without its exp(s2) -- each move at least one case by 0.015 or more.
  expect_true(all(abs(simulated - exact) < 3e-3))
})

test_that("weighted quantiles are the quantiles of the duplicated-row data", {
  v <- c(0.31, 0.12, 0.47, 0.12, 0.26, 0.40, 0.05)
  w <- c(2, 1, 3, 4, 1, 2, 5)
  probs <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  expect_equal(maihda_weighted_obs_quantile(v, w, probs),
               stats::quantile(rep(v, w), probs, names = FALSE))

  # Absent, all-1 or length-mismatched weights give the plain quantiles, on the
  # rule maihda_weighted_obs_mean() follows.
  plain <- stats::quantile(v, probs, names = FALSE)
  expect_identical(maihda_weighted_obs_quantile(v, NULL, probs), plain)
  expect_identical(maihda_weighted_obs_quantile(v, rep(1, 7), probs), plain)
  expect_identical(maihda_weighted_obs_quantile(v, c(1, 2), probs), plain)

  # A non-finite or non-positive weight drops its observation, as for the mean.
  expect_equal(maihda_weighted_obs_quantile(c(1, 2, 3, 4), c(2, NA, 0, 1), 0.5),
               stats::quantile(c(1, 1, 4), 0.5, names = FALSE))
  # One observation is every quantile.
  expect_equal(maihda_weighted_obs_quantile(0.3, 4, probs), rep(0.3, 7))
})

# ---- through the public function ---------------------------------------------

test_that("a null Poisson fit returns the single closed-form value", {
  skip_on_cran()
  m <- count_vpc_fit("null")
  r <- maihda_vpc_response(m)

  expect_s3_class(r, "maihda_vpc_response")
  expect_identical(r$scale, "response")
  expect_identical(r$method, "closed-form")
  expect_identical(r$family, "poisson")

  # Rebuilt from lme4's own estimates: lambda = exp(b0 + s2 / 2).
  s2 <- as.numeric(lme4::VarCorr(m$model)$stratum[1, 1])
  lambda <- exp(lme4::fixef(m$model)[["(Intercept)"]] + s2 / 2)
  expect_equal(r$estimate, count_vpc_by_hand(lambda, s2))
  expect_equal(r$var_between, s2)
  expect_identical(r$var_other, 0)
  expect_identical(r$theta, Inf)
  expect_equal(r$lambda, lambda)
  # One lambda, so one VPC: every summary of it is that value.
  expect_equal(r$lambda_range, c(lambda, lambda))
  expect_equal(r$median, r$estimate)
  expect_equal(r$quartiles, rep(r$estimate, 2))
  expect_equal(r$range, rep(r$estimate, 2))

  # It is a different quantity from the latent-scale VPC that summary() reports
  # for the same fit: s2 over s2 plus the lognormal level-1 variance.
  latent <- count_vpc_quiet(summary(m))$vpc$estimate
  expect_equal(latent, s2 / (s2 + log1p(1 / lambda)))
  expect_gt(abs(r$estimate - latent), 0.01)

  expect_output(print(r), "Count-scale VPC")
  expect_output(print(r), "Marginal expected count")
})

test_that("the count path draws no random number and reads neither n_sim nor seed", {
  skip_on_cran()
  m <- count_vpc_fit("null")
  state <- function() get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  before <- state()
  r <- maihda_vpc_response(m)
  expect_identical(state(), before)
  # Both belong to the binomial simulation: neither changes the value, and an
  # n_sim the simulation would refuse is not an error here.
  expect_identical(maihda_vpc_response(m, n_sim = 10, seed = 99), r)
  expect_identical(state(), before)
})

test_that("covariates and an offset make the VPC observation-specific", {
  skip_on_cran()
  m <- count_vpc_fit("adjusted")
  d <- count_vpc_data()
  r <- maihda_vpc_response(m)

  # Each observation's marginal expected count: its fixed part, its own offset and
  # half the stratum variance.
  s2 <- as.numeric(lme4::VarCorr(m$model)$stratum[1, 1])
  b <- lme4::fixef(m$model)
  lambda <- d$expo * exp(b[["(Intercept)"]] + b[["age"]] * d$age + s2 / 2)
  vpc_i <- count_vpc_by_hand(lambda, s2)

  expect_equal(r$estimate, mean(vpc_i))
  expect_equal(r$median, stats::median(vpc_i))
  expect_equal(r$quartiles, stats::quantile(vpc_i, c(0.25, 0.75), names = FALSE))
  expect_equal(r$range, range(vpc_i))
  expect_equal(r$lambda, mean(lambda))
  expect_equal(r$lambda_range, range(lambda))
  # The VPC rises with the expected count, so it spreads with it, and its mean is
  # not the VPC at the mean count.
  expect_gt(diff(r$range), 0.1)
  expect_gt(abs(r$estimate - count_vpc_by_hand(mean(lambda), s2)), 1e-3)

  # The offset passed as an argument is the same model and the same VPC.
  m_arg <- count_vpc_quiet(
    fit_maihda(yx ~ age + (1 | gender:race:edu), data = d, family = "poisson",
               offset = log(expo)))
  r_arg <- maihda_vpc_response(m_arg)
  expect_equal(r_arg$estimate, r$estimate, tolerance = 1e-5)
  expect_equal(r_arg$lambda_range, r$lambda_range, tolerance = 1e-5)

  expect_output(print(r), "observation-level VPCs")
})

test_that("a negative-binomial fit uses the NB2 expressions with alpha = 1 / theta", {
  skip_on_cran()
  m <- count_vpc_fit("negbin")
  d <- count_vpc_data()
  r <- maihda_vpc_response(m)

  # glmer.nb() refits with MASS::negative.binomial(theta), whose variance function
  # is mu + mu^2 / theta: the NB2 model with alpha = 1 / theta.
  theta <- lme4::getME(m$model, "glmer.nb.theta")
  expect_equal(stats::family(m$model)$variance(3), 3 + 3^2 / theta)
  expect_identical(r$family, "negbinomial")
  expect_equal(r$theta, theta)

  s2 <- as.numeric(lme4::VarCorr(m$model)$stratum[1, 1])
  b <- lme4::fixef(m$model)
  lambda <- exp(b[["(Intercept)"]] + b[["age"]] * d$age + s2 / 2)
  expect_equal(r$estimate, mean(count_vpc_by_hand(lambda, s2, alpha = 1 / theta)))
  expect_equal(r$var_between, s2)
  # The overdispersion is level-1 variance, so the share is below the Poisson one
  # at the same counts and stratum variance.
  expect_lt(r$estimate, mean(count_vpc_by_hand(lambda, s2)))

  expect_output(print(r), "Negative-binomial theta")
})

test_that("a fixed negative-binomial theta is used as supplied", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  d <- count_vpc_data()
  m <- count_vpc_quiet(
    fit_maihda(ynb ~ age + (1 | gender:race:edu), data = d,
               family = MASS::negative.binomial(2)))
  r <- maihda_vpc_response(m)

  expect_identical(r$family, "negbinomial")
  expect_equal(r$theta, 2)
  s2 <- as.numeric(lme4::VarCorr(m$model)$stratum[1, 1])
  b <- lme4::fixef(m$model)
  lambda <- exp(b[["(Intercept)"]] + b[["age"]] * d$age + s2 / 2)
  expect_equal(r$estimate, mean(count_vpc_by_hand(lambda, s2, alpha = 0.5)))
})

test_that("prior weights enter the summaries as frequencies", {
  skip_on_cran()
  # A weighted fit and the same rows repeated by their weights are one model, so
  # they return one count-scale VPC. The weights follow the exposure, hence the
  # expected count, so that ignoring them would move every summary.
  d <- count_vpc_data()[seq_len(300), ]
  d$w <- ifelse(d$expo > 2.25, 3L, 1L)
  expanded <- d[rep(seq_len(nrow(d)), d$w), ]
  rownames(expanded) <- NULL
  f <- yx ~ age + offset(log(expo)) + (1 | gender:race:edu)
  m_w <- count_vpc_quiet(fit_maihda(f, data = d, family = "poisson", weights = w))
  m_e <- count_vpc_quiet(fit_maihda(f, data = expanded, family = "poisson"))
  r_w <- maihda_vpc_response(m_w)
  r_e <- maihda_vpc_response(m_e)

  # Exact on the weighted fit's own rows: the summaries of the per-row VPCs, each
  # row counted as many times as its weight.
  s2 <- as.numeric(lme4::VarCorr(m_w$model)$stratum[1, 1])
  b <- lme4::fixef(m_w$model)
  lambda <- d$expo * exp(b[["(Intercept)"]] + b[["age"]] * d$age + s2 / 2)
  vpc_rep <- rep(count_vpc_by_hand(lambda, s2), d$w)
  expect_equal(r_w$estimate, mean(vpc_rep))
  expect_equal(r_w$median, stats::median(vpc_rep))
  expect_equal(r_w$quartiles, stats::quantile(vpc_rep, c(0.25, 0.75), names = FALSE))
  expect_equal(r_w$lambda, mean(rep(lambda, d$w)))
  # ... which the unweighted summary of those rows is not.
  expect_gt(abs(r_w$estimate - mean(count_vpc_by_hand(lambda, s2))), 0.01)

  # And equal to the duplicated-row fit, up to what two optimisations of the same
  # likelihood differ by.
  for (el in c("estimate", "median", "quartiles", "range", "lambda",
               "lambda_range", "var_between")) {
    expect_equal(r_w[[el]], r_e[[el]], tolerance = 1e-3, info = el)
  }
})

test_that("a binomial fit keeps the simulation and its fields", {
  skip_on_cran()
  d <- count_vpc_data()
  d$high <- as.integer(d$y > 3)
  m <- count_vpc_quiet(
    fit_maihda(high ~ (1 | gender:race:edu), data = d, family = "binomial"))
  r <- maihda_vpc_response(m, n_sim = 2000, seed = 1)

  expect_identical(r$method, "simulation")
  expect_identical(names(r), c("estimate", "scale", "method", "n_sim", "var_between",
                               "var_other", "lp_fixed", "inner_method", "inner_nodes"))
  expect_error(maihda_vpc_response(m, n_sim = 10), "n_sim")
  expect_output(print(r), "Response-scale VPC \\(simulation method\\)")
})

# ---- refusals ----------------------------------------------------------------

test_that("fits beyond the single stratum random intercept are refused", {
  skip_on_cran()
  d <- count_vpc_data()
  only <- "only random effect is the stratum random intercept"

  # A context is a second, crossed random intercept ...
  expect_error(maihda_vpc_response(count_vpc_fit("context")),
               paste0(only, ".*context random intercept\\(s\\) school"))

  # ... and so is any other random intercept written into the formula.
  d_st <- d
  d_st$stratum <- make_strata(d, vars = c("gender", "race", "edu"))$data$stratum
  m_site <- count_vpc_quiet(
    fit_maihda(y ~ (1 | stratum) + (1 | school), data = d_st, family = "poisson"))
  expect_error(maihda_vpc_response(m_site),
               paste0(only, ".*random effect\\(s\\) school"))

  # A crossed-dimensions fit has a random intercept for each dimension as well.
  a_cc <- count_vpc_quiet(
    maihda(y ~ (1 | gender:race:edu), data = d, family = "poisson",
           decomposition = "crossed-dimensions"))
  expect_error(maihda_vpc_response(a_cc$model), paste0(only, ".*crossed-dimensions"))

  # A random slope, with the message the other VPC paths give.
  m_slope <- count_vpc_quiet(
    fit_maihda(y ~ age + (age | stratum), data = d_st, family = "poisson"))
  expect_error(maihda_vpc_response(m_slope), "intercept-only random effects")

  # A longitudinal fit is refused as one, by its tag.
  m_lng <- count_vpc_fit("null")
  m_lng$longitudinal_info <- list(id = "id", time = "wave")
  expect_error(maihda_vpc_response(m_lng), paste0(only, ".*longitudinal"))
})

test_that("other links, families and engines are refused", {
  skip_on_cran()
  d <- count_vpc_data()

  # The expressions are those of the log link.
  m_sqrt <- count_vpc_quiet(
    fit_maihda(y ~ (1 | gender:race:edu), data = d,
               family = stats::poisson(link = "sqrt")))
  expect_error(maihda_vpc_response(m_sqrt), "defined for the log link.*sqrt")

  m_gauss <- count_vpc_quiet(fit_maihda(age ~ (1 | gender:race:edu), data = d))
  expect_error(maihda_vpc_response(m_gauss),
               "binomial.*Poisson and negative-binomial.*gaussian")

  m_brms <- count_vpc_fit("null")
  m_brms$engine <- "brms"
  expect_error(maihda_vpc_response(m_brms), "lme4 engine")

  expect_error(maihda_vpc_response(count_vpc_fit("null")$model), "maihda_model")
})

# ---- summary() ---------------------------------------------------------------

test_that("summary(response_vpc = TRUE) attaches the count-scale VPC, and only then", {
  skip_on_cran()
  m <- count_vpc_fit("adjusted")
  s0 <- count_vpc_quiet(summary(m))
  s1 <- count_vpc_quiet(summary(m, response_vpc = TRUE))

  # Off by default: neither computed nor printed.
  expect_null(s0$vpc_response)
  expect_false(any(grepl("Count-scale VPC", utils::capture.output(print(s0)),
                         fixed = TRUE)))

  expect_identical(s1$vpc_response, maihda_vpc_response(m))
  expect_output(print(s1), "Count-scale VPC")
  # Both VPCs of the fit refer to the same marginal expected count.
  expect_equal(s1$vpc_response$lambda, s1$count_vpc$lambda)
  # Nothing else in the summary moves, the latent-scale VPC least of all.
  rest <- c("vpc", "count_vpc", "variance_components", "decomposition", "context",
            "longitudinal", "discriminatory_accuracy", "stratum_estimates",
            "fixed_effects", "conf_level", "df_method")
  expect_identical(s1[rest], s0[rest])

  # There is nothing to simulate, so a seed changes nothing.
  expect_identical(
    count_vpc_quiet(summary(m, response_vpc = TRUE, seed = 3))$vpc_response,
    s1$vpc_response)
})

test_that("summary(response_vpc = TRUE) skips a count fit the closed form does not cover", {
  skip_on_cran()
  warned <- character()
  s <- withCallingHandlers(
    suppressMessages(summary(count_vpc_fit("context"), response_vpc = TRUE)),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_null(s$vpc_response)
  # Skipped like any family or engine it does not cover, not reported as a
  # calculation that failed.
  expect_false(any(grepl("could not be computed", warned, fixed = TRUE)))
})

test_that("maihda(response_vpc = TRUE) carries it on both model summaries", {
  skip_on_cran()
  d <- count_vpc_data()
  a <- count_vpc_quiet(
    maihda(y ~ age + (1 | gender:race:edu), data = d, family = "poisson",
           response_vpc = TRUE))
  expect_identical(a$summary$vpc_response, maihda_vpc_response(a$model))
  expect_identical(a$summary_adjusted$vpc_response,
                   maihda_vpc_response(a$model_adjusted))
  expect_identical(a$summary$vpc_response$method, "closed-form")
})
