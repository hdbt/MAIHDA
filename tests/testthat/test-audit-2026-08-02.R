# Regression tests for the 2026-08-02 audit findings.
#
#   1 [High] fit_maihda() accepted and fitted a binomial complementary-log-log
#            (cloglog) model, but summary()/VPC then stop()ed with "VPC residual
#            variance is not implemented for family 'binomial' with link 'cloglog'".
#            The latent-scale residual-variance helpers had logit (pi^2/3) and probit
#            (1) branches but no cloglog branch, even though the response-scale VPC is
#            deliberately link-agnostic (logit/probit/cloglog). The cloglog latent
#            residual follows a Gumbel (extreme-value type I) distribution, variance
#            pi^2/6; a branch was added to all three latent-variance helpers so the
#            fitted object is usable. Cumulative (ordinal) cloglog stays rejected at
#            fit time (maihda_ordinal_check_family), so only binomial-type families
#            reach the new branch.
#
#   2 [High, related to finding on make_strata/decompose] The numeric-stratum-dimension
#            linear-term warning maihda() emits -- a category code (education 1..5)
#            entering the adjusted model as a single linear slope rather than
#            categorical main effects -- did NOT fire from pcv_importance() or
#            stepwise_pcv(), which build the same models through the shared attribution
#            setup. maihda() already documents, warns, and (for >10-value numerics)
#            messages on auto-binning, so the only real gap was this missing warning at
#            those two entry points; it is now emitted there too. No estimand change
#            (models are bit-identical; only a warning was added), per the maintainer's
#            decision to keep numeric dimensions as-is rather than auto-factor them.

# Collect every warning an expression emits (expect_warning inspects only one, and
# these fits can also emit unrelated singular-fit warnings).
collect_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }
  )
  w
}

# ---- Finding 1: cloglog binomial latent VPC is usable -----------------------

# Smallest binomial MAIHDA fit on a chosen inverse link. The latent residual
# variance is a CONSTANT of the link, so convergence of the fit is irrelevant to the
# assertions -- suppress the incidental optimizer chatter.
cloglog_fit <- function(seed = 1, n = 1200, link = "cloglog") {
  set.seed(seed)
  d <- data.frame(gender = sample(c("F", "M"), n, replace = TRUE),
                  race   = sample(c("A", "B", "C"), n, replace = TRUE))
  sk <- interaction(d$gender, d$race, drop = TRUE)
  lp <- stats::rnorm(nlevels(sk), sd = 0.8)[sk] - 0.5
  p  <- switch(link,
               cloglog = 1 - exp(-exp(lp)),
               logit   = stats::plogis(lp),
               probit  = stats::pnorm(lp))
  d$y <- stats::rbinom(n, 1, p)
  d$stratum <- make_strata(d, vars = c("gender", "race"))$data$stratum
  suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = stats::binomial(link = link))))
}

test_that("binomial cloglog latent residual variance is the Gumbel pi^2/6", {
  skip_on_cran()
  fit <- cloglog_fit()
  # Regression: this call stop()ed with "not implemented ... 'cloglog'" before the fix.
  rv <- MAIHDA:::maihda_residual_variance_lme4(fit$model)
  expect_equal(rv, pi^2 / 6)
})

test_that("summary() of a binomial cloglog model is usable (was an error)", {
  skip_on_cran()
  fit <- cloglog_fit()
  # The exact symptom the finding reported: summary() used to stop().
  expect_error(summary(fit), NA)
  s <- summary(fit)
  vb <- as.numeric(lme4::VarCorr(fit$model)[["stratum"]][1])
  # The latent VPC is now formed on the pi^2/6 scale.
  expect_equal(s$vpc$estimate, vb / (vb + pi^2 / 6), tolerance = 1e-8)
  # The response-scale VPC was already link-agnostic and still resolves.
  vr <- maihda_vpc_response(fit, n_sim = 2000, seed = 1)
  expect_true(is.finite(vr$estimate) && vr$estimate > 0 && vr$estimate < 1)
})

test_that("adding the cloglog branch did not perturb logit / probit variances", {
  skip_on_cran()
  # The three link branches are mutually exclusive; assert the neighbours are intact.
  expect_equal(MAIHDA:::maihda_residual_variance_lme4(cloglog_fit(link = "logit")$model),
               pi^2 / 3)
  expect_equal(MAIHDA:::maihda_residual_variance_lme4(cloglog_fit(link = "probit")$model),
               1)
})

# ---- Finding 2: linear-dim warning fires from the PCV attribution entry points ---

# education 1..5 is a numeric category code with <= 10 unique values, so make_strata()
# keeps it numeric; used as a stratum DIMENSION (not a covariate) it enters the
# adjusted/subset models as a single linear slope. maihda() warns about exactly this;
# pcv_importance()/stepwise_pcv() were silent before the fix.
make_numeric_dim_strata <- function(seed = 1, n_per = 120) {
  set.seed(seed)
  grid <- expand.grid(gender = c("F", "M"), education = 1:5, stringsAsFactors = FALSE)
  d <- grid[rep(seq_len(nrow(grid)), each = n_per), ]
  d$y <- -0.5 * (d$gender == "M") +
    c(2, 0.2, -0.6, 0.1, 1.9)[d$education] + stats::rnorm(nrow(d))
  rownames(d) <- NULL
  make_strata(d, vars = c("gender", "education"))$data
}

test_that("pcv_importance() warns when a numeric stratum dimension enters linearly", {
  skip_on_cran()
  d <- make_numeric_dim_strata()
  ws <- collect_warnings(
    pcv_importance(d, "y", c("gender", "education"), family = "gaussian",
                   method = "shapley", approx = "exact"))
  expect_true(any(grepl("linear fixed effect", ws)))
})

test_that("stepwise_pcv() warns when a numeric stratum dimension enters linearly", {
  skip_on_cran()
  d <- make_numeric_dim_strata()
  ws <- collect_warnings(
    stepwise_pcv(d, "y", c("gender", "education"), family = "gaussian"))
  expect_true(any(grepl("linear fixed effect", ws)))
})

test_that("no linear-dim warning when the stratum dimension is a factor (no false positive)", {
  skip_on_cran()
  base <- make_numeric_dim_strata()
  base$education <- factor(base$education)
  d <- make_strata(base, vars = c("gender", "education"))$data  # factor is now the recorded type
  ws <- collect_warnings(
    stepwise_pcv(d, "y", c("gender", "education"), family = "gaussian"))
  expect_false(any(grepl("linear fixed effect", ws)))
})

test_that("a numeric COVARIATE (not a stratum dimension) does not trigger the warning", {
  skip_on_cran()
  # Confirms the fix is confined via intersect(strata_vars, vars): age is numeric and
  # in vars, but is not a stratum dimension, so it must NOT be flagged (this mirrors
  # the pre-existing expect_no_warning tests the fix must not break).
  set.seed(2)
  n <- 900
  d <- data.frame(gender = sample(c("m", "f"), n, TRUE),
                  race   = sample(c("a", "b", "c"), n, TRUE),
                  age    = stats::rnorm(n))
  d <- make_strata(d, c("gender", "race"))$data
  d$y <- 1 + ifelse(d$gender == "m", 0.8, 0) +
    c(a = 0, b = 0.9, c = -0.4)[d$race] + 0.3 * d$age + stats::rnorm(n)
  ws <- collect_warnings(
    pcv_importance(d, "y", c("gender", "race", "age"), method = "shapley",
                   approx = "exact"))
  expect_false(any(grepl("linear fixed effect", ws)))
})
