test_that("calculate_pcv works with basic models", {
  # Create test data with actual stratum effects
  set.seed(123)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 2)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit two models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum), 
                       data = data, 
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PCV
  pcv_result <- calculate_pcv(model1, model2)
  
  # Check structure
  expect_true(inherits(pcv_result, "pcv_result"))
  expect_true(is.numeric(pcv_result$pcv))
  expect_true(is.numeric(pcv_result$var_model1))
  expect_true(is.numeric(pcv_result$var_model2))
  expect_false(pcv_result$bootstrap)
  
  # Check that variances are positive
  expect_true(pcv_result$var_model1 > 0)
  expect_true(pcv_result$var_model2 > 0)
})

test_that("calculate_pcv works with bootstrap", {
  # Create test data with actual stratum effects
  set.seed(456)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1.5)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit two models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PCV with bootstrap (small number for testing; the low-n_boot
  # stability warning is expected and suppressed)
  pcv_result <- suppressWarnings(
    calculate_pcv(model1, model2, bootstrap = TRUE, n_boot = 50))
  
  # Check structure
  expect_true(inherits(pcv_result, "pcv_result"))
  expect_true(pcv_result$bootstrap)
  expect_true(!is.null(pcv_result$ci_lower))
  expect_true(!is.null(pcv_result$ci_upper))
  expect_true(is.numeric(pcv_result$ci_lower))
  expect_true(is.numeric(pcv_result$ci_upper))
  
  # A percentile bootstrap interval is well-ordered and finite, but it need NOT
  # contain the point estimate (the original-fit PCV): the bootstrap distribution can
  # be skewed -- especially at a small n_boot -- so containment is not a valid
  # invariant to assert.
  expect_true(is.finite(pcv_result$ci_lower) && is.finite(pcv_result$ci_upper))
  expect_true(pcv_result$ci_lower <= pcv_result$ci_upper)
})

test_that("calculate_pcv validates inputs", {
  # Create test data with actual stratum effects
  set.seed(789)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Invalid first argument
  expect_error(calculate_pcv("not a model", model1),
               "must be a maihda_model")
  
  # Invalid second argument
  expect_error(calculate_pcv(model1, "not a model"),
               "must be a maihda_model")
  
  # Both arguments invalid
  expect_error(calculate_pcv(data, data),
               "must be a maihda_model")
})

test_that("calculate_pcv rejects row-wise different stratum assignments", {
  set.seed(790)
  n <- 80
  d1 <- data.frame(
    stratum = factor(rep(seq_len(4), each = n / 4)),
    x = rnorm(n)
  )
  d1$y <- 2 + 0.4 * d1$x + rnorm(4, sd = 0.8)[d1$stratum] + rnorm(n, sd = 0.2)

  d2 <- d1
  d2$stratum <- factor(rep(seq_len(4), times = n / 4))

  model1 <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | stratum), data = d1)
  ))
  model2 <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | stratum), data = d2)
  ))

  expect_error(
    calculate_pcv(model1, model2),
    "assign each analytic row to the same stratum",
    fixed = TRUE
  )
})

test_that("calculate_pcv handles same model comparison", {
  # Create test data with actual stratum effects
  set.seed(111)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit same model twice
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PCV
  pcv_result <- calculate_pcv(model1, model2)
  
  # PCV should be very close to 0 (same model structure)
  expect_true(abs(pcv_result$pcv) < 0.1)
})

test_that("calculate_pcv calculates correct direction", {
  # Create test data with known structure
  set.seed(222)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 2)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Model without gender (should have more between-stratum variance)
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Model with gender (might reduce between-stratum variance if gender explains some)
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PCV
  pcv_result <- calculate_pcv(model1, model2)
  
  # Variances should be positive
  expect_true(pcv_result$var_model1 > 0)
  expect_true(pcv_result$var_model2 > 0)
  
  # PCV formula: (var1 - var2) / var1
  expected_pcv <- (pcv_result$var_model1 - pcv_result$var_model2) / pcv_result$var_model1
  expect_equal(pcv_result$pcv, expected_pcv)
})

test_that("calculate_pcv print method works", {
  # Create test data with actual stratum effects
  set.seed(333)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1.5)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Add stratum effects to outcome
  data$outcome <- 5 + 0.5 * data$age + stratum_effects[data$stratum] + rnorm(nrow(data), sd = 1)
  
  # Fit models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4")
  
  # Calculate PCV
  pcv_result <- calculate_pcv(model1, model2)
  
  # Print should work without error
  expect_output(print(pcv_result), "Proportional Change in Variance")
  expect_output(print(pcv_result), "PCV:")
  expect_output(print(pcv_result), "Between-stratum variance:")
})

test_that("calculate_pcv handles binomial models", {
  # Create test data for binomial with stratum effects on logit scale
  set.seed(444)
  n_strata <- 10
  n_per_stratum <- 10
  
  # Generate stratum-level random effects on logit scale
  stratum_effects <- rnorm(n_strata, mean = 0, sd = 1)
  
  data <- data.frame(
    stratum = rep(1:n_strata, each = n_per_stratum),
    age = rnorm(n_strata * n_per_stratum),
    gender = sample(c(0, 1), n_strata * n_per_stratum, replace = TRUE)
  )
  
  # Generate outcome on logit scale with stratum effects
  logit_p <- -0.5 + 0.3 * data$age + stratum_effects[data$stratum]
  prob <- plogis(logit_p)
  data$outcome <- rbinom(nrow(data), 1, prob)
  
  # Fit binomial models
  model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                       data = data,
                       engine = "lme4",
                       family = "binomial")
  model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum),
                       data = data,
                       engine = "lme4",
                       family = "binomial")
  
  # Calculate PCV
  pcv_result <- calculate_pcv(model1, model2)
  
  # Check structure
  expect_true(inherits(pcv_result, "pcv_result"))
  expect_true(is.numeric(pcv_result$pcv))
  expect_true(pcv_result$var_model1 > 0)
  expect_true(pcv_result$var_model2 > 0)
})

test_that("calculate_pcv handles zero variance error", {
  # Create test data with NO stratum effects (will result in singular fit)
  set.seed(555)
  data <- data.frame(
    stratum = rep(1:10, each = 10),
    age = rnorm(100),
    outcome = 5 + 0.5 * rnorm(100)  # No stratum effects
  )
  
  # Fit models (may have singular fit warnings)
  suppressWarnings({
    model1 <- fit_maihda(outcome ~ age + (1 | stratum),
                         data = data,
                         engine = "lme4")
    model2 <- fit_maihda(outcome ~ 1 + (1 | stratum),
                         data = data,
                         engine = "lme4")
  })
  
  # Should error due to zero/negative variance
  expect_error(calculate_pcv(model1, model2),
               "Between-stratum variance")
})

test_that("maihda_bootstrap_ci requires a minimum number of successful draws", {
  expect_error(
    MAIHDA:::maihda_bootstrap_ci(numeric(0), n_boot = 100, conf_level = 0.95, what = "VPC"),
    "All"
  )
  expect_error(
    MAIHDA:::maihda_bootstrap_ci(rep(0.5, 5), n_boot = 100, conf_level = 0.95, what = "VPC"),
    "at least"
  )
  ci <- MAIHDA:::maihda_bootstrap_ci(seq(0, 1, length.out = 50), n_boot = 100, conf_level = 0.95)
  expect_length(ci, 2)
  expect_true(ci[1] < ci[2])
})

test_that("calculate_pcv rejects models fitted to unrelated data of the same shape", {
  mk <- function(seed) {
    set.seed(seed)
    d <- data.frame(stratum = factor(rep(seq_len(8), each = 20)), x = rnorm(160))
    d$y <- 1 + 0.3 * d$x + rnorm(8, sd = 0.7)[d$stratum] + rnorm(160, sd = 0.4)
    fit_maihda(y ~ x + (1 | stratum), data = d)
  }
  # Same n, strata and default 1:160 row names, but different outcome values.
  m1 <- mk(101)
  m2 <- mk(202)
  expect_error(calculate_pcv(m1, m2), "outcome values differ")
})

test_that("calculate_pcv validates bootstrap arguments before model comparison", {
  fake_model <- structure(
    list(engine = "lme4"),
    class = "maihda_model"
  )

  expect_error(
    calculate_pcv(fake_model, fake_model, bootstrap = c(TRUE, FALSE)),
    "'bootstrap' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    calculate_pcv(fake_model, fake_model, bootstrap = TRUE, n_boot = 0),
    "'n_boot' must be a single whole number >= 10",
    fixed = TRUE
  )
  expect_error(
    calculate_pcv(fake_model, fake_model, bootstrap = TRUE, conf_level = 1),
    "'conf_level' must be a single number between 0 and 1",
    fixed = TRUE
  )
})

test_that("bootstrap_pcv rejects non-lme4 engines with honest per-engine guidance", {
  fake <- function(engine) {
    structure(list(engine = engine), class = "maihda_model")
  }
  # brms: the old message circularly told brms users to "refit with brms";
  # the honest guidance is that no PCV interval exists (point estimate only).
  expect_error(
    MAIHDA:::bootstrap_pcv(fake("brms"), fake("brms"), 10, 0.95),
    "point estimate"
  )
  brms_msg <- tryCatch(MAIHDA:::bootstrap_pcv(fake("brms"), fake("brms"), 10, 0.95),
                       error = conditionMessage)
  expect_false(grepl("refit with engine", brms_msg, fixed = TRUE))
  expect_error(
    MAIHDA:::bootstrap_pcv(fake("wemix"), fake("wemix"), 10, 0.95),
    "replicate weights"
  )
  expect_error(
    MAIHDA:::bootstrap_pcv(fake("ordinal"), fake("ordinal"), 10, 0.95),
    "point estimate"
  )
})

test_that("calculate_pvc is a deprecated alias of calculate_pcv", {
  set.seed(42)
  data <- data.frame(
    stratum = factor(rep(1:8, each = 20)),
    age = rnorm(160)
  )
  data$outcome <- 2 + 0.3 * data$age +
    rnorm(8, sd = 0.8)[data$stratum] + rnorm(160, sd = 0.5)

  model1 <- fit_maihda(outcome ~ 1 + (1 | stratum), data = data)
  model2 <- fit_maihda(outcome ~ age + (1 | stratum), data = data)

  new <- calculate_pcv(model1, model2)
  # Match on the condition class, not the message: .Deprecated()'s text is localized.
  expect_warning(old <- calculate_pvc(model1, model2), class = "deprecatedWarning")
  expect_identical(old, new)

  # The result answers to both the new and the historical class/element names.
  expect_s3_class(new, "pcv_result")
  expect_s3_class(new, "pvc_result")
  expect_identical(new$pvc, new$pcv)
})

test_that("PCV bootstrap reports zero-null-variance boundary draws", {
  # Regression: draws whose null model estimated a zero between-stratum
  # variance were silently converted to NA and dropped -- the advertised
  # percentile interval was conditional on a positive null variance with no
  # warning unless over half of ALL draws failed. Any boundary mass now warns
  # and is reported on the result.
  skip_on_cran()
  set.seed(8)
  n <- 320
  d <- data.frame(
    g = sample(c("F", "M"), n, replace = TRUE),
    r = sample(c("A", "B"), n, replace = TRUE),
    x = rnorm(n)
  )
  sk <- interaction(d$g, d$r, drop = TRUE)
  d$y <- 0.3 * d$x + rnorm(nlevels(sk), sd = 0.16)[sk] + rnorm(n)
  st <- make_strata(d, c("g", "r"))
  d$stratum <- st$data$stratum
  m1 <- suppressWarnings(suppressMessages(fit_maihda(y ~ (1 | stratum), data = d)))
  m2 <- suppressWarnings(suppressMessages(fit_maihda(y ~ x + (1 | stratum), data = d)))

  # capture_warnings: the simulated refits can also emit optimizer-convergence
  # warnings, which are irrelevant here.
  set.seed(1008)
  w <- capture_warnings(
    r <- suppressMessages(calculate_pcv(m1, m2, bootstrap = TRUE, n_boot = 40)))
  expect_true(any(grepl("conditional on a positive null variance", w)))
  expect_true(is.finite(r$n_boot_boundary) && r$n_boot_boundary > 0)
  expect_equal(r$n_boot_ok + r$n_boot_boundary, 40)
  expect_output(print(r), "boundary")
  # The surviving draws form a sane conditional interval (near-boundary
  # denominators that produced intervals like [-6.6e16, 1] are excluded).
  expect_true(r$ci_lower > -100 && r$ci_upper <= 1 + 1e-8)
})

test_that("calculate_pcv rejects an effectively singular null model", {
  # A null variance like 1e-9 passed the strict var1 <= 0 guard and produced
  # PCVs in the thousands (audit follow-up: PCV = -120427 with interval
  # [-6.6e16, 1]); the lme4 relative singularity tolerance now rejects the
  # degenerate denominator with a clear message.
  skip_on_cran()
  set.seed(4)
  n <- 320
  d <- data.frame(
    g = sample(c("F", "M"), n, replace = TRUE),
    r = sample(c("A", "B"), n, replace = TRUE),
    x = rnorm(n)
  )
  sk <- interaction(d$g, d$r, drop = TRUE)
  d$y <- 0.3 * d$x + rnorm(nlevels(sk), sd = 0.12)[sk] + rnorm(n)
  st <- make_strata(d, c("g", "r"))
  d$stratum <- st$data$stratum
  m1 <- suppressWarnings(suppressMessages(fit_maihda(y ~ (1 | stratum), data = d)))
  m2 <- suppressWarnings(suppressMessages(fit_maihda(y ~ x + (1 | stratum), data = d)))

  expect_error(suppressWarnings(suppressMessages(calculate_pcv(m1, m2))),
               "zero or negative|zero boundary")
})

test_that("calculate_pcv estimation basis: 'fitted' (default) uses REML, 'ML' refits", {
  # Audit finding 1: for a Gaussian lmer PCV the forced ML refit materially changes
  # the reported PCV with few strata. The estimation argument exposes both bases;
  # the default is now "fitted" (each fit's own REML variance), matching summary().
  skip_on_cran()
  set.seed(20260711)
  n_strata <- 24; n_per <- 12
  smean <- rnorm(n_strata, sd = 2)                    # stratum means
  x1_stratum <- smean + rnorm(n_strata, sd = 1)       # stratum-level covariate
  d <- data.frame(
    stratum = rep(seq_len(n_strata), each = n_per),
    age = rnorm(n_strata * n_per),
    x1 = rep(x1_stratum, each = n_per))
  d$outcome <- 5 + 0.5 * d$age + smean[d$stratum] + rnorm(nrow(d), sd = 1)
  # x1 is a stratum-level covariate that explains part of the between-stratum
  # variance, so the adjusted model has a clearly lower stratum variance.
  m1 <- fit_maihda(outcome ~ age + (1 | stratum), data = d, engine = "lme4")
  m2 <- fit_maihda(outcome ~ age + x1 + (1 | stratum), data = d, engine = "lme4")
  expect_true(lme4::isREML(m1$model))

  # Default = "fitted": PCV variances equal the fitted REML VarCorr variances.
  p_fit <- calculate_pcv(m1, m2)
  expect_identical(p_fit$estimation, "fitted")
  v1_reml <- as.numeric(lme4::VarCorr(m1$model)$stratum[1])
  v2_reml <- as.numeric(lme4::VarCorr(m2$model)$stratum[1])
  expect_equal(p_fit$var_model1, v1_reml)
  expect_equal(p_fit$var_model2, v2_reml)
  expect_equal(p_fit$pcv, (v1_reml - v2_reml) / v1_reml)

  # "ML": PCV variances equal the ML-refit VarCorr variances.
  p_ml <- calculate_pcv(m1, m2, estimation = "ML")
  expect_identical(p_ml$estimation, "ML")
  v1_ml <- as.numeric(lme4::VarCorr(lme4::refitML(m1$model))$stratum[1])
  v2_ml <- as.numeric(lme4::VarCorr(lme4::refitML(m2$model))$stratum[1])
  expect_equal(p_ml$var_model1, v1_ml)
  expect_equal(p_ml$var_model2, v2_ml)

  # The REML and ML adjusted variances differ (the finding), so the two bases give
  # a different PCV. (Direction is not asserted -- it depends on the design.)
  expect_false(isTRUE(all.equal(v2_reml, v2_ml)))
  expect_false(isTRUE(all.equal(p_fit$pcv, p_ml$pcv)))

  # print() states the variance basis; a bad value is rejected (message is locale-
  # dependent, so only the error itself is asserted).
  expect_output(print(p_fit), "as fitted")
  expect_output(print(p_ml), "ML-refit")
  expect_error(calculate_pcv(m1, m2, estimation = "reml"))
})

test_that("calculate_pcv estimation is a no-op for non-REML (glmer) fits", {
  # For a binomial glmer the fit is already ML, so "fitted" and "ML" coincide.
  skip_on_cran()
  set.seed(21)
  n_strata <- 20; n_per <- 40
  eff <- rnorm(n_strata, sd = 1)
  d <- data.frame(stratum = rep(seq_len(n_strata), each = n_per),
                  age = rnorm(n_strata * n_per),
                  z = rnorm(n_strata * n_per))
  d$y <- rbinom(nrow(d), 1, plogis(-0.2 + 0.5 * d$age + eff[d$stratum]))
  m1 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | stratum), data = d, family = "binomial")))
  m2 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + z + (1 | stratum), data = d, family = "binomial")))
  p_fit <- calculate_pcv(m1, m2)
  p_ml  <- calculate_pcv(m1, m2, estimation = "ML")
  expect_equal(p_fit$pcv, p_ml$pcv)
  expect_equal(p_fit$var_model1, p_ml$var_model1)
})

test_that("stepwise_pcv and pcv_importance accept the estimation argument", {
  # The estimation basis threads through the model-series PCV attributions too.
  skip_on_cran()
  set.seed(99)
  N <- 18 * 22
  d <- data.frame(
    g = factor(sample(letters[1:3], N, replace = TRUE)),   # 3 x 6 = up to 18 strata
    r = factor(sample(LETTERS[1:6], N, replace = TRUE)),
    age = rnorm(N))
  st <- make_strata(d, c("g", "r"))
  d <- st$data
  slev <- unique(as.character(d$stratum))
  eff <- stats::setNames(rnorm(length(slev), sd = 2), slev)  # per-stratum (incl. interaction)
  d$y <- 3 + 0.5 * d$age + eff[as.character(d$stratum)] + rnorm(nrow(d))

  sw_fit <- stepwise_pcv(d, "y", c("g", "r", "age"))                 # default fitted
  sw_ml  <- stepwise_pcv(d, "y", c("g", "r", "age"), estimation = "ML")
  expect_s3_class(sw_fit, "maihda_stepwise")
  # The two bases give different Total_PCV trajectories for the Gaussian lmer fits.
  expect_false(isTRUE(all.equal(sw_fit$Total_PCV, sw_ml$Total_PCV)))

  imp_fit <- pcv_importance(d, "y", c("g", "r", "age"))
  imp_ml  <- pcv_importance(d, "y", c("g", "r", "age"), estimation = "ML")
  expect_false(isTRUE(all.equal(imp_fit$total_pcv, imp_ml$total_pcv)))
})

test_that("stepwise_pcv rejects an effectively singular null model", {
  # Audit follow-up: calculate_pcv() rejects a boundary-level null denominator,
  # but stepwise_pcv() divided Total_PCV / Step_PCV by it, turning a degenerate
  # ~0 null variance into a spurious 100% (or a huge negative) PCV with no
  # warning. The same lme4 relative singularity guard now rejects it up front.
  # Reuses the boundary construction from the calculate_pcv boundary test above.
  skip_on_cran()
  set.seed(4)
  n <- 320
  d <- data.frame(
    g = sample(c("F", "M"), n, replace = TRUE),
    r = sample(c("A", "B"), n, replace = TRUE),
    x = rnorm(n)
  )
  sk <- interaction(d$g, d$r, drop = TRUE)
  d$y <- 0.3 * d$x + rnorm(nlevels(sk), sd = 0.12)[sk] + rnorm(n)
  st <- make_strata(d, c("g", "r"))
  d$stratum <- st$data$stratum

  expect_error(
    suppressWarnings(suppressMessages(stepwise_pcv(d, "y", "x"))),
    "zero boundary"
  )
})
