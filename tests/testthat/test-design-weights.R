# Design-weighted MAIHDA (sampling/survey weights): the wemix pseudo-ML engine,
# the brms pseudo-posterior weight plumbing (Stan-free), and the guard rails that
# keep sampling weights away from lme4's precision weights.

make_dw_data <- function(seed = 8101, n = 1200) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("A", "B", "C"), n, replace = TRUE),
    edu = sample(c("low", "high"), n, replace = TRUE),
    age = rnorm(n, 45, 10),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$gender, d$race, d$edu, drop = TRUE)
  u <- rnorm(nlevels(stratum), sd = 0.5)[stratum]
  d$y <- 2 + 0.3 * (d$gender == "M") + 0.5 * (d$race == "B") + 0.02 * d$age +
    u + rnorm(n, sd = 1.1)
  d$ybin <- rbinom(n, 1, stats::plogis(-0.5 + 0.4 * (d$gender == "M") + u))
  d$w <- runif(n, 0.5, 4)
  d
}

# ---- sampling-weight validation (no WeMix required) -------------------------

test_that("maihda_validate_sampling_weights validates the specification", {
  d <- data.frame(y = 1:5, w = c(1, 2, 3, 4, 5), chr = letters[1:5])

  expect_identical(maihda_validate_sampling_weights("w", d), "w")

  expect_error(maihda_validate_sampling_weights(c("a", "b"), d), "single column name")
  expect_error(maihda_validate_sampling_weights(1, d), "single column name")
  expect_error(maihda_validate_sampling_weights("", d), "single column name")
  expect_error(maihda_validate_sampling_weights("nope", d), "not found")
  expect_error(maihda_validate_sampling_weights("chr", d), "must be numeric")
  expect_error(maihda_validate_sampling_weights(".maihda_l2wt", d), "reserved")
  expect_error(maihda_validate_sampling_weights(".maihda_sw", d), "reserved")

  d$bad <- c(0, -1, NA, NaN, -Inf)
  expect_error(maihda_validate_sampling_weights("bad", d),
               "no positive finite values")
})

test_that("fit_maihda rejects sampling weights with the lme4 engine", {
  d <- make_dw_data()
  expect_error(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "lme4", sampling_weights = "w"),
    "precision weights"
  )
})

test_that("fit_maihda rejects sampling_weights together with precision weights", {
  d <- make_dw_data()
  d$pw <- 1
  expect_error(
    suppressMessages(fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
                                sampling_weights = "w", weights = pw)),
    "not both"
  )
})

test_that("engine = wemix requires sampling weights and the canonical structure", {
  d <- make_dw_data()

  expect_error(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d, engine = "wemix"),
    "requires 'sampling_weights'"
  )
  d$cnt <- rpois(nrow(d), 2)
  expect_error(
    fit_maihda(cnt ~ age + (1 | gender:race:edu), data = d, engine = "wemix",
               family = "poisson", sampling_weights = "w"),
    "gaussian\\(identity\\) and binomial\\(logit\\)"
  )
  # Extra random effects: not the canonical single (1 | stratum) structure.
  s <- make_strata(d, vars = c("gender", "race", "edu"))
  d2 <- s$data
  d2$site <- sample(c("s1", "s2"), nrow(d2), replace = TRUE)
  expect_error(
    fit_maihda(y ~ age + (1 | stratum) + (1 | site), data = d2,
               engine = "wemix", sampling_weights = "w"),
    "single intercept-only random effect"
  )
  # Contextual cross-classification needs crossed REs, which WeMix cannot fit.
  expect_error(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d, engine = "wemix",
               sampling_weights = "w", context = "edu"),
    "does not support 'context'"
  )
  # lme4-style data-masked engine arguments have no WeMix counterpart.
  expect_error(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d, engine = "wemix",
               sampling_weights = "w", subset = age > 40),
    "not supported by engine"
  )
})

test_that("maihda() rejects wemix-incompatible workflow options early", {
  d <- make_dw_data()
  expect_error(
    suppressMessages(maihda(y ~ age + (1 | gender:race:edu), data = d,
                            sampling_weights = "w",
                            decomposition = "crossed-dimensions")),
    "crossed random effects"
  )
  expect_error(
    suppressMessages(maihda(y ~ age + (1 | gender:race:edu), data = d,
                            sampling_weights = "w", bootstrap = TRUE)),
    "replicate weights"
  )
  expect_error(
    maihda(y ~ age + (1 | gender:race:edu), data = d,
           engine = "lme4", sampling_weights = "w"),
    "precision weights"
  )
})

test_that("stepwise_pcv and compare_maihda_groups reject lme4 + sampling weights", {
  d <- make_dw_data()
  s <- make_strata(d, vars = c("gender", "race", "edu"))
  expect_error(
    stepwise_pcv(s$data, "y", c("gender", "race"), engine = "lme4",
                 sampling_weights = "w"),
    "precision weights"
  )
  d$country <- sample(c("X", "Y"), nrow(d), replace = TRUE)
  expect_error(
    compare_maihda_groups(y ~ age + (1 | gender:race), d, group = "country",
                          engine = "lme4", sampling_weights = "w"),
    "precision weights"
  )
  expect_error(
    suppressMessages(
      compare_maihda_groups(y ~ age + (1 | gender:race), d, group = "country",
                            sampling_weights = "w",
                            decomposition = "crossed-dimensions")),
    "crossed random effects"
  )
})

# ---- brms weight plumbing (Stan-free) ---------------------------------------

test_that("maihda_brms_weights_formula rewrites the LHS addition term", {
  f1 <- maihda_brms_weights_formula(y ~ x + (1 | stratum), ".maihda_sw")
  expect_identical(
    paste(deparse(f1), collapse = " "),
    "y | weights(.maihda_sw) ~ x + (1 | stratum)"
  )

  # An existing addition term is extended, not replaced.
  f2 <- maihda_brms_weights_formula(y | trials(n) ~ x + (1 | stratum), ".maihda_sw")
  expect_identical(
    paste(deparse(f2), collapse = " "),
    "y | trials(n) + weights(.maihda_sw) ~ x + (1 | stratum)"
  )

  # A formula that already carries weights() conflicts.
  expect_error(
    maihda_brms_weights_formula(y | weights(w0) ~ x + (1 | stratum), ".maihda_sw"),
    "already carries a weights"
  )
})

test_that("maihda_prepare_brms_sampling_weights normalizes and drops bad rows", {
  d <- data.frame(y = rnorm(10), stratum = rep(c("a", "b"), 5),
                  w = c(2, 4, 6, 8, 10, 2, 4, NA, 0, -1))

  expect_warning(
    prep <- maihda_prepare_brms_sampling_weights(d, y ~ (1 | stratum), "w"),
    "dropped 3 row"
  )
  expect_equal(nrow(prep$data), 7)
  # Normalized to mean 1 so expansion weights do not inflate the effective n.
  expect_equal(mean(prep$data$.maihda_sw), 1)
  # Relative weights preserved.
  expect_equal(prep$data$.maihda_sw[2] / prep$data$.maihda_sw[1], 2)
  expect_match(paste(deparse(prep$formula), collapse = " "),
               "weights(.maihda_sw)", fixed = TRUE)

  d_bad <- data.frame(y = rnorm(3), stratum = c("a", "b", "a"), w = c(NA, 0, -2))
  expect_error(
    maihda_prepare_brms_sampling_weights(d_bad, y ~ (1 | stratum), "w"),
    "No usable rows"
  )
})

# ---- sampling-weight fingerprint ---------------------------------------------

test_that("maihda_sampling_weight_fingerprint distinguishes weight specifications", {
  base <- list(sampling_weights = NULL, data = data.frame(w = c(1, 2)))
  w1 <- list(sampling_weights = "w", data = data.frame(w = c(1, 2)))
  w2 <- list(sampling_weights = "w", data = data.frame(w = c(2, 1)))
  # brms analytic frames carry the normalized column instead of the original.
  wb <- list(sampling_weights = "w",
             data = stats::setNames(data.frame(c(2 / 3, 4 / 3)), ".maihda_sw"))

  expect_identical(maihda_sampling_weight_fingerprint(base), "none")
  expect_identical(maihda_sampling_weight_fingerprint(w1),
                   maihda_sampling_weight_fingerprint(w1))
  expect_false(identical(maihda_sampling_weight_fingerprint(w1),
                         maihda_sampling_weight_fingerprint(w2)))
  expect_false(identical(maihda_sampling_weight_fingerprint(w1),
                         maihda_sampling_weight_fingerprint(base)))
  expect_match(maihda_sampling_weight_fingerprint(wb), "^w:")
})

# ---- design-weighted AUC -----------------------------------------------------

test_that("the design-weighted AUC equals the expanded-data AUC", {
  # Integer weights so the weighted AUC has an exact expanded-data counterpart.
  prob <- c(0.1, 0.4, 0.35, 0.8, 0.5)
  y <- c(0, 0, 1, 1, 0)
  w <- c(1, 2, 3, 1, 2)
  weighted <- maihda_auc_weighted(prob, successes = w * y, trials = w)
  expanded <- maihda_auc(rep(prob, w), rep(y, w))
  expect_equal(weighted, expanded)
})

# ---- wemix engine (WeMix required) -------------------------------------------

test_that("fit_maihda auto-switches to wemix and fits a weighted gaussian model", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  expect_message(
    m <- fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
                    sampling_weights = "w"),
    "using engine \"wemix\"|wemix"
  )
  expect_s3_class(m, "maihda_model")
  expect_identical(m$engine, "wemix")
  expect_identical(m$sampling_weights, "w")
  expect_true(inherits(m$model, "WeMixResults"))
  expect_true("stratum" %in% names(m$data))
  expect_identical(m$diagnostics$engine, "wemix")
  expect_true(isTRUE(m$diagnostics$converged))
  expect_false(isTRUE(m$diagnostics$singular))

  # print() mentions the design-weighted fit.
  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(out, "Sampling weights: w")
  expect_match(out, "pseudo-maximum-likelihood")

  s <- summary(m)
  expect_s3_class(s, "maihda_summary")
  expect_gt(s$vpc$estimate, 0)
  expect_lt(s$vpc$estimate, 1)
  # Design-consistent (sandwich) fixed-effect standard errors are reported.
  expect_true("se" %in% names(s$fixed_effects))
  expect_true(all(is.finite(s$fixed_effects$se)))
  expect_true(all(s$fixed_effects$se > 0))
  # One stratum estimate per populated stratum, with finite conditional SEs.
  expect_equal(nrow(s$stratum_estimates),
               length(unique(as.character(m$data$stratum))))
  expect_true(all(is.finite(s$stratum_estimates$se)))

  expect_output(print(s), "Variance Partition Coefficient")
})

test_that("a unit-weight wemix fit reproduces the lme4 ML fit", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  d$wu <- 1
  mw <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "wu"))
  ml <- fit_maihda(y ~ age + (1 | gender:race:edu), data = d, REML = FALSE)

  expect_equal(summary(mw)$vpc$estimate, summary(ml)$vpc$estimate,
               tolerance = 1e-6)
  expect_equal(
    unname(mw$model$coef[c("(Intercept)", "age")]),
    unname(lme4::fixef(ml$model)[c("(Intercept)", "age")]),
    tolerance = 1e-6
  )
})

test_that("wemix rows with missing or non-positive weights are dropped", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  d$w[1] <- NA
  d$w[2] <- 0
  d$w[3] <- -1
  expect_warning(
    m <- suppressMessages(
      fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
                 engine = "wemix", sampling_weights = "w")),
    "dropped 3 row"
  )
  expect_equal(nrow(m$data), nrow(d) - 3)
})

test_that("predict_maihda works for wemix fits (individual and strata)", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  m <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w"))

  p_ind <- predict_maihda(m, type = "individual")
  expect_length(p_ind, nrow(m$data))
  expect_true(all(is.finite(p_ind)))
  # Gaussian identity: link and response scales agree.
  expect_equal(predict_maihda(m, type = "individual", scale = "link"), p_ind)

  p_str <- predict_maihda(m, type = "strata")
  expect_true(all(c("stratum", "predicted", "se", "lower_95", "upper_95")
                  %in% names(p_str)))
  expect_equal(nrow(p_str), length(unique(as.character(m$data$stratum))))

  # newdata restricts the strata and errors on an unseen stratum.
  nd <- m$data[1:10, ]
  p_sub <- predict_maihda(m, newdata = nd, type = "strata")
  expect_setequal(as.character(p_sub$stratum), unique(as.character(nd$stratum)))
  nd_bad <- nd
  nd_bad$stratum <- "no:such:stratum"
  expect_error(predict_maihda(m, newdata = nd_bad, type = "strata"),
               "not present in the fitted model")

  # Predictions including the stratum effect track the conditional means.
  expect_gt(stats::cor(p_ind, m$data$y), 0.3)
})

test_that("calculate_pvc works across wemix fits and guards mismatched weights", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  m0 <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w"))
  madj <- suppressMessages(
    fit_maihda(y ~ age + gender + race + edu + (1 | stratum),
               data = m0$original_data, engine = "wemix",
               sampling_weights = "w"))

  pv <- calculate_pvc(m0, madj)
  expect_s3_class(pv, "pvc_result")
  expect_true(is.finite(pv$pvc))
  expect_gt(pv$pvc, 0)
  expect_lt(pv$pvc, 1)

  # Different sampling weights on the same rows are not comparable.
  d2 <- m0$original_data
  d2$w2 <- d2$w * c(2, 0.5)
  m_other <- suppressMessages(
    fit_maihda(y ~ age + gender + race + edu + (1 | stratum), data = d2,
               engine = "wemix", sampling_weights = "w2"))
  expect_error(calculate_pvc(m0, m_other), "same sampling weights")
})

test_that("maihda() runs the design-weighted two-model decomposition", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  a <- suppressMessages(
    maihda(y ~ age + (1 | gender:race:edu), data = d, sampling_weights = "w"))
  expect_s3_class(a, "maihda_analysis")
  expect_identical(a$mode, "two-model")
  expect_identical(a$model$engine, "wemix")
  expect_identical(a$model_adjusted$engine, "wemix")
  expect_false(is.null(a$pcv))
  expect_true(is.finite(a$pcv$pvc))
  expect_output(print(a), "PCV")
})

test_that("wemix fits a design-weighted logistic MAIHDA with weighted DA", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  m <- suppressMessages(suppressWarnings(
    fit_maihda(ybin ~ (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w")))
  expect_identical(m$family$family, "binomial")

  s <- summary(m)
  expect_gt(s$vpc$estimate, 0)
  expect_lt(s$vpc$estimate, 1)
  # Latent-scale level-1 variance pi^2/3, matching the other engines.
  vc <- s$variance_components
  expect_equal(vc$variance[vc$component == "Within-stratum (residual)"],
               pi^2 / 3, tolerance = 1e-8)

  da <- s$discriminatory_accuracy
  expect_false(is.null(da))
  expect_true(isTRUE(da$weighted))
  expect_gt(da$auc, 0.5)
  expect_gte(da$mor, 1)
  out <- paste(capture.output(print(da)), collapse = "\n")
  expect_match(out, "design-weighted")

  expect_true(is.finite(maihda_mor(m)))
})

test_that("plot types build on a wemix fit", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  skip_if_not_installed("ggplot2")

  d <- make_dw_data()
  m <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w"))
  s <- summary(m)
  for (t in c("vpc", "predicted", "obs_vs_shrunken", "risk_vs_effect",
              "effect_decomp")) {
    p <- plot(m, type = t, summary_obj = s)
    expect_s3_class(p, "ggplot")
  }
})

test_that("stepwise_pcv runs design-weighted steps on one analytic sample", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  s <- make_strata(d, vars = c("gender", "race", "edu"))
  out <- suppressMessages(
    stepwise_pcv(s$data, "y", c("gender", "race", "edu"),
                 sampling_weights = "w"))
  expect_s3_class(out, "maihda_stepwise")
  expect_equal(nrow(out), 4)
  expect_true(all(is.finite(out$Variance)))
  # The dimensions genuinely explain between-stratum variance here.
  expect_gt(out$Total_PCV[4], 0)
})

test_that("compare_maihda_groups fits design-weighted per-group models", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  d$country <- sample(c("X", "Y"), nrow(d), replace = TRUE)
  g <- suppressMessages(
    compare_maihda_groups(y ~ age + (1 | gender:race:edu), d,
                          group = "country", sampling_weights = "w"))
  expect_s3_class(g, "maihda_group_comparison")
  expect_setequal(g$group, c("X", "Y"))
  expect_true(all(is.finite(g$vpc)))
  expect_true(all(is.finite(g$pcv)))
})

test_that("compare_maihda warns when sampling weights differ across models", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  m_w <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w"))
  d2 <- m_w$original_data
  d2$w2 <- d2$w * c(2, 0.5)
  m_w2 <- suppressMessages(
    fit_maihda(y ~ age + (1 | stratum), data = d2,
               engine = "wemix", sampling_weights = "w2"))

  expect_warning(compare_maihda(m_w, m_w2), "sampling weights")
  # Same weights: no sampling-weight warning.
  expect_no_warning(compare_maihda(m_w, m_w))
})

test_that("summary(bootstrap = TRUE) is rejected for the wemix engine", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  m <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w"))
  expect_error(summary(m, bootstrap = TRUE, n_boot = 10), "replicate weights")
})
