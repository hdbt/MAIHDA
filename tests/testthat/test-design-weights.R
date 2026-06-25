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

test_that("maihda_guard_reserved_weight_col rejects formula collisions only", {
  d <- data.frame(a = 1:3, .maihda_sw = 4:6, check.names = FALSE)

  # Present in data AND referenced by the formula -> would be overwritten -> error.
  expect_error(
    maihda_guard_reserved_weight_col(".maihda_sw", d, y ~ .maihda_sw, "brms"),
    "reserved internal column"
  )
  # Present in data but NOT referenced by the formula (the carry-along case, e.g.
  # maihda() refitting from a prior fit's $original_data) -> allowed.
  expect_true(maihda_guard_reserved_weight_col(".maihda_sw", d, y ~ a, "brms"))
  # Referenced by the formula but absent from the data -> no collision to guard
  # (an ordinary missing variable, caught later by the engine).
  expect_true(maihda_guard_reserved_weight_col(
    ".maihda_sw", data.frame(a = 1:3), y ~ .maihda_sw, "brms"))
})

test_that("the weighted engines reject a reserved column used in the formula", {
  d <- data.frame(y = rnorm(4), stratum = rep(c("a", "b"), 2),
                  w = c(1, 2, 3, 4),
                  .maihda_sw = c(10, 20, 30, 40),
                  .maihda_l2wt = c(50, 60, 70, 80),
                  check.names = FALSE)

  # brms path: a covariate named .maihda_sw would be clobbered by the normalized
  # likelihood weights, so fitting errors instead.
  expect_error(
    maihda_prepare_brms_sampling_weights(d, y ~ .maihda_sw + (1 | stratum), "w"),
    "reserved internal column"
  )
  # wemix path: the guard fires before WeMix::mix(), so no WeMix install needed.
  expect_error(
    maihda_fit_wemix(y ~ .maihda_l2wt + (1 | stratum), d,
                     stats::gaussian(), "w", list()),
    "reserved internal column"
  )

  # A reserved column merely present (not in the formula) is overwritten, not
  # rejected -- this keeps maihda()'s null/adjusted refits working when the column
  # rides along in $original_data.
  prep <- maihda_prepare_brms_sampling_weights(d, y ~ (1 | stratum), "w")
  expect_equal(mean(prep$data$.maihda_sw), 1)
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

test_that("wemix prefilters on the evaluated frame, not raw complete.cases", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data(n = 400)
  # Two rows whose transformed predictor log(age) is NaN. age itself is non-missing, so
  # a raw complete.cases() on the formula columns keeps them -- leaving the sampling
  # weight vector longer than the model matrix and erroring in WeMix::mix() with a
  # weight/X row-count mismatch. The evaluated analytic frame drops them up front.
  d$age[c(5, 9)] <- -1

  expect_warning(
    m <- suppressMessages(
      fit_maihda(y ~ log(age) + (1 | gender:race:edu), data = d,
                 engine = "wemix", sampling_weights = "w")),
    "dropped 2 row"
  )
  expect_equal(nrow(m$data), nrow(d) - 2)
  expect_true(all(is.finite(log(m$data$age))))   # no NaN-producing rows survived
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

test_that("wemix predictions add the formula offset", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  d <- make_dw_data()
  d$expo <- runif(nrow(d), 0.5, 4)
  m <- suppressMessages(
    fit_maihda(y ~ age + offset(log(expo)) + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w"))

  # Rebuild the link-scale prediction independently as X %*% beta + offset + u.
  off <- log(m$data$expo)
  beta <- m$model$coef
  tt <- stats::delete.response(stats::terms(reformulas::nobars(m$formula)))
  mf <- stats::model.frame(tt, m$data, na.action = stats::na.pass)
  X <- stats::model.matrix(tt, mf)
  re <- maihda_wemix_ranef_vector(m)
  u <- re[as.character(m$data$stratum)]; u[is.na(u)] <- 0
  ref <- drop(X[, names(beta), drop = FALSE] %*% beta) + off + unname(u)

  eta <- predict_maihda(m, type = "link")
  expect_equal(unname(eta), unname(ref), tolerance = 1e-8)
  # A fit that dropped the offset would differ by max|offset| (> 0.5 here).
  expect_gt(max(abs(off)), 0.5)
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

test_that("calculate_pvc/compare_maihda catch DIFFERENT wemix outcomes", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  # WeMixResults exposes no model frame, so the n/row/outcome checks must read the
  # wrapper's stored analytic $data. Two fits sharing formula/n/strata/weights/rows
  # but holding different OUTCOME values must not pass validation silently.
  d <- make_dw_data()
  m0 <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               engine = "wemix", sampling_weights = "w"))
  d_diff <- d
  d_diff$y <- d_diff$y + 100 + rnorm(nrow(d_diff), sd = 5)
  m_diff <- suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d_diff,
               engine = "wemix", sampling_weights = "w"))

  expect_error(calculate_pvc(m0, m_diff), "outcome values differ")
  expect_warning(compare_maihda(m0, m_diff), "analytic sample")
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

  # Stratum-level predictions on both scales; for a gaussian identity model the
  # link and response scales coincide.
  sp_resp <- maihda_stratum_predictions_wemix(m, s, scale = "response")
  sp_link <- maihda_stratum_predictions_wemix(m, s, scale = "link")
  expect_equal(sp_link$predicted_row, sp_resp$predicted_row)
  expect_true(all(c("stratum", "predicted_row", "n", "w_sum") %in% names(sp_link)))
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

  # Crossed-dimensions / contextual partitions are undefined for wemix; a model
  # that somehow carries those tags is rejected by summary().
  m_cc <- m
  m_cc$cc_info <- list(dim_groups = c(g = "g"), interaction_group = "stratum")
  expect_error(summary(m_cc), "crossed random effects")
  m_ctx <- m
  m_ctx$context_info <- list(groups = "site")
  expect_error(summary(m_ctx), "crossed random effects")
})

test_that("engine = wemix rejects aggregated binomial responses", {
  d <- make_dw_data()
  set.seed(1)
  d$succ <- rbinom(nrow(d), 5, 0.4)
  d$fail <- 5 - d$succ
  s <- make_strata(d, vars = c("gender", "race", "edu"))
  expect_error(
    fit_maihda(cbind(succ, fail) ~ (1 | stratum), data = s$data,
               engine = "wemix", family = "binomial", sampling_weights = "w"),
    "aggregated binomial"
  )
})

# ---- internal helper error branches (no WeMix fit required) -------------------

# A minimal stand-in for a fitted wemix maihda_model, enough for the variance /
# random-effect / prediction helpers to read.
make_fake_wemix_model <- function(varDF = NULL, ranefMat = NULL, coef = NULL,
                                  data = NULL, family = stats::gaussian(),
                                  formula = y ~ x + (1 | stratum)) {
  structure(
    list(
      model = structure(list(varDF = varDF, ranefMat = ranefMat, coef = coef),
                        class = "WeMixResults"),
      engine = "wemix",
      formula = formula,
      data = data,
      family = family,
      sampling_weights = "w"
    ),
    class = "maihda_model"
  )
}

test_that("maihda_wemix_variances validates the WeMix variance table", {
  expect_error(maihda_wemix_variances(make_fake_wemix_model(varDF = NULL)),
               "Could not read the variance components")

  no_stratum <- data.frame(grp = "Residual", var1 = NA, vcov = 1)
  expect_error(maihda_wemix_variances(make_fake_wemix_model(varDF = no_stratum)),
               "No 'stratum' random-effect variance")

  no_resid <- data.frame(grp = "stratum", var1 = "(Intercept)", vcov = 0.5)
  expect_error(maihda_wemix_variances(make_fake_wemix_model(varDF = no_resid)),
               "No residual variance")

  # Binomial-logit: the level-1 variance is the latent pi^2/3, no Residual row needed.
  m_bin <- make_fake_wemix_model(varDF = no_resid,
                                 family = stats::binomial(link = "logit"))
  v <- maihda_wemix_variances(m_bin)
  expect_equal(v$stratum, 0.5)
  expect_equal(v$residual, pi^2 / 3)
})

test_that("maihda_wemix_ranef_vector validates the random-effect table", {
  expect_error(maihda_wemix_ranef_vector(make_fake_wemix_model(ranefMat = NULL)),
               "No 'stratum' random effects")

  no_intercept <- list(stratum = data.frame(slope = c(0.1, -0.1),
                                            row.names = c("s1", "s2")))
  expect_error(
    maihda_wemix_ranef_vector(make_fake_wemix_model(ranefMat = no_intercept)),
    "must include an intercept"
  )
})

test_that("maihda_wemix_stratum_ranef collapses the SE at a boundary fit", {
  # Zero between-stratum variance: the conditional distribution of every stratum
  # effect collapses on 0, so the SE is 0 rather than undefined.
  ranefMat <- list(stratum = stats::setNames(
    data.frame(c(0, 0), row.names = c("s1", "s2")), "(Intercept)"))
  varDF <- data.frame(grp = c("stratum", "Residual"),
                      var1 = c("(Intercept)", NA), vcov = c(0, 1))
  d <- data.frame(stratum = c("s1", "s1", "s2"), w = c(1, 2, 1))
  m <- make_fake_wemix_model(varDF = varDF, ranefMat = ranefMat, data = d)

  out <- maihda_wemix_stratum_ranef(m)
  expect_equal(out$se, c(0, 0))
  expect_equal(out$random_effect, c(0, 0))
})

test_that("maihda_wemix_linpred errors when the design matrix cannot be rebuilt", {
  d <- data.frame(y = rnorm(4), x = rnorm(4), stratum = c("s1", "s2", "s1", "s2"))
  m <- make_fake_wemix_model(coef = c("(Intercept)" = 1, zzz = 2), data = d)
  expect_error(maihda_wemix_linpred(m, include_re = FALSE),
               "missing column")
})

test_that("maihda_fit_wemix errors when no usable rows remain", {
  skip_if_not_installed("WeMix")
  d <- data.frame(y = rnorm(4), stratum = c("s1", "s2", "s1", "s2"),
                  w = c(0, -1, NA, 0))
  expect_error(
    maihda_fit_wemix(y ~ (1 | stratum), d, stats::gaussian(), "w", list()),
    "No usable rows remain"
  )
})

test_that("maihda_stratum_predictions_wemix validates its inputs", {
  d <- data.frame(y = rnorm(4), x = rnorm(4), w = 1,
                  stratum = c("s1", "s2", "s1", "s2"))
  m <- make_fake_wemix_model(coef = c("(Intercept)" = 1, x = 0.5), data = d)
  no_stratum <- m
  no_stratum$data$stratum <- NULL
  expect_error(
    maihda_stratum_predictions_wemix(no_stratum, list(stratum_estimates = NULL)),
    "'stratum' variable not found"
  )
  expect_error(
    maihda_stratum_predictions_wemix(m, list(stratum_estimates = NULL)),
    "No stratum estimates"
  )
})

test_that("maihda_fit_diagnostics flags a boundary (singular) wemix fit", {
  singular <- structure(
    list(varDF = data.frame(grp = c("stratum", "Residual"),
                            var1 = c("(Intercept)", NA), vcov = c(0, 1))),
    class = "WeMixResults")
  diag_s <- maihda_fit_diagnostics(singular)
  expect_identical(diag_s$engine, "wemix")
  expect_true(isTRUE(diag_s$converged))
  expect_true(isTRUE(diag_s$singular))

  healthy <- structure(
    list(varDF = data.frame(grp = c("stratum", "Residual"),
                            var1 = c("(Intercept)", NA), vcov = c(0.4, 1))),
    class = "WeMixResults")
  expect_false(isTRUE(maihda_fit_diagnostics(healthy)$singular))
})

test_that("print.maihda_model labels brms sampling weights as pseudo-posterior", {
  d <- data.frame(y = rnorm(10), x = rnorm(10))
  m <- structure(
    list(
      model = stats::lm(y ~ x, data = d),
      engine = "brms",
      formula = y ~ x + (1 | stratum),
      family = list(family = "gaussian", link = "identity"),
      sampling_weights = "w",
      diagnostics = NULL
    ),
    class = "maihda_model"
  )
  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(out, "Sampling weights: w")
  expect_match(out, "pseudo-posterior")
})

test_that("maihda_prior_weights resolves design weights from either column", {
  d <- data.frame(y = rnorm(4), stratum = c("s1", "s2", "s1", "s2"),
                  w = c(2, 3, 4, 5))
  m <- make_fake_wemix_model(data = d)

  expect_equal(maihda_prior_weights(m), c(2, 3, 4, 5))

  # For a non-aggregated-binomial response the prediction weights coincide with the
  # prior weights (trial-weighting only kicks in for a cbind/trials() response).
  expect_equal(maihda_prediction_weights(m), maihda_prior_weights(m))

  # A brms analytic frame carries the normalized .maihda_sw column instead.
  m_brms <- m
  m_brms$data$w <- NULL
  m_brms$data$.maihda_sw <- c(1, 1.5, 0.5, 1)
  expect_equal(maihda_prior_weights(m_brms), c(1, 1.5, 0.5, 1))

  # Neither column resolvable: degrade to unit weights.
  m_none <- m
  m_none$data$w <- NULL
  expect_equal(maihda_prior_weights(m_none), rep(1, 4))
})

test_that("the sampling-weight fingerprint degrades when weights are unrecoverable", {
  d <- data.frame(y = rnorm(3), stratum = c("s1", "s2", "s1"))
  m <- make_fake_wemix_model(data = d)
  # sampling_weights is recorded as "w" but no weight column survives in the
  # analytic frame: fall back to the column name.
  expect_identical(maihda_sampling_weight_fingerprint(m), "col:w")

  m_unweighted <- m
  m_unweighted$sampling_weights <- NULL
  expect_identical(maihda_sampling_weight_fingerprint(m_unweighted), "none")
})

test_that("compare_maihda_groups rejects an unknown engine", {
  d <- make_dw_data()
  d$country <- sample(c("X", "Y"), nrow(d), replace = TRUE)
  expect_error(
    compare_maihda_groups(y ~ age + (1 | gender:race), d, group = "country",
                          engine = "nope"),
    "lme4, brms, wemix"
  )
})

test_that("fit_maihda routes sampling weights into a brms fit (Stan-free)", {
  skip_if_not_installed("brms")

  d <- make_dw_data()
  s <- make_strata(d, vars = c("gender", "race", "edu"))
  captured <- NULL
  local_mocked_bindings(
    brm = function(formula, data, family, ...) {
      captured <<- list(formula = formula, data = data, family = family)
      structure(list(), class = "brmsfit")
    },
    .package = "brms"
  )

  expect_message(
    m <- fit_maihda(y ~ age + (1 | stratum), data = s$data,
                    engine = "brms", sampling_weights = "w"),
    "pseudo-posterior"
  )
  expect_identical(m$engine, "brms")
  expect_identical(m$sampling_weights, "w")
  # The formula gained the weights() addition term, and the data the normalized
  # (mean-1) weight column the term references.
  expect_match(paste(deparse(captured$formula), collapse = " "),
               "weights(.maihda_sw)", fixed = TRUE)
  expect_true(".maihda_sw" %in% names(captured$data))
  expect_equal(mean(captured$data$.maihda_sw), 1)
  expect_equal(captured$data$.maihda_sw / captured$data$.maihda_sw[1],
               captured$data$w / captured$data$w[1])
})

test_that("maihda_reslice_dot_args slices only row-aligned forwarded values", {
  env <- new.env()
  dot_vals <- list(subset = 1:5, chains = 4L, mat = matrix(1:10, nrow = 5))
  for (nm in names(dot_vals)) {
    assign(paste0(".maihda_arg_", nm), dot_vals[[nm]], envir = env)
  }
  keep <- c(TRUE, FALSE, TRUE, FALSE, TRUE)

  maihda_reslice_dot_args(dot_vals, keep, 5L, env)

  # Row-aligned vector and matrix are sliced to the kept rows...
  expect_equal(get(".maihda_arg_subset", env), c(1L, 3L, 5L))
  expect_equal(get(".maihda_arg_mat", env),
               matrix(1:10, nrow = 5)[keep, , drop = FALSE])
  # ...a scalar (length != n) is left untouched.
  expect_equal(get(".maihda_arg_chains", env), 4L)
})

test_that("fit_maihda re-slices row-aligned dots after the brms weight drop", {
  skip_if_not_installed("brms")

  d <- make_dw_data(n = 12)
  # Force two rows out of the fit via a zero and a missing weight.
  d$w[3] <- 0
  d$w[7] <- NA
  s <- make_strata(d, vars = c("gender", "race", "edu"))

  captured <- NULL
  local_mocked_bindings(
    brm = function(formula, data, family, ...) {
      captured <<- list(data = data, dots = list(...))
      structure(list(), class = "brmsfit")
    },
    .package = "brms"
  )

  # `subset` is forwarded and row-aligned to the pre-drop data (12 rows); after the
  # weight drop the data is 10 rows, so the bound subset must be re-sliced to match
  # -- otherwise brms gets a 12-long subset against a 10-row data frame.
  sub <- rep(TRUE, nrow(s$data))
  suppressMessages(suppressWarnings(
    fit_maihda(y ~ age + (1 | stratum), data = s$data,
               engine = "brms", sampling_weights = "w", subset = sub)
  ))
  expect_equal(nrow(captured$data), 10L)
  expect_length(captured$dots$subset, 10L)
})

test_that("wemix unseen stratum: helper maps to zero, public path gates it", {
  ranefMat <- list(stratum = stats::setNames(
    data.frame(c(0.5, -0.5), row.names = c("s1", "s2")), "(Intercept)"))
  varDF <- data.frame(grp = c("stratum", "Residual"),
                      var1 = c("(Intercept)", NA), vcov = c(0.3, 1))
  d <- data.frame(y = rnorm(4), x = c(0, 1, 0, 1), w = 1,
                  stratum = c("s1", "s2", "s1", "s2"))
  m <- make_fake_wemix_model(varDF = varDF, ranefMat = ranefMat,
                             coef = c("(Intercept)" = 2, x = 1), data = d)

  nd <- data.frame(x = c(0, 0), stratum = c("s1", "unseen"))

  # The internal linpred primitive maps an unseen stratum to a zero random effect
  # (eta = fixed part). This is the population-average fallback, NOT the public
  # default: the seen stratum keeps its 0.5 effect, the unseen one drops to 0.
  eta <- maihda_wemix_linpred(m, newdata = nd, include_re = TRUE)
  expect_equal(unname(eta), c(2.5, 2))

  # Public individual predictions reject the unseen stratum by default, instead of
  # silently returning that fixed-only prediction.
  expect_error(
    predict_maihda(m, newdata = nd, type = "individual"),
    "not present in the fitted model"
  )

  # allow_new_levels = TRUE opts into the population average. The gaussian identity
  # link makes the response scale equal the latent eta computed above.
  pa <- predict_maihda(m, newdata = nd, type = "individual",
                       allow_new_levels = TRUE)
  expect_equal(unname(pa), c(2.5, 2))
})

test_that("predict_maihda supplies the brms weight column for user newdata", {
  skip_if_not_installed("brms")

  ranef_table <- data.frame(
    stratum = c("s1", "s2"),
    stratum_id = c(1L, 2L),
    random_effect = c(0.2, -0.2),
    se = c(0.1, 0.1),
    lower_95 = c(0, -0.4),
    upper_95 = c(0.4, 0),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    maihda_stratum_ranef_brms = function(model, group = "stratum") ranef_table
  )

  d <- data.frame(y = rnorm(4), stratum = c("s1", "s2", "s1", "s2"),
                  .maihda_sw = 1, check.names = FALSE)
  m <- structure(
    list(
      model = structure(list(), class = "brmsfit"),
      engine = "brms",
      formula = y | weights(.maihda_sw) ~ (1 | stratum),
      data = d,
      family = list(family = "gaussian", link = "identity"),
      sampling_weights = "w"
    ),
    class = "maihda_model"
  )

  # newdata lacking the internal weight column gets a unit weight injected and
  # the strata table is filtered to the requested stratum.
  out <- predict_maihda(m, newdata = data.frame(stratum = "s1"), type = "strata")
  expect_equal(as.character(out$stratum), "s1")
  expect_equal(out$predicted, 0.2)
})
