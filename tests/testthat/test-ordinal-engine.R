# Ordinal (cumulative) MAIHDA: the clmm-based "ordinal" engine, the pure
# cumulative-probability helpers shared with the brms path, the latent-scale
# VPC (pi^2/3 logit / 1 probit), and the family <-> engine handshakes.

make_ord_data <- function(seed = 7321, n = 800, sd_u = 0.6) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("A", "B", "C"), n, replace = TRUE),
    edu = sample(c("low", "high"), n, replace = TRUE),
    x = rnorm(n),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$gender, d$race, d$edu, drop = TRUE)
  u <- rnorm(nlevels(stratum), sd = sd_u)[stratum]
  lat <- u + 0.3 * d$x + rlogis(n)
  d$y <- factor(cut(lat, c(-Inf, -1, 0.5, 2, Inf), labels = 1:4), ordered = TRUE)
  d
}

fit_ord <- function(d, formula = y ~ x + (1 | gender:race:edu), ...) {
  suppressMessages(suppressWarnings(
    fit_maihda(formula, data = d, family = "ordinal", ...)
  ))
}

# ---- pure cumulative-probability helpers (no ordinal package needed) ----------

test_that("maihda_ordinal_category_probs matches hand-computed plogis differences", {
  eta <- c(-1, 0, 0.7)
  alpha <- c(-0.5, 0.4, 1.2)

  probs <- maihda_ordinal_category_probs(eta, alpha, link = "logit")
  expect_equal(dim(probs), c(3L, 4L))
  expect_equal(rowSums(probs), rep(1, 3))

  cum <- sapply(alpha, function(a) stats::plogis(a - eta))
  hand <- cbind(cum[, 1], cum[, 2] - cum[, 1], cum[, 3] - cum[, 2], 1 - cum[, 3])
  expect_equal(unname(probs), unname(hand))

  # Probit variant uses pnorm.
  probs_p <- maihda_ordinal_category_probs(eta, alpha, link = "probit")
  cum_p <- sapply(alpha, function(a) stats::pnorm(a - eta))
  expect_equal(unname(probs_p[, 1]), unname(cum_p[, 1]))
  expect_equal(rowSums(probs_p), rep(1, 3))

  expect_error(maihda_ordinal_category_probs(eta, c(1, 0), "logit"),
               "non-decreasing")
  expect_error(maihda_ordinal_category_probs(eta, alpha, "cloglog"),
               "Unsupported cumulative link")
})

test_that("maihda_ordinal_expected_score scores categories 1..K", {
  probs <- rbind(c(1, 0, 0), c(0, 0, 1), c(0.5, 0.25, 0.25))
  expect_equal(maihda_ordinal_expected_score(probs), c(1, 3, 1.75))

  # eta -> score composition, bounded in [1, K] and increasing in eta.
  alpha <- c(-1, 0, 1)
  sc <- maihda_ordinal_eta_to_score(c(-5, 0, 5), alpha, "logit")
  expect_true(all(sc >= 1 & sc <= 4))
  expect_true(all(diff(sc) > 0))
})

test_that("maihda_cumulative and maihda_family_is_ordinal agree on the marker forms", {
  expect_identical(maihda_cumulative(), list(family = "cumulative", link = "logit"))
  expect_identical(maihda_cumulative("probit")$link, "probit")
  expect_error(maihda_cumulative("cloglog"))

  expect_true(maihda_family_is_ordinal("ordinal"))
  expect_true(maihda_family_is_ordinal("cumulative"))
  expect_true(maihda_family_is_ordinal(maihda_cumulative("probit")))
  expect_false(maihda_family_is_ordinal("gaussian"))
  expect_false(maihda_family_is_ordinal(stats::gaussian()))
  expect_false(maihda_family_is_ordinal(NULL))
})

test_that("maihda_response_is_ordinal requires an ordered factor with 3+ levels", {
  d <- data.frame(
    yord = factor(sample(1:4, 40, TRUE), ordered = TRUE),
    yfac = factor(sample(letters[1:4], 40, TRUE)),
    y2 = factor(sample(1:2, 40, TRUE), ordered = TRUE),
    ynum = rnorm(40)
  )
  expect_true(maihda_response_is_ordinal(yord ~ 1, d))
  expect_false(maihda_response_is_ordinal(yfac ~ 1, d))  # unordered
  expect_false(maihda_response_is_ordinal(y2 ~ 1, d))    # 2 levels -> binomial
  expect_false(maihda_response_is_ordinal(ynum ~ 1, d))
})

# ---- validation / handshake paths (no fit) ------------------------------------

test_that("the ordinal family <-> engine handshake rejects impossible combinations", {
  d <- make_ord_data(n = 60)

  expect_error(
    fit_maihda(y ~ x + (1 | gender:race), data = d, family = "ordinal",
               engine = "lme4"),
    "lme4 cannot fit a cumulative"
  )
  expect_error(
    fit_maihda(x ~ y + (1 | gender:race), data = d, engine = "ordinal"),
    "supply\\s+family"
  )
})

test_that("ordinal-engine guards reject context, sampling weights, and engine dots", {
  d <- make_ord_data(n = 60)
  d$site <- sample(c("s1", "s2"), nrow(d), replace = TRUE)
  d$w <- runif(nrow(d), 0.5, 2)

  expect_error(
    suppressMessages(fit_maihda(y ~ x + (1 | gender:race), data = d,
                                family = "ordinal", context = "site")),
    "does not support 'context'"
  )
  expect_error(
    fit_maihda(y ~ x + (1 | gender:race), data = d, family = "ordinal",
               engine = "ordinal", sampling_weights = "w"),
    "does not support 'sampling_weights'"
  )
  expect_error(
    suppressMessages(fit_maihda(y ~ x + (1 | gender:race), data = d,
                                family = "ordinal", subset = x > 0)),
    "not supported by engine = \"ordinal\""
  )
  # Explicit wemix falls through to the wemix family rejection.
  expect_error(
    fit_maihda(y ~ x + (1 | gender:race), data = d, family = "ordinal",
               engine = "wemix", sampling_weights = "w"),
    "gaussian\\(identity\\) and binomial\\(logit\\)"
  )
  # A non-logit/probit cumulative link is rejected up front.
  expect_error(
    fit_maihda(y ~ x + (1 | gender:race), data = d,
               family = list(family = "cumulative", link = "cloglog")),
    "logit and probit links"
  )
})

test_that("a 2-level ordered factor takes the binomial auto-detect path", {
  d <- make_ord_data(n = 200)
  d$y2 <- factor(ifelse(as.integer(d$y) <= 2, "low", "high"),
                 levels = c("low", "high"), ordered = TRUE)
  expect_warning(
    m <- suppressMessages(fit_maihda(y2 ~ (1 | gender:race), data = d)),
    "appears to be binary"
  )
  expect_identical(maihda_model_family_name(m), "binomial")
})

# ---- clmm fit, summary, predict ------------------------------------------------

test_that("fit_maihda fits a cumulative model via clmm with auto-switch and contract", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()

  # Explicit family + default engine -> auto-switch message.
  expect_message(
    m <- suppressWarnings(fit_maihda(y ~ x + (1 | gender:race:edu), data = d,
                                     family = "ordinal")),
    "engine = \"ordinal\""
  )
  expect_s3_class(m$model, "clmm")
  expect_identical(m$engine, "ordinal")
  expect_identical(m$family, list(family = "cumulative", link = "logit"))
  expect_true("stratum" %in% names(m$data))
  expect_identical(m$diagnostics$engine, "ordinal")
  expect_true(isTRUE(m$diagnostics$converged))

  # All defaults: the ordered factor is detected and selects the model.
  expect_warning(
    m2 <- suppressMessages(fit_maihda(y ~ x + (1 | gender:race), data = d)),
    "ordered factor"
  )
  expect_s3_class(m2$model, "clmm")

  # Numeric response with the ordinal family errors helpfully.
  expect_error(
    suppressMessages(fit_maihda(x ~ 1 + (1 | gender:race), data = d,
                                family = "ordinal")),
    "ordered-factor\\s+response"
  )
  # An unordered factor is coerced (declared order) with a message.
  d3 <- d
  d3$y <- factor(as.character(d3$y), levels = levels(d3$y))
  expect_message(
    suppressWarnings(fit_maihda(y ~ x + (1 | gender:race), data = d3,
                                family = "ordinal", engine = "ordinal")),
    "coercing the response"
  )
})

test_that("summary of a clmm MAIHDA reports the latent-scale VPC and thresholds", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  m <- fit_ord(d)
  s <- summary(m)

  expect_true(s$vpc$estimate > 0 && s$vpc$estimate < 1)
  resid_row <- s$variance_components$variance[
    s$variance_components$component == "Within-stratum (residual)"]
  expect_equal(resid_row, (pi^2) / 3)

  # VPC = stratum / (stratum + pi^2/3) from the clmm variance.
  v <- maihda_clmm_variances(m)
  expect_equal(s$vpc$estimate, v$stratum / (v$stratum + v$residual))

  # Thresholds: K - 1 = 3 rows with finite SEs, shown by print().
  expect_s3_class(s$thresholds, "data.frame")
  expect_identical(nrow(s$thresholds), 3L)
  expect_true(all(is.finite(s$thresholds$se)))
  expect_output(print(s), "Thresholds")

  # Location coefficient table carries Hessian SEs.
  expect_identical(s$fixed_effects$term, "x")
  expect_true(is.finite(s$fixed_effects$se))

  # DA stays binomial-only; bootstrap is rejected with the brms recommendation.
  expect_null(s$discriminatory_accuracy)
  expect_error(summary(m, bootstrap = TRUE, n_boot = 10), "engine = \"brms\"")

  # A probit fit uses latent residual variance 1.
  mp <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | gender:race:edu), data = d,
               family = maihda_cumulative("probit"))))
  sp <- summary(mp)
  resid_p <- sp$variance_components$variance[
    sp$variance_components$component == "Within-stratum (residual)"]
  expect_equal(resid_p, 1)
})

test_that("clmm predictions work on both scales and respect newdata strata", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  m <- fit_ord(d)
  n <- nrow(m$data)

  eta <- predict_maihda(m, type = "individual", scale = "link")
  sc <- predict_maihda(m, type = "individual", scale = "response")
  expect_length(eta, n)
  expect_length(sc, n)
  expect_true(all(sc >= 1 & sc <= 4))

  # The category probabilities behind the score reproduce clmm's own fitted()
  # values (the probability of each observed category).
  probs <- maihda_ordinal_category_probs(eta, m$model$alpha, "logit")
  resp <- m$data$y
  expect_equal(unname(probs[cbind(seq_len(n), as.integer(resp))]),
               unname(as.numeric(fitted(m$model))), tolerance = 1e-6)

  p_str <- predict_maihda(m, type = "strata")
  expect_identical(nrow(p_str), 12L)
  expect_true(all(c("predicted", "se", "lower_95", "upper_95") %in% names(p_str)))

  nd <- m$data[1:4, ]
  expect_length(predict_maihda(m, newdata = nd, type = "individual"), 4L)
  one <- predict_maihda(m, newdata = nd[1, , drop = FALSE], type = "strata")
  expect_identical(nrow(one), 1L)
  expect_error(
    predict_maihda(m, newdata = data.frame(stratum = "no-such"), type = "strata"),
    "not present in the fitted model"
  )
})

test_that("maihda_mor returns the median cumulative odds ratio for a logit fit", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  m <- fit_ord(d)
  v <- maihda_clmm_variances(m)$stratum
  expect_equal(maihda_mor(m), exp(sqrt(2 * v) * stats::qnorm(0.75)))

  mp <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | gender:race:edu), data = d,
               family = maihda_cumulative("probit"))))
  expect_error(maihda_mor(mp), "logit link")

  expect_error(maihda_vpc_response(m), "lme4 engine|binomial")
  expect_error(maihda_discriminatory_accuracy(m), "cumulative")
})

# ---- PCV / workflows -----------------------------------------------------------

test_that("calculate_pvc and maihda() run the two-model ordinal decomposition", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ x + gender + race + edu + (1 | gender:race:edu),
           data = d, family = "ordinal")
  ))
  expect_s3_class(a, "maihda_analysis")
  expect_identical(a$model$engine, "ordinal")
  expect_true(is.finite(a$pcv$pvc))
  expect_output(print(a), "cumulative")

  # The standalone PVC over the pair agrees and the bootstrap is rejected.
  pcv <- calculate_pvc(a$model, a$model_adjusted)
  expect_equal(pcv$pvc, a$pcv$pvc)
  expect_error(
    calculate_pvc(a$model, a$model_adjusted, bootstrap = TRUE, n_boot = 10),
    "engine = \"brms\""
  )

  # All-defaults path: the ordered factor selects the ordinal model end-to-end.
  a2 <- suppressMessages(suppressWarnings(maihda(y ~ x + (1 | gender:race), data = d)))
  expect_identical(a2$model$engine, "ordinal")
  expect_true(is.finite(a2$pcv$pvc))

  expect_error(
    suppressMessages(maihda(y ~ x + (1 | gender:race), data = d,
                            family = "ordinal",
                            decomposition = "crossed-dimensions")),
    "crossed random effects"
  )
  expect_error(
    suppressMessages(maihda(y ~ x + (1 | gender:race), data = d,
                            family = "ordinal", bootstrap = TRUE)),
    "Bootstrap intervals are not available for engine = \"ordinal\""
  )
})

test_that("compare_maihda_groups and stepwise_pcv support the ordinal engine", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  d$country <- sample(c("X", "Y"), nrow(d), replace = TRUE)

  gc <- suppressMessages(suppressWarnings(
    compare_maihda_groups(y ~ x + (1 | gender:race), data = d,
                          group = "country", min_group_n = 30)
  ))
  expect_s3_class(gc, "maihda_group_comparison")
  expect_identical(attr(gc, "engine"), "ordinal")
  expect_true(all(is.finite(gc$vpc[gc$status == "ok"])))

  expect_error(
    suppressMessages(compare_maihda_groups(
      y ~ x + (1 | gender:race), data = d, group = "country",
      family = "ordinal", decomposition = "crossed-dimensions")),
    "crossed random effects"
  )

  st <- make_strata(d, vars = c("gender", "race"))
  sw <- suppressMessages(suppressWarnings(stepwise_pcv(st$data, "y", c("gender", "race"))))
  expect_s3_class(sw, "maihda_stepwise")
  expect_true(all(is.finite(sw$Variance)))
})

# ---- plots ---------------------------------------------------------------------

test_that("the plot layer renders a clmm MAIHDA (and the ternary errors cleanly)", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  m <- fit_ord(d)
  s <- summary(m)

  expect_identical(maihda_prediction_panel_auto_type(m$model), "ordinal")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_s3_class(plot(m, type = "vpc", summary_obj = s), "ggplot")
  expect_s3_class(plot(m, type = "predicted", summary_obj = s), "ggplot")
  expect_s3_class(plot(m, type = "predicted", summary_obj = s, scale = "link"), "ggplot")
  expect_s3_class(plot(m, type = "obs_vs_shrunken", summary_obj = s), "ggplot")
  expect_s3_class(plot(m, type = "risk_vs_effect", summary_obj = s), "ggplot")
  expect_s3_class(plot(m, type = "effect_decomp", summary_obj = s), "ggplot")
  expect_no_error(suppressWarnings(plot_prediction_deviation_panels(m)))
  expect_no_error(suppressWarnings(
    plot_prediction_deviation_panels(m, ordinal_mode = "expected_score")))

  expect_error(compute_maihda_ternary_data(m), "not yet supported for the ordinal")
})

# ---- fake-fixture accessor branches --------------------------------------------

test_that("clmm accessors raise targeted errors on malformed fits", {
  skip_if_not_installed("ordinal")

  fake <- structure(
    list(model = structure(list(), class = "clmm"),
         family = maihda_cumulative(),
         data = data.frame(stratum = c("s1", "s2"))),
    class = "maihda_model"
  )
  expect_error(maihda_clmm_variances(fake), "Could not read the 'stratum'")
  expect_error(maihda_clmm_stratum_ranef(fake), "No 'stratum' random effects")
  expect_error(maihda_clmm_thresholds(fake), "No thresholds")
})

# ---- brms plumbing (Stan-free) --------------------------------------------------

test_that("fit_maihda routes the ordinal family to brms::cumulative()", {
  skip_if_not_installed("brms")

  d <- make_ord_data(n = 80)
  captured <- NULL
  local_mocked_bindings(
    brm = function(formula, data, family, ...) {
      captured <<- list(formula = formula, data = data, family = family)
      structure(list(), class = "brmsfit")
    },
    .package = "brms"
  )

  m <- fit_maihda(y ~ x + (1 | gender:race), data = d,
                  family = "ordinal", engine = "brms")
  expect_identical(captured$family$family, "cumulative")
  expect_identical(captured$family$link, "logit")
  expect_true(is.ordered(captured$data$y))
  expect_identical(m$engine, "brms")

  m2 <- fit_maihda(y ~ x + (1 | gender:race), data = d,
                   family = maihda_cumulative("probit"), engine = "brms")
  expect_identical(captured$family$link, "probit")
  expect_identical(maihda_normalize_family_name(m2$family$family), "cumulative")
})

test_that("brms cumulative summary returns a draws-based latent VPC", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). The latent
  # residual stubs and probability arithmetic are covered Stan-free above.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- make_ord_data(n = 600)
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | gender:race:edu), data = d, family = "ordinal",
               engine = "brms", chains = 2, iter = 500, refresh = 0, seed = 1)))
  s <- summary(m)
  expect_true(s$vpc$estimate > 0 && s$vpc$estimate < 1)
  expect_true(is.finite(s$vpc$ci_lower) && is.finite(s$vpc$ci_upper))

  # Response-scale predictions collapse the fitted() probability array to the
  # expected category score.
  sc <- predict_maihda(m, type = "individual", scale = "response")
  expect_true(all(sc >= 1 & sc <= 4))
})
