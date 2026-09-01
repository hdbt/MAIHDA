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
               "strictly increasing")
  # Equal adjacent thresholds (a zero-probability category) are rejected too.
  expect_error(maihda_ordinal_category_probs(eta, c(0.4, 0.4), "logit"),
               "strictly increasing")
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

test_that("clmm prefilters on the evaluated frame, not raw complete.cases", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data(n = 400)
  d$pos <- stats::runif(nrow(d), 1, 5)
  # Two rows whose transformed predictor log(pos) is NaN. pos itself is non-missing, so
  # a raw complete.cases() on the formula columns keeps them -- the stored wrapper data
  # would then hold more rows than clmm actually fits, breaking downstream row alignment
  # (predictions / plots keyed off m$data). The evaluated analytic frame drops them.
  d$pos[c(4, 8)] <- -1

  expect_warning(
    m <- suppressMessages(fit_maihda(y ~ log(pos) + (1 | gender:race:edu),
                                     data = d, family = "ordinal")),
    "dropped 2 row"
  )
  expect_equal(nrow(m$data), nrow(d) - 2)
  # Stored analytic frame matches the rows clmm actually fit.
  expect_equal(nrow(m$data), as.integer(stats::nobs(m$model)))
  expect_true(all(is.finite(log(m$data$pos))))
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
  expect_error(suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 10)),
               "engine = \"brms\"")

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

  # Public individual predictions reject an unseen stratum by default rather than
  # silently returning a fixed-only (0-random-effect) prediction. The stratum is
  # supplied DIRECTLY, so the dimension columns are dropped: leaving them in place
  # would keep identifying the copied row's own stratum, which a supplied 'stratum'
  # must agree with (see the contradiction check in maihda_prepare_prediction_data).
  nd_unseen <- m$data[1, , drop = FALSE]
  nd_unseen[c("gender", "race", "edu")] <- NULL
  nd_unseen$stratum <- "999"
  expect_error(
    predict_maihda(m, newdata = nd_unseen, type = "individual", scale = "link"),
    "not present in the fitted model"
  )

  # allow_new_levels = TRUE opts into the population average: the latent location
  # equals the fixed part with the stratum random effect dropped (mapped to 0).
  pa_link <- predict_maihda(m, newdata = nd_unseen, type = "individual",
                            scale = "link", allow_new_levels = TRUE)
  fixed_only <- maihda_clmm_linpred(m, newdata = nd_unseen, include_re = FALSE)
  expect_true(is.finite(pa_link))
  expect_equal(unname(pa_link), unname(fixed_only))
})

test_that("clmm predictions add the formula offset (incl. an offset-only model)", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  d$expo <- runif(nrow(d), 0.5, 4)

  m <- fit_ord(d, formula = y ~ x + offset(log(expo)) + (1 | gender:race:edu))

  # The latent-location prediction must include the offset: rebuild it
  # independently as X %*% beta + offset + u and compare.
  off <- log(m$data$expo)
  beta <- m$model$beta
  tt <- stats::delete.response(stats::terms(reformulas::nobars(m$formula)))
  mf <- stats::model.frame(tt, m$data, na.action = stats::na.pass)
  X <- stats::model.matrix(tt, mf)
  ret <- maihda_clmm_stratum_ranef(m)
  re <- stats::setNames(ret$random_effect, ret$stratum)
  u <- re[as.character(m$data$stratum)]; u[is.na(u)] <- 0
  ref <- drop(X[, names(beta), drop = FALSE] %*% beta) + off + unname(u)

  eta <- predict_maihda(m, type = "individual", scale = "link")
  expect_equal(unname(eta), unname(ref), tolerance = 1e-8)
  # The offset is non-trivial, so a fit that dropped it would differ by max|offset|.
  expect_gt(max(abs(off)), 0.5)

  # An offset-only NULL model (no location coefficients) must still predict its
  # offset: the latent location is offset + u, not 0 + u.
  m0 <- fit_ord(d, formula = y ~ offset(log(expo)) + (1 | gender:race:edu))
  expect_length(m0$model$beta, 0L)
  ret0 <- maihda_clmm_stratum_ranef(m0)
  re0 <- stats::setNames(ret0$random_effect, ret0$stratum)
  u0 <- re0[as.character(m0$data$stratum)]; u0[is.na(u0)] <- 0
  eta0 <- predict_maihda(m0, type = "individual", scale = "link")
  expect_equal(unname(eta0), unname(log(m0$data$expo) + unname(u0)), tolerance = 1e-8)
  expect_gt(max(abs(eta0 - unname(u0))), 0.5)

  # The prediction-deviation panel rebuilds the clmm probabilities independently;
  # it must include the offset too.
  panel <- maihda_prediction_panel_ordinal_probs(m$model, m$data)
  exp_probs <- maihda_ordinal_category_probs(ref, m$model$alpha, m$model$link)
  panel_m <- as.matrix(panel); dimnames(panel_m) <- NULL
  exp_m <- exp_probs; dimnames(exp_m) <- NULL
  expect_equal(panel_m, exp_m, tolerance = 1e-8)
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

test_that("calculate_pcv and maihda() run the two-model ordinal decomposition", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_ord_data()
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ x + gender + race + edu + (1 | gender:race:edu),
           data = d, family = "ordinal")
  ))
  expect_s3_class(a, "maihda_analysis")
  expect_identical(a$model$engine, "ordinal")
  expect_true(is.finite(a$pcv$pcv))
  expect_output(print(a), "cumulative")

  # The standalone PCV over the pair agrees and the bootstrap is rejected with
  # the honest point-estimate guidance (not an engine-switch suggestion: no
  # engine computes a PCV interval).
  pcv <- calculate_pcv(a$model, a$model_adjusted)
  expect_equal(pcv$pcv, a$pcv$pcv)
  expect_error(
    suppressWarnings(
      calculate_pcv(a$model, a$model_adjusted, bootstrap = TRUE, n_boot = 10)),
    "point estimate"
  )

  # All-defaults path: the ordered factor selects the ordinal model end-to-end.
  a2 <- suppressMessages(suppressWarnings(maihda(y ~ x + (1 | gender:race), data = d)))
  expect_identical(a2$model$engine, "ordinal")
  expect_true(is.finite(a2$pcv$pcv))

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

test_that("the plot layer renders a clmm MAIHDA", {
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
  expect_s3_class(plot(m, type = "effect_decomp", summary_obj = s), "ggplot")
  expect_no_error(suppressWarnings(plot_prediction_deviation_panels(m)))
  expect_no_error(suppressWarnings(
    plot_prediction_deviation_panels(m, ordinal_mode = "expected_score")))
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

test_that("brms cumulative stratum predictions are expected category scores (Stan-free)", {
  skip_if_not_installed("brms")

  # Regression for the audit finding: the brms stratum-prediction helper applied
  # the SCALAR inverse link on the response scale, so for a cumulative (ordinal)
  # fit it returned a single cumulative probability in [0, 1] instead of the
  # EXPECTED CATEGORY SCORE sum_k k * P(Y = k) in [1, K] that maihda_table() and
  # the stratum plots document (and that the individual predict_maihda() path
  # already returns). It now maps the latent location to the score via the shared
  # cumulative helpers and posterior-mean thresholds.
  eta_fixed  <- c(0.10, 0.30, -0.20, 0.50)   # location mu (thresholds excluded)
  thresholds <- c(-0.5, 0.8)                  # K = 3 categories -> scores in [1, 3]
  strata     <- c("s1", "s1", "s2", "s2")

  # posterior_linpred(re_formula = NA) returns the location mu; fixef() carries
  # the thresholds as Intercept[1], Intercept[2] PLUS a location predictor row
  # 'x' that the threshold extractor must ignore. The mock returns what real
  # brms returns: the ndraws x nobs DRAWS matrix (no dimnames, no summary
  # columns; a summary= argument is ignored), which the package collapses via
  # maihda_brms_linpred_mean(). Symmetric +/- 0.05 jitter keeps the posterior
  # mean exactly eta_fixed while failing any non-averaging read of the matrix.
  local_mocked_bindings(
    posterior_linpred = function(object, ...) {
      rbind(
        matrix(rep(eta_fixed + 0.05, each = 20), nrow = 20),
        matrix(rep(eta_fixed - 0.05, each = 20), nrow = 20)
      )
    },
    fixef = function(object, ...) {
      matrix(c(thresholds, 0.4), ncol = 1,
             dimnames = list(c("Intercept[1]", "Intercept[2]", "x"), "Estimate"))
    },
    .package = "brms"
  )

  m <- structure(
    list(
      model = structure(list(family = list(family = "cumulative", link = "logit")),
                         class = "brmsfit"),
      engine = "brms",
      formula = y ~ x + (1 | stratum),
      data = data.frame(stratum = strata, stringsAsFactors = FALSE),
      family = list(family = "cumulative", link = "logit"),
      sampling_weights = NULL
    ),
    class = "maihda_model"
  )
  summ <- list(stratum_estimates = data.frame(
    stratum       = c("s1", "s2"),
    random_effect = c(0.2, -0.3),
    lower_95      = c(-0.1, -0.7),
    upper_95      = c(0.5,  0.1),
    stringsAsFactors = FALSE
  ))

  resp <- MAIHDA:::maihda_stratum_predictions_brms(m, summ, scale = "response")

  # Every reported value is an expected category score in [1, K = 3], NOT a
  # cumulative probability in [0, 1] (the bug). With these locations the scores
  # all exceed 1, which directly rules out the old plogis() output (always < 1).
  for (col in c("predicted_row", "lower_row", "upper_row", "fixed_row")) {
    expect_true(all(resp[[col]] >= 1 & resp[[col]] <= 3))
  }
  expect_true(all(resp$predicted_row > 1))

  # ...and they equal the per-stratum mean of the shared eta->score helper applied
  # row-wise (unit weights here) -- the same quantity the clmm path returns.
  idx <- match(strata, summ$stratum_estimates$stratum)
  u   <- summ$stratum_estimates$random_effect[idx]
  want_score <- tapply(
    MAIHDA:::maihda_ordinal_eta_to_score(eta_fixed + u, thresholds, "logit"),
    strata, mean
  )
  got_score <- stats::setNames(resp$predicted_row, as.character(resp$stratum))
  expect_equal(as.numeric(got_score[names(want_score)]),
               as.numeric(want_score), tolerance = 1e-10)

  # The link scale is untouched: the latent location mu + stratum effect, NOT a
  # score -- so the fix is scoped to scale = "response".
  lnk <- MAIHDA:::maihda_stratum_predictions_brms(m, summ, scale = "link")
  want_link <- tapply(eta_fixed + u, strata, mean)
  got_link  <- stats::setNames(lnk$predicted_row, as.character(lnk$stratum))
  expect_equal(as.numeric(got_link[names(want_link)]),
               as.numeric(want_link), tolerance = 1e-10)
})

test_that("brms cumulative summary returns a draws-based latent VPC", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). The latent
  # residual stubs and probability arithmetic are covered Stan-free above.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- make_ord_data(n = 600)
  # Sampler settings are load-bearing, not incidental. At the old chains = 2,
  # iter = 500 this fit reached max Rhat 1.08 with a bulk ESS of 26 on the
  # thresholds -- brms warned "some Rhats are > 1.05" out of the summary() call
  # below, and every assertion still passed, so the block proved only that the
  # code ran. iter = 2000 (1000 warmup) converges; keep the seed so the chains
  # are reproducible and the convergence assertion cannot flake.
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | gender:race:edu), data = d, family = "ordinal",
               engine = "brms", chains = 2, iter = 2000, warmup = 1000,
               refresh = 0, seed = 1)))
  # Guard the fixture: if a future brms/Stan release degrades the mixing, this
  # fails loudly instead of letting the weak assertions below pass on noise.
  # Also exercises maihda_fit_diagnostics()'s brms branch on a real Stan fit; the
  # verdict's direction and its 1.01 cut are pinned deterministically against mock
  # fits in test-audit-2026-07-26.R and test-audit-2026-08-31.R, which is what
  # keeps this guard from being vacuous. isTRUE() so a dropped or NULL
  # $diagnostics fails here rather than passing quietly.
  expect_true(isTRUE(m$diagnostics$converged))

  s <- summary(m)
  expect_true(s$vpc$estimate > 0 && s$vpc$estimate < 1)
  expect_true(is.finite(s$vpc$ci_lower) && is.finite(s$vpc$ci_upper))
  # Recovery, not just finiteness. Use the REALIZED between-stratum variance, not
  # the nominal sd_u: make_ord_data draws 12 effects at sd_u = 0.6 and this seed
  # realizes sd 0.8183 (var 0.6695), so against the logistic level-1 variance
  # pi^2/3 the latent VPC is 0.6695 / (0.6695 + 3.2899) = 0.169. A converged
  # posterior must cover it; the old chains = 2, iter = 500 fit could not support
  # a claim this sharp.
  expect_true(s$vpc$ci_lower < 0.169 && 0.169 < s$vpc$ci_upper)
  expect_true(s$vpc$estimate >= s$vpc$ci_lower && s$vpc$estimate <= s$vpc$ci_upper)

  # Response-scale predictions collapse the fitted() probability array to the
  # expected category score.
  sc <- predict_maihda(m, type = "individual", scale = "response")
  expect_true(all(sc >= 1 & sc <= 4))

  # A SINGLE prediction row must work too: the fitted() array's Estimate slice
  # drops to a bare vector for one row, which used to error (ncol() = NULL).
  sc1 <- predict_maihda(m, newdata = d[1, , drop = FALSE],
                        type = "individual", scale = "response")
  expect_length(sc1, 1L)
  expect_true(sc1 >= 1 && sc1 <= 4)
})
