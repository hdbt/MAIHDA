# Audit 2026-09-09 -- A3: structured (restricted) ordinal thresholds were read as
# the wrong number of categories.
#
# ordinal::clmm() stores the FREE threshold parameters in $alpha. Those equal the
# cut points only under the default threshold = "flexible". A structured request
# -- threshold = "equidistant" / "symmetric" / "symmetric2", part of clmm's
# documented API and reachable through fit_maihda()'s ... -- stores a shorter
# reparameterisation instead: an equidistant 5-category fit holds threshold.1 and
# spacing, two numbers, while the four cut points live in $Theta.
#
# Five sites read $alpha as if it were the cut points (individual predictions on
# the response scale, the summary threshold table, stratum response predictions,
# the proportional-odds bootstrap, and the prediction-deviation panels). Because
# threshold.1 < spacing often holds, the increasing-thresholds guard did NOT fire:
# a 5-category fit silently produced a 3-column probability matrix, an expected
# category score capped near 3 instead of 5, and a parametric bootstrap that
# simulated 3-category responses to calibrate a 5-category statistic.
#
# Fixed by maihda_clmm_cutpoints(), which expands the free parameters to the K-1
# cut points ($Theta, equivalently $tJac %*% $alpha) and asserts K-1 cut points for
# the K fitted categories. $tJac is the IDENTITY under "flexible", so the default
# path is bit-identical -- asserted below, since that is the only path most users
# take and a silent shift there would be a new defect.

# A 5-category ordinal fit; K = 5 makes the equidistant free-parameter count (2)
# differ from the cut-point count (4) by enough to be unmistakable.
maihda_a3_ord_data <- function(seed = 11, n = 900) {
  set.seed(seed)
  d <- data.frame(
    x  = stats::rnorm(n),
    d1 = sample(c("a", "b"), n, replace = TRUE),
    d2 = sample(c("p", "q"), n, replace = TRUE),
    d3 = sample(c("u", "v", "w"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  st  <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  eta <- 0.7 * d$x + stats::rnorm(nlevels(st), 0, 0.6)[st]
  cuts <- c(-0.8, 0.15, 1.10, 2.05)
  cum  <- vapply(cuts, function(a) stats::plogis(a - eta), numeric(n))
  cum  <- cbind(cum, 1)
  cum  <- t(apply(cum - cbind(0, cum[, -5, drop = FALSE]), 1, cumsum))
  cum[, 5] <- 1
  d$y <- factor(max.col(stats::runif(n) <= cum, ties.method = "first"),
                levels = 1:5, ordered = TRUE)
  d
}

maihda_a3_fit <- function(d, ...) {
  suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | d1:d2:d3), data = d, family = "ordinal", ...)
  ))
}

# ---- the defect: a structured threshold understated the category count ---------

test_that("structured thresholds expand to K-1 cut points, not the free parameters", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()

  for (th in c("equidistant", "symmetric", "symmetric2")) {
    m <- maihda_a3_fit(d, threshold = th)
    K <- nlevels(m$data$y)
    expect_identical(K, 5L)

    # The premise: clmm really does store fewer free parameters than cut points
    # here, so this is not a hypothetical branch.
    expect_lt(length(m$model$alpha), K - 1L)

    cut <- MAIHDA:::maihda_clmm_cutpoints(m$model)
    expect_length(cut, K - 1L)
    expect_false(is.unsorted(as.numeric(cut), strictly = TRUE))

    # They ARE the fit's own expanded cut points, both spellings.
    expect_equal(as.numeric(cut), as.numeric(m$model$Theta), tolerance = 0)
    expect_equal(as.numeric(cut),
                 as.numeric(m$model$tJac %*% as.numeric(m$model$alpha)),
                 tolerance = 1e-12)
  }
})

test_that("rebuilt category probabilities reproduce ordinal's own fitted values", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()

  for (th in c("flexible", "equidistant", "symmetric", "symmetric2")) {
    m   <- maihda_a3_fit(d, threshold = th)
    K   <- nlevels(m$data$y)
    eta <- MAIHDA:::maihda_clmm_linpred(m, include_re = TRUE)
    P   <- MAIHDA:::maihda_ordinal_category_probs(
      eta, MAIHDA:::maihda_clmm_cutpoints(m$model), m$family$link)

    # The invariant the finding asked for: one column per fitted category.
    expect_identical(ncol(P), K)
    expect_equal(rowSums(P), rep(1, nrow(P)), tolerance = 1e-12)

    # Independent oracle: ordinal computes fitted() = P(Y = y_i) internally from
    # the true cut points, whatever parameterisation was fitted. Pre-fix this was
    # off by up to 0.51 in probability under a structured threshold, and 254 of
    # 900 rows had an observed category with no column at all.
    idx <- cbind(seq_len(nrow(P)), as.integer(m$data$y))
    expect_equal(unname(P[idx]), unname(as.numeric(stats::fitted(m$model))),
                 tolerance = 1e-10)
  }
})

test_that("the response scale spans the full category range under a restricted fit", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()
  m <- maihda_a3_fit(d, threshold = "equidistant")
  K <- nlevels(m$data$y)

  score <- predict_maihda(m, type = "individual", scale = "response")
  expect_length(score, nrow(m$data))
  expect_true(all(score >= 1 & score <= K))

  # Pre-fix the score was capped near 3 (a 3-category expectation). It must now
  # reach well past the old ceiling.
  expect_gt(max(score), 3.5)

  # Independent check via the tail-sum identity E[Y] = 1 + sum_k P(Y > k), which
  # never forms the probability matrix at all.
  eta  <- MAIHDA:::maihda_clmm_linpred(m, include_re = TRUE)
  cut  <- MAIHDA:::maihda_clmm_cutpoints(m$model)
  tail <- 1 + rowSums(vapply(as.numeric(cut),
                             function(a) 1 - stats::plogis(a - eta),
                             numeric(length(eta))))
  expect_equal(unname(score), unname(tail), tolerance = 1e-12)

  # Stratum-level response predictions land on the same scale.
  sp <- MAIHDA:::maihda_stratum_predictions_ordinal(m, summary(m),
                                                    scale = "response")
  expect_true(all(sp$predicted >= 1 & sp$predicted <= K))
  expect_gt(max(sp$predicted), 2.5)
})

test_that("the summary threshold table reports the K-1 cut points", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()
  m <- maihda_a3_fit(d, threshold = "equidistant")
  K <- nlevels(m$data$y)

  tab <- summary(m)$thresholds
  expect_identical(nrow(tab), K - 1L)
  expect_identical(tab$term, c("1|2", "2|3", "3|4", "4|5"))
  expect_false(any(c("threshold.1", "spacing") %in% tab$term))
  expect_equal(tab$estimate, as.numeric(m$model$Theta), tolerance = 0)

  # Delta-method SEs for the expanded cut points: J V J' for the fit's Jacobian.
  V   <- stats::vcov(m$model)
  a   <- m$model$alpha
  ref <- sqrt(diag(m$model$tJac %*% V[names(a), names(a), drop = FALSE] %*%
                     t(m$model$tJac)))
  expect_equal(tab$se, as.numeric(ref), tolerance = 1e-12)
  expect_true(all(is.finite(tab$se) & tab$se > 0))

  # An equidistant fit constrains the spacing, so its cut-point SEs should be no
  # looser than the unconstrained fit's. A sanity check on the transform rather
  # than a discriminating one -- assert the lengths first so the comparison can
  # never quietly recycle a shorter vector into a passing result.
  flex_se <- summary(maihda_a3_fit(d))$thresholds$se
  expect_length(flex_se, length(tab$se))
  expect_true(all(tab$se <= flex_se + 1e-8))
})

test_that("the proportional-odds bootstrap simulates the fitted category count", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()
  m <- maihda_a3_fit(d, threshold = "equidistant")
  K <- nlevels(m$data$y)

  # The bootstrap draws responses from the category probabilities built on these
  # cut points; pre-fix that matrix had 3 columns, so every simulated statistic
  # was a 3-category LRT calibrating a 5-category observed one.
  cut <- MAIHDA:::maihda_clmm_cutpoints(m$model)
  expect_identical(length(cut) + 1L, K)

  po <- maihda_proportional_odds_test(m, n_sim = 9, seed = 3)
  expect_s3_class(po, "maihda_po_test")
  expect_identical(po$n_sim + po$n_failed, 9L)
  expect_true(all(is.finite(po$null_lrt)))
  expect_true(po$p_value > 0 && po$p_value <= 1)

  # The observed statistic tests K-2 extra threshold-specific slopes per term, so
  # a simulated 3-category draw could never reach it. The reference distribution
  # must now sit on the observed statistic's scale rather than an order of
  # magnitude below it.
  expect_identical(as.numeric(po$df), (K - 2) * po$n_terms)
  expect_gt(stats::median(po$null_lrt), 1)
})

# ---- confinement: the default flexible path must not move ---------------------

test_that("flexible thresholds are bit-identical to reading $alpha", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()

  for (link in c("logit", "probit")) {
    m <- suppressMessages(suppressWarnings(
      fit_maihda(y ~ x + (1 | d1:d2:d3), data = d,
                 family = maihda_cumulative(link))))

    alpha <- m$model$alpha
    cut   <- MAIHDA:::maihda_clmm_cutpoints(m$model)

    # tJac is the identity here, so the expansion is a no-op -- values AND names.
    expect_identical(as.numeric(cut), as.numeric(alpha))
    expect_identical(names(cut), names(alpha))
    expect_equal(m$model$tJac, diag(length(alpha)), ignore_attr = TRUE)

    eta <- MAIHDA:::maihda_clmm_linpred(m, include_re = TRUE)
    expect_identical(MAIHDA:::maihda_ordinal_category_probs(eta, cut, link),
                     MAIHDA:::maihda_ordinal_category_probs(eta, alpha, link))
    expect_identical(MAIHDA:::maihda_ordinal_eta_to_score(eta, cut, link),
                     MAIHDA:::maihda_ordinal_eta_to_score(eta, alpha, link))

    # The threshold table keeps the old labels, estimates and SEs exactly: the
    # delta method with J = I reproduces the untransformed vcov() diagonal.
    tab <- MAIHDA:::maihda_clmm_thresholds(m)
    expect_identical(tab$term, names(alpha))
    expect_identical(tab$estimate, as.numeric(alpha))
    expect_identical(tab$se,
                     as.numeric(sqrt(pmax(diag(stats::vcov(m$model)[
                       names(alpha), names(alpha), drop = FALSE]), 0))))
  }
})

test_that("cut-point recovery holds for a null model and a 3-category outcome", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()

  # Null (covariate-free) model: thresholds only, no beta.
  mn <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ 1 + (1 | d1:d2:d3), data = d, family = "ordinal")))
  expect_length(MAIHDA:::maihda_clmm_cutpoints(mn$model), nlevels(d$y) - 1L)
  expect_identical(nrow(MAIHDA:::maihda_clmm_thresholds(mn)), nlevels(d$y) - 1L)

  # 3 categories -> 2 cut points, equidistant (1 spacing + 1 anchor = 2 free
  # parameters, so alpha and Theta happen to agree in LENGTH but not in value).
  d3 <- d
  d3$y <- factor(pmin(as.integer(d$y), 3L), levels = 1:3, ordered = TRUE)
  m3 <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | d1:d2:d3), data = d3, family = "ordinal",
               threshold = "equidistant")))
  cut3 <- MAIHDA:::maihda_clmm_cutpoints(m3$model)
  expect_length(cut3, 2L)
  expect_equal(as.numeric(cut3), as.numeric(m3$model$Theta), tolerance = 0)
  # The second cut point is the anchor plus the spacing, NOT the spacing itself.
  expect_equal(as.numeric(cut3)[2],
               as.numeric(m3$model$alpha)[1] + as.numeric(m3$model$alpha)[2],
               tolerance = 1e-10)
})

# ---- the invariant ------------------------------------------------------------

test_that("a cut-point count inconsistent with the fitted categories is rejected", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- maihda_a3_ord_data()
  m <- maihda_a3_fit(d)

  bad <- m$model
  bad$Theta <- bad$Theta[, 1:2, drop = FALSE]
  expect_error(MAIHDA:::maihda_clmm_cutpoints(bad),
               "Cumulative fit inconsistency")

  # A fit carrying no thresholds at all still reports the original message.
  empty <- structure(list(), class = "clmm")
  expect_error(MAIHDA:::maihda_clmm_cutpoints(empty), "No thresholds")
})

test_that("an unnamed alpha degrades to NA SEs rather than a swallowed error", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # Hardening, not a reproduced field defect: ordinal names alpha on every fit
  # under all four threshold structures, so this state is not reachable from a
  # real fit. It is guarded because the natural guard is vacuously true --
  # names() of an unnamed vector is character(0), and all(logical(0)) is TRUE --
  # which let a 0 x 0 covariance block through to a non-conformable error that
  # summary()'s tryCatch swallows, silently dropping the threshold table.
  d <- maihda_a3_ord_data()
  m <- maihda_a3_fit(d, threshold = "equidistant")
  expect_false(is.null(names(m$model$alpha)))   # the premise: real fits DO name it

  bare <- m
  names(bare$model$alpha) <- NULL
  tab <- MAIHDA:::maihda_clmm_thresholds(bare)

  # The cut points themselves come from Theta and are unaffected...
  expect_identical(tab$estimate, MAIHDA:::maihda_clmm_thresholds(m)$estimate)
  expect_identical(tab$term, c("1|2", "2|3", "3|4", "4|5"))
  # ...only the standard errors, which need alpha's names to index vcov(), go NA.
  expect_true(all(is.na(tab$se)))

  # And summary() still renders the table rather than dropping it.
  s <- summary(bare)
  expect_identical(nrow(s$thresholds), 4L)
})
