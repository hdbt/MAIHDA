test_that("maihda_auc matches hand-computed rank-based AUC", {
  # Perfect separation -> 1; reversed -> 0; all-ties -> 0.5.
  expect_equal(maihda_auc(c(0.1, 0.2, 0.8, 0.9), c(0, 0, 1, 1)), 1)
  expect_equal(maihda_auc(c(0.9, 0.8, 0.2, 0.1), c(0, 0, 1, 1)), 0)
  expect_equal(maihda_auc(c(0.5, 0.5, 0.5, 0.5), c(0, 1, 0, 1)), 0.5)

  # A worked example: P(case score > control score), ties = 0.5.
  # cases {0.35, 0.8}, controls {0.1, 0.4}: pairs (0.35>0.1)T, (0.35>0.4)F,
  # (0.8>0.1)T, (0.8>0.4)T => 3/4 = 0.75.
  expect_equal(maihda_auc(c(0.1, 0.4, 0.35, 0.8), c(0, 0, 1, 1)), 0.75)

  # Logical y is accepted; equals the 0/1 version.
  expect_equal(maihda_auc(c(0.1, 0.4, 0.35, 0.8), c(FALSE, FALSE, TRUE, TRUE)), 0.75)
})

test_that("maihda_auc validates inputs and handles degenerate classes", {
  expect_error(maihda_auc(c(0.1, 0.2), c(1, 0, 1)), "same length")
  expect_error(maihda_auc("a", c(0, 1)), "numeric")
  expect_error(maihda_auc(c(0.1, 0.2), c(1, 2)), "binary")
  # Only one class present -> undefined AUC.
  expect_true(is.na(maihda_auc(c(0.1, 0.2, 0.3), c(1, 1, 1))))
})

test_that("maihda_auc does not integer-overflow on large samples", {
  # Regression: n1 * n0 here is 5.625e9, past .Machine$integer.max (~2.15e9). With
  # integer arithmetic (sum() of a logical is integer) the denominator overflowed to
  # NA; double arithmetic fixes it. Perfect separation gives AUC exactly 1.
  n_each <- 75000L
  expect_gt(as.double(n_each)^2, .Machine$integer.max)  # the test's own premise
  y <- rep(c(0L, 1L), each = n_each)
  prob <- rep(c(0.2, 0.8), each = n_each)
  expect_equal(maihda_auc(prob, y), 1)
})

test_that("maihda_auc_weighted does not integer-overflow on large integer counts", {
  # Same overflow guard for the count-weighted path: 75000 cases vs 75000 controls,
  # perfectly separated, as large integer counts (n1 * n0 = 5.625e9).
  expect_equal(
    maihda_auc_weighted(prob = c(0.2, 0.8),
                        successes = c(0L, 75000L),
                        trials    = c(75000L, 75000L)),
    1)
})

# A binomial MAIHDA model on synthetic data with a real between-stratum signal.
maihda_da_test_model <- function(seed = 123, n = 900, family = "binomial") {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  lp <- stats::rnorm(nlevels(sk), sd = 0.8)[sk]
  d$y <- stats::rbinom(n, 1, stats::plogis(lp))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = family)
  ))
}

test_that("maihda_mor equals exp(sqrt(2 * V_A) * qnorm(0.75)) for a logistic fit", {
  m <- maihda_da_test_model()
  v_a <- MAIHDA:::extract_between_variance(m)
  expect_equal(maihda_mor(m), exp(sqrt(2 * v_a) * stats::qnorm(0.75)))
  expect_gte(maihda_mor(m), 1)  # MOR is >= 1 by construction
})

test_that("maihda_mor errors for a non-logit binomial link (probit)", {
  # The MOR is an odds-ratio-scale quantity defined only for the logit link; a probit
  # fit must be rejected rather than returning exp(...) off the model's scale.
  m <- maihda_da_test_model(family = stats::binomial("probit"))
  expect_identical(MAIHDA:::maihda_model_link_name(m), "probit")
  expect_error(maihda_mor(m), "logit link")
})

test_that("maihda_discriminatory_accuracy reports AUC with mor = NA for a probit fit", {
  m <- maihda_da_test_model(family = stats::binomial("probit"))
  da <- maihda_discriminatory_accuracy(m)

  expect_s3_class(da, "maihda_da")
  expect_true(is.finite(da$auc) && da$auc >= 0 && da$auc <= 1)  # AUC is link-agnostic
  expect_true(is.na(da$mor))                                    # MOR undefined for probit
  expect_identical(da$link, "probit")
  expect_identical(da$family, "binomial")
  # print() explains WHY the MOR is NA rather than showing a bare NA
  expect_output(print(da), "requires the logit link")
})

test_that("maihda_mor errors for non-binomial models", {
  g <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(120), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    family = "gaussian"
  )))
  expect_error(maihda_mor(g$model), "binomial")
})

test_that("maihda_discriminatory_accuracy bundles AUC + MOR and reproduces the vignette computation", {
  m <- maihda_da_test_model()
  da <- maihda_discriminatory_accuracy(m)

  expect_s3_class(da, "maihda_da")
  expect_true(da$auc >= 0 && da$auc <= 1)
  expect_equal(da$mor, maihda_mor(m))
  expect_equal(da$n_case + da$n_control, nrow(m$data))
  expect_identical(da$family, "binomial")

  # Equivalence to the binary_outcomes vignette's hand-rolled AUC: the exported
  # function must return exactly what the documented one-liner produces.
  prob <- predict_maihda(m, type = "individual", scale = "response")
  y_obs <- as.numeric(lme4::getME(m$model, "y"))
  expect_equal(da$auc, maihda_auc(prob, y_obs))
})

test_that("discriminatory accuracy flags the AUC as apparent / in-sample", {
  m <- maihda_da_test_model()
  da <- maihda_discriminatory_accuracy(m)
  # The AUC is scored on the fitting rows, so it is apparent (optimistically biased).
  expect_true(isTRUE(da$apparent))
  # print() says so, so the value is not mistaken for out-of-sample discrimination.
  expect_output(print(da), "apparent")
})

test_that("maihda_discriminatory_accuracy rejects non-binomial models", {
  g <- suppressWarnings(suppressMessages(MAIHDA:::maihda_app_fit_models(
    MAIHDA::maihda_sim_data[seq_len(120), ],
    outcome_var = "health_outcome", grouping_vars = c("gender", "race"),
    family = "gaussian"
  )))
  expect_error(maihda_discriminatory_accuracy(g$model), "binomial")
})

test_that("maihda_auc_weighted equals the rank AUC on the expanded 0/1 data", {
  # Three probability levels with shared cases/controls; ties at equal probability
  # are counted as one half, exactly as maihda_auc() does on the expanded vectors.
  prob <- c(0.2, 0.5, 0.8)
  successes <- c(1, 2, 4)
  trials <- c(3, 4, 5)
  failures <- trials - successes

  expanded_prob <- unlist(Map(function(p, s, f) rep(p, s + f), prob, successes, failures))
  expanded_y <- unlist(Map(function(s, f) c(rep(1, s), rep(0, f)), successes, failures))

  expect_equal(MAIHDA:::maihda_auc_weighted(prob, successes, trials),
               maihda_auc(expanded_prob, expanded_y))

  # A degenerate class (no failures anywhere) yields NA, like maihda_auc().
  expect_true(is.na(MAIHDA:::maihda_auc_weighted(prob, trials, trials)))

  # Negative implied mass (successes > trials) is an invariant violation -- it
  # is how the fractional-weight defect produced an AUC above 1 -- and errors
  # rather than returning an out-of-bounds concordance.
  expect_error(MAIHDA:::maihda_auc_weighted(c(0.2, 0.8), c(2, 1), c(1.5, 1)),
               "negative case/control mass")
})

test_that("lme4 precision weights are ignored for the AUC (ordinary observation-level)", {
  # Audit finding 2: lme4 prior weights on a Bernoulli fit are PRECISION weights --
  # they scale the observation's likelihood/dispersion, not its population frequency
  # (a weight of 1.5 is not 1.5 population members). They must NOT be folded into
  # case/control mass, which reports a weighted concordance with no population-AUC
  # interpretation and silently changes the estimand based on a fitting control. The
  # AUC is now the ordinary observation-level concordance; the fit is flagged
  # precision_weights_ignored but reported as unweighted.
  skip_on_cran()
  set.seed(404)
  n <- 400
  d <- data.frame(gender = sample(c("F", "M"), n, replace = TRUE),
                  race = sample(c("A", "B"), n, replace = TRUE))
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- stats::rbinom(n, 1, stats::plogis(stats::rnorm(nlevels(sk), sd = 1.2)[sk]))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum

  # Non-uniform fractional precision weights: the AUC must equal the UNWEIGHTED rank
  # AUC, and must NOT be the (old) precision-weighted concordance.
  d$w <- stats::runif(n, 0.5, 2.5)
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "binomial", weights = w)))
  da <- maihda_discriminatory_accuracy(m)
  prob <- predict_maihda(m, type = "individual", scale = "response")
  expect_equal(da$auc, maihda_auc(prob, d$y))
  # The corrected estimand differs from the old precision-weighted concordance.
  w_fit <- as.numeric(stats::weights(m$model, type = "prior"))
  weighted_auc <- MAIHDA:::maihda_auc_weighted(prob, w_fit * d$y, w_fit)
  expect_false(isTRUE(all.equal(da$auc, weighted_auc)))
  # Reported as unweighted, with the precision-weights flag set.
  expect_false(isTRUE(da$weighted))
  expect_null(da$weight_type)
  expect_true(isTRUE(da$precision_weights_ignored))
  expect_equal(da$n_case, sum(d$y == 1))
  expect_equal(da$n_control, sum(d$y == 0))
  expect_true(is.finite(da$auc) && da$auc >= 0 && da$auc <= 1)
  expect_output(print(da), "precision weights")

  # Uniform precision weights: unweighted and weighted concordances coincide, so the
  # AUC is unchanged; still flagged and reported as unweighted.
  d1 <- d; d1$w <- 1.5
  m1 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d1, family = "binomial", weights = w)))
  da1 <- maihda_discriminatory_accuracy(m1)
  prob1 <- predict_maihda(m1, type = "individual", scale = "response")
  expect_equal(da1$auc, maihda_auc(prob1, d1$y))
  expect_false(isTRUE(da1$weighted))
  expect_true(isTRUE(da1$precision_weights_ignored))

  # Integer frequency weights still take the aggregated path unchanged (they ARE
  # replicated Bernoulli trials, so case/control mass is meaningful there).
  d3 <- d
  d3$w <- sample(1:3, n, replace = TRUE)
  m3 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d3, family = "binomial", weights = w)))
  da3 <- maihda_discriminatory_accuracy(m3)
  expect_equal(da3$n_case + da3$n_control, sum(d3$w))
  expect_true(is.finite(da3$auc) && da3$auc >= 0 && da3$auc <= 1)
  expect_false(isTRUE(da3$weighted))
  expect_false(isTRUE(da3$precision_weights_ignored))
})

test_that("AUC and MOR share the intersectional scope with extra random effects", {
  # Regression: for a site + stratum model the AUC used FULL predictions
  # (including the site effect) while the MOR read only the stratum variance --
  # a strong site effect carried a high AUC (0.90 in the audit repro) over a
  # negligible stratum effect (stratum-only 0.57), so the two summarized
  # different models. The headline AUC is now the intersectional-scope
  # concordance, with the full-model AUC reported alongside.
  skip_on_cran()
  set.seed(505)
  d <- expand.grid(stratum = factor(1:6), site = factor(1:8), rep = 1:30)
  stratum_u <- stats::rnorm(6, sd = 0.3)[d$stratum]
  site_u <- stats::rnorm(8, sd = 2.0)[d$site]
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.2 + stratum_u + site_u))

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | site) + (1 | stratum), data = d, family = "binomial")))
  da <- maihda_discriminatory_accuracy(m)

  # Headline AUC is the intersectional-scope concordance (site effect excluded)...
  eta_strata <- as.numeric(stats::predict(m$model, re.form = ~ (1 | stratum),
                                          type = "link"))
  expect_equal(da$auc, maihda_auc(eta_strata, d$y))
  expect_identical(da$auc_scope, "intersectional")
  # ...with the full-model AUC alongside; the strong site effect makes it larger.
  prob <- predict_maihda(m, type = "individual", scale = "response")
  expect_equal(da$auc_full, maihda_auc(prob, d$y))
  expect_gt(da$auc_full, da$auc)
  # The MOR is the same stratum-scope quantity.
  vc <- lme4::VarCorr(m$model)
  v_str <- as.numeric(as.matrix(vc$stratum)[1, 1])
  expect_equal(da$mor, exp(sqrt(2 * v_str) * stats::qnorm(0.75)))
  expect_output(print(da), "full model")

  # The canonical single-stratum fit is unchanged: one AUC, no full-model twin.
  m0 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "binomial")))
  da0 <- maihda_discriminatory_accuracy(m0)
  expect_identical(da0$auc_scope, "model")
  expect_null(da0$auc_full)
})

test_that("intersectional-scope AUC includes individual-level covariates (not strata-only)", {
  # Audit finding 3: with a non-stratum random effect present, auc_scope switches to
  # the intersectional scope, but maihda_da_scope_scores() retains the WHOLE
  # fixed-effects predictor. So an individual-level covariate in the fixed part (one
  # that varies WITHIN strata) enters this AUC -- it is an ADJUSTED intersectional
  # concordance, not strata-only discrimination, and matches the between-stratum MOR
  # scope only when the fixed part is intercept-only. The label is "intersectional"
  # (not "strata"), and the scope score is not constant within a stratum.
  skip_on_cran()
  set.seed(717)
  n <- 1200
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("A", "B"), n, replace = TRUE),
    site = factor(sample(1:6, n, replace = TRUE)),
    age = stats::rnorm(n))
  sk <- interaction(d$gender, d$race, drop = TRUE)
  site_u <- stats::rnorm(6, sd = 0.5)[d$site]
  strat_u <- stats::rnorm(nlevels(sk), sd = 0.8)[sk]
  d$y <- stats::rbinom(n, 1, stats::plogis(0.9 * d$age + strat_u + site_u))
  d$stratum <- make_strata(d, vars = c("gender", "race"))$data$stratum

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | stratum) + (1 | site), data = d, family = "binomial")))
  da <- maihda_discriminatory_accuracy(m)

  # Labelled the intersectional scope, NOT "strata".
  expect_identical(da$auc_scope, "intersectional")
  expect_false(identical(da$auc_scope, "strata"))

  # The scope score retains the age fixed effect, so it is NOT constant within a
  # stratum -- the decisive test that it is not strata-only discrimination.
  score <- MAIHDA:::maihda_da_scope_scores(
    m, MAIHDA:::maihda_da_re_scopes(m)$intersectional)
  within_sd <- tapply(score, d$stratum, function(z) stats::sd(z))
  expect_gt(max(within_sd, na.rm = TRUE), 1e-6)

  # print() no longer claims strata-only discrimination.
  out <- paste(utils::capture.output(print(da)), collapse = "\n")
  expect_false(grepl("intersectional strata", out, fixed = TRUE))
})

test_that("crossed-dimensions MOR uses the correlated-strata mixture", {
  # Regression, two layers deep. First: maihda_mor() on a crossed-dimensions fit
  # read only the interaction ("stratum") variance, ignoring the additive
  # dimension REs that are part of the between-stratum effect at the intersection
  # level. Then: having summed them, it fed the total to the INDEPENDENT-strata
  # closed form exp(sqrt(2 V) qnorm(.75)) -- but two crossed strata sharing a
  # dimension share that dimension's random effect, so the shared part cancels
  # from their difference instead of contributing 2 tau^2 to it. The MOR now
  # comes from the resulting mixture over strata pairs.
  skip_on_cran()
  set.seed(606)
  n <- 900
  d <- data.frame(g = factor(sample(c("F", "M"), n, replace = TRUE)),
                  r = factor(sample(c("A", "B", "C"), n, replace = TRUE)))
  d$stratum <- interaction(d$g, d$r, sep = "_")
  g_u <- c(F = 0.6, M = -0.6)[as.character(d$g)]
  r_u <- c(A = 0.5, B = 0, C = -0.5)[as.character(d$r)]
  d$y <- stats::rbinom(n, 1, stats::plogis(g_u + r_u))

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | g) + (1 | r) + (1 | stratum), data = d,
               family = "binomial")))
  m$cc_info <- list(dim_groups = c(g = "g", r = "r"),
                    interaction_group = "stratum")

  vars <- MAIHDA:::maihda_random_variances_lme4(m$model)
  got <- maihda_mor(m)

  # All three variances still enter (the first defect): dropping the dimension
  # REs would leave only the interaction, which is ~0 here, hence an MOR of ~1.
  expect_gt(got, 1.2)

  # Reference computed independently by enumerating every unordered pair of the
  # six observed strata: a pair differing on BOTH dimensions has difference
  # variance 2(tau_g^2 + tau_r^2 + tau_I^2), one sharing g only 2(tau_r^2 +
  # tau_I^2), one sharing r only 2(tau_g^2 + tau_I^2).
  cells <- unique(m$data[, c("g", "r")])
  vpair <- unlist(lapply(seq_len(nrow(cells) - 1), function(i) {
    vapply((i + 1):nrow(cells), function(j) {
      2 * (vars[["stratum"]] +
             (cells$g[i] != cells$g[j]) * vars[["g"]] +
             (cells$r[i] != cells$r[j]) * vars[["r"]])
    }, numeric(1))
  }))
  cdf <- function(x) mean(2 * stats::pnorm(x / sqrt(vpair)) - 1)
  ref <- exp(stats::uniroot(function(x) cdf(x) - 0.5,
                            c(1e-10, 10 * sqrt(max(vpair))),
                            tol = .Machine$double.eps^0.5)$root)
  expect_equal(got, ref, tolerance = 1e-6)

  # And it is STRICTLY BELOW the independent-strata closed form, because these
  # strata share dimensions. That inequality is the whole finding.
  expect_lt(got, exp(sqrt(2 * sum(vars[c("g", "r", "stratum")])) *
                       stats::qnorm(0.75)))
  # The dimension REs are intersectional: the DA reports a single-scope AUC.
  da <- maihda_discriminatory_accuracy(m)
  expect_identical(da$auc_scope, "model")
  expect_null(da$auc_full)
})

test_that("maihda_discriminatory_accuracy computes a count-weighted AUC for aggregated binomial", {
  skip_on_cran()
  set.seed(202)
  n <- 1500
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- stats::rbinom(n, 1, stats::plogis(stats::rnorm(nlevels(sk), sd = 0.8)[sk]))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  d <- d[!is.na(d$stratum), , drop = FALSE]

  # Aggregate the Bernoulli rows to per-stratum success/failure counts and fit the
  # SAME logistic MAIHDA as a cbind(success, failure) binomial model.
  agg <- stats::aggregate(y ~ stratum, data = d,
                          FUN = function(z) c(success = sum(z), failure = sum(1 - z)))
  agg <- data.frame(stratum = agg$stratum,
                    success = agg$y[, "success"],
                    failure = agg$y[, "failure"])
  m <- suppressWarnings(suppressMessages(
    fit_maihda(cbind(success, failure) ~ (1 | stratum), data = agg, family = "binomial")
  ))

  da <- maihda_discriminatory_accuracy(m)
  expect_s3_class(da, "maihda_da")
  expect_true(is.finite(da$auc) && da$auc >= 0 && da$auc <= 1)
  # Cases / controls are the TOTAL successes / failures, not the row count.
  expect_equal(da$n_case, sum(agg$success))
  expect_equal(da$n_control, sum(agg$failure))
  # The MOR is still defined for an aggregated logit fit.
  expect_true(is.finite(da$mor) && da$mor >= 1)

  # The reported AUC matches the rank AUC on the implied individual-level 0/1 data.
  prob <- predict_maihda(m, type = "individual", scale = "response")
  expanded_prob <- unlist(Map(function(p, s, f) rep(p, s + f), prob, agg$success, agg$failure))
  expanded_y <- unlist(Map(function(s, f) c(rep(1, s), rep(0, f)), agg$success, agg$failure))
  expect_equal(da$auc, maihda_auc(expanded_prob, expanded_y))
})

test_that("aggregated-binomial AUC stays trial-weighted when every row is all-case or all-control", {
  # Regression: each cbind(success, failure) cell is all-success OR all-failure, so the
  # response PROPORTIONS are exactly 0/1 and the old "response has values outside {0,1}"
  # heuristic misread the fit as individual-level Bernoulli. With several cells per
  # stratum, (1 | stratum) gives ONE fitted probability per stratum, so the unweighted
  # row-level path saw tied probabilities with mixed 0/1 labels -- collapsing the AUC to
  # 0.5 and reporting one pseudo-observation per row. Detection is now structural (the
  # cbind response is a two-column matrix), so the trial-weighted C-statistic and the
  # true success/failure totals are reported.
  agg <- data.frame(
    stratum = factor(c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s2")),
    success = c(90, 0, 70, 0,   0, 8, 0, 6),
    failure = c(0, 12, 0, 10,   85, 0, 60, 0)
  )
  m <- suppressWarnings(suppressMessages(
    fit_maihda(cbind(success, failure) ~ (1 | stratum), data = agg, family = "binomial")
  ))

  da <- maihda_discriminatory_accuracy(m)
  expect_s3_class(da, "maihda_da")
  # Totals are the trial-level success / failure counts, not the 8-row count.
  expect_equal(da$n_case, sum(agg$success))      # 174
  expect_equal(da$n_control, sum(agg$failure))   # 167
  expect_false(da$n_case + da$n_control == nrow(agg))

  # AUC is the trial-weighted C-statistic over the expanded individual-level 0/1 data,
  # NOT the unweighted rank over the stratum rows (which the tied probabilities would
  # have driven to 0.5).
  prob <- predict_maihda(m, type = "individual", scale = "response")
  expanded_prob <- unlist(Map(function(p, s, f) rep(p, s + f), prob, agg$success, agg$failure))
  expanded_y    <- unlist(Map(function(s, f) c(rep(1, s), rep(0, f)), agg$success, agg$failure))
  expect_equal(da$auc, maihda_auc(expanded_prob, expanded_y))
  expect_gt(da$auc, 0.7)                                 # the real separation survives
  expect_gt(da$auc, maihda_auc(prob, as.numeric(lme4::getME(m$model, "y"))))  # > row-level
})

test_that("maihda_discriminatory_accuracy computes a count-weighted AUC for a brms y | trials(n) fit", {
  # The brms analogue of the lme4 cbind(success, failure) test above. brms exposes no
  # weights.brmsfit, so the trial counts come from the trials() addition term parsed off
  # the formula, and the response-scale prediction is the expected COUNT (trials * p),
  # so it is divided by the trial counts to rank by probability. Compiles a Stan model,
  # so OPT-IN like the other brms tests.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  set.seed(404)
  n <- 1200
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- stats::rbinom(n, 1, stats::plogis(stats::rnorm(nlevels(sk), sd = 0.8)[sk]))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  d <- d[!is.na(d$stratum), , drop = FALSE]

  # Aggregate to per-stratum success / trial counts and fit a brms `success | trials(total)`
  # binomial MAIHDA (response stays binomial -- not rewritten to bernoulli).
  agg <- stats::aggregate(y ~ stratum, data = d,
                          FUN = function(z) c(success = sum(z), total = length(z)))
  agg <- data.frame(stratum = agg$stratum,
                    success = agg$y[, "success"],
                    total   = agg$y[, "total"])
  m <- suppressWarnings(suppressMessages(
    fit_maihda(success | trials(total) ~ (1 | stratum), data = agg,
               engine = "brms", family = "binomial",
               chains = 1, iter = 200, refresh = 0)
  ))
  expect_identical(m$family$family, "binomial")  # NOT rewritten to bernoulli

  da <- maihda_discriminatory_accuracy(m)
  expect_s3_class(da, "maihda_da")
  expect_true(is.finite(da$auc) && da$auc >= 0 && da$auc <= 1)
  # Cases / controls are the TOTAL successes / failures, not the row count.
  expect_equal(da$n_case, sum(agg$success))
  expect_equal(da$n_control, sum(agg$total - agg$success))
  expect_false(da$n_case + da$n_control == nrow(agg))
  # The MOR is still defined for an aggregated logit fit.
  expect_true(is.finite(da$mor) && da$mor >= 1)

  # The reported AUC is the count-weighted C-statistic over the implied 0/1 data,
  # ranked by the per-trial probability. predict_maihda(scale = "response") returns
  # that per-trial probability directly for a brms `y | trials(n)` fit -- it
  # normalises brms's expected success COUNT (trials * p) by the trial counts, so it
  # is a probability in [0, 1] like lme4's cbind() fit, NOT an expected count.
  prob <- predict_maihda(m, type = "individual", scale = "response")
  expect_true(all(prob >= 0 & prob <= 1))
  expanded_prob <- unlist(Map(function(p, s, f) rep(p, s + f), prob,
                              agg$success, agg$total - agg$success))
  expanded_y <- unlist(Map(function(s, f) c(rep(1, s), rep(0, f)),
                           agg$success, agg$total - agg$success))
  expect_equal(da$auc, maihda_auc(expanded_prob, expanded_y))
})

test_that("aggregated-binomial AUC folds sampling weights into the case/control mass", {
  # Regression guard: a brms `y | trials(n)` fit with sampling_weights took the
  # aggregated branch, which previously computed an UNWEIGHTED count AUC while still
  # reporting weighted = TRUE (a design-weighted AUC the README/print() claim). The
  # weights must fold into the per-row case/control mass, exactly as the Bernoulli
  # design-weighted branch does, so the reported AUC is the population concordance.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  set.seed(909)
  K <- 12
  agg <- data.frame(
    stratum = factor(sprintf("s%02d", seq_len(K))),
    total   = sample(20:40, K, replace = TRUE)
  )
  agg$success <- stats::rbinom(K, agg$total, stats::plogis(seq(-1.6, 1.6, length.out = K)))
  # Weights negatively associated with the risk gradient, so folding them shifts the
  # ranking mass enough to move the AUC (a no-op weighting would not detect the bug).
  agg$w <- rev(seq(0.2, 3, length.out = K))

  m <- suppressWarnings(suppressMessages(
    fit_maihda(success | trials(total) ~ (1 | stratum), data = agg,
               engine = "brms", family = "binomial", sampling_weights = "w",
               chains = 1, iter = 300, refresh = 0)
  ))

  # predict_maihda(scale = "response") already returns the per-trial probability for
  # a brms `y | trials(n)` fit (it normalises the expected success count internally).
  prob_row <- predict_maihda(m, type = "individual", scale = "response")
  success  <- as.numeric(maihda_da_observed_response(m))
  trials   <- maihda_brms_trial_counts(m)
  sw       <- maihda_prior_weights(m)

  auc_unwt <- maihda_auc_weighted(prob_row, success, trials)
  auc_wt   <- maihda_auc_weighted(prob_row, sw * success, sw * trials)

  da <- maihda_discriminatory_accuracy(m)
  expect_true(da$weighted)                               # claimed design-weighted ...
  expect_equal(da$auc, auc_wt)                           # ... and actually IS weighted
  expect_false(isTRUE(all.equal(da$auc, auc_unwt)))      # not the old unweighted value
  # Reported case/control totals stay UNWEIGHTED observation counts.
  expect_equal(da$n_case, sum(agg$success))
  expect_equal(da$n_control, sum(agg$total - agg$success))
})

test_that("DA helpers accept a brms Bernoulli family (not only 'binomial')", {
  # fit_maihda(engine = "brms") fits a binary 0/1 outcome with bernoulli(); the DA
  # helpers must treat that as a logistic MAIHDA model. Relabel a real lme4 binomial
  # fit's stored family to the bernoulli/logit object brms would carry, so the family
  # gate is exercised without compiling a Stan model.
  m <- maihda_da_test_model()
  m$family <- stats::binomial()      # baseline: a finite MOR / AUC under "binomial"
  ref_mor <- maihda_mor(m)
  ref_da <- maihda_discriminatory_accuracy(m)

  m$family <- structure(list(family = "bernoulli", link = "logit"), class = "family")
  expect_equal(maihda_mor(m), ref_mor)

  da <- maihda_discriminatory_accuracy(m)
  expect_s3_class(da, "maihda_da")
  expect_identical(da$family, "bernoulli")
  expect_equal(da$auc, ref_da$auc)
  expect_equal(da$mor, ref_mor)
})

# ---- summary() integration: DA + response-scale VPC as model summary slots --------

test_that("summary() attaches discriminatory accuracy for a binomial model", {
  m <- maihda_da_test_model()
  s <- summary(m)

  expect_s3_class(s$discriminatory_accuracy, "maihda_da")
  # Identical to calling the helper directly -- summary() just bundles it, no refit.
  expect_equal(s$discriminatory_accuracy$auc, maihda_discriminatory_accuracy(m)$auc)
  expect_equal(s$discriminatory_accuracy$mor, maihda_mor(m))
  # The response-scale VPC is opt-in: off by default.
  expect_null(s$vpc_response)
  expect_output(print(s), "Discriminatory accuracy")
})

test_that("summary(response_vpc = TRUE) attaches a reproducible response-scale VPC", {
  m <- maihda_da_test_model()
  s <- summary(m, response_vpc = TRUE, seed = 1)

  expect_s3_class(s$vpc_response, "maihda_vpc_response")
  expect_true(is.finite(s$vpc_response$estimate))
  # Seeded, so it reproduces the standalone helper exactly.
  expect_equal(s$vpc_response$estimate, maihda_vpc_response(m, seed = 1)$estimate)
  expect_output(print(s), "Response-scale VPC")
})

test_that("summary() attaches no DA / response VPC for a gaussian model", {
  set.seed(11)
  n <- 400
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  d$y <- stats::rnorm(n)
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  g <- fit_maihda(y ~ (1 | stratum), data = d)

  s <- summary(g)
  expect_null(s$discriminatory_accuracy)
  expect_null(s$vpc_response)
  # Even when explicitly requested, response VPC is silently skipped off-family.
  expect_null(summary(g, response_vpc = TRUE, seed = 1)$vpc_response)
})

test_that("maihda() surfaces discriminatory accuracy on the null and adjusted summaries", {
  skip_on_cran()
  set.seed(7)
  n <- 1500
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE),
    age    = stats::rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  lp <- stats::rnorm(nlevels(sk), sd = 0.8)[sk] + 0.3 * d$age
  d$y <- stats::rbinom(n, 1, stats::plogis(lp))

  a <- suppressWarnings(suppressMessages(
    maihda(y ~ age + (1 | gender:race), data = d, family = "binomial")
  ))

  expect_s3_class(a$summary$discriminatory_accuracy, "maihda_da")
  expect_s3_class(a$summary_adjusted$discriminatory_accuracy, "maihda_da")
  # The headline print shows the null-model DA line (with the adjusted AUC alongside).
  expect_output(print(a), "Discriminatory accuracy \\(null model\\)")
})

test_that("maihda(response_vpc = TRUE) attaches the response-scale VPC to the null summary", {
  skip_on_cran()
  set.seed(8)
  n <- 1500
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- stats::rbinom(n, 1, stats::plogis(stats::rnorm(nlevels(sk), sd = 0.8)[sk]))

  a <- suppressWarnings(suppressMessages(
    maihda(y ~ (1 | gender:race), data = d, family = "binomial",
           response_vpc = TRUE, seed = 1)
  ))
  expect_s3_class(a$summary$vpc_response, "maihda_vpc_response")
  expect_true(is.finite(a$summary$vpc_response$estimate))
  # The headline print() renders the response-scale VPC line.
  expect_output(print(a), "Response-scale VPC \\(null\\)")
})

test_that("maihda() headline DA shows MOR = NA for a non-logit (probit) binomial fit", {
  skip_on_cran()
  set.seed(9)
  n <- 1500
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- stats::rbinom(n, 1, stats::pnorm(stats::rnorm(nlevels(sk), sd = 0.5)[sk]))

  a <- suppressWarnings(suppressMessages(
    maihda(y ~ (1 | gender:race), data = d, family = stats::binomial("probit"))
  ))
  da <- a$summary$discriminatory_accuracy
  expect_s3_class(da, "maihda_da")
  expect_true(is.finite(da$auc))   # AUC is link-agnostic
  expect_true(is.na(da$mor))       # MOR undefined for a probit link
  # The headline print renders the NA MOR via the fmt() NA branch.
  expect_output(print(a), "MOR: NA")
})

# Binary stepwise data with a real between-stratum signal (and a gender main effect so
# adding it actually moves the discriminatory accuracy), exposing the data frame.
maihda_da_stepwise_data <- function(seed = 321, n = 900) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race   = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  lp <- stats::rnorm(nlevels(sk), sd = 0.8)[sk] + 0.6 * (d$gender == "M")
  d$y <- stats::rbinom(n, 1, stats::plogis(lp))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  d
}

test_that("stepwise_pcv carries the DA trajectory for a binary outcome", {
  d <- maihda_da_stepwise_data()
  out <- suppressWarnings(suppressMessages(
    stepwise_pcv(d, "y", c("gender", "race"), family = "binomial")
  ))

  expect_s3_class(out, "maihda_stepwise")
  expect_true(all(c("AUC", "Step_AUC", "Total_AUC", "MOR") %in% names(out)))
  expect_equal(nrow(out), 3L)  # null + 2 steps

  # AUC at each step equals maihda_discriminatory_accuracy() on the same fit, on the
  # same (here complete) analytic sample -- no extra fits, just read off each step.
  null_mod <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ 1 + (1 | stratum), data = d, family = "binomial")))
  m1 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ gender + (1 | stratum), data = d, family = "binomial")))
  m2 <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ gender + race + (1 | stratum), data = d, family = "binomial")))
  expect_equal(out$AUC[1], maihda_discriminatory_accuracy(null_mod)$auc)
  expect_equal(out$AUC[2], maihda_discriminatory_accuracy(m1)$auc)
  expect_equal(out$AUC[3], maihda_discriminatory_accuracy(m2)$auc)

  # Step_AUC / Total_AUC are ABSOLUTE deltas; the null row anchors both at 0.
  expect_equal(out$Step_AUC, c(0, diff(out$AUC)))
  expect_equal(out$Total_AUC, out$AUC - out$AUC[1])
  expect_equal(out$Step_AUC[1], 0)
  expect_equal(out$Total_AUC[1], 0)

  # MOR per step is exp(sqrt(2 * V_A) * qnorm(0.75)) on the between-stratum variance
  # already in the Variance column (logit link).
  expect_equal(out$MOR, exp(sqrt(2 * out$Variance) * stats::qnorm(0.75)))

  # print() surfaces the proportional-PCV vs absolute-delta-AUC legend.
  expect_output(print(out), "absolute changes in AUC")
})
