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
  # ranked by the per-trial probability: brms's response prediction is the expected
  # COUNT (trials * p), so it is divided by the trial counts first.
  pred_count <- predict_maihda(m, type = "individual", scale = "response")
  prob <- pred_count / agg$total
  expanded_prob <- unlist(Map(function(p, s, f) rep(p, s + f), prob,
                              agg$success, agg$total - agg$success))
  expanded_y <- unlist(Map(function(s, f) c(rep(1, s), rep(0, f)),
                           agg$success, agg$total - agg$success))
  expect_equal(da$auc, maihda_auc(expanded_prob, expanded_y))
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
