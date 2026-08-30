# Audit 2026-08-30 (fourth pass of the day): a brms aggregated-binomial outcome
# written as `success | trials(n)` keeps its denominator in the trials() addition
# term, not in the response. maihda_describe() extracted only the leftmost leaf of
# the `|` expression (maihda_describe_response_expr) and plot_obs_vs_shrunken()
# only stats::model.response()/all.vars()[1], so both handed the observed-outcome
# extractor a bare vector of SUCCESS counts, which defaulted the denominator to 1.
# A four-row 26-of-62 sample was reported as "26 events / 4 trials (650.0%)", and
# the observed-vs-shrunken plot put mean success COUNTS on the x-axis against
# trial-weighted predicted probabilities on the y-axis. The lme4
# cbind(success, failure) twin was always correct -- its denominator is structural,
# in the matrix response itself.
#
# Scope correction to the report: the existing trials parser was not merely
# "unused". brms also allows a CONSTANT trial count (`y | trials(20)`, which
# make_standata() recycles to one value per row), and maihda_brms_trial_counts()
# returned that length-1 vector unrecycled, so every length-guarded caller treated
# a constant-trials fit as having no trial counts at all: unit prediction weights,
# predict_maihda(scale = "response") returning an expected COUNT where it promises
# a probability, and maihda_discriminatory_accuracy() falling through to the
# observation-level path. Recycling now happens in maihda_trials_from_formula(),
# the shared formula-level extractor both the description and the plots use.

# Four rows, 26 successes out of 62 trials -- the report's reproduction.
audit_0830d_data <- function() {
  data.frame(
    y      = c(3, 8, 5, 10),
    n      = c(10, 20, 12, 20),
    gender = c("F", "F", "M", "M"),
    race   = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
}

# A maihda_model shell around a mock brmsfit: the observed-outcome path reads only
# $formula / $data / $family, so it needs no Stan.
audit_0830d_model <- function(d, formula) {
  structure(
    list(
      model = structure(list(data = d), class = "brmsfit"),
      engine = "brms",
      formula = formula,
      data = d,
      original_data = d,
      family = list(family = "binomial", link = "logit"),
      strata_vars = c("gender", "race"),
      strata_sep = " x ",
      strata_info = NULL,
      strata_autobin_info = NULL,
      context_vars = NULL,
      sampling_weights = NULL,
      longitudinal_info = NULL,
      response_recoding = NULL
    ),
    class = "maihda_model"
  )
}

test_that("maihda_trials_from_formula reads a trials() addition term off any formula", {
  d <- data.frame(y = 1:4, n = c(2, 4, 6, 8), w = 1)

  expect_equal(MAIHDA:::maihda_trials_from_formula(y | trials(n) ~ x, d), c(2, 4, 6, 8))

  # brms allows a CONSTANT trial count and recycles it to one value per row; so
  # must this, or every length-guarded caller sees "no trial counts".
  expect_equal(MAIHDA:::maihda_trials_from_formula(y | trials(5) ~ x, d), rep(5, 4))

  # Combined addition terms, in either order.
  expect_equal(MAIHDA:::maihda_trials_from_formula(y | trials(n) + weights(w) ~ x, d),
               c(2, 4, 6, 8))
  expect_equal(MAIHDA:::maihda_trials_from_formula(y | weights(w) + trials(n) ~ x, d),
               c(2, 4, 6, 8))

  # No trials() term anywhere: NULL, so every non-aggregated response is untouched.
  expect_null(MAIHDA:::maihda_trials_from_formula(y | weights(w) ~ x, d))
  expect_null(MAIHDA:::maihda_trials_from_formula(y ~ x, d))
  expect_null(MAIHDA:::maihda_trials_from_formula(~ x, d))

  # An unevaluable trials() term returns NULL rather than reaching outside the
  # data for a stray same-named object -- `n` is a common variable name, and
  # enclos = baseenv() keeps the search off the caller's frame.
  no_such_col <- c(99, 99, 99, 99)
  expect_null(MAIHDA:::maihda_trials_from_formula(y | trials(no_such_col) ~ x, d))

  # A length mismatch that is not a scalar is unusable, not silently recycled.
  expect_null(MAIHDA:::maihda_trials_from_formula(y | trials(n) ~ x, d, n = 2))

  # A non-numeric trials column is not coerced through factor level codes.
  d_f <- d
  d_f$n <- factor(d$n)
  expect_null(MAIHDA:::maihda_trials_from_formula(y | trials(n) ~ x, d_f))
})

test_that("maihda_describe reports successes/trials for a brms `y | trials(n)` outcome", {
  d <- audit_0830d_data()
  m <- audit_0830d_model(d, y | trials(n) ~ (1 | stratum))
  m$data$stratum <- m$original_data$stratum <-
    factor(c("F x A", "F x B", "M x A", "M x B"))

  desc <- maihda_describe(m)
  oo <- desc$outcome_overall

  # Was: 26 events / 4 trials, proportion 6.5.
  expect_equal(oo$outcome_events, 26)
  expect_equal(oo$outcome_trials, 62)
  expect_equal(oo$outcome_proportion, 26 / 62)
  expect_true(oo$outcome_proportion >= 0 && oo$outcome_proportion <= 1)

  # Per stratum, and in the slim $observations frame the describe plots read.
  expect_equal(desc$strata$outcome_trials, c(10, 20, 12, 20))
  expect_equal(desc$strata$outcome_proportion, c(3, 8, 5, 10) / c(10, 20, 12, 20))
  expect_equal(desc$observations$denominator, c(10, 20, 12, 20))

  # The printed line quotes the trials, not the row count.
  out <- paste(utils::capture.output(print(desc)), collapse = "\n")
  expect_match(out, "26 events / 62 trials")
  expect_false(grepl("650.0%", out, fixed = TRUE))
})

test_that("maihda_describe handles a constant trials() and the formula entry point", {
  d <- audit_0830d_data()

  # Formula path (no fitted model at all): the same extraction must apply, since
  # this is the pre-fit description of data destined for a brms trials() fit.
  desc <- maihda_describe(y | trials(n) ~ (1 | gender:race), data = d,
                          family = "binomial", flag_stratum_n = 0)
  expect_equal(desc$outcome_overall$outcome_events, 26)
  expect_equal(desc$outcome_overall$outcome_trials, 62)
  expect_equal(desc$outcome_overall$outcome_proportion, 26 / 62)

  # Constant trial count: 4 rows x 20 trials.
  desc_s <- maihda_describe(y | trials(20) ~ (1 | gender:race), data = d,
                            family = "binomial", flag_stratum_n = 0)
  expect_equal(desc_s$outcome_overall$outcome_trials, 80)
  expect_equal(desc_s$outcome_overall$outcome_proportion, 26 / 80)

  # A trials() term under a non-binomial family is not a denominator and is
  # dropped: the outcome stays a plain continuous summary.
  d$z <- c(1.5, 2.5, 3.5, 4.5)
  desc_g <- maihda_describe(z | trials(n) ~ (1 | gender:race), data = d,
                            family = "gaussian", flag_stratum_n = 0)
  expect_equal(desc_g$outcome_overall$outcome_mean, mean(d$z))
})

test_that("rows with unusable trial counts are excluded, and malformed rows warn", {
  d <- audit_0830d_data()

  # A missing or zero trial count is no observed outcome -- the rule the matrix
  # (cbind) branch already applied to a non-finite or zero row total.
  d_na <- d
  d_na$n[2] <- NA
  d_na$n[3] <- 0
  desc <- maihda_describe(y | trials(n) ~ (1 | gender:race), data = d_na,
                          family = "binomial", flag_stratum_n = 0)
  expect_equal(desc$overview$n_missing_outcome, 2)
  expect_equal(desc$outcome_overall$outcome_events, 3 + 10)
  expect_equal(desc$outcome_overall$outcome_trials, 10 + 20)

  # brms refuses to fit successes outside [0, trials] ("Number of trials is
  # smaller than the number of events"), so only the pre-fit description can see
  # such a row. Warn and drop it rather than report a proportion above 1.
  d_bad <- d
  d_bad$n[1] <- 2      # 3 successes out of 2 trials
  expect_warning(
    desc_bad <- maihda_describe(y | trials(n) ~ (1 | gender:race), data = d_bad,
                                family = "binomial", flag_stratum_n = 0),
    "successes outside \\[0, trials\\]")
  expect_equal(desc_bad$outcome_overall$outcome_events, 8 + 5 + 10)
  expect_equal(desc_bad$outcome_overall$outcome_trials, 20 + 12 + 20)
  expect_true(desc_bad$outcome_overall$outcome_proportion <= 1)

  # A dropped row must also be MISSING, or the report does not reconcile: it
  # would read "23 events / 52 trials" over 4 non-missing rows and 0 missing
  # outcomes, while only 3 rows contributed.
  expect_equal(desc_bad$overview$n_missing_outcome, 1)
  expect_equal(desc_bad$outcome_overall$n_nonmissing, 3)
  expect_equal(sum(desc_bad$strata$n_missing_outcome), 1)
  expect_equal(sum(desc_bad$strata$outcome_events),
               desc_bad$outcome_overall$outcome_events)
  expect_equal(sum(desc_bad$strata$outcome_trials),
               desc_bad$outcome_overall$outcome_trials)
  out_bad <- paste(utils::capture.output(print(desc_bad)), collapse = "\n")
  expect_match(out_bad, "3 non-missing")

  # The well-formed case is untouched: no row is spuriously marked missing.
  desc_ok <- maihda_describe(y | trials(n) ~ (1 | gender:race), data = d,
                             family = "binomial", flag_stratum_n = 0)
  expect_equal(desc_ok$overview$n_missing_outcome, 0)
  expect_equal(desc_ok$outcome_overall$n_nonmissing, 4)

  # Negative successes are malformed the same way.
  d_neg <- d
  d_neg$y[1] <- -1
  expect_warning(
    maihda_describe(y | trials(n) ~ (1 | gender:race), data = d_neg,
                    family = "binomial", flag_stratum_n = 0),
    "successes outside \\[0, trials\\]")
})

test_that("the observed-outcome extractor puts a trials() outcome on the proportion scale", {
  d <- audit_0830d_data()

  # Successes with the denominator supplied separately == the cbind() matrix form.
  od <- MAIHDA:::maihda_observed_outcome_for_plot(
    d$y, list(family = "binomial", link = "logit"), trials = d$n)
  od_cbind <- MAIHDA:::maihda_observed_outcome_for_plot(
    cbind(d$y, d$n - d$y), list(family = "binomial", link = "logit"))
  expect_equal(od$numerator, od_cbind$numerator)
  expect_equal(od$denominator, od_cbind$denominator)
  expect_equal(MAIHDA:::maihda_observed_weighted_mean(od$numerator, od$denominator),
               26 / 62)
  expect_equal(MAIHDA:::maihda_observed_sample_size(od$numerator, od$denominator), 62)

  # Denominator 1 without trials is the CORRECT contract for this helper -- a
  # numeric response with no separate denominator is one observation per row (a
  # Bernoulli 0/1 outcome). The defect was never here: it was the two call sites
  # not passing the trials they had. Pinned so the argument stays load-bearing
  # and the fix cannot be "moved down" into the extractor's default.
  od_none <- MAIHDA:::maihda_observed_outcome_for_plot(
    d$y, list(family = "binomial", link = "logit"))
  expect_equal(od_none$denominator, rep(1, 4))

  # Trials are ignored under a non-binomial family (they are not a denominator
  # there), and a non-numeric response with trials is rejected outright.
  od_g <- MAIHDA:::maihda_observed_outcome_for_plot(
    d$y, list(family = "gaussian", link = "identity"), trials = d$n)
  expect_equal(od_g$denominator, rep(1, 4))
  expect_error(
    MAIHDA:::maihda_observed_outcome_for_plot(
      factor(letters[1:4]), list(family = "binomial", link = "logit"), trials = d$n),
    "numeric success counts")
  expect_error(
    MAIHDA:::maihda_observed_outcome_for_plot(
      d$y, list(family = "binomial", link = "logit"), trials = d$n[1:2]),
    "one trial count per observation")
})

test_that("plot_obs_vs_shrunken sends a trials() outcome to the extractor (Stan-free)", {
  # The plot's own y-axis needs a real brms fit, but the defect is entirely on the
  # x-axis: the observed response and the trial counts the plot hands the
  # extractor. Exercise exactly that pairing off a mock fit, so the guard runs on
  # every platform.
  d <- audit_0830d_data()
  d$stratum <- factor(c("F x A", "F x B", "M x A", "M x B"))
  m <- audit_0830d_model(d, y | trials(n) ~ (1 | stratum))

  resp <- MAIHDA:::maihda_observed_response_from_model_frame(m$data, m$formula)
  trials <- MAIHDA:::maihda_trials_from_formula(m$formula, m$data)
  expect_equal(trials, d$n)

  od <- MAIHDA:::maihda_observed_outcome_for_plot(resp, m$family, trials = trials)
  obs <- vapply(seq_len(4), function(i)
    MAIHDA:::maihda_observed_weighted_mean(od$numerator[i], od$denominator[i]),
    numeric(1))
  # Every observed value is a proportion in [0, 1], on the same scale as the
  # predicted per-trial probabilities the y-axis carries -- not a success count.
  expect_equal(obs, d$y / d$n)
  expect_true(all(obs >= 0 & obs <= 1))
  expect_equal(MAIHDA:::maihda_observed_sample_size(od$numerator, od$denominator), 62)
})

test_that("a constant trials() fit is trial-weighted everywhere (Stan-free)", {
  d <- audit_0830d_data()
  d$stratum <- factor(c("s1", "s1", "s2", "s2"))
  m <- audit_0830d_model(d, y | trials(20) ~ (1 | stratum))

  # Was: a bare 20 of length 1, which every length guard below rejected.
  expect_equal(MAIHDA:::maihda_brms_trial_counts(m), rep(20, 4))
  expect_equal(MAIHDA:::maihda_prediction_weights(m), rep(20, 4))

  # The response-scale normalisation now applies: brms returns the expected COUNT
  # (trials * p), so 20 * 0.3 must come back as the probability 0.3.
  expect_equal(MAIHDA:::maihda_brms_response_to_prob(m, rep(6, 4), d), rep(0.3, 4))

  # And the discriminatory-accuracy path sees real aggregated counts.
  counts <- MAIHDA:::maihda_da_brms_aggregated_counts(m)
  expect_equal(counts$successes, d$y)
  expect_equal(counts$trials, rep(20, 4))

  # A vector trials() fit is unchanged, and a fit with no trials() term still has
  # no trial counts (never spuriously trial-weighted).
  m_v <- audit_0830d_model(d, y | trials(n) ~ (1 | stratum))
  expect_equal(MAIHDA:::maihda_prediction_weights(m_v), d$n)
  m_none <- audit_0830d_model(d, y ~ (1 | stratum))
  expect_null(MAIHDA:::maihda_brms_trial_counts(m_none))
  expect_equal(MAIHDA:::maihda_prediction_weights(m_none), rep(1, 4))
})

test_that("non-aggregated outcomes are untouched by the trials threading", {
  # The regression guard for the fix itself: every response that carries no
  # trials() term must summarise exactly as before.
  set.seed(830)
  n <- 80L
  d <- data.frame(
    g = rep(c("F", "M"), n / 2),
    r = rep(c("A", "B"), each = n / 2),
    w = stats::runif(n, 1, 3)
  )
  d$bin <- stats::rbinom(n, 1, 0.4)
  d$cont <- stats::rnorm(n)
  d$cnt <- stats::rpois(n, 2)

  bern <- maihda_describe(bin ~ (1 | g:r), data = d, flag_stratum_n = 0)
  expect_equal(bern$outcome_overall$outcome_proportion, mean(d$bin))
  expect_equal(bern$outcome_overall$outcome_trials, n)

  gauss <- maihda_describe(cont ~ (1 | g:r), data = d, flag_stratum_n = 0)
  expect_equal(gauss$outcome_overall$outcome_mean, mean(d$cont))

  cnt <- maihda_describe(cnt ~ (1 | g:r), data = d, family = "poisson",
                         flag_stratum_n = 0)
  expect_equal(cnt$outcome_overall$outcome_mean, mean(d$cnt))

  # A `weights()`-only addition term has no denominator to recover.
  wonly <- maihda_describe(bin | weights(w) ~ (1 | g:r), data = d,
                           family = "binomial", flag_stratum_n = 0)
  expect_equal(wonly$outcome_overall$outcome_proportion, mean(d$bin))

  # The lme4 cbind() form -- always correct -- and the trials form agree on every
  # summarised quantity. They are not identical objects: the outcome's NAME
  # legitimately differs ("cbind(succ, fail)" vs "succ"), which shows up in
  # $outcome and in $missingness$variable. The numbers are what must match.
  agg <- data.frame(succ = c(3, 8, 5, 10), tot = c(10, 20, 12, 20),
                    g = c("F", "F", "M", "M"), r = c("A", "B", "A", "B"))
  agg$fail <- agg$tot - agg$succ
  cb <- maihda_describe(cbind(succ, fail) ~ (1 | g:r), data = agg,
                        family = "binomial", flag_stratum_n = 0)
  tr <- maihda_describe(succ | trials(tot) ~ (1 | g:r), data = agg,
                        family = "binomial", flag_stratum_n = 0)
  expect_equal(cb$outcome_overall$outcome_events, tr$outcome_overall$outcome_events)
  expect_equal(cb$outcome_overall$outcome_trials, tr$outcome_overall$outcome_trials)
  expect_equal(cb$strata$outcome_proportion, tr$strata$outcome_proportion)
  expect_equal(cb$observations$denominator, tr$observations$denominator)

  # Sharper than the four columns above: enumerate EVERY component that differs,
  # so a future divergence anywhere in the object has to be declared here rather
  # than slipping past a hand-picked list. Only the outcome's name and the
  # recorded call may differ; $missingness differs solely in its `variable` label.
  differing <- names(cb)[!vapply(names(cb),
                                 function(k) isTRUE(all.equal(cb[[k]], tr[[k]])),
                                 logical(1))]
  expect_setequal(differing, c("outcome", "call", "missingness"))
  expect_equal(cb$missingness[, c("role", "n_missing", "pct_missing")],
               tr$missingness[, c("role", "n_missing", "pct_missing")])
  expect_equal(cb$missingness$variable[-1], tr$missingness$variable[-1])
})
