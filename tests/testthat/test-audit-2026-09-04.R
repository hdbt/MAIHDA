# Regression tests for the 2026-09-01c statistical audit pass.
#
# FINDING: maihda_describe() under-reported an aggregated-binomial outcome whenever
# the trial counts were supplied the way ?glm documents them -- a PROPORTION response
# with the denominators passed as `weights =` -- instead of as
# cbind(successes, failures).
#
# The two spellings are the same model (identical fixed effects to ~1e-10), but only
# the cbind() form puts the denominator where describe could see it. describe reads
# the outcome off the formula, and a proportion response leaves no trace of its
# trials there, so the denominator was taken to be 1 per ROW: a 12-row, 340-of-617
# sample was reported as "6.420795 events / 12 trials (53.5%)" -- a FRACTIONAL event
# count over a row count -- where the cbind() spelling of the same data reported
# "340 events / 617 trials (55.1%)". maihda_discriminatory_accuracy() read that same
# fit correctly as 340 cases / 277 controls, so describe and the AUC reported two
# different sample sizes for one model.
#
# The fix does not re-implement the detection rule: maihda_describe_model_agg() calls
# maihda_da_aggregated_counts(), the rule the AUC already uses, so the two cannot
# drift apart again. Deferring rather than duplicating is the whole design, and it is
# what makes this pass compose with the FIFTY-FOURTH pass (uncommitted in the main
# checkout when this was written), which rewrites that same rule to key on the
# WEIGHTS rather than the response, so that a 0/1 response with integral weights
# becomes trial counts too. describe follows the rule wherever it goes: these tests
# therefore pin the AGREEMENT between describe and the AUC, and pin exact counts only
# for the proportion-response spelling, whose reading both versions of the rule share.
# Anything asserting the 0/1 spelling's specific totals would pin a decision the user
# has already reversed.
#
# MERGE NOTE (this pass and the 54th both edit R/discriminatory_accuracy.R):
#   1. This pass adds a `context` argument to maihda_da_aggregated_counts() and
#      maihda_da_proportion_successes(), so the shared fractional-successes warning is
#      not worded "AUC" when maihda_describe() is what raised it. Keep it alongside
#      the 54th's `binomial_weights`: maihda_da_aggregated_counts(model,
#      binomial_weights = ..., context = ...).
#   2. This pass also FACTORS the (y, w) rule out of that wrapper into
#      maihda_agg_counts_from_weights(y, w, context), which the pre-fit describe path
#      calls with no model in hand. The 54th's rule change -- dropping the
#      any(y > 0 & y < 1) requirement so integral weights are trial counts whatever
#      the response -- therefore belongs INSIDE that helper now, not in the wrapper it
#      was originally written against. Applied there, both describe paths inherit it
#      for free and these tests still pass.
#
# PART 2 of the same finding: the pre-fit maihda_describe(formula, data) path had no
# `weights` argument at all, so it could not see the trial counts even in principle.
# It now takes one, with fit_maihda()'s non-standard-evaluation idiom, so the same
# call describes and fits the same sample -- and, being precision weights, they also
# exclude their own invalid rows from the analytic sample. It is the LAST formal
# rather than sitting beside `sampling_weights`, so no positional call written against
# 0.2.1 shifts.

# The audit's fixture: one aggregated table, three spellings, 340 of 617 events.
agg_fixture <- function() {
  a <- data.frame(
    stratum = factor(rep(1:6, each = 2)),
    x    = rep(c(0, 1), 6),
    succ = c(28, 41, 12, 33, 9, 26, 35, 47, 18, 31, 22, 38),
    tot  = c(50, 60, 40, 55, 45, 52, 48, 62, 44, 58, 47, 56)
  )
  a$fail <- a$tot - a$succ
  a$p <- a$succ / a$tot
  a
}

test_that("a proportion response with `weights =` trials is described as its trials", {
  skip_on_cran()
  a <- agg_fixture()
  expect_identical(sum(a$succ), 340)
  expect_identical(sum(a$tot), 617)
  # The rule's discriminator: a genuine (strictly fractional) proportion response.
  expect_true(any(a$p > 0 & a$p < 1))

  m_pr <- suppressWarnings(suppressMessages(fit_maihda(
    p ~ x + (1 | stratum), data = a, family = "binomial", weights = tot)))
  m_cb <- suppressMessages(fit_maihda(
    cbind(succ, fail) ~ x + (1 | stratum), data = a, family = "binomial"))
  # Same model, so any disagreement below is a reporting bug and nothing else.
  expect_equal(lme4::fixef(m_pr$model), lme4::fixef(m_cb$model), tolerance = 1e-6)

  d_pr <- maihda_describe(m_pr)
  d_cb <- maihda_describe(m_cb)

  # The finding: this was 6.420795 events / 12 trials.
  expect_equal(d_pr$outcome_overall$outcome_events, 340)
  expect_equal(d_pr$outcome_overall$outcome_trials, 617)
  expect_equal(d_pr$outcome_overall$outcome_proportion, 340 / 617)
  # An event count is a COUNT; the old value was not even a whole number.
  expect_equal(d_pr$outcome_overall$outcome_events,
               round(d_pr$outcome_overall$outcome_events))

  # Pinned against the cbind() spelling of identical data, not just against 340/617.
  expect_equal(d_pr$outcome_overall$outcome_events, d_cb$outcome_overall$outcome_events)
  expect_equal(d_pr$outcome_overall$outcome_trials, d_cb$outcome_overall$outcome_trials)
  expect_equal(d_pr$outcome_overall$outcome_proportion,
               d_cb$outcome_overall$outcome_proportion)
  # Per stratum too -- the overall total can agree while the strata do not.
  for (cl in c("outcome_events", "outcome_trials", "outcome_proportion")) {
    expect_equal(d_pr$strata[[cl]], d_cb$strata[[cl]], info = cl)
  }
  expect_equal(sum(d_pr$strata$outcome_events), 340)
  expect_equal(sum(d_pr$strata$outcome_trials), 617)

  # describe and the AUC now agree about the size of the sample.
  da <- maihda_discriminatory_accuracy(m_pr)
  expect_equal(da$n_case, d_pr$outcome_overall$outcome_events)
  expect_equal(da$n_case + da$n_control, d_pr$outcome_overall$outcome_trials)
  expect_output(print(d_pr), "340 events / 617 trials")
})

test_that("describe and the AUC agree on the frequency-cell spelling too", {
  skip_on_cran()
  a <- agg_fixture()
  # The frequency-weighted spelling of the SAME table: a 0/1 response whose counts
  # ride in `weights =`. Whether this reads as an aggregated binomial (340 of 617) or
  # as a precision-weighted Bernoulli (12 of 24) is a decision of the SHARED rule, and
  # the 54th pass moves it from the second reading to the first. This test therefore
  # asserts only what the fix itself owns and what holds under both readings: describe
  # and the AUC must give ONE answer for one model. Before the fix they could not --
  # describe had no way to see the weights at all.
  cell <- rbind(
    data.frame(stratum = a$stratum, x = a$x, y = 1, n = a$succ),
    data.frame(stratum = a$stratum, x = a$x, y = 0, n = a$fail))
  m <- suppressWarnings(suppressMessages(fit_maihda(
    y ~ x + (1 | stratum), data = cell, family = "binomial", weights = n)))

  d <- maihda_describe(m)
  da <- suppressWarnings(maihda_discriminatory_accuracy(m))
  expect_equal(da$n_case, d$outcome_overall$outcome_events)
  expect_equal(da$n_case + da$n_control, d$outcome_overall$outcome_trials)

  # And whichever reading is in force, it is one of exactly these two -- not a third
  # number arrived at by describe and the AUC disagreeing about the rows.
  expect_true(
    isTRUE(all.equal(c(d$outcome_overall$outcome_events,
                       d$outcome_overall$outcome_trials), c(340, 617))) ||
    isTRUE(all.equal(c(d$outcome_overall$outcome_events,
                       d$outcome_overall$outcome_trials), c(12, 24))))

  # describe tracks the rule rather than carrying its own copy: whatever the shared
  # detector says about this fit, describe's trial counts say the same.
  agg <- MAIHDA:::maihda_da_aggregated_counts(m)
  ours <- MAIHDA:::maihda_describe_model_agg(m, m$original_data)
  expect_identical(is.null(agg), is.null(ours))
  if (!is.null(agg)) {
    expect_equal(sum(agg$trials), sum(ours$trials, na.rm = TRUE))
    expect_equal(sum(agg$successes), sum(ours$successes, na.rm = TRUE))
  }
})

test_that("the weights-derived trials are aligned by row, not by position", {
  skip_on_cran()
  a <- agg_fixture()
  a$x[3] <- NA_real_          # lme4 drops row 3 for a missing COVARIATE
  m_pr <- suppressWarnings(suppressMessages(fit_maihda(
    p ~ x + (1 | stratum), data = a, family = "binomial", weights = tot)))
  # The hazard: the counts live on the analytic frame (11 rows), while describe
  # evaluates the outcome on the pre-fit frame (12 rows). Taking them as positionally
  # aligned would shift every count after row 3 into the wrong stratum.
  expect_identical(nrow(m_pr$original_data), 12L)
  expect_identical(nrow(m_pr$data), 11L)

  d <- maihda_describe(m_pr)
  expect_equal(d$outcome_overall$outcome_events, sum(a$succ[-3]))
  expect_equal(d$outcome_overall$outcome_trials, sum(a$tot[-3]))
  expect_equal(d$outcome_overall$n_nonmissing, 11L)

  # Row 3 is in stratum 2, and it is stratum 2 -- no other -- that loses it.
  expect_identical(as.character(a$stratum[3]), "2")
  expect_identical(d$strata$n_analytic, c(2L, 1L, 2L, 2L, 2L, 2L))
  expect_identical(d$strata$n_missing_outcome, c(0L, 1L, 0L, 0L, 0L, 0L))
  keep <- seq_len(12)[-3]
  expect_equal(d$strata$outcome_events,
               as.numeric(tapply(a$succ[keep], a$stratum[keep], sum)))
  expect_equal(d$strata$outcome_trials,
               as.numeric(tapply(a$tot[keep], a$stratum[keep], sum)))
})

test_that("the detector does not fire on models that are not aggregated binomials", {
  skip_on_cran()
  set.seed(3)
  st <- factor(rep(seq_len(8), each = 25))
  u <- stats::rnorm(8, sd = 0.8)[st]
  d <- data.frame(stratum = st, z = stats::rnorm(200),
                  y = stats::rbinom(200, 1, stats::plogis(u)))

  # (a) a plain Bernoulli fit: unchanged, and no trials anywhere.
  m0 <- suppressMessages(fit_maihda(
    y ~ z + (1 | stratum), data = d, family = "binomial"))
  expect_null(MAIHDA:::maihda_describe_model_agg(m0, m0$original_data))
  expect_equal(maihda_describe(m0)$outcome_overall$outcome_events, sum(d$y))
  expect_equal(maihda_describe(m0)$outcome_overall$outcome_trials, 200)

  # (b) a Gaussian fit whose outcome happens to lie in [0, 1], carrying integral
  # precision weights. It matches the aggregated detector's SHAPE, so the family
  # must be checked before the detector runs -- otherwise a Gaussian model warns
  # about "successes" it does not have. The reported summary stays mean/SD.
  g <- data.frame(stratum = st, y = stats::runif(200),
                  w = sample(1:5, 200, replace = TRUE))
  mg <- suppressMessages(fit_maihda(
    y ~ (1 | stratum), data = g, family = "gaussian", weights = w))
  expect_null(MAIHDA:::maihda_describe_model_agg(mg, mg$original_data))
  # expect_silent is the actual regression guard here: gating on the family only
  # AFTER running the detector still returns the right summary, and would differ
  # from this line alone -- the detector would already have warned.
  expect_silent(dg <- maihda_describe(mg))
  expect_true("outcome_mean" %in% names(dg$outcome_overall))
  expect_false("outcome_trials" %in% names(dg$outcome_overall))
  expect_equal(dg$outcome_overall$outcome_mean, mean(g$y))
})

# ---- part 2: the pre-fit formula path ---------------------------------------

# The fixture again, with dimension columns for the (1 | g:r) shorthand and no
# pre-built `stratum` column for make_strata() to collide with.
agg_fixture_dims <- function() {
  a <- agg_fixture()
  a$stratum <- NULL
  a$g <- factor(rep(c("M", "F"), 6))
  a$r <- factor(rep(c("A", "B", "C"), 4))
  a
}

test_that("maihda_describe(weights = ) reads pre-fit trial counts", {
  skip_on_cran()
  a <- agg_fixture_dims()

  d_w <- maihda_describe(p ~ x + (1 | g:r), data = a, family = "binomial",
                         weights = tot)
  d_cb <- maihda_describe(cbind(succ, fail) ~ x + (1 | g:r), data = a,
                          family = "binomial")

  expect_equal(d_w$outcome_overall$outcome_events, 340)
  expect_equal(d_w$outcome_overall$outcome_trials, 617)
  # Pinned against the cbind() spelling of identical data, per stratum as well.
  for (cl in c("outcome_events", "outcome_trials", "outcome_proportion")) {
    expect_equal(d_w$strata[[cl]], d_cb$strata[[cl]], info = cl)
  }

  # The whole point of matching fit_maihda()'s idiom: the same call, described
  # before the fit and after it, gives the same numbers.
  m <- suppressWarnings(suppressMessages(fit_maihda(
    p ~ x + (1 | g:r), data = a, family = "binomial", weights = tot)))
  d_fit <- maihda_describe(m)
  expect_equal(d_w$outcome_overall$outcome_events, d_fit$outcome_overall$outcome_events)
  expect_equal(d_w$outcome_overall$outcome_trials, d_fit$outcome_overall$outcome_trials)
})

test_that("pre-fit precision weights exclude their own invalid rows", {
  skip_on_cran()
  a <- agg_fixture_dims()
  a$tot0 <- a$tot
  a$tot0[3] <- 0          # lme4 drops a zero-weight row; so must the description
  d <- suppressWarnings(maihda_describe(
    p ~ x + (1 | g:r), data = a, family = "binomial", weights = tot0))
  expect_equal(d$outcome_overall$outcome_events, sum(a$succ[-3]))
  expect_equal(d$outcome_overall$outcome_trials, sum(a$tot[-3]))
  expect_equal(d$overview$n_analytic, 11)
})

test_that("`weights` is rejected where it cannot mean anything", {
  skip_on_cran()
  a <- agg_fixture_dims()
  m <- suppressWarnings(suppressMessages(fit_maihda(
    p ~ x + (1 | g:r), data = a, family = "binomial", weights = tot)))
  # A fitted model already carries its weights.
  expect_error(maihda_describe(m, weights = tot), "already carries")
  # Design weights and precision weights mean different things.
  expect_error(
    maihda_describe(p ~ x + (1 | g:r), data = a, family = "binomial",
                    weights = tot, sampling_weights = "tot"),
    "not both")
  expect_error(
    maihda_describe(p ~ x + (1 | g:r), data = a, family = "binomial",
                    weights = no_such_column),
    "could not evaluate 'weights'")
  expect_error(
    maihda_describe(p ~ x + (1 | g:r), data = a, family = "binomial",
                    weights = c(1, 2)),
    "one value per row")
})

test_that("pre-fit `weights` does not turn other families into events/trials", {
  skip_on_cran()
  set.seed(11)
  d <- data.frame(a1 = factor(rep(c("p", "q"), 30)),
                  a2 = factor(rep(c("u", "v", "w"), 20)),
                  w = sample(1:5, 60, replace = TRUE))
  # (a) a Gaussian outcome that happens to lie in [0, 1]: mean/SD, and silent.
  d$y <- stats::runif(60)
  expect_silent(
    dg <- maihda_describe(y ~ (1 | a1:a2), data = d, family = "gaussian", weights = w))
  expect_true("outcome_mean" %in% names(dg$outcome_overall))
  expect_false("outcome_trials" %in% names(dg$outcome_overall))

  # (b) a 0/1 binomial outcome. Which reading this gets belongs to the shared rule
  # (the 54th pass moves it), so pin the property this pass owns: the pre-fit
  # description and the description of the fit made by the SAME call agree.
  d$b <- stats::rbinom(60, 1, 0.5)
  db <- maihda_describe(b ~ (1 | a1:a2), data = d, family = "binomial", weights = w)
  m <- suppressWarnings(suppressMessages(fit_maihda(
    b ~ (1 | a1:a2), data = d, family = "binomial", weights = w)))
  d_fit <- maihda_describe(m)
  expect_equal(db$outcome_overall$outcome_events, d_fit$outcome_overall$outcome_events)
  expect_equal(db$outcome_overall$outcome_trials, d_fit$outcome_overall$outcome_trials)
})
