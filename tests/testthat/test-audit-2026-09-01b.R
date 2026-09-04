# Audit 2026-09-01 (second pass of the day; 2026-09-01 was taken by the ordinal
# proportional-odds pass, hence the "b" suffix).
#
# Finding [High]: equivalent binomial response encodings produced radically different
# AUCs. maihda_da_aggregated_counts() decided aggregation from the RESPONSE -- a row
# strictly inside (0, 1) -- so a fit whose every row is all-success or all-failure
# (y in {0, 1} with the counts in weights =) was read as a Bernoulli fit carrying
# "precision weights", which the AUC ignored.
#
# That shape is not exotic: it is what you get by collapsing individual 0/1 records to
# (covariate pattern x outcome) frequency cells. Because each such cell contributes one
# case and one control at an IDENTICAL score, every pair tied and the AUC was EXACTLY
# 0.5 -- no information -- printed with a confident note explaining why it was right.
#
# Measured on 3000 individuals in 12 strata: the individual 0/1 data, the cbind()
# spelling, and the coarse proportion+weights spelling all give AUC 0.786929 with
# 1427 / 1573; the same individuals collapsed one level finer gave 0.500000 with
# 24 / 24 -- from a fit whose fixed effects match to 4e-10. The reported AUC depended
# on collapse granularity, which is not a property of the model.
#
# The premise was verified rather than assumed: for the binomial family a prior weight
# IS a trial count -- a 600-row Bernoulli fit with weights 1:5 and the row-EXPANDED fit
# agree to 8e-14 in the fixed effects and 3e-13 in the fitted probabilities -- and
# ?glm says so outright. Non-integral weights cannot be counts and keep the
# observation-level path; binomial_weights overrides either way.

# Hand-built aggregated data: fixed counts, so the fits are deterministic and the
# cross-spelling assertions do not depend on a seed or on where the optimizer stops.
make_agg <- function() {
  succ <- c(28, 41, 12, 33,  9, 26, 35, 47, 18, 31, 22, 38)
  tot  <- c(50, 60, 40, 55, 45, 52, 48, 62, 44, 58, 47, 56)
  data.frame(stratum = factor(rep(1:6, each = 2)),
             x = rep(c(0, 1), times = 6),
             succ = succ, fail = tot - succ, tot = tot,
             p = succ / tot)
}

# The same data as individual 0/1 rows -- the ground truth the AUC must reproduce.
expand_agg <- function(a) {
  do.call(rbind, lapply(seq_len(nrow(a)), function(i) {
    data.frame(stratum = a$stratum[i], x = a$x[i],
               y = c(rep(1, a$succ[i]), rep(0, a$fail[i])))
  }))
}

# The same data collapsed to (stratum, x, outcome) frequency cells: y is exactly 0/1
# on every row and the count rides in weights =. This is the shape that broke.
cells_agg <- function(a) {
  out <- rbind(data.frame(stratum = a$stratum, x = a$x, y = 1, n = a$succ),
               data.frame(stratum = a$stratum, x = a$x, y = 0, n = a$fail))
  out[out$n > 0, ]
}

fit_b <- function(formula, data, ...) {
  suppressWarnings(suppressMessages(
    fit_maihda(formula, data = data, family = "binomial", ...)))
}

test_that("all four binomial response encodings give the same AUC and totals", {
  skip_on_cran()
  a    <- make_agg()
  ind  <- expand_agg(a)
  cell <- cells_agg(a)

  m_ind  <- fit_b(y ~ x + (1 | stratum), ind)
  m_cb   <- fit_b(cbind(succ, fail) ~ x + (1 | stratum), a)
  m_prop <- fit_b(p ~ x + (1 | stratum), a, weights = tot)
  m_cell <- fit_b(y ~ x + (1 | stratum), cell, weights = n)

  das   <- lapply(list(m_ind, m_cb, m_prop, m_cell), maihda_discriminatory_accuracy)
  aucs  <- vapply(das, function(d) d$auc, numeric(1))
  cases <- vapply(das, function(d) d$n_case, numeric(1))
  ctrls <- vapply(das, function(d) d$n_control, numeric(1))

  # A statistical summary must not change under likelihood-equivalent response
  # spelling. Before the fix the fourth entry was 0.5 with 12 / 12.
  expect_equal(max(abs(aucs - aucs[1])), 0, tolerance = 1e-6)
  expect_equal(cases, rep(sum(a$succ), 4))
  expect_equal(ctrls, rep(sum(a$fail), 4))
  # ... and the broken value really was the degenerate one, so this is not a
  # vacuous "they all agree because they are all 0.5" pass.
  expect_true(aucs[1] > 0.6)
  # None is flagged: all four are trial-count readings of the same data, so no
  # weights were discarded.
  expect_false(any(vapply(das, function(d) isTRUE(d$precision_weights_ignored),
                          logical(1))))

  # The fits themselves were always identical -- only the AUC diverged. Pinned so a
  # future change cannot "fix" the AUC by perturbing the model.
  expect_equal(unname(lme4::fixef(m_cell$model)), unname(lme4::fixef(m_cb$model)),
               tolerance = 1e-6)
  expect_equal(das[[4]]$mor, das[[2]]$mor, tolerance = 1e-6)
})

test_that("the frequency-cell AUC equals the row-expanded concordance exactly", {
  skip_on_cran()
  a    <- make_agg()
  cell <- cells_agg(a)
  m    <- fit_b(y ~ x + (1 | stratum), cell, weights = n)

  # Independent reconstruction: expand the cell rows by their counts and take the
  # ordinary rank AUC over the expanded 0/1 outcome. No package aggregation code in
  # the reference path.
  prob <- predict_maihda(m, type = "individual", scale = "response")
  idx  <- rep(seq_len(nrow(cell)), cell$n)
  auc  <- maihda_discriminatory_accuracy(m)$auc
  expect_equal(auc, maihda_auc(prob[idx], cell$y[idx]))

  # The exact failure mode: the observation-level AUC over these rows is 0.5, because
  # each (stratum, x) cell contributes one case and one control at the same score.
  expect_equal(maihda_auc(prob, cell$y), 0.5)
  expect_true(abs(auc - 0.5) > 0.1)
})

test_that("binomial_weights overrides the trial-count reading in both directions", {
  skip_on_cran()
  a    <- make_agg()
  cell <- cells_agg(a)
  m    <- fit_b(y ~ x + (1 | stratum), cell, weights = n)
  prob <- predict_maihda(m, type = "individual", scale = "response")

  auto   <- maihda_discriminatory_accuracy(m)
  forced <- maihda_discriminatory_accuracy(m, binomial_weights = "trials")
  expect_equal(auto$auc, forced$auc)
  expect_equal(auto$n_case, forced$n_case)

  # "analytic" restores the pre-fix observation-level reading, flagged as discarding
  # the weights, and says so in print().
  an <- maihda_discriminatory_accuracy(m, binomial_weights = "analytic")
  expect_equal(an$auc, maihda_auc(prob, cell$y))
  expect_equal(an$n_case, sum(cell$y == 1))
  expect_equal(an$n_control, sum(cell$y == 0))
  expect_true(isTRUE(an$precision_weights_ignored))
  expect_false(isTRUE(an$weighted))
  expect_output(print(an), "not read as trial counts")

  expect_error(maihda_discriminatory_accuracy(m, binomial_weights = "nope"), "'arg'")
})

test_that("non-integral weights cannot be trial counts and keep the observation path", {
  skip_on_cran()
  a    <- make_agg()
  cell <- cells_agg(a)
  cell$n <- cell$n + 0.5
  m    <- fit_b(y ~ x + (1 | stratum), cell, weights = n)
  prob <- predict_maihda(m, type = "individual", scale = "response")

  da <- maihda_discriminatory_accuracy(m)
  expect_equal(da$auc, maihda_auc(prob, cell$y))
  expect_equal(da$n_case, sum(cell$y == 1))
  expect_true(isTRUE(da$precision_weights_ignored))
  expect_null(MAIHDA:::maihda_da_aggregated_counts(m))

  # Forcing "trials" carries the fractional mass rather than rounding it into
  # observations that were never collected, and warns (2026-08-31 precedent).
  expect_warning(fda <- maihda_discriminatory_accuracy(m, binomial_weights = "trials"),
                 "not a whole number of successes")
  expect_equal(fda$n_case, sum(cell$n * cell$y))
  expect_false(isTRUE(fda$precision_weights_ignored))
})

test_that("binomial_weights does not touch cbind, unit-weight or proportion fits", {
  skip_on_cran()
  a   <- make_agg()
  ind <- expand_agg(a)

  # A cbind() response carries its denominator structurally; the argument is a no-op.
  m_cb <- fit_b(cbind(succ, fail) ~ x + (1 | stratum), a)
  base_cb <- maihda_discriminatory_accuracy(m_cb)
  for (bw in c("auto", "trials", "analytic")) {
    da <- maihda_discriminatory_accuracy(m_cb, binomial_weights = bw)
    expect_equal(da$n_case, sum(a$succ))
    expect_equal(da$auc, base_cb$auc)
    expect_false(isTRUE(da$precision_weights_ignored))
  }

  # An ordinary unweighted Bernoulli fit is unchanged in every mode.
  m_ind <- fit_b(y ~ x + (1 | stratum), ind)
  base  <- maihda_discriminatory_accuracy(m_ind)
  for (bw in c("auto", "trials", "analytic")) {
    da <- maihda_discriminatory_accuracy(m_ind, binomial_weights = bw)
    expect_equal(da$auc, base$auc)
    expect_equal(da$n_case, sum(ind$y == 1))
    expect_false(isTRUE(da$precision_weights_ignored))
  }

  # "analytic" on a genuine proportion response is incoherent -- there is no
  # observation-level case/control reading of 0.56 of a success -- and says so
  # instead of falling into maihda_auc()'s generic "must be binary" message.
  m_prop <- fit_b(p ~ x + (1 | stratum), a, weights = tot)
  expect_error(maihda_discriminatory_accuracy(m_prop, binomial_weights = "analytic"),
               "not available for a proportion response")
  expect_equal(maihda_discriminatory_accuracy(m_prop)$n_case, sum(a$succ))
})

test_that("an all-success / all-failure aggregated fit is read as trials", {
  skip_on_cran()
  # The reporter's explicit ask: every row all-success or all-failure. Here the
  # response is exactly 0/1 on EVERY row, so nothing but the weights can signal
  # aggregation.
  d <- data.frame(stratum = factor(rep(1:8, each = 2)),
                  x = rep(c(0, 1), times = 8),
                  y = c(1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1),
                  n = c(30, 44, 26, 38, 41, 22, 35, 47, 29, 33, 25, 40, 36, 28, 43, 31))
  m  <- fit_b(y ~ x + (1 | stratum), d, weights = n)
  da <- maihda_discriminatory_accuracy(m)
  expect_equal(da$n_case, sum(d$n[d$y == 1]))
  expect_equal(da$n_control, sum(d$n[d$y == 0]))
  expect_false(isTRUE(da$precision_weights_ignored))

  cnt <- MAIHDA:::maihda_da_aggregated_counts(m)
  expect_equal(cnt$trials, d$n)
  expect_equal(cnt$successes, d$n * d$y)

  # Structural equivalence with the cbind() spelling of the very same data.
  d2 <- data.frame(stratum = d$stratum, x = d$x,
                   succ = d$n * d$y, fail = d$n * (1 - d$y))
  m2 <- fit_b(cbind(succ, fail) ~ x + (1 | stratum), d2)
  expect_equal(da$auc, maihda_discriminatory_accuracy(m2)$auc, tolerance = 1e-6)
})
