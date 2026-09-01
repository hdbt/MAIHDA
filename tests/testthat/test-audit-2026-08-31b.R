# Audit 2026-08-31 (second pass) -- "aggregated-binomial AUC silently invents cases
# by rounding".
#
# maihda_da_aggregated_counts() recognises R's second aggregated-binomial idiom -- a
# PROPORTION response with the trial counts supplied as prior weights (?glm) -- and
# rebuilt the per-row success counts as round(y * w). It checked that the WEIGHTS
# were integral but never that y * w was, so a proportion that is not k / n for any
# integer k was rounded into whole observations that were never collected.
#
# CONFIRMED. 240 rows of "3.5 successes out of 7" carry a success mass of 840; the
# reconstruction returned 960. On a fitted model with half-integral successes the
# reported totals moved from 1140 / 540 to 1160 / 520 and the AUC from 0.6979 to
# 0.6625. Every such fit already draws glmer's "non-integer #successes in a binomial
# glm!" warning, so the data are known-malformed at fit time -- but the AUC path
# rounded them without a word.
#
# FIX: whole-number products are still snapped to the integer (the documented
# successes/trials spelling stays bit-identical to cbind(successes, failures)); a
# genuinely fractional product is carried into the weighted concordance AS-IS, with a
# warning. maihda_auc_weighted() already takes real-valued case/control mass, so
# nothing has to be fabricated to compute the AUC.

# ---- the pure boundary: maihda_da_proportion_successes() --------------------

test_that("a well-formed successes/trials proportion is reconstructed exactly and silently", {
  k <- c(0, 3, 7, 12, 20)
  n <- c(5, 8, 7, 40, 20)
  expect_silent(s <- MAIHDA:::maihda_da_proportion_successes(k / n, n))
  expect_identical(s, as.numeric(k))

  # The snap tolerance is RELATIVE: the representation error of (k / n) * n grows
  # with k, so a fixed absolute epsilon would stop recognising large counts.
  k_big <- c(1e6, 4e6, 12345678)
  n_big <- c(3e6, 7e6, 87654321)
  expect_silent(s_big <- MAIHDA:::maihda_da_proportion_successes(k_big / n_big, n_big))
  expect_identical(s_big, k_big)
})

test_that("a fractional success count is kept, not rounded into invented cases", {
  # The reported reproduction: 240 rows of 3.5 successes out of 7 trials.
  y <- rep(3.5 / 7, 240)
  w <- rep(7, 240)
  expect_warning(s <- MAIHDA:::maihda_da_proportion_successes(y, w),
                 "not a whole number of successes")
  expect_equal(sum(s), 840)
  # ... and specifically NOT the rounded reconstruction, which invented 120 cases.
  expect_false(isTRUE(all.equal(sum(s), sum(round(y * w)))))
  expect_equal(sum(round(y * w)), 960)

  # The mass is preserved row by row, not just in total.
  expect_equal(s, y * w)
})

test_that("only the fractional rows are left fractional; whole rows still snap", {
  y <- c(3 / 6, 3.5 / 7, 10 / 20, 1.5 / 5)
  w <- c(6, 7, 20, 5)
  expect_warning(s <- MAIHDA:::maihda_da_proportion_successes(y, w),
                 "in 2 of 4 rows")
  # rows 1 and 3 are exact integers, and are snapped (no floating-point residue)
  expect_identical(s[c(1, 3)], c(3, 10))
  # rows 2 and 4 keep their fractional mass
  expect_equal(s[c(2, 4)], c(3.5, 1.5))

  # The warning names the remedy, so a user is not left with a bare number.
  expect_warning(MAIHDA:::maihda_da_proportion_successes(y, w),
                 "cbind\\(successes, failures\\)")
})

# ---- the concordance takes fractional mass ---------------------------------

test_that("maihda_auc_weighted accepts fractional case/control mass", {
  # Two probability levels, half-integral mass in each: 1.5 cases of 3 at the low
  # score, 2.5 of 3 at the high one. Concordance by hand:
  #   n1 = 4, n0 = 2; controls strictly below the high level = 1.5;
  #   concordant = 1.5 * 0 + 2.5 * 1.5 + 0.5 * (1.5 * 1.5 + 2.5 * 0.5) = 5.5
  auc <- MAIHDA:::maihda_auc_weighted(c(0.2, 0.8), c(1.5, 2.5), c(3, 3))
  expect_equal(auc, 5.5 / (4 * 2))
  expect_true(is.finite(auc) && auc >= 0 && auc <= 1)
})

# ---- fractional totals must still print ------------------------------------

test_that("case/control totals render whether whole or fractional", {
  expect_identical(MAIHDA:::maihda_format_mass(840), "840")
  expect_identical(MAIHDA:::maihda_format_mass(5107.0000000001), "5107")
  expect_identical(MAIHDA:::maihda_format_mass(279.2496869), "279.25")
  expect_identical(MAIHDA:::maihda_format_mass(NA_real_), "NA")
  # Large whole totals must not go scientific ("1e+06 cases").
  expect_identical(MAIHDA:::maihda_format_mass(1e6), "1000000")

  # print.maihda_da() formatted the totals with "%d", which is an ERROR on a
  # fractional double -- the fix above would otherwise have broken printing.
  da <- structure(
    list(auc = 0.7, auc_scope = "model", auc_full = NULL, mor = 1.5,
         n_case = 5108.67, n_control = 4557.33, family = "binomial",
         link = "logit", engine = "lme4", weighted = FALSE, weight_type = NULL,
         precision_weights_ignored = FALSE, apparent = TRUE),
    class = "maihda_da")
  expect_output(print(da), "5108\\.67 / 4557\\.33")

  # maihda()'s own analysis print has a second copy of the same line, which the
  # self-check caught still on "%d".
  expect_output(
    MAIHDA:::maihda_print_analysis_da(list(discriminatory_accuracy = da)),
    "cases/controls: 5108\\.67/4557\\.33")
  # A whole-number fit still prints as counts there, not as "5107.00".
  da_whole <- da
  da_whole$n_case <- 5107
  da_whole$n_control <- 4559
  expect_output(
    MAIHDA:::maihda_print_analysis_da(list(discriminatory_accuracy = da_whole)),
    "cases/controls: 5107/4559")
})

# ---- end to end on fitted models -------------------------------------------

test_that("a half-integral aggregated fit reports its own mass, not a rounded one", {
  skip_on_cran()
  set.seed(42)
  strat <- factor(rep(seq_len(12), each = 20))
  u <- stats::rnorm(12, sd = 0.8)[strat]
  # k + 0.5 successes out of 7: no integer count produces these proportions.
  k_half <- pmin(pmax(round(7 * stats::plogis(u)), 0), 6) + 0.5
  d <- data.frame(stratum = strat, p = k_half / 7, n = rep(7L, 240))

  m <- suppressWarnings(suppressMessages(
    fit_maihda(p ~ (1 | stratum), data = d, family = "binomial", weights = n)))

  expect_warning(cnt <- MAIHDA:::maihda_da_aggregated_counts(m),
                 "not a whole number of successes")
  expect_equal(sum(cnt$successes), sum(d$p * d$n))
  expect_equal(sum(cnt$trials), sum(d$n))

  da <- suppressWarnings(maihda_discriminatory_accuracy(m))
  expect_equal(da$n_case, 1140)
  expect_equal(da$n_control, 540)
  # The rounded reconstruction reported 1160 / 520 -- 20 cases that do not exist.
  expect_false(isTRUE(all.equal(da$n_case, sum(round(d$p * d$n)))))

  # The warning is the user-facing half of the fix, and summary() reaches the DA
  # through maihda_try_optional(), which catches errors only. A suppressWarnings()
  # added there later would leave the fractional mass unannounced.
  w <- character(0)
  withCallingHandlers(
    summary(m),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("not a whole number of successes", w)))

  # The AUC is the concordance under the model's own binomial weighting.
  score <- predict_maihda(m, type = "individual", scale = "response")
  expect_equal(da$auc,
               MAIHDA:::maihda_auc_weighted(score, d$p * d$n, as.numeric(d$n)))
  # ... and differs materially from what rounding gave (0.6625).
  expect_gt(abs(da$auc - MAIHDA:::maihda_auc_weighted(score, round(d$p * d$n),
                                                      as.numeric(d$n))),
            0.01)
})

test_that("a well-formed proportion fit is unchanged and silent", {
  skip_on_cran()
  set.seed(7)
  strat <- factor(rep(seq_len(12), each = 20))
  u <- stats::rnorm(12, sd = 0.8)[strat]
  n <- sample(20:60, 240, replace = TRUE)
  k <- stats::rbinom(240, n, stats::plogis(u))
  d <- data.frame(stratum = strat, p = k / n, n = n, s = k, f = n - k)

  m_pr <- suppressMessages(
    fit_maihda(p ~ (1 | stratum), data = d, family = "binomial", weights = n))
  expect_silent(cnt <- MAIHDA:::maihda_da_aggregated_counts(m_pr))
  expect_identical(cnt$successes, as.numeric(k))

  da_pr <- maihda_discriminatory_accuracy(m_pr)
  expect_identical(da_pr$n_case, sum(as.numeric(k)))
  expect_identical(da_pr$n_control, sum(as.numeric(n - k)))
  expect_output(print(da_pr), paste0(sum(k), " / ", sum(n - k)))

  # The two spellings of the same aggregated data still agree exactly.
  m_cb <- suppressMessages(
    fit_maihda(cbind(s, f) ~ (1 | stratum), data = d, family = "binomial"))
  da_cb <- maihda_discriminatory_accuracy(m_cb)
  expect_equal(da_pr$auc, da_cb$auc)
  expect_equal(da_pr$n_case, da_cb$n_case)
  expect_equal(da_pr$n_control, da_cb$n_control)
})

test_that("Bernoulli fits are untouched by the aggregated path", {
  skip_on_cran()
  set.seed(3)
  strat <- factor(rep(seq_len(12), each = 20))
  u <- stats::rnorm(12, sd = 0.8)[strat]
  d <- data.frame(stratum = strat, y = stats::rbinom(240, 1, stats::plogis(u)))
  m <- suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, family = "binomial"))
  expect_null(MAIHDA:::maihda_da_aggregated_counts(m))
  da <- maihda_discriminatory_accuracy(m)
  expect_identical(da$n_case, sum(d$y == 1))
  expect_identical(da$n_control, sum(d$y == 0))
})
