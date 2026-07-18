# Regression tests for the 2026-07-27 audit findings.
#
#   1 [P1] maihda_prepare_prediction_data() range-checked auto-binned numeric
#          dimensions only on the REBUILD path. A caller-supplied 'stratum'
#          returned early before that check, so a row whose numeric dimension lay
#          outside every training bin was predicted SILENTLY, wearing the supplied
#          stratum's random effect -- while the identical row WITHOUT the stratum
#          column was correctly rejected. maihda_check_stratum_matches_dims() could
#          not catch it either: an out-of-range value cut()s to NA, so the implied
#          stratum is NA and the row is exempt there as "unresolvable". The check
#          now sits ahead of the early return, covering both paths: it still ERRORS
#          on the rebuild path (no stratum can be derived at all) but only WARNS
#          when a stratum was supplied, so callers relying on the old behaviour
#          keep their prediction until the escalation in a future release.
#   2 [P3] test-pcv_importance.R's print/plot fixture ran a seeded n_perm = 60
#          Monte-Carlo attribution under suppressMessages() only. That fixture warns
#          (largest MC_SE ~ 0.022 > 0.01) for every seed, so the package's own
#          convergence warning escaped as a stray suite warning; it is asserted with
#          expect_warning() there now.

# ---- Finding 1: supplied stratum bypassed the auto-bin bounds ------------------

# Tertile-binned `age` over [18, 80] crossed with gender. The stratum effects are
# large relative to the residual, so lending one to a row the bins cannot contain
# is a materially wrong prediction, not a rounding difference.
maihda_autobin_bounds_fixture <- function(n = 120) {
  set.seed(1015)
  d <- data.frame(age = runif(n, 18, 80),
                  gender = sample(c("F", "M"), n, replace = TRUE),
                  stringsAsFactors = FALSE)
  strata <- interaction(cut(d$age, breaks = 3), d$gender, drop = TRUE)
  d$y <- 2 + 0.1 * d$age + rnorm(nlevels(strata), sd = 0.6)[strata] +
    rnorm(n, sd = 0.4)
  d
}

test_that("maihda_check_autobin_in_range() flags only genuine out-of-range values", {
  f <- MAIHDA:::maihda_check_autobin_in_range
  info <- list(income = list(breaks = c(0, 50, 100), labels = c("lo", "hi")))

  # In range, including both closed boundaries (cut(include.lowest = TRUE)).
  expect_null(f(data.frame(income = c(0, 25, 100)), info))
  # NA is missing, not out of range -- it is the rebuild path's own concern.
  expect_null(f(data.frame(income = NA_real_), info))
  # No recipe, or the source column absent, leaves nothing to check.
  expect_null(f(data.frame(income = 1e6), NULL))
  expect_null(f(data.frame(other = 1e6), info))

  # No stratum supplied -> hard error: with no bin there is no stratum to predict
  # with, so there is nothing to return. Long-standing behaviour.
  expect_error(f(data.frame(income = 1e6), info),
               "outside the training auto-bin ranges", fixed = TRUE)
  expect_error(f(data.frame(income = -1), info), "income outside \\[0, 100\\]")
  # One bad row taints the batch rather than being dropped silently.
  expect_error(f(data.frame(income = c(25, 1e6)), info), "income outside")
  expect_error(f(data.frame(income = 1e6), info, has_stratum = FALSE),
               "Cannot rebuild 'stratum'", fixed = TRUE)

  # A supplied stratum WARNS instead, so callers who relied on the previously
  # silent behaviour keep their prediction; the message says it will escalate.
  expect_warning(f(data.frame(income = 1e6), info, has_stratum = TRUE),
                 "outside the training auto-bin ranges", fixed = TRUE)
  expect_warning(f(data.frame(income = 1e6), info, has_stratum = TRUE),
                 "will become an error in a future release", fixed = TRUE)
  # ... and it is a warning, not an error: the call returns rather than aborting.
  expect_null(suppressWarnings(f(data.frame(income = 1e6), info,
                                 has_stratum = TRUE)))
})

test_that("a supplied stratum does not bypass the numeric auto-bin bounds", {
  skip_on_cran()
  d <- maihda_autobin_bounds_fixture()
  model <- fit_maihda(y ~ age + gender + (1 | age:gender), data = d)
  known <- as.character(model$original_data$stratum[1])

  # The finding: identical rows, differing only in whether 'stratum' is supplied.
  # age = 1e9 falls in no training bin, so the intersection the row names cannot
  # exist in the fitted model. Neither row may pass silently -- but the severities
  # differ: the rebuild path has no stratum to predict with and errors, while the
  # supplied path warns and still returns, so existing callers are not broken.
  raw <- data.frame(age = 1e9, gender = "F", stringsAsFactors = FALSE)
  supplied <- data.frame(age = 1e9, gender = "F", stratum = known,
                         stringsAsFactors = FALSE)
  expect_error(predict_maihda(model, newdata = raw),
               "outside the training auto-bin ranges", fixed = TRUE)
  expect_warning(predict_maihda(model, newdata = supplied),
                 "outside the training auto-bin ranges", fixed = TRUE)
  expect_warning(predict_maihda(model, newdata = supplied),
                 "will become an error in a future release", fixed = TRUE)

  # The prediction is still produced -- the point of warning rather than erroring.
  # Note this holds because THIS formula uses the raw numeric `age` as its fixed
  # effect, so nothing references the (NA) binned factor. A fit that uses the
  # binned dimension as a fixed effect still fails downstream; see the adjusted
  # model in the null-model test below.
  expect_true(is.finite(as.numeric(
    suppressWarnings(predict_maihda(model, newdata = supplied)))))

  # Below the range as well as above it.
  expect_warning(
    predict_maihda(model, newdata = data.frame(age = 1, gender = "F",
                                               stratum = known,
                                               stringsAsFactors = FALSE)),
    "outside the training auto-bin ranges", fixed = TRUE)

  # Not relaxed by allow_new_levels, and reported for stratum-level predictions
  # too: a row in no bin is not a "new level" the population average stands in for.
  expect_warning(
    predict_maihda(model, newdata = supplied, allow_new_levels = TRUE),
    "outside the training auto-bin ranges", fixed = TRUE)
  expect_warning(
    predict_maihda(model, newdata = supplied, type = "strata"),
    "outside the training auto-bin ranges", fixed = TRUE)

  # A batch is flagged on account of its one bad row, not quietly predicted.
  row1 <- model$original_data[1, , drop = FALSE]
  expect_warning(
    predict_maihda(model, newdata = data.frame(
      age = c(row1$age, 1e9), gender = c(as.character(row1$gender), "F"),
      stratum = c(known, known), stringsAsFactors = FALSE)),
    "outside the training auto-bin ranges", fixed = TRUE)
})

test_that("the auto-bin bounds check does not reject legitimate predictions", {
  skip_on_cran()
  d <- maihda_autobin_bounds_fixture()
  model <- fit_maihda(y ~ age + gender + (1 | age:gender), data = d)
  row1 <- model$original_data[1, , drop = FALSE]
  known <- as.character(row1$stratum)

  # A row that agrees with its supplied stratum still predicts, and identically to
  # the same row with the stratum rebuilt from the dimensions.
  agreeing <- data.frame(age = row1$age, gender = as.character(row1$gender),
                         stratum = known, stringsAsFactors = FALSE)
  expect_equal(
    as.numeric(predict_maihda(model, newdata = agreeing)),
    as.numeric(predict_maihda(model, newdata = agreeing[, c("age", "gender")])),
    tolerance = 1e-8)

  # The closed boundaries of the training range are in range, not out of it.
  bounds <- range(model$strata_autobin_info$age$breaks)
  expect_length(
    predict_maihda(model, newdata = data.frame(age = bounds,
                                               gender = c("F", "M"),
                                               stringsAsFactors = FALSE)), 2L)

  # Predicting on the stored frames is unaffected.
  expect_length(predict_maihda(model, newdata = model$original_data), nrow(d))
  expect_length(predict_maihda(model), nrow(d))

  # A model with no auto-binned dimension never reaches the check.
  set.seed(88)
  d2 <- data.frame(gender = sample(c("F", "M"), 120, TRUE),
                   race = sample(c("a", "b"), 120, TRUE),
                   stringsAsFactors = FALSE)
  d2$y <- rnorm(120)
  m2 <- suppressMessages(fit_maihda(y ~ gender + race + (1 | gender:race), data = d2))
  expect_length(predict_maihda(m2, newdata = d2[1:10, c("gender", "race")]), 10L)
})

test_that("the bypass is reported on the null model, where it looked benign", {
  skip_on_cran()
  d <- maihda_autobin_bounds_fixture()
  ana <- suppressMessages(maihda(y ~ (1 | age:gender), data = d))
  nullm <- ana$model

  # The null model's fixed part is an intercept, so an out-of-range row returns a
  # perfectly plausible-looking value (intercept + the supplied stratum's random
  # effect) rather than an obviously absurd extrapolation -- which is exactly why
  # it needs to be flagged rather than left to the reader to notice.
  kn <- as.character(nullm$original_data$stratum[1])
  expect_warning(
    predict_maihda(nullm, newdata = data.frame(age = 1e9, gender = "F",
                                               stratum = kn,
                                               stringsAsFactors = FALSE)),
    "outside the training auto-bin ranges", fixed = TRUE)

  # The ADJUSTED model also carries the binned dimension as a FIXED effect, so an
  # out-of-range row has an NA .maihda_dim_age and the call cannot complete: the
  # warning fires, then lme4 fails on the unbuildable model matrix. That downstream
  # failure is not new (it is what happened before this check existed at all), but
  # it is why warning here only defers the problem -- the eventual escalation to an
  # error is what replaces an opaque lme4 message with an explanatory one.
  adj_warning <- NULL
  expect_error(
    withCallingHandlers(
      predict_maihda(ana$model_adjusted,
                     newdata = data.frame(age = 1e9, gender = "F",
                                          stratum = as.character(
                                            ana$model_adjusted$original_data$stratum[1]),
                                          stringsAsFactors = FALSE)),
      warning = function(cond) {
        adj_warning <<- conditionMessage(cond)
        invokeRestart("muffleWarning")
      }))
  expect_match(adj_warning, "outside the training auto-bin ranges", fixed = TRUE)

  # With the dimension columns absent there is nothing to range-check (just as
  # there is nothing to cross-check), so a bare supplied stratum is still trusted.
  expect_length(predict_maihda(nullm, newdata = data.frame(
    stratum = kn, stringsAsFactors = FALSE)), 1L)
})

test_that("an incomplete combination stays exempt, unlike an out-of-range one", {
  skip_on_cran()
  d <- maihda_autobin_bounds_fixture()
  ana <- suppressMessages(maihda(y ~ (1 | age:gender), data = d))
  nullm <- ana$model
  kn <- as.character(nullm$original_data$stratum[1])

  # An NA dimension genuinely says nothing, so a supplied stratum resolves the row
  # and it predicts SILENTLY. This is the distinction the fix turns on: out-of-range
  # is not unresolvable, it is resolvably outside every stratum the bins can
  # contain, so it is a contradiction and gets warned about (asserted above).
  na_row <- data.frame(age = NA_real_, gender = "F", stratum = kn,
                       stringsAsFactors = FALSE)
  expect_length(predict_maihda(nullm, newdata = na_row), 1L)
  expect_warning(predict_maihda(nullm, newdata = na_row), NA)

  # And the within-helper exemption is unchanged -- the warning comes from the
  # caller's separate maihda_check_autobin_in_range() call, not from here.
  expect_null(MAIHDA:::maihda_check_stratum_matches_dims(
    nullm, data.frame(age = 1e9, gender = "F", stratum = kn,
                      stringsAsFactors = FALSE)))
})
