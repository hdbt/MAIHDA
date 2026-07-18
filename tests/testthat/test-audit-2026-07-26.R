# Regression tests for the 2026-07-26 audit findings.
#
#   1 [P1] maihda_check_stratum_matches_dims() caught a contradicting supplied
#          'stratum' only where the dimensions mapped to an EXISTING training
#          stratum, and returned early for the WHOLE newdata as soon as any one
#          numeric dimension fell outside its training auto-bin range. So a
#          complete-but-unseen combination (rebuilt id NA) kept a supplied known
#          stratum's random effect, and a single out-of-range row disabled checking
#          for every other row beside it. The implied stratum now falls back to the
#          combination's own label, and the unresolvable-row exemption is per row.
#   2 [P1] test-ordinal-engine.R asserted that a copied training row with only its
#          'stratum' changed predicts under allow_new_levels = TRUE. Its dimension
#          columns still identified the original stratum, so the check (correctly)
#          rejected it and the non-brms integration job failed. The fixture now drops
#          the dimension columns; the contradiction itself is asserted here.
#   3 [P2] maihda_brms_individual_prediction() REPLACED dots$re_formula whenever a
#          grouping level was unseen, reintroducing random effects the caller had
#          explicitly excluded (re_formula = NA) or swapping a different term in for
#          the requested one. Zeroing an unseen level may now only NARROW the
#          caller's scope, and the lme4 spelling re.form is honoured in place.
#   4 [P2] maihda_fit_diagnostics() took a finite divergence count as sufficient
#          evidence of convergence. sum(na.rm = TRUE) is 0 for an absent or all-NA
#          divergence column, so a fit with no Rhat was reported converged on no
#          evidence at all. Only Rhat can now support a positive verdict.

# ---- Finding 1: supplied-stratum contradiction bypasses ----------------------

# A label-table object (no per-dimension columns), exercising the label-map branch
# of maihda_stratum_lookup(). "m x hi" is deliberately NOT a training stratum.
maihda_dimcheck_obj <- function() {
  list(
    strata_vars = c("d1", "d2"),
    strata_info = data.frame(stratum = c("1", "2"),
                             label = c("m x lo", "f x hi"),
                             stringsAsFactors = FALSE),
    strata_sep = " x ",
    strata_autobin_info = NULL)
}

# The same, with an auto-binned numeric dimension (breaks 0/50/100).
maihda_dimcheck_autobin_obj <- function() {
  list(
    strata_vars = c("gender", "income"),
    strata_info = data.frame(
      stratum = c("1", "2", "3", "4"),
      label = c("F x lo", "F x hi", "M x lo", "M x hi"),
      stringsAsFactors = FALSE),
    strata_sep = " x ",
    strata_autobin_info = list(
      income = list(breaks = c(0, 50, 100), labels = c("lo", "hi"))))
}

test_that("a supplied stratum is checked against an UNSEEN dimension combination", {
  obj <- maihda_dimcheck_obj()
  f <- MAIHDA:::maihda_check_stratum_matches_dims

  # "m x hi" is a complete combination the model never saw, so it maps to no
  # training stratum id. Supplying a KNOWN stratum for it used to pass unchecked.
  bad <- data.frame(stratum = "1", d1 = "m", d2 = "hi", stringsAsFactors = FALSE)
  expect_error(f(obj, bad), "not present when the model was fit")
  expect_error(f(obj, bad), "must match the dimension columns")

  # Naming that unseen combination by its own label is consistent, so it passes;
  # the downstream allow_new_levels logic then decides whether to allow it.
  ok <- data.frame(stratum = "m x hi", d1 = "m", d2 = "hi", stringsAsFactors = FALSE)
  expect_null(f(obj, ok))

  # A seen combination still reports the training stratum id, not a label.
  seen_bad <- data.frame(stratum = "2", d1 = "m", d2 = "lo", stringsAsFactors = FALSE)
  expect_error(f(obj, seen_bad), "identify stratum '1'")
})

test_that("one out-of-range row does not disable the check for the rest of newdata", {
  obj <- maihda_dimcheck_autobin_obj()
  f <- MAIHDA:::maihda_check_stratum_matches_dims

  # Row 1 contradicts (income 25 -> "lo" -> stratum 1, but stratum 3 supplied).
  contradiction <- data.frame(gender = "F", income = 25, stratum = "3",
                              stringsAsFactors = FALSE)
  expect_error(f(obj, contradiction), "identify stratum '1'")

  # Batched with an out-of-range row it must STILL error: the exemption for an
  # unresolvable row applies to that row only.
  out_of_range <- data.frame(gender = "F", income = 1e6, stratum = "1",
                             stringsAsFactors = FALSE)
  expect_error(f(obj, rbind(contradiction, out_of_range)), "identify stratum '1'")
  expect_error(f(obj, rbind(out_of_range, contradiction)), "identify stratum '1'")

  # The out-of-range row on its own stays exempt WITHIN this helper (its bin, and
  # so the stratum it implies, cannot be determined here) -- as does an incomplete
  # combination. It is not unremarked, though: since the 2026-07-27 audit the
  # caller separately warns about it via maihda_check_autobin_in_range(), ahead of
  # the supplied-'stratum' early return (see test-audit-2026-07-27.R).
  expect_null(f(obj, out_of_range))
  expect_null(f(obj, data.frame(gender = NA_character_, income = 25, stratum = "3",
                                stringsAsFactors = FALSE)))
  # ... and an in-range row that agrees is still accepted.
  expect_null(f(obj, data.frame(gender = "F", income = 25, stratum = "1",
                                stringsAsFactors = FALSE)))
})

# A fit whose F/B intersection is never observed, so that combination is complete
# but unseen at prediction time. The stratum effects are strong enough that lending
# one to the wrong intersection is a materially wrong prediction.
maihda_missing_cell_fixture <- function() {
  set.seed(4242)
  n <- 900
  d <- data.frame(g = sample(c("F", "M"), n, TRUE),
                  r = sample(c("A", "B", "C"), n, TRUE),
                  stringsAsFactors = FALSE)
  d <- d[!(d$g == "F" & d$r == "B"), ]
  sid <- droplevels(interaction(d$g, d$r))
  u <- rnorm(nlevels(sid), 0, 1.5)[as.integer(sid)]
  d$y <- 1 + 0.4 * (d$g == "M") + u + rnorm(nrow(d), sd = 0.6)
  d
}

test_that("prediction rejects a known stratum supplied for an unseen combination", {
  skip_on_cran()
  d <- maihda_missing_cell_fixture()
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ g + r + (1 | g:r), data = d)))
  expect_false("F x B" %in% maihda_known_strata(m))

  unseen <- data.frame(g = "F", r = "B", stringsAsFactors = FALSE)
  bad <- unseen
  bad$stratum <- as.character(m$strata_info$stratum[1])

  # Rejected with AND without the population-average opt-in: allow_new_levels
  # chooses how an unseen stratum is handled, it does not license pairing one
  # intersection's fixed effects with another intersection's random effect.
  expect_error(predict_maihda(m, newdata = bad, type = "link"),
               "not present when the model was fit")
  expect_error(predict_maihda(m, newdata = bad, type = "link",
                              allow_new_levels = TRUE),
               "not present when the model was fit")

  # The supported route is unchanged and still gives the zero-effect prediction.
  pa <- predict_maihda(m, newdata = unseen, type = "link", allow_new_levels = TRUE)
  expect_true(is.finite(as.numeric(pa)))

  # The two paths agree: the label the rebuild path assigns to an unseen
  # combination is exactly what the check accepts when it is supplied explicitly,
  # and it predicts identically. (This is why the implied stratum falls back to the
  # label rather than to NA.)
  round_trip <- unseen
  round_trip$stratum <- maihda_stratum_labels(unseen, m$strata_vars, m$strata_sep,
                                              m$strata_autobin_info)
  expect_equal(
    as.numeric(predict_maihda(m, newdata = round_trip, type = "link",
                              allow_new_levels = TRUE)),
    as.numeric(pa))

  # What the bypass was worth: the borrowed random effect was added to a
  # combination it does not belong to.
  re <- lme4::ranef(m$model)$stratum[, 1]
  expect_gt(abs(re[1]), 0.3)

  # No false positives: the fitted data, and a consistent supplied stratum, work.
  expect_length(predict_maihda(m, newdata = m$data, type = "link"), nrow(m$data))
  expect_true(is.finite(as.numeric(
    predict_maihda(m, newdata = m$data[1, , drop = FALSE], type = "link"))))
})

# ---- Finding 2: the ordinal fixture's contradiction ---------------------------

test_that("changing only the stratum on a copied training row is a contradiction", {
  skip_on_cran()
  d <- maihda_missing_cell_fixture()
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ g + r + (1 | g:r), data = d)))

  # This is the shape test-ordinal-engine.R used to assert would predict: the
  # dimension columns still identify the ORIGINAL stratum.
  nd <- m$data[1, , drop = FALSE]
  nd$stratum <- "999"
  expect_error(predict_maihda(m, newdata = nd, type = "link",
                              allow_new_levels = TRUE),
               "must match the dimension columns")

  # Dropping the dimension columns is the supported way to name a stratum directly
  # -- available when they are not also fixed effects, as in the null model below
  # (and in the ordinal fixture, whose fixed part is ~ x). The unseen stratum is
  # still rejected by default and yields the zero-effect prediction when opted into.
  null_m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ 1 + (1 | g:r), data = d)))
  nd2 <- null_m$data[1, , drop = FALSE]
  nd2[c("g", "r")] <- NULL
  nd2$stratum <- "999"
  expect_error(predict_maihda(null_m, newdata = nd2, type = "link"),
               "not present in the fitted model")
  expect_true(is.finite(as.numeric(
    predict_maihda(null_m, newdata = nd2, type = "link", allow_new_levels = TRUE))))
})

# ---- Finding 3: brms unseen-level handling vs the caller's scope --------------

# Stan-free: maihda_brms_individual_prediction() only needs $formula and $data, and
# the scope it hands to brms is captured by mocking the block-prediction helper.
maihda_scope_obj <- function() {
  list(
    formula = y ~ x + (1 | stratum) + (1 | school),
    data = data.frame(y = c(0, 1, 0, 1), x = c(0, 1, 0, 1),
                      stratum = c("1", "2", "1", "2"),
                      school = c("s1", "s2", "s1", "s2"),
                      stringsAsFactors = FALSE))
}

# The dots the helper would pass to brms for one row with an UNSEEN stratum and a
# SEEN school, given the dots the caller supplied.
maihda_scope_seen <- function(dots) {
  seen <- NULL
  testthat::local_mocked_bindings(
    maihda_brms_predict_rows = function(object, nd, scale, dots) {
      seen <<- dots
      rep(0, nrow(nd))
    })
  nd <- data.frame(x = 0, stratum = "999", school = "s1", stringsAsFactors = FALSE)
  MAIHDA:::maihda_brms_individual_prediction(
    maihda_scope_obj(), nd, "link", TRUE, dots)
  seen
}

test_that("unseen-level zeroing narrows the caller's re_formula instead of replacing it", {
  # No caller restriction: drop the unseen stratum, keep the seen school (unchanged).
  none <- maihda_scope_seen(list(allow_new_levels = TRUE))
  expect_s3_class(none$re_formula, "formula")
  expect_identical(all.vars(none$re_formula), "school")

  # re_formula = NA is a fixed-effects-only prediction. It must survive untouched;
  # it used to come back as "~ (1 | school)", reinstating an excluded effect.
  na_scope <- maihda_scope_seen(list(allow_new_levels = TRUE, re_formula = NA))
  expect_true(is.logical(na_scope$re_formula) && is.na(na_scope$re_formula))

  # Only the stratum requested, and the stratum is the unseen one -> nothing is
  # left, i.e. fixed-only. school must NOT be substituted in.
  only_str <- maihda_scope_seen(
    list(allow_new_levels = TRUE, re_formula = ~ (1 | stratum)))
  expect_true(is.logical(only_str$re_formula) && is.na(only_str$re_formula))

  # A requested term that is NOT the unseen one is passed through as written.
  only_sch <- maihda_scope_seen(
    list(allow_new_levels = TRUE, re_formula = ~ (1 | school)))
  expect_s3_class(only_sch$re_formula, "formula")
  expect_identical(all.vars(only_sch$re_formula), "school")

  # Both requested -> the unseen one is dropped, the seen one kept.
  both <- maihda_scope_seen(
    list(allow_new_levels = TRUE, re_formula = ~ (1 | stratum) + (1 | school)))
  expect_s3_class(both$re_formula, "formula")
  expect_identical(all.vars(both$re_formula), "school")
})

test_that("the lme4 spelling re.form is honoured in place, not shadowed", {
  # re.form = NA is fixed-only. A conflicting re_formula must not be injected
  # alongside it (brms lets re.form win on some paths and ignores it on others, so
  # the two spellings disagreeing would make the scale change the estimand).
  d <- maihda_scope_seen(list(allow_new_levels = TRUE, re.form = NA))
  expect_false("re_formula" %in% names(d))
  expect_true(is.logical(d$re.form) && is.na(d$re.form))

  # A narrowed re.form is written back under re.form, still never both.
  d2 <- maihda_scope_seen(
    list(allow_new_levels = TRUE, re.form = ~ (1 | stratum) + (1 | school)))
  expect_false("re_formula" %in% names(d2))
  expect_s3_class(d2$re.form, "formula")
  expect_identical(all.vars(d2$re.form), "school")
})

test_that("an uninterpretable re_formula is left to the caller, not overridden", {
  d <- maihda_scope_seen(list(allow_new_levels = TRUE, re_formula = "nonsense"))
  expect_identical(d$re_formula, "nonsense")
})

test_that("maihda_brms_requested_re resolves brms's re_formula contract", {
  f <- MAIHDA:::maihda_brms_requested_re
  expect_null(f(list())$bars)                       # absent    -> every term
  expect_identical(f(list())$name, "re_formula")
  expect_null(f(list(re_formula = NULL))$bars)      # NULL      -> every term
  expect_identical(f(list(re_formula = NA))$bars, list())   # NA -> no terms
  expect_length(f(list(re_formula = ~ (1 | a) + (1 | b)))$bars, 2L)
  # re.form wins when both are given, matching brms's own aliasing.
  expect_identical(f(list(re.form = NA, re_formula = ~ (1 | a)))$name, "re.form")
  expect_false(f(list(re_formula = 42))$understood)
})

# ---- Finding 4: brms convergence verdict needs Rhat ---------------------------

# The verdict maihda_fit_diagnostics() reaches for a brmsfit with the given Rhat
# values and divergence column, without fitting anything.
maihda_brms_verdict <- function(rhat_val, div_val) {
  testthat::local_mocked_bindings(
    rhat = function(x, ...) rhat_val,
    nuts_params = function(x, ...) {
      if (is.null(div_val)) {
        data.frame(Parameter = "accept_stat__", Value = 0.9)
      } else {
        data.frame(Parameter = "divergent__", Value = div_val)
      }
    },
    .package = "brms")
  MAIHDA:::maihda_fit_diagnostics(structure(list(), class = "brmsfit"))
}

test_that("a clean divergence count alone is not evidence of convergence", {
  skip_if_not_installed("brms")

  # Rhat present and fine -> converged; Rhat over the threshold -> not converged.
  expect_true(maihda_brms_verdict(c(a = 1.001, b = 1.004), 0)$converged)
  expect_false(maihda_brms_verdict(c(a = 1.05), 0)$converged)

  # No Rhat, zero divergences: nothing establishes that the chains mixed, so the
  # verdict is unknown -- this used to report a clean converged fit.
  no_rhat <- maihda_brms_verdict(NA_real_, 0)
  expect_true(is.na(no_rhat$converged))
  expect_length(no_rhat$messages, 0L)

  # An all-NA or absent divergence column sums to 0 under na.rm = TRUE; neither may
  # be read as "no divergences observed".
  expect_true(is.na(maihda_brms_verdict(NA_real_, NA_real_)$converged))
  expect_true(is.na(maihda_brms_verdict(NA_real_, NULL)$converged))

  # A diagnostic that DID flag a problem still yields FALSE, with the message.
  div <- maihda_brms_verdict(NA_real_, 3)
  expect_false(div$converged)
  expect_match(div$messages, "divergent transition", all = FALSE)

  # Divergences are still reported alongside a fine Rhat.
  both <- maihda_brms_verdict(c(a = 1.001), 5)
  expect_false(both$converged)
  expect_length(both$messages, 1L)
})
