# Audit 2026-09-10 (third pass): the ordinal offset relocation gave up on a
# formula whose right-hand side ends in a subtraction.
#
# maihda_offset_before_bars() (added 2026-09-09) moves every offset() term ahead
# of the random-effect terms, because ordinal::clmm() reads the offset by POSITION
# in the bar-free variables list while indexing the frame it built from the barred
# one. It split the right-hand side on top-level `+` only -- deliberately, so that
# nothing is promoted out of a subtraction -- which meant a right-hand side whose
# OUTERMOST operator is `-` presented itself as a single operand and no relocation
# happened at all.
#
# `- 1` is exactly what stats::update() leaves behind for a no-intercept formula,
# so a user writing
#
#     maihda(ord ~ 0 + x + a + b + offset(log(expo)) + (1 | a:b), engine = "ordinal")
#
# reached the engine as `ord ~ x + a + b + (1 | stratum) + offset(log(expo)) - 1`
# and was REFUSED by the alignment backstop -- a fit that is perfectly
# relocatable. The refusal was not gratuitous: handed to clmm() as written, that
# formula fits the stratum id as its offset, and on the data below the
# log-likelihood falls from -596.76 to -619.83 and the between-stratum variance
# inflates from 0.154 to 3.341. The backstop was right; the relocator was
# incomplete.
#
# The relocation now walks through a trailing subtraction, reordering only the
# minuend and only when no subtracted operand carries an offset -- an offset
# written INSIDE the subtraction must not be lifted out, because that would change
# the fixed design, and it stays for the backstop.

make_ordminus_data <- function(seed = 7, n = 600) {
  set.seed(seed)
  d <- data.frame(
    x    = stats::rnorm(n),
    expo = stats::runif(n, 0.5, 2),
    a    = factor(sample(letters[1:3], n, TRUE)),
    b    = factor(sample(LETTERS[1:3], n, TRUE)))
  d$st <- interaction(d$a, d$b, drop = TRUE)
  u <- stats::rnorm(nlevels(d$st), 0, 0.8)[as.integer(d$st)]
  lp <- 0.7 * d$x + u + log(d$expo)
  cuts <- stats::quantile(lp + stats::rlogis(n), c(0, 0.3, 0.6, 1))
  d$ord <- factor(cut(lp + stats::rlogis(n), breaks = cuts, include.lowest = TRUE),
                  ordered = TRUE)
  levels(d$ord) <- c("lo", "mid", "hi")
  d
}

# The design a formula implies once the bars are stripped, plus its offset --
# the two things a relocation must leave untouched.
ordminus_design <- function(f, data) {
  fixed <- maihda_nobars(f)
  mf <- stats::model.frame(fixed, data)
  list(x = stats::model.matrix(fixed, mf), offset = stats::model.offset(mf))
}


test_that("a trailing subtraction hides the whole right-hand side from the splitter", {
  # The splitter's contract is unchanged -- it still refuses to split a `-` -- so
  # a right-hand side ending in one is a single operand. That is the mechanism.
  f <- ord ~ x + a + (1 | st) + offset(log(expo)) - 1
  expect_length(maihda_rhs_plus_terms(f[[length(f)]]), 1L)
  expect_length(maihda_rhs_plus_terms(quote(x + a + (1 | st) + offset(log(expo)))), 4L)

  # ... and the relocation now reaches into it anyway.
  moved <- maihda_offset_before_bars(f)
  expect_identical(deparse1(moved),
                   "ord ~ x + a + offset(log(expo)) + (1 | st) - 1")
  expect_false(maihda_clmm_offset_aligned(f))
  expect_true(maihda_clmm_offset_aligned(moved))
})


test_that("relocation through a subtraction preserves the fixed design exactly", {
  d <- make_ordminus_data()
  spellings <- list(
    ord ~ x + a + (1 | st) + offset(log(expo)) - 1,
    ord ~ x + a + (1 | st) + offset(log(expo)) - b,
    ord ~ x + a + (1 | st) + offset(log(expo)) - b - 1,
    ord ~ x + a - b + (1 | st) + offset(log(expo)),
    ord ~ x + a + (1 | st) + offset(log(expo)))
  for (f in spellings) {
    moved <- maihda_offset_before_bars(f)
    before <- ordminus_design(f, d)
    after <- ordminus_design(moved, d)
    expect_identical(colnames(before$x), colnames(after$x))
    expect_equal(max(abs(before$x - after$x)), 0)
    expect_equal(max(abs(before$offset - after$offset)), 0)
    expect_true(maihda_clmm_offset_aligned(moved))
  }
})


test_that("an offset inside the subtraction is never lifted out", {
  # Moving it would change the fixed design, so the relocation declines and the
  # backstop keeps its teeth. Same for an offset buried in a `*` subtree.
  #
  # A HELPER-LEVEL invariant, not a spelling a user meets: measured end to end,
  # maihda_resolve_strata_formula() rewrites both of these into the canonical
  # `ord ~ x + a + (1 | stratum) + offset(log(expo))` before the engine is called,
  # so fit_maihda() fits them (and did so before this change too, byte-identically).
  # The guard matters for any formula that reaches maihda_fit_clmm() un-normalised.
  keep <- list(
    ord ~ x + a + (1 | st) - offset(log(expo)),
    ord ~ x + (1 | st) + a * offset(log(expo)))
  for (f in keep) {
    expect_identical(deparse1(maihda_offset_before_bars(f)), deparse1(f))
    expect_false(maihda_clmm_offset_aligned(maihda_offset_before_bars(f)))
  }
  d <- make_ordminus_data()
  for (f in keep) {
    resolved <- suppressWarnings(suppressMessages(
      maihda_resolve_strata_formula(f, d, autobin = FALSE)))$formula
    expect_true(maihda_clmm_offset_aligned(maihda_offset_before_bars(resolved)))
  }
})


test_that("a formula with nothing to move is returned unchanged", {
  # The bit-identical guarantee the no-offset path relies on: the same object,
  # not merely an equal one.
  unchanged <- list(
    ord ~ x + a + (1 | st),
    ord ~ x + a + offset(log(expo)) + (1 | st),
    ord ~ x + a + offset(log(expo)),
    ord ~ x + a + (1 | st) - 1,
    ord ~ 0 + x + a + offset(log(expo)) + (1 | st))
  for (f in unchanged) {
    expect_identical(maihda_offset_before_bars(f), f)
  }
})


test_that("the refused spelling now fits, and fits the SAME model", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_ordminus_data()
  canonical <- suppressWarnings(suppressMessages(
    maihda(ord ~ x + a + b + offset(log(expo)) + (1 | a:b), data = d,
           autobin = FALSE, engine = "ordinal", family = "ordinal")))
  # Before the fix this call errored: "The offset() term in this formula cannot be
  # moved ahead of the random-effect term automatically".
  no_int <- suppressWarnings(suppressMessages(
    maihda(ord ~ 0 + x + a + b + offset(log(expo)) + (1 | a:b), data = d,
           autobin = FALSE, engine = "ordinal", family = "ordinal")))

  ll <- function(m) as.numeric(stats::logLik(m$model))
  sv <- function(m) as.numeric(m$model$ST[[1]])^2
  # clmm needs an intercept and assumes one -- its thresholds ARE the intercepts --
  # so the `0 +` spelling is the same model, and must give the same fit.
  #
  # Tolerance 1e-3 is set from clmm's OWN floor, measured rather than guessed: two
  # clmm parameterisations of one model, fitted directly with no package code in
  # between, differ by 3.5e-09 in the log-likelihood and 3.9e-05 RELATIVE in the
  # stratum variance, and maihda() reproduces those figures to the digit. clmm
  # converges its variance components far more loosely than lmer does; the defect
  # this pins moves the same variance by a factor of 21.
  expect_equal(ll(canonical$model), ll(no_int$model), tolerance = 1e-6)
  expect_equal(ll(canonical$model_adjusted), ll(no_int$model_adjusted),
               tolerance = 1e-6)
  expect_equal(sv(canonical$model_adjusted), sv(no_int$model_adjusted),
               tolerance = 1e-3)
  expect_equal(canonical$pcv$pcv, no_int$pcv$pcv, tolerance = 1e-3)
  # The null models are the SAME formula, so they agree exactly, not merely closely.
  expect_identical(deparse1(canonical$formula), deparse1(no_int$formula))
  expect_equal(sv(canonical$model), sv(no_int$model), tolerance = 1e-12)
})


test_that("the mis-fit the backstop was catching is real", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Not a formality: handed the un-relocated formula, clmm fits the stratum id as
  # the offset. This is what would have happened had the relocation been made to
  # "succeed" by simply weakening the guard.
  d <- make_ordminus_data()
  res <- suppressWarnings(suppressMessages(
    maihda_resolve_strata_formula(
      ord ~ 0 + x + a + b + offset(log(expo)) + (1 | a:b), d, autobin = FALSE)))
  bad_f <- res$formula
  expect_false(maihda_clmm_offset_aligned(bad_f))
  good_f <- maihda_offset_before_bars(bad_f)
  expect_true(maihda_clmm_offset_aligned(good_f))

  bad <- suppressWarnings(ordinal::clmm(bad_f, data = res$data, Hess = TRUE))
  good <- suppressWarnings(ordinal::clmm(good_f, data = res$data, Hess = TRUE))
  expect_lt(as.numeric(stats::logLik(bad)), as.numeric(stats::logLik(good)) - 15)
  expect_gt(as.numeric(bad$ST[[1]])^2, 10 * as.numeric(good$ST[[1]])^2)
})


test_that("a cumulative model is not warned about a missing grand mean", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # The 2026-09-10 grand-mean guard warns when the fixed design cannot represent
  # the outcome mean. A cumulative model has no intercept COLUMN -- its free
  # thresholds are the intercepts -- so `0 +` is a no-op there and the warning
  # would be a false alarm. Proven against clmm directly rather than argued.
  d <- make_ordminus_data()
  m_no <- suppressWarnings(ordinal::clmm(ord ~ 0 + x + (1 | st), data = d, Hess = TRUE))
  m_yes <- suppressWarnings(ordinal::clmm(ord ~ x + (1 | st), data = d, Hess = TRUE))
  expect_equal(as.numeric(stats::logLik(m_no)), as.numeric(stats::logLik(m_yes)),
               tolerance = 1e-8)
  expect_equal(as.numeric(m_no$ST[[1]])^2, as.numeric(m_yes$ST[[1]])^2,
               tolerance = 1e-6)

  warned <- function(expr) {
    ws <- character(0)
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        ws <<- c(ws, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) invokeRestart("muffleMessage"))
    sum(grepl("cannot represent the outcome", ws))
  }
  expect_identical(warned(fit_maihda(ord ~ 0 + x + (1 | st), data = d,
                                     engine = "ordinal", family = "ordinal",
                                     autobin = FALSE)), 0L)
  # The Gaussian fit on the same data still warns -- the gate is the family, not
  # a blanket disabling.
  d$y <- as.numeric(d$x) + stats::rnorm(nrow(d))
  expect_identical(warned(fit_maihda(y ~ 0 + x + (1 | st), data = d,
                                     autobin = FALSE)), 1L)
})
