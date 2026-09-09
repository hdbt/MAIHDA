# Audit 2026-09-09 (c) -- a formula offset() placed AFTER (1 | stratum) made the
# ordinal engine fit the stratum column as its offset.
#
# ordinal::clmm() builds its model frame from subbars(formula) -- the bars turned
# into `+` -- then OVERWRITES that frame's terms attribute with the terms of
# nobars(formula) and reads the offset with model.offset(), i.e. by POSITION in
# the bar-free variables list. The two lists only line up while every
# random-effect term follows every offset term.
#
# The package produced the opposite order at every entry point. The strata
# shorthand resolves through stats::update(fixed, . ~ . + (1 | stratum)), which
# appends the random effect and pushes offset() to the end; the adjusted formula
# built by maihda() does the same; and a user writing
# y ~ x + (1 | stratum) + offset(off) supplies it directly. The mis-index then
# lands on the grouping column:
#
#   * stratum is an integer id column, so the offset became the STRATUM ID --
#     numeric, hence no error. On the reference fit below the thresholds moved
#     from (-0.83, -0.01, 0.89, 1.79) to (4.85, 5.65, 6.54, 7.42), the
#     between-stratum variance inflated from 0.47 to 13.50 (= var(1:12), the
#     variance of the ids the random intercept had to cancel), and the
#     log-likelihood fell 41.3 units. Silently.
#   * when the mis-indexed column is a factor or character instead -- the
#     adjusted formula puts a dimension column there -- it errors outright, so
#     the two-model ordinal path could not be run with an offset at all.
#
# maihda_offset_before_bars() relocates the offset terms ahead of the bars for
# the clmm call only. It is a no-op -- returning the formula UNCHANGED -- when
# there is no offset, no random effect, or the offsets already lead, so no-offset
# fits stay bit-identical; lme4/brms/WeMix read their offset from their own frame
# and keep the formula they always had.

make_offord_data <- function(seed = 11, n = 900) {
  set.seed(seed)
  d <- data.frame(
    x = stats::rnorm(n),
    z = stats::runif(n, 1, 4),
    d1 = sample(c("a", "b"), n, TRUE),
    d2 = sample(c("p", "q"), n, TRUE),
    d3 = sample(c("u", "v", "w"), n, TRUE),
    stringsAsFactors = FALSE
  )
  st <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  d$off <- log(d$z)
  u <- stats::rnorm(nlevels(st), 0, 0.6)[st]
  eta <- 0.7 * d$x + d$off + u
  cu <- vapply(c(-0.8, 0.15, 1.1, 2.05), function(a) stats::plogis(a - eta),
               numeric(n))
  cu <- cbind(cu, 1)
  cu <- t(apply(cu - cbind(0, cu[, -5, drop = FALSE]), 1, cumsum))
  cu[, 5] <- 1
  d$y <- factor(max.col(stats::runif(n) <= cu, ties.method = "first"),
                levels = 1:5, ordered = TRUE)
  # A Gaussian twin on the same strata/offset for the lme4 confinement check.
  d$g <- 3 + 0.7 * d$x + d$off + u + stats::rnorm(n)
  d
}

fit_off <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(fit_maihda(formula, data = data, ...)))
}

# ---- the pure relocation helper ---------------------------------------------

test_that("maihda_offset_before_bars leaves a formula that needs no move IDENTICAL", {
  # Identity (not merely equality) is the contract: an untouched formula is what
  # keeps every no-offset ordinal fit bit-identical to before the fix.
  unchanged <- list(
    y ~ x + (1 | stratum),
    y ~ x + offset(off) + (1 | stratum),
    y ~ offset(off) + (1 | stratum),
    y ~ x + offset(off),
    y ~ (1 | stratum),
    y ~ x,
    y ~ 1
  )
  for (f in unchanged) {
    expect_identical(maihda_offset_before_bars(f), f)
  }

  # A non-formula argument passes straight through.
  expect_identical(maihda_offset_before_bars("not a formula"), "not a formula")
})

test_that("maihda_offset_before_bars moves offsets ahead of every bar term", {
  expect_equal(maihda_offset_before_bars(y ~ x + (1 | stratum) + offset(off)),
               y ~ x + offset(off) + (1 | stratum))

  # The adjusted formula shape: the dimension main effects sit between the bar
  # and the offset, and must keep their order relative to each other.
  expect_equal(
    maihda_offset_before_bars(y ~ x + (1 | stratum) + d1 + d2 + d3 + offset(off)),
    y ~ x + offset(off) + (1 | stratum) + d1 + d2 + d3)

  # An offset-only null model.
  expect_equal(maihda_offset_before_bars(y ~ (1 | stratum) + offset(off)),
               y ~ offset(off) + (1 | stratum))

  # Several offsets: relative order preserved, one already ahead of the bar.
  expect_equal(maihda_offset_before_bars(y ~ x + offset(a) + (1 | s) + offset(b)),
               y ~ x + offset(a) + offset(b) + (1 | s))

  # Several bars: the offset must clear the FIRST one.
  expect_equal(maihda_offset_before_bars(y ~ x + (1 | s) + (1 | t) + offset(a)),
               y ~ x + offset(a) + (1 | s) + (1 | t))

  # Double bars are random-effect terms too.
  expect_equal(maihda_offset_before_bars(y ~ x + (1 || s) + offset(a)),
               y ~ x + offset(a) + (1 || s))

  # A PARENTHESIZED offset counts. terms() finds an offset through `(` and `+`
  # (attr(terms(y ~ x + (1|s) + (offset(a))), "offset") is 4), so these spellings
  # reach clmm's positional read exactly like a bare one and must move too. Found
  # by a self-check after the first version of this fix looked only at bare
  # offset() operands and left `(offset(off))` fitting the stratum column.
  expect_equal(maihda_offset_before_bars(y ~ x + (1 | s) + (offset(a))),
               y ~ x + (offset(a)) + (1 | s))
  expect_equal(maihda_offset_before_bars(y ~ x + (1 | s) + ((offset(a)))),
               y ~ x + ((offset(a))) + (1 | s))
  # A `+` group carrying an offset moves whole. nobars() is unaffected by which
  # side of the bar the group sits on, so the fixed formula is unchanged.
  expect_equal(maihda_offset_before_bars(y ~ x + (1 | s) + (z + offset(a))),
               y ~ x + (z + offset(a)) + (1 | s))
  expect_equal(maihda_offset_before_bars(y ~ (1 | s) + (x + offset(a))),
               y ~ (x + offset(a)) + (1 | s))

  # NOT an offset: terms() does not look through an ordinary call, so these are
  # plain variables and must not be moved.
  expect_identical(maihda_offset_before_bars(y ~ x + (1 | s) + I(offset(a))),
                   y ~ x + (1 | s) + I(offset(a)))
  expect_identical(maihda_offset_before_bars(y ~ x + (1 | s) + log(offset(a))),
                   y ~ x + (1 | s) + log(offset(a)))

  # The formula environment survives the rebuild.
  e <- new.env()
  f <- y ~ x + (1 | stratum) + offset(off)
  environment(f) <- e
  expect_identical(environment(maihda_offset_before_bars(f)), e)
})

test_that("maihda_clmm_offset_aligned accepts what clmm can read and rejects the rest", {
  # Safe: no offset, no bar, or the offset already ahead of every bar (including
  # after the relocation above).
  aligned <- list(
    y ~ x + (1 | stratum),
    y ~ x + offset(off),
    y ~ x + offset(off) + (1 | stratum),
    maihda_offset_before_bars(y ~ x + (1 | stratum) + offset(off)),
    maihda_offset_before_bars(y ~ x + (1 | stratum) + (z + offset(off)))
  )
  for (f in aligned) {
    expect_true(maihda_clmm_offset_aligned(f))
  }

  # An offset inside a `-` subtree AFTER the random effect cannot be relocated
  # without changing the fixed design, so the backstop refuses it rather than
  # letting clmm read the grouping column.
  bad <- y ~ x + (1 | stratum) - offset(off)
  expect_identical(maihda_offset_before_bars(bad), bad)
  expect_false(maihda_clmm_offset_aligned(bad))

  # A `-` offset BEFORE the random effect is still fine and must not be rejected.
  expect_true(maihda_clmm_offset_aligned(y ~ x - offset(off) + (1 | stratum)))
})

# ---- the fit ----------------------------------------------------------------

test_that("the ordinal fit uses the formula offset, not the grouping column", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_offord_data()
  m <- fit_off(y ~ x + offset(off) + (1 | d1:d2:d3), d, family = "ordinal")

  # Reference: the same model written with the offset BEFORE the random effect,
  # which is the spelling clmm parses correctly.
  ref <- ordinal::clmm(y ~ x + offset(off) + (1 | stratum), data = m$data,
                       link = "logit", Hess = TRUE)

  expect_equal(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(ref)))
  expect_equal(unname(m$model$Theta), unname(ref$Theta))
  expect_equal(unname(m$model$beta), unname(ref$beta))
  expect_equal(unname(stats::fitted(m$model)), unname(stats::fitted(ref)))
  expect_equal(as.numeric(ordinal::VarCorr(m$model)[[1]]),
               as.numeric(ordinal::VarCorr(ref)[[1]]))

  # The formula clmm actually recorded puts the offset ahead of the bar.
  expect_false(is.null(
    attr(stats::terms(maihda_nobars(stats::formula(m$model))), "offset")))
  expect_identical(stats::formula(m$model),
                   maihda_offset_before_bars(stats::formula(m$model)))

  # Negative control on the exact old failure: pre-fix the fit was BIT-IDENTICAL
  # to one whose offset is the integer stratum id. It must no longer be, and the
  # gap is large (41.3 log-likelihood units on this data), so a partial revert
  # cannot slip through on tolerance.
  wrong <- ordinal::clmm(y ~ x + offset(stratum) + (1 | stratum), data = m$data,
                         link = "logit", Hess = TRUE)
  expect_gt(abs(as.numeric(stats::logLik(m$model)) -
                  as.numeric(stats::logLik(wrong))), 10)
  expect_gt(min(abs(m$model$Theta - wrong$Theta)), 1)
  # ... and the stratum variance no longer absorbs var(1:12).
  expect_lt(as.numeric(ordinal::VarCorr(m$model)[[1]]), 2)
})

test_that("a user-written trailing offset is normalized for clmm too", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_offord_data()
  ds <- make_strata(d, vars = c("d1", "d2", "d3"))$data

  # The strata shorthand is not involved here: the formula reaches fit_maihda()
  # already carrying (1 | stratum), with the offset written last.
  trailing <- fit_off(y ~ x + (1 | stratum) + offset(off), ds, family = "ordinal")
  leading <- fit_off(y ~ x + offset(off) + (1 | stratum), ds, family = "ordinal")

  expect_equal(as.numeric(stats::logLik(trailing$model)),
               as.numeric(stats::logLik(leading$model)))
  expect_equal(unname(trailing$model$Theta), unname(leading$model$Theta))
  expect_equal(unname(trailing$model$beta), unname(leading$model$beta))
})

test_that("a parenthesized trailing offset fits correctly too", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_offord_data()
  ds <- make_strata(d, vars = c("d1", "d2", "d3"))$data
  ref <- ordinal::clmm(y ~ x + offset(off) + (1 | stratum), data = ds,
                       link = "logit", Hess = TRUE)

  # (offset(off)) after the bar silently fitted the stratum column; (x + offset(off))
  # after the bar silently DROPPED the offset (it matched the no-offset fit).
  for (f in list(y ~ x + (1 | stratum) + (offset(off)),
                 y ~ (1 | stratum) + (x + offset(off)))) {
    m <- fit_off(f, ds, family = "ordinal")
    expect_equal(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(ref)))
    expect_equal(unname(m$model$Theta), unname(ref$Theta))
  }
})

test_that("an offset clmm cannot read is refused, not fitted against the wrong column", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_offord_data()
  ds <- make_strata(d, vars = c("d1", "d2", "d3"))$data

  # An offset inside a `-`, `*` or `:` subtree AFTER the random effect cannot be
  # relocated without changing the fixed design, so the fit is refused. Before the
  # backstop each of these silently fitted the stratum column as the offset.
  for (f in list(y ~ x + (1 | stratum) - offset(off),
                 y ~ x + (1 | stratum) + z * offset(off),
                 y ~ x + (1 | stratum) + z:offset(off))) {
    expect_error(
      suppressMessages(suppressWarnings(
        fit_maihda(f, data = ds, family = "ordinal"))),
      "cannot be moved ahead of the random-effect term")
  }

  # The same models written with the offset first are accepted and fit.
  ok <- fit_off(y ~ x + offset(off) + (1 | stratum), ds, family = "ordinal")
  expect_s3_class(ok$model, "clmm")
})

test_that("an offset-free ordinal fit is untouched by the normalization", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_offord_data()
  m <- fit_off(y ~ x + (1 | d1:d2:d3), d, family = "ordinal")

  # The resolved formula is passed through unchanged, so the fit is the one the
  # package always produced.
  expect_identical(maihda_offset_before_bars(m$formula), m$formula)
  ref <- ordinal::clmm(y ~ x + (1 | stratum), data = m$data, link = "logit",
                       Hess = TRUE)
  expect_equal(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(ref)))
  expect_equal(unname(m$model$Theta), unname(ref$Theta))
})

test_that("maihda() runs the two-model ordinal decomposition with an offset", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  d <- make_offord_data()
  # Pre-fix this errored: the adjusted formula's mis-indexed offset landed on the
  # character dimension column d3 ("non-numeric argument to binary operator").
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ x + offset(off) + (1 | d1:d2:d3), data = d, family = "ordinal")))

  null_ref <- ordinal::clmm(y ~ x + offset(off) + (1 | stratum),
                            data = a$model$data, link = "logit", Hess = TRUE)
  adj_ref <- ordinal::clmm(y ~ x + offset(off) + d1 + d2 + d3 + (1 | stratum),
                           data = a$model_adjusted$data, link = "logit",
                           Hess = TRUE)

  expect_equal(as.numeric(stats::logLik(a$model$model)),
               as.numeric(stats::logLik(null_ref)))
  expect_equal(as.numeric(stats::logLik(a$model_adjusted$model)),
               as.numeric(stats::logLik(adj_ref)))
  expect_true(is.finite(a$pcv$pcv))
})

test_that("the clmm-only normalization leaves the lme4 engine alone", {
  skip_on_cran()

  d <- make_offord_data()
  m <- fit_off(g ~ x + offset(off) + (1 | d1:d2:d3), d)

  # The shared formula assembly is untouched: the resolved formula still carries
  # the trailing offset, which lme4 reads correctly from its own model frame.
  rhs <- paste(deparse(m$formula[[3]]), collapse = " ")
  expect_true(grepl("(1 | stratum) + offset(off)", rhs, fixed = TRUE))

  ds <- make_strata(d, vars = c("d1", "d2", "d3"))$data
  ref <- lme4::lmer(g ~ x + offset(off) + (1 | stratum), data = ds, REML = TRUE)
  expect_equal(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(ref)))
  expect_equal(unname(lme4::fixef(m$model)), unname(lme4::fixef(ref)))
  expect_equal(as.numeric(lme4::VarCorr(m$model)$stratum),
               as.numeric(lme4::VarCorr(ref)$stratum))
})
