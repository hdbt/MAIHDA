# Audit 2026-09-26 -- the proportional-odds statistic must not depend on how a
# covariate was SPELLED in the formula.
#
# FINDING (external audit, F07). maihda_ordinal_po_stat() rebuilt the fitted model's
# fixed-effects formula and fitted it against the clmm MODEL FRAME:
#
#   dat <- stats::model.frame(model)            # holds a column named "log(x)"
#   maihda_po_lrt(dat, resp, fixed_terms)       # fits y ~ log(x) against it
#
# A model frame stores every variable ALREADY EVALUATED under its deparsed name, so a
# fit written y ~ log(x) + (1 | stratum) leaves a column literally called "log(x)"
# and no column x at all. Asking clm() for y ~ log(x) on that frame makes R look up
# the symbol x, which the frame does not have, so evaluation escapes into the
# formula's enclosing environment. That went wrong two different ways:
#
#   * usually x was nowhere to be found, the refit errored, maihda_po_lrt()'s
#     tryCatch turned the error into NULL, and the check was silently LOST -- while
#     the identical model written y ~ lx + (1 | stratum) with lx <- log(x)
#     precomputed kept it. maihda_proportional_odds_test() then stopped with a
#     message blaming "a null, covariate-free model", which was not the reason.
#
#   * when a same-named object of the same length DID happen to be reachable -- the
#     chain from the package namespace runs through imports and base to the global
#     environment -- the statistic was computed from THAT object instead of from the
#     data the model was fitted to, and returned with no complaint at all.
#
# Scope: every term that is not a bare column name. Measured on HEAD, log(x),
# I(x^2), scale(x), poly(x, 2), factor(g) and log(x) + g all returned NULL, while
# bare columns, bare factors and their interactions were unaffected.
#
# FIX. maihda_po_refit_design() takes the fixed-effects columns from
# stats::model.matrix() against the model frame -- which matches each variable by its
# deparsed NAME instead of evaluating it -- and hands them on under generated
# syntactic names. The statistic is a likelihood ratio between a proportional and a
# nominal-effects fit, and both are invariant to any full-rank recoding of the same
# column space, so nothing moves for the spellings that already worked: all four were
# identical() to HEAD.

make_po_spell_data <- function(seed = 4, n = 360L) {
  set.seed(seed)
  d <- data.frame(g1 = rep(c("a", "b", "c"), length.out = n),
                  g2 = rep(c("p", "q", "r"), each = n / 3L))
  d <- make_strata(d, vars = c("g1", "g2"))$data
  d$xraw <- rlnorm(n, 0, 0.6)
  d$lx <- log(d$xraw)                 # the PRE-COMPUTED twin of log(xraw)
  d$grp <- factor(rep(c("lo", "mid", "hi"), length.out = n),
                  levels = c("lo", "mid", "hi"))
  u <- rnorm(nlevels(factor(d$stratum)), 0, 0.5)
  lat <- 0.8 * d$lx + 0.5 * (d$grp == "hi") +
    u[as.integer(factor(d$stratum))] + rlogis(n)
  d$y <- factor(cut(lat, c(-Inf, -0.8, 0.6, Inf)), labels = 1:3, ordered = TRUE)
  d
}

fit_po_spell <- function(formula, data) {
  suppressMessages(suppressWarnings(
    fit_maihda(formula, data = data, family = "ordinal")))
}

# ---- the premise: the two spellings are ONE model -----------------------------

test_that("log(x) and a precomputed log column are the same cumulative fit", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_spell_data()
  m_pre <- fit_po_spell(y ~ lx + (1 | stratum), d)
  m_form <- fit_po_spell(y ~ log(xraw) + (1 | stratum), d)

  # Same likelihood, same slope, same thresholds: any difference in a statistic
  # computed from these two fits is purely representational.
  expect_equal(as.numeric(stats::logLik(m_form$model)),
               as.numeric(stats::logLik(m_pre$model)), tolerance = 1e-10)
  expect_equal(unname(m_form$model$beta), unname(m_pre$model$beta),
               tolerance = 1e-10)
  expect_equal(unname(maihda_clmm_cutpoints(m_form$model)),
               unname(maihda_clmm_cutpoints(m_pre$model)), tolerance = 1e-10)

  # ... and the mechanism: the frame carries the EVALUATED column, not the raw one.
  fr <- stats::model.frame(m_form$model)
  expect_true("log(xraw)" %in% names(fr))
  expect_false("xraw" %in% names(fr))
  expect_equal(as.numeric(fr[["log(xraw)"]]), log(d$xraw))
})

# ---- the statistic no longer depends on the spelling -------------------------

test_that("the proportional-odds statistic is invariant to how a term is written", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_spell_data()

  ref <- maihda_ordinal_po_stat(fit_po_spell(y ~ lx + (1 | stratum), d)$model)
  expect_false(is.null(ref))

  # log(xraw) IS lx: same fit, so the same statistic, exactly.
  got <- maihda_ordinal_po_stat(fit_po_spell(y ~ log(xraw) + (1 | stratum), d)$model)
  expect_false(is.null(got))
  expect_equal(got$lrt, ref$lrt, tolerance = 1e-12)
  expect_identical(got$df, ref$df)

  # scale() is an affine recoding of the same column, so the LRT cannot move.
  sc <- maihda_ordinal_po_stat(fit_po_spell(y ~ scale(lx) + (1 | stratum), d)$model)
  expect_false(is.null(sc))
  expect_equal(sc$lrt, ref$lrt, tolerance = 1e-12)

  # factor(grp) IS grp, which was already a factor.
  rf <- maihda_ordinal_po_stat(fit_po_spell(y ~ grp + (1 | stratum), d)$model)
  ff <- maihda_ordinal_po_stat(fit_po_spell(y ~ factor(grp) + (1 | stratum), d)$model)
  expect_false(is.null(ff))
  expect_equal(ff$lrt, rf$lrt, tolerance = 1e-12)
  expect_identical(ff$df, rf$df)

  # and a transformed term alongside a bare one.
  both_pre <- maihda_ordinal_po_stat(fit_po_spell(y ~ lx + grp + (1 | stratum), d)$model)
  both_form <- maihda_ordinal_po_stat(
    fit_po_spell(y ~ log(xraw) + grp + (1 | stratum), d)$model)
  expect_false(is.null(both_form))
  expect_equal(both_form$lrt, both_pre$lrt, tolerance = 1e-12)
  expect_identical(both_form$df, both_pre$df)

  # I() against a precomputed square: an outside oracle, not just "not NULL".
  d$lxsq <- d$lx^2
  sq_form <- maihda_ordinal_po_stat(fit_po_spell(y ~ I(lx^2) + (1 | stratum), d)$model)
  sq_pre <- maihda_ordinal_po_stat(fit_po_spell(y ~ lxsq + (1 | stratum), d)$model)
  expect_false(is.null(sq_form))
  expect_equal(sq_form$lrt, sq_pre$lrt, tolerance = 1e-12)
  expect_identical(sq_form$df, sq_pre$df)

  # A MULTI-COLUMN term needs an oracle of its own: poly(lx, 2) is an orthogonal
  # basis for span{1, lx, lx^2}, and the LRT is invariant to a full-rank recoding
  # of the same space, so it must equal the raw {lx, lx^2} basis exactly.
  po <- maihda_ordinal_po_stat(fit_po_spell(y ~ poly(lx, 2) + (1 | stratum), d)$model)
  raw <- maihda_ordinal_po_stat(fit_po_spell(y ~ lx + lxsq + (1 | stratum), d)$model)
  expect_false(is.null(po))
  expect_equal(po$lrt, raw$lrt, tolerance = 1e-10)
  expect_identical(po$df, raw$df)
  # and the design's columns really do span that basis, both directions
  des <- maihda_po_refit_design(
    fit_po_spell(y ~ poly(lx, 2) + (1 | stratum), d)$model)
  X <- as.matrix(des$data[, des$terms, drop = FALSE])
  expect_lt(max(abs(stats::lm.fit(cbind(1, X), cbind(d$lx, d$lxsq))$residuals)), 1e-8)
  expect_lt(max(abs(stats::lm.fit(cbind(1, d$lx, d$lxsq), X)$residuals)), 1e-8)
})

test_that("every non-bare term spelling still yields a usable statistic", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_spell_data()
  specs <- list(
    "log(xraw)"   = y ~ log(xraw) + (1 | stratum),
    "I(lx^2)"     = y ~ I(lx^2) + (1 | stratum),
    "scale(lx)"   = y ~ scale(lx) + (1 | stratum),
    "poly(lx, 2)" = y ~ poly(lx, 2) + (1 | stratum),
    "factor(grp)" = y ~ factor(grp) + (1 | stratum)
  )
  for (nm in names(specs)) {
    s <- maihda_ordinal_po_stat(fit_po_spell(specs[[nm]], d)$model)
    expect_false(is.null(s), info = nm)
    expect_true(is.finite(s$lrt) && s$lrt >= 0, info = nm)
    expect_true(s$df >= 1, info = nm)
    # n_terms counts TERMS, not design columns: poly(lx, 2) is one term of two.
    expect_identical(s$n_terms, 1L, info = nm)
  }
  # a multi-column transformation contributes its columns to the df
  expect_identical(
    maihda_ordinal_po_stat(fit_po_spell(y ~ poly(lx, 2) + (1 | stratum), d)$model)$df,
    2L)
})

# ---- the refit design is the FITTED design -----------------------------------

test_that("the refit design carries the fitted columns, not re-evaluated ones", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_spell_data()
  m <- fit_po_spell(y ~ log(xraw) + grp + (1 | stratum), d)
  des <- maihda_po_refit_design(m$model)
  expect_false(is.null(des))

  # every design column has a syntactic, generated name a formula can carry
  expect_identical(des$terms, make.names(des$terms))
  expect_identical(des$resp, ".po_y")
  expect_identical(des$n_terms, 2L)
  expect_identical(nrow(des$data), nrow(stats::model.frame(m$model)))

  # the transformed column equals the FITTED values, and there is no intercept
  expect_false("(Intercept)" %in% names(des$data))
  expect_equal(as.numeric(des$data[[des$terms[1]]]), log(d$xraw))
  expect_identical(des$data[[des$resp]], stats::model.frame(m$model)[[1L]])

  # a covariate-free fit has no design at all
  m0 <- fit_po_spell(y ~ (1 | stratum), d)
  expect_null(maihda_po_refit_design(m0$model))
  expect_null(maihda_ordinal_po_stat(m0$model))
})

# ---- a rank-deficient design still yields its statistic ----------------------

test_that("aliased design columns are dropped rather than losing the statistic", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Taking the columns from model.matrix() gives up a rescue clm() performs for a
  # FORMULA: it drops an aliased factor TERM and warns, but will not do that for
  # plain numeric columns. y ~ grp + factor(grp) is the same variable twice, and it
  # produced a statistic through the old formula route -- so the design has to be
  # reduced to full rank itself, or the fix would have been a regression here.
  d <- make_po_spell_data()
  d$dup <- d$lx                                  # an exact copy
  d$konst <- 2                                   # zero variance

  ref_grp <- maihda_ordinal_po_stat(fit_po_spell(y ~ grp + (1 | stratum), d)$model)
  ref_lx <- maihda_ordinal_po_stat(fit_po_spell(y ~ lx + (1 | stratum), d)$model)
  expect_false(is.null(ref_grp))
  expect_false(is.null(ref_lx))

  # the same variable spelled twice spans the same space, so same LRT and same df
  twice <- maihda_ordinal_po_stat(
    fit_po_spell(y ~ grp + factor(grp) + (1 | stratum), d)$model)
  expect_false(is.null(twice))
  expect_equal(twice$lrt, ref_grp$lrt, tolerance = 1e-10)
  expect_identical(twice$df, ref_grp$df)

  # a duplicated numeric column, and a constant one, carry no information either
  for (f in list(y ~ lx + dup + (1 | stratum), y ~ lx + konst + (1 | stratum))) {
    got <- maihda_ordinal_po_stat(fit_po_spell(f, d)$model)
    expect_false(is.null(got), info = deparse(f))
    expect_equal(got$lrt, ref_lx$lrt, tolerance = 1e-10, info = deparse(f))
    expect_identical(got$df, ref_lx$df, info = deparse(f))
  }

  # the design keeps only the independent columns ...
  des <- maihda_po_refit_design(
    fit_po_spell(y ~ grp + factor(grp) + (1 | stratum), d)$model)
  expect_identical(length(des$terms), 2L)       # not the 4 model.matrix emits
  expect_identical(des$n_terms, 2L)             # n_terms still counts TERMS

  # ... and a covariate carrying no information at all is treated as none
  expect_null(maihda_po_refit_design(
    fit_po_spell(y ~ konst + (1 | stratum), d)$model))
})

# ---- the descriptive proxy is no longer silently dropped ----------------------

test_that("a transformed term keeps the descriptive proportional-odds proxy", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_spell_data()
  for (f in list(y ~ log(xraw) + (1 | stratum), y ~ factor(grp) + (1 | stratum))) {
    ad <- maihda_adequacy_checks(fit_po_spell(f, d)$model)
    expect_false(is.null(ad), info = deparse(f))
    expect_true("marginal_po_proxy" %in% names(ad), info = deparse(f))
    # still descriptive only: no flag, no p-value (audit 2026-08-03)
    expect_null(ad$marginal_po_proxy$flag)
    expect_null(ad$marginal_po_proxy$p)
  }
})

# ---- the calibrated test runs, and runs the SAME bootstrap -------------------

test_that("the bootstrap test is identical across two spellings of one model", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_spell_data()
  m_pre <- fit_po_spell(y ~ lx + (1 | stratum), d)
  m_form <- fit_po_spell(y ~ log(xraw) + (1 | stratum), d)

  r_pre <- suppressWarnings(
    maihda_proportional_odds_test(m_pre, n_sim = 19, seed = 7))
  r_form <- suppressWarnings(
    maihda_proportional_odds_test(m_form, n_sim = 19, seed = 7))

  expect_s3_class(r_form, "maihda_po_test")
  expect_identical(r_form$n_failed, 0L)
  expect_identical(r_form$n_sim, 19L)
  # the whole null distribution, not just the summary: same fit, same seed, so the
  # replicates must agree draw for draw.
  expect_equal(r_form$lrt, r_pre$lrt, tolerance = 1e-12)
  expect_equal(r_form$null_lrt, r_pre$null_lrt, tolerance = 1e-12)
  expect_identical(r_form$p_value, r_pre$p_value)
  expect_identical(r_form$n_terms, r_pre$n_terms)

  # a multi-column transformation bootstraps too
  r_poly <- suppressWarnings(maihda_proportional_odds_test(
    fit_po_spell(y ~ poly(lx, 2) + (1 | stratum), d), n_sim = 9, seed = 3))
  expect_s3_class(r_poly, "maihda_po_test")
  expect_identical(r_poly$df, 2L)
  expect_true(r_poly$n_sim >= 1L)
})

# ---- the two failure modes are told apart -----------------------------------

test_that("a covariate-free fit and a failed reconstruction get different messages", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_spell_data()

  m0 <- fit_po_spell(y ~ (1 | stratum), d)
  expect_error(maihda_proportional_odds_test(m0, n_sim = 5),
               "covariate-free")
  # the covariate-free message must NOT claim a convergence failure ...
  expect_error(maihda_proportional_odds_test(m0, n_sim = 5),
               "no covariate slopes to test")
  # ... and a fit WITH covariates must not be told it is covariate-free.
  m <- fit_po_spell(y ~ log(xraw) + (1 | stratum), d)
  expect_error(maihda_proportional_odds_test(m, n_sim = 19, seed = 1), NA)

  # The other branch IS reachable from a package fit -- a rank-deficient design
  # such as y ~ lx + dup, with dup an exact copy of lx, returns NULL from
  # maihda_po_lrt() -- but such a fit depends on the optimiser stopping in a
  # particular state, so the failing state is INJECTED here instead. That is what
  # makes the message assertion deterministic rather than optimiser-dependent.
  local_mocked_bindings(maihda_po_lrt = function(...) NULL, .package = "MAIHDA")
  expect_error(maihda_proportional_odds_test(m, n_sim = 3),
               "failed to converge")
  expect_error(maihda_proportional_odds_test(m, n_sim = 3),
               "could not be computed")
  # and that branch must NOT reach for the covariate-free wording
  expect_error(maihda_proportional_odds_test(m, n_sim = 3), "^(?!.*covariate-free).*$",
               perl = TRUE)
})

# ---- no leakage out of the fitted data --------------------------------------

test_that("a same-named object in the global environment is never used", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # The broken lookup escaped the model frame into the enclosing environment, whose
  # chain from this namespace reaches globalenv(); so the object has to live THERE,
  # not in this test's frame, for the defect to be reachable at all.
  skip_if(exists("a7_leak_v", envir = globalenv(), inherits = FALSE),
          "a7_leak_v already exists in the global environment")

  d <- make_po_spell_data()
  d$a7_leak_v <- d$xraw
  m <- fit_po_spell(y ~ log(a7_leak_v) + (1 | stratum), d)

  # what the FITTED data implies, computed without any formula transformation
  truth <- maihda_po_lrt(data.frame(y = d$y, v = log(d$a7_leak_v)), "y", "v")
  expect_false(is.null(truth))

  decoy <- rev(d$a7_leak_v)
  # control: the decoy really would give a different answer, so the assertion below
  # is not vacuous
  wrong <- maihda_po_lrt(data.frame(y = d$y, v = log(decoy)), "y", "v")
  expect_false(is.null(wrong))
  expect_gt(abs(wrong$lrt - truth$lrt), 1e-6)

  assign("a7_leak_v", decoy, envir = globalenv())
  on.exit(if (exists("a7_leak_v", envir = globalenv(), inherits = FALSE)) {
    rm("a7_leak_v", envir = globalenv())
  }, add = TRUE)

  got <- maihda_ordinal_po_stat(m$model)
  expect_false(is.null(got))
  expect_equal(got$lrt, truth$lrt, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(got$lrt, wrong$lrt, tolerance = 1e-6)))
})
