# Regression tests for the 2026-07-29 audit findings.
#
#   1 [High] maihda_fit_diagnostics() reported converged = TRUE for EVERY returned
#            WeMixResults, on the premise that WeMix errors when optimisation
#            fails. It does not: its Newton loop exits at `max_iteration` into an
#            empty `if (iteration >= max_iteration) {}` block, and the linear path
#            discards bobyqa's return code, so a fit that never reached a maximum
#            comes back silently. A max_iteration = 0 binomial fit returned the
#            untouched starting values (lnl -253.292 vs -253.144, stratum variance
#            5.122 vs 5.729) and was labelled converged. The verdict now rests on
#            evidence -- WeMix's own gradient test, re-derived from the likelihood
#            function the fit carries -- and is NA when none is readable. Nonsense
#            iteration limits (0, -5, 2.7, "ten") are rejected up front.
#   2 [High] The longitudinal count VPC moved only the internally CENTERED time
#            column when building its reporting grid, updating the user's original
#            time column solely when a formula offset existed. But
#            maihda_longitudinal_formula() removes only the BARE raw-time
#            polynomial, so terms like x:wave survive on the raw axis and were then
#            evaluated at each row's OWN wave while the result was reported as "VPC
#            at time t". Both columns now move together.
#   3 [High] The crossed-dimensions MOR applied the independent-strata closed form
#            exp(sqrt(2 V) qnorm(.75)) to the SUM of the dimension and interaction
#            variances. Two crossed strata sharing a dimension share that
#            dimension's random effect, which cancels from their difference, so the
#            difference is a mixture rather than N(0, 2V). On a 2x2 design carrying
#            only dimension-A variance 1 the package reported 2.596 against an
#            exact 1.569 -- 65% high. The MOR now comes from the mixture under a
#            stated sampling scheme (two distinct observed strata, uniform).
#   4 [Med]  maihda_ic()'s delta guard checked outcome, family, sample ids and
#            response values but not weights, so the same Gaussian model fitted
#            unweighted and with weights alternating 0.5/2 was ranked against
#            itself with a delta of 63.7 AIC and no warning. It now compares both
#            weight fingerprints, as calculate_pcv() and compare_maihda() do.
#   5 [Med]  Documentation only: fixed-effect SEs from a single person-weight
#            column were called "design-consistent" for complex surveys such as
#            NHANES and PISA, which a wrapper with no PSUs, sampling strata,
#            higher-stage weights, FPCs or replicate weights cannot support.

# ---- Finding 1: WeMix convergence must be evidence-based ---------------------

test_that("a nonsense wemix max_iteration is rejected before fitting", {
  bad <- list(0, -5L, 2.7, "ten", c(1, 2), NA_integer_, numeric(0))
  for (mi in bad) {
    expect_error(
      MAIHDA:::maihda_validate_wemix_max_iteration(list(max_iteration = mi)),
      "whole number",
      info = paste("max_iteration =", paste(deparse(mi), collapse = ""))
    )
  }
  # Valid limits pass through, coerced to integer.
  ok <- MAIHDA:::maihda_validate_wemix_max_iteration(list(max_iteration = 5))
  expect_identical(ok$max_iteration, 5L)
  # Other dot arguments are untouched, and an absent limit is a no-op.
  passthru <- MAIHDA:::maihda_validate_wemix_max_iteration(list(nQuad = 13))
  expect_identical(passthru, list(nQuad = 13))
})

test_that("wemix convergence is NA without evidence and FALSE on bad evidence", {
  # The old code said TRUE for this object, which carries nothing at all about how
  # the optimisation ended.
  bare <- structure(list(varDF = data.frame(grp = c("stratum", "Residual"),
                                            var1 = c("(Intercept)", NA),
                                            vcov = c(0.4, 1))),
                    class = "WeMixResults")
  expect_true(is.na(MAIHDA:::maihda_wemix_convergence(bare)$converged))

  # A likelihood that is PRESENT but not finite is positive evidence of failure.
  broken <- structure(list(lnl = NaN), class = "WeMixResults")
  verdict <- MAIHDA:::maihda_wemix_convergence(broken)
  expect_false(verdict$converged)
  expect_match(verdict$messages, "not finite")

  # So are non-finite coefficients or standard errors.
  expect_false(MAIHDA:::maihda_wemix_convergence(
    structure(list(lnl = -10, coef = c(1, Inf)), class = "WeMixResults"))$converged)
  expect_false(MAIHDA:::maihda_wemix_convergence(
    structure(list(lnl = -10, coef = c(1, 2), SE = c(0.1, NaN)),
              class = "WeMixResults"))$converged)

  # An unreadable likelihood function leaves the verdict unknown rather than good.
  opaque <- structure(list(lnl = -10, coef = c(1, 2), SE = c(0.1, 0.2),
                           is_adaptive = TRUE, vars = 1,
                           lnlf = function(...) stop("unfamiliar WeMix version")),
                      class = "WeMixResults")
  expect_true(is.na(MAIHDA:::maihda_wemix_convergence(opaque)$converged))
})

test_that("wemix convergence tracks the gradient of a real fit", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  set.seed(7)
  d <- expand.grid(gender = paste0("g", 1:3), race = paste0("r", 1:4), rep = 1:30)
  cell <- interaction(d$gender, d$race, drop = TRUE)
  u <- stats::rnorm(nlevels(cell), 0, 1.5)
  d$x1 <- stats::rnorm(nrow(d))
  d$y <- stats::rbinom(nrow(d), 1,
                       stats::plogis(-0.3 + 1.2 * d$x1 + u[as.integer(cell)]))
  d$w <- stats::runif(nrow(d), 0.5, 3)

  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x1 + (1 | gender:race), data = d, engine = "wemix",
               sampling_weights = "w")))

  # A fit that reached the maximum is still positively confirmed -- the fix must
  # not degrade every weighted fit to "unknown".
  expect_true(isTRUE(fit$diagnostics$converged))
  expect_length(fit$diagnostics$messages, 0)
  crit_ok <- MAIHDA:::maihda_wemix_gradient_criterion(fit$model)
  expect_true(is.finite(crit_ok))
  expect_lt(crit_ok, MAIHDA:::.maihda_wemix_grad_tol)

  # Now the abandoned-optimisation case, built deterministically rather than by
  # hoping an optimizer stops early: displace the reported estimates off the
  # maximum and restate the likelihood AT those estimates, which is exactly what a
  # fit returned from an exhausted iteration budget looks like. The self-check
  # therefore passes and the gradient -- correctly -- does not.
  off <- fit$model
  off$coef <- off$coef + c(0.4, -0.3)
  off$lnl <- sum(as.numeric(off$lnlf(c(as.numeric(off$coef),
                                       as.numeric(off$vars)))))
  expect_gt(MAIHDA:::maihda_wemix_gradient_criterion(off),
            MAIHDA:::.maihda_wemix_grad_tol)
  diag_off <- maihda_fit_diagnostics(off)
  expect_false(diag_off$converged)
  expect_match(diag_off$messages, "away from a maximum")
  # And the user actually sees it.
  expect_true(any(grepl("away from a maximum",
                        MAIHDA:::maihda_format_fit_diagnostics(diag_off))))

  # The self-check must veto a verdict when the rebuilt objective disagrees with
  # the fit's own reported likelihood (an unfamiliar parameterisation).
  mismatched <- fit$model
  mismatched$lnl <- fit$model$lnl + 500
  expect_true(is.na(maihda_fit_diagnostics(mismatched)$converged))
})

# ---- Finding 2: the reporting grid must move BOTH time columns ---------------

test_that("maihda_longitudinal_set_time moves the original time column too", {
  nd <- data.frame(.maihda_ctime = c(0, 1, 2), wave = c(10, 11, 12), x = 1:3)
  moved <- MAIHDA:::maihda_longitudinal_set_time(nd, ".maihda_ctime", 3,
                                                 orig_time = "wave", center = 10)
  expect_equal(moved$.maihda_ctime, rep(3, 3))
  expect_equal(moved$wave, rep(13, 3))     # 3 + 10, in lockstep
  expect_equal(moved$x, 1:3)               # everything else untouched

  # Without centering the two names coincide and only one column exists.
  same <- MAIHDA:::maihda_longitudinal_set_time(nd, "wave", 4)
  expect_equal(same$wave, rep(4, 3))
  expect_equal(same$.maihda_ctime, c(0, 1, 2))
})

test_that("longitudinal count VPC(t) evaluates raw-time terms at the report time", {
  skip_on_cran()
  set.seed(11)
  nid <- 160
  g1 <- sample(c("F", "M"), nid, TRUE)
  g2 <- sample(c("A", "B", "C"), nid, TRUE)
  xv <- stats::rnorm(nid)
  d <- do.call(rbind, lapply(seq_len(nid), function(i)
    data.frame(id = i, g1 = g1[i], g2 = g2[i], x = xv[i], wave = 10:14)))
  cell <- interaction(d$g1, d$g2, drop = TRUE)
  u <- stats::rnorm(nlevels(cell), 0, 0.45)
  ri <- stats::rnorm(nid, 0, 0.35)
  d$y <- stats::rpois(nrow(d), exp(-0.4 + 0.10 * (d$wave - 10) + 0.5 * d$x -
                                     0.09 * d$x * (d$wave - 10) +
                                     u[as.integer(cell)] + ri[d$id]))

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x * wave + (1 | g1:g2), data = d, family = "poisson",
               id = "id", time = "wave")))
  lng <- m$longitudinal_info
  time_term <- MAIHDA:::maihda_lng_time_term(lng)
  center <- MAIHDA:::maihda_lng_time_center(lng)

  # The premise: centering is active AND a raw-time term survived in the fixed
  # part. Without both, this test proves nothing.
  expect_identical(time_term, ".maihda_ctime")
  expect_equal(center, 10)
  expect_true("x:wave" %in% attr(stats::terms(
    MAIHDA:::maihda_nobars(stats::formula(m$model))), "term.labels"))

  vt <- suppressWarnings(suppressMessages(summary(m)))$longitudinal$vpc_t

  # Rebuild VPC(t) independently, moving BOTH time columns; and, for contrast, the
  # way the defect did it -- centered column only.
  Sigma_s <- MAIHDA:::maihda_re_block_lme4(m$model, "stratum", time_term,
                                           lng$time_degree)
  Sigma_i <- MAIHDA:::maihda_re_block_lme4(m$model, lng$id, time_term,
                                           lng$time_degree)
  frame <- MAIHDA:::maihda_model_frame(m$model)
  wts <- MAIHDA:::maihda_fit_prior_weights(m$model)
  vpc_at <- function(t_orig, move_raw) {
    t_c <- t_orig - center
    vs <- MAIHDA:::maihda_var_at_time(Sigma_s, t_c)
    vi <- MAIHDA:::maihda_var_at_time(Sigma_i, t_c)
    nd <- frame
    nd[[time_term]] <- t_c
    if (move_raw) nd[["wave"]] <- t_orig
    lambda <- pmax(exp(MAIHDA:::maihda_lme4_fixed_link(m$model, nd) + (vs + vi) / 2),
                   .Machine$double.eps)
    vs / (vs + vi + MAIHDA:::maihda_weighted_obs_mean(log1p(1 / lambda), wts))
  }

  correct <- vapply(vt$time, vpc_at, numeric(1), move_raw = TRUE)
  buggy   <- vapply(vt$time, vpc_at, numeric(1), move_raw = FALSE)
  expect_equal(vt$estimate, correct, tolerance = 1e-10)
  # The two really do differ, so the assertion above has teeth: reverting the fix
  # would make the reported trajectory match `buggy` instead.
  expect_gt(max(abs(correct - buggy)), 1e-4)

  # The headline reference-time VPC is on the same footing.
  s <- suppressWarnings(suppressMessages(summary(m)))
  expect_equal(s$vpc$estimate, vpc_at(s$vpc$ref_time, TRUE), tolerance = 1e-10)
})

# ---- Finding 3: crossed MOR is a mixture, not the independent closed form ----

# Drive maihda_mor_crossed() with a chosen variance vector and cell grid, with no
# fitting: the arithmetic under test is entirely in the random-effect variances
# and the stratum grid.
maihda_a729_crossed_mor <- function(cells, tau_dim, tau_int) {
  nms <- colnames(cells)
  vn <- stats::setNames(c(tau_dim, tau_int), c(nms, "stratum"))
  stub <- list(engine = "lme4", data = as.data.frame(cells),
               cc_info = list(dim_groups = stats::setNames(nms, nms),
                              interaction_group = "stratum"))
  testthat::local_mocked_bindings(
    maihda_random_variances_lme4 = function(...) vn, .package = "MAIHDA")
  MAIHDA:::maihda_mor_crossed(stub)
}

test_that("crossed MOR reduces to the closed form when only the interaction varies", {
  # Two DISTINCT cells always have independent interaction effects, so with no
  # additive dimension variance there is nothing to correlate and the mixture must
  # collapse onto exp(sqrt(2 V) qnorm(.75)) exactly. This pins the fix as a
  # generalisation rather than a different quantity.
  cells <- expand.grid(d1 = 1:2, d2 = 1:3)
  expect_equal(maihda_a729_crossed_mor(cells, c(0, 0), 0.8),
               exp(sqrt(2 * 0.8) * stats::qnorm(0.75)), tolerance = 1e-6)
  # No heterogeneity anywhere -> the two strata never differ -> MOR is exactly 1.
  expect_equal(maihda_a729_crossed_mor(cells, c(0, 0), 0), 1)
})

test_that("crossed MOR matches the audit's exact 2x2 counterexample", {
  # 2x2, only dimension A varies. Two of the six pairs share A (difference 0), the
  # other four differ (difference ~ N(0, 2)), so the median |difference| solves
  # 1/3 + 2/3 P(|N(0,2)| <= x) = 1/2, i.e. x = sqrt(2) qnorm(0.625).
  cells <- expand.grid(d1 = 1:2, d2 = 1:2)
  got <- maihda_a729_crossed_mor(cells, c(1, 0), 0)
  expect_equal(got, exp(sqrt(2) * stats::qnorm(0.625)), tolerance = 1e-6)
  expect_equal(round(got, 3), 1.569)
  # The old variance-sum formula, for the record.
  expect_equal(round(exp(sqrt(2 * 1) * stats::qnorm(0.75)), 3), 2.596)
})

test_that("crossed MOR is below the closed form whenever dimensions carry variance", {
  grids <- list(expand.grid(d1 = 1:2, d2 = 1:2),
                expand.grid(d1 = 1:3, d2 = 1:4),
                expand.grid(d1 = 1:2, d2 = 1:3, d3 = 1:2))
  for (cells in grids) {
    D <- ncol(cells)
    tau_dim <- rep(0.3, D)
    got <- maihda_a729_crossed_mor(cells, tau_dim, 0.1)
    closed <- exp(sqrt(2 * (sum(tau_dim) + 0.1)) * stats::qnorm(0.75))
    expect_gt(got, 1)
    expect_lt(got, closed)

    # Independent reference: enumerate every unordered pair of cells.
    v <- unlist(lapply(seq_len(nrow(cells) - 1), function(i)
      vapply((i + 1):nrow(cells), function(j)
        2 * (0.1 + sum(as.numeric(cells[i, ] != cells[j, ]) * tau_dim)),
        numeric(1))))
    ref <- exp(stats::uniroot(
      function(x) mean(2 * stats::pnorm(x / sqrt(v)) - 1) - 0.5,
      c(1e-10, 10 * sqrt(max(v))), tol = .Machine$double.eps^0.5)$root)
    expect_equal(got, ref, tolerance = 1e-6)
  }
})

test_that("an incomplete stratum grid is handled, and canonical fits are unchanged", {
  skip_on_cran()
  # Cells missing from the data must simply not enter the pair set.
  full <- expand.grid(d1 = 1:3, d2 = 1:3)
  partial <- full[-c(2, 6), ]
  expect_false(isTRUE(all.equal(maihda_a729_crossed_mor(full, c(0.4, 0.25), 0.2),
                                maihda_a729_crossed_mor(partial, c(0.4, 0.25), 0.2))))

  # A canonical (single-stratum) fit keeps the closed form: no cc_info, no change.
  set.seed(606)
  n <- 900
  d <- data.frame(g = factor(sample(c("F", "M"), n, TRUE)),
                  r = factor(sample(c("A", "B", "C"), n, TRUE)))
  d$y <- stats::rbinom(n, 1, stats::plogis(
    c(F = 0.6, M = -0.6)[as.character(d$g)] +
      c(A = 0.5, B = 0, C = -0.5)[as.character(d$r)]))
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | g:r), data = d, family = "binomial")))
  expect_null(m$cc_info)
  v_a <- MAIHDA:::extract_between_variance(m)
  expect_equal(maihda_mor(m), exp(sqrt(2 * v_a) * stats::qnorm(0.75)))
})

# ---- Finding 4: maihda_ic() must not rank across different weights -----------

test_that("maihda_ic() withholds the delta across differing weights", {
  skip_on_cran()
  set.seed(3)
  n <- 300
  d <- data.frame(g = factor(sample(c("F", "M"), n, TRUE)),
                  r = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = stats::rnorm(n))
  cell <- interaction(d$g, d$r, drop = TRUE)
  d$y <- stats::rnorm(n, 1 + 0.5 * d$x +
                        stats::rnorm(nlevels(cell), 0, 0.6)[cell], 1)
  d$w <- rep(c(0.5, 2), length.out = n)

  unw <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g:r), data = d)))
  wtd <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g:r), data = d, weights = w)))

  # Same model, same rows, same outcome -- only the likelihood being maximised
  # differs, which every other comparability guard in the package already catches.
  expect_warning(res <- maihda_ic(unw, wtd, model_names = c("unweighted", "weighted")),
                 "prior weights")
  expect_false("delta" %in% names(res))
  # The per-model criteria are still reported.
  expect_true(all(is.finite(res$AIC)))

  # Identical weights must KEEP the delta: this is the canonical null-vs-adjusted
  # comparison, and the guard must not break it.
  adj <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g:r), data = d, weights = w)))
  nul <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ 1 + (1 | g:r), data = d, weights = w)))
  kept <- maihda_ic(adj, nul, model_names = c("adjusted", "null"))
  expect_true("delta" %in% names(kept))
  expect_equal(min(kept$delta), 0)
})

test_that("maihda_ic() withholds the delta across differing sampling weights", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  set.seed(3)
  n <- 300
  d <- data.frame(g = factor(sample(c("F", "M"), n, TRUE)),
                  r = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = stats::rnorm(n))
  cell <- interaction(d$g, d$r, drop = TRUE)
  d$y <- stats::rnorm(n, 1 + 0.5 * d$x +
                        stats::rnorm(nlevels(cell), 0, 0.6)[cell], 1)
  d$sw <- stats::runif(n, 0.5, 3)

  unw <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g:r), data = d)))
  dsn <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g:r), data = d, engine = "wemix",
               sampling_weights = "sw")))
  expect_warning(res <- maihda_ic(unw, dsn, model_names = c("unweighted", "design")),
                 "sampling weights")
  expect_false("delta" %in% names(res))
})
