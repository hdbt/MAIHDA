# Regression tests for the 2026-07-20 audit findings on commit f0231ed.
#
#   #1 P1  maihda_longitudinal_pcv() gated its PCV denominators with only
#          `variance > 0`, so a null growth block on the singularity boundary
#          (between-stratum variance ~0) yielded PCV(t) in the hundreds/thousands
#          from a near-zero denominator. Now every denominator is tolerance-gated
#          against the null residual/latent variance and a fully-boundary null
#          block is flagged (null_at_boundary), so degenerate cells are NA.
#   #2 P1  maihda_resolve_strata_formula() detected duplicate intercept-only
#          stratum terms by string-comparing the LHS to "1", so the equivalent
#          spelling (1 | stratum) + (0 + 1 | stratum) bypassed the guard and lme4
#          split the between-stratum variance across 'stratum'/'stratum.1'. Now the
#          random-effect design columns are inspected instead.
#   #3 P1  The lme4 longitudinal paths froze a formula offset at its stored per-row
#          value / mean, so a time-dependent offset such as offset(0.5 * time) did
#          not track the trajectory time (fixed trajectory off by up to 1 link unit;
#          count VPC(t) biased). Now a formula offset is re-evaluated on the grid.
#   #4 P2  estimation_used ignored ml_refit_failed, so a failed refitML() (which
#          keeps REML) was still labelled a pure, correction-free "ML" comparison.
#          Now a skipped OR failed refit is "mixed".
#   #5 P2  estimation = "ML" on a brms fit is a no-op (no refit), yet the basis
#          resolved to "ML" and printed "ML-refit (correction-free...)". brms is a
#          Bayesian posterior, not ML; the basis is now reported as "posterior".

# ---- #1 P1: longitudinal PCV tolerance-gates the between-stratum denominator ----

test_that("maihda_longitudinal_pcv gates a boundary between-stratum denominator", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  a <- suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal"))
  null_m <- a$model; adj_m <- a$model_adjusted
  # Force the null between-stratum growth block onto the boundary (both variances ~0)
  # while the residual scale stays healthy (read from the real fit): every PCV
  # denominator a(t)' Sn a(t) is then boundary-level and MUST be gated to NA rather
  # than divided into an exploded ratio. Deterministic -- it does not rely on lme4
  # converging to a singular fit. (estimation = "fitted" so the models are not refit
  # and the mock's identity comparison holds.)
  testthat::local_mocked_bindings(
    maihda_re_block = function(object, group) {
      if (identical(object, adj_m)) matrix(c(1e-14, 0, 0, 1e-14), 2)
      else matrix(c(1e-13, 0, 0, 1e-13), 2)
    })
  p <- maihda_longitudinal_pcv(null_m, adj_m, estimation = "fitted")
  expect_true(is.na(p$pcv_intercept))
  expect_true(is.na(p$pcv_slope))
  expect_true(all(is.na(p$pcv_t$pcv)))
  # (The old code divided by these ~1e-13 denominators and returned finite, wildly
  # out-of-range ratios instead of NA.)
})

test_that("a no-signal longitudinal null flags null_at_boundary and NAs every PCV", {
  skip_on_cran()
  set.seed(42)
  combos <- expand.grid(gender = c("M", "F"), race = c("A", "B", "C", "D"),
                        edu = c("lo", "mid", "hi"), stringsAsFactors = FALSE)
  npg <- 8L; waves <- 0:3
  persons <- do.call(rbind, lapply(seq_len(nrow(combos)), function(k)
    data.frame(gender = combos$gender[k], race = combos$race[k], edu = combos$edu[k],
               pid = paste0("s", k, "_", seq_len(npg)), stringsAsFactors = FALSE)))
  long <- do.call(rbind, lapply(waves, function(w) { p <- persons; p$wave <- w; p }))
  pre <- stats::setNames(rnorm(nrow(persons), 0, 0.5), persons$pid)
  # A person random intercept + noise but NO stratum trajectory signal drives the null
  # (wave | stratum) growth block to the singularity boundary.
  long$y <- 1 + 0.3 * long$wave + pre[long$pid] + rnorm(nrow(long), 0, 0.8)
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ wave + (1 | gender:race:edu), data = long,
           id = "pid", time = "wave", decomposition = "longitudinal")))
  # lme4's singular-fit optimum is platform-dependent; the gate itself is covered
  # deterministically above, so only assert the end-to-end flag/print when this fit
  # actually reached the boundary (it does on the reference platform).
  skip_if_not(isTRUE(a$pcv$null_at_boundary),
              "null growth block did not reach the boundary on this platform")
  expect_true(is.na(a$pcv$pcv_intercept))
  expect_true(is.na(a$pcv$pcv_slope))
  fin <- a$pcv$pcv_t$pcv[is.finite(a$pcv$pcv_t$pcv)]
  expect_length(fin, 0L)
  expect_output(print(a$pcv), "singularity\\s+boundary")
})

test_that("a healthy longitudinal PCV is unaffected (finite, inside (0, 1))", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  a <- suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal"))
  expect_false(isTRUE(a$pcv$null_at_boundary))
  expect_gt(a$pcv$pcv_intercept, 0); expect_lt(a$pcv$pcv_intercept, 1)
  expect_gt(a$pcv$pcv_slope, 0);     expect_lt(a$pcv$pcv_slope, 1)
  expect_true(all(is.finite(a$pcv$pcv_t$pcv)))
  expect_true(all(a$pcv$pcv_t$pcv >= 0 & a$pcv$pcv_t$pcv <= 1))
})

# ---- #2 P1: equivalent duplicate stratum intercepts are rejected ------------

test_that("fit_maihda() rejects (1 | stratum) + (0 + 1 | stratum)", {
  sd <- make_strata(maihda_sim_data, c("gender", "race"))
  d <- sd$data
  # (0 + 1 | stratum) builds the same constant intercept column as (1 | stratum),
  # so the pair is non-identifiable -- previously accepted by the string guard.
  expect_error(
    suppressMessages(fit_maihda(health_outcome ~ (1 | stratum) + (0 + 1 | stratum), data = d)),
    "non-identifiable")
  # The literal duplicate is still rejected.
  expect_error(
    fit_maihda(health_outcome ~ (1 | stratum) + (1 | stratum), data = d),
    "non-identifiable")
  # A single stratum intercept is still accepted.
  expect_no_error(suppressMessages(fit_maihda(health_outcome ~ (1 | stratum), data = d)))
})

test_that("maihda_re_lhs_is_constant_intercept detects intercept columns by design", {
  d <- data.frame(x = 1:5, time = c(0, 1, 2, 3, 4))
  expect_true(MAIHDA:::maihda_re_lhs_is_constant_intercept(quote(1), d))
  expect_true(MAIHDA:::maihda_re_lhs_is_constant_intercept(quote(0 + 1), d))
  expect_false(MAIHDA:::maihda_re_lhs_is_constant_intercept(quote(time), d))     # slope
  expect_false(MAIHDA:::maihda_re_lhs_is_constant_intercept(quote(1 + time), d)) # int + slope
  expect_false(MAIHDA:::maihda_re_lhs_is_constant_intercept(quote(0 + time), d)) # slope only
})

# ---- #3 P1: time-dependent formula offsets are re-evaluated on the grid -----

test_that("maihda_lme4_formula_offset_at re-evaluates a formula offset on newdata", {
  set.seed(11)
  d <- data.frame(y = rnorm(40), x = rep(1:4, 10), g = factor(rep(1:8, 5)))
  m <- suppressMessages(lme4::lmer(y ~ x + offset(0.5 * x) + (1 | g), data = d))
  off <- MAIHDA:::maihda_lme4_formula_offset_at(m, data.frame(x = c(0, 2, 4)))
  expect_equal(off, c(0, 1, 2))
  # A fit with no formula offset returns NULL.
  m0 <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d))
  expect_null(MAIHDA:::maihda_lme4_formula_offset_at(m0, data.frame(x = c(0, 2, 4))))
})

test_that("maihda_lme4_formula_offset_at falls back to stored values for an offset-only variable", {
  set.seed(12)
  d <- data.frame(y = rnorm(40), x = rep(1:4, 10), g = factor(rep(1:8, 5)),
                  logE = log(runif(40, 1, 5)))
  m <- suppressMessages(lme4::lmer(y ~ x + offset(logE) + (1 | g), data = d))
  # The model frame stores only the derived "offset(logE)" column (logE is not otherwise
  # a predictor), so a grid built from it cannot re-evaluate the term; such a term is
  # necessarily time-invariant and must fall back to the stored per-row values (or their
  # mean for a representative-profile grid) rather than error.
  mf <- stats::model.frame(m)
  expect_false("logE" %in% names(mf))
  expect_equal(MAIHDA:::maihda_lme4_formula_offset_at(m, mf), d$logE)
  expect_equal(
    MAIHDA:::maihda_lme4_formula_offset_at(m, mf[rep(1L, 3), , drop = FALSE],
                                           fallback = "mean"),
    rep(mean(d$logE), 3))
})

test_that("longitudinal fixed trajectory tracks a time-dependent offset", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  m <- suppressMessages(
    fit_maihda(wellbeing ~ wave + offset(0.5 * wave) +
                 (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave"))
  grid <- sort(unique(maihda_long_data$wave))
  traj <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid)
  # Rebuild the expected fixed trajectory with the offset EVALUATED at each grid time.
  nd <- maihda_long_data[rep(1L, length(grid)), , drop = FALSE]
  nd$wave <- grid
  ct <- MAIHDA:::maihda_lng_time_term(m$longitudinal_info)
  if (!identical(ct, "wave")) {
    nd[[ct]] <- grid - MAIHDA:::maihda_lng_time_center(m$longitudinal_info)
  }
  expected <- MAIHDA:::maihda_lme4_fixed_link(m$model, nd, offset = NULL) + 0.5 * grid
  expect_equal(traj, expected, tolerance = 1e-8)
  # It must NOT equal the old frozen-at-mean behaviour (offset held at mean(0.5*wave)).
  frozen <- MAIHDA:::maihda_lme4_fixed_link(m$model, nd, offset = NULL) +
    mean(0.5 * maihda_long_data$wave)
  expect_gt(max(abs(traj - frozen)), 0.1)
})

test_that("a no-offset longitudinal trajectory is unchanged", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  m <- suppressMessages(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave"))
  grid <- sort(unique(maihda_long_data$wave))
  traj <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid)
  ref <- as.numeric(stats::predict(
    m$model, re.form = NA,
    newdata = transform(maihda_long_data[rep(1L, length(grid)), ], wave = grid)))
  expect_equal(traj, ref, tolerance = 1e-8)
})

# ---- #4 P2: a failed ML refit is reported as "mixed", not pure ML -----------

test_that("maihda_pcv_estimation_used folds in failed refits and the engine", {
  # skipped OR failed at a boundary -> mixed under ML; fitted is never mixed.
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("ML", TRUE), "mixed")
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("ML", FALSE), "ML")
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("fitted", TRUE), "fitted")
  # the "mixed" label no longer asserts a pure, correction-free comparison.
  expect_false(grepl("correction-free", MAIHDA:::maihda_pcv_basis_label("mixed")))
})

test_that("calculate_pcv() reports a failed ML refit as a mixed basis", {
  sd <- make_strata(maihda_sim_data, c("gender", "race"))
  d <- sd$data
  m1 <- suppressMessages(fit_maihda(health_outcome ~ (1 | stratum), data = d))
  m2 <- suppressMessages(fit_maihda(health_outcome ~ age + (1 | stratum), data = d))
  # Force refitML() to fail: the ML refit is warranted but cannot complete, so both
  # models stay on REML -- a mixed basis, not the pure ML the label used to claim. The
  # real refitML is restored via tryCatch(finally = ) immediately after the call (no
  # helper package: an undeclared withr:: here is an R CMD check unstated-test-
  # dependency WARNING).
  orig <- lme4::refitML
  assignInNamespace("refitML", function(object, ...) stop("forced refitML failure"),
                    ns = "lme4")
  res <- tryCatch(
    suppressWarnings(calculate_pcv(m1, m2, estimation = "ML")),
    finally = assignInNamespace("refitML", orig, ns = "lme4"))
  expect_true(isTRUE(res$ml_refit_failed))
  expect_identical(res$estimation_used, "mixed")
  expect_false(grepl("correction-free", MAIHDA:::maihda_pcv_basis_label(res$estimation_used)))
})

# ---- #5 P2: brms PCV is reported as a posterior basis, not an ML-refit -------

test_that("a brms PCV basis is 'posterior', never an ML-refit", {
  # estimation = "ML" is a no-op for brms (a Bayesian posterior, not ML), so the basis
  # is the as-fitted posterior scale regardless of the requested estimation.
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("ML", FALSE, engine = "brms"),
                   "posterior")
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("fitted", FALSE, engine = "brms"),
                   "posterior")
  # lme4 is unaffected (still a genuine ML refit basis).
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("ML", FALSE, engine = "lme4"), "ML")
  lbl <- MAIHDA:::maihda_pcv_basis_label("posterior")
  expect_true(grepl("posterior", lbl))
  expect_false(grepl("ML-refit", lbl))
})
