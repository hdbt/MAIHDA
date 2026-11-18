# Audit 2026-11-17, one finding from the independent verifier of the growth-grad-warning pass
# (test-audit-2026-11-16-growth-grad-warning.R): a growth fit lme4 flags was kept beside a
# singular fit of the same criterion.
#
# fit_maihda() fits a Gaussian growth model from lme4's start and from further starts and
# keeps lmer's own fit unless another is lower by more than 1e-4 in the criterion. Where the
# optimum lies on a variance bound, one start ends on the bound, where lme4's test of a
# singular fit sees it, and another a little short of it, at an entry of 1e-4 to 4e-2 of the
# residual SD, where it does not. The two differ little in what they estimate (VPC(t) within
# 0.002), but the fit short of the bound was reported as if the variance there were a small
# number and not zero: with no other bound reached, without the caveat of a singular fit and
# with lme4's "Model failed to converge with max|grad|" where it had stopped more than about
# 2e-6 above the optimum; and on any bound with a PCV that divides by that variance (a
# PCV_slope of 1.000 or -62.75, a PCV_intercept of 0.943 or -15,600,000, where the fit on
# the bound gives NA). On simulated designs of 6 to 48 strata the fit of the lowest criterion
# was short of a bound that another fit within 1e-4 was on in 4.8% of linear growth fits with
# a stratum slope and in 18.1% of quadratic ones without; a multi-start search of the
# criterion found the optimum on that bound in 349 of the 354 fits it was run on.
#
# lme4 runs no convergence check on a singular fit. The package's own measure of a fit's
# distance from a stationary point of its criterion (maihda_lmer_stationary_distance()) can
# be taken on one: on 344 such fits it was 0.85 to 1.04 times the square root of the fit's
# criterion above the searched optimum, and on each of the 5 singular fits on a bound the
# optimum was not on it could not be taken (a Hessian that is not positive definite).
#
# FIX: maihda_fit_lmer_growth() chooses as before (maihda_pick_lowest()) and then, when
# another fit within 1e-4 of the chosen one's criterion is on every variance bound the chosen
# fit is on and on at least one more (maihda_lme4_bounds_reached()), and stands at an optimum
# (maihda_lme4_singular_at_optimum(): lme4 calls it singular, its optimiser reported no
# failure, and that measure puts it within 1e-4 of a stationary point), keeps the one of
# those on the most bounds, of several the lowest (maihda_pick_singular()), whether or not
# lme4 flagged the fit chosen. The rule of the 2026-11-16 pass (maihda_pick_accepted()) then
# acts on what is left.
#
# This file draws no random numbers.

fbs17_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression with the text of every warning and message it raised, in order.
fbs17_heard <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(expr,
    warning = function(cnd) {
      seen <<- c(seen, paste("warning:", conditionMessage(cnd)))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) {
      seen <<- c(seen, paste("message:", conditionMessage(cnd)))
      invokeRestart("muffleMessage")
    })
  list(value = value, seen = seen)
}

# lme4's warning for a gradient above its tolerance and its note of a singular fit.
fbs17_flag <- "Model failed to converge with max|grad| = 0.00372085 (tol = 0.002, component 1)"
fbs17_note <- "boundary (singular) fit: see help('isSingular')"

# What print() shows, as one line.
fbs17_printed <- function(x) gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))

# A growth panel made without random numbers: 6 strata (2 x 3) of 4 persons at 5 waves. The
# persons' slopes take the same four values in every stratum, so the strata do not differ in
# slope and the REML optimum has the stratum block at rank one, its slope entry on the bound
# ("rank one"). Without the strata's levels (`levels = FALSE`) they do not differ in level
# either, and the optimum has the whole stratum block at zero, both its entries on their
# bounds ("flat").
fbs17_data <- function(levels = TRUE) {
  d <- expand.grid(wave = 0:4, person = 1:4, b = 1:3, a = 1:2)
  d$id <- factor(with(d, (a - 1) * 12 + (b - 1) * 4 + person))
  s <- with(d, (a - 1) * 3 + b)
  level_s <- if (levels) c(0.6, -0.5, 0.2, -0.3, 0.7, -0.7)[s] else 0
  level_p <- c(-1.1, 0.4, 0.9, -0.2)[d$person] * c(1, -1, 1, -1, 1, -1)[s]
  slope_p <- c(-0.45, 0.15, 0.45, -0.15)[d$person]
  noise <- sin(1.7 * as.integer(d$id) + 2.3 * d$wave) * 0.8 +
    cos(0.9 * as.integer(d$id) * (d$wave + 1)) * 0.5
  d$y <- 10 + level_s + level_p + (0.5 + slope_p) * d$wave + noise
  d$a <- factor(d$a, labels = c("a1", "a2"))
  d$b <- factor(d$b, labels = c("b1", "b2", "b3"))
  d[, c("id", "wave", "a", "b", "y")]
}

# The fit of the lmer call `call` (evaluated in `env`) at exactly `theta`, made by lme4's
# modular functions without an optimiser, with lme4's convergence record `lme4conv` and the
# optimiser's code `code` written by hand.
fbs17_placed <- function(call, env, theta, lme4conv = list(), code = 0L) {
  parts <- call
  parts[[1L]] <- quote(lme4::lFormula)
  parts$start <- NULL
  lmod <- eval(parts, env)
  devfun <- do.call(lme4::mkLmerDevfun, lmod)
  opt <- list(par = theta, fval = devfun(theta), feval = 1L, conv = code, message = NULL)
  lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr, mc = call,
                 lme4conv = lme4conv)
}

# The REML optimum of that call on the bounds it lies on: theta is (person: intercept,
# covariance, slope; stratum: the same), the entries that are zero in `start` are held at
# zero and the others found by lme4's own optimiser, run twice at tight tolerances from
# `start`. That the point is an optimum is checked where it is used.
fbs17_optimum <- function(call, env, start) {
  parts <- call
  parts[[1L]] <- quote(lme4::lFormula)
  parts$start <- NULL
  lmod <- eval(parts, env)
  devfun <- do.call(lme4::mkLmerDevfun, lmod)
  free <- which(start != 0)
  on_bounds <- function(p) devfun(replace(numeric(length(start)), free, p))
  lower <- environment(devfun)$lower[free]
  tight <- list(xtol_abs = 1e-10, ftol_abs = 1e-13, maxeval = 5000)
  end <- lme4::nloptwrap(start[free], on_bounds, lower = lower, upper = rep(Inf, length(free)),
                         control = tight)
  end <- lme4::nloptwrap(end$par, on_bounds, lower = lower, upper = rep(Inf, length(free)),
                         control = tight)
  replace(numeric(length(start)), free, end$par)
}

# Where the search of each panel's optimum starts: beside it, with the entries on a bound,
# and for the flat panel the stratum block's covariance entry, at zero.
fbs17_start <- list("rank one" = c(0.59, 0.28, 0.30, 0.31, 0.014, 0),
                    flat = c(0.41, 0.36, 0.20, 0, 0, 0))

# A value / conditions pair of maihda_eval_holding_conditions() for one kind of fit of a
# panel, placed at a theta made from the panel's optimum, with lme4's record written by hand.
# The rank-one panel: "singular" (the optimum), "short flagged" and "short accepted" (the
# stratum slope entry at 5e-4 instead of 0, a fit lme4 does not call singular, with lme4's
# gradient flag or an empty record), "failed" (the optimum under an optimiser's code of 5) and
# "elsewhere" (the stratum slope at zero, a singular fit that is not the optimum). The flat
# panel: "both" (the optimum, on the stratum intercept's bound and the slope's), "slope only"
# (the stratum intercept entry at 5e-4: on the slope's bound, short of the intercept's) and
# "neither flagged" / "neither accepted" (both entries at 5e-4).
fbs17_pair <- function(kind, call, env, optimum) {
  noted <- list(simpleMessage(paste0(fbs17_note, "\n")))
  on_bound <- function(theta, code = 0L) {
    list(value = fbs17_placed(call, env, theta, list(messages = fbs17_note), code = code),
         conditions = noted)
  }
  off_bound <- function(theta, flagged) {
    if (flagged) {
      list(value = fbs17_placed(call, env, theta, list(code = -1L, messages = fbs17_flag)),
           conditions = list(simpleWarning(fbs17_flag)))
    } else {
      list(value = fbs17_placed(call, env, theta), conditions = list())
    }
  }
  switch(kind,
    singular = on_bound(optimum),
    "short flagged" = off_bound(replace(optimum, 6L, 5e-4), TRUE),
    "short accepted" = off_bound(replace(optimum, 6L, 5e-4), FALSE),
    failed = on_bound(optimum, code = 5L),
    elsewhere = on_bound(replace(optimum, 5:6, 0)),
    both = on_bound(optimum),
    "slope only" = on_bound(replace(optimum, 4L, 5e-4)),
    "neither flagged" = off_bound(replace(optimum, c(4L, 6L), 5e-4), TRUE),
    "neither accepted" = off_bound(replace(optimum, c(4L, 6L), 5e-4), FALSE),
    stop("unknown kind of fit: ", kind))
}

# The optimum of the model fit_maihda() builds for each panel, found once for the file.
fbs17_found <- new.env()

# fit_maihda() on a panel with every candidate fit placed by hand (fbs17_pair()): lmer's own
# fit is `first`, and each fit from a further start one of `restarts`, in the order the
# fitter makes them (recycled). Which real fit ends on a bound, and which one lme4 flags,
# differs between machines, so none is left to an optimiser.
fbs17_fit_with <- function(first, restarts, panel = "rank one") {
  made <- NULL
  heard <- fbs17_heard(with_mocked_bindings(
    fit_maihda(y ~ wave + (1 | a:b), data = fbs17_data(levels = panel == "rank one"),
               id = "id", time = "wave"),
    maihda_eval_holding_conditions = function(call, env) {
      if (is.null(made)) {
        if (is.null(fbs17_found[[panel]])) {
          fbs17_found[[panel]] <- fbs17_optimum(call, env, fbs17_start[[panel]])
        }
        own <- fbs17_pair(first, call, env, fbs17_found[[panel]])
        made <<- list(own = own$value, restarts = 0L)
        return(own)
      }
      made$restarts <<- made$restarts + 1L
      fbs17_pair(restarts[(made$restarts - 1L) %% length(restarts) + 1L], call, env,
                 fbs17_found[[panel]])
    }))
  c(heard, list(own = made$own, optimum = fbs17_found[[panel]]))
}

# ---- the choice ---------------------------------------------------------------------------

test_that("a fit gives way to the fit at an optimum on the most variance bounds within 1e-4", {
  pick <- maihda_pick_singular
  y <- TRUE
  n <- FALSE
  # No fit beside the choice that may take its place: it stays, whatever bounds the others
  # are on.
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(0L, 1L, 1L), c(n, n, n), 1L), 1L)
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(1L, 1L, 1L), c(n, n, n), 1L), 1L)
  # Fits that may, within 1e-4 above or below: the one on the most bounds ...
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(0L, 1L, 2L), c(n, y, y), 1L), 3L)
  expect_identical(pick(c(10, 10 + 5e-5, 10 - 5e-5), c(0L, 2L, 1L), c(n, y, y), 1L), 2L)
  expect_identical(pick(c(10, 10 + 5e-5, 10 + 2e-5), c(1L, 3L, 2L), c(n, y, y), 1L), 2L)
  # ... and of several on as many, the lowest.
  expect_identical(pick(c(10, 10 + 5e-5, 10 - 5e-5), c(0L, 1L, 1L), c(n, y, y), 1L), 3L)
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(0L, 1L, 1L), c(n, y, y), 1L), 2L)
  expect_identical(pick(c(10, 10 + 5e-5, 10 + 2e-5), c(1L, 2L, 2L), c(n, y, y), 1L), 3L)
  expect_identical(pick(c(10, 10 + 5e-5), c(0L, 1L), c(n, y), 1L), 2L)
  # A fit on more bounds that may not take the choice's place does not, lower or not; one
  # beside it that may does.
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(0L, 2L, 1L), c(n, n, y), 1L), 3L)
  # A fit more than 1e-4 above is another optimum: the choice stays.
  expect_identical(pick(c(10, 10 + 2e-4, 11), c(0L, 1L, 1L), c(n, y, y), 1L), 1L)
  expect_identical(pick(c(10, 10 + 2e-4, 10 + 5e-5), c(0L, 2L, 1L), c(n, y, y), 1L), 3L)
  # A restart chosen for a criterion lower by more than 1e-4: a fit within 1e-4 of it takes
  # its place, and lmer's own fit, far above it, does not.
  expect_identical(pick(c(10, 9, 9 + 4e-5), c(0L, 0L, 1L), c(n, n, y), 2L), 3L)
  expect_identical(pick(c(10, 9, 9.5), c(1L, 0L, 1L), c(y, n, y), 2L), 2L)
  # A restart that failed has no criterion and is passed over.
  expect_identical(pick(c(10, NA, 10 + 5e-5), c(0L, 1L, 1L), c(n, y, y), 1L), 3L)
  expect_identical(pick(c(10, NA, NA), c(0L, 1L, 1L), c(n, y, y), 1L), 1L)
  # The tolerance is the caller's.
  expect_identical(pick(c(10, 10.05), c(0L, 1L), c(n, y), 1L, tol = 0.1), 2L)
  expect_identical(pick(c(10, 10.05), c(0L, 1L), c(n, y), 1L, tol = 0.01), 1L)
  expect_identical(formals(maihda_pick_singular)$tol, 1e-4)
  # A choice without a criterion, or no choice: returned as given.
  expect_identical(pick(c(NA, 12, 11), c(0L, 1L, 1L), c(n, y, y), 1L), 1L)
  expect_identical(pick(numeric(0), integer(0), logical(0), NA_integer_), NA_integer_)
  # A fit of which it is not known whether it may take the place does not; a count of bounds
  # that is not known counts as none.
  expect_identical(pick(c(10, 10 + 5e-5, 10 + 6e-5), c(0L, 1L, 1L), c(n, NA, y), 1L), 3L)
  expect_identical(pick(c(10, 10 + 5e-5), c(0L, 1L), c(n, NA), 1L), 1L)
  expect_identical(pick(c(10, 10 + 5e-5, 10 + 6e-5), c(0L, NA, 1L), c(n, y, y), 1L), 3L)
  expect_identical(pick(c(10, 10 + 5e-5), c(0L, NA), c(n, y), 1L), 2L)
})

test_that("fits on more bounds than the choice are measured, and only those within 1e-4 of it", {
  measured <- integer(0)
  tols <- numeric(0)
  local_mocked_bindings(maihda_lme4_singular_at_optimum = function(model, tol = 1e-4) {
    measured <<- c(measured, model$k)
    tols <<- c(tols, tol)
    isTRUE(model$opt)
  })
  fits <- lapply(1:4, function(k) list(value = list(k = k, opt = k != 3L), conditions = list()))
  beside <- function(crit, bounds, kept, ...) {
    measured <<- integer(0)
    tols <<- numeric(0)
    maihda_singular_fits_beside(fits, crit, bounds, kept, ...)
  }
  none <- c(FALSE, FALSE)
  first <- c(TRUE, FALSE)
  second <- c(FALSE, TRUE)
  both <- c(TRUE, TRUE)
  crit <- c(10, 10 + 5e-5, 10 - 5e-5, 10 + 2e-4)
  # The choice is on no bound: the fits within 1e-4 of it that are on one are measured, each
  # once, with the tolerance in force; one further off, or on none, is not.
  expect_identical(beside(crit, list(none, first, second, both), 1L), c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(measured, 2:3)
  expect_identical(tols, c(1e-4, 1e-4))
  expect_identical(beside(crit, list(none, first, none, both), 1L), c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(measured, 2L)
  expect_identical(beside(crit, list(none, first, second, both), 1L, tol = 1e-3),
                   c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(measured, 2:4)
  expect_identical(tols, rep(1e-3, 3L))
  # The choice is on a bound: a fit on that bound and another is measured; one on the other
  # bound alone, or on the same bound alone, is not.
  expect_identical(beside(crit, list(first, both, second, first), 1L), c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(measured, 2L)
  expect_identical(beside(crit, list(first, first, second, both), 1L), rep(FALSE, 4L))
  expect_identical(measured, integer(0))
  expect_identical(beside(crit, list(second, first, both, none), 1L), rep(FALSE, 4L))
  expect_identical(measured, 3L)
  # The choice is on every bound: nothing is measured.
  expect_identical(beside(crit, list(both, both, both, both), 1L), rep(FALSE, 4L))
  expect_identical(measured, integer(0))
  # The bounds of the choice are not known, it has no criterion, or there is none: nothing is
  # measured.
  expect_identical(beside(crit, list(NULL, first, both, both), 1L), rep(FALSE, 4L))
  expect_identical(beside(crit, list(c(NA, FALSE), first, both, both), 1L), rep(FALSE, 4L))
  expect_identical(beside(c(NA, 10, 10, 10), list(none, first, both, both), 1L), rep(FALSE, 4L))
  expect_identical(beside(crit, list(none, first, both, both), NA_integer_), rep(FALSE, 4L))
  expect_identical(measured, integer(0))
  # A fit whose bounds are not known, or do not match the choice's, is passed over.
  expect_identical(beside(crit, list(none, NULL, c(TRUE, NA), TRUE), 1L, tol = 1e-3),
                   rep(FALSE, 4L))
  expect_identical(measured, integer(0))
  # A fit that failed (NULL, no criterion) is passed over.
  measured <- integer(0)
  expect_identical(maihda_singular_fits_beside(list(fits[[1]], NULL, fits[[2]]),
                                               c(10, NA, 10 + 5e-5), list(none, NULL, first), 1L),
                   c(FALSE, FALSE, TRUE))
  expect_identical(measured, 2L)
})

# ---- lme4's verdict and the package's own check, read off the fit ---------------------------

test_that("a singular fit at its optimum is told from one short of the bound and from one elsewhere", {
  skip_on_cran()
  d <- fbs17_data()
  call <- quote(lme4::lmer(formula = y ~ wave + (wave | id) + (wave | stratum), data = d))
  env <- environment()
  d$stratum <- interaction(d$a, d$b)
  optimum <- fbs17_optimum(call, env, fbs17_start[["rank one"]])
  # The stratum block there has rank one: the slope entry on its bound beside an entry that is
  # not zero.
  expect_identical(optimum[6], 0)
  expect_gt(abs(optimum[5]), 0.005)
  at_optimum <- fbs17_placed(call, env, optimum, list(messages = fbs17_note))
  short <- fbs17_placed(call, env, replace(optimum, 6L, 5e-4))
  elsewhere <- fbs17_placed(call, env, replace(optimum, 5:6, 0), list(messages = fbs17_note))
  crit <- vapply(list(at_optimum, short, elsewhere), lme4::REMLcrit, numeric(1))
  # That point on the bound is an optimum of the panel's criterion, shown without asking an
  # optimiser where it stops: the criterion rises off the bound on either side of it (below),
  # and the package's measure puts the point on a stationary point (further below). Left to
  # itself lme4's optimiser ends up to 7e-5 short of the bound in that entry, by the last bits
  # of the data.
  expect_gt(lme4::REMLcrit(fbs17_placed(call, env, replace(optimum, 6L, -5e-4))), crit[1])
  # The fit short of the bound is within 1e-4 of the optimum, and more than lme4's gradient
  # tolerance allows above it; the singular fit elsewhere is far above.
  expect_gt(crit[2] - crit[1], 2e-6)
  expect_lt(crit[2] - crit[1], 1e-4)
  expect_gt(crit[3] - crit[1], 0.01)
  # lme4 calls the two fits on the bound singular and not the one short of it; the bound is
  # the fourth of the model's four (persons: intercept, slope; strata: intercept, slope).
  expect_true(maihda_lme4_calls_singular(at_optimum))
  expect_true(maihda_lme4_calls_singular(elsewhere))
  expect_false(maihda_lme4_calls_singular(short))
  expect_identical(maihda_lme4_bounds_reached(at_optimum), c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(maihda_lme4_bounds_reached(elsewhere), c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(maihda_lme4_bounds_reached(short), c(FALSE, FALSE, FALSE, FALSE))
  # The package's measure puts the fit at the optimum on a stationary point, and the others
  # at the square root of their criterion above it.
  expect_lt(maihda_lmer_stationary_distance(at_optimum), 1e-3)
  expect_lt(abs(maihda_lmer_stationary_distance(short) / sqrt(crit[2] - crit[1]) - 1), 0.05)
  expect_lt(abs(maihda_lmer_stationary_distance(elsewhere) / sqrt(crit[3] - crit[1]) - 1), 0.05)
  # So only the fit at the optimum is a singular fit at an optimum.
  expect_true(maihda_lme4_singular_at_optimum(at_optimum))
  expect_false(maihda_lme4_singular_at_optimum(elsewhere))
  expect_false(maihda_lme4_singular_at_optimum(short))
  # Not under an optimiser's code other than 0.
  expect_false(maihda_lme4_singular_at_optimum(
    fbs17_placed(call, env, optimum, list(messages = fbs17_note), code = 5L)))
  # A singular fit where the criterion still falls off the bound: the slope variance of the
  # sleep study held at zero. The measure cannot be taken there or is far from zero.
  sleep <- lme4::sleepstudy
  held <- fbs17_placed(quote(lme4::lmer(formula = Reaction ~ Days + (Days | Subject),
                                        data = sleep)), env, c(1.5, 0, 0))
  expect_true(maihda_lme4_calls_singular(held))
  expect_identical(maihda_lme4_bounds_reached(held), c(FALSE, TRUE))
  expect_false(maihda_lme4_singular_at_optimum(held))
})

test_that("a fit on one variance bound is told from the fit at the optimum on that bound and another", {
  skip_on_cran()
  d <- fbs17_data(levels = FALSE)
  call <- quote(lme4::lmer(formula = y ~ wave + (wave | id) + (wave | stratum), data = d))
  env <- environment()
  d$stratum <- interaction(d$a, d$b)
  optimum <- fbs17_optimum(call, env, fbs17_start[["flat"]])
  expect_identical(optimum[4:6], c(0, 0, 0))
  pairs <- lapply(c("slope only", "both", "neither accepted"), fbs17_pair, call = call,
                  env = env, optimum = optimum)
  fits <- lapply(pairs, `[[`, "value")
  crit <- vapply(fits, lme4::REMLcrit, numeric(1))
  # The fit with the whole stratum block at zero is the lowest; the fit short of the stratum
  # intercept's bound and the fit short of both are within 1e-4 of it.
  expect_gt(crit[1] - crit[2], 2e-6)
  expect_gt(crit[3] - crit[1], 2e-6)
  expect_lt(crit[3] - crit[2], 1e-4)
  # lme4 calls the first two singular alike; they differ in the bounds they are on.
  expect_identical(vapply(fits, maihda_lme4_calls_singular, logical(1)), c(TRUE, TRUE, FALSE))
  bounds <- lapply(fits, maihda_lme4_bounds_reached)
  expect_identical(bounds, list(c(FALSE, FALSE, FALSE, TRUE), c(FALSE, FALSE, TRUE, TRUE),
                                c(FALSE, FALSE, FALSE, FALSE)))
  # The fit on both bounds is on a stationary point: the criterion rises off either bound.
  expect_lt(maihda_lmer_stationary_distance(fits[[2]]), 1e-3)
  expect_true(maihda_lme4_singular_at_optimum(fits[[2]]))
  # From the fit on the slope's bound alone, the fit on both is the one beside it on more
  # bounds, and the one to keep; the fit on neither is not.
  beside <- maihda_singular_fits_beside(pairs, crit, bounds, 1L)
  expect_identical(beside, c(FALSE, TRUE, FALSE))
  expect_identical(maihda_pick_singular(crit, c(1L, 2L, 0L), beside, 1L), 2L)
  # From the fit on neither bound, both others are on more bounds and within 1e-4 of an
  # optimum; the one on the most bounds is kept.
  beside <- maihda_singular_fits_beside(pairs, crit, bounds, 3L)
  expect_identical(beside, c(TRUE, TRUE, FALSE))
  expect_identical(maihda_pick_singular(crit, c(1L, 2L, 0L), beside, 3L), 2L)
  # From the fit on both, nothing is on more.
  expect_identical(maihda_singular_fits_beside(pairs, crit, bounds, 2L), c(FALSE, FALSE, FALSE))
})

test_that("a singular fit is at an optimum within the caller's tolerance of a stationary point", {
  skip_on_cran()
  flat <- data.frame(y = rep(c(1, 2, 3, 4), 10), g = factor(rep(1:10, each = 4)))
  singular <- fbs17_quiet(lme4::lmer(y ~ 1 + (1 | g), data = flat))
  expect_true(maihda_lme4_calls_singular(singular))
  expect_identical(maihda_lme4_bounds_reached(singular), TRUE)
  distance <- NA_real_
  local_mocked_bindings(maihda_lmer_stationary_distance = function(model) distance)
  # Within 1e-4 of a stationary point is a distance below 0.01; within 4e-4, below 0.02.
  distance <- 0.009
  expect_true(maihda_lme4_singular_at_optimum(singular))
  distance <- 0.011
  expect_false(maihda_lme4_singular_at_optimum(singular))
  expect_true(maihda_lme4_singular_at_optimum(singular, tol = 4e-4))
  expect_identical(formals(maihda_lme4_singular_at_optimum)$tol, 1e-4)
  # A measure that cannot be taken vouches for nothing.
  distance <- NA_real_
  expect_false(maihda_lme4_singular_at_optimum(singular))
  # A fit lme4 does not call singular is not asked.
  distance <- 0
  sleep <- fbs17_quiet(lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy))
  expect_false(maihda_lme4_singular_at_optimum(sleep))
  expect_identical(maihda_lme4_bounds_reached(sleep), c(FALSE, FALSE))
  expect_true(maihda_lme4_singular_at_optimum(singular))
})

test_that("anything that is not an lmer fit is on no known bound and at no optimum", {
  for (x in list(NULL, list(crit = 1), stats::lm(dist ~ speed, data = datasets::cars))) {
    expect_false(maihda_lme4_calls_singular(x))
    expect_null(maihda_lme4_bounds_reached(x))
    expect_false(maihda_lme4_singular_at_optimum(x))
  }
})

# ---- the fitter -----------------------------------------------------------------------------

test_that("the growth fitter keeps the fit at an optimum on more bounds in place of a fit short of one", {
  env <- new.env()
  env$measured <- integer(0)
  # lmer's own fit and two restarts as the rows of `candidates`: criterion, lme4's flag, its
  # acceptance, the bounds the fit is on ("00" none of two, "01" the second, "11" both) and
  # whether a singular fit is at an optimum. A fit lme4 flags warns as lme4 does, and every
  # fit says which it is.
  env$fake_fit <- function(start = NULL) {
    k <- if (is.null(start)) 1L else as.integer(start[["k"]])
    row <- env$candidates[k, ]
    if (row$flag) warning(sprintf("lme4 flags fit %d", k))
    message(sprintf("fit %d", k))
    list(crit = row$crit, flag = row$flag, ok = row$ok, on = row$on, opt = row$opt, k = k)
  }
  local_mocked_bindings(
    maihda_growth_restart_starts = function(theta, cnms, id, scale = 0.5) list(c(k = 2), c(k = 3)),
    maihda_growth_centred_start = function(fit_call, fit_env, time, time_degree) list(c(k = 2)),
    maihda_growth_diagonal_start = function(theta, cnms, id, scale = 0.5) list(c(k = 3)),
    maihda_lmer_criterion = function(model) model$crit,
    maihda_lme4_check_flagged = function(model) isTRUE(model$flag),
    maihda_lme4_check_accepted = function(model) isTRUE(model$ok),
    maihda_lme4_bounds_reached = function(model) {
      if (is.null(model$on)) NULL else strsplit(model$on, "")[[1]] == "1"
    },
    maihda_lme4_singular_at_optimum = function(model, tol = 1e-4) {
      env$measured <- c(env$measured, model$k)
      env$tol <- tol
      isTRUE(model$opt)
    }
  )
  kept <- function(crit, flag, on, opt = on != "00", ok = !flag & on == "00", ...) {
    env$candidates <- data.frame(crit = crit, flag = flag, ok = ok, on = on, opt = opt,
                                 stringsAsFactors = FALSE)
    env$measured <- integer(0)
    fbs17_heard(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid", ...))
  }
  y <- TRUE
  n <- FALSE
  # lmer's own fit flagged and on no bound, a restart on one, at an optimum, within 1e-4: that
  # restart, and only what it raised; lme4's warning about lmer's fit is not passed on.
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(y, y, n), c("00", "00", "01"))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  expect_identical(env$measured, 3L)
  # Both restarts on the same bound and at an optimum: the lower one.
  h <- kept(c(10, 10 + 3e-5, 10 - 2e-5), c(y, n, n), c("00", "01", "01"))
  expect_identical(h$value$k, 3L)
  h <- kept(c(10, 10 + 3e-5, 10 + 6e-5), c(y, n, n), c("00", "01", "01"))
  expect_identical(h$value$k, 2L)
  expect_identical(env$measured, 2:3)
  # One restart on one bound and one on two, both at an optimum: the one on two, lower or not.
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(y, n, n), c("00", "01", "11"))
  expect_identical(h$value$k, 3L)
  h <- kept(c(10, 10 + 3e-5, 10 - 2e-5), c(n, n, n), c("00", "11", "10"))
  expect_identical(h$value$k, 2L)
  # lmer's own fit accepted by lme4: the restart on the bound all the same. lme4's flag is not
  # consulted.
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(n, n, n), c("00", "00", "01"))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  # lmer's own fit on one bound, as lme4's singular fits are, and a restart on that bound and
  # another, at an optimum: the restart.
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(n, n, n), c("01", "00", "11"))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  expect_identical(env$measured, 3L)
  # A restart on another bound alone, or on the same one alone, takes no such fit's place and
  # is not measured.
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(n, n, n), c("01", "10", "01"))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, "message: fit 1\n")
  expect_identical(env$measured, integer(0))
  # A restart lme4 accepted and one on a bound, at an optimum, beside a flagged fit: the one
  # on the bound, above the accepted one or below it.
  h <- kept(c(10, 10 - 5e-5, 10 + 3e-5), c(y, n, n), c("00", "00", "01"))
  expect_identical(h$value$k, 3L)
  h <- kept(c(10, 10 + 3e-5, 10 - 5e-5), c(y, n, n), c("00", "00", "01"))
  expect_identical(h$value$k, 3L)
  # So too when those two lie more than 1e-4 apart, each within 1e-4 of the flagged fit: the
  # fit on the bound is looked for beside the fit of the lowest criterion, before the accepted
  # one.
  h <- kept(c(10, 10 - 8e-5, 10 + 8e-5), c(y, n, n), c("00", "00", "01"))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  # The restart on the bound is not at an optimum: the fit lme4 accepted beside the flagged
  # one, as before; with none, lmer's own fit and lme4's warning.
  h <- kept(c(10, 10 - 5e-5, 10 + 3e-5), c(y, n, n), c("00", "00", "01"), opt = c(n, n, n))
  expect_identical(h$value$k, 2L)
  expect_identical(h$seen, "message: fit 2\n")
  expect_identical(env$measured, 3L)
  h <- kept(c(10, 10 - 5e-5, 10 + 3e-5), c(y, y, n), c("00", "00", "01"), opt = c(n, n, n))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, c("warning: lme4 flags fit 1", "message: fit 1\n"))
  # ... and a fit lme4 accepted stays beside it.
  h <- kept(c(10, 10 - 5e-5, 10 + 3e-5), c(n, n, n), c("00", "00", "01"), opt = c(n, n, n))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, "message: fit 1\n")
  # The fit on the bound is more than 1e-4 above: another optimum. It is not measured, and
  # lmer's own fit stays.
  h <- kept(c(10, 10 + 1e-3, 12), c(n, n, n), c("00", "01", "01"))
  expect_identical(h$value$k, 1L)
  expect_identical(env$measured, integer(0))
  # lmer's own fit on every bound: it stays beside fits lme4 accepted, and nothing is
  # measured.
  h <- kept(c(10, 10 - 5e-5, 10 + 3e-5), c(n, n, n), c("11", "00", "00"))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, "message: fit 1\n")
  expect_identical(env$measured, integer(0))
  h <- kept(c(10, 10 - 5e-5, 10 + 3e-5), c(n, n, n), c("11", "11", "01"))
  expect_identical(h$value$k, 1L)
  # A restart lower by more than 1e-4 is chosen; on no bound, it gives way to the other
  # restart within 1e-4 of it that is on one, and not to lmer's own fit far above.
  h <- kept(c(10, 9, 9 + 4e-5), c(n, y, n), c("01", "00", "01"))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  expect_identical(env$measured, 3L)
  h <- kept(c(10, 9, 9.5), c(n, y, n), c("01", "00", "01"))
  expect_identical(h$value$k, 2L)
  expect_identical(h$seen, c("warning: lme4 flags fit 2", "message: fit 2\n"))
  # The fitter's own tolerance is the rule's, and the measure's.
  h <- kept(c(10, 10 + 5e-3, 12), c(n, n, n), c("00", "01", "00"), tol = 0.01)
  expect_identical(h$value$k, 2L)
  expect_identical(env$tol, 0.01)
  h <- kept(c(10, 10 + 5e-3, 12), c(n, n, n), c("00", "01", "00"))
  expect_identical(h$value$k, 1L)
  # Quadratic growth takes the same rule: the centred-time restart beside lmer's fit, and the
  # diagonal start, made when those two differ by more than 0.01.
  h <- kept(c(10, 10 + 3e-5, NA), c(y, n, n), c("00", "01", "00"), time = "wave", time_degree = 2)
  expect_identical(h$value$k, 2L)
  expect_identical(h$seen, "message: fit 2\n")
  h <- kept(c(10, 10.5, 10 - 2e-5), c(n, n, n), c("01", "00", "11"), time = "wave", time_degree = 2)
  expect_identical(h$value$k, 3L)
  expect_identical(env$tol, 1e-4)
  # A fit handed in as the first one, as a bootstrap draw's refit is, is treated as lmer's own.
  env$candidates <- data.frame(crit = c(10, 10 + 3e-5, 10 - 2e-5), flag = c(y, n, n),
                               ok = c(n, n, n), on = c("00", "01", "01"), opt = c(n, y, y),
                               stringsAsFactors = FALSE)
  first <- list(value = list(crit = 10, flag = TRUE, ok = FALSE, on = "00", opt = FALSE, k = 1L),
                conditions = list(simpleWarning("lme4 flags the fit handed in")))
  h <- fbs17_heard(maihda_fit_lmer_growth(quote(fake_fit()), env, "pid", first = first))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  env$candidates$opt <- c(n, n, n)
  h <- fbs17_heard(maihda_fit_lmer_growth(quote(fake_fit()), env, "pid", first = first))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, "warning: lme4 flags the fit handed in")
})

# ---- fit_maihda() on panels whose optimum is on a bound ---------------------------------------

test_that("fit_maihda() returns the singular fit when lmer's own fit stopped short of the bound", {
  skip_on_cran()
  for (first in c("short flagged", "short accepted")) {
    h <- fbs17_fit_with(first = first, restarts = "singular")
    m <- h$value
    # The fit kept is the restart on the bound, not lmer's own.
    expect_false(lme4::isSingular(h$own))
    expect_true(lme4::isSingular(m$model))
    expect_true(is.numeric(stats::getCall(m$model)$start))
    expect_identical(unname(lme4::getME(m$model, "theta")), h$optimum)
    expect_lt(lme4::REMLcrit(m$model), lme4::REMLcrit(h$own))
    # lme4's note of the singular fit reaches the caller; its warning about the fit short of
    # the bound does not, nor the fit's record.
    expect_identical(h$seen, paste0("message: ", fbs17_note, "\n"))
    expect_true(m$diagnostics$singular)
    expect_identical(m$diagnostics$singular_terms, "stratum")
    expect_true(m$diagnostics$converged)
    expect_null(m$diagnostics$lme4_overruled)
    expect_null(m$diagnostics$stationary_distance)
    expect_match(fbs17_printed(m), "Singular fit: the random-effects block for 'stratum'", fixed = TRUE)
    expect_false(grepl("lme4 warned|failed to converge", fbs17_printed(m)))
    # The stratum slope is at the boundary for the fit kept and was not for lmer's own, whose
    # stratum intercept-slope correlation is short of one.
    expect_true(maihda_stratum_slope_at_boundary(m))
    expect_false(maihda_stratum_slope_at_boundary_lme4(h$own))
  }
})

test_that("fit_maihda() returns the fit on both bounds when lmer's own fit is short of one of them", {
  skip_on_cran()
  # lmer's own fit on the stratum slope's bound and short of the intercept's, a fit lme4 calls
  # singular already: the restart on both.
  h <- fbs17_fit_with(first = "slope only", restarts = "both", panel = "flat")
  m <- h$value
  expect_identical(maihda_lme4_bounds_reached(h$own), c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(maihda_lme4_bounds_reached(m$model), c(FALSE, FALSE, TRUE, TRUE))
  expect_true(is.numeric(stats::getCall(m$model)$start))
  expect_identical(unname(lme4::getME(m$model, "theta")), h$optimum)
  expect_identical(h$seen, paste0("message: ", fbs17_note, "\n"))
  expect_true(m$diagnostics$singular)
  # The whole stratum block is at the boundary for the fit kept, which a PCV reads as no
  # between-stratum variance to explain, and was not for lmer's own.
  expect_true(maihda_stratum_growth_at_boundary_lme4(m$model))
  expect_false(maihda_stratum_growth_at_boundary_lme4(h$own))
  # From a fit on neither bound, with a restart on one and a restart on both: the one on both,
  # whichever comes first.
  for (order in list(c("slope only", "both"), c("both", "slope only"))) {
    h <- fbs17_fit_with(first = "neither flagged", restarts = order, panel = "flat")
    expect_identical(maihda_lme4_bounds_reached(h$value$model), c(FALSE, FALSE, TRUE, TRUE))
    expect_identical(h$seen, paste0("message: ", fbs17_note, "\n"))
  }
})

test_that("a singular restart that is not at an optimum leaves the growth fit as it was", {
  skip_on_cran()
  # Under an optimiser's code other than 0, or a singular fit that is not the optimum: lmer's
  # own flagged fit and the record that overrules lme4. Nothing of the singular restart
  # reaches the caller; whether lme4's warning about the fit kept does is not this file's to
  # pin.
  for (left in c("failed", "elsewhere")) {
    h <- fbs17_fit_with(first = "short flagged", restarts = left)
    m <- h$value
    expect_identical(setdiff(h$seen, paste("warning:", fbs17_flag)), character(0))
    expect_false(lme4::isSingular(m$model))
    expect_null(stats::getCall(m$model)$start)
    expect_false(m$diagnostics$singular)
    expect_identical(m$diagnostics$lme4_overruled, fbs17_flag)
    expect_match(fbs17_printed(m), "lme4 warned that the model failed to converge", fixed = TRUE)
  }
  # A singular restart on the bound that the package's measure cannot vouch for.
  local_mocked_bindings(maihda_lme4_singular_at_optimum = function(model, tol = 1e-4) FALSE)
  h <- fbs17_fit_with(first = "short accepted", restarts = "singular")
  expect_identical(h$seen, character(0))
  expect_false(lme4::isSingular(h$value$model))
  expect_null(stats::getCall(h$value$model)$start)
  h <- fbs17_fit_with(first = "slope only", restarts = "both", panel = "flat")
  expect_null(stats::getCall(h$value$model)$start)
  expect_identical(maihda_lme4_bounds_reached(h$value$model), c(FALSE, FALSE, FALSE, TRUE))
})

# The blocks below pin what the remedy leaves alone; they pass without it too.

test_that("a growth fit on its optimum's bounds stays lmer's own beside fits short of them", {
  skip_on_cran()
  for (beside in c("short accepted", "short flagged")) {
    h <- fbs17_fit_with(first = "singular", restarts = beside)
    m <- h$value
    expect_identical(h$seen, paste0("message: ", fbs17_note, "\n"))
    expect_true(lme4::isSingular(m$model))
    expect_null(stats::getCall(m$model)$start)
    expect_true(m$diagnostics$singular)
  }
  h <- fbs17_fit_with(first = "both", restarts = c("slope only", "neither accepted"), panel = "flat")
  expect_identical(h$seen, paste0("message: ", fbs17_note, "\n"))
  expect_null(stats::getCall(h$value$model)$start)
  expect_identical(unname(lme4::getME(h$value$model, "theta")), h$optimum)
})

test_that("a growth fit short of the bound from every start is kept as before", {
  skip_on_cran()
  h <- fbs17_fit_with(first = "short accepted", restarts = "short accepted")
  expect_identical(h$seen, character(0))
  expect_false(lme4::isSingular(h$value$model))
  expect_null(stats::getCall(h$value$model)$start)
  expect_false(h$value$diagnostics$singular)
  h <- fbs17_fit_with(first = "short flagged", restarts = "short flagged")
  expect_identical(setdiff(h$seen, paste("warning:", fbs17_flag)), character(0))
  expect_false(lme4::isSingular(h$value$model))
  expect_null(stats::getCall(h$value$model)$start)
  expect_identical(h$value$diagnostics$lme4_overruled, fbs17_flag)
})
