# Audit 2026-10-13 (packet 2026-10-13-ordinal-ucminf-convergence): an ordinal fit made
# with ordinal's ucminf optimiser was always recorded as not converged -- CONFIRMED.
#
# fit_maihda(engine = "ordinal", control = ordinal::clmm.control(method = "ucminf")) on
# 2 x 3 strata of 30 rows gave a variance of 0.804437 (nlminb: 0.804439) with
# diagnostics$converged FALSE and the message "Stopped by small gradient (grtol).";
# print() and summary() showed it as a convergence problem and compare_maihda_groups()
# named every group. The record read the optimiser's return code as nlminb()'s, where 0
# is convergence. ucminf() has no 0: it returns 1 (gradient tolerance met), 2 (step
# tolerance met), 3 (evaluation limit), 4 (zero step from the line search) and negative
# codes for a search that did not start.
#
# Measured on 2,216 ucminf fits of simulated data, each beside its nlminb twin (8
# designs x 4 true SDs, logit and probit, 3 to 7 classes, Laplace and quadrature,
# flexible and equidistant thresholds, one covariate on five scales, a ladder of
# covariate scales, and 1,200 to 120,000 rows): 1,050 ended with code 1, 1,164 with
# code 4, 2 with code 2 and none with 3. Code 4 grows with the data in that survey: 1%
# of the fits that made a step up to 100 rows, 76% at 700 to 1,500, every one from
# 12,000 (the share depends on the design: 28% at those middle sizes in a second
# survey, 5% under quadrature). The 1,166 fits with 2 or 4 are of three kinds:
#   - at a stationary point: 772 -- 607 within 1e-4 of the best log-likelihood of their
#     model on their data, 164 on their way to an SD of zero, and 1 that is 1.9e-3
#     below an nlminb fit which itself reported false convergence;
#   - returned at the starting values, no step made: 392, an SD of exactly 1 and a
#     log-likelihood 0.3 to 36 below the maximum, with a covariate far from zero or on
#     a large scale (birth year, currency units; age in years at 12,000 rows);
#   - stopped short after some steps: 2, at 36,000 and 120,000 rows, 1.65 and 0.027
#     log-likelihood units below nlminb's fit (SD 0.985 against 0.809).
# ucminf() itself reports 1, 2 and 4 as convergence, all three kinds included. What
# tells them apart is the fit: with g the gradient and H the Hessian clmm() computes
# where it returns, g' H^-1 g / 2 is the log-likelihood the fit is short of a
# stationary point by. It equalled the measured shortfall wherever that could be
# measured (median ratio 1.01). Where the Hessian is positive definite, in 1,064 of
# the 1,166, it was at most 3.7e-4 for the first kind and at least 0.026 for the other
# two.
#
# Per the user: 1 counts; 3 and the negative codes do not; 2 and 4 are judged by the
# fit, converged within a tenth of a standard error of a stationary point
# (g' H^-1 g < 0.01), not converged beyond it with the distance in the record, and
# unknown where the fit cannot say -- except a search that made no step, which is then
# not converged, its estimates being the starting values. Which optimiser made a
# return is read off the return: clmm()'s general routine runs nlminb() when asked for
# ucminf().
#
# When the fit cannot say. Its Hessian is not positive definite (a fit on its way to
# an SD of zero, mostly). Or the step of 1e-4 by which clmm() differences every
# parameter (ordinal's deriv12()) is too coarse for one of them: what the stored
# gradient and Hessian then give is the error of the differences, not the fit's
# distance. There are two ways into that, each measured at a point that IS the maximum
# (the fit on the standardised covariate, carried exactly to another spelling of the
# covariate and handed to clmm() as its start):
#   - a coefficient whose standard error is below the step, a covariate in the
#     thousands on a few thousand rows (found in a second survey, made to check this
#     fix). The distance read was 4e-4 to 2e-3 times the squared ratio of the step to
#     the standard error, past the line of 0.1 once the ratio is 8 to 16, and sound
#     fits were recorded as not converged with a distance that was not theirs;
#   - a covariate far from zero, which ties its coefficient to the thresholds: a step
#     along the coefficient moves every observation's predictor by the step times the
#     covariate. With a birth year and classes of unequal size (70 / 15 / 10 / 5%), 7
#     of 16 maxima read 0.13 to 0.35 standard errors out with every standard error far
#     above the step, and of 134 fits from clmm()'s own start that had reached their
#     maximum three were recorded as not converged (two with a birth year).
# Both are one quantity, the size the third-order error of a central difference takes
# in the distance: step^2 * H_jj^1.5 * sqrt((H^-1)_jj) for parameter j, the squared
# step in conditional standard errors times the ratio of the marginal standard error
# to the conditional one. On 346 maxima of 37 designs the distance read was at most
# 0.003 times it wherever it was above 30; where it was at most 1 the largest distance
# read was 0.013 (172 maxima), and no reading passed 0.1 below 62. Past 1 the fit is
# not asked: of those 134 fits, 115 are then converged, 19 unknown and none not
# converged. On ordinary scales the 1,318 sound fits of the second survey were at most
# 0.052 standard errors out. At covariate SDs of 1e-6 and less the stored derivatives
# are noise the other way and nothing in the fit shows it; nlminb's record is no better
# there.
#
# No block draws a random number. The live blocks assert what holds whichever of its
# stops ucminf() takes on the platform.

uc13_observe <- function(expr) {
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(res = res, warns = warns)
}

# gender x ses strata of 30 rows: a stratum's outcomes are its shift plus the same
# fixed spread, cut into four classes, with the same thirty ages in every stratum.
uc13_years <- c(53, 49, 55, 51, 55, 43, 46, 55, 71, 56, 19, 53, 23, 41, 63, 45, 53, 54,
                54, 55, 38, 61, 23, 36, 55, 47, 51, 53, 41, 66)
uc13_strata <- function(shift = rep(0, 6)) {
  g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 30), , drop = FALSE]
  d$y <- rep(shift, each = 30) + rep(seq(-1.45, 1.45, length.out = 30), times = nrow(g))
  d$yo <- factor(findInterval(d$y, c(-0.8, 0, 0.8)) + 1L, levels = 1:4, ordered = TRUE)
  d$age <- rep(uc13_years, times = nrow(g))
  # One covariate on three scales: standardised, in years and in currency units.
  d$z <- (d$age - 45) / 12
  d$income <- 40000 + 1250 * (d$age - 45)
  rownames(d) <- NULL
  d
}
uc13_shift <- c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2)

uc13_fit <- function(d, formula = yo ~ (1 | gender:ses), ...) {
  suppressMessages(suppressWarnings(
    fit_maihda(formula, data = d, engine = "ordinal", ...)))
}
uc13_variance <- function(m) maihda_clmm_variances(m)$stratum
uc13_loglik <- function(m) as.numeric(stats::logLik(m$model))
uc13_printed <- function(x) utils::capture.output(suppressWarnings(print(x)))

# What ucminf() returns, as far as the package reads it: a search that did not start
# (a negative code) reports nothing in $info, and older versions of ucminf set no class.
uc13_ucminf <- function(code, evaluations = 20, classed = TRUE) {
  said <- c("1" = "Stopped by small gradient (grtol).",
            "2" = "Stopped by small step (xtol).",
            "3" = "Stopped by function evaluation limit (maxeval)",
            "4" = "Stopped by zero step from line search",
            "-5" = "Computation did not start: grtol or xtol <= 0.",
            "-6" = "Computation did not start: maxeval <= 0.")
  out <- list(par = c(-1.2, 0, 1.2, -0.9), value = 250, convergence = as.integer(code),
              message = said[[as.character(code)]])
  if (code > 0) {
    out$invhessian.lt <- rep(0.1, 10)
    out$info <- c(maxgradient = 2e-5, laststep = 1e-6, stepmax = 0.015,
                  neval = evaluations)
  }
  if (classed) {
    class(out) <- "ucminf"
  }
  out
}
# What nlminb() returns.
uc13_nlminb <- function(code = 0L, message = "relative convergence (4)") {
  list(par = c(-1.2, 0, 1.2, -0.9), objective = 250, convergence = as.integer(code),
       iterations = 12L, evaluations = c("function" = 15L, gradient = 60L),
       message = message)
}
# A stand-in for a clmm fit holding what the helpers read. `distance` gives it a
# gradient and a Hessian that put it that many standard errors from a stationary
# point; `curvature` holds the Hessian's eigenvalues, the first of which must be 4.
# Without `distance` the fit holds neither.
uc13_some <- c(0.3, -0.2, 0.1, 0, -0.4, 0.2)
uc13_zero <- rep(0, 6)
uc13_record <- function(returned, distance = NULL, curvature = c(4, 9, 1, 2), sd = 0.4,
                        effects = uc13_some, general = FALSE, method = "ucminf",
                        n_quad = 1L, loglik = -250) {
  fit <- structure(list(
    ST = list(matrix(sd, 1, 1, dimnames = list("(Intercept)", "(Intercept)"))),
    dims = list(nlev.re = c(stratum = length(effects))),
    coefficients = c("1|2" = -1.2, "2|3" = 0, "3|4" = 1.2),
    ranef = effects, control = list(useMatrix = general, method = method),
    nAGQ = n_quad, link = "logit", optRes = returned, logLik = loglik),
    class = "clmm")
  if (!is.null(distance)) {
    fit$gradient <- c(2 * distance, 0, 0, 0)
    fit$Hessian <- diag(curvature)
  }
  fit
}
uc13_no_step <- paste0("ucminf() made no step from the starting values (Stopped by zero ",
                       "step from line search): the estimates returned are those values.")
uc13_stopped <- function(distance, said = "Stopped by zero step from line search") {
  paste0("ucminf() stopped ", distance, " standard errors from a stationary point by ",
         "the fit's own gradient and Hessian (", said, ").")
}
# What is recorded of a live ucminf() fit whose gradient and Hessian cannot say:
# converged on ucminf()'s own gradient tolerance, not converged when the search made no
# step, unknown otherwise -- and never a distance.
uc13_expect_fit_silent <- function(m) {
  returned <- m$model$optRes
  expect_true(is.na(maihda_clmm_stationary_distance(m$model)))
  expect_false(any(grepl("standard errors from a stationary point",
                         m$diagnostics$messages, fixed = TRUE)))
  if (returned$convergence == 1) {
    expect_true(m$diagnostics$converged)
  } else if (maihda_ucminf_made_no_step(returned)) {
    expect_false(m$diagnostics$converged)
    expect_length(m$diagnostics$messages, 1L)
    expect_match(m$diagnostics$messages, "ucminf() made no step from the starting values",
                 fixed = TRUE)
  } else {
    expect_true(is.na(m$diagnostics$converged))
    expect_length(m$diagnostics$messages, 0L)
  }
}

test_that("a return is read by the rules of the optimiser that made it", {
  skip_if_not_installed("ordinal")
  expect_true(maihda_is_ucminf_return(uc13_ucminf(1)))
  expect_true(maihda_is_ucminf_return(uc13_ucminf(4)))
  # Without the class, the return is known by what ucminf() reports of its search.
  expect_true(maihda_is_ucminf_return(uc13_ucminf(1, classed = FALSE)))
  expect_true(maihda_is_ucminf_return(uc13_ucminf(4, classed = FALSE)))
  expect_true(maihda_is_ucminf_return(uc13_ucminf(-5)))
  expect_false(maihda_is_ucminf_return(uc13_ucminf(-5, classed = FALSE)))
  expect_false(maihda_is_ucminf_return(uc13_nlminb()))
  expect_false(maihda_is_ucminf_return(uc13_nlminb(1L, "false convergence (8)")))
  # A message alone is no evidence of the optimiser.
  expect_false(maihda_is_ucminf_return(
    list(convergence = 1L, message = "Stopped by small gradient (grtol).")))
  half <- uc13_nlminb()
  half$info <- c(maxgradient = 1e-5)
  expect_false(maihda_is_ucminf_return(half))
  expect_false(maihda_is_ucminf_return(NULL))
  expect_false(maihda_is_ucminf_return("ucminf"))

  # The same code, two verdicts.
  expect_true(maihda_clmm_converged(uc13_record(uc13_ucminf(1))))
  expect_true(maihda_clmm_converged(uc13_record(uc13_ucminf(1, classed = FALSE))))
  expect_false(maihda_clmm_converged(
    uc13_record(uc13_nlminb(1L, "false convergence (8)"), method = "nlminb")))
  expect_true(maihda_clmm_converged(uc13_record(uc13_nlminb(0L), method = "nlminb")))
  # ucminf() has no code 0.
  zero <- uc13_ucminf(1)
  zero$convergence <- 0L
  expect_false(maihda_clmm_converged(uc13_record(zero, distance = 0.001)))
  # The fit's gradient and Hessian are asked about ucminf()'s 2 and 4 only: nlminb()'s
  # return is its own verdict, whatever they say.
  expect_true(maihda_clmm_converged(
    uc13_record(uc13_nlminb(0L), distance = 3, method = "nlminb")))
  expect_false(maihda_clmm_converged(
    uc13_record(uc13_nlminb(1L, "false convergence (8)"), distance = 0.001,
                method = "nlminb")))
  expect_true(maihda_clmm_converged(uc13_record(uc13_ucminf(4), distance = 0.001)))
  expect_false(maihda_clmm_converged(uc13_record(uc13_ucminf(4), distance = 3)))

  # clmm()'s general routine asked for ucminf() runs nlminb(): the settings name
  # ucminf, the return is nlminb()'s and is read as that.
  general <- function(returned, sd = 0.4) {
    uc13_record(returned, sd = sd, effects = if (sd == 0) uc13_zero else uc13_some,
                general = TRUE, method = "ucminf")
  }
  expect_true(maihda_clmm_converged(general(uc13_nlminb(0L))))
  expect_false(maihda_clmm_converged(general(uc13_nlminb(1L, "false convergence (8)"))))
  expect_false(maihda_clmm_converged(
    general(uc13_nlminb(1L, "iteration limit reached without convergence (10)"))))
  expect_true(maihda_clmm_converged(
    general(uc13_nlminb(1L, "singular convergence (7)"), sd = 0)))
  expect_false(maihda_clmm_converged(general(uc13_nlminb(1L, "singular convergence (7)"))))
  # nlminb()'s boundary return has no ucminf() counterpart.
  at_zero <- uc13_ucminf(3)
  at_zero$message <- "singular convergence (7)"
  expect_false(maihda_clmm_converged(uc13_record(at_zero, sd = 0, effects = uc13_zero)))

  # No code to read: no verdict.
  none <- uc13_ucminf(1)
  none$convergence <- NULL
  expect_true(is.na(maihda_clmm_converged(uc13_record(none))))
  unknown <- uc13_ucminf(1)
  unknown$convergence <- NA_integer_
  expect_true(is.na(maihda_clmm_converged(uc13_record(unknown))))
  # A fit kept on the first fit's convergence carries that verdict, whatever made it.
  marked <- uc13_record(uc13_ucminf(3))
  expect_false(maihda_clmm_converged(marked))
  attr(marked, "maihda_converged") <- TRUE
  expect_true(maihda_clmm_converged(marked))
})

test_that("how far a fit is from a stationary point, by its gradient and Hessian", {
  expect_equal(maihda_clmm_stationary_distance(uc13_record(NULL, distance = 0.25)), 0.25)
  expect_equal(maihda_clmm_stationary_distance(uc13_record(NULL, distance = 3)), 3)
  expect_identical(maihda_clmm_stationary_distance(uc13_record(NULL, distance = 0)), 0)
  # sqrt(g' H^-1 g), worked by hand: H^-1 g = (4, -3) / 5 and g' H^-1 g = 7 / 5.
  fit <- list(gradient = c(a = 1, b = -1), Hessian = matrix(c(2, 1, 1, 3), 2))
  expect_equal(maihda_clmm_stationary_distance(fit), sqrt(1.4))
  # Half its square is the log-likelihood a quadratic is short of its maximum by.
  top <- c(0.4, -1.1)
  at <- c(0.7, -0.6)
  curve <- matrix(c(2, 1, 1, 3), 2)
  short <- 0.5 * drop(t(at - top) %*% curve %*% (at - top))
  quadratic <- list(gradient = drop(curve %*% (at - top)), Hessian = curve)
  expect_equal(maihda_clmm_stationary_distance(quadratic)^2 / 2, short)
  # The algebra does not depend on the units of a parameter: a covariate in years in
  # place of standard deviations divides its coefficient, and that coefficient's
  # standard error, by 12.
  rescaled <- function(by) {
    units <- diag(c(1, 1 / by))
    list(gradient = drop(solve(units) %*% fit$gradient),
         Hessian = solve(units) %*% fit$Hessian %*% solve(units))
  }
  expect_equal(maihda_clmm_stationary_distance(rescaled(12)), sqrt(1.4))
  expect_equal(maihda_clmm_stationary_distance(rescaled(1 / 300)), sqrt(1.4))
  # The derivatives clmm() stores do: they are central differences with a step of 1e-4
  # in every parameter. How coarse that step is for a parameter is the squared step in
  # its conditional standard errors times the ratio of its marginal standard error to
  # its conditional one, and past 1 the fit cannot say.
  expect_identical(.maihda_clmm_derivative_step, 1e-4)
  # For a parameter tied to no other it is the squared ratio of the step to its
  # standard error.
  expect_equal(maihda_clmm_step_coarseness(diag(c(4, 1e10)), diag(c(1 / 4, 1e-10))),
               c(4e-8, 100))
  # By hand for H = (2, 1; 1, 3) and H^-1 = (3, -1; -1, 2) / 5.
  expect_equal(maihda_clmm_step_coarseness(fit$Hessian, solve(fit$Hessian)),
               1e-8 * c(2^1.5 * sqrt(3 / 5), 3^1.5 * sqrt(2 / 5)))
  # In currency units the coefficient's standard error, sqrt(2 / 5) / 15000, is below
  # the step: the difference spans more than a standard error to either side.
  expect_true(is.na(maihda_clmm_stationary_distance(rescaled(15000))))
  fine <- function(se, gradient = c(0.02, 0)) {
    list(gradient = gradient, Hessian = diag(c(4, 1 / se^2)))
  }
  expect_equal(maihda_clmm_stationary_distance(fine(1)), 0.01)
  expect_equal(maihda_clmm_stationary_distance(fine(1.1e-4)), 0.01)
  expect_true(is.na(maihda_clmm_stationary_distance(fine(0.9e-4))))
  expect_true(is.na(maihda_clmm_stationary_distance(fine(7e-6))))
  # ... however far out the stored gradient puts the fit, and whichever parameter it is.
  expect_true(is.na(maihda_clmm_stationary_distance(fine(7e-6, c(5, 3e5)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(0, 0.02), Hessian = diag(c(1 / 7e-6^2, 4))))))
  # A standard error far above the step is no obstacle.
  expect_equal(maihda_clmm_stationary_distance(fine(1e4)), 0.01)
  # A covariate far from zero ties its coefficient to the threshold. Here a threshold
  # with a standard error of 1 / 10 and a coefficient with one of 1 / 20; moving the
  # covariate by `by` leaves the fit's distance and the coefficient's standard error
  # where they were and makes the threshold the one at a covariate of zero. A step
  # along the coefficient then moves the predictor by `by` steps.
  moved <- function(by, distance = 0.01) {
    old_from_new <- matrix(c(1, 0, -by, 1), 2)
    list(gradient = drop(t(old_from_new) %*% c(10 * distance, 0)),
         Hessian = t(old_from_new) %*% diag(c(100, 400)) %*% old_from_new)
  }
  for (by in c(0, 45, 100)) {
    expect_equal(maihda_clmm_stationary_distance(moved(by)), 0.01, tolerance = 1e-6)
    expect_equal(maihda_clmm_stationary_distance(moved(by, 2.3)), 2.3, tolerance = 1e-6)
  }
  # In years the step is fine for the coefficient (0.046). At a birth year it is not
  # (3 850), though no standard error is anywhere near as small as the step.
  in_years <- moved(45)
  expect_equal(max(maihda_clmm_step_coarseness(in_years$Hessian, solve(in_years$Hessian))),
               1e-8 * (100 * 45^2 + 400)^1.5 / 20, tolerance = 1e-6)
  born <- moved(1975)
  expect_gt(min(sqrt(diag(solve(born$Hessian)))), 0.04)
  expect_gt(max(maihda_clmm_step_coarseness(born$Hessian, solve(born$Hessian))), 3000)
  expect_true(is.na(maihda_clmm_stationary_distance(born)))
  expect_true(is.na(maihda_clmm_stationary_distance(moved(1975, 2.3))))
  # The line is at 1: 1e-8 * (100 by^2 + 400)^1.5 / 20 is 0.95 at 124 and 1.02 at 127.
  expect_equal(maihda_clmm_stationary_distance(moved(124)), 0.01, tolerance = 1e-6)
  expect_true(is.na(maihda_clmm_stationary_distance(moved(127))))

  # No pair of a gradient and a Hessian that go together: nothing to say.
  expect_true(is.na(maihda_clmm_stationary_distance(uc13_record(NULL))))
  expect_true(is.na(maihda_clmm_stationary_distance(list(gradient = c(1, -1)))))
  expect_true(is.na(maihda_clmm_stationary_distance(list(Hessian = diag(2)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(1, -1, 0), Hessian = diag(2)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(1, -1), Hessian = c(2, 3)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(1, NA), Hessian = diag(2)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(1, -1), Hessian = matrix(c(2, 1, 1, Inf), 2)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(1, -1), Hessian = matrix(list(2, 0, 0, 3), 2)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = function() 1, Hessian = diag(2)))))
  expect_true(is.na(maihda_clmm_stationary_distance(NULL)))
  # A Hessian that is not positive definite: a saddle, a flat direction.
  expect_true(is.na(maihda_clmm_stationary_distance(
    uc13_record(NULL, distance = 0.01, curvature = c(4, 9, 1, -2)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    uc13_record(NULL, distance = 0.01, curvature = c(4, 9, 1, 0)))))
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(0, 0), Hessian = matrix(c(1, 2, 2, 1), 2)))))
  # ... also one whose inverse has a positive diagonal.
  expect_true(is.na(maihda_clmm_stationary_distance(
    list(gradient = c(1, 0), Hessian = matrix(c(-1, 3, 3, -1), 2)))))
})

test_that("which of ucminf()'s stops count as convergence", {
  near <- uc13_record(NULL, distance = 0.01)
  far <- uc13_record(NULL, distance = 2.3)
  blind <- uc13_record(NULL)
  saddle <- uc13_record(NULL, distance = 0.01, curvature = c(4, 9, 1, -2))
  # A parameter with a standard error of 3e-5, below the step of clmm()'s derivatives,
  # beside a stored gradient that would put the fit 2.3 standard errors out, or 0.01.
  fine <- uc13_record(NULL, distance = 2.3, curvature = c(4, 9, 1, 1e9))
  fine_near <- uc13_record(NULL, distance = 0.01, curvature = c(4, 9, 1, 1e9))

  # Its gradient tolerance met: convergence, whatever the fit says and wherever.
  for (fit in list(near, far, blind, saddle, fine)) {
    expect_true(maihda_ucminf_converged(uc13_ucminf(1), fit))
    expect_true(maihda_ucminf_converged(uc13_ucminf(1, evaluations = 1), fit))
    # Out of evaluations, or not started: not, whatever the fit says.
    expect_false(maihda_ucminf_converged(uc13_ucminf(3), fit))
    expect_false(maihda_ucminf_converged(uc13_ucminf(3, evaluations = 500), fit))
    expect_false(maihda_ucminf_converged(uc13_ucminf(-5), fit))
    expect_false(maihda_ucminf_converged(uc13_ucminf(-6), fit))
  }
  for (code in c(2, 4)) {
    for (classed in c(TRUE, FALSE)) {
      stepped <- uc13_ucminf(code, classed = classed)
      still <- uc13_ucminf(code, evaluations = 1, classed = classed)
      # No further step could be made: the fit says whether that is a maximum ...
      expect_true(maihda_ucminf_converged(stepped, near))
      expect_false(maihda_ucminf_converged(stepped, far))
      # ... also of a search that made no step: one started on a maximum has none to
      # make.
      expect_true(maihda_ucminf_converged(still, near))
      expect_false(maihda_ucminf_converged(still, far))
      # Where the fit cannot say, a search that made no step is not converged and any
      # other is unknown.
      for (fit in list(blind, saddle, fine, fine_near)) {
        expect_true(is.na(maihda_ucminf_converged(stepped, fit)))
        expect_false(maihda_ucminf_converged(still, fit))
        expect_true(is.na(maihda_ucminf_converged(
          uc13_ucminf(code, evaluations = 2, classed = classed), fit)))
      }
    }
    # A tenth of a standard error is the line.
    expect_true(maihda_ucminf_converged(uc13_ucminf(code),
                                        uc13_record(NULL, distance = 0.099)))
    expect_false(maihda_ucminf_converged(uc13_ucminf(code),
                                         uc13_record(NULL, distance = 0.101)))
    expect_true(maihda_ucminf_converged(uc13_ucminf(code), uc13_record(NULL, distance = 0)))
    # Without a count of evaluations the search is not taken to have stood still.
    uncounted <- uc13_ucminf(code)
    uncounted$info <- c(maxgradient = 2e-5, laststep = 1e-6)
    expect_true(is.na(maihda_ucminf_converged(uncounted, blind)))
    expect_true(maihda_ucminf_converged(uncounted, near))
    expect_false(maihda_ucminf_converged(uncounted, far))
  }
  expect_identical(.maihda_clmm_stationary_within, 0.1)

  expect_true(maihda_ucminf_made_no_step(uc13_ucminf(4, evaluations = 1)))
  expect_true(maihda_ucminf_made_no_step(uc13_ucminf(1, evaluations = 1L)))
  expect_false(maihda_ucminf_made_no_step(uc13_ucminf(4, evaluations = 2)))
  expect_false(maihda_ucminf_made_no_step(uc13_ucminf(4, evaluations = 24)))
  # No count to read: not known to have stood still.
  expect_false(maihda_ucminf_made_no_step(uc13_ucminf(-5)))
  uncounted <- uc13_ucminf(4)
  uncounted$info <- c(maxgradient = 2e-5, neval = NA_real_)
  expect_false(maihda_ucminf_made_no_step(uncounted))
  uncounted$info <- c(maxgradient = 2e-5, laststep = 0)
  expect_false(maihda_ucminf_made_no_step(uncounted))
  # nlminb()'s `evaluations` is not ucminf()'s count.
  expect_false(maihda_ucminf_made_no_step(uc13_nlminb()))
  expect_false(maihda_ucminf_made_no_step(NULL))
})

test_that("the fit record of a ucminf() fit", {
  skip_if_not_installed("ordinal")
  converged <- list(uc13_record(uc13_ucminf(1)),
                    uc13_record(uc13_ucminf(1), distance = 2.3),
                    uc13_record(uc13_ucminf(4), distance = 0.01),
                    uc13_record(uc13_ucminf(2), distance = 0.01),
                    uc13_record(uc13_ucminf(4, classed = FALSE), distance = 0.01),
                    uc13_record(uc13_ucminf(4, evaluations = 1), distance = 0.01))
  for (fit in converged) {
    rec <- maihda_fit_diagnostics(fit)
    expect_identical(rec$engine, "ordinal")
    expect_true(rec$converged)
    expect_false(rec$singular)
    expect_length(rec$messages, 0L)
    # Nothing to print.
    expect_length(maihda_format_fit_diagnostics(rec), 0L)
  }
  # A search that was cut short or did not start: not converged, in ucminf()'s words.
  for (returned in list(uc13_ucminf(3), uc13_ucminf(-5), uc13_ucminf(-6))) {
    rec <- maihda_fit_diagnostics(uc13_record(returned, distance = 0.01))
    expect_false(rec$converged)
    expect_identical(rec$messages, returned$message)
    expect_identical(maihda_format_fit_diagnostics(rec),
                     c("Convergence problems of the ordinal::clmm fit:",
                       paste0("  - ", returned$message)))
  }
  # A search that stopped short after some steps: the record says how far.
  rec <- maihda_fit_diagnostics(uc13_record(uc13_ucminf(4), distance = 2.3087))
  expect_false(rec$converged)
  expect_false(rec$singular)
  expect_identical(rec$messages, uc13_stopped("2.3"))
  expect_identical(maihda_format_fit_diagnostics(rec),
                   c("Convergence problems of the ordinal::clmm fit:",
                     paste0("  - ", uc13_stopped("2.3"))))
  rec <- maihda_fit_diagnostics(uc13_record(uc13_ucminf(2, classed = FALSE),
                                            distance = 0.2273))
  expect_false(rec$converged)
  expect_identical(rec$messages, uc13_stopped("0.23", "Stopped by small step (xtol)"))
  # A search that made no step, away from a stationary point or with a fit that cannot
  # say: the record says what the estimates are.
  for (fit in list(uc13_record(uc13_ucminf(4, evaluations = 1), distance = 1.6, sd = 1),
                   uc13_record(uc13_ucminf(4, evaluations = 1), sd = 1),
                   uc13_record(uc13_ucminf(4, evaluations = 1), distance = 0.01,
                               curvature = c(4, 9, 1, -2), sd = 1),
                   uc13_record(uc13_ucminf(4, evaluations = 1), distance = 0.01,
                               curvature = c(4, 9, 1, 1e9), sd = 1))) {
    rec <- maihda_fit_diagnostics(fit)
    expect_false(rec$converged)
    expect_false(rec$singular)
    expect_identical(rec$messages, uc13_no_step)
    expect_identical(maihda_format_fit_diagnostics(rec),
                     c("Convergence problems of the ordinal::clmm fit:",
                       paste0("  - ", uc13_no_step)))
  }
  small_step <- maihda_fit_diagnostics(uc13_record(uc13_ucminf(2, evaluations = 1), sd = 1))
  expect_false(small_step$converged)
  expect_identical(small_step$messages,
                   paste0("ucminf() made no step from the starting values (Stopped by ",
                          "small step (xtol)): the estimates returned are those values."))
  # A search that made steps and a fit that cannot say: unknown, and nothing printed --
  # also when the stored gradient of a parameter too fine for it reads as far out.
  for (fit in list(uc13_record(uc13_ucminf(4)),
                   uc13_record(uc13_ucminf(2), distance = 0.01,
                               curvature = c(4, 9, 1, -2)),
                   uc13_record(uc13_ucminf(4), distance = 2.3,
                               curvature = c(4, 9, 1, 1e9)))) {
    rec <- maihda_fit_diagnostics(fit)
    expect_true(is.na(rec$converged))
    expect_length(rec$messages, 0L)
    expect_length(maihda_format_fit_diagnostics(rec), 0L)
  }

  # What the record says of a fit, optimiser by optimiser.
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(4, 1), distance = 1.6)),
    uc13_no_step)
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(4, 1, classed = FALSE))),
    uc13_no_step)
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(4), distance = 16.2)),
    uc13_stopped("16"))
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(4, 2), distance = 0.101)),
    uc13_stopped("0.1"))
  wordless <- uc13_ucminf(4, 1)
  wordless$message <- NULL
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(wordless)),
    paste0("ucminf() made no step from the starting values: the estimates returned are ",
           "those values."))
  wordless <- uc13_ucminf(4)
  wordless$message <- NULL
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(wordless, distance = 0.5)),
    paste0("ucminf() stopped 0.5 standard errors from a stationary point by the fit's ",
           "own gradient and Hessian."))
  # Only a search that could go no further is described by where it left the fit, and
  # one that made steps only by a distance the fit can give and that is past the line.
  expect_identical(maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(4))),
                   "Stopped by zero step from line search")
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(4), distance = 0.099)),
    "Stopped by zero step from line search")
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(4), distance = 2.3,
                                              curvature = c(4, 9, 1, 1e9))),
    "Stopped by zero step from line search")
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(3, 1), distance = 2.3)),
    "Stopped by function evaluation limit (maxeval)")
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(3), distance = 2.3)),
    "Stopped by function evaluation limit (maxeval)")
  expect_identical(
    maihda_clmm_optimiser_message(uc13_record(uc13_ucminf(1, 1), distance = 2.3)),
    "Stopped by small gradient (grtol).")
  expect_identical(
    maihda_clmm_optimiser_message(
      uc13_record(uc13_nlminb(1L, "false convergence (8)"), distance = 2.3)),
    "false convergence (8)")
  # A return that is not ucminf()'s is not described as its search.
  foreign <- list(convergence = 4L, message = "another optimiser's words",
                  info = c(neval = 1))
  expect_false(maihda_is_ucminf_return(foreign))
  expect_identical(maihda_clmm_optimiser_message(uc13_record(foreign, distance = 2.3)),
                   "another optimiser's words")
  nameless <- uc13_record(uc13_nlminb(1L, "false convergence (8)"))
  nameless$optRes$message <- ""
  expect_length(maihda_clmm_optimiser_message(nameless), 0L)
  nameless$optRes <- NULL
  expect_length(maihda_clmm_optimiser_message(nameless), 0L)
  # nlminb()'s failure keeps its own words in the record.
  rec <- maihda_fit_diagnostics(
    uc13_record(uc13_nlminb(1L, "false convergence (8)"), method = "nlminb"))
  expect_false(rec$converged)
  expect_identical(rec$messages, "false convergence (8)")

  # A ucminf() fit on the boundary is singular and converged.
  bound <- maihda_fit_diagnostics(uc13_record(uc13_ucminf(1), sd = 5e-5, n_quad = 7L))
  expect_true(bound$singular)
  expect_true(bound$converged)
  expect_length(bound$messages, 0L)
  # One left with every stratum effect at zero beside an SD that is not zero stays
  # flagged, whatever ucminf() and the fit say of its gradient.
  for (returned in list(uc13_ucminf(1, evaluations = 1), uc13_ucminf(4))) {
    flat <- maihda_fit_diagnostics(
      uc13_record(returned, distance = 0.001, sd = 1, effects = uc13_zero, n_quad = 3L))
    expect_false(flat$converged)
    expect_length(flat$messages, 1L)
    expect_match(flat$messages, "every stratum effect at exactly zero", fixed = TRUE)
  }
})

test_that("a boundary fit takes the convergence of a first fit made by ucminf()", {
  skip_if_not_installed("ordinal")
  args <- list(formula = yo ~ (1 | stratum), data = quote(data), link = "logit",
               Hess = TRUE)
  answer <- NULL
  local_mocked_bindings(clmm = function(...) answer, .package = "ordinal")
  run <- function(first) {
    maihda_clmm_general_refit(first, args, list(method = "ucminf"), environment(),
                              character(0))
  }
  first <- function(returned, ...) {
    uc13_record(returned, ..., sd = 1, effects = uc13_zero)
  }
  boundary <- function(returned) {
    uc13_record(returned, sd = 0, effects = uc13_zero, general = TRUE, method = "nlminb")
  }

  # nlminb() could not vouch for the boundary fit it was started on.
  answer <- boundary(uc13_nlminb(1L, "false convergence (8)"))
  for (vouching in list(first(uc13_ucminf(1, evaluations = 1)),
                        first(uc13_ucminf(1)),
                        first(uc13_ucminf(4), distance = 0.01),
                        first(uc13_ucminf(2, classed = FALSE), distance = 0.01),
                        first(uc13_ucminf(4, evaluations = 1), distance = 0.01))) {
    kept <- run(vouching)
    expect_true(attr(kept, "maihda_converged"))
    expect_identical(kept$optRes, answer$optRes)
    rec <- maihda_fit_diagnostics(kept)
    expect_true(rec$converged)
    expect_true(rec$singular)
    expect_length(rec$messages, 0L)
  }
  # A first fit ucminf() and its own gradient do not vouch for leaves the boundary fit
  # with nlminb()'s return: one cut short, one that stopped short, one at its starting
  # values, and two of which nothing is known.
  for (silent in list(first(uc13_ucminf(3), distance = 0.01),
                      first(uc13_ucminf(4), distance = 2.3),
                      first(uc13_ucminf(4, evaluations = 1)),
                      first(uc13_ucminf(-5)),
                      first(uc13_ucminf(4)),
                      first(uc13_ucminf(4), distance = 0.01,
                            curvature = c(4, 9, 1, 1e9)))) {
    kept <- run(silent)
    expect_null(attr(kept, "maihda_converged"))
    expect_identical(kept, answer)
    rec <- maihda_fit_diagnostics(kept)
    expect_false(rec$converged)
    expect_true(rec$singular)
    expect_identical(rec$messages, "false convergence (8)")
  }
  # A boundary fit that converged needs no mark.
  answer <- boundary(uc13_nlminb(0L))
  kept <- run(first(uc13_ucminf(1, evaluations = 1)))
  expect_null(attr(kept, "maihda_converged"))
  expect_true(maihda_fit_diagnostics(kept)$converged)
})

test_that("a ucminf() fit that is nlminb()'s fit is recorded as converged", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- uc13_strata(uc13_shift)
  by_ucminf <- ordinal::clmm.control(method = "ucminf")
  cases <- list(
    list(formula = yo ~ (1 | gender:ses), family = "ordinal"),
    list(formula = yo ~ (1 | gender:ses), family = maihda_cumulative("probit")),
    list(formula = yo ~ z + (1 | gender:ses), family = "ordinal"),
    list(formula = yo ~ age + (1 | gender:ses), family = "ordinal"),
    list(formula = yo ~ gender + ses + (1 | gender:ses), family = "ordinal"))
  for (case in cases) {
    plain <- uc13_fit(d, case$formula, family = case$family)
    expect_true(plain$diagnostics$converged)
    expect_false(maihda_is_ucminf_return(plain$model$optRes))
    spellings <- list(
      uc13_fit(d, case$formula, family = case$family, control = by_ucminf),
      uc13_fit(d, case$formula, family = case$family, method = "ucminf"))
    expect_identical(spellings[[1]]$model$optRes, spellings[[2]]$model$optRes)
    for (m in spellings) {
      returned <- m$model$optRes
      expect_true(maihda_is_ucminf_return(returned))
      expect_identical(m$model$control$method, "ucminf")
      expect_false(m$model$control$useMatrix)
      # It took steps and ended at its gradient tolerance or where none was left ...
      expect_true(returned$convergence %in% c(1, 2, 4))
      expect_false(maihda_ucminf_made_no_step(returned))
      # ... which the fit's own gradient and Hessian put at a stationary point.
      expect_lt(maihda_clmm_stationary_distance(m$model), 0.01)
      expect_equal(uc13_variance(m), uc13_variance(plain), tolerance = 1e-3)
      expect_equal(uc13_loglik(m), uc13_loglik(plain), tolerance = 1e-6)
      expect_true(m$diagnostics$converged)
      expect_false(m$diagnostics$singular)
      expect_length(m$diagnostics$messages, 0L)
      expect_false(any(grepl("Convergence problems", uc13_printed(m), fixed = TRUE)))
    }
  }

  # The summary and the two-model analysis report no problem either.
  m <- uc13_fit(d, control = by_ucminf)
  expect_gt(uc13_variance(m), 0.5)
  obs <- uc13_observe(summary(m))
  expect_length(obs$warns, 0L)
  expect_false(any(grepl("Convergence problems", uc13_printed(obs$res), fixed = TRUE)))
  both <- uc13_observe(maihda(yo ~ (1 | gender:ses), data = d, engine = "ordinal",
                              control = by_ucminf))
  expect_true(both$res$model$diagnostics$converged)
  expect_true(both$res$model_adjusted$diagnostics$converged)
  expect_false(any(grepl("Convergence problems", uc13_printed(both$res), fixed = TRUE)))
  expect_false(any(grepl("converge", both$warns, fixed = TRUE)))
})

test_that("a search started on a fit that is at its maximum leaves it converged", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- uc13_strata(uc13_shift)
  by_ucminf <- ordinal::clmm.control(method = "ucminf")
  for (formula in list(yo ~ (1 | gender:ses), yo ~ age + (1 | gender:ses))) {
    first <- uc13_fit(d, formula, control = by_ucminf)
    expect_true(first$diagnostics$converged)
    again <- uc13_fit(d, formula, control = by_ucminf,
                      start = as.numeric(first$model$optRes$par))
    returned <- again$model$optRes
    expect_true(maihda_is_ucminf_return(returned))
    # ucminf() has no step to make from there, and the fit is the one it was given.
    if (maihda_ucminf_made_no_step(returned)) {
      expect_identical(as.numeric(returned$par), as.numeric(first$model$optRes$par))
      expect_identical(uc13_variance(again), uc13_variance(first))
    }
    expect_equal(uc13_variance(again), uc13_variance(first), tolerance = 1e-6)
    expect_lt(maihda_clmm_stationary_distance(again$model), 0.01)
    expect_true(again$diagnostics$converged)
    expect_length(again$diagnostics$messages, 0L)
  }
})

test_that("a search cut short, or not started, is not converged", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- uc13_strata(uc13_shift)
  plain <- uc13_fit(d)

  short <- uc13_fit(d, control = ordinal::clmm.control(method = "ucminf", maxeval = 5))
  expect_equal(short$model$optRes$convergence, 3)
  expect_false(short$diagnostics$converged)
  expect_identical(short$diagnostics$messages,
                   "Stopped by function evaluation limit (maxeval)")
  # It is not the fit.
  expect_gt(abs(uc13_variance(short) - uc13_variance(plain)), 1e-3)
  printed <- uc13_printed(short)
  at <- which(printed == "  Convergence problems of the ordinal::clmm fit:")
  expect_length(at, 1L)
  expect_identical(printed[at + 1L],
                   "    - Stopped by function evaluation limit (maxeval)")

  for (control in list(ordinal::clmm.control(method = "ucminf", grtol = 0),
                       ordinal::clmm.control(method = "ucminf", maxeval = 0))) {
    idle <- uc13_fit(d, control = control)
    expect_lt(idle$model$optRes$convergence, 0)
    expect_false(idle$diagnostics$converged)
    expect_length(idle$diagnostics$messages, 1L)
    expect_match(idle$diagnostics$messages, "Computation did not start", fixed = TRUE)
    # What comes back is clmm()'s starting SD.
    expect_equal(uc13_variance(idle), 1)
  }
})

test_that("a search that stops short of a maximum is not called converged", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- uc13_strata(uc13_shift)
  by_ucminf <- function(...) ordinal::clmm.control(method = "ucminf", ...)
  reference <- uc13_fit(d)
  expect_true(reference$diagnostics$converged)
  expect_gt(uc13_variance(reference), 0.5)

  # Stopped after its first step by a step tolerance far too loose, and after steps too
  # short to matter by a trust region of next to nothing: code 2 with steps made.
  for (control in list(by_ucminf(xtol = 0.3), by_ucminf(stepmax = 1e-20))) {
    m <- uc13_fit(d, control = control)
    returned <- m$model$optRes
    expect_equal(returned$convergence, 2)
    expect_false(maihda_ucminf_made_no_step(returned))
    distance <- maihda_clmm_stationary_distance(m$model)
    expect_gt(distance, 0.2)
    expect_gt(uc13_loglik(reference) - uc13_loglik(m), 0.02)
    expect_gt(abs(uc13_variance(m) - uc13_variance(reference)), 0.05)
    expect_false(m$diagnostics$converged)
    expect_false(m$diagnostics$singular)
    expect_identical(m$diagnostics$messages,
                     uc13_stopped(format(signif(distance, 2)),
                                  "Stopped by small step (xtol)"))
    printed <- uc13_printed(m)
    expect_true("  Convergence problems of the ordinal::clmm fit:" %in% printed)
    expect_true(paste0("    - ", m$diagnostics$messages) %in% printed)
  }
  # A step tolerance that stops the search early but within a tenth of a standard
  # error: converged, and the fit is the reference's to that precision.
  loose <- uc13_fit(d, control = by_ucminf(xtol = 0.001))
  expect_equal(loose$model$optRes$convergence, 2)
  expect_lt(maihda_clmm_stationary_distance(loose$model), 0.05)
  expect_true(loose$diagnostics$converged)
  expect_length(loose$diagnostics$messages, 0L)
  expect_lt(uc13_loglik(reference) - uc13_loglik(loose), 0.005)

  # One model, the covariate standardised and in currency units: the same maximum.
  scaled <- uc13_fit(d, yo ~ z + (1 | gender:ses))
  expect_true(scaled$diagnostics$converged)
  raw <- uc13_fit(d, yo ~ income + (1 | gender:ses), control = by_ucminf())
  returned <- raw$model$optRes
  expect_true(maihda_is_ucminf_return(returned))
  # Its coefficient is too fine for clmm()'s derivatives: the fit cannot say.
  expect_true(is.na(maihda_clmm_stationary_distance(raw$model)))
  if (maihda_ucminf_made_no_step(returned)) {
    # The starting values come back, with the starting SD of 1.
    expect_equal(returned$convergence, 4)
    expect_identical(unname(as.numeric(returned$par)),
                     unname(as.numeric(raw$model$start)))
    expect_equal(uc13_variance(raw), 1)
    expect_lt(uc13_loglik(raw), uc13_loglik(scaled) - 0.5)
    expect_false(raw$diagnostics$converged)
    expect_false(raw$diagnostics$singular)
    expect_identical(raw$diagnostics$messages, uc13_no_step)
    printed <- uc13_printed(raw)
    expect_true("  Convergence problems of the ordinal::clmm fit:" %in% printed)
    expect_true(paste0("    - ", uc13_no_step) %in% printed)
  }
  # Whatever ucminf() does with it: never converged beside a log-likelihood more than a
  # tenth of a standard error's worth below the maximum.
  for (m in list(raw, loose, uc13_fit(d, control = by_ucminf(xtol = 0.01)),
                 uc13_fit(d, control = by_ucminf(xtol = 0.02)),
                 uc13_fit(d, yo ~ age + (1 | gender:ses), control = by_ucminf()))) {
    top <- if ("age" %in% all.vars(stats::formula(m$model)) ||
                 "income" %in% all.vars(stats::formula(m$model))) scaled else reference
    expect_true(!isTRUE(m$diagnostics$converged) ||
                  uc13_loglik(top) - uc13_loglik(m) < 0.005)
    expect_true(!isFALSE(m$diagnostics$converged) || length(m$diagnostics$messages) > 0L)
  }
})

test_that("a coefficient finer than the step of clmm()'s derivatives: the fit cannot say", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # The step is ordinal's own.
  expect_identical(.maihda_clmm_derivative_step,
                   formals(get("deriv12", envir = asNamespace("ordinal")))$delta)
  d <- uc13_strata(uc13_shift)
  by_ucminf <- ordinal::clmm.control(method = "ucminf")
  top <- uc13_fit(d, yo ~ z + (1 | gender:ses), control = by_ucminf)
  expect_true(top$diagnostics$converged)
  expect_lt(maihda_clmm_stationary_distance(top$model), 0.01)
  # The same model with the covariate in units k times finer, handed its maximum: the
  # thresholds and the SD as they are, the coefficient divided by k.
  carried <- function(k) {
    d$zk <- d$z * k
    at <- as.numeric(top$model$optRes$par)
    at[length(at) - 1L] <- at[length(at) - 1L] / k
    uc13_fit(d, yo ~ zk + (1 | gender:ses), control = by_ucminf, start = at)
  }
  stored <- function(m) {
    g <- as.numeric(m$model$gradient)
    H <- m$model$Hessian
    list(distance = sqrt(sum(g * solve(H, g))), finest = min(sqrt(diag(solve(H)))))
  }
  # With a standard error above the step, the derivatives say what is so: this is the
  # maximum, and the fit is recorded as converged though ucminf() had no step to make.
  for (k in c(100, 500)) {
    m <- carried(k)
    expect_equal(uc13_loglik(m), uc13_loglik(top), tolerance = 1e-8)
    expect_gt(stored(m)$finest, .maihda_clmm_derivative_step)
    expect_lt(maihda_clmm_stationary_distance(m$model), 0.01)
    expect_equal(maihda_clmm_stationary_distance(m$model), stored(m)$distance)
    expect_true(m$diagnostics$converged)
  }
  # With one far below it, they put the same maximum well past a tenth of a standard
  # error from a stationary point. The package reads no distance off them and records
  # none.
  m <- carried(40000)
  expect_equal(uc13_loglik(m), uc13_loglik(top), tolerance = 1e-8)
  expect_equal(uc13_variance(m), uc13_variance(top), tolerance = 1e-6)
  expect_lt(stored(m)$finest, .maihda_clmm_derivative_step / 5)
  expect_gt(stored(m)$distance, 0.2)
  # What is left is whether a step was made: none here, so the record says that the
  # estimates are the values the search was started on.
  uc13_expect_fit_silent(m)
})

test_that("a covariate far from zero: the fit cannot say", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Classes of unequal size: 120, 31, 20 and 9 rows.
  d <- uc13_strata(uc13_shift)
  d$yo <- factor(findInterval(d$y, c(0.5, 1.1, 1.7)) + 1L, levels = 1:4, ordered = TRUE)
  by_ucminf <- ordinal::clmm.control(method = "ucminf")
  top <- uc13_fit(d, yo ~ z + (1 | gender:ses), control = by_ucminf)
  expect_true(top$diagnostics$converged)
  expect_lt(maihda_clmm_stationary_distance(top$model), 0.01)
  # The same model on the covariate `from` + `per` z, handed its maximum: the
  # coefficient divided by `per`, the thresholds moved by `from` times it.
  carried <- function(from, per = 12) {
    d$x <- from + per * d$z
    at <- as.numeric(top$model$optRes$par)
    slope <- length(at) - 1L
    at[slope] <- at[slope] / per
    at[seq_len(slope - 1L)] <- at[seq_len(slope - 1L)] + from * at[slope]
    uc13_fit(d, yo ~ x + (1 | gender:ses), control = by_ucminf, start = at)
  }
  stored <- function(m) {
    g <- as.numeric(m$model$gradient)
    H <- m$model$Hessian
    list(distance = sqrt(sum(g * solve(H, g))), finest = min(sqrt(diag(solve(H)))),
         coarseness = max(maihda_clmm_step_coarseness(H, solve(H))))
  }
  # In years the step is fine for every parameter: the derivatives say what is so, and
  # the fit is recorded as converged.
  in_years <- carried(45)
  expect_equal(uc13_loglik(in_years), uc13_loglik(top), tolerance = 1e-8)
  expect_lt(stored(in_years)$coarseness, 0.01)
  expect_equal(maihda_clmm_stationary_distance(in_years$model), stored(in_years)$distance)
  expect_lt(maihda_clmm_stationary_distance(in_years$model), 0.01)
  expect_true(in_years$diagnostics$converged)
  # At a birth year, and 20,000 units out, no standard error is below the step, and the
  # step along the coefficient is too coarse all the same. The package reads no
  # distance off the fit and records none.
  far <- list(carried(1975, -12), carried(20045))
  for (m in far) {
    expect_equal(uc13_loglik(m), uc13_loglik(top), tolerance = 1e-8)
    expect_equal(uc13_variance(m), uc13_variance(top), tolerance = 1e-6)
    expect_gt(stored(m)$finest, .maihda_clmm_derivative_step)
    expect_gt(stored(m)$coarseness, 2)
    uc13_expect_fit_silent(m)
  }
  # 20,000 units out the stored derivatives put that maximum well past a tenth of a
  # standard error from a stationary point.
  expect_gt(stored(far[[2]])$distance, 0.2)
})

test_that("clmm()'s general routine asked for ucminf() is read as nlminb()'s fit", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- uc13_strata(uc13_shift)
  plain <- uc13_fit(d)
  obs <- uc13_observe(fit_maihda(
    yo ~ (1 | gender:ses), data = d, engine = "ordinal",
    control = ordinal::clmm.control(method = "ucminf", useMatrix = TRUE)))
  m <- obs$res
  # The settings name ucminf; the return is nlminb()'s.
  expect_identical(m$model$control$method, "ucminf")
  expect_true(m$model$control$useMatrix)
  returned <- m$model$optRes
  expect_false(maihda_is_ucminf_return(returned))
  expect_false(is.null(returned$objective))
  expect_identical(m$diagnostics$converged, returned$convergence == 0)
  expect_true(m$diagnostics$converged)
  expect_length(m$diagnostics$messages, 0L)
  expect_equal(uc13_variance(m), uc13_variance(plain), tolerance = 1e-3)
})

test_that("at a variance of zero with ucminf() as the first optimiser", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  by_ucminf <- ordinal::clmm.control(method = "ucminf")
  flat <- uc13_strata()
  for (formula in list(yo ~ (1 | gender:ses), yo ~ age + (1 | gender:ses))) {
    m <- uc13_fit(flat, formula, control = by_ucminf)
    # The boundary fit, made by the general routine with nlminb().
    expect_lt(uc13_variance(m), 1e-8)
    expect_true(m$diagnostics$singular)
    expect_true(m$model$control$useMatrix)
    expect_false(maihda_is_ucminf_return(m$model$optRes))
    # clmm()'s own first fit by ucminf() is one the package makes again; where it is
    # read as converged, the boundary fit is converged whatever nlminb() says of a fit
    # it was started on.
    first <- suppressWarnings(ordinal::clmm(
      stats::formula(m$model), data = m$data, link = "logit", Hess = TRUE,
      control = by_ucminf))
    expect_true(maihda_is_ucminf_return(first$optRes))
    expect_true(maihda_clmm_refit_wanted(first))
    if (isTRUE(maihda_clmm_converged(first))) {
      expect_true(m$diagnostics$converged)
      expect_length(m$diagnostics$messages, 0L)
    } else {
      expect_true(isTRUE(m$diagnostics$converged) || length(m$diagnostics$messages) > 0L)
    }
  }
})

test_that("a group comparison made with ucminf() names the groups that did not converge", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- rbind(cbind(uc13_strata(uc13_shift), country = "A"),
             cbind(uc13_strata(rev(uc13_shift) * 1.2), country = "B"),
             cbind(uc13_strata(c(0.6, -0.7, 0.1, 0.9, -0.5, -0.3)), country = "C"))
  compare <- function(...) {
    uc13_observe(compare_maihda_groups(yo ~ (1 | gender:ses), data = d,
                                       group = "country", engine = "ordinal", ...))
  }
  by_nlminb <- compare()
  by_ucminf <- compare(control = ordinal::clmm.control(method = "ucminf"))
  expect_false(any(grepl("fit problems", by_nlminb$warns, fixed = TRUE)))
  # No group is named, and the table is the nlminb() one.
  expect_false(any(grepl("fit problems", by_ucminf$warns, fixed = TRUE)))
  expect_identical(by_ucminf$warns, by_nlminb$warns)
  ours <- as.data.frame(by_ucminf$res)
  theirs <- as.data.frame(by_nlminb$res)
  expect_identical(ours$status, rep("ok", 3))
  expect_identical(ours$pcv_status, rep("ok", 3))
  expect_equal(ours$var_between, theirs$var_between, tolerance = 1e-3)
  expect_equal(ours$pcv, theirs$pcv, tolerance = 1e-3)

  # With every search cut short, every group and every adjusted model is.
  cut_short <- compare(control = ordinal::clmm.control(method = "ucminf", maxeval = 5))
  quality <- cut_short$warns[grepl("fit problems", cut_short$warns, fixed = TRUE)]
  expect_identical(quality,
                   maihda_group_fit_quality_message("ordinal", character(0),
                                                    c("A", "B", "C"), c("A", "B", "C")))
})
