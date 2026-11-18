# Audit 2026-11-16, one finding from CI: the default growth fit of the bundled panel passed on a
# false alarm of lme4's convergence check on some machines and not on others.
#
# fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = maihda_long_data,
# id = "id", time = "wave") fits lmer's model and two restarts and kept lmer's own fit, the
# three being within 1e-4 of one REML criterion. lme4's gradient check compares a scaled gradient
# with 0.002, which for this model is a criterion 2e-6 above the optimum, and lmer's default
# optimiser stops between 1e-8 and 2e-4 above it: 15 to 26% of numerically equivalent fits end
# above lme4's line. On which side a fit ends follows the last bits of the criterion, and Eigen,
# lme4's linear algebra, sizes the blocks of its products by the processor's L1 data cache. On
# processors with a 48 KB L1 data cache and on Apple Silicon lmer's own fit ends 3e-6 to 6e-6
# above the optimum and lme4 warned "Model failed to converge with max|grad| = 0.00213664 (tol =
# 0.002, component 1)"; on processors with 32 KB it ends 1e-6 above it and lme4 says nothing.
# Told either size, Eigen gave either fit on one machine, in each of 17 Linux jobs on GitHub's
# runners. The estimates agree to the optimiser's precision either way, and a restart that
# lme4 checked and accepted was there for each flagged fit.
#
# FIX: maihda_fit_lmer_growth() chooses as before (maihda_pick_lowest()) and then, when lme4's
# convergence checks flag the chosen fit (maihda_lme4_check_flagged()) and lme4 checked and
# accepted another whose criterion is within 1e-4 of it (maihda_lme4_check_accepted()), keeps
# the lowest of those (maihda_pick_accepted()). A chosen fit that lme4 does not flag stays, and
# so does a flagged one without such a neighbour, which the fitter returns with lme4's warning.
# A fit lme4 has not checked (a singular one, one too large for its checks) or whose optimiser
# reports a failure takes no flagged fit's place. What fit_maihda() does with a fit the fitter
# returns flagged, and with its warning, is in test-audit-2026-11-17-grad-warning-no-restarts.R.
#
# This file draws no random numbers.
#
# A singular fit that stands at an optimum takes the place of any fit that is not singular
# before the rule above acts (maihda_pick_singular(),
# test-audit-2026-11-17-flagged-beside-singular.R). The restarts "noted singular" below carry
# lme4's note on a fit it does not call singular, so that rule leaves them alone.

gg16_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression with the text of every warning and message it raised, in order.
gg16_heard <- function(expr) {
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

# lme4's warning for a gradient above its tolerance, as it worded the one of the finding, and
# its note of a singular fit.
gg16_message <- "Model failed to converge with max|grad| = 0.00213664 (tol = 0.002, component 1)"
gg16_singular_note <- "boundary (singular) fit: see help('isSingular')"

# Whether the record of the lme4 fit `model` holds that kind of message.
gg16_gradient_flag <- function(model) {
  any(grepl("failed to converge", unlist(model@optinfo$conv$lme4$messages), fixed = TRUE))
}

# What print() shows, as one line.
gg16_printed <- function(x) gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))

# A value / conditions pair of maihda_eval_holding_conditions() as lme4 leaves a fit its
# gradient check flags: the message in the fit's record, and the warning among the conditions.
gg16_as_flagged <- function(pair) {
  pair$value@optinfo$conv$lme4 <- list(code = -1L, messages = gg16_message)
  pair$conditions <- c(pair$conditions, list(simpleWarning(gg16_message)))
  pair
}

# The pair as lme4 leaves a fit its checks accept: an empty record and no such warning.
gg16_as_accepted <- function(pair) {
  pair$value@optinfo$conv$lme4 <- list()
  pair$conditions <- Filter(function(cnd) !grepl("failed to converge", conditionMessage(cnd), fixed = TRUE),
                            pair$conditions)
  pair
}

# The default growth fit of the bundled panel with lme4's verdict on each of its three fits
# written by hand. lmer's own fit is made as fit_maihda() makes it; each restart is that same
# fit with the restart's start in its call, the same optimum reached from another start. Which
# of the real fits lme4 flags differs between machines, so none of them is left to it.
# `first`: whether lme4's check flags lmer's own fit. `restarts`: what lme4 left in the
# restarts: "accepted", "flagged", "unchecked" (no record: a fit lme4 ran no check on),
# "noted singular" (its note alone, as on a singular fit) or "failed" (an accepted record
# under an optimiser's code other than 0).
gg16_fit_with <- function(first, restarts) {
  real <- maihda_eval_holding_conditions
  own <- NULL
  gg16_heard(with_mocked_bindings(
    suppressMessages(fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                                data = MAIHDA::maihda_long_data, id = "id", time = "wave")),
    maihda_eval_holding_conditions = function(call, env) {
      if (is.null(own)) {
        own <<- gg16_as_accepted(real(call, env))
        return(if (first) gg16_as_flagged(own) else own)
      }
      restart <- own
      restart$value@call$start <- call$start
      switch(restarts,
             accepted = restart,
             flagged = gg16_as_flagged(restart),
             unchecked = {
               restart$value@optinfo$conv["lme4"] <- list(NULL)
               restart
             },
             "noted singular" = {
               restart$value@optinfo$conv$lme4 <- list(messages = gg16_singular_note)
               restart
             },
             failed = {
               restart$value@optinfo$conv$opt <- 5L
               restart
             })
    }))
}

# ---- the choice ---------------------------------------------------------------------------

test_that("a flagged choice gives way to the lowest fit within 1e-4 that lme4 accepted", {
  pick <- maihda_pick_accepted
  y <- TRUE
  n <- FALSE
  # A choice lme4 does not flag stays, whatever it says of the others.
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(n, n, n), c(y, y, y), 1L), 1L)
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(n, y, y), c(y, n, n), 1L), 1L)
  # A flagged choice with accepted fits within 1e-4, above or below: the lowest of them.
  expect_identical(pick(c(10, 10 + 5e-5, 10 - 5e-5), c(y, n, n), c(n, y, y), 1L), 3L)
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(y, n, n), c(n, y, y), 1L), 2L)
  expect_identical(pick(c(10, 10 + 5e-5, 10 - 5e-5), c(y, n, y), c(n, y, n), 1L), 2L)
  expect_identical(pick(c(10, 10 + 5e-5, 10 + 2e-5), c(y, n, n), c(n, y, y), 1L), 3L)
  # An accepted fit more than 1e-4 above is not the same optimum: the choice stays.
  expect_identical(pick(c(10, 10 + 2e-4, 11), c(y, n, n), c(n, y, y), 1L), 1L)
  # Every fit flagged: the choice stays.
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(y, y, y), c(n, n, n), 1L), 1L)
  # A fit that is neither flagged nor accepted (lme4 ran no check on it, or its optimiser
  # failed) takes no flagged fit's place, lower or not; an accepted one beside it does.
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(y, n, n), c(n, n, n), 1L), 1L)
  expect_identical(pick(c(10, 10 - 5e-5, 10 + 5e-5), c(y, n, n), c(n, n, y), 1L), 3L)
  # A restart chosen for a criterion lower by more than 1e-4, and flagged: lmer's own fit, far
  # above it, does not take its place, accepted or not; another restart within 1e-4 does.
  expect_identical(pick(c(10, 9, 9 + 5e-5), c(n, y, y), c(y, n, n), 2L), 2L)
  expect_identical(pick(c(10, 9, 9 + 5e-5), c(n, y, n), c(y, n, y), 2L), 3L)
  # A restart that failed has no criterion and is passed over.
  expect_identical(pick(c(10, NA, 10 + 5e-5), c(y, n, n), c(n, y, y), 1L), 3L)
  expect_identical(pick(c(10, NA, NA), c(y, n, n), c(n, y, y), 1L), 1L)
  # The tolerance is the caller's.
  expect_identical(pick(c(10, 10.05), c(y, n), c(n, y), 1L, tol = 0.1), 2L)
  expect_identical(pick(c(10, 10.05), c(y, n), c(n, y), 1L, tol = 0.01), 1L)
  expect_identical(formals(maihda_pick_accepted)$tol, 1e-4)
  # A choice without a criterion, or no choice: returned as given.
  expect_identical(pick(c(NA, 12, 11), c(y, n, n), c(n, y, y), 1L), 1L)
  expect_identical(pick(numeric(0), logical(0), logical(0), NA_integer_), NA_integer_)
  # A fit whose verdict is not known takes no flagged fit's place, and keeps its own.
  expect_identical(pick(c(10, 10 + 5e-5, 10 + 6e-5), c(y, n, n), c(n, NA, y), 1L), 3L)
  expect_identical(pick(c(10, 10 + 5e-5), c(y, n), c(n, NA), 1L), 1L)
  expect_identical(pick(c(10, 10 + 5e-5), c(NA, n), c(n, y), 1L), 1L)
})

# ---- lme4's verdict, read off the fit -------------------------------------------------------

test_that("lme4's convergence checks are read off the record lme4 leaves in a fit", {
  skip_on_cran()
  sleep <- lme4::sleepstudy
  verdict <- function(model) c(flagged = maihda_lme4_check_flagged(model),
                               accepted = maihda_lme4_check_accepted(model))
  checked <- gg16_quiet(lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleep))
  expect_identical(checked@optinfo$conv$lme4, list())
  expect_identical(verdict(checked), c(flagged = FALSE, accepted = TRUE))
  # lme4's own gradient message: a fit stopped nine evaluations in by a step tolerance of 0.3,
  # with the optimiser's code at 0, as in the finding.
  loose <- gg16_quiet(lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleep,
                                 control = lme4::lmerControl(optCtrl = list(xtol_rel = 0.3))))
  expect_true(loose@optinfo$conv$opt == 0)
  expect_match(unlist(loose@optinfo$conv$lme4$messages), "Model failed to converge with max|grad|",
               fixed = TRUE)
  expect_identical(verdict(loose), c(flagged = TRUE, accepted = FALSE))
  # Its Hessian messages: a fit stopped after three evaluations, and each message alone.
  short <- gg16_quiet(lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleep,
                                 control = lme4::lmerControl(optCtrl = list(maxeval = 3))))
  expect_match(paste(unlist(short@optinfo$conv$lme4$messages), collapse = " "), "Hessian",
               fixed = TRUE)
  expect_identical(verdict(short), c(flagged = TRUE, accepted = FALSE))
  for (text in c("Model failed to converge: degenerate  Hessian with 1 negative eigenvalues",
                 "Model is nearly unidentifiable: very large eigenvalue\n - Rescale variables?")) {
    hessian <- checked
    hessian@optinfo$conv$lme4 <- list(code = -3L, messages = text)
    expect_identical(verdict(hessian), c(flagged = TRUE, accepted = FALSE))
  }
  # A singular fit is not flagged, and not accepted either: lme4 leaves its note and runs no
  # gradient or Hessian check on it. Ten groups with one mean.
  flat <- data.frame(y = rep(c(1, 2, 3, 4), 10), g = factor(rep(1:10, each = 4)))
  singular <- gg16_quiet(lme4::lmer(y ~ 1 + (1 | g), data = flat))
  expect_identical(unlist(singular@optinfo$conv$lme4$messages), gg16_singular_note)
  expect_identical(verdict(singular), c(flagged = FALSE, accepted = FALSE))
  # The same without lme4's note, as with its singular check switched off.
  unnoted <- singular
  unnoted@optinfo$conv$lme4 <- list()
  expect_identical(verdict(unnoted), c(flagged = FALSE, accepted = FALSE))
  # A singular fit that a check flags as well is flagged.
  both <- singular
  both@optinfo$conv$lme4$messages <- c(gg16_singular_note, gg16_message)
  expect_identical(verdict(both), c(flagged = TRUE, accepted = FALSE))
  # A fit lme4 ran no check on (here: asked for no derivatives) has no record.
  unchecked <- gg16_quiet(lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleep,
                                     control = lme4::lmerControl(calc.derivs = FALSE)))
  expect_null(unchecked@optinfo$conv$lme4)
  expect_identical(verdict(unchecked), c(flagged = FALSE, accepted = FALSE))
  # An accepted record under an optimiser's code other than 0 is not accepted.
  failed <- checked
  failed@optinfo$conv$opt <- 5L
  expect_identical(verdict(failed), c(flagged = FALSE, accepted = FALSE))
  # The records written by hand for the blocks below read as lme4's own.
  pair <- list(value = checked, conditions = list())
  expect_identical(verdict(gg16_as_flagged(pair)$value), c(flagged = TRUE, accepted = FALSE))
  expect_identical(verdict(gg16_as_accepted(gg16_as_flagged(pair))$value),
                   c(flagged = FALSE, accepted = TRUE))
  expect_length(gg16_as_accepted(gg16_as_flagged(pair))$conditions, 0L)
})

test_that("anything that is not an lme4 fit is neither flagged nor accepted", {
  for (x in list(NULL, list(crit = 1), stats::lm(dist ~ speed, data = datasets::cars))) {
    expect_false(maihda_lme4_check_flagged(x))
    expect_false(maihda_lme4_check_accepted(x))
  }
})

# ---- the fitter -----------------------------------------------------------------------------

test_that("the growth fitter keeps a restart lme4 accepted in place of a flagged fit beside it", {
  env <- new.env()
  # lmer's own fit and two restarts as the criteria and lme4's verdicts of `candidates`, one
  # row per fit; a fit that lme4 flags warns as lme4 does, and every fit says which it is.
  env$fake_fit <- function(start = NULL) {
    k <- if (is.null(start)) 1L else as.integer(start[["k"]])
    row <- env$candidates[k, ]
    if (row$flag) warning(sprintf("lme4 flags fit %d", k))
    message(sprintf("fit %d", k))
    list(crit = row$crit, flag = row$flag, ok = row$ok, k = k)
  }
  local_mocked_bindings(
    maihda_growth_restart_starts = function(theta, cnms, id, scale = 0.5) list(c(k = 2), c(k = 3)),
    maihda_growth_centred_start = function(fit_call, fit_env, time, time_degree) list(c(k = 2)),
    maihda_growth_diagonal_start = function(theta, cnms, id, scale = 0.5) list(c(k = 3)),
    maihda_lmer_criterion = function(model) model$crit,
    maihda_lme4_check_flagged = function(model) isTRUE(model$flag),
    maihda_lme4_check_accepted = function(model) isTRUE(model$ok)
  )
  kept <- function(crit, flag, ok = !flag, ...) {
    env$candidates <- data.frame(crit = crit, flag = flag, ok = ok)
    gg16_heard(maihda_fit_lmer_growth(quote(fake_fit()), env, id = "pid", ...))
  }
  # lmer's own fit flagged, both restarts accepted within 1e-4: the lower restart, and only its
  # own conditions; lme4's warning about lmer's fit is not passed on.
  h <- kept(c(10, 10 + 3e-5, 10 - 2e-5), c(TRUE, FALSE, FALSE))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  h <- kept(c(10, 10 + 3e-5, 10 + 6e-5), c(TRUE, FALSE, FALSE))
  expect_identical(h$value$k, 2L)
  expect_identical(h$seen, "message: fit 2\n")
  # One restart accepted, the other flagged too.
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(TRUE, TRUE, FALSE))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  # No restart accepted: lmer's own fit, with lme4's warning. So too when the restarts are not
  # flagged but were not checked either, as fits too large for lme4's checks are.
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(TRUE, TRUE, TRUE))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, c("warning: lme4 flags fit 1", "message: fit 1\n"))
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(TRUE, FALSE, FALSE), ok = c(FALSE, FALSE, FALSE))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, c("warning: lme4 flags fit 1", "message: fit 1\n"))
  h <- kept(c(10, 10 - 2e-5, 10 + 3e-5), c(TRUE, FALSE, FALSE), ok = c(FALSE, FALSE, TRUE))
  expect_identical(h$value$k, 3L)
  # The accepted restart is another optimum, more than 1e-4 above: lmer's own fit, with the
  # warning.
  h <- kept(c(10, 10 + 1e-3, 12), c(TRUE, FALSE, FALSE))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, c("warning: lme4 flags fit 1", "message: fit 1\n"))
  # lmer's own fit accepted: it stays, beside lower restarts within 1e-4, flagged or not.
  h <- kept(c(10, 10 - 2e-5, 10 - 5e-5), c(FALSE, FALSE, TRUE))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, "message: fit 1\n")
  # A restart lower by more than 1e-4 is kept, flagged or not, when the other restart is not
  # beside it; flagged, it gives way to the other restart within 1e-4 that lme4 accepted.
  h <- kept(c(10, 9, 9.5), c(FALSE, TRUE, FALSE))
  expect_identical(h$value$k, 2L)
  expect_identical(h$seen, c("warning: lme4 flags fit 2", "message: fit 2\n"))
  h <- kept(c(10, 9, 9 + 4e-5), c(FALSE, TRUE, FALSE))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  # The fitter's own tolerance is the rule's: 0.005 above is beside the flagged fit at 0.01
  # and another optimum at the default.
  h <- kept(c(10, 10 + 5e-3, 12), c(TRUE, FALSE, FALSE), tol = 0.01)
  expect_identical(h$value$k, 2L)
  h <- kept(c(10, 10 + 5e-3, 12), c(TRUE, FALSE, FALSE))
  expect_identical(h$value$k, 1L)
  # Quadratic growth takes the same rule: the centred-time restart beside lmer's flagged fit,
  # and the diagonal start, made when those two differ by more than 0.01.
  h <- kept(c(10, 10 + 3e-5, NA), c(TRUE, FALSE, FALSE), time = "wave", time_degree = 2)
  expect_identical(h$value$k, 2L)
  expect_identical(h$seen, "message: fit 2\n")
  h <- kept(c(10, 10.5, 10 - 2e-5), c(TRUE, FALSE, FALSE), time = "wave", time_degree = 2)
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  h <- kept(c(10, 10.5, 10 - 2e-5), c(TRUE, FALSE, TRUE), time = "wave", time_degree = 2)
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, c("warning: lme4 flags fit 1", "message: fit 1\n"))
  # A fit handed in as the first one, as a bootstrap draw's refit is, is treated as lmer's own.
  env$candidates <- data.frame(crit = c(10, 10 + 3e-5, 10 - 2e-5), flag = c(TRUE, FALSE, FALSE),
                               ok = c(FALSE, TRUE, TRUE))
  first <- list(value = list(crit = 10, flag = TRUE, ok = FALSE, k = 1L),
                conditions = list(simpleWarning("lme4 flags the fit handed in")))
  h <- gg16_heard(maihda_fit_lmer_growth(quote(fake_fit()), env, "pid", first = first))
  expect_identical(h$value$k, 3L)
  expect_identical(h$seen, "message: fit 3\n")
  env$candidates$flag <- c(TRUE, TRUE, TRUE)
  env$candidates$ok <- c(FALSE, FALSE, FALSE)
  h <- gg16_heard(maihda_fit_lmer_growth(quote(fake_fit()), env, "pid", first = first))
  expect_identical(h$value$k, 1L)
  expect_identical(h$seen, "warning: lme4 flags the fit handed in")
})

# ---- fit_maihda() on the bundled panel --------------------------------------------------------

test_that("fit_maihda() returns a restart lme4 accepted when lme4 flags lmer's own growth fit", {
  skip_on_cran()
  h <- gg16_fit_with(first = TRUE, restarts = "accepted")
  m <- h$value
  # Nothing of lme4's warning reaches the caller, the fit or its record.
  expect_identical(h$seen, character(0))
  expect_false(gg16_gradient_flag(m$model))
  expect_true(is.numeric(stats::getCall(m$model)$start))
  expect_true(m$diagnostics$converged)
  expect_null(m$diagnostics$lme4_overruled)
  expect_null(m$diagnostics$stationary_distance)
  expect_false(grepl("lme4 warned|failed to converge", gg16_printed(m)))
})

# The blocks below pin what the growth fitter leaves alone. Each later fit of the call, the one
# fit_maihda() makes from a flagged fit's estimates included, is written as its restarts are.

test_that("a growth fit flagged from every start keeps lmer's own fit and its record", {
  skip_on_cran()
  h <- gg16_fit_with(first = TRUE, restarts = "flagged")
  m <- h$value
  # The record calls the fit converged, so lme4's warning is not passed on.
  expect_identical(h$seen, character(0))
  expect_true(gg16_gradient_flag(m$model))
  expect_null(stats::getCall(m$model)$start)
  # The record measures the fit and overrules lme4, and print() says so.
  expect_true(m$diagnostics$converged)
  expect_identical(m$diagnostics$lme4_overruled, gg16_message)
  expect_lt(m$diagnostics$stationary_distance, 0.1)
  expect_match(gg16_printed(m),
               "lme4 warned that the model failed to converge; in centred and scaled units the fit is",
               fixed = TRUE)
})

test_that("restarts lme4 did not check, or whose optimiser failed, take no flagged fit's place", {
  skip_on_cran()
  for (left in c("unchecked", "noted singular", "failed")) {
    h <- gg16_fit_with(first = TRUE, restarts = left)
    expect_identical(h$seen, character(0))
    expect_true(gg16_gradient_flag(h$value$model))
    expect_null(stats::getCall(h$value$model)$start)
    expect_identical(h$value$diagnostics$lme4_overruled, gg16_message)
  }
})

test_that("a growth fit lme4 accepts from every start is lmer's own, as before", {
  skip_on_cran()
  h <- gg16_fit_with(first = FALSE, restarts = "accepted")
  m <- h$value
  expect_identical(h$seen, character(0))
  expect_null(stats::getCall(m$model)$start)
  expect_false(gg16_gradient_flag(m$model))
  expect_null(m$diagnostics$lme4_overruled)
})
