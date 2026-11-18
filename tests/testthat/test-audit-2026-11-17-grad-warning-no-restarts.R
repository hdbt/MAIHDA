# Audit 2026-11-17, one finding left open by the pass of 2026-11-16: fit_maihda() passed on
# lme4's "Model failed to converge with max|grad| = ... (tol = 0.002, component 1)" for fits
# whose own record then called them converged, wherever it had no other fit to offer.
#
# CONFIRMED for the three kinds of fit named. (a) A linear growth model with stratum_slope =
# FALSE is fitted once: lme4 flagged 13 of 600 such fits on 576 rows (16 strata of 6 persons, 6
# waves), 28 of 300 numerically equivalent copies of the bundled panel, 23 of 200 on 3,000 rows
# and 7 of 60 on 8,000, every one 1.9e-6 to 2.0e-5 above its REML optimum, with the warning
# passed on and the record converged. (b) Quadratic growth without a stratum slope came back
# flagged from every start in 26 of 120 fits, 14 more beside a singular fit. (c) Outside growth
# models lme4 flagged no lmer fit without a random slope (0 of 1,500), but the logistic model
# with the strata's main effects on the bundled maihda_health_data in 12 of 13 orders of its
# rows (lme4 2.1-0), 5e-6 to 1.2e-4 above its maximum, and, at the optimum itself, fits with a
# covariate far from zero or in large units (a random slope on an age in years: 52 of 100).
# What stops lmer short of lme4's own line is nloptr's default xtol_rel of 1e-4 (10 of 13
# flagged fits) and lme4's ftol_abs of 1e-8 (3 of 13).
#
# Per the user, (1) a fit that lme4's convergence checks flag is refitted once from its
# estimates, an lmer fit with xtol_rel 1e-5 and ftol_abs 1e-10 unless the call has a control of
# the user's own, a glmer fit by its second stage, and the refit is returned when lme4 checks
# and accepts it (maihda_lme4_accepted_refit(), maihda_lme4_flagged_refit()); and (2) when the
# fit returned is still flagged and its record calls it converged, lme4's warning is not passed
# on: the record, print() and summary() keep it (maihda_lme4_conditions_passed_on()). A fit the
# record does not call converged keeps its warning.
#
# Which real fit lme4 flags differs between machines, so lme4's verdict on each fit of the
# fit_maihda() blocks is written by hand. This file draws no random numbers.

gw17_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The value of an expression with the text of every warning and message it raised, in order.
gw17_heard <- function(expr) {
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

# lme4's messages as it words them: a gradient above its tolerance (the finding's), a large
# eigenvalue of the Hessian and a large ratio of its eigenvalues, and its note of a singular
# fit.
gw17_gradient <- "Model failed to converge with max|grad| = 0.00319134 (tol = 0.002, component 1)"
gw17_eigenvalue <- "Model is nearly unidentifiable: very large eigenvalue\n - Rescale variables?"
gw17_ratio <- "Model is nearly unidentifiable: large eigenvalue ratio\n - Rescale variables?"
gw17_singular_note <- "boundary (singular) fit: see help('isSingular')"

# What print() shows, as one line.
gw17_printed <- function(x) gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))

# The warnings lme4 raises for the messages `said` of its checks: one for each message of the
# gradient check, and one for all those of the Hessian check, joined by ";".
gw17_warnings <- function(said) {
  hessian <- grepl("eigenvalue|Hessian", said)
  texts <- c(said[!hessian], if (any(hessian)) paste(said[hessian], collapse = ";"))
  lapply(texts, simpleWarning)
}

# Whether a condition is one of lme4's convergence-check warnings as written here.
gw17_is_check <- function(cnd) {
  inherits(cnd, "warning") &&
    grepl("Model failed to converge|Model is nearly unidentifiable", conditionMessage(cnd))
}

# A value / conditions pair as lme4 leaves a fit its checks flag: the messages in the fit's
# record, and the warnings among the conditions.
gw17_as_flagged <- function(pair, said = gw17_gradient) {
  pair$value@optinfo$conv$lme4 <- list(code = -1L, messages = said)
  pair$conditions <- c(pair$conditions, gw17_warnings(said))
  pair
}

# The pair as lme4 leaves a fit its checks accept: an empty record and no such warning.
gw17_as_accepted <- function(pair) {
  pair$value@optinfo$conv$lme4 <- list()
  pair$conditions <- Filter(Negate(gw17_is_check), pair$conditions)
  pair
}

# 120 persons of the bundled panel (600 rows, 12 strata): a growth fit without a stratum
# slope that is neither singular nor near a bound.
gw17_panel <- function() {
  d <- MAIHDA::maihda_long_data
  d[d$id %in% unique(d$id)[1:120], ]
}

# Six strata of 30 rows with a 0 / 1 outcome and a covariate on a scale of about 1.
gw17_binary <- function() {
  years <- c(53, 49, 55, 51, 55, 43, 46, 55, 71, 56, 19, 53, 23, 41, 63, 45, 53, 54, 54, 55,
             38, 61, 23, 36, 55, 47, 51, 53, 41, 66)
  g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"), stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 30L), , drop = FALSE]
  d$z <- (rep(years, times = nrow(g)) - 45) / 12
  latent <- rep(c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2), each = 30L) +
    rep(seq(-1.45, 1.45, length.out = 30L), times = nrow(g)) + 0.4 * d$z
  d$yb <- as.integer(latent > 0)
  rownames(d) <- NULL
  d
}

# fit_maihda() with lme4's verdicts written by hand. The first fit is made as fit_maihda()
# makes it and then reads as flagged (`first`, the messages lme4 left) or accepted (NULL); the
# refit from its estimates, the only later fit of these models, is made for real and then reads
# as `refit` says: "accepted", "flagged", "unchecked" (no record), "noted singular" (lme4's
# note alone) or "failed" (an accepted record under an optimiser's code other than 0).
# Returned with what was heard: `refits`, the number of refits made, and `first`, the first
# fit as it was made.
gw17_fit_with <- function(fit, first = gw17_gradient, refit = "accepted") {
  real_hold <- maihda_hold_fit_conditions
  real_eval <- maihda_eval_holding_conditions
  refits <- 0L
  made <- NULL
  heard <- gw17_heard(with_mocked_bindings(
    suppressMessages(fit()),
    maihda_hold_fit_conditions = function(expr) {
      pair <- gw17_as_accepted(real_hold(expr))
      made <<- pair$value
      if (is.null(first)) pair else gw17_as_flagged(pair, first)
    },
    maihda_eval_holding_conditions = function(call, env) {
      pair <- real_eval(call, env)
      if (is.null(call$start)) {
        return(pair)
      }
      refits <<- refits + 1L
      pair <- gw17_as_accepted(pair)
      switch(refit,
             accepted = pair,
             flagged = gw17_as_flagged(pair),
             unchecked = {
               pair$value@optinfo$conv["lme4"] <- list(NULL)
               pair
             },
             "noted singular" = {
               pair$value@optinfo$conv$lme4 <- list(messages = gw17_singular_note)
               pair
             },
             failed = {
               pair$value@optinfo$conv$opt <- 5L
               pair
             })
    }))
  c(heard, list(refits = refits, first = made))
}

gw17_growth <- function(...) {
  fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = gw17_panel(),
             id = "id", time = "wave", stratum_slope = FALSE, ...)
}
gw17_logistic <- function(...) {
  fit_maihda(yb ~ z + (1 | gender:ses), data = gw17_binary(), family = "binomial", ...)
}

# The stop of an lmer refit, as its call carries it.
gw17_stop <- quote(lme4::lmerControl(optCtrl = list(xtol_rel = 1e-5, ftol_abs = 1e-10)))

# ---- holding what a fit raises ----------------------------------------------------------------

test_that("what a fit raises is held in order, and released when the fit stops with an error", {
  held <- gw17_heard(maihda_hold_fit_conditions({
    message("one")
    warning("two")
    message("three")
    42
  }))
  expect_identical(held$seen, character(0))
  expect_identical(held$value$value, 42)
  expect_identical(vapply(held$value$conditions, conditionMessage, character(1)),
                   c("one\n", "two", "three\n"))
  expect_identical(vapply(held$value$conditions, function(cnd) inherits(cnd, "warning"), logical(1)),
                   c(FALSE, TRUE, FALSE))
  # The fit is evaluated where it is written.
  made <- FALSE
  pair <- maihda_hold_fit_conditions({
    made <- TRUE
    "fit"
  })
  expect_true(made)
  expect_identical(pair, list(value = "fit", conditions = list()))
  # An error: what the fit had raised until then is signalled, then the error passes.
  stopped <- gw17_heard(tryCatch(
    maihda_hold_fit_conditions({
      warning("before the error")
      message("a note")
      stop("the fit failed")
    }),
    error = function(e) conditionMessage(e)))
  expect_identical(stopped$value, "the fit failed")
  expect_identical(stopped$seen, c("warning: before the error", "message: a note\n"))
  # An error caught inside the fit, as a restart's is, releases nothing.
  inner <- gw17_heard(maihda_hold_fit_conditions({
    warning("kept back")
    tryCatch(stop("a restart failed"), error = function(e) NULL)
    "fit"
  }))
  expect_identical(inner$seen, character(0))
  expect_identical(inner$value$value, "fit")
  expect_length(inner$value$conditions, 1L)
})

# ---- which warnings are passed on -------------------------------------------------------------

test_that("lme4's convergence warnings are passed on unless the record calls the fit converged", {
  passed <- function(conditions, record) {
    vapply(maihda_lme4_conditions_passed_on(conditions, record), conditionMessage, character(1))
  }
  raised <- list(simpleMessage("a note\n"), simpleWarning(gw17_gradient),
                 simpleWarning("another warning"), simpleWarning(gw17_eigenvalue))
  all_of_them <- c("a note\n", gw17_gradient, "another warning", gw17_eigenvalue)
  # The record overrules both of lme4's messages: its two warnings go, the rest stays in order.
  both <- list(converged = TRUE, lme4_overruled = c(gw17_gradient, gw17_eigenvalue))
  expect_identical(passed(raised, both), c("a note\n", "another warning"))
  # lme4 leaves a list when two checks fire.
  expect_identical(passed(raised, list(converged = TRUE,
                                       lme4_overruled = list(gw17_gradient, gw17_eigenvalue))),
                   c("a note\n", "another warning"))
  # Only what the record holds as overruled goes.
  expect_identical(passed(raised, list(converged = TRUE, lme4_overruled = gw17_gradient)),
                   c("a note\n", "another warning", gw17_eigenvalue))
  # lme4's Hessian check raises one warning for all of its messages, joined by ";": it goes
  # when the record overrules them, beside the gradient check's or alone, and stays when the
  # record overrules only one of them.
  joined <- paste(gw17_eigenvalue, gw17_ratio, sep = ";")
  hessian <- list(simpleWarning(gw17_gradient), simpleWarning(joined), simpleWarning("another"))
  expect_identical(passed(hessian, list(converged = TRUE,
                                        lme4_overruled = c(gw17_gradient, gw17_eigenvalue, gw17_ratio))),
                   "another")
  expect_identical(passed(hessian[2:3], list(converged = TRUE,
                                             lme4_overruled = c(gw17_eigenvalue, gw17_ratio))),
                   "another")
  expect_identical(passed(hessian[2:3], list(converged = TRUE, lme4_overruled = gw17_eigenvalue)),
                   c(joined, "another"))
  # Messages joined in another order than lme4 left them are not lme4's warning.
  expect_identical(passed(list(simpleWarning(paste(gw17_ratio, gw17_eigenvalue, sep = ";"))),
                          list(converged = TRUE, lme4_overruled = c(gw17_eigenvalue, gw17_ratio))),
                   paste(gw17_ratio, gw17_eigenvalue, sep = ";"))
  # A fit the record does not call converged keeps every warning: one it did not overrule,
  # one it could not measure, and one it overruled and then found at a boundary trap.
  expect_identical(passed(raised, list(converged = FALSE, messages = gw17_gradient)), all_of_them)
  expect_identical(passed(raised, list(converged = NA)), all_of_them)
  expect_identical(passed(raised, list(converged = FALSE, lme4_overruled = gw17_gradient)),
                   all_of_them)
  # A converged record that overrules nothing passes everything on.
  expect_identical(passed(raised, list(converged = TRUE)), all_of_them)
  expect_identical(passed(raised, list(converged = TRUE, lme4_overruled = character(0))),
                   all_of_them)
  # A message with the text of an overruled warning is not a warning and stays; a warning that
  # differs from lme4's by its value stays; space around the text does not count.
  expect_identical(passed(list(simpleMessage(gw17_gradient)), both), gw17_gradient)
  other <- sub("0.00319134", "0.00213664", gw17_gradient, fixed = TRUE)
  expect_identical(passed(list(simpleWarning(other)), both), other)
  expect_identical(passed(list(simpleWarning(paste0(gw17_gradient, "\n"))), both), character(0))
  # The condition objects are passed on as they were raised, and nothing gives nothing.
  expect_identical(maihda_lme4_conditions_passed_on(raised, both), raised[c(1L, 3L)])
  expect_null(maihda_lme4_conditions_passed_on(NULL, both))
  expect_identical(maihda_lme4_conditions_passed_on(list(), both), list())
})

# ---- the refit --------------------------------------------------------------------------------

test_that("an lmer fit lme4 flags is refitted from its estimates under a tighter stop", {
  skip_on_cran()
  env <- new.env()
  env$sleep <- lme4::sleepstudy
  call <- quote(lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleep))
  # A fit stopped nine evaluations in by a step tolerance of 0.3: lme4's gradient check flags
  # it, and its optimiser reports no failure.
  loose_call <- call
  loose_call$control <- quote(lme4::lmerControl(optCtrl = list(xtol_rel = 0.3)))
  loose <- gw17_quiet(eval(loose_call, env))
  expect_true(maihda_lme4_check_flagged(loose))
  expect_true(loose@optinfo$conv$opt == 0)
  best <- gw17_quiet(eval(call, env))
  expect_gt(lme4::REMLcrit(loose) - lme4::REMLcrit(best), 0.01)

  refit <- maihda_lme4_flagged_refit(loose, call, env)
  expect_named(refit, c("value", "conditions"))
  again <- refit$value
  # It starts at the fit's variance parameters and carries the tighter stop in its call.
  expect_identical(stats::getCall(again)$start, unname(lme4::getME(loose, "theta")))
  expect_identical(stats::getCall(again)$control, gw17_stop)
  expect_identical(again@optinfo$control[c("xtol_rel", "ftol_abs")],
                   list(xtol_rel = 1e-5, ftol_abs = 1e-10))
  expect_identical(.maihda_flagged_refit_control, gw17_stop)
  # It ends at the optimum, where lmer's own fit ends, and lme4 accepts it.
  expect_lt(abs(lme4::REMLcrit(again) - lme4::REMLcrit(best)), 1e-5)
  expect_true(maihda_lme4_check_accepted(again))
  expect_length(refit$conditions, 0L)

  # A control of the user's own in the call is the refit's too: the same loose stop.
  env$.maihda_arg_control <- lme4::lmerControl(optCtrl = list(xtol_rel = 0.3))
  user_call <- call
  user_call$control <- quote(.maihda_arg_control)
  user <- maihda_lme4_flagged_refit(loose, user_call, env)
  expect_identical(stats::getCall(user$value)$control, quote(.maihda_arg_control))
  expect_identical(user$value@optinfo$control$xtol_rel, 0.3)
  expect_identical(stats::getCall(user$value)$start, unname(lme4::getME(loose, "theta")))
  expect_lte(lme4::REMLcrit(user$value), lme4::REMLcrit(loose))
  # So is one under a shorter spelling of the argument: the call gets no second control (with
  # two, lmer() would stop), and lmer() stores the user's under the full name.
  short_call <- call
  short_call$contro <- quote(.maihda_arg_control)
  short <- maihda_lme4_flagged_refit(loose, short_call, env)
  expect_named(short, c("value", "conditions"))
  expect_identical(stats::getCall(short$value)$control, quote(.maihda_arg_control))
  expect_identical(short$value@optinfo$control$xtol_rel, 0.3)

  # A refit that ends more than `tol` above the fit is not offered: none can end 1 below the
  # optimum, and at the default any refit that does not rise is.
  expect_null(maihda_lme4_flagged_refit(best, call, env, tol = -1))
  expect_named(maihda_lme4_flagged_refit(best, call, env), c("value", "conditions"))
  expect_identical(formals(maihda_lme4_flagged_refit)$tol, 1e-4)
  # A maximum-likelihood fit is compared by its deviance.
  ml_call <- call
  ml_call$REML <- FALSE
  ml <- gw17_quiet(eval(ml_call, env))
  ml_again <- maihda_lme4_flagged_refit(ml, ml_call, env)$value
  expect_false(lme4::isREML(ml_again))
  expect_lt(abs(stats::deviance(ml_again) - stats::deviance(ml)), 1e-5)

  # A glmer.nb fit is not refitted, nor a fit of any other function.
  expect_null(maihda_lme4_flagged_refit(best, quote(lme4::glmer.nb(y ~ 1 + (1 | g), data = d)), env))
  expect_null(maihda_lme4_flagged_refit(best, quote(brms::brm(y ~ 1 + (1 | g), data = d)), env))
})

test_that("a glmer fit lme4 flags is refitted by glmer's second stage from its estimates", {
  skip_on_cran()
  env <- new.env()
  env$d <- gw17_binary()
  call <- quote(lme4::glmer(yb ~ z + (1 | gender:ses), data = d, family = stats::binomial()))
  fit <- gw17_quiet(eval(call, env))
  expect_identical(fit@optinfo$optimizer, "Nelder_Mead")
  # The last pass of a refitted fit leaves a fit of that optimiser alone ...
  expect_null(maihda_glmer_finish(list(value = fit, conditions = list()), call, env))
  # ... and the refit of a flagged fit runs it once more, from the fit's estimates.
  refit <- maihda_lme4_flagged_refit(fit, call, env)
  expect_named(refit, c("value", "conditions"))
  again <- refit$value
  start <- stats::getCall(again)$start
  expect_identical(start$theta, unname(lme4::getME(fit, "theta")))
  expect_identical(start$fixef, lme4::fixef(fit))
  expect_identical(stats::getCall(again)$control, quote(.maihda_arg_finish_control))
  expect_false(get(".maihda_arg_finish_control", envir = env)$nAGQ0initStep)
  expect_identical(again@optinfo$optimizer, "Nelder_Mead")
  expect_true(again@optinfo$conv$opt == 0)
  expect_lte(maihda_glmer_rise(again, fit), 1e-4)
  # lme4's own log-likelihoods of two fits at one point differ by the precision of its
  # penalised least squares (1.3e-4 here before lme4 2.1-0).
  expect_lt(abs(as.numeric(stats::logLik(again)) - as.numeric(stats::logLik(fit))), 1e-3)

  # That pass is asked for a glmer() fit and not for a glmer.nb() fit, whose dispersion it
  # would search again.
  asked <- character(0)
  local_mocked_bindings(maihda_glmer_finish = function(kept, fit_call, fit_env, tol = 1e-4,
                                                       again = FALSE) {
    asked <<- c(asked, paste(deparse(fit_call[[1L]]), again))
    NULL
  })
  nb_call <- call
  nb_call[[1L]] <- quote(lme4::glmer.nb)
  expect_null(maihda_lme4_flagged_refit(fit, nb_call, env))
  expect_identical(asked, character(0))
  expect_null(maihda_lme4_flagged_refit(fit, call, env))
  expect_identical(asked, "lme4::glmer TRUE")
})

test_that("a flagged fit gives way to its refit only when lme4 checked and accepted the refit", {
  env <- new.env()
  env$made <- 0L
  # The fit's environment, where the refit of a glmer fit binds its control.
  fit_env <- new.env()
  control_name <- ".maihda_arg_finish_control"
  local_mocked_bindings(
    maihda_lme4_check_flagged = function(model) isTRUE(model$flag),
    maihda_lme4_check_accepted = function(model) isTRUE(model$ok),
    maihda_lme4_flagged_refit = function(model, fit_call, fit_env, tol = 1e-4) {
      env$made <- env$made + 1L
      env$tol <- tol
      assign(control_name, "the refit's", envir = fit_env)
      if (is.function(env$refit)) env$refit() else env$refit
    })
  first <- list(value = list(flag = TRUE, k = 1L),
                conditions = list(simpleWarning("lme4 flags fit 1")))
  kept <- function(refit, fitted = first, ...) {
    env$refit <- refit
    maihda_lme4_accepted_refit(fitted, quote(fit()), fit_env, ...)
  }
  # A refit lme4 accepted is returned, with what it raised and nothing of the first fit, and
  # its control stays bound for its call.
  accepted <- list(value = list(ok = TRUE, k = 2L), conditions = list(simpleMessage("fit 2\n")))
  expect_identical(kept(accepted), accepted)
  expect_identical(env$made, 1L)
  expect_identical(get(control_name, envir = fit_env), "the refit's")
  # A refit lme4 flagged too, one it did not check, none at all, and one that failed: the
  # first fit, with its warning, and the fit's environment as it was: a control bound there
  # for the fit's own call is that control again ...
  not_returned <- list(list(value = list(flag = TRUE, ok = FALSE), conditions = list()),
                       list(value = list(k = 3L), conditions = list()),
                       NULL,
                       function() stop("the refit failed"))
  for (refit in not_returned) {
    assign(control_name, "the fit's", envir = fit_env)
    expect_identical(kept(refit), first)
    expect_identical(get(control_name, envir = fit_env), "the fit's")
  }
  # ... and where none was bound, none is left.
  for (refit in not_returned) {
    rm(list = control_name, envir = fit_env)
    expect_identical(kept(refit), first)
    expect_false(exists(control_name, envir = fit_env, inherits = FALSE))
    assign(control_name, "the fit's", envir = fit_env)
  }
  # A refit that is returned keeps its own control bound, in place of the first fit's.
  expect_identical(kept(accepted), accepted)
  expect_identical(get(control_name, envir = fit_env), "the refit's")
  # A fit lme4 does not flag is returned as it is, and no refit is made.
  env$made <- 0L
  clean <- list(value = list(flag = FALSE, k = 1L), conditions = list(simpleMessage("fit 1\n")))
  expect_identical(kept(accepted, fitted = clean), clean)
  expect_identical(env$made, 0L)
  # The tolerance is the caller's.
  kept(accepted)
  expect_identical(env$tol, 1e-4)
  kept(accepted, tol = 0.01)
  expect_identical(env$tol, 0.01)
})

# ---- fit_maihda(): a growth model fitted once ---------------------------------------------------

test_that("fit_maihda() returns the refit of a flagged growth fit when lme4 accepts it", {
  skip_on_cran()
  h <- gw17_fit_with(gw17_growth, refit = "accepted")
  m <- h$value
  expect_identical(h$refits, 1L)
  # Nothing of lme4's warning reaches the caller, the fit or its record.
  expect_identical(h$seen, character(0))
  expect_false(maihda_lme4_check_flagged(m$model))
  expect_true(m$diagnostics$converged)
  expect_null(m$diagnostics$lme4_overruled)
  expect_null(m$diagnostics$stationary_distance)
  expect_false(grepl("lme4 warned|failed to converge", gw17_printed(m)))
  # The fit is the refit: started at the first fit's estimates, under the tighter stop, and at
  # the optimum the first fit is at.
  expect_identical(stats::getCall(m$model)$start, unname(lme4::getME(h$first, "theta")))
  expect_identical(stats::getCall(m$model)$control, gw17_stop)
  expect_lt(abs(lme4::REMLcrit(m$model) - lme4::REMLcrit(h$first)), 1e-4)
  expect_identical(m$longitudinal_info$stratum_slope, FALSE)
})

test_that("a flagged growth fit whose refit lme4 does not accept stays, without lme4's warning", {
  skip_on_cran()
  for (left in c("flagged", "unchecked", "noted singular", "failed")) {
    h <- gw17_fit_with(gw17_growth, refit = left)
    m <- h$value
    expect_identical(h$refits, 1L)
    # lmer's own fit, as lme4 left it ...
    expect_null(stats::getCall(m$model)$start)
    expect_null(stats::getCall(m$model)$control)
    expect_identical(unlist(m$model@optinfo$conv$lme4$messages), gw17_gradient)
    # ... which the record measures and calls converged, so lme4's warning is not passed on;
    # the record and print() keep it.
    expect_identical(h$seen, character(0))
    expect_true(m$diagnostics$converged)
    expect_identical(m$diagnostics$lme4_overruled, gw17_gradient)
    expect_lt(m$diagnostics$stationary_distance, 0.1)
    expect_match(gw17_printed(m),
                 "lme4 warned that the model failed to converge; in centred and scaled units the fit is",
                 fixed = TRUE)
  }
})

test_that("a growth fit lme4 accepts is returned as before, without a refit", {
  skip_on_cran()
  h <- gw17_fit_with(gw17_growth, first = NULL)
  m <- h$value
  expect_identical(h$refits, 0L)
  expect_identical(h$seen, character(0))
  expect_null(stats::getCall(m$model)$start)
  expect_null(stats::getCall(m$model)$control)
  expect_null(m$diagnostics$lme4_overruled)
  # It is lmer's fit of the model, to the last digit.
  direct <- gw17_quiet(lme4::lmer(wellbeing ~ wave + (wave | id) + (1 | stratum), data = m$data))
  expect_identical(lme4::getME(m$model, "theta"), lme4::getME(direct, "theta"))
  expect_identical(lme4::REMLcrit(m$model), lme4::REMLcrit(direct))
})

test_that("a fit from starting values of the user's own is not refitted", {
  skip_on_cran()
  # No refit is asked for, however lme4 reads the fit.
  real <- maihda_lme4_accepted_refit
  asked <- 0L
  local_mocked_bindings(maihda_lme4_accepted_refit = function(fitted, fit_call, fit_env, tol = 1e-4) {
    asked <<- asked + 1L
    real(fitted, fit_call, fit_env, tol)
  })
  gw17_fit_with(gw17_growth, refit = "accepted")
  expect_identical(asked, 1L)
  theta <- c(0.9, 0.05, 0.35, 1)
  h <- gw17_fit_with(function() gw17_growth(start = theta), refit = "accepted")
  m <- h$value
  expect_identical(asked, 1L)
  expect_identical(h$refits, 0L)
  expect_true(is.name(stats::getCall(m$model)$start))
  expect_identical(unlist(m$model@optinfo$conv$lme4$messages), gw17_gradient)
  # The record still decides on the warning.
  expect_identical(h$seen, character(0))
  expect_identical(m$diagnostics$lme4_overruled, gw17_gradient)
  # A shorter spelling of the argument is the user's start as well.
  h <- gw17_fit_with(function() gw17_growth(sta = theta), refit = "accepted")
  expect_identical(asked, 1L)
  expect_identical(h$refits, 0L)
  # So it is for a glmer fit, whose start is a list.
  gw17_fit_with(gw17_logistic, refit = "accepted")
  expect_identical(asked, 2L)
  h <- gw17_fit_with(function() gw17_logistic(start = list(theta = 0.5)), refit = "accepted")
  expect_identical(asked, 2L)
  expect_identical(h$refits, 0L)
})

test_that("a fit the record does not call converged keeps lme4's warning", {
  skip_on_cran()
  # A fit stopped early by a step tolerance of 0.3, far from its optimum, with no refit on
  # offer: lme4's checks flag it, the optimiser reports no failure, and the record does not
  # overrule lme4.
  local_mocked_bindings(maihda_lme4_flagged_refit = function(model, fit_call, fit_env, tol = 1e-4) NULL)
  h <- gw17_heard(suppressMessages(
    gw17_growth(control = lme4::lmerControl(optCtrl = list(xtol_rel = 0.3)))))
  m <- h$value
  said <- unlist(m$model@optinfo$conv$lme4$messages)
  expect_gt(length(said), 0L)
  expect_true(m$model@optinfo$conv$opt == 0)
  expect_false(m$diagnostics$converged)
  expect_null(m$diagnostics$lme4_overruled)
  # Every warning lme4 raised for it is passed on, as lme4 worded it.
  expect_identical(h$seen, paste("warning:", said))
  expect_match(gw17_printed(m), "Convergence warnings reported by lme4", fixed = TRUE)
})

test_that("warnings of an lme4 fit that are not its convergence checks' are passed on", {
  skip_on_cran()
  real_hold <- maihda_hold_fit_conditions
  h <- gw17_heard(with_mocked_bindings(
    suppressMessages(gw17_growth()),
    maihda_hold_fit_conditions = function(expr) {
      pair <- gw17_as_flagged(gw17_as_accepted(real_hold(expr)))
      pair$conditions <- c(list(simpleWarning("a warning of another kind")), pair$conditions)
      pair
    },
    maihda_lme4_flagged_refit = function(model, fit_call, fit_env, tol = 1e-4) NULL))
  expect_identical(h$seen, "warning: a warning of another kind")
  expect_identical(h$value$diagnostics$lme4_overruled, gw17_gradient)
})

# ---- fit_maihda(): outside growth models --------------------------------------------------------

test_that("fit_maihda() returns the refit of a flagged glmer fit when lme4 accepts it", {
  skip_on_cran()
  h <- gw17_fit_with(gw17_logistic, refit = "accepted")
  m <- h$value
  expect_identical(h$refits, 1L)
  expect_identical(h$seen, character(0))
  expect_false(maihda_lme4_check_flagged(m$model))
  expect_true(m$diagnostics$converged)
  expect_null(m$diagnostics$lme4_overruled)
  # The refit: glmer's second stage from the first fit's estimates, its control written into
  # the call as lme4 writes one.
  start <- stats::getCall(m$model)$start
  expect_identical(start$theta, unname(lme4::getME(h$first, "theta")))
  expect_identical(start$fixef, lme4::fixef(h$first))
  expect_identical(stats::getCall(m$model)$control, quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  fit_env <- environment(stats::formula(m$model))
  expect_false(get(".maihda_arg_finish_control", envir = fit_env)$nAGQ0initStep)
  expect_lte(maihda_glmer_rise(m$model, h$first), 1e-4)
  expect_lt(abs(as.numeric(stats::logLik(m$model)) - as.numeric(stats::logLik(h$first))), 1e-3)
})

test_that("a flagged glmer fit whose refit lme4 flags too stays, without lme4's warnings", {
  skip_on_cran()
  # The gradient check's warning and the Hessian check's, which joins its two messages.
  all_three <- c(gw17_gradient, gw17_eigenvalue, gw17_ratio)
  h <- gw17_fit_with(gw17_logistic, first = all_three, refit = "flagged")
  m <- h$value
  expect_identical(h$refits, 1L)
  expect_null(stats::getCall(m$model)$start)
  expect_identical(h$seen, character(0))
  expect_true(m$diagnostics$converged)
  expect_identical(m$diagnostics$lme4_overruled, all_three)
  expect_lt(m$diagnostics$stationary_distance, 0.1)
  expect_match(gw17_printed(m), "lme4 warned that the model failed to converge", fixed = TRUE)
  # The refit that was not returned leaves nothing bound in the fit's environment.
  expect_false(exists(".maihda_arg_finish_control", envir = environment(stats::formula(m$model)),
                      inherits = FALSE))
  # The Hessian check's alone, as with a covariate in large units.
  h <- gw17_fit_with(gw17_logistic, first = c(gw17_eigenvalue, gw17_ratio), refit = "flagged")
  expect_identical(h$seen, character(0))
  expect_identical(h$value$diagnostics$lme4_overruled, c(gw17_eigenvalue, gw17_ratio))
  expect_match(gw17_printed(h$value), "lme4 warned that the model is nearly unidentifiable",
               fixed = TRUE)
})

test_that("the record is still told which fits their fitter has searched", {
  skip_on_cran()
  real <- maihda_fit_diagnostics
  told <- NULL
  local_mocked_bindings(maihda_fit_diagnostics = function(model, longitudinal_info = NULL,
                                                          searched = FALSE) {
    told <<- searched
    real(model, longitudinal_info, searched = searched)
  })
  d <- gw17_binary()
  d$stratum <- interaction(d$gender, d$ses, drop = TRUE)
  d$y <- d$z + rep(c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2), each = 30L) +
    rep(seq(-1.45, 1.45, length.out = 30L), times = 6L)
  # An lmer fit and a glmer fit with a random slope are made from their fitters' searches ...
  gw17_quiet(fit_maihda(y ~ z + (1 + z | stratum), data = d))
  expect_true(told)
  gw17_quiet(fit_maihda(yb ~ z + (1 + z | stratum), data = d, family = "binomial"))
  expect_true(told)
  # ... and fits without one, or of a growth model, are not.
  gw17_quiet(fit_maihda(y ~ z + (1 | stratum), data = d))
  expect_false(told)
  gw17_quiet(fit_maihda(yb ~ z + (1 | stratum), data = d, family = "binomial"))
  expect_false(told)
  gw17_quiet(gw17_growth())
  expect_false(told)
})

test_that("an intercept-only lmer fit lme4 accepts raises and returns what it did before", {
  skip_on_cran()
  d <- gw17_binary()
  d$y <- d$z + rep(c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2), each = 30L) +
    rep(seq(-1.45, 1.45, length.out = 30L), times = 6L)
  h <- gw17_heard(fit_maihda(y ~ z + (1 | gender:ses), data = d))
  m <- h$value
  expect_identical(h$seen, character(0))
  expect_null(stats::getCall(m$model)$start)
  expect_null(stats::getCall(m$model)$control)
  direct <- gw17_quiet(lme4::lmer(y ~ z + (1 | stratum), data = m$data))
  expect_identical(lme4::getME(m$model, "theta"), lme4::getME(direct, "theta"))
  # lme4's note of a singular fit is a message, and is passed on as one: ten strata of one
  # mean.
  flat <- data.frame(y = rep(c(1, 2, 3, 4), 10), g = factor(rep(1:10, each = 4)))
  h <- gw17_heard(suppressWarnings(fit_maihda(y ~ 1 + (1 | g), data = flat)))
  expect_true(any(grepl("message: boundary (singular) fit", h$seen, fixed = TRUE)))
})
