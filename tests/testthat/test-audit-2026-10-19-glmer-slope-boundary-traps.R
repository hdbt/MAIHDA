# Audit pass 2026-10-19: glmer fits with a random slope, outside growth models, that stop at
# a bound below the maximum of their likelihood.
#
# A random-effect block with a slope, such as (1 + x | stratum), can stop with its intercept
# Cholesky entry at its bound of 0 beside a non-zero intercept-slope entry: glmer's first
# stage takes it there and its second stage cannot cross the bound, though the covariance
# there is the same with that entry's sign reversed and the likelihood is higher on that
# other side. It happens where the intercept variance is small next to the slope variance,
# as in an adjusted model on additive strata or strata that fan out from x = 0. The refits
# and the convergence record of the glmer growth fits did not cover these fits. Now a glmer
# or glmer.nb fit that is not a growth model is refitted once from the other side of such a
# bound, and the record flags any glmer fit still stopped at one, naming the block's terms.
#
# Random numbers are drawn only inside with_seed(), which puts the suite's stream back.

with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  force(expr)
}

quiet <- function(e) suppressWarnings(suppressMessages(e))

heard <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(expr,
    message = function(m) { seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { seen <<- c(seen, paste("warning:", conditionMessage(w))); invokeRestart("muffleWarning") })
  list(value = value, seen = seen)
}

# Eight additive strata (A x B x C) with a stratum slope on x and no stratum effect at x = 0,
# so that the adjusted model's intercept variance is small next to its slope variance.
slope_design <- function(seed, n = 600) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1)
    st <- sample(8, n, TRUE)
    x <- round(stats::rnorm(n), 2)
    u1 <- stats::rnorm(8, 0, 0.25)
    add <- 0.4 * (g$A - 0.5) - 0.3 * (g$B - 0.5) + 0.2 * (g$C - 0.5)
    eta <- 0.5 + 0.2 * x + add[st] + u1[st] * x
    data.frame(y = stats::rpois(n, exp(eta)), x = x, A = factor(g$A[st]), B = factor(g$B[st]),
               C = factor(g$C[st]), stratum = factor(st))
  })
}

slope_formula <- y ~ x + A + B + C + (1 + x | stratum)

# Seed 7: lme4's own fit stops with the stratum intercept entry at 0 (theta 0, 0.2119,
# 0.0101) at a log-likelihood of -982.13197; the maximum, from the best of lme4's fits and a
# search over its deviance function with no bound, is -981.32970.
trap_7 <- c(0, 0.2119, 0.0101)

# ---- which fits ---------------------------------------------------------------------------

test_that("glmer and glmer.nb fits outside growth models are refitted, unless the user gives a start", {
  expect_true(maihda_glmer_bound_applies(quote(lme4::glmer), NULL, list()))
  expect_true(maihda_glmer_bound_applies(quote(lme4::glmer.nb), NULL, list()))
  expect_true(maihda_glmer_bound_applies(quote(lme4::glmer), NULL, list(control = 1, nAGQ = 2)))
  expect_false(maihda_glmer_bound_applies(quote(lme4::glmer), NULL, list(start = 1)))
  expect_false(maihda_glmer_bound_applies(quote(lme4::lmer), NULL, list()))
  expect_false(maihda_glmer_bound_applies(quote(brms::brm), NULL, list()))
  # A growth model keeps its own refits.
  lng <- list(id = "id", time = "wave", time_degree = 1L, stratum_slope = TRUE)
  expect_false(maihda_glmer_bound_applies(quote(lme4::glmer), lng, list()))
  expect_true(maihda_glmer_growth_applies(quote(lme4::glmer), lng, list()))
})

test_that("fit_maihda() sends every glmer fit outside growth models through the refit", {
  skip_on_cran()
  calls <- 0L
  real <- maihda_fit_glmer_bound
  local_mocked_bindings(maihda_fit_glmer_bound = function(...) { calls <<- calls + 1L; real(...) })
  d <- slope_design(7)
  quiet(fit_maihda(y ~ x + (1 | stratum), data = d, family = "poisson"))
  quiet(fit_maihda(slope_formula, data = d, family = "poisson"))
  quiet(fit_maihda(y ~ x + (1 | stratum), data = d, family = "negbinomial"))
  expect_identical(calls, 3L)
  # Not lmer, and not a start of the user's own.
  quiet(fit_maihda(log1p(y) ~ x + (1 + x | stratum), data = d))
  quiet(fit_maihda(slope_formula, data = d, family = "poisson", start = list(theta = trap_7)))
  expect_identical(calls, 3L)
})

# ---- the refit, with stand-in fits --------------------------------------------------------

test_that("a fit is refitted once across a bound, and the refit kept only when it is returned", {
  across <- NULL
  local_mocked_bindings(
    maihda_glmer_mirror_refit = function(kept, fit_call, fit_env, tol = 1e-4, optimizer = NULL) {
      across <<- list(kept$value$crit, tol, optimizer)
      if (kept$value$crit > 9) list(value = list(crit = 6), conditions = list(simpleMessage("across\n"))) else NULL
    }
  )
  env <- new.env()
  env$fake_fit <- function(crit) {
    message("first fit says")
    warning("first fit warns")
    list(crit = crit)
  }
  # The refit is made from the first fit, by bobyqa with the tolerance given, and speaks
  # alone.
  h <- heard(maihda_fit_glmer_bound(quote(fake_fit(10)), env, tol = 0.5))
  expect_identical(across, list(10, 0.5, "bobyqa"))
  expect_identical(h$value$crit, 6)
  expect_identical(h$seen, "across\n")
  expect_identical(formals(maihda_fit_glmer_bound)$tol, 1e-4)
  # No refit returned: the first fit stands, with its own message and warning.
  h <- heard(maihda_fit_glmer_bound(quote(fake_fit(8)), env))
  expect_identical(h$value$crit, 8)
  expect_identical(h$seen, c("first fit says\n", "warning: first fit warns"))
  # A refit that fails is passed over; a first fit that fails fails.
  local_mocked_bindings(maihda_glmer_mirror_refit = function(...) stop("the refit fails"))
  h <- heard(maihda_fit_glmer_bound(quote(fake_fit(10)), env))
  expect_identical(h$value$crit, 10)
  expect_identical(h$seen, c("first fit says\n", "warning: first fit warns"))
  env$fake_fit <- function(crit) stop("the first fit fails")
  expect_error(maihda_fit_glmer_bound(quote(fake_fit(10)), env), "the first fit fails")
})

# ---- the record's sentence ----------------------------------------------------------------

test_that("the record's sentence names the terms of a block outside growth models", {
  tail <- ": the fit stopped at a bound of lme4's parametrisation, not at the maximum."
  expect_identical(
    maihda_lme4_trap_message("stratum", "(Intercept)", 2L, c("(Intercept)", "x")),
    paste0("The intercept variance of the 'stratum' block is at zero, where the likelihood ",
           "still rises with the sign of its covariance with x reversed", tail))
  expect_match(maihda_lme4_trap_message("stratum", "(Intercept)", 3L, c("(Intercept)", "x", "z")),
               "rises with the signs of its covariances with x and z reversed:", fixed = TRUE)
  # A later column at its bound: what is left of its variance after the terms before it.
  expect_identical(
    maihda_lme4_trap_message("stratum", "x", 3L, c("(Intercept)", "x", "z")),
    paste0("In the 'stratum' block the variance of x left after the intercept is zero, where ",
           "the likelihood still rises with the sign of its covariance with z reversed", tail))
  expect_match(maihda_lme4_trap_message("g", "z", 4L, c("(Intercept)", "x", "z", "w")),
               "the variance of z left after the intercept and x is zero", fixed = TRUE)
  # A block without an intercept.
  expect_identical(
    maihda_lme4_trap_message("site", "x", 2L, c("x", "z")),
    paste0("The variance of x in the 'site' block is at zero, where the likelihood still rises ",
           "with the sign of its covariance with z reversed", tail))
  # Without the terms, or with a term not among them, a growth block's sentence.
  growth <- maihda_lme4_trap_message("id", "(Intercept)", 2L)
  expect_match(growth, "covariance with the slope reversed", fixed = TRUE)
  expect_identical(maihda_lme4_trap_message("id", "(Intercept)", 2L, c("x", "z")), growth)
  expect_match(maihda_lme4_trap_message("id", ".maihda_ctime", 3L),
               "the slope is perfectly correlated with the intercept", fixed = TRUE)
})

test_that("a block's terms are found by the position of its entries in theta", {
  # lme4 can give two blocks the same name.
  cnms <- list(stratum = "z", stratum = c("(Intercept)", "x"), site = c("(Intercept)", "w", "v"))
  expect_identical(maihda_lme4_block_terms(cnms, 1L), "z")
  expect_identical(maihda_lme4_block_terms(cnms, 2L), c("(Intercept)", "x"))
  expect_identical(maihda_lme4_block_terms(cnms, 4L), c("(Intercept)", "x"))
  expect_identical(maihda_lme4_block_terms(cnms, 5L), c("(Intercept)", "w", "v"))
  expect_identical(maihda_lme4_block_terms(cnms, 10L), c("(Intercept)", "w", "v"))
})

# ---- live fits ----------------------------------------------------------------------------

test_that("a cross-sectional glmer fit stopped at a bound is refitted to the maximum", {
  skip_on_cran()
  for (s in c(7, 1)) {
    # Seed 1: lme4's own fit stops the same way, at -969.99706; the maximum is -969.35826.
    m <- quiet(fit_maihda(slope_formula, data = slope_design(s), family = "poisson"))
    target <- if (s == 7) -981.32970 else -969.35826
    expect_lt(abs(as.numeric(stats::logLik(m$model)) - target), 0.005, label = paste("seed", s))
    expect_gt(lme4::getME(m$model, "theta")[[1]], 0.01)
    expect_true(m$diagnostics$converged)
    expect_null(m$diagnostics$boundary_traps)
  }
})

test_that("the refit starts across the bound and recovers a fit stopped there", {
  skip_on_cran()
  env <- new.env()
  env$data <- slope_design(7)
  call <- quote(lme4::glmer(y ~ x + A + B + C + (1 + x | stratum), data = data,
                            family = stats::poisson))
  trapped_call <- call
  trapped_call$start <- list(theta = trap_7)
  seen <- list()
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen[[length(seen) + 1L]] <<- call
    real(call, env)
  })
  m <- maihda_fit_glmer_bound(trapped_call, env)
  # The first fit stays at the bound, 0.80 below the maximum ...
  expect_length(seen, 2L)
  first <- quiet(real(trapped_call, env))$value
  expect_lt(lme4::getME(first, "theta")[[1]], 1e-4)
  expect_lt(as.numeric(stats::logLik(first)), -982.1)
  # ... and the refit, from its estimates with the intercept entry across the bound by
  # glmer's second stage alone, by bobyqa, reaches it.
  th <- unname(lme4::getME(first, "theta"))
  expect_equal(seen[[2]]$start$theta, c(0.05, -th[2], th[3]))
  expect_identical(seen[[2]]$start$fixef, lme4::fixef(first))
  expect_false(eval(seen[[2]]$control, env)$nAGQ0initStep)
  expect_identical(eval(seen[[2]]$control, env)$optimizer, c("bobyqa", "bobyqa"))
  expect_identical(eval(seen[[2]]$control, env)$optCtrl, list(maxfun = 1e5))
  expect_lt(abs(as.numeric(stats::logLik(m)) - (-981.32970)), 0.005)
  # A fit not at a bound is not refitted.
  seen <- list()
  inside <- new.env()
  inside$data <- slope_design(3)
  m2 <- maihda_fit_glmer_bound(call, inside)
  expect_length(seen, 1L)
  expect_gt(lme4::getME(m2, "theta")[[1]], 1e-4)
})

test_that("the refit is made by bobyqa with the settings it takes, or else by the control's own", {
  skip_on_cran()
  env <- new.env()
  env$data <- slope_design(7)
  trapped_call <- quote(lme4::glmer(y ~ x + A + B + C + (1 + x | stratum), data = data,
                                    family = stats::poisson))
  trapped_call$start <- list(theta = trap_7)
  kept <- quiet(maihda_eval_holding_conditions(trapped_call, env))
  expect_lt(lme4::getME(kept$value, "theta")[[1]], 1e-4)
  seen <- NULL
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen <<- call
    real(call, env)
  })
  refit_with <- function(ctrl, optimizer = "bobyqa") {
    env$ctrl <- ctrl
    call <- quote(lme4::glmer(y ~ x + A + B + C + (1 + x | stratum), data = data,
                              family = stats::poisson, control = ctrl))
    out <- maihda_glmer_mirror_refit(kept, call, env, optimizer = optimizer)
    list(out = out, control = eval(seen$control, env))
  }
  # A Nelder-Mead control of the user's own: its first stage and its evaluation limit kept,
  # a setting bobyqa does not take dropped, so bobyqa does not warn about it.
  r <- refit_with(lme4::glmerControl(optimizer = "Nelder_Mead",
                                     optCtrl = list(FtolAbs = 1e-8, maxfun = 3e4)))
  expect_identical(r$control$optimizer, c("Nelder_Mead", "bobyqa"))
  expect_identical(r$control$optCtrl, list(maxfun = 3e4))
  expect_false(r$control$nAGQ0initStep)
  expect_lt(abs(as.numeric(stats::logLik(r$out$value)) - (-981.32970)), 0.005)
  said <- vapply(r$out$conditions, conditionMessage, character(1))
  expect_false(any(grepl("unused control arguments", said, fixed = TRUE)))
  # Without a limit set, bobyqa's is 1e5 evaluations; a bobyqa control keeps its settings.
  r <- refit_with(lme4::glmerControl(optimizer = "Nelder_Mead", optCtrl = list(FtolAbs = 1e-8)))
  expect_identical(r$control$optCtrl, list(maxfun = 1e5))
  r <- refit_with(lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-8)))
  expect_identical(r$control$optCtrl, list(rhoend = 1e-8, maxfun = 1e5))
  # Without an optimiser named, the refit keeps the control's own.
  r <- refit_with(lme4::glmerControl(optimizer = "Nelder_Mead", optCtrl = list(FtolAbs = 1e-8)),
                  optimizer = NULL)
  expect_identical(r$control$optimizer, c("Nelder_Mead", "Nelder_Mead"))
  expect_identical(r$control$optCtrl, list(FtolAbs = 1e-8))
})

test_that("a theta that is not a full Cholesky factor per block has no columns at a bound", {
  two <- list(stratum = c("(Intercept)", "x"))
  expect_length(maihda_lme4_trap_columns(c(0, 0.2, 0.1), two), 1L)
  # lme4's diag(1 + x | stratum): two entries for a two-term block.
  expect_length(maihda_lme4_trap_columns(c(0, 0.2), two), 0L)
  both <- list(stratum = c("(Intercept)", "x"), A = c("(Intercept)", "z"))
  expect_length(maihda_lme4_trap_columns(c(0, 0.2, 0.1, 0.3, 0, 0.1), both), 1L)
  expect_length(maihda_lme4_trap_columns(c(0, 0.2, 0.1, 0.3, 0.1), both), 0L)
})

test_that("the record flags a cross-sectional fit stopped at a bound where the likelihood still rises", {
  skip_on_cran()
  d <- slope_design(7)
  # The user's own start is fitted as given: lme4's trap.
  m <- quiet(fit_maihda(slope_formula, data = d, family = "poisson", start = list(theta = trap_7)))
  expect_lt(as.numeric(stats::logLik(m$model)), -982.1)
  dg <- m$diagnostics
  expect_false(dg$converged)
  expect_identical(dg$boundary_traps$block, "stratum")
  expect_identical(dg$boundary_traps$term, "(Intercept)")
  expect_gt(dg$boundary_traps$drop, 0.01)
  sentence <- paste0("The intercept variance of the 'stratum' block is at zero, where the ",
                     "likelihood still rises with the sign of its covariance with x reversed: ",
                     "the fit stopped at a bound of lme4's parametrisation, not at the maximum.")
  expect_true(sentence %in% dg$messages)
  out <- maihda_format_fit_diagnostics(dg)
  expect_true("Convergence problems of the lme4 fit:" %in% out)
  expect_true(paste0("  - ", sentence) %in% out)
  # print() shows it.
  printed <- utils::capture.output(print(m))
  expect_true(any(grepl("Convergence problems of the lme4 fit:", printed, fixed = TRUE)))
  # The same fit read as a growth fit is worded as one.
  growth <- maihda_fit_diagnostics(m$model, list(id = "id", time = "wave"))
  expect_true(any(grepl("covariance with the slope reversed", growth$messages, fixed = TRUE)))
  # At the maximum there is nothing to flag.
  ok <- quiet(fit_maihda(slope_formula, data = d, family = "poisson"))
  expect_identical(nrow(maihda_lme4_boundary_traps(ok$model)), 0L)
})

test_that("an intercept-only glmer fit is made as before", {
  skip_on_cran()
  d <- slope_design(7)
  seen <- 0L
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen <<- seen + 1L
    real(call, env)
  })
  m <- quiet(fit_maihda(y ~ x + A + B + C + (1 | stratum), data = d, family = "poisson"))
  expect_identical(seen, 1L)
  direct <- quiet(lme4::glmer(y ~ x + A + B + C + (1 | stratum), data = d, family = stats::poisson))
  expect_identical(lme4::getME(m$model, "theta"), lme4::getME(direct, "theta"))
  expect_identical(lme4::fixef(m$model), lme4::fixef(direct))
  expect_identical(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(direct)))
  # lme4's own message on a fit still reaches the caller: strata of A alone, with A among
  # the fixed effects, have no variance left.
  h <- heard(fit_maihda(y ~ A + (1 | A), data = d[names(d) != "stratum"], family = "poisson"))
  expect_true(any(grepl("boundary (singular) fit", h$seen, fixed = TRUE)))
})

test_that("the refit keeps an offset given as an argument and the record leaves it alone", {
  skip_on_cran()
  d <- slope_design(7)
  d$expo <- with_seed(70, stats::runif(nrow(d), 0.5, 2))
  d$off <- log(d$expo)
  off <- d$off
  m <- quiet(fit_maihda(slope_formula, data = d, family = "poisson", offset = log(expo)))
  # The kept fit is the one lme4 makes with the same offset from its estimates.
  direct <- quiet(lme4::glmer(slope_formula, data = d, family = stats::poisson, offset = off,
                              start = list(theta = unname(lme4::getME(m$model, "theta")),
                                           fixef = lme4::fixef(m$model)),
                              control = lme4::glmerControl(nAGQ0initStep = FALSE)))
  expect_lt(abs(as.numeric(stats::logLik(m$model)) - as.numeric(stats::logLik(direct))), 1e-6)
  expect_identical(unname(m$model@frame[["(offset)"]]), off)
  # A fit at a bound, given as the user's start, is measured without moving its offset.
  t0 <- quiet(fit_maihda(slope_formula, data = d, family = "poisson", offset = log(expo),
                         start = list(theta = c(0, unname(lme4::getME(m$model, "theta"))[2:3]))))
  expect_identical(unname(t0$model@frame[["(offset)"]]), off)
  expect_identical(unname(lme4::getME(t0$model, "offset")), off)
})
