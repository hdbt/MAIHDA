# Audit 2026-10-14 (packet 2026-10-14-ordinal-covariate-scale): with a covariate far
# from zero or on a large or tiny scale, fit_maihda(engine = "ordinal") returned
# ordinal::clmm()'s starting values or a fit short of the maximum, and its standard
# errors were wrong even at the maximum -- CONFIRMED, wider than reported.
#
# On 2 x 3 strata of 30 rows with one covariate in exact affine spellings (z, age in
# years, a birth year, income in thousands and in currency units), nlminb() left a birth
# year and currency units at clmm()'s starting values with "false convergence (8)": SD 1,
# VPC 0.2331, and maihda() printed a PCV of 0.0000, the ratio of two starting SDs. In a
# survey of 96 simulated data sets (6 to 48 strata, 48 to 1,200 rows, true SD 0 to 1; 960
# fits) 335 of 339 "false convergence" fits were the starting values, and 12 fits that
# nlminb() reported as converged (11 with a birth year, 1 with income in thousands) were
# 5e-4 to 0.65 log-likelihood units short and recorded as converged; age in years failed
# the same way on 12,000 and 120,000 rows, and age + I(age^2) on about 1,000. A covariate
# with an SD of 1e-4 to 1e-7 left 11 of 24 fits short and recorded as converged.
# ordinal's "Rescale variables?"
# warning comes from its fixed-effects start fit (clm()'s largest Hessian eigenvalue above
# 1 / sqrt(.Machine$double.eps)), not from the clmm fit, and was absent on 12,000 rows.
#
# Not in the finding, same cause: clmm() takes the gradient and Hessian by central
# differences with a fixed step of 1e-4 in every parameter, too coarse for the
# coefficient of such a covariate, so even at the maximum the standard errors were wrong:
# a birth-year coefficient's SE 0.013 to 1.77 times the right one over 700 interior fits
# of that survey, and 7.5 times too small at a variance of zero on this file's fixture,
# against clm()'s analytic SEs.
#
# Per the user, the package now fits the model on its fixed-effects design with every
# column centred and scaled (only scaled under threshold = "symmetric2"), the boundary
# fit included, and returns the clmm() fit of the formula built at the carried-back
# maximum without a further search, with the gradient and Hessian of the rescaled fit
# carried over exactly, so that its standard errors are that fit's. A caller's start is
# in the formula's own terms, and a fit given one is made as the formula stands.
#
# No block draws a random number.

cs14_observe <- function(expr) {
  warns <- character(0)
  notes <- character(0)
  res <- withCallingHandlers(
    expr,
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      notes <<- c(notes, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(res = res, warns = warns, notes = notes)
}

# gender x ses strata of 30 rows: a stratum's outcomes are its shift plus the same
# fixed spread, cut into four classes, with the same thirty ages in every stratum, and
# the age in five exact affine spellings: x = from + per * z.
cs14_years <- c(53, 49, 55, 51, 55, 43, 46, 55, 71, 56, 19, 53, 23, 41, 63, 45, 53, 54,
                54, 55, 38, 61, 23, 36, 55, 47, 51, 53, 41, 66)
cs14_spellings <- list(z = c(0, 1), age = c(45, 12), birth = c(1975, -12),
                       income_k = c(40, 15), income = c(40000, 15000))
cs14_strata <- function(shift = rep(0, 6), cuts = c(-0.8, 0, 0.8)) {
  g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 30), , drop = FALSE]
  d$y <- rep(shift, each = 30) + rep(seq(-1.45, 1.45, length.out = 30), times = nrow(g))
  d$yo <- factor(findInterval(d$y, cuts) + 1L, levels = seq_len(length(cuts) + 1L),
                 ordered = TRUE)
  d$age <- rep(cs14_years, times = nrow(g))
  d$z <- (d$age - 45) / 12
  for (s in names(cs14_spellings)) {
    d[[s]] <- cs14_spellings[[s]][1] + cs14_spellings[[s]][2] * d$z
  }
  rownames(d) <- NULL
  d
}
cs14_shift <- c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2)

cs14_fit <- function(d, formula, ...) {
  suppressMessages(suppressWarnings(
    fit_maihda(formula, data = d, engine = "ordinal", ...)))
}
cs14_formula <- function(x) stats::as.formula(paste("yo ~", x, "+ (1 | gender:ses)"))
cs14_variance <- function(m) maihda_clmm_variances(m)$stratum
cs14_loglik <- function(m) as.numeric(stats::logLik(m$model))
cs14_se <- function(m) {
  nm <- c(names(m$model$alpha), names(m$model$beta))
  sqrt(diag(stats::vcov(m$model)))[nm]
}
# The cumulative predictor theta_j - x' beta at the rows of `x` (one column per
# coefficient), one column per cut point.
cs14_predictor <- function(theta, beta, x) {
  outer(-drop(x %*% beta), theta, "+")
}

test_that("the map from rescaled to the formula's own parameters", {
  # By hand, flexible thresholds: three cut points, two columns centred on 45 and 0.5
  # and scaled by 12 and 0.5.
  map <- maihda_clmm_rescaled_map(diag(3), centre = c(45, 0.5), spread = c(12, 0.5))
  expect_equal(map, rbind(cbind(diag(3), matrix(c(45 / 12, 0.5 / 0.5), 3, 2, byrow = TRUE)),
                          cbind(matrix(0, 2, 3), diag(c(1 / 12, 2)))))
  # The defining property: the rescaled parameters on the rescaled columns and the
  # mapped ones on the columns as they are give the same predictor at every row, under
  # each threshold structure that can take up a shift.
  x <- cbind(c(19, 45, 71, 33, 60), c(0, 1, 1, 0, 1))
  centre <- c(45, 0.6)
  spread <- c(12, 0.49)
  xs <- sweep(sweep(x, 2, centre), 2, spread, "/")
  structures <- list(flexible = diag(3),
                     equidistant = cbind(1, 0:2),
                     symmetric = rbind(c(1, -1), c(1, 0), c(1, 1)))
  for (tJac in structures) {
    alpha_s <- if (ncol(tJac) == 3L) c(-1.2, 0.1, 1.3) else c(0.2, 0.9)
    beta_s <- c(0.7, -0.4)
    par <- drop(maihda_clmm_rescaled_map(tJac, centre, spread) %*% c(alpha_s, beta_s))
    alpha <- par[seq_len(ncol(tJac))]
    beta <- par[ncol(tJac) + 1:2]
    expect_equal(cs14_predictor(drop(tJac %*% alpha), beta, x),
                 cs14_predictor(drop(tJac %*% alpha_s), beta_s, xs))
    expect_equal(beta, beta_s / spread)
  }
  # Thresholds that cannot take up a shift (threshold = "symmetric2"): scaled only.
  symmetric2 <- matrix(c(-1, 0, 1), 3, 1)
  par <- drop(maihda_clmm_rescaled_map(symmetric2, c(0, 0), spread) %*% c(0.8, 0.7, -0.4))
  expect_equal(par, c(0.8, 0.7 / 12, -0.4 / 0.49))
  expect_error(maihda_clmm_rescaled_map(symmetric2, centre, spread),
               "cannot take up a shift")
  # No column: nothing to map.
  expect_identical(maihda_clmm_rescaled_map(diag(3), numeric(0), numeric(0)), diag(3))
})

test_that("the call that rebuilds a fit at a point without a search", {
  args <- list(formula = yo ~ x + (1 | stratum), data = quote(data), link = "logit",
               Hess = TRUE)
  dots <- list(start = c(1, 2), meth = "ucminf", useMatrix = FALSE, Hess = TRUE,
               iter.max = 50L, nAGQ = 5L, threshold = "equidistant",
               control = list(maxIter = 100))
  plain <- maihda_clmm_rebuild_call(args, dots, general = FALSE)
  expect_identical(plain[[1L]], quote(ordinal::clmm))
  expect_identical(plain$start, quote(.maihda_clmm_start))
  expect_identical(plain$method, "nlminb")
  expect_identical(plain$iter.max, 0L)
  expect_false(plain$Hess)
  expect_identical(sum(names(plain) == "Hess"), 1L)
  expect_identical(sum(names(plain) == "iter.max"), 1L)
  expect_false(any(c("meth", "useMatrix") %in% names(plain)))
  expect_identical(plain$nAGQ, 5L)
  expect_identical(plain$threshold, "equidistant")
  expect_identical(plain$control, list(maxIter = 100))
  expect_identical(plain$formula, args$formula)
  general <- maihda_clmm_rebuild_call(args, dots, general = TRUE)
  expect_true(general$useMatrix)
  expect_null(general$nAGQ)
  expect_identical(general$start, quote(.maihda_clmm_start))
  # Unnamed dots are passed on as they are.
  expect_identical(maihda_clmm_rebuild_call(args, list(), general = FALSE)$link, "logit")
})

test_that("the design a model is fitted on", {
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  d$stratum <- factor(paste(d$gender, d$ses))
  d$gender <- factor(d$gender)
  d$expo <- rep(c(1, 2, 4), length.out = nrow(d))
  env <- new.env()
  env$data <- d
  formula <- yo ~ birth + gender + offset(log(expo)) + (1 | stratum)
  environment(formula) <- env
  args <- list(formula = formula, data = quote(data), link = "logit", Hess = TRUE)
  design <- maihda_clmm_rescaled_design(args, list(contrasts = list(gender = "contr.sum"),
                                                   meth = "ucminf"), env)
  expect_identical(design$columns, c("birth", "gender1"))
  expect_identical(design$generated, c(".maihda_x1", ".maihda_x2"))
  sum_coded <- ifelse(d$gender == "F", 1, -1)
  expect_equal(design$centre, c(mean(d$birth), mean(sum_coded)))
  expect_equal(design$spread, c(stats::sd(d$birth), stats::sd(sum_coded)))
  fitted_on <- design$env$data
  expect_equal(fitted_on$.maihda_x1, (d$birth - mean(d$birth)) / stats::sd(d$birth))
  expect_equal(fitted_on$.maihda_offset, log(d$expo))
  expect_identical(
    deparse(design$args$formula),
    deparse(yo ~ .maihda_x1 + .maihda_x2 + offset(.maihda_offset) + (1 | stratum)))
  expect_identical(environment(design$args$formula), design$env)
  # The design has no factor left for `contrasts` to code; the rest is passed on.
  expect_null(design$dot_vals$contrasts)
  expect_identical(design$dot_vals$meth, "ucminf")
  expect_identical(design$call$formula, design$args$formula)
  # The caller's optimiser, however it is spelled, goes to the rescaled fit and does not
  # keep clmm() from handing out its design.
  for (dots in list(list(method = "ucminf"), list(meth = "ucminf"),
                    list(control = ordinal::clmm.control(method = "ucminf")))) {
    spelled <- maihda_clmm_rescaled_design(args, dots, env)
    expect_false(is.null(spelled))
    expect_identical(spelled$dot_vals, dots)
  }
  # The data the caller handed over is not touched.
  expect_false(any(startsWith(names(env$data), ".maihda")))

  # Thresholds that cannot take up a shift: the columns are scaled only.
  two <- maihda_clmm_rescaled_design(args, list(threshold = "symmetric2"), env)
  expect_identical(two$centre, c(0, 0))
  expect_equal(two$env$data$.maihda_x1, d$birth / stats::sd(d$birth))
  # A caller's start, a model without fixed effects, or a fixed part that does not vary:
  # nothing to rescale.
  expect_null(maihda_clmm_rescaled_design(args, list(start = 1:5), env))
  null_formula <- yo ~ 1 + (1 | stratum)
  environment(null_formula) <- env
  expect_null(maihda_clmm_rescaled_design(
    list(formula = null_formula, data = quote(data), link = "logit", Hess = TRUE), list(),
    env))
  # A constant column is left as it is, for clmm() to drop as it would.
  env$data$level <- 3
  level_formula <- yo ~ birth + level + (1 | stratum)
  environment(level_formula) <- env
  level <- maihda_clmm_rescaled_design(
    list(formula = level_formula, data = quote(data), link = "logit", Hess = TRUE), list(),
    env)
  expect_identical(level$centre[2], 0)
  expect_identical(level$spread[2], 1)
  expect_identical(level$env$data$.maihda_x2, rep(3, nrow(d)))
  # Data that already hold a name the rescaled design would use: not rescaled.
  env$data$.maihda_x1 <- 0
  expect_null(maihda_clmm_rescaled_design(args, list(), env))
})

test_that("a fit made on rescaled columns is asked by the derivatives taken there", {
  # The derivatives as clmm() took them on rescaled columns put the fit 0.25 standard
  # errors from a stationary point; carried over to a coefficient far below the step of
  # clmm()'s differences they are the same distance, but such derivatives taken there
  # could not be read.
  fine <- list(gradient = c(0.5, 0), Hessian = diag(c(4, 1)))
  carried <- list(gradient = c(0.5, 0), Hessian = diag(c(4, 1 / 7e-6^2)))
  expect_true(is.na(maihda_clmm_stationary_distance(carried)))
  attr(carried, "maihda_rescaled") <- fine
  expect_equal(maihda_clmm_stationary_distance(carried), 0.25)
  expect_equal(maihda_clmm_stationary_distance(fine), 0.25)
})

test_that("one covariate in five units gives one fit", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  fits <- lapply(names(cs14_spellings), function(x) cs14_fit(d, cs14_formula(x)))
  names(fits) <- names(cs14_spellings)
  top <- fits$z
  expect_gt(cs14_variance(top), 0.5)
  v_top <- stats::vcov(top$model)[c("1|2", "2|3", "3|4", "z"), c("1|2", "2|3", "3|4", "z")]
  # One maximum reached by searches on the five spellings' rescaled columns, which are
  # equal to rounding: the same log-likelihood, the parameters and their standard errors
  # to the optimiser's tolerance (3.4e-6 here).
  for (x in names(cs14_spellings)) {
    m <- fits[[x]]
    from <- cs14_spellings[[x]][1]
    per <- cs14_spellings[[x]][2]
    expect_equal(cs14_loglik(m), cs14_loglik(top), tolerance = 1e-8)
    expect_equal(cs14_variance(m), cs14_variance(top), tolerance = 1e-4)
    expect_true(m$diagnostics$converged)
    expect_false(m$diagnostics$singular)
    expect_length(m$diagnostics$messages, 0L)
    expect_false(is.null(attr(m$model, "maihda_rescaled")))
    # The coefficient and the cut points are the standardised fit's, in these units ...
    expect_equal(unname(m$model$beta), unname(top$model$beta) / per, tolerance = 1e-4)
    expect_equal(as.numeric(m$model$Theta),
                 as.numeric(top$model$Theta) + from * unname(top$model$beta) / per,
                 tolerance = 1e-4)
    # ... and so are their standard errors, by the exact linear map.
    jac <- rbind(cbind(diag(3), from / per), c(0, 0, 0, 1 / per))
    expect_equal(unname(cs14_se(m)), sqrt(diag(jac %*% v_top %*% t(jac))),
                 tolerance = 1e-4)
  }
  # The analysis: the standardised fit's VPC and PCV, and no convergence problem.
  analyse <- function(x) {
    suppressMessages(suppressWarnings(
      maihda(cs14_formula(x), data = d, engine = "ordinal")))
  }
  reference <- analyse("z")
  for (x in c("birth", "income")) {
    a <- analyse(x)
    expect_true(a$model$diagnostics$converged)
    expect_true(a$model_adjusted$diagnostics$converged)
    expect_equal(a$pcv$pcv, reference$pcv$pcv, tolerance = 1e-6)
    expect_gt(a$pcv$pcv, 0.1)
    printed <- utils::capture.output(print(a))
    expect_false(any(grepl("Convergence problems", printed, fixed = TRUE)))
  }
})

test_that("standard errors at a variance of zero are those of the model without it", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Identical strata: the between-stratum variance is zero, and the fit is the model
  # without the random effect, whose standard errors clm() gives from analytic
  # derivatives.
  flat <- cs14_strata()
  for (x in names(cs14_spellings)) {
    m <- cs14_fit(flat, cs14_formula(x))
    expect_lt(cs14_variance(m), 1e-8)
    expect_true(m$diagnostics$singular)
    expect_true(m$diagnostics$converged)
    # clm() warns that a birth year or currency units make its Hessian's eigenvalues
    # large; its derivatives are analytic all the same.
    without <- suppressWarnings(ordinal::clm(stats::as.formula(paste("yo ~", x)),
                                             data = flat))
    expect_equal(cs14_loglik(m), as.numeric(stats::logLik(without)), tolerance = 1e-8)
    se_clm <- sqrt(diag(stats::vcov(without)))
    expect_equal(unname(cs14_se(m)), unname(se_clm[names(cs14_se(m))]), tolerance = 1e-5)
  }
})

test_that("the record of a fit made on rescaled columns", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  m <- cs14_fit(d, cs14_formula("birth"))
  model <- m$model
  # The optimiser's return is the rescaled search's, its point carried back to the
  # formula's own terms; the stored derivatives are that search's, carried over.
  n_fixed <- length(model$coefficients)
  expect_equal(unname(model$optRes$par[seq_len(n_fixed)]), unname(model$coefficients))
  expect_equal(exp(model$optRes$par[n_fixed + 1L]), sqrt(cs14_variance(m)))
  expect_identical(names(model$gradient), c(names(model$coefficients), "ST1"))
  expect_identical(dimnames(model$Hessian), list(names(model$gradient), names(model$gradient)))
  taken <- attr(model, "maihda_rescaled")
  expect_identical(taken$columns, "birth")
  expect_equal(taken$centre, mean(d$birth))
  expect_equal(taken$spread, stats::sd(d$birth))
  carry <- diag(5)
  carry[1:3, 4] <- taken$centre / taken$spread
  carry[4, 4] <- 1 / taken$spread
  expect_equal(unname(model$Hessian),
               unname(t(solve(carry)) %*% taken$Hessian %*% solve(carry)))
  expect_equal(unname(model$gradient), unname(drop(t(solve(carry)) %*% taken$gradient)))
  # The search was nlminb()'s, as asked, from clmm()'s own start (the fixed-effects fit
  # and an SD of 1) on the rescaled columns, which is that fit in these units; its
  # counts and largest gradient are the record's.
  expect_identical(model$control$method, "nlminb")
  expect_false(model$control$useMatrix)
  without <- suppressWarnings(ordinal::clm(yo ~ birth, data = d))
  expect_equal(unname(model$start[seq_len(n_fixed)]), unname(stats::coef(without)),
               tolerance = 1e-4)
  expect_identical(model$start[n_fixed + 1L], 0)
  expect_gt(model$neval, 50)
  expect_identical(model$info$niter, paste0(model$neval, "(", model$Niter, ")"))
  expect_identical(model$info$max.grad,
                   formatC(max(abs(model$gradient)), digits = 2, format = "e"))
  # The call is the caller's, as clmm() records a call.
  expect_null(model$call$start)
  expect_null(model$call$iter.max)
  expect_null(model$call$useMatrix)
  expect_true(model$call$Hess)
  expect_identical(model$call$data, quote(data))
  expect_identical(names(model$call)[2:5], c("formula", "data", "Hess", "link"))
})

test_that("the call stored with a fit is the one the fit stands for", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  # The caller asked for clmm()'s general routine: the caller's call, which update()
  # can make again.
  general <- cs14_fit(d, cs14_formula("birth"),
                      control = ordinal::clmm.control(useMatrix = TRUE))
  expect_true(general$model$control$useMatrix)
  expect_false(is.null(attr(general$model, "maihda_rescaled")))
  expect_true(general$diagnostics$converged)
  for (made_here in c("start", "iter.max", "useMatrix", "method")) {
    expect_null(general$model$call[[made_here]])
  }
  again <- stats::update(general$model, evaluate = FALSE)
  expect_s3_class(suppressWarnings(eval(again, list(data = general$data))), "clmm")
  # The boundary fit the package makes with that routine stands for its call.
  flat <- cs14_fit(cs14_strata(), cs14_formula("birth"))
  expect_identical(flat$model$call[["start"]], quote(.maihda_clmm_start))
  expect_true(flat$model$call[["useMatrix"]])
})

test_that("a boundary fit vouched for by the first fit keeps that verdict", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata()
  d$stratum <- factor(paste(d$gender, d$ses))
  env <- new.env()
  env$data <- d
  formula <- yo ~ birth + (1 | stratum)
  environment(formula) <- env
  args <- list(formula = formula, data = quote(data), link = "logit", Hess = TRUE)
  design <- maihda_clmm_rescaled_design(args, list(), env)
  rescaled <- suppressWarnings(eval(design$call, design$env))
  attr(rescaled, "maihda_converged") <- TRUE
  kept <- maihda_clmm_carry_back(rescaled, design, args, list(), env,
                                 as.call(c(list(quote(ordinal::clmm)), args)))
  expect_true(attr(kept, "maihda_converged"))
  attr(rescaled, "maihda_converged") <- NULL
  plain <- maihda_clmm_carry_back(rescaled, design, args, list(), env,
                                  as.call(c(list(quote(ordinal::clmm)), args)))
  expect_null(attr(plain, "maihda_converged"))
})

test_that("a fit that cannot be carried back is made as the formula stands", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  local_mocked_bindings(maihda_clmm_carry_back = function(...) stop("no carry"))
  heard <- cs14_observe(fit_maihda(cs14_formula("birth"), data = d, engine = "ordinal",
                                   family = "ordinal"))
  m <- heard$res
  expect_null(attr(m$model, "maihda_rescaled"))
  # clmm()'s own search in birth years, as before: its start fit's warning, once.
  expect_identical(sum(grepl("Rescale variables", heard$warns, fixed = TRUE)), 1L)
  expect_false(m$diagnostics$converged)
})

test_that("ucminf(), quadrature and structured thresholds on rescaled columns", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  same <- function(m, reference) {
    expect_true(m$diagnostics$converged)
    expect_equal(cs14_loglik(m), cs14_loglik(reference), tolerance = 1e-8)
    expect_equal(cs14_variance(m), cs14_variance(reference), tolerance = 1e-4)
  }
  by_ucminf <- ordinal::clmm.control(method = "ucminf")
  m <- cs14_fit(d, cs14_formula("birth"), control = by_ucminf)
  same(m, cs14_fit(d, cs14_formula("z"), control = by_ucminf))
  expect_true(maihda_is_ucminf_return(m$model$optRes))
  expect_identical(m$model$control$method, "ucminf")
  expect_lt(maihda_clmm_stationary_distance(m$model), 0.01)

  same(cs14_fit(d, cs14_formula("income"), nAGQ = 7),
       cs14_fit(d, cs14_formula("z"), nAGQ = 7))
  same(cs14_fit(d, cs14_formula("birth"), family = maihda_cumulative("probit")),
       cs14_fit(d, cs14_formula("z"), family = maihda_cumulative("probit")))

  five <- cs14_strata(cs14_shift, cuts = c(-1, -0.4, 0.2, 0.8))
  for (structure in c("equidistant", "symmetric")) {
    same(cs14_fit(five, cs14_formula("birth"), threshold = structure),
         cs14_fit(five, cs14_formula("z"), threshold = structure))
  }
  # "symmetric2" fixes the latent mean, so a shift of the covariate is another model; a
  # change of its scale is not.
  five$birth_12 <- five$birth / 12
  two <- cs14_fit(five, cs14_formula("birth"), threshold = "symmetric2")
  same(two, cs14_fit(five, cs14_formula("birth_12"), threshold = "symmetric2"))
  expect_identical(attr(two$model, "maihda_rescaled")$centre, 0)
})

test_that("a squared covariate, a second covariate and an offset", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  d$expo <- rep(c(1, 2, 4), length.out = nrow(d))
  squared <- cs14_fit(d, yo ~ gender + ses + age + I(age^2) + (1 | gender:ses))
  orthogonal <- cs14_fit(d, yo ~ gender + ses + poly(z, 2) + (1 | gender:ses))
  expect_true(squared$diagnostics$converged)
  expect_equal(cs14_loglik(squared), cs14_loglik(orthogonal), tolerance = 1e-8)
  with_offset <- cs14_fit(d, yo ~ income + offset(log(expo)) + (1 | gender:ses))
  expect_true(with_offset$diagnostics$converged)
  expect_false(is.null(attr(with_offset$model, "maihda_rescaled")))
  direct <- ordinal::clmm(yo ~ z + offset(log(expo)) + (1 | stratum), data = with_offset$data,
                          link = "logit", Hess = TRUE)
  expect_equal(cs14_loglik(with_offset), as.numeric(stats::logLik(direct)), tolerance = 1e-8)
  expect_gt(abs(cs14_loglik(with_offset) - cs14_loglik(cs14_fit(d, cs14_formula("z")))), 1)
  # An aliased column is dropped from both fits alike.
  d$age_twice <- 2 * d$age + 1
  aliased <- cs14_fit(d, yo ~ age + age_twice + (1 | gender:ses))
  expect_identical(names(aliased$model$beta), "age")
  expect_equal(cs14_loglik(aliased), cs14_loglik(cs14_fit(d, cs14_formula("z"))),
               tolerance = 1e-8)
})

test_that("a fit given a start is made as the formula stands", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  top <- cs14_fit(d, cs14_formula("z"))
  at <- as.numeric(top$model$optRes$par)
  m <- cs14_fit(d, cs14_formula("z"), start = at)
  expect_null(attr(m$model, "maihda_rescaled"))
  expect_identical(as.numeric(m$model$start), at)
  # A model with no fixed effect has nothing to rescale.
  plain <- cs14_fit(d, yo ~ (1 | gender:ses))
  expect_null(attr(plain$model, "maihda_rescaled"))
})

test_that("a rescaled search cut short is recorded in its own words", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  reference <- cs14_fit(d, cs14_formula("birth"))
  for (control in list(ordinal::clmm.control(iter.max = 2),
                       ordinal::clmm.control(method = "ucminf", maxeval = 5))) {
    m <- cs14_fit(d, cs14_formula("birth"), control = control)
    expect_false(m$diagnostics$converged)
    expect_identical(m$diagnostics$messages, as.character(m$model$optRes$message))
    expect_gt(cs14_loglik(reference) - cs14_loglik(m), 1e-3)
  }
  expect_identical(
    cs14_fit(d, cs14_formula("birth"), control = ordinal::clmm.control(iter.max = 2))$
      diagnostics$messages,
    "iteration limit reached without convergence (10)")
})

test_that("what a fit on rescaled columns says, it says once", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs14_strata(cs14_shift)
  heard <- cs14_observe(fit_maihda(cs14_formula("birth"), data = d, engine = "ordinal",
                                   family = "ordinal"))
  # ordinal's "Rescale variables?" warning was about its fixed-effects start fit, which
  # is now made on the rescaled columns.
  expect_false(any(grepl("Rescale variables", heard$warns, fixed = TRUE)))
  expect_length(heard$warns, 0L)
  # A column clmm() drops is announced by the fit on rescaled columns, once.
  d$age_twice <- 2 * d$age + 1
  twice <- cs14_observe(fit_maihda(yo ~ age + age_twice + (1 | gender:ses), data = d,
                                   engine = "ordinal", family = "ordinal"))
  dropped <- grepl("dropping", c(twice$warns, twice$notes), fixed = TRUE)
  expect_identical(sum(dropped), 1L)
})
