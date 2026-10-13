# Audit 2026-10-13 (packet 2026-10-13-ordinal-clmm-start-value): engine = "ordinal"
# at and near a between-stratum variance of zero -- CONFIRMED, both parts.
#
#   1. On strata with the same outcome proportions the fit reported a between-stratum
#      variance of exactly 1 and a VPC of 0.2331 (0.5 under probit), not flagged: six
#      identical strata of 30 gave variance 1 after one iteration with convergence
#      code 0, a group comparison showed the flat group at a VPC of 0.2331 beside
#      0.1965 for a group whose strata differ, and the two-model PCV was 0. The 1 is
#      ordinal::clmm()'s starting value. Its routine for a single random intercept
#      tests the stratum effects' convergence before it computes their Hessian;
#      effects that pass at their starting value of zero come back with the Hessian as
#      it was handed in, and the Laplace approximation is then the likelihood of the
#      model without the random effect, flat in the SD. It happens whenever every
#      stratum's score is zero at the pooled estimates: identical strata of any size,
#      strata that differ but balance, and any data once the fixed part spans the
#      strata (yo ~ gender * ses + (1 | gender:ses)). The ML variance is 0 each time.
#   2. A fit that reaches the boundary by itself stops where the optimiser's relative
#      tolerance is met, between 9e-10 and 2.4e-07 on 108 simulated boundary fits, so
#      14 of them were not flagged singular and 4 of 56 two-model analyses with a
#      boundary null printed a PCV of 0.88 to 0.97 from two boundary variances.
#
# Per the user: such a fit -- every stratum effect exactly zero, or an SD below 1% of
# the latent residual SD -- is made again by clmm()'s general routine,
# clmm.control(useMatrix = TRUE), which computes that Hessian and bounds the SD at
# zero, and that fit is kept; no threshold changes and no gradient test. Where the
# second fit cannot be made (quadrature) or fails, the fit is recorded as not
# converged. nlminb()'s "singular convergence" return from a variance of zero is the
# boundary fit, not a failure.
#
# How the second fit is made and when it is kept came out of measuring that routine.
# Left to its own start (SD 1) it comes to rest on its bound of zero for small
# variances that are not zero as well (ML SDs of 0.013 to 0.05 on six strata of 200
# rows), and with a covariate on a raw scale (age in years) it uses up nlminb()'s 150
# iterations in 17 of 23 boundary fits. So it is started from the first fit's fixed
# part and an SD of zero, where it returns the boundary fit in one to four iterations,
# and that fit is kept only with a log-likelihood not below the first fit's by more
# than 1e-9 of its size. Started on the fit it is to return, nlminb() reports "false
# convergence" for most models with such a covariate; a boundary fit is then kept on
# the first fit's own convergence.
#
# The live blocks assert what holds whether or not clmm()'s first routine misbehaves:
# the package's fit, never clmm()'s own return. No block draws a random number except
# through cs13_keep_seed(), which puts the caller's seed back.

cs13_observe <- function(expr) {
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(res = res, warns = warns)
}

cs13_keep_seed <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  force(expr)
}

cs13_cut <- function(y) {
  factor(findInterval(y, c(-0.8, 0, 0.8)) + 1L, levels = 1:4, ordered = TRUE)
}

# gender x ses strata of `n` rows: a stratum's outcomes are its shift plus the same
# fixed spread, cut into four classes, so shifts of zero give identical strata.
cs13_strata <- function(shift = rep(0, 6), n = 30, gender = c("F", "M"),
                        ses = c("low", "mid", "high")) {
  g <- expand.grid(gender = gender, ses = ses, stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = n), , drop = FALSE]
  d$y <- rep(shift, each = n) + rep(seq(-1.45, 1.45, length.out = n), times = nrow(g))
  d$yo <- cs13_cut(d$y)
  rownames(d) <- NULL
  d
}
cs13_shift <- c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2)

cs13_fit <- function(d, formula = yo ~ (1 | gender:ses), ...) {
  suppressMessages(suppressWarnings(
    fit_maihda(formula, data = d, engine = "ordinal", ...)))
}
cs13_variance <- function(m) maihda_clmm_variances(m)$stratum

# A stand-in for a clmm fit holding what the helpers read: the SD, the stratum
# effects, which routine made it, the quadrature setting, the link and the
# optimiser's return.
cs13_record <- function(sd, effects, general = FALSE, n_quad = 1L, link = "logit",
                        code = 0L, message = "relative convergence (4)",
                        loglik = -250) {
  structure(list(
    ST = list(matrix(sd, 1, 1, dimnames = list("(Intercept)", "(Intercept)"))),
    dims = list(nlev.re = c(stratum = length(effects))),
    coefficients = c("1|2" = -1.2, "2|3" = 0, "3|4" = 1.2),
    ranef = effects, control = list(useMatrix = general), nAGQ = n_quad, link = link,
    optRes = list(convergence = code, message = message), logLik = loglik),
    class = "clmm")
}
cs13_zero <- rep(0, 6)
cs13_some <- c(0.3, -0.2, 0.1, 0, -0.4, 0.2)

test_that("which fits are made again", {
  skip_if_not_installed("ordinal")
  latent <- pi / sqrt(3)

  expect_equal(maihda_clmm_stratum_sd(cs13_record(0.3, cs13_some)), 0.3)
  expect_identical(maihda_clmm_stratum_sd(cs13_record(0, cs13_zero, general = TRUE)), 0)
  expect_true(is.na(maihda_clmm_stratum_sd(structure(list(), class = "clmm"))))

  # Every effect exactly zero, from the single-intercept routine.
  expect_true(maihda_clmm_effects_at_start(cs13_record(1, cs13_zero)))
  expect_false(maihda_clmm_effects_at_start(cs13_record(1, c(cs13_zero[-1], 1e-300))))
  expect_false(maihda_clmm_effects_at_start(cs13_record(0, cs13_zero, general = TRUE)))
  expect_false(maihda_clmm_effects_at_start(cs13_record(1, numeric(0))))
  no_control <- cs13_record(1, cs13_zero)
  no_control$control <- NULL
  expect_false(maihda_clmm_effects_at_start(no_control))
  expect_false(maihda_clmm_effects_at_start(list(ranef = cs13_zero)))

  # The stuck fit, whatever its SD.
  expect_true(maihda_clmm_refit_wanted(cs13_record(1, cs13_zero)))
  expect_true(maihda_clmm_refit_wanted(cs13_record(0.37, cs13_zero)))
  # A fit with effects and an SD well away from zero is left alone.
  expect_false(maihda_clmm_refit_wanted(cs13_record(0.5, cs13_some)))
  # Below 1% of the latent residual SD, on the fit's own latent scale.
  expect_true(maihda_clmm_refit_wanted(cs13_record(0.0099 * latent, cs13_some)))
  expect_false(maihda_clmm_refit_wanted(cs13_record(0.0101 * latent, cs13_some)))
  expect_true(maihda_clmm_refit_wanted(cs13_record(0.0099, cs13_some, link = "probit")))
  expect_false(maihda_clmm_refit_wanted(cs13_record(0.0101, cs13_some, link = "probit")))
  expect_true(maihda_clmm_refit_wanted(cs13_record(0.015, cs13_some)))
  expect_false(maihda_clmm_refit_wanted(cs13_record(0.015, cs13_some, link = "probit")))
  # A fit of the general routine is not made a third time.
  expect_false(maihda_clmm_refit_wanted(cs13_record(0, cs13_zero, general = TRUE)))
  expect_false(maihda_clmm_refit_wanted(cs13_record(1e-5, cs13_some, general = TRUE)))
  # Quadrature belongs to the single-intercept routine alone.
  expect_false(maihda_clmm_refit_wanted(cs13_record(1, cs13_zero, n_quad = 7L)))
  expect_false(maihda_clmm_refit_wanted(cs13_record(1, cs13_zero, n_quad = -5L)))
  expect_false(maihda_clmm_refit_wanted(cs13_record(1e-5, cs13_some, n_quad = 3L)))
  expect_true(maihda_clmm_refit_wanted(cs13_record(1, cs13_zero, n_quad = 0L)))
  expect_false(maihda_clmm_refit_wanted(list(ranef = cs13_zero)))
})

test_that("the second call changes the routine and nothing the caller chose", {
  args <- list(formula = yo ~ (1 | stratum), data = quote(data), link = "probit",
               Hess = TRUE)
  cl <- maihda_clmm_general_call(args, list(
    nAGQ = 1, start = c(0, 0, 0, 0), threshold = "equidistant", model = TRUE,
    control = list(maxIter = 200), useMatrix = FALSE, method = "ucminf",
    meth = "ucminf", gradTol = 1e-6))
  expect_identical(cl[[1]], quote(ordinal::clmm))
  expect_identical(names(cl)[-1],
                   c("formula", "data", "link", "Hess", "threshold", "model", "control",
                     "gradTol", "start", "useMatrix", "method"))
  expect_true(cl$useMatrix)
  expect_identical(cl$method, "nlminb")
  expect_identical(cl$link, "probit")
  expect_identical(cl$data, quote(data))
  expect_identical(cl$threshold, "equidistant")
  expect_identical(cl$control, list(maxIter = 200))
  expect_identical(cl$gradTol, 1e-6)
  # The caller's start is gone; the one that takes its place is named, not written out.
  expect_identical(cl$start, quote(.maihda_clmm_start))

  bare <- maihda_clmm_general_call(args, list())
  expect_identical(names(bare)[-1],
                   c("formula", "data", "link", "Hess", "start", "useMatrix", "method"))
  # Unnamed dots are passed on as they came.
  unnamed <- maihda_clmm_general_call(args, list(3))
  expect_identical(names(unnamed)[-1],
                   c("formula", "data", "link", "Hess", "", "start", "useMatrix",
                     "method"))

  # The second fit starts from the first fit's fixed part, on the boundary.
  expect_identical(maihda_clmm_general_start(cs13_record(1, cs13_zero)),
                   list(c(-1.2, 0, 1.2), 0))
  expect_identical(
    maihda_clmm_general_start(cs13_record(3e-5, cs13_some, link = "probit")),
    list(c(-1.2, 0, 1.2), 0))
})

test_that("the optimiser's return is read as a verdict", {
  skip_if_not_installed("ordinal")
  expect_identical(maihda_nlminb_port_code("singular convergence (7)"), 7L)
  expect_identical(maihda_nlminb_port_code("relative convergence (4)"), 4L)
  expect_identical(
    maihda_nlminb_port_code("both X-convergence and relative convergence (5)"), 5L)
  expect_identical(
    maihda_nlminb_port_code("iteration limit reached without convergence (10)"), 10L)
  expect_identical(maihda_nlminb_port_code("Stopped by small gradient (grtol)."),
                   NA_integer_)
  expect_identical(maihda_nlminb_port_code("(7) singular convergence"), NA_integer_)
  expect_identical(maihda_nlminb_port_code(NA_character_), NA_integer_)
  expect_identical(maihda_nlminb_port_code(c("a (7)", "b (7)")), NA_integer_)
  expect_identical(maihda_nlminb_port_code(7), NA_integer_)
  expect_identical(maihda_nlminb_port_code(NULL), NA_integer_)

  expect_true(maihda_clmm_converged(cs13_record(0.4, cs13_some)))
  # nlminb()'s "singular convergence" from an SD of exactly zero is the boundary fit.
  at_zero <- function(message, sd = 0) {
    cs13_record(sd, cs13_zero, general = TRUE, code = 1L, message = message)
  }
  expect_true(maihda_clmm_converged(at_zero("singular convergence (7)")))
  # The same below the variance of 1e-8 that marks a singular fit ...
  expect_true(maihda_clmm_converged(at_zero("singular convergence (7)", sd = 5e-5)))
  # ... but not above it, and no other return.
  expect_false(maihda_clmm_converged(at_zero("singular convergence (7)", sd = 2e-4)))
  expect_false(maihda_clmm_converged(at_zero("singular convergence (7)", sd = 0.3)))
  expect_false(maihda_clmm_converged(at_zero("false convergence (8)")))
  expect_false(maihda_clmm_converged(
    at_zero("iteration limit reached without convergence (10)")))
  expect_false(maihda_clmm_converged(at_zero("Stopped by small gradient (grtol).")))
  no_return <- cs13_record(0.4, cs13_some)
  no_return$optRes <- NULL
  expect_true(is.na(maihda_clmm_converged(no_return)))
  expect_true(is.na(maihda_clmm_converged(cs13_record(0.4, cs13_some,
                                                      code = NA_integer_))))

  # A fit kept on the first fit's convergence carries that verdict as a mark.
  vouched <- at_zero("false convergence (8)")
  attr(vouched, "maihda_converged") <- TRUE
  expect_true(maihda_clmm_converged(vouched))

  expect_true(maihda_clmm_at_boundary(cs13_record(0, cs13_zero, general = TRUE)))
  expect_true(maihda_clmm_at_boundary(cs13_record(5e-5, cs13_some)))
  expect_false(maihda_clmm_at_boundary(cs13_record(2e-4, cs13_some)))
  expect_false(maihda_clmm_at_boundary(structure(list(), class = "clmm")))

  # A second fit can take the first one's place when it converged ...
  expect_true(maihda_clmm_refit_usable(cs13_record(0, cs13_zero, general = TRUE)))
  expect_true(maihda_clmm_refit_usable(cs13_record(0.2, cs13_some, general = TRUE)))
  expect_true(maihda_clmm_refit_usable(at_zero("singular convergence (7)")))
  # ... or sits on the boundary: nlminb() cannot always say of a fit it was started on
  # that it converged.
  expect_true(maihda_clmm_refit_usable(at_zero("false convergence (8)")))
  expect_true(maihda_clmm_refit_usable(
    at_zero("iteration limit reached without convergence (10)")))
  expect_true(maihda_clmm_refit_usable(at_zero("false convergence (8)", sd = 5e-5)))
  # Off the boundary only its own convergence counts.
  expect_false(maihda_clmm_refit_usable(at_zero("false convergence (8)", sd = 0.2)))
  expect_false(maihda_clmm_refit_usable(at_zero("false convergence (8)", sd = 2e-4)))
  expect_false(maihda_clmm_refit_usable(
    cs13_record(0, cs13_zero, general = TRUE, loglik = -Inf)))
  expect_false(maihda_clmm_refit_usable(cs13_record(NaN, cs13_zero, general = TRUE)))
  expect_false(maihda_clmm_refit_usable(NULL))
  expect_false(maihda_clmm_refit_usable(list(ST = list(matrix(0)))))

  # ... and only with a log-likelihood that is not below the first fit's, to 1e-9 of
  # its size.
  second <- function(loglik) cs13_record(0, cs13_zero, general = TRUE, loglik = loglik)
  first <- function(loglik) cs13_record(0.01, cs13_some, loglik = loglik)
  expect_true(maihda_clmm_refit_no_worse(second(-250), first(-250)))
  expect_true(maihda_clmm_refit_no_worse(second(-249), first(-250)))
  expect_true(maihda_clmm_refit_no_worse(second(-250.0000002), first(-250)))
  expect_false(maihda_clmm_refit_no_worse(second(-250.0000003), first(-250)))
  expect_true(maihda_clmm_refit_no_worse(second(-25000.00002), first(-25000)))
  expect_false(maihda_clmm_refit_no_worse(second(-25000.00003), first(-25000)))
  expect_true(maihda_clmm_refit_no_worse(second(-0.5 - 8e-10), first(-0.5)))
  expect_false(maihda_clmm_refit_no_worse(second(-0.5 - 2e-9), first(-0.5)))
  # Nothing to compare with: the second fit stands. Nothing to compare: it does not.
  expect_true(maihda_clmm_refit_no_worse(second(-250), first(NA_real_)))
  no_loglik <- first(-250)
  no_loglik$logLik <- NULL
  expect_true(maihda_clmm_refit_no_worse(second(-250), no_loglik))
  expect_false(maihda_clmm_refit_no_worse(second(NA_real_), first(-250)))
})

test_that("the fit record of a fit the general routine could not take over", {
  skip_if_not_installed("ordinal")
  # A quadrature fit left at clmm()'s starting SD: not converged, and it says why.
  rec <- maihda_fit_diagnostics(cs13_record(1, cs13_zero, n_quad = 3L))
  expect_identical(rec$engine, "ordinal")
  expect_false(rec$converged)
  expect_false(rec$singular)
  expect_length(rec$messages, 1L)
  expect_match(rec$messages, "every stratum effect at exactly zero", fixed = TRUE)
  expect_match(rec$messages, "between-stratum SD of 1.", fixed = TRUE)
  expect_match(rec$messages, "With nAGQ = 1", fixed = TRUE)
  # The same state under the Laplace approximation (the second fit failed) does not
  # send the reader to a setting already in use.
  rec1 <- maihda_fit_diagnostics(cs13_record(0.37, cs13_zero))
  expect_false(rec1$converged)
  expect_match(rec1$messages, "between-stratum SD of 0.37.", fixed = TRUE)
  expect_false(grepl("nAGQ", rec1$messages, fixed = TRUE))
  # With the optimiser's own complaint, both are kept.
  rec2 <- maihda_fit_diagnostics(cs13_record(1, cs13_zero, n_quad = 3L, code = 1L,
                                             message = "false convergence (8)"))
  expect_false(rec2$converged)
  expect_length(rec2$messages, 2L)
  expect_identical(rec2$messages[1], "false convergence (8)")

  # Effects of zero beside an SD of zero are a boundary fit and nothing else.
  flat <- maihda_fit_diagnostics(cs13_record(5e-5, cs13_zero))
  expect_true(flat$singular)
  expect_true(flat$converged)
  expect_length(flat$messages, 0L)
  # The general routine's boundary return.
  bound <- maihda_fit_diagnostics(cs13_record(0, cs13_zero, general = TRUE, code = 1L,
                                              message = "singular convergence (7)"))
  expect_true(bound$singular)
  expect_true(bound$converged)
  expect_length(bound$messages, 0L)
  # Another non-zero return stays a failure, in the optimiser's words ...
  unvouched <- cs13_record(0, cs13_zero, general = TRUE, code = 1L,
                           message = "false convergence (8)")
  failed <- maihda_fit_diagnostics(unvouched)
  expect_false(failed$converged)
  expect_identical(failed$messages, "false convergence (8)")
  # ... unless the fit was kept on the first fit's convergence.
  attr(unvouched, "maihda_converged") <- TRUE
  kept <- maihda_fit_diagnostics(unvouched)
  expect_true(kept$converged)
  expect_true(kept$singular)
  expect_length(kept$messages, 0L)
  # A fit with effects and an SD away from zero is clean.
  clean <- maihda_fit_diagnostics(cs13_record(0.4, cs13_some))
  expect_true(clean$converged)
  expect_false(clean$singular)
  expect_length(clean$messages, 0L)

  # What print() puts above those lines claims no source for them.
  lines <- maihda_format_fit_diagnostics(rec)
  expect_identical(lines[1], "Convergence problems of the ordinal::clmm fit:")
  expect_identical(lines[2], paste0("  - ", rec$messages))
  expect_false(any(grepl("reported by", lines, fixed = TRUE)))
})

test_that("a second fit speaks only when it is kept, and never twice", {
  skip_if_not_installed("ordinal")
  first <- cs13_record(1, cs13_zero)
  args <- list(formula = yo ~ (1 | stratum), data = quote(data), link = "logit",
               Hess = TRUE)
  dots <- list(nAGQ = 1, start = c(0, 0, 0, 0), threshold = "equidistant")
  seen <- NULL
  answer <- NULL
  local_mocked_bindings(
    clmm = function(...) {
      seen <<- list(...)
      warning("said by both fits")
      message("noted by both fits")
      warning("said by the second fit alone")
      if (is.function(answer)) answer() else answer
    },
    .package = "ordinal")
  run <- function(first_fit = first) {
    warns <- character(0)
    notes <- character(0)
    out <- withCallingHandlers(
      maihda_clmm_general_refit(first_fit, args, dots, environment(),
                                c("said by both fits", "noted by both fits\n")),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        notes <<- c(notes, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    list(out = out, warns = warns, notes = notes)
  }

  # A converged second fit is kept as it is, with what it alone had to say.
  answer <- cs13_record(0, cs13_zero, general = TRUE)
  kept <- run()
  expect_identical(kept$out, answer)
  expect_identical(kept$warns, "said by the second fit alone")
  expect_length(kept$notes, 0L)
  # It was asked for the general routine, without the quadrature setting, and from the
  # first fit's fixed part and an SD of zero in place of the caller's start.
  expect_identical(names(seen), c("formula", "data", "link", "Hess", "threshold",
                                  "start", "useMatrix", "method"))
  expect_true(seen$useMatrix)
  expect_identical(seen$method, "nlminb")
  expect_identical(seen$threshold, "equidistant")
  expect_identical(seen$start, list(c(-1.2, 0, 1.2), 0))

  # A boundary fit nlminb() could not vouch for is kept on the first fit's
  # convergence, and tells the fit record so.
  answer <- cs13_record(0, cs13_zero, general = TRUE, code = 1L,
                        message = "false convergence (8)")
  vouched <- run()
  expect_true(attr(vouched$out, "maihda_converged"))
  expect_identical(vouched$out$optRes, answer$optRes)
  expect_identical(vouched$out$ST, answer$ST)
  rec <- maihda_fit_diagnostics(vouched$out)
  expect_true(rec$converged)
  expect_true(rec$singular)
  expect_length(rec$messages, 0L)
  # When the first fit did not converge either, the boundary fit is kept as it came,
  # and its record gives the optimiser's own return.
  astray <- cs13_record(1, cs13_zero, code = 1L, message = "false convergence (8)")
  unvouched <- run(astray)
  expect_identical(unvouched$out, answer)
  rec <- maihda_fit_diagnostics(unvouched$out)
  expect_false(rec$converged)
  expect_true(rec$singular)
  expect_identical(rec$messages, "false convergence (8)")

  # One that did not converge away from the boundary is dropped, and so is what it
  # said.
  answer <- cs13_record(0.2, cs13_some, general = TRUE, code = 1L,
                        message = "false convergence (8)")
  dropped <- run()
  expect_identical(dropped$out, first)
  expect_length(dropped$warns, 0L)
  expect_length(dropped$notes, 0L)
  # So is one that converged to a lower log-likelihood than the first fit's.
  answer <- cs13_record(0, cs13_zero, general = TRUE, loglik = -250.01)
  lower <- run()
  expect_identical(lower$out, first)
  expect_length(lower$warns, 0L)
  # A difference the optimisers cannot tell apart does not count.
  answer <- cs13_record(0, cs13_zero, general = TRUE, loglik = -250.0000001)
  expect_identical(run()$out, answer)
  # One that fails leaves the first fit in place, silently.
  answer <- function() stop("the second fit failed")
  failed <- run()
  expect_identical(failed$out, first)
  expect_length(failed$warns, 0L)
})

test_that("strata with the same outcomes: a variance of zero, not clmm()'s start", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  cases <- list(
    list(d = cs13_strata(), family = "ordinal", latent = pi^2 / 3),
    list(d = cs13_strata(n = 8), family = "ordinal", latent = pi^2 / 3),
    list(d = cs13_strata(shift = rep(0, 12), gender = c("a", "b", "c"),
                         ses = c("p", "q", "r", "s")),
         family = "ordinal", latent = pi^2 / 3),
    list(d = cs13_strata(), family = maihda_cumulative("probit"), latent = 1))
  for (case in cases) {
    m <- cs13_fit(case$d, family = case$family)
    v <- maihda_clmm_variances(m)
    expect_lt(v$stratum, 1e-8)
    expect_equal(v$residual, case$latent)
    expect_true(m$diagnostics$singular)
    expect_true(m$diagnostics$converged)
    expect_length(m$diagnostics$messages, 0L)
    # The fit that is kept is the general routine's.
    expect_true(m$model$control$useMatrix)
    obs <- cs13_observe(summary(m))
    expect_length(obs$warns, 0L)
    expect_lt(obs$res$vpc$estimate, 1e-8)
    expect_true(all(is.finite(obs$res$thresholds$se)))
    expect_true(all(obs$res$stratum_estimates$random_effect == 0))
    expect_true(all(obs$res$stratum_estimates$se == 0))
    # It is the fit of the model without a random effect.
    plain <- ordinal::clm(yo ~ 1, data = case$d, link = m$model$link)
    expect_equal(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(plain)),
                 tolerance = 1e-6)
    expect_equal(unname(maihda_clmm_cutpoints(m$model)), unname(as.numeric(plain$alpha)),
                 tolerance = 1e-3)
  }
  # clmm()'s own fit of the first case is one the package makes again.
  d <- cs13_strata()
  d$stratum <- interaction(d$gender, d$ses, drop = TRUE)
  expect_true(maihda_clmm_refit_wanted(
    suppressWarnings(ordinal::clmm(yo ~ 1 + (1 | stratum), data = d))))
})

test_that("strata that balance, and a fixed part that spans the strata", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  counts <- function(n, gender, ses) {
    data.frame(gender = gender, ses = ses,
               yo = factor(rep(1:4, times = n), levels = 1:4, ordered = TRUE))
  }
  # Three different outcome distributions, every one with a score of zero at the
  # pooled thresholds.
  balanced <- rbind(counts(c(7, 8, 8, 7), "F", "low"), counts(c(5, 10, 10, 5), "M", "low"),
                    counts(c(10, 5, 5, 10), "F", "mid"), counts(c(7, 8, 8, 7), "M", "mid"),
                    counts(c(5, 10, 10, 5), "F", "high"),
                    counts(c(10, 5, 5, 10), "M", "high"))
  m <- cs13_fit(balanced)
  expect_lt(cs13_variance(m), 1e-8)
  expect_true(m$diagnostics$singular)
  expect_true(m$diagnostics$converged)

  # Strata that do differ: their variance is there until the fixed part takes it all.
  d <- cs13_strata(cs13_shift)
  null <- cs13_fit(d)
  expect_gt(cs13_variance(null), 0.5)
  spanned <- cs13_fit(d, yo ~ gender * ses + (1 | gender:ses))
  expect_lt(cs13_variance(spanned), 1e-8)
  expect_true(spanned$diagnostics$singular)
  expect_true(spanned$diagnostics$converged)
  # The main effects alone leave the interaction to the strata.
  additive <- cs13_fit(d, yo ~ gender + ses + (1 | gender:ses))
  expect_gt(cs13_variance(additive), 0.1)
  expect_false(additive$diagnostics$singular)
})

test_that("a fit that reaches the boundary by itself is flagged and has no PCV", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Identical strata with one row moved up a class: clmm()'s first routine stops at a
  # variance of 3.4e-08, above the singular flag and the PCV's boundary test alike.
  d <- cs13_strata()
  d$yo[which(d$yo == "2")[1]] <- "3"
  m <- cs13_fit(d)
  expect_lt(cs13_variance(m), 1e-8)
  expect_true(m$diagnostics$singular)
  expect_true(m$diagnostics$converged)
  expect_length(cs13_observe(summary(m))$warns, 0L)

  obs <- cs13_observe(maihda(yo ~ (1 | gender:ses), data = d, engine = "ordinal"))
  expect_null(obs$res$pcv)
  expect_true(any(grepl("the PCV could not be computed", obs$warns, fixed = TRUE)))
  expect_true(obs$res$model$diagnostics$singular)

  # Simulated strata with little between them: clmm()'s first routine stops at a
  # variance of 1.04e-08, just above the singular flag.
  noise <- cs13_keep_seed({
    set.seed(5)
    x <- cs13_strata(shift = rep(0, 12), n = 20, gender = c("a", "b", "c"),
                     ses = c("p", "q", "r", "s"))
    y <- stats::rnorm(nrow(x), sd = 0.84)
    x$yo <- cs13_cut(y + rep(stats::rnorm(12, sd = 0.1), each = 20))
    x
  })
  mn <- cs13_fit(noise)
  expect_lt(cs13_variance(mn), 1e-8)
  expect_true(mn$diagnostics$singular)
  expect_true(mn$diagnostics$converged)
  expect_length(mn$diagnostics$messages, 0L)
})

test_that("a covariate on a scale far from 1 does not keep the boundary fit away", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Ages in years, the same thirty in every stratum. With the first set clmm()'s
  # general routine uses up nlminb()'s 150 iterations from its own start; with the
  # second, started on the boundary fit, it reports "false convergence" from the start
  # itself, as does the first routine.
  years <- list(
    c(53, 49, 55, 51, 55, 43, 46, 55, 71, 56, 19, 53, 23, 41, 63, 45, 53, 54, 54, 55,
      38, 61, 23, 36, 55, 47, 51, 53, 41, 66),
    c(61, 38, 35, 55, 55, 63, 62, 59, 57, 33, 45, 63, 41, 60, 54, 51, 45, 48, 58, 44,
      56, 40, 50, 54, 52, 44, 45, 40, 68, 47))
  with_age <- function(age, nudge = FALSE) {
    d <- cs13_strata()
    d$age <- rep(age, times = 6)
    if (nudge) {
      d$yo[which(d$yo == "2")[1]] <- "3"
    }
    d
  }
  for (d in list(with_age(years[[1]]), with_age(years[[2]]),
                 with_age(years[[1]], nudge = TRUE),
                 with_age(years[[2]], nudge = TRUE))) {
    m <- cs13_fit(d, yo ~ age + (1 | gender:ses))
    expect_lt(cs13_variance(m), 1e-8)
    expect_true(m$diagnostics$singular)
    expect_true(m$model$control$useMatrix)
    # Converged, or said not to be in the optimiser's words: never silently neither.
    expect_true(isTRUE(m$diagnostics$converged) ||
                  length(m$diagnostics$messages) > 0L)
    expect_false(any(grepl("every stratum effect at exactly zero",
                           m$diagnostics$messages, fixed = TRUE)))
    obs <- cs13_observe(summary(m))
    expect_length(obs$warns, 0L)
    expect_true(all(is.finite(obs$res$thresholds$se)))
    expect_true(all(is.finite(obs$res$fixed_effects$se)))
    # The fixed part is that of the model without a random effect.
    plain <- ordinal::clm(yo ~ age, data = d)
    expect_equal(as.numeric(stats::logLik(m$model)), as.numeric(stats::logLik(plain)),
                 tolerance = 1e-8)
    expect_equal(unname(m$model$beta), unname(plain$beta), tolerance = 1e-4)
  }
})

test_that("a small variance that is not zero stays as the first routine fitted it", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Six strata of 200 rows with class counts 45 / 55 / 55 / 45. In each pair of strata
  # `moved` rows go from class 2 to class 3 in one and from 3 to 2 in the other, so
  # the pooled counts never change and the strata differ by as little as wanted.
  pairs <- function(moved) {
    base <- c(45, 55, 55, 45)
    up <- function(k) c(0, -k, k, 0)
    counts <- rbind(base + up(moved[1]), base - up(moved[1]), base + up(moved[2]),
                    base - up(moved[2]), base + up(moved[3]), base - up(moved[3]))
    g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
                     stringsAsFactors = FALSE)
    do.call(rbind, lapply(seq_len(6), function(j) {
      data.frame(gender = g$gender[j], ses = g$ses[j],
                 yo = factor(rep(1:4, times = counts[j, ]), levels = 1:4, ordered = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }
  # An ML SD of 0.0129 (201-node quadrature), below 1% of the latent SD: the boundary
  # fit is 1.6e-04 lower in log-likelihood than the first routine's fit.
  m <- cs13_fit(pairs(c(16, 15, 12)))
  direct <- ordinal::clmm(yo ~ (1 | stratum), data = m$data, link = "logit", Hess = TRUE)
  expect_true(maihda_clmm_refit_wanted(direct))
  expect_equal(sqrt(cs13_variance(m)), 0.0129, tolerance = 0.05)
  expect_false(m$diagnostics$singular)
  expect_true(m$diagnostics$converged)
  expect_gte(as.numeric(stats::logLik(m$model)),
             as.numeric(stats::logLik(direct)) - 1e-5)

  # The same design with an ML SD of zero: the first routine stops short of it, above
  # the 1e-8 that flags a singular fit, and the boundary fit is the better of the two.
  m0 <- cs13_fit(pairs(c(20, 13, 7)))
  expect_lt(cs13_variance(m0), 1e-8)
  expect_true(m0$diagnostics$singular)
  expect_true(m0$diagnostics$converged)
  expect_true(m0$model$control$useMatrix)
})

test_that("what clmm() signals reaches the caller once, however many fits are made", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  real <- ordinal::clmm
  local_mocked_bindings(
    clmm = function(...) {
      warning("clmm speaking")
      message("clmm noting")
      cl <- match.call()
      cl[[1]] <- real
      eval(cl, parent.frame())
    },
    .package = "ordinal")
  heard <- function(d) {
    warns <- character(0)
    notes <- character(0)
    m <- withCallingHandlers(
      fit_maihda(yo ~ (1 | gender:ses), data = d, engine = "ordinal",
                 family = "ordinal"),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(msg) {
        notes <<- c(notes, conditionMessage(msg))
        invokeRestart("muffleMessage")
      })
    list(m = m, warns = sum(warns == "clmm speaking"),
         notes = sum(notes == "clmm noting\n"))
  }
  # Two fits: identical strata.
  flat <- heard(cs13_strata())
  expect_true(flat$m$model$control$useMatrix)
  expect_lt(cs13_variance(flat$m), 1e-8)
  expect_identical(flat$warns, 1L)
  expect_identical(flat$notes, 1L)
  # One fit: strata that differ.
  one <- heard(cs13_strata(cs13_shift))
  expect_false(one$m$model$control$useMatrix)
  expect_identical(one$warns, 1L)
  expect_identical(one$notes, 1L)
})

test_that("a fit away from the boundary is clmm()'s own, as before", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs13_strata(cs13_shift)
  m <- cs13_fit(d)
  expect_false(m$model$control$useMatrix)
  expect_false(m$diagnostics$singular)
  expect_true(m$diagnostics$converged)
  direct <- ordinal::clmm(yo ~ (1 | stratum), data = m$data, link = "logit", Hess = TRUE)
  for (part in c("coefficients", "ST", "logLik", "ranef", "condVar", "Hessian",
                 "gradient")) {
    expect_identical(m$model[[part]], direct[[part]])
  }
  expect_identical(m$model$optRes$par, direct$optRes$par)
  # The stored call is the one the caller's fit was made with.
  expect_null(m$model$call$useMatrix)
})

test_that("a flat ordinal group is compared as a flat lme4 group is", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  flat <- cs13_strata()
  flat$country <- "flat"
  differ <- cs13_strata(cs13_shift)
  differ$country <- "differ"
  d <- rbind(flat, differ)
  obs <- cs13_observe(compare_maihda_groups(yo ~ (1 | gender:ses), data = d,
                                            group = "country", engine = "ordinal"))
  res <- as.data.frame(obs$res)
  at <- function(group, column) res[[column]][res$group == group]
  expect_lt(at("flat", "var_between"), 1e-8)
  expect_lt(at("flat", "vpc"), 1e-8)
  expect_true(is.na(at("flat", "pcv")))
  expect_identical(at("flat", "pcv_status"), "failed")
  expect_identical(res$status, c("ok", "ok"))
  # The group whose strata differ keeps the numbers of its own fit.
  own <- cs13_fit(differ)
  expect_identical(at("differ", "var_between"), cs13_variance(own))
  expect_gt(at("differ", "vpc"), 0.15)
  expect_gt(at("differ", "vpc"), at("flat", "vpc"))
  expect_identical(at("differ", "pcv_status"), "ok")
  # The flat group is the one the fit-quality warning names.
  quality <- obs$warns[grepl("fit problems", obs$warns, fixed = TRUE)]
  expect_identical(quality,
                   maihda_group_fit_quality_message("ordinal", "flat", character(0)))
  expect_true(any(grepl("Between-stratum variance in model1 is zero or negative",
                        obs$warns, fixed = TRUE)))

  # The same rows on lme4, with the outcome as it was before it was cut.
  lin <- as.data.frame(cs13_observe(compare_maihda_groups(
    y ~ (1 | gender:ses), data = d, group = "country"))$res)
  expect_identical(lin$pcv_status[lin$group == "flat"], at("flat", "pcv_status"))
  expect_lt(lin$vpc[lin$group == "flat"], 1e-8)
})

test_that("what the caller passed reaches the second fit", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs13_strata()
  flat_fit <- function(...) {
    cs13_observe(fit_maihda(yo ~ (1 | gender:ses), data = d, engine = "ordinal", ...))
  }
  is_boundary <- function(m) {
    cs13_variance(m) < 1e-8 && isTRUE(m$diagnostics$singular) &&
      isTRUE(m$diagnostics$converged) && isTRUE(m$model$control$useMatrix)
  }
  quiet <- function(warns) warns[!grepl("ordered factor", warns, fixed = TRUE)]

  plain <- flat_fit()
  expect_true(is_boundary(plain$res))
  expect_length(quiet(plain$warns), 0L)
  # The stored call names the data and its start, and stays short.
  expect_identical(plain$res$model$call$data, quote(data))
  expect_identical(plain$res$model$call$start, quote(.maihda_clmm_start))
  expect_lt(length(deparse(plain$res$model$call)), 5)
  expect_true(plain$res$model$call$useMatrix)
  expect_identical(plain$res$model$start[[2]], 0)

  tuned <- flat_fit(control = ordinal::clmm.control(gradTol = 1e-6, maxIter = 200,
                                                    innerCtrl = "noWarn"))
  expect_true(is_boundary(tuned$res))
  expect_identical(tuned$res$model$control$ctrl$gradTol, 1e-6)
  expect_identical(tuned$res$model$control$ctrl$maxIter, 200)
  expect_identical(tuned$res$model$control$ctrl$innerCtrl, "noWarn")
  expect_identical(tuned$res$model$control$method, "nlminb")

  loose <- flat_fit(gradTol = 1e-6)
  expect_true(is_boundary(loose$res))
  expect_identical(loose$res$model$control$ctrl$gradTol, 1e-6)

  # Another optimiser for the first routine: the second uses its own, without a word.
  for (other in list(flat_fit(control = ordinal::clmm.control(method = "ucminf")),
                     flat_fit(method = "ucminf"))) {
    expect_true(is_boundary(other$res))
    expect_identical(other$res$model$control$method, "nlminb")
    expect_length(quiet(other$warns), 0L)
  }

  expect_true(is_boundary(flat_fit(start = c(-1.19, 0, 1.19, 0))$res))
  shaped <- flat_fit(threshold = "equidistant")
  expect_true(is_boundary(shaped$res))
  expect_identical(shaped$res$model$threshold, "equidistant")
  probit <- flat_fit(family = maihda_cumulative("probit"))
  expect_true(is_boundary(probit$res))
  expect_identical(probit$res$model$link, "probit")
  # A caller who asks for the general routine gets it in one fit.
  asked <- flat_fit(control = ordinal::clmm.control(useMatrix = TRUE))
  expect_true(is_boundary(asked$res))
  expect_null(asked$res$model$call$method)
})

test_that("a quadrature fit left at zero effects is never called converged", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- cs13_strata()
  for (n_quad in c(3, 5, 7)) {
    m <- cs13_fit(d, nAGQ = n_quad)
    # Quadrature stays with the routine that offers it.
    expect_false(m$model$control$useMatrix)
    expect_identical(m$model$nAGQ, as.integer(n_quad))
    stuck <- maihda_clmm_effects_at_start(m$model) && !isTRUE(m$diagnostics$singular)
    if (stuck) {
      expect_false(m$diagnostics$converged)
      expect_true(any(grepl("every stratum effect at exactly zero",
                            m$diagnostics$messages, fixed = TRUE)))
    } else {
      expect_true(m$diagnostics$singular)
    }
  }
})
