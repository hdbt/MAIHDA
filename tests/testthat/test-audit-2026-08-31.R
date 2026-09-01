# Audit 2026-08-31 -- "some opt-in Bayesian tests accept non-converged fits".
#
# The finding was CONFIRMED: the opt-in Stan blocks fitted models that had not
# mixed (max Rhat 1.08 with a bulk ESS of 26 on the ordinal thresholds; 1.075
# with an ESS of 36 on the tidiers fit) and asserted only that numbers came back
# finite and in [0, 1]. Those blocks are now tuned to convergence and each one
# asserts isTRUE(model$diagnostics$converged).
#
# That makes maihda_fit_diagnostics()'s brms verdict load-bearing for the whole
# opt-in suite: if it ever returned TRUE unconditionally, every one of those new
# guards would pass vacuously and the suite would silently go back to proving
# only that the code runs. test-audit-2026-07-26.R already pins the verdict's
# direction (clean Rhat -> TRUE, 1.05 -> FALSE, divergences -> FALSE, missing
# Rhat -> NA). What it does not pin is the THRESHOLD itself, and the threshold is
# the whole reason the guard has teeth: the package cuts at 1.01, while brms only
# warns above 1.05, so there is a band -- exactly where these fixtures used to
# sit -- in which brms says nothing and the package must still say "not
# converged". Pin the boundary and that band here.

# The verdict maihda_fit_diagnostics() reaches for a brmsfit with the given Rhat
# values and divergence count, without fitting anything.
brms_verdict_at <- function(rhat_val, div_val = 0) {
  testthat::local_mocked_bindings(
    rhat = function(x, ...) rhat_val,
    nuts_params = function(x, ...) data.frame(Parameter = "divergent__", Value = div_val),
    .package = "brms")
  MAIHDA:::maihda_fit_diagnostics(structure(list(), class = "brmsfit"))
}

test_that("the brms convergence verdict cuts at Rhat 1.01, not brms's 1.05", {
  skip_if_not_installed("brms")

  # The comparison is strictly greater-than, so 1.01 exactly is still converged
  # and anything above it is not. Neither side of this boundary was covered.
  expect_true(brms_verdict_at(c(a = 1.0100))$converged)
  expect_false(brms_verdict_at(c(a = 1.0101))$converged)

  # The band that made this audit possible: brms itself warns only above 1.05, so
  # a fit in [1.01, 1.05] draws no Stan warning at all. The package must still
  # refuse to call it converged -- that is what the opt-in fixtures now rely on.
  for (r in c(1.02, 1.03, 1.04, 1.05)) {
    v <- brms_verdict_at(c(a = r))
    expect_false(v$converged)
    expect_match(v$messages, "may not have converged", all = FALSE)
  }

  # The verdict is over ALL parameters, not the first: a single bad stratum
  # effect among otherwise clean ones is still a failure. The old fixtures failed
  # exactly this way -- the thresholds and r_stratum[] entries were the bad ones
  # while sd_stratum__Intercept sat at 1.009.
  mixed <- c(b_Intercept = 1.001, sd_stratum__Intercept = 1.009,
             `r_stratum[1,Intercept]` = 1.053)
  expect_false(brms_verdict_at(mixed)$converged)
})

test_that("divergences alone can withhold the converged verdict", {
  skip_if_not_installed("brms")

  # Rhat clean but the sampler diverged. This is the half of the bar that forced
  # adapt_delta = 0.95 onto the 6-stratum tidiers fixture: its Rhat was already
  # 1.006 with a bulk ESS of 677, and 5 divergences in 4000 draws still (rightly)
  # withheld the verdict. A test asserting Rhat/ESS alone would have passed it.
  clean_rhat <- c(b_Intercept = 1.006, sd_stratum__Intercept = 1.005)
  expect_true(brms_verdict_at(clean_rhat, div_val = 0)$converged)

  div <- brms_verdict_at(clean_rhat, div_val = c(rep(0, 3995), rep(1, 5)))
  expect_false(div$converged)
  expect_match(div$messages, "5 divergent transition", all = FALSE)
})

# Found while retuning the opt-in fixtures above, not in the audit report: the
# longitudinal decomposition gave different verdicts run to run on IDENTICAL
# settings with seed = 1. The cause is an argument-name collision, not the
# sampler -- maihda() has a `seed` FORMAL (the response-scale VPC simulation
# seed), so seed = never lands in `...` and never reaches brms::brm().
# fit_maihda() has no such formal, so its seed does reach brm(): the same
# spelling means different things in the two entry points, and a user reading
# `seed` as "make this Bayesian fit reproducible" silently gets a different
# posterior every run. maihda() now warns; these pin that warning's scope.

test_that("maihda() warns that 'seed' does not seed the brms sampler", {
  skip_if_not_installed("brms")

  d <- data.frame(
    g1 = rep(c("a", "b"), each = 60),
    g2 = rep(c("x", "y", "z"), times = 40),
    x  = seq_len(120) / 120)
  d$y <- 1 + 0.3 * d$x + as.numeric(factor(paste(d$g1, d$g2))) * 0.1

  # Collect the warnings an expression emits before it finishes OR errors.
  # warning() does not unwind, so a bare maihda(engine = "brms") call here would
  # go on to compile and run a Stan model; `offset` is rejected by fit_maihda()'s
  # brms handshake immediately after the warning point, which aborts the call
  # cheaply and keeps this block Stan-free.
  warns_of <- function(expr) {
    w <- character(0)
    tryCatch(
      withCallingHandlers(
        expr,
        warning = function(x) {
          w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
        }),
      error = function(e) NULL)
    w
  }
  seeded <- function(w) any(grepl("forwarded to brms::brm", w, fixed = TRUE))

  # brms + seed -> warn.
  expect_true(seeded(warns_of(
    maihda(y ~ x + (1 | g1:g2), data = d, engine = "brms", seed = 1,
           offset = rep(0, nrow(d))))))

  # SCOPE, both directions. brms WITHOUT a seed must stay silent: the warning is
  # about the argument collision, not about the engine.
  expect_false(seeded(warns_of(
    maihda(y ~ x + (1 | g1:g2), data = d, engine = "brms",
           offset = rep(0, nrow(d))))))

  # And a seed on a likelihood engine is the documented response-scale VPC
  # simulation seed, honoured there, so it must stay silent too --
  # test-discriminatory-accuracy.R and the introduction vignette both rely on
  # maihda(response_vpc = TRUE, seed = ) meaning what it says.
  expect_false(seeded(warns_of(
    maihda(y ~ x + (1 | g1:g2), data = d, engine = "lme4", seed = 1))))
})
