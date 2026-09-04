# Regression tests for the 2026-09-03 audit finding.
#
#   [Med] "Non-converged bootstrap fits are retained in intervals." The parametric
#         PCV/VPC bootstrap RETAINS draws whose refit optimiser did not converge
#         (optinfo$conv$opt != 0), counting and warning about them but never
#         invalidating the interval no matter how large their share -- unlike the
#         boundary and hard-failure shares, both of which are gated. maihda_bootstrap_ci()
#         is blind to non-convergence (those draws are finite, so they count as
#         n_ok/successes), so a mostly-non-converged bootstrap cleared the gate
#         unremarked.
#
#   VERDICT: PARTIAL. The mechanism is real (confirmed by injection below) but the
#   severity premise is overstated: optinfo$conv$opt != 0 fires on ~0% of refits even
#   for deliberately hard GLMM fits (the frequent near-boundary/singular draws are a
#   SEPARATE, already-gated path). FIX (user's call: warn + flag, still return; no
#   converged-only sensitivity interval): a shared maihda_report_nonconvergence()
#   helper emits the existing retained-note at/below a documented 50% share and an
#   escalated warning above it, returning interval_reliable, which is attached to the
#   PCV, VPC, crossed-dimensions, contextual, longitudinal and PCV-importance results
#   and surfaced by print(). Retention itself (the 2026-07-31 finding #4 decision) is
#   preserved.

# ---- The reliability helper: threshold logic, deterministic (no fit) ----------

test_that("maihda_report_nonconvergence flags above the threshold, not at or below", {
  f <- MAIHDA:::maihda_report_nonconvergence

  # Above 50%: unreliable, escalated warning.
  expect_warning(hi <- f(60, 100, "PCV"), "reliability threshold")
  expect_false(hi)

  # At/below 50%: reliable, standard retained-note.
  expect_warning(lo <- f(40, 100, "PCV"), "retained in the interval")
  expect_true(lo)

  # Exactly 50% is NOT above "more than half" -> reliable.
  expect_warning(mid <- f(50, 100, "PCV"), "retained in the interval")
  expect_true(mid)

  # Zero non-converged: reliable and SILENT (no warning at all).
  expect_no_warning(z <- f(0, 100, "PCV"))
  expect_true(z)

  # Degenerate inputs are treated as "nothing to report" (reliable, silent).
  expect_no_warning(expect_true(f(5, 0, "PCV")))
  expect_no_warning(expect_true(f(NA_integer_, 100, "PCV")))

  # The threshold is a parameter: 30% is unreliable when the ceiling is 20%.
  expect_warning(strict <- f(30, 100, "PCV", threshold = 0.2), "reliability threshold")
  expect_false(strict)
})

# ---- PCV bootstrap: high non-convergence flags but does NOT refuse the interval ---

# Force EVERY refit to report non-convergence via the low-false-positive detector the
# bootstrap already consults. This is the deterministic way to reach the bad branch:
# optinfo$conv$opt != 0 fires on ~0% of real refits, so we inject it rather than hope
# an optimiser stops early.
make_pcv_pair <- function(seed = 11, n_strata = 12, n_per = 12) {
  set.seed(seed)
  se <- stats::rnorm(n_strata, 0, 1.2)
  d <- data.frame(stratum = rep(seq_len(n_strata), each = n_per),
                  age = stats::rnorm(n_strata * n_per),
                  gender = sample(0:1, n_strata * n_per, TRUE))
  d$outcome <- 5 + 0.5 * d$age + se[d$stratum] + stats::rnorm(nrow(d))
  list(
    m1 = fit_maihda(outcome ~ age + (1 | stratum), data = d, engine = "lme4"),
    m2 = fit_maihda(outcome ~ age + gender + (1 | stratum), data = d, engine = "lme4")
  )
}

test_that("a fully non-converged PCV bootstrap is flagged unreliable but still returned", {
  skip_on_cran()
  mm <- make_pcv_pair()

  testthat::local_mocked_bindings(maihda_lme4_optimizer_failed = function(model) TRUE)

  # The interval is RETURNED (the user's decision: warn + flag, do not refuse) with an
  # escalated warning -- the gate is deliberately not turned into a hard error here.
  expect_warning(
    res <- calculate_pcv(mm$m1, mm$m2, bootstrap = TRUE, n_boot = 60),
    "reliability threshold")

  expect_true(is.finite(res$ci_lower) && is.finite(res$ci_upper))
  expect_false(res$interval_reliable)
  # Every contributing draw is counted non-converged, and it never exceeds n_boot_ok.
  expect_identical(res$n_boot_nonconverged, res$n_boot_ok)
  expect_true(res$n_boot_nonconverged > 0)

  # print() discloses both the retained count and the unreliable flag (the count was
  # stored but invisible before this fix).
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "did not converge")
  expect_match(out, "UNRELIABLE")
})

test_that("a converged PCV bootstrap is rated reliable and prints no unreliable flag", {
  skip_on_cran()
  mm <- make_pcv_pair()
  # No mock: the real detector reports conv$opt, which does not fire here.
  res <- suppressWarnings(calculate_pcv(mm$m1, mm$m2, bootstrap = TRUE, n_boot = 60))
  expect_true(isTRUE(res$interval_reliable))
  expect_identical(res$n_boot_nonconverged, 0L)
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_false(grepl("UNRELIABLE", out))
})

# ---- VPC summary path carries the same flag -----------------------------------

test_that("a fully non-converged VPC bootstrap summary is flagged unreliable", {
  skip_on_cran()
  set.seed(21)
  n_strata <- 12; n_per <- 12
  se <- stats::rnorm(n_strata, 0, 1.2)
  d <- data.frame(stratum = rep(seq_len(n_strata), each = n_per),
                  x = stats::rnorm(n_strata * n_per))
  d$y <- 2 + 0.4 * d$x + se[d$stratum] + stats::rnorm(nrow(d))
  m <- fit_maihda(y ~ x + (1 | stratum), data = d, engine = "lme4")

  testthat::local_mocked_bindings(maihda_lme4_optimizer_failed = function(model) TRUE)
  s <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 40))

  expect_false(s$vpc$interval_reliable)
  expect_true(is.finite(s$vpc$ci_lower) && is.finite(s$vpc$ci_upper))
  expect_true(s$vpc$n_boot_nonconverged > 0)
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "UNRELIABLE")
})
