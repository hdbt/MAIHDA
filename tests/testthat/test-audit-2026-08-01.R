# Regression tests for the 2026-08-01 audit findings.
#
#   1 [High] pcv_importance(bootstrap = TRUE) excluded null-boundary / degenerate
#            draws from every attribution interval (keeping them out of the reducer's
#            success-fraction denominator via n_excluded) but disclosed NOTHING: no
#            warning and no returned boundary count, unlike calculate_pcv(). A heavy
#            boundary share could therefore clear the majority-of-eligible gate and
#            return a 95% interval from a handful of survivors, unremarked. It now
#            warns and reports n_boot_boundary, mirroring calculate_pcv().
#   2 [Med]  Optimizer-failure counting (NEWS: "bootstrap draws whose refit optimizer
#            did not converge are counted and reported") was wired into the standard
#            VPC, PCV, and PCV-importance bootstraps only. The crossed-dimensions,
#            contextual, and longitudinal bootstraps never called
#            maihda_lme4_optimizer_failed(), so their n_boot_ok could include
#            non-converged refits with no n_boot_nonconverged and no warning. The
#            longitudinal per-time band also used a bare length(col) >= 10L, bypassing
#            maihda_bootstrap_ci()'s majority-of-eligible rule for time-specific
#            numerical failures. All three now count non-convergence, and the band is
#            gated by the same floor + majority rule as the headline.

# ---- Finding 1: pcv_importance() discloses boundary exclusions ----------------

# Near-boundary fixture: the null model's between-stratum variance sits just above
# the zero boundary, so a sizeable share of bootstrap null refits collapse onto it
# and are legitimately excluded from the attribution. Deterministic given the seed.
make_boundary_pcv_data <- function(seed = 1, sd_u = 0.25) {
  set.seed(seed)
  strata <- expand.grid(d1 = c("a", "b", "c"), d2 = c("p", "q", "r", "s"),
                        stringsAsFactors = FALSE)
  dat <- strata[rep(seq_len(nrow(strata)), each = 10), ]
  dat$x <- rnorm(nrow(dat))
  rownames(dat) <- NULL
  dat <- make_strata(dat, c("d1", "d2"))$data
  u <- rnorm(max(dat$stratum), sd = sd_u)
  dat$y <- 0.5 * dat$x + u[dat$stratum] + rnorm(nrow(dat), sd = 1)
  dat
}

# Collect every warning an expression emits (expect_warning only inspects one).
collect_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr,
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  w
}

test_that("pcv_importance() reports n_boot_boundary and warns when boundary draws are excluded", {
  skip_on_cran()
  dat <- make_boundary_pcv_data()
  ws <- character(0)
  res <- withCallingHandlers(
    suppressMessages(pcv_importance(dat, "y", c("d1", "d2"), method = "shapley",
                                    bootstrap = TRUE, n_boot = 300)),
    warning = function(c) { ws <<- c(ws, conditionMessage(c)); invokeRestart("muffleWarning") })

  # The count is now part of the bootstrap result's contract (was absent before).
  expect_true("n_boot_boundary" %in% names(res))
  # This construction reliably drives many null refits onto the boundary.
  expect_gt(res$n_boot_boundary, 0)
  # ... and the exclusion is now surfaced as a warning, as calculate_pcv() does.
  expect_true(any(grepl("boundary|conditional on a positive null variance",
                        ws, ignore.case = TRUE)))
})

test_that("a healthy pcv_importance() bootstrap carries n_boot_boundary = 0 and does not over-warn", {
  skip_on_cran()
  set.seed(3)
  strata <- expand.grid(d1 = c("a", "b", "c"), d2 = c("p", "q", "r", "s"),
                        stringsAsFactors = FALSE)
  dat <- strata[rep(seq_len(nrow(strata)), each = 40), ]
  dat$x <- rnorm(nrow(dat)); rownames(dat) <- NULL
  dat <- make_strata(dat, c("d1", "d2"))$data
  # a clear between-stratum signal keeps every null refit off the boundary
  u <- rnorm(max(dat$stratum), sd = 1.0)
  dat$y <- 0.5 * dat$x + u[dat$stratum] + rnorm(nrow(dat), sd = 1)

  ws <- collect_warnings(
    res <<- suppressMessages(pcv_importance(dat, "y", c("d1", "d2"), method = "shapley",
                                            bootstrap = TRUE, n_boot = 300)))
  expect_identical(res$n_boot_boundary, 0L)
  expect_false(any(grepl("boundary|conditional on a positive null variance",
                         ws, ignore.case = TRUE)))
})

# ---- Finding 2a: longitudinal per-time band obeys the majority-of-eligible rule ---

test_that("maihda_longitudinal_vpc_band applies both the floor and the majority rule", {
  band <- MAIHDA:::maihda_longitudinal_vpc_band
  set.seed(1)

  # 10 finite of 1000 clears the old bare `length(col) >= 10L` gate but is only 1% of
  # the eligible draws -- a biased handful. It must now yield NA, not a band.
  expect_identical(band(c(rnorm(10), rep(NA_real_, 990)), n_boot = 1000, conf_level = 0.95),
                   c(NA_real_, NA_real_))
  # A majority of finite draws yields a proper, ordered band.
  b <- band(c(rnorm(600), rep(NA_real_, 400)), n_boot = 1000, conf_level = 0.95)
  expect_true(all(is.finite(b)) && b[1] <= b[2])
  # The absolute floor of 10 still bites even when a majority is finite (9 of 10).
  expect_identical(band(rnorm(9), n_boot = 10, conf_level = 0.95), c(NA_real_, NA_real_))
  # Exactly at both thresholds (10 of 10) forms a band.
  expect_true(all(is.finite(band(rnorm(10), n_boot = 10, conf_level = 0.95))))
  # Majority boundary at n_boot = 200: 60 fails (< 100), 120 passes.
  expect_identical(band(c(rnorm(60), rep(NA_real_, 140)), n_boot = 200, conf_level = 0.95),
                   c(NA_real_, NA_real_))
  expect_true(all(is.finite(band(c(rnorm(120), rep(NA_real_, 80)), n_boot = 200, conf_level = 0.95))))
})

# ---- Finding 2b: crossed / contextual / longitudinal report n_boot_nonconverged ---

# Recursively test whether a key appears anywhere in a nested list (the field lives at
# different depths across the summary variants).
has_key <- function(x, key) {
  if (is.list(x)) {
    if (!is.null(names(x)) && key %in% names(x)) return(TRUE)
    return(any(vapply(x, has_key, logical(1), key = key)))
  }
  FALSE
}

test_that("crossed, contextual, and longitudinal bootstrap summaries disclose convergence", {
  skip_on_cran()
  set.seed(11)
  N <- 900
  ds <- data.frame(a = sample(c("m", "f"), N, TRUE),
                   b = sample(c("lo", "mid", "hi"), N, TRUE),
                   site = sample(paste0("S", 1:20), N, TRUE), x = rnorm(N))
  strat <- interaction(ds$a, ds$b, drop = TRUE)
  ds$y <- rnorm(nlevels(strat), sd = 1)[strat] +
          rnorm(20, sd = 0.8)[factor(ds$site)] + 0.4 * ds$x + rnorm(N, sd = 1.2)

  # standard VPC already disclosed this; it is the reference the others must match.
  m_std <- suppressMessages(fit_maihda(y ~ a + b + (1 | a:b), data = ds))
  s_std <- suppressWarnings(suppressMessages(summary(m_std, bootstrap = TRUE, n_boot = 15)))
  expect_true(has_key(s_std, "n_boot_nonconverged"))

  # contextual cross-classified (was absent before the fix)
  m_ctx <- suppressMessages(fit_maihda(y ~ x + (1 | a:b), data = ds, context = "site"))
  s_ctx <- suppressWarnings(suppressMessages(summary(m_ctx, bootstrap = TRUE, n_boot = 15)))
  expect_true(has_key(s_ctx, "n_boot_nonconverged"))

  # crossed-dimensions decomposition (was absent before the fix)
  a_cc <- suppressWarnings(suppressMessages(
    maihda(y ~ a + b + (1 | a:b), data = ds, decomposition = "crossed-dimensions")))
  cc_model <- NULL
  for (nm in names(a_cc)) {
    el <- a_cc[[nm]]
    if (inherits(el, "maihda_model") && !is.null(el$cc_info)) { cc_model <- el; break }
  }
  if (is.null(cc_model) && !is.null(a_cc$cc_info)) cc_model <- a_cc
  skip_if(is.null(cc_model), "could not locate the crossed-dimensions model")
  s_cc <- suppressWarnings(suppressMessages(summary(cc_model, bootstrap = TRUE, n_boot = 15)))
  expect_true(has_key(s_cc, "n_boot_nonconverged"))

  # longitudinal (was absent before the fix)
  set.seed(7)
  ni <- 160; waves <- 0:4
  ld <- data.frame(id = factor(rep(seq_len(ni), each = length(waves))),
                   wave = rep(waves, ni),
                   g = rep(sample(c("m", "f"), ni, TRUE), each = length(waves)),
                   e = rep(sample(c("lo", "hi"), ni, TRUE), each = length(waves)))
  ld$y <- 1 + 0.2 * ld$wave + rnorm(ni, sd = 0.6)[as.integer(ld$id)] + rnorm(nrow(ld), sd = 1)
  m_lng <- suppressMessages(fit_maihda(y ~ wave + (1 | g:e), data = ld, id = "id", time = "wave"))
  s_lng <- suppressWarnings(suppressMessages(summary(m_lng, bootstrap = TRUE, n_boot = 15)))
  expect_true(has_key(s_lng, "n_boot_nonconverged"))
})
