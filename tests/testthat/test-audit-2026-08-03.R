# Audit pass 2026-08-03: model-adequacy checks.
#
# The package checked only COMPUTATIONAL adequacy (singular / non-converged fits);
# it never checked whether the chosen LIKELIHOOD was adequate, so an overdispersed
# or zero-inflated count, a non-normal stratum random effect, autocorrelated
# longitudinal residuals, or a violated proportional-odds assumption all produced a
# confidently-reported VPC/PCV with no caveat. maihda_fit_diagnostics() now runs a
# suite of adequacy probes that surface in the print()/summary() "Fit diagnostics"
# block. These tests pin each pure probe on synthetic inputs (deterministic, no
# reliance on an optimiser stopping in a particular place) plus a few live-fit
# integration checks, including the false-positive control that a well-specified
# model stays silent.

# ---- overdispersion statistic (pure) ----------------------------------------

test_that("maihda_overdispersion_stat is the Pearson chisq over resid df", {
  st <- maihda_overdispersion_stat(rep(2, 100), 50)
  expect_equal(st$chisq, 400)          # sum(2^2) over 100 = 400
  expect_equal(st$rdf, 50)
  expect_equal(st$ratio, 8)            # 400 / 50
  expect_lt(st$p, 1e-6)                # chisq 400 on 50 df

  # a well-specified count has chisq ~ rdf -> ratio ~ 1
  clean <- maihda_overdispersion_stat(rep(1, 50), 50)
  expect_equal(clean$ratio, 1)
  expect_gt(clean$p, 0.1)

  # guards: no residuals, non-positive df
  expect_null(maihda_overdispersion_stat(numeric(0), 50))
  expect_null(maihda_overdispersion_stat(rep(1, 10), 0))
  expect_null(maihda_overdispersion_stat(rep(1, 10), NA_real_))
})

# ---- zero-inflation statistic (pure) ----------------------------------------

test_that("maihda_zeroinflation_stat is expected/observed zeros", {
  st <- maihda_zeroinflation_stat(observed_zeros = 100, expected_zeros = 25, n = 500)
  expect_equal(st$observed, 100)
  expect_equal(st$expected, 25)
  expect_equal(st$ratio, 0.25)         # model predicts a quarter of the zeros seen

  # more expected than observed -> ratio > 1 (no zero inflation)
  expect_equal(maihda_zeroinflation_stat(20, 30, 200)$ratio, 1.5)

  # guard: no observed zeros -> nothing to judge
  expect_null(maihda_zeroinflation_stat(0, 10, 100))
})

# ---- random-effect normality statistic (pure) -------------------------------

test_that("maihda_re_normality_stat flags a skewed BLUP vector and passes a normal one", {
  # deterministic exponential-quantile grid: strongly right-skewed
  skewed <- stats::qexp(stats::ppoints(60))
  st_s <- maihda_re_normality_stat(skewed)
  expect_gt(st_s$skew, 0.75)
  expect_lt(st_s$shapiro_p, 0.01)

  # deterministic normal-quantile grid: symmetric, Shapiro does not reject
  normalish <- stats::qnorm(stats::ppoints(60))
  st_n <- maihda_re_normality_stat(normalish)
  expect_lt(abs(st_n$skew), 0.2)
  expect_gt(st_n$shapiro_p, 0.1)

  # guards: too few values, or no variance
  expect_null(maihda_re_normality_stat(c(1, 2, 3)))
  expect_null(maihda_re_normality_stat(rep(4, 40)))
})

test_that("maihda_re_standardize rescales by the marginal BLUP sd and drops degenerate groups", {
  z <- maihda_re_standardize(c(1, 2), c(0.5, 0.75), tau2 = 1)
  expect_equal(z, c(1 / sqrt(0.5), 2 / sqrt(0.25)))

  # a group whose condvar reaches tau2 has no marginal BLUP variance: drop, not divide
  z2 <- maihda_re_standardize(c(1, 2, 3), c(0.5, 1, 2), tau2 = 1)
  expect_equal(z2, 1 / sqrt(0.5))

  # guards: boundary tau2, mismatched lengths
  expect_identical(maihda_re_standardize(c(1, 2), c(0.1, 0.1), tau2 = 0), numeric(0))
  expect_identical(maihda_re_standardize(c(1, 2), c(0.1), tau2 = 1), numeric(0))
})

test_that("unequal-information scale mixture no longer masquerades as RE non-normality", {
  # Audit regression: Var(u_hat_j) = tau2 - v_j depends on cell size, so raw BLUPs
  # from unequal strata are a scale mixture of normals -- leptokurtic with exactly
  # Gaussian true effects. The shape gate must trip on the raw mixture and must NOT
  # trip once each half is standardized by its model-implied marginal sd.
  th <- maihda_adequacy_thresholds()
  base <- stats::qnorm(stats::ppoints(40))
  raw <- c(base * 1, base * 4)          # sd 1 (heavily shrunk) vs sd 4 (data-rich)
  raw_st <- maihda_re_normality_stat(raw)
  expect_gt(raw_st$excess_kurtosis, th$re_excess_kurtosis)
  expect_lt(raw_st$shapiro_p, th$re_shapiro_p)

  # tau2 = 16 with v = 15 (tiny cells) / 0 (huge cells) reproduces those sds exactly
  z <- maihda_re_standardize(raw, c(rep(15, 40), rep(0, 40)), tau2 = 16)
  std_st <- maihda_re_normality_stat(z)
  expect_lt(abs(std_st$excess_kurtosis), th$re_excess_kurtosis)
  expect_gt(std_st$shapiro_p, th$re_shapiro_p)
})

# ---- residual autocorrelation statistic (pure) ------------------------------

test_that("maihda_resid_autocorr_stat computes pooled lag-1 correlation", {
  # perfectly alternating residuals within one unit -> lag-1 correlation -1
  alt <- maihda_resid_autocorr_stat(rep(c(1, -1), 5), rep(1L, 10), 1:10)
  expect_equal(alt$acf1, -1)
  expect_equal(alt$n_pairs, 9)
  expect_equal(alt$gap, 1)

  # two units, cur = prev + 1 within each -> perfect positive correlation, pooled
  st <- maihda_resid_autocorr_stat(
    resid = c(1, 2, 3, 1, 2, 3),
    id = c(1L, 1L, 1L, 2L, 2L, 2L),
    time = c(1, 2, 3, 1, 2, 3))
  expect_equal(st$acf1, 1)
  expect_equal(st$n_pairs, 4)          # 2 consecutive pairs per unit
})

test_that("maihda_resid_autocorr_stat pairs only the modal time gap", {
  # gaps 1, 1, 2 within one unit -> modal gap 1, the wide pair is dropped
  st <- maihda_resid_autocorr_stat(c(1, 2, 1, 9), rep(1L, 4), c(0, 1, 2, 4))
  expect_equal(st$gap, 1)
  expect_equal(st$n_pairs, 2)          # (t0,t1) and (t1,t2) only

  # tie between gap 1 and gap 2 -> smallest wins
  st2 <- maihda_resid_autocorr_stat(
    c(1, 2, 9, 5, 3, 7), c(1L, 1L, 1L, 2L, 2L, 2L),
    c(0, 1, 3, 0, 1, 3))               # gaps 1,2 per unit -> 2 of each
  expect_equal(st2$gap, 1)
  expect_equal(st2$n_pairs, 2)
})

test_that("maihda_resid_autocorr_stat orders by time within unit and skips singletons", {
  ordered <- maihda_resid_autocorr_stat(c(1, 5, 2, 6), rep(1L, 4), c(1, 2, 3, 4))
  # same (time, residual) pairs supplied out of order must give the same result
  shuffled <- maihda_resid_autocorr_stat(c(2, 1, 6, 5), rep(1L, 4), c(3, 1, 4, 2))
  expect_equal(shuffled$acf1, ordered$acf1)

  # a single-observation unit contributes no consecutive pair
  with_singleton <- maihda_resid_autocorr_stat(
    c(1, 2, 3, 99), c(1L, 1L, 1L, 2L), c(1, 2, 3, 1))
  expect_equal(with_singleton$n_pairs, 2)

  expect_null(maihda_resid_autocorr_stat(c(1, 2), c(1L, 1L), c(1, 2)))  # < 3 obs
})

# ---- thresholds are pinned --------------------------------------------------

test_that("adequacy thresholds are the documented conservative values", {
  th <- maihda_adequacy_thresholds()
  expect_equal(th$overdispersion_ratio, 1.5)
  expect_equal(th$zeroinflation_ratio, 0.9)
  expect_equal(th$re_min_levels, 20L)
  expect_equal(th$autocorr_min, 0.3)
})

# ---- surfacing: one caveat line per FLAGGED check (pure) ---------------------

test_that("maihda_format_adequacy renders exactly the flagged checks", {
  expect_identical(maihda_format_adequacy(NULL), character(0))

  # nothing flagged -> no lines
  unflagged <- list(overdispersion = list(ratio = 1.1, p = 0.4,
                                          family = "poisson", flag = FALSE))
  expect_identical(maihda_format_adequacy(unflagged), character(0))

  full <- list(
    overdispersion = list(ratio = 6.8, p = 1e-9, family = "poisson", flag = TRUE),
    zeroinflation = list(observed = 583, expected = 157, ratio = 0.27, flag = TRUE),
    re_normality = list(group = "stratum", term = "(Intercept)", skew = 1.1,
                        excess_kurtosis = 0.6, shapiro_p = 0.002, flag = TRUE),
    autocorrelation = list(acf1 = 0.51, n_pairs = 2100, gap = 1, flag = TRUE),
    # descriptive only since 2026-09-01: stored, never rendered as a caveat
    marginal_po_proxy = list(lrt = 40, df = 1, n_terms = 1))
  lines <- maihda_format_adequacy(full)
  expect_match(lines, "Overdispersion", all = FALSE)
  expect_match(lines, "negbinomial", all = FALSE)          # remedy for a poisson fit
  expect_match(lines, "Zero inflation", all = FALSE)
  expect_match(lines, "Random-effect distribution", all = FALSE)
  expect_match(lines, "standardized BLUPs", all = FALSE)
  expect_match(lines, "Residual autocorrelation", all = FALSE)
  expect_match(lines, "at time gap 1", all = FALSE)
  # the fixed-only proportional-odds proxy is confounded with stratum
  # heterogeneity and must NEVER surface as an automatic caveat
  expect_false(any(grepl("Proportional odds", lines, fixed = TRUE)))

  # the negbinomial overdispersion message points beyond the NB, not back to it
  nb <- maihda_format_adequacy(list(
    overdispersion = list(ratio = 3, p = 1e-4, family = "negbinomial", flag = TRUE)))
  expect_match(nb, "beyond the negative binomial", all = FALSE)
  expect_false(any(grepl("consider family", nb)))

  # a flagged random SLOPE names the term instead of claiming intercepts
  slope <- maihda_format_adequacy(list(
    re_normality = list(group = "pid", term = "t", skew = 1.2,
                        excess_kurtosis = 2.0, shapiro_p = 0.001, flag = TRUE)))
  expect_match(slope, "random 't' slopes", all = FALSE)
  expect_false(any(grepl("intercepts", slope)))
})

# ---- integration: the diagnostics slot carries adequacy ---------------------

test_that("a well-specified Gaussian fit raises no adequacy caveat", {
  skip_on_cran()
  set.seed(303)
  strata <- expand.grid(a = factor(1:6), b = factor(1:5))
  d <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i)
    data.frame(a = strata$a[i], b = strata$b[i], x = rnorm(40))))
  idx <- interaction(d$a, d$b, drop = TRUE)
  d$y <- 1 + 0.5 * d$x + rnorm(nlevels(idx), 0, 0.6)[as.integer(idx)] + rnorm(nrow(d))
  m <- suppressWarnings(fit_maihda(y ~ x + (1 | a:b), data = d))
  # the field exists; nothing is flagged
  flagged <- Filter(function(x) isTRUE(x$flag), m$diagnostics$adequacy)
  expect_length(flagged, 0)
  expect_identical(maihda_format_adequacy(m$diagnostics$adequacy), character(0))
})

test_that("an overdispersed Poisson fit flags overdispersion and zero inflation", {
  skip_on_cran()
  set.seed(304)
  strata <- expand.grid(a = factor(1:4), b = factor(1:3), c = factor(1:2))
  d <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i)
    data.frame(a = strata$a[i], b = strata$b[i], c = strata$c[i], x = rnorm(60))))
  idx <- interaction(d$a, d$b, d$c, drop = TRUE)
  mu <- exp(1 + 0.3 * d$x + rnorm(nlevels(idx), 0, 0.4)[as.integer(idx)])
  d$y <- rpois(nrow(d), mu * rgamma(nrow(d), 0.5, 0.5))   # gamma frailty -> overdispersed
  m <- suppressWarnings(fit_maihda(y ~ x + (1 | a:b:c), data = d, family = "poisson"))
  expect_true(isTRUE(m$diagnostics$adequacy$overdispersion$flag))
  expect_gt(m$diagnostics$adequacy$overdispersion$ratio, 1.5)
  expect_true(isTRUE(m$diagnostics$adequacy$zeroinflation$flag))
  expect_match(maihda_format_adequacy(m$diagnostics$adequacy),
               "Overdispersion", all = FALSE)
})

test_that("a longitudinal fit with AR(1) residuals flags autocorrelation", {
  skip_on_cran()
  set.seed(305)
  np <- 150
  base <- expand.grid(pid = 1:np, t = 0:14)
  strat <- data.frame(pid = 1:np, a = factor(sample(1:4, np, TRUE)),
                      b = factor(sample(1:3, np, TRUE)))
  d <- merge(base, strat, by = "pid")
  d <- d[order(d$pid, d$t), ]
  si <- interaction(d$a, d$b, drop = TRUE)
  e <- numeric(nrow(d)); rho <- 0.85
  for (i in seq_len(nrow(d))) {
    e[i] <- if (d$t[i] == 0) rnorm(1) else rho * e[i - 1] + rnorm(1, 0, sqrt(1 - rho^2))
  }
  d$y <- 2 + 0.2 * d$t + rnorm(nlevels(si), 0, 0.7)[as.integer(si)] +
    rnorm(np, 0, 0.5)[d$pid] + rnorm(np, 0, 0.05)[d$pid] * d$t + e
  m <- suppressWarnings(fit_maihda(y ~ (1 | a:b), data = d, id = "pid", time = "t"))
  expect_true(isTRUE(m$diagnostics$adequacy$autocorrelation$flag))
  expect_gt(m$diagnostics$adequacy$autocorrelation$acf1, 0.3)
})

test_that("a CORRECT short panel is not flagged despite the negative detrending artifact", {
  # Audit regression: conditional residuals subtract the estimated per-id
  # intercept+slope, which mechanically induces NEGATIVE lag-1 correlation on a
  # correctly specified model -- around -0.4 at 4 waves with well-resolved random
  # effects. The one-sided flag must not read that artifact as misspecification.
  skip_on_cran()
  set.seed(306)
  np <- 150
  base <- expand.grid(pid = 1:np, t = 0:3)                 # 4 waves
  strat <- data.frame(pid = 1:np, a = factor(sample(1:4, np, TRUE)),
                      b = factor(sample(1:3, np, TRUE)))
  d <- merge(base, strat, by = "pid")
  d <- d[order(d$pid, d$t), ]
  si <- interaction(d$a, d$b, drop = TRUE)
  d$y <- 2 + 0.3 * d$t +
    rnorm(nlevels(si), 0, 0.5)[as.integer(si)] +
    0.10 * d$t * rnorm(nlevels(si))[as.integer(si)] +
    rnorm(np, 0, 1.2)[d$pid] +                             # light shrinkage
    rnorm(np, 0, 0.4)[d$pid] * d$t +
    rnorm(nrow(d), 0, 0.7)                                 # INDEPENDENT errors
  m <- suppressWarnings(fit_maihda(y ~ (1 | a:b), data = d, id = "pid", time = "t"))
  ac <- m$diagnostics$adequacy$autocorrelation
  expect_false(is.null(ac))
  expect_lt(ac$acf1, -0.3)             # the artifact is there, clearly
  expect_false(isTRUE(ac$flag))        # and it must NOT flag
})

test_that("a correct model with strongly unequal strata is not RE-flagged", {
  # Audit regression: Gaussian true effects, cells of 5 vs 800. The raw-BLUP
  # scale-mixture artifact is demonstrably present in this fixture; the
  # standardized pipeline must not read it as non-normality.
  skip_on_cran()
  set.seed(302)
  sizes <- c(rep(5L, 20), rep(800L, 20))
  ab <- expand.grid(a = factor(1:8), b = factor(1:5))       # 40 strata
  u <- rnorm(40, 0, 0.15)                                   # GAUSSIAN truth
  d <- do.call(rbind, lapply(1:40, function(j)
    data.frame(a = ab$a[j], b = ab$b[j], x = rnorm(sizes[j]), j = j)))
  d$y <- 1 + 0.5 * d$x + u[d$j] + rnorm(nrow(d), 0, 1)
  m <- suppressWarnings(fit_maihda(y ~ x + (1 | a:b), data = d[, c("a", "b", "x", "y")]))
  # the raw-BLUP artifact would have flagged this fit under BOTH gate conditions
  th <- maihda_adequacy_thresholds()
  raw <- maihda_re_normality_stat(lme4::ranef(m$model)$stratum[["(Intercept)"]])
  expect_gt(raw$excess_kurtosis, th$re_excess_kurtosis)
  expect_lt(raw$shapiro_p, th$re_shapiro_p)
  # the standardized pipeline does not read it as misspecification
  expect_false(isTRUE(m$diagnostics$adequacy$re_normality$flag))
})

test_that("the standardized RE check still detects a genuinely skewed random effect", {
  skip_on_cran()
  set.seed(308)
  sizes <- c(rep(5L, 20), rep(800L, 20))                    # same unequal design
  ab <- expand.grid(a = factor(1:8), b = factor(1:5))
  u <- (rexp(40) - 1) * 0.6                                 # strongly skewed truth
  d <- do.call(rbind, lapply(1:40, function(j)
    data.frame(a = ab$a[j], b = ab$b[j], x = rnorm(sizes[j]), j = j)))
  d$y <- 1 + 0.5 * d$x + u[d$j] + rnorm(nrow(d), 0, 1)
  m <- suppressWarnings(fit_maihda(y ~ x + (1 | a:b), data = d[, c("a", "b", "x", "y")]))
  re <- m$diagnostics$adequacy$re_normality
  expect_true(isTRUE(re$flag))
  expect_identical(re$term, "(Intercept)")
})

# ---- ordinal proportional-odds (live clmm) ----------------------------------

make_po_data <- function(seed, violate) {
  set.seed(seed)
  strata <- expand.grid(a = factor(1:4), b = factor(1:3))
  d <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i)
    data.frame(a = strata$a[i], b = strata$b[i], x = rnorm(120))))
  idx <- interaction(d$a, d$b, drop = TRUE)
  u <- rnorm(nlevels(idx), 0, 0.5)[as.integer(idx)]
  if (violate) {
    lat <- u + rlogis(nrow(d))
    d$y <- factor(ifelse(lat < -1 - 1.5 * d$x, 1L,
                         ifelse(lat < 1 - 0.1 * d$x, 2L, 3L)), 1:3, ordered = TRUE)
  } else {
    lat <- 0.8 * d$x + u + rlogis(nrow(d))
    d$y <- factor(ifelse(lat < -1, 1L, ifelse(lat < 1, 2L, 3L)), 1:3, ordered = TRUE)
  }
  d
}

# REWRITTEN 2026-09-01: this test used to assert that the fixed-only screen
# FLAGGED the violated fit via a raw chi-squared p-value. That screen was
# confounded with stratum heterogeneity (see test-audit-2026-09-01.R), so the
# flag was removed. What survives is the separation itself, now judged by the
# calibrated bootstrap test rather than by an automatic caveat.
test_that("the calibrated proportional-odds test separates violated from proportional", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  m_bad <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = make_po_data(311, TRUE), family = "ordinal")))
  t_bad <- maihda_proportional_odds_test(m_bad, n_sim = 99, seed = 311)
  expect_lt(t_bad$p_value, 0.05)

  m_ok <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = make_po_data(312, FALSE), family = "ordinal")))
  t_ok <- maihda_proportional_odds_test(m_ok, n_sim = 99, seed = 312)
  expect_gt(t_ok$p_value, 0.05)

  # neither fit may raise an automatic proportional-odds caveat any more
  for (m in list(m_bad, m_ok)) {
    expect_null(m$diagnostics$adequacy$proportional_odds)
    expect_false(any(grepl("Proportional odds",
                           maihda_format_adequacy(m$diagnostics$adequacy),
                           fixed = TRUE)))
  }
})

test_that("a null ordinal model has no covariate slopes to test", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  m0 <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ (1 | a:b), data = make_po_data(313, TRUE), family = "ordinal")))
  expect_null(maihda_ordinal_po_stat(m0$model))
})
