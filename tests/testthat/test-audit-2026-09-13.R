# Audit 2026-09-13: "A8. Monte Carlo correction is presented as exact
# finite-sample calibration" -- CONFIRMED, together with a floating-point defect
# in the interval that the same help page said agrees exactly with the p-value.
#
# (1) The wording. ?summary.maihda_model said the df_method = "bootstrap" p-value
# "is exact at any n_boot for which (n_boot + 1) alpha is a whole number", and
# ?maihda_proportional_odds_test that proportional-odds data "rejects at the
# nominal rate by construction". Both simulate from a fit whose parameters were
# estimated from the same data. The plus-one Monte Carlo p-value is exact only
# when the simulated reference IS the null distribution of the statistic, which
# a plug-in parametric bootstrap's is not (Dufour 2006, J. Econometrics 133,
# 443-477). Measured on the shipped functions, at replicate counts for which the
# old text promised exactness:
#   binomial MAIHDA, 4 strata of 120, stratum SD 0.5, no dimension effects:
#     P(p <= 0.05) = 0.137 (350 tests, n_boot = 19) and 0.141 (320, n_boot = 99)
#     -- the same 14% at both, against 36-40% for the Wald z; 23% on the singular
#     56% of fits and 3% on the rest; p-value distribution uniform at p = 1e-20;
#   proportional-odds test, 4 strata of 1,000, stratum SD 1.5, n_sim = 19:
#     rejection 0.067 at a nominal 5% and 0.155 at 10% over 600 datasets,
#     uniformity rejected at p = 3e-11 (the chi-squared reference it replaces
#     rejected 41% of the same correctly specified fits).
# The two help pages now describe parametric-bootstrap approximations with
# finite Monte Carlo resolution. A level is a long simulation, not a unit test,
# so the help-page block below pins the wording and the audit log keeps the runs.
#
# (2) The interval. The critical value's rank was
# ceiling((1 - conf_level) * (B + 1)) - 1. 1 - 0.95 is stored as
# 0.05000000000000004 and 1 - 0.90 as 0.09999999999999998, so wherever
# (B + 1) * alpha is a whole number the product landed just above it at some
# levels and at or below it at others. The interval therefore excluded zero when
# p <= alpha at the first kind and only when p < alpha at the second. Live, with
# n_boot = 199, a p-value of exactly 0.05 excluded zero at conf_level = 0.95
# while one of exactly 0.10 did not at 0.90. The fix applies p <= alpha, the
# usual Monte Carlo test rule, at every level through
# maihda_boot_critical_rank(). Which levels that moves is decided by storage and
# has to be measured, not assumed: 0.50, 0.60, 0.66, 0.68, 0.75, 0.80 and 0.90
# move by one order statistic on whole-number grids, while 0.70, 0.85, 0.95,
# 0.975, 0.98, 0.99, 0.995 and 0.999 are bit-identical.

test_that("the bootstrap critical rank excludes zero exactly when p <= alpha", {
  # Every attainable count c = #{|t*| >= |t|}: zero falls outside the interval
  # when the observed statistic beats the r-th largest draw, i.e. when
  # c <= r - 1, and that must be exactly p = (1 + c) / (B + 1) <= alpha with
  # alpha typed as the literal the user means by the level.
  # 0.50 and 0.75 are in the set because their complements are EXACTLY
  # representable, which is one of the two ways the old index took the strict
  # rule; 0.90 and 0.80 are the other (stored a hair below).
  lev <- c(0.50, 0.75, 0.80, 0.90, 0.95, 0.99)
  alpha <- c(0.50, 0.25, 0.20, 0.10, 0.05, 0.01)
  for (i in seq_along(lev)) {
    for (B in c(10:60, 99L, 199L, 999L, 1999L)) {
      r <- maihda_boot_critical_rank(B, lev[i])
      cc <- 0:B
      expect_identical(cc <= r - 1L, (1 + cc) / (B + 1) <= alpha[i],
                       info = sprintf("conf_level %.2f, B = %d", lev[i], B))
    }
  }
  # The fixed points the help page and the code comments rely on.
  expect_identical(maihda_boot_critical_rank(199L, 0.95), 10L)
  expect_identical(maihda_boot_critical_rank(999L, 0.95), 50L)
  expect_identical(maihda_boot_critical_rank(199L, 0.90), 20L)
  expect_identical(maihda_boot_critical_rank(40L, 0.95), 2L)
  expect_identical(maihda_boot_critical_rank(19L, 0.95), 1L)
  # With 18 draws no attainable p reaches 0.05, so the interval is unbounded.
  expect_identical(maihda_boot_critical_rank(18L, 0.95), 0L)
})

test_that("the default levels keep the pre-fix critical value", {
  # Negative control. At conf_level 0.95 and 0.99 the old index already applied
  # p <= alpha -- the representation error happened to land above the whole
  # number -- so nothing moves there, for any n. Where a level does move it moves
  # by exactly one, and only where (B + 1) * alpha is whole.
  old <- function(B, conf_level) as.integer(ceiling((1 - conf_level) * (B + 1)) - 1)
  B <- 10:5000
  for (cl in c(0.95, 0.99)) {
    expect_identical(maihda_boot_critical_rank(B, cl), old(B, cl),
                     info = sprintf("conf_level %.2f", cl))
  }
  moved <- maihda_boot_critical_rank(B, 0.90) - old(B, 0.90)
  expect_identical(which(moved != 0L), which((B + 1L) %% 10L == 0L))
  expect_true(all(moved[moved != 0L] == 1L))

  # Which levels move is an accident of storage, not a property of 0.90 and
  # 0.80: a complement stored at or below its decimal value took the strict rule
  # and moves, one stored above it did not. Pinned because the first write-up of
  # this fix claimed the change was confined to 0.90 and 0.80, and it is not.
  moves <- function(cl) any(maihda_boot_critical_rank(B, cl) != old(B, cl))
  expect_true(all(vapply(c(0.50, 0.60, 0.75, 0.80, 0.90), moves, logical(1))))
  expect_false(any(vapply(c(0.70, 0.85, 0.95, 0.975, 0.99), moves, logical(1))))
})

# A Gaussian fixture whose one tested coefficient sits between the bootstrap
# reference's 90th and 95th percentiles, so two bootstrap seeds land its p-value
# EXACTLY on a level -- the only case in which the two rules disagree.
a8_tie_data <- function(seed = 6, beta = 0.45) {
  set.seed(seed)
  g <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:2))
  d <- g[rep(seq_len(8), each = 40), , drop = FALSE]
  d$stratum <- factor(rep(seq_len(8), each = 40))
  u <- stats::rnorm(8, 0, 0.5)
  d$y <- u[d$stratum] + stats::rnorm(nrow(d)) + beta * (d$d1 == "2")
  rownames(d) <- NULL
  d
}

test_that("a bootstrap p-value on the level excludes zero at 0.90 as at 0.95", {
  skip_on_cran()
  d <- a8_tie_data()
  fit <- fit_maihda(y ~ d1 + (1 | stratum), data = d)
  excluded <- function(fe) fe$lower[2] > 0 || fe$upper[2] < 0

  # conf_level = 0.90 through the public path: 19 of the 199 draws reach the
  # observed statistic, so p = 20 / 200 = 0.10 exactly. Before the fix this
  # interval contained zero, [-0.0049, 1.4164]; the p-value was the same.
  set.seed(15)
  fe90 <- suppressWarnings(summary(fit, df_method = "bootstrap", n_boot = 199,
                                   conf_level = 0.90))$fixed_effects
  expect_identical(fe90$p_value[2], 0.1,
                   info = "fixture premise: this seed must land p exactly on 0.10")
  expect_true(excluded(fe90))

  # conf_level = 0.95: p = 10 / 200 = 0.05 exactly already excluded zero, and
  # the fix leaves this level bit-identical.
  set.seed(17)
  fe95 <- suppressWarnings(maihda_bootstrap_fixef(fit$model, n_boot = 199L,
                                                  conf_level = 0.95))
  expect_identical(fe95$p_value[2], 0.05,
                   info = "fixture premise: this seed must land p exactly on 0.05")
  expect_true(excluded(fe95))
})

test_that("the bootstrap help pages call both tests approximations, not exact", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree (test_local / load_all);
  # R CMD check copies only tests/, so there is nothing to read there.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  # Whitespace-normalised, because Rd re-wraps and a phrase can straddle a line.
  rd <- function(f) {
    gsub("[[:space:]]+", " ",
         paste(readLines(file.path(man, f), warn = FALSE), collapse = " "))
  }
  # Spelling-agnostic: remove the legitimate uses of "exact", then no form of the
  # word may remain. Banning only the old sentences would let a reworded relapse
  # through.
  no_exact_left <- function(txt, allowed) {
    for (a in allowed) txt <- gsub(a, "", txt, fixed = TRUE)
    !grepl("exact", txt, ignore.case = TRUE)
  }

  sm <- rd("summary.maihda_model.Rd")
  expect_true(grepl("\\section{Fixed-effect reference distribution}", sm, fixed = TRUE))
  sec <- sub(".*\\\\section\\{Fixed-effect reference distribution\\}\\{", "", sm)
  sec <- sub("\\\\section\\{Two VPCs.*", "", sec)
  expect_lt(nchar(sec), nchar(sm))            # the section was actually isolated
  expect_true(grepl("approximation, not an exact test", sec, fixed = TRUE))
  expect_true(grepl("sets the Monte Carlo resolution, not that approximation", sec,
                    fixed = TRUE))
  expect_true(grepl("No \\code{n_boot} makes the test exact", sec, fixed = TRUE))
  expect_true(grepl("exactly when the p-value is at most \\code{1 - conf_level}",
                    sec, fixed = TRUE))
  expect_true(grepl("algebraic, not a coverage guarantee", sec, fixed = TRUE))
  expect_true(no_exact_left(sec, c("approximation, not an exact test",
                                   "No \\code{n_boot} makes the test exact",
                                   "exactly when the p-value is at most",
                                   "estimated at exactly zero")))
  expect_false(grepl("tighten", sec, ignore.case = TRUE))
  expect_false(grepl("true null", sec, fixed = TRUE))
  expect_false(grepl("significant", sec, fixed = TRUE))

  po <- rd("maihda_proportional_odds_test.Rd")
  expect_true(grepl("parametric-bootstrap approximation, not an exact test", po,
                    fixed = TRUE))
  expect_true(grepl("Monte Carlo resolution", po, fixed = TRUE))
  # (The symmetric-threshold sentence this list also allowed was removed as false
  # by the third 2026-09-13 pass; see test-audit-2026-09-13c.R.)
  expect_true(no_exact_left(po, "approximation, not an exact test"))
  expect_false(grepl("by construction", po, fixed = TRUE))
  expect_false(grepl("nominal rate", po, fixed = TRUE))
  expect_false(grepl("removes that confounding", po, fixed = TRUE))

  bf <- rd("maihda_bootstrap_fixef.Rd")
  expect_false(grepl("true null", bf, fixed = TRUE))
  expect_true(grepl("parametric-bootstrap approximation", bf, fixed = TRUE))
})
