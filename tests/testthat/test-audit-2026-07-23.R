# Regression tests for the 2026-07-23 audit finding: reordered-but-identical
# analytic samples were incorrectly rejected by the model-comparison guards.
#
#   [P2] calculate_pcv() required identical ROW ORDER (validate_pcv_models) and
#        compared reorder-dependent stratum labels, while the response / prior-weight
#        / sampling-weight fingerprints pasted values in row order despite a comment
#        claiming order-invariance. The same order-sensitivity produced false
#        "analytic sample differs" warnings in compare_maihda() and false delta
#        issues in maihda_ic(). Two fits to the SAME observations (same persistent
#        row ids) in a different order describe the same analytic sample -- their
#        between-stratum variances, and hence the PCV, are order-invariant -- so the
#        comparison now aligns on the row ids (matching response, stratum PARTITION
#        and weights after alignment) and accepts them, while still rejecting
#        genuinely different data.

# ---- unit: order-invariant id alignment ---------------------------------------

test_that("maihda_order_by_ids canonicalises by row id and falls back safely", {
  f <- MAIHDA:::maihda_order_by_ids
  # Values travel with their ids into a canonical (id-sorted) order.
  expect_identical(f(c(10, 20, 30), c("c", "a", "b")), c(20, 30, 10))
  # The SAME (id -> value) map supplied in a different order yields the SAME
  # canonical sequence -- the property the fingerprints rely on.
  expect_identical(f(c(10, 20, 30), c("c", "a", "b")),
                   f(c(20, 30, 10), c("a", "b", "c")))
  # Unusable ids (NULL / length mismatch / NA / duplicated) leave the order alone,
  # i.e. degrade to the conservative order-sensitive behaviour.
  expect_identical(f(c(1, 2, 3), NULL), c(1, 2, 3))
  expect_identical(f(c(1, 2, 3), c("a", "b")), c(1, 2, 3))
  expect_identical(f(c(1, 2, 3), c("a", NA, "b")), c(1, 2, 3))
  expect_identical(f(c(1, 2, 3), c("a", "a", "b")), c(1, 2, 3))
})

test_that("maihda_order_by_ids reorders MATRIX rows by id, preserving dimensions", {
  f <- MAIHDA:::maihda_order_by_ids
  # An aggregated-binomial cbind(success, failure) response reaches the helper as an
  # n x 2 matrix. It must key on the number of ROWS (n), not length() (== 2n), which
  # previously failed the `length(ids) == length(values)` guard and skipped alignment.
  m <- matrix(c(10, 20, 30, 1, 2, 3), ncol = 2)  # 3 rows x 2 cols
  out <- f(m, c("c", "a", "b"))
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(3L, 2L))
  # Whole rows travel with their ids into id-sorted order; columns stay put.
  expect_identical(out, matrix(c(20, 30, 10, 2, 3, 1), ncol = 2))
  # The SAME (id -> row) map supplied in a different row order yields the SAME
  # canonical matrix -- the order-invariance the matrix fingerprint now relies on.
  m2 <- matrix(c(20, 30, 10, 2, 3, 1), ncol = 2)
  expect_identical(f(m, c("c", "a", "b")), f(m2, c("a", "b", "c")))
})

# ---- unit: stratum partition key (label- and order-independent) ---------------

test_that("maihda_stratum_partition_key ignores labelling and row order", {
  f <- MAIHDA:::maihda_stratum_partition_key
  # Same partition {r1,r2 | r3 | r4}, different LABELS and different ROW ORDER.
  k1 <- f(c("A", "A", "B", "C"), c("r1", "r2", "r3", "r4"))
  k2 <- f(c("z", "y", "x", "x"), c("r4", "r3", "r2", "r1"))
  expect_identical(k1, c(1L, 1L, 2L, 3L))
  expect_identical(k1, k2)
  # A genuinely different partition {r1,r3 | r2 | r4} is NOT equal.
  k3 <- f(c("A", "B", "A", "C"), c("r1", "r2", "r3", "r4"))
  expect_false(identical(k1, k3))
})

# ---- integration: reordered-but-identical fits are accepted --------------------

maihda_build_reorder_fixture <- function() {
  set.seed(42)
  N  <- 160
  d1 <- factor(sample(c("m", "f"), N, replace = TRUE))
  d2 <- factor(sample(c("lo", "mid", "hi"), N, replace = TRUE))
  sid <- droplevels(interaction(d1, d2))
  u  <- rnorm(nlevels(sid), 0, 0.6)[as.integer(sid)]
  y  <- 1 + 0.5 * (d1 == "m") + 0.3 * (d2 == "hi") + u + rnorm(N)
  dat <- data.frame(y = y, d1 = d1, d2 = d2)
  rownames(dat) <- paste0("obs", seq_len(N))
  dat
}

test_that("PCV/compare/IC accept a reordered-but-identical fit and agree numerically", {
  skip_on_cran()
  dat <- maihda_build_reorder_fixture()
  null_fit <- suppressMessages(suppressWarnings(fit_maihda(y ~ 1 + (1 | d1:d2), data = dat)))
  adj_fit  <- suppressMessages(suppressWarnings(fit_maihda(y ~ d1 + d2 + (1 | d1:d2), data = dat)))

  set.seed(7)
  shuf     <- dat[sample(nrow(dat)), ]
  adj_shuf <- suppressMessages(suppressWarnings(fit_maihda(y ~ d1 + d2 + (1 | d1:d2), data = shuf)))

  # Row ids are the same SET but a different ORDER.
  expect_true(setequal(MAIHDA:::maihda_wrapper_row_ids(adj_fit),
                       MAIHDA:::maihda_wrapper_row_ids(adj_shuf)))
  expect_false(identical(MAIHDA:::maihda_wrapper_row_ids(adj_fit),
                         MAIHDA:::maihda_wrapper_row_ids(adj_shuf)))

  # The response fingerprint is now genuinely order-invariant (its comment claim).
  expect_identical(MAIHDA:::maihda_wrapper_response_fingerprint(adj_fit),
                   MAIHDA:::maihda_wrapper_response_fingerprint(adj_shuf))

  # calculate_pcv no longer aborts, and returns the same PCV as the in-order fit.
  base  <- calculate_pcv(null_fit, adj_fit)
  reord <- calculate_pcv(null_fit, adj_shuf)
  expect_equal(reord$pcv, base$pcv, tolerance = 1e-4)

  # compare_maihda emits NO false "analytic sample" warning.
  w <- character(0)
  withCallingHandlers(
    compare_maihda(null_fit, adj_shuf, ic = FALSE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("analytic sample", w, fixed = TRUE)))

  # maihda_ic reports no false delta issue.
  expect_length(MAIHDA:::maihda_ic_delta_issues(list(null_fit, adj_shuf)), 0L)
})

# ---- integration: genuinely different samples are STILL rejected ---------------

test_that("PCV still rejects genuinely different analytic samples after the fix", {
  skip_on_cran()
  mk <- function(seed) {
    set.seed(seed)
    d <- data.frame(stratum = factor(rep(seq_len(8), each = 20)), x = rnorm(160))
    d$y <- 1 + 0.3 * d$x + rnorm(8, sd = 0.7)[d$stratum] + rnorm(160, sd = 0.4)
    suppressMessages(suppressWarnings(fit_maihda(y ~ x + (1 | stratum), data = d)))
  }
  # Same n, strata and default 1:160 names, but different OUTCOME values.
  expect_error(calculate_pcv(mk(101), mk(202)), "outcome values differ")

  # Same rows and outcome, but a different stratum PARTITION.
  set.seed(5)
  n  <- 80
  d1 <- data.frame(stratum = factor(rep(seq_len(4), each = n / 4)), x = rnorm(n))
  d1$y <- 2 + 0.4 * d1$x + rnorm(4, sd = 0.8)[d1$stratum] + rnorm(n, sd = 0.2)
  d2 <- d1
  d2$stratum <- factor(rep(seq_len(4), times = n / 4))
  m1 <- suppressMessages(suppressWarnings(fit_maihda(y ~ x + (1 | stratum), data = d1)))
  m2 <- suppressMessages(suppressWarnings(fit_maihda(y ~ x + (1 | stratum), data = d2)))
  expect_error(calculate_pcv(m1, m2), "same stratum")
})

# ---- integration: reordered AGGREGATED-BINOMIAL (cbind matrix) fits -------------
# The vector-response fixture above never exercised a matrix response, which is the
# form the response fingerprint mishandled: model.response() returns an n x 2
# cbind(success, failure) matrix (length 2n), so row alignment was skipped and two
# identical fits in different row order got different fingerprints -- a false
# "outcome values differ" in calculate_pcv(), a false "analytic sample" warning in
# compare_maihda(), and a suppressed delta in maihda_ic().

maihda_build_binom_fixture <- function() {
  set.seed(11)
  N  <- 150
  d1 <- factor(sample(c("m", "f"), N, replace = TRUE))
  d2 <- factor(sample(c("lo", "mid", "hi"), N, replace = TRUE))
  trials  <- sample(6:18, N, replace = TRUE)
  p       <- stats::plogis(-0.2 + 0.5 * (d1 == "m") + 0.3 * (d2 == "hi"))
  success <- stats::rbinom(N, trials, p)
  dat <- data.frame(success = success, failure = trials - success, d1 = d1, d2 = d2)
  rownames(dat) <- paste0("row", seq_len(N))
  dat
}

test_that("PCV/compare/IC accept a reordered aggregated-binomial (cbind) fit", {
  skip_on_cran()
  dat     <- maihda_build_binom_fixture()
  fit_adj <- function(d) suppressMessages(suppressWarnings(
    fit_maihda(cbind(success, failure) ~ d1 + d2 + (1 | d1:d2),
               data = d, family = stats::binomial())))
  null_fit <- suppressMessages(suppressWarnings(
    fit_maihda(cbind(success, failure) ~ 1 + (1 | d1:d2),
               data = dat, family = stats::binomial())))
  adj_fit  <- fit_adj(dat)
  set.seed(3)
  adj_shuf <- fit_adj(dat[sample(nrow(dat)), ])

  # Same row-id SET, different ORDER -- the analytic sample is identical.
  expect_true(setequal(MAIHDA:::maihda_wrapper_row_ids(adj_fit),
                       MAIHDA:::maihda_wrapper_row_ids(adj_shuf)))

  # The two-column (n x 2) matrix fingerprint is now order-invariant.
  expect_identical(MAIHDA:::maihda_wrapper_response_fingerprint(adj_fit),
                   MAIHDA:::maihda_wrapper_response_fingerprint(adj_shuf))

  # calculate_pcv no longer aborts and returns the same PCV as the in-order fit.
  base  <- calculate_pcv(null_fit, adj_fit)
  reord <- calculate_pcv(null_fit, adj_shuf)
  expect_equal(reord$pcv, base$pcv, tolerance = 1e-4)

  # compare_maihda emits NO false "analytic sample" warning.
  w <- character(0)
  withCallingHandlers(
    compare_maihda(null_fit, adj_shuf, ic = FALSE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("analytic sample", w, fixed = TRUE)))

  # maihda_ic reports no false delta issue.
  expect_length(MAIHDA:::maihda_ic_delta_issues(list(null_fit, adj_shuf)), 0L)
})

test_that("PCV still rejects a genuinely different aggregated-binomial sample", {
  skip_on_cran()
  # Hold the strata (d1, d2) and trials fixed so the SAME rows/strata/default names
  # are shared; only the success/failure counts differ, isolating the matrix
  # response fingerprint guard.
  set.seed(9)
  N      <- 120
  d1     <- factor(sample(c("m", "f"), N, replace = TRUE))
  d2     <- factor(sample(c("lo", "hi"), N, replace = TRUE))
  trials <- sample(6:18, N, replace = TRUE)
  mk <- function(seed) {
    set.seed(seed)
    success <- stats::rbinom(N, trials, stats::plogis(-0.2 + 0.5 * (d1 == "m")))
    d <- data.frame(success = success, failure = trials - success, d1 = d1, d2 = d2)
    suppressMessages(suppressWarnings(
      fit_maihda(cbind(success, failure) ~ d1 + d2 + (1 | d1:d2),
                 data = d, family = stats::binomial())))
  }
  expect_error(calculate_pcv(mk(51), mk(52)), "outcome values differ")
})
