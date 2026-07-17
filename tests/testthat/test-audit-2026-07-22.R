# Regression tests for the 2026-07-22 audit findings on commit f993586.
#
#   #1 P1  A character (row-name) `subset` was left unchanged by
#          maihda_normalize_subset(): its length != n_full, so in
#          compare_maihda_groups() slice_full() dropped it from the min_group_n
#          guard (guard ignored the subset) while slice_dots_for_group() still
#          forwarded the whole vector to each per-group fit. The guard therefore
#          counted rows the fit excludes -- undersized groups slipped past
#          min_group_n -- and lme4 applied the raw character subset via base `['s
#          character row-indexing, which PARTIAL-matches and inserts phantom NA rows
#          once per-group slicing has de-duplicated a name (r40 -> r40.1), pulling in
#          the wrong observations. Now every subset form (numeric / logical /
#          character) is resolved to a full-length positional logical mask up front
#          (character via rownames(data), by EXACT %in%), so the guard and the fit
#          see identical rows and the engine never sees a raw character subset.
#   #2 P2  maihda_resolve_strata_formula() detected duplicate stratum intercepts
#          from the formula intercept ATTRIBUTE only, so (1 | stratum) + (0 + one |
#          stratum) with a constant `one` (a slope column identical to the intercept)
#          reported "no intercept" for the second term and was accepted -- lme4 then
#          split the between-stratum variance across 'stratum' / 'stratum.1'. Now each
#          stratum term's random-effect design matrix is built on the ANALYTIC rows
#          and its column space is tested for the intercept vector
#          (maihda_re_lhs_spans_intercept), catching a constant slope column and full
#          dummy encodings while still allowing a genuinely varying (0 + x | stratum).

# ---- #1 P1: character subset unit resolution (exact, no partial matching) ------

test_that("maihda_normalize_subset resolves character names by EXACT match", {
  f <- MAIHDA:::maihda_normalize_subset
  rn <- c("r1", "r2", "r3", "r40.1")
  # A character subset becomes a full-length logical mask over the named rows.
  expect_identical(f(c("r2", "r3"), 4L, rn), c(FALSE, TRUE, TRUE, FALSE))
  # "r40" is absent exactly and must NOT partial-match "r40.1" (the base-`[` bug).
  expect_identical(f("r40", 4L, rn), c(FALSE, FALSE, FALSE, FALSE))
  # Without row names available, the character vector is returned unchanged.
  expect_identical(f("r2", 4L), "r2")
  # Numeric / logical forms keep their existing positional-mask behaviour.
  expect_identical(f(c(1L, 3L), 4L, rn), c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(f(c(TRUE, FALSE, TRUE, FALSE), 4L, rn), c(TRUE, FALSE, TRUE, FALSE))
})

# ---- #1 P1: character subset does not bypass the min_group_n guard -------------

test_that("compare_maihda_groups() applies a character subset to the min_group_n guard", {
  skip_on_cran()
  set.seed(1)
  N <- 120
  dat <- data.frame(
    y  = rnorm(N),
    d1 = factor(rep(c("m", "f"), length.out = N)),
    d2 = factor(rep(c("lo", "hi"), each = 2, length.out = N)),
    grp = factor(c(rep("A", 60), rep("B", 60))))
  rownames(dat) <- paste0("r", seq_len(N))
  # Group A rows are r1..r60, group B rows r61..r120. Eligible under the subset:
  # A = 8 (< min_group_n = 10, must be skipped), B = 25.
  subset_names <- c(paste0("r", 1:8), paste0("r", 61:85))

  res <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ 1 + (1 | d1:d2), data = dat, group = "grp",
                          family = "gaussian", subset = subset_names,
                          min_group_n = 10, shared_strata = TRUE)))
  a <- res[res$group == "A", ]
  b <- res[res$group == "B", ]
  # Previously the guard ignored the character subset, counted all 30 analytic rows,
  # and fitted group A on only 8 rows with status "ok".
  expect_match(a$status, "skipped")
  expect_equal(a$n, 8)
  expect_equal(b$n, 25)

  # A character subset and the equivalent logical mask must give identical outcomes.
  res_log <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ 1 + (1 | d1:d2), data = dat, group = "grp",
                          family = "gaussian", subset = rownames(dat) %in% subset_names,
                          min_group_n = 10, shared_strata = TRUE)))
  expect_identical(res$status, res_log$status)
  expect_identical(res$n, res_log$n)
})

# ---- #1 P1: fitted ROW IDENTITIES (not only counts) match the exact names ------

test_that("fit_maihda() fits exactly the character-named rows, not base-[ partial matches", {
  skip_on_cran()
  set.seed(505)
  ns <- 20; per <- 10; N <- ns * per
  d <- data.frame(y = rnorm(N),
                  stratum = factor(rep(seq_len(ns), each = per)))
  rn <- paste0("r", seq_len(N))
  # Drop the exact name "r40" and leave only its ".1" completion present, so base
  # `['s character indexing would UNIQUELY partial-match "r40" -> "r40.1" (row 40).
  rn[40] <- "r40.1"
  rownames(d) <- rn
  sel <- paste0("r", 1:40)            # names r1..r40; "r40" is now absent exactly

  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = d, subset = sel)))
  fr <- fit$model@frame
  # Exact resolution keeps r1..r39 only (39 rows); base `[` would have added the
  # wrong "r40.1" row (40 rows).
  expect_equal(nrow(fr), 39)
  # The row-40 observation ("r40.1") is NOT among the fitted rows.
  expect_false(d$y[40] %in% fr$y)
})

# ---- #2 P2: intercept-span detection on the design matrix ----------------------

test_that("maihda_re_lhs_spans_intercept detects numerically-realized intercepts", {
  g <- MAIHDA:::maihda_re_lhs_spans_intercept
  set.seed(9)
  dd <- data.frame(one = 1, x = rnorm(60),
                   f = factor(sample(c("a", "b", "c"), 60, replace = TRUE)))
  expect_true(g(quote(1), dd))
  expect_true(g(quote(1 + x), dd))
  expect_true(g(quote(x), dd))          # implicit intercept: lme4 reads (x|g) as (1+x|g)
  expect_true(g(quote(0 + one), dd))    # constant column IS the intercept -- the fix
  expect_true(g(quote(0 + f), dd))      # full dummy set: columns sum to the intercept
  expect_false(g(quote(0 + x), dd))     # genuinely varying slope, no intercept
})

test_that("fit_maihda() rejects a constant slope column masquerading as a stratum intercept", {
  skip_on_cran()
  set.seed(7)
  ns <- 40; per <- 12
  d <- data.frame(
    y = rnorm(ns * per),
    stratum = factor(rep(seq_len(ns), each = per)),
    one = 1,
    x = rnorm(ns * per),
    f = factor(sample(c("a", "b", "c"), ns * per, replace = TRUE)))
  d$y <- d$y + rnorm(ns, sd = 1)[as.integer(d$stratum)]

  # (0 + one | stratum) with constant `one` is a second stratum intercept.
  expect_error(
    suppressMessages(fit_maihda(y ~ 1 + (1 | stratum) + (0 + one | stratum), data = d)),
    "non-identifiable")
  # A full dummy encoding whose columns span the intercept, alongside (1 | stratum).
  expect_error(
    suppressMessages(fit_maihda(y ~ 1 + (1 | stratum) + (0 + f | stratum), data = d)),
    "non-identifiable")
  # A genuinely varying slope-only term stays identifiable and accepted.
  expect_no_error(
    suppressWarnings(suppressMessages(
      fit_maihda(y ~ 1 + (1 | stratum) + (0 + x | stratum), data = d))))
  # A single compound intercept+slope term is one intercept -- accepted.
  expect_no_error(
    suppressWarnings(suppressMessages(
      fit_maihda(y ~ 1 + (1 + x | stratum), data = d))))
})

test_that("intercept-span detection uses the ANALYTIC rows (constant only after subsetting)", {
  skip_on_cran()
  set.seed(13)
  ns <- 40; per <- 12; N <- ns * per
  d <- data.frame(
    y = rnorm(N),
    stratum = factor(rep(seq_len(ns), each = per)))
  d$y <- d$y + rnorm(ns, sd = 1)[as.integer(d$stratum)]
  # `z` varies in the raw data but is constant (== 5) on every row after the first
  # stratum, so on the analytic sample selected below it is a second intercept.
  d$z <- ifelse(seq_len(N) <= per, rnorm(N), 5)
  keep <- seq_len(N) > per

  expect_error(
    suppressMessages(fit_maihda(y ~ 1 + (1 | stratum) + (0 + z | stratum),
                                data = d, subset = keep)),
    "non-identifiable")
})
