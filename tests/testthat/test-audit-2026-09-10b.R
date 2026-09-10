# Audit 2026-09-10 (second pass): an equivalent coding of the adjusted model must
# not change the derived null model.
#
# `y ~ x + a + b + (1 | a:b)` and `y ~ 0 + x + a + b + (1 | a:b)` are the SAME
# adjusted model -- with no intercept R codes the first factor as cell means, so the
# two designs span one space and fit identical values -- but the null was derived
# with update(. ~ . - a - b), which took the grand mean out with the dimensions only
# in the second spelling. The null's stratum random intercept then absorbed the
# outcome mean: on the data below the null between-stratum variance went from 3.01
# to 1941.78, the null VPC from 0.40 to 0.998 and the PCV from 0.600 to 0.999, with
# the adjusted fit bit-identical either way. The same reduction feeds the
# crossed-dimensions base formula, the longitudinal null and the per-group
# formulas, so all four decompositions moved.

maihda_a5_data <- function(n = 900, seed = 11) {
  set.seed(seed)
  a <- factor(sample(letters[1:4], n, TRUE))
  b <- factor(sample(LETTERS[1:3], n, TRUE))
  x <- stats::rnorm(n)
  f <- factor(sample(c("p", "q"), n, TRUE))
  st <- interaction(a, b, drop = TRUE)
  u <- stats::rnorm(nlevels(st), 0, 1.2)[as.integer(st)]
  # A LARGE outcome mean is what makes the defect unmistakable: it is the mean the
  # random intercept is forced to absorb when the fixed part cannot carry it.
  y <- 40 + 0.8 * x + 0.5 * (f == "q") + as.numeric(a) * 1.1 +
    as.numeric(b) * 0.6 + u + stats::rnorm(n, 0, 2)
  data.frame(y = y, x = x, a = a, b = b, f = f,
             g = factor(sample(c("G1", "G2"), n, TRUE)))
}

# Cross-parameterisation comparisons below use tolerance = 1e-6: the two spellings
# span one column space but are two separate lme4 optimisations, so they agree to
# optimiser precision (~1e-8 relative), not to machine precision. The defect being
# pinned moves these numbers by ~0.4 absolute, so the margin is ample.
maihda_a5_int <- function(f) attr(stats::terms(maihda_nobars(f)), "intercept")
maihda_a5_var <- function(m) as.data.frame(lme4::VarCorr(m$model))$vcov[1]

# Count the conditions an expression signals, muffling them. A handler established
# INSIDE a suppressWarnings() never sees the condition (the suppressor muffles
# first), so the counting handler must be the innermost one.
maihda_a5_conditions <- function(expr) {
  msgs <- character(0)
  warns <- character(0)
  withCallingHandlers(
    force(expr),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(messages = msgs, warnings = warns)
}


test_that("the two adjusted codings really are the same model", {
  d <- maihda_a5_data()
  x_int <- stats::model.matrix(y ~ x + a + b, d)
  x_cell <- stats::model.matrix(y ~ 0 + x + a + b, d)
  # Same rank, and a combined rank that adds nothing: one column space, two spellings.
  expect_identical(qr(x_int)$rank, qr(x_cell)$rank)
  expect_identical(qr(cbind(x_int, x_cell))$rank, qr(x_int)$rank)
  m_int <- lme4::lmer(y ~ x + a + b + (1 | a:b), data = d, REML = TRUE)
  m_cell <- lme4::lmer(y ~ 0 + x + a + b + (1 | a:b), data = d, REML = TRUE)
  expect_equal(unname(stats::fitted(m_int)), unname(stats::fitted(m_cell)),
               tolerance = 1e-6)
  expect_equal(as.numeric(stats::logLik(m_int)), as.numeric(stats::logLik(m_cell)),
               tolerance = 1e-6)
})


test_that("the derived null keeps the grand mean under cell-means coding", {
  d <- maihda_a5_data()
  r_int <- maihda(y ~ x + a + b + (1 | a:b), data = d, autobin = FALSE)
  r_cell <- suppressMessages(
    maihda(y ~ 0 + x + a + b + (1 | a:b), data = d, autobin = FALSE))

  # The null formulas now agree, and both carry an intercept.
  expect_identical(maihda_a5_int(r_int$formula), 1L)
  expect_identical(maihda_a5_int(r_cell$formula), 1L)

  # The adjusted fits were always identical -- only the null moved.
  expect_equal(maihda_a5_var(r_int$model_adjusted),
               maihda_a5_var(r_cell$model_adjusted), tolerance = 1e-6)
  # ... and now the null does not.
  expect_equal(maihda_a5_var(r_int$model), maihda_a5_var(r_cell$model),
               tolerance = 1e-6)
  expect_equal(r_int$summary$vpc, r_cell$summary$vpc, tolerance = 1e-6)
  expect_equal(r_int$pcv$pcv, r_cell$pcv$pcv, tolerance = 1e-6)
  # Absolute anchors: the correct values, not merely equal ones. Before the fix the
  # cell-means arm reported 1941.78 / 0.9977 / 0.9994.
  expect_equal(maihda_a5_var(r_cell$model), 3.005578, tolerance = 1e-4)
  expect_equal(r_cell$pcv$pcv, 0.600439, tolerance = 1e-4)
  expect_lt(maihda_a5_var(r_cell$model), 10)
})


test_that("the partial branch keeps the grand mean too", {
  # Only one dimension written in the formula: maihda() fits a clean null and then
  # builds the adjusted from it, so the null's lost intercept propagates.
  d <- maihda_a5_data()
  p_int <- suppressMessages(maihda(y ~ x + a + (1 | a:b), data = d, autobin = FALSE))
  p_cell <- suppressMessages(
    maihda(y ~ 0 + x + a + (1 | a:b), data = d, autobin = FALSE))
  expect_identical(maihda_a5_int(p_cell$formula), 1L)
  expect_equal(p_int$pcv$pcv, p_cell$pcv$pcv, tolerance = 1e-6)
  expect_equal(p_cell$pcv$pcv, 0.600439, tolerance = 1e-4)
})


test_that("the crossed-dimensions partition is invariant to the coding", {
  skip_on_cran()
  d <- maihda_a5_data()
  c_int <- suppressWarnings(maihda(y ~ x + a + b + (1 | a:b), data = d,
                                   autobin = FALSE,
                                   decomposition = "crossed-dimensions"))
  c_cell <- suppressWarnings(suppressMessages(
    maihda(y ~ 0 + x + a + b + (1 | a:b), data = d, autobin = FALSE,
           decomposition = "crossed-dimensions")))
  expect_equal(c_int$decomposition$additive_share,
               c_cell$decomposition$additive_share, tolerance = 1e-6)
  expect_equal(unname(c_int$decomposition$per_dim),
               unname(c_cell$decomposition$per_dim), tolerance = 1e-6)
  # Before the fix dimension `a`'s additive random intercept absorbed the mean:
  # per_dim a was 1940.6 against the correct 1.94, and additive_share 0.9994.
  expect_lt(c_cell$decomposition$per_dim[["a"]], 10)
  expect_lt(c_cell$decomposition$additive_share, 0.9)
})


test_that("the per-group decomposition is invariant to the coding", {
  skip_on_cran()
  d <- maihda_a5_data()
  g_int <- as.data.frame(compare_maihda_groups(
    y ~ x + a + b + (1 | a:b), data = d, group = "g", autobin = FALSE))
  g_cell <- as.data.frame(suppressMessages(compare_maihda_groups(
    y ~ 0 + x + a + b + (1 | a:b), data = d, group = "g", autobin = FALSE)))
  expect_equal(g_int$var_between, g_cell$var_between, tolerance = 1e-6)
  expect_equal(g_int$vpc, g_cell$vpc, tolerance = 1e-6)
  expect_equal(g_int$pcv, g_cell$pcv, tolerance = 1e-6)
  # Every group's between-stratum variance was ~1900 before the fix.
  expect_true(all(g_cell$var_between < 10))
})


test_that("the longitudinal PCV is invariant to the coding", {
  skip_on_cran()
  set.seed(3)
  nid <- 240
  waves <- 4
  gd <- factor(sample(c("W", "M"), nid, TRUE))
  ed <- factor(sample(c("Lo", "Hi", "Md"), nid, TRUE))
  sti <- interaction(gd, ed, drop = TRUE)
  ui <- stats::rnorm(nlevels(sti), 0, 1)[as.integer(sti)]
  si <- stats::rnorm(nlevels(sti), 0, 0.3)[as.integer(sti)]
  pid <- rep(seq_len(nid), each = waves)
  wv <- rep(0:(waves - 1), nid)
  yl <- 30 + 2 * wv + (as.numeric(gd) * 1.2 + as.numeric(ed) * 0.8 + ui)[pid] +
    si[pid] * wv + stats::rnorm(nid, 0, 1)[pid] + stats::rnorm(nid * waves, 0, 1)
  ld <- data.frame(id = pid, wave = wv, gender = gd[pid], edu = ed[pid],
                   wb = yl, cov = stats::rnorm(nid * waves))

  l_int <- suppressWarnings(maihda(
    wb ~ wave + cov + gender + edu + (1 | gender:edu), data = ld,
    id = "id", time = "wave", autobin = FALSE))
  l_cell <- suppressWarnings(suppressMessages(maihda(
    wb ~ 0 + wave + cov + gender + edu + (1 | gender:edu), data = ld,
    id = "id", time = "wave", autobin = FALSE)))
  expect_equal(l_int$pcv$pcv_intercept, l_cell$pcv$pcv_intercept, tolerance = 1e-6)
  expect_equal(l_int$pcv$pcv_slope, l_cell$pcv$pcv_slope, tolerance = 1e-6)
  # The slope PCV is the one that moved furthest: 0.2019 -> 0.9934 before the fix.
  expect_lt(l_cell$pcv$pcv_slope, 0.5)
})


test_that("maihda_fixed_spans_intercept() reads the design, not the formula", {
  d <- maihda_a5_data(n = 300, seed = 4)
  expect_true(maihda_fixed_spans_intercept(y ~ x + (1 | a:b), d))
  # A bare numeric with no intercept provably does not span it ...
  expect_false(maihda_fixed_spans_intercept(y ~ x + (1 | a:b) - 1, d))
  # ... but a factor coded as cell means does, though the formula says `- 1`.
  expect_true(maihda_fixed_spans_intercept(y ~ f + (1 | a:b) - 1, d))
  expect_true(maihda_fixed_spans_intercept(y ~ a + (1 | a:b) - 1, d))
  # Stripping the bars can leave nothing but `1`, which lme4 fits with an intercept.
  expect_identical(maihda_a5_int(y ~ (1 | a:b) - 1), 1L)
  expect_true(maihda_fixed_spans_intercept(y ~ (1 | a:b) - 1, d))
  # An unbuildable design is NA -- never a guess from the written `0 +`.
  expect_true(is.na(maihda_fixed_spans_intercept(y ~ 0 + not_a_column, d)))
  expect_true(is.na(maihda_fixed_spans_intercept(y ~ x - 1, d[0, ])))
})


test_that("maihda_drop_fixed_terms() restores the mean only where it is lost", {
  d <- maihda_a5_data(n = 300, seed = 4)
  dims <- c("a", "b")

  # (1) Already-intercept source: the reduction is exactly what update() gives.
  keep <- maihda_drop_fixed_terms(y ~ x + a + b + (1 | a:b), dims, d)
  expect_identical(maihda_a5_int(keep), 1L)
  expect_identical(deparse1(keep), deparse1(
    stats::update(y ~ x + a + b + (1 | a:b), . ~ . - a - b)))

  # (2) Cell-means source, numeric survivor: the grand mean comes back.
  fixed <- maihda_drop_fixed_terms(y ~ 0 + x + a + b + (1 | a:b), dims, d)
  expect_identical(maihda_a5_int(fixed), 1L)

  # (3) Factor survivor: the reduction still spans the mean, so it is left alone --
  # `- 1` is preserved and the fixed-effect table keeps its cell-means coding.
  left <- maihda_drop_fixed_terms(y ~ 0 + f + a + b + (1 | a:b), dims, d)
  expect_identical(maihda_a5_int(left), 0L)

  # (4) A source that genuinely has no grand mean (numeric dimensions, a regression
  # through the origin): adding one would put a column in the null the adjusted
  # model does not have, so the parameterization is kept.
  dn <- d
  dn$a <- as.numeric(dn$a)
  dn$b <- as.numeric(dn$b)
  origin <- maihda_drop_fixed_terms(y ~ 0 + x + a + b + (1 | a:b), dims, dn)
  expect_identical(maihda_a5_int(origin), 0L)

  # (5) Nothing to remove is the identity, including on a `0 +` formula.
  expect_identical(
    deparse1(maihda_drop_fixed_terms(y ~ 0 + x + (1 | a:b), character(0), d)),
    deparse1(y ~ 0 + x + (1 | a:b)))
})


test_that("an intercept-bearing decomposition is untouched", {
  # Negative control: the formulas that were always right must not move, and the fix
  # must not blanket-add an intercept wherever `- 1` appears.
  d <- maihda_a5_data()
  k_int <- maihda(y ~ f + a + b + (1 | a:b), data = d, autobin = FALSE)
  k_cell <- maihda(y ~ 0 + f + a + b + (1 | a:b), data = d, autobin = FALSE)
  expect_identical(maihda_a5_int(k_cell$formula), 0L)
  expect_equal(maihda_a5_var(k_int$model), maihda_a5_var(k_cell$model),
               tolerance = 1e-6)
  expect_equal(k_int$pcv$pcv, k_cell$pcv$pcv, tolerance = 1e-6)

  z_int <- maihda(y ~ a + b + (1 | a:b), data = d, autobin = FALSE)
  z_cell <- maihda(y ~ 0 + a + b + (1 | a:b), data = d, autobin = FALSE)
  expect_identical(maihda_a5_int(z_cell$formula), 1L)
  expect_equal(maihda_a5_var(z_int$model), maihda_a5_var(z_cell$model),
               tolerance = 1e-6)
})


test_that("the grand-mean message fires once, and only where it applies", {
  skip_on_cran()
  d <- maihda_a5_data(n = 500, seed = 6)
  hit <- function(x) sum(grepl("keeps the grand mean", x))

  one <- maihda_a5_conditions(
    maihda(y ~ 0 + x + a + b + (1 | a:b), data = d, autobin = FALSE))
  expect_identical(hit(one$messages), 1L)

  # A group comparison repeats the same reduction per group; it must not repeat the
  # message.
  grp <- maihda_a5_conditions(
    maihda(y ~ 0 + x + a + b + (1 | a:b), data = d, autobin = FALSE, group = "g"))
  expect_identical(hit(grp$messages), 1L)

  # compare_maihda_groups() reports it itself when called directly.
  direct <- maihda_a5_conditions(compare_maihda_groups(
    y ~ 0 + x + a + b + (1 | a:b), data = d, group = "g", autobin = FALSE))
  expect_identical(hit(direct$messages), 1L)

  # Silent where nothing was restored.
  plain <- maihda_a5_conditions(
    maihda(y ~ x + a + b + (1 | a:b), data = d, autobin = FALSE))
  expect_identical(hit(plain$messages), 0L)
  factor_cov <- maihda_a5_conditions(
    maihda(y ~ 0 + f + a + b + (1 | a:b), data = d, autobin = FALSE))
  expect_identical(hit(factor_cov$messages), 0L)
})


test_that("a fit whose fixed part cannot carry the mean is warned about", {
  skip_on_cran()
  d <- maihda_a5_data(n = 500, seed = 6)
  hit <- function(x) sum(grepl("cannot represent the outcome", x))

  # A null the user supplied without a grand mean: not an equivalent recoding, so it
  # is fitted as written -- but the stratum variance absorbs the mean and the
  # VPC/MOR/PCV stop describing the strata, so it is flagged.
  supplied <- maihda_a5_conditions(
    maihda(y ~ 0 + x + (1 | a:b), data = d, autobin = FALSE))
  expect_identical(hit(supplied$warnings), 1L)

  expect_identical(
    hit(maihda_a5_conditions(fit_maihda(y ~ 0 + x + (1 | a:b),
                                        data = d, autobin = FALSE))$warnings), 1L)

  # Silent on every sound fit, including the `- 1` spellings that still span.
  expect_identical(
    hit(maihda_a5_conditions(fit_maihda(y ~ x + (1 | a:b),
                                        data = d, autobin = FALSE))$warnings), 0L)
  expect_identical(
    hit(maihda_a5_conditions(fit_maihda(y ~ 0 + f + (1 | a:b),
                                        data = d, autobin = FALSE))$warnings), 0L)
  expect_identical(
    hit(maihda_a5_conditions(
      maihda(y ~ 0 + x + a + b + (1 | a:b), data = d, autobin = FALSE))$warnings), 0L)
})
