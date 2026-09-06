# Audit 2026-09-06
#
#   A1 [High] CONFIRMED -- summary(df_method = "bootstrap") did not impose the
#   null it advertises. For each fixed-effect term the reference distribution was
#   simulated from stats::update(formula, . ~ . - <term>), but removing a term
#   from a formula does not always remove it from the MODEL. R's marginality
#   rules recode a surviving higher-order term to absorb a dropped marginal one:
#   for `y ~ x * f`, `. ~ . - x` gives `f + x:f`, which R codes with a full dummy
#   expansion (fa:x, fb:x) instead of a contrast, so the "reduced" design spans
#   exactly the same column space as the full one -- same rank, identical fitted
#   values, identical logLik. Simulating from it draws data that still carry the
#   estimated effect, so the reference distribution of |t*| centres on the
#   OBSERVED statistic and the p-value collapses towards 0.5 no matter how large
#   the effect is. On the fit below (x estimate 2.013, SE 0.059, statistic 34.0)
#   the bootstrap reported p = 0.54 and an interval of [-0.170, 4.196] covering
#   zero. The null distribution it built had median |t*| = 34.17.
#
#   The failure is degeneracy, not anticonservatism. Because the "reduced" fit
#   reproduces the full model, its draws carry the estimated effect and |t*|
#   lands on |t_obs| whatever the truth, so the p-value is pinned near 0.5. Over
#   300 replicates (n = 400, 12 strata, B = 99) with beta_x truly zero and x:f
#   non-zero, the formula reduction rejected at 0.0000 against a nominal 5%
#   (median p = 0.590). Its power was 0.0000 at beta_x = 0.30 (200 reps, median
#   p 0.505) AND at beta_x = 0.60 (150 reps, median p 0.510). A median p of
#   0.59 / 0.51 / 0.51 across the three truths is the signature: the statistic is
#   insensitive to the effect, so the reported non-significance was uninformative
#   rather than merely conservative -- the test could not reject at any size.
#
#   The repaired reference, scored through the SHIPPED entry point
#   (summary(fit, df_method = "bootstrap")) rather than a prototype loop:
#   rejection 0.0650 under the null (200 reps, MC SE 0.0174, so the interval
#   covers the nominal 5%; median p 0.490) and power 0.9750 at beta_x = 0.30
#   (120 reps). A prototype of the restriction alone gave 0.0533 and 0.9750 on
#   the same design, which is the same answer within Monte Carlo error.
#
#   The affected shapes are every term that is marginal to a higher-order term
#   still in the model: the continuous main effect of `x * f`, BOTH factor main
#   effects of `f * g`, all main effects and all two-way terms under `f * g * h`,
#   a nested `f / g`, and a transformed main effect such as poly(x, 2) under an
#   interaction. `x1 * x2` (two continuous) and any additive fixed part are
#   sound, which is why the 2026-09-02 pass's GLMM calibration -- measured on
#   `ev ~ d1 + d2 + d3` -- is unaffected and stands. A NO-INTERCEPT fit is wrong
#   for a second reason: update() turns `y ~ 0 + x` into `y ~ (1 | st) - 1`,
#   which lme4 fits WITH an intercept the full model never had.
#
#   Fixed by verifying the reduction rather than assuming it: the refitted design
#   must not span the full model's column space (maihda_same_column_space). Where
#   it does, the constraint is imposed on the fitted design columns directly
#   (maihda_restrict_fixef), which needs no formula algebra and so handles
#   transformed terms, weights, offsets and matrix responses alike. Where the
#   formula reduction is already correct it is still used, so every previously
#   valid p-value is unchanged.
#
#   A1b [High] CONFIRMED (found while validating A1, unrelated mechanism) --
#   every parametric bootstrap in the package failed outright on any lme4 fit
#   whose data carried an NA. lme4::refit() treats a `newresp` with no
#   "na.action" attribute as being on the ORIGINAL data scale and drops the rows
#   the fit's na.action removed; a simulated draw is already on the model frame's
#   rows, so it was truncated a second time and the refit errored with
#   "replacement has <n - k> rows, data has <n>". Every draw failed, so
#   summary(bootstrap = TRUE), summary(df_method = "bootstrap"),
#   calculate_pcv(bootstrap = TRUE) and the longitudinal VPC band all raised
#   "All ... refits failed". The bundled maihda_health_data has 208 such rows.
#   Fixed by tagging each draw with the na.action of the model being refitted
#   INTO (maihda_refit_draw), which is the "already aligned" signal refit()
#   tests for. FOUR of the five refit sites were affected, not five:
#   pcv_importance() routes through the same helper but reduces to ONE
#   complete-case analytic sample before fitting, so its subset models carry no
#   na.action and it worked before the fix as well -- measured, 0 of 96 draws
#   needed the tag. Its routing is defensive consistency, not a repair.

maihda_a1_data <- function(seed = 1, n = 600) {
  set.seed(seed)
  f <- factor(sample(c("a", "b"), n, TRUE))
  x <- stats::rnorm(n)
  st <- factor(sample(paste0("s", 1:6), n, TRUE))
  u <- stats::rnorm(6, 0, 0.5)[as.integer(st)]
  data.frame(
    y = 1 + 2 * x + 0.5 * (f == "b") + 0.3 * x * (f == "b") + u + stats::rnorm(n),
    x = x, f = f, stratum = st)
}

# ---------------------------------------------------------------- A1: mechanism

test_that("dropping a marginal term from a formula does not impose its null", {
  d <- maihda_a1_data()
  X_full <- stats::model.matrix(~ x * f, d)
  # This is exactly what stats::update(. ~ . - x) hands lmer.
  X_red <- stats::model.matrix(~ f + x:f, d)
  # The defect in one line: the "reduced" design has the SAME rank as the full
  # one, because x:f is expanded to fa:x AND fb:x once the x main effect is gone.
  expect_identical(qr(X_red)$rank, qr(X_full)$rank)
  expect_true(maihda_same_column_space(X_red, X_full))
  # The null it was supposed to impose -- the x column removed -- is one smaller.
  X_null <- X_full[, colnames(X_full) != "x", drop = FALSE]
  expect_identical(qr(X_null)$rank, qr(X_full)$rank - 1L)
  expect_false(maihda_same_column_space(X_red, X_null))
})

test_that("maihda_same_column_space() separates a real reduction from a recoding", {
  d <- maihda_a1_data()
  # Additive: dropping a term genuinely removes a column.
  X <- stats::model.matrix(~ x + f, d)
  expect_true(maihda_same_column_space(
    stats::model.matrix(~ f, d), X[, colnames(X) != "x", drop = FALSE]))
  # Equal rank is not enough on its own: a same-sized design spanning a
  # DIFFERENT direction must be rejected.
  A <- stats::model.matrix(~ x, d)
  B <- stats::model.matrix(~ f, d)
  expect_identical(qr(A)$rank, qr(B)$rank)
  expect_false(maihda_same_column_space(A, B))
  expect_false(maihda_same_column_space(A, NULL))
})

test_that("the formula reduction fails for every marginal-term shape, and only those", {
  set.seed(2)
  n <- 300
  d <- data.frame(
    y = stats::rnorm(n),
    x1 = stats::rnorm(n), x2 = stats::rnorm(n),
    f = factor(sample(c("a", "b"), n, TRUE)),
    g = factor(sample(c("p", "q", "r"), n, TRUE)),
    h = factor(sample(c("u", "v"), n, TRUE)))
  # TRUE when update(. ~ . - term) really does impose that term's null.
  reduces <- function(fml, term) {
    tt <- stats::terms(fml, data = d)
    X <- stats::model.matrix(tt, d)
    cols <- which(attr(X, "assign") == match(term, attr(tt, "term.labels")))
    X_null <- X[, setdiff(seq_len(ncol(X)), cols), drop = FALSE]
    X_red <- stats::model.matrix(
      stats::terms(stats::update(fml, paste(". ~ . -", term)), data = d), d)
    maihda_same_column_space(X_red, X_null)
  }
  # Broken: a term marginal to a surviving higher-order term.
  expect_false(reduces(y ~ x1 * f, "x1"))
  expect_false(reduces(y ~ f * g, "f"))
  expect_false(reduces(y ~ f * g, "g"))
  expect_false(reduces(y ~ f * g * h, "f"))
  expect_false(reduces(y ~ f * g * h, "f:g"))
  expect_false(reduces(y ~ f / g, "f"))
  expect_false(reduces(y ~ poly(x1, 2) * f, "poly(x1, 2)"))
  # Sound: nothing higher-order survives to absorb the dropped term.
  expect_true(reduces(y ~ x1 * f, "f"))
  expect_true(reduces(y ~ x1 * f, "x1:f"))
  expect_true(reduces(y ~ x1 * x2, "x1"))
  expect_true(reduces(y ~ x1 + f, "x1"))
  expect_true(reduces(y ~ f * g, "f:g"))
})

# ------------------------------------------------------------------- A1: repair

test_that("maihda_restrict_fixef() imposes exactly the tested term's null", {
  skip_on_cran()
  d <- maihda_a1_data()
  m <- lme4::lmer(y ~ x * f + (1 | stratum), data = d, REML = TRUE)
  X <- lme4::getME(m, "X")
  j <- which(colnames(X) == "x")
  red <- maihda_restrict_fixef(m, j)
  expect_s4_class(red, "lmerMod")
  X_red <- lme4::getME(red, "X")
  # Exactly the full design minus the x column -- not merely the right rank.
  expect_true(maihda_same_column_space(
    X_red, X[, setdiff(seq_len(ncol(X)), j), drop = FALSE]))
  expect_identical(qr(X_red)$rank, qr(X)$rank - 1L)
  # It is a genuinely different fit, unlike the formula reduction, whose
  # log-likelihood matches the full model's to the optimizer's tolerance.
  expect_lt(as.numeric(stats::logLik(red)), as.numeric(stats::logLik(m)) - 1)
  expect_identical(stats::nobs(red), stats::nobs(m))
  # Simulating from it is what the bootstrap needs.
  expect_length(maihda_simulate_lme4(red, nsim = 2), 2L)
})

test_that("maihda_restrict_fixef() carries weights, offsets and a matrix response", {
  skip_on_cran()
  set.seed(4)
  n <- 500
  d <- data.frame(
    x = stats::rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
    st = factor(sample(paste0("s", 1:8), n, TRUE)),
    w = stats::runif(n, 0.5, 2), expo = stats::runif(n, 1, 5), trials = 20L)
  u <- stats::rnorm(8, 0, 0.5)[as.integer(d$st)]
  eta <- -0.3 + 0.8 * d$x + 0.4 * (d$f == "b") + 0.5 * d$x * (d$f == "b") + u
  d$y <- 1 + 2 * d$x + u + stats::rnorm(n)
  d$succ <- stats::rbinom(n, d$trials, stats::plogis(eta))
  d$cnt <- stats::rpois(n, exp(eta - 1) * d$expo)

  ok <- function(m) {
    X <- lme4::getME(m, "X")
    j <- which(colnames(X) == "x")
    red <- maihda_restrict_fixef(m, j)
    !is.null(red) && maihda_same_column_space(
      lme4::getME(red, "X"), X[, setdiff(seq_len(ncol(X)), j), drop = FALSE])
  }
  expect_true(ok(lme4::lmer(y ~ x * f + (1 | st), d, REML = TRUE, weights = w)))
  expect_true(ok(suppressWarnings(lme4::glmer(
    cbind(succ, trials - succ) ~ x * f + (1 | st), d, family = stats::binomial()))))
  expect_true(ok(suppressWarnings(lme4::glmer(
    cnt ~ x * f + offset(log(expo)) + (1 | st), d, family = stats::poisson()))))
  expect_true(ok(lme4::lmer(y ~ x * f + (1 + x | st), d, REML = TRUE)))
})

test_that("maihda_restrict_fixef() handles a model with nothing left to keep", {
  skip_on_cran()
  # A no-intercept fit whose only fixed term is the one under test leaves NO
  # columns to retain, and the correct null is a zero-column fixed design (lme4
  # fits one). Two traps here, both hit during this audit's self-check:
  #   - paste0(stem, "x", integer(0)) is the single string "<stem>x", NOT
  #     character(0), so the restricted formula named a column the loop never
  #     created and the whole restriction failed to build;
  #   - the formula reduction is wrong here for a SECOND reason -- update() turns
  #     `y ~ 0 + x` into `y ~ (1 | st) - 1`, which lme4 fits WITH an intercept
  #     the full model never had, so its draws are from a different model
  #     entirely rather than from beta_x = 0.
  set.seed(5)
  n <- 400
  d <- data.frame(x = stats::rnorm(n),
                  st = factor(sample(paste0("s", 1:8), n, TRUE)))
  u <- stats::rnorm(8, 0, 0.5)[as.integer(d$st)]
  d$y <- 0.6 * d$x + u + stats::rnorm(n)
  m <- lme4::lmer(y ~ 0 + x + (1 | st), data = d, REML = TRUE)
  expect_identical(ncol(lme4::getME(m, "X")), 1L)

  red <- maihda_restrict_fixef(m, 1L)
  expect_s4_class(red, "lmerMod")
  expect_identical(ncol(lme4::getME(red, "X")), 0L)
  expect_length(maihda_simulate_lme4(red, nsim = 2), 2L)

  # The formula reduction must be REJECTED here, not silently used.
  old <- maihda_refit_reduced(m, stats::update(stats::formula(m), ". ~ . - x"))
  expect_true("(Intercept)" %in% names(lme4::fixef(old)))
  expect_false(maihda_same_column_space(
    lme4::getME(old, "X"), lme4::getME(m, "X")[, 0, drop = FALSE]))

  fit <- fit_maihda(y ~ 0 + x + (1 | st), data = d)
  set.seed(2)
  fe <- suppressWarnings(summary(fit, df_method = "bootstrap", n_boot = 24))$fixed_effects
  expect_gt(abs(fe$statistic[fe$term == "x"]), 5)
  expect_lte(fe$p_value[fe$term == "x"], 1 / 24)
  expect_gt(fe$lower[fe$term == "x"], 0)
})

test_that("df_method = 'bootstrap' rejects a huge effect under an interaction", {
  skip_on_cran()
  d <- maihda_a1_data()
  fit <- fit_maihda(y ~ x * f + (1 | stratum), data = d)
  set.seed(11)
  s <- suppressWarnings(summary(fit, df_method = "bootstrap", n_boot = 24))
  fe <- s$fixed_effects
  x_row <- fe$term == "x"
  # The statistic is enormous; before the fix the bootstrap returned p = 0.54 and
  # an interval spanning zero, because its null draws carried the effect itself.
  expect_gt(abs(fe$statistic[x_row]), 20)
  expect_lte(fe$p_value[x_row], 1 / 24)
  expect_gt(fe$lower[x_row], 0)
  # The interval and the p-value still agree, as they do on the sound path.
  expect_true(all(fe$lower[fe$term != "(Intercept)"] > 0))
})

test_that("a sound formula reduction is still the path taken", {
  skip_on_cran()
  set.seed(7)
  n <- 400
  d <- data.frame(x = stats::rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                  st = factor(sample(paste0("s", 1:8), n, TRUE)))
  u <- stats::rnorm(8, 0, 0.5)[as.integer(d$st)]
  d$y <- 1 + 0.4 * d$x + 0.3 * (d$f == "b") + u + stats::rnorm(n)
  m <- lme4::lmer(y ~ x + f + (1 | st), data = d, REML = TRUE)
  X <- lme4::getME(m, "X")
  labs <- attr(stats::terms(stats::formula(m, fixed.only = TRUE)), "term.labels")
  for (k in seq_along(labs)) {
    cols <- which(attr(X, "assign") == k)
    red <- maihda_refit_reduced(
      m, stats::update(stats::formula(m), paste(". ~ . -", labs[k])))
    # An additive model reaches maihda_restrict_fixef() never: the formula
    # reduction already is the null, so these results are untouched by the fix.
    expect_false(maihda_same_column_space(lme4::getME(red, "X"), X))
    expect_true(maihda_same_column_space(
      lme4::getME(red, "X"), X[, setdiff(seq_len(ncol(X)), cols), drop = FALSE]))
  }
})

# ------------------------------------------------------------------------- A1b

maihda_a1b_data <- function(seed = 1, n = 400, n_na = 20) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n),
                  st = factor(sample(paste0("s", 1:8), n, TRUE)))
  u <- stats::rnorm(8, 0, 0.6)[as.integer(d$st)]
  d$y <- 1 + 0.5 * d$x + u + stats::rnorm(n)
  d$y[sample(n, n_na)] <- NA
  d
}

test_that("maihda_refit_draw() survives a fit whose data carried NAs", {
  skip_on_cran()
  d <- maihda_a1b_data()
  m <- lme4::lmer(y ~ x + (1 | st), data = d, REML = TRUE)
  expect_false(is.null(attr(m@frame, "na.action")))
  sim <- maihda_simulate_lme4(m, nsim = 1)
  expect_identical(length(sim[[1]]), nrow(m@frame))
  # The bug, pinned: lme4 truncates an untagged draw by the fit's own na.action.
  expect_error(lme4::refit(m, newresp = sim[[1]]))
  # The fix: tagged with that na.action, the draw is taken as already aligned.
  expect_s4_class(maihda_refit_draw(m, sim[[1]]), "lmerMod")
})

test_that("maihda_refit_draw() leaves a complete-data refit bit-identical", {
  skip_on_cran()
  d <- maihda_a1b_data(n_na = 0)
  m <- lme4::lmer(y ~ x + (1 | st), data = d, REML = TRUE)
  expect_null(attr(m@frame, "na.action"))
  sim <- maihda_simulate_lme4(m, nsim = 1)
  expect_equal(lme4::fixef(maihda_refit_draw(m, sim[[1]])),
               lme4::fixef(lme4::refit(m, newresp = sim[[1]])))
})

test_that("every bootstrap path runs on data with missing values", {
  skip_on_cran()
  d <- maihda_a1b_data()
  fit <- fit_maihda(y ~ x + (1 | st), data = d)
  null_fit <- fit_maihda(y ~ 1 + (1 | st), data = d)
  expect_false(is.null(attr(fit$model@frame, "na.action")))

  s <- suppressWarnings(summary(fit, bootstrap = TRUE, n_boot = 15))
  expect_true(all(is.finite(c(s$vpc$ci_lower, s$vpc$ci_upper))))

  sb <- suppressWarnings(summary(fit, df_method = "bootstrap", n_boot = 15))
  tested <- sb$fixed_effects$term != "(Intercept)"
  expect_true(all(is.finite(sb$fixed_effects$p_value[tested])))

  p <- suppressWarnings(calculate_pcv(null_fit, fit, bootstrap = TRUE, n_boot = 15))
  expect_true(all(is.finite(c(p$ci_lower, p$ci_upper))))
})

test_that("the longitudinal VPC band also survives missing values", {
  skip_on_cran()
  # The fourth affected refit site (longitudinal.R). Before the fix this raised
  # "All VPC bootstrap refits failed"; every one of its draws needs the tag.
  set.seed(3)
  nid <- 120
  d <- do.call(rbind, lapply(seq_len(nid), function(i)
    data.frame(id = i, time = 0:3,
               g = sample(c("a", "b"), 1), h = sample(c("p", "q"), 1))))
  grp <- interaction(d$g, d$h, drop = TRUE)
  us <- stats::rnorm(nlevels(grp), 0, 0.4)[as.integer(grp)]
  ui <- stats::rnorm(nid, 0, 0.5)[d$id]
  d$y <- 1 + 0.3 * d$time + us + ui + stats::rnorm(nrow(d))
  d$y[sample(nrow(d), 30)] <- NA

  a <- suppressWarnings(suppressMessages(
    maihda(y ~ 1 + (1 | g:h), data = d, id = "id", time = "time",
           bootstrap = TRUE, n_boot = 12)))
  lng <- a$summary_null$longitudinal
  if (is.null(lng)) lng <- a$summary_adjusted$longitudinal
  expect_false(is.null(lng))
  expect_true(any(is.finite(unlist(lng$vpc_t))))
})

test_that("pcv_importance() never carried the na.action defect", {
  skip_on_cran()
  # Scope correction found in this pass's second self-check: pcv_importance()
  # routes through maihda_refit_draw() but its shared preamble reduces to ONE
  # complete-case analytic sample, so its subset models carry no na.action and
  # the tag is a no-op. Pinned so a future change to that preamble -- one that
  # let an NA-bearing frame reach the refit -- shows up here rather than as a
  # silently reintroduced failure.
  set.seed(2)
  n <- 600
  d <- data.frame(d1 = factor(sample(c("a", "b"), n, TRUE)),
                  d2 = factor(sample(c("p", "q"), n, TRUE)),
                  d3 = factor(sample(c("u", "v"), n, TRUE)))
  sk <- interaction(d$d1, d$d2, d$d3, drop = TRUE)
  d$y <- 1 + stats::rnorm(nlevels(sk), 0, 0.6)[as.integer(sk)] + stats::rnorm(n)
  d$y[sample(n, 40)] <- NA
  ds <- make_strata(d, c("d1", "d2", "d3"))$data

  seen <- integer(0)
  res <- testthat::with_mocked_bindings(
    suppressWarnings(
      pcv_importance(ds, "y", c("d1", "d2", "d3"), bootstrap = TRUE, n_boot = 12)),
    maihda_refit_draw = function(model, newresp) {
      seen <<- c(seen, as.integer(!is.null(attr(model@frame, "na.action"))))
      lme4::refit(model, newresp = newresp)   # the PRE-FIX call
    },
    .package = "MAIHDA")
  expect_gt(length(seen), 0L)
  expect_identical(sum(seen), 0L)             # no refit target carries an na.action
  expect_false(is.null(res))                  # and the untagged refit still works
})
