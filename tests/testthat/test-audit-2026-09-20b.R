# Audit 2026-09-20 (b): the fixed-effect bootstrap constrained an entire TERM
# while reporting a p-value and interval for each of that term's COEFFICIENTS.
#
# For a three-level factor `f` the fitted design carries two columns, `fb` and
# `fc`. maihda_bootstrap_fixef() removed BOTH, simulated from `fb = fc = 0`, and
# then referred each coefficient's own Wald statistic to those draws. The row
# labelled `fc` states `fc = 0` with `fb` free; the draws came from a model in
# which `fb` had been deleted as well. Measured on the reported design (binomial,
# 12 strata of 60, beta_fb = 3, beta_fc = 0.4):
#
#   full fit, either spelling   fb 2.5576353  fc 0.5791719 (se 0.4033327)
#                               between-stratum variance 0.2402725
#   null the package used       (Intercept) only, tau^2 1.4584371
#   null the row claims         (Intercept) + fb = 2.27323, tau^2 0.2990457
#
# The deleted sibling lands in the stratum variance -- a factor of 4.9 here --
# so the reference distribution is built from data the fit does not describe. Nor
# is it stable: `y ~ f` and `y ~ fb + fc` are the SAME fitted model (identical
# fixef, identical logLik, identical between-stratum variance above), and the
# second spelling puts each coefficient in a term of its own, so it already got
# the null the row claims. On maihda_f4_data() below at B = 99, HEAD gave the
# `fc` contrast p = 0.13, [-0.2103895, 1.2477903] spelled as a factor and
# p = 0.14, [-0.1116736, 1.1490744] spelled as indicators.
#
# How much that moves a p-value is a separate question from which hypothesis is
# being tested, and the honest answer is that it depends on the design. The
# statistic is studentized, so it is close to pivotal: on the binomial design
# above, 1000 draws from each null gave 95% critical values of 2.5548 (MC se
# 0.126) and 2.6242 (MC se 0.118) -- indistinguishable, although the means and
# SDs of |t*| were 2.46/20.0 against 1.47/10.7, the wrong null's far tail being
# much heavier. On a Gaussian design of the same shape (12 strata of 20,
# beta_fb = 3, tau = 0.25; 250 replicates x 99 draws) the wrong null never
# reached the variance boundary where the right one did 8% of the time; its mean
# 95% critical value was 2.3727 against 2.1748 (paired over the same datasets,
# so that 9% inflation is measured well), and it rejected a true null at 3.2%
# (MC se 1.1) against 4.4% (1.3) for a nominal 5%, the two disagreeing on 2.8%
# of datasets; at a true beta_fc = 0.35 it rejected at 27.6% (2.8) against 30.0%
# (2.9), disagreeing on 9.6%. So the practical effect runs from negligible to
# systematically conservative; the reason to fix it is that the row's hypothesis
# and the row's reference distribution have to be the same one.
#
# CONFIRMED. The fix restricts one design COLUMN at a time -- the siblings are
# nuisance parameters, kept and re-estimated -- which is the restriction the
# coefficient names under whatever contrasts the model was fitted with.
# maihda_restrict_fixef() already built exactly that fit; the loop never asked
# for it. n_boot is now spent per tested coefficient rather than per term.
#
# CONFINED: a model whose terms each span one design column is untouched, since
# there the two restrictions coincide and the formula route is still taken first.
# Measured identical() against a HEAD worktree for the indicator spelling and for
# a `y ~ x + h` fit; only the factor spelling moved, and it moved onto the
# indicator spelling's own numbers -- exactly, at a matching column order. What
# holds for every order is the null; the draws are taken one coefficient at a
# time, so reordering the columns reorders the blocks and gives another Monte
# Carlo realisation of that same null. Every pre-existing bootstrap test uses
# two-level factors or continuous covariates, which is why this went unseen.

maihda_f4_data <- function(seed = 4, n_per = 25) {
  set.seed(seed)
  dg <- expand.grid(f = factor(c("a", "b", "c")), g = factor(1:4))
  d <- dg[rep(seq_len(nrow(dg)), each = n_per), ]
  d$stratum <- interaction(d$f, d$g)
  # Named so that the indicator spelling reproduces the factor spelling's
  # coefficient names AND its column order: (Intercept), fb, fc.
  d$fb <- as.numeric(d$f == "b")
  d$fc <- as.numeric(d$f == "c")
  u <- stats::rnorm(nlevels(d$stratum), 0, 0.4)
  # A large sibling effect on b and a small one on c: deleting b is what inflates
  # the stratum variance of the null the package used to simulate from.
  d$y <- 2 * d$fb + 0.3 * d$fc + u[as.integer(d$stratum)] + stats::rnorm(nrow(d))
  rownames(d) <- NULL
  d
}

# The loop's own order -- terms outermost, columns within a term -- and the RNG
# it consumes: nothing before the first simulate() call, then one block of
# n_boot draws per tested COLUMN.
maihda_f4_replay <- function(m, B, restrict, conf_level = 0.95) {
  X <- lme4::getME(m, "X")
  asg <- attr(X, "assign")
  tl <- attr(stats::terms(stats::formula(m, fixed.only = TRUE)), "term.labels")
  est <- unname(lme4::fixef(m))
  se <- unname(sqrt(diag(as.matrix(stats::vcov(m)))))
  p <- lo <- hi <- rep(NA_real_, length(est))
  for (k in seq_along(tl)) {
    cols <- which(asg == k)
    for (j in cols) {
      red <- restrict(m, j, cols, tl[k])
      sim <- maihda_simulate_lme4(red, nsim = B)
      ts <- vapply(seq_len(B), function(i) {
        bm <- lme4::refit(m, newresp = sim[[i]])
        abs(lme4::fixef(bm)[[j]] / sqrt(diag(as.matrix(stats::vcov(bm))))[[j]])
      }, numeric(1))
      p[j] <- min(1, (1 + sum(ts >= abs(est[j] / se[j]))) / (B + 1))
      crit <- sort(ts, decreasing = TRUE)[maihda_boot_critical_rank(B, conf_level)]
      lo[j] <- est[j] - crit * se[j]
      hi[j] <- est[j] + crit * se[j]
    }
  }
  list(p_value = p, lower = lo, upper = hi)
}

# --------------------------------------------------------------- the mechanism

test_that("a coefficient's null keeps the other columns of its own term", {
  skip_on_cran()
  d <- maihda_f4_data()
  m <- lme4::lmer(y ~ f + (1 | stratum), data = d, REML = TRUE)
  X <- lme4::getME(m, "X")
  expect_identical(colnames(X), c("(Intercept)", "fb", "fc"))
  expect_identical(attr(X, "assign"), c(0L, 1L, 1L))   # ONE term, two columns

  j_fc <- match("fc", colnames(X))
  per_coef <- maihda_restrict_fixef(m, j_fc)
  per_term <- maihda_restrict_fixef(m, which(attr(X, "assign") == 1L))

  # The coefficient's own null still carries its sibling, and its variance
  # components are those of the model the row describes -- pinned against an
  # independently spelled fit, so this does not merely restate the helper.
  sib <- lme4::lmer(y ~ fb + (1 | stratum), data = d, REML = TRUE)
  expect_length(lme4::fixef(per_coef), 2L)
  expect_equal(unname(lme4::fixef(per_coef)), unname(lme4::fixef(sib)))
  expect_equal(as.numeric(lme4::VarCorr(per_coef)$stratum),
               as.numeric(lme4::VarCorr(sib)$stratum))

  # The term's null deletes the sibling, and the sibling's effect reappears as
  # between-stratum variance. That was the reference HEAD used for BOTH rows.
  expect_length(lme4::fixef(per_term), 1L)
  expect_gt(as.numeric(lme4::VarCorr(per_term)$stratum),
            3 * as.numeric(lme4::VarCorr(per_coef)$stratum))
})

# ------------------------------------------------------------------ the fix

test_that("the reference is built from the model without THAT COEFFICIENT", {
  skip_on_cran()
  d <- maihda_f4_data()
  m <- lme4::lmer(y ~ f + (1 | stratum), data = d, REML = TRUE)
  B <- 39L

  set.seed(88)
  fx <- suppressWarnings(maihda_bootstrap_fixef(m, n_boot = B, conf_level = 0.95))
  set.seed(88)
  want <- maihda_f4_replay(m, B, function(m, j, cols, lab)
    maihda_restrict_fixef(m, j))
  set.seed(88)
  head_way <- maihda_f4_replay(m, B, function(m, j, cols, lab)
    maihda_restrict_fixef(m, cols))

  tested <- 2:3
  # The replay reads its critical value at rank maihda_boot_critical_rank(B),
  # where the shipped code uses the count of finite draws. They agree only while
  # every refit succeeds, so say so: a future mismatch should point at a failed
  # refit rather than at the restriction.
  expect_identical(attr(fx, "n_boot_ok")[tested], c(B, B))
  expect_equal(fx$p_value[tested], want$p_value[tested])
  expect_equal(fx$lower[tested], want$lower[tested])
  expect_equal(fx$upper[tested], want$upper[tested])
  # ... and NOT what constraining the whole term gives, or the test has no teeth.
  expect_false(isTRUE(all.equal(fx$upper[tested], head_way$upper[tested])))
  # The intercept still has no reduced model to simulate from.
  expect_true(is.na(fx$p_value[1]) && is.na(fx$lower[1]) && is.na(fx$upper[1]))
  # Estimate and standard error are untouched: only the reference moved.
  expect_equal(fx$estimate, unname(lme4::fixef(m)))
  expect_equal(fx$se, unname(sqrt(diag(as.matrix(stats::vcov(m))))))
})

test_that("a factor term and its indicator spelling now share the null", {
  skip_on_cran()
  d <- maihda_f4_data()
  mf <- lme4::lmer(y ~ f + (1 | stratum), data = d, REML = TRUE)
  mi <- lme4::lmer(y ~ fb + fc + (1 | stratum), data = d, REML = TRUE)
  # The same fitted model, reached two ways: one term of two columns against two
  # terms of one. Same columns in the same order, so the draws line up.
  expect_equal(unname(lme4::fixef(mf)), unname(lme4::fixef(mi)))
  expect_equal(as.numeric(lme4::VarCorr(mf)$stratum),
               as.numeric(lme4::VarCorr(mi)$stratum))
  expect_identical(attr(lme4::getME(mf, "X"), "assign"), c(0L, 1L, 1L))
  expect_identical(attr(lme4::getME(mi, "X"), "assign"), c(0L, 1L, 2L))

  B <- 39L
  set.seed(31); a <- suppressWarnings(maihda_bootstrap_fixef(mf, B, 0.95))
  set.seed(31); b <- suppressWarnings(maihda_bootstrap_fixef(mi, B, 0.95))
  expect_equal(a$p_value, b$p_value)
  expect_equal(a$lower, b$lower)
  expect_equal(a$upper, b$upper)

  # What holds for EVERY spelling is the restriction; the numbers coincide above
  # because the columns are also in the same order, and the blocks are drawn one
  # coefficient at a time. Spelled y ~ fc + fb -- the order the report used --
  # the `fc` block is drawn first instead, which is a different Monte Carlo
  # realisation of the SAME null: at this B and seed its `fc` interval is
  # [-0.3193886, 1.3567890] where the two above give [-0.6258361, 1.6632370],
  # both at p = 0.15. So pin the null, which is what the fix makes invariant.
  mr <- lme4::lmer(y ~ fc + fb + (1 | stratum), data = d, REML = TRUE)
  tau2_null <- function(m, nm) {
    X <- lme4::getME(m, "X")
    red <- maihda_restrict_fixef(m, match(nm, colnames(X)))
    c(as.numeric(lme4::VarCorr(red)$stratum), length(lme4::fixef(red)))
  }
  expect_equal(tau2_null(mf, "fc"), tau2_null(mi, "fc"))
  expect_equal(tau2_null(mf, "fc"), tau2_null(mr, "fc"))
})

test_that("the restriction is the one the coefficient names under its contrasts", {
  skip_on_cran()
  d <- maihda_f4_data()
  # Sum coding, where a coefficient is a level's deviation from the mean of the
  # level means rather than a difference from a reference. The help page says the
  # restriction follows the coding, so check the identity that claim implies.
  # (A control: maihda_restrict_fixef() itself is unchanged by this pass.)
  m <- lme4::lmer(y ~ f + (1 | stratum), data = d, REML = TRUE,
                  contrasts = list(f = stats::contr.sum))
  X <- lme4::getME(m, "X")
  expect_identical(colnames(X), c("(Intercept)", "f1", "f2"))
  expect_identical(attr(X, "assign"), c(0L, 1L, 1L))

  b <- lme4::fixef(maihda_restrict_fixef(m, match("f1", colnames(X))))
  # With f1 constrained the three level means are mu, mu + f2 and mu - f2, so
  # level 1's mean IS their mean and its deviation is exactly zero.
  means <- c(b[[1]], b[[1]] + b[[2]], b[[1]] - b[[2]])
  expect_equal(means[1] - mean(means), 0)

  # The bootstrap runs on a non-default contrast at all: the null is built from
  # the fitted design columns, which carry whatever coding was used.
  set.seed(3)
  fx <- suppressWarnings(maihda_bootstrap_fixef(m, n_boot = 19L, conf_level = 0.90))
  expect_true(all(is.finite(fx$p_value[2:3])))
  expect_true(all(is.finite(fx$lower[2:3])))
})

test_that("an interaction's columns are restricted one at a time", {
  skip_on_cran()
  d <- maihda_f4_data()
  set.seed(12)
  d$x <- stats::rnorm(nrow(d))
  d$y <- d$y + 0.5 * d$x + 0.4 * d$x * d$fb
  m <- lme4::lmer(y ~ x * f + (1 | stratum), data = d, REML = TRUE)
  X <- lme4::getME(m, "X")
  asg <- attr(X, "assign")
  # x:f spans two columns, and so does f.
  expect_identical(sum(asg == 3L), 2L)

  j <- which(asg == 3L)[1]
  red <- maihda_restrict_fixef(m, j)
  # One column gone, not the pair: the sibling interaction column survives.
  expect_length(lme4::fixef(red), ncol(X) - 1L)

  B <- 25L
  set.seed(6)
  fx <- suppressWarnings(maihda_bootstrap_fixef(m, n_boot = B, conf_level = 0.95))
  set.seed(6)
  want <- maihda_f4_replay(m, B, function(m, j, cols, lab)
    maihda_restrict_fixef(m, j))
  set.seed(6)
  head_way <- maihda_f4_replay(m, B, function(m, j, cols, lab)
    maihda_restrict_fixef(m, cols))
  tested <- which(asg > 0L)
  # The p-values alone do not discriminate here: at B = 25 four of the five sit
  # on the floor 1 / 26 and the fifth agrees by coincidence, so they are the same
  # under either restriction. The interval endpoints are continuous and do
  # separate them -- assert on those, and on the contrast with the term-wide
  # restriction, or this block would pass on the unfixed code.
  expect_equal(fx$p_value[tested], want$p_value[tested])
  expect_equal(fx$lower[tested], want$lower[tested])
  expect_equal(fx$upper[tested], want$upper[tested])
  expect_false(isTRUE(all.equal(fx$lower[tested], head_way$lower[tested])))
  # The one-column term `x` comes first, so its block is the same draws either
  # way: only the multi-column terms move.
  jx <- which(asg == 1L)
  expect_equal(fx$lower[jx], head_way$lower[jx])
})

# ------------------------------------------------------------- confinement

test_that("a one-column term still takes the formula reduction, unchanged", {
  skip_on_cran()
  set.seed(9)
  n <- 400
  d <- data.frame(x = stats::rnorm(n),
                  h = factor(sample(c("lo", "hi"), n, TRUE)),
                  st = factor(sample(paste0("s", 1:8), n, TRUE)))
  d$y <- 0.5 * d$x + 0.4 * (d$h == "hi") +
    stats::rnorm(8, 0, 0.5)[as.integer(d$st)] + stats::rnorm(n)
  m <- lme4::lmer(y ~ x + h + (1 | st), data = d, REML = TRUE)
  expect_identical(attr(lme4::getME(m, "X"), "assign"), c(0L, 1L, 2L))

  B <- 39L
  set.seed(5)
  fx <- suppressWarnings(maihda_bootstrap_fixef(m, n_boot = B, conf_level = 0.95))
  # Replayed through the FORMULA, not the design matrix: where a term is one
  # column the two are the same restriction and the formula route is still first.
  set.seed(5)
  want <- maihda_f4_replay(m, B, function(m, j, cols, lab)
    stats::update(m, formula. = stats::update(stats::formula(m),
                                              paste(". ~ . -", lab))))
  expect_equal(fx$p_value[2:3], want$p_value[2:3])
  expect_equal(fx$lower[2:3], want$lower[2:3])
  expect_equal(fx$upper[2:3], want$upper[2:3])
})

test_that("the binomial case from the report is spelling-invariant too", {
  skip_on_cran()
  set.seed(19)
  dg <- expand.grid(f = factor(c("a", "b", "c")), g = factor(1:4))
  d <- dg[rep(seq_len(nrow(dg)), each = 30), ]
  d$stratum <- interaction(d$f, d$g)
  d$fb <- as.numeric(d$f == "b")
  d$fc <- as.numeric(d$f == "c")
  u <- stats::rnorm(nlevels(d$stratum), 0, 0.45)
  d$y <- stats::rbinom(nrow(d), 1,
                       stats::plogis(-1 + 2.5 * d$fb + 0.4 * d$fc +
                                       u[as.integer(d$stratum)]))
  mf <- suppressWarnings(lme4::glmer(y ~ f + (1 | stratum), data = d,
                                     family = stats::binomial()))
  mi <- suppressWarnings(lme4::glmer(y ~ fb + fc + (1 | stratum), data = d,
                                     family = stats::binomial()))
  expect_equal(unname(lme4::fixef(mf)), unname(lme4::fixef(mi)))

  # A glmer refit is not cheap, so few draws -- and at 90% rather than 95%, where
  # floor(0.10 * (B + 1)) is still 2 and the critical value stays finite.
  B <- 19L
  set.seed(77); a <- suppressWarnings(maihda_bootstrap_fixef(mf, B, 0.90))
  set.seed(77); b <- suppressWarnings(maihda_bootstrap_fixef(mi, B, 0.90))
  expect_true(all(is.finite(a$lower[2:3])))
  expect_equal(a$p_value, b$p_value)
  expect_equal(a$lower, b$lower)
  expect_equal(a$upper, b$upper)
})

# ------------------------------------------------------------------- the docs

test_that("the help pages state the restriction per coefficient", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  rd <- function(f) {
    gsub("[[:space:]]+", " ",
         paste(readLines(file.path(man, f), warn = FALSE), collapse = " "))
  }
  sm <- rd("summary.maihda_model.Rd")
  expect_true(grepl("each fixed-effect coefficient the model is refitted with that coefficient",
                    sm, fixed = TRUE))
  expect_true(grepl("The restriction is on the coefficient, not on its term", sm,
                    fixed = TRUE))
  expect_true(grepl("the siblings kept and re-estimated", sm, fixed = TRUE))
  expect_true(grepl("refits \\emph{per tested coefficient}", sm, fixed = TRUE))
  # ... and the df_method argument's own line, where the phrase wraps.
  expect_true(grepl("refits \\emph{per tested fixed-effect coefficient}", sm,
                    fixed = TRUE))
  expect_false(grepl("refits \\emph{per term}", sm, fixed = TRUE))
  expect_false(grepl("refits \\emph{per fixed-effect term}", sm, fixed = TRUE))

  bf <- rd("maihda_bootstrap_fixef.Rd")
  expect_true(grepl("For each fixed-effect \\emph{coefficient}", bf, fixed = TRUE))
  expect_true(grepl("One coefficient at a time, not one term at a time", bf,
                    fixed = TRUE))
  expect_false(grepl("replicates \\emph{per term}", bf, fixed = TRUE))
})
