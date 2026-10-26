# Audit pass 2026-10-26: random-effect blocks whose grouping factors tie in number of levels.
#
# lme4::getME(fit, "mmList") and lme4::getME(fit, "cnms") put the blocks in different orders
# where lme4 sorts them by their number of levels (a block written after one with fewer) and
# two of them have as many levels: lme4 keeps such tied blocks in the formula's order and
# reformulas, which orders cnms, reverses them. The measures that re-express each block by its
# term columns paired the two lists by position, found them disagree and fell back silently.
# The standardised refit start was skipped, and fit_maihda() kept lme4's own stop of a fit
# with a slope on a birth year, far below the optimum and often recorded converged; the lmer
# re-check, the search records of singular lmer and glmer fits, the singular verdict in scaled
# terms and the trap record's Newton step were skipped; the glmer re-check was taken in lme4's
# own units, where a slope in large units gives none. The
# fixed-effect containment degrees of freedom paired a block's columns with another block's
# grouping factor. Each block's columns are now found by its grouping factor and terms.
#
# Random numbers are drawn only inside with_seed(), which puts the suite's stream back.

with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  force(expr)
}

quiet <- function(e) suppressWarnings(suppressMessages(e))

stay_lmer <- function(fn, par, lower, upper, control) {
  list(par = par, fval = fn(par), conv = 0, feval = 1)
}

# Two grouping factors of 8 levels behind one of 30, a slope on x in g1's block, a covariate z
# constant within g1, and a Poisson count with a slope on a birth year in g1's block.
tied_data <- function(seed = 11, n = 1500) {
  with_seed(seed, {
    d <- data.frame(g1 = factor(sample(8, n, TRUE)), g2 = factor(sample(8, n, TRUE)),
                    gbig = factor(sample(30, n, TRUE)), x = round(stats::rnorm(n), 2))
    zlev <- stats::rnorm(8)
    d$z <- zlev[d$g1]
    d$yr <- round(stats::runif(n, 1940, 2000))
    u1 <- stats::rnorm(8, 0, 0.5)
    v1 <- stats::rnorm(8, 0, 0.3)
    u2 <- stats::rnorm(8, 0, 0.4)
    ub <- stats::rnorm(30, 0, 0.3)
    d$y <- 1 + 0.5 * d$x + 0.3 * d$z + u1[d$g1] + v1[d$g1] * d$x + u2[d$g2] + ub[d$gbig] +
      stats::rnorm(n)
    d$cnt <- stats::rpois(n, exp(0.2 + 0.2 * d$z + 0.4 * u1[d$g1] + 0.004 * (d$yr - 1970) +
                                   stats::rnorm(8, 0, 0.003)[d$g1] * (d$yr - 1970) +
                                   0.3 * u2[d$g2] + 0.3 * ub[d$gbig]))
    d
  })
}

# The same model written with the tied blocks in the formula's order and in the order cnms has.
tied <- y ~ x + z + (1 + x | g1) + (1 | g2) + (1 | gbig)
sorted <- y ~ x + z + (1 | gbig) + (1 | g2) + (1 + x | g1)

test_that("a block's columns are found by its grouping factor and terms, not its position", {
  cnms <- list(gbig = "(Intercept)", g2 = "(Intercept)", g1 = c("(Intercept)", "x"))
  one <- matrix(1, 3, 1, dimnames = list(NULL, "(Intercept)"))
  two <- cbind("(Intercept)" = 1, x = c(0.5, -1, 2))
  big <- one * 2
  # g2's intercept block listed before gbig's: the columns alone would pair them wrongly.
  values <- list("1 | g2" = one, "1 + x | g1" = two, "1 | gbig" = big)
  out <- maihda_lme4_match_blocks(values, cnms)
  expect_identical(names(out), names(cnms))
  expect_identical(unname(out), list(big, one, two))
  # Each matrix is taken once.
  expect_identical(unname(maihda_lme4_match_blocks(list("1 | g" = one, "1 | g" = big),
                                                   list(g = "(Intercept)", g = "(Intercept)"))),
                   list(one, big))
  # Two blocks on one factor, as (x || g) gives, are told apart by their terms.
  split <- list("1 | g" = one, "0 + x | g" = two[, "x", drop = FALSE])
  expect_identical(unname(maihda_lme4_match_blocks(split[2:1],
                                                   list(g = "(Intercept)", g = "x"))),
                   unname(split))
  # A factor whose name is not syntactic: cnms names it bare and an interaction of it with its
  # backticks, as deparse1() deparses the bar's right-hand side.
  odd <- list("1 | `my g`:h" = one, "1 + x | `my g`" = two)
  expect_identical(maihda_lme4_match_blocks(odd, list("my g" = c("(Intercept)", "x"),
                                                      "`my g`:h" = "(Intercept)")),
                   list("my g" = two, "`my g`:h" = one))
  # Unnamed values are taken in their own order; a block with nothing to match finds none.
  expect_identical(maihda_lme4_match_blocks(unname(values), cnms), unname(values))
  expect_null(maihda_lme4_match_blocks(values[-1], cnms))
  expect_null(maihda_lme4_match_blocks(c(values[1:2], "1 | g3" = list(one)), cnms))
})

test_that("each measure gives the tied spelling of a model what the sorted spelling gets", {
  skip_on_cran()
  d <- tied_data()
  m <- quiet(lme4::lmer(tied, data = d))
  s <- quiet(lme4::lmer(sorted, data = d))
  expect_equal(lme4::REMLcrit(m), lme4::REMLcrit(s), tolerance = 1e-8)
  # The lmer re-check, the search's re-expression and the standardised start.
  expect_false(is.na(maihda_lmer_stationary_distance(m)))
  expect_equal(maihda_lmer_stationary_distance(m), maihda_lmer_stationary_distance(s),
               tolerance = 1e-6)
  expect_false(is.null(maihda_lmer_reexpressed(m)))
  expect_length(maihda_lmer_standardised_start(m), 1L)
  # The singular verdict in centred and scaled terms.
  expect_false(is.null(maihda_re_boundary_components_scaled_lme4(m)))
  expect_identical(maihda_re_boundary_components_scaled_lme4(m),
                   maihda_re_boundary_components_scaled_lme4(s))
  # The containment degrees of freedom: x and z are contained in g1's block (8 levels, three
  # columns absorbed), not in g2's.
  expect_identical(maihda_containment_df(m), maihda_containment_df(s))
  expect_identical(unname(maihda_containment_df(m)), c(5, 5, 5))
  # lme4's two orders of the blocks disagree here; the columns follow cnms.
  expect_false(identical(sub(".*\\| ", "", names(lme4::getME(m, "mmList"))),
                         names(lme4::getME(m, "cnms"))))
  cols <- maihda_lme4_block_columns(m)
  expect_identical(names(cols), names(lme4::getME(m, "cnms")))
  expect_identical(lapply(cols, colnames), lapply(lme4::getME(m, "cnms"), unname))
  expect_identical(maihda_lmer_block_scales(cols, lme4::getME(m, "cnms")),
                   maihda_lmer_block_scales(maihda_lme4_block_columns(s), lme4::getME(s, "cnms")))
})

# A slope on a birth year (1940 to 2000) in the stratum's block, beside a factor with as many
# levels and one with more: lme4's own fit stops 22 log-likelihood units below the optimum,
# which the refit from the standardised start reaches.
birth_year_data <- function() {
  with_seed(6, {
    n <- 1500
    d <- data.frame(stratum = factor(sample(8, n, TRUE)), g2 = factor(sample(8, n, TRUE)),
                    gbig = factor(sample(30, n, TRUE)), yr = round(stats::runif(n, 1940, 2000)))
    s <- as.integer(d$stratum)
    d$y <- 1 + 0.01 * (d$yr - 1970) + stats::rnorm(8, 0, 0.3)[s] +
      stats::rnorm(8, 0, 0.01)[s] * (d$yr - 1970) + stats::rnorm(8, 0, 0.3)[d$g2] +
      stats::rnorm(30, 0, 0.3)[d$gbig] + stats::rnorm(n)
    d
  })
}

test_that("fit_maihda() refits a tied birth-year slope fit to the optimum of the sorted spelling", {
  skip_on_cran()
  d <- birth_year_data()
  tied_fit <- quiet(fit_maihda(y ~ yr + (1 + yr | stratum) + (1 | g2) + (1 | gbig), data = d))
  sorted_fit <- quiet(fit_maihda(y ~ yr + (1 | gbig) + (1 | g2) + (1 + yr | stratum), data = d))
  own <- quiet(lme4::lmer(y ~ yr + (1 + yr | stratum) + (1 | g2) + (1 | gbig), data = d))
  ll <- function(m) as.numeric(stats::logLik(m))
  expect_equal(ll(tied_fit$model), ll(sorted_fit$model), tolerance = 1e-3)
  expect_gt(ll(tied_fit$model) - ll(own), 20)
  expect_true(tied_fit$diagnostics$converged)
})

# Age in years (20 to 70) and a binary z with stratum slopes in g1's block, the two 8-level
# factors behind the 30-level one: lme4's own fit of the tied spelling stops with g1's
# intercept and age columns at their bounds; across the bound the age column's criterion
# rises at every step, and only the Newton step shows it falling.
age_tied_data <- function() {
  with_seed(105, {
    n <- 1500
    d <- data.frame(g1 = factor(sample(8, n, TRUE)), g2 = factor(sample(8, n, TRUE)),
                    gbig = factor(sample(30, n, TRUE)), age = round(stats::runif(n, 20, 70), 1),
                    z = stats::rbinom(n, 1, 0.5))
    ua <- stats::rnorm(8, 0, 0.02)
    uz <- stats::rnorm(8, 0, 0.3)
    d$y <- 1 + 0.02 * d$age + 0.3 * d$z + ua[d$g1] * d$age + uz[d$g1] * d$z +
      stats::rnorm(8, 0, 0.3)[d$g2] + stats::rnorm(30, 0, 0.3)[d$gbig] + stats::rnorm(n)
    d
  })
}
trap_age_tied <- c(0.2867115, 0.2696478, 0, 0.01840528, 0.0556246, 0, -0.07930906, 0.202071)

test_that("the trap record takes its Newton step in a tied block", {
  skip_on_cran()
  d <- age_tied_data()
  hold <- lme4::lmerControl(optimizer = stay_lmer, calc.derivs = FALSE, boundary.tol = 0)
  tied_fit <- quiet(lme4::lmer(y ~ age + z + (1 + age + z | g1) + (1 | g2) + (1 | gbig), data = d,
                               start = trap_age_tied, control = hold))
  sorted_fit <- quiet(lme4::lmer(y ~ age + z + (1 | gbig) + (1 | g2) + (1 + age + z | g1), data = d,
                                 start = trap_age_tied, control = hold))
  expect_identical(unname(lme4::getME(tied_fit, "theta")), trap_age_tied)
  tr <- maihda_lme4_boundary_traps(tied_fit)
  expect_identical(tr, maihda_lme4_boundary_traps(sorted_fit))
  age <- tr[tr$term == "age", ]
  expect_lt(age$across, 0)
  expect_gt(age$drop, 0.5)
})

test_that("a singular tied fit is searched as the sorted spelling is", {
  skip_on_cran()
  d <- tied_data()
  m <- quiet(lme4::lmer(tied, data = d))
  trap <- unname(lme4::getME(m, "theta"))
  trap[3] <- 0
  held <- quiet(lme4::lmer(tied, data = d, start = trap,
                           control = lme4::lmerControl(optimizer = stay_lmer, calc.derivs = FALSE)))
  held_sorted <- quiet(lme4::lmer(sorted, data = d, start = trap,
                                  control = lme4::lmerControl(optimizer = stay_lmer,
                                                              calc.derivs = FALSE)))
  expect_true(lme4::isSingular(held))
  gain <- maihda_lmer_search_gain(held)
  expect_false(is.na(gain))
  expect_equal(gain, maihda_lmer_search_gain(held_sorted), tolerance = 1e-6)
  expect_identical(maihda_lme4_boundary_traps(held), maihda_lme4_boundary_traps(held_sorted))
})

test_that("a tied glmer fit with a birth-year slope is re-checked as the sorted spelling is", {
  skip_on_cran()
  d <- tied_data(n = 600)
  g <- quiet(lme4::glmer(cnt ~ yr + z + (1 + yr | g1) + (1 | g2) + (1 | gbig), data = d,
                         family = stats::poisson))
  gs <- quiet(lme4::glmer(cnt ~ yr + z + (1 | gbig) + (1 | g2) + (1 + yr | g1), data = d,
                          family = stats::poisson))
  expect_false(is.na(maihda_lme4_stationary_distance(g)))
  expect_equal(maihda_lme4_stationary_distance(g), maihda_lme4_stationary_distance(gs),
               tolerance = 1e-6)
  # The singular verdict in centred and scaled terms, and the search of a fit lme4 calls
  # singular.
  expect_false(is.null(maihda_re_boundary_components_scaled_lme4(g)))
  expect_identical(maihda_re_boundary_components_scaled_lme4(g),
                   maihda_re_boundary_components_scaled_lme4(gs))
  expect_false(is.na(maihda_glmer_search_gain(g)))
  expect_equal(maihda_glmer_search_gain(g), maihda_glmer_search_gain(gs), tolerance = 1e-6)
})

test_that("blocks lme4 orders alike are taken as they come, and (x || g) ties with itself", {
  skip_on_cran()
  d <- tied_data()
  # No tie: the reordered blocks agree, and an already sorted formula is not reordered.
  for (f in list(y ~ x + (1 + x | g1) + (1 | gbig), y ~ x + (1 | gbig) + (1 + x | g1))) {
    m <- quiet(lme4::lmer(f, data = d))
    expect_identical(unname(maihda_lme4_block_columns(m)), unname(lme4::getME(m, "mmList")))
  }
  # (x || g1) is two blocks on g1, sorted behind gbig's: lme4 lists (1 | g1) before
  # (0 + x | g1), cnms the other way round.
  m <- quiet(lme4::lmer(y ~ x + (x || g1) + (1 | gbig), data = d))
  cnms <- lme4::getME(m, "cnms")
  expect_null(maihda_lmer_block_scales(lme4::getME(m, "mmList"), cnms))
  cols <- maihda_lme4_block_columns(m)
  expect_identical(lapply(cols, colnames), lapply(cnms, unname))
  expect_false(is.null(maihda_lmer_block_scales(cols, cnms)))
  # A tied factor whose name is not syntactic.
  d$`my g` <- d$g1
  m <- quiet(lme4::lmer(y ~ x + (1 + x | `my g`) + (1 | g2) + (1 | gbig), data = d))
  cnms <- lme4::getME(m, "cnms")
  cols <- maihda_lme4_block_columns(m)
  expect_identical(names(cols), names(cnms))
  expect_identical(lapply(cols, colnames), lapply(cnms, unname))
})
