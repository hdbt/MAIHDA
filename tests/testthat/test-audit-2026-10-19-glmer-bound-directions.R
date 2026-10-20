# Audit pass 2026-10-19: glmer fits stopped at a variance bound in a random-effect block of
# three or more terms.
#
# Where a column of a block's Cholesky factor has its diagonal entry at the bound of 0,
# every re-factoring of the block's rows below that column gives the same covariance
# matrix, and leaving the bound moves that term's covariances with the later terms in a
# direction R q (R those rows, q a unit vector) that depends on the factoring. With one
# later term the directions are the two signs of their covariance; with two or more the
# reflection of the entries below (q = -e1) is one direction of many, and the likelihood
# can rise only in others. A glmer fit outside growth models stopped there is now refitted
# from across the bound along the direction in which the deviance falls fastest at first
# order, and the convergence record measures every such stop along that direction.
#
# Random numbers are drawn only inside with_seed(), which puts the suite's stream back.

with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  force(expr)
}

quiet <- function(e) suppressWarnings(suppressMessages(e))

# Twelve additive strata (A x B) with stratum slopes on x and z and no stratum effect at
# x = z = 0, so that the intercept variance of (1 + x + z | stratum) is small next to the
# slope variances.
three_term_design <- function(seed) {
  with_seed(seed, {
    n <- 1500
    A <- sample(0:2, n, TRUE)
    B <- sample(0:3, n, TRUE)
    st <- 4L * A + B + 1L
    x <- round(stats::rnorm(n), 2)
    z <- round(stats::rnorm(n), 2)
    ux <- stats::rnorm(12, 0, 0.2)
    uz <- stats::rnorm(12, 0, 0.15)
    eta <- 0.3 + 0.2 * x - 0.1 * z + c(0, 0.3, 0.5)[A + 1] + c(0, -0.2, 0.1, 0.3)[B + 1] +
      ux[st] * x + uz[st] * z
    data.frame(y = stats::rpois(n, exp(eta)), x = x, z = z, A = factor(A), B = factor(B),
               stratum = factor(st))
  })
}

three_term_formula <- y ~ x + z + A + B + (1 + x + z | stratum)

# lme4's own fits of seeds 3 and 22 stop with the intercept entry at 0. Seed 3: at
# -2483.04958, where a refit from the reflection stopped at -2483.01334, recorded
# converged; the maximum, from a search over the deviance function with no bound and from
# glmmTMB, is -2482.34487. Seed 22: at -2516.64742, the maximum -2516.25888; across the
# bound the deviance rises by 0.0028 at a step of 0.001 by the reflection and falls by 0.047
# along the first-order direction.
trap_3 <- c(0, 0.1521901, -0.06153039, 0.01991965, -0.03108059, 0.2094773)
fixef_3 <- c(0.2542837, 1.050656e-06, -0.164203, 0.3269555, 0.5209243, -0.1968458,
             0.1190303, 0.2784862)
trap_22 <- c(0, 0.01742157, 0.07183298, 0.2065725, -0.001465654, 0.1711376)
fixef_22 <- c(0.3466641, 0.1433912, -0.06572225, 0.2776074, 0.49189, -0.2455793,
              0.102721, 0.2643339)

# An optimiser for glmerControl() that keeps its start: a fit held at a given point.
stay <- function(par, fn, lower, upper, control = list(), ...) {
  list(par = par, fval = fn(par), conv = 0L, message = "start kept")
}

# The 136th's two-term design and lme4's stop on its seed 7 (see that pass's test file).
slope_design <- function(seed, n = 600) {
  with_seed(seed, {
    g <- expand.grid(A = 0:1, B = 0:1, C = 0:1)
    st <- sample(8, n, TRUE)
    x <- round(stats::rnorm(n), 2)
    u1 <- stats::rnorm(8, 0, 0.25)
    add <- 0.4 * (g$A - 0.5) - 0.3 * (g$B - 0.5) + 0.2 * (g$C - 0.5)
    eta <- 0.5 + 0.2 * x + add[st] + u1[st] * x
    data.frame(y = stats::rpois(n, exp(eta)), x = x, A = factor(g$A[st]), B = factor(g$B[st]),
               C = factor(g$C[st]), stratum = factor(st))
  })
}
trap_7 <- c(0, 0.2119, 0.0101)

# The covariance matrix of each block of a theta.
block_covs <- function(theta, cnms) {
  sizes <- lengths(cnms) * (lengths(cnms) + 1) / 2
  lapply(seq_along(cnms), function(b) {
    n <- length(cnms[[b]])
    L <- matrix(0, n, n)
    L[lower.tri(L, diag = TRUE)] <- theta[sum(sizes[seq_len(b - 1L)]) + seq_len(sizes[b])]
    tcrossprod(L)
  })
}

# ---- across the bound ---------------------------------------------------------------------

test_that("a column at a zero diagonal is re-factored without moving the covariance", {
  cnms <- list(id = c("(Intercept)", "t", "t2"), stratum = c("(Intercept)", "x", "z", "w"))
  theta <- c(0, 0.3, -0.2, 0.5, 0.6, 0.7, 0, 0.2, -0.1, 0.3, 0.4, 0.05, -0.2, 0.6, 0.1, 0.5)
  qs <- with_seed(5, matrix(stats::rnorm(40), 10))
  for (index in c(1L, 7L)) {
    k <- if (index == 1L) 3L else 4L
    for (r in 1:10) {
      q <- qs[r, seq_len(k)] / sqrt(sum(qs[r, seq_len(k)]^2))
      expect_equal(block_covs(maihda_lme4_refactor_theta(theta, cnms, index, q, 0), cnms),
                   block_covs(theta, cnms), tolerance = 1e-12)
    }
    # q = (-1, 0, ...) is the reflection.
    expect_equal(maihda_lme4_refactor_theta(theta, cnms, index, c(-1, rep(0, k - 1L)), 0.05),
                 maihda_lme4_mirror_theta(theta, cnms, index, 0.05), tolerance = 1e-12)
  }
  # Off the bound the term's covariances with the later terms are the diagonal entry times
  # R q, and those among the later terms stay; the other block is left as it is.
  q <- c(0.6, 0, -0.8, 0)
  across <- maihda_lme4_refactor_theta(theta, cnms, 7L, q, 0.05)
  S <- block_covs(across, cnms)[[2]]
  L <- matrix(0, 4, 4)
  L[lower.tri(L, diag = TRUE)] <- theta[7:16]
  expect_equal(S[1, 1], 0.05^2)
  expect_equal(S[2:4, 1], 0.05 * drop(L[2:4, ] %*% q), tolerance = 1e-12)
  expect_equal(S[2:4, 2:4], tcrossprod(L)[2:4, 2:4], tolerance = 1e-12)
  expect_identical(across[1:6], theta[1:6])
})

test_that("a positive semi-definite matrix gets a lower-triangular factor", {
  P <- matrix(c(4, 2, 2, 3), 2)
  expect_equal(maihda_lower_chol_psd(P), t(chol(P)), tolerance = 1e-12)
  # Of rank one: the second pivot is 0 and its column stays 0.
  S <- matrix(c(4, 2, 2, 1), 2)
  M <- maihda_lower_chol_psd(S)
  expect_equal(tcrossprod(M), S, tolerance = 1e-12)
  expect_identical(M[, 2], c(0, 0))
  # A zero first pivot.
  expect_equal(maihda_lower_chol_psd(matrix(c(0, 0, 0, 2), 2)), matrix(c(0, 0, 0, sqrt(2)), 2))
})

test_that("the first-order direction is the one in which the deviance falls fastest", {
  cnms <- list(stratum = c("(Intercept)", "x", "z", "w"))
  theta <- c(0, 0.2, -0.1, 0.3, 0.4, 0.05, -0.2, 0.6, 0.1, 0.5)
  L <- matrix(0, 4, 4)
  L[lower.tri(L, diag = TRUE)] <- theta
  R <- L[2:4, ]
  # A deviance linear in the intercept's covariances with the later terms, with gradient g
  # in them, plus a term in its variance.
  g <- c(1.5, -2, 0.7)
  dev <- function(th) {
    S <- block_covs(th, cnms)[[1]]
    sum(g * S[2:4, 1]) + 3 * S[1, 1]
  }
  across <- maihda_lme4_first_order_theta(theta, cnms, 1L, 0.05, dev)
  M <- matrix(0, 4, 4)
  M[lower.tri(M, diag = TRUE)] <- across
  h <- drop(crossprod(R, g))
  expect_equal(M[2:4, 1], -drop(R %*% h) / sqrt(sum(h^2)), tolerance = 1e-6)
  expect_identical(M[1, 1], 0.05)
  expect_equal(block_covs(replace(across, 1L, 0), cnms), block_covs(theta, cnms), tolerance = 1e-12)
  # A fraction of the length scales the column and leaves the later block non-singular.
  short <- maihda_lme4_first_order_theta(theta, cnms, 1L, 0.05, dev, fraction = 0.9)
  M9 <- matrix(0, 4, 4)
  M9[lower.tri(M9, diag = TRUE)] <- short
  expect_equal(M9[2:4, 1], 0.9 * M[2:4, 1], tolerance = 1e-12)
  expect_equal(block_covs(replace(short, 1L, 0), cnms), block_covs(theta, cnms), tolerance = 1e-12)
  expect_gt(min(diag(M9)[2:4]), 1e-3)
  expect_lt(min(abs(diag(M)[2:4])), 1e-8)
  # Along it the deviance falls at a step of 1e-3 by more than the record counts; by the
  # reflection it does not.
  expect_gt(dev(theta) - dev(maihda_lme4_first_order_theta(theta, cnms, 1L, 1e-3, dev)),
            .maihda_lme4_trap_drop)
  expect_lt(dev(theta) - dev(maihda_lme4_mirror_theta(theta, cnms, 1L, 1e-3)),
            .maihda_lme4_trap_drop)
  # One entry below: the reflection's side, with that entry as large as the block allows.
  two <- list(stratum = c("(Intercept)", "x"))
  across2 <- maihda_lme4_first_order_theta(c(0, 0.3, 0.4), two, 1L, 0.05,
                                           function(th) 2 * th[1] * th[2] + 5 * th[1]^2)
  expect_equal(across2, c(0.05, -0.5, 0), tolerance = 1e-9)
  # No direction from a deviance that is flat or cannot be evaluated, or for an index that
  # is not a column's diagonal entry.
  expect_null(maihda_lme4_first_order_theta(theta, cnms, 1L, 0.05, function(th) 7))
  expect_null(maihda_lme4_first_order_theta(theta, cnms, 1L, 0.05, function(th) NA_real_))
  expect_null(maihda_lme4_first_order_theta(theta, cnms, 2L, 0.05, dev))
})

# ---- the refit ----------------------------------------------------------------------------

test_that("the refit's start takes each column across in turn", {
  # A two-term block at its bound, and a four-term block with its first two columns at 0.
  cnms <- list(id = c("(Intercept)", "t"), stratum = c("(Intercept)", "x", "z", "w"))
  theta <- c(0, 0.3, 0.4, 0, 0.2, -0.1, 0.3, 0, 0.05, -0.2, 0.6, 0.1, 0.5)
  cols <- maihda_lme4_trap_columns(theta, cnms)
  expect_identical(vapply(cols, `[[`, integer(1), "index"), c(1L, 4L, 8L))
  g <- c(1.5, -2, 0.7)
  dev <- function(th) {
    S <- block_covs(th, cnms)
    2 * S[[1]][2, 1] + sum(g * S[[2]][2:4, 1]) + 3 * S[[2]][1, 1]
  }
  start <- maihda_lme4_cross_start(theta, cnms, cols, 0.05, dev)
  # The two-term block is reflected, though its first-order direction would differ; the
  # four-term block's first column goes along the first-order direction to 0.9 of the
  # length the block allows, which takes its second column off the bound, so that one is
  # left as it then is.
  expect_identical(start[1:3], c(0.05, -0.3, 0.4))
  expect_equal(maihda_lme4_first_order_theta(theta, cnms, 1L, 0.05, dev)[1:3], c(0.05, -0.5, 0),
               tolerance = 1e-9)
  reflected <- maihda_lme4_mirror_theta(theta, cnms, 1L, 0.05)
  want <- maihda_lme4_first_order_theta(reflected, cnms, 4L, 0.05, dev, fraction = 0.9)
  expect_identical(start, want)
  expect_gt(start[8], 1e-4)
  # The later terms' block starts non-singular; at the full length it would not.
  expect_gt(min(start[c(8, 11, 13)]), 1e-3)
  full <- maihda_lme4_first_order_theta(reflected, cnms, 4L, 0.05, dev)
  expect_lt(min(abs(full[c(8, 11, 13)])), 1e-8)
  # Without a deviance function every column is reflected.
  expect_identical(maihda_lme4_cross_start(theta, cnms, cols, 0.05),
                   maihda_lme4_mirror_theta(theta, cnms, c(1L, 4L, 8L), 0.05))
  # A direction that cannot be found falls back to the reflection.
  expect_identical(maihda_lme4_cross_start(theta, cnms, cols, 0.05, function(th) NA_real_),
                   maihda_lme4_mirror_theta(theta, cnms, c(1L, 4L, 8L), 0.05))
})

test_that("the refit takes a column followed by two or more terms along the first-order direction", {
  skip_on_cran()
  skip_slow_test()
  env <- new.env()
  env$data <- three_term_design(3)
  env$stay <- stay
  held <- quote(lme4::glmer(y ~ x + z + A + B + (1 + x + z | stratum), data = data,
                            family = stats::poisson,
                            control = lme4::glmerControl(nAGQ0initStep = FALSE, optimizer = stay)))
  held$start <- list(theta = trap_3, fixef = fixef_3)
  kept <- quiet(maihda_eval_holding_conditions(held, env))
  expect_identical(unname(lme4::getME(kept$value, "theta")), trap_3)
  call <- quote(lme4::glmer(y ~ x + z + A + B + (1 + x + z | stratum), data = data,
                            family = stats::poisson))
  seen <- NULL
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(maihda_eval_holding_conditions = function(call, env) {
    seen <<- call
    real(call, env)
  })
  cnms <- lme4::getME(kept$value, "cnms")
  out <- quiet(maihda_glmer_mirror_refit(kept, call, env, optimizer = "bobyqa", first_order = TRUE))
  want <- maihda_lme4_first_order_theta(trap_3, cnms, 1L, 0.05,
                                        maihda_lme4_theta_deviance(kept$value), fraction = 0.9)
  expect_equal(seen$start$theta, unname(want))
  # The reflection's start would have been another one, and it went back to the bound.
  mirrored <- unname(maihda_lme4_mirror_theta(trap_3, cnms, 1L, 0.05))
  expect_gt(max(abs(seen$start$theta - mirrored)), 0.05)
  expect_lt(abs(as.numeric(stats::logLik(out$value)) - (-2482.34487)), 0.005)
  # Without first_order, as for growth fits, the reflection.
  quiet(maihda_glmer_mirror_refit(kept, call, env, optimizer = "bobyqa"))
  expect_identical(seen$start$theta, mirrored)
})

test_that("a column followed by one term is reflected, without a deviance function", {
  skip_on_cran()
  env <- new.env()
  env$data <- slope_design(7)
  call <- quote(lme4::glmer(y ~ x + A + B + C + (1 + x | stratum), data = data,
                            family = stats::poisson))
  trapped_call <- call
  trapped_call$start <- list(theta = trap_7)
  kept <- quiet(maihda_eval_holding_conditions(trapped_call, env))
  th <- lme4::getME(kept$value, "theta")
  expect_lt(th[[1]], 1e-4)
  seen <- NULL
  built <- 0L
  real <- maihda_eval_holding_conditions
  local_mocked_bindings(
    maihda_eval_holding_conditions = function(call, env) {
      seen <<- call
      real(call, env)
    },
    maihda_lme4_theta_deviance = function(model) {
      built <<- built + 1L
      NULL
    })
  quiet(maihda_glmer_mirror_refit(kept, call, env, optimizer = "bobyqa", first_order = TRUE))
  expect_identical(seen$start$theta,
                   unname(maihda_lme4_mirror_theta(th, lme4::getME(kept$value, "cnms"), 1L, 0.05)))
  expect_identical(built, 0L)
})

test_that("a three-term fit the reflection left short is refitted to its maximum", {
  skip_on_cran()
  skip_slow_test()
  m <- quiet(fit_maihda(three_term_formula, data = three_term_design(3), family = "poisson"))
  expect_lt(abs(as.numeric(stats::logLik(m$model)) - (-2482.34487)), 0.005)
  expect_gt(lme4::getME(m$model, "theta")[[1]], 0.02)
  expect_true(m$diagnostics$converged)
  expect_null(m$diagnostics$boundary_traps)
})

# ---- the record ---------------------------------------------------------------------------

test_that("the record's drop is the larger of the first-order and the reflection's", {
  # A two-term block at its bound; the deviance falls with the covariance and curves fast in
  # it, so the first-order direction at its full length (0.5) overshoots and the
  # reflection's shorter step (0.3) falls further.
  two <- list(stratum = c("(Intercept)", "x"))
  theta <- c(0, 0.3, 0.4)
  dev <- function(th) {
    s <- th[1] * th[2]
    12 * s + 2.4e4 * s^2
  }
  first <- dev(theta) - dev(maihda_lme4_first_order_theta(theta, two, 1L, 1e-3, dev))
  reflected <- dev(theta) - dev(maihda_lme4_mirror_theta(theta, two, 1L, 1e-3))
  expect_lt(first, .maihda_lme4_trap_drop)
  expect_gt(reflected, .maihda_lme4_trap_drop)
  expect_identical(maihda_lme4_column_drop(theta, two, 1L, dev), reflected)
  # Where the first-order direction falls further, its fall; NA where nothing can be
  # evaluated.
  cnms <- list(stratum = c("(Intercept)", "x", "z", "w"))
  theta4 <- c(0, 0.2, -0.1, 0.3, 0.4, 0.05, -0.2, 0.6, 0.1, 0.5)
  g <- c(1.5, -2, 0.7)
  dev4 <- function(th) {
    S <- block_covs(th, cnms)[[1]]
    sum(g * S[2:4, 1]) + 3 * S[1, 1]
  }
  expect_identical(maihda_lme4_column_drop(theta4, cnms, 1L, dev4),
                   dev4(theta4) - dev4(maihda_lme4_first_order_theta(theta4, cnms, 1L, 1e-3, dev4)))
  expect_identical(maihda_lme4_column_drop(theta4, cnms, 1L, function(th) NA_real_), NA_real_)
})

test_that("the record measures a stop along the direction in which the likelihood rises", {
  skip_on_cran()
  # Seed 22's stop, held there by the user's start and an optimiser that keeps it.
  m <- quiet(fit_maihda(three_term_formula, data = three_term_design(22), family = "poisson",
                        start = list(theta = trap_22, fixef = fixef_22),
                        control = lme4::glmerControl(nAGQ0initStep = FALSE, optimizer = stay)))
  expect_identical(unname(lme4::getME(m$model, "theta")), trap_22)
  dg <- m$diagnostics
  expect_false(dg$converged)
  expect_identical(dg$boundary_traps$block, "stratum")
  expect_identical(dg$boundary_traps$term, "(Intercept)")
  expect_identical(dg$boundary_traps$n_terms, 3L)
  expect_gt(dg$boundary_traps$drop, 0.04)
  # By the reflection the deviance rises across the bound.
  at <- maihda_lme4_theta_deviance(m$model)
  cnms <- lme4::getME(m$model, "cnms")
  expect_gt(at(maihda_lme4_mirror_theta(trap_22, cnms, 1L, 1e-3)) - at(trap_22), 0.002)
  sentence <- paste0("The intercept variance of the 'stratum' block is at zero, where the ",
                     "likelihood still rises when that variance and its covariances with x ",
                     "and z move off zero: the fit stopped at a bound of lme4's ",
                     "parametrisation, not at the maximum.")
  expect_true(sentence %in% dg$messages)
})

test_that("the record's sentence for a column followed by two or more terms", {
  tail <- ": the fit stopped at a bound of lme4's parametrisation, not at the maximum."
  off <- function(what, later) {
    paste0(what, ", where the likelihood still rises when that variance and its covariances ",
           "with ", later, " move off zero", tail)
  }
  expect_identical(maihda_lme4_trap_message("stratum", "(Intercept)", 3L, c("(Intercept)", "x", "z")),
                   off("The intercept variance of the 'stratum' block is at zero", "x and z"))
  expect_identical(maihda_lme4_trap_message("stratum", "x", 4L, c("(Intercept)", "x", "z", "w")),
                   off("In the 'stratum' block the variance of x left after the intercept is zero",
                       "z and w"))
  expect_identical(maihda_lme4_trap_message("g", "x", 3L, c("x", "z", "w")),
                   off("The variance of x in the 'g' block is at zero", "z and w"))
  # A growth block of three or more terms.
  expect_identical(maihda_lme4_trap_message("id", "(Intercept)", 3L),
                   off("The intercept variance of the 'id' block is at zero", "the time terms"))
  # Followed by one term: the sign of their covariance reversed.
  expect_identical(maihda_lme4_trap_message("stratum", "x", 3L, c("(Intercept)", "x", "z")),
                   paste0("In the 'stratum' block the variance of x left after the intercept is ",
                          "zero, where the likelihood still rises with the sign of its covariance ",
                          "with z reversed", tail))
  expect_identical(maihda_lme4_trap_message("id", "(Intercept)", 2L),
                   paste0("The intercept variance of the 'id' block is at zero, where the ",
                          "likelihood still rises with the sign of its covariance with the slope ",
                          "reversed", tail))
})
