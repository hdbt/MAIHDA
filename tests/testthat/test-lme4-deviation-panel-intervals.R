# The 95% bars that plot_prediction_deviation_panels() draws for lme4 fits.
#
# A stratum panel plots m_j = sum_i a_ji v_i, the weighted mean of the plotted values
# v_i over the stratum's rows (a_ji = w_i / sum_{i in j} w_i, the panel's prior
# weights); v_i is the linear predictor eta_i on a Gaussian panel and g^-1(eta_i) on a
# binomial or count panel. Var(m_j) = d_j' C d_j, with d_j = sum_i a_ji v_i'(eta_i)
# [z_i; x_i] and C the joint covariance of (b, beta) given the variance parameters:
# C = sigma^2 T H^-1 T', T = diag(Lambda, I),
# H = [Lambda'Z'WZ Lambda + I, Lambda'Z'WX; X'WZ Lambda, X'WX] (Henderson 1975, in
# lme4's parametrisation, Bates et al. 2015), W the prior weights of an lmer fit or the
# converged working weights of a glmer fit. The bar is m_j +/- z SE on a Gaussian
# panel, and the delta-method interval on the logit (binomial) or log (count) scale of
# m_j, mapped back. A case-level panel has one bar per row. The references below build
# C densely from that definition; the package solves with lme4's sparse factor, whose
# fill-reducing permutation is not the identity once a fit has several grouping terms.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
ldp_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

# Twelve unbalanced strata (5 to 70 rows), a covariate that varies within strata,
# precision weights and an exposure.
ldp_data <- function() {
  ldp_with_seed(20261027, {
    d <- expand.grid(gender = c("f", "m"), race = c("a", "b", "c"), edu = c("lo", "hi"),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    d <- d[rep(seq_len(nrow(d)), c(60, 8, 45, 12, 70, 5, 30, 15, 50, 9, 40, 20)), ]
    cell <- interaction(d$gender, d$race, d$edu, drop = TRUE)
    d$x <- stats::rnorm(nrow(d))
    d$w <- stats::runif(nrow(d), 0.5, 2)
    d$expo <- stats::runif(nrow(d), 0.5, 3)
    u <- stats::rnorm(nlevels(cell), 0, 0.7)[cell]
    dims <- c(f = 0, m = 0.4)[d$gender] + c(a = 0, b = -0.5, c = 0.6)[d$race] +
      c(lo = 0, hi = 0.3)[d$edu]
    d$y <- 1 + 0.5 * d$x + dims + u + stats::rnorm(nrow(d))
    d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(-0.4 + 0.6 * d$x + u))
    d$cnt <- stats::rpois(nrow(d), d$expo * exp(-0.2 + 0.3 * d$x + u))
    rownames(d) <- NULL
    d
  })
}

# C from its definition, dense.
ldp_ref_cov <- function(fit) {
  X <- as.matrix(lme4::getME(fit, "X"))
  Z <- as.matrix(lme4::getME(fit, "Z"))
  Lam <- as.matrix(lme4::getME(fit, "Lambda"))
  fam <- stats::family(fit)
  eta <- as.numeric(stats::predict(fit, type = "link"))
  pw <- as.numeric(stats::weights(fit))
  W <- if (lme4::isLMM(fit)) pw else pw * fam$mu.eta(eta)^2 / fam$variance(fam$linkinv(eta))
  q <- ncol(Z)
  p <- ncol(X)
  ZL <- Z %*% Lam
  H <- rbind(cbind(crossprod(ZL, W * ZL) + diag(q), crossprod(ZL, W * X)),
             cbind(crossprod(X, W * ZL), crossprod(X, W * X)))
  Tm <- diag(q + p)
  Tm[seq_len(q), seq_len(q)] <- Lam
  stats::sigma(fit)^2 * Tm %*% solve(H) %*% t(Tm)
}

# Reference bars on the fitted rows `rows`: one per stratum, or one per row when
# `strata` is NULL.
ldp_ref_bars <- function(fit, strata = NULL, w = NULL, type = "gaussian",
                         rows = seq_len(stats::nobs(fit))) {
  C <- ldp_ref_cov(fit)
  ZX <- cbind(as.matrix(lme4::getME(fit, "Z")), as.matrix(lme4::getME(fit, "X")))[rows, , drop = FALSE]
  fam <- stats::family(fit)
  eta <- as.numeric(stats::predict(fit, type = "link"))[rows]
  on_response <- type != "gaussian"
  v <- if (on_response) fam$linkinv(eta) else eta
  sl <- if (on_response) fam$mu.eta(eta) else rep(1, length(eta))
  if (is.null(w)) w <- rep(1, length(eta))
  key <- if (is.null(strata)) as.character(seq_along(eta)) else as.character(strata)
  g <- if (is.null(strata)) key else sort(unique(key))
  A <- vapply(g, function(k) ifelse(key == k, w, 0) / sum(w[key == k]), numeric(length(key)))
  m <- colSums(A * v)
  D <- t(A * sl) %*% ZX
  se <- sqrt(rowSums((D %*% C) * D))
  z <- stats::qnorm(0.975)
  ci <- switch(type,
    gaussian = list(m - z * se, m + z * se),
    binomial = list(stats::plogis(stats::qlogis(m) - z * se / (m * (1 - m))),
                    stats::plogis(stats::qlogis(m) + z * se / (m * (1 - m)))),
    poisson = list(exp(log(m) - z * se / m), exp(log(m) + z * se / m)))
  data.frame(key = g, fitted = m, lower = ci[[1]], upper = ci[[2]], stringsAsFactors = FALSE)
}

ldp_warnings <- function(expr) {
  n <- 0L
  value <- withCallingHandlers(expr, warning = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  list(value = value, n = n)
}

# Largest difference between a panel's bars and the reference, matched by stratum (or,
# for a case-level panel, by the row id the panel carries): absolute, or relative for
# a glmer fit, whose stored sparse factor comes from the last weight update of its
# fitting iterations and so agrees with the converged weights to about 1e-6.
ldp_bar_gap <- function(pd, ref, by = c("stratum", "id"), relative = FALSE) {
  by <- match.arg(by)
  idx <- match(as.character(pd[[by]]), ref$key)
  got <- c(pd$fitted, pd$ci_lower, pd$ci_upper)
  want <- c(ref$fitted[idx], ref$lower[idx], ref$upper[idx])
  if (relative) max(abs(got / want - 1)) else max(abs(got - want))
}

test_that("maihda_delta_interval() maps the interval back from the identity, logit and log scales", {
  z <- stats::qnorm(0.975)
  id <- MAIHDA:::maihda_delta_interval(2, 0.5, "identity")
  expect_equal(c(id$lower, id$upper), c(2 - z * 0.5, 2 + z * 0.5))
  lg <- MAIHDA:::maihda_delta_interval(0.2, 0.05, "logit")
  expect_equal(c(lg$lower, lg$upper),
               stats::plogis(stats::qlogis(0.2) + c(-1, 1) * z * 0.05 / 0.16))
  lo <- MAIHDA:::maihda_delta_interval(3, 0.6, "log")
  expect_equal(c(lo$lower, lo$upper), exp(log(3) + c(-1, 1) * z * 0.2))
  # asymmetric about the mean near the edge of the range, and inside it
  expect_gt(lg$upper - 0.2, 0.2 - lg$lower)
  expect_true(lg$lower > 0 && lo$lower > 0)
  # a mean on the edge of its range has no interval
  expect_true(all(is.na(unlist(MAIHDA:::maihda_delta_interval(c(0, 1), c(0.1, 0.1), "logit")))))
  expect_true(all(is.na(unlist(MAIHDA:::maihda_delta_interval(0, 0.1, "log")))))
  half <- MAIHDA:::maihda_delta_interval(0.3, 0.02, "logit", level = 0.5)
  expect_equal(half$upper, stats::plogis(stats::qlogis(0.3) + stats::qnorm(0.75) * 0.02 / 0.21))
})

test_that("maihda_lme4_combination_var() equals the joint covariance of the mixed-model equations", {
  d <- ldp_data()
  fits <- list(
    single = lme4::lmer(y ~ x + (1 | race), data = d),
    crossed = suppressMessages(lme4::lmer(
      y ~ x + (1 | gender) + (1 | race) + (1 | edu) + (1 | gender:race:edu), data = d)),
    binomial = lme4::glmer(yb ~ x + (1 | race) + (1 | gender:race:edu), data = d,
                           family = stats::binomial),
    # a correlated random slope: Lambda is not symmetric, so Lambda and Lambda' differ
    slope = suppressMessages(lme4::lmer(y ~ x + (1 + x | race), data = d))
  )
  perm <- function(f) as.integer(lme4::getME(f, "L")@perm)
  expect_identical(perm(fits$single), seq_along(perm(fits$single)) - 1L)
  expect_false(identical(perm(fits$crossed), seq_along(perm(fits$crossed)) - 1L))
  expect_false(identical(perm(fits$binomial), seq_along(perm(fits$binomial)) - 1L))
  expect_false(isSymmetric(as.matrix(lme4::getME(fits$slope, "Lambda"))))
  tol <- c(single = 1e-9, crossed = 1e-9, binomial = 1e-5, slope = 1e-9)
  for (nm in names(fits)) {
    f <- fits[[nm]]
    C <- ldp_ref_cov(f)
    q <- ncol(lme4::getME(f, "Z"))
    D <- ldp_with_seed(7, matrix(stats::rnorm(nrow(C) * 9), nrow(C), 9))
    v_ref <- colSums(D * (C %*% D))
    v <- MAIHDA:::maihda_lme4_combination_var(f, D[seq_len(q), , drop = FALSE],
                                              D[-seq_len(q), , drop = FALSE], block = 4L)
    expect_lt(max(abs(v / v_ref - 1)), tol[[nm]])
  }
  expect_identical(MAIHDA:::maihda_lme4_combination_var(fits$single, matrix(0, 3, 0),
                                                        matrix(0, 2, 0)), numeric(0))
})

test_that("lme4 stratum bars are the joint-covariance interval of the plotted stratum mean", {
  d <- ldp_data()
  mg <- suppressMessages(fit_maihda(y ~ x + (1 | gender:race:edu), data = d))
  pg <- ldp_warnings(plot_prediction_deviation_panels(mg))
  expect_identical(pg$n, 0L)
  pd <- pg$value[[2]]$data
  expect_lt(ldp_bar_gap(pd, ldp_ref_bars(mg$model, mg$data$stratum)), 1e-8)
  # what the error-bar layer draws is these bars, row for row
  lay <- ggplot2::layer_data(pg$value[[2]], 2)
  expect_equal(lay$ymin, pd$ci_lower[match(lay$x, pd$rank)])
  expect_equal(lay$ymax, pd$ci_upper[match(lay$x, pd$rank)])

  mb <- suppressMessages(fit_maihda(yb ~ x + (1 | gender:race:edu), data = d,
                                    family = "binomial"))
  pb <- ldp_warnings(plot_prediction_deviation_panels(mb))
  expect_identical(pb$n, 0L)
  pdb <- pb$value[[2]]$data
  expect_lt(ldp_bar_gap(pdb, ldp_ref_bars(mb$model, mb$data$stratum, type = "binomial"),
                        relative = TRUE), 1e-5)
  expect_true(all(pdb$ci_lower > 0 & pdb$ci_upper < 1))
  # logit-scale bars: longer towards 0.5
  low <- pdb$fitted < 0.5
  expect_true(any(low))
  expect_true(all((pdb$ci_upper - pdb$fitted)[low] > (pdb$fitted - pdb$ci_lower)[low]))

  mp <- suppressMessages(fit_maihda(cnt ~ x + offset(log(expo)) + (1 | gender:race:edu),
                                    data = d, family = "poisson"))
  pp <- ldp_warnings(plot_prediction_deviation_panels(mp))
  expect_identical(pp$n, 0L)
  pdp <- pp$value[[2]]$data
  expect_lt(ldp_bar_gap(pdp, ldp_ref_bars(mp$model, mp$data$stratum, type = "poisson"),
                        relative = TRUE), 1e-5)
  expect_true(all(pdp$ci_lower > 0))
})

test_that("prior weights weight the stratum mean and its interval alike", {
  d <- ldp_data()
  mw <- suppressMessages(fit_maihda(y ~ x + (1 | gender:race:edu), data = d, weights = w))
  pd <- plot_prediction_deviation_panels(mw)[[2]]$data
  ref <- ldp_ref_bars(mw$model, mw$data$stratum, w = d$w)
  expect_lt(ldp_bar_gap(pd, ref), 1e-8)
  unweighted <- ldp_ref_bars(mw$model, mw$data$stratum)
  expect_gt(max(abs(unweighted$fitted - ref$fitted)), 1e-3)
})

test_that("a crossed-dimensions analysis gets bars from the permuted factor", {
  d <- ldp_data()
  an <- suppressMessages(suppressWarnings(
    maihda(y ~ x + (1 | gender:race:edu), data = d, decomposition = "crossed-dimensions")))
  fit <- an$model$model
  expect_false(identical(as.integer(lme4::getME(fit, "L")@perm),
                         seq_len(ncol(lme4::getME(fit, "Z"))) - 1L))
  p <- ldp_warnings(plot(an, type = "prediction_deviation"))
  expect_identical(p$n, 0L)
  pd <- p$value[[2]]$data
  expect_lt(ldp_bar_gap(pd, ldp_ref_bars(fit, an$model$data$stratum)), 1e-8)
})

test_that("bars for supplied data follow the rows supplied", {
  d <- ldp_data()
  mg <- suppressMessages(fit_maihda(y ~ x + (1 | gender:race:edu), data = d))
  mb <- suppressMessages(fit_maihda(yb ~ x + (1 | gender:race:edu), data = d,
                                    family = "binomial"))
  for (m in list(mg, mb)) {
    type <- if (identical(m, mb)) "binomial" else "gaussian"
    glmm <- type == "binomial"
    tol <- if (glmm) 1e-5 else 1e-8
    all_rows <- plot_prediction_deviation_panels(m, data = m$data)[[2]]$data
    expect_lt(ldp_bar_gap(all_rows, ldp_ref_bars(m$model, m$data$stratum, type = type),
                          relative = glmm), tol)
    rows <- seq(1, nrow(m$data), by = 3)
    sub <- plot_prediction_deviation_panels(m, data = m$data[rows, ])[[2]]$data
    expect_lt(ldp_bar_gap(sub, ldp_ref_bars(m$model, m$data$stratum[rows], type = type,
                                            rows = rows), relative = glmm), tol)
    # rows from half of the strata: the random effects of the supplied rows must be
    # placed among all the fitted ones, not in the order they come
    half <- which(as.integer(factor(m$data$stratum)) %% 2 == 0)
    part <- plot_prediction_deviation_panels(m, data = m$data[half, ])[[2]]$data
    expect_lt(ldp_bar_gap(part, ldp_ref_bars(m$model, m$data$stratum[half], type = type,
                                             rows = half), relative = glmm), tol)
  }
})

test_that("supplied rows are rebuilt with the fit's levels, coefficients and grouping codes", {
  d <- ldp_data()
  # a random slope: each level carries two coefficients, placed level by level
  fs <- suppressMessages(lme4::lmer(y ~ x + (1 + x | race), data = d))
  keep <- which(d$race != "b")
  ps <- plot_prediction_deviation_panels(fs, data = d[keep, ])[[2]]$data
  expect_lt(ldp_bar_gap(ps, ldp_ref_bars(fs, rows = keep), by = "id"), 1e-8)
  # a fixed factor with a level absent from the supplied rows keeps the fit's coding
  ff <- lme4::lmer(y ~ x + edu + (1 | race), data = d)
  hi <- which(d$edu == "hi")
  pf <- plot_prediction_deviation_panels(ff, data = d[hi, ])[[2]]$data
  expect_true(all(is.finite(pf$ci_lower)))
  expect_lt(ldp_bar_gap(pf, ldp_ref_bars(ff, rows = hi), by = "id"), 1e-8)
  # a character copy of an integer-coded stratum: "10" sorts before "2" as text
  mg <- suppressMessages(fit_maihda(y ~ x + (1 | gender:race:edu), data = d))
  expect_true(is.numeric(mg$data$stratum) && length(unique(mg$data$stratum)) >= 10)
  sup <- mg$data
  sup$stratum <- as.character(sup$stratum)
  pc <- plot_prediction_deviation_panels(mg, data = sup)[[2]]$data
  expect_lt(ldp_bar_gap(pc, ldp_ref_bars(mg$model, mg$data$stratum)), 1e-8)
})

test_that("case-level bars of a bare lme4 fit come from the joint covariance, fitted or supplied rows", {
  d <- ldp_data()
  f <- lme4::lmer(y ~ x + (1 | race) + (1 | edu), data = d)
  expect_false(identical(as.integer(lme4::getME(f, "L")@perm),
                         seq_len(ncol(lme4::getME(f, "Z"))) - 1L))
  ref <- ldp_ref_bars(f)
  p0 <- ldp_warnings(plot_prediction_deviation_panels(f))
  p1 <- ldp_warnings(plot_prediction_deviation_panels(f, data = d))
  expect_identical(c(p0$n, p1$n), c(0L, 0L))
  expect_lt(ldp_bar_gap(p0$value[[2]]$data, ref, by = "id"), 1e-8)
  expect_lt(ldp_bar_gap(p1$value[[2]]$data, ref, by = "id"), 1e-8)
  # supplied rows without one race: its random effect drops out of the new design
  keep <- which(d$race != "b")
  p2 <- plot_prediction_deviation_panels(f, data = d[keep, ])
  expect_lt(ldp_bar_gap(p2[[2]]$data, ldp_ref_bars(f, rows = keep), by = "id"), 1e-8)

  fb <- lme4::glmer(yb ~ x + (1 | race), data = d, family = stats::binomial)
  pb <- ldp_warnings(plot_prediction_deviation_panels(fb, data = d))
  expect_identical(pb$n, 0L)
  pdb <- pb$value[[2]]$data
  expect_lt(ldp_bar_gap(pdb, ldp_ref_bars(fb, type = "binomial"), by = "id", relative = TRUE), 1e-5)
  expect_true(all(pdb$ci_lower > 0 & pdb$ci_upper < 1))
})

test_that("supplied rows whose random effects cannot be placed draw no bars", {
  d <- ldp_data()
  # one grouping factor carrying two random-effect terms
  f <- suppressMessages(lme4::lmer(y ~ x + (1 | race) + (0 + x | race), data = d))
  fitted_rows <- plot_prediction_deviation_panels(f)[[2]]$data
  expect_lt(ldp_bar_gap(fitted_rows, ldp_ref_bars(f), by = "id"), 1e-8)
  supplied <- ldp_warnings(plot_prediction_deviation_panels(f, data = d))
  expect_identical(supplied$n, 0L)
  expect_true(all(is.na(supplied$value[[2]]$data$ci_lower)))
  expect_true(all(is.na(supplied$value[[2]]$data$ci_upper)))
})
