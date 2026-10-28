# Audit 2026-10-06, packet glmm-pcv-ml-bias.
#
# The PCV of a glmer MAIHDA compares two maximum-likelihood variances. With few or
# small strata these run low, the adjusted model's most, so the PCV runs high.
# calculate_pcv() and maihda() take estimation = "REML" for binomial and Poisson lme4
# fits: REML-type variances, the Laplace approximation with the fixed effects
# integrated out as well as the random effects. print() says when the compared
# variances are maximum-likelihood estimates. No default changed.

# Twelve strata (2 x 2 x 3) with a fixed table of successes: additive effects on the
# logit plus a fixed set of stratum departures. No random numbers.
audit_glmm_table <- function(departures = TRUE) {
  grid <- expand.grid(a = factor(1:2), b = factor(1:2), c = factor(1:3))
  lin <- -0.6 + 0.5 * (grid$a == "2") + 0.4 * (grid$b == "2") + c(0, 0.3, 0.6)[grid$c]
  if (departures) {
    lin <- lin + c(0.35, -0.30, 0.10, -0.25, 0.40, -0.15, -0.35, 0.30, 0.05, 0.20,
                   -0.40, 0.05)
    grid$n <- c(60, 45, 80, 30, 55, 70, 40, 65, 50, 75, 35, 60)
  } else {
    grid$n <- 200
  }
  grid$s <- round(grid$n * stats::plogis(lin))
  grid
}

# The table as one Bernoulli row per observation.
audit_glmm_rows <- function(grid) {
  d <- grid[rep(seq_len(nrow(grid)), grid$n), c("a", "b", "c")]
  d$y <- unlist(lapply(seq_len(nrow(grid)), function(j) {
    rep(c(1L, 0L), c(grid$s[j], grid$n[j] - grid$s[j]))
  }))
  rownames(d) <- NULL
  d
}

audit_quiet <- function(expr) suppressMessages(suppressWarnings(expr))

# Null and adjusted logistic fits on the Bernoulli rows.
audit_glmm_fits <- function(grid = audit_glmm_table()) {
  data <- make_strata(audit_glmm_rows(grid), c("a", "b", "c"))$data
  list(
    data = data,
    null = audit_quiet(fit_maihda(y ~ 1 + (1 | stratum), data = data,
                                  family = "binomial")),
    adjusted = audit_quiet(fit_maihda(y ~ a + b + c + (1 | stratum), data = data,
                                      family = "binomial"))
  )
}

# The REML-type variance from its definition, with no lme4 object involved: minus
# twice the log of the Laplace approximation of the likelihood with BOTH the fixed
# effects (flat prior) and the standard-normal random effects integrated out, taken
# at their joint mode, minimised over the random-effect standard deviation theta.
# One row per stratum; `trials` for a binomial count, `offset` for a Poisson rate.
audit_laplace_reml_variance <- function(X, y, trials = NULL, offset = 0,
                                        family = c("binomial", "poisson")) {
  family <- match.arg(family)
  p <- ncol(X)
  J <- nrow(X)
  criterion <- function(theta) {
    pieces <- function(par) {
      u <- par[p + seq_len(J)]
      lin <- as.numeric(X %*% par[seq_len(p)] + theta * u + offset)
      if (family == "binomial") {
        mu <- trials * stats::plogis(lin)
        w <- mu * (1 - stats::plogis(lin))
        loglik <- sum(y * lin - trials * log1p(exp(lin)))
      } else {
        mu <- exp(lin)
        w <- mu
        loglik <- sum(y * lin - mu)
      }
      list(u = u, loglik = loglik,
           grad = c(crossprod(X, y - mu), theta * (y - mu) - u),
           hess = rbind(cbind(crossprod(X, w * X), theta * t(w * X)),
                        cbind(theta * (w * X), diag(theta^2 * w + 1, J))))
    }
    # Newton steps to the joint mode, halved while they fail to raise the penalised
    # log-likelihood, from the intercept that matches the overall rate.
    penalised <- function(k) k$loglik - 0.5 * sum(k$u^2)
    start <- if (family == "binomial") {
      stats::qlogis(sum(y) / sum(trials))
    } else {
      log(sum(y) / sum(exp(rep_len(offset, J))))
    }
    par <- c(start, rep(0, p - 1 + J))
    for (iteration in 1:500) {
      k <- pieces(par)
      step <- solve(k$hess, k$grad)
      repeat {
        value <- penalised(pieces(par + step))
        if ((is.finite(value) && value >= penalised(k) - 1e-12) ||
            max(abs(step)) < 1e-14) {
          break
        }
        step <- step / 2
      }
      par <- par + step
      if (max(abs(step)) < 1e-12) break
    }
    k <- pieces(par)
    -2 * (k$loglik - 0.5 * sum(k$u^2)) +
      as.numeric(determinant(k$hess, logarithm = TRUE)$modulus)
  }
  opt <- stats::optimize(criterion, c(0, 3), tol = 1e-10)
  if (criterion(0) <= opt$objective) 0 else opt$minimum^2
}

test_that("the REML-type variance is the Laplace estimate with the fixed effects integrated out", {
  grid <- audit_glmm_table()
  expect_equal(grid$s, c(26, 18, 38, 15, 29, 36, 18, 46, 26, 50, 18, 43))
  fits <- audit_glmm_fits(grid)

  got_null <- MAIHDA:::maihda_restricted_variance_lme4(fits$null$model)
  got_adj <- MAIHDA:::maihda_restricted_variance_lme4(fits$adjusted$model)

  # Against the definition, computed from the 12-row table without lme4.
  want_null <- audit_laplace_reml_variance(matrix(1, nrow(grid), 1), grid$s, grid$n)
  want_adj <- audit_laplace_reml_variance(stats::model.matrix(~ a + b + c, grid),
                                          grid$s, grid$n)
  expect_gt(want_adj, 0.01)
  expect_lt(abs(got_null - want_null), 1e-5)
  expect_lt(abs(got_adj - want_adj), 1e-5)

  # Against glmmTMB(REML = TRUE) 1.1.14 on the same table (values recorded here so
  # the package is not needed).
  expect_lt(abs(got_null - 0.12775011), 1e-5)
  expect_lt(abs(got_adj - 0.01482264), 1e-5)
})

test_that("the one-parameter search finds a minimum next to zero, at zero, or far out", {
  search <- MAIHDA:::maihda_minimise_nonneg

  # An interior minimum well inside the first step of a scan in equal steps, on a
  # criterion that is flat to rounding next to zero.
  near <- search(function(t) 50 * (t - 0.058)^2 + 1e-11 * sin(1e9 * t))
  expect_identical(near$status, "ok")
  expect_lt(abs(near$theta - 0.058), 1e-5)
  expect_identical(near$n_failed, 0L)
  # A small interior minimum is found too.
  expect_lt(abs(search(function(t) 1000 * (t - 0.002)^2)$theta - 0.002), 1e-6)

  # A minimum at zero comes back as exactly zero.
  expect_identical(search(function(t) 3 + 40 * t^2)$theta, 0)
  expect_identical(search(function(t) 3 + 0.2 * t)$theta, 0)
  # A dip below zero's value that only rounding could make is not an interior minimum.
  expect_identical(search(function(t) 3 - 1e-10 * (t > 1e-3 && t < 2e-3))$theta, 0)

  # Beyond the first upper limit the scan widens.
  far <- search(function(t) (t - 5)^2, upper = 2)
  expect_identical(far$status, "ok")
  expect_lt(abs(far$theta - 5), 1e-5)
  expect_gte(far$upper, 8)

  # No minimum, and nothing to evaluate.
  expect_identical(search(function(t) -t)$status, "unbounded")
  expect_identical(search(function(t) stop("no"))$status, "unevaluable")
  # Failed evaluations are counted and rank below every finite value.
  some <- search(function(t) if (t > 1) stop("no") else 50 * (t - 0.3)^2)
  expect_lt(abs(some$theta - 0.3), 1e-6)
  expect_gt(some$n_failed, 0L)
})

test_that("an interior minimum is found whatever maximum likelihood made of the variance", {
  # Twelve strata of 300 with a small interaction. glmer puts the adjusted variance
  # at or next to zero here (a relative standard deviation of 5e-8 on the machine
  # this was written on), while the REML-type criterion has its minimum well inside:
  # 0.0034058, which glmmTMB(REML = TRUE) 1.1.14 also returns.
  set.seed(26)
  grid <- expand.grid(a = factor(1:3), b = factor(1:4))
  lin <- -0.3 + seq(0, 0.6, length.out = 3)[grid$a] +
    seq(0, 0.5, length.out = 4)[grid$b] + stats::rnorm(12, 0, 0.07)
  rows <- grid[rep(1:12, each = 300), ]
  x <- stats::rnorm(3600)
  rows$y <- stats::rbinom(3600, 1, stats::plogis(rep(lin, each = 300) + 0.4 * x))
  null <- audit_quiet(fit_maihda(y ~ 1 + (1 | a:b), data = rows, family = "binomial"))
  adj <- audit_quiet(fit_maihda(y ~ a + b + (1 | a:b), data = rows,
                                family = "binomial"))
  res <- calculate_pcv(null, adj, estimation = "REML")
  expect_lt(res$var_model2_fitted, 1e-6)
  expect_lt(abs(res$var_model2 - 0.0034058), 1e-4)
  expect_false(res$adjusted_at_boundary)
  expect_lt(res$pcv, 0.99)
})

test_that("estimation = \"REML\" compares the REML-type variances and keeps the fitted ones", {
  fits <- audit_glmm_fits()
  before <- list(lme4::getME(fits$null$model, "theta"), lme4::fixef(fits$adjusted$model),
                 lme4::VarCorr(fits$adjusted$model))

  def <- calculate_pcv(fits$null, fits$adjusted)
  reml <- calculate_pcv(fits$null, fits$adjusted, estimation = "REML")

  # The default is lme4's own two variances, as before.
  expect_identical(def$var_model1,
                   as.numeric(lme4::VarCorr(fits$null$model)$stratum[1, 1]))
  expect_identical(def$var_model2,
                   as.numeric(lme4::VarCorr(fits$adjusted$model)$stratum[1, 1]))
  expect_identical(def$estimation, "fitted")
  expect_identical(def$estimation_used, "fitted")
  expect_true(def$ml_variances)
  expect_true(def$reml_available)
  expect_null(def$pcv_fitted)
  # On this table maximum likelihood puts the adjusted variance on the boundary.
  expect_true(def$adjusted_at_boundary)
  expect_gt(def$pcv, 1 - 1e-6)

  expect_s3_class(reml, "pcv_result")
  expect_identical(reml$estimation, "REML")
  expect_identical(reml$estimation_used, "REML")
  expect_identical(reml$var_model1,
                   MAIHDA:::maihda_restricted_variance_lme4(fits$null$model))
  expect_identical(reml$var_model2,
                   MAIHDA:::maihda_restricted_variance_lme4(fits$adjusted$model))
  expect_identical(reml$pcv, (reml$var_model1 - reml$var_model2) / reml$var_model1)
  expect_s3_class(tryCatch(reml$pvc, warning = function(w) w), "deprecatedWarning")
  expect_identical(unclass(reml)$pvc, reml$pcv)
  expect_lt(abs(reml$pcv - 0.883971), 1e-4)
  expect_false(reml$adjusted_at_boundary)
  expect_false(reml$ml_variances)
  # Both variances are larger than their maximum-likelihood counterparts, the
  # adjusted one by proportionally more, so the PCV is lower.
  expect_gt(reml$var_model1, def$var_model1)
  expect_gt(reml$var_model2, def$var_model2)
  expect_lt(def$var_model2 / reml$var_model2, def$var_model1 / reml$var_model1 - 0.5)
  expect_lt(reml$pcv, def$pcv - 0.05)
  # The fitted values stay on the result.
  expect_identical(reml$var_model1_fitted, def$var_model1)
  expect_identical(reml$var_model2_fitted, def$var_model2)
  expect_identical(reml$pcv_fitted, def$pcv)

  # The two fits are read, not changed.
  expect_identical(lme4::getME(fits$null$model, "theta"), before[[1]])
  expect_identical(lme4::fixef(fits$adjusted$model), before[[2]])
  expect_identical(lme4::VarCorr(fits$adjusted$model), before[[3]])
})

test_that("an aggregated cbind() fit gives the Bernoulli fit's REML-type variances", {
  grid <- audit_glmm_table()
  fits <- audit_glmm_fits(grid)
  agg <- make_strata(grid, c("a", "b", "c"))$data
  agg$f <- agg$n - agg$s
  null <- audit_quiet(fit_maihda(cbind(s, f) ~ 1 + (1 | stratum), data = agg,
                                 family = "binomial"))
  adj <- audit_quiet(fit_maihda(cbind(s, f) ~ a + b + c + (1 | stratum), data = agg,
                                family = "binomial"))
  from_counts <- calculate_pcv(null, adj, estimation = "REML")
  from_rows <- calculate_pcv(fits$null, fits$adjusted, estimation = "REML")
  expect_identical(from_counts$estimation_used, "REML")
  expect_lt(abs(from_counts$var_model1 - from_rows$var_model1), 1e-6)
  expect_lt(abs(from_counts$var_model2 - from_rows$var_model2), 1e-6)

  # The same table as proportions with the trial counts as prior weights: the weights
  # are read from the fit, so the answer is again the table's.
  agg$p <- agg$s / agg$n
  null_w <- audit_quiet(fit_maihda(p ~ 1 + (1 | stratum), data = agg,
                                   family = "binomial", weights = n))
  adj_w <- audit_quiet(fit_maihda(p ~ a + b + c + (1 | stratum), data = agg,
                                  family = "binomial", weights = n))
  from_weights <- calculate_pcv(null_w, adj_w, estimation = "REML")
  expect_lt(abs(from_weights$var_model1 - from_rows$var_model1), 1e-6)
  expect_lt(abs(from_weights$var_model2 - from_rows$var_model2), 1e-6)
  # Without the weights the twelve proportions would be twelve single trials.
  unweighted <- audit_laplace_reml_variance(matrix(1, nrow(grid), 1), agg$p,
                                            rep(1, nrow(grid)))
  expect_gt(abs(unweighted - from_rows$var_model1), 0.01)
})

test_that("a Poisson fit with an exposure offset is handled on its own rows", {
  grid <- audit_glmm_table()
  grid$expo <- grid$n
  grid$cnt <- round(grid$expo * exp(
    -1.2 + 0.4 * (grid$a == "2") + c(0, 0.2, 0.5)[grid$c] +
      c(0.45, -0.30, 0.15, -0.45, 0.40, -0.15, -0.40, 0.30, 0.10, 0.25, -0.45, 0.15)))
  grid$lexpo <- log(grid$expo)
  data <- make_strata(grid, c("a", "b", "c"))$data
  null <- audit_quiet(fit_maihda(cnt ~ 1 + offset(lexpo) + (1 | stratum), data = data,
                                 family = "poisson"))
  adj <- audit_quiet(fit_maihda(cnt ~ a + b + c + offset(lexpo) + (1 | stratum),
                                data = data, family = "poisson"))
  res <- calculate_pcv(null, adj, estimation = "REML")
  expect_identical(res$estimation_used, "REML")

  want_null <- audit_laplace_reml_variance(matrix(1, nrow(grid), 1), grid$cnt,
                                           offset = grid$lexpo, family = "poisson")
  want_adj <- audit_laplace_reml_variance(stats::model.matrix(~ a + b + c, grid),
                                          grid$cnt, offset = grid$lexpo,
                                          family = "poisson")
  expect_gt(want_adj, 0.01)
  # lme4's inner iterations stop at their own tolerance, which leaves the minimiser
  # within about 1e-5 of the oracle's on these counts.
  expect_lt(abs(res$var_model1 - want_null), 1e-4)
  expect_lt(abs(res$var_model2 - want_adj), 1e-4)
  # Dropping the offset from the oracle moves it, so the offset is being used.
  no_offset <- audit_laplace_reml_variance(matrix(1, nrow(grid), 1), grid$cnt,
                                           family = "poisson")
  expect_gt(abs(no_offset - want_null), 0.01)
  expect_gt(res$var_model2, res$var_model2_fitted)
})

test_that("an exactly additive table gives a REML-type adjusted variance of zero", {
  fits <- audit_glmm_fits(audit_glmm_table(departures = FALSE))
  res <- audit_quiet(calculate_pcv(fits$null, fits$adjusted, estimation = "REML"))
  expect_identical(res$var_model2, 0)
  expect_true(res$adjusted_at_boundary)
  expect_identical(res$pcv, 1)
  expect_gt(res$var_model1, res$var_model1_fitted)
  expect_output(print(res), "singularity boundary")
})

test_that("REML lmer fits are compared as fitted under estimation = \"REML\"", {
  data <- audit_glmm_fits()$data
  data$yg <- 0.5 * (data$a == "2") + 0.3 * as.integer(data$c) +
    0.4 * sin(as.integer(data$stratum) * 2.3) + sin(seq_len(nrow(data)) * 1.7)
  null <- audit_quiet(fit_maihda(yg ~ 1 + (1 | stratum), data = data))
  adj <- audit_quiet(fit_maihda(yg ~ a + b + c + (1 | stratum), data = data))
  def <- calculate_pcv(null, adj)
  reml <- calculate_pcv(null, adj, estimation = "REML")
  expect_identical(reml$pcv, def$pcv)
  expect_identical(reml$var_model1, def$var_model1)
  expect_identical(reml$var_model2, def$var_model2)
  expect_identical(reml$estimation, "REML")
  expect_identical(reml$estimation_used, "fitted")
  expect_false(def$ml_variances)
  expect_false(def$reml_available)
  expect_null(reml$pcv_fitted)
  out <- paste(utils::capture.output(print(def)), collapse = "\n")
  expect_match(out, "Variance basis: as fitted (REML for Gaussian lmer, matching summary())",
               fixed = TRUE)
  expect_false(grepl("maximum-likelihood variances", out, fixed = TRUE))

  # The same pair fitted by maximum likelihood: the note appears, without pointing to
  # an option that does not apply, and estimation = "REML" is refused.
  null_ml <- audit_quiet(fit_maihda(yg ~ 1 + (1 | stratum), data = data, REML = FALSE))
  adj_ml <- audit_quiet(fit_maihda(yg ~ a + b + c + (1 | stratum), data = data,
                                   REML = FALSE))
  def_ml <- calculate_pcv(null_ml, adj_ml)
  expect_true(def_ml$ml_variances)
  expect_false(def_ml$reml_available)
  out_ml <- paste(utils::capture.output(print(def_ml)), collapse = "\n")
  expect_match(out_ml, "Variance basis: as fitted (maximum likelihood)", fixed = TRUE)
  expect_match(out_ml, "Note: these are maximum-likelihood variances", fixed = TRUE)
  expect_false(grepl("estimation = \"REML\"", out_ml, fixed = TRUE))
  expect_error(calculate_pcv(null_ml, adj_ml, estimation = "REML"),
               "fitted by maximum likelihood", fixed = TRUE)
})

test_that("estimation = \"REML\" is refused where no REML-type variance exists", {
  fits <- audit_glmm_fits()

  # A bootstrap interval would be on the maximum-likelihood refits.
  expect_error(calculate_pcv(fits$null, fits$adjusted, estimation = "REML",
                             bootstrap = TRUE),
               "bootstrap = TRUE is not available with estimation = \"REML\"",
               fixed = TRUE)

  # A further random effect.
  data <- fits$data
  data$site <- factor(rep_len(c("s1", "s2", "s3", "s4", "s5"), nrow(data)))
  null_ctx <- audit_quiet(fit_maihda(y ~ 1 + (1 | stratum), data = data,
                                     family = "binomial", context = "site"))
  adj_ctx <- audit_quiet(fit_maihda(y ~ a + b + c + (1 | stratum), data = data,
                                    family = "binomial", context = "site"))
  expect_error(calculate_pcv(null_ctx, adj_ctx, estimation = "REML"),
               "only random effect is the stratum intercept", fixed = TRUE)
  expect_false(audit_quiet(calculate_pcv(null_ctx, adj_ctx))$reml_available)

  # The engines that estimate by maximum likelihood only, and a brms posterior.
  expect_error(MAIHDA:::maihda_pcv_reml_check(list(engine = "wemix")),
               "not available for engine = \"wemix\"", fixed = TRUE)
  expect_error(MAIHDA:::maihda_pcv_reml_check(list(engine = "ordinal")),
               "not available for engine = \"ordinal\"", fixed = TRUE)
  expect_true(MAIHDA:::maihda_pcv_reml_check(list(engine = "brms")))
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("REML", FALSE, engine = "brms"),
                   "posterior")
  expect_identical(MAIHDA:::maihda_pcv_estimation_used("REML", FALSE, engine = "lme4"),
                   "fitted")

  # A glmer family other than binomial and Poisson.
  gd <- data.frame(g = factor(rep(1:12, each = 6)))
  gd$y <- exp(0.3 * sin(as.integer(gd$g) * 1.9) + 0.4 * cos(seq_len(nrow(gd)) * 2.1))
  gamma_fit <- audit_quiet(lme4::glmer(y ~ 1 + (1 | g), data = gd,
                                       family = stats::Gamma(link = "log")))
  expect_error(MAIHDA:::maihda_pcv_reml_check(list(engine = "lme4", model = gamma_fit)),
               "binomial and Poisson glmer fits", fixed = TRUE)
})

test_that("print() names a maximum-likelihood basis and the REML-type one", {
  fits <- audit_glmm_fits()
  def <- calculate_pcv(fits$null, fits$adjusted)
  reml <- calculate_pcv(fits$null, fits$adjusted, estimation = "REML")

  lines <- utils::capture.output(print(def))
  out <- paste(lines, collapse = "\n")
  expect_match(out, "Variance basis: as fitted (maximum likelihood)", fixed = TRUE)
  expect_match(out, "Note: these are maximum-likelihood variances", fixed = TRUE)
  # The argument is never broken across two lines.
  expect_true(any(grepl("estimation = \"REML\" compares", lines, fixed = TRUE)))
  expect_lte(max(nchar(cli::ansi_strip(lines))), 84)

  lines_reml <- utils::capture.output(print(reml))
  out_reml <- paste(lines_reml, collapse = "\n")
  expect_match(out_reml, "Variance basis: REML-type", fixed = TRUE)
  expect_match(out_reml, sprintf("As fitted (maximum likelihood): %.6f -> %.6f, PCV %.4f",
                                 def$var_model1, def$var_model2, def$pcv), fixed = TRUE)
  expect_false(grepl("Note: these are maximum-likelihood variances", out_reml,
                     fixed = TRUE))

  # Where the option does not apply the note does not name it.
  other <- def
  other$reml_available <- FALSE
  out_other <- paste(utils::capture.output(print(other)), collapse = "\n")
  expect_match(out_other, "Note: these are maximum-likelihood variances", fixed = TRUE)
  expect_false(grepl("estimation = \"REML\"", out_other, fixed = TRUE))
  # Nor beside a bootstrap interval, which the option does not give for a glmer pair.
  boot <- def
  boot$bootstrap <- TRUE
  boot$ci_lower <- 0.5
  boot$ci_upper <- 1
  boot$conf_level <- 0.95
  out_boot <- paste(utils::capture.output(print(boot)), collapse = "\n")
  expect_match(out_boot, "Note: these are maximum-likelihood variances", fixed = TRUE)
  expect_false(grepl("estimation = \"REML\"", out_boot, fixed = TRUE))
  # A result saved before these fields existed prints as it did.
  old <- def
  old$ml_variances <- NULL
  old$reml_available <- NULL
  out_old <- paste(utils::capture.output(print(old)), collapse = "\n")
  expect_match(out_old, "Variance basis: as fitted (REML for Gaussian lmer",
               fixed = TRUE)
  expect_false(grepl("Note: these are maximum-likelihood variances", out_old,
                     fixed = TRUE))
})

test_that("maihda() takes estimation = \"REML\" in the two-model decomposition only", {
  grid <- audit_glmm_table()
  rows <- audit_glmm_rows(grid)
  fits <- audit_glmm_fits(grid)
  want <- calculate_pcv(fits$null, fits$adjusted, estimation = "REML")

  def <- audit_quiet(maihda(y ~ a + b + c + (1 | a:b:c), data = rows,
                            family = "binomial"))
  expect_identical(def$pcv$estimation_used, "fitted")
  out_def <- paste(utils::capture.output(print(def)), collapse = "\n")
  expect_match(out_def, "Note: these are maximum-likelihood variances", fixed = TRUE)
  expect_match(out_def, "estimation = \"REML\" compares", fixed = TRUE)
  # maihda() refuses the option with a group comparison, so that note does not name it.
  grouped <- def
  grouped$groups <- data.frame(group = c("g1", "g2"), n = c(10L, 12L))
  grouped$group_var <- "c"
  out_grouped <- paste(utils::capture.output(print(grouped)), collapse = "\n")
  expect_match(out_grouped, "Note: these are maximum-likelihood variances", fixed = TRUE)
  expect_false(grepl("estimation = \"REML\"", out_grouped, fixed = TRUE))

  reml <- audit_quiet(maihda(y ~ a + b + c + (1 | a:b:c), data = rows,
                             family = "binomial", estimation = "REML"))
  expect_identical(reml$pcv$estimation_used, "REML")
  expect_lt(abs(reml$pcv$var_model1 - want$var_model1), 1e-6)
  expect_lt(abs(reml$pcv$var_model2 - want$var_model2), 1e-6)
  expect_lt(abs(reml$pcv$pcv - want$pcv), 1e-6)
  out_reml <- paste(utils::capture.output(print(reml)), collapse = "\n")
  expect_match(out_reml, "REML-type variances (estimation = \"REML\")", fixed = TRUE)
  expect_match(out_reml, "As fitted (maximum likelihood):", fixed = TRUE)
  expect_false(grepl("Note: these are maximum-likelihood variances", out_reml,
                     fixed = TRUE))
  # The single-model summaries stay the fits' own.
  expect_identical(reml$summary$vpc$estimate, def$summary$vpc$estimate)
  # The results table says which basis its PCV row is on; the default needs no note.
  note <- maihda_table(reml)$models_note
  expect_match(note, "REML-type between-stratum variances (estimation = \"REML\"",
               fixed = TRUE)
  expect_match(note, "own maximum-likelihood estimate", fixed = TRUE)
  expect_false(grepl("maximum-likelihood-refitted", note, fixed = TRUE))
  expect_null(maihda_table(def)$models_note)

  # Refused before any model is fitted.
  expect_error(maihda(y ~ a + b + (1 | a:b), data = rows, family = "binomial",
                      group = "c", estimation = "REML"),
               "not available with 'group'", fixed = TRUE)
  rows$id <- seq_len(nrow(rows))
  rows$wave <- rep_len(0:2, nrow(rows))
  expect_error(maihda(y ~ wave + (1 | a:b), data = rows, id = "id", time = "wave",
                      estimation = "REML"),
               "not available for a longitudinal MAIHDA", fixed = TRUE)
  expect_error(MAIHDA:::maihda_reml_workflow_check("ordinal", "two-model", NULL),
               "not available for engine = \"ordinal\"", fixed = TRUE)
  expect_error(MAIHDA:::maihda_reml_workflow_check("wemix", "two-model", NULL),
               "not available for engine = \"wemix\"", fixed = TRUE)
  # The crossed-dimensions decomposition computes no PCV, so there is nothing to refuse.
  expect_true(MAIHDA:::maihda_reml_workflow_check("lme4", "crossed-dimensions", "c"))
  # Refused after the first fit, before the bootstrap summaries.
  expect_error(audit_quiet(maihda(y ~ a + b + c + (1 | a:b:c), data = rows,
                                  family = "binomial", estimation = "REML",
                                  bootstrap = TRUE, n_boot = 200)),
               "bootstrap = TRUE is not available with estimation = \"REML\"",
               fixed = TRUE)
})
