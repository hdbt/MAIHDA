# Audit 2026-10-19 (sparse-realised-truth): attr(maihda_sparse_data, "truth") held the
# generator's targets -- a between-stratum variance of 0.35 (VPC 0.26) on y and 1.5
# (latent VPC 0.31) on event, an interaction share of 40% -- and ?maihda_sparse_data
# quoted them as the truth. data-raw/maihda_sparse_data.R standardises its additive and
# interaction parts separately, and over 36 strata they come out correlated (-0.12):
# the 36 drawn stratum effects have a variance of 0.308 (VPC 0.235) on y and 1.319
# (0.286) on the latent scale of event. Their interaction share is 40.6% by the
# mean-square split of the realised-strata block (33.6% by sums of squares). The truth
# attribute now holds what the drawn effects carry, with the effects themselves.
#
# The vignette held lme4's interaction share against that truth and read the whole gap
# as the sparse cells. lme4 reports a ratio of variance parameters, and the dimensions'
# variances, each estimated from two or three levels, add up to more than the spread
# of the main effects among the strata: with 1,000 individuals in every stratum lme4
# reads 27.6% (y, REML; noise-free limit 28.3%) and 34.7% (event, ML; limit 35.1%) on
# the same effects, whose two sets differ only in scale. The vignette now gives those
# shares, from its cache.
#
# The first block rebuilds the generator under its own seed and puts the caller's
# random-number state back; the others draw nothing.

# Evaluate `expr` under a fixed seed and restore the caller's RNG state, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
srt_with_seed <- function(seed, expr) {
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

# The draws of data-raw/maihda_sparse_data.R, in its order and arithmetic.
srt_generate <- function() {
  gender    <- c("Women", "Men")
  ethnicity <- c("White", "Black", "Asian")
  education <- c("Low", "High")
  age_group <- c("Young", "Mid", "Older")
  grid <- expand.grid(gender = gender, ethnicity = ethnicity,
                      education = education, age_group = age_group,
                      stringsAsFactors = FALSE)
  S <- nrow(grid)
  key <- apply(grid, 1, paste, collapse = ":")
  dim_fx <- function(levels, sd) stats::setNames(stats::rnorm(length(levels), 0, sd), levels)
  b_gender <- dim_fx(gender, 0.30); b_eth <- dim_fx(ethnicity, 0.30)
  b_edu    <- dim_fx(education, 0.30); b_age <- dim_fx(age_group, 0.30)
  raw_add  <- b_gender[grid$gender] + b_eth[grid$ethnicity] +
              b_edu[grid$education] + b_age[grid$age_group]
  raw_int  <- stats::rnorm(S, 0, 0.45)
  std <- function(x) (x - mean(x)) / stats::sd(x)
  add0 <- std(raw_add)
  int0 <- std(raw_int)
  target_share <- 0.40
  between_g <- 0.35; between_b <- 1.5
  g_stratum <- 50 + sqrt(between_g * (1 - target_share)) * add0 +
    sqrt(between_g * target_share) * int0
  b_eta <- 0.0 + sqrt(between_b * (1 - target_share)) * add0 +
    sqrt(between_b * target_share) * int0
  names(g_stratum) <- names(b_eta) <- key
  N <- 240L
  p <- stats::runif(S, 0.25, 1)^1.5
  p <- p / sum(p)
  who <- sample(seq_len(S), N, replace = TRUE, prob = p)
  d <- data.frame(
    gender    = factor(grid$gender[who],    levels = gender),
    ethnicity = factor(grid$ethnicity[who], levels = ethnicity),
    education = factor(grid$education[who],  levels = education),
    age_group = factor(grid$age_group[who], levels = age_group),
    y         = g_stratum[who] + stats::rnorm(N, 0, 1.0),
    event     = factor(ifelse(stats::runif(N) < stats::plogis(b_eta[who]), "Yes", "No"),
                       levels = c("No", "Yes")),
    stringsAsFactors = FALSE
  )
  rownames(d) <- NULL
  list(data = d, grid = grid, y = g_stratum - 50, event = b_eta - 0.0)
}

# The share an lme4 crossed-dimensions fit converges to as every stratum grows, for
# noise-free effects on the full 2 x 3 x 2 x 3 grid: the REML (or ML) estimate from the
# likelihood of the 36 effects, e ~ N(mu 1, sum_d s_d Z_d Z_d' + s_I I).
srt_limit_share <- function(effects, grid, REML) {
  Z <- lapply(names(grid), function(k) {
    m <- stats::model.matrix(~ 0 + factor(grid[[k]]))
    m %*% t(m)
  })
  J <- length(effects); one <- rep(1, J)
  nll <- function(lsd) {
    v <- exp(2 * pmin(pmax(lsd, -12), 4))
    V <- v[5] * diag(J)
    for (i in 1:4) V <- V + v[i] * Z[[i]]
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    solve_v <- function(x) backsolve(R, forwardsolve(t(R), x))
    vi1 <- solve_v(one)
    r <- effects - sum(solve_v(effects)) / sum(vi1)
    (2 * sum(log(diag(R))) + sum(r * solve_v(r)) + if (REML) log(sum(vi1)) else 0) / 2
  }
  fits <- lapply(list(rep(log(0.3), 5), c(rep(log(0.1), 4), log(0.3))), function(start)
    stats::optim(start, nll, method = "BFGS", control = list(reltol = 1e-12, maxit = 500)))
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  v <- exp(2 * pmin(pmax(best$par, -12), 4))
  v[5] / sum(v)
}

test_that("the truth attribute holds what the 36 drawn stratum effects carry", {
  data(maihda_sparse_data, envir = environment())
  tr <- attr(maihda_sparse_data, "truth")
  rebuilt <- srt_with_seed(2026, srt_generate())

  # The shipped rows and stored effects are the generator's (to rounding, which can
  # differ in the last bit on another platform).
  plain <- maihda_sparse_data
  attr(plain, "truth") <- NULL
  expect_equal(plain, rebuilt$data, tolerance = 1e-10)
  expect_identical(names(tr$gaussian$stratum_effects), names(rebuilt$y))
  expect_identical(names(tr$binary_latent$stratum_effects), names(rebuilt$event))
  expect_equal(unname(tr$gaussian$stratum_effects), unname(rebuilt$y), tolerance = 1e-10)
  expect_equal(unname(tr$binary_latent$stratum_effects), unname(rebuilt$event),
               tolerance = 1e-10)
  expect_identical(tr$gaussian$intercept, 50)
  expect_identical(tr$binary_latent$intercept, 0)

  # The realised split of summary()'s realised-strata block, the effects as one draw.
  cells <- rebuilt$grid
  for (o in c("gaussian", "binary_latent")) {
    eff <- unname(tr[[o]]$stratum_effects)
    U <- lapply(cells, function(level) matrix(tapply(eff, level, mean), 1L))
    rr <- maihda_cc_realised_from_draws(Z = matrix(eff, 1L), U = U, cells = cells,
                                        within = tr[[o]]$residual_var)
    expect_equal(rr$df_interaction, 29L, info = o)
    expect_equal(tr[[o]]$between_var, rr$between_var, tolerance = 1e-10, info = o)
    expect_equal(tr[[o]]$interaction_var, rr$interaction_var, tolerance = 1e-10, info = o)
    expect_equal(tr[[o]]$additive_var, rr$between_var - rr$interaction_var,
                 tolerance = 1e-10, info = o)
    expect_equal(tr[[o]]$interaction_share, rr$interaction_share, tolerance = 1e-10,
                 info = o)
    expect_equal(tr[[o]]$vpc, rr$vpc, tolerance = 1e-10, info = o)
  }
  expect_equal(tr$gaussian$residual_var, 1)
  expect_equal(tr$binary_latent$residual_var, pi^2 / 3)

  # Not the targets: the drawn effects fall short of 0.35 and 1.5, and the share is
  # the mean-square split, not the sum-of-squares one.
  expect_lt(tr$gaussian$between_var, 0.32)
  expect_lt(tr$binary_latent$between_var, 1.35)
  expect_gt(tr$gaussian$interaction_share, 0.40)
  expect_lt(tr$gaussian$interaction_share, 0.41)
  expect_equal(tr$target_interaction_share, 0.40)
})

test_that("the vignette's cache follows the data and holds lme4's share with ample data", {
  skip_on_cran()
  cache <- file.path(testthat::test_path("..", ".."), "vignettes", "sparse_precomputed.rds")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(file.exists(cache), "the package tree is not available")
  data(maihda_sparse_data, envir = environment())
  tr <- attr(maihda_sparse_data, "truth")
  pc <- readRDS(cache)
  expect_identical(pc$truth, tr)
  expect_equal(pc$lme4_large$n_per_stratum, 1000)

  # The fits with 1,000 individuals in every stratum read lme4's own large-sample share:
  # REML for the Gaussian outcome, ML (glmer) for the binary one.
  grid <- expand.grid(gender = c("Women", "Men"), ethnicity = c("White", "Black", "Asian"),
                      education = c("Low", "High"), age_group = c("Young", "Mid", "Older"),
                      stringsAsFactors = FALSE)
  expect_identical(names(tr$gaussian$stratum_effects), apply(grid, 1, paste, collapse = ":"))
  lim_g <- srt_limit_share(unname(tr$gaussian$stratum_effects), grid, REML = TRUE)
  lim_b <- srt_limit_share(unname(tr$binary_latent$stratum_effects), grid, REML = FALSE)
  expect_lt(abs(lim_g - 0.2829), 0.001)
  expect_lt(abs(lim_b - 0.3509), 0.001)
  expect_lt(abs(pc$lme4_large$gaussian$share - lim_g), 0.02)
  expect_lt(abs(pc$lme4_large$binary$share - lim_b), 0.02)

  # What the vignette says of them: below the truth; the sparse lme4 and brms points
  # below them; the binary lme4 point near zero; the brms intervals reach the truth.
  share <- tr$gaussian$interaction_share
  for (o in c("gaussian", "binary")) {
    expect_lt(pc$lme4_large[[o]]$share, share)
    expect_lt(pc$lme4[[o]]$share, pc$lme4_large[[o]]$share)
    expect_lt(pc$brms[[o]]$share, pc$lme4_large[[o]]$share)
    expect_gt(pc$brms[[o]]$ci[2], share)
  }
  expect_lt(pc$lme4$binary$share, 0.1)
})

test_that("the help page, the vignette and the data script quote the drawn truth", {
  skip_on_cran()
  root <- testthat::test_path("..", "..")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(file.path(root, "man")), "the package tree is not available")
  flat <- function(...) {
    p <- file.path(root, ...)
    if (!file.exists(p)) return(NA_character_)
    gsub("[[:space:]]+", " ", paste(readLines(p, warn = FALSE), collapse = " "))
  }
  data(maihda_sparse_data, envir = environment())
  tr <- attr(maihda_sparse_data, "truth")
  pct <- function(x) sprintf("%.0f\\%%", 100 * x)

  rd <- flat("man", "maihda_sparse_data.Rd")
  expect_true(grepl(sprintf("a variance of %.2f (a VPC of %.2f), of which %s is the",
                            tr$gaussian$between_var, tr$gaussian$vpc,
                            pct(tr$gaussian$interaction_share)), rd, fixed = TRUE))
  expect_true(grepl(sprintf("a variance of %.2f (a VPC of %.2f), again %s interaction",
                            tr$binary_latent$between_var, tr$binary_latent$vpc,
                            pct(tr$binary_latent$interaction_share)), rd, fixed = TRUE))
  for (old in c("VPC 0.26", "VPC is 0.31", "exact generative truth", "again 40\\%")) {
    expect_false(grepl(old, rd, fixed = TRUE), info = old)
  }
  # The shares differ by estimator, not by outcome: both sets of effects have one shape.
  expect_true(grepl("add up to more than the spread of the main effects", rd, fixed = TRUE))
  expect_true(grepl("the gap between the two being that of REML and maximum likelihood",
                    rd, fixed = TRUE))
  cache <- file.path(root, "vignettes", "sparse_precomputed.rds")
  if (file.exists(cache)) {
    pc <- readRDS(cache)
    expect_true(grepl(sprintf("about %s on \\code{y} and %s on \\code{event}, the gap",
                              pct(pc$lme4_large$gaussian$share),
                              pct(pc$lme4_large$binary$share)), rd, fixed = TRUE))
    expect_true(grepl(sprintf("about %s on \\code{y} and %s on \\code{event}. The brms",
                              pct(pc$lme4$gaussian$share), pct(pc$lme4$binary$share)),
                      rd, fixed = TRUE))
  }

  v <- flat("vignettes", "bayesian_sparse_maihda.Rmd")
  if (!is.na(v)) {
    for (old in c("true 40%", "the 40% either", "weakly identified at this density",
                  "constructed to account for", "neither engine can locate",
                  "# the true interaction share: 0.40")) {
      expect_false(grepl(old, v, fixed = TRUE), info = old)
    }
    expect_true(grepl("Only part of that gap is the sparse cells.", v, fixed = TRUE))
    expect_true(grepl("lme4 puts the share at `r pct(pc$lme4_large$gaussian$share)`", v,
                      fixed = TRUE))
    expect_true(grepl("geom_point(data = large, shape = 1", v, fixed = TRUE))
    expect_true(grepl("that of REML (Gaussian) and maximum likelihood (binary)", v,
                      fixed = TRUE))
    # The lme4 column holds a REML fit (Gaussian) beside the ML ones.
    expect_false(grepl("\"lme4 (ML)\"", v, fixed = TRUE))
  }

  raw <- flat("data-raw", "maihda_sparse_data.R")
  if (!is.na(raw)) {
    for (old in c("showcase", "falsely precise", "regularise the variance off the boundary",
                  "claim *recovery*", "realised true variances")) {
      expect_false(grepl(old, raw, fixed = TRUE), info = old)
    }
    expect_true(grepl("stratum_effects = effects", raw, fixed = TRUE))
  }
})
