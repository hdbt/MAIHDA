# Audit 2026-10-19 (sparse-vignette-lme4-boundary): vignette("bayesian_sparse_maihda")
# said that both crossed-dimensions lme4 fits of maihda_sparse_data are singular because
# maximum likelihood "has placed the interaction variance at the boundary" -- CONFIRMED
# wrong. The variance at zero is a dimension's: gender's in the Gaussian fit (4.0e-10),
# gender's (0) and ethnicity's (1.3e-10) in the binary fit. The interaction variance is
# 0.105 and 0.022, interaction shares of 23.4% and 3.1% for a true 40% (of the
# variance of the 36 stratum effects the generator drew, 33.6% is non-additive). The
# vignette now says which variance is at zero and that the interaction share is
# understated, and no longer reads the lme4 point as possibly "none at all". Words
# only; no number changes.
#
# No block draws a random number: data() and maihda() with lme4 draw nothing.

test_that("the sparse vignette's lme4 fits have a dimension at zero, not the interaction", {
  skip_on_cran()
  data(maihda_sparse_data, envir = environment())
  quiet <- function(expr) suppressMessages(suppressWarnings(expr))
  # The vignette's `lme4` chunk.
  m_g <- quiet(maihda(y ~ 1 + (1 | gender:ethnicity:education:age_group),
                      data = maihda_sparse_data, decomposition = "crossed-dimensions",
                      engine = "lme4"))
  m_b <- quiet(maihda(event ~ 1 + (1 | gender:ethnicity:education:age_group),
                      data = maihda_sparse_data, decomposition = "crossed-dimensions",
                      engine = "lme4", family = "binomial"))
  truth <- attr(maihda_sparse_data, "truth")

  # The named dimensions' variances are at zero and no other term is. The REML optimum
  # of the Gaussian fit has gender's variance exactly at zero, but whether lme4 calls
  # the fit singular depends on where its optimiser stops near it (a Cholesky factor
  # below 1e-4; the row order alone moves the stop up to about 2e-3), so the variances
  # are read, and the record only for the terms it names.
  check <- function(m, zero, true_share) {
    dg <- m$model$diagnostics
    dec <- m$decomposition
    expect_true(all(dec$per_dim[zero] < 1e-4))
    expect_true(all(dec$per_dim[setdiff(names(dec$per_dim), zero)] > 0.01))
    expect_true(all(dg$singular_terms %in% zero))
    # The interaction variance is away from zero; its share is below the truth.
    expect_gt(dec$interaction_var, 0.01)
    expect_lt(dec$interaction_share, true_share)
  }
  check(m_g, "gender", truth$gaussian$interaction_share)
  check(m_b, c("gender", "ethnicity"), truth$binary_latent$interaction_share)
  # The binary fit's two zeros sit well inside lme4's threshold.
  expect_true(isTRUE(m_b$model$diagnostics$singular))
  expect_setequal(m_b$model$diagnostics$singular_terms, c("gender", "ethnicity"))

  # The vignette's text quotes the shares from its cache: they are these fits' shares.
  cache <- file.path(testthat::test_path("..", ".."), "vignettes", "sparse_precomputed.rds")
  if (file.exists(cache)) {
    pc <- readRDS(cache)
    pct <- function(x) sprintf("%.0f%%", 100 * x)
    expect_identical(pct(pc$lme4$gaussian$share), pct(m_g$decomposition$interaction_share))
    expect_identical(pct(pc$lme4$binary$share), pct(m_b$decomposition$interaction_share))
  }
})

test_that("the sparse vignette says which variance the lme4 fits put at zero", {
  skip_on_cran()
  root <- testthat::test_path("..", "..")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  p <- file.path(root, "vignettes", "bayesian_sparse_maihda.Rmd")
  skip_if_not(file.exists(p), "the package tree is not available")
  v <- gsub("[[:space:]]+", " ", paste(readLines(p, warn = FALSE), collapse = " "))
  for (old in c("has placed the interaction variance at",
                "interaction share is pulled down",
                "or none at all")) {
    expect_false(grepl(old, v, fixed = TRUE), info = old)
  }
  expect_true(grepl(paste0(
    "Both fits are singular: maximum likelihood has placed a dimension's variance at ",
    "zero (gender's in the Gaussian fit; gender's and ethnicity's in the binary one), ",
    "and the interaction variance is positive in both."), v, fixed = TRUE))
  expect_true(grepl("The interaction share is understated all the same", v, fixed = TRUE))
})
