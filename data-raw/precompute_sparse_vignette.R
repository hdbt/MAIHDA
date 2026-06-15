# data-raw/precompute_sparse_vignette.R
# Precompute the brms (and reference lme4) fits for the "Bayesian MAIHDA for sparse
# intersections" vignette, and cache *lightweight* summaries to
# vignettes/sparse_precomputed.rds.
#
# WHY: brms requires Stan + a compiler and takes minutes; CRAN's and pkgdown's
# builders cannot run it. The vignette therefore shows the brm() call in a non-
# evaluated chunk and renders its results from this cache. We store only small
# numbers / data frames -- NOT the multi-MB brmsfit objects -- to keep the package
# light.
#
# Re-run from the package root whenever maihda_sparse_data or the model settings
# change:  Rscript data-raw/precompute_sparse_vignette.R

stopifnot(requireNamespace("brms", quietly = TRUE))
suppressMessages(devtools::load_all(quiet = TRUE))
library(brms)

load("data/maihda_sparse_data.rda")
d  <- maihda_sparse_data
tr <- attr(d, "truth")
form <- ~ 1 + (1 | gender:ethnicity:education:age_group)

## ---- helpers ----------------------------------------------------------------
brms_diag <- function(maihda_obj) {
  bf <- maihda_obj$model$model
  np <- brms::nuts_params(bf)
  ndraws <- (length(unique(np$Iteration))) * length(unique(np$Chain))
  list(
    max_rhat   = max(brms::rhat(bf), na.rm = TRUE),
    min_ess    = min(brms::neff_ratio(bf), na.rm = TRUE) * ndraws,
    divergences = sum(np$Parameter == "divergent__" & np$Value == 1)
  )
}

# Build a fit's summary entry explicitly (avoids c()'s list-splicing name mangling).
pack <- function(maihda_obj, with_diag = FALSE) {
  dec <- maihda_obj$decomposition
  out <- list(share = dec$interaction_share,
              ci = dec$interaction_share_ci,
              additive_share = dec$additive_share,
              singular = isTRUE(maihda_obj$model$diagnostics$singular))
  if (with_diag) out$diag <- brms_diag(maihda_obj)
  out
}

## ---- lme4 reference fits (fast, deterministic) ------------------------------
lme4_g <- maihda(stats::update(form, y ~ .), data = d,
                 decomposition = "crossed-dimensions", engine = "lme4")
lme4_b <- maihda(stats::update(form, event ~ .), data = d,
                 decomposition = "crossed-dimensions", engine = "lme4", family = "binomial")

## ---- brms fits --------------------------------------------------------------
brms_g <- maihda(
  stats::update(form, y ~ .), data = d,
  decomposition = "crossed-dimensions", engine = "brms",
  prior = set_prior("normal(0, 0.5)", class = "sd"),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 1,
  refresh = 0, silent = 2, control = list(adapt_delta = 0.97, max_treedepth = 12)
)
brms_b <- maihda(
  stats::update(form, event ~ .), data = d,
  decomposition = "crossed-dimensions", engine = "brms", family = "binomial",
  prior = set_prior("normal(0, 1)", class = "sd"),
  chains = 4, iter = 2500, warmup = 1000, cores = 4, seed = 1,
  refresh = 0, silent = 2, control = list(adapt_delta = 0.99, max_treedepth = 12)
)

## ---- assemble lightweight cache --------------------------------------------
precomp <- list(
  truth = tr,
  lme4  = list(
    gaussian = pack(lme4_g),
    binary   = pack(lme4_b)
  ),
  brms = list(
    gaussian = pack(brms_g, with_diag = TRUE),
    binary   = pack(brms_b, with_diag = TRUE)
  ),
  prior = list(gaussian = "normal(0, 0.5)", binary = "normal(0, 1)"),
  settings = list(gaussian = "4 chains x 2000, adapt_delta = 0.97",
                  binary   = "4 chains x 2500, adapt_delta = 0.99")
)

saveRDS(precomp, "vignettes/sparse_precomputed.rds", version = 2)

cat("\n=== cached summary ===\n")
cat(sprintf("TRUE interaction share: gaussian=%.1f%%  binary=%.1f%%\n",
            100 * tr$gaussian$interaction_share, 100 * tr$binary_latent$interaction_share))
cat(sprintf("lme4 gaussian=%.1f%% (sing=%s) | binary=%.1f%% (sing=%s)\n",
            100 * precomp$lme4$gaussian$share, precomp$lme4$gaussian$singular,
            100 * precomp$lme4$binary$share, precomp$lme4$binary$singular))
cat(sprintf("brms gaussian=%.1f%% CI[%.1f,%.1f] rhat=%.3f div=%d\n",
            100 * precomp$brms$gaussian$share, 100 * precomp$brms$gaussian$ci[1],
            100 * precomp$brms$gaussian$ci[2], precomp$brms$gaussian$diag$max_rhat,
            precomp$brms$gaussian$diag$divergences))
cat(sprintf("brms binary  =%.1f%% CI[%.1f,%.1f] rhat=%.3f div=%d\n",
            100 * precomp$brms$binary$share, 100 * precomp$brms$binary$ci[1],
            100 * precomp$brms$binary$ci[2], precomp$brms$binary$diag$max_rhat,
            precomp$brms$binary$diag$divergences))
cat("saved vignettes/sparse_precomputed.rds\n")
