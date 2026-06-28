# Precompute the design-weighted (WeMix) BRFSS MAIHDA for the vignette.
#
# The vignette's headline analysis is the design-based fit:
#   maihda(frequent_distress ~ <dims> + (1 | <strata>),
#          family = "binomial", sampling_weights = "survey_weight")
# which uses WeMix (design-consistent pseudo-maximum-likelihood). WeMix's adaptive
# quadrature over the 432 intersectional strata makes a single fit take minutes and
# the full two-model maihda() far longer -- impossible at package-build time. So the
# fit is run here ONCE on the full individual data and its results are cached for the
# article to display (the vignette shows the real call but does not evaluate it).
#
# Outputs:
#   vignettes/brfss_precomputed.rds   -- small cache the vignette reads (numbers + tables)
#   vignettes/figures/brfss_*.png     -- the three case-study figures
#
# Run from the package root (full data):
#   Rscript data-raw/brfss_wemix_precompute.R
# Smoke test on a subsample (fast, for validating the pipeline only -- NOT the
# published numbers):
#   Rscript data-raw/brfss_wemix_precompute.R 2000 1
#
# Quadrature: defaults to nQuad = 7 with fast = TRUE (WeMix's C++ path). WeMix's
# own defaults are nQuad = 13, fast = FALSE; for this single random-intercept model
# adaptive quadrature is essentially converged well below nQuad = 13, so nQuad = 7
# should match the nQuad = 13 estimates very closely while running several times
# faster. The originally shipped cache was built with WeMix's full default
# (nQuad = 13, fast = FALSE); pass nquad = NA to reproduce that exact (slower)
# setting. If you regenerate at nQuad = 7, sanity-check the headline numbers against
# the shipped values before relying on the refreshed cache.

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(MAIHDA)
  }
  library(dplyr); library(ggplot2)
})

run_brfss_wemix_precompute <- function(n_sub = NA_integer_, nquad = 7L,
                                       data_dir = file.path("data", "brfss-2024"),
                                       vignette_dir = "vignettes",
                                       figures_dir = file.path("vignettes", "figures")) {
  analytic_path <- file.path(data_dir, "brfss_reduced_analytic.rds")
  if (!file.exists(analytic_path)) {
    stop("Analytic individual file not found: ", analytic_path,
         "\nRun data-raw/brfss_reduced_case_study.R first.", call. = FALSE)
  }
  d <- as.data.frame(readRDS(analytic_path))
  if (!is.na(n_sub) && n_sub < nrow(d)) {
    set.seed(2024); d <- d[sample.int(nrow(d), n_sub), , drop = FALSE]
    d[] <- lapply(d, function(x) if (is.factor(x)) droplevels(x) else x)
    message("SMOKE TEST: subsampled to ", nrow(d), " rows (numbers are NOT publishable)")
  }

  f <- frequent_distress ~ sex + race_ethnicity + age_group + education + income + disability +
    (1 | sex:race_ethnicity:age_group:education:income:disability)
  dots <- list(formula = f, data = d, family = "binomial",
               sampling_weights = "survey_weight", interactions = "BH")
  if (!is.na(nquad)) { dots$nQuad <- nquad; dots$fast <- TRUE }

  message("Fitting WeMix maihda() on ", format(nrow(d), big.mark = ","),
          " rows", if (!is.na(nquad)) paste0(" (nQuad=", nquad, ")") else "", " ...")
  t0 <- Sys.time()
  fit <- suppressWarnings(do.call(maihda, dots))
  message(sprintf("  fit done in %.1f min (engine=%s)",
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), fit$model$engine))

  g  <- as.data.frame(generics::glance(fit))
  tab <- maihda_table(fit)
  mt <- tab$models
  mstat <- function(s, col) {
    v <- mt[mt$statistic == s, col]; if (length(v) == 0) NA_real_ else as.numeric(v[1])
  }

  # Per-stratum observed counts from the (null model's) individual analytic frame,
  # joined to the ranked table by the stratum id.
  md <- fit$model$data
  ranked <- tab$strata
  if (!is.null(md) && "stratum" %in% names(md) && "stratum" %in% names(ranked)) {
    counts <- md |>
      mutate(stratum = as.character(stratum)) |>
      group_by(stratum) |>
      summarise(raw_n = n(),
                raw_cases = sum(frequent_distress),
                raw_prevalence = mean(frequent_distress),
                weighted_prevalence = stats::weighted.mean(frequent_distress, survey_weight),
                .groups = "drop")
    ranked <- ranked |>
      mutate(stratum = as.character(stratum)) |>
      left_join(counts, by = "stratum")
  } else {
    ranked$raw_n <- NA_integer_; ranked$weighted_prevalence <- NA_real_
  }
  rcols <- intersect(c("rank", "label", "predicted", "predicted_lower",
                       "predicted_upper", "weighted_prevalence", "raw_n"), names(ranked))
  ranked_small <- ranked[, rcols, drop = FALSE]

  inter <- as.data.frame(fit$interactions)
  inter <- inter[order(-abs(inter$interaction)), , drop = FALSE]
  icols <- intersect(c("label", "interaction", "lower", "upper", "p_adjusted",
                       "direction", "flagged"), names(inter))

  # Sparsity counts from the full per-stratum sizes.
  sizes <- if (exists("counts")) counts$raw_n else NA_integer_

  pc <- list(
    smoke = !is.na(n_sub),
    analytic_n = nrow(d),
    n_strata = nrow(ranked),
    n_sparse_lt20 = if (all(!is.na(sizes))) sum(sizes < 20) else NA_integer_,
    n_sparse_lt50 = if (all(!is.na(sizes))) sum(sizes < 50) else NA_integer_,
    vpc = as.numeric(g$vpc), vpc_adjusted = mstat("VPC/ICC", "adjusted"),
    pcv = as.numeric(g$pcv), auc = as.numeric(g$auc),
    auc_adjusted = as.numeric(g$auc.adjusted), mor = as.numeric(g$mor),
    mor_adjusted = mstat("MOR", "adjusted"),
    between_var_null = mstat("Between-stratum variance", "null"),
    between_var_adjusted = mstat("Between-stratum variance", "adjusted"),
    n_flagged = sum(inter$flagged, na.rm = TRUE), n_interactions = nrow(inter),
    glance = g, model_table = mt,
    top_strata = utils::head(ranked_small, 10),
    bottom_strata = utils::tail(ranked_small, 10),
    top_interactions = utils::head(inter[, icols, drop = FALSE], 20)
  )

  cache_path <- file.path(vignette_dir, "brfss_precomputed.rds")
  saveRDS(pc, cache_path, version = 2)
  message("Wrote ", cache_path, " (", format(file.size(cache_path), big.mark = ","), " bytes)")

  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  save_fig <- function(p, name, h = 4.8) {
    ggplot2::ggsave(file.path(figures_dir, name), plot = p, width = 7, height = h, dpi = 150)
    message("Wrote ", file.path(figures_dir, name))
  }
  save_fig(plot(fit, type = "vpc"), "brfss_vpc.png", h = 4.2)
  save_fig(plot(fit, type = "predicted", n_strata = 30, highlight_interactions = TRUE),
           "brfss_predicted_top30.png", h = 6.5)
  save_fig(plot(fit, type = "effect_decomp", highlight_interactions = TRUE),
           "brfss_effect_decomp.png", h = 6)

  message("DONE.")
  invisible(pc)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  n_sub <- if (length(a) >= 1) as.integer(a[1]) else NA_integer_
  nquad <- if (length(a) >= 2) as.integer(a[2]) else 7L
  run_brfss_wemix_precompute(n_sub = n_sub, nquad = nquad)
}
