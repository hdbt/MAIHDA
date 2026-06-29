# Precompute the UNWEIGHTED full-data MAIHDA for the vignette.
#
# The vignette's analysis is the unweighted MAIHDA on all ~352,714 complete-case
# records. It is fitted ONCE here and the results/figures are cached for the article
# to display (the vignette shows the individual-level call but does not evaluate it).
#
# HOW IT IS FITTED: on the grouped one-row-per-stratum table, cbind(cases, controls)
# ~ dims + (1 | stratum). Because the outcome is binomial and every predictor is a
# stratum-defining dimension, this is EXACTLY the unweighted individual-level fit on
# all 352,714 records -- the cell counts are sufficient statistics, so nothing is
# lost (verified: identical VPC/PCV/BLUPs/ROPE to the individual fit, to ~1e-8). We
# use the grouped form purely because it is instant; the literal 352k-row individual
# glmer takes ~10 min for identical results. (Historical note: that individual fit
# used to return AUC = NA via an integer overflow in maihda_auc's n1 * n0; fixed in
# R/discriminatory_accuracy.R, with a regression test, so either form now gives the
# same AUC.)
#
# Outputs:
#   vignettes/brfss_precomputed.rds   -- small cache the vignette reads
#   vignettes/figures/brfss_*.png     -- the three case-study figures
#
# Run from the package root:  Rscript data-raw/brfss_precompute.R

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE)
  else library(MAIHDA)
  library(dplyr); library(ggplot2)
})

run_brfss_precompute <- function(data_dir = file.path("data", "brfss-2024"),
                                 vignette_dir = "vignettes",
                                 figures_dir = file.path("vignettes", "figures")) {
  analytic_path <- file.path(data_dir, "brfss_reduced_analytic.rds")
  if (!file.exists(analytic_path)) {
    stop("Analytic individual file not found: ", analytic_path,
         "\nRun the recode in data-raw/brfss_reduced_case_study.R first.", call. = FALSE)
  }
  sv <- c("sex", "race_ethnicity", "age_group", "education", "income", "disability")
  d <- as.data.frame(readRDS(analytic_path))[, c("frequent_distress", sv)]
  n_individuals <- nrow(d)

  # Aggregate the full individual data to its sufficient statistics (one row per
  # stratum). This loses nothing for a binomial model with stratum-level predictors.
  strata <- d |>
    group_by(across(all_of(sv))) |>
    summarise(raw_n = n(), raw_cases = sum(frequent_distress),
              raw_controls = raw_n - raw_cases, raw_prevalence = mean(frequent_distress),
              .groups = "drop")
  strata$label <- do.call(paste, c(strata[sv], sep = " × "))
  message("Fitting unweighted MAIHDA on ", format(n_individuals, big.mark = ","),
          " records (", nrow(strata), " strata) ...")

  fit <- suppressWarnings(maihda(
    cbind(raw_cases, raw_controls) ~
      sex + race_ethnicity + age_group + education + income + disability +
      (1 | sex:race_ethnicity:age_group:education:income:disability),
    data = strata, family = "binomial", interactions = "BH"))

  g  <- as.data.frame(generics::glance(fit))
  tab <- maihda_table(fit)
  mt <- tab$models
  mstat <- function(s, col) { v <- mt[mt$statistic == s, col]; if (length(v)) as.numeric(v[1]) else NA_real_ }

  # Observed cell sizes / prevalences for the ranked table, attached from the
  # aggregated table by the stratum label (maihda_table uses the same " × " labels).
  ranked <- tab$strata |>
    left_join(dplyr::select(strata, label, raw_n, observed = raw_prevalence), by = "label")
  rcols <- intersect(c("rank", "label", "predicted", "predicted_lower",
                       "predicted_upper", "observed", "raw_n"), names(ranked))
  ranked_small <- ranked[, rcols, drop = FALSE]

  rope_summ <- function(dd) {
    ir <- maihda_interactions(fit, rope = dd)
    tb <- table(factor(ir$decision, c("relevant", "negligible", "inconclusive")))
    list(d = dd, or = exp(dd), n = nrow(ir),
         relevant = as.integer(tb[["relevant"]]),
         negligible = as.integer(tb[["negligible"]]),
         inconclusive = as.integer(tb[["inconclusive"]]), ints = ir)
  }
  rp <- rope_summ(0.4); rs <- rope_summ(0.2)
  inter <- as.data.frame(rp$ints)
  inter <- inter[order(-abs(inter$interaction)), , drop = FALSE]
  icols <- intersect(c("label", "interaction", "lower", "upper", "p_adjusted",
                       "direction", "flagged", "decision"), names(inter))

  pc <- list(
    analytic_n = n_individuals, n_strata = nrow(ranked),
    n_sparse_lt20 = sum(strata$raw_n < 20), n_sparse_lt50 = sum(strata$raw_n < 50),
    vpc = as.numeric(g$vpc), vpc_adjusted = mstat("VPC/ICC", "adjusted"),
    pcv = as.numeric(g$pcv), auc = as.numeric(g$auc),
    auc_adjusted = as.numeric(g$auc.adjusted), mor = as.numeric(g$mor),
    mor_adjusted = mstat("MOR", "adjusted"),
    between_var_null = mstat("Between-stratum variance", "null"),
    between_var_adjusted = mstat("Between-stratum variance", "adjusted"),
    n_flagged = sum(fit$interactions$flagged, na.rm = TRUE),
    n_interactions = nrow(fit$interactions),
    max_abs_int = max(abs(inter$interaction), na.rm = TRUE),
    rope_primary = rp[c("d", "or", "n", "relevant", "negligible", "inconclusive")],
    rope_strict  = rs[c("d", "or", "n", "relevant", "negligible", "inconclusive")],
    glance = g, model_table = mt,
    top_strata = utils::head(ranked_small, 10),
    bottom_strata = utils::tail(ranked_small, 10),
    top_interactions = utils::head(inter[, icols, drop = FALSE], 20)
  )

  cache_path <- file.path(vignette_dir, "brfss_precomputed.rds")
  saveRDS(pc, cache_path, version = 2)
  message("Wrote ", cache_path, " (", format(file.size(cache_path), big.mark = ","), " bytes); ",
          sprintf("VPC=%.3f PCV=%.3f AUC=%.3f MOR=%.3f flagged=%d ROPE0.4 rel/neg/inc=%d/%d/%d",
                  pc$vpc, pc$pcv, pc$auc, pc$mor, pc$n_flagged,
                  pc$rope_primary$relevant, pc$rope_primary$negligible, pc$rope_primary$inconclusive))

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

if (sys.nframe() == 0L) run_brfss_precompute()
