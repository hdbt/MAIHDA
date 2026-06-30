# Precompute the UNWEIGHTED full-data MAIHDA for the vignette.
#
# The vignette's analysis is the unweighted MAIHDA on all ~352,714 complete-case
# records. It is fitted ONCE here and the results/figures are cached for the article
# to display (the vignette shows the individual-level call but does not evaluate it).
#
# HOW IT IS FITTED: on the full INDIVIDUAL records,
# frequent_distress ~ dims + (1 | stratum). This is slow (~10 min) but it carries the
# true per-stratum sample sizes, which the UpSet plot's intersection-size bars need.
# (The equivalent grouped cbind(cases, controls) fit is instant and gives identical
# VPC/PCV/BLUPs/ROPE/AUC -- the cell counts are sufficient statistics -- but it reports
# a stratum "size" of 1 for every cell, so its UpSet size bars come out all equal.
# Hence the individual fit here. AUC on the full individual fit relies on the integer-
# overflow fix in R/discriminatory_accuracy.R.)
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

  message("Fitting unweighted MAIHDA on ", format(n_individuals, big.mark = ","),
          " individual records (slow ~10 min) ...")
  fit <- suppressWarnings(maihda(
    frequent_distress ~
      sex + race_ethnicity + age_group + education + income + disability +
      (1 | sex:race_ethnicity:age_group:education:income:disability),
    data = d, family = "binomial", interactions = "BH"))

  g  <- as.data.frame(generics::glance(fit))
  tab <- maihda_table(fit)
  mt <- tab$models
  mstat <- function(s, col) { v <- mt[mt$statistic == s, col]; if (length(v)) as.numeric(v[1]) else NA_real_ }

  # Observed cell sizes / prevalences for the ranked table + sparsity counts,
  # aggregated from the (individual) model frame by stratum id.
  counts <- fit$model$data |>
    mutate(stratum = as.character(stratum)) |>
    group_by(stratum) |>
    summarise(raw_n = n(), observed = mean(frequent_distress), .groups = "drop")
  ranked <- tab$strata |>
    mutate(stratum = as.character(stratum)) |>
    left_join(counts, by = "stratum")
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
    n_sparse_lt20 = sum(counts$raw_n < 20), n_sparse_lt50 = sum(counts$raw_n < 50),
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
  save_fig <- function(p, name, h = 4.8, w = 7) {
    ggplot2::ggsave(file.path(figures_dir, name), plot = p, width = w, height = h, dpi = 150)
    message("Wrote ", file.path(figures_dir, name))
  }
  # Highlight by the ROPE equivalence decision (the practically-relevant strata at
  # rope = 0.4), not the zero-centred FDR flag -- the substantive "which intersections
  # matter" set rather than the much larger "differs from zero" set.
  save_fig(plot(fit, type = "vpc"), "brfss_vpc.png", h = 4.2)
  # UpSet-style alternative to the text "predicted" view: a category matrix replaces
  # the long intersectional axis labels (legible at 432 strata).
  save_fig(plot(fit, type = "upset", n_strata = 20, select = "deviation",
                highlight_interactions = TRUE, highlight_by = "rope", rope = 0.4),
           "brfss_upset.png", h = 8.5, w = 8)
  save_fig(plot(fit, type = "effect_decomp",
                highlight_interactions = TRUE, highlight_by = "rope", rope = 0.4),
           "brfss_effect_decomp.png", h = 6)
  message("DONE.")
  invisible(pc)
}

if (sys.nframe() == 0L) run_brfss_precompute()
