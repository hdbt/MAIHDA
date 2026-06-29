# Full reproducible pipeline for the BRFSS 2024 case study: download the public-use
# file, build the collapsed analytic data, fit the grouped-binomial MAIHDA, and
# export results under data/brfss-2024/.
#
# This script is the SINGLE SOURCE OF TRUTH for the case study. The recode/strata
# logic in vignettes/real_data_brfss.Rmd mirrors the block below for illustration;
# if you change the recoding here, mirror it there. Running this script to the end
# also runs data-raw/brfss_precompute.R, which fits the unweighted MAIHDA on the
# full data and caches the numbers/figures the vignette shows
# (vignettes/brfss_precomputed.rds + vignettes/figures/brfss_*.png), so the rendered
# article cannot drift from this run. That full-data fit is the slow step (minutes).

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(MAIHDA)
  }
  library(dplyr)
  library(ggplot2)
  library(haven)
})

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", sprintf(...), "\n")
  flush.console()
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

safe_plot <- function(expr, path, width = 7, height = 4.8) {
  tryCatch({
    p <- force(expr)
    ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 160)
    TRUE
  }, error = function(e) {
    log_msg("Plot failed for %s: %s", basename(path), conditionMessage(e))
    FALSE
  })
}

export_analysis <- function(analysis, prefix, output_dir, strata_lookup,
                            highlight_interactions = TRUE) {
  log_msg("Exporting %s results", prefix)

  glance <- generics::glance(analysis)
  tab <- maihda_table(analysis)
  ranked <- tab$strata
  count_cols <- c("raw_n", "raw_cases", "raw_controls", "raw_prevalence",
                  "weighted_n", "weighted_cases", "weighted_controls",
                  "weighted_prevalence")
  counts_joined <- FALSE
  model_data <- tryCatch(analysis$model$data, error = function(e) NULL)
  if (!is.null(model_data) && "stratum" %in% names(model_data) &&
      "stratum" %in% names(ranked)) {
    model_count_cols <- intersect(
      c("stratum", "raw_n", "raw_cases", "raw_controls", "raw_prevalence",
        "weighted_n", "weighted_cases", "weighted_controls",
        "weighted_prevalence"),
      names(model_data)
    )
    raw_counts <- model_data |>
      select(all_of(model_count_cols)) |>
      distinct()
    raw_counts$stratum <- as.character(raw_counts$stratum)
    ranked$stratum <- as.character(ranked$stratum)
    ranked <- dplyr::left_join(ranked, raw_counts, by = "stratum")
    counts_joined <- "raw_n" %in% names(ranked) && any(!is.na(ranked$raw_n))
  }
  if (!counts_joined) {
    ranked <- ranked |>
      select(-any_of(count_cols))
    ranked <- dplyr::left_join(ranked, strata_lookup, by = "label")
  }

  write_csv(glance, file.path(output_dir, paste0(prefix, "_glance.csv")))
  write_csv(tab$models, file.path(output_dir, paste0(prefix, "_model_table.csv")))
  write_csv(ranked, file.path(output_dir, paste0(prefix, "_ranked_strata.csv")))

  if (!is.null(analysis$interactions)) {
    int_cols <- intersect(
      c("stratum", "label", "interaction", "lower", "upper",
        "se", "p_value", "p_adjusted", "direction", "flagged"),
      names(analysis$interactions)
    )
    int_out <- analysis$interactions |>
      arrange(desc(abs(.data$interaction))) |>
      select(all_of(int_cols))
    write_csv(int_out, file.path(output_dir, paste0(prefix, "_interactions.csv")))
  }

  suffix <- if (isTRUE(highlight_interactions)) "_highlighted" else "_plain"
  safe_plot(
    plot(analysis, type = "vpc"),
    file.path(output_dir, paste0(prefix, "_vpc.png"))
  )
  safe_plot(
    plot(
      analysis,
      type = "predicted",
      n_strata = 30,
      highlight_interactions = highlight_interactions
    ),
    file.path(output_dir, paste0(prefix, "_predicted_top30", suffix, ".png")),
    height = 6
  )
  safe_plot(
    plot(
      analysis,
      type = "effect_decomp",
      highlight_interactions = highlight_interactions
    ),
    file.path(output_dir, paste0(prefix, "_effect_decomp", suffix, ".png")),
    height = 6
  )

  list(
    glance = glance,
    model_table = tab$models,
    ranked_strata = ranked,
    interactions = analysis$interactions
  )
}

data_dir <- file.path("data", "brfss-2024")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

highlight_interactions <- !tolower(Sys.getenv(
  "MAIHDA_BRFSS_HIGHLIGHT_INTERACTIONS", "true"
)) %in% c("false", "0", "no")

brfss_url <- "https://www.cdc.gov/brfss/annual_data/2024/files/LLCP2024XPT.zip"
zip_file <- file.path(data_dir, "LLCP2024XPT.zip")

if (!file.exists(zip_file)) {
  log_msg("Downloading BRFSS 2024 XPT zip")
  download.file(brfss_url, zip_file, mode = "wb", quiet = FALSE)
} else {
  log_msg("Using cached zip: %s", zip_file)
}

xpt_listing <- unzip(zip_file, list = TRUE)
xpt_name <- xpt_listing$Name[grepl("\\.xpt\\s*$", xpt_listing$Name, ignore.case = TRUE)][1]
xpt_file <- file.path(data_dir, trimws(xpt_name))

if (!file.exists(xpt_file)) {
  log_msg("Unzipping %s", trimws(xpt_name))
  unzip(zip_file, exdir = data_dir, overwrite = TRUE)
  extracted_xpt <- list.files(data_dir, full.names = TRUE)
  extracted_xpt <- extracted_xpt[grepl("\\.xpt\\s*$", basename(extracted_xpt), ignore.case = TRUE)]
  if (length(extracted_xpt) == 0L) {
    stop("No XPT file found after unzipping ", zip_file)
  }
  xpt_file <- extracted_xpt[1]
} else {
  log_msg("Using cached XPT: %s", xpt_file)
}

raw_cols <- c(
  "_MENT14D", "SEXVAR", "_IMPRACE", "_AGEG5YR", "EDUCA", "_INCOMG1",
  "DEAF", "BLIND", "DECIDE", "DIFFWALK", "DIFFDRES", "DIFFALON",
  "_STATE", "_LLCPWT"
)

log_msg("Reading only relevant BRFSS columns")
brfss_raw <- haven::read_xpt(xpt_file, col_select = dplyr::all_of(raw_cols))
log_msg("Reduced raw columns: %s rows x %s columns", nrow(brfss_raw), ncol(brfss_raw))

disability_vars <- c("DEAF", "BLIND", "DECIDE", "DIFFWALK", "DIFFDRES", "DIFFALON")
disability_items <- as.data.frame(lapply(brfss_raw[disability_vars], as.numeric))
any_disability <- rowSums(disability_items == 1, na.rm = TRUE) > 0
all_answered_no <- rowSums(disability_items == 2, na.rm = TRUE) == length(disability_vars)

log_msg("Creating collapsed analytic dataset")
brfss_reduced <- brfss_raw |>
  transmute(
    frequent_distress = case_when(
      .data[["_MENT14D"]] == 3 ~ 1L,
      .data[["_MENT14D"]] %in% c(1, 2) ~ 0L,
      TRUE ~ NA_integer_
    ),
    sex = factor(
      case_when(
        SEXVAR == 1 ~ "Male",
        SEXVAR == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      levels = c("Male", "Female")
    ),
    race_ethnicity = factor(
      case_when(
        .data[["_IMPRACE"]] == 1 ~ "White, non-Hispanic",
        .data[["_IMPRACE"]] == 2 ~ "Black, non-Hispanic",
        .data[["_IMPRACE"]] == 5 ~ "Hispanic",
        .data[["_IMPRACE"]] %in% c(3, 4, 6) ~ "Other race/ethnicity",
        TRUE ~ NA_character_
      ),
      levels = c("White, non-Hispanic", "Black, non-Hispanic",
                 "Hispanic", "Other race/ethnicity")
    ),
    age_group = factor(
      case_when(
        .data[["_AGEG5YR"]] %in% 1:3 ~ "18-34",
        .data[["_AGEG5YR"]] %in% 4:9 ~ "35-64",
        .data[["_AGEG5YR"]] %in% 10:13 ~ "65+",
        TRUE ~ NA_character_
      ),
      levels = c("18-34", "35-64", "65+")
    ),
    education = factor(
      case_when(
        EDUCA %in% 1:4 ~ "HS or less",
        EDUCA == 5 ~ "Some college",
        EDUCA == 6 ~ "College graduate",
        TRUE ~ NA_character_
      ),
      levels = c("HS or less", "Some college", "College graduate")
    ),
    income = factor(
      case_when(
        .data[["_INCOMG1"]] %in% 1:2 ~ "<$25k",
        .data[["_INCOMG1"]] %in% 3:4 ~ "$25k-<$50k",
        .data[["_INCOMG1"]] %in% 5:7 ~ "$50k+",
        TRUE ~ NA_character_
      ),
      levels = c("<$25k", "$25k-<$50k", "$50k+")
    ),
    disability = factor(
      case_when(
        any_disability ~ "Any disability",
        all_answered_no ~ "No disability",
        TRUE ~ NA_character_
      ),
      levels = c("No disability", "Any disability")
    ),
    state = factor(.data[["_STATE"]]),
    survey_weight = as.numeric(.data[["_LLCPWT"]])
  ) |>
  filter(
    !is.na(frequent_distress),
    !is.na(sex),
    !is.na(race_ethnicity),
    !is.na(age_group),
    !is.na(education),
    !is.na(income),
    !is.na(disability),
    !is.na(state),
    !is.na(survey_weight),
    survey_weight > 0
  ) |>
  mutate(analysis_weight = survey_weight / mean(survey_weight))

write_csv(brfss_reduced, file.path(data_dir, "brfss_reduced_analytic.csv"))
saveRDS(brfss_reduced, file.path(data_dir, "brfss_reduced_analytic.rds"),
        compress = "xz")

stratum_vars <- c(
  "sex", "race_ethnicity", "age_group", "education", "income", "disability"
)

log_msg("Aggregating collapsed weighted strata")
strata_reduced <- brfss_reduced |>
  group_by(across(all_of(stratum_vars))) |>
  summarise(
    raw_n = n(),
    raw_cases = sum(frequent_distress),
    raw_controls = raw_n - raw_cases,
    raw_prevalence = mean(frequent_distress),
    weighted_cases = sum(analysis_weight * frequent_distress),
    weighted_controls = sum(analysis_weight * (1 - frequent_distress)),
    weighted_n = weighted_cases + weighted_controls,
    weighted_prevalence = weighted_cases / weighted_n,
    .groups = "drop"
  ) |>
  arrange(desc(weighted_prevalence))

strata_reduced$label <- do.call(
  paste,
  c(strata_reduced[stratum_vars], sep = " \u00d7 ")
)

write_csv(strata_reduced, file.path(data_dir, "brfss_reduced_strata.csv"))
saveRDS(strata_reduced, file.path(data_dir, "brfss_reduced_strata.rds"),
        compress = "xz")

log_msg("Analytic rows: %s", format(nrow(brfss_reduced), big.mark = ","))
log_msg("Weighted prevalence: %.3f",
        sum(brfss_reduced$frequent_distress * brfss_reduced$analysis_weight) /
          sum(brfss_reduced$analysis_weight))
log_msg("Collapsed strata: %s; n < 20: %s; n < 50: %s",
        nrow(strata_reduced), sum(strata_reduced$raw_n < 20),
        sum(strata_reduced$raw_n < 50))

weighted_formula <- cbind(weighted_cases, weighted_controls) ~
  sex + race_ethnicity + age_group + education + income + disability +
  (1 | sex:race_ethnicity:age_group:education:income:disability)

log_msg("Fitting collapsed weighted grouped-binomial MAIHDA model")
weighted <- suppressWarnings(
  maihda(
    weighted_formula,
    data = strata_reduced,
    family = "binomial",
    interactions = "BH",
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
)

lookup <- strata_reduced |>
  select(label, raw_n, raw_cases, raw_controls, raw_prevalence,
         weighted_n, weighted_cases, weighted_controls, weighted_prevalence)

results <- list(
  analytic_n = nrow(brfss_reduced),
  n_strata = nrow(strata_reduced),
  n_sparse_lt20 = sum(strata_reduced$raw_n < 20),
  n_sparse_lt50 = sum(strata_reduced$raw_n < 50),
  highlight_interactions = highlight_interactions,
  weighted = export_analysis(
    weighted,
    "brfss_reduced_weighted",
    data_dir,
    lookup,
    highlight_interactions = highlight_interactions
  )
)

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0) return("NA")
  x <- as.numeric(x)
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

fmt_p <- function(x) {
  if (length(x) == 0) return("NA")
  x <- as.numeric(x)
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "e", digits = 2)
  out
}

model_stat <- function(model_table, statistic, column) {
  value <- model_table[model_table$statistic == statistic, column, drop = TRUE]
  if (length(value) == 0) NA_real_ else as.numeric(value[1])
}

top_interaction_lines <- function(interactions, n = 5) {
  if (is.null(interactions) || !"flagged" %in% names(interactions)) {
    return("  Strongest flagged interactions: not requested")
  }
  flagged <- interactions |>
    filter(.data$flagged) |>
    arrange(desc(abs(.data$interaction))) |>
    head(n)
  if (nrow(flagged) == 0) {
    return("  Strongest flagged interactions: none")
  }
  c(
    "  Strongest flagged interactions:",
    sprintf(
      "    %s: %s [%s, %s], FDR p=%s (%s)",
      flagged$label,
      fmt_num(flagged$interaction),
      fmt_num(flagged$lower),
      fmt_num(flagged$upper),
      fmt_p(flagged$p_adjusted),
      flagged$direction
    )
  )
}

top_ranked_lines <- function(ranked, n = 5) {
  top <- head(ranked, n)
  c(
    "  Highest predicted-risk strata:",
    sprintf(
      "    %s: predicted=%s, weighted prevalence=%s, raw n=%s",
      top$label,
      fmt_num(top$predicted),
      fmt_num(top$weighted_prevalence),
      format(top$raw_n, big.mark = ",", trim = TRUE)
    )
  )
}

model_summary_lines <- function(label, result) {
  glance <- result$glance
  model_table <- result$model_table
  interactions <- result$interactions
  flagged_n <- if (is.null(interactions)) NA_integer_ else
    sum(interactions$flagged, na.rm = TRUE)
  interaction_n <- if (is.null(interactions)) NA_integer_ else nrow(interactions)

  c(
    sprintf("%s:", label),
    sprintf(
      "  VPC null=%s; adjusted=%s; PCV=%s",
      fmt_num(glance$vpc),
      fmt_num(model_stat(model_table, "VPC/ICC", "adjusted")),
      fmt_num(glance$pcv)
    ),
    sprintf(
      "  AUC null=%s; adjusted=%s; MOR null=%s; adjusted=%s",
      fmt_num(glance$auc),
      fmt_num(glance$auc.adjusted),
      fmt_num(glance$mor),
      fmt_num(model_stat(model_table, "MOR", "adjusted"))
    ),
    sprintf(
      "  BH-flagged interactions: %s of %s",
      flagged_n,
      interaction_n
    ),
    top_interaction_lines(interactions),
    top_ranked_lines(result$ranked_strata)
  )
}

summary_lines <- c(
  "BRFSS 2024 reduced/collapsed MAIHDA case study",
  "",
  sprintf("Analytic rows: %s", format(results$analytic_n, big.mark = ",")),
  sprintf("Collapsed strata: %s", results$n_strata),
  sprintf("Strata with raw n < 20: %s", results$n_sparse_lt20),
  sprintf("Strata with raw n < 50: %s", results$n_sparse_lt50),
  sprintf("Interaction highlighting in plots: %s", results$highlight_interactions),
  "",
  "Weighting: BRFSS final adult weights normalized to mean 1 in the analytic sample.",
  "",
  model_summary_lines(
    "Weighted grouped-binomial model",
    results$weighted
  ),
  "",
  "Key output files:",
  "  brfss_reduced_analytic.csv / .rds",
  "  brfss_reduced_strata.csv / .rds",
  "  brfss_reduced_weighted_* CSV/PNG outputs"
)

saveRDS(results, file.path(data_dir, "brfss_reduced_case_study_results.rds"),
        compress = "xz")
summary_path <- file.path(data_dir, "brfss_reduced_case_study_summary.txt")
summary_con <- file(summary_path, open = "w", encoding = "UTF-8")
writeLines(summary_lines, summary_con, useBytes = FALSE)
close(summary_con)

log_msg("Done. Reduced results written under %s",
        normalizePath(data_dir, winslash = "/"))

# Precompute the unweighted full-data MAIHDA for the vignette (cache + figures).
# This is the slow step (full-data glmer over the intersectional strata, two
# models: minutes), and it is why the vignette ships a cache rather than fitting at
# build time.
log_msg("Precomputing unweighted full-data MAIHDA for the vignette (slow)...")
source(file.path("data-raw", "brfss_precompute.R"))
run_brfss_precompute(data_dir = data_dir)
log_msg("Vignette cache + figures refreshed.")
