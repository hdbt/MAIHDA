# data-raw/fetch_ess_youth.R
# Builds the youth analysis frames for the ESS case-study vignette from the
# European Social Survey integrated files (Rounds 9-11). The ESS microdata is
# license-restricted and NOT shipped with the package: download the integrated
# .dta files yourself into data-raw/ess/ (see data-raw/ess/README.md), then run
# this script to regenerate the analysis frames. Only derived, aggregate outputs
# (the precomputed summaries + figures) are distributed; see
# data-raw/precompute_youth_vignette.R.
#
# Source: European Social Survey ESS9-ESS11, ESS Data Portal (Sikt).

if (!requireNamespace("haven", quietly = TRUE)) {
  stop("Install 'haven' to read the ESS .dta files: install.packages('haven')")
}
suppressMessages({library(haven); library(dplyr)})

# ---- knobs ----------------------------------------------------------------
ESS_DIR  <- "data-raw/ess"
ROUND_FILES <- c(
  file.path(ESS_DIR, "ESS9e03_3",  "ESS9e03_3.dta"),
  file.path(ESS_DIR, "ESS10e03_3", "ESS10e03_3.dta"),
  file.path(ESS_DIR, "ESS11e04_1", "ESS11e04_1.dta")
)
AGE_MIN <- 15L
AGE_MAX <- 29L
# ---------------------------------------------------------------------------

need <- c("essround", "idno", "cntry", "gndr", "agea", "mnactic", "eisced",
          "brncntr", "facntr", "mocntr", "anweight", "pspwght")

read_round <- function(f) {
  if (!file.exists(f)) stop("Missing ESS file: ", f, ".\nSee data-raw/ess/README.md.")
  haven::read_dta(f, col_select = dplyr::any_of(need)) |> haven::zap_labels()
}

# Base youth frame: the three social dimensions, country, and the raw main-activity
# code. `mnactic`: 1 paid work, 2 education, 3 unemployed-looking, 4 unemployed-not-
# looking, 5 sick/disabled, 6 retired, 7 community/military service, 8 housework/
# care, 9 other.
base <- dplyr::bind_rows(lapply(ROUND_FILES, read_round)) |>
  filter(!is.na(agea), agea >= AGE_MIN, agea <= AGE_MAX) |>
  transmute(
    essround = as.integer(essround),
    country  = as.character(cntry),
    age      = as.integer(agea),
    mnactic  = as.integer(mnactic),
    gender = factor(gndr, levels = c(1, 2), labels = c("male", "female")),
    migration = factor(case_when(
      brncntr == 2                               ~ "first_gen",
      brncntr == 1 & (facntr == 2 | mocntr == 2) ~ "second_gen",
      brncntr == 1 & facntr == 1 & mocntr == 1   ~ "native",
      TRUE                                       ~ NA_character_
    ), levels = c("native", "second_gen", "first_gen")),
    education = factor(case_when(
      eisced %in% c(1, 2)    ~ "low",
      eisced %in% c(3, 4, 5) ~ "mid",
      eisced %in% c(6, 7)    ~ "high",
      TRUE                   ~ NA_character_
    ), levels = c("high", "mid", "low")),
    anweight = anweight
  ) |>
  filter(!is.na(gender), !is.na(migration), !is.na(education))
# NB: no 'stratum' column is stored -- maihda()'s (1 | gender:migration:education)
# shorthand builds the strata itself and records the defining dimensions.

# Frame A -- UNEMPLOYMENT among the economically active (ILO labour-force
# denominator): employed (mnactic 1) vs unemployed-and-seeking (mnactic 3).
ess_youth <- base |>
  filter(mnactic %in% c(1, 3)) |>
  mutate(unemployed = as.integer(mnactic == 3)) |>
  select(essround, country, age, unemployed, gender, migration, education, anweight)

# Frame B -- NEET (not in employment, education or training) on ALL youth; pulls in
# the inactive/home-care channel. NEET numerator = {3,4,5,8,9}; denominator excludes
# the ambiguous retired (6) and community/military service (7).
ess_youth_neet <- base |>
  filter(mnactic %in% c(1, 2, 3, 4, 5, 8, 9)) |>
  mutate(neet = as.integer(mnactic %in% c(3, 4, 5, 8, 9))) |>
  select(essround, country, age, neet, gender, migration, education, anweight)

# ---- diagnostics ----------------------------------------------------------
nstrata <- function(d) nlevels(interaction(d$gender, d$migration, d$education, drop = TRUE))
cat("ESS youth frames (rounds", paste(sort(unique(base$essround)), collapse = "/"), ")\n")
cat(sprintf("  UNEMPLOYMENT: N=%d, %d countries, %d strata (min cell %d), rate %.3f\n",
            nrow(ess_youth), dplyr::n_distinct(ess_youth$country), nstrata(ess_youth),
            min(table(interaction(ess_youth$gender, ess_youth$migration, ess_youth$education, drop = TRUE))),
            mean(ess_youth$unemployed)))
cat(sprintf("  NEET        : N=%d, %d countries, %d strata (min cell %d), rate %.3f\n",
            nrow(ess_youth_neet), dplyr::n_distinct(ess_youth_neet$country), nstrata(ess_youth_neet),
            min(table(interaction(ess_youth_neet$gender, ess_youth_neet$migration, ess_youth_neet$education, drop = TRUE))),
            mean(ess_youth_neet$neet)))

saveRDS(ess_youth,      file.path(ESS_DIR, "ess_youth_analysis.rds"))
saveRDS(ess_youth_neet, file.path(ESS_DIR, "ess_youth_neet.rds"))
cat("  saved ess_youth_analysis.rds + ess_youth_neet.rds\n")
