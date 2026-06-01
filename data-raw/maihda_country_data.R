# data-raw/maihda_country_data.R
# Builds `maihda_country_data`: a cross-national educational-achievement dataset
# used to demonstrate comparing intersectional inequality (VPC/ICC) ACROSS a
# higher-level grouping variable (country) via compare_maihda_groups() / maihda().
#
# Source: OECD Programme for International Student Assessment (PISA) 2018, accessed
# through the 'learningtower' package (MIT licensed), which provides a cleaned,
# reproducible extract of the public OECD PISA microdata. 'learningtower' is only
# required to *regenerate* this dataset (like 'NHANES' for maihda_health_data); it
# is not a dependency of the MAIHDA package itself.
#
# Strata variables (intersectional): gender x ses (socioeconomic tertile).
# Grouping variable: country. Outcome: PISA mathematics score.

if (!requireNamespace("learningtower", quietly = TRUE)) {
  stop("Install 'learningtower' to regenerate this dataset: install.packages('learningtower')")
}

library(dplyr)

set.seed(2026)

# Six countries chosen a priori for geographic/economic diversity (NOT selected by
# their VPC). Stored with readable names.
country_map <- c(
  FIN = "Finland",
  DEU = "Germany",
  GBR = "United Kingdom",
  ITA = "Italy",
  JPN = "Japan",
  MEX = "Mexico"
)

n_per_country <- 600

raw <- learningtower::load_student("2018")

clean <- raw %>%
  filter(as.character(.data$country) %in% names(country_map)) %>%
  filter(!is.na(.data$gender), !is.na(.data$escs),
         !is.na(.data$math), !is.na(.data$read)) %>%
  mutate(
    country = factor(unname(country_map[as.character(.data$country)]),
                     levels = unname(country_map)),
    gender = factor(as.character(.data$gender), levels = c("female", "male")),
    # drop SPSS label/format attributes carried on the score columns
    math = as.numeric(.data$math),
    reading = as.numeric(.data$read),
    escs = as.numeric(.data$escs)
  )

# Balanced subsample per country for a compact, evenly-weighted teaching dataset.
sampled <- clean %>%
  group_by(.data$country) %>%
  slice_sample(n = n_per_country) %>%
  ungroup()

# Global socioeconomic tertiles on the pooled sample, so a stratum denotes the
# same SES band in every country (matches shared_strata = TRUE).
ses_breaks <- quantile(sampled$escs, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
ses_breaks[1] <- -Inf
ses_breaks[length(ses_breaks)] <- Inf

maihda_country_data <- sampled %>%
  mutate(
    ses = cut(.data$escs, breaks = ses_breaks,
              labels = c("Low", "Medium", "High"), include.lowest = TRUE),
    # PISA proficiency Level 2 baseline in mathematics is ~420 points.
    low_math = factor(ifelse(.data$math < 420, "Yes", "No"), levels = c("No", "Yes"))
  ) %>%
  transmute(
    country = .data$country,
    gender = .data$gender,
    ses = .data$ses,
    escs = round(.data$escs, 3),
    math = round(.data$math, 1),
    reading = round(.data$reading, 1),
    low_math = .data$low_math
  ) %>%
  as.data.frame()

stopifnot(
  nrow(maihda_country_data) == n_per_country * length(country_map),
  !anyNA(maihda_country_data[, c("country", "gender", "ses", "math")])
)

usethis::use_data(maihda_country_data, overwrite = TRUE, compress = "xz")
