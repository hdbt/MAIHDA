# data-raw/maihda_health_data.R
# Prepares a subset of the National Health and Nutrition Examination Survey (NHANES)
# for pedagogical use in the MAIHDA R package.

if (!requireNamespace("NHANES", quietly = TRUE)) {
  install.packages("NHANES")
}

library(dplyr)

set.seed(42)

# Select a meaningful subset of the NHANES data
# - BMI: Continuous outcome variable
# - Obese: Binary outcome variable (BMI >= 30)
# - Age: Continuous covariate
# - Gender: Dichotomous strata grouping
# - Race1: Polytomous strata grouping
# - Education: Polytomous strata grouping
# - Poverty: Income-to-poverty ratio (covariate)
maihda_health_data <- NHANES::NHANES %>%
  select(BMI, Age, Gender, Race1, Education, Poverty) %>%
  filter(!is.na(BMI), !is.na(Gender), !is.na(Race1), !is.na(Education)) %>%
  # Sample 3000 rows to keep the package lightweight but statistically meaningful
  sample_n(3000) %>%
  mutate(
    Obese = factor(ifelse(BMI >= 30, "Yes", "No")),
    Gender = as.factor(Gender),
    Race = as.factor(Race1),
    Education = as.factor(Education)
  ) %>%
  select(BMI, Obese, Age, Gender, Race, Education, Poverty)

usethis::use_data(maihda_health_data, overwrite = TRUE)
