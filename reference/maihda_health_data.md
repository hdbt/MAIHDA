# NHANES Health Data Subset for MAIHDA Use

A pedagogical subset of the National Health and Nutrition Examination
Survey (NHANES) dataset, serving as a real-world example for Multilevel
Analysis of Individual Heterogeneity and Discriminatory Accuracy
(MAIHDA). Contains selected records demonstrating intersectional
demographic health inequalities.

## Usage

``` r
maihda_health_data
```

## Format

A data frame with 3,000 rows and 7 variables:

- BMI:

  Body Mass Index (kg/m^2), a continuous outcome variable.

- Obese:

  Factor indicating obesity status (No/Yes).

- Age:

  Age in years at screening, a continuous covariate.

- Gender:

  Gender of the participant (male/female).

- Race:

  Self-reported race/ethnicity.

- Education:

  Educational attainment level.

- Poverty:

  Poverty to income ratio, a continuous covariate. Some values may be
  missing.

## Source

Derived from the `NHANES` R package. Original data collected by the
Centers for Disease Control and Prevention (CDC).

## Examples

``` r
data(maihda_health_data)

# Example usage:
# strata_result <- make_strata(maihda_health_data, vars = c("Gender", "Race", "Education"))
# model <- fit_maihda(BMI ~ Age + (1 | stratum), data = strata_result$data)
```
