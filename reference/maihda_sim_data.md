# Simulated Health Data for MAIHDA Use

A simulated dataset for demonstrating Multilevel Analysis of Individual
Heterogeneity and Discriminatory Accuracy (MAIHDA).

## Usage

``` r
maihda_sim_data
```

## Format

A data frame with 500 rows and 6 variables:

- id:

  Unique participant identifier.

- gender:

  Gender of the participant.

- race:

  Simulated race/ethnicity category.

- education:

  Educational attainment level.

- age:

  Age in years, a continuous covariate.

- health_outcome:

  A continuous simulated health outcome.

## Source

Simulated for the purpose of the MAIHDA package.

## Examples

``` r
data(maihda_sim_data)
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race", "education"))
```
