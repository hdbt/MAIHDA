# Plot Model Comparison

Creates a plot comparing VPC/ICC across multiple models.

## Usage

``` r
plot_comparison(comparison_df)
```

## Arguments

- comparison_df:

  A data frame from
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md).

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
# Create strata and models using simulated data
strata_1 <- make_strata(maihda_sim_data, vars = c("gender", "race"))
strata_2 <- make_strata(maihda_sim_data, vars = c("gender", "race", "education"))

model1 <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_1$data)
model2 <- fit_maihda(health_outcome ~ age + gender + (1 | stratum), data = strata_2$data)

comparison <- compare_maihda(model1, model2, bootstrap = TRUE)
plot_comparison(comparison)

# }
```
