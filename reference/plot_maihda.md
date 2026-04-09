# Plot MAIHDA Model Results

Creates various plots for visualizing MAIHDA model results including
variance partition coefficient comparisons, observed vs. shrunken
estimates, and predicted subgroup values with confidence intervals.

## Usage

``` r
plot_maihda(
  object,
  type = c("vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect", "effect_decomp"),
  summary_obj = NULL,
  n_strata = 50,
  ...
)
```

## Arguments

- object:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- type:

  Character string specifying plot type:

  - "vpc": Variance partition coefficient visualization

  - "obs_vs_shrunken": Observed vs. shrunken stratum means

  - "predicted": Predicted values for each stratum with confidence
    intervals

- summary_obj:

  Optional maihda_summary object from
  [`summary_maihda()`](https://hdbt.github.io/MAIHDA/reference/summary_maihda.md).
  If NULL, will be computed.

- n_strata:

  Maximum number of strata to display in predicted plot. Default is 50.
  Use NULL for all strata.

- ...:

  Additional arguments (not currently used).

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)

# VPC plot
plot_maihda(model, type = "vpc")


# Observed vs shrunken plot
plot_maihda(model, type = "obs_vs_shrunken")
#> function (x, y, ...) 
#> UseMethod("plot")
#> <bytecode: 0x55d1418e7c60>
#> <environment: namespace:base>

# Predicted values with confidence intervals
plot_maihda(model, type = "predicted")

# }
```
