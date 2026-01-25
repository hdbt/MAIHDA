# Plot MAIHDA Model Results

Creates various plots for visualizing MAIHDA model results including
caterpillar plots, variance partition coefficient comparisons, observed
vs. shrunken estimates, and predicted subgroup values with confidence
intervals.

## Usage

``` r
plot_maihda(object, type = c("caterpillar", "vpc", "obs_vs_shrunken", "predicted"),
            summary_obj = NULL, n_strata = 50, ...)
```

## Arguments

- object:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- type:

  Character string specifying plot type: "caterpillar" for caterpillar
  plot of stratum random effects, "vpc" for variance partition
  coefficient visualization, "obs_vs_shrunken" for observed vs. shrunken
  stratum means, or "predicted" for predicted values for each stratum
  with confidence intervals.

- summary_obj:

  Optional maihda_summary object from
  [`summary_maihda()`](https://hdbt.github.io/MAIHDA/reference/summary_maihda.md).
  If NULL, will be computed.

- n_strata:

  Maximum number of strata to display in caterpillar plot or predicted
  plot. Default is 50. Use NULL for all strata.

- ...:

  Additional arguments (not currently used).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)

# Caterpillar plot
plot_maihda(model, type = "caterpillar")

# VPC plot
plot_maihda(model, type = "vpc")

# Observed vs shrunken plot
plot_maihda(model, type = "obs_vs_shrunken")

# Predicted values with confidence intervals
plot_maihda(model, type = "predicted")
} # }
```
