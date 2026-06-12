# Plot MAIHDA Model Results

Creates various plots for visualizing MAIHDA model results including
variance partition coefficient comparisons, observed vs. shrunken
estimates, and predicted subgroup values with confidence intervals.

## Usage

``` r
# S3 method for class 'maihda_model'
plot(
  x,
  type = c("all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect",
    "effect_decomp", "ternary", "prediction_deviation", "context_vpc"),
  summary_obj = NULL,
  n_strata = 50,
  ...
)
```

## Arguments

- x:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- type:

  Character string specifying plot type:

  - "vpc": Variance partition coefficient visualization

  - "obs_vs_shrunken": Observed vs. shrunken stratum means. The y-axis
    (model-based estimate) includes the fixed effects, so for a
    covariate-adjusted model the distance from the diagonal reflects
    both shrinkage *and* covariate adjustment, not shrinkage alone; it
    is a pure shrinkage view only for an intercept-only (null) model

  - "predicted": Predicted values for each stratum with confidence
    intervals

  - "risk_vs_effect": Quadrant scatterplot of each stratum's mean
    predicted outcome against its random effect

  - "effect_decomp": Visualizes additive vs intersectional deviation
    from global mean

  - "ternary": Ternary diagnostic of the relative additive,
    intersectional, and uncertainty signals per stratum (a
    normalized-magnitude diagnostic, not a variance decomposition)

  - "prediction_deviation": Detailed deviation panels for individuals or
    strata

  - "context_vpc": Stratum vs. context variance bars for a contextual
    cross-classified fit (`fit_maihda(context = )`); errors otherwise

  - "all": Generate all available plots (default if not specified)

- summary_obj:

  Optional maihda_summary object from
  [`summary()`](https://rdrr.io/r/base/summary.html). If NULL, will be
  computed.

- n_strata:

  Maximum number of strata to display in the predicted plot. When there
  are more strata than this, the first `n_strata` (in stratum order) are
  shown and the plot caption notes how many were omitted. Default is 50.
  Use NULL for all strata.

- ...:

  Additional arguments (not currently used).

## Value

A ggplot2 object, or a list of ggplot2 objects if type = "all".

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)

# VPC plot
plot(model, type = "vpc")


# Generate all plots
plots <- plot(model)





#> Warning: Removing Layer 2 ('PositionNudge'), as it is not an approved position (for ternary plots) under the present ggtern package.


# }
```
