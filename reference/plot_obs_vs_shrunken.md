# Observed vs. Shrunken Estimates Plot

Observed vs. Shrunken Estimates Plot

## Usage

``` r
plot_obs_vs_shrunken(object, summary_obj)
```

## Arguments

- object:

  A maihda_model object

- summary_obj:

  A maihda_summary object

## Value

A ggplot2 object

## Details

The x-axis is each stratum's raw observed mean; the y-axis is the
model-based stratum estimate, which includes the fixed-effect
contribution. For an intercept-only (null) model the vertical distance
from the diagonal is pure shrinkage toward the grand mean. For a
covariate-adjusted model the model estimate also moves with the
stratum's covariate profile, so distance from the diagonal reflects
*both* shrinkage and covariate adjustment and should not be read as
shrinkage alone. The caption notes which case applies.
