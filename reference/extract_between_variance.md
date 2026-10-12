# Extract Between-Stratum Variance

Internal function to extract between-stratum variance from a MAIHDA
model.

## Usage

``` r
extract_between_variance(model, brms_point = c("mean", "median"))
```

## Arguments

- model:

  A maihda_model object

- brms_point:

  For a `brms` fit, the posterior summary of the variance draws to
  return: `"mean"` (default) or `"median"`. Ignored by the other
  engines.

## Value

Numeric value of between-stratum variance
