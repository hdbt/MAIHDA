# Extract Between-Stratum Variance

Internal function to extract between-stratum variance from a MAIHDA
model. For a `brms` fit it is the posterior median of the variance
draws.

## Usage

``` r
extract_between_variance(model)
```

## Arguments

- model:

  A maihda_model object

## Value

Numeric value of between-stratum variance
