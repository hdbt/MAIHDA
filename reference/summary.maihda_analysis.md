# Summarize a MAIHDA Analysis

Returns the variance summary (VPC/ICC, variance components, stratum
estimates) of the fitted model. The per-group comparison, when present,
is attached as the `"groups"` attribute.

## Usage

``` r
# S3 method for class 'maihda_analysis'
summary(object, ...)
```

## Arguments

- object:

  A `maihda_analysis` object from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

- ...:

  Additional arguments (not used).

## Value

The `maihda_summary` for the fitted model.
