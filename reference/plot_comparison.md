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
if (FALSE) { # \dontrun{
comparison <- compare_maihda(model1, model2, bootstrap = TRUE)
plot_comparison(comparison)
} # }
```
