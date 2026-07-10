# Plot a MAIHDA Model Comparison

Plots VPC/ICC across the models compared by
[`compare_maihda`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md).
Dispatched via [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
on the classed result.

## Usage

``` r
# S3 method for class 'maihda_comparison'
plot(x, ...)
```

## Arguments

- x:

  A `maihda_comparison` object from
  [`compare_maihda`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md).

- ...:

  Additional arguments (not used).

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))

null_model <- fit_maihda(health_outcome ~ 1 + (1 | stratum), data = strata$data)
adj_model  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata$data)

comparison <- compare_maihda(null_model, adj_model, bootstrap = TRUE)
#> boundary (singular) fit: see help('isSingular')
plot(comparison)

# }
```
