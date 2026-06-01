# Plot a MAIHDA Group Comparison

Visualises the output of
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
either as a point/forest plot of the VPC/ICC by group, or as stacked
variance-composition bars (between- vs within-stratum share) by group.
Dispatched via [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
on the classed result.

## Usage

``` r
# S3 method for class 'maihda_group_comparison'
plot(x, type = c("vpc", "components"), ...)
```

## Arguments

- x:

  A `maihda_group_comparison` object from
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).

- type:

  Either "vpc" (default) for VPC by group with optional bootstrap
  confidence intervals, or "components" for stacked variance
  proportions.

- ...:

  Additional arguments (not used).

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
data(maihda_health_data)
cmp <- compare_maihda_groups(BMI ~ Age + (1 | Gender:Race),
                             data = maihda_health_data, group = "Education")
plot(cmp, type = "vpc")

plot(cmp, type = "components")

# }
```
