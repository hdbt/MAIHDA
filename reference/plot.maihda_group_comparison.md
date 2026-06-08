# Plot a MAIHDA Group Comparison

Visualises the output of
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
as a point/forest plot of the VPC/ICC by group, as stacked
variance-composition bars (between- vs within-stratum share) by group,
or as bars of the absolute between-stratum (intersectional) variance by
group. Dispatched via
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the classed
result.

## Usage

``` r
# S3 method for class 'maihda_group_comparison'
plot(x, type = c("vpc", "components", "between_variance"), ...)
```

## Arguments

- x:

  A `maihda_group_comparison` object from
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).

- type:

  One of "vpc" (default) for VPC by group with optional bootstrap
  confidence intervals, "components" for stacked variance proportions,
  or "between_variance" for the absolute between-stratum variance by
  group. The VPC is a *share* of the unexplained variance;
  "between_variance" shows the *magnitude* the ratio cannot convey (two
  groups with very different VPCs can share a between-stratum variance,
  and vice versa).

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

plot(cmp, type = "between_variance")

# }
```
