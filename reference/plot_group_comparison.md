# Plot a MAIHDA Group Comparison

Visualises the output of
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
either as a point/forest plot of the VPC/ICC by group, or as stacked
variance-composition bars (between- vs within-stratum share) by group.

## Usage

``` r
plot_group_comparison(x, type = c("vpc", "components"))
```

## Arguments

- x:

  A `maihda_group_comparison` object from
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).

- type:

  Either "vpc" (default) for VPC by group with optional bootstrap
  confidence intervals, or "components" for stacked variance
  proportions.

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
data(maihda_health_data)
cmp <- compare_maihda_groups(BMI ~ Age + (1 | Gender:Race),
                             data = maihda_health_data, group = "Education")
plot_group_comparison(cmp, type = "vpc")

plot_group_comparison(cmp, type = "components")

# }
```
