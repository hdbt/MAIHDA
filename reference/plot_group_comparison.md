# Plot a MAIHDA Group Comparison (deprecated)

Deprecated. Use [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
on the
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
result instead, e.g. `plot(cmp, type = "vpc")`.

## Usage

``` r
plot_group_comparison(x, type = c("vpc", "components"))
```

## Arguments

- x:

  A `maihda_group_comparison` object.

- type:

  Either "vpc" (default) or "components".

## Value

A ggplot2 object.
