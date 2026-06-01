# Plot MAIHDA Ternary Diagram

Renders the ternary decomposition produced by
[`compute_maihda_ternary_data`](https://hdbt.github.io/MAIHDA/reference/compute_maihda_ternary_data.md).
Dispatched via [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
on the classed result.

## Usage

``` r
# S3 method for class 'maihda_ternary'
plot(
  x,
  size_var = "n",
  color_var = "label",
  label_top_n = 5,
  label_by = c("interaction_signal", "uncertainty", "n"),
  alpha = 0.7,
  ...
)
```

## Arguments

- x:

  A `maihda_ternary` object from `compute_maihda_ternary_data`.

- size_var:

  Column name for point sizing.

- color_var:

  Column name for point colors.

- label_top_n:

  Number of top strata to label.

- label_by:

  Variable used to determine top strata.

- alpha:

  Point transparency.

- ...:

  Additional arguments (not used).

## Value

A plot object.
