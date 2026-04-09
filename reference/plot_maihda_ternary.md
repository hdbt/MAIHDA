# Plot MAIHDA Ternary Diagram

Plot MAIHDA Ternary Diagram

## Usage

``` r
plot_maihda_ternary(
  ternary_data,
  size_var = "n",
  color_var = "label",
  label_top_n = 5,
  label_by = c("interaction_signal", "uncertainty", "n"),
  alpha = 0.7
)
```

## Arguments

- ternary_data:

  Data output from `compute_maihda_ternary_data`.

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

## Value

A plot object.
