# Plot a PCV importance attribution

Bar chart of the per-variable contributions to the total PCV, ordered by
size, with a dashed zero line (negative bars flag suppression) and, when
available, parametric-bootstrap confidence intervals. The bars sum to
the full-model Total PCV (the efficiency property).

## Usage

``` r
# S3 method for class 'maihda_pcv_importance'
plot(x, ...)
```

## Arguments

- x:

  A `maihda_pcv_importance` object from
  [`pcv_importance`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md).

- ...:

  Additional arguments (not used).

## Value

A ggplot object.
