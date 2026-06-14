# Time-specific PCV plot (longitudinal MAIHDA)

The additive share – the proportional change in the between-stratum
(trajectory) variance from the null to the adjusted model – as a
function of time. A high, flat curve means intersectional trajectory
inequalities are "mostly additive".

## Usage

``` r
plot_pcv_trajectory(pcv)
```

## Arguments

- pcv:

  A `maihda_long_pcv` from a longitudinal
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) pair.

## Value

A ggplot2 object.
