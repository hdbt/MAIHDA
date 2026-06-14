# Longitudinal MAIHDA proportional change in variance (PCV)

Compares the stratum-level random-effect covariance block of the null
growth model with that of the adjusted model (null + dimension main
effects + their `dim:time` interactions). Reports the PCV in the
baseline (intercept) variance and in the slope variance – the
additive-vs-multiplicative split of the intersectional trajectory
inequality (Bell, Evans, Holman & Leckie 2024) – and the time-specific
PCV over the supplied times.

## Usage

``` r
maihda_longitudinal_pcv(null_model, adjusted_model, times = NULL)
```

## Arguments

- null_model, adjusted_model:

  Longitudinal `maihda_model`s from a
  `maihda(decomposition = "longitudinal")` pair.

- times:

  Optional numeric times for the time-specific PCV; defaults to the null
  model's reporting grid.

## Value

An object of class `maihda_long_pcv`.
