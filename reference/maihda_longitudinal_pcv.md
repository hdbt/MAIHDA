# Longitudinal MAIHDA proportional change in variance (PCV)

Compares the stratum-level random-effect covariance block of the null
growth model with that of the adjusted model (null + dimension main
effects + their `dim:time` interactions). Reports the PCV in the
baseline variance and in the instantaneous-slope variance at baseline –
the additive-vs-multiplicative split of the intersectional trajectory
inequality (Bell, Evans, Holman & Leckie 2024) – and the time-specific
PCV over the supplied times. Both are evaluated at the observed baseline
time (`ref_time`), so they are invariant to how the time axis is coded
(for linear growth the slope variance is the same at every time, so this
reduces to the slope-variance cell).

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
