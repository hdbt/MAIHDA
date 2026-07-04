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

## Details

As in
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
REML `lmer` growth fits are refitted with maximum likelihood
([`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html)) before the
comparison: the null and adjusted models differ in fixed effects (the
dimensions' main effects and their `dim:time` interactions), across
which REML variance estimates are not comparable – using them biases
both PCVs downward, overstating the multiplicative/interaction share.
The stored models (and the single-model summaries computed from them,
e.g. the time-varying VPC) keep their REML fit; `ml_refit` on the result
records whether the refit applied. glmer (GLMM) and brms fits are
already on the ML / posterior scale.
