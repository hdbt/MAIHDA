# Deprecated: use calculate_pcv()

`calculate_pvc()` is the former name of
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md):
the statistic is the PCV (proportional change in variance), but the
historical function name transposed the acronym. `calculate_pvc()` now
forwards to
[`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
with a deprecation warning and will be removed in a future release.

## Usage

``` r
calculate_pvc(
  model1,
  model2,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95
)
```

## Arguments

- model1:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  This is the reference model (typically a simpler or baseline model).

- model2:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  This is the comparison model (typically a more complex model with
  additional predictors).

- bootstrap:

  Logical indicating whether to compute bootstrap confidence intervals
  for the PCV. Default is FALSE. **lme4 engine only**: the parametric
  bootstrap relies on lme4's
  [`simulate()`](https://rdrr.io/r/stats/simulate.html)/`refit()`, so
  for the brms, wemix, and ordinal engines the PCV is reported as a
  point estimate and `bootstrap = TRUE` is an error (see Details).

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default is 0.95.

## Value

See
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).
