# Per-stratum trajectory parameters for a longitudinal MAIHDA

The stratum-level random-effect estimates as a wide table: the random
intercept (baseline deviation) and the random slope(s) on time
(trajectory deviation), one row per stratum. This is the longitudinal
shape of `predict_maihda(type = "strata")` – a stratum is now a
*trajectory*, not a single value.

## Usage

``` r
maihda_longitudinal_strata_predictions(object)
```

## Arguments

- object:

  A longitudinal `maihda_model`.

## Value

A data frame: `stratum`, `stratum_id`, optional `label`, `intercept`,
`slope`(, `slope2`, ...).
