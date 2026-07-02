# Per-stratum trajectory parameters for a longitudinal MAIHDA

The stratum-level random-effect estimates as a wide table, one row per
stratum: the stratum's deviation at the baseline time (`baseline`, the
longitudinal analogue of a cross-sectional stratum BLUP), the random
intercept (`intercept`; the deviation at the model's coefficient origin
– the internal centering offset, i.e. the observed baseline, when
centering applied, or raw time 0 otherwise) and the random slope(s) on
time (`slope`, ...). This is the longitudinal shape of
`predict_maihda(type = "strata")` – a stratum is now a *trajectory*, not
a single value.

## Usage

``` r
maihda_longitudinal_strata_predictions(object)
```

## Arguments

- object:

  A longitudinal `maihda_model`.

## Value

A data frame: `stratum`, `stratum_id`, optional `label`, `baseline`,
`intercept`, `slope`(, `slope2`, ...).

## Details

`baseline` is \\a(t_0 - c)' coef\\ with \\a(t) = (1, t, t^2, ...)\\,
\\t_0 = \\ the reference (baseline) time `ref_time = min(time)` and
\\c\\ the internal centering offset; it equals `intercept` whenever
`ref_time` coincides with the centering origin (the usual case for a
centered fit, and the zero-anchored case for a raw fit).
