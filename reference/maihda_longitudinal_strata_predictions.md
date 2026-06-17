# Per-stratum trajectory parameters for a longitudinal MAIHDA

The stratum-level random-effect estimates as a wide table, one row per
stratum: the stratum's deviation at the baseline time (`baseline`, the
longitudinal analogue of a cross-sectional stratum BLUP), the raw random
intercept at time 0 (`intercept`) and the random slope(s) on time
(`slope`, ...). This is the longitudinal shape of
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

`baseline` is \\a(t_0)' coef\\ with \\a(t) = (1, t, t^2, ...)\\ and
\\t_0 = \\ the reference (baseline) time `ref_time = min(time)`; the
package defines the baseline at `ref_time`, so it equals the raw
`intercept` (deviation at time 0) only when time is zero-referenced.
