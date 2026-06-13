# Stratum random-effect table for a cumulative (clmm) fit

Mirrors `maihda_stratum_ranef_lme4()`: one row per stratum with the
conditional mode, its conditional standard error (from
[`ordinal::condVar()`](https://rdrr.io/pkg/ordinal/man/ranef.html),
which returns conditional *variances*), and a 95% interval. At a
boundary fit (zero between-stratum variance) the conditional
distribution collapses on 0, so the SE is 0.

## Usage

``` r
maihda_clmm_stratum_ranef(object)
```

## Arguments

- object:

  A `maihda_model` with engine `"ordinal"`.

## Value

A data frame with `stratum`, `stratum_id`, `random_effect`, `se`,
`lower_95`, `upper_95`.
