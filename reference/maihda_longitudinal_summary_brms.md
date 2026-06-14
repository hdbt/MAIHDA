# Time-varying VPC summary for a longitudinal MAIHDA (brms, linear growth)

Time-varying VPC summary for a longitudinal MAIHDA (brms, linear growth)

## Usage

``` r
maihda_longitudinal_summary_brms(object, conf_level = 0.95)
```

## Arguments

- object:

  A longitudinal `maihda_model` (brms engine, time_degree 1).

- conf_level:

  Credible-interval level.

## Value

As `maihda_longitudinal_summary_lme4`, with posterior bands.
