# Time-varying VPC summary for a longitudinal MAIHDA (lme4)

Time-varying VPC summary for a longitudinal MAIHDA (lme4)

## Usage

``` r
maihda_longitudinal_summary_lme4(
  object,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95
)
```

## Arguments

- object:

  A longitudinal `maihda_model` (lme4 engine).

- bootstrap, n_boot, conf_level:

  Parametric-bootstrap controls for the VPC(t) band.

## Value

A list with `vpc_result` (the reference-time VPC, for the headline
print), `variance_components`, and `longitudinal` (the trajectory).
