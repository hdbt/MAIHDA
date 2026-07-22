# Guard against unidentified growth terms (distinct-times checks)

A degree-`d` growth curve is only estimable when the observed times
actually span it. Two conditions, both required: the polynomial basis
needs at least `d + 1` distinct observed time values (on a single time
even the linear slope does not exist; on two times a quadratic is
collinear with the linear term), and the person-level slope block
`(time | id)` needs at least one id measured at two distinct times (ids
repeating only at a single time are replicate measurements, not a growth
design). Data failing either check can still FIT: lme4 drops
rank-deficient fixed-effect columns but keeps the random-effect columns,
leaving the likelihood flat (or near-flat) in the affected covariance
parameters, so the "slope variances" it returns are arbitrary –
observed: optimizer starting values reported back verbatim on all-equal
times, and a ~1300-fold optimizer-dependent swing in the between-stratum
variance at unobserved times for a quadratic on two waves – with at most
a "boundary (singular) fit" note. Reject with guidance instead.

## Usage

``` r
maihda_check_longitudinal_times(data, id, time, time_degree)
```

## Arguments

- data:

  The model data (full input or analytic sample).

- id, time:

  Column names of the person identifier and time variable.

- time_degree:

  Integer polynomial degree of the growth curve.

## Value

`NULL`, invisibly; stops when the growth terms are unidentified.

## Details

Called on the full input by
[`maihda_validate_longitudinal()`](https://hdbt.github.io/MAIHDA/reference/maihda_validate_longitudinal.md)
and again on the analytic sample by
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md):
dropped rows (subset, missing outcome/covariate/weight) can only REDUCE
the distinct times, so the early check never rejects a fittable model,
and the re-check catches data whose identifying occasions sit on rows
the fit drops.
