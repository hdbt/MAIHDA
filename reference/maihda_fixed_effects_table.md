# Fixed-effect table with Wald statistics

Internal helper that assembles the `fixed_effects` slot of a
`maihda_summary` for the frequentist engines from point estimates and
standard errors: the Wald statistic (`estimate / se`), its two-sided
normal-approximation p-value, and the matching Wald interval at
`conf_level`. The normal approximation is the same one
[`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
uses for the stratum BLUPs – no denominator degrees of freedom are
estimated, so for a Gaussian fit these are z-based, not
Satterthwaite/Kenward-Roger. That is essentially exact for a covariate
varying within stratum, but anticonservative for terms constant within a
stratum (the intercept and an adjusted model's dimension main effects),
whose effective sample size is the number of strata rather than \\n\\;
see the “Fixed-effect statistics” section of
[`maihda_tidiers`](https://hdbt.github.io/MAIHDA/reference/maihda_tidiers.md).

## Usage

``` r
maihda_fixed_effects_table(term, estimate, se, conf_level = 0.95)
```

## Arguments

- term:

  Character vector of coefficient names.

- estimate:

  Numeric vector of point estimates.

- se:

  Numeric vector of standard errors (`NULL` for none).

- conf_level:

  Interval level.

## Value

A data frame with `term`, `estimate`, `se`, `statistic`, `p_value`,
`lower` and `upper`.
