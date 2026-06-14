# Build the 3-level growth formula for a longitudinal MAIHDA

Given a base formula already carrying the covariates and the resolved
stratum grouping (`y ~ covars + (1 | stratum)`), returns the growth
formula
`y ~ covars + time(+ I(time^2)...) + (time... | id) + (time... | stratum)`:
the time polynomial enters the fixed part (if absent) and a random
intercept+slope block is placed at both the individual and stratum
levels. Any random effects in the base formula are replaced by this
canonical structure.

## Usage

``` r
maihda_longitudinal_formula(base_formula, id, time, time_degree)
```

## Arguments

- base_formula:

  The resolved formula (fixed part + stratum grouping).

- id, time, time_degree:

  The longitudinal specification.

## Value

The growth formula (same environment as `base_formula`).
