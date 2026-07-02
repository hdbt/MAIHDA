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
maihda_longitudinal_formula(
  base_formula,
  id,
  time,
  time_degree,
  orig_time = time
)
```

## Arguments

- base_formula:

  The resolved formula (fixed part + stratum grouping).

- id, time, time_degree:

  The longitudinal specification; `time` is the model variable the
  growth terms are built on (the centered column when centering
  applies).

- orig_time:

  The user's original time column name; differs from `time` only when
  centering applies.

## Value

The growth formula (same environment as `base_formula`).

## Details

When the growth terms are built on internally centered time (`time` is
the derived `.maihda_ctime` column; see the file header), any bare
raw-time polynomial the user wrote in the fixed part (`orig_time`,
`I(orig_time^2)`, ...) is *replaced* by the centered terms rather than
kept alongside them, which would be perfectly collinear.
