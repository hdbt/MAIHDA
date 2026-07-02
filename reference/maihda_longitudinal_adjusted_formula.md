# Build the adjusted-model formula for a longitudinal MAIHDA decomposition

The longitudinal analogue of
[`maihda_adjusted_formula`](https://hdbt.github.io/MAIHDA/reference/maihda_adjusted_formula.md):
the null growth model plus the stratum dimensions' additive main effects
AND their interactions with the time polynomial (`dim:time`), so the
remaining stratum-level intercept/slope variance is the intersectional
interaction beyond additive. Auto-binned numeric dimensions reuse their
reconstructed tertile factor (`.maihda_dim_*`, via
[`maihda_adjusted_terms`](https://hdbt.github.io/MAIHDA/reference/maihda_adjusted_terms.md)).

## Usage

``` r
maihda_longitudinal_adjusted_formula(
  null_formula,
  strata_vars,
  autobin_info,
  data,
  time,
  time_degree
)
```

## Arguments

- null_formula:

  The fitted null growth formula.

- strata_vars, autobin_info, data:

  Stratum metadata (as for `maihda_adjusted_formula`).

- time, time_degree:

  The longitudinal specification; `time` must be the variable the growth
  terms were built on (the internally centered column for a centered
  fit, `maihda_lng_time_term()`), so the `dim:time` interactions
  reference the same terms as the null formula.

## Value

A list with `formula` and `data`, or `NULL` if fewer than two dimensions
are available.
