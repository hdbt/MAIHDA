# Fit a design-weighted MAIHDA model via WeMix

Internal engine call for `fit_maihda(engine = "wemix")`. Builds the
analytic sample (complete cases on the model variables and the weight
column, positive weights only) so the stored `data` matches the rows
WeMix fits, attaches the constant level-2 weight column (strata are
exhaustive population cells, sampled with certainty), and calls
[`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
with the unconditional weights `c(level1, level2)`.

## Usage

``` r
maihda_fit_wemix(formula, data, family, sampling_weights, dot_vals)
```

## Arguments

- formula:

  The resolved model formula (with `(1 | stratum)`).

- data:

  The data (after strata creation / response recoding).

- family:

  The resolved family object (gaussian-identity or binomial-logit).

- sampling_weights:

  Name of the level-1 sampling-weight column.

- dot_vals:

  Named list of evaluated `...` arguments forwarded to
  [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
  (e.g. `nQuad`, `verbose`, `fast`).

## Value

A list with `model` (the `WeMixResults`) and `data` (the analytic data
frame actually fitted, including the weight columns).
