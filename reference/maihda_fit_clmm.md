# Fit a cumulative MAIHDA model via ordinal::clmm

Internal engine call for `fit_maihda(engine = "ordinal")`. Builds the
analytic sample (complete cases on the model variables) so the stored
`data` matches the rows clmm fits, then calls
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) with
`Hess = TRUE` (needed for the threshold standard errors). The analytic
frame is passed by NAME (bound in a private environment) so the call
clmm stores stays one line: ordinal's `print`/`summary` methods deparse
`call$data`, and embedding the frame there made printing an ordinal fit
dump the whole data set. Any formula
[`offset()`](https://rdrr.io/r/stats/offset.html) term is moved ahead of
the `(1 | stratum)` term first, because clmm reads the offset by
position in the bar-free variables list (see
`maihda_offset_before_bars`).

## Usage

``` r
maihda_fit_clmm(formula, data, family, dot_vals)
```

## Arguments

- formula:

  The resolved model formula (with `(1 | stratum)`).

- data:

  The data (after strata creation / response preparation).

- family:

  The cumulative family marker (link "logit" or "probit").

- dot_vals:

  Named list of evaluated `...` arguments forwarded to
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) (e.g.
  `nAGQ`, `control`).

## Value

A list with `model` (the `clmm` fit) and `data` (the analytic data frame
actually fitted).
