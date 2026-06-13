# Fit a cumulative MAIHDA model via ordinal::clmm

Internal engine call for `fit_maihda(engine = "ordinal")`. Builds the
analytic sample (complete cases on the model variables) so the stored
`data` matches the rows clmm fits, then calls
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) with
`Hess = TRUE` (needed for the threshold standard errors).

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
