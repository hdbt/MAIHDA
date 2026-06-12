# Inject sampling weights into a brms formula

Rewrites `y ~ ...` as `y | weights(w) ~ ...`. An existing addition term
(e.g. an aggregated-binomial `y | trials(n)`) is extended with
`+ weights(w)`; a formula that already carries a
[`weights()`](https://rdrr.io/r/stats/weights.html) addition term is
rejected (the two weight specifications would conflict).

## Usage

``` r
maihda_brms_weights_formula(formula, wcol)
```

## Arguments

- formula:

  The model formula.

- wcol:

  Name of the (normalized) weight column.

## Value

The rewritten formula (same environment).
