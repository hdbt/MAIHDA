# Rows complete on every variable a brms model will use

Mirrors the rows brms retains after its own NA exclusion: complete on
the response and its addition-term variables (`y | trials(n)`), on the
fixed-effect terms – evaluated through the model frame, so a transformed
predictor whose transformation yields `NA`/`NaN` counts as incomplete –
and on every variable in the random-effect terms.

## Usage

``` r
maihda_brms_complete_rows(formula, data)
```

## Arguments

- formula:

  The (pre-weights-injection) model formula.

- data:

  The model data.

## Value

A logical vector over the rows of `data`.
