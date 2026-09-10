# Expanded cut points of a cumulative (clmm) fit

The \\K-1\\ cut points of a \\K\\-category cumulative model, whatever
threshold structure was fitted.

## Usage

``` r
maihda_clmm_cutpoints(model)
```

## Arguments

- model:

  A fitted `clmm`.

## Value

A numeric vector of \\K-1\\ increasing cut points, named `"1|2"`,
`"2|3"`, ... when the fit supplies those labels.

## Details

`clmm`'s `$alpha` holds the *free* threshold parameters, which equal the
cut points only under the default `threshold = "flexible"`. A structured
threshold request (`"equidistant"`, `"symmetric"`, `"symmetric2"`, all
part of `clmm`'s documented API and reachable through `fit_maihda`'s
`...`) leaves `$alpha` holding a shorter, differently-parameterised
vector – `threshold.1` and `spacing` for the equidistant case – so
reading it as the cut points understates the number of categories and
mislocates every one of them. The expanded cut points are stored on the
fit as `$Theta` (a \\1 \times (K-1)\\ matrix), equivalently
`$tJac %*% $alpha`; `$tJac` is the identity under `"flexible"`, so this
returns `$alpha` unchanged there.
