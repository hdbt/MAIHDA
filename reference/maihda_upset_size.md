# Recommended Figure Size for the UpSet Stratum Plot

Computes sensible `width` and `height` (in inches) for
`plot(object, type = "upset")`, so a knitr chunk or
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html) call
can size the figure to its content. The UpSet composite grows *taller*
with the number of matrix rows (one per binary 0/1 dimension, one per
level of a multi-level factor) and *wider* with the number of strata
columns shown, so a single fixed size tends to crop or stretch it –
particularly for multi-level designs (many rows) or a large `n_strata`
(many columns; UpSet is an inherently wide format).

## Usage

``` r
maihda_upset_size(object, n_strata = 50)
```

## Arguments

- object:

  A `maihda_model` from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  or a `maihda` analysis from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md); it must
  carry the per-dimension stratum table from
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).

- n_strata:

  Maximum number of strata the plot will show – pass the same value you
  give [`plot()`](https://rdrr.io/r/graphics/plot.default.html). `NULL`
  means all strata. Default 50.

## Value

A list with numeric `width` and `height` (inches) plus the `rows`
(matrix rows) and `cols` (strata shown) they derive from.

## See also

[`plot.maihda_model`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md)

## Examples

``` r
# \donttest{
strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
model <- fit_maihda(health_outcome ~ (1 | stratum), data = strata$data)
sz <- maihda_upset_size(model, n_strata = 30)
ggplot2::ggsave(
  tempfile(fileext = ".png"),
  plot(model, type = "upset", n_strata = 30),
  width = sz$width, height = sz$height)
# }
```
