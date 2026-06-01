# Generate Ternary Plot from MAIHDA Model

Generate Ternary Plot from MAIHDA Model

## Usage

``` r
maihda_ternary_plot(model, summary_obj = NULL, ...)
```

## Arguments

- model:

  A fitted MAIHDA model.

- summary_obj:

  Optional output from `summary_maihda`.

- ...:

  Additional arguments passed to `compute_maihda_ternary_data` and
  [`plot.maihda_ternary`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_ternary.md).

## Value

A list containing `data` and `plot`.
