# Plot a MAIHDA Analysis

Dispatches to the model's plots (see
[`plot.maihda_model`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md))
for the model-level `type`s, and to the group comparison for
`"group_vpc"` and `"group_components"` when
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) was called
with a `group`.

## Usage

``` r
# S3 method for class 'maihda_analysis'
plot(x, type = "all", ...)
```

## Arguments

- x:

  A `maihda_analysis` object from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

- type:

  One of the
  [`plot.maihda_model`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md)
  types ("all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect",
  "effect_decomp", "ternary", "prediction_deviation") or a group type
  ("group_vpc", "group_components"). Default "all".

- ...:

  Additional arguments passed to the underlying plot method.

## Value

A ggplot2 object, or (for `type = "all"`) an invisible list of them.
