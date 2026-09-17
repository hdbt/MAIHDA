# Effect Decomposition Plot

Decomposes the total deviation from the overall mean into the additive
(fixed) component and the intersectional (random) component for each
stratum.

## Usage

``` r
plot_effect_decomposition(
  object,
  summary_obj,
  top_n_labels = 10,
  highlight = NULL
)
```

## Arguments

- object:

  A maihda_model object

- summary_obj:

  A maihda_summary object

- top_n_labels:

  Number of most extreme strata to label

- highlight:

  Optional character vector of stratum ids to highlight. When supplied,
  labels are restricted to these strata rather than the most extreme
  overall deviations.

## Value

A ggplot2 object

## Details

Only the intersectional random effects enter the deviations and the
overall mean: the stratum effect and, for a crossed-dimensions fit, the
dimension effects. A contextual random effect, or any other grouping, is
part of neither component and is excluded.
