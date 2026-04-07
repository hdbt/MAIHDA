# Effect Decomposition Plot

Decomposes the total deviation from the overall mean into the additive
(fixed) component and the intersectional (random) component for each
stratum.

## Usage

``` r
plot_effect_decomposition(object, summary_obj, top_n_labels = 10)
```

## Arguments

- object:

  A maihda_model object

- summary_obj:

  A maihda_summary object

- top_n_labels:

  Number of most extreme strata to label

## Value

A ggplot2 object
