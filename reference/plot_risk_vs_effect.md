# Risk vs. Intersectional Effect Plot

Creates a quadrant scatterplot comparing overall marginal predicted risk
against pure intersectional effects (shrunken residuals). Points
represent strata.

## Usage

``` r
plot_risk_vs_effect(object, summary_obj, top_n_labels = 10)
```

## Arguments

- object:

  A maihda_model object

- summary_obj:

  A maihda_summary object

- top_n_labels:

  Number of most extreme strata to label (by absolute effect size)

## Value

A ggplot2 object
