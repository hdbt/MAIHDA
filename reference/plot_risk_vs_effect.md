# Mean Prediction vs. Stratum Random Effect Plot

Creates a quadrant scatterplot comparing each stratum's mean predicted
outcome against its stratum random effect (shrunken between-stratum
deviation). Points represent strata. Whether a higher predicted value is
"worse" or "better" depends on the outcome, so the axes are not framed
as risk. The random effect equals the *pure* intersectional
(interaction) component only when the additive main effects of the
strata variables are included in the model; otherwise it also absorbs
those omitted main effects.

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
