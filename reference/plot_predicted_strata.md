# Plot Predicted Stratum Values with Confidence Intervals

Plot Predicted Stratum Values with Confidence Intervals

## Usage

``` r
plot_predicted_strata(
  object,
  summary_obj,
  n_strata,
  scale = c("response", "link")
)
```

## Arguments

- object:

  A maihda_model object

- summary_obj:

  A maihda_summary object

- n_strata:

  Maximum number of strata to display

- scale:

  Prediction scale: "response" (default) or "link"

## Value

A ggplot2 object
