# Plot Predicted Stratum Values with Confidence Intervals

Plot Predicted Stratum Values with Confidence Intervals

## Usage

``` r
plot_predicted_strata(
  object,
  summary_obj,
  n_strata,
  scale = c("response", "link"),
  highlight = NULL,
  only_flagged = FALSE,
  select = c("order", "deviation")
)
```

## Arguments

- object:

  A maihda_model object

- summary_obj:

  A maihda_summary object

- n_strata:

  Maximum number of strata to display (the first n_strata, in stratum
  order)

- scale:

  Prediction scale: "response" (default) or "link"

- highlight:

  Optional character vector of highlighted stratum ids – the flagged
  strata, or the ROPE-relevant strata under `highlight_by = "rope"` –
  with the interaction-screen parameters attached as attributes
  (including `highlight_by` and `rope`).

- only_flagged:

  When TRUE, show only the highlighted strata (those in `highlight`); a
  captioned empty panel is returned if none are.

- select:

  When the `n_strata` cap drops strata, which to keep: `"order"`
  (default) the first n_strata in stratum order, or `"deviation"` the
  n_strata furthest from the reference line (largest
  `|predicted - reference|`, so both tails). Flagged strata are kept
  regardless; this governs the fill / the unflagged case. The displayed
  x-axis stays in stratum order either way.

## Value

A ggplot2 object
