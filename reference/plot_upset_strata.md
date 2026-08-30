# UpSet-style Predicted Stratum Plot

Composite alternative to the text-labelled `"predicted"` view that
replaces the long intersectional axis labels with an UpSet-style
category matrix. Three panels share one column order: a top bar of
intersection (stratum) sizes, a middle matrix encoding each stratum's
category on every dimension, and a bottom panel of predicted values with
conditional intervals. Columns are ordered by intersection size (largest
first) by default, or by any other `order_by` rule – `"predicted_desc"`
turns the view into the ranked caterpillar of
[`plot_predicted_strata()`](https://hdbt.github.io/MAIHDA/reference/plot_predicted_strata.md)
drawn against the category matrix instead of long text labels. Binary
0/1 (or logical) dimensions collapse to a single present/absent row;
multi-level factors get one row per level, and each column lights
exactly one dot per dimension.

## Usage

``` r
plot_upset_strata(
  object,
  summary_obj,
  n_strata,
  scale = c("response", "link"),
  highlight = NULL,
  only_flagged = FALSE,
  select = c("order", "deviation"),
  order_by = c("size", "predicted_desc", "stratum", "predicted_asc", "deviation"),
  quantity = c("predicted", "interaction")
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
  regardless; this governs the fill / the unflagged case. It controls
  *which* strata are shown, separately from how they are ordered for
  display.

- order_by:

  Left-to-right column order (**display-only**; does not change which
  strata are shown, nor any value): `"size"` (default, largest
  intersection first – the UpSet convention), `"stratum"` native stratum
  order, or `"predicted_desc"` / `"predicted_asc"` / `"deviation"`,
  which sort on the quantity the estimate panel actually shows and so
  follow `quantity`. `NULL` takes this view's default.

## Value

A patchwork object stacking the three panels.
