# Plot a MAIHDA Analysis

Dispatches each `type` to the model it is valid on. The VPC and
shrinkage views (`"vpc"`, `"obs_vs_shrunken"`, `"predicted"`) use the
**null** model. The additive-vs-intersectional views (`"effect_decomp"`,
`"prediction_deviation"`) use the **adjusted** model, whose fixed
effects carry the dimensions' additive part so the stratum random effect
is the pure interaction; with fewer than two dimensions (no adjusted
model) they fall back to the null model. Group types (`"group_vpc"`,
`"group_components"`, `"group_between_variance"`, `"group_pcv"`,
`"group_additive_share"`) use the group comparison when
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) was called
with a `group`.

## Usage

``` r
# S3 method for class 'maihda_analysis'
plot(
  x,
  type = "all",
  highlight_interactions = FALSE,
  only_flagged = FALSE,
  highlight_by = c("flag", "rope"),
  rope = NULL,
  select = c("order", "deviation"),
  order_by = c("predicted_desc", "stratum", "predicted_asc", "deviation"),
  ...
)
```

## Arguments

- x:

  A `maihda_analysis` object from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

- type:

  One of the model types ("all", "vpc", "obs_vs_shrunken", "predicted",
  "upset" (the UpSet-style alternative to "predicted"; forwards
  `quantity` via `...`), "effect_decomp", "prediction_deviation"), the
  contextual type ("context_vpc", a stratum-vs-context variance bar;
  requires `maihda(context = )`), a longitudinal type ("vpc_trajectory",
  "trajectories", "pcv_trajectory"; requires
  `decomposition = "longitudinal"`), or a group type ("group_vpc",
  "group_components", "group_between_variance", "group_pcv",
  "group_additive_share"). Default "all". For a longitudinal analysis
  "all" shows the VPC-over-time, the stratum trajectories, and the
  time-specific PCV.

- highlight_interactions:

  Highlight strata with a credibly non-zero intersectional interaction
  on the BLUP-based views (see
  [`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  and
  [`plot`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md)).
  `FALSE` (default), `TRUE` (computed from this analysis's adjusted /
  crossed-dimensions model), a multiple-testing method such as `"BH"`,
  or a `maihda_interactions` object. The flags are computed once from
  the correct (adjusted) model and reused across views.

- only_flagged:

  Show only the highlighted strata on the `"predicted"` and
  `"obs_vs_shrunken"` views instead of dimming the rest (see
  [`plot`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md)).
  `FALSE` (default) keeps every stratum; `TRUE` restricts those views to
  the highlighted strata (those carrying a credibly non-zero
  interaction, or – under `highlight_by = "rope"` – those classified
  ROPE-`"relevant"`) and, when `highlight_interactions` is left `FALSE`,
  turns the highlight on with this analysis's stored diagnostic.
  `"effect_decomp"` ignores it (it stays highlighted in context).

- highlight_by:

  Which interaction-screen column defines the highlighted strata:
  `"flag"` (default, the zero-centred `flagged` column) or `"rope"` (the
  equivalence `decision` column, highlighting the strata classified
  `"relevant"`). See
  [`plot`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md).
  `"rope"` requires `rope` (or a `highlight_interactions` object built
  with one).

- rope:

  Equivalence region forwarded to
  [`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  when computing the screen, so `highlight_by = "rope"` has a `decision`
  column to read: a single positive `d` for the symmetric region
  `c(-d, d)` on the latent (link) scale, or `c(lower, upper)`. `NULL`
  (default) adds no equivalence classification.

- select:

  When the `n_strata` cap drops strata on the `"predicted"` (or
  longitudinal `"trajectories"`) view, which to keep: `"order"`
  (default, first n_strata in stratum order) or `"deviation"` (the
  n_strata furthest from the reference, both tails). See
  [`plot`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md).

- order_by:

  Display order of the strata on the `"predicted"` view (display-only;
  does not change which strata are shown, nor any value):
  `"predicted_desc"` (default) highest predicted first, `"stratum"`
  native order, `"predicted_asc"` lowest first, or `"deviation"` largest
  `|predicted - reference|` first. See
  [`plot`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_model.md).

- ...:

  Additional arguments passed to the underlying plot method.

## Value

A ggplot2 object, or (for `type = "all"`) an invisible list of them.
