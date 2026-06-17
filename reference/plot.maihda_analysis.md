# Plot a MAIHDA Analysis

Dispatches each `type` to the model it is valid on. The VPC and
shrinkage views (`"vpc"`, `"obs_vs_shrunken"`, `"predicted"`) use the
**null** model. The additive-vs-intersectional views
(`"risk_vs_effect"`, `"effect_decomp"`, `"ternary"`,
`"prediction_deviation"`) use the **adjusted** model, whose fixed
effects carry the dimensions' additive part so the stratum random effect
is the pure interaction; with fewer than two dimensions (no adjusted
model) they fall back to the null model. Group types (`"group_vpc"`,
`"group_components"`, `"group_between_variance"`, `"group_pcv"`) use the
group comparison when
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) was called
with a `group`.

## Usage

``` r
# S3 method for class 'maihda_analysis'
plot(x, type = "all", highlight_interactions = FALSE, ...)
```

## Arguments

- x:

  A `maihda_analysis` object from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

- type:

  One of the model types ("all", "vpc", "obs_vs_shrunken", "predicted",
  "risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation"),
  the contextual type ("context_vpc", a stratum-vs-context variance bar;
  requires `maihda(context = )`), a longitudinal type ("vpc_trajectory",
  "trajectories", "pcv_trajectory"; requires
  `decomposition = "longitudinal"`), or a group type ("group_vpc",
  "group_components", "group_between_variance", "group_pcv"). Default
  "all". For a longitudinal analysis "all" shows the VPC-over-time, the
  stratum trajectories, and the time-specific PCV.

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

- ...:

  Additional arguments passed to the underlying plot method.

## Value

A ggplot2 object, or (for `type = "all"`) an invisible list of them.
