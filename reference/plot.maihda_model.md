# Plot MAIHDA Model Results

Creates various plots for visualizing MAIHDA model results including
variance partition coefficient comparisons, observed vs. shrunken
estimates, and predicted subgroup values with confidence intervals.

## Usage

``` r
# S3 method for class 'maihda_model'
plot(
  x,
  type = c("all", "vpc", "obs_vs_shrunken", "predicted", "upset", "effect_decomp",
    "prediction_deviation", "context_vpc", "vpc_trajectory", "trajectories"),
  summary_obj = NULL,
  n_strata = 50,
  highlight_interactions = FALSE,
  only_flagged = FALSE,
  highlight_by = c("flag", "rope"),
  rope = NULL,
  select = c("order", "deviation"),
  order_by = c("predicted_desc", "stratum", "predicted_asc", "deviation"),
  quantity = c("predicted", "interaction"),
  ...
)
```

## Arguments

- x:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- type:

  Character string specifying plot type:

  - "vpc": Variance partition coefficient visualization

  - "obs_vs_shrunken": Observed vs. shrunken stratum means. The y-axis
    (model-based estimate) includes the fixed effects, so for a
    covariate-adjusted model the distance from the diagonal reflects
    both shrinkage *and* covariate adjustment, not shrinkage alone; it
    is a pure shrinkage view only for an intercept-only (null) model

  - "predicted": Predicted values for each stratum with confidence
    intervals. By default the strata are ordered highest-to-lowest
    predicted value (a ranked caterpillar plot); see `order_by` to
    change or disable the ordering

  - "upset": UpSet-style alternative to `"predicted"` – an
    intersection-size bar, a category matrix encoding each stratum's
    level on every dimension, and the predicted-value panel, all sharing
    one column order. Replaces the long intersectional axis labels with
    the matrix. Binary 0/1 dimensions show as a single present/absent
    row; multi-level factors get one row per level

  - "effect_decomp": Visualizes additive vs intersectional deviation
    from global mean

  - "prediction_deviation": Detailed deviation panels for individuals or
    strata

  - "context_vpc": Stratum vs. context variance bars for a contextual
    cross-classified fit (`fit_maihda(context = )`); errors otherwise

  - "vpc_trajectory": Time-varying VPC/ICC curve for a **longitudinal**
    fit (`fit_maihda(id =, time =)`); errors otherwise. For a
    longitudinal model `"vpc"` and `"all"` also route here

  - "trajectories": Predicted per-stratum mean trajectories over time
    (longitudinal fits only)

  - "all": Generate all available plots (default if not specified)

- summary_obj:

  Optional maihda_summary object from
  [`summary()`](https://rdrr.io/r/base/summary.html). If NULL, will be
  computed.

- n_strata:

  Maximum number of strata to display in the predicted plot. When there
  are more strata than this, the first `n_strata` (in stratum order) are
  shown and the plot caption notes how many were omitted. Default is 50.
  Use NULL for all strata.

- highlight_interactions:

  Highlight the strata that carry a credibly non-zero intersectional
  interaction (from
  [`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md))
  on the BLUP-based views (`"effect_decomp"`, `"predicted"`,
  `"obs_vs_shrunken"`); other views ignore it. `FALSE` (default) off;
  `TRUE` computes the flags with
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  defaults; or pass a multiple-testing method such as `"BH"` or a
  `maihda_interactions` object to reuse a specific `conf_level`/
  `adjust`. For the pure-interaction reading the model should be the
  adjusted (or crossed-dimensions) model – e.g. via
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis, which routes these views to the adjusted model
  automatically. Which column of the screen drives the highlight set is
  governed by `highlight_by`.

- only_flagged:

  Show *only* the flagged strata rather than dimming the rest. `FALSE`
  (default) keeps every stratum (flagged ones highlighted); `TRUE`
  restricts the `"predicted"` and `"obs_vs_shrunken"` views to the
  highlighted strata (those carrying a credibly non-zero interaction, or
  – under `highlight_by = "rope"` – those classified ROPE-`"relevant"`),
  so a highlighted stratum is never hidden by the `n_strata` cap. When
  `TRUE` and `highlight_interactions` is left `FALSE`, the flags are
  computed with
  [`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  defaults; pass `highlight_interactions` to choose the
  `conf_level`/`adjust`. A captioned empty panel is returned when no
  stratum is flagged. It does not apply to `"effect_decomp"` (whose
  waterfall exists to show each flagged stratum's place in the *full*
  distribution); that view stays highlighted. Independently of this
  argument, whenever interactions are highlighted the `n_strata` cap on
  `"predicted"` becomes flag-aware: every flagged stratum is kept and
  the remaining slots are filled according to `select`.

- highlight_by:

  Which interaction-screen column defines the highlighted strata:
  `"flag"` (default) uses the zero-centred `flagged` column (credibly
  non-zero interaction), preserving the historical behaviour; `"rope"`
  uses the equivalence `decision` column, highlighting the strata
  classified `"relevant"` (interaction interval entirely outside the
  region of practical equivalence). `"rope"` requires a screen carrying
  a `decision` column: either pass `rope`, or supply a
  `maihda_interactions` object built with `rope`; otherwise it errors.

- rope:

  Equivalence region (a "smallest interaction of interest") forwarded to
  [`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  when the screen is computed here (i.e. when `highlight_interactions`
  is `TRUE` or a `p.adjust` method name). `NULL` (default) adds no
  equivalence classification; a single positive `d` means the symmetric
  region `c(-d, d)` on the latent (link) scale, or supply
  `c(lower, upper)`. Needed for `highlight_by = "rope"` unless the
  supplied `maihda_interactions` object already carries a `decision`
  column. Ignored when a precomputed `maihda_interactions` object is
  passed (its own `rope` is used).

- select:

  When the `n_strata` cap must drop strata, which to keep: `"order"`
  (default; the first n_strata in stratum order, the historical
  behaviour) or `"deviation"` (the n_strata furthest from the reference
  line – largest `|predicted - reference|`, so the most extreme strata
  in *both* directions). Applies to `"predicted"` and, for a
  longitudinal fit, `"trajectories"` (where it keeps the strata whose
  trajectories swing furthest from the population curve). Flagged strata
  are always kept; this governs the fill and the unflagged case.
  `select` changes *which* strata appear; their left-to-right display
  order is a separate choice governed by `order_by`.

- order_by:

  For `type = "predicted"`, how to order the strata that are displayed
  (**display-only** – it does not change *which* strata are shown, that
  is `n_strata`/`select`, nor the predicted values, intervals, reference
  line, or highlighted set): `"predicted_desc"` (default) orders labels
  from the highest predicted value to the lowest, `"stratum"` keeps the
  native stratum order (the previous behaviour), `"predicted_asc"`
  orders from lowest to highest, and `"deviation"` orders by largest
  absolute deviation from the reference line
  (`|predicted - reference|`). Ignored by the other plot types
  (`"upset"` orders by intersection size).

- quantity:

  For `type = "upset"`, which quantity the bottom panel shows:
  `"predicted"` (default) the stratum's predicted value (fixed + random
  effect) against the across-strata reference line, or `"interaction"`
  the stratum random effect (the BLUP) against zero – the deviation from
  the model's fixed prediction, which is the *pure* intersectional
  interaction when the dimension main effects are in the model (the
  adjusted model). Ignored by the other plot types.

- ...:

  Additional arguments (not currently used).

## Value

For a single `type`, a ggplot2 object that you can extend with the usual
`+` grammar (themes,
[`labs()`](https://ggplot2.tidyverse.org/reference/labs.html), added
layers, or a replacement fill/colour scale). Two types return a richer
object: `"prediction_deviation"` returns a patchwork of two panels
(theme every panel at once with `& theme_*()`). `type = "all"` returns a
named list of ggplot2 objects.

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)

# VPC plot
plot(model, type = "vpc")


# Single-type plots are ggplot objects -- restyle them with ggplot2:
plot(model, type = "vpc") +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "Variance partition, restyled")


# Generate all plots (a named list); pick one out to restyle it:
plots <- plot(model)





plots$predicted + ggplot2::theme_bw()

# }
```
