# Plot MAIHDA Model Results

Creates various plots for visualizing MAIHDA model results including
variance partition coefficient comparisons, observed vs. shrunken
estimates, and predicted subgroup values with confidence intervals.

## Usage

``` r
# S3 method for class 'maihda_model'
plot(
  x,
  type = c("all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect",
    "effect_decomp", "ternary", "prediction_deviation", "context_vpc", "vpc_trajectory",
    "trajectories"),
  summary_obj = NULL,
  n_strata = 50,
  highlight_interactions = FALSE,
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
    intervals

  - "risk_vs_effect": Quadrant scatterplot of each stratum's mean
    predicted outcome against its random effect

  - "effect_decomp": Visualizes additive vs intersectional deviation
    from global mean

  - "ternary": Ternary diagnostic of the relative additive,
    intersectional, and uncertainty signals per stratum (a
    normalized-magnitude diagnostic, not a variance decomposition)

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
  automatically.

- ...:

  Additional arguments (not currently used).

## Value

For a single `type`, a ggplot2 object that you can extend with the usual
`+` grammar (themes,
[`labs()`](https://ggplot2.tidyverse.org/reference/labs.html), added
layers, or a replacement fill/colour scale). Two types return a richer
object: `"prediction_deviation"` returns a patchwork of two panels
(theme every panel at once with `& theme_*()`), and `"ternary"` returns
a ggtern object (use the `ggtern::theme_*()` family rather than the
standard ggplot2 themes). `type = "all"` returns a named list of ggplot2
objects.

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





#> Warning: Removing Layer 2 ('PositionNudge'), as it is not an approved position (for ternary plots) under the present ggtern package.


plots$predicted + ggplot2::theme_bw()

# }
```
