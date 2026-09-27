# Observed vs. Shrunken Estimates Plot

Observed vs. Shrunken Estimates Plot

## Usage

``` r
plot_obs_vs_shrunken(
  object,
  summary_obj,
  highlight = NULL,
  only_flagged = FALSE
)
```

## Arguments

- object:

  A maihda_model object

- summary_obj:

  A maihda_summary object

- highlight:

  Optional character vector of highlighted stratum ids (flagged, or
  ROPE-relevant under `highlight_by = "rope"`), with the
  interaction-screen parameters attached as attributes.

- only_flagged:

  When TRUE, show only the highlighted strata; a captioned empty panel
  is returned if none are.

## Value

A ggplot2 object

## Details

The x-axis is each stratum's raw observed mean; the y-axis is the
model-based stratum estimate, which includes the fixed-effect
contribution. Both are on the scale the response was *fitted* on, which
is what makes the \\y = x\\ diagonal meaningful: for a transformed
response such as `log(y) ~ x` the observed values are stratum means of
`log(y)`, not of `y`. To read the panel on the original scale, fit the
model on that scale; back-transforming only one axis would not give
comparable quantities. For an intercept-only (null) model the vertical
distance from the diagonal is pure shrinkage toward the grand mean. For
a covariate-adjusted model the model estimate also moves with the
stratum's covariate profile, so distance from the diagonal reflects
*both* shrinkage and covariate adjustment and should not be read as
shrinkage alone.
