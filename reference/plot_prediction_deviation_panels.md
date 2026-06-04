# Plot Prediction Deviation Panels

Creates an advanced, publication-ready two-panel dashboard for
visualizing predicted values and highlighting the most notable cases or
strata. What "notable" means depends on the model type, and the labelled
points are *not* statistical outliers in the regression-diagnostic
sense:

- Gaussian and Poisson (and the ordinal `"expected_score"` mode): the
  cases/strata whose prediction sits furthest from the mean prediction
  (largest deviation), ranked by absolute deviation.

- Binomial: the cases/strata with the largest absolute deviance
  residual, i.e. where the observed 0/1 outcome is least consistent with
  the fitted probability (worst-fit points), ranked by \\\|deviance
  residual\|\\.

- Ordinal `"surprise"` mode: the cases/strata with the highest surprise
  \\-\log P(\text{observed category})\\, i.e. the least probable
  observations under the model.

## Usage

``` r
plot_prediction_deviation_panels(
  model,
  data = NULL,
  type = c("auto", "gaussian", "poisson", "binomial", "ordinal"),
  ordinal_mode = c("surprise", "expected_score"),
  top_n_labels = 5,
  strata_info = NULL
)
```

## Arguments

- model:

  A fitted model object (e.g., from \`lm()\`, \`glm()\`,
  \`MASS::polr()\`, or \`lme4::glmer()\`).

- data:

  The original data frame used to fit the model. If \`NULL\`, attempts
  to extract from the model.

- type:

  Model type: "auto" (default), "gaussian", "poisson", "binomial", or
  "ordinal".

- ordinal_mode:

  For ordinal models: "surprise" (default, based on observation
  probability) or "expected_score".

- top_n_labels:

  Number of points to label on the plot. The ranking metric depends on
  the model type (see Description): deviation from the mean prediction
  for Gaussian/Poisson and the ordinal expected-score mode, absolute
  deviance residual for binomial, and surprise for the ordinal surprise
  mode. Default is 5.

- strata_info:

  Optional data frame of strata labels, generally extracted from
  \`maihda_model\` objects.

## Value

A \`patchwork\` object containing two \`ggplot2\` panels.
