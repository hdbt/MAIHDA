# Plot Prediction Deviation Panels

Creates an advanced, publication-ready two-panel dashboard for
visualizing predicted values and identifying deviant cases in linear,
binomial, or ordinal models.

## Usage

``` r
plot_prediction_deviation_panels(
  model,
  data = NULL,
  type = c("auto", "gaussian", "binomial", "ordinal"),
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

  Model type: "auto" (default), "gaussian", "binomial", or "ordinal".

- ordinal_mode:

  For ordinal models: "surprise" (default, based on observation
  probability) or "expected_score".

- top_n_labels:

  Number of extreme/deviant cases to label on the plot. Default is 5.

- strata_info:

  Optional data frame of strata labels, generally extracted from
  \`maihda_model\` objects.

## Value

A \`patchwork\` object containing two \`ggplot2\` panels.
