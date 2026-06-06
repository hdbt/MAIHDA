# Predict from MAIHDA Model

Makes predictions from a fitted MAIHDA model, either at the stratum
level or individual level.

## Usage

``` r
predict_maihda(
  object,
  newdata = NULL,
  type = c("individual", "strata", "response", "link"),
  scale = c("response", "link"),
  ...
)
```

## Arguments

- object:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- newdata:

  Optional data frame for making predictions. If NULL, uses the original
  data from model fitting.

- type:

  Character string specifying prediction type:

  - "individual": Individual-level predictions including random effects

  - "strata": Stratum-level predictions (random effects only)

  For backward compatibility, "link" or "response" may also be passed
  here and will be interpreted as individual-level predictions on that
  scale.

- scale:

  Character string specifying the prediction scale for individual-level
  predictions: "response" (default) or "link".

- ...:

  Additional arguments passed to predict method of underlying model.

## Value

Depending on type:

- For "individual": A numeric vector of predicted values on the
  requested scale

- For "strata": A data frame with stratum ID and predicted random
  effect. When `newdata` is supplied, the result is restricted to the
  strata present in `newdata` (and a stratum the model never saw is an
  error, as for "individual"); when `newdata` is `NULL`, every training
  stratum is returned.

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)

# Individual predictions
pred_ind <- predict_maihda(model, type = "individual")

# Stratum predictions
pred_strata <- predict_maihda(model, type = "strata")
# }
```
