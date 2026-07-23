# Predict method for maihda_model objects

S3 [`predict`](https://rdrr.io/r/stats/predict.html) method for MAIHDA
models: a thin alias of
[`predict_maihda`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md),
so `predict(model, type = "strata")` and
`predict_maihda(model, type = "strata")` are interchangeable. See
[`predict_maihda`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
for the full documentation of the arguments and return value.

## Usage

``` r
# S3 method for class 'maihda_model'
predict(
  object,
  newdata = NULL,
  type = c("individual", "strata", "response", "link"),
  scale = c("response", "link"),
  allow_new_levels = FALSE,
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

  - "strata": Stratum-level predictions (random effects only). For a
    longitudinal (growth-curve) fit a stratum is a *trajectory*, so this
    returns the per-stratum trajectory parameters (baseline deviation,
    random intercept and random slope(s)) rather than a single random
    effect. A *non-longitudinal* fit whose stratum random effects
    include random slopes (a hand-written `(1 + x | stratum)`) is an
    error here, as in
    [`summary()`](https://rdrr.io/r/base/summary.html): the single value
    this would return is the intercept alone – the stratum effect only
    where the slope variables are zero.

  For backward compatibility, "link" or "response" may also be passed
  here and will be interpreted as individual-level predictions on that
  scale.

- scale:

  Character string specifying the prediction scale for individual-level
  predictions: "response" (default) or "link". For a cumulative
  (ordinal) model the "link" scale is the latent location \\\eta\\ and
  the "response" scale is the *expected category score* \\\sum_k k P(Y =
  k)\\ (categories scored 1..K in their declared order). For an
  aggregated-binomial fit (an lme4 `cbind(success, failure)` or a brms
  `success | trials(n)`) the "response" scale is the per-trial
  *probability* on both engines (not the expected success count).

- allow_new_levels:

  Logical. By default (`FALSE`) a stratum in `newdata` that the model
  never saw – whether supplied directly as a `stratum` column or rebuilt
  from the grouping variables – is an error, for every engine, matching
  lme4's default. Set `TRUE` to instead predict unseen strata with the
  stratum random effect dropped (treated as zero), while keeping any
  *other* random effect the row participates in (e.g. a contextual
  `(1 | school)` intercept from `fit_maihda(context = )`, or a
  longitudinal growth term) – the same behaviour as lme4's
  `allow.new.levels`, which zeroes only the unseen level's effect and
  keeps seen ones. For the usual single-stratum model the stratum is the
  only random effect, so this is the *population-average*
  (fixed-effects-only) prediction. This affects `type = "individual"`
  only: a stratum-level prediction (`type = "strata"`) has no random
  effect to report for an unseen stratum, so unseen strata remain an
  error there regardless.

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

## See also

[`predict_maihda`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
