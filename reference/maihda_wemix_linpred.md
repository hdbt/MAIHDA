# Fixed-part (and optionally full) linear predictor of a wemix fit

WeMix's own [`predict()`](https://rdrr.io/r/stats/predict.html) method
needs the grouping structure re-resolved and offers no fixed-only form,
so predictions are built directly from the coefficient vector and the
stored stratum effects: the fixed design matrix is constructed with the
training data's factor levels and multiplied by `coef`, and `include_re`
adds each row's stratum effect (conditional mode; an unseen stratum
contributes 0 – the population-average fallback that
[`predict_maihda`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
only reaches when `allow_new_levels = TRUE`, having otherwise rejected
unseen strata upstream). Everything is on the link scale.

## Usage

``` r
maihda_wemix_linpred(object, newdata = NULL, include_re = TRUE)
```

## Arguments

- object:

  A `maihda_model` with engine `"wemix"`.

- newdata:

  Data to predict for; defaults to the analytic data.

- include_re:

  Add the stratum random effect (conditional mode)?

## Value

A numeric vector of link-scale predictions.
