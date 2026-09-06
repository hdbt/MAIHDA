# Location linear predictor of a cumulative (clmm) fit

`predict.clmm` does not exist, so the location part \\\eta = x'\beta (+
u)\\ is built directly: the fixed design matrix is constructed with the
training data's factor levels AND transformation basis (so a
data-dependent term such as `scale(x)` uses the fit's centre and scale
rather than recomputing them from `newdata`) and multiplied by the
location coefficients `beta` (a clmm has *no* intercept column – it is
absorbed by the thresholds – so `beta`'s names select the right
columns), any formula offset term is evaluated on `newdata` and added
(so an offset-only null model still predicts its offset), and
`include_re` adds each row's stratum conditional mode (an unseen stratum
contributes 0 – the population-average fallback that
[`predict_maihda`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
only reaches when `allow_new_levels = TRUE`, having otherwise rejected
unseen strata upstream). Everything is on the latent (link) scale; map
through
[`maihda_ordinal_eta_to_score`](https://hdbt.github.io/MAIHDA/reference/maihda_ordinal_eta_to_score.md)
for the response-scale expected category score.

## Usage

``` r
maihda_clmm_linpred(object, newdata = NULL, include_re = TRUE)
```

## Arguments

- object:

  A `maihda_model` with engine `"ordinal"`.

- newdata:

  Data to predict for; defaults to the analytic data.

- include_re:

  Add the stratum random effect (conditional mode)?

## Value

A numeric vector of latent-scale location predictions.
