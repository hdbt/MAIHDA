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
  residual, i.e. where the observed 0/1 outcome – for an aggregated
  binomial (`cbind(successes, failures)`, a proportion with its trial
  counts as prior weights, or a brms `y | trials(n)` fit), the observed
  successes out of trials – is least consistent with the fitted
  probability (worst-fit points), ranked by \\\|deviance residual\|\\.
  As in R's binomial family, a prior weight counts as that many trials,
  so a 0/1 outcome with weights is scored as that many successes or
  failures. Every row's outcome is coded as the model coded it – by
  label, whatever the order of a factor's levels in `data` – and a
  factor with more than two levels, as in R's binomial family, as its
  first level against all the others. Predictions are per-trial
  probabilities for every binomial fit. A row whose outcome or trial
  count is unknown has no residual: it is left out of the stratum
  residual means, is never labelled, and is drawn as a hollow point.

- Ordinal `"surprise"` mode: the cases/strata with the highest surprise
  \\-\log P(\text{observed category})\\, i.e. the least probable
  observations under the model. It needs the observed response in
  `data`: a row whose response is missing is left out, and one whose
  category is not among the model's fitted categories is left out with a
  warning.

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

  A fitted model object (e.g., from
  [`lm()`](https://rdrr.io/r/stats/lm.html),
  [`glm()`](https://rdrr.io/r/stats/glm.html),
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html), or
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)).

- data:

  The original data frame used to fit the model. If `NULL`, attempts to
  extract from the model. When supplied, everything is computed for its
  own rows: their predictions, their binomial residuals (from their own
  outcomes and trial counts) and their weights in the stratum summaries.
  A row's weight is read from its column in `data` – the sampling-weight
  column, or the column a bare model's `weights =` argument names
  ([`glm()`](https://rdrr.io/r/stats/glm.html), `lmer()`,
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html), ...) – or
  else from the fitted row with the same row name and prediction (a
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  model stores its weights as values). If any plotted row's weight
  cannot be found, every row is weighted equally, with a warning.

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
  `maihda_model` objects.

## Value

A patchwork object containing two ggplot2 panels.

## Details

The vertical bars are 95% intervals for the plotted values. On a stratum
panel each bar is an interval for the stratum's weighted mean
prediction, the point drawn: posterior quantiles of that mean for a brms
fit, and for an lme4 fit a delta-method interval from the joint
covariance of the fixed and random effects (Henderson's mixed-model
equations), formed on the logit scale for a probability and the log
scale for a count, so that a bar stays inside the range of the response.
A case-level panel of an lme4 fit draws the same interval for each row.
The lme4 interval treats the variance parameters as known: it runs
narrow when strata are few or small or an outcome is rare, and a random
effect whose variance is estimated at zero adds no uncertainty to it.
Other models draw case-level bars from the standard errors of their
[`predict()`](https://rdrr.io/r/stats/predict.html) method and no
stratum bars; WeMix and ordinal fits draw none.

## References

Henderson, C. R. (1975). Best linear unbiased estimation and prediction
under a selection model. *Biometrics*, 31(2), 423-447.
[doi:10.2307/2529430](https://doi.org/10.2307/2529430)

Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*,
67(1), 1-48.
[doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01)

Kackar, R. N., & Harville, D. A. (1984). Approximations for standard
errors of estimators of fixed and random effects in mixed linear models.
*Journal of the American Statistical Association*, 79(388), 853-862.
[doi:10.1080/01621459.1984.10477102](https://doi.org/10.1080/01621459.1984.10477102)
