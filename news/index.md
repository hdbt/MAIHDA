# Changelog

## MAIHDA 0.1.8

### General Updates & New Features

- Added
  [`plot_prediction_deviation_panels()`](https://hdbt.github.io/MAIHDA/reference/plot_prediction_deviation_panels.md)
  function for visualizing predicted values and identifying deviant
  cases.
- Added
  [`plot_risk_vs_effect()`](https://hdbt.github.io/MAIHDA/reference/plot_risk_vs_effect.md)
  function to create a quadrant scatterplot comparing overall marginal
  predicted risk against pure intersectional effects.
- Added
  [`plot_effect_decomposition()`](https://hdbt.github.io/MAIHDA/reference/plot_effect_decomposition.md)
  function to visually decompose the total deviation from the overall
  mean into additive and intersectional components.
- Replaced the redundant “caterpillar” plot with the “predicted” plot in
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the
  interactive dashboard.
- Added automatic tertile binning (via an `autobin` parameter) for
  numeric grouping variables with more than 10 unique values in
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).
- Updated the interactive Shiny Dashboard
  ([`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
  to include the new visualizations and a toggle for auto-binning
  continuous strata variables.
- Added detection for binomial data.
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  will now automatically detect binomial outcomes and switch to the
  appropriate family.

### Bug Fixes

- **VPC/ICC Calculation Fix**: Corrected the residual variance
  estimation for binomial and ordinal models. The package now accurately
  applies the theoretical level-1 variance ($`\pi^2 / 3`$ for `"logit"`
  links and $`1`$ for `"probit"` links) internally when summarizing
  models or bootstrapping the variance partition coefficient, avoiding
  deflated VPC/ICC metrics.

## MAIHDA 0.1.7

CRAN release: 2026-04-05

### General Updates & New Features

- Added
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  function to sequentially estimate proportional change in variance
  (PCV) by adding predictors one-by-one.
- Added a fully-featured interactive Shiny Dashboard (via
  [`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
  for visual data exploration, model fitting, and performance
  visualization.
- Improved bootstrap methods for more efficient confidence interval
  estimation.
- Added missing documentation block for the `maihda_sim_data` dataset to
  resolve `R CMD check` warnings.
- Updated test suite setup: `tests/testthat.R` was modified to correctly
  use `test_check("MAIHDA")` instead of `shinytest2`.
- Added `importFrom(stats, as.formula)` for the `stepwise_pcv` function
  to prevent undefined warnings.
- Updated `introduction.Rmd` vignette: added standard CRAN installation
  instructions, and improved text clarity.

## MAIHDA 0.1.0

CRAN release: 2026-04-03

### Initial Release

- Initial CRAN submission
- Added
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  function for creating intersectional strata
- Added
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  function for fitting multilevel models with lme4 (default) or brms
  engines
- Added [`summary()`](https://rdrr.io/r/base/summary.html) function for
  variance partition and stratum estimates
- Added
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  function for individual and stratum-level predictions
- Added [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  function with three plot types:
  - Caterpillar plots of stratum random effects
  - Variance partition coefficient visualization
  - Observed vs. shrunken estimates comparison
- Added
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  function for comparing models with bootstrap confidence intervals
- Added comprehensive documentation and vignettes
- Added unit tests for core functionality

### Bug Fixes and Improvements

- Enhanced
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  to properly handle missing values (NA) in input variables:
  - Observations with missing values in any stratum variable are now
    assigned NA stratum
  - Missing values are no longer included as valid stratum categories
  - Added comprehensive tests for missing value handling
