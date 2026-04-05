# MAIHDA 0.1.7

## General Updates & New Features

* Added `stepwise_pcv()` function to sequentially estimate proportional change in variance (PCV) by adding predictors one-by-one.
* Added a fully-featured interactive Shiny Dashboard (via `run_maihda_app()`) for visual data exploration, model fitting, and performance visualization.
* Improved bootstrap methods for more efficient confidence interval estimation.
* Added missing documentation block for the `maihda_sim_data` dataset to resolve `R CMD check` warnings.
* Updated test suite setup: `tests/testthat.R` was modified to correctly use `test_check("MAIHDA")` instead of `shinytest2`.
* Added `importFrom(stats, as.formula)` for the `stepwise_pcv` function to prevent undefined warnings.
* Updated `introduction.Rmd` vignette: added standard CRAN installation instructions, and improved text clarity.

# MAIHDA 0.1.0

## Initial Release

* Initial CRAN submission
* Added `make_strata()` function for creating intersectional strata
* Added `fit_maihda()` function for fitting multilevel models with lme4 (default) or brms engines
* Added `summary_maihda()` function for variance partition and stratum estimates
* Added `predict_maihda()` function for individual and stratum-level predictions
* Added `plot_maihda()` function with three plot types:
  * Caterpillar plots of stratum random effects
  * Variance partition coefficient visualization
  * Observed vs. shrunken estimates comparison
* Added `compare_maihda()` function for comparing models with bootstrap confidence intervals
* Added comprehensive documentation and vignettes
* Added unit tests for core functionality

## Bug Fixes and Improvements

* Enhanced `make_strata()` to properly handle missing values (NA) in input variables:
  * Observations with missing values in any stratum variable are now assigned NA stratum
  * Missing values are no longer included as valid stratum categories
  * Added comprehensive tests for missing value handling
