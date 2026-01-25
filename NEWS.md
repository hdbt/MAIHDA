# MAIHDA 0.1.0

## Initial Release

* Initial CRAN submission
* Added `make_strata()` function for creating intersectional strata
* Added `fit_maihda()` function for fitting multilevel models with lme4 (default) or brms engines
* Added `summary_maihda()` function for variance partition and stratum estimates
* Added `predict_maihda()` function for individual and stratum-level predictions
* Added `plot_maihda()` function with three plot types:
  - Caterpillar plots of stratum random effects
  - Variance partition coefficient visualization
  - Observed vs. shrunken estimates comparison
* Added `compare_maihda()` function for comparing models with bootstrap confidence intervals
* Added comprehensive documentation and vignettes
* Added unit tests for core functionality

## Bug Fixes and Improvements

* Enhanced `make_strata()` to properly handle missing values (NA) in input variables:
  - Observations with missing values in any stratum variable are now assigned NA stratum
  - Missing values are no longer included as valid stratum categories
  - Added comprehensive tests for missing value handling
