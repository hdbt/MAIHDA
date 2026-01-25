# Changelog

## MAIHDA 0.1.0

### Initial Release

- Initial CRAN submission
- Added
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  function for creating intersectional strata
- Added
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  function for fitting multilevel models with lme4 (default) or brms
  engines
- Added
  [`summary_maihda()`](https://hdbt.github.io/MAIHDA/reference/summary_maihda.md)
  function for variance partition and stratum estimates
- Added
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  function for individual and stratum-level predictions
- Added
  [`plot_maihda()`](https://hdbt.github.io/MAIHDA/reference/plot_maihda.md)
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
