# MAIHDA 0.1.10

## New Features

* Added `compare_maihda_groups()` to compare intersectional inequality (VPC/ICC and between-/within-stratum variance) across levels of a higher-level grouping variable such as country, region, or survey wave. It fits a stratified MAIHDA model per group, by default using shared/global strata so VPCs are directly comparable, with optional per-group bootstrap confidence intervals.
* Added `plot_group_comparison()` to visualize the result either as a VPC-by-group forest plot or as stacked variance-composition bars.

## Bug Fixes

* Fixed parametric-bootstrap confidence intervals for VPC (`summary()`) and PVC (`calculate_pvc()`): failed `refit()` iterations were silently recorded as `0` instead of being dropped, biasing intervals toward zero and suppressing the high-failure-rate warning. Failed iterations are now excluded correctly.
* Corrected the Poisson VPC residual-variance approximation to `log(1 + 1/mu)` (Stryhn et al. 2006); the previous `1/mu` linearization biased the VPC downward for low-count outcomes.
* `plot_prediction_deviation_panels()` no longer draws zero-width "95% CI" bars when the underlying model does not supply standard errors (e.g. `lme4::merMod`); intervals are omitted instead of collapsed.


# MAIHDA 0.1.9

## Bug Fixes

* Clarified the Shiny dashboard PVC/HUD interpretation so negative PVC values are shown as variance unmasking rather than as unexplained interaction variance.
* Fixed the coverage workflow failure-artifact upload configuration.


# MAIHDA 0.1.8

## General Updates & New Features

* Added `plot_prediction_deviation_panels()` function for visualizing predicted values and identifying deviant cases.
* Added `plot_risk_vs_effect()` function to create a quadrant scatterplot comparing overall marginal predicted risk against pure intersectional effects.
* Added `plot_effect_decomposition()` function to visually decompose the total deviation from the overall mean into additive and intersectional components.
* Replaced the redundant "caterpillar" plot with the "predicted" plot in `plot()` and the interactive dashboard.
* Added automatic tertile binning (via an `autobin` parameter) for numeric grouping variables with more than 10 unique values in `make_strata()`.
* Updated the interactive Shiny Dashboard (`run_maihda_app()`) to include the new visualizations and a toggle for auto-binning continuous strata variables.
* Added detection for binomial data. `fit_maihda()` will now automatically detect binomial outcomes and switch to the appropriate family.

## Bug Fixes

* **VPC/ICC Calculation Fix**: Corrected the residual variance estimation for binomial and ordinal models. The package now accurately applies the theoretical level-1 variance ($\pi^2 / 3$ for `"logit"` links and $1$ for `"probit"` links) internally when summarizing models or bootstrapping the variance partition coefficient, avoiding deflated VPC/ICC metrics.

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
* Added `summary()` function for variance partition and stratum estimates
* Added `predict_maihda()` function for individual and stratum-level predictions
* Added `plot()` function with three plot types:
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
