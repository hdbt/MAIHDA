# MAIHDA 0.1.11

## New Features

* Added `maihda()`, a single high-level entry point that runs the standard workflow in one call: it fits the model, summarises the VPC/ICC and variance components, and -- when a `group` is supplied -- also compares intersectional inequality across that grouping variable. It returns one consistent `maihda_analysis` object (with `print`, `summary`, and `plot` methods); the `groups` slot is simply `NULL` when no grouping variable is given.
* Added the `maihda_country_data` dataset (OECD PISA 2018, accessed via the `learningtower` package): 3,600 students across six countries with gender x socioeconomic-status strata and mathematics-achievement outcomes. It showcases `compare_maihda_groups()` / `maihda(group = "country")`, since intersectional inequality (VPC/ICC) genuinely differs across the countries.

## Improvements

* Plotting is now unified under the base `plot()` generic. `compare_maihda()` and `compute_maihda_ternary_data()` now return classed objects, so `plot()` dispatches automatically:
  * `plot()` on a `compare_maihda()` result (was `plot_comparison()`)
  * `plot()` on a `compare_maihda_groups()` result, with `type = "vpc"`/"components" (was `plot_group_comparison()`)
  * `plot()` on a `compute_maihda_ternary_data()` result (was `plot_maihda_ternary()`)
* The old `plot_comparison()`, `plot_group_comparison()`, and `plot_maihda_ternary()` functions still work but are **deprecated** and emit a one-time warning pointing to `plot()`.
* `fit_maihda()` now records lme4 fit-quality diagnostics (singular fit and convergence warnings) on the model object; `print()` on the model and `summary()` surface a "Fit diagnostics" note so a boundary/non-converged fit -- which makes the VPC/PCV unreliable -- is no longer silent.
* `compare_maihda_groups()` now raises a single aggregated warning naming any group whose lme4 fit was singular or failed to converge (per-group fits on small strata are where this is most likely), so an unreliable per-group VPC -- often pinned at 0 by a boundary fit -- is no longer silent.
* The `"risk_vs_effect"` plot no longer frames the outcome axis as "risk". A higher predicted value is not universally "bad" (it depends on the outcome), so the axis and title now read as a neutral "mean predicted value/probability", with a note that the direction depends on the outcome.
* Clarified the documentation of `plot_prediction_deviation_panels()` to match the implementation: the labelled points use a per-type metric -- deviation from the mean prediction for Gaussian/Poisson (and the ordinal expected-score mode), the absolute deviance residual for binomial, and surprise for the ordinal surprise mode -- and are not statistical outliers or model-misfit "deviants".
* Clarified that the per-group VPC/ICC in `compare_maihda_groups()` is the between-stratum *share* of variance, which can differ across groups because of the residual variance as well as the between-stratum variance. The documentation now points to the `var_between` column for comparing the absolute amount of intersectional variation, and notes that overlap of separate per-group intervals is not a valid test of whether two groups' VPCs differ.
* Clarified the PCV documentation (`calculate_pvc()`, `stepwise_pcv()`, and the print method): the PCV is a model-dependent change in between-stratum variance and equals variance "explained" only when the second model nests the first; the stepwise PCV is order-dependent and not a variable's unique contribution.

## Bug Fixes

* VPC/ICC is no longer reported for a Gaussian model fit with a non-identity link (e.g. `gaussian(link = "log")`). The residual variance is on the response scale while the between-stratum variance is on the link scale, so `summary()`, `compare_maihda()`, and `calculate_pvc()` now raise a clear error instead of silently returning an invalid variance partition.
* Binary-outcome auto-detection now keys off the analytic model frame -- after applying covariate transformations (e.g. `log(x)`), dropping rows with missing values, and applying any `subset=` -- instead of the raw outcome column. An outcome that is only 0/1 once excluded rows are removed is now correctly fit with `family = "binomial"`.
* `fit_maihda()` now builds the `lmer()`/`glmer()`/`brm()` call explicitly when forwarding `...`, so data-masked engine arguments such as `weights=`, `subset=`, and `offset=` work instead of failing with "..1 used in an incorrect context".
* `plot_prediction_deviation_panels()` now plots Poisson/count models on the response (expected-count) scale with count labels, rather than routing them through the Gaussian link-scale branch.
* `compare_maihda_groups(min_group_n = ...)` now guards the analytic sample size (the rows the model actually fits) rather than the raw group row count, so a group with enough raw rows but a tiny usable sample is skipped instead of being fit on a handful of observations.
* `n_boot` for bootstrap intervals must now be at least 10 (the minimum number of successful refits an interval requires); an unusably small `n_boot` fails immediately with a clear message instead of only erroring after the bootstrap runs.


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
