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
* Clarified the PCV documentation (`calculate_pvc()`, `stepwise_pcv()`, and the print method): the PCV is a model-dependent change in between-stratum variance and equals variance "explained" only when the second model nests the first; the stepwise PCV is order-dependent and not a variable's unique contribution. The vignette and Shiny app no longer describe PCV as variance causally "explained" by main effects or treat a negative PCV as evidence of hidden structural inequality.
* Corrected the `summary()` VPC/ICC documentation: the denominator is the total unexplained variance, which includes the variance of any additional random effects (not just between-stratum + residual), and a note on the weighted-Gaussian level-1 variance was added.
* Documented which families the MAIHDA variance summaries support (`gaussian("identity")`, binomial/Bernoulli with logit/probit, `poisson("log")`); other families such as `Gamma(link = "log")` will fit but `summary()`/VPC/PCV will stop with a clear "not implemented" error rather than returning an invalid number.
* README clarifications: a note that the CRAN release can lag this repository (so the newest features may require the GitHub development version), and the interactive-dashboard dependency list now includes `future`, `promises`, and `haven` (for SPSS/Stata uploads).
* Completed the PCV wording cleanup in the remaining vignette and Shiny-app text (the prior pass missed several spots), and softened over-strong app labels: the quadrant plot is "Mean Prediction vs. Stratum Effect" (not "Risk vs. Intersectional Effect"), the cumulative-PCV chart is "Change in Between-Stratum Variance" (not "Variance Explained"), and "deviant strata" is now "most extreme strata".
* Removed the stale checked-in rendered vignette HTML (`vignettes/*.html`); these are build artifacts generated from the `.Rmd` sources and had drifted out of sync with the corrected text. They are now git-ignored and added to `.Rbuildignore`, so a locally rendered HTML is never shipped in the tarball and R CMD build regenerates `inst/doc` from the `.Rmd`.
* Fixed an invalid documented URL flagged by `R CMD check --as-cran`: the `maihda_country_data` `@source` linked to `https://www.oecd.org/pisa/data/`, which returns HTTP 403. It now links to the CRAN page of the `learningtower` package (the reproducible access route used to build the dataset), keeping the OECD PISA 2018 attribution in the text.
* The `data-raw/maihda_health_data.R` regeneration script no longer calls `install.packages("NHANES")` as a side effect; like the country-data script it now stops with a clear message asking the developer to install the dependency.

## Performance

* `make_strata()` (and the prediction-time stratum lookup) now matches rows to strata with a single vectorised, collision-safe key match instead of an O(rows x strata x variables) row-by-row scan, so it scales to large survey datasets. Behaviour is unchanged, including the correct handling of values that contain the stratum-label separator.

## Bug Fixes

* VPC/ICC is no longer reported for a Gaussian model fit with a non-identity link (e.g. `gaussian(link = "log")`). The residual variance is on the response scale while the between-stratum variance is on the link scale, so `summary()`, `compare_maihda()`, and `calculate_pvc()` now raise a clear error instead of silently returning an invalid variance partition.
* Binary-outcome auto-detection now keys off the analytic model frame -- after applying covariate transformations (e.g. `log(x)`), dropping rows with missing values, dropping rows with a missing prior weight, and applying any `subset=` (including negative/character row indices) -- instead of the raw outcome column. An outcome that is only 0/1 once excluded rows are removed is now correctly fit with `family = "binomial"`.
* Data-masked engine arguments forwarded through `...` (e.g. `weights=`, `subset=`, `offset=`) now work at any nesting depth, including through the `maihda()` and `compare_maihda_groups()` wrappers. Arguments are captured as quosures and resolved with the data mask (`rlang::eval_tidy`), fixing the previous "object not found" / "..1 used in an incorrect context" failures (a direct `fit_maihda()` call worked, but the wrappers did not).
* A binary outcome is now recoded to 0/1 in a way that no longer breaks a `subset=` expression referencing the original response labels (e.g. `subset = y %in% c("no", "yes")`): the subset is evaluated against the original response before recoding.
* `compare_maihda_groups()` now slices forwarded `weights=`/`subset=`/`offset=` to each group's rows before fitting, not just for the row-count guard. An external (non-column) `weights` vector previously failed every group with a length mismatch, and an external `subset` vector could be recycled onto the wrong rows of later groups; both now align correctly per group.
* The Gaussian VPC/ICC now accounts for prior `weights`. With weights the per-observation residual variance is `sigma^2 / w_i`, so the level-1 variance reported is the mean conditional residual variance `mean(sigma^2 / w_i)` rather than a single `sigma^2`; unweighted models are unchanged.
* `plot_effect_decomposition()` now uses the stratum random effect (BLUP) itself as the intersectional component instead of (full prediction - fixed prediction). With additional random effects such as `(1 | site)` the latter wrongly absorbed those effects into the stratum component.
* The stratum-level "surprise" in `plot_prediction_deviation_panels()` (ordinal, surprise mode) is now the average per-observation surprise, `mean(-log(p))` (log loss), instead of `-log(mean(p))`, which could change stratum rankings.
* The Shiny app (`run_maihda_app()`) now also auto-detects a numeric 0/1 outcome and fits it as binomial under the default family, matching the core API, instead of silently fitting a linear probability model.
* `compare_maihda_groups()` now warns when groups end up with different populated strata even under `shared_strata = TRUE`, since their VPCs are then estimated over different stratum support and are not strictly directly comparable.
* `plot_prediction_deviation_panels()` now plots Poisson/count models on the response (expected-count) scale with count labels, rather than routing them through the Gaussian link-scale branch.
* `compare_maihda_groups(min_group_n = ...)` now guards the analytic sample size (the rows the model actually fits) rather than the raw group row count, so a group with enough raw rows but a tiny usable sample is skipped instead of being fit on a handful of observations.
* `predict_maihda(type = "strata")` now respects `newdata`: it returns only the strata present in `newdata` (and errors on a stratum the model never saw, as `type = "individual"` does) instead of always returning every training stratum. With `newdata = NULL` it still returns all strata.
* `compare_maihda_groups()` now counts populated strata for the pre-fit guard on the analytic model frame, not the raw subgroup. A group with two raw strata but only one stratum left after missing-row removal is now cleanly skipped as VPC-undefined instead of failing during fitting with "grouping factors must have > 1 sampled level".
* `n_boot` for bootstrap intervals must now be at least 10 (the minimum number of successful refits an interval requires); an unusably small `n_boot` fails immediately with a clear message instead of only erroring after the bootstrap runs.
* `calculate_pvc()` and `compare_maihda()` now detect differing prior `weights`: previously two models fit to the same rows/outcome/strata but with different weights compared as if equivalent (PCV returned silently, no warning), even though weights change the variance estimates. `calculate_pvc()` now errors and `compare_maihda()` warns; an unweighted fit and an explicit unit-weight fit are still treated as equal.
* When a two-level non-0/1 outcome is recoded to 0/1, the chosen mapping (which level becomes the modeled event = 1) is now reported via a `message()` and stored on the model as `$response_recoding`. Previously the mapping followed alphabetical (character) or declared (factor) level order silently, with no signal at all when `family = "binomial"` was passed explicitly, so the modeled event could be inverted unnoticed. The 0/1 assignment rule is unchanged (it matches base `glm`).
* Stratum-level plot aggregations (`plot(type = "predicted")`, `"risk_vs_effect"`, `"obs_vs_shrunken"`, `"effect_decomp"`) now honour prior `weights`, collapsing per-stratum predictions with a prior-weight-weighted mean (and weighting the reference lines by the summed weights). This makes the plots consistent with the weighted Gaussian VPC for survey/weighted fits; unweighted models are unaffected, and aggregated-binomial (`cbind`) fits, whose `weights(type = "prior")` are the trials, are left unweighted to avoid double-counting.
* The Shiny app's "Compute Bootstrap CIs" control now actually produces the VPC/ICC bootstrap intervals it (and the vignette) advertise. The bootstrap was previously applied only to the PCV, so the headline VPC/ICC was shown as a point estimate with no interval. The interval is now computed in the background worker (keeping the UI responsive) and displayed alongside the VPC/ICC in the Model Summary and Interactive Explorer tabs; the PCV interval in the PCV Results tab is unchanged.
* The Shiny app (`run_maihda_app()`) no longer aborts the whole analysis when the baseline (null) model has zero or negative between-stratum variance (a singular / no-between-variation fit), which makes `calculate_pvc()` error by design. The fitted model, VPC/ICC, summaries, and plots are now shown as usual and the PCV is reported as unavailable (with the underlying reason) instead of failing the entire workflow.

## Diagnostics

* `brms` fits now record MCMC convergence diagnostics (maximum Rhat and the number of divergent transitions) alongside the engine, surfaced in the "Fit diagnostics" block of `print()`/`summary()` so a non-converged or divergent Bayesian fit is no longer silent.
* Bootstrap VPC/PCV intervals now report Monte Carlo error: `summary()` and `calculate_pvc()` print the number of successful bootstrap draws and the Monte Carlo standard error so the precision of an interval can be judged.


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
