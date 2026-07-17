# Changelog

## MAIHDA (development version)

### API changes

- Gaussian PCV calculations now use each model’s fitted (REML)
  between-stratum variance by default. Use `estimation = "ML"` in
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md),
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md),
  or [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) to
  restore the previous ML-refit behaviour.
- `pcv_importance(method = "sequential")` is soft-deprecated. Use
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  for an order-dependent path or `method = "shapley"` for
  order-invariant attribution.

### Bug fixes

- Longitudinal PCV calculations now return `NA` when the null growth
  variance is effectively zero. The result records this in
  `null_at_boundary`.
- [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  now rejects duplicate stratum intercepts, including equivalent
  spellings such as `(1 | stratum) + (0 + 1 | stratum)` and compound
  terms that repeat the intercept such as
  `(1 | stratum) + (1 + x | stratum)`, while still allowing an
  uncorrelated slope such as `(1 | stratum) + (0 + x | stratum)`.
- Time-dependent formula offsets are now re-evaluated at each time point
  in longitudinal predictions and count VPC trajectories.
- Failed or skipped ML refits are reported as
  `estimation_used = "mixed"`; `brms` comparisons are reported as
  `"posterior"`.
- PCV boundary detection now covers `lme4`, `wemix`, and `ordinal` fits
  and is shared by
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md).
- Singular adjusted fits are marked in PCV results and print output
  without treating a near-100% PCV as automatically erroneous.
- ML-refit failures now warn, retain their actual estimation basis, and
  no longer produce information-criterion deltas across mixed REML/ML
  fits.
- External offsets are retained in fitted-row predictions and included
  in family detection, response recoding, automatic binning, and
  longitudinal time centering.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now computes shared numeric-strata cut points from the pooled analytic
  sample.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  and [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  now exclude rows with missing offsets when selecting the family and
  calculating analytic sample sizes.
- Failed per-group PCV decompositions now set `pcv_status = "failed"`
  and report the affected groups.
- Subsetting a `maihda_ic` object now preserves the metadata required by
  its print method.
- [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  now warns and omits deltas when likelihood and Bayesian information
  criteria are mixed.
- Ordinal models now recheck the number of observed response categories
  after analytic-sample filtering.
- Longitudinal validation and time centering now use the
  transformation-aware analytic model frame.
- PCV result objects now retain and print the requested and actual
  estimation basis.
- `fit_maihda(engine = "brms")` now rejects unsupported lme4-style
  `weights`, `subset`, and `offset` arguments.
- Design-weighted `brms` workflows now complete derived null and
  adjusted fits without conflicting with the internal weight column.
- Design-weighted `brms` fits now keep the full pre-fit data as
  `original_data`, so
  [`maihda_describe()`](https://hdbt.github.io/MAIHDA/reference/maihda_describe.md)
  reports the original total, missingness, and excluded-weight counts
  rather than the analytic sample alone.
- Non-positive or non-finite Gaussian precision weights are dropped
  before fitting.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  and [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  now apply the same non-positive/non-finite precision-weight exclusion
  before family/engine detection, shared-strata binning, and analytic
  sample-size checks.
- Individual predictions now reject missing strata unless
  `allow_new_levels = TRUE`.
- Count-family longitudinal VPCs now evaluate residual variance at each
  time point; the summary includes `var_resid_t`.
- Crossed-dimensions predictions now exclude contextual random effects
  from the stratum baseline.
- AUC calculations no longer treat lme4 precision weights as population
  frequencies.
- The AUC scope is now labelled `"intersectional"` when fixed effects
  contribute to the prediction.
- Crossed-dimensions response-scale VPCs now include additive dimension
  variances in the between-stratum component.
- Contextual binary models now report discriminatory accuracy and
  optional response-scale VPCs.
- `brms` convergence is reported as unknown when no MCMC diagnostic is
  available.
- Optional summaries and plots now warn when a requested component fails
  instead of silently omitting it.
- WAIC and PSIS-LOO reliability warnings are now passed through.
- Bootstrap requests below 200 replications now warn about unstable
  percentile intervals.
- Training predictions with an external offset now use the fitted
  model’s offset-aware prediction path.
- Crossed-dimensions models now reject additional random effects in the
  formula and direct users to `context` or the two-model decomposition.
- Sampling-weighted stepwise and importance analyses now exclude invalid
  weights before family detection.
- Automatic numeric-strata cut points are now based on the full analytic
  sample.
- PCV comparisons now require matching non-stratum random-effects
  structures.
- The `brms` count longitudinal VPC interval now propagates residual
  uncertainty by posterior draw.
- Contextual binary
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  results now include the AUC and MOR trajectory.
- Scalar PCV helpers now reject crossed-dimensions fits and direct users
  to the crossed-dimensions decomposition.
- Longitudinal ID/stratum checks now run on the analytic sample.
- The `brms` variance parser now supports grouping-variable names
  containing `__`.

### Testing and CI

- Added an `integration-tests` CI job for optional non-brms backends and
  expanded regression coverage for PCV, weighting, AUC, offsets,
  longitudinal models, and crossed-dimensions models.
- Added `test-audit-2026-07-13.R` for the latest regression fixes.

## MAIHDA 0.2.1

CRAN release: 2026-07-09

### New features

- Added
  [`maihda_describe()`](https://hdbt.github.io/MAIHDA/reference/maihda_describe.md)
  for pre-model sample, strata, missingness, weighting, and outcome
  summaries, with print and plot methods.
- [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now support stratified contextual analyses using `group` and `context`
  together.
- `plot(type = "predicted")` now orders strata by predicted value by
  default and accepts `order_by`.
- VPC plots for `maihda_analysis` objects can show the null model,
  adjusted model, or both.
- Added contextual models to
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md).
- Added
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  for Shapley and dominance-based PCV attribution, with optional
  bootstrap intervals.

### API changes

- Added a [`predict()`](https://rdrr.io/r/stats/predict.html) method for
  `maihda_model` objects.
- Renamed
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  to
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).
  The old name and `$pvc` result element remain as deprecated aliases.

### Bug fixes

- Preserved brms response addition terms such as `trials()` and
  [`weights()`](https://rdrr.io/r/stats/weights.html) when formulas are
  rebuilt.
- Normalized brms sampling weights after analytic-sample filtering.
- Rejected longitudinal IDs that are reused across strata.
- Corrected AUC handling for aggregated binomial data and fractional
  precision weights.
- Aligned AUC and MOR to the same intersectional model scope.
- Included non-stratum random effects in response-scale VPC
  calculations.
- Reported bootstrap boundary draws and rejected effectively zero PCV
  denominators.
- Fixed zero-row strata creation and aggregated-binomial formulas that
  use only the strata shorthand.
- Updated brms linear-predictor extraction for current brms releases.
- Fixed single-row ordinal predictions, longitudinal time centering,
  count residual variance, longitudinal PCV fitting, and higher-order
  slope PCV indexing.
- Removed whitespace padding from mixed-type stratum labels.

### Improvements

- Warned when numeric stratum dimensions enter adjusted models as linear
  effects.
- Improved guidance for unsupported bootstrap engines and documented
  latent-scale PCV rescaling.
- Corrected documentation for bootstrap Monte Carlo error, plotting,
  group comparisons, interaction estimates, and fitted-sample
  descriptions.

## MAIHDA 0.2.0

CRAN release: 2026-07-02

### New features

- Added the `"upset"` plot type and
  [`maihda_upset_size()`](https://hdbt.github.io/MAIHDA/reference/maihda_upset_size.md)
  for compact intersection displays.
- Set
  [`theme_maihda()`](https://hdbt.github.io/MAIHDA/reference/theme_maihda.md)
  as the default theme for package plots.
- Added `select = "deviation"` to retain the most extreme strata when
  plots are truncated.
- Added `only_flagged` to BLUP plots and ensured flagged strata are not
  hidden by display limits.
- Changed the default interaction adjustment to Benjamini-Hochberg and
  added ROPE-based decisions.
- Removed the redundant `"risk_vs_effect"` and `"ternary"` plot types
  and the optional `ggtern` dependency.

### Bug fixes

- Warned when model comparisons mix likelihood and Bayesian
  information-criterion scales.
- Rejected fixed interactions among stratum dimensions in standard
  MAIHDA decompositions.
- Restored discriminatory-accuracy summaries for brms
  aggregated-binomial fits.
- Limited ML-refit skipping to boundary stratum variances rather than
  any singular random effect.
- Weighted aggregated-binomial stratum predictions by trial count.
- Corrected brms ordinal stratum predictions to return expected category
  scores.
- Rejected scalar interaction diagnostics for longitudinal fits.
- Corrected Bayesian probability-of-direction values for negative
  effects.
- Standardized unseen-stratum prediction behaviour across engines and
  added `allow_new_levels`.
- Normalized brms aggregated-binomial predictions to per-trial
  probabilities.
- Returned `NA` rather than `NaN` for undefined crossed-dimensions
  shares.
- Included formula offsets in `wemix` and `ordinal` predictions.
- Selected families and engines from the analytic sample in high-level
  workflows.
- Reserved internal `.maihda_dim_` column names used for automatic
  binning.
- Added a note when
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)
  combines fitted REML variance rows with an ML-based PCV.

### Improvements

- Added colour-aware console output through `cli`, with automatic
  plain-text fallback.
- Added `var_between_adjusted_ml` to grouped comparisons.

## MAIHDA 0.1.11

CRAN release: 2026-06-18

### New features

- [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) now
  computes and reports intersectional interaction diagnostics by
  default.
- Added [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html) methods
  for models, summaries, and analyses.
- Added longitudinal growth-curve MAIHDA through the `id`, `time`, and
  `time_degree` arguments, including VPC and PCV trajectories.
- Added cumulative ordinal models through the `ordinal` and `brms`
  engines.
- Added negative-binomial models for overdispersed count outcomes.
- Added design-weighted models through `sampling_weights`, with `wemix`
  and brms support.
- Added contextual cross-classified MAIHDA through the `context`
  argument.
- Renamed the `"cross-classified"` decomposition to
  `"crossed-dimensions"`; the old name remains a deprecated alias.
- Added the high-level
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  workflow for null/adjusted model fitting and PCV reporting.
- Added the `maihda_country_data`, `maihda_sparse_data`, and
  `maihda_long_data` example datasets.
- Added
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)
  for model summaries and ranked-strata exports.
- Added automatic AUC and MOR summaries for binomial models and optional
  response-scale VPCs.
- Added
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  for AIC/BIC and WAIC/LOOIC comparisons.
- Added AUC and MOR trajectories to
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  for binary outcomes.
- Added
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  and plot highlighting for flagged strata.

### Improvements

- Restricted MOR to cumulative-logit and binary-logit models.
- Clarified the conditional interpretation of response-scale VPCs for
  adjusted models.
- Unified comparison and diagnostic plotting under the base
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic;
  older plotting helpers are deprecated.
- Added lme4 fit diagnostics to model output and grouped-comparison
  warnings.
- Added group plots for absolute between-stratum variance.
- Clarified PCV, VPC, weighting, supported-family, and plot terminology
  across the documentation and Shiny app.
- Removed checked-in rendered vignette HTML and corrected documentation
  dependencies and links.
- Removed package-install side effects from the health-data regeneration
  script.

### Performance

- Replaced row-by-row stratum matching with vectorized, collision-safe
  key matching.

### Bug fixes

- Applied analytic-sample identity checks to design-weighted fits.
- Evaluated longitudinal baseline PCV at the observed reference time.
- Rejected cross-sectional stratum summaries for longitudinal models.
- Reported analytic sample sizes for `wemix` information-criterion
  tables.
- Corrected Gaussian PCV comparisons by using ML refits for models with
  different fixed effects.
- Rejected VPC/ICC calculations for Gaussian models with non-identity
  links.
- Based binary-family detection and numeric-strata binning on the
  analytic sample.
- Fixed nested forwarding of data-masked `weights`, `subset`, and
  `offset` arguments.
- Preserved original response labels while evaluating subset
  expressions.
- Correctly sliced external weights, subsets, and offsets in grouped
  fits.
- Accounted for precision weights in Gaussian residual variance and
  stratum-level plot aggregation.
- Corrected effect-decomposition random effects and ordinal surprise
  scores.
- Aligned Shiny family detection and bootstrap VPC/ICC output with the
  core API.
- Added grouped-comparison checks for analytic sample size, populated
  strata, and shared stratum support.
- Corrected response-scale plotting for count models.
- Made `predict_maihda(type = "strata")` respect `newdata`.
- Required at least 10 bootstrap replications.
- Rejected PCV comparisons with different prior weights.
- Reported the event mapping when two-level outcomes are recoded to 0/1.
- Kept Shiny analyses available when the null between-stratum variance
  is zero.
- Accepted bare family functions in grouped comparisons.

### Diagnostics

- Added maximum R-hat and divergent-transition counts for brms fits.
- Added successful-draw counts and Monte Carlo standard errors to
  bootstrap output.

## MAIHDA 0.1.10

### New features

- Added
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  for VPC/ICC and variance comparisons across groups.
- Added
  [`plot_group_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_group_comparison.md)
  for forest and variance-composition plots.

### Bug fixes

- Excluded failed parametric-bootstrap refits from VPC and PCV
  intervals.
- Corrected the Poisson residual-variance approximation to
  `log(1 + 1 / mu)`.
- Omitted unsupported confidence intervals from prediction-deviation
  plots instead of drawing zero-width intervals.

## MAIHDA 0.1.9

### Bug fixes

- Clarified the interpretation of negative PCV values in the Shiny
  dashboard.
- Fixed the coverage workflow’s failure-artifact upload configuration.

## MAIHDA 0.1.8

CRAN release: 2026-05-16

### New features

- Added prediction-deviation, predicted-value, risk-versus-effect, and
  effect-decomposition plots.
- Added automatic tertile binning for numeric grouping variables with
  more than 10 unique values.
- Expanded the Shiny dashboard with the new plots and automatic-binning
  controls.
- Added automatic detection of binomial outcomes.

### Bug fixes

- Corrected latent residual variance for binomial and ordinal VPC/ICC
  calculations.

## MAIHDA 0.1.7

CRAN release: 2026-04-05

### New features

- Added
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md).
- Added the interactive Shiny dashboard.
- Improved bootstrap confidence-interval performance.
- Updated package tests, imports, dataset documentation, and the
  introductory vignette.

## MAIHDA 0.1.0

CRAN release: 2026-04-03

### Initial release

- Added intersectional-strata creation and multilevel model fitting with
  lme4 and brms.
- Added variance summaries, individual and stratum predictions, model
  comparison, and core plots.
- Added documentation, vignettes, and tests.
- Improved missing-value handling in
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).
