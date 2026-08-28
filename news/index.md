# Changelog

## MAIHDA (development version)

### New features

- [`summary()`](https://rdrr.io/r/base/summary.html) on a `maihda`
  analysis gained `which = "adjusted"`, matching
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html). It returns
  the null model by default, whose fixed effects are the intercept and
  covariates only because the strata dimensions are its random-effect
  grouping – easy to misread as the dimension main effects having gone
  missing. [`print()`](https://rdrr.io/r/base/print.html) now names
  which of the two models it is showing.
- [`summary()`](https://rdrr.io/r/base/summary.html) now reports
  fixed-effect standard errors, Wald statistics, two-sided p-values and
  intervals for the lme4, WeMix and ordinal engines, so
  `tidy(component = "fixed")` returns the full broom shape –
  `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high` – instead
  of estimates beside all-`NA` columns. The interval level follows
  `summary(conf_level = )`; no denominator degrees of freedom are
  estimated, so a Gaussian fit’s statistics are z-based rather than
  Satterthwaite/Kenward-Roger.
- A longitudinal [`summary()`](https://rdrr.io/r/base/summary.html) now
  also reports the intercept and slope VPCs of Bell et al. (2024, eq. 5)
  as `$longitudinal$vpc_intercept` and `$vpc_slope`. These exclude the
  occasion-level residual from the denominator, so they measure the
  between-stratum share of between-INDIVIDUAL variation in where a
  trajectory starts and how fast it changes, and are unaffected by
  measurement noise; the headline VPC keeps the residual and answers the
  discriminatory-accuracy question. `vpc_slope` is the only scale-free
  summary of how intersectional the rates of change are, and it pairs
  with the `PCV_slope` the decomposition already reported. Report the
  headline VPC unless you mean the trajectory question –
  [`?summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
  and the longitudinal vignette set out the contrast.
- [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  and [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  gained `stratum_slope` for longitudinal fits. `stratum_slope = FALSE`
  fits `(time | id) + (1 | stratum)`: individuals keep their growth
  curves but strata differ in level only, so the between-stratum
  variance is constant over time and the decomposition reports
  `PCV_intercept` alone. It is the remedy when the stratum slope
  variance sits at the singularity boundary, which few strata, few
  occasions per stratum, or irregular measurement times routinely
  produce. The VPC still varies with time, through the person-level
  slope variance and the residual.
- [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  models now carry likelihood-adequacy diagnostics that
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) surface alongside
  the singular/convergence caveats: Poisson/negative-binomial
  overdispersion (Pearson ratio) and zero inflation, stratum
  random-effect non-normality, longitudinal residual autocorrelation,
  and an approximate ordinal proportional-odds screen. A well-specified
  model stays silent.

### API changes

- The
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)
  intercept row now carries an interval (`*_lower`/`*_upper`, previously
  always `NA`): the summary’s Wald interval for the likelihood engines,
  the credible interval for brms. Being a between-stratum quantity it is
  too narrow when the strata are few, as the printed footnote now says.
  The variance and SD rows remain point estimates.
- The `fixed_effects` element of a
  [`summary()`](https://rdrr.io/r/base/summary.html) object gained
  `statistic`, `p_value`, `lower` and `upper` columns (and `se` for
  lme4), and `summary(conf_level = )` now sets the brms
  credible-interval quantiles as well as the Wald ones. Code reading
  `term`/`estimate`/`se` is unaffected; code assuming the exact column
  set is not.
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
- [`maihda_mor()`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md)
  on a crossed-dimensions fit now returns the mixture MOR over pairs of
  distinct strata instead of applying the independent-strata closed form
  to the summed variance. Reported values fall, most where the variance
  sits in the additive dimensions.
- `fit_maihda(engine = "wemix")` now rejects a `max_iteration` that is
  not a single whole number of at least 1.
- The Poisson and negative-binomial level-1 variance behind the VPC is
  now evaluated at a single mean count, `log(1 + 1/mean(lambda))`, as
  Stryhn et al. (2006) and Nakagawa, Johnson & Schielzeth (2017) define
  it, instead of averaging `log(1 + 1/lambda_i)` over rows. Null-model
  count VPCs are unchanged; adjusted, offset, weighted, and longitudinal
  count VPCs rise.

### Documentation

- [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  now states the predictive target of the Bayesian criteria: `brms`
  WAIC/LOOIC condition on the fitted random effects and so assess
  prediction of new observations within the represented strata, not
  generalisation to a new stratum (a leave-one-group-out
  cross-validation question), while the likelihood engines’ AIC/BIC are
  computed from the marginal likelihood.
- [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  now documents that its default BH flags are a conservative screen:
  partial pooling deflates a truly-null stratum’s Wald tail (null z
  variance is about the shrinkage fraction), so the flags under-flag
  rather than exceed the nominal false-discovery rate. The Bayesian path
  is described as already partially pooled instead of
  “multiplicity-free”, and the flag description drops causal phrasing.
- The finding-interactions vignette no longer calls the default BH flags
  “false-discovery-rate controlled”. It now matches
  [`?maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md):
  the flags are a conservative FDR screen that under-flags rather than
  over-flags, not an exact error-rate guarantee.
- The design-weighted documentation no longer describes the fixed-effect
  standard errors as design-consistent for complex surveys.
  `sampling_weights` takes one person-level weight column and cannot
  represent PSUs, sampling strata, higher-stage weights,
  finite-population corrections, or replicate weights, so it delivers
  population-weighted point estimates rather than general design-based
  inference.

### Performance

- `maihda(decomposition = "crossed-dimensions")` now fits the model once
  instead of twice. The preliminary pass that resolves the strata and
  family no longer refits the supplied formula only to discard it, which
  roughly halves the fitting cost (most noticeable for `brms`).

### Bug fixes

- A singular fit is no longer reported as a convergence failure. lme4
  files “boundary (singular) fit” among its own convergence messages, so
  a fit whose optimizer returned code 0 was recorded as
  `converged = FALSE` and printed twice – once as “Singular fit”, once
  under a “Convergence warnings reported by lme4” header – which also
  double-counted such groups in
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)’s
  singular and non-converged warnings.
- The singular-fit report now names the random-effects block that is at
  the boundary and scopes its VPC/PCV caveat accordingly. It previously
  said the between-stratum variance “may be unreliable” for every
  singular fit, which is false when the boundary sits in a non-stratum
  block – routine in a longitudinal fit whose `(time | id)` block has no
  person-level slope variation.
- `PCV_slope` is now `NA` when the null model’s between-stratum block is
  rank-deficient (a perfect intercept-slope correlation), with a note
  saying why. The stratum variation has collapsed onto a single
  direction in (intercept, slope) space, so the slope variance is not a
  free parameter and the null and adjusted fits need not have collapsed
  onto the same direction – the ratio does not compare the same quantity
  before and after adjustment. It was not merely uncertain but
  explosive: across replicates of an irregular design with no true
  stratum slope variance, rank-deficient fits returned values from -196
  to +1. The existing denominator guard could not catch it, because it
  scales a slope variance (outcome-units squared per time-unit squared)
  against the residual variance, so its relative threshold has no fixed
  meaning for that cell. `PCV_intercept` and `PCV(t)` are unaffected and
  still reported. [`print()`](https://rdrr.io/r/base/print.html) also
  gained the adjusted-model boundary note the cross-sectional
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  has carried all along.
- [`print()`](https://rdrr.io/r/base/print.html) on a longitudinal
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) now
  shows both growth fits’ diagnostics, labelled. Only the null model’s
  were shown, so a singular or non-converged adjusted fit – which pins
  the reported additive share – delivered its headline PCV in silence.
- [`print()`](https://rdrr.io/r/base/print.html) on a cross-sectional
  two-model
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) now
  shows fit diagnostics, which that branch printed for neither model: a
  non-converged fit, or one singular in a non-stratum block such as a
  `context` random intercept, produced its headline PCV in silence. A
  stratum-block singularity in the ADJUSTED fit is deliberately not
  banner-warned – it is the expected shape of an additive decomposition
  – and is reported instead by the boundary note below, so the
  diagnostics block stays a signal rather than firing on nearly every
  healthy analysis.
- The “PCV is pinned near 100%” boundary note now reaches
  [`print()`](https://rdrr.io/r/base/print.html) on a
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis.
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  has recorded `adjusted_at_boundary` all along and
  [`print.pcv_result()`](https://hdbt.github.io/MAIHDA/reference/print.pcv_result.md)
  showed the caveat, but the analysis print formats its PCV inline
  rather than calling `print(x$pcv)`, so a “100% additive” headline read
  off a singular adjusted fit carried no caveat anywhere the reader
  looks. Both methods now share one note.
- The per-group variance extractors behind the contextual and
  crossed-dimensions partitions now enforce the intercept-only contract
  the ordinary summary already applied. A contextual lme4 fit with
  `(1 + x | stratum)` previously returned a partition whose total
  silently dropped the slope variance and intercept-slope covariance,
  and a slope-only context `(0 + x | site)` produced an all-`NA`
  partition; the brms path accepted a slope-only `sd_<group>__<slope>`
  column, squaring slope draws as the group’s intercept variance.
- `maihda(decomposition = "crossed-dimensions")` now rejects a random
  slope on an allowed grouping factor, e.g. `(1 + x | stratum)`, instead
  of silently rewriting it to the canonical intercept-only crossed
  formula. The builder previously checked only the grouping side of each
  random-effect bar, so the slope vanished without a warning through
  both [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  and the per-group crossed workflow.
- `predict_maihda(type = "strata")` now rejects a non-longitudinal fit
  whose stratum random effects include random slopes, matching the
  summary guard. It previously returned the intercept column alone as a
  complete-looking table (estimate, SE, interval) – the stratum effect
  at zero of the slope variables, an extrapolation that can reorder or
  sign-flip strata. The brms extractor also no longer presents a
  slope-only `(0 + x | stratum)` block’s slope as the stratum effect.
  Longitudinal trajectory predictions and summaries are unchanged.
- [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  now rejects a longitudinal growth curve the observed times cannot
  identify: fewer than `time_degree + 1` distinct measurement times, or
  no person measured at two distinct times. Such models previously
  fitted and returned arbitrary, optimizer-dependent slope variances
  flagged only as a singular fit. The check runs on the input and again
  on the analytic sample after row exclusions.
- [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  now reports `WAIC`/`LOOIC` as `NA` for a sampling-weighted `brms` fit,
  with estimator `"Bayesian (weighted pseudo-posterior)"`, matching the
  `wemix` treatment: the weighted pointwise log-likelihoods of a
  pseudo-posterior are not log predictive densities and define no
  standard criterion.
- [`summary()`](https://rdrr.io/r/base/summary.html), the VPC, and the
  PCV helpers now support a binomial model fitted with the complementary
  log-log (`cloglog`) link, using the extreme-value latent level-1
  variance `pi^2/6`. Such a model previously fitted but then stopped
  with a “not implemented” error, even though the response-scale VPC
  already handled every binomial link.
- [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  and
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  now emit the same warning
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) does
  when a numeric stratum dimension (a category code) enters the models
  as a single linear term rather than categorical main effects. The
  fitted models and the reported estimates are unchanged; only the
  previously missing warning was added.
- `pcv_importance(bootstrap = TRUE)` now warns and records
  `n_boot_boundary` when bootstrap draws are excluded because the null
  model’s between-stratum variance hit the zero boundary, so an
  attribution interval conditional on only a handful of surviving draws
  is no longer returned silently.
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  already disclosed this.
- The crossed-dimensions, contextual, and longitudinal VPC bootstraps
  now count and report non-converged refits (`n_boot_nonconverged` and a
  warning), so their `n_boot_ok` no longer implies a convergence that
  was never checked. The longitudinal per-time bands also require a
  majority of finite draws to form, not only ten.
- `fit_maihda(family = "binomial")` no longer recodes a two-level
  proportion response (e.g. 0.25 and 0.5) to 0/1. A numeric response is
  Bernoulli only when its values are whole numbers, so a proportion of
  successes supplied with `weights =` trial counts stays an aggregated
  binomial fit instead of a silently different one.
- Bootstrap intervals (VPC, PCV, crossed and contextual decomposition,
  longitudinal VPC,
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md))
  now require a majority of the eligible refits to succeed, not only an
  absolute minimum of ten. An extreme failure rate makes the interval
  unavailable rather than returning one from a handful of draws;
  legitimately excluded draws (e.g. PCV boundary draws) are not counted
  against the success fraction.
- The crossed-dimensions MOR now returns 1, not `NA`, when at least half
  of the stratum pairs have zero contrast variance. The median odds
  ratio is exactly 1 there, but the mixture root search could not
  bracket it.
- Bootstrap draws whose refit optimizer did not converge are counted and
  reported (`n_boot_nonconverged` and a warning) instead of silently
  counted as successful refits, so `n_boot_ok` no longer implies
  convergence that was never checked.
- Integer `weights=` on a Bernoulli fit are no longer read as aggregated
  binomial trial counts by
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md).
  Aggregation is now inferred structurally from a
  [`cbind()`](https://rdrr.io/r/base/cbind.html) matrix response or a
  proportion response, so precision weights no longer change the AUC
  estimand or inflate the reported case/control totals.
- `wemix` fits are no longer reported as converged unconditionally.
  WeMix returns the last iterate without warning when its optimisation
  is abandoned, so convergence is now judged from its own gradient
  criterion and reported as `NA` when no evidence is readable.
- Longitudinal count VPC trajectories now evaluate raw-time fixed-effect
  terms such as `x:wave` at the reporting time. Under internal time
  centering only the derived centered column was moved, so those terms
  stayed at each row’s own observed time.
- [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  now warns and omits the delta when the models differ in prior or
  sampling weights, which change the likelihood being maximised.
- [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  now reports `Step_PCV` as `NA` when the preceding step’s
  between-stratum variance is at the singularity boundary, instead of
  dividing by optimizer noise. The affected steps are listed in an
  `"undefined_step_pcv"` attribute and noted by
  [`print()`](https://rdrr.io/r/base/print.html); `Total_PCV` is
  unchanged.
- The response-scale VPC for a model with non-stratum random effects now
  estimates the total probability variance on the same basis as its
  numerator, removing an `n_sim / (n_sim - 1)` inflation that was
  largest at small `n_sim`.
- `calculate_pcv(bootstrap = TRUE)` now permutes each simulated response
  into each model’s own row order before refitting, so a bootstrap
  across two fits of the same observations in a different row order no
  longer corrupts the interval.
- Individual predictions now reject a supplied `stratum` that
  contradicts the intersectional dimension columns in the same
  `newdata`, instead of pairing one intersection’s fixed effects with
  another’s random effect. The check also covers dimension combinations
  the model never saw, and a row whose dimensions cannot be resolved no
  longer exempts the rest of the `newdata`.
- Predictions now check the numeric auto-bin ranges whether or not
  `newdata` supplies a `stratum` column. A row whose auto-binned numeric
  dimension falls outside the training range previously went unremarked
  when a `stratum` was supplied, silently pairing that stratum’s random
  effect with a combination the training bins cannot contain; it now
  warns. Existing results are unchanged, and the warning will become an
  error in a future release.
- `lme4` convergence reporting now also consults the optimizer return
  code and message, so a fit whose optimizer stopped early
  (e.g. `bobyqa` hitting `maxfun`) is no longer reported as converged
  when `lme4`’s own gradient check happens to pass.
- `brms` individual predictions with `allow_new_levels = TRUE` now zero
  every unseen grouping level (context and longitudinal, not only
  stratum) to match `lme4`, instead of sampling unseen non-stratum
  effects from the random-effects distribution.
- Zeroing an unseen grouping level no longer overrides a caller-supplied
  `re_formula` or `re.form`. The requested scope is narrowed to drop the
  unseen terms rather than replaced, so `re_formula = NA` stays a
  fixed-effects-only prediction and a partial `re_formula` is not
  swapped for a different term.
- Longitudinal PCV now records `estimation_used` (`"fitted"`, `"ML"`, or
  `"mixed"`), so a boundary skip or failed ML refit that leaves a mixed
  REML/ML comparison is reported rather than mistaken for a clean fitted
  one.
- [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md),
  and
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  now treat two fits to the same observations supplied in a different
  row order as the same analytic sample, instead of rejecting the
  reordered fit or warning about a differing sample. The comparison
  aligns on the row identifiers and matches the response, stratum
  partition, and weights after alignment.
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
- `brms` convergence is reported as unknown unless an R-hat is
  available; a divergent-transition count alone is not evidence that the
  chains converged.
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
