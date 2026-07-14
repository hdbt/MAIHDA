# Changelog

## MAIHDA (development version)

### API Changes

- **The Gaussian PCV now defaults to each model’s own fitted (REML)
  between-stratum variance, via a new `estimation` argument.**
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md),
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md),
  and [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  gain `estimation = c("fitted", "ML")`. Previously the PCV *always*
  refitted a Gaussian `lmer` model with maximum likelihood before
  differencing the between-stratum variances. That is defensible — it
  removes REML’s fixed-effects-specific degrees-of-freedom correction,
  matching
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)/[`anova()`](https://rdrr.io/r/stats/anova.html)
  — but ML variance components are downward-biased in finite samples,
  most sharply with few strata (the usual MAIHDA regime) and more so for
  the adjusted model, which *inflates* the reported PCV. The new default
  `"fitted"` differences each model’s own REML variance (the values
  [`summary()`](https://rdrr.io/r/base/summary.html) reports, and
  conventional MAIHDA practice), avoiding that bias; `estimation = "ML"`
  restores the previous behaviour for a correction-free cross-model
  comparison. **This changes the Gaussian-`lmer` PCV reported by
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)/[`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  relative to 0.2.1** — in a 30-stratum example the PCV moved from
  0.94 (ML) to 0.85 (REML). Binomial/Poisson/ordinal and brms fits are
  already on the ML scale and are unaffected. The result object records
  the basis in `$estimation`, and
  [`print()`](https://rdrr.io/r/base/print.html) states it.

- **`pcv_importance(method = "sequential")` is soft-deprecated in favour
  of
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md).**
  The `"sequential"` method now emits a deprecation warning (it still
  runs and returns the same result) and will be removed in a future
  release; `"shapley"` and `"dominance"` are unaffected. It was a
  strictly weaker duplicate of
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md):
  the sequential *trajectory* quantities — the step-specific `Step_PCV`
  (normalised by the previous step, so it does not obey the efficiency
  identity
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  is built on) and, for a binary outcome, the discriminatory-accuracy
  path (`AUC`, `Step_AUC`/`Total_AUC`, `MOR`) — live naturally in
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)’s
  row-per-step table and have no place in the row-per-variable
  attribution object. Use
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  for the order-dependent path, or `method = "shapley"` for an
  order-invariant split of the same total PCV.

### Bug Fixes

- **An external `offset=` is now retained in every downstream `lme4`
  prediction, not just at fit time.** Although fitting honoured the
  offset, several prediction paths rebuilt `newdata` from the stored
  model frame and lost it: `predict.merMod()` silently drops an
  *external* offset (the `offset=` argument) on the `newdata` path, and
  *errors* on a *formula*
  [`offset()`](https://rdrr.io/r/stats/offset.html) term because the
  model frame stores the offset under the derived
  `offset(...)`/`(offset)` column name rather than the raw variable.
  This affected the ranked-strata table and `plot(type = "predicted")`
  (per-stratum expected counts came out too low for an external offset,
  and errored for a formula offset), the count longitudinal time-varying
  VPC (`VPC(t)` biased by the dropped offset), the
  mean-covariate-profile growth trajectory, and the prediction-deviation
  panels. Predictions on the fitted rows now reuse the stored linear
  predictor (which includes the offset), and predictions on derived
  grids rebuild the fixed-effects linear predictor directly and add the
  model’s offset (held at its per-row value, or its mean for the
  representative-profile trajectory); a genuine external prediction
  grid, where an external offset cannot be reconstructed, is rejected
  with a directed error. No-offset fits are numerically unchanged.

- **[`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now auto-bins shared numeric strata on the pooled analytic sample, so
  `subset=` matches prefiltered input.** With `shared_strata = TRUE` the
  shared tertile cut-points were computed from *all* input rows rather
  than the rows the per-group fits actually use, so a row excluded by
  `subset=`, a missing/zero weight, a missing offset, missingness, or a
  missing group still shifted the breaks — silently redefining other
  rows’ strata and changing the reported VPCs.
  `compare_maihda_groups(data, subset = keep)` now returns the same
  strata (and VPCs) as fitting the pre-filtered data, mirroring
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)’s
  existing analytic-sample binning.

- **[`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  and [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  now fold the external offset into family detection and the per-group
  analytic n.** As in
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
  an offset-`NA` row (which `lme4` drops) was still seen by the
  group/workflow family-and-engine selection and counted in the reported
  per-group `n`/`min_group_n` guard. A stray out-of-sample outcome value
  on such a row could flip the detected family — e.g. an ordered outcome
  whose third category sits only on the offset-`NA` row selected the
  ordinal engine, then failed against
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)’s
  binary analytic sample. Detection and the analytic row count now key
  off the same offset-aware sample the engine fits.

- **A failed
  [`lme4::refitML()`](https://rdrr.io/pkg/lme4/man/refitML.html) is no
  longer silently reported as an ML result.** When `estimation = "ML"`
  was requested but `refitML()` failed,
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  kept the original REML fit yet still labelled the result
  `estimation = "ML"`;
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  could likewise compute an AIC/BIC delta across a mixed REML/ML basis.
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  now warns on an ML-refit failure and records `ml_refit_failed = TRUE`,
  and
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  withholds the delta (with a warning) when an intended ML refit fell
  back to REML, leaving the criteria on incomparable bases.

- **A singular adjusted fit that saturates the PCV near 100% is flagged
  under any `estimation` basis.**
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)’s
  adjusted-model boundary warning was gated on `estimation = "ML"`, so
  under the default `estimation = "fitted"` a singular adjusted model
  produced an unqualified `PCV = 1` with no caveat. The warning now
  fires whenever model2 is on the singularity boundary regardless of the
  basis, and the result records `adjusted_at_boundary`. The status is
  propagated as silent attributes to
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  (`boundary_steps`) and
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  (`full_at_boundary`) — surfaced in their
  [`print()`](https://rdrr.io/r/base/print.html) methods rather than as
  a warning, since a saturated additive full model is the common,
  legitimate case there.

- **Subsetting a `maihda_ic` result no longer breaks its
  [`print()`](https://rdrr.io/r/base/print.html) method.** Column/row
  subsetting kept the `maihda_ic` class but dropped the `ic_primary`
  attribute, so [`print()`](https://rdrr.io/r/base/print.html) printed
  the table and then errored on a zero-length condition
  (`!is.na(NULL)`). A `[.maihda_ic` method now preserves the metadata
  across subsetting, and
  [`print.maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/print.maihda_ic.md)
  guards the attribute defensively.

- **An external `offset=` no longer skews automatic family detection,
  response recoding, strata auto-binning, or longitudinal time-centering
  through a row `lme4` later drops.** The evaluated top-level `offset`
  was forwarded to `lme4` (which puts it in its model frame, so
  `na.omit` removes offset-`NA` rows) but was omitted from the
  analytic-sample masks the preprocessing keys off. A row later dropped
  for a missing offset was therefore still seen upstream: it could flip
  the auto-detected family (a stray out-of-sample outcome value making
  an otherwise-binary outcome look continuous, so a linear model was
  fitted where a `binomial` one was intended), shift the response
  recoding / strata cut-points, or — in a longitudinal fit — drag the
  internal time centre to an occasion the fit never uses (biasing the
  growth basis and the time-varying VPC/PCV). The offset is now folded
  into the shared row mask (`maihda_row_mask()`) alongside
  `subset`/`weights`, so family detection, 0/1 recoding, auto-binning,
  and longitudinal validation/centering all see exactly the rows the
  engine fits. Only the `lme4` engine takes an offset (the others reject
  it), and a fit with no missing offset is unaffected.

- **A per-group PCV decomposition that fails is no longer silently
  reported as a success.** In
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  (two-model mode), the adjusted-model fit and
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  were wrapped in a
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html) that turned any
  error into `NULL`, leaving `pcv = NA` while the group’s `status`
  stayed `"ok"` and no warning explained the gap — an incomplete
  decomposition that looked complete (a common trigger: a stratum
  dimension constant within a group, whose adjusted main effect is a
  one-level factor `lmer` rejects). The result now carries a
  `pcv_status` column recording the decomposition outcome per group
  (`"ok"`, `"failed"`, or `"singular"`), the captured error is surfaced
  in an aggregated warning naming the affected groups, and an adjusted
  fit that is singular — pinning the adjusted variance near 0 so the PCV
  saturates near 100% — is flagged separately as a boundary artifact.
  The group’s own `status` still reflects its (successful) null VPC
  model.

- **[`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  now warns and withholds the delta when a comparison mixes likelihood
  and Bayesian information criteria.** A direct
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  call spanning the likelihood engines (`lme4`/`ordinal`, reporting
  AIC/BIC) and `brms` (reporting WAIC/LOOIC) checked only for differing
  outcome, family, or sample — not for a mixed criterion scale — so a
  same-family `lme4`-vs-`brms` table selected AIC as the primary
  criterion, ranked on it, and left the Bayesian row a bare `NA` delta
  with no caveat. It now applies the same scale guard
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  already uses: when the populated criteria span both scales it warns
  that AIC/BIC and WAIC/LOOIC are not comparable and omits the delta
  (the per-model criteria are still reported). Same-scale comparisons
  are unaffected.

- **An explicit ordinal (`family = "ordinal"`) model whose top category
  occurs only on rows with a missing predictor no longer silently
  collapses to a binary fit.** The response category count was validated
  on the full data (`nlevels >= 3`) *before* the analytic-sample drop,
  so a category confined to rows the fit later removes (a missing
  predictor or outcome) passed the check — yet
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html), and
  the brms `cumulative()` path, then fit only the observed categories
  without complaint, producing a one-threshold (binary) model while the
  wrapper still recorded three levels and response predictions used the
  wrong 1..K scale. The count is now re-checked on the analytic sample
  (after complete-case filtering, with any now-empty category dropped)
  for both the `ordinal` and `brms` engines, erroring with a clear
  message rather than degrading the model order — the brms check runs
  *before* the Stan fit. A clean fit with every category observed is
  unaffected.

- **A longitudinal MAIHDA with a non-finite transformed predictor no
  longer mis-anchors the internal time centre or raises a false
  cross-stratum id error.** The analytic-sample mask that the
  repeated-measures check, the globally-unique-id check, and the
  internal time-centering all key off was built with
  [`complete.cases()`](https://rdrr.io/r/stats/complete.cases.html) on
  the *raw* model columns, so a row whose fixed-effect transformation is
  non-finite (e.g. `log(x)` of `x <= 0`) — which `lmer` drops — was
  retained. That excluded row could drag the time centre to an
  out-of-sample occasion (biasing the growth basis and the time-varying
  VPC/PCV) and, if its stratum-defining values differed, trigger a
  spurious “id … appear in more than one stratum” error. The mask now
  uses the transformation-aware analytic model frame (as the
  `wemix`/`ordinal` engines already did), intersected with an explicit
  complete-case check on the id/time columns, with the previous raw
  check retained only as a fallback.

- **[`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now retain and print the chosen variance-estimation basis, and the PCV
  documentation no longer claims the basis is always ML.** These three
  gained `estimation = c("fitted", "ML")` alongside
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  but, unlike it, did not record the choice on their result objects, so
  a serialized result lost that (statistically material) provenance.
  Each now stores the basis —
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  in `$estimation`,
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  in an `"estimation"` attribute that survives subsetting — and states
  it in [`print()`](https://rdrr.io/r/base/print.html). Documentation
  passages in
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  that still said the PCV is *always* refitted with maximum likelihood —
  stale since `"fitted"` became the default — now describe the
  `estimation`-dependent behaviour.

- **`fit_maihda(engine = "brms")` now rejects the lme4-style
  `weights`/`subset`/`offset` arguments instead of silently ignoring
  them.** These are not
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  arguments — brms absorbs them into `...` and drops them — whereas
  family auto-detection and strata auto-binning upstream *did* honour
  them, so preprocessing described one analytic sample while brms fit
  another (silently altering coefficients, variance components, VPC, and
  PCV). They are now rejected with a message pointing to formula
  addition terms ([`weights()`](https://rdrr.io/r/stats/weights.html),
  [`offset()`](https://rdrr.io/r/stats/offset.html)), design weights via
  `sampling_weights`, or prefiltering `data` — mirroring the existing
  `wemix`/`ordinal` guards. Only the `lme4` engine, which applies them
  directly, still accepts them.

- **Design-weighted `brms` MAIHDA now completes its second
  (null/adjusted) fit.** A sampling-weighted brms fit stores a formula
  carrying the internal `weights(.maihda_sw)` addition term and a
  `.maihda_sw` data column.
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  (two-model and crossed-dimensions) and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  derive the null/adjusted models from that stored formula and refit
  with `sampling_weights=`, so the reserved-column guard saw
  `.maihda_sw` in both the formula and the data and aborted — failing
  weighted `maihda(..., engine = "brms")` outright, and letting
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  silently record `NA` PCVs. The weight-preparation step now strips the
  internal term (and drops the stale column) before re-normalizing from
  the original weight column, so the derived fits succeed. Single-model
  weighted
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  was unaffected.

- **Zero (or negative / non-finite) Gaussian precision `weights=` are
  now dropped before the `lmer` fit instead of producing a degenerate
  fit with a finite VPC.** lme4 keeps a zero-weight row and returns a
  degenerate fit (`logLik = -Inf`, `NA` gradient), after which the
  residual-variance helper silently discarded the zero weight and still
  reported a finite VPC — with variance estimates materially different
  from fitting after removing the row. Such rows are now excluded from
  the analytic sample up front (with a warning), so binary/ordinal
  detection, strata auto-binning, and the fit all use the same rows and
  the fit matches fitting the row-removed data exactly. (Only `lme4`
  takes precision weights; the other engines already reject them.)

- **[`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  no longer returns a silent population-average for a row whose stratum
  is missing.** For `type = "individual"` predictions with the default
  `allow_new_levels = FALSE`, a row with an `NA` stratum — supplied
  directly, or produced by a missing stratum-defining dimension —
  slipped through the unseen-stratum gate (which dropped `NA`s before
  checking), and the `wemix`/`ordinal` engines then mapped its absent
  random effect to zero, yielding a fixed-effects-only prediction where
  `lme4` correctly errors. The gate now rejects a missing stratum for
  individual predictions across all engines (unless
  `allow_new_levels = TRUE` opts into the population-average fallback);
  `type = "strata"`, where an `NA` stratum simply yields no row, is
  unchanged.

- **The
  [`summary.maihda_model()`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
  documentation now describes the count-family VPC correctly.** The note
  stated the Poisson / negative-binomial residual variance is built from
  posterior-*mean* fixed effects and random-effect variances “rather
  than per draw”; the implementation in fact propagates the marginal
  expected counts *per draw* for the intercept-only VPC structures
  (strata, crossed-dimensions, contextual), holding the posterior-mean
  plug-in only for the random-slope / longitudinal fallback (the
  negative-binomial `shape` draws are always propagated). Documentation
  only; results are unchanged.

- **The count longitudinal VPC trajectory now evaluates the level-1
  (residual) variance at each reporting time instead of reusing one
  sample-wide average.** For a Poisson or negative-binomial longitudinal
  (growth) MAIHDA, the latent-scale level-1 variance
  `log(1 + 1/λ[ + 1/θ])` depends on the marginal expected count `λ`,
  which changes over time in a growth model. The time-varying VPC
  summary computed a single average of that residual across all
  observations and reused it at every point on the reporting grid,
  biasing the whole `VPC(t)` trajectory and its bootstrap/posterior
  bands. The residual is now recomputed at each grid time — setting the
  model’s time term to that value and marginalizing the per-row level-1
  term over the observed covariate rows, with the log-normal mean
  correction from the growth-block variance at that time — on both the
  lme4 (point estimate + parametric bootstrap) and brms (posterior)
  engines. Gaussian and binomial longitudinal fits are unaffected (their
  level-1 variance is constant in time). The longitudinal summary object
  gains a `var_resid_t` vector (the residual over the grid); its scalar
  `var_resid` is now the reference-time value rather than the
  whole-sample average.

- **Crossed-dimensions stratum predictions no longer fold in a
  contextual random effect.** When a
  `decomposition = "crossed-dimensions"` fit also carried a contextual
  random intercept (`context =`), the per-stratum predicted-outcome
  baseline was built as `predict(<all random effects>) − u_stratum`,
  which left the contextual/site effect in the baseline. The ranked
  stratum predictions — `maihda_table()$strata` and
  `plot(type = "predicted")` — then depended on each stratum’s observed
  context composition rather than the intended intersectional scope. The
  baseline is now scoped (via `re.form`/`re_formula`) to the fixed
  effects plus the additive dimension random effects only, excluding any
  non-intersectional random effect; the interaction stratum RE is added
  back per row. Both the lme4 and brms engines. Canonical
  single-`(1 | stratum)` models and plain contextual fits (no
  crossed-dimensions `cc_info`) are unaffected.

- **Discriminatory-accuracy AUC no longer misreads lme4 precision
  weights as population frequencies.** For a Bernoulli `lme4` fit with
  non-integer `weights=` (precision weights, which scale the
  observation’s likelihood/dispersion, not its population frequency),
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
  previously folded the weights into case/control mass and reported a
  weighted Mann–Whitney concordance — a quantity with no population-AUC
  interpretation, and one that silently changed the estimand based on a
  fitting control. Such fits now report the ordinary observation-level
  AUC (`weighted = FALSE`), with a new
  `precision_weights_ignored = TRUE` flag and a
  [`print()`](https://rdrr.io/r/base/print.html) note. The
  sampling-weighted (design-based) and aggregated trial-count AUC paths
  — where case/control mass genuinely represents population/replication
  mass — are unchanged.

- **The intersectional-scope AUC is no longer mislabelled as strata-only
  discrimination.** When a model carries non-stratum random effects (a
  contextual `(1 | school)` or an explicit `(1 | site)`), the headline
  AUC excludes those but retains the *entire* fixed-effects predictor —
  so an adjusted model’s individual-level covariates (e.g. `age`) enter
  it. `auc_scope` is renamed from `"strata"` to `"intersectional"`, and
  [`print()`](https://rdrr.io/r/base/print.html)/the documentation now
  state that this is an adjusted intersectional concordance (which
  matches the between-stratum MOR’s scope only when the fixed part is
  intercept-only), not the discriminatory accuracy of the strata alone —
  for which you score the null (strata-only) model.

- **[`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  and
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  now reject a boundary-level (effectively singular) null denominator,
  as
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  already did.** Both tested only whether the null between-stratum
  variance was strictly positive, so a degenerate fit that lands just
  off zero (e.g. `1.6e-17` rather than exactly `0`) slipped through and
  every PCV was divided by it — turning a null with no usable
  between-stratum variation into a spurious 100% or a huge negative PCV
  (a reported Total PCV around `-2.8e14`) with no warning. Both now
  apply the same lme4 relative-singularity guard
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  uses, and
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  also drops any bootstrap draw whose null refit lands on the boundary;
  they error with a clear message rather than tabulating a degenerate
  ratio. The Shiny dashboard degrades a degenerate stepwise
  decomposition to an empty panel rather than aborting the fit.

- **The response-scale VPC now counts the additive dimension variances
  as between-stratum for a crossed-dimensions fit.**
  [`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  treated only the pure intersection (`stratum`) component as
  between-stratum variance and folded the additive dimension random
  effects into the non-stratum (“other”) variance it integrates over —
  so for a `decomposition = "crossed-dimensions"` model it understated
  the intersectional response-scale VPC by orders of magnitude (in one
  fit ~0.005 against the correct ~0.22). For a crossed-dimensions fit
  the between-stratum variance is now the sum of the additive dimension
  effects and the interaction — the total intersectional variance,
  exactly as the MOR (`maihda_mor_between_variance()`) and the
  latent-scale crossed-dimensions VPC read it — and only genuinely
  contextual/non-intersectional effects enter the integrated
  `var_other`. Canonical single-stratum and contextual fits are
  unchanged.

- **[`summary()`](https://rdrr.io/r/base/summary.html) no longer
  suppresses the discriminatory accuracy and response-scale VPC of a
  contextual binary fit.** The gate that attaches these binomial
  companions skipped every model carrying a `context =` random
  intercept, on the (now stale) grounds that their estimand would not
  match the stratum-vs-context partition. It does now: the headline AUC
  is the intersectional-scope concordance that *excludes* the context
  random effect (with the all-effects value reported separately as
  `auc_full`), and
  [`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  integrates the context variance into its denominator. A contextual
  binomial `lme4` fit therefore now reports its discriminatory accuracy,
  and `summary(model, response_vpc = TRUE)` its response-scale VPC,
  instead of silently returning `NULL`. Crossed-dimensions and
  longitudinal fits still skip them.

- **A brms fit is no longer recorded as converged when no MCMC
  diagnostic is available.** `maihda_fit_diagnostics()` derived a brms
  fit’s convergence solely from whether it had captured a warning
  message, so a fit for which *neither* the maximum R-hat *nor* the
  divergent-transition count could be computed – a single-chain fit, a
  variational/approximate algorithm, an unusual object – produced an
  empty message set and was reported as `converged = TRUE` on no
  evidence. Convergence is now `NA` unless at least one of those
  diagnostics was actually available.

- **Requested summary and plot outputs that fail now say so instead of
  vanishing.** The discriminatory accuracy, the response-scale VPC
  (`summary(response_vpc = TRUE)`), the automatic interaction
  diagnostics, and the individual panels of `plot(type = "all")` were
  each wrapped in a bare `tryCatch(., error = function(e) NULL)`, so a
  genuine error in one silently produced a partial object with no sign a
  calculation had failed. Each now still degrades to `NULL` (the core
  result never breaks) but re-emits the original error as a warning, via
  a shared `maihda_try_optional()` helper, so an explicitly requested
  output is never dropped without a trace. The analysis-level
  `plot(maihda_analysis, type = "all")` montage — a separate code path
  from the model-level plot — is now covered too, so its optional
  null/adjusted/group panels likewise warn on failure rather than
  silently vanishing (the deliberately mutually-exclusive
  `group_pcv`/`group_additive_share` pair, where exactly one applies by
  decomposition mode, stays quiet).

- **WAIC and PSIS-LOO reliability warnings are no longer suppressed.**
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  wrapped [`brms::waic()`](https://mc-stan.org/loo/reference/waic.html)
  / [`brms::loo()`](https://mc-stan.org/loo/reference/loo.html) in
  `suppressWarnings(suppressMessages())` and kept only the scalar
  criterion, hiding high-Pareto-k and low-ESS warnings while still
  letting models be ranked on those criteria. It now keeps brms’s
  progress *messages* quiet but captures its reliability *warnings* and
  re-emits them, so an information-criterion comparison never reports a
  criterion as if reliable when its own diagnostics say otherwise.

- **A low bootstrap-replication count now warns.** Requesting `n_boot`
  below ~200 (the hard floor stays 10, so the fast internal tests are
  unaffected; the default remains 1000) now warns that a percentile
  interval’s tail endpoints are order statistics from few draws and are
  unstable. Applies wherever an interval is bootstrapped – VPC, PCV,
  longitudinal, contextual, cross-classified, and the Shapley
  attribution.

- **`calculate_pcv(estimation = "ML")` warns when the ML refit makes the
  adjusted model singular.** An ML refit can push a small-but-positive
  REML between-stratum variance onto the zero boundary; the PCV then
  collapses toward 1 (“covariates explain all between-stratum variance”)
  as a boundary artefact rather than a substantive result. This case now
  warns and points back to the default `estimation = "fitted"` (REML).
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  inherit the guard through
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).

- **[`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  no longer drops an external offset from training predictions.** A
  model fit with `fit_maihda(..., offset = log(exposure))` stores the
  offset outside the formula. The wrapper always replaced a `NULL`
  `newdata` with the stored model frame and called
  `predict(newdata = )`, a path `predict.merMod` evaluates *without* the
  external offset — so individual predictions, and the discriminatory
  accuracy, tables, and plots built on them, were computed as if the
  offset were zero (errors tracking exposure exactly). Training
  predictions (no `newdata`) now use the engine’s own no-`newdata` path,
  which retains the offset; a formula
  [`offset()`](https://rdrr.io/r/stats/offset.html) term still works on
  `newdata`; and an external offset with `newdata` — which cannot be
  reconstructed from the stored fit — now errors with a pointer to
  writing the offset into the formula.

- **`decomposition = "crossed-dimensions"` no longer silently drops a
  random effect written in the formula.** The crossed-model builder
  strips the random part and re-adds one intercept per stratum dimension
  plus the intersection intercept, so an explicit `(1 | site)` in the
  formula vanished and its variance was misallocated to the strata or
  residual (the two-model path kept it). Such a term now raises a
  targeted error pointing to `context =` — which composes with the
  crossed model — or `decomposition = "two-model"`; a grouping supplied
  through `context =` is unaffected.

- **Longitudinal time-centering now uses the analytic sample.** The
  growth-term centering offset was computed from the full input before
  `subset`, missing-value, and weight filtering (only `ref_time` was
  later recomputed from the fitted frame). An analytic sample beginning
  well away from 0 — e.g. a `subset` that drops the early waves — could
  therefore keep `time_center = 0`, defeating the numerical-stability
  protection and risking a non-convergent or false-optimum growth fit.
  The centre is now taken over the rows that survive
  `subset`/missingness/weights, so `fit_maihda(..., subset = keep)`
  centres and fits identically to `fit_maihda(data[keep, ])`.

- **Weighted
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)/[`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  no longer let invalid sampling-weight rows drive family detection.**
  Their shared setup filtered to complete cases but not the
  finite-positive sampling-weight mask the design-weighted (WeMix)
  engine applies, so a zero/negative/`Inf`-weight row still fed the
  automatic binomial/ordinal family detection and inflated `n_obs`. A
  single zero-weight row carrying an out-of-sample outcome value (e.g. a
  `2` in an otherwise binary column) could silently flip the analysis to
  Gaussian. Those rows are now excluded before detection and counting,
  matching the analytic sample every fit uses.

- **Automatic numeric strata cut-points now use the full analytic
  sample.**
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)‘s
  tertile breaks were computed from the full input before
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  dropped rows, so `fit_maihda(..., subset = keep)` binned on different
  quantiles — and assigned rows to different strata — than fitting
  `data[keep, ]` (in one bimodal example the breaks reached the max of a
  cluster the subset had removed).
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  gains an internal `bin_rows` argument, and
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  now passes the rows that survive **every** engine-specific exclusion —
  `subset`, a missing outcome/covariate (via the model frame), a missing
  precision weight, and an invalid (non-finite or non-positive) sampling
  weight — computed with the same analytic-frame mask the fit uses. An
  excluded row carrying an extreme but non-missing value of the binning
  variable (a zero-weight or missing-outcome row, say) can no longer
  pull the quantiles out and silently redefine other rows’ strata; in a
  reproducer the top tertile break moved from 1060 to 120 once such rows
  were excluded.
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  called directly is unchanged.

- **[`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  now check the non-stratum random-effects structure.** Their
  comparability checks pinned the outcome, family, sample, weights, and
  stratum assignment, but not the *other* grouping terms — so a
  `(1 | stratum)` model and a `(1 | stratum) + (1 | site)` model were
  accepted, folding a changed variance decomposition into the claimed
  covariate-attributable PCV.
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  now errors when the models’ random-effects grouping structure differs
  (or a shared non-stratum grouping is reassigned);
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  adds a matching “random-effects structure” comparability warning.

- **The brms count longitudinal VPC band now propagates residual
  uncertainty draw-by-draw.** For a Poisson/negative-binomial growth
  fit, the time-varying level-1 residual was built from posterior-*mean*
  fixed-effect predictors and a posterior-mean random-effect variance
  correction — a single value (Poisson) or one varying only through the
  dispersion draws (NB) — while the between-level variances in the same
  VPC ratio varied per draw, understating the credible band for low or
  strongly time-varying counts. The residual is now computed per
  posterior draw (per-draw linear predictor and per-draw random-effect
  variance, reusing the helper the cross-sectional count VPC already
  uses), falling back to the previous plug-in only if the draw axes
  cannot be aligned. The lme4 engine and non-count families are
  unchanged.

- **Contextual binary
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  now reports the AUC/MOR trajectory.** It suppressed the
  discriminatory-accuracy columns for any `context =` fit on the stale
  assumption that the AUC would include the context effect; the headline
  AUC is in fact the intersectional-scope concordance that *excludes* it
  (and the MOR is the between-stratum quantity), so both share the scope
  of the net-of-context `Step_PCV`/`Total_PCV` columns — matching
  [`summary.maihda_model()`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md),
  which was already updated. The columns are no longer dropped.
  ([`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)’s
  documentation, which still described the old suppressed behaviour, has
  been corrected to match.)

- **The two-model PCV now rejects a crossed-dimensions fit instead of
  silently using only its interaction variance.**
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  (and
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)/[`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md),
  all via
  [`extract_between_variance()`](https://hdbt.github.io/MAIHDA/reference/extract_between_variance.md))
  read the between-stratum variance of a
  `decomposition = "crossed-dimensions"` fit as the single intersection
  (`stratum`) random-effect variance. But for a crossed-dimensions model
  the between-stratum variance is the *sum* of the additive dimension
  effects and the interaction (`maihda_cc_partition()`:
  `between = additive + interaction`), exactly as the response-scale VPC
  and the MOR already read it. Using only the interaction component
  silently dropped the additive part and could even reverse the PCV’s
  sign (a reproducer moved from a reported −0.0135 to a correct
  +0.0255). Such a fit is now rejected with a pointer to
  `maihda(decomposition = "crossed-dimensions")`, whose summary reports
  the additive/interaction shares (the crossed-dimensions analogue of
  the PCV); the canonical single-`(1 | stratum)` PCV is unchanged, and
  the cc-aware
  VPC/MOR/[`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  paths (which never routed through this scalar) are unaffected.

- **Longitudinal id/stratum validation now runs on the analytic
  sample.** The repeated-measures check and the cross-stratum id check
  operated on the full resolved data, before `subset`, missingness, and
  weight filtering. A row the fit drops could therefore inject a
  spurious `(id, stratum)` pairing — an excluded occasion assigning one
  person to a second stratum made `fit_maihda(..., subset = keep)` fail
  with an “id appears in more than one stratum” error while fitting the
  equivalent `data[keep, ]` succeeded. Both checks now key off the same
  analytic mask used for time-centering and fitting, so a fit and its
  pre-filtered equivalent validate identically; the full-data structural
  gate in
  [`maihda_validate_longitudinal()`](https://hdbt.github.io/MAIHDA/reference/maihda_validate_longitudinal.md)
  is retained (it can only be over-lenient, never over-strict).

- **The brms group-variance parser no longer corrupts grouping names
  containing `__`.** `maihda_group_variance_draws_brms()` recovered each
  grouping factor’s name by deleting everything from the *first* `__` in
  the `sd_<group>__<coef>` posterior column, so a group literally named
  `site__id` (a valid contextual or crossed grouping) was truncated to
  `site` — then reported missing by the by-name lookups in the
  crossed-dimensions/contextual brms summary. The group is now recovered
  by splitting on the *last* `__` (brms coefficient names never contain
  `__`), which keeps names with `__` intact and still lumps a group’s
  coefficients together so the random-slope rejection is preserved.

### Testing / CI

- **Added a routine `integration-tests` CI job plus regression tests for
  the estimand issues above.** The optional-backend tests (WeMix,
  ordinal, negative-binomial, response-scale VPC, weighted counts) are
  guarded by `skip_on_cran()`, so the ordinary `R CMD check` did not
  exercise them; a new `integration-tests` job runs the full non-brms
  suite with `NOT_CRAN=true` on every push/PR (brms keeps its own
  `brms-tests.yaml`). New regression tests pin the `estimation` basis,
  the precision-weight AUC, the intersectional AUC scope, the
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)/[`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  boundary-null rejection, the crossed-dimensions response-scale VPC
  partition, and the contextual-fit discriminatory accuracy / response
  VPC.

- **Added `test-audit-2026-07-13.R`** pinning the corrected behaviour
  for the nine defects above: the external-offset prediction path, the
  crossed-dimensions extra-random-effect rejection, longitudinal
  centering on the analytic sample, invalid-weight exclusion from
  weighted family detection, subset-consistent auto-bin cut-points, the
  PCV random-effects-structure guard, the draw-by-draw brms count VPC
  residual (brms-gated), and the contextual
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  AUC/MOR trajectory.

## MAIHDA 0.2.1

CRAN release: 2026-07-09

### New Features

- **New
  [`maihda_describe()`](https://hdbt.github.io/MAIHDA/reference/maihda_describe.md):
  pre-model “Table 1” sample descriptives, built from the same machinery
  as the model** ([\#60](https://github.com/hdbt/MAIHDA/issues/60)). The
  pre-model counterpart of
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md):
  one call reports the total and complete-case analytic sample, observed
  vs. expected intersectional strata (empty cells enumerated, small
  cells flagged — never dropped), per-dimension distributions,
  **family-aware** per-stratum outcome summaries (proportions for a
  binomial outcome, category scores for an ordinal one — never a
  Gaussian mean/SD on a 0/1 outcome), missingness accounting, per-unit
  context tables for a contextual design, weighted counts/means under
  `sampling_weights`, and concrete data-quality warnings (auto-binned /
  ID-like / linear-numeric dimensions, sparse strata, concentrated
  outcome missingness, weakly identified contexts). Strata IDs, labels,
  and counts are guaranteed to match a subsequent
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  /
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  / [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  because the strata shorthand, family detection, and analytic row
  selection are resolved by the same shared internals as the fitters.
  Accepts a formula + data, or a fitted `maihda_model` /
  `maihda_analysis` to describe the exact analytic sample post hoc;
  ships with [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  (`"stratum_size"`, `"outcome"`, `"missingness"`), and every table is
  an export-ready data frame.

- **[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now accept `group` and `context` together (stratified × contextual
  MAIHDA).** Supplying both now runs the stratified comparison where
  **each per-group fit is itself a contextual cross-classified model** —
  `context` is forwarded into every per-group
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  call.

- **`plot(type = "predicted")` gains an `order_by` argument and now
  orders the strata by predicted value by default (a ranked caterpillar
  plot).** The predicted-values plot previously drew the strata in
  native stratum order, while `maihda_table()$strata` already ranked
  them by predicted outcome.

- **`plot(<maihda_analysis>, type = "vpc")` gains a `model` argument to
  show the null VPC, the adjusted VPC, or both as one change plot.** A
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis carries both a null and an adjusted model, but
  `plot(a, type = "vpc")` only ever drew the null-model VPC, with no way
  to see the adjusted partition and no cue as to which model was shown.

- **[`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  gains a `context` argument for a contextual cross-classified stepwise
  PCV.**

- **New
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md):
  order-invariant PCV attribution via Shapley values and dominance
  analysis** ([\#63](https://github.com/hdbt/MAIHDA/issues/63)).
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)’s
  per-step PCV depends on the entry order of the variables;
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)
  treats the PCV as a value function over variable subsets,
  `v(S) = (V0 − V(S)) / V0`, and apportions the full-model Total PCV
  fairly across the predictors. `method = "shapley"` (the default)
  averages each variable’s marginal PCV over all entry orders — exact up
  to ~10 variables (every subset model fit once and cached), with a
  Monte-Carlo permutation approximation (`approx = "montecarlo"`,
  per-variable MC standard errors and a convergence warning) beyond;
  `method = "dominance"` adds Budescu’s conditional and complete
  dominance detail (general dominance coincides with the Shapley
  values); `method = "sequential"` keeps the order-dependent path for
  continuity. All methods satisfy the efficiency identity — the
  contributions sum exactly to the full-model Total PCV — and are
  reported *signed*, so a suppressor variable shows up as negative
  rather than being normalised away. Attribute among the stratum
  dimensions to split the additive share (“which dimension drives the
  additive between-stratum inequality”), or among individual-level
  covariates as the order-free counterpart of
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  (the latent-scale rescaling caveat of
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  carries over for non-Gaussian families). Full parity with
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  on `engine`/`family`/`context`/`sampling_weights`, the shared
  complete-case analytic sample, and auto-binned dimension
  reconstruction — the two functions now literally share their setup
  code. Optional parametric-bootstrap CIs per contribution
  (`bootstrap = TRUE`; lme4 + exact attribution only, costing `n_boot ×`
  the number of subset models). Ships with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

### API Changes

- **Added a [`predict()`](https://rdrr.io/r/stats/predict.html) S3
  method for `maihda_model` objects.** A thin alias of
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md),
  so the `predict(model, type = "strata")` call the documentation and
  summary output recommend now works directly.

- **[`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  has been renamed to
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).**
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  remains as a deprecated alias that warns and forwards, and will be
  removed in a future release. The result object carries the estimate in
  `$pcv` and is classed `c("pcv_result", "pvc_result")`; the old `$pvc`
  element is kept as a deprecated duplicate, so existing code and
  objects saved by earlier versions keep working (including
  [`print()`](https://rdrr.io/r/base/print.html)). All internal
  consumers
  ([`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md),
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md),
  the tidiers, and the Shiny app), the documentation, and the vignettes
  now use the corrected name.

### Bug Fixes

- **A brms addition-term response survived automatic strata creation
  with the wrong model.**
  [`reformulas::nobars()`](https://rdrr.io/pkg/reformulas/man/nobars.html)
  descends into a formula’s left-hand side, so
  `y | trials(n) ~ x + (1 | a:b)` was silently rewritten to
  `x ~ (1 | stratum)` — the response replaced by a predictor — in
  automatic strata creation, the crossed-dimensions builder, and the
  longitudinal growth-formula builder. Bar terms are now stripped from
  the right-hand side only (a shared internal helper used at every
  fixed-part extraction), preserving `| trials()`, `| weights()`, and
  [`cbind()`](https://rdrr.io/r/base/cbind.html) responses exactly.

- **brms sampling weights are now normalized on the analysis sample.**
  Weights were normalized to mean 1 *before* brms dropped rows
  incomplete on the model variables, so the weights actually fitted
  could have an arbitrary mean — rescaling the pseudo-posterior’s
  effective sample size (a two-row case left a single surviving weight
  of 0.02, a ~50× tempered likelihood). Rows incomplete on the outcome,
  the (possibly transformed) predictors, or the grouping variables are
  now dropped, with a warning, before normalization.

- **Longitudinal ids reused across strata are rejected.** `(time | id)`
  treats every row sharing an id as the same person, so ids numbered
  within a site or group (person “1” in every site) silently merged
  different people’s trajectories. An id appearing in more than one
  stratum now errors with guidance (build a globally unique id, or fix a
  stratum that changed between occasions at a reference occasion).

- **Fractional lme4 precision weights no longer produce an impossible
  AUC.** Non-integer prior weights (e.g. 1.5) tripped the
  aggregated-binomial heuristic — `round(y * 1.5) = 2` successes out of
  1.5 trials implies negative failures — reporting AUC values above 1
  and doubled case totals. Only integer trial counts mark aggregation
  now; fractional-weight fits get the weighted Mann–Whitney concordance
  (reported with `weight_type = "precision"`), and the weighted AUC
  internally rejects negative case/control mass.

- **AUC and MOR now summarize the same model scope.** For a model with
  random effects beyond the intersectional partition (a contextual
  `(1 | school)` or an explicit `(1 | site)`),
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
  computed the AUC from full predictions (site effects included) while
  the MOR used only the stratum variance — a strong site effect could
  carry a 0.90 AUC over a negligible stratum effect. The headline `auc`
  is now the intersectional-scope concordance (matching the MOR), with
  the full-model AUC reported alongside as `auc_full` (`auc_scope` names
  the scope).
  [`maihda_mor()`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md)
  on a crossed-dimensions fit now sums the additive dimension variances
  with the interaction component instead of reading the interaction
  alone.

- **The response-scale VPC integrates over non-stratum random effects.**
  [`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  simulated only the stratum effect and silently ignored any other
  random intercepts, overstating the stratum share (0.049 reported where
  a coherent nested Monte Carlo gives 0.013). The simulation now
  integrates over the combined non-stratum effects; the result reports
  their summed latent variance as `var_other`.

- **PCV bootstrap draws on the zero-variance boundary are reported, and
  a degenerate PCV denominator is rejected.** Bootstrap draws whose null
  model estimated a zero between-stratum variance were silently dropped,
  making the percentile interval conditional on a positive null variance
  with no signal; any boundary mass now warns, is returned as
  `n_boot_boundary`, and is repeated by
  [`print()`](https://rdrr.io/r/base/print.html). Boundary detection
  uses lme4’s relative singularity tolerance — a strictly positive but
  effectively singular variance (1e-9) previously passed a strict-zero
  check and produced PCVs in the thousands with intervals like
  \[-6.6e16, 1\];
  [`calculate_pcv()`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  now also rejects an effectively singular null model.

- [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  on a zero-row data frame now fails immediately with a clear message
  instead of a base replacement error.

- [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  (and the new
  [`maihda_describe()`](https://hdbt.github.io/MAIHDA/reference/maihda_describe.md))
  errored on an aggregated binomial outcome combined with the strata
  shorthand and no covariates —
  e.g. `cbind(successes, failures) ~ (1 | gender:race)` — because
  `nobars()` returns a bare call (not a formula) for a call-valued
  response with a bars-only right-hand side.

- The brms engine’s linear-predictor reads were broken under current
  brms (\>= ~2.16): `posterior_linpred()` ignores a `summary = TRUE`
  argument and returns the raw draws matrix, so several brms paths
  errored.

- A single-row `newdata` errored for response-scale predictions from a
  brms cumulative (ordinal) fit.

- Longitudinal (growth-curve) fits on a time axis that does not start at
  0 could silently converge to a wrong solution; the growth terms are
  now fit on internally centered time.

- The Poisson / negative-binomial latent-scale level-1 variance was
  evaluated at the conditional fitted means (BLUPs included); it is now
  evaluated at the marginal expected count, matching the cited
  references.

- The longitudinal PCV compared REML variances across models with
  different fixed effects; it is now computed from ML-refitted growth
  models.

- Stratum display labels are no longer whitespace-padded for mixed-type
  dimensions.

- The longitudinal slope PCV read the wrong covariance cell for
  `time_degree >= 2`

### Improvements

- [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) now
  warns when a numeric stratum-defining dimension would enter the
  adjusted model as a linear fixed effect instead of categorical main
  effects.

- calculate_pcv(bootstrap = TRUE)\` on a non-lme4 engine now explains
  what is available.

- Documented the latent-scale rescaling caveat on the PCV.

### Documentation

- The bootstrap Monte Carlo error is now labelled as what it is — the
  Monte Carlo standard error of the bootstrap **mean**
  (`sd(draws)/sqrt(n)`), a coarse gauge of bootstrap noise — in
  [`print()`](https://rdrr.io/r/base/print.html) output and the
  internals; it is not the sampling uncertainty of the percentile
  interval’s endpoints.
- A batch of documentation corrections from the 0.2.1 audit:
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  documents its full six-element return; the `"upset"` plot type is
  documented as the three-panel patchwork it returns; the
  group-comparison docs no longer deny cross-classification when
  `context =` is supplied;
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)’s
  description acknowledges its single-model crossed-dimensions mode;
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  documents that brms posterior intervals are returned without
  bootstrapping;
  [`pcv_importance()`](https://hdbt.github.io/MAIHDA/reference/pcv_importance.md)’s
  exact-cost examples count the null fit (256 fits at k = 8, 1024 at k =
  10); the README calls `Context_Variance` an absolute variance rather
  than a share; fit diagnostics are documented for all four engines;
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
  documents its `weighted`/`weight_type` fields; and the CITATION year
  is 2026.
- Fixed a mislabelled figure in the reporting vignette: the
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) caterpillar
  plotted the null model’s stratum estimates, the total between-stratum
  deviations, but the text called them “interaction estimates” (the pure
  interactions require `tidy(a, which = "adjusted")`; the vignette now
  says so).
- The plot-interpretation vignette no longer claims a bare
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  fit is “equivalent” to the
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis for its plots:
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the
  analysis routes the VPC/predicted/shrinkage views to the null model
  and the effect-decomposition views to the adjusted model, which a
  single fit cannot do.
- The group-comparison vignette now explains its own `group_pcv` figure:
  with the real PISA data every country’s PCV is essentially 1 (the
  gender-by-SES gaps are purely additive, so the per-country adjusted
  fits are singular-boundary), and the flat row of 100% bars is itself
  the substantive result, not a rendering problem.
- Restored the accidentally truncated opening paragraph of the “Planning
  a MAIHDA analysis” vignette.

## MAIHDA 0.2.0

CRAN release: 2026-07-02

### New Features

- **Removed the `"risk_vs_effect"` plot type (and the internal
  `plot_risk_vs_effect()`).** The quadrant view of each stratum’s mean
  prediction against its stratum random effect largely duplicated what
  `"effect_decomp"` and `"predicted"` already show more directly, so it
  has been dropped from
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for both
  `maihda_model` and `maihda_analysis` objects, from `type = "all"`, and
  from the Shiny app’s plot menu.
- Removed the `"ternary"` plot type and the
  `compute_maihda_ternary_data()`, `maihda_ternary_plot()`, and
  `plot_maihda_ternary()` functions. The ternary diagnostic normalised
  each stratum’s additive, intersection-specific, and uncertainty
  signals to sum to 1, which discards effect magnitude and puts two
  effect components and an estimation-error term on a single triangle.
  Its content is covered more directly by `"effect_decomp"` (the
  additive-vs-intersectional split with magnitudes intact) and the
  `"predicted"` / `"upset"` views (uncertainty shown as intervals), so
  it has been dropped from
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for both
  `maihda_model` and `maihda_analysis` objects, from `type = "all"`, and
  from the Shiny app, and the optional `ggtern` dependency removed.
- **Added the `"upset"` plot type and
  [`maihda_upset_size()`](https://hdbt.github.io/MAIHDA/reference/maihda_upset_size.md).**
  An UpSet-style alternative to `"predicted"` that replaces the long
  intersectional axis labels with a category matrix (an
  intersection-size bar, a dot matrix encoding each stratum’s level on
  every dimension, and a predicted-value panel, all sharing one column
  order). Binary 0/1 dimensions collapse to a present/absent row;
  multi-level factors get one row per level. A `quantity` argument
  switches the bottom panel between the predicted value (fixed + random)
  and the stratum random effect (interaction).
  [`maihda_upset_size()`](https://hdbt.github.io/MAIHDA/reference/maihda_upset_size.md)
  returns content-scaled figure dimensions.
- **[`theme_maihda()`](https://hdbt.github.io/MAIHDA/reference/theme_maihda.md)
  now is set as standard theme for ggplot objects**
- **The `predicted` (and longitudinal `trajectories`) plot can now keep
  the most extreme strata when it truncates, instead of the first by
  stratum order.** When there are more strata than the `n_strata` cap,
  the view dropped strata in stratum order, effectively arbitrary with
  respect to how far a subgroup sits from the population mean, so the
  most striking strata could be the ones omitted.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) gains an
  opt-in `select` argument: `"order"` (default, unchanged) keeps the
  first n_strata in stratum order; `"deviation"` keeps the n_strata
  furthest from the reference line (largest `|predicted - reference|`,
  so the most extreme strata in both directions, not just one tail). For
  a longitudinal `trajectories` plot it keeps the strata whose
  trajectories swing furthest from the population curve (peak
  `|random deviation|` over the time grid). Selection and display are
  kept separate: `select` changes which strata appear, but the x-axis
  stays in stratum order. It composes with the flag-aware cap (flagged
  strata are always kept; `select` governs the fill) and the caption
  names the rule used (“the 12 strata furthest from the reference, of
  200”).
- **The BLUP plots can now show only the flagged strata, and never hide
  a flagged stratum behind the display cap.** When there are many strata
  the `"predicted"` view caps the number drawn (`n_strata`, default 50)
  and kept them in stratum order, so a stratum flagged by
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  that fell past the cap was dropped from the figure entirely; the
  highlight could not reach it.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) gains
  `only_flagged`: `TRUE` restricts the `"predicted"` and
  `"obs_vs_shrunken"` views to the strata carrying a credibly non-zero
  interaction (and turns the highlight on with the stored diagnostic if
  it was off), with a caption naming the screen (e.g. “Showing the 7
  flagged strata (95% interval, BH-adjusted) of 200 total”) and a
  graceful captioned empty panel when none are flagged. Independently,
  whenever interactions are highlighted the `n_strata` cap on
  `"predicted"` is now flag-aware: every flagged stratum is kept and the
  remaining slots are filled in stratum order, so the signal the
  highlight exists to surface is never silently truncated away. The
  across-strata reference line is still computed from the full set, so
  filtering never shifts it. `"effect_decomp"` deliberately ignores
  `only_flagged` (its waterfall exists to show each flagged stratum’s
  place in the full distribution) and says so; the framing stays on
  “flagged”, not “significant”, consistent with the diagnostic’s
  exploratory, partial-pooling reading.
- **The interaction diagnostic now defaults to FDR control and gains an
  equivalence (ROPE) reading.**
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  now defaults to `adjust = "BH"` (Benjamini-Hochberg) rather than
  `"none"`: fitting and flagging every stratum in one call is a
  screening question, so the flags should control the false-discovery
  rate by default. Pass `adjust = "none"` (or `interactions = "none"` to
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)/[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md))
  for the uncorrected, per-stratum individual-testing view. A new `rope`
  argument adds an equivalence / smallest-interaction-of-interest
  reading (Schuirmann 1987; Kruschke 2018): supply a region of practical
  equivalence (a half-width `d` for `c(-d, d)`, or `c(lower, upper)`, on
  the link scale) and each stratum gains a `decision` of `"relevant"`
  (interval entirely outside the region), `"negligible"` (entirely
  inside), or `"inconclusive"` (straddling a bound), reported by
  [`print()`](https://rdrr.io/r/base/print.html). The documentation now
  also keeps two ideas the literature distinguishes apart – partial
  pooling regularises magnitude/sign (Gelman, Hill & Yajima 2012) while
  whether to correct depends on the inferential structure of the claim
  (Rubin 2021) – rather than treating shrinkage as itself a multiplicity
  correction.

### Bug Fixes

- **`compare_maihda(ic = TRUE)` could show information criteria on
  incompatible scales without a warning.** The comparability check only
  warned when the models differed in outcome, weights, family/link,
  analytic sample, or strata, then appended whichever
  information-criterion columns existed. Two same-family models fitted
  on different engines, e.g. a Gaussian `lme4` fit (reporting AIC/BIC)
  and a Gaussian `brms` fit (reporting WAIC/LOOIC) – agree on all the
  checked aspects, so no warning fired, yet the appended table placed
  the likelihood AIC/BIC and the Bayesian WAIC/LOOIC side by side. Those
  are different scales and are not comparable to each other (as
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)’s
  documentation already notes).
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  now emits an explicit warning whenever the appended criteria span both
  the likelihood (AIC/BIC) and the Bayesian (WAIC/LOOIC) scales, so a
  cross-engine comparison can no longer be read as if the four numbers
  were on one ruler. An all-`lme4`/all-`ordinal` (AIC/BIC only) or
  all-`brms` (WAIC/LOOIC only) comparison is unaffected.

- **A fixed interaction among the stratum dimensions
  (e.g. `Gender * Race`) silently corrupted the decomposition.** R
  expands `Gender * Race` to `Gender + Race + Gender:Race`, but the
  workflow only looked for the dimensions’ additive main effects when
  deciding whether the supplied formula was the fully-specified adjusted
  model. It found both `Gender` and `Race`, treated the fit as the
  adjusted model, and derived the null by removing only those main
  effects – leaving the fixed `Gender:Race` interaction in **both** the
  null and the adjusted formulas. That fixed cell-means term duplicates
  the intersectional `(1 | Gender:Race)` random intercept: it absorbs
  the between-stratum variance into fixed effects, pinning the stratum
  variance at a singular boundary so the PCV came back `NULL`/invalid
  (and every per-group PCV came back `NA` in
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  / `maihda(group = )`).
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now **reject** such a formula up front with an actionable error – the
  MAIHDA adjusted model is defined to carry only the dimensions’
  additive main effects, with the intersection estimated by the stratum
  random effect (write `Gender + Race`, and use
  `decomposition = "crossed-dimensions"` or
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  to quantify the interaction).
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  likewise now warns when handed a bare model that carries such a fixed
  dimension interaction (its stratum random effects are no longer the
  pure interaction the diagnostic reports). Detection is robust to
  variable order within the interaction and to non-syntactic names, and
  a legitimate covariate-by-dimension interaction (e.g. `age * Gender`)
  is left untouched.

- **Discriminatory accuracy silently dropped the AUC/MOR for a `brms`
  aggregated-binomial (`y | trials(n)`) fit.** The model layer fits such
  responses, and [`summary()`](https://rdrr.io/r/base/summary.html)
  attaches the discriminatory accuracy automatically for any
  binomial/Bernoulli fit, but
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
  only recognised the aggregated structure for `lme4`
  `cbind(success, failure)` fits. A `brms` `y | trials(n)` fit therefore
  reached
  [`maihda_auc()`](https://hdbt.github.io/MAIHDA/reference/maihda_auc.md)
  with the per-row success counts\* as the response, errored on the
  non-0/1 values, and the summary quietly omitted the DA. It now detects
  a `brms` aggregated binomial via the existing trial-count extraction
  path (`maihda_brms_trial_counts()`, which parses the `trials()`
  addition term off the stored formula – `brms` exposes no
  `weights.brmsfit`) and computes the same count-weighted C-statistic
  the `lme4` [`cbind()`](https://rdrr.io/r/base/cbind.html) path uses.
  The rows are ranked by the per-trial probability
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  returns (now normalised to a probability on both engines – see below),
  so the AUC ranks by probability and not by a count that confounds the
  probability with the number of trials; `n_case` / `n_control` are the
  total successes / failures. Bernoulli and `lme4` fits are unaffected.

- [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  (and therefore
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md))
  could **silently report a REML-based PCV for a singular
  multi-random-effect model**. The ML refit
  ([`lme4::refitML()`](https://rdrr.io/pkg/lme4/man/refitML.html),
  needed because REML between-stratum variances are not comparable
  across models with different fixed effects) was skipped whenever the
  fit was globally singular, but a fit can be singular because a
  **non-stratum** random effect (e.g. an extra `(1 | site)` grouping
  factor) is on the boundary while the stratum variance is comfortably
  nonzero. In that case the package returned the REML PCV instead of the
  intended ML PCV. The ML refit is now skipped only when the **stratum**
  variance itself is on the boundary (effectively zero, judged on
  [`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html)’s
  own relative tolerance), or when `refitML()` fails; a non-stratum
  boundary no longer suppresses it. The boundary skip for a
  zero-variance stratum is preserved, so the zero-variance guard in
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  is not masked.

- **Aggregated-binomial stratum predictions were row-weighted instead of
  trial-weighted.** For a `cbind(success, failure)` (or `y | trials(n)`)
  fit the rows of an intersectional stratum can carry very different
  numbers of binomial trials, but the per-stratum prediction means
  averaged the fitted probabilities / linear predictors with unit
  weights, so a 5-trial row counted as much as a 200-trial row. The unit
  weighting is correct for the observed stratum summaries (which sum
  successes over summed trials, the trials are already in the
  denominator), but wrong for the model predictions. Prediction
  aggregation now uses a dedicated `maihda_prediction_weights()` that
  weights each row by its binomial trial count – read from
  `stats::weights(model, type = "prior")` for an `lme4`
  [`cbind()`](https://rdrr.io/r/base/cbind.html) fit, or parsed from the
  formula’s `trials()` addition term for a `brms` `y | trials(n)` fit
  (which exposes `model.frame.brmsfit` but no `weights.brmsfit`, so the
  prior-weights route is unavailable and the counts would otherwise fall
  back to unit weights) – while the observed numerator/denominator path
  is unchanged. This corrects the predicted stratum means and the
  `w_sum` weights feeding
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md),
  the predicted-strata / effect-decomposition plots, and the
  prediction-deviation panels; unweighted and non-binomial fits are
  unaffected (the weights reduce exactly to the previous values).

- **`brms` cumulative (ordinal) stratum predictions were on the wrong
  response scale.** The stratum-level prediction helper applied the
  **scalar inverse link** for the `engine = "brms"`,
  `family = "ordinal"` response scale, returning a single cumulative
  probability in `[0, 1]` rather than the expected category score
  `sum_k k * P(Y = k)` in `[1, K]` that the rest of the package
  documents and computes – the individual
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  path already collapses the
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
  category-probability array to that score, and the `clmm`
  (`engine = "ordinal"`) stratum path already used it. This silently
  mis-scaled (and mis-ranked)
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)
  and the stratum plots for a Bayesian cumulative model. The brms
  stratum helper now maps the latent location to the expected category
  score with the same shared cumulative helpers (using the
  posterior-mean thresholds), so the brms and `clmm` cumulative paths
  agree and match the documented response scale; the link scale is
  unchanged.

- **[`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  returned a meaningless scalar diagnostic for a longitudinal fit.** A
  direct call on a longitudinal
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis (or a bare longitudinal model) fell through to the
  crossed-dimensions branch and reported a single per-stratum BLUP –
  which drops the random slope and so describes a cross-section, not the
  trajectory. It now errors with a pointer to the trajectory tools,
  matching the automatic-attachment path that already skips longitudinal
  objects.

- **The Bayesian `pd` column was mislabelled.**
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  reported `pd = P(BLUP > 0)` under the heading “probability of
  direction”, so a strongly **negative** interaction (whose direction is
  near-certain) showed `pd` near 0 rather than near 1. `pd` is now the
  conventional probability of direction `max(P(BLUP > 0), P(BLUP < 0))`
  (in `[0.5, 1]`, à la `bayestestR::p_direction`), with the sign
  reported separately in the existing `direction` column.

- **[`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  silently accepted an unseen stratum on the `wemix` and `ordinal`
  paths.** When `newdata` already carried a `stratum` column, the
  preparation step returned early and skipped the unseen-stratum check,
  so a misspelled or genuinely new stratum flowed through to the
  WeMix/`clmm` linear-predictor helpers, which map a missing random
  effect to 0, yielding a fixed-only prediction that looked valid and
  contradicted both the documented contract (unseen strata are an error,
  as for `type = "strata"`) and `lme4`’s default behaviour.
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  now rejects an unseen stratum by default for every engine and
  prediction type, whether the stratum is supplied directly or rebuilt
  from the grouping variables. A new `allow_new_levels` argument
  (default `FALSE`) opts into the previous behaviour explicitly: for
  `type = "individual"` it returns a
  population-average(fixed-effects-only) prediction for unseen strata,
  dropping the stratum random effect (forwarded as `allow.new.levels` to
  lme4 and `allow_new_levels` to brms). Stratum-level predictions have
  no random effect to report for an unseen stratum, so they remain an
  error regardless.

- **`predict_maihda(scale = "response")` returned expected counts, not
  probabilities, for a `brms` aggregated-binomial (`y | trials(n)`)
  fit.** `brms`’s
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) /
  `posterior_epred()` return the expected number of successes
  (`trials * p`) for a binomial fit with a `trials()` term, whereas
  `lme4`’s `type = "response"` for a `cbind(success, failure)` fit
  returns the per-trial probability. The same
  `predict_maihda(scale = "response")` call therefore produced two
  different scales depending on the engine, so anyone ranking or
  comparing predicted risks across `brms` strata was comparing expected
  counts (which confound the probability with the per-row number of
  trials) rather than probabilities. The `brms` response-scale
  prediction is now normalised by the per-row trial counts (evaluated on
  the prediction `newdata`), so both engines return a per-trial
  probability in `[0, 1]`; the discriminatory-accuracy AUC, the Shiny
  app’s absolute-prediction panel, and any downstream ranking now read a
  true probability. Bernoulli and non-binomial fits are unaffected (no
  `trials()` term -\> no normalisation).

- **`predict_maihda(allow_new_levels = TRUE)` did not give the
  documented zero-effect prediction for an unseen stratum on the `brms`
  engine.** The contract is to drop the stratum random effect (treat it
  as zero) for a stratum the model never saw. The `brms` path merely
  forwarded `allow_new_levels = TRUE`, but `brms`’s default
  `sample_new_levels = "uncertainty"` draws a new-stratum effect from
  the estimated random-effects distribution rather than zeroing it, so
  the prediction silently carried a random per-call group effect.
  Unseen-stratum rows are now split off and predicted with an
  `re_formula` that excludes the stratum term, while seen-stratum rows
  keep their estimated effect (a blanket exclusion would have dropped
  the seen strata’s effects too). Any other random effect the unseen row
  participates in – a contextual `(1 | school)` intercept from
  `fit_maihda(context = )`, a longitudinal `(time | id)` growth term –
  is kept , exactly as `lme4`’s `allow.new.levels` zeroes only the
  unseen level’s effect and retains seen ones; for the usual
  single-stratum model the excluding `re_formula` is `NA`, the
  fixed-effects-only population average. `lme4` and the
  `wemix`/`ordinal` paths already behaved this way and are unchanged, so
  the engines now agree on the same model.

- **The crossed-dimensions decomposition reported `NaN` additive /
  interaction shares when there was no between-stratum variance.** The
  shares split the between-strata variance, so for a degenerate fit
  whose additive and interaction variances are both estimated at exactly
  zero they are `0 / 0 = NaN`, which leaked through
  [`summary()`](https://rdrr.io/r/base/summary.html) and the
  comparison/tidier outputs. `maihda_cc_partition()` now returns `NA` (a
  defined “undefined”) for the shares when `between == 0`, and for the
  VPC when the total variance is `0`;
  [`print()`](https://rdrr.io/r/base/print.html) shows
  `NA (no between-strata variance to split)`. Non-degenerate fits are
  unchanged.

- **A formula [`offset()`](https://rdrr.io/r/stats/offset.html) term was
  silently dropped from `wemix` and `ordinal` predictions.**
  [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
  and [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html)
  both honour an `offset(.)` term in the model formula when fitting, but
  the manual linear-predictor helpers (no `predict.clmm` exists, and
  WeMix’s own [`predict()`](https://rdrr.io/r/stats/predict.html) offers
  no fixed-only/scale control) rebuilt the linear predictor as
  `X %*% beta (+ u)` from the design matrix only – and
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) never
  produces a column for an offset term, so the offset was omitted. Every
  link-scale prediction (and the response-scale value derived from it)
  was therefore off by exactly the per-row offset:
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  for both engines, and the ordinal
  [`plot_prediction_deviation_panels()`](https://hdbt.github.io/MAIHDA/reference/plot_prediction_deviation_panels.md)
  probabilities, which rebuild the `clmm` location independently. The
  helpers now evaluate the formula’s offset on the prediction data
  ([`stats::model.offset()`](https://rdrr.io/r/stats/model.extract.html)
  of the rebuilt model frame) and add it back, including for an
  offset-only null model whose location is otherwise identically zero.
  Fits with no offset are unaffected.

- **[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  chose the engine from the raw outcome, before `subset`/weights.** The
  wrappers pass `engine` explicitly to every sub-fit, so they
  pre-selected `engine = "ordinal"` from an ordered-factor outcome – but
  they checked the raw column, while
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  detects the family on the analytic sample (after `subset` and weight
  filtering). An ordered 3-level outcome subset to two observed levels
  is a binary analytic sample: a direct
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  fits it as `binomial`/`lme4`, but the wrappers pinned
  `engine = "ordinal"` and then errored
  ([`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)) or
  failed every group
  ([`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md))
  with an engine/family contradiction. The ordinal pre-check now runs on
  the same analytic sample
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  uses (the forwarded `subset`/weights are resolved against the data
  mask first), so the engine choice can no longer contradict the
  resolved family. Relatedly,
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) now
  forwards the evaluated `subset`/`weights` (not the raw expressions) to
  its derived null/adjusted fits, mirroring
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md):
  a response-referencing `subset`
  (e.g. `subset = y %in% c("lo", "mid")`) was otherwise re-evaluated by
  each derived fit against `$original_data`, whose response has already
  been recoded to 0/1, selecting zero rows.

- **Auto-binning a numeric stratum dimension could silently overwrite a
  user column.** When
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  (or the `(1 | a:b)` shorthand) auto-bins a numeric dimension `v`, the
  adjusted-model and prediction machinery add an internal factor column
  named `.maihda_dim_<v>` (referenced by the adjusted /
  crossed-dimensions formulae and rebuilt for prediction `newdata`).
  Neither writer checked whether the user’s data already carried a
  column of that name, so an existing `.maihda_dim_<v>` variable was
  clobbered – and, because that column then enters the model, the fit or
  prediction changed silently. The `.maihda_dim_` prefix is now
  **reserved**:
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  rejects a pre-existing `.maihda_dim_<v>` column for any dimension it
  is about to auto-bin, with an actionable rename hint (or set
  `autobin = FALSE`), mirroring the existing reserved-weight-column
  guard. The internal name is now generated through one shared helper so
  the fit-time and predict-time writers cannot drift, and predicting on
  a fitted model’s own data (which legitimately carries the package’s
  copy) is unaffected.

- **[`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)
  could show REML variance rows alongside an ML-based PCV without saying
  so.**
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  refits a Gaussian `lme4` fit with maximum likelihood before
  differencing the between-stratum variances (REML variances are not
  comparable across models with different fixed effects), but the
  table’s *Between-stratum variance* / *SD* and *VPC/ICC* rows are each
  model’s own (REML) estimate – the quantities
  [`summary()`](https://rdrr.io/r/base/summary.html) reports. The two
  are on different variance bases by design, so
  `PCV != (var_null - var_adj) / var_null` read off the table, which
  looks like an inconsistency.
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)
  now attaches a `models_note` (shown by
  [`print()`](https://rdrr.io/r/base/print.html)) explaining this
  whenever the PCV’s variance basis actually differs from the displayed
  rows; it stays silent for already-ML engines
  (`glmer`/`brms`/`wemix`/`ordinal`) and for a boundary null where
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  keeps the REML fit. The displayed numbers are unchanged, so the table
  still agrees exactly with
  [`summary()`](https://rdrr.io/r/base/summary.html).

### Improvements

- **Console output is now colour-coded** (via `cli`). The print methods
  across the package –
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analyses, model/summary, PCV, discriminatory accuracy, the interaction
  diagnostic, group comparison, tables, information criteria, and the
  rest – bold section titles, accent the headline values, and dim
  secondary notes, reusing the **plot accent colour** so the console and
  the figures agree. The palette is deliberately **valence-neutral**:
  colour marks emphasis (a finding to look at, a neutral conclusion, a
  de-emphasised note), never good-vs-bad – so, e.g., a `negligible`
  equivalence result is shown in a neutral colour, not green. It
  **degrades to plain text automatically** wherever ANSI is unsupported
  (knitr/vignettes, `R CMD check`, `testthat`, `NO_COLOR`), so rendered
  and captured output is byte-for-byte unchanged. `cli` joins `Imports`.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  (and `maihda(group = )`) gain a **`var_between_adjusted_ml`** column:
  the adjusted model’s *actual* between-stratum variance, read straight
  off the adjusted fit on the maximum-likelihood scale the PCV is
  computed on. The existing `var_between_adjusted` column is, by design,
  a derived coherence quantity (`var_between * (1 - pcv)`, on the REML
  `var_between`/`vpc` scale so the table satisfies
  `pcv = (var_between - var_between_adjusted) / var_between` exactly)
  and is **not** the adjusted fit’s own variance; the documentation now
  states this explicitly and `var_between_adjusted_ml` reports the
  literal adjusted variance for anyone who needs it. The two differ only
  by the small REML-vs-ML gap in the null variance.

## MAIHDA 0.1.11

CRAN release: 2026-06-18

### New Features

- **[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) now
  reports the intersectional interactions by default.** “Which strata
  genuinely interact” is the scientific payoff of MAIHDA, so it no
  longer needs a separate call:
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  computes
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
  on the adjusted (or crossed-dimensions) model, stores it as the
  `interactions` slot, and surfaces a one-line headline in
  [`print()`](https://rdrr.io/r/base/print.html) (flagged count, the
  strongest stratum, and – crucially – the multiplicity stance actually
  used, so an uncorrected scan is never silently read as corrected). The
  computation is essentially free (it reads the stratum estimates the
  summary already holds; no refit), and it is skipped for a longitudinal
  decomposition (whose interaction is a trajectory, not a scalar).
  Control it with the new `interactions` argument: `TRUE` (default)
  computes it with the diagnostic’s own default correction; `FALSE`
  skips it; or pass a `p.adjust` method name to choose the correction.
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  gains the same argument as an **opt-in** (`interactions = FALSE` by
  default, since a single fit is often a null model where the diagnostic
  does not apply). `plot(..., highlight_interactions = TRUE)` now reuses
  the stored diagnostic, so the plot highlights and the printed headline
  can no longer disagree.
- Added **`broom`-style
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html)
  methods** for `maihda_model`, `maihda_summary`, and `maihda_analysis`,
  so MAIHDA results drop straight into tidy data for `ggplot2`,
  `gt`/`flextable` tables, and other downstream tooling.
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
  stratum (intersection) random-effect estimates by default (with the
  human-readable intersectional label), or the variance-components table
  (`component = "variance"`) or fixed effects (`component = "fixed"`),
  all in broom’s `estimate`/`std.error`/`conf.low`/`conf.high` shape;
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) on a
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis takes `which = c("null", "adjusted")`.
  [`glance()`](https://generics.r-lib.org/reference/glance.html) returns
  the **MAIHDA headline as a one-row tibble** – the VPC/ICC (with its
  bootstrap/posterior interval), and for a
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis the **PCV**, the adjusted-model AUC, the additive/interaction
  shares, AUC and MOR for a binary outcome, plus
  `n_strata`/`nobs`/`engine`/`family` – a row no generic mixed-model
  tool emits, since the PCV needs the null+adjusted pair. The layout is
  uniform across the `lme4`, `brms`, `wemix`, and `ordinal` engines. The
  generics come from the lightweight **`generics`** package (the same
  ones `broom`/`broom.mixed` re-export), so
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html)
  dispatch whether the user has `broom`, `generics`, or just `MAIHDA`
  loaded, with no hard `broom` dependency; the raw fixed-effect/per-row
  tidying `broom.mixed` already does on the underlying fit is
  intentionally not duplicated.
- Added **longitudinal (growth-curve) MAIHDA** – the life-course
  extension of Bell, Evans, Holman & Leckie (2024, *Soc Sci Med*
  351:116955) – via new `id`, `time`, and `time_degree` arguments on
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  and [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md).
  Supplying `id` (the person/unit measured repeatedly) and `time` (a
  numeric measurement-time column) fits a **3-level growth curve** –
  occasions within individuals within intersectional strata – with a
  random intercept *and* slope on time at *both* the individual and
  stratum levels (`outcome ~ time + (time | id) + (time | stratum)`; the
  growth slopes are added automatically, so you write only the strata
  shorthand `(1 | var1:var2)`). The between-stratum variance, and
  therefore the **VPC, is then a function of time**
  (`VarS(t) = a(t)' Sigma_s a(t)`):
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the
  baseline (reference-time) VPC, the full time-varying VPC trajectory,
  and the stratum/individual intercept-slope-covariance blocks, and a
  new plot type `"vpc_trajectory"` draws the VPC-over-time curve (with
  `"trajectories"` for the predicted per-stratum mean trajectories;
  `plot(type = "vpc")`/`"all"` route there for a longitudinal fit).
  `maihda(decomposition = "longitudinal")` (selected automatically when
  `id`/`time` are supplied) fits a null and an adjusted growth model –
  the adjusted model adding the dimensions’ main effects **and their
  `dim:time` interactions** – and reports the **PCV separately for the
  baseline (intercept) and the slope variance**: the
  additive-vs-multiplicative split of the intersectional *trajectory*
  inequality (the paper’s “partly multiplicative but mostly additive”
  finding), returned as a `maihda_long_pcv` object with `pcv_intercept`,
  `pcv_slope`, and a time-specific `pcv_t` (and a `"pcv_trajectory"`
  plot). All families with a defined level-1 variance are supported
  (`gaussian`/`binomial`/`poisson`/`negbinomial`, latent scale for
  non-Gaussian) on `engine = "lme4"` (any degree) and `engine = "brms"`
  (linear growth; posterior credible bands on the VPC trajectory).
  `predict_maihda(type = "strata")` returns each stratum’s deviation at
  the baseline (reference) time plus its raw intercept and slope (a
  stratum is now a *trajectory*; the baseline deviation, evaluated at
  `ref_time = min(time)`, is the longitudinal analogue of a
  cross-sectional stratum BLUP and differs from the raw time-0 intercept
  whenever time is not zero-referenced). The intercept-only guards
  elsewhere are untouched – a longitudinal fit is tagged and routed to
  the time-varying path, while every other model still rejects random
  slopes – and
  [`extract_between_variance()`](https://hdbt.github.io/MAIHDA/reference/extract_between_variance.md)/[`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  reject a longitudinal model with a pointer to the time-varying
  decomposition. Design-weighted, contextual, `wemix`/`ordinal`, and
  group-comparison longitudinal models are out of scope (each errors). A
  new bundled dataset **`maihda_long_data`** (600 persons x 5 waves,
  gender x ethnicity x education strata, with constructed
  mostly-additive trajectory differences) demonstrates it.
- Added **ordinal (cumulative) MAIHDA** for ordered-factor outcomes –
  the multicategorical extension the MAIHDA tutorials flag as priority
  future work. `family = "ordinal"` (alias `"cumulative"`; probit via
  the new exported `maihda_cumulative("probit")` or
  `brms::cumulative("probit")`) fits a cumulative (proportional-odds)
  model: a new **`engine = "ordinal"`** wraps
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) for the
  frequentist path (selected automatically for an ordinal family under
  the default engine, with a message) and `engine = "brms"` uses
  [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html).
  An ordered-factor outcome with 3+ levels under the all-default
  family/engine selects the model automatically (with a warning; a
  2-level ordered factor still takes the binomial auto-detect path), and
  an unordered factor response is coerced to ordered in its declared
  level order with a message. The VPC/ICC lives on the latent scale –
  level-1 variance pi^2/3 (logit) or 1 (probit), the same latent
  treatment as binomial models, so cumulative VPCs are comparable to
  logistic ones – and [`summary()`](https://rdrr.io/r/base/summary.html)
  gains a `thresholds` slot (the cut points that take the intercept’s
  place, with Hessian SEs) printed alongside the location fixed effects.
  Because `predict.clmm` does not exist, predictions are built from the
  stored coefficients: the link scale is the latent location eta =
  x’beta + u, and the response scale is the **expected category score**
  sum_k k\*P(Y = k) (the quantity the plots label “Average Expected
  Category Score”), with P(Y \<= k) = g^-1(alpha_k - eta) differenced by
  pure, unit-tested helpers shared with the brms path (whose
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) probability
  array
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md)
  now collapses to the same score). The two-model
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  decomposition and PCV,
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  (engine handshakes mirror the `sampling_weights` precedent), the
  stratum plots, `obs_vs_shrunken` (observed mean category score),
  `risk_vs_effect`, `effect_decomp` (exact on the latent scale), and the
  deviation panels all work;
  [`maihda_mor()`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md)
  now also accepts cumulative-logit fits, returning the **median
  cumulative odds ratio**. Explicit limits with targeted errors:
  logit/probit links only, the canonical single `(1 | stratum)`
  structure (no `context`/crossed-dimensions – use brms), no
  `sampling_weights` on the clmm path, no parametric bootstrap (clmm has
  no simulate/refit; brms provides credible intervals),
  AUC/[`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  stay binomial-only, and the ternary diagnostic is not yet supported.
  `ordinal` joins `Suggests`.
- Added the **negative-binomial family** for overdispersed count
  outcomes: `family = "negbinomial"` on
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  (and therefore
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)).
  The dispersion parameter theta is *estimated* from the data –
  [`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html) under
  the default engine, the `shape` parameter under `engine = "brms"` –
  and a fixed-theta `MASS::negative.binomial(theta)` family object is
  also accepted with lme4, honouring the supplied theta. The VPC/ICC
  uses the lognormal latent-scale level-1 variance
  `log(1 + 1/mu + 1/theta)` (Nakagawa, Johnson & Schielzeth 2017, *J R
  Soc Interface*), the negative-binomial analogue of the Stryhn et
  al. (2006) Poisson approximation already used by the package (it
  reduces to it as theta grows); for brms the `shape` draws are
  propagated into the VPC credible interval. The two-model
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  decomposition,
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md),
  parametric-bootstrap intervals (via
  [`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html), which
  holds theta fixed at its estimate – the interval is conditional on
  theta, as documented), prediction, and the stratum plots (routed to
  the count branch) all work; the log link is required, and the wemix
  engine,
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md),
  and
  [`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  reject the family with targeted messages. Internally, engine-specific
  family labels are now canonicalised so the family/link comparability
  checks in
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  no longer depend on raw label strings: two fits whose theta is
  *estimated* (lme4’s theta-embedding `"Negative Binomial(<theta>)"`,
  brms’s `shape`) compare equal, since each label otherwise carried its
  own theta estimate and could never match. A *fixed*, user-specified
  theta (`MASS::negative.binomial(theta)`) is treated differently – it
  is part of the model specification, so it stays in the comparability
  key and two fits with different fixed thetas are (correctly) rejected
  as incomparable.
- Added **design-weighted MAIHDA** (survey/sampling weights) via a new
  `sampling_weights` argument on
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).
  Sampling weights are *not* the same thing as lme4’s `weights=`
  (precision weights, which rescale the residual variance), so feeding
  survey weights to lmer/glmer maximises the wrong objective and gives
  invalid population estimates; supplying `sampling_weights` with
  `engine = "lme4"` is therefore an **error**, and supplying it with the
  default engine switches (with a message) to the new
  **`engine = "wemix"`**: weighted pseudo-maximum-likelihood via
  [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
  (Rabe-Hesketh & Skrondal 2006), the estimator built for
  NAEP/PISA-style survey analysis. The individual weights enter at level
  1 unchanged and the level-2 (stratum) weights are 1, because
  intersectional strata are exhaustive population cells sampled with
  certainty. The wemix engine supports the canonical MAIHDA structure –
  `gaussian(identity)` or `binomial(logit)` with the single
  `(1 | stratum)` intercept – and flows through the whole toolkit:
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the
  design-weighted VPC/ICC (latent-scale pi^2/3 level-1 variance for
  logistic fits, matching the other engines) and design-consistent
  (sandwich) fixed-effect standard errors; stratum estimates carry
  analytic conditional SEs (the design-weighted analogue of lme4’s
  `condVar`, reducing to it at unit weights);
  [`predict_maihda()`](https://hdbt.github.io/MAIHDA/reference/predict_maihda.md),
  the stratum plots (aggregated with the sampling weights, so stratum
  summaries are population-representative),
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  (which now also refuses to compare fits with *different* sampling
  weights), the
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  two-model decomposition,
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  all work. For a binomial fit,
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
  computes the **design-weighted AUC** (each observation contributes its
  weight as case/control mass) and flags it in the print method.
  Alternatively `engine = "brms"` accepts `sampling_weights` as
  likelihood weights (`y | weights(w)`, normalized to mean 1), giving a
  *pseudo-posterior*: point estimates are design-consistent but credible
  intervals are not design-based (a message says so). Limitations are
  explicit rather than silent: no parametric bootstrap for wemix (a
  design-based interval would need replicate weights – a possible future
  extension), and no crossed random effects, so `context =` and
  `decomposition = "crossed-dimensions"` require lme4/brms. Rows with
  missing or non-positive weights are dropped with a warning. A
  unit-weight wemix fit reproduces the lme4 ML fit to numerical
  precision. `WeMix` joins `Suggests`.
- Added the **contextual cross-classified MAIHDA** – the
  “cross-classified MAIHDA” of the literature (e.g. patients
  cross-classified by intersectional stratum and *hospital*, or students
  by stratum and *school*) – via a new `context` argument on
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  and [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md).
  `context = "school"` (one or more column names) enters each context as
  a crossed random intercept alongside the intersectional stratum
  effect, `outcome ~ covars + (1 | stratum) + (1 | school)`, and the
  summaries then partition the unexplained variance into
  **between-stratum vs. between-context vs. residual**: the headline
  VPC/ICC becomes the between-stratum share *net of* the context, each
  context gets its own `Context: <name>` row in the variance-components
  table, and a new `$context` summary element carries the per-context
  variances and shares (with parametric-bootstrap intervals for lme4 via
  `bootstrap = TRUE`, and posterior credible intervals for brms). In
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) the
  context random intercept is carried by *both* the null and the
  adjusted model, so the PCV is computed with the context partialled
  out; `context` also composes with
  `decomposition = "crossed-dimensions"` (the context variance then
  enters the single fit and its VPC denominator). A new plot type
  `"context_vpc"` (on both `maihda_model` and `maihda_analysis`) bars
  the stratum vs. context variances with their shares, and
  `plot(type = "vpc")` renders the contextual split automatically.
  `context` cannot be combined with `group` (a *stratified* per-level
  comparison vs. a *joint* crossed model are different designs;
  supplying both errors), a context variable may not be a stratum
  dimension or appear as a fixed effect, and a context with few levels
  weakly identifies its variance (consider `engine = "brms"`). A
  manually written `... + (1 | school)` still fits and is summarised
  generically as “Other random effects”; only `context =` activates the
  labelled partition.
- **Renamed** the `decomposition` value `"cross-classified"` to
  **`"crossed-dimensions"`** in
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  (and in the Shiny app’s decomposition toggle), because that mode
  crosses the stratum *dimensions’* main effects as random intercepts –
  whereas “cross-classified MAIHDA” in the literature means the
  contextual stratum-by-place model now fitted via `context` (see
  above). The old value is accepted as a **deprecated alias** that warns
  and maps to `"crossed-dimensions"`. Note for code that inspects
  results: a `maihda_analysis` from this mode now has
  `mode = "crossed-dimensions"` (was `"cross-classified"`), and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  results carry `attr(, "decomposition") == "crossed-dimensions"`;
  printed output and plot titles say “crossed-dimensions” accordingly.
- Added [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
  a single high-level entry point that runs the standard two-model
  MAIHDA workflow in one call. It fits the **null** model (the formula
  you supply) and the **adjusted** model (the same model plus the
  additive main effects of the stratum dimensions), summarises the null
  model’s VPC/ICC, and reports the **PCV** – the proportional change in
  between-stratum variance from null to adjusted (the additive share of
  the intersectional inequality). When a `group` is supplied it also
  runs this decomposition within each group (the
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  result gains a `pcv` column). It returns one consistent
  `maihda_analysis` object with new
  `model_adjusted`/`summary_adjusted`/`pcv` slots (alongside the
  unchanged null-model `model`/`summary`), and `print`, `summary`, and
  `plot` methods;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) routes the
  VPC/shrinkage views to the null model and the
  additive-vs-intersectional views (`effect_decomp`, `risk_vs_effect`,
  `ternary`) to the adjusted model, and gains `type = "group_pcv"`.
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) is
  intrinsically a decomposition and has no single-model mode – use
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  for a single fit. A numeric dimension auto-binned for the strata
  enters the adjusted model as its reconstructed tertile factor (not a
  linear term). Because it cannot decompose without an intersection,
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) errors
  (rather than returning a half-result) when the stratum-defining
  variables are unidentifiable (e.g. a hand-built `stratum` column) or
  define only one dimension; the shorthand `(1 | var1:var2)` and
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  paths both record the dimensions and decompose normally.
- Added the `maihda_country_data` dataset (OECD PISA 2018, accessed via
  the `learningtower` package): 3,600 students across six countries with
  gender x socioeconomic-status strata and mathematics-achievement
  outcomes. It showcases
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  / `maihda(group = "country")`, since intersectional inequality
  (VPC/ICC) genuinely differs across the countries.
- Added the `maihda_sparse_data` dataset and a new vignette, **“Bayesian
  MAIHDA for sparse intersections”**, showing why `engine = "brms"` is
  the safer estimator when intersectional cells are small. The
  (simulated) data carry a *known* interaction – 40% of the
  between-stratum variance, on a Gaussian outcome `y` and a binary
  outcome `event` – across 36 strata with deliberately skewed sizes
  (median 6, half the cells below 5). Under that sparsity the
  maximum-likelihood (lme4) crossed-dimensions fit is singular and
  over-shrinks the interaction (to ~23% Gaussian, ~3% binary – i.e. a
  real interaction read as “fully additive”), with no uncertainty
  attached; a weakly-informative prior on the random-effect SDs
  regularises the variance off the boundary and returns a calibrated
  credible interval that covers the truth. The vignette also documents
  the precompute-and-cache pattern for shipping brms results in a build
  (Stan cannot run on CRAN’s / pkgdown’s builders).
- Added
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md),
  a one-call **publication-ready results export** that assembles the two
  canonical MAIHDA write-up deliverables (cf. Evans et al. 2024’s
  tutorial) from a fitted
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis (or a single
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  model): (a) a **model-results table** contrasting the null (Model 1)
  and adjusted (Model 2) fits – intercept, between-stratum variance and
  SD, VPC/ICC, the PCV, and, for a binary outcome, the AUC and Median
  Odds Ratio – and (b) a **ranked-strata table** ordering every
  intersectional stratum by its model-predicted outcome (their Table 4),
  with the predicted value’s conditional interval, the stratum size, and
  the stratum random effect. It introduces no new estimator – the
  model-results table reuses the quantities from
  [`summary()`](https://rdrr.io/r/base/summary.html) and the
  ranked-strata table reuses `plot(type = "predicted")`’s stratum
  predictions, so the table agrees exactly with
  [`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html).
  The `$models` data frame is numeric and export-ready (statistics in
  rows, an estimate + `*_lower`/`*_upper` interval columns per model:
  the VPC bootstrap/posterior interval and the bootstrap PCV interval
  are carried, other rows are point estimates), and the
  [`print()`](https://rdrr.io/r/base/print.html) method renders the
  familiar “estimate \[low, high\]” layout plus the top/bottom strata
  (`n_strata`). It adapts to every fit type: a crossed-dimensions
  analysis gets “Additive share”/“Interaction share” rows instead of the
  PCV, a contextual cross-classified analysis (`context =`) gets a
  “Context share (VPC)” row, an ordinal fit leaves the intercept row
  `NA` (its category thresholds are reported by
  [`summary()`](https://rdrr.io/r/base/summary.html), not the table),
  and `which = "adjusted"` ranks the strata by the adjusted rather than
  the null model. Works across the lme4, brms, wemix, and ordinal
  engines.
- [`summary()`](https://rdrr.io/r/base/summary.html) of a
  binomial/Bernoulli MAIHDA model – and therefore
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md), whose
  summaries it builds – now reports the **discriminatory accuracy**
  automatically: the AUC / C-statistic and Median Odds Ratio (the “DA”
  in MAIHDA), as a new `discriminatory_accuracy` slot shown by the
  [`print()`](https://rdrr.io/r/base/print.html) methods, so the
  strata’s discriminatory accuracy sits alongside the VPC without a
  separate call. For
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md), the
  headline [`print()`](https://rdrr.io/r/base/print.html) shows the null
  model’s AUC/MOR with the adjusted model’s AUC for comparison. The
  response-scale (probability-scale) VPC is attached on request via
  `summary(model, response_vpc = TRUE, seed = )` or
  `maihda(..., response_vpc = TRUE, seed = )` (it is simulation-based,
  hence opt-in and seeded). Both are computed from the already-fitted
  model with no refit, and are skipped for non-binomial outcomes and for
  the cross-classified fit (whose single-stratum between-variance the
  MOR/response-VPC require is not defined across crossed random
  effects). The standalone
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md),
  [`maihda_auc()`](https://hdbt.github.io/MAIHDA/reference/maihda_auc.md),
  [`maihda_mor()`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md),
  and
  [`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  are unchanged.
- Added
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md),
  which **surfaces relative-fit information criteria** for choosing
  between model *structures* (different covariate sets, strata
  definitions, or families) – the question the VPC/PCV do not answer. It
  reports `AIC`/`BIC` for the likelihood engines (lme4, ordinal `clmm`)
  and the Bayesian `WAIC`/`LOOIC` for brms, takes one or more
  `maihda_model`s (or a `maihda_analysis`, which expands into its null
  and adjusted rows), and adds a `delta` column (gap from the best model
  on the primary criterion). Crucially it handles the **REML pitfall**:
  `lmer` fits Gaussian models by REML, whose AIC/BIC are *not*
  comparable across models with different fixed effects (the canonical
  null-vs-adjusted MAIHDA case), so when more than one model is supplied
  any REML fit is refitted with maximum likelihood via
  [`lme4::refitML()`](https://rdrr.io/pkg/lme4/man/refitML.html) before
  the criteria are read – matching what
  [`anova()`](https://rdrr.io/r/stats/anova.html) does on `lme4` models
  – and the `estimator` column records it. Design-weighted (`wemix`)
  pseudo-likelihood fits report `NA` (no standard AIC/BIC is defined).
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  now appends these criteria alongside the VPC/ICC by default
  (`ic = TRUE`); set `ic = FALSE` for the lean VPC-only table.
- [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
  now also reports the **discriminatory-accuracy trajectory** for a
  binary (binomial/Bernoulli) outcome, alongside the
  between-stratum-variance / PCV trajectory it already produced – the
  stepwise discriminatory-accuracy analysis of Merlo et al. (2016). Each
  step gains an `AUC` (that model’s C-statistic; step 0 is the
  strata-only discriminatory accuracy), `Step_AUC` and `Total_AUC` (the
  **absolute** change in AUC – delta-AUC – versus the previous step and
  versus the null, in contrast to the **proportional**
  `Step_PCV`/`Total_PCV`), and `MOR` (the Median Odds Ratio, logit link
  only) column. No extra models are fit: the columns are read off each
  step’s already-fitted model via
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md),
  so a design-weighted stepwise (`sampling_weights`) reports the
  design-weighted AUC, and the columns are simply absent for non-binary
  outcomes (the gaussian/poisson/ordinal table is unchanged). A new
  [`print()`](https://rdrr.io/r/base/print.html) method for the
  `maihda_stepwise` table notes the proportional-vs-absolute
  distinction.
- Added
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md),
  a diagnostic that flags which strata carry a credibly non-zero
  **intersectional interaction** – the heart of “where is there genuine
  intersectionality”. It reads each stratum’s interaction BLUP (the
  stratum random effect of the **adjusted** / crossed-dimensions model,
  i.e. the departure from the additive main-effects prediction) and
  flags the strata whose effect is credibly different from zero,
  returning a ranked `maihda_interactions` table (flagged strata first,
  by interaction magnitude). For the frequentist engines
  (`lme4`/`wemix`/`ordinal`) it uses the BLUP’s conditional standard
  error with an optional multiplicity correction (`adjust`, default
  `"none"`; the docs steer to FDR `"BH"` for screening many strata, with
  the full set of `p.adjust` methods available for a reviewer who needs
  a specific one); for `brms` it uses the **exact posterior tail** – a
  credible interval and the probability of direction `pd = P(BLUP > 0)`
  – rather than a normal approximation, and `adjust` is inert (the
  Bayesian answer is multiplicity-free). Passing a
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis uses the right (adjusted) model automatically; passing a bare
  null model is caught with a warning, since its stratum effects mix the
  additive and interaction parts. The help explains why a correction is
  optional on already-shrunken BLUP estimates (Gelman, Hill & Yajima
  2012; Gelman & Carlin 2014).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) gains a
  `highlight_interactions` argument (on both a `maihda_model` and a
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis) that focuses and stars the
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)-flagged
  strata on the BLUP-based views (`effect_decomp`, `predicted`,
  `obs_vs_shrunken`). Pass `TRUE` (flags computed with defaults), a
  multiple-testing method such as `"BH"`, or a precomputed
  `maihda_interactions` object to reuse a specific
  `conf_level`/`adjust`; for an analysis the flags are computed once
  from the adjusted model and reused across views. In `effect_decomp`,
  labels follow the selected flagged set, so a BH screen labels only
  strata that survive the BH adjustment. `FALSE` (default) leaves every
  plot unchanged.

### Improvements

- [`maihda_mor()`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md)
  now requires the **logit** link, not merely the binomial family. The
  Median Odds Ratio is an odds-ratio-scale quantity derived from the
  logistic latent variance, so applying its
  `exp(sqrt(2 * V_A) * qnorm(0.75))` formula to a
  `binomial(link = "probit")` fit returned a number that was not on the
  model’s scale.
  [`maihda_mor()`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md)
  now errors for a non-logit binomial link, and
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
  reports the (link-agnostic) AUC with `mor = NA` for such fits,
  recording the link and explaining the `NA` in its
  [`print()`](https://rdrr.io/r/base/print.html) method.
- Clarified that
  [`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  collapses the fixed part to its mean linear predictor before
  simulating, so for an adjusted (covariate) model the response-scale
  VPC is a conditional-at-mean estimate (evaluated at the average
  covariate profile), not one marginalised over the covariate
  distribution. It is exact for the canonical null/strata-only model;
  the documentation now states this and recommends reading it from the
  null model.
- Updated the “MAIHDA for binary outcomes” vignette, which still claimed
  the package shipped no AUC/MOR helper and defined a local one. It now
  uses the exported
  [`maihda_discriminatory_accuracy()`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md),
  [`maihda_auc()`](https://hdbt.github.io/MAIHDA/reference/maihda_auc.md),
  and
  [`maihda_mor()`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md),
  and points to
  [`maihda_vpc_response()`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md)
  for the response-scale VPC.
- Added the missing `shinycssloaders` dependency to the README’s
  interactive-dashboard list (it is in `DESCRIPTION` Suggests and used
  by
  [`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md)).
- Plotting is now unified under the base
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic.
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  and `compute_maihda_ternary_data()` now return classed objects, so
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) dispatches
  automatically:
  - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
    [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
    result (was
    [`plot_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_comparison.md))
  - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
    [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
    result, with `type = "vpc"`/“components” (was
    [`plot_group_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_group_comparison.md))
  - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
    `compute_maihda_ternary_data()` result (was `plot_maihda_ternary()`)
- The old
  [`plot_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_comparison.md),
  [`plot_group_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_group_comparison.md),
  and `plot_maihda_ternary()` functions still work but are
  **deprecated** and emit a one-time warning pointing to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
- [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  now records lme4 fit-quality diagnostics (singular fit and convergence
  warnings) on the model object;
  [`print()`](https://rdrr.io/r/base/print.html) on the model and
  [`summary()`](https://rdrr.io/r/base/summary.html) surface a “Fit
  diagnostics” note so a boundary/non-converged fit – which makes the
  VPC/PCV unreliable – is no longer silent.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now raises a single aggregated warning naming any group whose lme4 fit
  was singular or failed to converge (per-group fits on small strata are
  where this is most likely), so an unreliable per-group VPC – often
  pinned at 0 by a boundary fit – is no longer silent.
- The `"risk_vs_effect"` plot no longer frames the outcome axis as
  “risk”. A higher predicted value is not universally “bad” (it depends
  on the outcome), so the axis and title now read as a neutral “mean
  predicted value/probability”, with a note that the direction depends
  on the outcome.
- Clarified the documentation of
  [`plot_prediction_deviation_panels()`](https://hdbt.github.io/MAIHDA/reference/plot_prediction_deviation_panels.md)
  to match the implementation: the labelled points use a per-type metric
  – deviation from the mean prediction for Gaussian/Poisson (and the
  ordinal expected-score mode), the absolute deviance residual for
  binomial, and surprise for the ordinal surprise mode – and are not
  statistical outliers or model-misfit “deviants”.
- Clarified that the per-group VPC/ICC in
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  is the between-stratum *share* of variance, which can differ across
  groups because of the residual variance as well as the between-stratum
  variance. The documentation now points to the `var_between` column for
  comparing the absolute amount of intersectional variation, and notes
  that overlap of separate per-group intervals is not a valid test of
  whether two groups’ VPCs differ.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  result gains `type = "between_variance"` (and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  `maihda(group = )` analysis gains `type = "group_between_variance"`),
  which bars the absolute between-stratum variance by group – the
  *magnitude* of intersectional variation that the VPC *share* cannot
  convey. All group plots now name any groups omitted because their VPC
  was not estimable instead of dropping them silently.
- Clarified the PCV documentation
  ([`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md),
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and the print method): the PCV is a model-dependent change in
  between-stratum variance and equals variance “explained” only when the
  second model nests the first; the stepwise PCV is order-dependent and
  not a variable’s unique contribution. The vignette and Shiny app no
  longer describe PCV as variance causally “explained” by main effects
  or treat a negative PCV as evidence of hidden structural inequality.
- Corrected the [`summary()`](https://rdrr.io/r/base/summary.html)
  VPC/ICC documentation: the denominator is the total unexplained
  variance, which includes the variance of any additional random effects
  (not just between-stratum + residual), and a note on the
  weighted-Gaussian level-1 variance was added.
- Documented which families the MAIHDA variance summaries support
  (`gaussian("identity")`, binomial/Bernoulli with logit/probit,
  `poisson("log")`); other families such as `Gamma(link = "log")` will
  fit but [`summary()`](https://rdrr.io/r/base/summary.html)/VPC/PCV
  will stop with a clear “not implemented” error rather than returning
  an invalid number.
- README clarifications: a note that the CRAN release can lag this
  repository (so the newest features may require the GitHub development
  version), and the interactive-dashboard dependency list now includes
  `future`, `promises`, and `haven` (for SPSS/Stata uploads).
- Completed the PCV wording cleanup in the remaining vignette and
  Shiny-app text (the prior pass missed several spots), and softened
  over-strong app labels: the quadrant plot is “Mean Prediction
  vs. Stratum Effect” (not “Risk vs. Intersectional Effect”), the
  cumulative-PCV chart is “Change in Between-Stratum Variance” (not
  “Variance Explained”), and “deviant strata” is now “most extreme
  strata”.
- Removed the stale checked-in rendered vignette HTML
  (`vignettes/*.html`); these are build artifacts generated from the
  `.Rmd` sources and had drifted out of sync with the corrected text.
  They are now git-ignored and added to `.Rbuildignore`, so a locally
  rendered HTML is never shipped in the tarball and R CMD build
  regenerates `inst/doc` from the `.Rmd`.
- Fixed an invalid documented URL flagged by `R CMD check --as-cran`:
  the `maihda_country_data` `@source` linked to
  `https://www.oecd.org/pisa/data/`, which returns HTTP 403. It now
  links to the CRAN page of the `learningtower` package (the
  reproducible access route used to build the dataset), keeping the OECD
  PISA 2018 attribution in the text.
- The `data-raw/maihda_health_data.R` regeneration script no longer
  calls `install.packages("NHANES")` as a side effect; like the
  country-data script it now stops with a clear message asking the
  developer to install the dependency.

### Performance

- [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  (and the prediction-time stratum lookup) now matches rows to strata
  with a single vectorised, collision-safe key match instead of an
  O(rows x strata x variables) row-by-row scan, so it scales to large
  survey datasets. Behaviour is unchanged, including the correct
  handling of values that contain the stratum-label separator.

### Bug Fixes

- [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  now apply their analytic-sample identity checks to **design-weighted
  (`wemix`) fits**. A `WeMixResults` object exposes no
  [`nobs()`](https://rdrr.io/r/stats/nobs.html)/[`model.frame()`](https://rdrr.io/r/stats/model.frame.html),
  so the row-count, row-identity and outcome-fingerprint guards
  previously fell through to a silent pass: two wemix fits sharing a
  formula, `n` and strata but fit to **completely different outcome
  values** compared as if they were the same analytic sample
  ([`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  returned a PCV,
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  reported VPCs, neither warned). The guards now fall back to the
  wrapper’s stored analytic `$data` (the rows the engine actually fit),
  so a mismatched outcome is caught –
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  errors and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  warns – exactly as for the lme4 path. Genuinely matched wemix fits are
  unaffected.
- The **longitudinal PCV baseline** is now the between-stratum variance
  at the *observed* baseline time (`ref_time = min(time)`), not the raw
  time-0 random-intercept variance.
  [`maihda_longitudinal_pcv()`](https://hdbt.github.io/MAIHDA/reference/maihda_longitudinal_pcv.md)
  (and the `maihda_long_pcv` print method) evaluated `Sigma[1, 1]`
  directly, which is the variance at `time = 0` – correct only when time
  is centred so the baseline is 0. For a model whose time does not start
  at zero (e.g. waves 10:12) this reported the PCV of an extrapolated
  time-0 variance, which is meaningless and can even be negative; the
  baseline PCV is now `a(t)'Sigma a(t)` evaluated at `ref_time`,
  matching how [`summary()`](https://rdrr.io/r/base/summary.html)
  reports the baseline VPC. The slope-variance PCV is unchanged (it is
  invariant to where time is zeroed).
- Cross-sectional, single-value-per-stratum summaries are now **refused
  for longitudinal (growth-curve) models** rather than silently
  returning a misleading scalar. A growth model’s stratum estimand is a
  *trajectory* (random intercept + slope), so the scalar BLUP that
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)’s
  ranked-strata table and
  `plot(type = "predicted" / "obs_vs_shrunken" / "risk_vs_effect" / "effect_decomp" / "prediction_deviation" / "ternary")`
  build (which adds only the random intercept and drops the slope)
  misrepresents it. These paths now error – or, for
  [`maihda_table()`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md),
  omit the ranking with an explanatory note – and point to the
  trajectory tools (`predict(type = "strata")`,
  `plot(type = "trajectories")`, `plot(type = "vpc_trajectory")`).
  [`summary()`](https://rdrr.io/r/base/summary.html)’s stratum-estimates
  print is relabelled “baseline (intercept) deviations” for a
  longitudinal fit, with the same pointer.
- [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  now reports the analytic sample size `n` for a **design-weighted
  (`wemix`) fit** instead of `NA`. `WeMixResults` has no
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) method, so the IC table
  fell back to `NA`; it now uses `nrow(model$data)` (the rows the fit
  used), matching the `nobs` that
  [`glance()`](https://generics.r-lib.org/reference/glance.html) already
  reports for the same model.
- **Corrected the PCV for Gaussian `lmer` models (the REML pitfall).**
  The proportional change in between-stratum variance compares the
  stratum variance of a *null* and an *adjusted* model that differ in
  their fixed effects (the adjusted model adds the dimensions’ additive
  main effects). `lmer` fits Gaussian models by REML, and a REML
  variance estimate is **not** comparable across models with different
  fixed effects, so the PCV was biased – it overstated the adjusted
  model’s residual between-stratum variance and therefore **understated
  the additive share (the PCV)**.
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  – and hence
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)’s PCV,
  [`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
  and the per-group PCV in
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  – now refits any REML `lmer` model with maximum likelihood
  ([`lme4::refitML()`](https://rdrr.io/pkg/lme4/man/refitML.html))
  before reading the variances (and before the parametric bootstrap, so
  the interval matches), exactly as
  [`maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
  already does for AIC/BIC and as
  [`anova()`](https://rdrr.io/r/stats/anova.html) does on `lme4` models.
  In simulations with a known 60% additive share (40% interaction), the
  reported PCV moved from ~50% (biased, REML) to ~58–60% (ML),
  recovering the truth. `glmer` (GLMM) fits, the brms/wemix/ordinal
  engines, single-model VPC/ICC summaries (which correctly stay REML,
  being comparison-free), and singular/boundary fits are unaffected. In
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md),
  `var_between_adjusted` is now reported as `var_between * (1 - pcv)` so
  it shares the REML scale of the VPC’s `var_between` and the table
  stays internally coherent.
- VPC/ICC is no longer reported for a Gaussian model fit with a
  non-identity link (e.g. `gaussian(link = "log")`). The residual
  variance is on the response scale while the between-stratum variance
  is on the link scale, so
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md),
  and
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  now raise a clear error instead of silently returning an invalid
  variance partition.
- Binary-outcome auto-detection now keys off the analytic model frame –
  after applying covariate transformations (e.g. `log(x)`), dropping
  rows with missing values, dropping rows with a missing prior weight,
  and applying any `subset=` (including negative/character row indices)
  – instead of the raw outcome column. An outcome that is only 0/1 once
  excluded rows are removed is now correctly fit with
  `family = "binomial"`.
- Data-masked engine arguments forwarded through `...` (e.g. `weights=`,
  `subset=`, `offset=`) now work at any nesting depth, including through
  the [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  and
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  wrappers. Arguments are captured as quosures and resolved with the
  data mask
  ([`rlang::eval_tidy`](https://rlang.r-lib.org/reference/eval_tidy.html)),
  fixing the previous “object not found” / “..1 used in an incorrect
  context” failures (a direct
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  call worked, but the wrappers did not).
- A binary outcome is now recoded to 0/1 in a way that no longer breaks
  a `subset=` expression referencing the original response labels
  (e.g. `subset = y %in% c("no", "yes")`): the subset is evaluated
  against the original response before recoding.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now slices forwarded `weights=`/`subset=`/`offset=` to each group’s
  rows before fitting, not just for the row-count guard. An external
  (non-column) `weights` vector previously failed every group with a
  length mismatch, and an external `subset` vector could be recycled
  onto the wrong rows of later groups; both now align correctly per
  group.
- The Gaussian VPC/ICC now accounts for prior `weights`. With weights
  the per-observation residual variance is `sigma^2 / w_i`, so the
  level-1 variance reported is the mean conditional residual variance
  `mean(sigma^2 / w_i)` rather than a single `sigma^2`; unweighted
  models are unchanged.
- [`plot_effect_decomposition()`](https://hdbt.github.io/MAIHDA/reference/plot_effect_decomposition.md)
  now uses the stratum random effect (BLUP) itself as the intersectional
  component instead of (full prediction - fixed prediction). With
  additional random effects such as `(1 | site)` the latter wrongly
  absorbed those effects into the stratum component.
- The stratum-level “surprise” in
  [`plot_prediction_deviation_panels()`](https://hdbt.github.io/MAIHDA/reference/plot_prediction_deviation_panels.md)
  (ordinal, surprise mode) is now the average per-observation surprise,
  `mean(-log(p))` (log loss), instead of `-log(mean(p))`, which could
  change stratum rankings.
- The Shiny app
  ([`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
  now also auto-detects a numeric 0/1 outcome and fits it as binomial
  under the default family, matching the core API, instead of silently
  fitting a linear probability model.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now warns when groups end up with different populated strata even
  under `shared_strata = TRUE`, since their VPCs are then estimated over
  different stratum support and are not strictly directly comparable.
- [`plot_prediction_deviation_panels()`](https://hdbt.github.io/MAIHDA/reference/plot_prediction_deviation_panels.md)
  now plots Poisson/count models on the response (expected-count) scale
  with count labels, rather than routing them through the Gaussian
  link-scale branch.
- `compare_maihda_groups(min_group_n = ...)` now guards the analytic
  sample size (the rows the model actually fits) rather than the raw
  group row count, so a group with enough raw rows but a tiny usable
  sample is skipped instead of being fit on a handful of observations.
- `predict_maihda(type = "strata")` now respects `newdata`: it returns
  only the strata present in `newdata` (and errors on a stratum the
  model never saw, as `type = "individual"` does) instead of always
  returning every training stratum. With `newdata = NULL` it still
  returns all strata.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now counts populated strata for the pre-fit guard on the analytic
  model frame, not the raw subgroup. A group with two raw strata but
  only one stratum left after missing-row removal is now cleanly skipped
  as VPC-undefined instead of failing during fitting with “grouping
  factors must have \> 1 sampled level”.
- `n_boot` for bootstrap intervals must now be at least 10 (the minimum
  number of successful refits an interval requires); an unusably small
  `n_boot` fails immediately with a clear message instead of only
  erroring after the bootstrap runs.
- [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  now detect differing prior `weights`: previously two models fit to the
  same rows/outcome/strata but with different weights compared as if
  equivalent (PCV returned silently, no warning), even though weights
  change the variance estimates.
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  now errors and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  warns; an unweighted fit and an explicit unit-weight fit are still
  treated as equal.
- When a two-level non-0/1 outcome is recoded to 0/1, the chosen mapping
  (which level becomes the modeled event = 1) is now reported via a
  [`message()`](https://rdrr.io/r/base/message.html) and stored on the
  model as `$response_recoding`. Previously the mapping followed
  alphabetical (character) or declared (factor) level order silently,
  with no signal at all when `family = "binomial"` was passed
  explicitly, so the modeled event could be inverted unnoticed. The 0/1
  assignment rule is unchanged (it matches base `glm`).
- Stratum-level plot aggregations (`plot(type = "predicted")`,
  `"risk_vs_effect"`, `"obs_vs_shrunken"`, `"effect_decomp"`) now honour
  prior `weights`, collapsing per-stratum predictions with a
  prior-weight-weighted mean (and weighting the reference lines by the
  summed weights). This makes the plots consistent with the weighted
  Gaussian VPC for survey/weighted fits; unweighted models are
  unaffected, and aggregated-binomial (`cbind`) fits, whose
  `weights(type = "prior")` are the trials, are left unweighted to avoid
  double-counting.
- The Shiny app’s “Compute Bootstrap CIs” control now actually produces
  the VPC/ICC bootstrap intervals it (and the vignette) advertise. The
  bootstrap was previously applied only to the PCV, so the headline
  VPC/ICC was shown as a point estimate with no interval. The interval
  is now computed in the background worker (keeping the UI responsive)
  and displayed alongside the VPC/ICC in the Model Summary and
  Interactive Explorer tabs; the PCV interval in the PCV Results tab is
  unchanged.
- The Shiny app
  ([`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
  no longer aborts the whole analysis when the baseline (null) model has
  zero or negative between-stratum variance (a singular /
  no-between-variation fit), which makes
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  error by design. The fitted model, VPC/ICC, summaries, and plots are
  now shown as usual and the PCV is reported as unavailable (with the
  underlying reason) instead of failing the entire workflow.
- [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  now accepts a bare family *function* (e.g. `family = stats::gaussian`)
  – one of the documented forms (“as in `fit_maihda`”), alongside a
  family name (`"gaussian"`) and a family object
  ([`gaussian()`](https://rdrr.io/r/stats/family.html)). The per-group
  fits already handled it; only the family metadata recorded on the
  result did not, so the call fit every group and then failed with
  “object of type ‘closure’ is not subsettable”. The family function is
  now resolved to a family object up front, exactly as
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  does.

### Diagnostics

- `brms` fits now record MCMC convergence diagnostics (maximum Rhat and
  the number of divergent transitions) alongside the engine, surfaced in
  the “Fit diagnostics” block of
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  so a non-converged or divergent Bayesian fit is no longer silent.
- Bootstrap VPC/PCV intervals now report Monte Carlo error:
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  print the number of successful bootstrap draws and the Monte Carlo
  standard error of the bootstrap mean (`sd(draws)/sqrt(n)`), a coarse
  gauge of how much the bootstrap distribution’s centre would move with
  a different seed. (It is not the sampling uncertainty of the
  percentile interval’s endpoints, which is a larger order-statistic
  quantity.)

## MAIHDA 0.1.10

### New Features

- Added
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  to compare intersectional inequality (VPC/ICC and
  between-/within-stratum variance) across levels of a higher-level
  grouping variable such as country, region, or survey wave. It fits a
  stratified MAIHDA model per group, by default using shared/global
  strata so VPCs are directly comparable, with optional per-group
  bootstrap confidence intervals.
- Added
  [`plot_group_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_group_comparison.md)
  to visualize the result either as a VPC-by-group forest plot or as
  stacked variance-composition bars.

### Bug Fixes

- Fixed parametric-bootstrap confidence intervals for VPC
  ([`summary()`](https://rdrr.io/r/base/summary.html)) and PVC
  ([`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)):
  failed `refit()` iterations were silently recorded as `0` instead of
  being dropped, biasing intervals toward zero and suppressing the
  high-failure-rate warning. Failed iterations are now excluded
  correctly.
- Corrected the Poisson VPC residual-variance approximation to
  `log(1 + 1/mu)` (Stryhn et al. 2006); the previous `1/mu`
  linearization biased the VPC downward for low-count outcomes.
- [`plot_prediction_deviation_panels()`](https://hdbt.github.io/MAIHDA/reference/plot_prediction_deviation_panels.md)
  no longer draws zero-width “95% CI” bars when the underlying model
  does not supply standard errors
  (e.g. [`lme4::merMod`](https://rdrr.io/pkg/lme4/man/merMod-class.html));
  intervals are omitted instead of collapsed.

## MAIHDA 0.1.9

### Bug Fixes

- Clarified the Shiny dashboard PVC/HUD interpretation so negative PVC
  values are shown as variance unmasking rather than as unexplained
  interaction variance.
- Fixed the coverage workflow failure-artifact upload configuration.

## MAIHDA 0.1.8

CRAN release: 2026-05-16

### General Updates & New Features

- Added
  [`plot_prediction_deviation_panels()`](https://hdbt.github.io/MAIHDA/reference/plot_prediction_deviation_panels.md)
  function for visualizing predicted values and identifying deviant
  cases.
- Added `plot_risk_vs_effect()` function to create a quadrant
  scatterplot comparing overall marginal predicted risk against pure
  intersectional effects.
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
