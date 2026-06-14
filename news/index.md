# Changelog

## MAIHDA 0.1.11

### New Features

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
  `predict_maihda(type = "strata")` returns each stratum’s intercept and
  slope (a stratum is now a *trajectory*). The intercept-only guards
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
  family labels are now canonicalised (notably lme4’s theta-embedding
  `"Negative Binomial(<theta>)"`), so the family/link comparability
  checks in
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  and
  [`compare_maihda()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
  no longer depend on raw label strings – previously two NB fits could
  never compare equal because each label carried its own theta estimate.
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
  the stratum random effect. It computes nothing new – every quantity is
  read from the summaries already attached to the analysis, so the table
  agrees exactly with
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
  “Context share (VPC)” row, an ordinal fit’s thresholds stand in for
  the intercept, and `which = "adjusted"` ranks the strata by the
  adjusted rather than the null model. Works across the lme4, brms,
  wemix, and ordinal engines.
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
  analysis) that rings and stars the
  [`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)-flagged
  strata on the BLUP-based views (`effect_decomp`, `predicted`,
  `obs_vs_shrunken`). Pass `TRUE` (flags computed with defaults) or a
  precomputed `maihda_interactions` object to reuse a specific
  `conf_level`/`adjust`; for an analysis the flags are computed once
  from the adjusted model and reused across views. `FALSE` (default)
  leaves every plot unchanged.

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
  and
  [`compute_maihda_ternary_data()`](https://hdbt.github.io/MAIHDA/reference/compute_maihda_ternary_data.md)
  now return classed objects, so
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
    [`compute_maihda_ternary_data()`](https://hdbt.github.io/MAIHDA/reference/compute_maihda_ternary_data.md)
    result (was
    [`plot_maihda_ternary()`](https://hdbt.github.io/MAIHDA/reference/plot_maihda_ternary.md))
- The old
  [`plot_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_comparison.md),
  [`plot_group_comparison()`](https://hdbt.github.io/MAIHDA/reference/plot_group_comparison.md),
  and
  [`plot_maihda_ternary()`](https://hdbt.github.io/MAIHDA/reference/plot_maihda_ternary.md)
  functions still work but are **deprecated** and emit a one-time
  warning pointing to
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
  convey. The VPC plot now also carries an interpretive caption (it is
  descriptive, and overlapping intervals are not a difference test), and
  all group plots now name any groups omitted because their VPC was not
  estimable instead of dropping them silently.
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
  standard error so the precision of an interval can be judged.

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
