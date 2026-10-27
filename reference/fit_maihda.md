# Fit MAIHDA Model

Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
Heterogeneity and Discriminatory Accuracy) using lme4, brms, WeMix (for
design-weighted (survey) data), or – for an ordered-factor outcome – a
cumulative link mixed model via
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html).

## Usage

``` r
fit_maihda(
  formula,
  data,
  engine = "lme4",
  family = "gaussian",
  autobin = TRUE,
  context = NULL,
  sampling_weights = NULL,
  id = NULL,
  time = NULL,
  time_degree = 1,
  stratum_slope = TRUE,
  interactions = FALSE,
  count_approximation = c("lognormal", "delta", "trigamma"),
  ...
)
```

## Arguments

- formula:

  A formula specifying the model. Can include a random effect for
  stratum (e.g., `outcome ~ fixed_vars + (1 | stratum)`) or can directly
  specify the intersection variables to be used for forming strata
  (e.g., `outcome ~ fixed_vars + (1 | var1:var2:var3)`). If variables
  other than "stratum" are provided in the random effect,
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  will be called internally to compute the strata and the formula will
  be updated. A random slope, such as `(1 + x | stratum)`, is fitted,
  but [`summary()`](https://rdrr.io/r/base/summary.html), the VPC, the
  PCV and the other MAIHDA summaries need intercept-only random effects
  outside a longitudinal model (see `id`). Fitted by `glmer()`
  (`glmer.nb()` for `family = "negbinomial"`) with `nAGQ` of 1 or more,
  a model with a random slope is refitted once when its fit stops with
  the variance of a block's first term, usually its intercept, at zero
  (in a block of three or more terms, also what is left of a later
  term's variance after the terms before it): glmer's first stage can
  leave it there, where the likelihood can still rise as that variance
  and the term's covariances with the later terms move off zero, and its
  second stage stops there too. The refit starts from the other side of
  that bound: with one later term, the sign of their covariance
  reversed; with two or more, along the direction in which the
  likelihood rises fastest. It is made by glmer's second stage with the
  `bobyqa` optimiser and replaces the fit when its -2 log-likelihood is
  lower by more than 1e-4. Such a model is then searched with each
  random-effect block's terms re-expressed as described for `lmer()`
  below and the fixed effects as the coefficients of centred and scaled
  columns, from lme4's starting values of the variance parameters with
  the fit's fixed effects and from the fit's own estimates, and refitted
  from the lower end by glmer's second stage with `bobyqa` when its -2
  log-likelihood is lower by more than 1e-4: with a slope covariate far
  from zero or in large units, such as a birth or calendar year, glmer's
  own second stage cannot resolve the slope's variance and can stop
  short of the maximum, often with what is left of that variance at
  zero, where lme4 calls the fit singular. Fitting then takes about one
  and a half times as long, and up to about four times. A fit refitted
  either way is refitted once more from its estimates by glmer's second
  stage with the optimiser your control names for that stage
  (`Nelder_Mead` by default; where Nelder-Mead's own first steps fail,
  or it stops without converging, its first steps are 0.1% of each
  estimate, at least 1e-6). This adds up to about three quarters to the
  time of a `glmer()` fit so refitted, and about doubles that of a
  `glmer.nb()` fit. Kept unless it ends higher than the refit or its
  optimiser reports a failure, the fit then records that optimiser, as
  glmer's own fits do;
  [`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html) of a fit
  without a control runs the optimiser the fit records, at its package's
  defaults. Fitted by `lmer()` outside a longitudinal model, a model
  with a random slope is refitted from the other side of such a bound
  when lmer's fit stops at one, and from a start carried back from the
  same model with each random-effect block's terms, other than an
  intercept, centred at their means and divided by their standard
  deviations in a block with an intercept, and divided by their root
  mean squares in a block without one, such as `(0 + x | ctx)` or the
  slope of `(1 + x || stratum)`. Besides stopping at the bound, lmer can
  stop short of the optimum when a slope covariate is far from zero or
  in large units, such as age in years or a calendar year, and where
  what is left of the variance of a block's last term after the terms
  before it is zero while the criterion still falls as that variance
  moves off zero; the fit kept and lmer's own fit are then refitted once
  more from that variance moved to the lowest criterion along it. A
  refit replaces lmer's fit when its REML criterion (deviance for a
  maximum-likelihood fit) is lower by more than 1e-4, the lowest of them
  kept; this takes about four times as long as one fit, and up to about
  nine times where lmer's own fit stops with a variance at zero. Pass
  `start` through `...` to fit from your own starting values only.

- data:

  A data frame containing the variables in the formula.

- engine:

  Character string specifying which engine to use: "lme4" (default),
  "brms", "wemix" (design-weighted pseudo-maximum-likelihood via
  [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html);
  requires `sampling_weights`), or "ordinal" (cumulative link mixed
  model via
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html);
  requires an ordinal family). When `sampling_weights` is supplied and
  `engine` is left at its default, the engine switches to "wemix"
  automatically (with a message); likewise an ordinal family (or an
  auto-detected ordered-factor outcome) switches the default engine to
  "ordinal".

- family:

  Character string, family object, or family function specifying the
  model family. Common options: "gaussian", "binomial", "poisson",
  "negbinomial". Default is "gaussian". `family = "negbinomial"` fits an
  overdispersed count model with the dispersion parameter theta
  *estimated* from the data: lme4 via
  [`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html) and
  brms via its `shape` parameter (log link only; not supported by the
  wemix engine). A fixed-theta `MASS::negative.binomial(theta)` family
  object is also accepted with `engine = "lme4"` and is fitted with
  `glmer()`, honouring the supplied theta. With lme4 an estimated theta
  needs `nAGQ` of 0 or 1: `glmer.nb()` mis-estimates theta under
  adaptive quadrature, so `nAGQ > 1` is an error there (a fixed theta is
  not affected). `family = "ordinal"` (alias `"cumulative"`; or
  [`maihda_cumulative`](https://hdbt.github.io/MAIHDA/reference/maihda_cumulative.md)`("probit")`
  /
  [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  for a non-logit link) fits a cumulative (proportional-odds) model for
  an *ordered-factor* outcome:
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) under
  the automatic "ordinal" engine,
  [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  under `engine = "brms"`. The VPC/ICC lives on the latent scale
  (level-1 variance \\\pi^2/3\\ logit / 1 probit, as for binomial
  models) and response-scale predictions are *expected category scores*
  (categories scored 1..K in order). An ordered-factor outcome with 3+
  levels under the default family selects this model automatically, with
  a warning. The logit and probit links are supported;
  `sampling_weights`, `context`, and lme4-style
  `weights`/`subset`/`offset` arguments are not available on the clmm
  path. If the outcome variable appears to be binary and the default
  family is used, the function will automatically switch to "binomial",
  recode two-level responses to 0/1 for `glmer()`, and issue a warning.
  When a two-level non-0/1 response is recoded (on either the
  auto-detected or an explicit `family = "binomial"` path), the mapping
  follows the usual convention – the first level becomes 0 (reference)
  and the second becomes 1 (the modeled event), where "first/second"
  means alphabetical order for a character outcome and the declared
  order for a factor. The chosen mapping is reported via a
  [`message()`](https://rdrr.io/r/base/message.html) and stored on the
  result as `$response_recoding`; set the factor levels (or supply a 0/1
  outcome) to control which level is the event. Although any valid
  family object is accepted for fitting, the MAIHDA variance summaries
  ([`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md),
  VPC/ICC, PCV) are only defined for `gaussian("identity")`, the
  binomial/Bernoulli families with a logit, probit, or complementary
  log-log (`cloglog`) link (latent level-1 variance \\\pi^2/3\\, 1, and
  \\\pi^2/6\\ respectively), `poisson("log")`, and the negative binomial
  with a log link (level-1 variance `log(1 + 1/lambda + 1/theta)` at the
  *marginal* expected count `lambda`; Nakagawa, Johnson & Schielzeth
  2017 – see
  [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
  for how `lambda` is computed). Other families (for example
  `Gamma(link = "log")`) will fit, but
  [`summary()`](https://rdrr.io/r/base/summary.html) and the VPC/PCV
  helpers will stop with an "not implemented" error because no level-1
  variance is defined for them.

- autobin:

  Logical indicating whether numeric variables used only for automatic
  strata creation should be binned by
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).
  Default is TRUE.

- context:

  Optional character vector naming one or more higher-level *context*
  columns in `data` (e.g. `"school"`, `"hospital"`, `"region"`). Each
  enters the model as a crossed intercept-only random effect alongside
  the intersectional stratum effect –
  `outcome ~ covars + (1 | stratum) + (1 | context)` – giving the
  *contextual cross-classified MAIHDA* of the literature (individuals
  cross-classified by stratum and place/institution).
  [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
  then partitions the unexplained variance into between-stratum vs.
  between-context vs. residual, and the headline VPC/ICC remains the
  between-stratum share (now net of the context). A context variable may
  not be a stratum dimension or `"stratum"` itself, and may not already
  appear as a fixed-effect term (its variance would then be absorbed by
  the fixed part). A context with few levels (say \< 10) identifies its
  variance weakly with either engine: `lme4` estimates it imprecisely,
  at zero (a singular fit) when the levels differ little, and with
  `engine = "brms"` its posterior depends on the prior on its standard
  deviation (under the default prior it comes out high when the context
  effects are smaller than that prior's scale, most with two or three
  levels), and the context share and the VPC depend on it. Prefer
  contexts with many levels; with few, report the context share with its
  interval (a `brms` summary also gives the realised partition, which
  depends little on that prior when each context and stratum has enough
  data; see the section Realised strata and contexts of
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)), or
  enter the context as a fixed effect instead
  (`outcome ~ context + (1 | ...)`, without `context =`) when only the
  strata net of it are wanted. Writing the random effect directly in the
  formula (`... + (1 | school)`) fits the same model but is summarised
  generically as "Other random effects"; only `context =` activates the
  labelled contextual partition. Not supported by the `wemix` engine.

- sampling_weights:

  Optional single character string naming a numeric column of `data`
  holding individual *sampling* (survey/design) weights, for a
  **weighted MAIHDA** on survey data. Sampling weights are not the same
  thing as lme4's `weights=` (precision weights, which rescale the
  residual variance), so combining `sampling_weights` with
  `engine = "lme4"` is an error.

  **What this does and does not buy you.** A single person-level weight
  column is all this argument accepts: there is no representation of the
  sampling hierarchy – no primary sampling units, no sampling strata, no
  higher-stage or conditional weights, no finite-population corrections,
  and no replicate weights. Weighting the likelihood makes the point
  estimates target the population-weighted estimand, which is the main
  reason to use it. The standard errors, however, are *not* general
  design-based standard errors for a multistage or stratified sample,
  because none of the design features that drive them are supplied.
  Treat the intervals as approximate unless your design genuinely
  matches the assumptions below; for full design-based inference on such
  data, use replicate weights with a survey-specific package.

  Two engines support them:

  - `engine = "wemix"` (chosen automatically when `engine` is left at
    its default): weighted pseudo-maximum-likelihood via
    [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
    (Rabe-Hesketh & Skrondal 2006), the estimator used for NAEP/PISA
    analysis. The individual weights enter at level 1 unchanged and the
    level-2 (stratum) weights are 1, because intersectional strata are
    treated as exhaustive population cells included with certainty – so
    this is a *single-stage* weighted model, and it is that assumption
    the inference rests on. Supports `gaussian(identity)` and
    `binomial(logit)` models with the canonical single `(1 | stratum)`
    random intercept and no
    [`offset()`](https://rdrr.io/r/stats/offset.html) term:
    [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
    leaves an offset out of the model it fits, so a formula with one is
    an error. For a Gaussian outcome, fit the response minus the offset
    instead (the same model) and add the offset back to its predictions.
    A `binomial(logit)` model whose between-stratum variance is
    estimated at zero is an error as well:
    [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
    cannot return that fit. It is a common result for the adjusted model
    when the strata are additive, and `engine = "brms"` returns a fit
    there. Fixed-effect standard errors are the sandwich (robust) errors
    WeMix reports, which account for the weighting and for dependence
    within the model's own grouping (the intersectional strata) – but
    not for clustering or stratification induced by the sample design,
    which the strata do not represent. The VPC/PCV are reported as point
    estimates (no bootstrap – see
    [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)).

  - `engine = "brms"`: the weights enter the model as likelihood weights
    (`y | weights(w)`), normalized to mean 1, giving a
    *pseudo-posterior*: point estimates target the population-weighted
    estimand but credible intervals are not design-based – interpret
    them cautiously. Full design-based (replicate-weight / linearised)
    variances for the variance components are not computed.

  Rows with a missing or non-positive sampling weight are dropped with a
  warning. The column names `.maihda_sw` (brms likelihood weight) and
  `.maihda_l2wt` (WeMix level-2 weight) are **reserved** for the
  design-weighted engines: they are written into the analytic data
  internally, so a model that supplies `sampling_weights` may not also
  reference a variable of either name in its formula (fitting would
  error). Default `NULL` (unweighted).

- id:

  Optional single character string naming a person/unit identifier
  column for a **longitudinal (growth-curve) MAIHDA** on long-format
  data (one row per measurement occasion). Id values must be *globally*
  unique to a person – ids numbered within a site or group (person "1"
  in every site) would merge different people's trajectories, and an id
  appearing in more than one stratum is rejected with an error. Supplied
  together with `time`, it makes the model a 3-level growth curve –
  occasions within individuals (`id`) within intersectional strata –
  with a random intercept and slope on `time` at *both* the individual
  and stratum levels. The growth random effects are added automatically:
  write the strata shorthand `(1 | var1:var2)` (or `(1 | stratum)`)
  only, not the slopes. The between-stratum variance (and hence the VPC)
  then becomes a function of time;
  [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
  reports the time-varying VPC. Longitudinal fits are supported by
  `engine = "lme4"`/`"brms"` only (not `wemix`/`ordinal`), and are
  incompatible with `context` and `sampling_weights`. Default `NULL`
  (cross-sectional). See Bell, Evans, Holman & Leckie (2024).

- time:

  Optional single character string naming a numeric measurement-time
  column (e.g. wave 0, 1, 2, ... or age), required for a longitudinal
  MAIHDA; see `id`. When the time axis does not start at 0 (age,
  calendar year, waves coded 10, 11, ...), the growth terms are fit on
  internally *centered* time (`time - min(time)`, with a message): the
  raw polynomial basis over an offset range is ill-conditioned and can
  silently converge to a wrong solution. With `engine = "lme4"` the
  model is also fit on time divided by the power of ten that puts its
  span between 1 and 10 (days over five years are divided by 1000; a
  span already between 1 and 10 is left as it is), with a message:
  counted in fine units (days, or months for a quadratic curve), lme4's
  optimisers stop far from the optimum, mostly without a convergence
  warning. A model fitted by `glmer()` (`glmer.nb()` for
  `family = "negbinomial"`) is refitted as well: glmer's first stage
  often leaves a growth block's intercept variance at zero, where the
  likelihood still rises with the sign of the intercept-slope covariance
  reversed, and its second stage stops there. Linear growth with a
  stratum slope is refitted from two more starts, which reset the
  person-level block and start the stratum block at half the
  person-level intercept SD and 0.71 times its slope SD, with an
  intercept-slope correlation of 0.71 or -0.71; a refit replaces the fit
  when its -2 log-likelihood is lower by more than 1e-4. Then any glmer
  growth fit that still has a block's intercept variance at zero is
  refitted from the other side of that bound by glmer's second stage
  with the `bobyqa` optimiser and, when that refit is kept, once more:
  across the bound again if it stops at one, otherwise from where it
  stops; each refit is kept on the same terms. Every glmer growth fit is
  then searched as a glmer model with a random slope is (see `formula`),
  with the penalised least squares inside each likelihood evaluation run
  to a relative tolerance of 1e-10 (lme4's is 1e-7), and refitted from
  the lower end when its -2 log-likelihood is lower by 0.001 or more:
  glmer's own second stage stops quadratic growth fits short of the
  maximum, often near a zero variance of the quadratic term, where lme4
  calls the fit singular. This takes about three times as long as one
  fit; the searches add about a fifth to that, and a refit from them
  about as much again. Pass `start` through `...` to fit from your own
  starting values only; they are on the internal time when `time` is
  scaled. All results (the time-varying VPC, the PCV, the fixed effects
  in [`summary()`](https://rdrr.io/r/base/summary.html), plots,
  predictions) are reported on the original `time` scale; the lme4 fit
  stored in the result's `model` element is on the internal time. The
  column name `.maihda_ctime` is reserved for the internal time
  variable. Default `NULL`.

- time_degree:

  Polynomial degree of the growth curve when `time` is supplied: 1
  (default) linear, 2 quadratic, etc. The brms engine supports degree 1
  only.

- stratum_slope:

  Longitudinal only: keep the stratum-level random slope(s) on `time`?
  `TRUE` (default) fits the canonical Bell et al. (2024) structure,
  `(time | id) + (time | stratum)`, in which the between-stratum
  variance is a function of time. `FALSE` fits
  `(time | id) + (1 | stratum)`: the individual level keeps its growth
  block, but strata differ in *level* only, so the between-stratum
  variance – and the numerator of the VPC – is constant over time and no
  `PCV_slope` is defined. It assumes that the strata do not diverge over
  time. A `(time | stratum)` block that collapses to a perfect
  intercept-slope correlation does not show that: few strata, few
  occasions per stratum or irregular measurement times produce such fits
  when the strata diverge too, and where they do, `FALSE` moves their
  slope differences into the person-level slope variance. The VPC still
  varies with time through the person-level slope variance and the
  residual, so this is a time-constant between-stratum *variance*, not a
  time-constant VPC. A Gaussian growth model fitted by `lmer()` with a
  stratum slope, or with quadratic or higher growth, is fitted from
  lme4's starting values and again from more starts, and a refit
  replaces lmer's fit when its REML criterion (deviance for a
  maximum-likelihood fit) is lower by more than 1e-4: lme4's default
  optimiser can stop at a boundary of a random-effect block away from
  the optimum. With linear growth the starts are scaled by the
  person-level block, with a positive and a negative stratum
  intercept-slope correlation; with quadratic or higher growth one start
  is lmer's fit of the same model on time centred at the middle of its
  range, and, when that refit and lmer's fit differ by more than 0.01,
  another has both covariance blocks diagonal. Pass `start` through
  `...` to fit from your own starting values only; they are on the
  internal time when `time` is scaled (see `time`).

- interactions:

  Opt-in per-stratum interaction diagnostic
  ([`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)),
  attached as the `interactions` slot and shown by
  [`print()`](https://rdrr.io/r/base/print.html). `FALSE` (default)
  skips it; `TRUE` computes it with the diagnostic's default correction
  (`adjust = "BH"`); or pass a
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) method name,
  including `"none"` for the uncorrected view. It is meaningful only on
  an *adjusted* model (the dimensions' main effects in the fixed part);
  on a null model `maihda_interactions` warns. This is the single-fit
  parallel to the default-on `interactions` of
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

- count_approximation:

  Which latent-scale level-1 (observation) variance approximation the
  VPC/ICC of a *log-link count* model uses, from table 1 of Nakagawa,
  Johnson & Schielzeth (2017): `"lognormal"` (default, \\\ln(1 +
  1/\lambda\\ \[+\\ 1/\theta\])\\), `"delta"` (\\1/\lambda\\ \[+\\
  1/\theta\]\\), or `"trigamma"` (\\\psi_1(\lambda)\\ for Poisson,
  \\\psi_1((1/\lambda + 1/\theta)^{-1})\\ for the negative binomial),
  with \\\lambda\\ the marginal expected count (see
  [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)).
  The three agree above about \\\lambda = 2\\ and diverge sharply below
  it, where the VPC depends heavily on the choice;
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the method
  and \\\lambda\\, and warns below that threshold. Nakagawa et al.
  recommend the trigamma form; at low counts it gives the smallest VPC
  of the three. The default is the approximation `performance::icc()`
  uses, but evaluated at a different \\\lambda\\ on an adjusted model,
  so the two VPCs agree exactly only for a null model. Used by
  [`summary()`](https://rdrr.io/r/base/summary.html), the bootstrap
  intervals and the longitudinal VPC(t); ignored for every other family.
  The count-scale VPC of Leckie et al. (2020), from
  [`maihda_vpc_response`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md),
  partitions the variance of the counts themselves and does not depend
  on this choice.

- ...:

  Additional arguments passed to `lmer`/`glmer` (lme4), `brm` (brms), or
  [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
  (wemix; e.g. `nQuad`, `fast`). The wemix engine rejects `mix()`'s
  `center_grand` and `center_group`, under any partial spelling: WeMix
  centres the covariates internally without keeping the centring
  constants, so predictions could not reproduce the fit. Centre
  covariates in `data` before fitting instead. It likewise rejects an
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the formula
  (see `sampling_weights`). The lme4-style `weights` (precision
  weights), `subset`, and `offset` arguments are honoured only by the
  `lme4` engine, which applies them directly. The `wemix`, `ordinal`,
  and `brms` engines reject them: none takes them as a top-level fitting
  argument (brms in particular expects weighting/offset as formula
  addition terms, `weights(.)` / `offset(.)`, and design weights via
  `sampling_weights`). Prefilter `data` instead of using `subset` on
  those engines. On `engine = "brms"` the response may carry the
  `trials()` and [`weights()`](https://rdrr.io/r/stats/weights.html)
  addition terms (`y | trials(n)`, `y | weights(w)`); any other –
  `se()`, `rate()`, `cens()`,
  [`trunc()`](https://rdrr.io/r/base/Round.html), `mi()`, ... – is
  refused, since the VPC, predictions and summaries do not model it.
  Write an exposure as `+ offset(log(expo))` rather than `rate(expo)`:
  for a Poisson model the same model, which they do (for a negative
  binomial, brms's `rate()` also scales the shape by the exposure, which
  the offset does not). An `engine = "brms"` fit without a `seed` takes
  its sampler seed from R's random-number stream, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  reproduces it, with the chains run in parallel (`cores > 1`) as well.

## Value

A maihda_model object containing:

- model:

  The fitted model object (lme4, brms, WeMix, or ordinal::clmm). An
  ordinal fit with a between-stratum variance of zero comes from
  `clmm()`'s general routine, `clmm.control(useMatrix = TRUE)`, the one
  that can return it. An ordinal model with fixed effects is fitted on
  its design with the columns centred and scaled (scaled only under
  `threshold = "symmetric2"`), since a covariate far from zero or on a
  large or tiny scale can stall the optimiser, and is returned in the
  formula's own terms with the gradient and Hessian of that fit, which
  give its standard errors; a fit given a `start` is made as the formula
  stands. An lme4 fit's call gives its control as a call to
  [`lme4::lmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html)
  or
  [`lme4::glmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html),
  so [`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html) and
  `confint(method = "boot")` refit it as they refit lme4's own fits. A
  glmer fit that `fit_maihda()` last refitted by bobyqa, or by
  Nelder-Mead with first steps of its own where Nelder-Mead's own fail,
  names that refit's control instead, bound in the formula's
  environment, and both stop on it: a draw refitted from the fit with
  those steps can stop short of its maximum. A control whose optimiser
  is a function is named the same way, as
  [`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html) does not
  take one

- engine:

  The engine used ("lme4", "brms", "wemix", or "ordinal")

- sampling_weights:

  The sampling-weight column name when supplied, NULL otherwise

- formula:

  The model formula

- data:

  The data used for fitting

- family:

  The family used

- strata_info:

  The strata information from make_strata() if available, NULL otherwise

- context_vars:

  The context variable name(s) when `context` was supplied, NULL
  otherwise

- interactions:

  The `maihda_interactions` diagnostic when `interactions` is not
  `FALSE`, NULL otherwise

- response_recoding:

  For a recoded two-level outcome, a data frame mapping each original
  level to its 0/1 value and role (reference/event); NULL when no
  recoding occurred

- diagnostics:

  Fit-quality diagnostics, surfaced by the print and summary methods:
  singular fit / convergence for lme4, WeMix and ordinal (clmm) fits,
  and MCMC convergence (maximum Rhat, divergent transitions) for brms.
  lme4 checks a fit in the units of its covariates (a `glmer()` fit in
  those of its fixed effects and random slopes, an `lmer()` fit in those
  of its random slopes), so its checks can fire at the maximum when such
  a covariate is far from zero or in large units, as a calendar year or
  an income is. When they fire and the optimiser reports success, the
  fit counts as converged if it lies within a tenth of a standard error
  of a stationary point of its likelihood (the restricted likelihood for
  an lmer fit by REML), measured with those covariates centred and
  scaled (a slope in a [`diag()`](https://rdrr.io/r/base/diag.html)
  block only scaled); the distance is `stationary_distance`, and lme4's
  messages are kept in `lme4_overruled` when the fit converged by it.
  lme4's singular-fit test, a diagonal entry of a random-effect block's
  relative Cholesky factor below 1e-4, is taken in the same units and
  can fire on a block with a slope on such a covariate that is not at
  the boundary: in an lmer fit, or a glmer fit with `nAGQ` of 1 or more,
  outside a longitudinal model, a block with a slope is at the boundary
  when a diagonal entry of its factor with the terms centred and scaled
  (in a block without an intercept, divided by their root mean square)
  is below 0.01, and lme4's own verdict is `lme4_singular` where the two
  differ. A fit with a [`diag()`](https://rdrr.io/r/base/diag.html),
  `cs()` or `ar1()` block keeps lme4's test. A glmer fit (`nAGQ` of 1 or
  more), or an lmer fit outside a longitudinal model, with the variance
  of a random-effect block's first term, usually its intercept, at zero
  (or, in a block of three or more terms, what is left of a later term's
  variance after the terms before it, such as a quadratic growth block's
  slope perfectly correlated with its intercept), where the likelihood
  still rises as the block's variances and covariances move, or as that
  variance and the term's covariances with the later terms move off
  zero, has stopped at a bound of lme4's parametrisation and counts as
  not converged; those blocks are in `boundary_traps`. A slope's entry
  can lie within lme4's 1e-4 of that bound at its natural size when its
  covariate is in large units, as a covariate in thousands is; such a
  column is not at a bound (`at_bound` is `FALSE`), and the fit counts
  as not converged when the likelihood still rises as the block's
  variances and covariances move. An lmer fit outside a longitudinal
  model that lme4 calls singular, which lme4 does not check further, or
  one with a random slope that lme4 does not check because of its size
  (by default one of 10,000 rows or more, or with 10 or more variance
  parameters), counts as not converged when a search of the same model,
  from lme4's starting values with each random-effect block's terms
  centred and scaled and from a variance at zero moved off it, reaches a
  REML criterion (deviance for a maximum-likelihood fit) 0.01 or more
  below the fit's; that amount is `search_gain`. A fit refitted as
  described under `formula` has been through those searches already, so
  this concerns a fit made from your own `start` or a singular one
  without a random slope. A glmer fit with `nAGQ` of 1 or more and a
  random slope, outside a longitudinal model and with no
  [`diag()`](https://rdrr.io/r/base/diag.html), `cs()` or `ar1()` block,
  that lme4 calls singular is searched from the fit with its
  random-effect terms and fixed effects centred and scaled, and counts
  as not converged when the search reaches a deviance 0.01 or more below
  the fit's (`search_gain`); such a fit refitted as described under
  `formula` has been through such searches already and is not searched
  again. A glmer fit with `nAGQ` of 1 or more and a random slope,
  outside a longitudinal model, that lme4 does not check because of its
  size (by default one of 10,000 rows or more, or with 20 or more
  parameters, fixed effects included) and does not call singular is
  measured as one whose checks fired, and counts as not converged, with
  `lme4_unchecked` TRUE, a tenth of a standard error or more from a
  stationary point; where that cannot be measured it is searched in the
  same way. An lme4 fit also carries likelihood-adequacy caveats – count
  overdispersion and zero inflation, stratum random-effect non-normality
  and longitudinal residual autocorrelation – reported only when a
  conservative threshold is crossed, since the VPC/PCV and interaction
  estimates are conditional on the likelihood holding. Every one of
  those checks is lme4-only, so a cumulative (clmm) fit raises no
  adequacy caveat at all: it stores a single descriptive fixed-only
  proportional-odds proxy statistic, never flagged because it cannot be
  separated from stratum heterogeneity; use
  [`maihda_proportional_odds_test`](https://hdbt.github.io/MAIHDA/reference/maihda_proportional_odds_test.md)
  to test that assumption

## Examples

``` r
# \donttest{
# Standard approach: manually create strata first
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race", "education"))
model <- fit_maihda(health_outcome ~ age + (1 | stratum),
                    data = strata_result$data,
                    engine = "lme4")

# Simplified approach: specify stratifying variables directly in the grouping structure
# The function internally calls make_strata() to create intersectionals
model2 <- fit_maihda(health_outcome ~ age + (1 | gender:race:education),
                     data = maihda_sim_data,
                     engine = "lme4")

# Contextual cross-classified MAIHDA: strata crossed with a higher-level
# context (here country) -- the literature's cross-classified MAIHDA.
data(maihda_country_data)
model3 <- fit_maihda(math ~ 1 + (1 | gender:ses),
                     data = maihda_country_data,
                     context = "country")
summary(model3)  # between-stratum vs. between-country vs. residual
#> MAIHDA Model Summary
#> ====================
#> 
#> Variance Partition Coefficient (VPC/ICC):
#>   Estimate: 0.1032
#> 
#> Variance Components:
#>                  component variance    sd proportion
#>   Between-stratum (random)    915.2 30.25     0.1032
#>           Context: country   1137.1 33.72     0.1283
#>  Within-stratum (residual)   6813.0 82.54     0.7685
#>                      Total   8865.3 94.16     1.0000
#> 
#> Contextual Cross-Classified Partition (stratum x context):
#>   Between-stratum (intersectional) variance: 915.2323 (share 10.3%)
#>   Context 'country' variance: 1137.1067 (share 12.8%)
#>   Note: the headline VPC/ICC is the between-stratum share conditional on
#>   the context random effect(s). The context share is the between-context
#>   component of the unexplained variance.
#> 
#> Fixed Effects (Wald t, 95% intervals):
#>         term estimate    se statistic df  p_value lower upper
#>  (Intercept)    492.3 18.55     26.54  5 1.42e-06 444.6 539.9
#>   df: containment (between-within).
#> 
#> Stratum Estimates (first 10):
#>  stratum stratum_id           label random_effect    se lower_95 upper_95
#>        1          1   male × Medium         7.404 9.689   -11.59   26.396
#>        2          2 female × Medium        -7.027 9.740   -26.12   12.063
#>        3          3   female × High        30.825 9.729    11.76   49.894
#>        4          4      male × Low       -28.600 9.741   -47.69   -9.508
#>        5          5    female × Low       -37.651 9.724   -56.71  -18.592
#>        6          6     male × High        35.048 9.713    16.01   54.085
# }
```
