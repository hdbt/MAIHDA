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
  interactions = FALSE,
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
  be updated.

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
  `glmer()`, honouring the supplied theta. `family = "ordinal"` (alias
  `"cumulative"`; or
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
  the fixed part). A context with few levels (say \< 10) weakly
  identifies its variance and often yields a singular lme4 fit; the
  `brms` engine handles this better. Writing the random effect directly
  in the formula (`... + (1 | school)`) fits the same model but is
  summarised generically as "Other random effects"; only `context =`
  activates the labelled contextual partition. Not supported by the
  `wemix` engine.

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
    random intercept. Fixed-effect standard errors are the sandwich
    (robust) errors WeMix reports, which account for the weighting and
    for dependence within the model's own grouping (the intersectional
    strata) – but not for clustering or stratification induced by the
    sample design, which the strata do not represent. The VPC/PCV are
    reported as point estimates (no bootstrap – see
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
  silently converge to a wrong solution. All results (the time-varying
  VPC, the PCV, plots, predictions) are reported on the original `time`
  scale; the column name `.maihda_ctime` is reserved for the internal
  centered variable. Default `NULL`.

- time_degree:

  Polynomial degree of the growth curve when `time` is supplied: 1
  (default) linear, 2 quadratic, etc. The brms engine supports degree 1
  only.

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

- ...:

  Additional arguments passed to `lmer`/`glmer` (lme4), `brm` (brms), or
  [`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
  (wemix; e.g. `nQuad`, `fast`). The lme4-style `weights` (precision
  weights), `subset`, and `offset` arguments are honoured only by the
  `lme4` engine, which applies them directly. The `wemix`, `ordinal`,
  and `brms` engines reject them: none takes them as a top-level fitting
  argument (brms in particular expects weighting/offset as formula
  addition terms, `weights(.)` / `offset(.)`, and design weights via
  `sampling_weights`). Prefilter `data` instead of using `subset` on
  those engines.

## Value

A maihda_model object containing:

- model:

  The fitted model object (lme4, brms, WeMix, or ordinal::clmm)

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
  singular fit / convergence for lme4 and WeMix, MCMC convergence
  (maximum Rhat, divergent transitions) for brms, and the optimizer
  convergence code for an ordinal (clmm) fit. For lme4 and clmm fits it
  also carries likelihood-adequacy caveats – count overdispersion and
  zero inflation, stratum random-effect non-normality, longitudinal
  residual autocorrelation, and an approximate ordinal proportional-odds
  screen – reported only when a conservative threshold is crossed, since
  the VPC/PCV and interaction estimates are conditional on the
  likelihood holding

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
#> Fixed Effects:
#>         term estimate
#>  (Intercept)    492.3
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
