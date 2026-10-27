# Run a Complete MAIHDA Analysis

A single high-level entry point that runs the standard two-model MAIHDA
workflow and returns one bundled object. It fits the **null** model
(covariates plus the intersectional random intercept, *excluding* the
stratum dimensions' main effects) and the **adjusted** model (the null
plus the additive main effects of the stratum-defining dimensions),
summarises the variance partition (VPC/ICC) of the null model, and
reports the **PCV** – the proportional change in between-stratum
variance from the null to the adjusted model, i.e. the additive share of
the intersectional inequality. When a higher-level grouping variable is
supplied it also compares this decomposition across that variable's
levels.

## Usage

``` r
maihda(
  formula,
  data,
  group = NULL,
  context = NULL,
  engine = "lme4",
  family = "gaussian",
  decomposition = c("two-model", "crossed-dimensions", "longitudinal"),
  autobin = TRUE,
  shared_strata = TRUE,
  min_group_n = 30,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95,
  response_vpc = FALSE,
  seed = NULL,
  sampling_weights = NULL,
  id = NULL,
  time = NULL,
  time_degree = 1,
  stratum_slope = TRUE,
  interactions = TRUE,
  estimation = c("fitted", "ML", "REML"),
  ...
)
```

## Arguments

- formula:

  A model formula, using either the intersectional shorthand
  `outcome ~ covars + (1 | var1:var2)` or `... + (1 | stratum)` when
  `data` already has a `stratum` column from
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).
  The dimensions' additive main effects may be listed (the
  fully-specified adjusted model) or omitted (added automatically, with
  a message); see Details. Any dimension main effect listed must be
  written as its bare column name: a transformed term such as
  `factor(edu)` is rejected, because only a bare name is recognised as a
  dimension's main effect and the transform would otherwise survive into
  the derived null model. Transform the column in `data` instead
  (`data$edu <- factor(data$edu)`). For the same reason a fixed
  interaction involving a dimension is rejected, whether with another
  dimension (`gender * race`) or with a covariate (`age * gender`): only
  the bare main effect is removed, so the interaction survives into the
  null model, which would then already adjust for that dimension.

- data:

  A data frame with the model variables (and the `group` variable, if
  used).

- group:

  Optional character string naming a higher-level grouping variable
  (e.g. `"country"`). When supplied,
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  is run and attached to the result. `group` runs a *stratified*
  analysis (one independent model per level); to instead model a higher
  level *jointly*, crossed with the strata, use `context`. The two are
  different designs but they *compose*: supplying both runs the
  stratified comparison where each per-group fit is itself a contextual
  cross-classified model (see `context`).

- context:

  Optional character vector naming higher-level *context* column(s) in
  `data` (e.g. `"school"`, `"hospital"`, `"region"`), forwarded to
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  Each context enters every fitted model as a crossed random intercept –
  `outcome ~ covars + (1 | stratum) + (1 | context)` – the *contextual
  cross-classified MAIHDA* of the literature. The summaries then
  partition the unexplained variance into between-stratum vs.
  between-context vs. residual, the headline VPC/ICC becomes the
  between-stratum share *net of* the context, and the PCV decomposition
  is computed with the context partialled out (the context random
  intercept is carried by both the null and the adjusted model). May be
  combined with `group`: the stratified comparison then fits a
  contextual cross-classified model within each group and reports each
  group's stratum-vs-context partition (a per-group subset shrinks the
  context level count, so weak identification is warned about per
  group). Not available for the `wemix`/`ordinal` engines (they fit no
  crossed random effect). A context with few levels identifies its
  variance weakly with either engine: `lme4` estimates it imprecisely,
  at zero (a singular fit) when the levels differ little, and with
  `engine = "brms"` its posterior depends on the prior on its standard
  deviation – with two levels roughly in proportion to the prior's scale
  – and under the default prior it comes out high when the context
  effects are smaller than that prior's scale, most with two or three
  levels; the context share and the VPC depend on it.
  [`print()`](https://rdrr.io/r/base/print.html) of the analysis and of
  its [`summary()`](https://rdrr.io/r/base/summary.html) names the
  contexts of a `brms` fit that have fewer than 10 levels, and the
  contextual partition of a `brms` summary holds the level counts as
  `context_levels`. A `brms` fit also reports the realised partition:
  the variances and shares of the effects of the strata and context
  levels in the data, which depend little on that prior when each
  context and stratum has enough data (see the section Realised strata
  and contexts). Prefer contexts with many levels (see `context` in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  for what to do with few).

- engine:

  Modeling engine, "lme4" (default), "brms", "wemix" (the
  design-weighted pseudo-maximum-likelihood fit; requires
  `sampling_weights` and is selected automatically when they are
  supplied with the default engine), or "ordinal" (cumulative link mixed
  model via
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html);
  selected automatically for an ordinal family or an ordered-factor
  outcome).

- family:

  Model family. Default "gaussian". As in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
  a binary outcome is auto-detected when `family` is left at the
  default, and the same resolved family is then used for the group
  comparison so all models agree.

- decomposition:

  How to decompose the intersectional inequality into additive and
  interaction parts. `"two-model"` (default) is the classic MAIHDA
  approach: a null model and an adjusted model (the dimensions' additive
  main effects as *fixed* effects), with the additive share read from
  the PCV (proportional change in between-stratum variance).
  `"crossed-dimensions"` fits a **single** model that enters each
  dimension's additive main effect as a *random* intercept –
  `outcome ~ covars + (1 | dim1) + ... + (1 | stratum)` – so each
  dimension's RE variance is its additive contribution and the
  intersection (`stratum`) RE variance is the interaction beyond
  additive; the additive and interaction *shares* of the total
  between-strata variance are read directly from that one fit. The two
  modes target the same scientific question with different estimators,
  so their additive shares are conceptually parallel but not numerically
  identical. The crossed-dimensions additive share is a partial-pooling
  estimate, and each dimension's variance is estimated from that
  dimension's levels alone, so a dimension with few levels (a binary sex
  variable has two) is poorly identified with either engine. `lme4`
  often returns a singular fit. With `engine = "brms"` the posterior of
  such a variance depends on the prior on its standard deviation – with
  two levels roughly in proportion to the prior's scale – and the VPC
  and the additive share follow it: under the default prior they are
  overstated, most with two- to four-level dimensions, and with two- or
  three-level dimensions the additive share comes out high even where
  the dimensions have no additive effect.
  [`print()`](https://rdrr.io/r/base/print.html) names the dimensions of
  a `brms` fit that have fewer than 10 levels. A `brms` fit also reports
  the realised strata: the VPC, shares and variances of the effects of
  the strata and dimension levels in the data, which depend little on
  that prior when each stratum has enough data, whatever the number of
  levels (see the section Realised strata). Otherwise prefer
  `"two-model"` (with few strata see `estimation`) unless every
  dimension has about ten or more levels. See Details.
  (`"cross-classified"` is accepted as a deprecated alias for
  `"crossed-dimensions"`, with a warning: in the MAIHDA literature
  “cross-classified” refers to the contextual stratum-by-place model,
  which this package fits via `context`.) `"longitudinal"` fits a
  3-level **growth-curve** MAIHDA (requires `id` and `time`; selected
  automatically when they are supplied): a null and an adjusted growth
  model, where the adjusted model adds the dimensions' main effects *and
  their interactions with time* (`dim:time`). The between-stratum
  variance is then time-varying and the PCV is reported separately for
  the baseline (intercept) and the slope variance – the
  additive-vs-multiplicative split of the intersectional trajectory
  inequality (Bell, Evans, Holman & Leckie 2024). As with the two-model
  PCV, the `estimation` argument sets the variance-estimation basis: the
  default `"fitted"` keeps each REML `lmer` growth fit's own variance,
  while `"ML"` refits them with maximum likelihood before the
  null-vs-adjusted comparison, which overstates the PCV when strata are
  few (see
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)).
  The reported time-varying VPC always keeps each fit's own (REML)
  estimate. See
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- autobin:

  Logical passed to
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md);
  tertile-bins numeric grouping variables. Default TRUE.

- shared_strata:

  Logical, forwarded to
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  when `group` is supplied: build strata once on the full data so VPCs
  are comparable across groups (TRUE, default) or rebuild them within
  each group.

- min_group_n:

  Minimum group size for the per-group comparison, forwarded to
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).
  Default 30.

- bootstrap:

  Logical; compute parametric-bootstrap VPC confidence intervals (lme4
  only) for both the overall summary and the per-group comparison.
  Default FALSE.

- n_boot:

  Number of bootstrap samples when `bootstrap = TRUE`.

- conf_level:

  Confidence level for bootstrap intervals. Default 0.95.

- response_vpc:

  Logical; for a binomial (lme4) outcome, also attach the response-scale
  VPC
  ([`maihda_vpc_response`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md))
  to the model summaries. It is estimated by simulation, so it is opt-in
  (default `FALSE`) and uses `seed`. For a Poisson or negative-binomial
  (lme4) outcome it attaches the count-scale VPC instead, a closed form
  (see
  [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)).
  The discriminatory accuracy (AUC + MOR) is attached automatically for
  a binomial/Bernoulli outcome regardless of this flag (see Details).

- seed:

  Optional integer seed for the response-scale VPC simulation. It is NOT
  forwarded to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html): an
  `engine = "brms"` fit's MCMC sampler is unseeded by it (a warning says
  so). Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before
  `maihda()` to make a Bayesian fit reproducible, or use
  `fit_maihda(seed = )`, which has no `seed` formal and so passes it
  through to `brm()`.

- sampling_weights:

  Optional name of a sampling-weight column for a **weighted MAIHDA** on
  survey data; see
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  for the full semantics (engine selection, the pseudo-likelihood
  weighting, and – importantly – why a single person weight does not
  deliver design-based inference for a multistage or stratified sample).
  Both the null and the adjusted model (and any per-group fits) use the
  same weights, so the PCV is a weighted decomposition. Not compatible
  with `engine = "lme4"`, `decomposition = "crossed-dimensions"` under
  the wemix engine, or `bootstrap = TRUE`.

- id, time, time_degree:

  For a **longitudinal** MAIHDA: the person/unit identifier column, the
  numeric measurement-time column, and the growth-curve polynomial
  degree (1 = linear). Supplying `id`/`time` selects
  `decomposition = "longitudinal"`. See
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  for the model structure; `group`, `context`, and `sampling_weights`
  are not supported alongside them. Default `NULL` (cross-sectional).

- stratum_slope:

  Longitudinal only: keep the stratum-level random slope(s) on `time`?
  `FALSE` fits a random intercept at the stratum level, so the
  between-stratum variance is constant over time and the decomposition
  reports `PCV_intercept` only. Applied to the null and adjusted growth
  models alike. See
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- interactions:

  Whether to compute the per-stratum interaction diagnostic
  ([`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md))
  on the adjusted / crossed-dimensions model and attach it as the
  `interactions` slot, surfaced in
  [`print()`](https://rdrr.io/r/base/print.html). `TRUE` (default) uses
  the diagnostic's own default correction (`adjust = "BH"`, FDR – since
  scanning all strata is a screening question); `FALSE` skips it; or
  pass a [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) method name
  – including `"none"` for the uncorrected, per-stratum
  individual-testing view. Uses `conf_level`. Not computed for a
  longitudinal decomposition. The computation is cheap (it reads the
  stratum estimates the summary already holds; no refit).

- estimation:

  Variance-estimation basis for the PCV (`"two-model"` and
  `"longitudinal"` decompositions), `"fitted"` (default), `"ML"` or
  `"REML"`; see
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).
  `"fitted"` differences each model's own between-stratum variance (the
  REML estimate for a Gaussian `lmer` fit, as
  [`summary()`](https://rdrr.io/r/base/summary.html) reports); `"ML"`
  refits REML `lmer` fits with maximum likelihood first. `"REML"`
  compares REML-type variances for a binomial or Poisson lme4 fit in the
  `"two-model"` decomposition, whose maximum-likelihood PCV runs high
  with few or small strata; it is an error with `group`, with `context`
  on such a fit, and in the `"longitudinal"` decomposition.

- ...:

  Additional arguments passed to
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  (and on to `lmer`/`glmer`). Includes `count_approximation`, which
  selects the log-link count level-1 variance approximation for the VPC
  and is applied to every derived fit.

## Value

An object of class `maihda_analysis`: a list with

- model:

  the fitted `maihda_model` (see
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md));
  the **null** model in `"two-model"` mode, or the single
  **crossed-dimensions** model in `"crossed-dimensions"` mode

- summary:

  the model's `maihda_summary` (VPC/ICC, variance components, stratum
  estimates; plus the additive/interaction `decomposition` in
  crossed-dimensions mode, and the stratum-vs-context `context`
  partition when `context` is supplied)

- model_adjusted:

  the fitted **adjusted** `maihda_model` (`"two-model"` and
  `"longitudinal"` modes; `NULL` otherwise)

- summary_adjusted:

  the adjusted model's `maihda_summary`, or `NULL`

- pcv:

  the proportional change in variance: the `pcv_result` from
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
  in `"two-model"` mode, or a `maihda_long_pcv` (the intercept/slope and
  time-specific PCV) in `"longitudinal"` mode; `NULL` otherwise. Both
  use the variance-estimation basis set by `estimation`: the default
  `"fitted"` keeps each Gaussian `lmer` fit's own REML variance
  (matching the single-model summaries), while `"ML"` refits the REML
  fits with maximum likelihood, which overstates the PCV when strata are
  few (see
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)).
  In `"longitudinal"` mode both PCVs, the slope's first, are understated
  on average when the strata differ little after adjustment next to the
  precision of one stratum's own estimate. For `engine = "brms"` in
  `"longitudinal"` mode the variances compared are posterior medians of
  the per-draw between-stratum variance – at the baseline, of the slope,
  and at each reporting time – and, like the two-model `brms` PCV, the
  result is understated on average when strata are few

- decomposition:

  the additive/interaction partition (additive and interaction variances
  and shares, per-dimension variances and, for a `brms` fit, each
  dimension's number of levels in `dim_levels` and the realised strata
  in `realised`, described in the section Realised strata;
  `"crossed-dimensions"` mode only, `NULL` otherwise)

- groups:

  a `maihda_group_comparison` when `group` is supplied, otherwise `NULL`

- interactions:

  the `maihda_interactions` diagnostic (per-stratum interaction BLUPs
  and flags) when `interactions` is not `FALSE`, otherwise `NULL`

- mode:

  `"two-model"`, `"crossed-dimensions"`, or `"longitudinal"`

- context_vars:

  the context variable name(s) when `context` was supplied, otherwise
  `NULL`

- formula, adjusted_formula, group_var, call:

  bookkeeping for printing

## Details

**Binomial companions.** For a binary outcome the model summaries also
carry the discriminatory accuracy (AUC / C-statistic and Median Odds
Ratio) – the "DA" in MAIHDA – automatically, so the null model's
strata-only AUC sits alongside its VPC; set `response_vpc = TRUE` to add
the (simulation-based) response-scale VPC as well. These are read from
`summary(x)` and the attached `summary_adjusted`, and the headline
[`print()`](https://rdrr.io/r/base/print.html) shows the null-model AUC.

This is a convenience wrapper around
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
and
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).
It *always* runs a decomposition, never a single plain fit: by default
the two-model null/adjusted pair, or – with
`decomposition = "crossed-dimensions"` – one crossed model whose
additive dimension random effects carry the additive share. For a single
undecomposed fit (e.g. just the null-model VPC / discriminatory
accuracy), call
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
directly.

**The dimensions' additive main effects.** You may write them in the
formula – the fully-specified, lme4-native adjusted model
`outcome ~ covars + var1 + var2 + (1 | var1:var2)` – or omit them.
Either way the null *excludes* the dimension main effects and the
adjusted *includes* them: when the formula already lists them it is
taken as the adjusted model and the null is derived by dropping them;
when they are missing `maihda()` adds them to the adjusted model and
emits a [`message()`](https://rdrr.io/r/base/message.html) so the
decomposition stays explicit. Only the *additive* main effects belong
here: a fixed interaction among the stratum dimensions – `var1 * var2`,
which R expands to `var1 + var2 + var1:var2` – duplicates the
intersectional stratum random intercept (it absorbs the between-stratum
variance into fixed cell means, which makes the PCV invalid), so it is
**rejected with an error**. Write `var1 + var2`; the intersection is
estimated by the stratum random effect (and quantified by
`decomposition = "crossed-dimensions"` or
[`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)).
The dimensions themselves are read from the random term: the shorthand
`(1 | var1:var2)` and
[`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
both record them, and a numeric dimension that
[`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
auto-binned enters the adjusted model as its reconstructed tertile
factor (matching the strata), not as a linear term. Because `maihda()`
is intrinsically a decomposition, it **errors** (rather than returning a
null-only result) when it cannot build the adjusted model – when the
dimensions cannot be recovered (a hand-built `stratum` column records
none) or there is only one dimension (no intersection to decompose). Use
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
for those single-model fits.

## Numeric stratum dimensions

The stratum-defining dimensions are categorical. A factor or character
dimension enters the adjusted model as a categorical main effect, and a
numeric dimension with more than 10 unique values is tertile-binned when
`autobin = TRUE` (see
[`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)).
A numeric dimension that is *not* auto-binned – a category code such as
`education = 1, 2, 3` (few unique values), or a many-valued numeric
fitted with `autobin = FALSE` – enters the adjusted model as a single
*linear* fixed effect rather than categorical main effects, which
changes the PCV/decomposition interpretation. `maihda()` therefore
**warns** in that case (for three or more distinct values – a binary
numeric dimension is exempt, since with two levels a linear term and a
factor are equivalent): wrap category codes in
[`factor()`](https://rdrr.io/r/base/factor.html) before fitting, or, for
a genuinely continuous variable, use it as a covariate rather than a
stratum dimension.

## Realised strata

For a `brms` crossed-dimensions fit, `decomposition$realised` and the
"Realised strata" block of
[`print()`](https://rdrr.io/r/base/print.html) summarise, over the
posterior draws, the effects of the \\J\\ strata in the data. In each
draw stratum \\j\\ has the effect \\\zeta_j\\, the sum of its
dimensions' level effects and its own (interaction) effect. These are
finite-population summaries (Gelman 2005): they describe the strata and
levels in the data, where the variance components describe the
distributions those were drawn from. A dimension with two levels says a
lot about the former and little about the latter, so the realised values
depend little on the prior on the standard deviations when each stratum
has enough data to pin down its own effect. With few observations per
stratum each draw also carries the uncertainty about every stratum's
effect: the realised variances then come out too high and move with the
prior. Compare the intervals of the stratum estimates with the spread
between the strata. Per draw:

- `between_var`:

  \\\sum_j (\zeta_j - \bar\zeta)^2 / (J - 1)\\, half the mean squared
  difference between two distinct strata.

- `vpc`:

  `between_var` over itself plus the within-stratum variance:
  \\\pi^2/3\\ for a logit link, 1 for a probit link, \\\pi^2/6\\ for a
  cloglog link, the residual variance of a Gaussian fit. Not given for a
  Poisson or negative-binomial fit, whose within-stratum variance is
  computed from the variance components. With a `context`, the realised
  variance of the context effects is added to the denominator (see the
  section Realised strata and contexts).

- `interaction_var`:

  the residual sum of squares of the least-squares fit of the
  \\\zeta_j\\ on an intercept and the dimensions' levels, divided by its
  residual degrees of freedom, `df_interaction`. Strata missing from the
  grid of dimension levels are allowed.

- `additive_share`:

  1 minus `interaction_var / between_var`, and `interaction_share` the
  rest: a mean-square split, like a PCV computed from REML variances. It
  is near zero when the dimensions explain none of the differences
  between the strata, and can then fall below zero.

- `per_dim`:

  the variance of each dimension's level effects, with its number of
  levels minus one as divisor.

- `mor`:

  for a logit link, the median odds ratio between two distinct strata of
  the data: the exponential of the median of \\\|\zeta_j -
  \zeta\_{j'}\|\\ over the pairs. With many strata it is computed from
  an evenly spaced subset of at least 1000 draws.

Each is the posterior median, with the equal-tailed credible interval at
`conf_level` in the matching `*_ci` element. Every stratum counts once,
whatever its size. `n_strata` is \\J\\. An `lme4` crossed-dimensions fit
reports its variance components only.

## Realised strata and contexts

For a `brms` fit with a `context`, `summary()$context$realised` and the
"Realised strata and contexts" block of
[`print()`](https://rdrr.io/r/base/print.html) give the realised
partition: over the posterior draws, the finite-population variances
(Gelman 2005) of the effects of the strata and of each context's levels
in the data. A context's variance component describes the distribution
its levels were drawn from, and with few levels it follows the prior on
its standard deviation; the variance of the effects of the levels in the
data describes those levels, which the data pin down when each level has
enough observations. With two contexts it is half the squared difference
between their effects. Per draw, for a grouping factor with \\C\\ levels
in the data and effects \\u_k\\:

- `var_stratum`:

  the realised variance of the strata: of the stratum intercepts, or of
  the \\\zeta_j\\ of a crossed-dimensions fit (its `between_var`).

- `per_context`:

  each context's \\\sum_k (u_k - \bar u)^2 / (C - 1)\\, and
  `context_var_total` their sum.

- `other_var`:

  the same, summed over any other random effect in the formula; given
  only when there is one.

- `vpc_stratum`:

  `var_stratum` over the sum of the realised variances and the
  within-stratum variance (see `vpc` in Realised strata), and
  `vpc_context` and `vpc_context_total` each context's and all contexts'
  share of that sum. Given for a Gaussian identity fit and for a
  binomial, Bernoulli or cumulative fit on its latent scale, not for a
  Poisson or negative-binomial fit. For a crossed-dimensions fit
  `vpc_stratum` is `decomposition$realised$vpc`.

Each is the posterior median, with the equal-tailed credible interval at
`conf_level` in the matching `*_ci` element; for `per_context` and
`vpc_context` a matrix with a row per context. Every level counts once,
whatever its size. `n_strata` is the number of strata, and the
partition's `context_levels` the number of each context's levels.
[`print()`](https://rdrr.io/r/base/print.html) of a two-model analysis
shows the realised VPC and context share of the null model. These values
depend little on the prior on the standard deviations when each context
and stratum has enough data. With few observations per context each draw
also carries the uncertainty about every context's effect, and they then
move with the prior. An `lme4` fit reports its variance components only.

## References

Gelman, A. (2005). Analysis of variance – why it is more important than
ever. *The Annals of Statistics*, 33(1), 1-53.
[doi:10.1214/009053604000001048](https://doi.org/10.1214/009053604000001048)

## See also

[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
for the single-model fitter,
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
for the group comparison, and
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
for the variance summary.

## Examples

``` r
# \donttest{
data(maihda_health_data)

# One call: null + adjusted fit, VPC summary, and PCV decomposition. Writing the
# dimensions' additive main effects (Gender + Race) gives the fully-specified
# adjusted model; maihda() derives the null by dropping them.
a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race), data = maihda_health_data)
a                              # VPC (null) and PCV (null -> adjusted)
#> MAIHDA Analysis
#> ===============
#> 
#> Null formula:    BMI ~ Age + (1 | stratum)
#> Adjusted formula:BMI ~ Age + Gender + Race + (1 | stratum)
#> Engine: lme4 | Family: gaussian
#> VPC/ICC (null): 0.0585
#> PCV (null -> adjusted): 0.4963
#> Between-stratum variance: 2.7383 (null) -> 1.3792 (adjusted)
#>   ~49.6% of the between-stratum variance is additive (the dimensions' main
#>   effects); the remainder is the between-stratum variance remaining after the
#>   additive main effects -- a model-dependent quantity
#> Strata: 10
#> Intersectional interactions: 4 of 10 strata flagged (95% interval, BH-adjusted)
#>   strongest: male × Black (-1.290, below)
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
#> 
a$pcv                          # proportional change in between-stratum variance
#> Proportional Change in Variance (PCV)
#> =====================================
#> 
#> PCV: 0.4963
#> 
#> Variance basis: as fitted (REML for Gaussian lmer, matching summary())
#> 
#> Between-stratum variance:
#>   Model 1: 2.738282
#>   Model 2: 1.379216
#>   Change:  1.359066 (49.63%)
#> 
#> Interpretation (PCV is the proportional change in between-stratum
#> variance between the models):
#>   Between-stratum variance is 49.6% lower in Model 2 than in Model 1.
a$formula                      # null:     BMI ~ Age + (1 | stratum)
#> BMI ~ Age + (1 | stratum)
#> <environment: 0x55abff27ef48>
a$adjusted_formula             # adjusted: null + Gender + Race main effects
#> BMI ~ Age + Gender + Race + (1 | stratum)
#> <environment: 0x55abffef3bb8>

# Omitting them is equivalent -- maihda() adds them to the adjusted model and
# emits a message; the null and PCV are identical to the explicit form above.
a0 <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data)
#> maihda(): added the additive main effect(s) of the stratum dimension(s) Gender, Race to the adjusted model; the null model excludes them. List them in the formula to specify the adjusted model explicitly.

plot(a, type = "vpc")          # null model

plot(a, type = "effect_decomp")# adjusted model (additive vs intersectional)


# Crossed-dimensions decomposition: one model, the dimensions' main effects entered
# as RANDOM intercepts. The additive and interaction shares of the between-strata
# variance are read directly from the single fit (no null/adjusted pair).
cc <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data,
             decomposition = "crossed-dimensions")
#> boundary (singular) fit: see help('isSingular')
cc                                    # VPC and additive/interaction shares
#> MAIHDA Analysis
#> ===============
#> 
#> Decomposition:   crossed-dimensions (single model)
#> Formula:         BMI ~ Age + (1 | Gender) + (1 | Race) + (1 | stratum)
#> Engine: lme4 | Family: gaussian
#> Fit diagnostics:
#>   Singular fit: the random-effects block for 'Gender' is at (or near) the boundary
#>     (a variance estimated at ~0, or a perfect correlation between an intercept and a slope).
#>     'Gender' is a stratum dimension of this crossed-dimensions model: the additive variance
#>     has nothing from it, so read the additive and interaction shares with caution.
#>     The between-stratum variance and the VPC come from the terms not at the boundary.
#> 
#> 
#> VPC/ICC: 0.0637
#> 
#> Additive vs. Intersectional Decomposition (crossed-dimensions):
#>   Additive (sum of dimension main effects) variance: 1.8415
#>   Intersectional interaction variance:               1.1593
#>   Total between-strata variance:                     3.0008
#>   Additive share of between-strata variance:    61.4%
#>   Interaction share of between-strata variance: 38.6%
#>   Per-dimension additive variance:
#>     Gender: 0.0000
#>     Race: 1.8415
#>   Note: the additive share is the crossed-dimensions analogue of the PCV but
#>   a different estimator; interpret the interaction share cautiously.
#> 
#> Strata: 10
#> Intersectional interactions: 0 of 10 strata flagged (95% interval, BH-adjusted)
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
#> 
cc$decomposition$additive_share       # crossed-dimensions analogue of the PCV
#> [1] 0.6136712
cc$formula                            # BMI ~ Age + (1|Gender) + (1|Race) + (1|stratum)
#> BMI ~ Age + (1 | Gender) + (1 | Race) + (1 | stratum)
#> <environment: 0x55abffd3a1b0>

# Add a higher-level grouping variable to also compare across its levels.
# maihda_country_data has a real country grouping (PISA achievement data):
data(maihda_country_data)
a2 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
             group = "country")
#> boundary (singular) fit: see help('isSingular')
#> maihda(): added the additive main effect(s) of the stratum dimension(s) gender, ses to the adjusted model; the null model excludes them. List them in the formula to specify the adjusted model explicitly.
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
a2
#> MAIHDA Analysis
#> ===============
#> 
#> Null formula:    math ~ (1 | stratum)
#> Adjusted formula:math ~ (1 | stratum) + gender + ses
#> Engine: lme4 | Family: gaussian
#> VPC/ICC (null): 0.1493
#> PCV (null -> adjusted): 1.0000
#> Between-stratum variance: 1352.2436 (null) -> 0.0000 (adjusted)
#>   ~100.0% of the between-stratum variance is additive (the dimensions' main
#>   effects); the remainder is the between-stratum variance remaining after the
#>   additive main effects -- a model-dependent quantity
#> 
#> Note: the adjusted model's between-stratum variance is at the singularity boundary,
#> so the PCV is pinned near 100% and its interval/SE are unreliable. A near-complete
#> PCV is consistent with genuinely additive strata (no interaction beyond the main
#> effects) as well as a degenerate fit; the two cannot be told apart from this fit.
#> 
#> Strata: 6
#> Intersectional interactions: 0 of 6 strata flagged (95% interval, BH-adjusted)
#> 
#> Group comparison by 'country':
#> MAIHDA Group Comparison
#> =======================
#> 
#> Group variable: country 
#> Engine: lme4  | Family: gaussian  | Strata: shared/global 
#> Variance basis: as fitted (REML for Gaussian lmer, matching summary()) 
#> 
#>           group   n n_strata     vpc var_between var_other var_residual    pcv
#>         Finland 600        6 0.10994       785.8         0         6361 1.0000
#>         Germany 600        6 0.14448      1271.6         0         7529 1.0000
#>           Italy 600        6 0.11890      1065.3         0         7895 1.0000
#>           Japan 600        6 0.13344      1032.3         0         6704 0.9266
#>          Mexico 600        6 0.13649       771.5         0         4881 1.0000
#>  United Kingdom 600        6 0.06011       470.5         0         7357 1.0000
#>  var_between_adjusted var_between_adjusted_ml pcv_status status
#>                  0.00                    0.00   singular     ok
#>                  0.00                    0.00   singular     ok
#>                  0.00                    0.00   singular     ok
#>                 75.78                   75.78         ok     ok
#>                  0.00                    0.00   singular     ok
#>                  0.00                    0.00   singular     ok
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
#> 
plot(a2, type = "group_vpc")

plot(a2, type = "group_pcv")


# Contextual cross-classified MAIHDA: instead of one model per country (group=),
# model the strata CROSSED with country in a single fit. The summary partitions
# the unexplained variance into between-stratum vs. between-country vs. residual,
# and the PCV is computed with country partialled out.
a3 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
             context = "country")
#> maihda(): added the additive main effect(s) of the stratum dimension(s) gender, ses to the adjusted model; the null model excludes them. List them in the formula to specify the adjusted model explicitly.
a3
#> MAIHDA Analysis
#> ===============
#> 
#> Null formula:    math ~ (1 | stratum) + (1 | country)
#> Adjusted formula:math ~ (1 | stratum) + (1 | country) + gender + ses
#> Engine: lme4 | Family: gaussian
#> Context: country (crossed contextual random intercept in the null and adjusted models)
#> VPC/ICC (null): 0.1032
#> Context share (null): 0.1283 (between-country share of unexplained variance)
#> PCV (null -> adjusted): 0.9978
#> Between-stratum variance: 915.2323 (null) -> 1.9877 (adjusted)
#>   ~99.8% of the between-stratum variance is additive (the dimensions' main
#>   effects); the remainder is the between-stratum variance remaining after the
#>   additive main effects -- a model-dependent quantity
#> Strata: 6
#> Intersectional interactions: 0 of 6 strata flagged (95% interval, BH-adjusted)
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
#> 
a3$summary$context$vpc_context_total  # the country (general contextual) share
#> [1] 0.1282643
plot(a3, type = "context_vpc")

# }
```
