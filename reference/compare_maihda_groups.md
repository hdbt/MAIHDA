# Compare MAIHDA Metrics Across Levels of a Grouping Variable

Fits a separate random-intercept MAIHDA model (intercept-only *random*
effects; any fixed-effect covariates in `formula` are still used) within
each level of a higher-level grouping variable (for example country,
region, or survey wave) and reports how the variance partition
coefficient (VPC/ICC) and the between-/within-stratum variance
components differ across those groups. When the strata are defined by at
least two dimensions it also fits the adjusted model (the dimensions'
additive main effects) within each group and reports the per-group `pcv`
– the proportional change in between-stratum variance, i.e. the additive
share of that group's intersectional inequality.

## Usage

``` r
compare_maihda_groups(
  formula,
  data,
  group,
  engine = "lme4",
  family = "gaussian",
  shared_strata = TRUE,
  min_group_n = 30,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95,
  autobin = TRUE,
  decomposition = c("two-model", "crossed-dimensions"),
  context = NULL,
  sampling_weights = NULL,
  warn_linear = TRUE,
  estimation = c("fitted", "ML"),
  ...
)
```

## Arguments

- formula:

  A model formula. Either the shorthand intersectional form
  `outcome ~ covars + (1 | var1:var2)` (strata are built automatically)
  or `outcome ~ covars + (1 | stratum)` when `data` already contains a
  `stratum` column from
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).

- data:

  A data frame containing the variables in `formula` and the grouping
  variable.

- group:

  Character string naming the grouping variable in `data` (e.g.
  `"country"`). A separate model is fitted for each non-missing level.

- engine:

  Modeling engine, "lme4" (default), "brms", "wemix" (the
  design-weighted fit; requires `sampling_weights` and is selected
  automatically when they are supplied with the default engine), or
  "ordinal" (cumulative link mixed model via
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html);
  selected automatically for an ordinal family or an ordered-factor
  outcome).

- family:

  Model family. Default "gaussian". As in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
  a binary outcome is auto-detected once on the full data and switched
  to "binomial" (with a warning) so every group uses the same family.

- shared_strata:

  Logical. When TRUE (default) intersectional strata are defined once on
  the full data so that a stratum denotes the same combination in every
  group; this makes the stratum *definitions* comparable across groups.
  Note that a group may still not contain every stratum, so two groups'
  VPCs can be estimated over different sets of populated strata – they
  are then not strictly directly comparable, and the function warns when
  this happens. When FALSE, strata are rebuilt independently within each
  group (stratum identities are then not comparable across groups at
  all).

- min_group_n:

  Minimum size of the *analytic* sample a group must have – the rows
  that survive the model frame (covariate transformations applied, rows
  with a missing outcome/covariate dropped) – to be modelled. Groups
  with a smaller usable sample are skipped with a warning, even if they
  have more raw rows. Default 30.

- bootstrap:

  Logical; compute per-group parametric-bootstrap VPC confidence
  intervals. lme4 engine only. Default FALSE.

- n_boot:

  Number of bootstrap samples when `bootstrap = TRUE`. Default 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default 0.95.

- autobin:

  Logical passed to
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  controlling tertile binning of numeric grouping variables. Default
  TRUE.

- decomposition:

  Per-group additive-vs-interaction decomposition: the two-model null
  -\> adjusted PCV (`"two-model"`, default) or the single
  crossed-dimensions model (`"crossed-dimensions"`; `"cross-classified"`
  is a deprecated alias that warns). The crossed-dimensions form
  requires `shared_strata = TRUE` and at least two stratum dimensions,
  and adds the `var_additive`, `var_interaction`, `additive_share` and
  `interaction_share` columns (in place of `pcv` /
  `var_between_adjusted`); `var_between` is then the total
  between-strata variance (additive + interaction). See
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) for the
  underlying model and its caveats.

- context:

  Optional character vector naming higher-level *context* column(s) in
  `data` (e.g. `"school"`, `"region"`). Each per-group fit then becomes
  a contextual cross-classified model –
  `outcome ~ covars + (1 | stratum) + (1 | context)` – so within every
  group the stratum random intercept is crossed with the context random
  intercept(s). `vpc`/`var_between` are then the between-stratum
  quantities *net of* the context, and two columns report the per-group
  context partition: `var_context` (the between-context variance, summed
  over contexts) and `vpc_context` (the contexts' share of the group's
  unexplained variance). The per-context split is kept on the
  `"context_per"` attribute. A per-group subset shrinks each context's
  level count, so groups with too few context levels (\< 10) to identify
  the context variance are named in a single warning. Forwarded to
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md);
  not available for the `wemix`/`ordinal` engines (they fit no crossed
  random effect), and the `context` may not name the `group` variable
  itself.

- sampling_weights:

  Optional name of a sampling-weight column in `data` for
  design-weighted per-group fits; see
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  The column is sliced with each group's rows, so every group is fitted
  with its own members' weights. Not compatible with `engine = "lme4"`,
  `bootstrap = TRUE`, or (under the wemix engine)
  `decomposition = "crossed-dimensions"`.

- warn_linear:

  Logical; when `TRUE` (default) a one-time warning is emitted if a
  numeric stratum-defining dimension enters the per-group *adjusted*
  models as a raw linear slope rather than as categorical main effects,
  which silently changes the PCV/decomposition interpretation (see
  [`maihda_adjusted_terms`](https://hdbt.github.io/MAIHDA/reference/maihda_adjusted_terms.md)).
  Set to `FALSE` internally when
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  delegates here, because it has already warned on the same data.

- estimation:

  Variance-estimation basis for each group's per-group PCV, `"fitted"`
  (default) or `"ML"`; see
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).
  Affects Gaussian `lmer` fits only. Recorded as an `"estimation"`
  attribute on the returned comparison; the basis actually used is
  recorded as `"estimation_used"` (`"mixed"` when a group's ML refit was
  skipped at the boundary, leaving it on REML). Both are shown by
  [`print()`](https://rdrr.io/r/base/print.html).

- ...:

  Additional arguments passed to
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  (and on to `lmer`/`glmer`).

## Value

A `data.frame` of class `maihda_group_comparison` with one row per group
and columns `group`, `n`, `n_strata`, `vpc`, `var_between`, `var_other`,
`var_residual`, `status` (and `ci_lower`/`ci_upper` whenever a group's
summary supplies an interval: an lme4 bootstrap CI with
`bootstrap = TRUE`, or a brms posterior credible interval, which is
returned without bootstrapping). When the strata are defined by at least
two dimensions, two further columns report the per-group null -\>
adjusted decomposition: `pcv` (proportional change in between-stratum
variance when the dimensions' additive main effects are added; computed
on the variance basis set by `estimation` – see
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)),
`var_between_adjusted` (a *derived* coherence quantity, reported as
`var_between * (1 - pcv)` so it shares the scale of the REML
`var_between`/`vpc` and the table satisfies
`pcv = (var_between - var_between_adjusted) / var_between` exactly – it
is **not** the adjusted fit's own variance), and
`var_between_adjusted_ml` (the adjusted model's *actual* between-stratum
variance, read straight off the adjusted fit on the same basis as the
PCV; under the default `estimation = "fitted"` this is the REML
variance, coinciding with `var_between_adjusted`, and under
`estimation = "ML"` it is the ML variance, differing from
`var_between_adjusted` only by the small REML-vs-ML gap in the null
variance; the `_ml` suffix is retained for output-schema continuity).
All three are `NA` for a group whose adjusted fit failed. A fourth
column, `pcv_status`, records the decomposition outcome per group:
`"ok"` when the PCV was computed, `"failed"` when the adjusted model or
PCV errored (`pcv` is then `NA` and the group is named in a warning –
the group's own `status` can still be `"ok"` because its null VPC model
succeeded), and `"singular"` when the adjusted fit was singular and the
PCV saturated near 100% – *not* warned about, since a singular adjusted
fit is indistinguishable from genuinely additive strata (a legitimate,
common result), so it is flagged only in this column. These four columns
are omitted entirely when the strata have a single dimension. When
`context` is supplied, two further columns report each group's
contextual partition: `var_context` (the between-context variance,
summed over contexts) and `vpc_context` (the contexts' share of the
group's unexplained variance); the per-context split is on the
`"context_per"` attribute and the context name(s) on `"context_var"`.
These are dropped when no context is supplied. `n` is the analytic
sample size used by the model (after dropping rows with a missing
outcome/covariate) for both fitted and skipped groups, falling back to
the raw row count only when the model frame cannot be built. `var_other`
is the variance of any additional random effects and is 0 for the
canonical single-stratum model. Groups that were skipped or failed have
`NA` metrics and an explanatory `status`.

## Details

It estimates one VPC per group as a stratified analysis: each group is
modelled independently – the grouping variable itself is not modelled as
a crossed random effect, and the strata are not adjusted for it. (A
`context` argument, if supplied, does add a crossed contextual random
intercept *within* each group's model, making each per-group fit
contextual cross-classified; the groups themselves remain independently
fitted.)

The VPC is the *share* of the unexplained variance that lies between
strata, not the absolute magnitude of intersectional inequality. Because
it is a ratio, a group's VPC can differ from another's because the
between-stratum variance differs, because the within-stratum (residual)
variance differs, or both – two groups with the same between-stratum
variance can have very different VPCs. To compare the absolute amount of
between-stratum (intersectional) variation across groups, read the
returned `var_between` column alongside the VPC rather than treating a
higher VPC as "more inequality".

It is **descriptive**: it reports each group's VPC (with an interval
when available – an lme4 bootstrap CI or a brms credible interval) for
side-by-side comparison, but does not test whether the VPCs differ
between groups. The per-group intervals describe each group's own
uncertainty; whether two intervals overlap is *not* a valid test of the
difference between their VPCs, which would require modelling that
difference directly.

Robustness: a group whose *analytic* sample (rows surviving the model
frame) has fewer than `min_group_n` observations is always skipped with
a warning. A group with fewer than two populated strata is also skipped
(VPC is undefined with a single stratum) when the stratum membership is
known before fitting – that is, when `shared_strata = TRUE` or `data`
already carries a `stratum` column. Under `shared_strata = FALSE` strata
are rebuilt inside each group, so a degenerate single-stratum group is
instead reported with a "fit failed" status rather than a pre-fit skip.
A singular fit yields a VPC of 0 rather than an error (unlike
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)).
A hard fit failure in one group records `NA` and a status note without
aborting the whole comparison.

Fit-quality diagnostics: for the `lme4` engine, groups whose model is
singular or fails to converge keep a `status` of `"ok"` (the fit did
complete) but are named in a single aggregated warning, because their
VPC/ICC may be unreliable – a singular fit usually pins the
between-stratum variance at the boundary, giving a VPC of 0.

## See also

[`compare_maihda`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
for comparing different models on the same data;
[`plot.maihda_group_comparison`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_group_comparison.md)
for visualising the result.

## Examples

``` r
# \donttest{
data(maihda_country_data)
# How does gender x SES inequality in PISA math scores differ across countries?
cmp <- compare_maihda_groups(
  math ~ 1 + (1 | gender:ses),
  data = maihda_country_data,
  group = "country"
)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
print(cmp)
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
plot(cmp, type = "vpc")

# }
```
