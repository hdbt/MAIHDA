# Stepwise Proportional Change in Variance (PCV)

Estimates the proportional change in variance (PCV) sequentially by
fitting intermediate (partially-adjusted) models, adding each predictor
one-by-one. The step-specific PCV is the change in between-stratum
variance contributed by a predictor *given the variables already in the
model*. Because the steps are sequential it is order-dependent: it
reflects each variable's marginal, model-dependent change, not an
order-invariant “unique” contribution.

## Usage

``` r
stepwise_pcv(
  data,
  outcome,
  vars,
  engine = "lme4",
  family = "gaussian",
  context = NULL,
  sampling_weights = NULL,
  estimation = c("fitted", "ML")
)
```

## Arguments

- data:

  Data frame with observations. Ensure \`make_strata()\` was run first
  so the \`stratum\` variable exists.

- outcome:

  Character string; the dependent variable.

- vars:

  Character vector; predictors (strata groupings & covariates) to add
  sequentially to the model.

- engine:

  Modeling engine ("lme4", "brms", "wemix", or "ordinal"). Default is
  "lme4"; switches to "wemix" automatically when `sampling_weights` is
  supplied, and to "ordinal" for an ordinal family or ordered-factor
  outcome.

- family:

  Error distribution and link function. Default is "gaussian".

- context:

  Optional character vector naming one or more higher-level *context*
  columns in `data` (e.g. `"school"`, `"region"`), forwarded to every
  step's
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  call so that each model carries a crossed contextual random intercept,
  `(1 | context)`, alongside the intersectional stratum effect – the
  stepwise analogue of `maihda(context = )` (the literature's contextual
  cross-classified MAIHDA). The reported `Step_PCV` / `Total_PCV` are
  then the between-stratum PCV **net of** the context (the context
  intercept is held in the null model and every step, so the stratum
  variance is isolated from it), and an extra `Context_Variance` column
  reports the between-context variance held at each step. Context
  columns join the complete-case filter so every step shares one
  analytic sample. A context variable may not be a stratum dimension,
  `"stratum"` itself, or one of `vars`. Supported by the `lme4` and
  `brms` engines only (`wemix` and `ordinal` fit no crossed random
  effects). See
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  for the full contextual-model semantics. Default `NULL` (no context).

- sampling_weights:

  Optional name of a sampling-weight column for design-weighted stepwise
  fits; see
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  The weight column joins the complete-case filter so every step uses
  the same analytic sample.

- estimation:

  Variance-estimation basis for the between-stratum variances compared
  across steps, `"fitted"` (default) or `"ML"`; see
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).
  Affects Gaussian `lmer` fits only.

## Value

A data.frame (class `maihda_stepwise`) showing the sequential models,
the between-stratum variance at each step, and both the step-specific
and total PCV. For a **binary** (binomial/Bernoulli) outcome it also
carries the discriminatory-accuracy trajectory: `AUC` (the C-statistic
of each step's model – step 0 is the strata-only discriminatory
accuracy), `Step_AUC` and `Total_AUC` (the *absolute* change in AUC,
delta-AUC, versus the previous step and versus the null), and `MOR` (the
Median Odds Ratio, logit link only). These columns are absent for
non-binary outcomes. When `context` is supplied, a `Context_Variance`
column reports the between-context variance held at each step, and the
discriminatory-accuracy trajectory is omitted even for a binary outcome
– the AUC would include the context random effect and so mismatch the
net-of-context `Step_PCV` / `Total_PCV`, exactly as
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
drops it for a contextual fit.

## Details

All models are fit on the complete cases for \`outcome\`, \`stratum\`,
and all variables in \`vars\` so that each sequential variance
comparison uses the same analytic sample.

For a binary (or ordinal) outcome the sequential between-stratum
variances live on the latent scale, whose level-1 variance is fixed by
the link: a step that adds an *individual-level* variable (one varying
within strata) rescales the latent metric itself, so its `Step_PCV`
mixes explained variance with rescaling and can be understated or
negative on that account – see the “Latent-scale families and rescaling”
note in
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).
Steps that add a stratum-constant dimension are largely unaffected.

For a binary outcome the table additionally tracks discriminatory
accuracy (Merlo et al. 2016): `AUC` is each model's C-statistic and
`Step_AUC` / `Total_AUC` are its *absolute* change (delta-AUC), in
contrast to the *proportional* `Step_PCV` / `Total_PCV`. The `MOR` is
reported for the logit link (`NA` otherwise) and is a monotone transform
of the between-stratum variance already in `Variance`. For a
design-weighted fit (`sampling_weights`) the AUC is the design-weighted
(population) C-statistic. Reuses
[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
on each step's fitted model, so no additional models are fit. Note that
adding a *stratum-defining* dimension (one already encoded by the
strata) typically leaves the AUC essentially unchanged: it re-partitions
the between-stratum variance (so the PCV and MOR move) but not the
per-stratum predicted ranking the rank-based AUC depends on. The AUC
trajectory is therefore most informative for individual-level covariates
that vary *within* strata.

## References

Merlo, J., Wagner, P., Ghith, N., & Leckie, G. (2016). An original
stepwise multilevel logistic regression analysis of discriminatory
accuracy: the case of neighbourhoods and health. *PLOS ONE*, 11(4),
e0153778.

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
stepwise_pcv(strata_result$data, "health_outcome", c("gender", "race", "age"))
#>  Step      Model        Added_Variable Variance Step_PCV Total_PCV
#>     0 Null Model None (Intercept only)   26.715   0.0000    0.0000
#>     1    Model 1                gender   30.863  -0.1553   -0.1553
#>     2    Model 2                  race    2.346   0.9240    0.9122
#>     3    Model 3                   age    3.032  -0.2922    0.8865

# Contextual cross-classified stepwise PCV: strata crossed with a higher-level
# context (country). Step_PCV / Total_PCV are then net of the country intercept,
# and a Context_Variance column reports the between-country variance per step.
cc <- make_strata(maihda_country_data, c("gender", "ses"))
stepwise_pcv(cc$data, "math", c("gender", "ses"), context = "country")
#>  Step      Model        Added_Variable Variance Context_Variance Step_PCV
#>     0 Null Model None (Intercept only)  915.232             1137   0.0000
#>     1    Model 1                gender 1114.595             1136  -0.2178
#>     2    Model 2                   ses    1.988             1128   0.9982
#>  Total_PCV
#>     0.0000
#>    -0.2178
#>     0.9978
#> 
#> Step_PCV / Total_PCV are the between-stratum PCV NET OF the context random
#> intercept held in every model; Context_Variance is that (summed)
#> between-context variance at each step.
#> 
# }
```
