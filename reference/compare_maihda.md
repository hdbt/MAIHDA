# Compare MAIHDA Models

Compares variance partition coefficients (VPC/ICC) across multiple
MAIHDA models, with optional bootstrap confidence intervals.

## Usage

``` r
compare_maihda(
  ...,
  model_names = NULL,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95,
  ic = TRUE
)
```

## Arguments

- ...:

  Multiple maihda_model objects to compare.

- model_names:

  Optional character vector of names for the models.

- bootstrap:

  Logical; for **lme4** models, compute parametric-bootstrap VPC
  confidence intervals. Default FALSE. It does not apply to **brms**
  models, which always return a posterior credible interval (so passing
  `bootstrap = TRUE` with brms models errors) – their interval is
  included regardless.

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for the VPC interval (lme4 bootstrap CI or brms
  credible interval). Default is 0.95.

- ic:

  Logical; append relative-fit information criteria to the table for
  comparing model *structures*: `AIC`/`BIC` for the likelihood engines
  (lme4, ordinal) and `WAIC`/`LOOIC` for brms (see
  [`maihda_ic`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)).
  Default TRUE. REML `lmer` fits are refitted with ML so AIC/BIC are
  comparable across different fixed effects. Set FALSE for the lean
  VPC-only table.

## Value

A `maihda_comparison` data frame of VPC/ICC by model. Interval columns
(`ci_lower`/`ci_upper`) are included when any model supplies an interval
– an lme4 bootstrap CI or a brms posterior credible interval. When
`ic = TRUE`, information-criteria columns (`AIC`/`BIC` or
`WAIC`/`LOOIC`, whichever apply) are appended.

## Details

VPCs are only directly comparable when the models share an outcome,
family/link, analytic sample, and strata – the canonical use is nested
models (e.g. null vs covariate-adjusted) on the *same* data and strata,
to show how the VPC attenuates. If the supplied models differ in any of
these, `compare_maihda()` still returns the table but issues a single
warning, because the VPCs are then not directly comparable. The same
comparability caveat applies to the appended information criteria (see
[`maihda_ic`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)). In
addition, when the appended criteria mix scales – likelihood `AIC`/`BIC`
(lme4/ordinal) shown alongside Bayesian `WAIC`/`LOOIC` (brms), which can
happen for a same-family lme4-vs-brms comparison that the family/link
check does not flag – `compare_maihda()` warns, because those criteria
are on different scales and are not comparable to each other.

## Examples

``` r
# \donttest{
# Canonical use: nested models on the SAME data and strata (null vs adjusted)
strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))

null_model <- fit_maihda(health_outcome ~ 1 + (1 | stratum), data = strata$data)
adj_model  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata$data)

# Compare without bootstrap
comparison <- compare_maihda(null_model, adj_model,
                            model_names = c("Null", "Adjusted"))

# Compare with bootstrap CI
comparison_boot <- compare_maihda(null_model, adj_model,
                                 model_names = c("Null", "Adjusted"),
                                 bootstrap = TRUE, n_boot = 500)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
# }
```
