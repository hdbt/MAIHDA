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
  conf_level = 0.95
)
```

## Arguments

- ...:

  Multiple maihda_model objects to compare.

- model_names:

  Optional character vector of names for the models.

- bootstrap:

  Logical indicating whether to compute bootstrap confidence intervals.
  Default is FALSE.

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default is 0.95.

## Value

A data frame comparing VPC/ICC across models with optional confidence
intervals.

## Details

VPCs are only directly comparable when the models share an outcome,
family/link, analytic sample, and strata – the canonical use is nested
models (e.g. null vs covariate-adjusted) on the *same* data and strata,
to show how the VPC attenuates. If the supplied models differ in any of
these, `compare_maihda()` still returns the table but issues a single
warning, because the VPCs are then not directly comparable.

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
