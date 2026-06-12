# Bootstrap a contextual cross-classified MAIHDA partition (lme4)

Parametric bootstrap (simulate from the fitted model, refit) of the
between-stratum VPC and the contexts' total share for a contextual
cross-classified fit, returning a percentile interval for each via
`maihda_bootstrap_ci`. lme4 only – brms returns posterior credible
intervals directly.

## Usage

``` r
bootstrap_context(model, ctx_vars, n_boot, conf_level)
```

## Arguments

- model:

  The underlying lme4 model object.

- ctx_vars:

  Character vector of context grouping factors.

- n_boot:

  Number of bootstrap samples.

- conf_level:

  Confidence level.

## Value

A list with `vpc` (between-stratum share) and `context_vpc` (contexts'
total share), each a length-2 interval carrying `n_ok`/`mc_se`
attributes.
