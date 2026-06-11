# Bootstrap a cross-classified MAIHDA partition (lme4)

Parametric bootstrap (simulate from the fitted model, refit) of the
cross-classified VPC and the additive / interaction shares, returning a
percentile interval for each via `maihda_bootstrap_ci`. lme4 only – brms
returns posterior credible intervals directly.

## Usage

``` r
bootstrap_cc(model, cc, n_boot, conf_level)
```

## Arguments

- model:

  The underlying lme4 model object.

- cc:

  The `cc_info` list.

- n_boot:

  Number of bootstrap samples.

- conf_level:

  Confidence level.

## Value

A list with `vpc`, `additive_share`, `interaction_share`, each a
length-2 interval carrying `n_ok`/`mc_se` attributes.
