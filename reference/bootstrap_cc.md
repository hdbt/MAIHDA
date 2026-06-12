# Bootstrap a crossed-dimensions MAIHDA partition (lme4)

Parametric bootstrap (simulate from the fitted model, refit) of the
crossed-dimensions VPC and the additive / interaction shares, returning
a percentile interval for each via `maihda_bootstrap_ci`. lme4 only –
brms returns posterior credible intervals directly. When `ctx_vars`
names contextual random intercepts, their variance enters each refit's
VPC denominator and a `context_vpc` interval (the contexts' total share)
is returned too.

## Usage

``` r
bootstrap_cc(model, cc, n_boot, conf_level, ctx_vars = character(0))
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

- ctx_vars:

  Character vector of contextual grouping factors (may be empty).

## Value

A list with `vpc`, `additive_share`, `interaction_share` (and
`context_vpc` when `ctx_vars` is non-empty), each a length-2 interval
carrying `n_ok`/`mc_se` attributes.
