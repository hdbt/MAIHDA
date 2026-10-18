# Crossed-dimensions variance summary (brms)

brms counterpart of
[`maihda_cc_summary_lme4`](https://hdbt.github.io/MAIHDA/reference/maihda_cc_summary_lme4.md):
computes the additive / interaction partition per posterior draw and
returns posterior point estimates with credible intervals for the VPC
and the shares (no bootstrap – the posterior already supplies the
interval). A contextual random intercept (`object$context_info` set)
enters the per-draw VPC denominator and is reported as its own component
row and `context` element. Each variance component is the posterior
median of its draws, and the additive, between-strata and context totals
are sums of those components. When the draws hold the level effects,
`decomposition$realised` carries the realised strata and
`context$realised` the realised partition (see the sections Realised
strata and Realised strata and contexts of
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)).

## Usage

``` r
maihda_cc_summary_brms(object, cc, conf_level, point = c("median", "mean"))
```

## Arguments

- object:

  A `maihda_model` (crossed-dimensions, brms engine).

- cc:

  The `cc_info` list.

- conf_level:

  Credible-interval level.

- point:

  Posterior point estimate of the VPC and the shares, "median" (default)
  or "mean".

## Value

A list with `variance_components`, `vpc_result`, `decomposition` and
`context` (NULL without a context).
