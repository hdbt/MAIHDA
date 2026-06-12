# Contextual cross-classified variance summary (brms)

brms counterpart of
[`maihda_context_summary_lme4`](https://hdbt.github.io/MAIHDA/reference/maihda_context_summary_lme4.md):
computes the stratum / context / residual partition per posterior draw
and returns posterior point estimates with credible intervals for the
between-stratum VPC and the contexts' total share (no bootstrap – the
posterior supplies the interval).

## Usage

``` r
maihda_context_summary_brms(
  object,
  ctx,
  conf_level,
  point = c("median", "mean")
)
```

## Arguments

- object:

  A `maihda_model` with `context_info` (brms engine).

- ctx:

  The `context_info` list.

- conf_level:

  Credible-interval level.

- point:

  Posterior point estimate, "median" (default) or "mean".

## Value

A list with `variance_components`, `vpc_result`, `context`.
