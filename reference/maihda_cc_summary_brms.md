# Cross-classified variance summary (brms)

brms counterpart of
[`maihda_cc_summary_lme4`](https://hdbt.github.io/MAIHDA/reference/maihda_cc_summary_lme4.md):
computes the additive / interaction partition per posterior draw and
returns posterior point estimates with credible intervals for the VPC
and the shares (no bootstrap – the posterior already supplies the
interval).

## Usage

``` r
maihda_cc_summary_brms(object, cc, conf_level, point = c("median", "mean"))
```

## Arguments

- object:

  A `maihda_model` (cross-classified, brms engine).

- cc:

  The `cc_info` list.

- conf_level:

  Credible-interval level.

- point:

  Posterior point estimate, "median" (default) or "mean".

## Value

A list with `variance_components`, `vpc_result`, `decomposition`.
