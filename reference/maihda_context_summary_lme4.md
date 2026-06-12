# Contextual cross-classified variance summary (lme4)

Internal helper for
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
when the model carries a contextual random intercept
(`object$context_info` set, `fit_maihda(context = )`) without the
crossed-dimensions decomposition. Partitions the unexplained variance
into between-stratum vs. between-context (one share per context
variable) vs. residual. The headline VPC stays the between-stratum share
of all unexplained variance – numerically identical to the generic
single-stratum summary, which folds the context into "Other random
effects" – but the context is now named, given its own component row(s),
and returned as a `context` element.

## Usage

``` r
maihda_context_summary_lme4(object, ctx, vc, bootstrap, n_boot, conf_level)
```

## Arguments

- object:

  A `maihda_model` with `context_info`.

- ctx:

  The `context_info` list (`context_vars`).

- vc:

  The model's `VarCorr`.

- bootstrap, n_boot, conf_level:

  Bootstrap controls.

## Value

A list with `variance_components`, `vpc_result`, `context`.
