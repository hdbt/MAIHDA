# Cross-classified variance summary (lme4)

Internal helper for
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
when the model is a cross-classified MAIHDA fit (`object$cc_info` set).
Partitions the crossed random-effect variances into the additive (sum of
the dimension REs) and interaction (intersection RE) components, builds
the variance-components table and the VPC, and – when `bootstrap = TRUE`
– adds parametric-bootstrap intervals for the VPC and the
additive/interaction shares.

## Usage

``` r
maihda_cc_summary_lme4(object, cc, vc, bootstrap, n_boot, conf_level)
```

## Arguments

- object:

  A `maihda_model` (cross-classified).

- cc:

  The `cc_info` list (`dim_groups`, `interaction_group`).

- vc:

  The model's `VarCorr`.

- bootstrap, n_boot, conf_level:

  Bootstrap controls.

## Value

A list with `variance_components`, `vpc_result`, `decomposition`.
