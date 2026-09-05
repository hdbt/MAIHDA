# Null-restricted parametric-bootstrap fixed effects for an lme4 fit

Internal helper. For each fixed-effect term, refits the model without
that term, simulates `n_boot` responses from the reduced fit, refits the
full model on each, and refers the observed Wald statistic to the
bootstrap distribution of \\\|t^\*\|\\ under a true null. Returns the
shape
[`maihda_fixed_effects_table`](https://hdbt.github.io/MAIHDA/reference/maihda_fixed_effects_table.md)
produces, with `df` `NA`: the reference is an empirical distribution,
not a \\t\\.

## Usage

``` r
maihda_bootstrap_fixef(model, n_boot, conf_level)
```

## Arguments

- model:

  An lme4 model object (`lmerMod` or `glmerMod`).

- n_boot:

  Number of bootstrap replicates *per term*.

- conf_level:

  Interval level.

## Value

A data frame with `term`, `estimate`, `se`, `statistic`, `df`,
`p_value`, `lower` and `upper`, carrying `n_boot_ok`,
`n_boot_nonconverged` and `interval_reliable` attributes – the last
`FALSE` when more than half the contributing draws failed to converge,
as for the VPC and PCV bootstraps (`maihda_report_nonconvergence`).

## Details

The intercept has no reduced model to simulate from – a MAIHDA intercept
is a reference-category level rather than a term that can be dropped –
so its p-value and interval are `NA`, not a Wald z that is miscalibrated
in exactly the regime this function exists for.
