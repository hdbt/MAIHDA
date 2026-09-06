# Null-restricted parametric-bootstrap fixed effects for an lme4 fit

Internal helper. For each fixed-effect term, refits the model with that
term's coefficients constrained to zero, simulates `n_boot` responses
from the restricted fit, refits the full model on each, and refers the
observed Wald statistic to the bootstrap distribution of \\\|t^\*\|\\
under a true null. Returns the shape
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

The restriction is imposed by dropping the term from the formula and
then *verifying* that the refitted design no longer spans the full
model's column space. It usually does not, but R's marginality rules
recode a surviving higher-order term to absorb a removed marginal one –
`. ~ . - x` applied to `y ~ x * f` gives `f + x:f`, whose full dummy
expansion spans exactly the original space – and such a "reduction"
constrains nothing. Where the spans agree the constraint is instead
imposed on the fitted design columns directly (`maihda_restrict_fixef`).

The intercept has no reduced model to simulate from – a MAIHDA intercept
is a reference-category level rather than a term that can be dropped –
so its p-value and interval are `NA`, not a Wald z that is miscalibrated
in exactly the regime this function exists for.
