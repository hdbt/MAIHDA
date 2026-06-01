# Resolve shared strata and the fitting formula for group comparison

Resolve shared strata and the fitting formula for group comparison

## Usage

``` r
maihda_prepare_group_strata(formula, data, shared_strata, autobin = TRUE)
```

## Arguments

- formula:

  User formula.

- data:

  Full data frame.

- shared_strata:

  Logical; build strata once on the full data.

- autobin:

  Logical passed to make_strata.

## Value

list(data, formula, strata_vars).
