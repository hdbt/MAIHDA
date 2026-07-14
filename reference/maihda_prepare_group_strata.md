# Resolve shared strata and the fitting formula for group comparison

Resolve shared strata and the fitting formula for group comparison

## Usage

``` r
maihda_prepare_group_strata(
  formula,
  data,
  shared_strata,
  autobin = TRUE,
  bin_rows = NULL
)
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

- bin_rows:

  Optional logical mask (over rows of `data`) restricting the rows used
  to compute shared numeric auto-bin cut-points to the pooled analytic
  sample; forwarded to
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).
  `NULL` bins on all rows.

## Value

list(data, formula, strata_vars).
