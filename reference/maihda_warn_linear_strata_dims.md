# Warn when a numeric stratum dimension enters the adjusted model as a linear term

In MAIHDA the stratum-defining dimensions are categorical, but the
adjusted / decomposition model adds an un-binned numeric dimension by
its raw column name – i.e. as a single linear slope, not categorical
main effects (see the `else` branch of
[`maihda_adjusted_terms`](https://hdbt.github.io/MAIHDA/reference/maihda_adjusted_terms.md)).
That silently changes the PCV/decomposition interpretation. This helper
flags the offending dimensions and warns; it does *not* alter the model
(the term is still entered as-is). The flagged case is a stratum
variable that is *numeric* in the data, *not* a key in `autobin_info`,
and has **three or more** distinct values – either a category code with
few unique values, or a many-valued numeric fitted with
`autobin = FALSE`. Categorical (factor/character) dimensions and
auto-binned numerics – reconstructed as the `.maihda_dim_*` tertile
factors – are excluded, and so is a binary (two-level) numeric
dimension: with only two levels a linear term and the two-level factor
span the same design space, so the fit and the PCV coincide and there is
nothing to warn about.

## Usage

``` r
maihda_warn_linear_strata_dims(strata_vars, autobin_info, data, fn = "maihda")
```

## Arguments

- strata_vars:

  Character vector of stratum-defining variable names.

- autobin_info:

  Auto-binning recipe (`strata_autobin_info`); its names are the
  auto-binned numeric dimensions, which are excluded from the warning.

- data:

  Data frame carrying the original stratum-defining columns.

- fn:

  Name of the calling function, used in the warning prefix.

## Value

Invisibly, the character vector of flagged dimension names (empty if
none).

## Details

Called once per decomposition at each public entry point (not inside
[`maihda_adjusted_terms()`](https://hdbt.github.io/MAIHDA/reference/maihda_adjusted_terms.md),
which runs several times per decomposition), so the warning fires
exactly once.
