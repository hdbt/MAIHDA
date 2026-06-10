# Reconstruct the adjusted-model main-effect terms for a MAIHDA decomposition

For each stratum-defining variable, returns the model term to add as an
additive fixed main effect in the adjusted model, plus the data
augmented with any reconstructed binned factors. A categorical dimension
is used directly; a numeric dimension that
[`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
tertile-binned is reconstructed as the *same* binned factor (using the
stored breaks/labels), because
[`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
bins only a temporary copy and leaves the original numeric column intact
– adding the raw numeric column would wrongly enter a linear term
instead of the factor that defines the strata.

## Usage

``` r
maihda_adjusted_terms(strata_vars, autobin_info, data)
```

## Arguments

- strata_vars:

  Character vector of stratum-defining variable names.

- autobin_info:

  Named list of `list(breaks, labels)` per auto-binned variable (the
  `strata_autobin_info` stored on a `maihda_model`).

- data:

  Data frame containing the original stratum-defining columns.

## Value

A list with `terms` (character vector of RHS term names) and `data` (the
input augmented with any `.maihda_dim_*` binned columns).
