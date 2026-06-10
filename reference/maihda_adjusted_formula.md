# Build the adjusted-model formula and data for a MAIHDA decomposition

Given a fitted null model's formula (in `(1 | stratum)` form) and its
stored strata metadata, returns the adjusted formula – the null formula
plus the additive main effects of the stratum dimensions – and the data
carrying any reconstructed binned factors. Returns `NULL` when fewer
than two dimensions are available, because there is no intersection to
decompose and the single dimension's main effect would render the
stratum random intercept redundant (singular).

## Usage

``` r
maihda_adjusted_formula(null_formula, strata_vars, autobin_info, data)
```

## Arguments

- null_formula:

  The null model formula using `(1 | stratum)`.

- strata_vars:

  Character vector of stratum-defining variables.

- autobin_info:

  Auto-binning recipe (`strata_autobin_info`).

- data:

  The null model's data (`original_data`) with the dimension columns.

## Value

A list with `formula` and `data`, or `NULL` if fewer than two dimensions
are available.
