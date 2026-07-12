# Build the crossed-dimensions-model formula and data for a MAIHDA decomposition

The crossed-dimensions alternative to the two-model (fixed-effects PCV)
decomposition (the function name keeps the historical "cross_classified"
spelling). Given a null model's formula (in `(1 | stratum)` form,
carrying only the covariates) and the stratum metadata, returns the
single crossed formula – the covariates plus an *additive random
intercept for each stratum dimension* plus the intersection (`stratum`)
random intercept – together with the data carrying any reconstructed
binned factors. In the fitted model each dimension's RE variance is that
dimension's additive main-effect variance and the `stratum` RE variance
is the interaction beyond additive; see
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

## Usage

``` r
maihda_cross_classified_formula(
  null_formula,
  strata_vars,
  autobin_info,
  data,
  interaction_group = "stratum",
  context = NULL
)
```

## Arguments

- null_formula:

  The null model formula using `(1 | stratum)` (covariates only – any
  dimension main effects written as fixed terms should be removed first,
  because they enter as random effects here).

- strata_vars:

  Character vector of stratum-defining variables.

- autobin_info:

  Auto-binning recipe (`strata_autobin_info`).

- data:

  The null model's data (`original_data`) with the `stratum` column and
  the dimension columns.

- interaction_group:

  Name of the intersection grouping factor (the column whose random
  intercept captures the interaction). Default `"stratum"`.

- context:

  Optional character vector of higher-level grouping variables that the
  caller re-appends as contextual random intercepts (via `context =` on
  the fit). Named here only so the extra-random-effect guard treats them
  as legitimate rather than flagging them; the builder itself does not
  add them.

## Value

A list with `formula`, `data`, `dim_groups` (a named character vector
mapping each `strata_var` to its random-effect grouping-factor name) and
`interaction_group` (`"stratum"`); or `NULL` if fewer than two
dimensions are available.

## Details

Returns `NULL` when fewer than two dimensions are available (there is no
intersection to decompose). The dimension grouping factor reuses the
dimension's own column for a categorical dimension and the reconstructed
`.maihda_dim_*` tertile factor for an auto-binned numeric dimension (via
[`maihda_adjusted_terms`](https://hdbt.github.io/MAIHDA/reference/maihda_adjusted_terms.md)),
so the additive REs are crossed on exactly the levels that define the
strata.
