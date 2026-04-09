# Compute Ternary Data for MAIHDA Models

Compute Ternary Data for MAIHDA Models

## Usage

``` r
compute_maihda_ternary_data(
  model,
  summary_obj = NULL,
  scale = c("link", "response"),
  reference_values = NULL,
  uncertainty_method = c("auto", "se", "ci_width", "posterior_sd"),
  include_na_strata = FALSE,
  verbose = TRUE
)
```

## Arguments

- model:

  A fitted MAIHDA model object from \`fit_maihda()\`.

- summary_obj:

  Optional output from \`summary_maihda()\`.

- scale:

  Character, either "link" or "response".

- reference_values:

  List or data.frame of reference values for covariates.

- uncertainty_method:

  Character indicating how to extract uncertainty.

- include_na_strata:

  Logical, whether to include strata with missing data.

- verbose:

  Logical, whether to print messages.

## Value

A tidy tibble with ternary coordinates.
