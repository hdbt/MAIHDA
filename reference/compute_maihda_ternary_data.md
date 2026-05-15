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

  Optional output from \`summary()\`.

- scale:

  Character, either "link" or "response".

- reference_values:

  List or data.frame of reference values for covariates.

- uncertainty_method:

  Character indicating how to extract uncertainty. "auto" uses
  conditional standard errors for lme4 models and posterior standard
  deviations for brms models. "ci_width" uses the 95% interval width.

- include_na_strata:

  Logical, whether to include strata with missing data.

- verbose:

  Logical, whether to print messages.

## Value

A tidy tibble with ternary coordinates.
