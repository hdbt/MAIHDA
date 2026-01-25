# Summarize MAIHDA Model

Provides a summary of a MAIHDA model including variance partition
coefficients (VPC/ICC) and stratum-specific estimates.

## Usage

``` r
summary_maihda(object, bootstrap = FALSE, n_boot = 1000, 
               conf_level = 0.95, ...)
```

## Arguments

- object:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- bootstrap:

  Logical indicating whether to compute bootstrap confidence intervals
  for VPC/ICC. Default is FALSE.

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default is 0.95.

- ...:

  Additional arguments (not currently used).

## Value

A maihda_summary object containing:

- vpc:

  Variance Partition Coefficient (ICC) with optional CI

- variance_components:

  Data frame of variance components

- stratum_estimates:

  Data frame of stratum-specific random effects

- fixed_effects:

  Fixed effects estimates

- model_summary:

  Original model summary

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)
summary_result <- summary_maihda(model)

# With bootstrap CI
summary_boot <- summary_maihda(model, bootstrap = TRUE, n_boot = 500)
} # }
```
