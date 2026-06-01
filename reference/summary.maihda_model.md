# Summarize MAIHDA Model

Provides a summary of a MAIHDA model including variance partition
coefficients (VPC/ICC) and stratum-specific estimates.

## Usage

``` r
# S3 method for class 'maihda_model'
summary(object, bootstrap = FALSE, n_boot = 1000, conf_level = 0.95, ...)
```

## Arguments

- object:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- bootstrap:

  Logical indicating whether to compute bootstrap confidence intervals
  for VPC/ICC. Default is FALSE. Currently supported for lme4 models
  only.

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

  Data frame of stratum-specific random effects with labels if available

- fixed_effects:

  Fixed effects estimates

- model_summary:

  Original model summary

## Note

For `brms` models the VPC/ICC is currently a point estimate computed
from the posterior *summary* of the random-effect standard deviations –
i.e. it squares the posterior mean SD (\\E\[\sigma\]^2\\) rather than
averaging the VPC over posterior draws (\\E\[\sigma^2\]\\), and it does
not return a credible interval. It is therefore slightly biased and
omits posterior uncertainty; a future release will compute the VPC from
posterior draws. For lme4 models, use `bootstrap = TRUE` for an
interval.

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
summary_result <- summary(model)

# With bootstrap CI
# summary_boot <- summary(model, bootstrap = TRUE, n_boot = 50)
# }
```
