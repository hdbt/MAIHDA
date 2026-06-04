# Calculate Proportional Change in Between-Stratum Variance (PCV)

Calculates the proportional change in between-stratum variance (PCV)
between two MAIHDA models. The PCV measures how much the between-stratum
variance changes when moving from one model to another, and is
calculated as: PCV = (Var_model1 - Var_model2) / Var_model1. (The
function and result object retain the historical "pvc" naming; “PVC” and
“PCV” refer to the same quantity.)

## Usage

``` r
calculate_pvc(
  model1,
  model2,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95
)
```

## Arguments

- model1:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  This is the reference model (typically a simpler or baseline model).

- model2:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  This is the comparison model (typically a more complex model with
  additional predictors).

- bootstrap:

  Logical indicating whether to compute bootstrap confidence intervals
  for the PCV. Default is FALSE.

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default is 0.95.

## Value

A list containing:

- pvc:

  The estimated proportional change in variance

- var_model1:

  Between-stratum variance from model1

- var_model2:

  Between-stratum variance from model2

- ci_lower:

  Lower bound of confidence interval (if bootstrap = TRUE)

- ci_upper:

  Upper bound of confidence interval (if bootstrap = TRUE)

- bootstrap:

  Logical indicating if bootstrap was used

## Details

The PVC is interpreted as the proportional reduction (or increase if
negative) in between-stratum variance when moving from model1 to model2.
A positive PVC indicates that model2 explains some of the
between-stratum variance present in model1, while a negative PVC
suggests that model2 has more unexplained between-stratum variance.

When bootstrap = TRUE, the function uses a parametric bootstrap: it
simulates new responses from model2 and refits both models with
[`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html) for each
simulated response to obtain confidence intervals for the PVC estimate.

## Examples

``` r
# \donttest{
# Create strata and fit two models
strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
model1 <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
model2 <- fit_maihda(health_outcome ~ age + gender + (1 | stratum), data = strata_result$data)

# Calculate PVC without bootstrap
pvc_result <- calculate_pvc(model1, model2)
print(pvc_result$pvc)
#> [1] -0.1626748

# Calculate PVC with bootstrap CI
# pvc_boot <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 500)
# print(pvc_boot)
# }
```
