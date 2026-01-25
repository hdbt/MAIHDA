# Calculate Proportional Change in Between-Stratum Variance (PVC)

Calculates the proportional change in between-stratum variance (PVC)
between two MAIHDA models. The PVC measures how much the between-stratum
variance changes when moving from one model to another, and is
calculated as: PVC = (Var_model1 - Var_model2) / Var_model1

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
  for PVC. Default is FALSE.

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

When bootstrap = TRUE, the function resamples the data with replacement
and refits both models for each bootstrap sample to obtain confidence
intervals for the PVC estimate.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit two models
model1 <- fit_maihda(outcome ~ age + (1 | stratum), data = data)
model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum), data = data)

# Calculate PVC without bootstrap
pvc_result <- calculate_pvc(model1, model2)
print(pvc_result$pvc)

# Calculate PVC with bootstrap CI
pvc_boot <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 500)
print(pvc_boot)
} # }
```
