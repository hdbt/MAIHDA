# Compare MAIHDA Models

Compares variance partition coefficients (VPC/ICC) across multiple
MAIHDA models, with optional bootstrap confidence intervals.

## Usage

``` r
compare_maihda(..., model_names = NULL, bootstrap = FALSE, 
               n_boot = 1000, conf_level = 0.95)
```

## Arguments

- ...:

  Multiple maihda_model objects to compare.

- model_names:

  Optional character vector of names for the models.

- bootstrap:

  Logical indicating whether to compute bootstrap confidence intervals.
  Default is FALSE.

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default is 0.95.

## Value

A data frame comparing VPC/ICC across models with optional confidence
intervals.

## Examples

``` r
if (FALSE) { # \dontrun{
model1 <- fit_maihda(outcome ~ age + (1 | stratum), data = data1)
model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum), data = data2)

# Compare without bootstrap
comparison <- compare_maihda(model1, model2, 
                            model_names = c("Base", "With Gender"))

# Compare with bootstrap CI
comparison_boot <- compare_maihda(model1, model2,
                                 model_names = c("Base", "With Gender"),
                                 bootstrap = TRUE, n_boot = 500)
} # }
```
