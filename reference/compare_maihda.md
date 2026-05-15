# Compare MAIHDA Models

Compares variance partition coefficients (VPC/ICC) across multiple
MAIHDA models, with optional bootstrap confidence intervals.

## Usage

``` r
compare_maihda(
  ...,
  model_names = NULL,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95
)
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
# \donttest{
# Create strata and models using simulated data
strata_1 <- make_strata(maihda_sim_data, vars = c("gender", "race"))
strata_2 <- make_strata(maihda_sim_data, vars = c("gender", "race", "education"))

model1 <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_1$data)
model2 <- fit_maihda(health_outcome ~ age + gender + (1 | stratum), data = strata_2$data)

# Compare without bootstrap
comparison <- compare_maihda(model1, model2,
                            model_names = c("Base", "With Gender"))

# Compare with bootstrap CI
comparison_boot <- compare_maihda(model1, model2,
                                 model_names = c("Base", "With Gender"),
                                 bootstrap = TRUE, n_boot = 500)
#> boundary (singular) fit: see help('isSingular')
# }
```
