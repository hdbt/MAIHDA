# Stepwise Proportional Change in Variance (PCV)

Estimates the proportional change in variance (PCV) sequentially by
fitting intermediate (partially-adjusted) models, adding each predictor
one-by-one. The step-specific PCV is the change in between-stratum
variance contributed by a predictor *given the variables already in the
model*. Because the steps are sequential it is order-dependent: it
reflects each variable's marginal, model-dependent change, not an
order-invariant “unique” contribution.

## Usage

``` r
stepwise_pcv(data, outcome, vars, engine = "lme4", family = "gaussian")
```

## Arguments

- data:

  Data frame with observations. Ensure \`make_strata()\` was run first
  so the \`stratum\` variable exists.

- outcome:

  Character string; the dependent variable.

- vars:

  Character vector; predictors (strata groupings & covariates) to add
  sequentially to the model.

- engine:

  Modeling engine ("lme4" or "brms"). Default is "lme4".

- family:

  Error distribution and link function. Default is "gaussian".

## Value

A data.frame showing the sequential models, the between-stratum variance
at each step, and both the step-specific and total PCV.

## Details

All models are fit on the complete cases for \`outcome\`, \`stratum\`,
and all variables in \`vars\` so that each sequential variance
comparison uses the same analytic sample.

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
stepwise_pcv(strata_result$data, "health_outcome", c("gender", "race", "age"))
#>   Step      Model        Added_Variable  Variance   Step_PCV  Total_PCV
#> 1    0 Null Model None (Intercept only) 26.714703  0.0000000  0.0000000
#> 2    1    Model 1                gender 30.862914 -0.1552782 -0.1552782
#> 3    2    Model 2                  race  2.346242  0.9239786  0.9121741
#> 4    3    Model 3                   age  3.031709 -0.2921554  0.8865153
# }
```
