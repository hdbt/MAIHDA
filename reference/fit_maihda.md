# Fit MAIHDA Model

Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
Heterogeneity and Discriminatory Accuracy) using either lme4 or brms.

## Usage

``` r
fit_maihda(
  formula,
  data,
  engine = "lme4",
  family = "gaussian",
  autobin = TRUE,
  ...
)
```

## Arguments

- formula:

  A formula specifying the model. Can include a random effect for
  stratum (e.g., `outcome ~ fixed_vars + (1 | stratum)`) or can directly
  specify the intersection variables to be used for forming strata
  (e.g., `outcome ~ fixed_vars + (1 | var1:var2:var3)`). If variables
  other than "stratum" are provided in the random effect,
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  will be called internally to compute the strata and the formula will
  be updated.

- data:

  A data frame containing the variables in the formula.

- engine:

  Character string specifying which engine to use: "lme4" (default) or
  "brms".

- family:

  Character string or family object specifying the model family. Common
  options: "gaussian", "binomial", "poisson". Default is "gaussian". If
  the outcome variable appears to be binary (0/1) and the default family
  is used, the function will automatically switch to "binomial" and
  issue a warning.

- autobin:

  Logical indicating whether numeric variables used only for automatic
  strata creation should be binned by
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).
  Default is TRUE.

- ...:

  Additional arguments passed to `lmer`/`glmer` (lme4) or `brm` (brms).

## Value

A maihda_model object containing:

- model:

  The fitted model object (lme4 or brms)

- engine:

  The engine used ("lme4" or "brms")

- formula:

  The model formula

- data:

  The data used for fitting

- family:

  The family used

- strata_info:

  The strata information from make_strata() if available, NULL otherwise

## Examples

``` r
# \donttest{
# Standard approach: manually create strata first
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race", "education"))
model <- fit_maihda(health_outcome ~ age + (1 | stratum),
                    data = strata_result$data,
                    engine = "lme4")

# Simplified approach: specify stratifying variables directly in the grouping structure
# The function internally calls make_strata() to create intersectionals
model2 <- fit_maihda(health_outcome ~ age + (1 | gender:race:education),
                     data = maihda_sim_data,
                     engine = "lme4")
# }
```
