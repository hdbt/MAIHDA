# Fit MAIHDA Model

Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
Heterogeneity and Discriminatory Accuracy) using either lme4 or brms.

## Usage

``` r
fit_maihda(formula, data, engine = "lme4", family = "gaussian", ...)
```

## Arguments

- formula:

  A formula specifying the model. Should include random effect for
  stratum (e.g., `outcome ~ fixed_vars + (1 | stratum)`).

- data:

  A data frame containing the variables in the formula.

- engine:

  Character string specifying which engine to use: "lme4" (default) or
  "brms".

- family:

  Character string or family object specifying the model family. Common
  options: "gaussian", "binomial", "poisson". Default is "gaussian".

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Create strata
strata_result <- make_strata(data, vars = c("gender", "race"))

# Fit model with lme4
model <- fit_maihda(outcome ~ age + (1 | stratum),
                    data = strata_result$data,
                    engine = "lme4")

# Fit model with brms (if installed)
model_brms <- fit_maihda(outcome ~ age + (1 | stratum),
                         data = strata_result$data,
                         engine = "brms")
} # }
```
