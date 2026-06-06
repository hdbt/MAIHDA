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

  Character string, family object, or family function specifying the
  model family. Common options: "gaussian", "binomial", "poisson".
  Default is "gaussian". If the outcome variable appears to be binary
  and the default family is used, the function will automatically switch
  to "binomial", recode two-level responses to 0/1 for `glmer()`, and
  issue a warning. When a two-level non-0/1 response is recoded (on
  either the auto-detected or an explicit `family = "binomial"` path),
  the mapping follows the usual convention – the first level becomes 0
  (reference) and the second becomes 1 (the modeled event), where
  "first/second" means alphabetical order for a character outcome and
  the declared order for a factor. The chosen mapping is reported via a
  [`message()`](https://rdrr.io/r/base/message.html) and stored on the
  result as `$response_recoding`; set the factor levels (or supply a 0/1
  outcome) to control which level is the event. Although any valid
  family object is accepted for fitting, the MAIHDA variance summaries
  ([`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md),
  VPC/ICC, PCV) are only defined for `gaussian("identity")`, the
  binomial/Bernoulli families with a logit or probit link, and
  `poisson("log")`. Other families (for example `Gamma(link = "log")`)
  will fit, but [`summary()`](https://rdrr.io/r/base/summary.html) and
  the VPC/PCV helpers will stop with an "not implemented" error because
  no level-1 variance is defined for them.

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

- response_recoding:

  For a recoded two-level outcome, a data frame mapping each original
  level to its 0/1 value and role (reference/event); NULL when no
  recoding occurred

- diagnostics:

  Fit-quality diagnostics (singular fit / convergence) for lme4 models,
  surfaced by the print and summary methods

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
