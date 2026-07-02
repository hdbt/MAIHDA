# Calculate Proportional Change in Between-Stratum Variance (PCV)

Calculates the proportional change in between-stratum variance (PCV)
between two MAIHDA models. The PCV measures how much the between-stratum
variance changes when moving from one model to another, and is
calculated as: PCV = (Var_model1 - Var_model2) / Var_model1.

## Usage

``` r
calculate_pcv(
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

- pcv:

  The estimated proportional change in variance

- pvc:

  Deprecated duplicate of `pcv`, kept so code written against the
  historical
  [`calculate_pvc()`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)
  spelling keeps working; it will be removed in a future release

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

The PCV is the proportional change in between-stratum variance when
moving from model1 to model2: a positive value means model2 has lower
between-stratum variance, a negative value means higher. It is the share
of model1's between-stratum variance *explained* by model2 only in the
canonical nested case, where model2 adds fixed-effect predictors to
model1 on the same outcome, analytic sample and strata. The function
does not require nesting, so for non-nested models the PCV is simply a
model-dependent difference in variance, not an explained proportion.

**REML vs ML.** `lmer` fits Gaussian models by REML, whose
between-stratum variance estimate is *not* comparable across models with
different fixed effects – exactly the canonical null-vs-adjusted PCV,
where the adjusted model adds the dimensions' main effects.
`calculate_pcv()` therefore refits any REML `lmer` model with maximum
likelihood ([`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html))
before reading the variances (and before the parametric bootstrap, so
the interval matches), matching
[`maihda_ic`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md) and
[`anova()`](https://rdrr.io/r/stats/anova.html) on `lme4` models. Using
REML estimates here biases the PCV (it overstates the residual
between-stratum variance of the adjusted model). GLMM fits (`glmer`) and
the brms/wemix/ordinal engines are already on the maximum-likelihood
scale and are unaffected; single-model VPC/ICC summaries keep their REML
fit, since that comparison-free quantity is not subject to the pitfall.

When bootstrap = TRUE, the function uses a parametric bootstrap: it
simulates new responses from model2 and refits both models with
[`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html) for each
simulated response to obtain confidence intervals for the PCV estimate.
For negative-binomial models (`glmer.nb`) `refit()` holds the dispersion
parameter theta fixed at its original estimate, so the interval is
conditional on the estimated theta.

## See also

[`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
for the sequential (one-variable-at-a-time) PCV, and
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) which
computes the canonical null-vs-adjusted PCV automatically.

## Examples

``` r
# \donttest{
# Create strata and fit two models
strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
model1 <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
model2 <- fit_maihda(health_outcome ~ age + gender + (1 | stratum), data = strata_result$data)

# Calculate PCV without bootstrap
pcv_result <- calculate_pcv(model1, model2)
print(pcv_result$pcv)
#> [1] 0.007185334

# Calculate PCV with bootstrap CI
# pcv_boot <- calculate_pcv(model1, model2, bootstrap = TRUE, n_boot = 500)
# print(pcv_boot)
# }
```
