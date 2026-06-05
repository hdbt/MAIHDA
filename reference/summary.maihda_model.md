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

  Logical indicating whether to compute parametric bootstrap confidence
  intervals for VPC/ICC. Default is FALSE. Supported for lme4 models
  only; `brms` models always return a posterior credible interval (see
  Details), so `bootstrap = TRUE` is rejected for them.

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for the VPC/ICC interval – the lme4 bootstrap CI or
  the brms posterior credible interval. Default is 0.95.

- ...:

  Additional arguments (not currently used).

## Value

A maihda_summary object containing:

- vpc:

  Variance Partition Coefficient (ICC); for lme4 with `bootstrap = TRUE`
  and for all brms models this includes
  `ci_lower`/`ci_upper`/`conf_level`

- variance_components:

  Data frame of variance components

- stratum_estimates:

  Data frame of stratum-specific random effects with labels if available

- fixed_effects:

  Fixed effects estimates

- model_summary:

  Original model summary

- diagnostics:

  Fit-quality diagnostics (singular fit / convergence) carried over from
  the fitted model and reported by the print method

## Note

For `lme4` models a VPC/ICC interval is obtained from a parametric
bootstrap (`bootstrap = TRUE`). For `brms` models the VPC/ICC is
summarised directly from the posterior draws: the reported estimate is
the posterior median of the per-draw VPC (\\E\[\sigma^2\]\\-based, not
the biased \\E\[\sigma\]^2\\) and the interval is a central credible
interval at `conf_level` (default 95%), so no `bootstrap` argument is
needed. The variance-components table reports the posterior-mean
variance components, so the stratum proportion shown there may differ
slightly from the headline VPC because the median of a ratio is not the
ratio of means. For non-Gaussian `brms` families the level-1 (residual)
variance uses the usual latent-scale approximation; for `poisson(log)`
it is evaluated at the posterior-mean fitted values rather than per draw
to avoid an expensive \\ndraws \times nobs\\ computation.

## Interpreting the VPC/ICC

The VPC is the between-stratum variance divided by the total
*unexplained* variance. For the canonical single-stratum model that
denominator is between-stratum + residual, but if the model includes
additional random effects (e.g. `(1 | site)`) their variance is included
in the denominator too (between-stratum + other random effects +
residual), so the VPC is the between-stratum *share* of all unexplained
variance. It is a conditional/residual ICC that excludes variance
captured by the fixed effects, so for models with covariates it is
conditional on them. It is most commonly read from the null model
`outcome ~ 1 + (1 | stratum)`, where it is the total between-stratum
share. For non-Gaussian families the level-1 (residual) variance uses a
latent/distributional approximation (e.g. \\\pi^2/3\\ for logistic), so
the VPC is on that latent scale; for a *weighted* Gaussian model the
level-1 variance is the mean conditional residual variance,
\\\bar{\sigma^2 / w_i}\\, since the per-observation residual variance is
\\\sigma^2 / w_i\\. The stratum random effects represent the total
between-stratum deviation; they equal the *pure* intersectional
(interaction) component only when the additive main effects of the
strata variables are included in the model.

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
