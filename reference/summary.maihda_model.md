# Summarize MAIHDA Model

Provides a summary of a MAIHDA model including variance partition
coefficients (VPC/ICC) and stratum-specific estimates.

## Usage

``` r
# S3 method for class 'maihda_model'
summary(
  object,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95,
  response_vpc = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A maihda_model object from
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- bootstrap:

  Logical indicating whether to compute parametric bootstrap confidence
  intervals for VPC/ICC. Default is FALSE. Supported for lme4 models
  only; `brms` models always return a posterior credible interval (see
  Details), so `bootstrap = TRUE` is rejected for them. For a
  negative-binomial model (`glmer.nb`) the bootstrap refits via
  [`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html), which
  holds the dispersion parameter theta fixed at its original estimate,
  so the interval is conditional on the estimated theta (theta's own
  sampling uncertainty is not propagated). The `ordinal` (clmm) engine
  has no simulate/refit machinery, so `bootstrap = TRUE` is rejected
  there (use `engine = "brms"` for interval estimates).

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000.

- conf_level:

  Confidence level for the VPC/ICC interval – the lme4 bootstrap CI or
  the brms posterior credible interval. Default is 0.95.

- response_vpc:

  Logical; for a binomial (lme4) model, also compute the response-scale
  VPC
  ([`maihda_vpc_response`](https://hdbt.github.io/MAIHDA/reference/maihda_vpc_response.md))
  and attach it as the `vpc_response` slot. It is estimated by
  simulation, so it is opt-in (default `FALSE`) and uses `seed` for
  reproducibility. Ignored for other families/engines.

- seed:

  Optional integer seed for the response-scale VPC simulation when
  `response_vpc = TRUE`.

- ...:

  Additional arguments (not currently used).

## Value

A maihda_summary object containing:

- vpc:

  Variance Partition Coefficient (ICC); for lme4 with `bootstrap = TRUE`
  and for all brms models this includes
  `ci_lower`/`ci_upper`/`conf_level`. For a contextual cross-classified
  fit this is the *between-stratum* share of all unexplained variance
  (net of the context)

- variance_components:

  Data frame of variance components. For a contextual cross-classified
  fit (`fit_maihda(context = )`) each context appears as its own
  `Context: <name>` row

- context:

  For a contextual cross-classified fit, the stratum vs. context
  partition: per-context variances and shares, the contexts' total share
  (`vpc_context_total`, with an interval when bootstrapped or for brms),
  and the between-stratum share (`vpc_stratum`); `NULL` otherwise

- discriminatory_accuracy:

  For a binomial/Bernoulli outcome, the `maihda_da` object (AUC + MOR)
  from
  [`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md);
  `NULL` otherwise. A contextual fit (`fit_maihda(context = )`) is
  included – its headline AUC is the intersectional-scope concordance
  that excludes the context random effect. `NULL` for a
  crossed-dimensions fit (whose headline here is the
  additive/interaction decomposition) and a longitudinal fit

- vpc_response:

  The response-scale VPC (`maihda_vpc_response`) when
  `response_vpc = TRUE` for a binomial lme4 model, including a
  contextual fit (the context variance enters the VPC denominator);
  `NULL` otherwise (including for crossed-dimensions and longitudinal
  fits)

- stratum_estimates:

  Data frame of stratum-specific random effects with labels if available

- fixed_effects:

  Fixed effects estimates

- thresholds:

  For a cumulative (ordinal) clmm fit, the threshold (cut point)
  estimates with standard errors – the cumulative model's "intercepts";
  NULL otherwise

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
variance uses the usual latent-scale approximation. For the count
families (Poisson, negative binomial) the marginal expected counts are
propagated *per draw* – from the fixed-part linear-predictor draws and
each draw's total random-intercept variance – for the intercept-only VPC
structures (strata, crossed-dimensions, contextual), so the credible
interval reflects fixed-effect and random-effect variance uncertainty;
the negative-binomial `shape` draws are always propagated. A
random-slope (longitudinal) structure, whose per-row random-effect
design that fast path does not carry, instead holds the marginal
expected counts at a posterior-mean plug-in (constant across draws) to
avoid an expensive \\ndraws \times nobs\\ computation.

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
latent/distributional approximation (\\\pi^2/3\\ for logistic;
\\\log(1 + 1/\lambda)\\ for Poisson per Stryhn et al. 2006 and
\\\log(1 + 1/\lambda + 1/\theta)\\ for the negative binomial per
Nakagawa, Johnson & Schielzeth 2017, each evaluated at a single
*marginal* expected count \\\lambda\\: the mean over the analytic sample
of the row-level \\\lambda_i = \exp(x_i'\beta + v_i/2)\\ – the
fixed-part prediction with the log-normal correction for the row's total
random-effect variance \\v_i\\. The counts are averaged *before* the
transform, which is where the cited \\\lambda\\ is defined and which
reduces to Nakagawa et al.'s \\\lambda = \exp(\beta_0 + \sigma^2/2)\\ in
the null model; *not* at the conditional fitted means, whose BLUPs would
tie the level-1 variance to the realized random effects), so the VPC is
on that latent scale; for a *weighted* Gaussian model the level-1
variance is the mean conditional residual variance, \\\bar{\sigma^2 /
w_i}\\, since the per-observation residual variance is \\\sigma^2 /
w_i\\. The stratum random effects represent the total between-stratum
deviation; they equal the *pure* intersectional (interaction) component
only when the additive main effects of the strata variables are included
in the model.

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
