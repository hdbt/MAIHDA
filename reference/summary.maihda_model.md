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
  df_method = c("between-within", "normal", "bootstrap"),
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
  there (use `engine = "brms"` for interval estimates). For a Gaussian
  model carrying lme4 precision `weights`, the simulated responses draw
  each residual at \\\sigma / \sqrt{w_i}\\, so the interval rests on the
  same \\\sigma^2 / w_i\\ semantics as the point estimate above.

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

- df_method:

  Reference distribution for the fixed-effect p-values and intervals of
  an `lme4` fit: `"between-within"` (default) a \\t\\ on containment
  degrees of freedom for a Gaussian fit and a z elsewhere, `"normal"` a
  z, `"bootstrap"` a null-restricted parametric bootstrap costing
  `n_boot` refits *per fixed-effect term*. `"bootstrap"` is the
  reference to use for a GLMM term that is constant within a stratum,
  such as an adjusted model's dimension main effects. Every other engine
  uses a z regardless.

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

- longitudinal:

  For a longitudinal (growth-curve) fit, the time-varying summary:
  `vpc_t` (the VPC over a reporting grid, with bootstrap or posterior
  bands), the per-level variances over that grid, the stratum and
  individual covariance blocks, and the two *trajectory VPCs*
  `vpc_intercept` and `vpc_slope` described below, each with an interval
  in `vpc_intercept_ci` / `vpc_slope_ci` (a posterior credible interval
  for brms, a bootstrap interval for a bootstrapped lme4 fit, `NA`
  otherwise) and the basis in `trajectory_vpc_method`. `NULL` for a
  cross-sectional fit

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

- count_vpc:

  For a log-link count model, the `approximation` the level-1 variance
  came from, the marginal count `lambda` (and `theta` /
  `lambda_effective` for the negative binomial) it was evaluated at, the
  `alternatives` all three approximations give at that `lambda`
  (`level1_variance` is the one used), and `low_count` – `TRUE` below
  the \\\lambda = 2\\ threshold above which Nakagawa et al. (2017)
  report the three agree. These are plug-in values at a single `lambda`:
  on the likelihood engines that is exactly the number in
  `variance_components`, but a `brms` summary works draw by draw and
  reports \\E\[\sigma^2_e\]\\ there, which differs slightly. `NULL` for
  every other family.

- vpc_response:

  The response-scale VPC (`maihda_vpc_response`) when
  `response_vpc = TRUE` for a binomial lme4 model, including a
  contextual fit (the context variance enters the VPC denominator);
  `NULL` otherwise (including for crossed-dimensions and longitudinal
  fits)

- stratum_estimates:

  Data frame of stratum-specific random effects with labels if available

- fixed_effects:

  Fixed-effect estimates. For the lme4, WeMix and ordinal engines a data
  frame with `term`, `estimate`, `se`, `statistic`, `df`, `p_value` and
  the Wald interval `lower`/`upper` at `conf_level`; `df` is `NA`
  wherever the reference is a z. The WeMix standard errors are its
  sandwich (robust) ones. For brms, the
  [`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
  matrix (posterior mean, `Est.Error` and the credible-interval
  quantiles at `conf_level`). Available in a tidy, engine-independent
  shape from `tidy(x, component = "fixed")`

- conf_level:

  The interval level used for the fixed effects (and, when bootstrapped
  or Bayesian, the VPC)

- df_method:

  The reference distribution the `fixed_effects` table used,
  `"between-within"`, `"normal"` or `"bootstrap"`

- thresholds:

  For a cumulative (ordinal) clmm fit, the threshold (cut point)
  estimates with standard errors – the cumulative model's "intercepts";
  NULL otherwise

- model_summary:

  Original model summary

- diagnostics:

  Fit-quality diagnostics (singular fit / convergence) carried over from
  the fitted model and reported by the print method

- strata_autobin_info:

  The auto-binning recipe carried over from the fitted model: for each
  numeric stratum dimension
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  discretised, its `breaks` and `labels`. The cut-points are quantiles
  of the analytic sample, so they define the strata (and hence the
  estimand); the print method reports them. An empty list when nothing
  was binned

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
Nakagawa, Johnson & Schielzeth 2017 – their `"delta"` and `"trigamma"`
alternatives are available via `fit_maihda(count_approximation = )`, and
the choice is reported in the printed summary and in `count_vpc` because
the three diverge materially below a marginal count of 2 – each
evaluated at a single *marginal* expected count \\\lambda\\: the mean
over the analytic sample of the row-level \\\lambda_i = \exp(x_i'\beta +
v_i/2)\\ – the fixed-part prediction with the log-normal correction for
the row's total random-effect variance \\v_i\\. The counts are averaged
*before* the transform, which is where the cited \\\lambda\\ is defined
and which reduces to Nakagawa et al.'s \\\lambda = \exp(\beta_0 +
\sigma^2/2)\\ in the null model; *not* at the conditional fitted means,
whose BLUPs would tie the level-1 variance to the realized random
effects), so the VPC is on that latent scale; for a *weighted* Gaussian
model the level-1 variance is the mean conditional residual variance,
\\\bar{\sigma^2 / w_i}\\, since the per-observation residual variance is
\\\sigma^2 / w_i\\. The stratum random effects represent the total
between-stratum deviation; they equal the *pure* intersectional
(interaction) component only when the additive main effects of the
strata variables are included in the model.

## Fixed-effect reference distribution

A Gaussian `lme4` fit refers each Wald statistic to a \\t\\ on
*containment* (between-within) degrees of freedom, reported in the `df`
column: a term absorbed by a random-effect grouping is tested against
that grouping's units minus the terms it absorbs, and a term absorbed by
none against \\n\\ minus the random-effect levels. A random slope counts
as absorbing. Set `df_method = "normal"` for a z instead.

A GLMM, a WeMix pseudo-ML fit and an
[`ordinal::clmm`](https://rdrr.io/pkg/ordinal/man/clmm.html) fit have no
finite-sample \\t\\ and use the Wald z; a brms summary reports the
posterior. For Kenward-Roger or Satterthwaite, apply pbkrtest or
lmerTest to `x$model`.

`df_method = "bootstrap"` replaces that reference for an `lme4` fit with
at least one fixed-effect term, Gaussian or not, and is the one to use
for a GLMM – whose z is anticonservative for a term constant within a
stratum, most severely when the strata are few. For each fixed-effect
term the model is refitted with that term's coefficients *constrained to
zero*, `n_boot` responses are simulated from the restricted fit, the
full model is refitted on each, and the observed Wald statistic is
referred to the resulting distribution of \\\|t^\*\|\\ under a true
null. The estimate and standard error are unchanged; the p-value and the
interval both come from that distribution and agree exactly, zero
falling outside the interval precisely when the p-value is significant.
`df` is `NA`, and so are the intercept's p-value and interval: a MAIHDA
intercept is a reference-category level rather than a term that can be
dropped, so it has no null model to simulate from.

The constraint is imposed on the fitted design and verified, not assumed
from the formula. Removing a term from a formula does not always remove
it from the model: R's marginality rules recode a surviving higher-order
term to absorb a dropped marginal one, so for `y ~ x * f` the formula
`. ~ . - x` still spans the original column space and leaves the
coefficient under test entirely unrestricted. The same holds for either
main effect of `f * g`, for every main effect and two-way term under a
three-way interaction, and for a nested `f / g`. Where that happens the
term's design columns are constrained directly instead. A model whose
fixed part is additive is unaffected: there, dropping the term from the
formula already is the null.

It costs `n_boot` refits *per term*, and is a separate bootstrap from
the `bootstrap = TRUE` VPC interval, which is not reused. The smallest
reportable p-value is \\1 / (n\\boot + 1)\\.

Budget for it. A Gaussian refit takes milliseconds, but a binomial one
takes about a second at \\n = 1000\\ and tens of seconds at \\n =
6000\\, so the default `n_boot = 1000` on a three-dimension GLMM is
roughly an hour at the smaller size and impractical at the larger. The
p-value is exact at any `n_boot` for which \\(n\\boot + 1)\alpha\\ is a
whole number – 199 and 999 at the 5% level – while the interval
endpoints, being order statistics, keep tightening with more draws;
`n_boot = 199` is the usual compromise for a GLMM.

## Two VPCs for a longitudinal fit

A longitudinal summary reports **two different variance partitions**,
and they are not interchangeable. They differ in one term of the
denominator: the level-1 (occasion) residual variance \\\sigma^2_e\\,
which in a growth model is *within-individual volatility* – how far a
single measurement falls from that person's own smooth trajectory. It
mixes measurement error, genuine short-term fluctuation, and any misfit
of the assumed functional form. A cross-sectional MAIHDA cannot separate
it from between-individual variance at all; repeated measures are what
split the two.

**The headline VPC** (`vpc`, and `longitudinal$vpc_t` over time) keeps
it: \$\$VPC_S(t) = \frac{Var_S(t)}{Var_S(t) + Var_I(t) +
\sigma^2_e}.\$\$ This is the discriminatory-accuracy question – how much
of an *observed measurement* at time \\t\\ a stratum accounts for – and
it is the quantity comparable to published cross-sectional MAIHDA VPCs.
Report this one unless you specifically mean the trajectory question.

**The trajectory VPCs** (`longitudinal$vpc_intercept` and
`longitudinal$vpc_slope`) drop it, following Bell et al. (2024),
equation (5): \$\$VPC\_{intercept} = \frac{Var_S(t_0)}{Var_S(t_0) +
Var_I(t_0)}, \qquad VPC\_{slope} =
\frac{SlopeVar_S(t_0)}{SlopeVar_S(t_0) + SlopeVar_I(t_0)}.\$\$ These ask
how intersectionally patterned people's *trajectories* are i.e., what
share of the between-individual variation in where a trajectory starts,
and in how fast it changes, lies between strata. Because \\\sigma^2_e\\
is absent neither is affected by how noisy the outcome measure is, which
makes them comparable across studies using different instruments.

Both are evaluated at the baseline \\t_0\\ (`ref_time`, the earliest
observed time), pairing with `PCV_intercept` and `PCV_slope` from
`maihda(decomposition = "longitudinal")`. The intercept VPC depends on
where time is zeroed and the slope VPC does not, as Bell et al. note;
their own examples centre on mean age rather than the baseline, so an
intercept VPC replicated from the paper will differ from the one
reported here unless the reference points are aligned. `vpc_slope` is
`NA` when the model was fit with `stratum_slope = FALSE` (no
between-stratum slope variance exists to take a share of).

Both come with an interval in `vpc_intercept_ci` / `vpc_slope_ci`, and
`trajectory_vpc_method` records its basis. For a `brms` fit the two
shares are computed *per posterior draw* and reported as the posterior
median with a credible interval, matching `vpc_t` and the headline VPC
on the same fit; for an `lme4` fit the point estimates are the plug-in
from the fitted covariance blocks and `summary(bootstrap = TRUE)` adds a
parametric-bootstrap interval (`NA` without one). Report the interval:
these shares are poorly determined when the strata are few, and one
spanning half the unit interval is an ordinary result rather than an
unusual one – the twelve strata of `maihda_long_data`, fitted with
`brms` over 150 individuals, give an intercept VPC of 0.58 running from
0.34 to 0.82.

## References

Bell, A., Evans, C., Holman, D., & Leckie, G. (2024). Extending
intersectional multilevel analysis of individual heterogeneity and
discriminatory accuracy (MAIHDA) to study individual longitudinal
trajectories, with application to mental health in the UK. *Social
Science & Medicine*, 351, 116955.
[doi:10.1016/j.socscimed.2024.116955](https://doi.org/10.1016/j.socscimed.2024.116955)

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
