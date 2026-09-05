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
  conf_level = 0.95,
  estimation = c("fitted", "ML")
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
  for the PCV. Default is FALSE. **lme4 engine only**: the parametric
  bootstrap relies on lme4's
  [`simulate()`](https://rdrr.io/r/stats/simulate.html)/`refit()`, so
  for the brms, wemix, and ordinal engines the PCV is reported as a
  point estimate and `bootstrap = TRUE` is an error (see Details).

- n_boot:

  Number of bootstrap samples if bootstrap = TRUE. Default is 1000. A
  value below about 200 warns that the interval's tail endpoints are
  unstable (the hard minimum is 10).

- conf_level:

  Confidence level for bootstrap intervals. Default is 0.95.

- estimation:

  Variance-estimation basis for the cross-model comparison, one of
  `"fitted"` (default) or `"ML"`. `"fitted"` differences each model's
  own between-stratum variance (the REML estimate for a Gaussian `lmer`
  fit); `"ML"` refits any REML `lmer` fit with maximum likelihood first,
  for a correction-free comparison. The choice affects Gaussian `lmer`
  fits only – `glmer` and the wemix/ordinal engines are already
  maximum-likelihood, and a `brms` fit is a Bayesian posterior (not ML),
  so `"ML"` is a no-op for all of them; a `brms` comparison is reported
  on the as-fitted *posterior* basis rather than as an ML-refit. See
  Details for the finite-sample tradeoff. Whenever model2 (the adjusted
  model) sits on the singularity boundary – under *any* `estimation`
  basis – the PCV is pinned near 1; this is recorded as
  `adjusted_at_boundary = TRUE` and noted by
  [`print()`](https://rdrr.io/r/base/print.html). It is not treated as
  an error or warned about: a singular fit is indistinguishable from
  genuinely additive strata (a common, legitimate result).

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

- estimation:

  The variance-estimation basis requested (`"fitted"` or `"ML"`)

- estimation_used:

  The basis actually used: `"fitted"`, `"ML"`, `"mixed"`, or
  `"posterior"`. `"mixed"` arises under `estimation = "ML"` when a model
  keeps its REML fit – its between-stratum variance is on the boundary
  so the ML refit is skipped, or
  [`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html) failed – so the
  comparison is partly REML rather than a pure, correction-free ML one.
  `"posterior"` is reported for a `brms` comparison (a Bayesian
  posterior, on which `"ML"` is a no-op).
  [`print()`](https://rdrr.io/r/base/print.html) states the basis

- adjusted_at_boundary:

  Logical; `TRUE` when model2's between-stratum variance is on the
  singularity boundary, so the PCV is pinned near 1 (100%) and its
  interval/SE are unreliable – consistent with genuinely additive strata
  as well as a degenerate fit
  ([`print()`](https://rdrr.io/r/base/print.html) states this)

- ml_refit_failed:

  Logical; `TRUE` when `estimation = "ML"` was requested but
  [`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html) failed for a
  model, so its REML fit was used instead (the comparison is then not on
  a pure ML basis)

- ci_lower:

  Lower bound of confidence interval (if bootstrap = TRUE)

- ci_upper:

  Upper bound of confidence interval (if bootstrap = TRUE)

- bootstrap:

  Logical indicating if bootstrap was used

- n_boot_nonconverged:

  Number of contributing bootstrap draws whose refit optimiser did not
  converge (`optinfo$conv$opt != 0`); these are retained in the interval
  and counted, so `n_boot_ok` does not imply convergence

- interval_reliable:

  Logical (bootstrap only); `FALSE` when more than half the contributing
  draws did not converge, in which case the interval is still returned
  but flagged unreliable (see Details)

## Details

The PCV is the proportional change in between-stratum variance when
moving from model1 to model2: a positive value means model2 has lower
between-stratum variance, a negative value means higher. It is the share
of model1's between-stratum variance *explained* by model2 only in the
canonical nested case, where model2 adds fixed-effect predictors to
model1 on the same outcome, analytic sample and strata. The function
does not require nesting, so for non-nested models the PCV is simply a
model-dependent difference in variance, not an explained proportion.

**REML vs ML (the `estimation` argument).** `lmer` fits Gaussian models
by REML, and two considerations pull in opposite directions when the PCV
differences two such fits. On one hand, the REML *likelihood* is not
comparable across models with different fixed effects, and REML applies
a model-specific degrees-of-freedom correction that differs between the
null and the adjusted fit; refitting both with maximum likelihood
([`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html)) puts the two
between-stratum variances on a common, correction-free basis, matching
[`maihda_ic`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md) and
[`anova()`](https://rdrr.io/r/stats/anova.html) on `lme4` models. On the
other hand, ML variance-component estimates are downward-biased in
finite samples – most sharply with few strata (the usual MAIHDA regime),
and more so for the adjusted model (more fixed effects, larger REML
correction) – which *inflates* the reported PCV relative to the REML
estimates. Both are defensible point estimates of each model's
between-stratum variance, so `estimation` selects between them:

- `"fitted"` (default):

  use each model's own fitted between-stratum variance – the REML
  estimate for an `lmer` Gaussian fit. This matches the variances
  [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
  reports and conventional MAIHDA practice, and avoids ML's
  finite-sample downward bias.

- `"ML"`:

  refit any REML `lmer` fit with maximum likelihood before reading the
  variances (and before the parametric bootstrap, so the interval
  matches), for a correction-free cross-model comparison.

The choice affects Gaussian `lmer` fits only: GLMM fits (`glmer`) and
the wemix/ordinal engines are already on the maximum-likelihood scale,
so `"fitted"` and `"ML"` coincide there. A `brms` fit is a Bayesian Stan
*posterior*, not a maximum-likelihood fit, so `"ML"` performs no refit
and the comparison is reported on the as-fitted posterior basis
(`estimation_used = "posterior"`), never as an ML-refit. Single-model
VPC/ICC summaries always keep their REML fit, since that comparison-free
quantity is not subject to the cross-model pitfall.

**Latent-scale families and rescaling.** For families whose level-1
variance is a fixed latent-scale constant – binomial/Bernoulli
(\\\pi^2/3\\ logit, 1 probit) and the cumulative (ordinal) model – the
linear predictor is identified only up to scale. Adding predictors that
explain *within-stratum* (individual-level) variation cannot shrink that
fixed level-1 variance; the latent scale stretches instead, inflating
the coefficients and the between-stratum variance alike (Bauer 2009;
Mood 2010). Part of a null-vs-adjusted change in the between-stratum
variance is then rescaling rather than genuinely explained variance, so
latent-scale PCVs tend to be understated and can turn negative on this
account alone. The canonical MAIHDA adjusted model – which adds the
stratum dimensions' main effects, constant *within* each stratum – is
largely unaffected, but the caveat is first-order whenever an added
predictor varies within strata (an individual-level covariate, as in the
[`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
steps that add one). The count families' level-1 variance is not a fixed
constant, but as with any non-identity link the same non-collapsibility
logic applies in attenuated form. Gaussian identity-link PCVs are not
subject to this.

When bootstrap = TRUE, the function uses a parametric bootstrap: it
simulates new responses from model2 and refits both models with
[`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html) for each
simulated response to obtain confidence intervals for the PCV estimate.
For negative-binomial models (`glmer.nb`) `refit()` holds the dispersion
parameter theta fixed at its original estimate, so the interval is
conditional on the estimated theta.

A bootstrap draw whose *null-model* between-stratum variance lands on
the zero boundary has no defined PCV (the denominator is zero); such
draws are excluded, so the percentile interval is *conditional on
estimating a positive null variance*. Whenever any draws hit the
boundary the function warns, reports the count as `n_boot_boundary` on
the result, and [`print()`](https://rdrr.io/r/base/print.html) repeats
the caveat – a sizeable boundary share signals weak between-stratum
variation, and the PCV itself is then fragile.

Bootstrap draws whose refit optimiser reports non-convergence
(`optinfo$conv$opt != 0`) are *retained* in the interval – lme4's
post-hoc relative-gradient flag is a frequent false positive on
simulated refits, and the optimiser's own return code fires on well
under 1% of refits in practice – but the count is reported as
`n_boot_nonconverged` and [`print()`](https://rdrr.io/r/base/print.html)
discloses it. As a documented ceiling, when more than half the
contributing draws fail to converge the interval is still returned but
flagged `interval_reliable = FALSE` with an escalated warning; treat
such an interval as indicative only and check for singular or failing
fits. This is distinct from the boundary exclusion above: non-converged
draws still carry a defined PCV.

The bootstrap is available for the `lme4` engine only. For the other
engines the PCV is a *point estimate*: a brms fit's posterior credible
interval (reported by
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md))
covers a single fit's VPC/ICC, not the PCV, which compares two
separately fitted models – no posterior interval for the PCV itself is
computed – and a design-based (wemix) interval would require replicate
weights.

## References

Bauer, D. J. (2009). A note on comparing the estimates of models for
cluster-correlated or longitudinal data with binary or ordinal outcomes.
*Psychometrika*, 74(1), 97-105.

Mood, C. (2010). Logistic regression: why we cannot do what we think we
can do, and what we can do about it. *European Sociological Review*,
26(1), 67-82.

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
#> [1] -0.1626748

# Calculate PCV with bootstrap CI
# pcv_boot <- calculate_pcv(model1, model2, bootstrap = TRUE, n_boot = 500)
# print(pcv_boot)
# }
```
