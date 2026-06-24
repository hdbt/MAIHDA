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
stepwise_pcv(
  data,
  outcome,
  vars,
  engine = "lme4",
  family = "gaussian",
  sampling_weights = NULL
)
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

  Modeling engine ("lme4", "brms", "wemix", or "ordinal"). Default is
  "lme4"; switches to "wemix" automatically when `sampling_weights` is
  supplied, and to "ordinal" for an ordinal family or ordered-factor
  outcome.

- family:

  Error distribution and link function. Default is "gaussian".

- sampling_weights:

  Optional name of a sampling-weight column for design-weighted stepwise
  fits; see
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  The weight column joins the complete-case filter so every step uses
  the same analytic sample.

## Value

A data.frame (class `maihda_stepwise`) showing the sequential models,
the between-stratum variance at each step, and both the step-specific
and total PCV. For a **binary** (binomial/Bernoulli) outcome it also
carries the discriminatory-accuracy trajectory: `AUC` (the C-statistic
of each step's model – step 0 is the strata-only discriminatory
accuracy), `Step_AUC` and `Total_AUC` (the *absolute* change in AUC,
delta-AUC, versus the previous step and versus the null), and `MOR` (the
Median Odds Ratio, logit link only). These columns are absent for
non-binary outcomes.

## Details

All models are fit on the complete cases for \`outcome\`, \`stratum\`,
and all variables in \`vars\` so that each sequential variance
comparison uses the same analytic sample.

For a binary outcome the table additionally tracks discriminatory
accuracy (Merlo et al. 2016): `AUC` is each model's C-statistic and
`Step_AUC` / `Total_AUC` are its *absolute* change (delta-AUC), in
contrast to the *proportional* `Step_PCV` / `Total_PCV`. The `MOR` is
reported for the logit link (`NA` otherwise) and is a monotone transform
of the between-stratum variance already in `Variance`. For a
design-weighted fit (`sampling_weights`) the AUC is the design-weighted
(population) C-statistic. Reuses
[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
on each step's fitted model, so no additional models are fit. Note that
adding a *stratum-defining* dimension (one already encoded by the
strata) typically leaves the AUC essentially unchanged: it re-partitions
the between-stratum variance (so the PCV and MOR move) but not the
per-stratum predicted ranking the rank-based AUC depends on. The AUC
trajectory is therefore most informative for individual-level covariates
that vary *within* strata.

## References

Merlo, J., Wagner, P., Ghith, N., & Leckie, G. (2016). An original
stepwise multilevel logistic regression analysis of discriminatory
accuracy: the case of neighbourhoods and health. *PLOS ONE*, 11(4),
e0153778.

## Examples

``` r
# \donttest{
strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
stepwise_pcv(strata_result$data, "health_outcome", c("gender", "race", "age"))
#>  Step      Model        Added_Variable  Variance Step_PCV Total_PCV
#>     0 Null Model None (Intercept only) 2.324e+01  0.00000   0.00000
#>     1    Model 1                gender 2.290e+01  0.01457   0.01457
#>     2    Model 2                  race 7.564e-14  1.00000   1.00000
#>     3    Model 3                   age 0.000e+00  1.00000   1.00000
# }
```
