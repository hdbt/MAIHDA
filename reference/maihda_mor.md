# Median Odds Ratio (MOR) for a logistic MAIHDA model

The Median Odds Ratio translates the between-stratum variance of a
logistic MAIHDA model onto the odds-ratio scale: the median relative
change in the odds of the outcome when comparing two individuals from
randomly chosen strata (higher- vs lower-risk).
`MOR = exp(sqrt(2 * V_A) * qnorm(0.75))`, where `V_A` is the
between-stratum (latent, logit-scale) variance. An MOR of 1 indicates no
between-stratum heterogeneity. The MOR is defined only for the **logit**
link (it is the median *odds* ratio); a non-logit binomial fit such as
`probit` is rejected, because its latent variance is on a different
scale and the `exp(...)` above would not be an odds ratio.

For a **cumulative-logit** (ordinal) MAIHDA model the same formula
applies to the latent logit-scale between-stratum variance and is the
*median cumulative odds ratio*: the median relative change in the odds
of being at or below any given outcome category between two randomly
chosen strata (under the model's proportional-odds assumption it is the
same for every category split).

**Scope.** `V_A` is the *between-stratum* variance. Contextual
(`context = `) and other non-stratum random effects are never included:
the MOR compares two individuals from randomly chosen strata *within the
same context*.

**Crossed-dimensions fits use a different calculation.** The closed form
above assumes the two strata's random effects are *independent*, which
is what makes their difference \\N(0, 2 V_A)\\. In a crossed-dimensions
fit (from `maihda(decomposition = "crossed-dimensions")`) a stratum's
effect is the sum of its dimension effects plus the intersection effect,
so two strata sharing a dimension – say two strata that are both
"female" – share that dimension's random effect and are **correlated**.
Their difference is then a *mixture* of normals, one component per
pattern of shared dimensions, and applying the closed form to the summed
variance overstates the MOR (substantially so when the variance sits
mainly in the additive dimensions, and not at all when it sits entirely
in the interaction).

The MOR reported for such a fit is therefore computed from that mixture,
under an explicit sampling scheme: **two distinct strata drawn uniformly
at random from the strata present in the analytic sample**. Writing
\\\tau^2_d\\ for dimension \\d\\'s variance and \\\tau^2_I\\ for the
intersection variance, a pair differing on the dimension set \\D^\*\\
has difference variance \\v = 2(\tau^2_I + \sum\_{d \in D^\*}
\tau^2_d)\\, and the MOR is `exp(x)` for the `x` solving \\\sum\_{pairs}
(2\Phi(x/\sqrt{v}) - 1) / n\_{pairs} = 0.5\\. For a canonical
single-stratum fit the two calculations coincide, and that closed form
is used.

## Usage

``` r
maihda_mor(model)
```

## Arguments

- model:

  A `maihda_model` from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  fitted with a `binomial` (lme4), `bernoulli` (brms), or `cumulative`
  (ordinal) family and a **logit** link.

## Value

A single number (the MOR, \\\ge 1\\), or `NA_real_` if the
between-stratum variance is unavailable – which for a crossed-dimensions
fit also covers the case where the stratum grid needed for the mixture
cannot be resolved (fewer than two strata, absent dimension columns, or
more than 12 dimensions).

## References

Larsen, K., & Merlo, J. (2005). Appropriate assessment of neighborhood
effects on individual health: integrating random and fixed effects in
multilevel logistic regression. *American Journal of Epidemiology*,
161(1), 81-88.

## See also

[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
