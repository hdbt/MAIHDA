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
between-stratum variance is unavailable.

## References

Larsen, K., & Merlo, J. (2005). Appropriate assessment of neighborhood
effects on individual health: integrating random and fixed effects in
multilevel logistic regression. *American Journal of Epidemiology*,
161(1), 81-88.

## See also

[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
