# Variance components of a wemix MAIHDA fit

Reads the between-stratum variance (and, for a linear model, the
residual variance) from the `WeMixResults` variance table. For a
binomial-logit model the level-1 variance is the usual latent-scale
\\\pi^2/3\\, matching the lme4/brms summaries.

## Usage

``` r
maihda_wemix_variances(object)
```

## Arguments

- object:

  A `maihda_model` with engine `"wemix"`.

## Value

A list with `stratum` and `residual` variances.
