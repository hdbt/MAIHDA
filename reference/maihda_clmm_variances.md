# Variance components of a cumulative (clmm) MAIHDA fit

Reads the between-stratum variance from
[`ordinal::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) and
pairs it with the latent-scale level-1 variance (\\\pi^2/3\\ for logit,
1 for probit), matching the latent treatment of binomial models in the
other engines.

## Usage

``` r
maihda_clmm_variances(object)
```

## Arguments

- object:

  A `maihda_model` with engine `"ordinal"`.

## Value

A list with `stratum` and `residual` variances.
