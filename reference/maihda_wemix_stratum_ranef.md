# Stratum random-effect table for a wemix fit

Mirrors `maihda_stratum_ranef_lme4()`: one row per stratum with the
random-effect estimate (conditional mode), a conditional standard error,
and a 95% interval. WeMix reports no conditional variances, so the SE is
computed analytically from the weighted pseudo-likelihood: for a
Gaussian model the conditional precision of \\u_j\\ is \\1/\tau^2 +
\sum_j w\_{ij}/\sigma^2\\ (the design-weighted analogue of lme4's
`condVar`, to which it reduces at unit weights), and for a
binomial-logit model the Laplace curvature at the conditional mode,
\\1/\tau^2 + \sum_j w\_{ij}\\\hat p\_{ij}(1-\hat p\_{ij})\\. These are
model-based approximations, not design-based (replicate-weight)
uncertainty.

## Usage

``` r
maihda_wemix_stratum_ranef(object)
```

## Arguments

- object:

  A `maihda_model` with engine `"wemix"`.

## Value

A data frame with `stratum`, `stratum_id`, `random_effect`, `se`,
`lower_95`, `upper_95`.
