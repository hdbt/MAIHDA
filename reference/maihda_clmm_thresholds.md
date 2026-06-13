# Threshold (cut-point) estimates of a cumulative (clmm) MAIHDA fit

The thresholds \\\alpha_k\\ take the place of the intercept in a
cumulative model: \\P(Y \le k) = g^{-1}(\alpha_k - \eta)\\. Standard
errors come from the Hessian-based
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (hence `Hess = TRUE` at
fit time) and degrade to `NA` when unavailable.

## Usage

``` r
maihda_clmm_thresholds(object)
```

## Arguments

- object:

  A `maihda_model` with engine `"ordinal"`.

## Value

A data frame with `term`, `estimate`, `se`.
