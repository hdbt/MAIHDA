# Threshold (cut-point) estimates of a cumulative (clmm) MAIHDA fit

The thresholds \\\alpha_k\\ take the place of the intercept in a
cumulative model: \\P(Y \le k) = g^{-1}(\alpha_k - \eta)\\. Reported as
the \\K-1\\ expanded cut points
([`maihda_clmm_cutpoints`](https://hdbt.github.io/MAIHDA/reference/maihda_clmm_cutpoints.md)),
so the table means the same thing under every threshold structure rather
than echoing the free parameters of a constrained fit. Standard errors
come from the Hessian-based
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (hence `Hess = TRUE` at
fit time) and degrade to `NA` when unavailable; under a structured
threshold they are the delta-method SEs of the cut points, \\J V J'\\
for the fit's threshold Jacobian `$tJac`. That Jacobian is the identity
under the default `threshold = "flexible"`, where the transform
reproduces the direct [`vcov()`](https://rdrr.io/r/stats/vcov.html) SEs
exactly.

## Usage

``` r
maihda_clmm_thresholds(object)
```

## Arguments

- object:

  A `maihda_model` with engine `"ordinal"`.

## Value

A data frame with `term`, `estimate`, `se`.
