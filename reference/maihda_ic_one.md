# Information criteria for a single MAIHDA model

Internal worker for
[`maihda_ic`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md):
returns a one-row data frame of the fit criteria for one `maihda_model`,
dispatched on the fitted object's class (mirroring
`maihda_fit_diagnostics`).

## Usage

``` r
maihda_ic_one(model, ml = FALSE)
```

## Arguments

- model:

  A `maihda_model`.

- ml:

  Logical; for a REML `lmer` fit, refit with ML via
  [`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html) before reading
  AIC/BIC (used when comparing models that may differ in fixed effects).

## Value

A one-row data frame with `n`, `estimator`, `df`, `logLik`, `AIC`,
`BIC`, `WAIC`, `LOOIC` (NA where not applicable to the engine).
