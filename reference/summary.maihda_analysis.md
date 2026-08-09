# Summarize a MAIHDA Analysis

Returns the variance summary (VPC/ICC, variance components, stratum
estimates) of one of the analysis's models. The per-group comparison,
when present, is attached as the `"groups"` attribute.

## Usage

``` r
# S3 method for class 'maihda_analysis'
summary(object, which = c("null", "adjusted"), ...)
```

## Arguments

- object:

  A `maihda_analysis` object from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

- which:

  Which model to summarize: `"null"` (default) or `"adjusted"`. A
  crossed-dimensions analysis has a single model and accepts only
  `"null"`.

- ...:

  Additional arguments (not used).

## Value

The `maihda_summary` for the requested model.

## Which model

A two-model
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) analysis
holds a **null** model (strata as the only predictor) and an
**adjusted** model (the dimensions' additive main effects added as fixed
effects), and this returns the null model's summary by default – the one
the headline VPC/ICC comes from. The strata dimensions are the
random-effect grouping there, not fixed effects, so the printed
fixed-effect table shows only the intercept and any covariates. Use
`which = "adjusted"` for the model whose fixed effects carry the
dimensions' additive main effects (the same object as
`object$summary_adjusted`), matching
[`tidy`](https://hdbt.github.io/MAIHDA/reference/maihda_tidiers.md)'s
`which` argument.
