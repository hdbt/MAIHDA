# Tag a summary with its role in a two-model analysis

Internal helper. The null and adjusted summaries of a
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) analysis
are indistinguishable once pulled out of the analysis object, which
makes a printed summary easy to misread – the null model's fixed effects
are the intercept and covariates only, because the strata dimensions are
its random-effect grouping rather than fixed effects. Stamping the role
lets
[`print.maihda_summary`](https://hdbt.github.io/MAIHDA/reference/print.maihda_summary.md)
say which model it is showing.

## Usage

``` r
maihda_tag_role(s, role)
```

## Arguments

- s:

  A `maihda_summary`, or `NULL`.

- role:

  `"null"` or `"adjusted"`.

## Value

`s` with a `"maihda_role"` attribute (`NULL` in, `NULL` out).
