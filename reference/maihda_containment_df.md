# Containment degrees of freedom for a Gaussian mixed-model fixed effect

Internal helper. Returns between-within (containment) denominator
degrees of freedom for each fixed-effect column of a Gaussian `lmerMod`:
a column absorbed by a random-effect term is tested against that term's
units minus the columns it absorbs, a column absorbed by none against
\\n\\ minus the random-effect levels. A column counts as absorbed when,
within every level of the grouping factor, it is a linear combination of
that term's random-effect covariates – constant within the group for a
random intercept, any group-specific linear function of time for a
random slope on time.

## Usage

``` r
maihda_containment_df(model)
```

## Arguments

- model:

  A fitted model. `NULL` for anything that is not a Gaussian `lmerMod`.

## Value

A named numeric vector of degrees of freedom, or `NULL`.

## Details

Where several terms absorb a column the smallest applies, taken over the
terms that materially contribute to that coefficient's variance
(`maihda_variance_shares`); a term whose variance component is estimated
at essentially zero cannot limit precision. That filter can narrow the
set but never empty it, so a singular fit keeps its stratum degrees of
freedom.

The per-source degrees of freedom are deliberately not
Satterthwaite-combined: in the balanced one-way case the stratum and
residual shares of a stratum-level contrast come from the same mean
square, so combining them as independent chi-squares inflates the result
several-fold. Cross-checked against pbkrtest in
`tests/testthat/test-fixed-effect-df.R`; pbkrtest is declared nowhere
and those checks skip without it.
