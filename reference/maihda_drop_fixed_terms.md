# Drop fixed terms from a formula without dropping the grand mean

The single place a MAIHDA derives a reduced model by removing the
stratum dimensions' fixed main effects: the two-model null, the
crossed-dimensions base formula, the longitudinal null and the per-group
formulas all come through here.

## Usage

``` r
maihda_drop_fixed_terms(formula, terms, data, notify = FALSE, fn = "maihda")
```

## Arguments

- formula:

  The source formula (the adjusted / supplied model).

- terms:

  Character vector of fixed terms to remove (raw, unquoted names).

- data:

  Data frame the fixed designs are built on.

- notify:

  Logical; message once when the grand mean is actually restored.
  `FALSE` at the per-group call sites, which repeat a reduction the
  overall decomposition has already reported.

- fn:

  Name of the calling function, used in the message prefix.

## Value

The reduced formula, in the same environment as `formula`.

## Details

`update(f, . ~ . - a - b)` removes the named terms, but on a
no-intercept formula it removes the grand mean with them.
`y ~ x + a + b + (1 | a:b)` and `y ~ 0 + x + a + b + (1 | a:b)` are the
SAME adjusted model – with no intercept R codes `a` as cell means, so
the two designs span one space and fit identical values – yet the plain
reduction turns the first into `y ~ x + (1 | stratum)` and the second
into `y ~ x + (1 | stratum) - 1`. The second null has no way to
represent the outcome's mean, so the stratum random intercept absorbs
it: the between-stratum variance, and with it the VPC, MOR and PCV,
become functions of the arbitrary origin of the response rather than of
the strata.

The grand mean is therefore restored whenever the reduction loses it AND
the source model carried it. Conditioning on the source keeps the
reduced model nested: a fit that genuinely has no intercept in its span
(`y ~ 0 + x + a + b` with numeric dimensions – a regression through the
origin, not an equivalent recoding) keeps its parameterization, because
adding a grand mean there would put a column in the null that the
adjusted model does not have. A reduction that still spans the intercept
some other way – a surviving factor covariate coded as cell means, or a
fixed part that collapses to `1` once the bars are stripped – is left
exactly as it was.
