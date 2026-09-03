# Fixed-effect table with Wald statistics

Internal helper that assembles the `fixed_effects` slot of a
`maihda_summary` for the frequentist engines: the Wald statistic
(`estimate / se`), its two-sided p-value and the matching Wald interval
at `conf_level`. With `df` the reference is a \\t\\ on those degrees of
freedom, otherwise the normal.

## Usage

``` r
maihda_fixed_effects_table(term, estimate, se, conf_level = 0.95, df = NULL)
```

## Arguments

- term:

  Character vector of coefficient names.

- estimate:

  Numeric vector of point estimates.

- se:

  Numeric vector of standard errors (`NULL` for none).

- conf_level:

  Interval level.

- df:

  Denominator degrees of freedom, named by `term` or in `term` order;
  `NULL` (default) for the normal approximation.

## Value

A data frame with `term`, `estimate`, `se`, `statistic`, `df`,
`p_value`, `lower` and `upper`. `df` is `NA` on the normal path.
