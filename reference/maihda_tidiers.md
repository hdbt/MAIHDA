# Tidy a MAIHDA summary, model, or analysis

[`tidy`](https://generics.r-lib.org/reference/tidy.html) methods that
return the MAIHDA estimates as a tidy `tibble`, ready for downstream
tables (`gt`, `flextable`) and `ggplot2`. They read the slots that
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
already computes and add no new statistics.

## Usage

``` r
# S3 method for class 'maihda_summary'
tidy(x, component = c("strata", "variance", "fixed"), ...)

# S3 method for class 'maihda_model'
tidy(x, component = c("strata", "variance", "fixed"), ...)

# S3 method for class 'maihda_analysis'
tidy(
  x,
  component = c("strata", "variance", "fixed"),
  which = c("null", "adjusted"),
  ...
)
```

## Arguments

- x:

  A `maihda_summary` (from
  [`summary`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)),
  a `maihda_model` (from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)),
  or a `maihda_analysis` (from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)).

- component:

  Which estimates to return:

  `"strata"`

  :   (default) the stratum (intersection) random-effect estimates – one
      row per stratum, with `estimate`, `std.error` and
      `conf.low`/`conf.high`, plus the human-readable intersectional
      `label` when available.

  `"variance"`

  :   the variance-components table (between-stratum, any other random
      effects, residual, and total) with each component's variance, SD
      and proportion.

  `"fixed"`

  :   the fixed-effect estimates in broom's
      `term`/`estimate`/`std.error`/`statistic`/
      `p.value`/`conf.low`/`conf.high` shape.

- ...:

  Unused, for S3 consistency.

- which:

  For a `maihda_analysis`, whether to tidy the `"null"` (default) or
  `"adjusted"` model's summary.

## Value

A `tibble`. For `component = "strata"`: columns `stratum`, `label`,
`estimate`, `std.error`, `conf.low`, `conf.high`. For `"variance"`:
`component`, `variance`, `sd`, `proportion`. For `"fixed"`: `term`,
`estimate`, `std.error`, `statistic`, `df`, `p.value`, `conf.low`,
`conf.high`.

## Fixed-effect statistics

For the lme4, WeMix and ordinal engines `statistic` is the Wald
statistic (`estimate / std.error`), `p.value` its two-sided p-value,
`conf.low`/`conf.high` the matching Wald interval at the summary's
`conf_level`, and `df` the denominator degrees of freedom they were
built on.

A Gaussian lme4 fit reports containment (between-within) `df` and a
\\t\\; every other engine leaves `df` `NA` and uses a z. Under
`summary(df_method = "bootstrap")` `p.value` and the interval come from
the parametric bootstrap and `df` is `NA`. See
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md).
For brms the estimate is the posterior mean with its `Est.Error` and
credible interval, and `statistic`/`df`/`p.value` are `NA`. Standard
errors are `NA` where they are undefined (a boundary fit).

## See also

[`glance.maihda_analysis`](https://hdbt.github.io/MAIHDA/reference/maihda_glance.md)
for the one-row model summary.

## Examples

``` r
data("maihda_health_data")
m <- fit_maihda(BMI ~ Age + (1 | Gender:Race:Education), data = maihda_health_data)
tidy(m)                       # stratum estimates
#> # A tibble: 50 × 6
#>    stratum label                           estimate std.error conf.low conf.high
#>    <chr>   <chr>                              <dbl>     <dbl>    <dbl>     <dbl>
#>  1 1       male × Hispanic × Some College  -0.245       1.04   -2.29      1.80  
#>  2 2       male × Black × College Grad      0.772       1.11   -1.41      2.96  
#>  3 3       female × White × College Grad   -1.82        0.352  -2.51     -1.13  
#>  4 4       male × Hispanic × 8th Grade      0.921       1.32   -1.66      3.51  
#>  5 5       female × Mexican × 8th Grade     1.82        0.951  -0.0409    3.69  
#>  6 6       male × White × College Grad     -0.743       0.357  -1.44     -0.0431
#>  7 7       female × White × High School     0.00494     0.419  -0.816     0.826 
#>  8 8       male × White × Some College      1.02        0.356   0.326     1.72  
#>  9 9       female × White × 9 - 11th Grade  0.685       0.604  -0.498     1.87  
#> 10 10      female × Hispanic × High School  0.287       1.07   -1.81      2.38  
#> # ℹ 40 more rows
tidy(m, component = "variance")
#> # A tibble: 3 × 4
#>   component                 variance    sd proportion
#>   <chr>                        <dbl> <dbl>      <dbl>
#> 1 Between-stratum (random)      2.90  1.70     0.0627
#> 2 Within-stratum (residual)    43.4   6.59     0.937 
#> 3 Total                        46.3   6.80     1     
tidy(m, component = "fixed")
#> # A tibble: 2 × 8
#>   term        estimate std.error statistic    df  p.value  conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 (Intercept)  28.2      0.453       62.3     49 2.56e-48 27.3        29.1   
#> 2 Age           0.0150   0.00736      2.04  2949 4.11e- 2  0.000610    0.0295
```
