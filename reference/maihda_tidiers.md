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
`estimate`, `std.error`, `statistic`, `p.value`, `conf.low`,
`conf.high`.

## Fixed-effect statistics

For the lme4, WeMix and ordinal engines `statistic` is the Wald \\z\\
(`estimate / std.error`), `p.value` its two-sided normal-approximation
p-value, and `conf.low`/`conf.high` the matching Wald interval at the
`conf_level` the summary was computed at (95% by default). No
denominator degrees of freedom are estimated, so these are z-based
rather than Satterthwaite or Kenward-Roger.

Which terms that approximation is safe for depends on how they vary. A
covariate that varies *within* a stratum has an effective sample size of
the number of observations, and z is essentially exact for it. The
intercept and the dimension main effects of an *adjusted* MAIHDA model
are constant within a stratum by construction, so their effective sample
size is the **number of strata**, however large \\n\\ is: with few
strata their intervals are too narrow and their p-values
anticonservative. As a rough guide, a Kenward-Roger check gives those
terms about 5% wider intervals at 50 strata but roughly 40% wider at 10.
When that matters, apply lmerTest or pbkrtest to the underlying fit
(`x$model`, an `lmerMod`), or use the brms engine.

For brms the estimate is the posterior mean with its `Est.Error` and
credible interval, and `statistic`/`p.value` are `NA`. Standard errors
are `NA` where they are undefined (a boundary fit).

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
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic p.value  conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>     <dbl>     <dbl>
#> 1 (Intercept)  28.2      0.453       62.3   0      27.3        29.1   
#> 2 Age           0.0150   0.00736      2.04  0.0410  0.000616    0.0295
```
