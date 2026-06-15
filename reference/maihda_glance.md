# Glance at a MAIHDA model or analysis

[`glance`](https://generics.r-lib.org/reference/glance.html) methods
that return the MAIHDA headline as a one-row `tibble`: the variance
partition coefficient (VPC/ICC), and – for a `maihda_analysis` – the
proportional change in variance (PCV), plus the additive/interaction
shares and the discriminatory accuracy (AUC, MOR) for a binomial
outcome. The layout is uniform across the lme4, brms, WeMix and ordinal
engines. No other package emits this row from the underlying fit,
because PCV needs the null+adjusted pair that only a `maihda_analysis`
holds.

## Usage

``` r
# S3 method for class 'maihda_summary'
glance(x, ...)

# S3 method for class 'maihda_model'
glance(x, ...)

# S3 method for class 'maihda_analysis'
glance(x, ...)
```

## Arguments

- x:

  A `maihda_summary`, `maihda_model`, or `maihda_analysis`.

- ...:

  Unused, for S3 consistency.

## Value

A one-row `tibble`. `glance.maihda_analysis` adds `pcv` (with
`pcv.conf.low`/`pcv.conf.high` when bootstrapped or from a brms
posterior), the adjusted-model `auc.adjusted`, `nobs`, `family` and
`mode` to the columns produced for a single summary.

## See also

[`maihda_tidiers`](https://hdbt.github.io/MAIHDA/reference/maihda_tidiers.md)
for the per-estimate
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) methods.

## Examples

``` r
data("maihda_health_data")
a <- maihda(BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
            data = maihda_health_data)
glance(a)
#> # A tibble: 1 × 16
#>      vpc vpc.conf.low vpc.conf.high   pcv pcv.conf.low pcv.conf.high
#>    <dbl>        <dbl>         <dbl> <dbl>        <dbl>         <dbl>
#> 1 0.0627           NA            NA 0.766           NA            NA
#> # ℹ 10 more variables: additive.share <dbl>, interaction.share <dbl>,
#> #   auc <dbl>, auc.adjusted <dbl>, mor <dbl>, n_strata <int>, nobs <int>,
#> #   engine <chr>, family <chr>, mode <chr>
```
