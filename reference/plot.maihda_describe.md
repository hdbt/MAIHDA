# Plot a MAIHDA sample description

Descriptive plots for a
[`maihda_describe`](https://hdbt.github.io/MAIHDA/reference/maihda_describe.md)
object: the stratum-size distribution (the sparsity of the strata
space), the family-aware outcome distribution, and per-variable
missingness.

## Usage

``` r
# S3 method for class 'maihda_describe'
plot(x, type = c("stratum_size", "outcome", "missingness"), ...)
```

## Arguments

- x:

  A `maihda_describe` object.

- type:

  One of `"stratum_size"` (default; histogram of stratum sizes with the
  small-stratum threshold marked, empty strata included as zeros when
  they were enumerated), `"outcome"` (histogram of a continuous/count
  outcome, or category bars for a binomial/ordinal one), or
  `"missingness"` (percent missing per variable).

- ...:

  Additional arguments (not used).

## Value

A ggplot2 object, extendable with the usual `+` grammar.

## Examples

``` r
desc <- maihda_describe(BMI ~ Age + (1 | Gender:Race:Education),
                        data = maihda_health_data)
plot(desc, type = "stratum_size")

plot(desc, type = "outcome")

plot(desc, type = "missingness")
```
