# Cumulative (ordinal) family marker for MAIHDA models

Specifies a cumulative (proportional-odds) model for an ordinal outcome
in [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
/ [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md), with a
choice of link: `maihda_cumulative("logit")` (the default, equivalent to
`family = "ordinal"`) or `maihda_cumulative("probit")`. It plays the
role a `stats` family object plays for the other families – there is no
cumulative family constructor in `stats`, and using
[`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
would require brms for a frequentist fit.

## Usage

``` r
maihda_cumulative(link = c("logit", "probit"))
```

## Arguments

- link:

  The cumulative link: `"logit"` (default) or `"probit"`. These are the
  links for which the latent-scale VPC is defined (level-1 variance
  \\\pi^2/3\\ and 1 respectively).

## Value

A family marker list with elements `family = "cumulative"` and `link`.

## See also

[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)

## Examples

``` r
maihda_cumulative()
#> $family
#> [1] "cumulative"
#> 
#> $link
#> [1] "logit"
#> 
maihda_cumulative("probit")
#> $family
#> [1] "cumulative"
#> 
#> $link
#> [1] "probit"
#> 
```
