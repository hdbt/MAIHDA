# Stratum vs. Context Variance Plot (contextual cross-classified MAIHDA)

One bar per variance component – the between-stratum (intersectional)
variance, each context's variance, any other random effects, and the
residual – on the variance scale, with each component's share of the
total printed above its bar. Complements
[`plot_vpc()`](https://hdbt.github.io/MAIHDA/reference/plot_vpc.md)'s
stacked proportion bar by showing the *magnitudes* the shares are
computed from.

## Usage

``` r
plot_context_vpc(summary_obj)
```

## Arguments

- summary_obj:

  A `maihda_summary` from a contextual cross-classified fit
  (`fit_maihda(context = )`).

## Value

A ggplot2 object.
