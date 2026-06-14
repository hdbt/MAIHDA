# Time-varying VPC trajectory plot (longitudinal MAIHDA)

The between-stratum share of variance (VPC/ICC) as a function of time,
with a confidence/credible ribbon when available. The headline
reference-time VPC is marked. For a longitudinal MAIHDA the VPC is not a
single number – the between-stratum variance is a random intercept +
slope on time – so this curve replaces the cross-sectional VPC bar.

## Usage

``` r
plot_vpc_trajectory(summary_obj)
```

## Arguments

- summary_obj:

  A `maihda_summary` from a longitudinal model.

## Value

A ggplot2 object.
