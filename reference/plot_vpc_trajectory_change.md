# Null-vs-Adjusted VPC-over-time Change Plot (longitudinal MAIHDA)

Overlays the null and adjusted models' time-varying VPC/ICC curves on
one axis, so the reduction in the between-stratum share after adjustment
is visible over the whole time range. Backs
`plot(<maihda_analysis>, type = "vpc", model = "both")` for a
longitudinal analysis; complements
[`plot_pcv_trajectory`](https://hdbt.github.io/MAIHDA/reference/plot_pcv_trajectory.md)
(the additive-share PCV(t)).

## Usage

``` r
plot_vpc_trajectory_change(null_summary, adjusted_summary)
```

## Arguments

- null_summary, adjusted_summary:

  Longitudinal `maihda_summary` objects.

## Value

A ggplot2 object.
