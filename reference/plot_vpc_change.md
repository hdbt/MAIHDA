# Null-vs-Adjusted Variance-Partition Change Plot

Draws the null and adjusted models' variance partitions as two stacked
bars on a shared axis and annotates the PCV, so the drop in the
between-stratum share after accounting for the dimensions' additive main
effects is visible in one figure. Backs
`plot(<maihda_analysis>, type = "vpc", model = "both")` for a
cross-sectional analysis.

## Usage

``` r
plot_vpc_change(null_summary, adjusted_summary, pcv = NULL)
```

## Arguments

- null_summary, adjusted_summary:

  The `maihda_summary` objects of the null and adjusted models of a
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis.

- pcv:

  Optional `pcv_result` (the analysis's `$pcv`); its estimate is shown
  in the subtitle. `NULL` or a non-finite estimate drops the PCV clause
  (e.g. a boundary fit where the PCV is undefined).

## Value

A ggplot2 object.
