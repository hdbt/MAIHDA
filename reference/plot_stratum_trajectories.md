# Stratum mean-trajectory plot (longitudinal MAIHDA)

One predicted line per stratum over time – the fixed-part trajectory
plus each stratum's random intercept and slope (BLUPs) – the
longitudinal analogue of the predicted-strata caterpillar. Shows how the
intersectional groups fan out (or converge) over time.

## Usage

``` r
plot_stratum_trajectories(object, summary_obj, n_strata = 50)
```

## Arguments

- object:

  A longitudinal `maihda_model`.

- summary_obj:

  Its `maihda_summary`.

- n_strata:

  Maximum number of strata to draw (by stratum order); the rest are
  noted in the caption.

## Value

A ggplot2 object.
