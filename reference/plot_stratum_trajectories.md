# Stratum mean-trajectory plot (longitudinal MAIHDA)

One predicted line per stratum over time – the fixed-part trajectory
plus each stratum's random intercept and slope (BLUPs) – the
longitudinal analogue of the predicted-strata caterpillar. Shows how the
intersectional groups fan out (or converge) over time.

## Usage

``` r
plot_stratum_trajectories(
  object,
  summary_obj,
  n_strata = 50,
  select = c("order", "deviation")
)
```

## Arguments

- object:

  A longitudinal `maihda_model`.

- summary_obj:

  Its `maihda_summary`.

- n_strata:

  Maximum number of strata to draw; the rest are noted in the caption.

- select:

  When the cap drops strata, which to keep: `"order"` (default) the
  first n_strata in stratum order, or `"deviation"` the n_strata whose
  trajectory swings furthest from the population curve (largest peak
  `|random deviation|` over the time grid, either direction).

## Value

A ggplot2 object.
