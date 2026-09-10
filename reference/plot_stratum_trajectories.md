# Stratum mean-trajectory plot (longitudinal MAIHDA)

One predicted line per stratum over time – that stratum's fixed-part
trajectory plus its own random intercept and slope (BLUPs) – the
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
  first n_strata in stratum order, or `"deviation"` the n_strata with
  the largest peak `|random deviation|` over the time grid, either
  direction. For a null model that is the departure from the population
  curve; for an adjusted model it is the departure from the stratum's
  own additive prediction, i.e. the strongest intersectional
  interactions.

## Value

A ggplot2 object.

## Details

The fixed part varies with the stratum drawn: the stratum-defining
dimensions take that stratum's values, so an adjusted growth model's
dimension main effects and `dim:time` interactions are included, while
every *other* covariate is held at a shared reference profile (mean for
a numeric, modal level for a factor). The lines are therefore predicted
stratum trajectories net of covariate composition, excluding the
individual-level random effects. For a null growth model the fixed part
carries no dimension terms and every line shares the one population
trajectory.
