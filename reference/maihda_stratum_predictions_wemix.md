# Per-stratum predictions for a wemix fit

wemix counterpart of `maihda_stratum_predictions_lme4()`: per-stratum
means of the fixed-part prediction plus the stratum effect, aggregated
with the SAMPLING weights so the stratum-level summaries are
design-weighted (population-representative under the weights), unlike
the lme4 prior-weight aggregation.

## Usage

``` r
maihda_stratum_predictions_wemix(
  object,
  summary_obj,
  scale = c("response", "link")
)
```

## Arguments

- object:

  A `maihda_model` with engine `"wemix"`.

- summary_obj:

  Its `maihda_summary` (for the stratum estimates).

- scale:

  "response" or "link".

## Value

A data frame as from `maihda_weighted_stratum_aggregate()`.
