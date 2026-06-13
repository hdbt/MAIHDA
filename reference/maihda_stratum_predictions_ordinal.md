# Per-stratum predictions for a cumulative (clmm) fit

Ordinal counterpart of
[`maihda_stratum_predictions_wemix()`](https://hdbt.github.io/MAIHDA/reference/maihda_stratum_predictions_wemix.md):
per-stratum aggregates of the location prediction plus the stratum
effect, on the latent (link) scale or as the expected category score
(response scale).

## Usage

``` r
maihda_stratum_predictions_ordinal(
  object,
  summary_obj,
  scale = c("response", "link")
)
```

## Arguments

- object:

  A `maihda_model` with engine `"ordinal"`.

- summary_obj:

  Its `maihda_summary` (for the stratum estimates).

- scale:

  "response" (expected category score) or "link" (latent).

## Value

A data frame as from `maihda_weighted_stratum_aggregate()`.
