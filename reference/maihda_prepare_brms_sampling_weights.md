# Prepare data and formula for a sampling-weighted brms fit

Drops rows with missing or non-positive sampling weights (with a
warning), normalizes the remaining weights to mean 1 – likelihood
weights scale the effective sample size, so unnormalized expansion
weights (summing to the population) would massively overstate the
information in the data – and rewrites the formula with a
[`weights()`](https://rdrr.io/r/stats/weights.html) addition term.

## Usage

``` r
maihda_prepare_brms_sampling_weights(data, formula, sampling_weights)
```

## Arguments

- data:

  The model data.

- formula:

  The model formula.

- sampling_weights:

  Name of the sampling-weight column.

## Value

A list with `data` (weights column `.maihda_sw` added), `formula`
(rewritten), and `keep` (a logical mask over the input rows marking
those retained, so the caller can re-slice any row-aligned forwarded
arguments to the same rows).
