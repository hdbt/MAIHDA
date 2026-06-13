# Expected category score from a probability matrix

The response-scale summary of a cumulative model used throughout the
package (the plot layer's "Average Expected Category Score"): \\\sum_k
k\\ p_k\\, with categories scored 1..K in order.

## Usage

``` r
maihda_ordinal_expected_score(probs)
```

## Arguments

- probs:

  A category-probability matrix (rows = observations).

## Value

A numeric vector of expected scores in \\\[1, K\]\\.
