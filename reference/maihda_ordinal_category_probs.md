# Category probabilities of a cumulative model

Pure function (no fit object): given latent locations `eta`, ordered
thresholds `alpha` and the link, returns the category probability matrix
via \\P(Y \le k) = g^{-1}(\alpha_k - \eta)\\ and differencing. Rows are
observations, columns categories `1..K` (`K = length(alpha) + 1`).

## Usage

``` r
maihda_ordinal_category_probs(eta, thresholds, link = "logit")
```

## Arguments

- eta:

  Numeric vector of latent locations.

- thresholds:

  Numeric vector of increasing thresholds \\\alpha_k\\.

- link:

  `"logit"` or `"probit"`.

## Value

A numeric matrix with `length(eta)` rows that sum to 1.
