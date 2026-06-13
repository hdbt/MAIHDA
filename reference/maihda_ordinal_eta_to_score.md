# Latent location to expected category score

Convenience composition of
[`maihda_ordinal_category_probs`](https://hdbt.github.io/MAIHDA/reference/maihda_ordinal_category_probs.md)
and
[`maihda_ordinal_expected_score`](https://hdbt.github.io/MAIHDA/reference/maihda_ordinal_expected_score.md).

## Usage

``` r
maihda_ordinal_eta_to_score(eta, thresholds, link = "logit")
```

## Arguments

- eta:

  Numeric vector of latent locations.

- thresholds:

  Numeric vector of increasing thresholds.

- link:

  `"logit"` or `"probit"`.

## Value

A numeric vector of expected category scores.
