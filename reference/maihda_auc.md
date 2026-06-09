# Area under the ROC curve (C-statistic), rank-based

Computes the AUC / C-statistic as the Mann-Whitney U statistic: the
probability that a randomly chosen case (`y == 1`) is assigned a higher
predicted value than a randomly chosen non-case (`y == 0`), with ties
counting as one half. This needs no external package. An AUC of 0.5 is
chance; 1 is perfect separation.

## Usage

``` r
maihda_auc(prob, y)
```

## Arguments

- prob:

  Numeric vector of predicted probabilities (or any score where larger
  means more case-like).

- y:

  Observed binary outcome as 0/1 numeric or logical, the same length as
  `prob`.

## Value

A single number in `[0, 1]`, or `NA_real_` if either class is absent.

## References

Merlo, J., Wagner, P., Ghith, N., & Leckie, G. (2016). An original
stepwise multilevel logistic regression analysis of discriminatory
accuracy: the case of neighbourhoods and health. *PLOS ONE*, 11(4),
e0153778.

## Examples

``` r
maihda_auc(c(0.1, 0.4, 0.35, 0.8), c(0, 0, 1, 1))
#> [1] 0.75
```
