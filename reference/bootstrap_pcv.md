# Bootstrap PCV

Internal function to compute bootstrap confidence intervals for PCV.

## Usage

``` r
bootstrap_pcv(model1, model2, n_boot, conf_level)
```

## Arguments

- model1:

  First maihda_model object

- model2:

  Second maihda_model object

- n_boot:

  Number of bootstrap samples

- conf_level:

  Confidence level

## Value

A vector with lower and upper confidence bounds
