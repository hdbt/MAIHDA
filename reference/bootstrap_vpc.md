# Bootstrap VPC/ICC

Internal function to compute bootstrap confidence intervals for VPC.

## Usage

``` r
bootstrap_vpc(model, data, formula, n_boot, conf_level)
```

## Arguments

- model:

  An lme4 model object

- data:

  The data used to fit the model

- formula:

  The model formula

- n_boot:

  Number of bootstrap samples

- conf_level:

  Confidence level

## Value

A vector with lower and upper confidence bounds
