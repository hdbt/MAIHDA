# Plot a MAIHDA Analysis

Dispatches each `type` to the model it is valid on. The VPC and
shrinkage views (`"vpc"`, `"obs_vs_shrunken"`, `"predicted"`) use the
**null** model. The additive-vs-intersectional views
(`"risk_vs_effect"`, `"effect_decomp"`, `"ternary"`,
`"prediction_deviation"`) use the **adjusted** model, whose fixed
effects carry the dimensions' additive part so the stratum random effect
is the pure interaction; with fewer than two dimensions (no adjusted
model) they fall back to the null model. Group types (`"group_vpc"`,
`"group_components"`, `"group_between_variance"`, `"group_pcv"`) use the
group comparison when
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) was called
with a `group`.

## Usage

``` r
# S3 method for class 'maihda_analysis'
plot(x, type = "all", ...)
```

## Arguments

- x:

  A `maihda_analysis` object from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md).

- type:

  One of the model types ("all", "vpc", "obs_vs_shrunken", "predicted",
  "risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation"),
  the contextual type ("context_vpc", a stratum-vs-context variance bar;
  requires `maihda(context = )`), or a group type ("group_vpc",
  "group_components", "group_between_variance", "group_pcv"). Default
  "all".

- ...:

  Additional arguments passed to the underlying plot method.

## Value

A ggplot2 object, or (for `type = "all"`) an invisible list of them.
